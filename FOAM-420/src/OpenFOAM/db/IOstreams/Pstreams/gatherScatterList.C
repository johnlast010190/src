/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM (R) <http://www.openfoam.org/>.

    FOAMcore is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    FOAMcore is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with FOAMcore.  If not, see <http://www.gnu.org/licenses/>.

Copyright
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Description
    Gather data from all processors onto single processor according to some
    communication schedule (usually linear-to-master or tree-to-master).
    The gathered data will be a list with element procID the data from processor
    procID. Before calling every processor should insert its value into
    Values[UPstream::myProcNo(comm)].
    Note: after gather every processor only knows its own data and that of the
    processors below it. Only the 'master' of the communication schedule holds
    a fully filled List. Use scatter to distribute the data.

\*---------------------------------------------------------------------------*/

#include "mpiFunctors.H"
#include "FastSerialiser.H"
#include "db/IOstreams/Pstreams/IPstream.H"
#include "db/IOstreams/Pstreams/OPstream.H"
#include "primitives/contiguous/contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<bool All, typename T>
void performGather(
    List<T>& Values,
    const label comm
) {
    if (!UPstream::parRun() || UPstream::nProcs(comm) == 1) {
        return;
    }

    using Gather = std::conditional_t<All, mpiAllGather, mpiGather>;
    using Gatherv = std::conditional_t<All, mpiAllGatherv, mpiGatherv>;

    label nProcs = UPstream::nProcs(comm);
    int pid = UPstream::myProcNo(comm);
    bool master = UPstream::master(comm);

    FOAM_ASSERT(Values.size() == nProcs) {
        FatalErrorInFunction
            << "Size of list:" << Values.size()
            << " does not equal the number of peers in the communicator for a gather operation:"
            << nProcs
            << Foam::abort(FatalError);
    }

    bool failure = false;
    if constexpr (contiguous<T>()) {
        failure |= Gather{}(
            // Send.
            // nullptr is used as a stand-in for MPI_IN_PLACE, to avoid having to include the MPI header here
            // because that would involve un-insane-ifying the dynamic code system.
            ((master || All) ? nullptr : &Values[pid]), sizeof(T),

            // Recv
            Values.data(), sizeof(T),

            comm
        );
    } else {
        // Serialise the message into the buffer.
        OStringStream ss{IOstream::BINARY};
        fastSerialise(ss, Values[pid]);
        string buf = ss.str();
        int length = static_cast<int>(buf.size());

        // TODO: It is possible to be smarter. We use this code path right now for every type that
        //       is not contiguous (ie. can't just be memcpy'd). Buuut: we only actually need to
        //       do this relatively inefficient gather-lengths-then-gather-values thing when we
        //       don't know the lengths upfront. A type that contains a heap alloaction but for which
        //       the length is known (either from some other part of the logic, or just because it is
        //       always the same) can be handled more efficiently - just use the known lengths/displacements
        //       directly.
        List<int> lengths(nProcs);
        lengths[pid] = length;

        if constexpr (All) {
            Pstream::allGatherList(lengths, comm);
        } else {
            Pstream::gatherList(lengths, comm);
        }

        // Really must be int, not label, since the MPI library works in ints.
        // Presumably we never send things that would cause an overflow...? :D
        List<int> displacements;
        string rcvBuf;

        if (All || master) {
            displacements.resize(nProcs);

            int acc = 0;
            for (int i = 0; i < displacements.size(); i++) {
                displacements[i] = acc;
                acc += lengths[i];
            }
            rcvBuf = string{static_cast<size_t>(acc), '\0'};
        }

        // Now the master knows how much to receive from every thread, and every thread
        // knows how much to send (and has serialised it). So we can do the thing.
        failure |= Gatherv{}(
            // Send
            buf.data(), length,

            // Recv
            rcvBuf.data(), lengths.data(), displacements.data(),

            comm
        );

        // Deserialise the results.
        if (All || master) {
            for (int i = 0; i < nProcs; i++) {
                string fragment
                {
                    rcvBuf.data() + displacements[i],
                    static_cast<size_t>(lengths[i])
                };
                IStringStream iss(fragment, IOstream::BINARY);
                fastDeserialise(iss, Values[i]);
            }
        }
    }

    FOAM_ASSERT(!failure) {
        FatalErrorInFunction
            << "MPI failed to gather."
            << Foam::abort(FatalError);
    }
}

template<class T>
void Pstream::allGatherList
(
    List<T>& Values,
    label comm
)
{
    performGather<true>(Values, comm);
}

template<class T>
void Pstream::gatherList(
    List<T>& Values,
    label comm
) {
    performGather<false>(Values, comm);
}

template<class T>
void Pstream::scatterList
(
    const List<UPstream::commsStruct>& comms,
    List<T>& Values,
    const int tag,
    const label comm
)
{
    if (!UPstream::parRun() || UPstream::nProcs(comm) == 1) {
        return;
    }

    FOAM_ASSERT(Values.size() == UPstream::nProcs(comm)) {
        FatalErrorInFunction
            << "Size of list:" << Values.size()
            << " does not equal the number of processors:"
            << UPstream::nProcs(comm)
            << Foam::abort(FatalError);
    }

    // Get my communication order
    const commsStruct& myComm = comms[UPstream::myProcNo(comm)];

    // Reveive from up
    if (myComm.above() != -1)
    {
        const labelList& notBelowLeaves = myComm.allNotBelow();

        if constexpr (contiguous<T>())
        {
            List<T> receivedValues(notBelowLeaves.size());

            UIPstream::read
            (
                UPstream::commsTypes::scheduled,
                myComm.above(),
                reinterpret_cast<char*>(receivedValues.begin()),
                receivedValues.byteSize(),
                tag,
                comm
            );

            forAll(notBelowLeaves, leafI)
            {
                Values[notBelowLeaves[leafI]] = receivedValues[leafI];
            }
        }
        else
        {
            IPstream fromAbove
            (
                UPstream::commsTypes::scheduled,
                myComm.above(),
                0,
                tag,
                comm
            );

            forAll(notBelowLeaves, leafI)
            {
                label leafID = notBelowLeaves[leafI];
                fromAbove >> Values[leafID];

                if (debug)
                {
                    Pout<< " received through "
                        << myComm.above() << " data for:" << leafID
                        << " data:" << Values[leafID] << endl;
                }
            }
        }
    }

    // Send to my downstairs neighbours
    forAllReverse(myComm.below(), belowI)
    {
        label belowID = myComm.below()[belowI];
        const labelList& notBelowLeaves = comms[belowID].allNotBelow();

        if constexpr (contiguous<T>())
        {
            List<T> sendingValues(notBelowLeaves.size());

            forAll(notBelowLeaves, leafI)
            {
                sendingValues[leafI] = Values[notBelowLeaves[leafI]];
            }

            OPstream::write
            (
                UPstream::commsTypes::scheduled,
                belowID,
                reinterpret_cast<const char*>(sendingValues.begin()),
                sendingValues.byteSize(),
                tag,
                comm
            );
        }
        else
        {
            OPstream toBelow
            (
                UPstream::commsTypes::scheduled,
                belowID,
                0,
                tag,
                comm
            );

            // Send data destined for all other processors below belowID
            forAll(notBelowLeaves, leafI)
            {
                label leafID = notBelowLeaves[leafI];
                toBelow << Values[leafID];

                if (debug)
                {
                    Pout<< " sent through "
                        << belowID << " data for:" << leafID
                        << " data:" << Values[leafID] << endl;
                }
            }
        }
    }
}


template<class T>
void Pstream::scatterList(List<T>& Values, const int tag, const label comm)
{
    if (UPstream::nProcs(comm) < UPstream::nProcsSimpleSum)
    {
        scatterList(UPstream::linearCommunication(comm), Values, tag, comm);
    }
    else
    {
        scatterList(UPstream::treeCommunication(comm), Values, tag, comm);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
