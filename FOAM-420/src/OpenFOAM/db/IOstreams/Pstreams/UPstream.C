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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/Pstreams/UPstream.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "db/IOstreams/Pstreams/PstreamReduceOps.H"
#include "global/debug/debug.H"
#include "global/debug/registerSwitch.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/IOstreams.H"
#include "containers/Lists/SubList/SubList.H"
#include "PstreamGlobals.H"
#include "primitives/ints/int/int.H"
#include "global/fileOperations/collatedFileOperation/collatedFileOperation.H"

#include <mpi.h>
#include <cstring>
#include <cstdlib>
#include <csignal>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// The min value and default for MPI buffers length
constexpr int minBufLen = 20000000;

// Track if we have attached MPI buffers
static bool ourBuffers = false;

// Track if we initialized MPI
static bool ourMpi = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::UPstream::setParRun(const label nProcs, const bool haveThreads)
{
    if (nProcs == 0)
    {
        parRun_ = false;
        haveThreads_ = haveThreads;

        freeCommunicator(UPstream::worldComm);
        label comm = allocateCommunicator(-1, labelList(1, Zero), false);
        if (comm != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::worldComm:" << UPstream::worldComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = "";
        Perr.prefix() = "";
    }
    else
    {
        parRun_ = true;
        haveThreads_ = haveThreads;

        // Redo worldComm communicator (this has been created at static
        // initialisation time)
        freeCommunicator(UPstream::worldComm);
        label comm = allocateCommunicator(-1, identity(nProcs), true);
        if (comm != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "problem : comm:" << comm
                << "  UPstream::worldComm:" << UPstream::worldComm
                << Foam::exit(FatalError);
        }

        Pout.prefix() = '[' +  name(myProcNo(Pstream::worldComm)) + "] ";
        Perr.prefix() = '[' +  name(myProcNo(Pstream::worldComm)) + "] ";
    }
}


Foam::label Foam::UPstream::allocateCommunicator
(
    const label parentIndex,
    const labelList& subRanks,
    const bool doPstream
)
{
    label index;
    if (!freeComms_.empty())
    {
        index = freeComms_.pop();
    }
    else
    {
        // Extend storage
        index = parentCommunicator_.size();

        myProcNo_.append(-1);
        procIDs_.append(List<int>(0));
        parentCommunicator_.append(-1);
        linearCommunication_.append(List<commsStruct>(0));
        treeCommunication_.append(List<commsStruct>(0));
    }

    if (debug)
    {
        Pout<< "Communicators : Allocating communicator " << index << endl
            << "    parent : " << parentIndex << endl
            << "    procs  : " << subRanks << endl
            << endl;
    }

    // Initialise; overwritten by allocatePstreamCommunicator
    myProcNo_[index] = 0;

    // Convert from label to int
    procIDs_[index].setSize(subRanks.size());
    forAll(procIDs_[index], i)
    {
        procIDs_[index][i] = subRanks[i];

        // Enforce incremental order (so index is rank in next communicator)
        if (i >= 1 && subRanks[i] <= subRanks[i-1])
        {
            FatalErrorInFunction
                << "subranks not sorted : " << subRanks
                << " when allocating subcommunicator from parent "
                << parentIndex
                << Foam::abort(FatalError);
        }
    }
    parentCommunicator_[index] = parentIndex;

    // Size but do not fill structure - this is done on-the-fly
    linearCommunication_[index] = List<commsStruct>(procIDs_[index].size());
    treeCommunication_[index] = List<commsStruct>(procIDs_[index].size());

    if (doPstream && parRun())
    {
        allocatePstreamCommunicator(parentIndex, index);
    }

    return index;
}


void Foam::UPstream::freeCommunicator
(
    const label communicator,
    const bool doPstream
)
{
    if (debug)
    {
        Pout<< "Communicators : Freeing communicator " << communicator << endl
            << "    parent   : " << parentCommunicator_[communicator] << endl
            << "    myProcNo : " << myProcNo_[communicator] << endl
            << endl;
    }

    if (doPstream && parRun())
    {
        freePstreamCommunicator(communicator);
    }
    myProcNo_[communicator] = -1;
    //procIDs_[communicator].clear();
    parentCommunicator_[communicator] = -1;
    linearCommunication_[communicator].clear();
    treeCommunication_[communicator].clear();

    freeComms_.push(communicator);
}


void Foam::UPstream::freeCommunicators(const bool doPstream)
{
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freeCommunicator(communicator, doPstream);
        }
    }
}


int Foam::UPstream::baseProcNo(const label myComm, const int myProcID)
{
    int procID = myProcID;
    label comm = myComm;

    while (parent(comm) != -1)
    {
        const List<int>& parentRanks = UPstream::procID(comm);
        procID = parentRanks[procID];
        comm = UPstream::parent(comm);
    }

    return procID;
}


Foam::label Foam::UPstream::procNo(const label myComm, const int baseProcID)
{
    const List<int>& parentRanks = procID(myComm);
    label parentComm = parent(myComm);

    if (parentComm == -1)
    {
        return findIndex(parentRanks, baseProcID);
    }
    else
    {
        const label parentRank = procNo(parentComm, baseProcID);
        return findIndex(parentRanks, parentRank);
    }
}


Foam::label Foam::UPstream::procNo
(
    const label myComm,
    const label currentComm,
    const int currentProcID
)
{
    label physProcID = UPstream::baseProcNo(currentComm, currentProcID);
    return procNo(myComm, physProcID);
}


template<>
Foam::UPstream::commsStruct&
Foam::UList<Foam::UPstream::commsStruct>::operator[](const label procID)
{
    UPstream::commsStruct& t = v_[procID];

    if (t.allBelow().size() + t.allNotBelow().size() + 1 != size())
    {
        // Not yet allocated

        label above(-1);
        labelList below(0);
        labelList allBelow(0);

        if (size() < UPstream::nProcsSimpleSum)
        {
            // Linear schedule

            if (procID == 0)
            {
                below.setSize(size()-1);
                for (label procI = 1; procI < size(); procI++)
                {
                    below[procI-1] = procI;
                }
            }
            else
            {
                above = 0;
            }
        }
        else
        {
            // Use tree like schedule. For 8 procs:
            // (level 0)
            //      0 receives from 1
            //      2 receives from 3
            //      4 receives from 5
            //      6 receives from 7
            // (level 1)
            //      0 receives from 2
            //      4 receives from 6
            // (level 2)
            //      0 receives from 4
            //
            // The sends/receives for all levels are collected per processor
            // (one send per processor; multiple receives possible) creating
            // a table:
            //
            // So per processor:
            // proc     receives from   sends to
            // ----     -------------   --------
            //  0       1,2,4           -
            //  1       -               0
            //  2       3               0
            //  3       -               2
            //  4       5               0
            //  5       -               4
            //  6       7               4
            //  7       -               6

            label mod = 0;

            for (label step = 1; step < size(); step = mod)
            {
                mod = step * 2;

                if (procID % mod)
                {
                    above = procID - (procID % mod);
                    break;
                }
                else
                {
                    for
                    (
                        label j = procID + step;
                        j < size() && j < procID + mod;
                        j += step
                    )
                    {
                        below.append(j);
                    }
                    for
                    (
                        label j = procID + step;
                        j < size() && j < procID + mod;
                        j++
                    )
                    {
                        allBelow.append(j);
                    }
                }
            }
        }
        t = UPstream::commsStruct(size(), procID, above, below, allBelow);
    }
    return t;
}


template<>
const Foam::UPstream::commsStruct&
Foam::UList<Foam::UPstream::commsStruct>::operator[](const label procID) const
{
    return const_cast<UList<UPstream::commsStruct>&>(*this).operator[](procID);
}


namespace Foam
{
    // Register re-reader
    class addcommsTypeToOpt
    :
        public ::Foam::simpleRegIOobject
    {
    public:

        addcommsTypeToOpt(const char* name)
        :
            ::Foam::simpleRegIOobject(Foam::debug::addOptimisationObject, name)
        {}

        virtual ~addcommsTypeToOpt() = default;

        virtual void readData(Foam::Istream& is)
        {
            UPstream::defaultCommsType =
                UPstream::commsTypeNames.read(is);
        }

        virtual void writeData(Foam::Ostream& os) const
        {
            os << UPstream::commsTypeNames[UPstream::defaultCommsType];
        }
    };

    addcommsTypeToOpt addcommsTypeToOpt_("commsType");
}

Foam::label Foam::UPstream::worldComm(0);

Foam::label Foam::UPstream::warnComm(-1);

Foam::FixedList<MPI_Datatype, Foam::MaxMpiCustomTypeLengthInBytes + 1> Foam::UPstream::dataTypes{};

int Foam::UPstream::nPollProcInterfaces
(
    Foam::debug::optimisationSwitch("nPollProcInterfaces", 0)
);
registerOptSwitch
(
    "nPollProcInterfaces",
    int,
    Foam::UPstream::nPollProcInterfaces
);


int Foam::UPstream::maxCommsSize
(
    Foam::debug::optimisationSwitch("maxCommsSize", 0)
);
registerOptSwitch
(
    "maxCommsSize",
    int,
    Foam::UPstream::maxCommsSize
);


const int Foam::UPstream::mpiBufferSize
(
    Foam::debug::optimisationSwitch("mpiBufferSize", 0)
);


// ************************************************************************* //


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

static void attachOurBuffers()
{
    if (ourBuffers)
    {
        return;  // Already attached
    }
    ourBuffers = true;

    // Use UPstream::mpiBufferSize (optimisationSwitch),
    // but allow override with MPI_BUFFER_SIZE env variable (int value)

#ifndef SGIMPI
    int len = 0;

    const std::string str(Foam::getEnv("MPI_BUFFER_SIZE"));
    if (str.empty() || !Foam::read(str.c_str(), len) || len <= 0)
    {
        len = Foam::UPstream::mpiBufferSize;
    }

    if (len < minBufLen)
    {
        len = minBufLen;
    }

    if (Foam::UPstream::debug)
    {
        Foam::Pout<< "UPstream::init : buffer-size " << len << '\n';
    }

    char* buf = new char[len];

    if (MPI_SUCCESS != MPI_Buffer_attach(buf, len))
    {
        delete[] buf;
        Foam::Pout<< "UPstream::init : could not attach buffer\n";
    }
#endif
}


static void detachOurBuffers()
{
    if (!ourBuffers)
    {
        return;  // Nothing to detach
    }
    ourBuffers = false;

    // Some MPI notes suggest that the return code is MPI_SUCCESS when
    // no buffer is attached.
    // Be extra careful and require a non-zero size as well.

#ifndef SGIMPI
    int len = 0;
    char* buf = nullptr;

    if (MPI_SUCCESS == MPI_Buffer_detach(&buf, &len) && len)
    {
        delete[] buf;
    }
#endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NOTE:
// valid parallel options vary between implementations, but flag common ones.
// if they are not removed by MPI_Init(), the subsequent argument processing
// will notice that they are wrong
void Foam::UPstream::addValidParOptions(HashTable<string>& validParOptions)
{
    validParOptions.insert("np", "");
    validParOptions.insert("p4pg", "PI file");
    validParOptions.insert("p4wd", "directory");
    validParOptions.insert("p4amslave", "");
    validParOptions.insert("p4yourname", "hostname");
    validParOptions.insert("machinefile", "machine file");
}


bool Foam::UPstream::initNull()
{
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init\n"
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        if (debug)
        {
            Pout<< "UPstream::initNull : was already initialized\n";
        }
    }
    else
    {
        // Not already initialized

        MPI_Init_thread
        (
            nullptr,    // argc
            nullptr,    // argv
            MPI_THREAD_SINGLE,
            &flag       // provided_thread_support
        );

        ourMpi = true;
    }

    // Could also attach buffers etc.

    return true;
}


bool Foam::UPstream::init(int& argc, char**& argv, const bool needsThread)
{
    int numprocs = 0, myRank = 0;
    int provided_thread_support = 0;
    int flag = 0;

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized - this is an error
        FatalErrorInFunction
            << "MPI was already finalized - cannot perform MPI_Init" << endl
            << Foam::abort(FatalError);

        return false;
    }

    MPI_Initialized(&flag);
    if (flag)
    {
        // Already initialized.
        // Warn if we've called twice, but skip if initialized externally

        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already initialized - cannot perform MPI_Init" << nl
                << "This could indicate an application programming error!"
                << endl;

            return true;
        }
        else if (debug)
        {
            Pout<< "UPstream::init : was already initialized\n";
        }
    }
    else
    {
        MPI_Init_thread
        (
            &argc,
            &argv,
            (
                needsThread
              ? MPI_THREAD_MULTIPLE
              : MPI_THREAD_SINGLE
            ),
            &provided_thread_support
        );

        ourMpi = true;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if (debug)
    {
        Pout<< "UPstream::init : procs=" << numprocs
            << " rank:" << myRank << endl;
    }

    if (numprocs <= 1)
    {
        FatalErrorInFunction
            << "attempt to run parallel on 1 processor"
            << Foam::abort(FatalError);
    }

    // Initialise parallel structure
    setParRun(numprocs, provided_thread_support == MPI_THREAD_MULTIPLE);

    generateCustomMpiTypes();
    attachOurBuffers();

    return true;
}

void Foam::UPstream::generateCustomMpiTypes() {
    for (int i = 1; i < MaxMpiCustomTypeLengthInBytes; i++) {
        MPI_Type_contiguous(i, MPI_BYTE, &dataTypes[i]);
        MPI_Type_commit(&dataTypes[i]);
    }
}

void Foam::UPstream::exit(int errnum)
{
    if (debug)
    {
        Pout<< "UPstream::exit\n";
    }

    int flag = 0;

    MPI_Initialized(&flag);
    if (!flag)
    {
        // Not initialized - just exit
        std::exit(errnum);
    }

    MPI_Finalized(&flag);
    if (flag)
    {
        // Already finalized elsewhere?
        if (ourMpi)
        {
            WarningInFunction
                << "MPI was already finalized (by a connected program?)\n";
        }
        else if (debug)
        {
            Pout<< "UPstream::exit : was already finalized\n";
        }
    }
    else
    {
        detachOurBuffers();
    }


    const label nOutstanding = PstreamGlobals::outstandingRequests_.size();
    if (nOutstanding)
    {
        PstreamGlobals::outstandingRequests_.clear();

        WarningInFunction
            << "There were still " << nOutstanding
            << " outstanding MPI_Requests." << nl
            << "Which means your code exited before doing a "
            << " UPstream::waitRequests()." << nl
            << "This should not happen for a normal code exit."
            << nl;
    }

    // Clean mpi communicators
    forAll(myProcNo_, communicator)
    {
        if (myProcNo_[communicator] != -1)
        {
            freePstreamCommunicator(communicator);
        }
    }

    if (!flag)
    {
        // MPI not already finalized

        if (!ourMpi)
        {
            WarningInFunction
                << "Finalizing MPI, but was initialized elsewhere\n";
        }

        if (errnum == 0)
        {
            MPI_Finalize();
        }
        else
        {
            MPI_Abort(MPI_COMM_WORLD, errnum);
        }
    }

    std::exit(errnum);
}


void Foam::UPstream::abort()
{
    MPI_Abort(MPI_COMM_WORLD, 1);

    // This will not be reached, but ensures the function can be marked
    // [[noreturn]]:
    std::abort();
}

void Foam::UPstream::allToAll
(
    const labelUList& sendData,
    labelUList& recvData,
    const label communicator
)
{
    label np = nProcs(communicator);

    if (sendData.size() != np || recvData.size() != np)
    {
        FatalErrorInFunction
            << "Size of sendData " << sendData.size()
            << " or size of recvData " << recvData.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        recvData.deepCopy(sendData);
    }
    else
    {
        if
        (
            MPI_Alltoall
            (
                // NOTE: const_cast is a temporary hack for
                // backward-compatibility with versions of OpenMPI < 1.7.4
                // OpenMPI 1.7.4 is more than ten years old.
                const_cast<label*>(sendData.begin()),
                sizeof(label),
                MPI_BYTE,
                recvData.begin(),
                sizeof(label),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoall failed for " << sendData
                << " on communicator " << communicator
                << Foam::abort(FatalError);
        }

    }
}


void Foam::UPstream::allToAll
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,

    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        sendSizes.size() != np
     || sendOffsets.size() != np
     || recvSizes.size() != np
     || recvOffsets.size() != np
    )
    {
        FatalErrorInFunction
            << "Size of sendSize " << sendSizes.size()
            << ", sendOffsets " << sendOffsets.size()
            << ", recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        if (recvSizes[0] != sendSizes[0])
        {
            FatalErrorInFunction
                << "Bytes to send " << sendSizes[0]
                << " does not equal bytes to receive " << recvSizes[0]
                << Foam::abort(FatalError);
        }
        std::memmove(recvData, &sendData[sendOffsets[0]], recvSizes[0]);
    }
    else
    {

        if
        (
            MPI_Alltoallv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                MPI_BYTE,
                PstreamGlobals::MPICommunicators_[communicator]
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Alltoallv failed for sendSizes " << sendSizes
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

    }
}


void Foam::UPstream::gather
(
    const char* sendData,
    int sendSize,

    char* recvData,
    const UList<int>& recvSizes,
    const UList<int>& recvOffsets,
    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        UPstream::master(communicator)
     && (recvSizes.size() != np || recvOffsets.size() < np)
    )
    {
        // Note: allow recvOffsets to be e.g. 1 larger than np so we
        // can easily loop over the result

        FatalErrorInFunction
            << "Size of recvSizes " << recvSizes.size()
            << " or recvOffsets " << recvOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, sendSize);
    }
    else
    {

        if
        (
            MPI_Gatherv
            (
                const_cast<char*>(sendData),
                sendSize,
                MPI_BYTE,
                recvData,
                const_cast<int*>(recvSizes.begin()),
                const_cast<int*>(recvOffsets.begin()),
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Gatherv failed for sendSize " << sendSize
                << " recvSizes " << recvSizes
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

    }
}


void Foam::UPstream::scatter
(
    const char* sendData,
    const UList<int>& sendSizes,
    const UList<int>& sendOffsets,

    char* recvData,
    int recvSize,
    const label communicator
)
{
    label np = nProcs(communicator);

    if
    (
        UPstream::master(communicator)
     && (sendSizes.size() != np || sendOffsets.size() != np)
    )
    {
        FatalErrorInFunction
            << "Size of sendSizes " << sendSizes.size()
            << " or sendOffsets " << sendOffsets.size()
            << " is not equal to the number of processors in the domain "
            << np
            << Foam::abort(FatalError);
    }

    if (!UPstream::parRun())
    {
        std::memmove(recvData, sendData, recvSize);
    }
    else
    {

        if
        (
            MPI_Scatterv
            (
                const_cast<char*>(sendData),
                const_cast<int*>(sendSizes.begin()),
                const_cast<int*>(sendOffsets.begin()),
                MPI_BYTE,
                recvData,
                recvSize,
                MPI_BYTE,
                0,
                MPI_Comm(PstreamGlobals::MPICommunicators_[communicator])
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Scatterv failed for sendSizes " << sendSizes
                << " sendOffsets " << sendOffsets
                << " communicator " << communicator
                << Foam::abort(FatalError);
        }

    }
}


void Foam::UPstream::allocatePstreamCommunicator
(
    const label parentIndex,
    const label index
)
{
    if (index == PstreamGlobals::MPIGroups_.size())
    {
        // Extend storage with dummy values
        MPI_Group newGroup = MPI_GROUP_NULL;
        PstreamGlobals::MPIGroups_.append(newGroup);
        MPI_Comm newComm = MPI_COMM_NULL;
        PstreamGlobals::MPICommunicators_.append(newComm);
    }
    else if (index > PstreamGlobals::MPIGroups_.size())
    {
        FatalErrorInFunction
            << "PstreamGlobals out of sync with UPstream data. Problem."
            << Foam::exit(FatalError);
    }


    if (parentIndex == -1)
    {
        // Allocate world communicator

        if (index != UPstream::worldComm)
        {
            FatalErrorInFunction
                << "world communicator should always be index "
                << UPstream::worldComm << Foam::exit(FatalError);
        }

        PstreamGlobals::MPICommunicators_[index] = MPI_COMM_WORLD;
        MPI_Comm_group(MPI_COMM_WORLD, &PstreamGlobals::MPIGroups_[index]);
        MPI_Comm_rank
        (
            PstreamGlobals::MPICommunicators_[index],
           &myProcNo_[index]
        );

        // Set the number of processes to the actual number
        int numProcs;
        MPI_Comm_size(PstreamGlobals::MPICommunicators_[index], &numProcs);

        //procIDs_[index] = identity(numProcs);
        procIDs_[index].setSize(numProcs);
        forAll(procIDs_[index], i)
        {
            procIDs_[index][i] = i;
        }
    }
    else
    {
        // Create new group
        MPI_Group_incl
        (
            PstreamGlobals::MPIGroups_[parentIndex],
            procIDs_[index].size(),
            procIDs_[index].begin(),
           &PstreamGlobals::MPIGroups_[index]
        );

        // Create new communicator
        MPI_Comm_create
        (
            PstreamGlobals::MPICommunicators_[parentIndex],
            PstreamGlobals::MPIGroups_[index],
           &PstreamGlobals::MPICommunicators_[index]
        );

        if (PstreamGlobals::MPICommunicators_[index] == MPI_COMM_NULL)
        {
            myProcNo_[index] = -1;
        }
        else
        {
            if
            (
                MPI_Comm_rank
                (
                    PstreamGlobals::MPICommunicators_[index],
                   &myProcNo_[index]
                )
            )
            {
                FatalErrorInFunction
                    << "Problem :"
                    << " when allocating communicator at " << index
                    << " from ranks " << procIDs_[index]
                    << " of parent " << parentIndex
                    << " cannot find my own rank"
                    << Foam::exit(FatalError);
            }
        }
    }
}


void Foam::UPstream::freePstreamCommunicator(const label communicator)
{
    if (communicator != UPstream::worldComm)
    {
        if (PstreamGlobals::MPICommunicators_[communicator] != MPI_COMM_NULL)
        {
            // Free communicator. Sets communicator to MPI_COMM_NULL
            MPI_Comm_free(&PstreamGlobals::MPICommunicators_[communicator]);
        }
        if (PstreamGlobals::MPIGroups_[communicator] != MPI_GROUP_NULL)
        {
            // Free greoup. Sets group to MPI_GROUP_NULL
            MPI_Group_free(&PstreamGlobals::MPIGroups_[communicator]);
        }
    }
}


Foam::label Foam::UPstream::nRequests()
{
    return PstreamGlobals::outstandingRequests_.size();
}


void Foam::UPstream::resetRequests(const label i)
{
    if (i < PstreamGlobals::outstandingRequests_.size())
    {
        PstreamGlobals::outstandingRequests_.setSize(i);
    }
}


void Foam::UPstream::waitRequests(const label start)
{
    if (UPstream::debug)
    {
        Pout<< "UPstream::waitRequests : starting wait for "
            << PstreamGlobals::outstandingRequests_.size()-start
            << " outstanding requests starting at " << start << endl;
    }

    if (PstreamGlobals::outstandingRequests_.size())
    {
        SubList<MPI_Request> waitRequests
        (
            PstreamGlobals::outstandingRequests_,
            PstreamGlobals::outstandingRequests_.size() - start,
            start
        );


        if
        (
            MPI_Waitall
            (
                waitRequests.size(),
                waitRequests.begin(),
                MPI_STATUSES_IGNORE
            )
        )
        {
            FatalErrorInFunction
                << "MPI_Waitall returned with error" << Foam::endl;
        }


        resetRequests(start);
    }

    if (debug)
    {
        Pout<< "UPstream::waitRequests : finished wait." << endl;
    }
}


void Foam::UPstream::waitRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::waitRequest : starting wait for request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }


    if
    (
        MPI_Wait
        (
           &PstreamGlobals::outstandingRequests_[i],
            MPI_STATUS_IGNORE
        )
    )
    {
        FatalErrorInFunction
            << "MPI_Wait returned with error" << Foam::endl;
    }


    if (debug)
    {
        Pout<< "UPstream::waitRequest : finished wait for request:" << i
            << endl;
    }
}


bool Foam::UPstream::finishedRequest(const label i)
{
    if (debug)
    {
        Pout<< "UPstream::finishedRequest : checking request:" << i
            << endl;
    }

    if (i >= PstreamGlobals::outstandingRequests_.size())
    {
        FatalErrorInFunction
            << "There are " << PstreamGlobals::outstandingRequests_.size()
            << " outstanding send requests and you are asking for i=" << i
            << nl
            << "Maybe you are mixing blocking/non-blocking comms?"
            << Foam::abort(FatalError);
    }

    int flag;
    MPI_Test
    (
       &PstreamGlobals::outstandingRequests_[i],
       &flag,
        MPI_STATUS_IGNORE
    );

    if (debug)
    {
        Pout<< "UPstream::finishedRequest : finished request:" << i
            << endl;
    }

    return flag != 0;
}


int Foam::UPstream::allocateTag(const char* s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


int Foam::UPstream::allocateTag(const word& s)
{
    int tag;
    if (PstreamGlobals::freedTags_.size())
    {
        tag = PstreamGlobals::freedTags_.remove();
    }
    else
    {
        tag = PstreamGlobals::nTags_++;
    }

    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = 'X';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::allocateTag " << s
            << " : tag:" << tag
            << endl;
    }

    return tag;
}


void Foam::UPstream::freeTag(const char* s, const int tag)
{
    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}


void Foam::UPstream::freeTag(const word& s, const int tag)
{
    if (debug)
    {
        //if (UPstream::lateBlocking > 0)
        //{
        //    string& poutp = Pout.prefix();
        //    poutp[poutp.size()-2*(UPstream::lateBlocking+2)+tag] = ' ';
        //    Perr.prefix() = Pout.prefix();
        //}
        Pout<< "UPstream::freeTag " << s << " tag:" << tag << endl;
    }
    PstreamGlobals::freedTags_.append(tag);
}

void Foam::Pstream::barrier(const Foam::label comm)
{
    if (Pstream::parRun())
    {
        MPI_Barrier(PstreamGlobals::MPICommunicators_[comm]);
    }
}
