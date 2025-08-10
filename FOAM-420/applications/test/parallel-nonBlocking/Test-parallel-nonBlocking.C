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
    (c) 2012-2017 OpenFOAM Foundation

Application
    Test-parallel-nonBlocking

Description
    Test for various non-blocking parallel routines.

\*---------------------------------------------------------------------------*/

#include "containers/Lists/List/List.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistribute.H"
#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "db/IOstreams/Pstreams/IPstream.H"
#include "db/IOstreams/Pstreams/OPstream.H"
#include "primitives/Vector/vector/vector.H"
#include "db/IOstreams/IOstreams.H"
#include "primitives/random/Random/Random.H"
#include "primitives/Tuple2/Tuple2.H"
#include "db/IOstreams/Pstreams/PstreamBuffers.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"


    // Test PstreamBuffers
    // ~~~~~~~~~~~~~~~~~~~
    if (false)
    {
        Perr<< "\nStarting transfers\n" << endl;

        vector data
        (
            Pstream::myProcNo(),
            Pstream::myProcNo(),
            Pstream::myProcNo()
        );

        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            Perr<< "slave sending to master "
                << Pstream::masterNo() << endl;
            UOPstream toMaster(Pstream::masterNo(), pBufs);
            toMaster << data;
        }

        // Start sending and receiving and block
        pBufs.finishedSends();

        // Consume
        DynamicList<vector> allData;
        if (Pstream::myProcNo() == Pstream::masterNo())
        {
            // Collect my own data
            allData.append(data);

            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Perr << "master receiving from slave " << slave << endl;
                UIPstream fromSlave(slave, pBufs);
                allData.append(vector(fromSlave));
            }
        }


        // Send allData back
        PstreamBuffers pBufs2(Pstream::commsTypes::nonBlocking);
        if (Pstream::myProcNo() == Pstream::masterNo())
        {
            for
            (
                int slave=Pstream::firstSlave();
                slave<=Pstream::lastSlave();
                slave++
            )
            {
                Perr << "master sending to slave " << slave << endl;
                UOPstream toSlave(slave, pBufs2);
                toSlave << allData;
            }
        }

        // Start sending and receiving and block
        pBufs2.finishedSends();

        // Consume
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            Perr<< "slave receiving from master "
                << Pstream::masterNo() << endl;
            UIPstream fromMaster(Pstream::masterNo(), pBufs2);
            fromMaster >> allData;
            Perr<< allData << endl;
        }
    }


    // Test non-blocking reductions
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    scalar data1 = 1.0;
    label request1 = -1;
    {
        Foam::reduce(data1, sumOp<scalar>(), request1);
    }

    scalar data2 = 0.1;
    label request2 = -1;
    {
        Foam::reduce(data2, sumOp<scalar>(), request2);
    }


    // Do a non-blocking send inbetween
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UOPstream toProc(proci, pBufs);
            toProc << Pstream::myProcNo();
        }

        // Start sending and receiving and block
        pBufs.finishedSends();

        // Consume
        for (label proci = 0; proci < Pstream::nProcs(); proci++)
        {
            UIPstream fromProc(proci, pBufs);
            label data;
            fromProc >> data;

            if (data != proci)
            {
                FatalErrorInFunction
                    << "From processor " << proci << " received " << data
                    << " but expected " << proci
                    << exit(FatalError);
            }
        }
    }


    if (request1 != -1)
    {
        Pout<< "Waiting for non-blocking reduce with request " << request1
            << endl;
        Pstream::waitRequest(request1);
    }
    Info<< "Reduced data1:" << data1 << endl;

    if (request2 != -1)
    {
        Pout<< "Waiting for non-blocking reduce with request " << request1
            << endl;
        Pstream::waitRequest(request2);
    }
    Info<< "Reduced data2:" << data2 << endl;


    // Clear any outstanding requests
    Pstream::resetRequests(0);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
