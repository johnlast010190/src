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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

Application
    parallelTest

Description
    Test for various parallel routines.

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

using namespace Foam;


void testMapDistribute()
{
    Random rndGen(43544*Pstream::myProcNo());

    // Generate random data.
    List<Tuple2<label, List<scalar>>> complexData(100);
    forAll(complexData, i)
    {
    #if defined(WIN64) || defined(WIN32)
        complexData[i].first() = rndGen.bit();
        #else
        complexData[i].first() = rndGen.position<label>(0, Pstream::nProcs()-1);
        #endif
        complexData[i].second().setSize(3);
        complexData[i].second()[0] = 1;
        complexData[i].second()[1] = 2;
        complexData[i].second()[2] = 3;
    }

    // Send all ones to processor indicated by .first()

    // Count how many to send
    labelList nSend(Pstream::nProcs(), 0);
    forAll(complexData, i)
    {
        label procI = complexData[i].first();
        nSend[procI]++;
    }

    // Collect items to be sent
    labelListList sendMap(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendMap[procI].setSize(nSend[procI]);
    }
    nSend = 0;
    forAll(complexData, i)
    {
        label procI = complexData[i].first();
        sendMap[procI][nSend[procI]++] = i;
    }

    // Sync how many to send
    labelList nRecv;
    Pstream::exchangeSizes(sendMap, nRecv);

    // Collect items to be received
    labelListList recvMap(Pstream::nProcs());
    forAll(recvMap, procI)
    {
        recvMap[procI].setSize(nRecv[procI]);
    }

    label constructSize = 0;
    // Construct with my own elements first
    forAll(recvMap[Pstream::myProcNo()], i)
    {
        recvMap[Pstream::myProcNo()][i] = constructSize++;
    }
    // Construct from other processors
    forAll(recvMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            forAll(recvMap[procI], i)
            {
                recvMap[procI][i] = constructSize++;
            }
        }
    }

    // Construct distribute map (destructively)
    mapDistribute map(constructSize, sendMap.xfer(), recvMap.xfer());

    // Distribute complexData
    map.distribute(complexData);

    Pout<< "complexData:" << complexData << endl;
}


template<class T>
void testTransfer(const T& input)
{
    T data = input;

    if (Pstream::master())
    {
        Perr<<"test transfer (" << (typeid(T).name()) << "): " << data << nl << endl;
    }

    if (Pstream::myProcNo() != Pstream::masterNo())
    {
        {
            Perr<< "slave sending to master " << Pstream::masterNo() << endl;
            OPstream toMaster(Pstream::commsTypes::blocking, Pstream::masterNo());
            toMaster << data;
        }

        Perr<< "slave receiving from master " << Pstream::masterNo() << endl;
        IPstream fromMaster(Pstream::commsTypes::blocking, Pstream::masterNo());
        fromMaster >> data;
        Perr<< data << endl;
    }
    else
    {
        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            Perr<< "master receiving from slave " << slave << endl;
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);
            fromSlave >> data;
            Perr<< data << endl;
        }

        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            Perr<< "master sending to slave " << slave << endl;
            OPstream toSlave(Pstream::commsTypes::blocking, slave);
            toSlave << data;
        }
    }
}


template<class T>
void testTokenized(const T& data)
{
    token tok;

    if (Pstream::master())
    {
        Perr<<"test tokenized \"" << data << "\"" << nl << endl;
    }

    if (Pstream::myProcNo() != Pstream::masterNo())
    {
        {
            Perr<< "slave sending to master " << Pstream::masterNo() << endl;
            OPstream toMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            toMaster << data;
        }

        Perr<< "slave receiving from master " << Pstream::masterNo() << endl;
        IPstream fromMaster
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo()
        );

        fromMaster >> tok;
        Perr<< tok.info() << endl;
    }
    else
    {
        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            Perr<< "master receiving from slave " << slave << endl;
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);
            fromSlave >> tok;
            Perr<< tok.info() << endl;
        }

        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            Perr<< "master sending to slave " << slave << endl;
            OPstream toSlave(Pstream::commsTypes::blocking, slave);
            toSlave << data;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    testMapDistribute();

    if (!Pstream::parRun())
    {
        Info<< "\nWarning: not parallel - skipping further tests\n" << endl;
        return 0;
    }

    Info<< "\nStarting transfers\n\n" << endl;

    testTransfer(vector(0, 1, 2));
    testTransfer(label(1234));
    testTransfer(scalar(3.14159));
    testTransfer(string("test   string"));
    testTransfer(string("  x "));
    testTransfer(word("3.141 59"));  // bad word, but transfer doesn't care

    testTokenized(label(1234));
    testTokenized(scalar(3.14159));
    testTokenized('a');
    testTokenized('$');  // will not tokenize well

    testTokenized(string("test   string1"));
    testTokenized("test   string1");
    testTokenized(word("3.141 59"));  // bad word, but transfer doesn't care

    testTokenized(string("  a "));
    testTokenized("  a ");

    testTokenized(string("  $ "));
    testTokenized("  $ ");  // reduces to 'char' and will not tokenize well

    testTokenized(string("  $$ "));
    testTokenized("  $$ "); // reduces to 'word' and is tagged as such


    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
