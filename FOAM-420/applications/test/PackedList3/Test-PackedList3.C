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
    (c) 2011 OpenFOAM Foundation

Application

Description

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "primitives/bools/lists/boolList.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "containers/HashTables/StaticHashTable/StaticHashTable.H"
#include "cpuTime/cpuTime.H"
#include <vector>
#include "containers/Lists/PackedList/PackedBoolList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
    const label n = 100000000;
    const label nReport = 1000000;

    cpuTime timer;

    // test inserts
    // PackedBoolList
    PackedBoolList packed;
    for (label i = 0; i < n; i++)
    {
        if ((i % nReport) == 0 && i)
        {
            Info<< "i:" << i << " in " << timer.cpuTimeIncrement() << " s"
                <<endl;
        }
        packed[i] = 1;
    }
    Info<< "insert test: " << n << " elements in "
        << timer.cpuTimeIncrement() << " s\n\n";

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
