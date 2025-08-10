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

Description
    Miscellaneous tests for HashTable

\*---------------------------------------------------------------------------*/

#include "containers/HashTables/HashTable/HashTable.H"
#include "containers/HashTables/HashPtrTable/HashPtrTable.H"
#include "containers/HashTables/Map/Map.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    HashTable<label, Foam::string> table1
    {
        {"kjhk", 10},
        {"kjhk2", 12}
    };

    Info<< "table1: " << table1 << nl
        << "toc: " << table1.toc() << endl;

    HashTable<label, label, Hash<label>> table2
    {
        {3, 10},
        {5, 12},
        {7, 16}
    };

    Info<< "table2: " << table2 << nl
        << "toc: " << table2.toc() << endl;

    Map<label> table3(1);
    table3.transfer(table2);

    Info<< "table2: " << table2 << nl
        << "toc: " << table2.toc() << endl;

    Info<< "table3: " << table3 << nl
        << "toc: " << table3.toc() << endl;

    Map<label> table4(table3.xfer());

    Info<< "table3: " << table3 << nl
        << "toc: " << table3.toc() << endl;

    Info<< "table4: " << table4 << nl
        << "toc: " << table4.toc() << endl;

    HashPtrTable<label, Foam::string> ptable1(0);
    ptable1.insert("kjhkjh", new label(10));

    Info<< "PtrTable toc: " << ptable1.toc() << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
