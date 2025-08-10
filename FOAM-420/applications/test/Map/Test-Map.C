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
    Test-Map

Description

\*---------------------------------------------------------------------------*/

#include "containers/HashTables/Map/Map.H"
#include <map>
#include "db/IOstreams/IOstreams.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    Map<bool> banana{{5, true}};

    // Taking a const iterator from find does not work!
    // Also, fails later on op==
    Map<bool>::const_iterator bananaIter = banana.find(5);

    // This works but now I can change the value.
    //Map<bool>::iterator bananaIter = banana.find(5);

    if (!bananaIter.found()) // same as  (bananaIter == banana.end())
    {
        Info<< "not found" << endl;
    }
    else
    {
        Info<< "5 is " << bananaIter() << endl;
    }

    // Same with STL
    Info<< "Same with STL" << endl;

    std::map<label, bool> STLbanana{{5, true}};
    std::map<label, bool>::const_iterator STLbananaIter = STLbanana.find(5);

    if (STLbananaIter == STLbanana.end())
    {
        Info<< "not found" << endl;
    }
    else
    {
        Info<< "5 is " << STLbananaIter->second << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
