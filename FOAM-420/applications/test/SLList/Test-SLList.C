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

Application

Description

\*---------------------------------------------------------------------------*/

#include "include/OSspecific.H"

#include "db/IOstreams/IOstreams.H"
#include "SLList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    SLList<scalar> myList{2.1, 3.4};
    myList = {2.1, 3.4, 4.3};

    for (int i = 0; i<10; i++)
    {
        myList.append(1.3*i);
    }

    myList.append(100.3);
    myList.append(500.3);

    Info<< nl << "And again using STL iterator: " << nl << endl;

    for (const auto& val : myList)
    {
        Info<< "element:" << val << endl;
    }

    Info<< nl << "And again using STL const_iterator: " << nl << endl;

    const SLList<scalar>& const_myList = myList;

    forAllConstIters(const_myList, iter)
    {
        Info<< "element:" << *iter << endl;
    }

    forAllIters(myList, iter)
    {
        Info<< "Removing element:" << *iter << endl;
        myList.remove(iter);
    }

    for (const auto& val : const_myList)
    {
        Info<< "element:" << val << endl;
    }


    for (int i = 0; i<10; i++)
    {
        myList.append(1.3*i);
    }

    myList.append(100.3);
    myList.append(500.3);

    Info<< nl << "Testing transfer: " << nl << endl;
    Info<< "original: " << myList << endl;

    SLList<scalar> newList;
    newList.transfer(myList);

    Info<< nl << "source: " << myList << nl
        << nl << "target: " << newList << endl;


    Info<< nl << "Done." << endl;
    return 0;
}


// ************************************************************************* //
