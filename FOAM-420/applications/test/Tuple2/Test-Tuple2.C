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
    Tuple2Test

Description

\*---------------------------------------------------------------------------*/

#include "primitives/Tuple2/Tuple2.H"
#include "primitives/ints/label/label.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "containers/Lists/List/List.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main()
{
    typedef Tuple2<label, scalar> indexedScalar;

    indexedScalar t2(1, 3.2);

    Info<< "tuple: "
        << t2 << " "
        << t2.first() << " " << t2.second() << endl;

    List<indexedScalar> list1(10);
    forAll(list1, i)
    {
        list1[i] = indexedScalar(-i, i*i);
    }

    sort(list1);

    Info<< "tuples:" << nl
        << list1
        << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
