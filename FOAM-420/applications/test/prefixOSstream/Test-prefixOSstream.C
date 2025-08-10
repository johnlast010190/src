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
    prefixOSstreamTest

Description

\*---------------------------------------------------------------------------*/

#include "include/OSspecific.H"

#include "db/IOstreams/IOstreams.H"
#include "db/IOstreams/StringStreams/IStringStream.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/Vector/vector/vector.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "global/argList/argList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList args(argc, argv);

    //Pout.prefix() = '[' + name(Pstream::myProcNo()) + "] ";

    List<vector> list(IStringStream("1 ((0 1 2))")());
    Info<< list << endl;

    List<vector> list2
    (
        IStringStream
        (
            "(\
                 (0 1 2)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
                 (3 4 5)\
             )"
        )()
    );
    Pout<< list2 << endl;

    Info<< findIndex(list2, vector(3, 4, 5)) << endl;

    return 0;
}


// ************************************************************************* //
