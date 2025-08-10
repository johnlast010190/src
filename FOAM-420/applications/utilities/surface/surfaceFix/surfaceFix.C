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
    (c) 2010-2011, Esi Ltd

Description
    A general purpose surface fixing tool. Please extend as appropriate.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "triSurface/triSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input file");
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName surfName(args[1]);

    Info<< "Reading surf from " << surfName << " ..." << nl << endl;

    fileName outFileName(args[2]);

    triSurface surf(surfName);

    //Fix 1
    //Check surface region names and make sure
    geometricSurfacePatchList& patches = surf.patches();

    forAll(patches, pI)
    {
        Info<< patches[pI].index() << " " << patches[pI].name() << endl;
    }

    triSurface surf2(surf, surf.patches(), surf.points());

    Info<< nl << "Writing surf to " << outFileName << " ..." << endl;

    surf2.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
