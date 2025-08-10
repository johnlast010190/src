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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description
    Finds feature edges and corners of a triangulated surface

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "primitives/strings/fileName/fileName.H"
#include "utilities/meshes/triSurf/triSurf.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "include/OSspecific.H"
#include "include/demandDrivenData.H"
#include "utilities/triSurfaceTools/triSurfaceImportSurfaceAsSubset/triSurfaceImportSurfaceAsSubset.H"
#include <cstdlib>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("master surface file");
    argList::validArgs.append("import surface file");

    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName importFileName(args.args()[2]);

    triSurf originalSurface(inFileName);

    triSurf importedSurface(importFileName);

    triSurfaceImportSurfaceAsSubset importSurf(originalSurface);

    importSurf.addSurfaceAsSubset(importedSurface, importFileName.lessExt());

    if (inFileName.ext() == "fms")
    {
        originalSurface.writeSurface(inFileName);
    }
    else
    {
        fileName newName = inFileName.lessExt();
        newName.append(".fms");
        Warning << "Writting surface as " << newName
            << " to preserve the subset!!" << endl;

        originalSurface.writeSurface(newName);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
