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
    Reads the surface mesh, remove the selected facets
    and writes the modified mesh into a new file

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/triSurfaceTools/triSurfaceExtrude2DEdges/triSurfaceExtrude2DEdges.H"
#include "include/demandDrivenData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");

    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    //- read the input surface
    triSurf origSurf(inFileName);

    //- remove the selected facets
    triSurfaceExtrude2DEdges extruder(origSurf);

    const triSurf* newSurfacePtr = extruder.extrudeSurface();

    if (!newSurfacePtr)
        FatalError << "Extruding of the edge mesh failed!" << exit(FatalError);

    //- write the modified surface mesh
    newSurfacePtr->writeSurface(outFileName);

    deleteDemandDrivenData(newSurfacePtr);

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
