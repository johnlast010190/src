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
#include <cstdlib>
#include <sstream>

#include "utilities/triSurfaceTools/triSurfaceDetectFeatureEdges/triSurfaceDetectFeatureEdges.H"
#include "utilities/triSurfaceTools/triSurfacePatchManipulator/triSurfacePatchManipulator.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList::validOptions.insert("angle", "scalar");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    if (outFileName == inFileName)
    {
        FatalErrorIn(args.executable())
            << "Output file " << outFileName
            << " would overwrite the input file."
            << exit(FatalError);
    }

    scalar tol(45.0);
    if (args.options().found("angle"))
    {
        const scalar ang = readScalar(IStringStream(args.options()["angle"])());
        tol = ang;
    }
    else
    {
        Info<< "Using 45 deg as default angle!" << endl;
    }

    triSurf originalSurface(inFileName);

    triSurfaceDetectFeatureEdges edgeDetector(originalSurface, tol);
    edgeDetector.detectFeatureEdges();

    if (outFileName.ext() == "fms" || outFileName.ext() == "FMS")
    {
        Info<< "Writing : " << outFileName << endl;
        originalSurface.writeSurface(outFileName);
    }
    else
    {
        triSurfacePatchManipulator manipulator(originalSurface);
        const triSurf* newSurfPtr = manipulator.surfaceWithPatches();

        Info<< "Writing : " << outFileName << endl;
        newSurfPtr->writeSurface(outFileName);

        deleteDemandDrivenData(newSurfPtr);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
