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
    Scales the mesh into other units.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "utilities/meshes/polyMeshGen/polyMeshGen.H"
#include "utilities/helperFunctions/helperFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::validArgs.append("scalingFactor");

#include "include/setRootCase.H"
#include "include/createTime.H"

    const scalar scalingFactor(help::textToScalar(args.args()[1]));

    Info<< "Scaling mesh vertices by a factor " << scalingFactor << endl;

    //- read the mesh from disk
    polyMeshGen pmg(runTime);

    Info<< "Reading mesh" << endl;
    pmg.read();

    //- scale the points
    pointFieldPMG& pts = pmg.points();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 100)
    # endif
    forAll(pts, pointI)
        pts[pointI] *= scalingFactor;

    //- write the mesh back on disk
    Info<< "Writting scaled mesh" << endl;
    pmg.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
