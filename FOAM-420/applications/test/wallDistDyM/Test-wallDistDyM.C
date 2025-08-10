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
    (c) 2017 OpenCFD Ltd.
    (c) 2017 OpenFOAM Foundation

Description
    Calculate and write the distance-to-wall field for a moving mesh.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "fvMesh/fvMesh.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    Info<< "Mesh read in = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;

    // Calculate initial mesh-to-mesh mapping. Note that this should be
    // done under the hood, e.g. as a MeshObject
    mesh.update();

    Info<< "Time now = " << runTime.timeName() << endl;

    // Wall-reflection vectors
    const volVectorField& n = wallDist::New(mesh).n();
    n.write();

    // Wall distance
    const volScalarField& y = wallDist::New(mesh).y();
    y.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
