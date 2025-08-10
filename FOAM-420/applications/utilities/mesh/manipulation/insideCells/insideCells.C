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
    (c) 2016 OpenCFD Ltd.

Application
    insideCells

Group
    grpMeshManipulationUtilities

Description
    Picks up cells with cell centre 'inside' of surface.
    Requires surface to be closed and singly connected.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "triSurface/triSurface.H"
#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include "sets/topoSets/cellSet.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create a cellSet for cells with their centres inside the defined "
        "surface.\n"
        "Surface must be closed and singly connected."
    );

    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("cellSet");

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createPolyMesh.H"

    const fileName surfName = args[1];
    const fileName setName  = args[2];

    // Read surface
    Info<< "Reading surface from " << surfName << endl;
    triSurface surf(surfName);

    // Destination cellSet.
    cellSet insideCells(mesh, setName, IOobject::NO_READ);


    // Construct search engine on surface
    triSurfaceSearch querySurf(surf);

    boolList inside(querySurf.calcInside(mesh.cellCentres()));

    forAll(inside, celli)
    {
        if (inside[celli])
        {
            insideCells.insert(celli);
        }
    }


    Info<< "Selected " << returnReduce(insideCells.size(), sumOp<label>())
        << " of " << mesh.globalData().nTotalCells()
        << " cells" << nl << nl
        << "Writing selected cells to cellSet " << insideCells.name()
        << nl << nl
        << "Use this cellSet e.g. with subsetMesh : " << nl << nl
        << "    subsetMesh " << insideCells.name()
        << nl << endl;

    insideCells.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
