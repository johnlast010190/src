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

\*---------------------------------------------------------------------------*/

#include "include/demandDrivenData.H"
#include "utilities/smoothers/topology/topologicalCleaner/topologicalCleaner.H"
#include "utilities/decomposeCells/decomposeCells.H"
#include "utilities/smoothers/topology/checkCellConnectionsOverFaces/checkCellConnectionsOverFaces.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void topologicalCleaner::decomposeCells()
{
    if (!changed_)
         return;

    Foam::decomposeCells dc(mesh_);
    dc.decomposeMesh(decomposeCell_);

}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from points, cells, boundary faces, and octree
topologicalCleaner::topologicalCleaner
(
    polyMeshGen& mesh
)
:
    mesh_(mesh),
    changed_(false),
    decomposeCell_(mesh.cells().size(), false)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

topologicalCleaner::~topologicalCleaner()
{
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool topologicalCleaner::cleanTopology()
{
    checkInvalidConnectionsForVerticesCells();

    checkInvalidConnectionsForVerticesFaces();

    checkNonConsecutiveBoundaryVertices();

    checkNonMappableCells();

    checkNonMappableFaces();

    decomposeCells();

    if (checkCellConnectionsOverFaces(mesh_).checkCellGroups())
        changed_ = true;

    return changed_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
