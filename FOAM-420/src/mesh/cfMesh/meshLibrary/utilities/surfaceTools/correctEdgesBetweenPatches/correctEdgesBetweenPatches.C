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

#include "utilities/surfaceTools/correctEdgesBetweenPatches/correctEdgesBetweenPatches.H"
#include "include/demandDrivenData.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "utilities/decomposeCells/decomposeCells.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const meshSurfaceEngine& correctEdgesBetweenPatches::meshSurface() const
{
    if (!msePtr_)
        msePtr_ = new meshSurfaceEngine(mesh_);

    return *msePtr_;
}

//- delete mesh surface
void correctEdgesBetweenPatches::clearMeshSurface()
{
    deleteDemandDrivenData(msePtr_);
}

void correctEdgesBetweenPatches::replaceBoundary()
{
    clearMeshSurface();

    polyMeshGenModifier(mesh_).replaceBoundary
    (
        patchNames_,
        newBoundaryFaces_,
        newBoundaryOwners_,
        newBoundaryPatches_
    );
}

void correctEdgesBetweenPatches::decomposeCorrectedCells()
{
    if (decompose_)
    {
        clearMeshSurface();

        decomposeCells dc(mesh_);
        dc.decomposeMesh(decomposeCell_);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
correctEdgesBetweenPatches::correctEdgesBetweenPatches(polyMeshGen& mesh)
:
    mesh_(mesh),
    msePtr_(nullptr),
    patchNames_(mesh.boundaries().size()),
    patchTypes_(mesh.boundaries().size()),
    newBoundaryFaces_(),
    newBoundaryOwners_(),
    newBoundaryPatches_(),
    decomposeCell_(mesh_.cells().size(), false),
    decompose_(false)
{
    const PtrList<boundaryPatch>& boundaries = mesh_.boundaries();
    forAll(boundaries, patchI)
    {
        patchNames_[patchI] = boundaries[patchI].patchName();
        patchTypes_[patchI] = boundaries[patchI].patchType();
    }

    //decomposeProblematicFaces();

    decomposeConcaveFaces();

    patchCorrection();

    decomposeCorrectedCells();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

correctEdgesBetweenPatches::~correctEdgesBetweenPatches()
{
    deleteDemandDrivenData(msePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
