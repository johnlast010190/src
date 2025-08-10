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

#include "utilities/surfaceTools/meshSurfaceMapper/meshSurfaceMapper.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "utilities/surfaceTools/meshSurfacePartitioner/meshSurfacePartitioner.H"
#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/triSurfaceTools/triSurfacePartitioner/triSurfacePartitioner.H"
#include "include/demandDrivenData.H"
#include "utilities/octrees/meshOctree/meshOctree.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper::createMeshSurfacePartitioner() const
{
    surfaceEnginePartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine_);
}

void meshSurfaceMapper::createTriSurfacePartitioner() const
{
    surfPartitionerPtr_ = new triSurfacePartitioner(meshOctree_.surface());
}

void meshSurfaceMapper::clearOut()
{
    if (deletePartitioner_)
        deleteDemandDrivenData(surfaceEnginePartitionerPtr_);
    deleteDemandDrivenData(surfPartitionerPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceMapper::meshSurfaceMapper
(
    const meshSurfaceEngine& mse,
    const meshOctree& octree
)
:
    surfaceEngine_(mse),
    meshOctree_(octree),
    surfaceEnginePartitionerPtr_(nullptr),
    deletePartitioner_(true),
    surfPartitionerPtr_(nullptr)
{
    if (Pstream::parRun())
    {
        //- allocate bpAtProcs and other addressing
        //- this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }
}

meshSurfaceMapper::meshSurfaceMapper
(
    const meshSurfacePartitioner& mPart,
    const meshOctree& octree
)
:
    surfaceEngine_(mPart.surfaceEngine()),
    meshOctree_(octree),
    surfaceEnginePartitionerPtr_(&mPart),
    deletePartitioner_(false),
    surfPartitionerPtr_(nullptr)
{
    if (Pstream::parRun())
    {
        //- allocate bpAtProcs and other addressing
        //- this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceMapper::~meshSurfaceMapper()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
