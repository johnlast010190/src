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

#include "utilities/surfaceTools/meshSurfaceMapper2D/meshSurfaceMapper2D.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "utilities/surfaceTools/meshSurfacePartitioner/meshSurfacePartitioner.H"
#include "utilities/meshes/polyMeshGen2DEngine/polyMeshGen2DEngine.H"
#include "utilities/meshes/triSurf/triSurf.H"
#include "utilities/triSurfaceTools/triSurfacePartitioner/triSurfacePartitioner.H"
#include "include/demandDrivenData.H"
#include "utilities/octrees/meshOctree/meshOctree.H"
#include "utilities/helperFunctions/helperFunctions.H"

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshSurfaceMapper2D::findActiveBoundaryEdges()
{
    const VRWGraph& edgeFaces = surfaceEngine_.edgeFaces();

    const polyMeshGen2DEngine& mesh2DEngine = this->mesh2DEngine();
    const boolList& activeFace = mesh2DEngine.activeFace();

    const label startFace = surfaceEngine_.mesh().boundaries()[0].patchStart();

    activeBoundaryEdges_.clear();

    //- check which edges are at the boundary
    forAll(edgeFaces, edgeI)
    {
        if (edgeFaces.sizeOfRow(edgeI) == 2)
        {
            const bool active0 = activeFace[startFace + edgeFaces(edgeI, 0)];
            const bool active1 = activeFace[startFace + edgeFaces(edgeI, 1)];

            if (active0 && active1)
                activeBoundaryEdges_.append(edgeI);
        }
    }

    if (Pstream::parRun())
    {
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndEdgeAddressing();
        const Map<label>& otherProc = surfaceEngine_.otherEdgeFaceAtProc();

        const DynList<label>& neiProcs = surfaceEngine_.beNeiProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(neiProcs, i)
                exchangeData[neiProcs[i]].clear();

        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label beI = it();

            if (activeFace[startFace + edgeFaces(beI, 0)])
            {
                exchangeData[otherProc[beI]].append(it.key());
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const label beI = globalToLocal[receivedData[i]];

            if (activeFace[startFace + edgeFaces(beI, 0)])
                activeBoundaryEdges_.append(beI);
        }
    }
}

void meshSurfaceMapper2D::create2DEngine() const
{
    polyMeshGen& mesh = const_cast<polyMeshGen&>(surfaceEngine_.mesh());
    mesh2DEnginePtr_ = new polyMeshGen2DEngine(mesh);
}

void meshSurfaceMapper2D::createTriSurfacePartitioner() const
{
    surfPartitionerPtr_ = new triSurfacePartitioner(meshOctree_.surface());
}

void meshSurfaceMapper2D::createMeshSurfacePartitioner() const
{
    meshPartitionerPtr_ = new meshSurfacePartitioner(surfaceEngine_);
}

void meshSurfaceMapper2D::clearOut()
{
    deleteDemandDrivenData(mesh2DEnginePtr_);
    deleteDemandDrivenData(surfPartitionerPtr_);
    deleteDemandDrivenData(meshPartitionerPtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

meshSurfaceMapper2D::meshSurfaceMapper2D
(
    const meshSurfaceEngine& mse,
    const meshOctree& octree
)
:
    surfaceEngine_(mse),
    meshOctree_(octree),
    mesh2DEnginePtr_(nullptr),
    surfPartitionerPtr_(nullptr),
    meshPartitionerPtr_(nullptr)
{
    if (Pstream::parRun())
    {
        //- allocate bpAtProcs and other addressing
        //- this is done here to prevent possible deadlocks
        surfaceEngine_.bpAtProcs();
    }

    findActiveBoundaryEdges();

    createMeshSurfacePartitioner();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceMapper2D::~meshSurfaceMapper2D()
{
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
