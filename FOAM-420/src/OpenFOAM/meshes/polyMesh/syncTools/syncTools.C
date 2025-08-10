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

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/syncTools/syncTools.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

#ifdef FOAM_USE_TBB
  #include <tbb/parallel_for.h>
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::syncTools::swapBoundaryCellPositions
(
    const polyMesh& mesh,
    const UList<point>& cellData,
    List<point>& neighbourCellData
)
{
    if (cellData.size() != mesh.nCells())
    {
        FatalErrorInFunction
            << "Number of cell values " << cellData.size()
            << " is not equal to the number of cells in the mesh "
            << mesh.nCells() << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label nBnd = mesh.nFaces()-mesh.nInternalFaces();

    neighbourCellData.setSize(nBnd);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        if (!isA<indirectPolyPatch>(pp))
        {
            const labelUList& faceCells = pp.faceCells();
            forAll(faceCells, i)
            {
                label bFacei = pp.localToGlobal(i) - mesh.nInternalFaces();
                neighbourCellData[bFacei] = cellData[faceCells[i]];
            }
        }
    }
    syncTools::swapBoundaryFacePositions(mesh, neighbourCellData);
}


Foam::PackedBoolList Foam::syncTools::getMasterPoints(const polyMesh& mesh)
{
    PackedBoolList isMasterPoint(mesh.nPoints(), true);

    const globalMeshData& globalData = mesh.globalData();
    const labelList& meshPoints = globalData.coupledPatch().meshPoints();
    const labelListList& slaves = globalData.globalPointSlaves();
    const labelListList& transformedSlaves =
            globalData.globalPointTransformedSlaves();

    forAll(meshPoints, coupledPointi)
    {
        label meshPointi = meshPoints[coupledPointi];
        isMasterPoint[meshPointi] = (
            slaves[coupledPointi].size()
          + transformedSlaves[coupledPointi].size()
        ) > 0;
    }

    return isMasterPoint;
}


Foam::PackedBoolList Foam::syncTools::getMasterEdges(const polyMesh& mesh)
{
    PackedBoolList isMasterEdge(mesh.nEdges(), true);

    const globalMeshData& globalData = mesh.globalData();
    const labelList& meshEdges = globalData.coupledPatchMeshEdges();
    const labelListList& slaves = globalData.globalEdgeSlaves();
    const labelListList& transformedSlaves =
        globalData.globalEdgeTransformedSlaves();

    forAll(meshEdges, coupledEdgeI)
    {
        label meshEdgeI = meshEdges[coupledEdgeI];
        isMasterEdge[meshEdgeI] = (
              slaves[coupledEdgeI].size()
            + transformedSlaves[coupledEdgeI].size()
        ) > 0;
    }

    return isMasterEdge;
}


Foam::PackedBoolList Foam::syncTools::getMasterFaces
(
    const polyMesh& mesh,
    const bool threaded /*=false*/
)
{
    PackedBoolList isMasterFace(mesh.nFaces(), 1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        if (patches[patchi].coupled())
        {
            const coupledPolyPatch& pp =
                refCast<const coupledPolyPatch>(patches[patchi]);

            if (!pp.owner())
            {
                if (threaded && pp.size() > 2048)
                {
#ifndef FOAM_USE_TBB
                    WarningInFunction
                        << "FOAM not linked against TBB! Cannot run multithreaded."
                        << endl;
#else
                    const auto grainSize =
                        std::max<size_t>
                        (
                            pp.size() / tbb::this_task_arena::max_concurrency(),
                            1024
                        );

                    tbb::parallel_for
                    (
                        tbb::blocked_range<size_t>(0, pp.size(), grainSize),
                        [&](const tbb::blocked_range<size_t>& r)
                        {
                            for (size_t i = r.begin(); i < r.end(); ++i)
                            {
                                isMasterFace.unset(pp.start()+i);
                            }
                        },
                        tbb::simple_partitioner()
                    );
#endif
                }
                else
                {
                    forAll(pp, i)
                    {
                        isMasterFace.unset(pp.start()+i);
                    }
                }
            }
        }
    }

    return isMasterFace;
}


Foam::PackedBoolList Foam::syncTools::getInternalOrMasterFaces
(
    const polyMesh& mesh
)
{
    PackedBoolList isMasterFace(mesh.nFaces(), 1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            if (!refCast<const coupledPolyPatch>(pp).owner())
            {
                forAll(pp, i)
                {
                    isMasterFace.unset(pp.start()+i);
                }
            }
        }
        else
        {
            forAll(pp, i)
            {
                isMasterFace.unset(pp.start()+i);
            }
        }
    }

    return isMasterFace;
}


Foam::PackedBoolList Foam::syncTools::getInternalOrCoupledFaces
(
    const polyMesh& mesh
)
{
    PackedBoolList isMasterFace(mesh.nFaces(), 1);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (!pp.coupled())
        {
            forAll(pp, i)
            {
                isMasterFace.unset(pp.start()+i);
            }
        }
    }

    return isMasterFace;
}


// ************************************************************************* //
