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
    (c) 2014-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "layerManipulate/layerManipulate.H"
#include "externalDisplacementMeshMover/fieldSmoother/fieldSmoother.H"
#include "algorithms/PointEdgeWave/PointEdgeWave.H"
#include "sets/topoSets/pointSet.H"
#include "sets/topoSets/faceSet.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "hessianMeshOptimization/leastSquaresCurveFit/leastSquaresCurveFit.H"
#include "meshes/polyMesh/syncTools/dummyTransform.H"
#include "edgeClassification/edgeClassification.H"
#include "meshes/polyMesh/polyMeshCheck/polyMeshTools.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "motionSmoother/polyMeshGeometry/polyMeshGeometry.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

namespace Foam
{

defineTypeNameAndDebug(layerManipulate, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::layerManipulate::makeFaceZonePatch()
{
    labelList addressing;

    // Count faces.
    label nFaces = 0;
    forAll(mesh_.faceZones(), zoneID)
    {
        nFaces += mesh_.faceZones()[zoneID].size();
    }
    addressing.setSize(nFaces);

    nFaces = 0;
    forAll(mesh_.faceZones(), zoneID)
    {
        const labelList& zoneFaces = mesh_.faceZones()[zoneID];
        forAll(zoneFaces, i)
        {
            label faceI = zoneFaces[i];
            addressing[nFaces++] = faceI;
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh_.faces(), addressing),
            mesh_.points()
        )
    );
}


Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::layerManipulate::makeLayerPatch()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& patchToNLayers = layerParams_.numLayers();

    DynamicList<label> adaptPatches(patches.size());

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && patchToNLayers[patchI] != -1)
        {
            adaptPatches.append(patchI);
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        meshRefinement::makePatch
        (
            mesh_,
            adaptPatches
        )
    );
}

Foam::autoPtr<Foam::indirectPrimitivePatch>
Foam::layerManipulate::makeGrownUpPatch(const bool onlySnap)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& patchToNLayers = layerParams_.numLayers();
    const List<Switch>& reSnap = layerParams_.reSnap();

    DynamicList<label> grownUpPatches(patches.size());
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && patchToNLayers[patchI] == -1)
        {
            if (onlySnap)
            {
                if (reSnap[patchI])
                {
                    grownUpPatches.append(patchI);
                }
            }
            else
            {
                grownUpPatches.append(patchI);
            }
        }
    }

    return autoPtr<indirectPrimitivePatch>
    (
        meshRefinement::makePatch
        (
            mesh_,
            grownUpPatches
        )
    );
}

// Tries and find a medial axis point. Done by comparing vectors to nearest
// wall point for both vertices of edge.
bool Foam::layerManipulate::isMaxEdge
(
    const fvMesh& mesh,
    const List<pointData>& pointWallDist,
    const label edgeI,
    const scalar minCos
)
{
    int dummyTrackData = 0;
    const edge& e = mesh.edges()[edgeI];
    if
    (
        !pointWallDist[e[0]].valid(dummyTrackData)
        || !pointWallDist[e[1]].valid(dummyTrackData)
    )
    {
        return false;
    }

    const pointField& points = mesh.points();

    // Do not mark edges with one side on moving wall.

    vector v0(points[e[0]] - pointWallDist[e[0]].origin());
    scalar magV0(mag(v0));

    if (magV0 < SMALL)
    {
        return false;
    }

    vector v1(points[e[1]] - pointWallDist[e[1]].origin());
    scalar magV1(mag(v1));

    if (magV1 < SMALL)
    {
        return false;
    }

    //- Detect based on extrusion vector differing for both endpoints
    //  the idea is that e.g. a sawtooth wall can still be extruded
    //  successfully as long as it is done all to the same direction.
    if ((pointWallDist[e[0]].v() & pointWallDist[e[1]].v()) < minCos)
    {
        vector origVec =
            pointWallDist[e[1]].origin() - pointWallDist[e[0]].origin();
        scalar magOrigVec = mag(origVec);

        if (magOrigVec > SMALL)
        {
            origVec /= magOrigVec;
            if
            (
                (origVec & pointWallDist[e[1]].v()) > 0
                && (origVec & pointWallDist[e[0]].v()) < 0
             )
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


void Foam::layerManipulate::calculateLayerStacks
(
    const indirectPrimitivePatch& pp,
    const boolList& boundaryFaces,
    const boolList& boundaryEdges,
    labelList& cellLayerNumber,
    labelList& layerFaceType,
    labelList& layerPointType,
    labelList& stackEdgeLayer
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& owners = mesh_.faceOwner();

    const volScalarField& layerCells =
        mesh_.lookupObject<volScalarField>("layerStacks");

    forAll(mesh_.cells(), cellI)
    {
        if (layerCells[cellI] == -1)
        {
            cellLayerNumber[cellI] = -1;
        }
    }

    boolList bPoints(mesh_.nPoints(), false);
    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        if (layerCells[owners[meshFaceI]] < 0)
        {
            continue;
        }

        layerFaceType[meshFaceI] = 2;

        label cellI = owners[meshFaceI];

        if (cellLayerNumber[cellI] == 0)
        {
            cellLayerNumber[cellI]++;
        }

        face f = pp[i];
        forAll(f, fp)
        {
            label meshPointI = f[fp];
            bPoints[meshPointI] = true;
            const labelList& pEdges = mesh_.pointEdges()[meshPointI];

            forAll(pEdges, pEI)
            {
                labelHashSet cEdges(mesh_.cellEdges()[cellI]);
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if
                    (
                        !boundaryEdges[edgeI] && stackEdgeLayer[edgeI] == 0
                        && cEdges.found(edgeI)
                     )
                    {
                        stackEdgeLayer[edgeI]++;
                    }
                }
            }
        }
        const labelList& fEdges = mesh_.faceEdges()[meshFaceI];
        labelHashSet cFaces(mesh_.cells()[cellI]);
        forAll(fEdges, fEI)
        {
            label edgeI = fEdges[fEI];

            const labelList& eFaces = mesh_.edgeFaces()[edgeI];
            forAll(eFaces, eFI)
            {
                label faceI = eFaces[eFI];
                if
                (
                    faceI != meshFaceI && layerFaceType[faceI] == -1
                    && !boundaryFaces[faceI] && cFaces.found(faceI)
                )
                {
                    layerFaceType[faceI] = 1;
                }
            }
        }

        forAll(f, fp)
        {
            label meshPointI = f[fp];
            const labelList& pFaces = mesh_.pointFaces()[meshPointI];
            forAll(pFaces, pFI)
            {
                label faceI = pFaces[pFI];
                if
                (
                    layerFaceType[faceI] == -1
                    && !boundaryFaces[faceI] && cFaces.found(faceI)
                 )
                {
                    layerFaceType[faceI] = 1;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        bPoints,
        orEqOp<bool>(),
        false
    );

    forAll(bPoints, meshPointI)
    {
        if (!bPoints[meshPointI])
        {
            continue;
        }

        const labelList& pCells = mesh_.pointCells()[meshPointI];
        forAll(pCells, pCI)
        {
            label cellI = pCells[pCI];

            if (cellLayerNumber[cellI] == 0)
            {
                cellLayerNumber[cellI]++;
            }

            const labelList& pEdges = mesh_.pointEdges()[meshPointI];

            labelHashSet cEdges(mesh_.cellEdges()[cellI]);
            forAll(pEdges, pEI)
            {
                label edgeI = pEdges[pEI];
                if
                (
                    !boundaryEdges[edgeI] && stackEdgeLayer[edgeI] == 0
                    && cEdges.found(edgeI)
                 )
                {
                    stackEdgeLayer[edgeI]++;
                }
            }
            const labelList& pFaces = mesh_.pointFaces()[meshPointI];
            labelHashSet cFaces(mesh_.cells()[cellI]);

            forAll(pFaces, pFI)
            {
                label faceI = pFaces[pFI];
                if
                (
                    layerFaceType[faceI] == -1
                    &&!boundaryFaces[faceI] && cFaces.found(faceI)
                 )
                {
                    layerFaceType[faceI] = 1;
                }
            }
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncEdgeList
        (
            mesh_,
            stackEdgeLayer,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncFaceList
        (
            mesh_,
            layerFaceType,
            maxEqOp<label>()
        );

        labelList neiLayerCells(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiCellLayerNumber(mesh_.nFaces()-mesh_.nInternalFaces());
        for
        (
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            neiLayerCells[faceI-mesh_.nInternalFaces()] =
                layerCells[owners[faceI]];
            neiCellLayerNumber[faceI-mesh_.nInternalFaces()] =
                cellLayerNumber[owners[faceI]];
        }
        syncTools::swapBoundaryFaceList(mesh_, neiLayerCells);
        syncTools::swapBoundaryFaceList(mesh_, neiCellLayerNumber);

        forAll(mesh_.faces(), meshFaceI)
        {
            label own = owners[meshFaceI];

            label patchI = patches.whichPatch(meshFaceI);

            if (layerFaceType[meshFaceI] != -1)
            {
                continue;
            }

            if (patchI == -1 || patches[patchI].coupled())
            {
                label neiCLN = -1;
                label nLayerCell = -1;
                label nei = -1;

                if (patchI == -1)
                {
                    nei = mesh_.faceNeighbour()[meshFaceI];
                    neiCLN = cellLayerNumber[nei];
                    nLayerCell = layerCells[nei];
                }
                else
                {
                    neiCLN =
                        neiCellLayerNumber[meshFaceI-mesh_.nInternalFaces()];
                    nLayerCell =
                        neiLayerCells[meshFaceI-mesh_.nInternalFaces()];
                }

                label cellI = -1;
                label prevLayer = -1;

                if (layerCells[own] == nLayerCell)
                {
                    if (cellLayerNumber[own] > 0 && neiCLN == 0)
                    {
                        if (patchI == -1)
                        {
                            cellI = nei;
                            prevLayer = cellLayerNumber[own];
                        }
                        else
                        {
                            nSet++;
                        }
                    }
                    else if
                    (
                        neiCLN > 0 && cellLayerNumber[own] == 0
                    )
                    {
                        cellI = own;
                        prevLayer = neiCLN;
                    }
                }

                if (cellI != -1)
                {
                    nSet++;
                    cellLayerNumber[cellI] = prevLayer + 1;
                    layerFaceType[meshFaceI] = 2;

                    face f = mesh_.faces()[meshFaceI];
                    const labelList& fEdges = mesh_.faceEdges()[meshFaceI];
                    labelHashSet fEdgeSet(fEdges);

                    forAll(f, fp)
                    {
                        label meshPointI = f[fp];
                        const labelList& pEdges =
                            mesh_.pointEdges()[meshPointI];

                        forAll(pEdges, pEI)
                        {
                            labelHashSet cEdges(mesh_.cellEdges()[cellI]);
                            forAll(pEdges, pEI)
                            {
                                label edgeI = pEdges[pEI];
                                if
                                (
                                    stackEdgeLayer[edgeI] == 0
                                    && cEdges.found(edgeI)
                                    && !fEdgeSet.found(edgeI)
                                )
                                {
                                    stackEdgeLayer[edgeI]++;
                                }
                            }
                        }
                    }

                    labelHashSet cFaces(mesh_.cells()[cellI]);
                    forAll(fEdges, fEI)
                    {
                        label edgeI = fEdges[fEI];

                        const labelList& eFaces = mesh_.edgeFaces()[edgeI];
                        forAll(eFaces, eFI)
                        {
                            label faceI = eFaces[eFI];
                            if
                            (
                                faceI != meshFaceI && layerFaceType[faceI] == -1
                                && cFaces.found(faceI)
                             )
                            {
                                layerFaceType[faceI] = 1;
                            }
                        }
                    }
                    forAll(f, fp)
                    {
                        label meshPointI = f[fp];
                        const labelList& pFaces =
                            mesh_.pointFaces()[meshPointI];
                        forAll(pFaces, pFI)
                        {
                            label faceI = pFaces[pFI];
                            if
                            (
                                layerFaceType[faceI] == -1
                                && cFaces.found(faceI)
                            )
                            {
                                layerFaceType[faceI] = 1;
                            }
                        }
                    }
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    labelList neiLayerCells(mesh_.nFaces()-mesh_.nInternalFaces());
    for
    (
        label faceI = mesh_.nInternalFaces();
        faceI < mesh_.nFaces();
        faceI++
    )
    {
        neiLayerCells[faceI-mesh_.nInternalFaces()] =
            layerCells[owners[faceI]];
    }
    syncTools::swapBoundaryFaceList(mesh_, neiLayerCells);

    //Mark outer stack
    forAll(mesh_.faces(), faceI)
    {
        label patchI = patches.whichPatch(faceI);

        if (patchI == -1 || patches[patchI].coupled())
        {
            label own = owners[faceI];
            label nLayerCell = -1;

            if (patchI == -1)
            {
                label nei = mesh_.faceNeighbour()[faceI];
                nLayerCell = layerCells[nei];
            }
            else
            {
                nLayerCell = neiLayerCells[faceI-mesh_.nInternalFaces()];
            }

            if
            (
                (layerCells[own] != nLayerCell)
                && layerFaceType[faceI] == -1
            )
            {
                layerFaceType[faceI] = 3;
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        layerFaceType,
        maxEqOp<label>()
    );

    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        layerPointType[meshPointI] = 1;
    }

    syncTools::syncPointList
    (
        mesh_,
        layerPointType,
        maxEqOp<label>(),
        label(-1)
    );

    forAll(mesh_.faces(), faceI)
    {
        label patchI = patches.whichPatch(faceI);

        if (patchI == -1 || patches[patchI].coupled())
        {
            label own = owners[faceI];
            label nLayerCell = -1;

            if (patchI == -1)
            {
                label nei = mesh_.faceNeighbour()[faceI];
                nLayerCell = layerCells[nei];
            }
            else
            {
                nLayerCell = neiLayerCells[faceI-mesh_.nInternalFaces()];
            }

            if (layerFaceType[faceI] == 2)
            {
                const face f = mesh_.faces()[faceI];
                forAll(f,fp)
                {
                    if (layerPointType[f[fp]] == -1)
                    {
                        layerPointType[f[fp]] = 0;
                    }
                }
            }
            else if (layerFaceType[faceI] == 3)
            {
                const face f = mesh_.faces()[faceI];
                if (layerCells[own] == -1 || nLayerCell == -1)
                {
                    forAll(f,fp)
                    {
                        layerPointType[f[fp]] = 2;
                    }
                }
            }
        }
    }

    forAll(mesh_.faces(), faceI)
    {
        label patchI = patches.whichPatch(faceI);

        if (patchI == -1 || patches[patchI].coupled())
        {
            label own = owners[faceI];
            label nLayerCell = -1;
            if (patchI == -1)
            {
                label nei = mesh_.faceNeighbour()[faceI];
                nLayerCell = layerCells[nei];
            }
            else
            {
                nLayerCell = neiLayerCells[faceI-mesh_.nInternalFaces()];
            }

            if (layerFaceType[faceI] == 3)
            {
                const face f = mesh_.faces()[faceI];
                if (layerCells[own] != -1 && nLayerCell != -1)
                {
                    forAll(f,fp)
                    {
                        layerPointType[f[fp]] = 3;
                    }
                }
            }
        }
    }
    syncTools::syncPointList(mesh_, layerPointType, maxEqOp<label>(), label(-1));


    forAll(mesh_.cells(), cellI)
    {
        if (layerCells[cellI] > -1)
        {
            const labelList& cellPoints = mesh_.cellPoints()[cellI];
            forAll(cellPoints, cPtI)
            {
                label pointi = cellPoints[cPtI];

                if (layerPointType[pointi] == -1)
                {
                    layerPointType[cellPoints[cPtI]] = 3;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        layerPointType,
        maxEqOp<label>(),
        label(-1)
    );

    if (debug)
    {
        Time& runTime = const_cast<Time&>(mesh_.time());
        runTime++;
        mesh_.write();

        volScalarField layerCount
        (
            IOobject
            (
                "layerCount",
                runTime.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );
        forAll(layerCount, cellI)
        {
            layerCount[cellI] = cellLayerNumber[cellI];
        }
        layerCount.write();
    }
}


void Foam::layerManipulate::writeLayerInfo()
{
    Info<<"Writing layer information "<<endl;

    const volScalarField& layerCells =
        mesh_.lookupObject<volScalarField>("layerStacks");

    autoPtr<indirectPrimitivePatch> ppPtr = makeLayerPatch();
    indirectPrimitivePatch& pp = ppPtr();

    const labelList& meshPoints = pp.meshPoints();

    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh_.edges(),
            mesh_.pointEdges()
        )
    );

    boolList boundaryFaces(mesh_.nFaces(), false);
    boolList boundaryEdges(mesh_.nEdges(), false);

    //Boundary point types
    // -1 : not on boundary
    // 0 : grown patch
    // 1 : grown up
    // 2 : grown up feature pt
    // 3 : stationary corner points
    labelList boundaryPoints(mesh_.nPoints(), -1);

    forAll(pp, i)
    {
        boundaryFaces[pp.addressing()[i]] = true;
    }

    forAll(meshPoints, ptI)
    {
        boundaryPoints[meshPoints[ptI]] = 0;
    }

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh_));
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    forAll(pp.edges(), edgeI)
    {
        boundaryEdges[meshEdges[edgeI]] = true;
    }

    syncTools::syncEdgeList
    (
        mesh_,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh_,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    labelList cellLayerNumber(mesh_.nCells(), 0);
    labelList stackEdgeLayer(mesh_.nEdges(), 0);

    //layer face type
    // 1 : side faces
    // 2 : inner stack faces
    // 3 : outer stack faces
    labelList layerFaceType(mesh_.nFaces(), -1);
    labelList layerPointType(mesh_.nPoints(), -1);

    calculateLayerStacks
    (
        pp,
        boundaryFaces,
        boundaryEdges,
        cellLayerNumber,
        layerFaceType,
        layerPointType,
        stackEdgeLayer
    );

    boolList stackEdges(mesh_.nEdges(), false);
    forAll(mesh_.edges(), edgeI)
    {
        if (stackEdgeLayer[edgeI] != 0)
        {
            stackEdges[edgeI] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        stackEdges,
        orEqOp<bool>(),
        false
    );

    boolList excludedFaces(pp.size(), false);
    edgeClassification eClass
    (
        mesh_,
        mesh_.points(),
        pp,
        meshEdges,
        excludedFaces,
        0.707,
        0.707
    );
    pointField pointNormals = eClass.calculatePointNormals
    (
        excludedFaces,
        0,
        true
    );

    labelList layerCount(mesh_.nPoints(), -1);
    pointField pointOrigin(mesh_.nPoints(), vector::zero);

    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        layerCount[meshPointI] = 0;
    }

    boolList edgeSet(mesh_.nEdges(), false);
    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh_,
            layerCount,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncEdgeList
        (
            mesh_,
            edgeSet,
            orEqOp<bool>(),
            false              // null value
        );

        forAll(mesh_.edges(), edgeI)
        {
            if (!isMasterEdge[edgeI])
            {
                continue;
            }

            if (!boundaryEdges[edgeI] && stackEdges[edgeI] && !edgeSet[edgeI])
            {
                edge e = mesh_.edges()[edgeI];

                if (layerCount[e[0]] != -1)
                {
                    if (layerPointType[e[0]] < 2)
                    {
                        if (layerCount[e[1]] == -1)
                        {
                            layerCount[e[1]] = layerCount[e[0]] + 1;
                        }
                        edgeSet[edgeI] = true;
                        nSet++;
                    }
                }

                if (!edgeSet[edgeI] && layerCount[e[1]] != -1)
                {
                    if (layerPointType[e[1]] < 2)
                    {
                        if (layerCount[e[0]] == -1)
                        {
                            layerCount[e[0]] = layerCount[e[1]] + 1;
                        }
                        edgeSet[edgeI] = true;
                        nSet++;
                    }
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    const labelList& owners = mesh_.faceOwner();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList totalNumLayers(mesh_.nPoints(), -1);
    forAll(mesh_.faces(), facei)
    {
        if (layerFaceType[facei] == 3)
        {
            label own = owners[facei];
            label layerCount = -1;

            if (layerCells[own] != -1)
            {
                layerCount = cellLayerNumber[own];
            }
            else
            {
                label patchI = patches.whichPatch(facei);
                if (patchI == -1)
                {
                    label nei = mesh_.faceNeighbour()[facei];
                    if (layerCells[nei] != -1)
                    {
                        layerCount = cellLayerNumber[nei];
                    }
                }
            }

            if (layerCount != -1)
            {
                face f = mesh_.faces()[facei];
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    totalNumLayers[pointi] = max
                    (
                        totalNumLayers[pointi],layerCount
                    );
                }
            }
        }
    }

    pointField outerPt(mesh_.nPoints(), vector::zero);
    forAll(mesh_.points(), pointi)
    {
        if (layerPointType[pointi] >= 2)
        {
            outerPt[pointi] = mesh_.points()[pointi];
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh_,
            totalNumLayers,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh_,
            outerPt,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        forAll(mesh_.edges(), edgeI)
        {
            if (!boundaryEdges[edgeI] && stackEdges[edgeI])
            {
                edge e = mesh_.edges()[edgeI];

                if
                (
                    totalNumLayers[e[0]] != -1 && totalNumLayers[e[1]] == -1
                    && boundaryPoints[e[0]] != 0
                )
                {
                    totalNumLayers[e[1]] = totalNumLayers[e[0]];
                    outerPt[e[1]] = outerPt[e[0]];
                    nSet++;
                }
                else if
                (
                    totalNumLayers[e[1]] != -1 && totalNumLayers[e[0]] == -1
                    && boundaryPoints[e[1]] != 0
                )
                {
                    totalNumLayers[e[0]] = totalNumLayers[e[1]];
                    outerPt[e[0]] = outerPt[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    labelList patchNumLayers(pp.meshPoints().size(), 0);
    labelList patchLevel(pp.meshPoints().size(), 0);
    scalarField patchFCH(pp.meshPoints().size(), 0);
    scalarField patchLayerHeight(pp.meshPoints().size(), 0);

    forAll(pp.meshPoints(), pointi)
    {
        label meshPointI = pp.meshPoints()[pointi];
        patchNumLayers[pointi] = totalNumLayers[meshPointI];
        patchLevel[pointi] = pointLevel_[meshPointI];

        const labelList& pEdges = mesh_.pointEdges()[meshPointI];

        scalar firstCellHeight = GREAT;
        forAll(pEdges, pEI)
        {
            label edgeI = pEdges[pEI];

            if (!boundaryEdges[edgeI] && stackEdges[edgeI])
            {
                vector eVec = mesh_.edges()[edgeI].vec(mesh_.points());
                scalar projEdgeDist = mag(eVec & pointNormals[pointi]);
                firstCellHeight = min(firstCellHeight, projEdgeDist);
            }
        }
        patchFCH[pointi] = firstCellHeight;

        vector layerVec = outerPt[meshPointI] - mesh_.points()[meshPointI];
        patchLayerHeight[pointi] = mag(layerVec & pointNormals[pointi]);
    }

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        patchLayerHeight,
        minEqOp<scalar>(),
        scalar(0.)          // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        patchFCH,
        minEqOp<scalar>(),
        scalar(0.)          // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        patchFCH,
        minEqOp<scalar>(),
        scalar(0.)          // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        patchLevel,
        maxEqOp<label>(),
        label(0)          // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        patchNumLayers,
        maxEqOp<label>(),
        label(0)          // null value
    );

    //Convert to face data for output
    labelList facePatchNumLayers(pp.size(), labelMax);
    labelList facePatchLevel(pp.size(), 0);
    scalarField facePatchFCH(pp.size(), 0);
    scalarField facePatchLayerHeight(pp.size(), 0);

    label nPatches = patches.size();
    scalarField avePatchNumLayer(nPatches, scalar(0));
    scalarField avePatchFCH(nPatches, scalar(0));
    scalarField avePatchThickness(nPatches, scalar(0));

    forAll(pp, facei)
    {
        label own = owners[pp.addressing()[facei]];
        if (layerCells[own] < 0)
        {
           continue;
        }

        const face& f = pp.localFaces()[facei];

        forAll(f, fp)
        {
            facePatchFCH[facei] += patchFCH[f[fp]];
            facePatchLayerHeight[facei] += patchLayerHeight[f[fp]];
            facePatchNumLayers[facei] = min
                (patchNumLayers[f[fp]],facePatchNumLayers[facei]);
            facePatchLevel[facei] = max
                (patchLevel[f[fp]],facePatchLevel[facei]);
        }
        facePatchFCH[facei] /= f.size();
        facePatchLayerHeight[facei] /= f.size();

        label patchi = patches.whichPatch(pp.addressing()[facei]);
        if (patchi != -1)
        {
            avePatchNumLayer[patchi] += facePatchNumLayers[facei];
            avePatchFCH[patchi] += facePatchFCH[facei];
            avePatchThickness[patchi] += facePatchLayerHeight[facei];
        }
    }

    boolList validPatches(patches.size(), false);
    forAll(patches, patchi)
    {
        if (!patches[patchi].coupled())
        {
            label sz = patches[patchi].size();
            reduce(sz, sumOp<label>());
            if (sz > 0)
            {
                scalar avgNumLayers =
                    returnReduce(avePatchNumLayer[patchi], sumOp<scalar>())/sz;
                if (avgNumLayers > 0)
                {
                    validPatches[patchi] = true;
                }
            }
        }
    }

    // Find maximum length of a patch name, for a nicer output
    label maxPatchNameLen = 0;
    forAll(patches, patchi)
    {
        if (validPatches[patchi])
        {
            word patchName = patches[patchi].name();
            maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
        }
    }

    Info<< nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
        << setw(0) << " faces    layers avg thickness[m]" << nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << " "
        << setw(0) << "                 near-wall overall" << nl
        << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
        << setw(0) << " -----    ------ --------- -------" << endl;

    forAll(patches, patchi)
    {
        if (validPatches[patchi])
        {
            label sz = patches[patchi].size();
            reduce(sz, sumOp<label>());

            if (sz > 0)
            {
                scalar avgThickness =
                    returnReduce(avePatchThickness[patchi], sumOp<scalar>())/sz;
                scalar avgNearWallThickness =
                    returnReduce(avePatchFCH[patchi], sumOp<scalar>())/sz;
                scalar avgNumLayers =
                    returnReduce(avePatchNumLayer[patchi], sumOp<scalar>())/sz;
                Info<< setf(ios_base::left) << setw(maxPatchNameLen)
                    << patches[patchi].name() << setprecision(3)
                    << " " << setw(8) << sz
                    << " " << setw(6) << avgNumLayers
                    << " " << setw(8) << avgNearWallThickness
                    << "  " << setw(8) << avgThickness
                    << endl;
            }
        }
    }

    if (layerParams_.writeVTK())
    {
        simpleVTKWriter layerVTK
        (
            pp.localFaces(),
            pp.localPoints()
        );
        layerVTK.addFaceData("numLayers", facePatchNumLayers);
        layerVTK.addFaceData("fch", facePatchFCH);
        layerVTK.addFaceData("level", facePatchLevel);
        layerVTK.addFaceData("layerHeight", facePatchLayerHeight);

        fileName linfoPath("layerInfo"/mesh_.name());
        if (Pstream::master())
        {
            if (!isDir(linfoPath))
            {
                mkDir(linfoPath);
            }
        }
        layerVTK.write(linfoPath/"layerInfo.vtk");
    }
}


void Foam::layerManipulate::setLayerMethod
(
    const indirectPrimitivePatch& pp,
    const labelList& meshPoints,
    const labelList& ppPtLevel,
    const labelList& ppPtLevelMin,
    List<Tuple2<word,scalar>>& layerMethod,
    scalarField& ppMaxLayerThickness
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const List<wordList>& layerSpec = layerParams_.layerSpec();
    const scalarField& patchFCH = layerParams_.fch();
    const scalarField& patchExpansionRatio = layerParams_.expansionRatio();
    const scalarField& patchMaxLayerThickness =
        layerParams_.maxLayerThickness();
    const scalarField& patchFinalLayerThickness =
        layerParams_.finalLayerThickness();
    const labelList& patchNumLayers = layerParams_.numLayers();
    const List<Switch>& relativeSizes = layerParams_.relativeSizes();

    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        label patchI = patches.whichPatch(meshFaceI);

        face f = pp.localFaces()[i];

        word method = layerSpec[patchI][1];
        if (method == "fch" || method == "rfch")
        {
            bool relative(method == "rfch" ? true : false);
            scalar relativeSclaing = 1;
            forAll(f,fp)
            {
                if (relative)
                {
                    relativeSclaing = edge0Length_ / (1<<ppPtLevel[f[fp]]);
                }
                layerMethod[f[fp]] = Tuple2<word,scalar>
                    ("fch",relativeSclaing*patchFCH[patchI]);
            }
        }
        else if (method == "expansionRatio")
        {
            forAll(f,fp)
            {
                layerMethod[f[fp]] = Tuple2<word,scalar>
                    ("expansionRatio",patchExpansionRatio[patchI]);
            }
        }
        else
        {
            WarningInFunction
                << "Cannot find correct layer method for patch "
                << patches[patchI].name()
                << " requires fch, rfch or expansionRatio to be set"
                << endl;
        }

        word optionalConstraint = layerSpec[patchI][2];
        Switch relSizePatch = relativeSizes[patchI];
        if (optionalConstraint == "maxLayerThickness")
        {
            scalar plt = patchMaxLayerThickness[patchI];
            if (relSizePatch)
            {
                plt *= edge0Length_;
            }

            forAll(f,fp)
            {
                label pti = f[fp];
                scalar lLen = plt;
                if (relSizePatch)
                {
                    label lowerLev(1<<ppPtLevelMin[pti]);
                    label upperLev(1<<ppPtLevel[pti]);

                    scalar levelWeight = 0.5*((1.0/lowerLev)+(1.0/upperLev));
                    lLen *= levelWeight;
                }
                ppMaxLayerThickness[pti] = min
                (
                    ppMaxLayerThickness[pti],
                    lLen
                );
            }
        }
        else if (optionalConstraint == "finalLayerThickness")
        {
            if (method == "fch" || method == "rfch")
            {
                scalar pflt = patchFinalLayerThickness[patchI];
                scalar pfch = patchFCH[patchI];
                label pnl = patchNumLayers[patchI];

                bool relative(method == "rfch" ? true : false);
                forAll(f,fp)
                {
                    label pti = f[fp];
                    scalar ptFCH = pfch;
                    scalar ptFLT = pflt;

                    label lowerLev(1<<ppPtLevelMin[pti]);
                    label upperLev(1<<ppPtLevel[pti]);
                    scalar levelWeight = 0.5*((1.0/lowerLev)+(1.0/upperLev));

                    if (relative)
                    {
                        ptFCH *= edge0Length_*levelWeight;
                    }
                    if (relSizePatch)
                    {
                        ptFLT *= edge0Length_*levelWeight;
                    }

                    scalar str = 1.0;
                    if (pnl > 1)
                    {
                        str = pow((ptFLT/ptFCH),(1.0/(pnl-1)));
                    }
                    scalar lLen = pnl*ptFCH;
                    if (str < 1.0 - SMALL || str > 1.0 + SMALL)
                    {
                        lLen = ptFCH*(1.-pow(str,pnl))
                            / (1.-str);
                    }
                    ppMaxLayerThickness[pti] = min
                    (
                        ppMaxLayerThickness[pti],
                        lLen
                    );
                }
            }
            else
            {
                scalar pflt = patchFinalLayerThickness[patchI];
                label pnl = patchNumLayers[patchI];
                scalar invstr = 1.0/patchExpansionRatio[patchI];

                forAll(f,fp)
                {
                    label pti = f[fp];
                    scalar ptFLT = pflt;

                    label lowerLev(1<<ppPtLevelMin[pti]);
                    label upperLev(1<<ppPtLevel[pti]);
                    scalar levelWeight = 0.5*((1.0/lowerLev)+(1.0/upperLev));

                    if (relSizePatch)
                    {
                        ptFLT *= edge0Length_*levelWeight;
                    }
                    scalar lLen = pnl*ptFLT;
                    if (invstr < 1.0 - SMALL || invstr > 1.0 + SMALL)
                    {
                        lLen = ptFLT*(1.-pow(invstr,pnl))
                            / (1.-invstr);
                    }
                    ppMaxLayerThickness[pti] = min
                    (
                        ppMaxLayerThickness[pti],
                        lLen
                    );
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        ppMaxLayerThickness,
        minEqOp<scalar>(),
        GREAT         // null value
    );

    return;
}


void Foam::layerManipulate::updateGeometry
(
    const pointField& newPoints,
    pointField& newFaceCentres,
    vectorField& newFaceAreas,
    pointField& newCellCentres,
    scalarField& newCellVolumes
)
{
    newFaceCentres = vector::zero;
    newFaceAreas = vector::zero;
    newCellCentres = vector::zero;
    newCellVolumes = 0;

    scalarField newMagFaceAreas(mesh_.nFaces(), Zero);

    if (layerParams_.fastGeomUpdate())
    {
        scalar a1 = (1.0/3.0);
        point avePt = Zero;
        label nextFp = 0;
        vector sumN = Zero;
        scalar sumA = 0.0;
        vector sumAc = Zero;

        vector n;
        vector c;
        scalar a;

        const faceList& fs = mesh_.faces();

        forAll(fs, facei)
        {
            const face& f = fs[facei];
            label nPoints = f.size();
            if (nPoints == 3)
            {
                newFaceCentres[facei] = (1.0/3.0)*
                (
                    newPoints[f[0]] + newPoints[f[1]] + newPoints[f[2]]
                );
                newFaceAreas[facei] = 0.5*
                (
                    (newPoints[f[1]] - newPoints[f[0]])
                    ^(newPoints[f[2]] - newPoints[f[0]])
                );
            }
            else
            {
                sumN = Zero;
                sumA = 0.0;
                sumAc = Zero;
                const point& pt0 = newPoints[f[0]];
                nextFp = 1;
                for (label fp = 0; fp < nPoints-2; fp++)
                {
                    const point& pt1 = newPoints[f[nextFp]];
                    nextFp = f.fcIndex(nextFp);
                    const point& pt2 = newPoints[f[nextFp]];
                    n = (pt1 - pt0)^(pt2 -pt0);
                    c = pt0 + pt1 +pt2;
                    a = mag(n);
                    sumA += a;
                    sumAc += a*c;
                    sumN += n;
                }
                if (sumA > SMALL)
                {
                    newFaceCentres[facei] = a1*sumAc/sumA;
                }
                else
                {
                    avePt = Zero;
                    forAll(f,fp)
                    {
                        avePt += newPoints[f[fp]];
                    }
                    newFaceCentres[facei] = (avePt/nPoints);
                }
                newFaceAreas[facei] = 0.5*sumN;
            }
        }

        const labelList& owners = mesh_.faceOwner();
        const labelList& neighbours = mesh_.faceNeighbour();
        pointField cEst(mesh_.nCells(), Zero);
        forAll(owners, facei)
        {
            label own = owners[facei];
            cEst[own] += newFaceCentres[facei];
        }
        forAll(neighbours, facei)
        {
            label nei = neighbours[facei];
            cEst[nei] += newFaceCentres[facei];
        }
        forAll(cEst, celli)
        {
            cEst[celli] /= mesh_.cells()[celli].size();
        }

        scalarField cellCentreMag(mesh_.nCells(), 0);
        forAll(owners, facei)
        {
            label own = owners[facei];
            const point& cc = cEst[own];
            const point& fC = newFaceCentres[facei];
            const point& fA = newFaceAreas[facei];
            scalar  pyr3Vol = fA & (fC - cc);
            vector pc = 0.75*fC + 0.25*cc;
            scalar pyr3VolMag = mag(pyr3Vol);
            newCellCentres[own] += pyr3VolMag*pc;
            // Accumulate face-pyramid volume
            newCellVolumes[own] += pyr3Vol;
            cellCentreMag[own] += pyr3VolMag;
         }

        forAll(neighbours, facei)
        {
            label nei = neighbours[facei];
            const point& cc = cEst[nei];
            const point& fC = newFaceCentres[facei];
            const point& fA = newFaceAreas[facei];
            scalar  pyr3Vol = fA & (cc - fC);
            vector pc = 0.75*fC + 0.25*cc;
            scalar pyr3VolMag = mag(pyr3Vol);
            newCellCentres[nei] += pyr3VolMag*pc;
            // Accumulate face-pyramid volume
            newCellVolumes[nei] += pyr3Vol;
            cellCentreMag[nei] += pyr3VolMag;
        }

        forAll(newCellCentres, celli)
        {
            scalar cellVol = cellCentreMag[celli];
            if (cellVol > VSMALL)
            {
                newCellCentres[celli] /= cellVol;
            }
            else
            {
                newCellCentres[celli] = cEst[celli];
            }
        }

        newCellVolumes  *= (1.0/3.0);
    }
    else
    {
        // Geometrical calculations
        mesh_.makeFaceCentresAndAreas
        (
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newMagFaceAreas
        );

        mesh_.makeCellCentresAndVols
        (
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes
        );
    }

    return;
}


void Foam::layerManipulate::deactivateSquishCells
(
    const labelList& layerPointType,
    const labelList& boundaryPoints,
    const pointField& newFaceCentres,
    const vectorField& newFaceAreas,
    const pointField& newCellCentres,
    labelList& stationaryPts
)
{
    const volScalarField& layerCells =
        mesh_.lookupObject<volScalarField>("layerStacks");
    const scalar& squishTol = layerParams_.squishTol();

    forAll(stationaryPts, pointi)
    {
        if (stationaryPts[pointi] == 5)
        {
            stationaryPts[pointi] = -2;
        }
    }

    forAll(mesh_.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            continue;
        }

        label nLayerPts = 0;
        label nZonePts = 0;
        label nBoundaryPts = 0;
        const labelList& cPts = mesh_.cellPoints()[celli];
        forAll(cPts, cPtI)
        {
            label pointi = cPts[cPtI];
            if (layerPointType[pointi] != -1)
            {
                nLayerPts++;
            }
            //zone pts
            if
            (
                stationaryPts[pointi] == -1
                || stationaryPts[pointi] == 1
            )
            {
                nZonePts++;
            }
            //boundary pts
            if (boundaryPoints[pointi] > -1)
            {
                nBoundaryPts++;
            }
        }
        if (nLayerPts > 0 || nZonePts > 0 || nBoundaryPts > 0)
        {
            continue;
        }

        bool squished = false;
        const cell& c = mesh_.cells()[celli];
        scalar squishVal = scalar(0);
        scalar totFaceArea = scalar(0);
        forAll(c, cFI)
        {
            label facei = c[cFI];
            vector fcToCC = newFaceCentres[facei] - newCellCentres[celli];
            scalar magfcToCC = mag(fcToCC);
            vector faceArea = newFaceAreas[facei];
            scalar magFaceArea = mag(faceArea);
            totFaceArea += magFaceArea;
            if (magFaceArea < SMALL || magfcToCC < SMALL)
            {
                squished = true;
            }
            else
            {
                fcToCC /= magfcToCC;
                squishVal += mag(faceArea & fcToCC);
            }
        }
        if (totFaceArea > SMALL)
        {
            squishVal /= totFaceArea;
            if (squished || squishVal < squishTol)
            {
                forAll(cPts, cPtI)
                {
                    label pointi = cPts[cPtI];
                    if (stationaryPts[pointi] == -2)
                    {
                        stationaryPts[pointi] = 5;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    DynamicList<label> nbrCells(mesh_.nCells()/10);
    forAll(mesh_.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            continue;
        }

        label nLayerPts = 0;
        label nStaticPts = 0;
        const labelList& cPts = mesh_.cellPoints()[celli];
        forAll(cPts, cPtI)
        {
            label pointi = cPts[cPtI];

            if (stationaryPts[pointi] == 5)
            {
                nStaticPts++;
            }
            if (layerPointType[pointi] != -1)
            {
                nLayerPts++;
            }
        }
        if (nLayerPts > 0 && nStaticPts > 0)
        {
            nbrCells.append(celli);
        }
    }

    forAll(nbrCells, nbri)
    {
        label celli = nbrCells[nbri];
        const labelList& cPts = mesh_.cellPoints()[celli];
        forAll(cPts, cPtI)
        {
            label pointi = cPts[cPtI];
            if (stationaryPts[pointi] == -2)
            {
                stationaryPts[pointi] = 5;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );
}

void Foam::layerManipulate::smoothLayerStack
(
    const meshControl& controller,
    const label nSmoothIter,
    const dictionary& grownUpGeometryDict,
    const dictionary& grownUpZoneGeometryDict,
    const labelList& layerOffset,
    labelList& adjustedNLayers,
    const bool adjustFinalNLayers
)
{
    Info<<"Smoothing layer cells "<<endl;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const volScalarField& layerCells =
        mesh_.lookupObject<volScalarField>("layerStacks");

    autoPtr<indirectPrimitivePatch> ppPtr = makeLayerPatch();
    indirectPrimitivePatch& pp = ppPtr();

    boolList layerEdges(mesh_.nEdges(), false);
    forAll(mesh_.cells(), cellI)
    {
        if (layerCells[cellI] != -1)
        {
            const labelList& cellEdges = mesh_.cellEdges()[cellI];
            forAll(cellEdges, cEI)
            {
                layerEdges[cellEdges[cEI]] = true;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh_,
        layerEdges,
        orEqOp<bool>(),
        false
    );

    const labelList& meshPoints = pp.meshPoints();
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh_.edges(),
            mesh_.pointEdges()
        )
    );

    //Set layer growth method
    labelList ppPtLevel(meshPoints.size(), -1);
    labelList ppPtLevelMin(meshPoints.size(), labelMax);
    forAll(meshPoints, pti)
    {
        const labelList& pFaces = pp.pointFaces()[pti];
        forAll(pFaces, pfi)
        {
            label meshFaceI = pp.addressing()[pFaces[pfi]];
            label own = mesh_.faceOwner()[meshFaceI];
            ppPtLevel[pti] = max(cellLevel_[own],ppPtLevel[pti]);
            ppPtLevelMin[pti] = min(cellLevel_[own],ppPtLevelMin[pti]);
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        ppPtLevel,
        maxEqOp<label>(),
        label(-1)         // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        ppPtLevelMin,
        minEqOp<label>(),
        labelMax         // null value
    );

    List<Tuple2<word,scalar>> layerMethod(meshPoints.size());
    scalarField ppMaxLayerThickness(meshPoints.size(), GREAT);
    setLayerMethod
    (
        pp,
        meshPoints,
        ppPtLevel,
        ppPtLevelMin,
        layerMethod,
        ppMaxLayerThickness
    );

    boolList boundaryFaces(mesh_.nFaces(), false);
    boolList boundaryEdges(mesh_.nEdges(), false);

    //Boundary point types
    // -1 : not on boundary
    // 0 : grown patch
    // 1 : grown up
    // 2 : grown up feature pt
    // 3 : stationary corner points
    labelList boundaryPoints(mesh_.nPoints(), -1);

    forAll(pp, i)
    {
        boundaryFaces[pp.addressing()[i]] = true;
    }

    forAll(meshPoints, ptI)
    {
        boundaryPoints[meshPoints[ptI]] = 0;
    }

    autoPtr<indirectPrimitivePatch> grownUpSnapPPPtr = makeGrownUpPatch(true);
    indirectPrimitivePatch& grownUpSnapPP = grownUpSnapPPPtr();

    autoPtr<searchableSurfaces> grownUpGeometryPtr
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                      // dummy name
                mesh_.time().constant(),     // directory
                "triSurface",               // instance
                mesh_.time(),                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
             ),
             grownUpGeometryDict
         )
    );

    autoPtr<searchableSurfaces> grownUpZoneGeometryPtr
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                      // dummy name
                mesh_.time().constant(),     // directory
                "triSurface",               // instance
                mesh_.time(),                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
             ),
             grownUpZoneGeometryDict
         )
    );

    autoPtr<indirectPrimitivePatch> grownUpPPPtr = makeGrownUpPatch();
    indirectPrimitivePatch& grownUpPP = grownUpPPPtr();

    autoPtr<indirectPrimitivePatch> fzonePPPtr = makeFaceZonePatch();
    indirectPrimitivePatch& fzonePP = fzonePPPtr();

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh_));
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));
    boolList featurePts(mesh_.nPoints(), false);
    boolList featureEdges(mesh_.nEdges(), false);

    //Mark zone feature edges
    if
    (
        layerParams_.dualReSnapZones()
        && returnReduce(fzonePP.size(), sumOp<label>()) != 0
    )
    {
        boolList excludedFaces(fzonePP.size(), false);

        const labelList zoneEdges
        (
            fzonePP.meshEdges
            (
                mesh_.edges(),
                mesh_.pointEdges()
            )
        );

        edgeClassification eClass
        (
            mesh_,
            mesh_.points(),
            fzonePP,
            zoneEdges,
            excludedFaces,
            0.8191,
            0.8191,
            -GREAT,
            true //flip face based on owner zone
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();
        forAll(zoneEdges, edgeI)
        {
            if
            (
                eType[edgeI].first() == edgeClassification::CONCAVE
                || eType[edgeI].first() == edgeClassification::CONVEX
            )
            {
                label meshEdgeI = zoneEdges[edgeI];
                edge e = mesh_.edges()[meshEdgeI];
                featureEdges[meshEdgeI] = true;
                featurePts[e[0]] = true;
                featurePts[e[1]] = true;
            }
        }
    }

    //Calculate feature edges grown up patch
    {
        boolList excludedFaces(grownUpPP.size(), false);
        const labelList meshEdgesGrownUp
        (
            grownUpPP.meshEdges
            (
                mesh_.edges(),
                mesh_.pointEdges()
            )
        );

        edgeClassification eClass
        (
            mesh_,
            mesh_.points(),
            grownUpPP,
            meshEdgesGrownUp,
            excludedFaces,
            0.8191,
            0.8191
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        forAll(meshEdgesGrownUp, edgeI)
        {
            if
            (
                eType[edgeI].first() == edgeClassification::CONCAVE
                || eType[edgeI].first() == edgeClassification::CONVEX
            )
            {
                label meshEdgeI = meshEdgesGrownUp[edgeI];
                edge e = mesh_.edges()[meshEdgeI];
                featureEdges[meshEdgeI] = true;
                featurePts[e[0]] = true;
                featurePts[e[1]] = true;
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            featureEdges,
            orEqOp<bool>(),
            false
        );

        syncTools::syncPointList
        (
            mesh_,
            featurePts,
            orEqOp<bool>(),
            false
        );

        labelList nFeatureEdges(mesh_.nPoints(), 0);
        forAll(mesh_.points(), pointi)
        {
            labelList pEdges = mesh_.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                if (isMasterEdge[pEdges[pEI]] && featureEdges[pEdges[pEI]])
                {
                    nFeatureEdges[pointi]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            nFeatureEdges,
            plusEqOp<label>(),
            label(0)
        );

        forAll(mesh_.points(), pointi)
        {
            if (boundaryPoints[pointi] == -1)
            {
                if (featurePts[pointi])
                {
                    if (nFeatureEdges[pointi] > 2)
                    {
                        boundaryPoints[pointi] = 3;
                    }
                    else
                    {
                        boundaryPoints[pointi] = 2;
                    }
                }
            }
        }

        forAll(grownUpPP.meshPoints(), ptI)
        {
            label pointi = grownUpPP.meshPoints()[ptI];

            if (boundaryPoints[pointi] == -1)
            {
                boundaryPoints[pointi] = 1;
            }
        }
    }

    forAll(pp.edges(), edgeI)
    {
        boundaryEdges[meshEdges[edgeI]] = true;
    }

    syncTools::syncEdgeList
    (
        mesh_,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh_,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    labelList cellLayerNumber(mesh_.nCells(), 0);
    labelList stackEdgeLayer(mesh_.nEdges(), 0);

    //layer face type
    // 1 : side faces
    // 2 : inner stack faces
    // 3 : outer stack faces
    labelList layerFaceType(mesh_.nFaces(), -1);
    labelList layerPointType(mesh_.nPoints(), -1);

    calculateLayerStacks
    (
        pp,
        boundaryFaces,
        boundaryEdges,
        cellLayerNumber,
        layerFaceType,
        layerPointType,
        stackEdgeLayer
    );

    autoPtr<indirectPrimitivePatch> outerShellPtr;
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        DynamicList<label> extFaces(mesh_.nFaces()/10);
        forAll(mesh_.faces(), facei)
        {
            if (layerFaceType[facei] == 3)
            {
                extFaces.append(facei);
            }
        }
        outerShellPtr.reset
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>(mesh_.faces(), extFaces),
                mesh_.points()
             )
        );
    }

    //Calculate layer interface weights based on number of layer
    // and non-layers cells in order to reduce jagged layer interface
    scalarField nonLayerWeights(mesh_.nPoints(),0);
    scalarField layerWeights(mesh_.nPoints(),0);
    forAll(mesh_.points(), pointi)
    {
        if (layerPointType[pointi] == 2)
        {
            if (boundaryPoints[pointi] == 1)
            {
                const labelList& pFaces = mesh_.pointFaces()[pointi];
                forAll(pFaces, pFI)
                {
                    label faceI = pFaces[pFI];
                    label patchI = patches.whichPatch(faceI);

                    if (patchI != -1 && !patches[patchI].coupled())
                    {
                        label cellI = mesh_.faceOwner()[faceI];
                        if (layerCells[cellI] != -1)
                        {
                            layerWeights[pointi] += scalar(1);
                        }
                        else
                        {
                            nonLayerWeights[pointi] += scalar(1);
                        }
                    }
                }
            }
            else
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    label cellI = pCells[pCI];
                    if (layerCells[cellI] != -1)
                    {
                        layerWeights[pointi] += scalar(1);
                    }
                    else
                    {
                        nonLayerWeights[pointi] += scalar(1);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        layerWeights,
        plusEqOp<scalar>(),
        scalar(0.)          // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        nonLayerWeights,
        plusEqOp<scalar>(),
        scalar(0.)          // null value
    );

    //Calculate weights at interface of layer and non-layer cells
    scalar nliw = layerParams_.layerInterfaceWeights();
    scalar liw = scalar(1) - nliw;
    forAll(mesh_.points(), pointi)
    {
        if (layerPointType[pointi] == 2)
        {
            scalar nPointCells =  nonLayerWeights[pointi]
                + layerWeights[pointi];

            if (nonLayerWeights[pointi] > 0 && layerWeights[pointi] > 0)
            {
                nonLayerWeights[pointi] = nliw/nonLayerWeights[pointi];
                layerWeights[pointi] = liw/layerWeights[pointi];
            }
            else if (nPointCells > 0)
            {
                nonLayerWeights[pointi] = 1.0 / nPointCells;
                layerWeights[pointi] = 1.0 / nPointCells;
            }
        }
    }

    boolList stackEdges(mesh_.nEdges(), false);
    forAll(mesh_.edges(), edgeI)
    {
        if (stackEdgeLayer[edgeI] != 0)
        {
            stackEdges[edgeI] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        stackEdges,
        orEqOp<bool>(),
        false
    );

    pointField pointNormals,unsmoothedPointNormals;
    labelList ftrPointOrigin(mesh_.nPoints(), -1);
    //Mark feature points
    {
        boolList excludedFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh_,
            mesh_.points(),
            pp,
            meshEdges,
            excludedFaces,
            0.707,//convex
            0.93969,//concave
            -0.707//-0.9848 //baffle
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        forAll(meshEdges, edgeI)
        {
            edge e = mesh_.edges()[meshEdges[edgeI]];
            if (eType[edgeI].first() == edgeClassification::BAFFLE)
            {
                ftrPointOrigin[e[0]] = max(ftrPointOrigin[e[0]],label(2));
                ftrPointOrigin[e[1]] = max(ftrPointOrigin[e[1]],label(2));
            }
            else if
            (
                eType[edgeI].first() == edgeClassification::CONVEX
            )
            {
                ftrPointOrigin[e[0]] = max(ftrPointOrigin[e[0]],label(1));
                ftrPointOrigin[e[1]] = max(ftrPointOrigin[e[1]],label(1));
            }
            else if
            (
                eType[edgeI].first() == edgeClassification::CONCAVE
            )
            {
                ftrPointOrigin[e[0]] = max(ftrPointOrigin[e[0]],label(0));
                ftrPointOrigin[e[1]] = max(ftrPointOrigin[e[1]],label(0));
            }
        }
    }

    //Calculate area weighted point normals (feature points fixed)
    {
        boolList excludedFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh_,
            mesh_.points(),
            pp,
            meshEdges,
            excludedFaces,
            0.707,//-0.5,
            0.707//-0.5
        );
        pointNormals = eClass.calculatePointNormals
        (
            excludedFaces,
            4,
            true
        );
        unsmoothedPointNormals =
            eClass.calculatePointNormals(excludedFaces, 0, true);

    }

    syncTools::syncPointList
    (
        mesh_,
        ftrPointOrigin,
        maxEqOp<label>(),
        label(-1)
    );

    forAll(ftrPointOrigin, pointi)
    {
        if (ftrPointOrigin[pointi] == 2)
        {
            const labelList& pFaces = mesh_.pointFaces()[pointi];
            forAll(pFaces, pFI)
            {
                label faceI = pFaces[pFI];
                if (!boundaryFaces[faceI] && mesh_.faces()[faceI].size() == 3)
                {
                    ftrPointOrigin[pointi] = -1;
                    break;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        ftrPointOrigin,
        minEqOp<label>(),
        label(-1)
    );

    mesh_.clearOut();
    label nMedialAxisIter = mesh_.globalData().nTotalPoints();

    // Distance to wall
    List<pointData> pointWallDist(mesh_.nPoints());
    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;

    // 1. Calculate distance to points where displacement is specified.
    {
        labelList nSurfLayer(mesh_.nPoints(), -1);
        forAll(mesh_.points(), pointi)
        {
            if (layerPointType[pointi] >= 2)
            {
                nSurfLayer[pointi] = 0;
            }
        }

        while (true)
        {
            label nSet = 0;

            syncTools::syncPointList
            (
                mesh_,
                nSurfLayer,
                maxEqOp<label>(),
                label(-1)              // null value
            );

            forAll(mesh_.edges(), edgeI)
            {
                if (!boundaryEdges[edgeI] && stackEdges[edgeI])
                {
                    edge e = mesh_.edges()[edgeI];

                    if
                    (
                        nSurfLayer[e[0]] != -1 && nSurfLayer[e[1]] == -1
                        && boundaryPoints[e[0]] != 0
                    )
                    {
                        nSurfLayer[e[1]] = nSurfLayer[e[0]]+1;
                        nSet++;
                    }
                    else if
                    (
                        nSurfLayer[e[1]] != -1 && nSurfLayer[e[0]] == -1
                        && boundaryPoints[e[1]] != 0
                    )
                    {
                        nSurfLayer[e[0]] = nSurfLayer[e[1]]+1;
                        nSet++;
                    }
                }
            }

            if (returnReduce(nSet, sumOp<label>()) == 0)
            {
                break;
            }
        }

        scalarField estimatedLayerHeight(meshPoints.size(), 0);

        scalarField len(meshPoints.size(), 0);
        labelList nPointFaces(meshPoints.size(), 0);

        forAll(meshPoints, ptI)
        {
            const labelList& pFaces = pp.pointFaces()[ptI];

            forAll(pFaces, pFI)
            {
                label meshFaceI = pp.addressing()[pFaces[pFI]];
                label own = mesh_.faceOwner()[meshFaceI];
                len[ptI] += edge0Length_ / (1<<cellLevel_[own]);
            }
            nPointFaces[ptI] += pFaces.size();
        }

        syncTools::syncPointList
        (
            mesh_,
            meshPoints,
            len,
            plusEqOp<scalar>(),
            scalar(0.)          // null value
        );

        syncTools::syncPointList
        (
            mesh_,
            meshPoints,
            nPointFaces,
            plusEqOp<label>(),
            label(0)          // null value
        );

        //assume stretching for total layer thickness estimate
        // needs to be made exact
        scalar str = 0.769;//0.869;
        forAll(meshPoints, ptI)
        {
            scalar ratio = (1-pow(str, nSurfLayer[meshPoints[ptI]]))
                / (scalar(1.)-str);
            estimatedLayerHeight[ptI] = ratio*(len[ptI]/nPointFaces[ptI]);
        }

        // Seed data.
        DynamicList<pointData> wallInfo(meshPoints.size());
        DynamicList<label> wallPoints(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label meshPointI = meshPoints[patchPointI];
            wallPoints.append(meshPointI);
            wallInfo.append
            (
                pointData
                (
                    mesh_.points()[meshPointI],
                    0.0,
                    estimatedLayerHeight[patchPointI],
                    unsmoothedPointNormals[patchPointI]
                )
            );
        }
        wallInfo.shrink();
        wallPoints.shrink();

        List<pointData> edgeWallDist(mesh_.nEdges());
        // Do all calculations
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh_,
            wallPoints,
            wallInfo,

            pointWallDist,
            edgeWallDist,
            0,
            dummyTrackData
        );
        wallDistCalc.iterate(nMedialAxisIter);
    }

    scalarField medialDist(mesh_.nPoints(), scalar(0));
    scalarField medialLayerHeight(mesh_.nPoints(), scalar(0));

    // 1. Medial axis points
    {
        const edgeList& edges = mesh_.edges();
        List<pointData> pointMedialDist(mesh_.nPoints());
        List<pointData> edgeMedialDist(mesh_.nEdges());

        // Seed point data.
        DynamicList<pointData> maxInfo(meshPoints.size());
        DynamicList<label> maxPoints(meshPoints.size());

        forAll(edges, edgeI)
        {
            if (isMaxEdge(mesh_, pointWallDist, edgeI, -0.5))
            {
                // Both end points of edge have very different nearest wall
                // point. Mark both points as medial axis points.

                // Approximate medial axis location on edge.
                //const point medialAxisPt = e.centre(points);
                const edge& e = edges[edgeI];
                vector eVec = e.vec(mesh_.points());
                scalar eMag = mag(eVec);
                if (eMag > VSMALL)
                {
                    eVec /= eMag;

                    // Calculate distance along edge
                    const point& p0 = mesh_.points()[e[0]];
                    const point& p1 = mesh_.points()[e[1]];
                    scalar dist0 = (p0-pointWallDist[e[0]].origin()) & eVec;
                    scalar dist1 = (pointWallDist[e[1]].origin()-p1) & eVec;
                    scalar s = 0.5*(dist1+eMag+dist0);

                    point medialAxisPt;
                    if (s <= dist0)
                    {
                        medialAxisPt = p0;
                    }
                    else if (s >= dist0+eMag)
                    {
                        medialAxisPt = p1;
                    }
                    else
                    {
                        medialAxisPt = p0+(s-dist0)*eVec;
                    }

                    scalar totalLayerHeight = pointWallDist[e[0]].s() +
                        pointWallDist[e[1]].s();
                    vector origToOrig = pointWallDist[e[1]].origin()
                        - pointWallDist[e[0]].origin();

                    forAll(e, ep)
                    {
                        label pointi = e[ep];
                        maxPoints.append(pointi);
                        maxInfo.append
                        (
                            pointData
                            (
                                medialAxisPt,   //points[pointi],
                                magSqr(mesh_.points()[pointi]-medialAxisPt),
                                totalLayerHeight,
                                origToOrig
                            )
                        );
                    }
                }
            }
        }
        maxInfo.shrink();
        maxPoints.shrink();

        // Do all calculations
        PointEdgeWave<pointData> medialDistCalc
        (
            mesh_,
            maxPoints,
            maxInfo,

            pointMedialDist,
            edgeMedialDist,
            0,
            dummyTrackData
        );
        medialDistCalc.iterate(2*nMedialAxisIter);

        // Extract medial axis distance as pointScalarField
        forAll(meshPoints, ptI)
        {
            label meshPointI = meshPoints[ptI];
            if (pointMedialDist[meshPointI].valid(dummyTrackData))
            {
                scalar pLen = 2*edge0Length_ / (1<<ppPtLevel[ptI]);
                scalar maxDist = max(mag(pointMedialDist[meshPointI].v()),pLen);
                scalar mDist =
                    Foam::sqrt(pointMedialDist[meshPointI].distSqr());
                if (mDist < maxDist)
                {
                    medialDist[meshPointI] = mDist;
                }
                else
                {
                    medialDist[meshPointI] = GREAT;
                }
            }
            else
            {
                medialDist[meshPointI] = GREAT;
            }
        }

        forAll(mesh_.points(), pointi)
        {
            if (pointMedialDist[pointi].valid(dummyTrackData))
            {
                medialLayerHeight[pointi] = pointMedialDist[pointi].s();
            }
            else
            {
                medialLayerHeight[pointi] = GREAT;
            }
        }
    }

    if (debug)
    {
        simpleVTKWriter medialVTK
        (
            pp.localFaces(),
            pp.localPoints()
        );
        scalarField medialHeight(meshPoints.size(), 0);

        forAll(meshPoints, ptI)
        {
            label meshPointI = meshPoints[ptI];
            medialHeight[ptI] = medialDist[meshPointI];
        }
        medialVTK.addPointData("medialHeight", medialHeight);

        medialVTK.write("medialVTK.vtk");
    }

    labelList layerCount(mesh_.nPoints(), -1);

    pointField pointOrigin(mesh_.nPoints(), vector::zero);
    pointField pointSurfNormal(mesh_.nPoints(), vector::zero);

    scalarField maxLayerThickness(mesh_.nPoints(), GREAT);

    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        layerCount[meshPointI] = layerOffset[ptI];
        pointOrigin[meshPointI] = pp.localPoints()[ptI];
        pointSurfNormal[meshPointI] = -pointNormals[ptI];
        maxLayerThickness[meshPointI] = ppMaxLayerThickness[ptI];
    }

    labelList nVisited(mesh_.nPoints(), 0);

    boolList edgeSet(mesh_.nEdges(), false);

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh_,
            ftrPointOrigin,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh_,
            pointOrigin,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            pointSurfNormal,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            layerCount,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh_,
            medialDist,
            maxEqOp<scalar>(),
            scalar(0)
        );

        syncTools::syncEdgeList
        (
            mesh_,
            edgeSet,
            orEqOp<bool>(),
            false              // null value
        );

        syncTools::syncPointList
        (
            mesh_,
            maxLayerThickness,
            minEqOp<scalar>(),
            GREAT
        );

        forAll(mesh_.edges(), edgeI)
        {
            if (!isMasterEdge[edgeI])
            {
                continue;
            }

            if (!boundaryEdges[edgeI] && stackEdges[edgeI] && !edgeSet[edgeI])
            {
                edge e = mesh_.edges()[edgeI];

                if (layerCount[e[0]] != -1)
                {
                    if (layerPointType[e[0]] < 2)
                    {
                        if (layerCount[e[1]] == -1)
                        {
                            layerCount[e[1]] = layerCount[e[0]] + 1;
                        }
                        edgeSet[edgeI] = true;

                        nVisited[e[1]]++;

                        maxLayerThickness[e[1]] = maxLayerThickness[e[0]];
                        ftrPointOrigin[e[1]] = ftrPointOrigin[e[0]];
                        pointOrigin[e[1]] += pointOrigin[e[0]];
                        medialDist[e[1]] += medialDist[e[0]];
                        pointSurfNormal[e[1]] += pointSurfNormal[e[0]];
                        nSet++;
                    }
                }

                if (!edgeSet[edgeI] && layerCount[e[1]] != -1)
                {
                    if (layerPointType[e[1]] < 2)
                    {
                        if (layerCount[e[0]] == -1)
                        {
                            layerCount[e[0]] = layerCount[e[1]] + 1;
                        }
                        edgeSet[edgeI] = true;

                        nVisited[e[0]]++;

                        maxLayerThickness[e[0]] = maxLayerThickness[e[1]];
                        ftrPointOrigin[e[0]] = ftrPointOrigin[e[1]];
                        pointOrigin[e[0]] += pointOrigin[e[1]];
                        medialDist[e[0]] += medialDist[e[1]];
                        pointSurfNormal[e[0]] += pointSurfNormal[e[1]];
                        nSet++;
                    }
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        nVisited,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh_.points(), pointi)
    {
        if (nVisited[pointi] > 1)
        {
            pointOrigin[pointi] /= nVisited[pointi];
            pointSurfNormal[pointi] /= nVisited[pointi];
            medialDist[pointi] /= nVisited[pointi];
        }
    }

    //Check for unitinitilised hanging points
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        forAll(mesh_.points(), pointi)
        {
            if (layerPointType[pointi] == 2 && nVisited[pointi] == 0)
            {
                const labelList& pEdges = mesh_.pointEdges()[pointi];
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    edge e = mesh_.edges()[edgei];
                    label otherPt(e[0] ==  pointi ? e[1] : e[0]);
                    if (layerEdges[edgei] && layerPointType[otherPt] == 2)
                    {
                        ftrPointOrigin[pointi] =
                            max(ftrPointOrigin[otherPt],ftrPointOrigin[pointi]);
                        maxLayerThickness[pointi] = min
                        (
                            maxLayerThickness[otherPt],
                            maxLayerThickness[pointi]
                        );
                        pointOrigin[pointi] += pointOrigin[otherPt];
                        pointSurfNormal[pointi] += pointSurfNormal[otherPt];
                        medialDist[pointi] += medialDist[otherPt];
                        nVisited[pointi]++;
                    }
                }
            }
            else
            {
                nVisited[pointi] = 1;
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            nVisited,
            plusEqOp<label>(),
            label(0)
        );
        syncTools::syncPointList
        (
            mesh_,
            ftrPointOrigin,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh_,
            maxLayerThickness,
            minEqOp<scalar>(),
            GREAT
        );

        syncTools::syncPointList
        (
            mesh_,
            pointOrigin,
            plusEqOp<vector>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            pointSurfNormal,
            plusEqOp<vector>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            medialDist,
            plusEqOp<scalar>(),
            scalar(0)
        );

        forAll(mesh_.points(), pointi)
        {
            if (nVisited[pointi] > 1)
            {
                pointOrigin[pointi] /= nVisited[pointi];
                pointSurfNormal[pointi] /= nVisited[pointi];
                medialDist[pointi] /= nVisited[pointi];
            }
        }
    }

    labelList totalNumLayers(mesh_.nPoints(), -1);
    forAll(mesh_.faces(), facei)
    {
        if (layerFaceType[facei] == 3)
        {
            face f = mesh_.faces()[facei];
            forAll(f,fp)
            {
                label pointi = f[fp];
                totalNumLayers[pointi] = max
                (
                    totalNumLayers[pointi],layerCount[pointi]
                );
            }
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh_,
            totalNumLayers,
            maxEqOp<label>(),
            label(-1)
        );

        forAll(mesh_.edges(), edgeI)
        {
            if (!boundaryEdges[edgeI] && stackEdges[edgeI])
            {
                edge e = mesh_.edges()[edgeI];

                if
                (
                    totalNumLayers[e[0]] != -1 && totalNumLayers[e[1]] == -1
                    && boundaryPoints[e[0]] != 0
                )
                {
                    totalNumLayers[e[1]] = totalNumLayers[e[0]];
                    nSet++;
                }
                else if
                (
                    totalNumLayers[e[1]] != -1 && totalNumLayers[e[0]] == -1
                    && boundaryPoints[e[1]] != 0
                )
                {
                    totalNumLayers[e[0]] = totalNumLayers[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    pointField newPoints = mesh_.points();
    pointField tPts(mesh_.points());

    pointField newCellCentres(mesh_.nCells(), vector::zero);
    scalarField newCellVolumes(mesh_.nCells(), 0);
    pointField newFaceCentres(mesh_.nFaces(), vector::zero);
    vectorField newFaceAreas(mesh_.nFaces(), vector::zero);

    labelList resetPts(mesh_.nPoints(), -1);

    //Prevent faceZones from moving if further than estimated
    // layer depth from surface
    // startionaryPts: -2 (moving) -1 (moving zone),
    // 1 (stationary), 2 (gap), 3 (non layer medial pts)
    labelList stationaryPts(mesh_.nPoints(), -2);
    scalar zoneLayersScaling
    (
        layerParams_.dualReSnapZones()
        ? GREAT : layerParams_.zoneLayersScaling()
    );
    forAll(fzonePP.meshPoints(), pointi)
    {
        label meshPointI = fzonePP.meshPoints()[pointi];
        if (!pointWallDist[meshPointI].valid(dummyTrackData))
        {
            stationaryPts[meshPointI] = -1;
            continue;
        }

        scalar maxLayerHeightSqr = zoneLayersScaling
            *sqr(pointWallDist[meshPointI].s());
        if
        (
            pointWallDist[meshPointI].distSqr() == GREAT
            || pointWallDist[meshPointI].distSqr() >= maxLayerHeightSqr
        )
        {
            stationaryPts[meshPointI] = 1;
        }
        else
        {
            stationaryPts[meshPointI] = -1;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    //Check for cells with all points zoned and set stationary
    forAll(mesh_.cells(), celli)
    {
        const labelList& cPts = mesh_.cellPoints()[celli];
        label nZonePts = 0;
        forAll(cPts, i)
        {
            label meshPointI = cPts[i];
            if
            (
                stationaryPts[meshPointI] == -1
                || stationaryPts[meshPointI] == 1
            )
            {
                nZonePts++;
            }
        }
        if (nZonePts == cPts.size())
        {
            forAll(cPts, i)
            {
                label meshPointI = cPts[i];
                stationaryPts[meshPointI] = 1;
            }
        }
    }

    //Block based on medial axis edge cells
    {
        forAll(mesh_.edges(), edgei)
        {
            if (isMaxEdge(mesh_, pointWallDist, edgei, -0.5))
            {
                const labelList& eCells = mesh_.edgeCells()[edgei];
                forAll(eCells, eCI)
                {
                    label celli = eCells[eCI];
                    if (layerCells[celli] > -1)
                    {
                        continue;
                    }

                    const labelList& cPts = mesh_.cellPoints()[celli];
                    label nLayerPts = 0;

                    forAll(cPts, i)
                    {
                        label pointi = cPts[i];
                        if (layerPointType[pointi] == -1)
                        {
                            nLayerPts++;
                        }
                    }

                    if (nLayerPts == cPts.size())
                    {
                        forAll(cPts, i)
                        {
                            label pointi = cPts[i];
                            if (stationaryPts[pointi] != 1)
                            {
                                stationaryPts[pointi] = 3;
                            }
                        }
                    }
                }
            }
        }
    }

    //Prevent motion of gap cell points
    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh_.lookupObject<volScalarField>("gapCells");

        forAll(mesh_.cells(), celli)
        {
            if (gapCells[celli] > -1)
            {
                const labelList& cPts = mesh_.cellPoints()[celli];
                forAll(cPts, i)
                {
                    stationaryPts[cPts[i]] = 2;
                }
            }
        }
    }

    //Prevent motion at points further than the maximum projection distance
    const scalar maxProjDist = layerParams_.maxProjectionDist();
    if (maxProjDist < GREAT)
    {
        forAll(stationaryPts,pointi)
        {
            if
            (
                pointWallDist[pointi].valid(dummyTrackData)
                && pointWallDist[pointi].distSqr() > maxProjDist
            )
            {
                stationaryPts[pointi] = 4;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    boolList hangingNodes(mesh_.nPoints(), false);
    List<labelList> hangingNbrs(mesh_.nPoints(), labelList());
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        forAll(mesh_.cells(), celli)
        {
            const labelList& cEdges = mesh_.cellEdges()[celli];
            const labelList& cPts = mesh_.cellPoints()[celli];
            labelHashSet cEdgesSet(cEdges);
            label cLevel = cellLevel_[celli];

            forAll(cPts, cPtI)
            {
                label pointi = cPts[cPtI];
                label pLevel = pointLevel_[pointi];
                if (pLevel > cLevel)
                {
                    const labelList& pEdges = mesh_.pointEdges()[pointi];
                    DynamicList<label> oPts(pEdges.size());
                    DynamicList<label> hPts(pEdges.size());

                    if (boundaryPoints[pointi] > 1)
                    {
                        continue;
                    }

                    bool grownUpPt(boundaryPoints[pointi] == 1 ? true : false);

                    label nBoundaryNbr = 0;
                    label nGrownUpNbr = 0;
                    label nZoneNbr = 0;
                    forAll(pEdges, pEI)
                    {
                        label edgei = pEdges[pEI];
                        if (cEdgesSet.found(edgei))
                        {
                            edge e = mesh_.edges()[edgei];
                            label otherPt(e[0] == pointi ? e[1] : e[0]);

                            if (boundaryPoints[otherPt] == 1)
                            {
                                nGrownUpNbr++;
                            }
                            else if
                            (
                                boundaryPoints[otherPt] != -1
                            )
                            {
                                nBoundaryNbr++;
                            }
                            else if
                            (
                                stationaryPts[otherPt] == -1
                                || stationaryPts[otherPt] == 1
                            )
                            {
                                nZoneNbr++;
                            }

                            if (pointLevel_[otherPt] > cLevel)
                            {
                                if (grownUpPt)
                                {
                                    if (boundaryPoints[otherPt] == 1)
                                    {
                                        hPts.append(otherPt);
                                    }
                                }
                                else
                                {
                                    hPts.append(otherPt);
                                }
                            }
                            else
                            {
                                if (grownUpPt)
                                {
                                    if (boundaryPoints[otherPt] == 1)
                                    {
                                        oPts.append(otherPt);
                                    }
                                }
                                else
                                {
                                    oPts.append(otherPt);
                                }
                            }
                        }
                    }

                    if
                    (
                        boundaryPoints[pointi] == -1
                        &&
                        (
                            (
                                nGrownUpNbr > 0
                                && (nBoundaryNbr+nGrownUpNbr) > 1
                            )
                            || nZoneNbr > 1
                        )

                    )
                    {
                        continue;
                    }

                    if (hPts.size() == 4)
                    {
                        hangingNodes[pointi] = true;
                        hangingNbrs[pointi] = hPts;
                    }
                    else if (oPts.size() == 2)
                    {
                        hangingNodes[pointi] = true;
                        hangingNbrs[pointi] = oPts;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            hangingNodes,
            orEqOp<bool>(),
            false
        );
    }

    scalarField initExpansion(meshPoints.size(), 20);
    pointField outerPos = newPoints;

    for (label iter = 0; iter < nSmoothIter; iter++)
    {
        if (iter % 10 == 0)
        {
            Info<<"Iter : "<<iter<<endl;
        }

        resetPts = -1;

        //Used for carying out stack points and old points field
        // for error reversal
        tPts = newPoints;
        outerPos = newPoints;

        if (debug)
        {
            Time& runTime = const_cast<Time&>(mesh_.time());
            runTime++;
            pointField origPts(mesh_.points());
            mesh_.movePoints(newPoints);
            mesh_.write();
            mesh_.movePoints(origPts);
        }

        updateGeometry
        (
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes
        );

        if (layerParams_.squishTol() > -1 && iter > 10 && (iter % 2 == 0))
        {
            deactivateSquishCells
            (
                layerPointType,
                boundaryPoints,
                newFaceCentres,
                newFaceAreas,
                newCellCentres,
                stationaryPts
            );
        }

        newPoints = vector::zero;

        projectOuterLayerFaces
        (
            outerShellPtr,
            layerPointType,
            layerFaceType,
            ftrPointOrigin,
            stationaryPts,
            stackEdges,
            isMasterFace,
            newFaceAreas,
            newCellCentres,
            pointOrigin,
            pointSurfNormal,
            tPts,
            newPoints,
            resetPts
        );

        limitLayerHeight
        (
            ftrPointOrigin,
            layerPointType,
            maxLayerThickness,
            pointSurfNormal,
            pointOrigin,
            newPoints
        );

        resetPoints
        (
            iter,
            controller,
            layerPointType,
            layerFaceType,
            layerEdges,
            medialDist,
            medialLayerHeight,
            pointOrigin,
            pointSurfNormal,
            tPts,
            outerPos,
            newCellVolumes,
            newFaceAreas,
            newFaceCentres,
            newCellCentres,
            newPoints,
            resetPts,
            false
        );

        tPts = newPoints;

        updateGeometry
        (
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes
        );

        volumeSmoothPts
        (
            controller,
            isMasterEdge,
            hangingNodes,
            hangingNbrs,
            layerPointType,
            ftrPointOrigin,
            boundaryPoints,
            stationaryPts,
            layerWeights,
            nonLayerWeights,
            stackEdges,
            featureEdges,
            newCellVolumes,
            newFaceAreas,
            newFaceCentres,
            newCellCentres,
            pointOrigin,
            pointSurfNormal,
            tPts,
            newPoints,
            resetPts
        );

        resetPoints
        (
            iter,
            controller,
            layerPointType,
            layerFaceType,
            layerEdges,
            medialDist,
            medialLayerHeight,
            pointOrigin,
            pointSurfNormal,
            tPts,
            outerPos,
            newCellVolumes,
            newFaceAreas,
            newFaceCentres,
            newCellCentres,
            newPoints,
            resetPts,
            true
        );

        tPts = vector::zero;

        //Final iteration try to correct concave outer cells
        if
        (
            (controller.algorithm() == meshControl::EXTRUDE)
            && (iter == nSmoothIter-1)
        )
        {
            correctConcaveCells
            (
                layerPointType,
                layerFaceType,
                newFaceCentres,
                newCellCentres,
                newFaceAreas,
                tPts,
                newPoints
            );
        }

        tPts = vector::zero;

        reSnap
        (
            grownUpSnapPP,
            fzonePP,
            grownUpGeometryPtr,
            grownUpZoneGeometryPtr,
            boundaryPoints,
            stationaryPts,
            newFaceAreas,
            newFaceCentres,
            tPts,
            newPoints
        );

        outerPos = vector::zero;

        reProjectOuter
        (
            ftrPointOrigin,
            layerPointType,
            maxLayerThickness,
            boundaryEdges,
            stackEdges,
            pointSurfNormal,
            pointOrigin,
            newPoints,
            outerPos
        );

        stretchLayerMesh
        (
            controller,
            layerPointType,
            boundaryEdges,
            stackEdges,
            boundaryPoints,
            stationaryPts,
            pp,
            meshPoints,
            layerMethod,
            totalNumLayers,
            layerCount,
            pointSurfNormal,
            pointOrigin,
            unsmoothedPointNormals,
            outerPos,
            initExpansion,
            newPoints,
            tPts
        );

        if
        (
            adjustFinalNLayers
            && (iter == nSmoothIter-1)
        )
        {
            adjustNumLayers
            (
                pp,
                layerMethod,
                initExpansion,
                pointOrigin,
                outerPos,
                totalNumLayers,
                adjustedNLayers
            );
        }
    }

    //reuse tPts field for calculating point displacement
    tPts = (newPoints-mesh_.points());

    syncTools::syncPointList
    (
        mesh_,
        tPts,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    newPoints = mesh_.points() + tPts;

    mesh_.movePoints(newPoints);
}


void Foam::layerManipulate::stretchLayers
(
    const meshControl& controller
)
{
    Info<<"Stretching layer cells "<<endl;

    const volScalarField& layerCells =
        mesh_.lookupObject<volScalarField>("layerStacks");

    autoPtr<indirectPrimitivePatch> ppPtr = makeLayerPatch();
    indirectPrimitivePatch& pp = ppPtr();

    boolList layerEdges(mesh_.nEdges(), false);
    forAll(mesh_.cells(), cellI)
    {
        if (layerCells[cellI] != -1)
        {
            const labelList& cellEdges = mesh_.cellEdges()[cellI];
            forAll(cellEdges, cEI)
            {
                layerEdges[cellEdges[cEI]] = true;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh_,
        layerEdges,
        orEqOp<bool>(),
        false
    );

    const labelList& meshPoints = pp.meshPoints();
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh_.edges(),
            mesh_.pointEdges()
        )
    );

    //Set layer growth method
    labelList ppPtLevel(meshPoints.size(), -1);
    labelList ppPtLevelMin(meshPoints.size(), labelMax);
    forAll(meshPoints, pti)
    {
        const labelList& pFaces = pp.pointFaces()[pti];
        forAll(pFaces, pfi)
        {
            label meshFaceI = pp.addressing()[pFaces[pfi]];
            label own = mesh_.faceOwner()[meshFaceI];
            ppPtLevel[pti] = max(cellLevel_[own],ppPtLevel[pti]);
            ppPtLevelMin[pti] = min(cellLevel_[own],ppPtLevelMin[pti]);
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        ppPtLevel,
        maxEqOp<label>(),
        label(-1)         // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        ppPtLevelMin,
        minEqOp<label>(),
        labelMax         // null value
    );

    List<Tuple2<word,scalar>> layerMethod(meshPoints.size());
    scalarField ppMaxLayerThickness(meshPoints.size(), GREAT);
    setLayerMethod
    (
        pp,
        meshPoints,
        ppPtLevel,
        ppPtLevelMin,
        layerMethod,
        ppMaxLayerThickness
    );

    boolList boundaryFaces(mesh_.nFaces(), false);
    boolList boundaryEdges(mesh_.nEdges(), false);

    //Boundary point types
    // -1 : not on boundary
    // 0 : grown patch
    // 1 : grown up
    // 2 : grown up feature pt
    // 3 : stationary corner points
    labelList boundaryPoints(mesh_.nPoints(), -1);

    forAll(pp, i)
    {
        boundaryFaces[pp.addressing()[i]] = true;
    }

    forAll(meshPoints, ptI)
    {
        boundaryPoints[meshPoints[ptI]] = 0;
    }

    autoPtr<indirectPrimitivePatch> grownUpPPPtr = makeGrownUpPatch();
    indirectPrimitivePatch& grownUpPP = grownUpPPPtr();

    autoPtr<indirectPrimitivePatch> fzonePPPtr = makeFaceZonePatch();
    indirectPrimitivePatch& fzonePP = fzonePPPtr();

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh_));
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));
    boolList featurePts(mesh_.nPoints(), false);
    boolList featureEdges(mesh_.nEdges(), false);

    //Mark zone feature edges
    if
    (
        layerParams_.dualReSnapZones()
        && returnReduce(fzonePP.size(), sumOp<label>()) != 0
    )
    {
        boolList excludedFaces(fzonePP.size(), false);

        const labelList zoneEdges
        (
            fzonePP.meshEdges
            (
                mesh_.edges(),
                mesh_.pointEdges()
            )
        );

        edgeClassification eClass
        (
            mesh_,
            mesh_.points(),
            fzonePP,
            zoneEdges,
            excludedFaces,
            0.8191,
            0.8191,
            -GREAT,
            true //flip face based on owner zone
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        forAll(zoneEdges, edgeI)
        {
            if
            (
                eType[edgeI].first() == edgeClassification::CONCAVE
                || eType[edgeI].first() == edgeClassification::CONVEX
            )
            {
                label meshEdgeI = zoneEdges[edgeI];
                edge e = mesh_.edges()[meshEdgeI];
                featureEdges[meshEdgeI] = true;
                featurePts[e[0]] = true;
                featurePts[e[1]] = true;
            }
        }
    }

    //Calculate feature edges grown up patch
    {
        boolList excludedFaces(grownUpPP.size(), false);

        const labelList meshEdgesGrownUp
        (
            grownUpPP.meshEdges
            (
                mesh_.edges(),
                mesh_.pointEdges()
            )
        );

        edgeClassification eClass
        (
            mesh_,
            mesh_.points(),
            grownUpPP,
            meshEdgesGrownUp,
            excludedFaces,
            0.8191,
            0.8191
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        forAll(meshEdgesGrownUp, edgeI)
        {
            if
            (
                eType[edgeI].first() == edgeClassification::CONCAVE
                || eType[edgeI].first() == edgeClassification::CONVEX
            )
            {
                label meshEdgeI = meshEdgesGrownUp[edgeI];
                edge e = mesh_.edges()[meshEdgeI];
                featureEdges[meshEdgeI] = true;
                featurePts[e[0]] = true;
                featurePts[e[1]] = true;
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            featureEdges,
            orEqOp<bool>(),
            false
        );

        syncTools::syncPointList
        (
            mesh_,
            featurePts,
            orEqOp<bool>(),
            false
        );

        labelList nFeatureEdges(mesh_.nPoints(), 0);
        forAll(mesh_.points(), pointi)
        {
            labelList pEdges = mesh_.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                if (isMasterEdge[pEdges[pEI]] && featureEdges[pEdges[pEI]])
                {
                    nFeatureEdges[pointi]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            nFeatureEdges,
            plusEqOp<label>(),
            label(0)
        );

        forAll(mesh_.points(), pointi)
        {
            if (boundaryPoints[pointi] == -1)
            {
                if (featurePts[pointi])
                {
                    if (nFeatureEdges[pointi] > 2)
                    {
                        boundaryPoints[pointi] = 3;
                    }
                    else
                    {
                        boundaryPoints[pointi] = 2;
                    }
                }
            }
        }

        forAll(grownUpPP.meshPoints(), ptI)
        {
            label pointi = grownUpPP.meshPoints()[ptI];

            if (boundaryPoints[pointi] == -1)
            {
                boundaryPoints[pointi] = 1;
            }
        }
    }

    forAll(pp.edges(), edgeI)
    {
        boundaryEdges[meshEdges[edgeI]] = true;
    }

    syncTools::syncEdgeList
    (
        mesh_,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh_,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    labelList cellLayerNumber(mesh_.nCells(), 0);
    labelList stackEdgeLayer(mesh_.nEdges(), 0);

    //layer face type
    // 1 : side faces
    // 2 : inner stack faces
    // 3 : outer stack faces
    labelList layerFaceType(mesh_.nFaces(), -1);
    labelList layerPointType(mesh_.nPoints(), -1);

    calculateLayerStacks
    (
        pp,
        boundaryFaces,
        boundaryEdges,
        cellLayerNumber,
        layerFaceType,
        layerPointType,
        stackEdgeLayer
    );

    boolList stackEdges(mesh_.nEdges(), false);
    forAll(mesh_.edges(), edgeI)
    {
        if (stackEdgeLayer[edgeI] != 0)
        {
            stackEdges[edgeI] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        stackEdges,
        orEqOp<bool>(),
        false
    );

    pointField pointNormals,unsmoothedPointNormals;
    //Calculate area weighted point normals (feature points fixed)
    {
        boolList excludedFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh_,
            mesh_.points(),
            pp,
            meshEdges,
            excludedFaces,
            0.707,//-0.5,
            0.707//-0.5
        );

        pointNormals = eClass.calculatePointNormals
        (
            excludedFaces,
            4,
            true
        );
        unsmoothedPointNormals =
            eClass.calculatePointNormals(excludedFaces, 0, true);
    }

    mesh_.clearOut();

    labelList layerCount(mesh_.nPoints(), -1);

    pointField pointOrigin(mesh_.nPoints(), vector::zero);
    pointField pointSurfNormal(mesh_.nPoints(), vector::zero);

    scalarField maxLayerThickness(mesh_.nPoints(), GREAT);

    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        layerCount[meshPointI] = label(0);
        pointOrigin[meshPointI] = pp.localPoints()[ptI];
        pointSurfNormal[meshPointI] = -pointNormals[ptI];
        maxLayerThickness[meshPointI] = ppMaxLayerThickness[ptI];
    }

    labelList nVisited(mesh_.nPoints(), 0);

    boolList edgeSet(mesh_.nEdges(), false);

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh_,
            pointOrigin,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            pointSurfNormal,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            layerCount,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncEdgeList
        (
            mesh_,
            edgeSet,
            orEqOp<bool>(),
            false              // null value
        );

        syncTools::syncPointList
        (
            mesh_,
            maxLayerThickness,
            minEqOp<scalar>(),
            GREAT
        );

        forAll(mesh_.edges(), edgeI)
        {
            if (!isMasterEdge[edgeI])
            {
                continue;
            }

            if (!boundaryEdges[edgeI] && stackEdges[edgeI] && !edgeSet[edgeI])
            {
                edge e = mesh_.edges()[edgeI];

                if (layerCount[e[0]] != -1)
                {
                    if (layerPointType[e[0]] < 2)
                    {
                        if (layerCount[e[1]] == -1)
                        {
                            layerCount[e[1]] = layerCount[e[0]] + 1;
                        }
                        edgeSet[edgeI] = true;

                        nVisited[e[1]]++;

                        maxLayerThickness[e[1]] = maxLayerThickness[e[0]];
                        pointOrigin[e[1]] += pointOrigin[e[0]];
                        pointSurfNormal[e[1]] += pointSurfNormal[e[0]];
                        nSet++;
                    }
                }

                if (!edgeSet[edgeI] && layerCount[e[1]] != -1)
                {
                    if (layerPointType[e[1]] < 2)
                    {
                        if (layerCount[e[0]] == -1)
                        {
                            layerCount[e[0]] = layerCount[e[1]] + 1;
                        }
                        edgeSet[edgeI] = true;

                        nVisited[e[0]]++;

                        maxLayerThickness[e[0]] = maxLayerThickness[e[1]];
                        pointOrigin[e[0]] += pointOrigin[e[1]];
                        pointSurfNormal[e[0]] += pointSurfNormal[e[1]];
                        nSet++;
                    }
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        nVisited,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh_.points(), pointi)
    {
        if (nVisited[pointi] > 1)
        {
            pointOrigin[pointi] /= nVisited[pointi];
            pointSurfNormal[pointi] /= nVisited[pointi];
        }
    }

    //Check for unitinitilised hanging points
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        forAll(mesh_.points(), pointi)
        {
            if (layerPointType[pointi] == 2 && nVisited[pointi] == 0)
            {
                const labelList& pEdges = mesh_.pointEdges()[pointi];
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    edge e = mesh_.edges()[edgei];
                    label otherPt(e[0] ==  pointi ? e[1] : e[0]);
                    if (layerEdges[edgei] && layerPointType[otherPt] == 2)
                    {
                        maxLayerThickness[pointi] = min
                        (
                            maxLayerThickness[otherPt],
                            maxLayerThickness[pointi]
                        );
                        pointOrigin[pointi] += pointOrigin[otherPt];
                        pointSurfNormal[pointi] += pointSurfNormal[otherPt];
                        nVisited[pointi]++;
                    }
                }
            }
            else
            {
                nVisited[pointi] = 1;
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            nVisited,
            plusEqOp<label>(),
            label(0)
        );

        syncTools::syncPointList
        (
            mesh_,
            maxLayerThickness,
            minEqOp<scalar>(),
            GREAT
        );

        syncTools::syncPointList
        (
            mesh_,
            pointOrigin,
            plusEqOp<vector>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            pointSurfNormal,
            plusEqOp<vector>(),
            vector::zero
        );

        forAll(mesh_.points(), pointi)
        {
            if (nVisited[pointi] > 1)
            {
                pointOrigin[pointi] /= nVisited[pointi];
                pointSurfNormal[pointi] /= nVisited[pointi];
            }
        }
    }

    labelList totalNumLayers(mesh_.nPoints(), -1);
    forAll(mesh_.faces(), facei)
    {
        if (layerFaceType[facei] == 3)
        {
            face f = mesh_.faces()[facei];
            forAll(f,fp)
            {
                label pointi = f[fp];
                totalNumLayers[pointi] = max
                (
                    totalNumLayers[pointi],layerCount[pointi]
                );
            }
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh_,
            totalNumLayers,
            maxEqOp<label>(),
            label(-1)
        );

        forAll(mesh_.edges(), edgeI)
        {
            if (!boundaryEdges[edgeI] && stackEdges[edgeI])
            {
                edge e = mesh_.edges()[edgeI];

                if
                (
                    totalNumLayers[e[0]] != -1 && totalNumLayers[e[1]] == -1
                    && boundaryPoints[e[0]] != 0
                )
                {
                    totalNumLayers[e[1]] = totalNumLayers[e[0]];
                    nSet++;
                }
                else if
                (
                    totalNumLayers[e[1]] != -1 && totalNumLayers[e[0]] == -1
                    && boundaryPoints[e[1]] != 0
                )
                {
                    totalNumLayers[e[0]] = totalNumLayers[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }


    //Prevent faceZones from having blended profile
    labelList stationaryPts(mesh_.nPoints(), -2);
    forAll(fzonePP.meshPoints(), pointi)
    {
        label meshPointI = fzonePP.meshPoints()[pointi];
        stationaryPts[meshPointI] = -1;
    }

    syncTools::syncPointList
    (
        mesh_,
        stationaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    scalarField initExpansion(meshPoints.size(), 20);
    pointField newPoints(mesh_.points());
    pointField tPts(mesh_.points());

    pointField outerPos(mesh_.points());
    boolList pointVisited(mesh_.nPoints(), false);
    forAll(mesh_.points(), pointi)
    {
        if (layerPointType[pointi] >= 2)
        {
            pointVisited[pointi] = true;
        }
    }

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh_,
            outerPos,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            pointVisited,
            orEqOp<bool>(),
            false
        );

        forAll(mesh_.edges(), edgeI)
        {
            if (!boundaryEdges[edgeI] && stackEdges[edgeI])
            {
                edge e = mesh_.edges()[edgeI];

                if
                (
                    pointVisited[e[0]] && !pointVisited[e[1]]
                )
                {
                    pointVisited[e[1]] = true;
                    outerPos[e[1]] = outerPos[e[0]];
                    nSet++;
                }
                else if
                (
                    pointVisited[e[1]] && !pointVisited[e[0]]
                )
                {
                    pointVisited[e[0]] = true;
                    outerPos[e[0]] = outerPos[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }

    stretchLayerMesh
    (
        controller,
        layerPointType,
        boundaryEdges,
        stackEdges,
        boundaryPoints,
        stationaryPts,
        pp,
        meshPoints,
        layerMethod,
        totalNumLayers,
        layerCount,
        pointSurfNormal,
        pointOrigin,
        unsmoothedPointNormals,
        outerPos,
        initExpansion,
        newPoints,
        tPts
    );

    mesh_.movePoints(newPoints);
}





void Foam::layerManipulate::volumeSmoothPts
(
    const meshControl& controller,
    const PackedBoolList& isMasterEdge,
    const boolList& hangingNodes,
    const List<labelList>& hangingNbrs,
    const labelList& layerPointType,
    const labelList& ftrPointOrigin,
    const labelList& boundaryPoints,
    const labelList& stationaryPts,
    const scalarField& layerWeights,
    const scalarField& nonLayerWeights,
    const boolList& stackEdges,
    const boolList& featureEdges,
    const scalarField& newCellVolumes,
    const vectorField& newFaceAreas,
    const pointField& newFaceCentres,
    const pointField& newCellCentres,
    const pointField& pointOrigin,
    const vectorField& pointSurfNormal,
    const pointField& tPts,
    pointField& newPoints,
    labelList& resetPts
)
{
    scalar weight = 0.5;
    scalar invWeight = 1-weight;
    scalar reducedWeight = 0.1;
    scalar invReducedWeight = 1-reducedWeight;

    scalarField sumVol(mesh_.nPoints(), 0);
    pointField avePt(mesh_.nPoints(), 0);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const volScalarField& layerCells =
        mesh_.lookupObject<volScalarField>("layerStacks");

    forAll(newPoints, pointi)
    {
        if
        (
            boundaryPoints[pointi] != 0
            && (stationaryPts[pointi] < 0 || stationaryPts[pointi] == 3)
            && (layerPointType[pointi] == -1 || layerPointType[pointi] >= 2)
         )
        {
            if (boundaryPoints[pointi] > 0)
            {
                if
                (
                    layerPointType[pointi] == 2
                    && boundaryPoints[pointi] == 1
                )
                {
                    labelList pFaces = mesh_.pointFaces()[pointi];
                    forAll(pFaces, pFI)
                    {
                        label faceI = pFaces[pFI];
                        label patchI = patches.whichPatch(faceI);

                        if (patchI != -1 && !patches[patchI].coupled())
                        {
                            label own = mesh_.faceOwner()[faceI];
                            if (layerCells[own] != -1)
                            {
                                scalar wf = -1;
                                if
                                (
                                    controller.algorithm()
                                    == meshControl::DUAL
                                 )
                                {
                                    wf = layerWeights[pointi]
                                        *mag(newFaceAreas[faceI]);
                                }
                                else
                                {
                                    label  level =
                                        cellLevel_[mesh_.faceOwner()[faceI]];
                                    scalar len =
                                        edge0Length_ / pow(2., level);
                                    wf = layerWeights[pointi]
                                        *sqrt(mag(newFaceAreas[faceI]))/len;
                                }
                                avePt[pointi] +=
                                    (wf * newFaceCentres[faceI]);
                                sumVol[pointi] += wf;
                            }
                            else
                            {
                                scalar wf = -1;
                                if
                                (
                                    controller.algorithm()
                                    == meshControl::DUAL
                                 )
                                {
                                    wf = nonLayerWeights[pointi]
                                        *mag(newFaceAreas[faceI]);
                                }
                                else
                                {
                                    label level =
                                        cellLevel_[mesh_.faceOwner()[faceI]];
                                    scalar len = edge0Length_
                                        / pow(2., level);
                                    wf = nonLayerWeights[pointi]
                                        *sqrt(mag(newFaceAreas[faceI]))/len;
                                }
                                avePt[pointi] +=
                                    (wf * newFaceCentres[faceI]);
                                sumVol[pointi] += wf;
                            }
                        }
                    }
                }
                else
                {
                    if (boundaryPoints[pointi] == 1)
                    {
                        if (hangingNodes[pointi])
                        {
                            const labelList& nbrs = hangingNbrs[pointi];
                            forAll(nbrs, nbri)
                            {
                                avePt[pointi] += newPoints[nbrs[nbri]];
                                sumVol[pointi] += scalar(1);
                            }
                        }
                        else
                        {
                            labelList pFaces = mesh_.pointFaces()[pointi];
                            forAll(pFaces, pFI)
                            {
                                label faceI = pFaces[pFI];
                                label patchI = patches.whichPatch(faceI);
                                if
                                (
                                    patchI != -1
                                    && !patches[patchI].coupled()
                                 )
                                {
                                    scalar wf = -1;
                                    if
                                    (
                                        controller.algorithm()
                                        == meshControl::DUAL
                                     )
                                    {
                                        wf = mag(newFaceAreas[faceI]);
                                    }
                                    else
                                    {
                                        label  level = cellLevel_
                                            [mesh_.faceOwner()[faceI]];
                                        scalar len = edge0Length_
                                            / pow(2., level);
                                        wf = sqrt(mag(newFaceAreas[faceI]))
                                            /len;
                                    }

                                    avePt[pointi] +=
                                        (wf * newFaceCentres[faceI]);
                                    sumVol[pointi] += wf;
                                }
                            }
                        }
                    }
                    else if (boundaryPoints[pointi] == 2)
                    {
                        labelList pEdges = mesh_.pointEdges()[pointi];
                        forAll(pEdges, pEI)
                        {
                            label edgeI = pEdges[pEI];
                            if (featureEdges[edgeI] && isMasterEdge[edgeI])
                            {
                                edge e = mesh_.edges()[edgeI];
                                scalar edgeWeight = e.mag(newPoints);
                                if
                                (
                                    controller.algorithm()
                                    != meshControl::DUAL
                                 )
                                {
                                    label level = max
                                    (
                                        pointLevel_[e[0]],
                                        pointLevel_[e[1]]
                                    );
                                    scalar len = edge0Length_
                                        / pow(2., level);
                                    edgeWeight /= len;
                                }
                                avePt[pointi] +=
                                    (edgeWeight * e.centre(newPoints));
                                sumVol[pointi] += edgeWeight;
                            }
                        }
                    }
                }
            }
            else
            {
                if (hangingNodes[pointi] && layerPointType[pointi] != 2)
                {
                    const labelList& nbrs = hangingNbrs[pointi];
                    forAll(nbrs, nbri)
                    {
                        avePt[pointi] += newPoints[nbrs[nbri]];
                        sumVol[pointi] += scalar(1);
                    }
                }
                else
                {
                    labelList pCells = mesh_.pointCells()[pointi];
                    forAll(pCells, pCI)
                    {
                        label cellI = pCells[pCI];

                        if (layerPointType[pointi] == 2)
                        {
                            scalar volWeight = -1;
                            if (controller.algorithm() == meshControl::DUAL)
                            {
                                volWeight = mag(newCellVolumes[cellI]);
                            }
                            else
                            {
                                label  level = cellLevel_[cellI];
                                scalar len = edge0Length_ / pow(2., level);
                                volWeight =
                                    cbrt(mag(newCellVolumes[cellI]))/len;
                            }

                            if (layerCells[cellI] != -1)
                            {
                                scalar wf = layerWeights[pointi]
                                    * volWeight;
                                avePt[pointi] += (wf*newCellCentres[cellI]);
                                sumVol[pointi] += wf;
                            }
                            else
                            {
                                scalar wf = nonLayerWeights[pointi]
                                    * volWeight;
                                avePt[pointi] += (wf*newCellCentres[cellI]);
                                sumVol[pointi] += wf;
                            }
                        }
                        else
                        {
                            scalar volWeight = -1;
                            if (controller.algorithm() == meshControl::DUAL)
                            {
                                volWeight = mag(newCellVolumes[cellI]);
                            }
                            else
                            {
                                label  level = cellLevel_[cellI];
                                scalar len = edge0Length_ / pow(2., level);
                                volWeight =
                                    cbrt(mag(newCellVolumes[cellI]))/len;
                            }
                            avePt[pointi] +=
                                (volWeight * newCellCentres[cellI]);
                            sumVol[pointi] += volWeight;
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        sumVol,
        plusEqOp<scalar>(),
        scalar(0)
     );

    syncTools::syncPointList
    (
        mesh_,
        avePt,
        plusEqOp<point>(),
        vector::zero
     );

    resetPts = -1;

    forAll(newPoints, pointi)
    {
        if
        (
            boundaryPoints[pointi] != 0
            && (layerPointType[pointi] == -1 || layerPointType[pointi] >= 2)
        )
        {
            if (sumVol[pointi] > VSMALL)
            {
                if (stationaryPts[pointi] < 0 && layerPointType[pointi] >= 2)
                {
                    point updatedPt;

                    if (ftrPointOrigin[pointi] > 0)
                    {
                        updatedPt = invReducedWeight*tPts[pointi]
                            +reducedWeight*(avePt[pointi]/sumVol[pointi]);
                    }
                    else
                    {
                        updatedPt = invWeight*tPts[pointi]
                            +weight*(avePt[pointi]/sumVol[pointi]);
                    }

                    newPoints[pointi] = updatedPt;

                    //Do not freeze boundary points
                    if (boundaryPoints[pointi] > 0)
                    {
                        continue;
                    }

                    labelList pEdges = mesh_.pointEdges()[pointi];

                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];

                        if (stackEdges[edgeI])
                        {
                            edge e = mesh_.edges()[edgeI];
                            label otherPt = (e[0] ==  pointi ? e[1] : e[0]);

                            vector toSurf = updatedPt
                                - pointOrigin[otherPt];
                            toSurf /= (mag(toSurf) + SMALL);
                            vector sNorm = pointSurfNormal[otherPt];
                            sNorm /= (mag(sNorm) + SMALL);

                            vector toSurfOrig =
                                tPts[pointi] - pointOrigin[otherPt];
                            toSurfOrig /= (mag(toSurfOrig) + SMALL);

                            scalar dP1 = (toSurf&sNorm);
                            scalar dP2 = (toSurfOrig&sNorm);

                            if
                            (
                                dP1 < dP2 && dP1 < 0.866
                                && ftrPointOrigin[otherPt] != -1
                             )
                            {
                                resetPts[pointi] = label(1);
                                break;
                            }
                        }
                    }
                }
                else
                {
                    if (stationaryPts[pointi] > -1)
                    {
                        if (stationaryPts[pointi] == 3)
                        {
                            newPoints[pointi] = invReducedWeight*tPts[pointi]
                                + reducedWeight*(avePt[pointi]/sumVol[pointi]);
                        }
                        else
                        {
                            newPoints[pointi] = tPts[pointi];
                        }
                    }
                    else
                    {
                        newPoints[pointi] = invWeight*tPts[pointi]
                            +weight*(avePt[pointi]/sumVol[pointi]);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        resetPts,
        maxEqOp<label>(),
        label(-1)
     );
}


void Foam::layerManipulate::stretchLayerMesh
(
    const meshControl& controller,
    const labelList& layerPointType,
    const boolList& boundaryEdges,
    const boolList& stackEdges,
    const labelList& boundaryPoints,
    const labelList& stationaryPts,
    const indirectPrimitivePatch& pp,
    const labelList& meshPoints,
    const List<Tuple2<word,scalar>>& layerMethod,
    const labelList& totalNumLayers,
    const labelList& layerCount,
    const vectorField& pointSurfNormal,
    const pointField& pointOrigin,
    const pointField& unsmoothedPointNormals,
    const pointField& outerPos,
    scalarField& initExpansion,
    pointField& newPoints,
    pointField& tPts
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    scalar tol = 0.0001;
    scalar maxExpansion = 4;
    scalar minStretch = layerParams_.minStretch();
    point greatPoint(GREAT, GREAT, GREAT);

    scalarField expansionRatio(mesh_.nPoints(), -1);

    const List<wordList>& layerSpec = layerParams_.layerSpec();

    const label fixedLayerCount = layerParams_.nFCHLayers();
    const List<Switch> fixedFCH = layerParams_.fixedFCH();
    bool setExtactFCH = false;
    forAll(fixedFCH, patchi)
    {
        if (fixedFCH[patchi])
        {
            word method = layerSpec[patchi][1];
            if (method == "fch" || method == "rfch")
            {
                setExtactFCH = true;
                break;
            }
        }
    }

    vectorField exactFCH;
    boolList fixedFCHPts;
    if (setExtactFCH)
    {
        exactFCH.setSize(mesh_.nPoints(), vector::zero);
        fixedFCHPts.setSize(meshPoints.size(), false);
        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];
            label patchi = patches.whichPatch(meshFaceI);
            if (fixedFCH[patchi])
            {
                face f = pp.localFaces()[i];
                forAll(f,fp)
                {
                    fixedFCHPts[f[fp]] = true;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh_,
            meshPoints,
            fixedFCHPts,
            orEqOp<bool>(),
            false         // null value
        );
    }

    //Calculate stretching at primitive patch points
    forAll(meshPoints, ptI)
    {
        label meshPointI = meshPoints[ptI];

        if (layerMethod[ptI].first() == "fch")
        {
            scalar firstCellHeight = layerMethod[ptI].second();

            scalar xn = initExpansion[ptI];
            scalar x = -1.0;

            label repeat = 0;

            label numLayers = totalNumLayers[meshPointI];
            vector layerVec = outerPos[meshPointI]-pointOrigin[meshPointI];
            scalar layerHeight = mag(layerVec);

            if (layerHeight < SMALL)
            {
                expansionRatio[meshPointI] = 1.0;
                continue;
            }

            vector unitLayerVec = layerVec / layerHeight;
            scalar dProd = max
            (
                0.5,
                mag(unitLayerVec & unsmoothedPointNormals[ptI])
            );

            //correct first cell height
            firstCellHeight /= dProd;

            if (firstCellHeight > layerHeight)
            {
                expansionRatio[meshPointI] = 1.0;
                continue;
            }

            if (xn < 1.0 + tol && xn > 1.0 - tol)
            {
                xn = maxExpansion;
                x = -1.0;
                repeat = 0;
                while (!(repeat++ >= 200 || mag(x - xn) < tol))
                {
                    x = xn;
                    scalar f1 = firstCellHeight * pow(x, numLayers)
                        + layerHeight*(1 - x) - firstCellHeight;
                    scalar f2 = firstCellHeight * numLayers
                        * pow(x, numLayers -1) - layerHeight;
                    xn = x - f1/(f2+SMALL);
                }
            }
            else
            {
                while (!(repeat++ >= 200 || mag(x - xn) < tol))
                {
                    x = xn;
                    scalar f1 = firstCellHeight * pow(x, numLayers)
                        + layerHeight*(1 - x) - firstCellHeight;
                    scalar f2 = firstCellHeight * numLayers
                        * pow(x, numLayers -1) - layerHeight;
                    xn = x - f1/(f2+SMALL);
                }
            }

            if (xn < 1.0 + tol && xn > 1.0 - tol)
            {
                xn = 0.0;
                x = -1.0;
                repeat = 0;
                while (!(repeat++ >= 200 || mag(x - xn) < tol))
                {
                    x = xn;
                    scalar f1 = firstCellHeight * pow(x, numLayers)
                        + layerHeight*(1 - x) - firstCellHeight;
                    scalar f2 = firstCellHeight * numLayers
                        * pow(x, numLayers -1) - layerHeight;
                    xn = x - f1/(f2+SMALL);
                }
            }
            scalar limitMin = max(minStretch,xn);
            initExpansion[ptI] = limitMin;
            expansionRatio[meshPointI] = min(limitMin,maxExpansion);
            if (setExtactFCH && fixedFCHPts[ptI])
            {
                xn = expansionRatio[meshPointI];
                scalar actualFCH;
                scalar layer2height;
                if (xn < 1.0 - tol || xn > 1.0 + tol)
                {
                    actualFCH = layerHeight*(1. - xn)
                        /(1.-pow(xn, numLayers));
                    layer2height = (1.+xn)*actualFCH;
                }
                else
                {
                    actualFCH = layerHeight/numLayers;
                    layer2height = 2*actualFCH;
                }
                scalar inputFCH = layerMethod[ptI].second();
                if (inputFCH < 0.8*layer2height)
                {
                    exactFCH[meshPointI] =
                        -inputFCH*unsmoothedPointNormals[ptI];
                }
            }
        }
        else
        {
            expansionRatio[meshPointI] = layerMethod[ptI].second();
        }
    }

    tPts = -greatPoint;

    //Check if layer blending required
    bool blend = false;
    if (controller.algorithm() == meshControl::EXTRUDE)
    {
        blend = layerParams_.extrudeBlend();
    }

    //Calculate layer stretched layer profile
    {
        boolList pointVisited(mesh_.nPoints(), false);
        forAll(meshPoints, ptI)
        {
            label meshPointI = meshPoints[ptI];
            pointVisited[meshPointI] = true;
        }

        forAll(mesh_.points(), pointi)
        {
            if (layerPointType[pointi] >= 2)
            {
                tPts[pointi] = outerPos[pointi];
            }
        }

        while (true)
        {
            label nSet = 0;

            syncTools::syncPointList
            (
                mesh_,
                expansionRatio,
                maxEqOp<scalar>(),
                maxExpansion
            );

            if (setExtactFCH)
            {
                syncTools::syncPointList
                (
                    mesh_,
                    exactFCH,
                    maxMagSqrEqOp<point>(),
                    vector::zero
                );
            }

            syncTools::syncPointList
            (
                mesh_,
                pointVisited,
                orEqOp<bool>(),
                false
             );

            forAll(mesh_.edges(), edgeI)
            {
                if (!boundaryEdges[edgeI] && stackEdges[edgeI])
                {
                    edge e = mesh_.edges()[edgeI];

                    label pt0 = -1;
                    label pt1 = -1;

                    if (pointVisited[e[0]] && !pointVisited[e[1]])
                    {
                        pt0 = e[0];
                        pt1 = e[1];
                    }
                    else if (pointVisited[e[1]] && !pointVisited[e[0]])
                    {
                        pt0 = e[1];
                        pt1 = e[0];
                    }

                    if
                    (
                        pt0 == -1
                        || boundaryPoints[pt1] == 0
                        || layerPointType[pt1] >=2
                    )
                    {
                        continue;
                    }

                    pointVisited[pt1] = true;
                    expansionRatio[pt1] = expansionRatio[pt0];
                    if (setExtactFCH)
                    {
                        exactFCH[pt1] = exactFCH[pt0];
                    }
                    vector dir = (outerPos[pt1]-pointOrigin[pt1]);
                    scalar str = expansionRatio[pt0];

                    label lcount = layerCount[pt1];

                    if
                    (
                        setExtactFCH && lcount <= fixedLayerCount
                        && mag(exactFCH[pt0]) > 0
                    )
                    {
                        tPts[pt1] = pointOrigin[pt1] + lcount*exactFCH[pt0];
                    }
                    else if (str > 1.0 + tol || str < 1.0 - tol)
                    {
                        scalar lh = mag(dir) + SMALL;
                        scalar fch = lh *
                        (
                            (1.-str)/(1.-pow(str,totalNumLayers[pt1]))
                        );
                        dir /= lh;
                        scalar ratio = fch*(1.-pow(str,lcount))
                            /(1.-str);

                        if
                        (
                            blend && boundaryPoints[pt1] != 1
                            && stationaryPts[pt1] != -1
                        )
                        {
                            scalar weight = min
                            (
                                scalar(1.0),
                                max
                                (
                                    scalar(0.0),ratio/lh
                                )
                            );

                            weight = sin
                            (
                                weight
                                * 0.5*Foam::constant::mathematical::pi
                            );
                            weight = Foam::pow(weight, scalar(0.25));

                            dir = (1.-weight)*pointSurfNormal[pt1]
                                + weight*dir;
                        }
                        tPts[pt1] = pointOrigin[pt1]
                            + ratio*dir;
                    }
                    else
                    {
                        scalar ratio = lcount / scalar(totalNumLayers[pt1]);
                        if
                        (
                            blend && boundaryPoints[pt1] != 1
                            && stationaryPts[pt1] != -1
                        )
                        {
                            scalar lh = mag(dir) + SMALL;
                            vector udir(dir/lh);

                            scalar weight = min
                            (
                                scalar(1.0),
                                max
                                (
                                    scalar(0.0),ratio
                                )
                            );

                            weight = sin
                            (
                                weight
                                * 0.5*Foam::constant::mathematical::pi
                            );
                            weight = Foam::pow(weight, scalar(0.25));

                            dir = lh*
                            (
                                (1.-weight)*pointSurfNormal[pt1]
                                + weight*udir
                            );
                        }
                        tPts[pt1] = pointOrigin[pt1]
                            + ratio*dir;
                    }
                    nSet++;
                }
            }

            if (returnReduce(nSet, sumOp<label>()) == 0)
            {
                break;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        tPts,
        maxEqOp<point>(),
        -greatPoint
    );

    //Reset layer points into newPoints field
    forAll(tPts, pointi)
    {
        if (tPts[pointi] != -greatPoint)
        {
            newPoints[pointi] = tPts[pointi];
        }
    }

    tPts = (newPoints-mesh_.points());

    syncTools::syncPointList
    (
        mesh_,
        tPts,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    newPoints = mesh_.points() + tPts;
}


void Foam::layerManipulate::adjustNumLayers
(
    const indirectPrimitivePatch& pp,
    const List<Tuple2<word,scalar>>& layerMethod,
    const scalarField& initExpansion,
    const pointField& pointOrigin,
    const pointField& outerPos,
    const labelList& totalNumLayers,
    labelList& adjustedNLayers
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& meshPoints = pp.meshPoints();
    scalarField minTargetExp(meshPoints.size(), GREAT);
    scalarField maxTargetExp(meshPoints.size(), -GREAT);

    const List<Tuple2<scalar,scalar>>& targetExpansions =
       layerParams_.targetExpansions();

    forAll(pp, i)
    {
        label meshfacei = pp.addressing()[i];
        label patchi = patches.whichPatch(meshfacei);
        const face& f = pp.localFaces()[i];
        const Tuple2<scalar,scalar>& targetExp = targetExpansions[patchi];
        scalar minExp = targetExp.first();
        scalar maxExp = targetExp.second();
        forAll(f,fp)
        {
            label pointi = f[fp];
            minTargetExp[pointi] = min(minTargetExp[pointi],minExp);
            maxTargetExp[pointi] = max(maxTargetExp[pointi],maxExp);
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        minTargetExp,
        minEqOp<scalar>(),
        -GREAT         // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        maxTargetExp,
        maxEqOp<scalar>(),
        GREAT         // null value
    );

    forAll(meshPoints, ptI)
    {
        label meshPointI = meshPoints[ptI];
        if (layerMethod[ptI].first() == "fch")
        {
            scalar eR = initExpansion[ptI];
            scalar expMin = minTargetExp[ptI];
            scalar expMax = maxTargetExp[ptI];
            scalar expTarget = -1;
            if (eR < expMin)
            {
                expTarget = expMin;
            }
            if (eR > expMax)
            {
                expTarget = expMax;
            }
            label numLayers = totalNumLayers[meshPointI];
            scalar firstCellHeight = layerMethod[ptI].second();
            vector layerVec = outerPos[meshPointI]-pointOrigin[meshPointI];
            scalar layerHeight = mag(layerVec);
            label wantedLayers = numLayers;
            if (wantedLayers > 1 && expTarget > -1)
            {
               if (expTarget < 1.0 - SMALL || expTarget > 1.0 + SMALL)
               {
                   scalar a1 = 1.
                       - (layerHeight/firstCellHeight)*(1. - expTarget);
                   wantedLayers = max(label(2),log(a1)/log(expTarget));
               }
               else
               {
                   wantedLayers = max(label(2),layerHeight/firstCellHeight);
               }
            }
            adjustedNLayers[ptI] = wantedLayers;
        }
    }

    return;
}


void Foam::layerManipulate::reProjectOuter
(
    const labelList& ftrPointOrigin,
    const labelList& layerPointType,
    const scalarField& maxLayerThickness,
    const boolList& boundaryEdges,
    const boolList& stackEdges,
    const vectorField& pointSurfNormal,
    const pointField& pointOrigin,
    const pointField& newPoints,
    pointField& outerPos
)
{
    boolList pointVisited(mesh_.nPoints(), false);
    forAll(mesh_.points(), pointi)
    {
        if (layerPointType[pointi] >= 2)
        {
            pointVisited[pointi] = true;

            if (ftrPointOrigin[pointi] == 2)
            {
                point startPt = pointOrigin[pointi];
                scalar layerHeight = mag(newPoints[pointi]-startPt);
                point endPt = startPt
                    + 2.*layerHeight*pointSurfNormal[pointi];
                pointHit lHit = linePointRef(startPt, endPt).nearestDist
                    (newPoints[pointi]);
                if (lHit.hit())
                {
                    outerPos[pointi] =
                        0.5*(lHit.hitPoint()+newPoints[pointi]);
                }
                else
                {
                    outerPos[pointi] = newPoints[pointi];
                }
            }
            else
            {
                outerPos[pointi] = newPoints[pointi];
            }
        }
    }

    limitLayerHeight
    (
        ftrPointOrigin,
        layerPointType,
        maxLayerThickness,
        pointSurfNormal,
        pointOrigin,
        outerPos
    );

    while (true)
    {
        label nSet = 0;

        syncTools::syncPointList
        (
            mesh_,
            outerPos,
            maxMagSqrEqOp<point>(),
            vector::zero
         );

        syncTools::syncPointList
        (
            mesh_,
            pointVisited,
            orEqOp<bool>(),
            false
         );

        forAll(mesh_.edges(), edgeI)
        {
            if (!boundaryEdges[edgeI] && stackEdges[edgeI])
            {
                edge e = mesh_.edges()[edgeI];

                if
                (
                    pointVisited[e[0]] && !pointVisited[e[1]]
                )
                {
                    pointVisited[e[1]] = true;
                    outerPos[e[1]] = outerPos[e[0]];
                    nSet++;
                }
                else if
                (
                    pointVisited[e[1]] && !pointVisited[e[0]]
                )
                {
                    pointVisited[e[0]] = true;
                    outerPos[e[0]] = outerPos[e[1]];
                    nSet++;
                }
            }
        }

        if (returnReduce(nSet, sumOp<label>()) == 0)
        {
            break;
        }
    }
}


void Foam::layerManipulate::reSnap
(
    const indirectPrimitivePatch& grownUpSnapPP,
    const indirectPrimitivePatch& fzonePP,
    const autoPtr<searchableSurfaces>& grownUpGeometryPtr,
    const autoPtr<searchableSurfaces>& grownUpZoneGeometryPtr,
    const labelList& boundaryPoints,
    const labelList& stationaryPts,
    const vectorField& newFaceAreas,
    const pointField& newFaceCentres,
    pointField& tPts,
    pointField& newPoints
)
{
    //Re-snap if any patches selected
    if (returnReduce(grownUpSnapPP.size(), sumOp<label>()) != 0)
    {
        const labelList& grownUpSnapPPMeshPts = grownUpSnapPP.meshPoints();

        labelList snapSurf = identity(grownUpGeometryPtr().size());
        pointField snapPts(grownUpSnapPPMeshPts.size());
        scalarField snapDist(grownUpSnapPPMeshPts.size());
        forAll(snapPts, ptI)
        {
            label meshPointI = grownUpSnapPPMeshPts[ptI];
            snapPts[ptI] = newPoints[meshPointI];
            snapDist[ptI] =
                sqr(edge0Length_ / (1<<pointLevel_[meshPointI]));
        }

        List<pointIndexHit> info;
        grownUpGeometryPtr().findNearest
        (
            snapPts,
            snapDist,
            snapSurf,
            info
        );

        forAll(snapPts, ptI)
        {
            if (info[ptI].hit())
            {
                label meshPointI = grownUpSnapPPMeshPts[ptI];
                if
                (
                    boundaryPoints[meshPointI] == 1
                    || boundaryPoints[meshPointI] == 2
                )
                {
                    tPts[meshPointI] += info[ptI].hitPoint()
                        - newPoints[meshPointI];
                }
            }
        }
    }

    if
    (
        grownUpZoneGeometryPtr().size()
        && returnReduce(fzonePP.meshPoints().size(), sumOp<label>()) != 0
    )
    {
        const labelList& zoneMeshPts = fzonePP.meshPoints();
        pointField snapPts(zoneMeshPts.size());
        scalarField snapDist(zoneMeshPts.size());

        labelList snapSurf = identity(grownUpZoneGeometryPtr().size());
        pointField ppNewZonePts(zoneMeshPts.size(), vector::zero);
        scalarField ppNewZoneWeights(zoneMeshPts.size(), Zero);

        //locally smooth zone points
        forAll(snapPts, ptI)
        {
            if
            (
                stationaryPts[zoneMeshPts[ptI]] > -1
                || boundaryPoints[zoneMeshPts[ptI]] > -1
             )
            {
                continue;
            }

            const labelList& pFaces = fzonePP.pointFaces()[ptI];
            forAll(pFaces, pfI)
            {
                label facei = fzonePP.addressing()[pFaces[pfI]];
                scalar fArea = mag(newFaceAreas[facei]);
                ppNewZonePts[ptI] += fArea*newFaceCentres[facei];
                ppNewZoneWeights[ptI] += fArea;
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            zoneMeshPts,
            ppNewZonePts,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh_,
            zoneMeshPts,
            ppNewZoneWeights,
            plusEqOp<scalar>(),
            scalar(0.0)         // null value
        );

        forAll(snapPts, ptI)
        {
            label meshPointI = zoneMeshPts[ptI];
            point avePt = newPoints[meshPointI];
            if (ppNewZoneWeights[ptI] > SMALL)
            {
                avePt = ppNewZonePts[ptI]/ppNewZoneWeights[ptI];
            }

            snapPts[ptI] = avePt;
            snapDist[ptI] =
                sqr(edge0Length_ / (1<<pointLevel_[meshPointI]));
        }

        List<pointIndexHit> info;
        grownUpZoneGeometryPtr().findNearest
        (
            snapPts,
            snapDist,
            snapSurf,
            info
         );

        forAll(snapPts, ptI)
        {
            if (info[ptI].hit())
            {
                label meshPointI = zoneMeshPts[ptI];
                if
                (
                    (
                       boundaryPoints[meshPointI] == -1
                       || boundaryPoints[meshPointI] == 1
                    )
                    &&
                    (
                        stationaryPts[meshPointI] < 0
                        || stationaryPts[meshPointI] == 3
                    )
                 )
                {
                    tPts[meshPointI] += info[ptI].hitPoint()
                        - newPoints[meshPointI];
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        tPts,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    newPoints = newPoints + tPts;
}


void Foam::layerManipulate::correctConcaveCells
(
    const labelList& layerPointType,
    const labelList& layerFaceType,
    const pointField& newFaceCentres,
    const pointField& newCellCentres,
    const vectorField& newFaceAreas,
    pointField& tPts,
    pointField& newPoints
)
{
    const volScalarField& layerCells =
        mesh_.lookupObject<volScalarField>("layerStacks");

    faceSet errorFaces(mesh_, "errorFaces", mesh_.nFaces()/100+1);
    labelList checkFaces(identity(mesh_.nFaces()));
    Foam::polyMeshGeometry::checkFacePyramids
    (
        false,
        -SMALL,
        mesh_,
        newFaceCentres,
        newCellCentres,
        newFaceAreas,
        newPoints,
        checkFaces,
        List<labelPair>(0),
        &errorFaces
    );

    boolList errorCells(mesh_.nCells(), false);
    forAllConstIter(labelHashSet, errorFaces, iter)
    {
        label facei = iter.key();
        label own = mesh_.faceOwner()[facei];
        errorCells[own] = true;
        if (mesh_.isInternalFace(facei))
        {
            label nei = mesh_.faceNeighbour()[facei];
            errorCells[nei] = true;
        }
    }

    labelList nMoved(mesh_.nPoints(), 0);
    forAll(mesh_.cells(), celli)
    {
        if (errorCells[celli] &&  layerCells[celli] < 0)
        {
            const cell& c = mesh_.cells()[celli];
            DynamicList<label> stackFaces(c.size());
            forAll(c, cFI)
            {
                label facei = c[cFI];
                if (layerFaceType[facei] == 3)
                {
                    stackFaces.append(facei);
                }
            }

            if (stackFaces.size() > 1)
            {
                const labelList& cEdges = mesh_.cellEdges()[celli];
                labelHashSet stackFaceSet(stackFaces);
                forAll(cEdges, cEI)
                {
                    label edgei = cEdges[cEI];
                    const edge& e = mesh_.edges()[edgei];
                    if
                    (
                        layerPointType[e[0]] == 2
                        && layerPointType[e[1]] == 2
                    )
                    {
                        const labelList& eFaces =
                            mesh_.edgeFaces()[edgei];
                        DynamicList<label> cellEdgeFaces(2);
                        forAll(eFaces, eFI)
                        {
                            label facei = eFaces[eFI];
                            if
                            (
                                mesh_.faces()[facei].size() == 3
                                && stackFaceSet.found(facei)
                             )
                            {
                                cellEdgeFaces.append(facei);
                            }
                        }
                        if (cellEdgeFaces.size() == 2)
                        {
                            scalarField fAreas(2);
                            pointField fCentres(2);
                            pointField fNormals(2);
                            DynamicList<label> otherPts(4);
                            forAll(cellEdgeFaces, cEFI)
                            {
                                label facei = cellEdgeFaces[cEFI];
                                const face& f = mesh_.faces()[facei];
                                label own = mesh_.faceOwner()[facei];
                                fNormals[cEFI] = f.areaNormal(newPoints);
                                if (own == celli)
                                {
                                    fNormals[cEFI] = -fNormals[cEFI];
                                }
                                fCentres[cEFI] = f.centre(newPoints);
                                fAreas[cEFI] = mag(fNormals[cEFI]);
                                fNormals[cEFI] /= (fAreas[cEFI]+SMALL);
                                forAll(f,fp)
                                {
                                    if (f[fp] != e[0] && f[fp] != e[1])
                                    {
                                        otherPts.append(f[fp]);
                                    }
                                }
                            }
                            scalar sumArea = fAreas[0]+fAreas[1];

                            vector eNorm
                            (
                                fAreas[0]*fNormals[0]
                                +fAreas[1]*fNormals[1]
                            );
                            eNorm /= (sumArea + SMALL);

                            if (mag(eNorm) < VSMALL)
                            {
                                continue;
                            }

                            point eC = e.centre(newPoints);
                            vector fC2eC = fCentres[0] - eC;
                            //check if concave cell
                            if ((fC2eC&eNorm) < 0)
                            {
                                otherPts.append(e[0]);
                                otherPts.append(e[1]);
                                point aveFC = fAreas[0]*fCentres[0]
                                    + fAreas[1]*fCentres[1];
                                aveFC /= sumArea;
                                plane fPlane(aveFC, eNorm);
                                forAll(otherPts, oPtI)
                                {
                                    label pointi = otherPts[oPtI];
                                    vector disp = fPlane.nearestPoint
                                    (
                                        newPoints[pointi]
                                    ) - newPoints[pointi];
                                    tPts[pointi] += disp;
                                    nMoved[pointi]++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        tPts,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        nMoved,
        plusEqOp<label>(),
        label(0)         // null value
     );

    forAll(nMoved, pointi)
    {
        if (nMoved[pointi] > 0)
        {
            point avePt = newPoints[pointi]
                + (tPts[pointi]/nMoved[pointi]);
            newPoints[pointi] = avePt;
        }
    }
}


void Foam::layerManipulate::resetPoints
(
    const label iter,
    const meshControl& controller,
    const labelList& layerPointType,
    const labelList& layerFaceType,
    const boolList& layerEdges,
    const scalarField& medialDist,
    const scalarField& medialLayerHeight,
    const pointField& pointOrigin,
    const pointField& pointSurfNormal,
    const pointField& tPts,
    const pointField& outerPos,
    scalarField& newCellVolumes,
    vectorField& newFaceAreas,
    pointField& newFaceCentres,
    pointField& newCellCentres,
    pointField& newPoints,
    labelList& resetPts,
    bool correctErrors
)
{
    if (iter > 0)
    {
        forAll(mesh_.points(), pointi)
        {
            if (layerPointType[pointi] == 2)
            {
                vector toSurf = newPoints[pointi] - pointOrigin[pointi];
                scalar dist0 = mag(toSurf);

                if
                (
                    dist0 > 0.8*medialDist[pointi]
                    && dist0 < 1.1*medialLayerHeight[pointi]
                 )
                {
                    newPoints[pointi] = tPts[pointi];
                }
            }
        }
    }

    forAll(mesh_.points(), pointi)
    {
        if (resetPts[pointi] == 1)
        {
            newPoints[pointi] = tPts[pointi];
            const labelList& pEdges = mesh_.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgeI = pEdges[pEI];
                if (!layerEdges[edgeI])
                {
                    edge e = mesh_.edges()[edgeI];
                    label otherPt = (e[0] ==  pointi ? e[1] : e[0]);
                    if (layerPointType[otherPt] >= 2)
                    {
                        newPoints[otherPt] = tPts[otherPt];
                    }
                }
            }
        }
    }

    //If errors introduced reset outer layer point
    if (iter > 0 && correctErrors)
    {
        scalar maxOrtho = layerParams_.dualOrtho();

        const volScalarField& layerCells =
            mesh_.lookupObject<volScalarField>("layerStacks");

        updateGeometry
        (
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes
        );

        faceSet errorFaces(mesh_, "errorFaces", mesh_.nFaces()/100+1);
        labelList checkFaces(identity(mesh_.nFaces()));

        Foam::polyMeshGeometry::checkFacePyramids
        (
            false,
            -SMALL,
            mesh_,
            newFaceCentres,
            newCellCentres,
            newFaceAreas,
            newPoints,
            checkFaces,
            List<labelPair>(0),
            &errorFaces
         );

        if (maxOrtho < 180.0-SMALL && iter > 5)
        {
            Foam::polyMeshGeometry::checkFaceDotProduct
            (
                false,
                maxOrtho,
                scalar(180),
                mesh_,
                newFaceCentres,
                newCellCentres,
                newFaceAreas,
                checkFaces,
                List<labelPair>(0),
                &errorFaces
            );
        }

        boolList errorCells(mesh_.nCells(), false);
        forAllConstIter(labelHashSet, errorFaces, iter)
        {
            label facei = iter.key();
            label own = mesh_.faceOwner()[facei];
            if (layerCells[own] > -1)
            {
                errorCells[own] = true;
            }

            if (mesh_.isInternalFace(facei))
            {
                label nei = mesh_.faceNeighbour()[facei];
                if (layerCells[nei] > -1)
                {
                    errorCells[nei] = true;
                }
            }
        }

        boolList errorPts(mesh_.nPoints(), false);
        forAll(errorCells, celli)
        {
            if (errorCells[celli])
            {
                const labelList& cellPts = mesh_.cellPoints()[celli];
                forAll(cellPts, cPI)
                {
                    errorPts[cellPts[cPI]] = true;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh_,
            errorPts,
            orEqOp<bool>(),
            false
         );

        forAll(errorPts, pointi)
        {
            if (errorPts[pointi])
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];
                    if (layerCells[celli] > -1)
                    {
                        errorCells[celli] = true;
                    }
                }
            }
        }

        scalarField faceLayerCells(mesh_.nFaces(), -1);
        forAll(errorCells, celli)
        {
            if (errorCells[celli])
            {
                scalar layerCellI = layerCells[celli];
                const cell& c = mesh_.cells()[celli];
                forAll(c, cFI)
                {
                    label facei = c[cFI];
                    faceLayerCells[facei] =
                        max(faceLayerCells[facei],layerCellI);
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh_,
            faceLayerCells,
            maxEqOp<scalar>()
        );

        while (true)
        {
            label nReset = 0;
            forAll(errorCells, celli)
            {
                if (!errorCells[celli] && layerCells[celli] > -1)
                {
                    scalar layerCellI = layerCells[celli];
                    const cell& c = mesh_.cells()[celli];
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        if (faceLayerCells[facei] == layerCellI)
                        {
                            errorCells[celli] = true;
                            break;
                        }
                    }
                }
            }

            if (returnReduce(nReset, sumOp<label>()) == 0)
            {
                break;
            }

            forAll(errorCells, celli)
            {
                if (errorCells[celli])
                {
                    scalar layerCellI = layerCells[celli];
                    const cell& c = mesh_.cells()[celli];
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        faceLayerCells[facei] =
                            max(faceLayerCells[facei],layerCellI);
                    }
                }
            }

            syncTools::syncFaceList
            (
                mesh_,
                faceLayerCells,
                maxEqOp<scalar>()
             );
        }

        forAll(mesh_.cells(), celli)
        {
            if (errorCells[celli] && layerCells[celli] > -1)
            {
                const labelList& cPts = mesh_.cellPoints()[celli];
                forAll(cPts, cPtI)
                {
                    label pointi = cPts[cPtI];
                    resetPts[pointi] = label(0);
                }
            }
        }

        //Check first non layer cells becoming too deformed (acute angles)
        if (controller.algorithm() == meshControl::EXTRUDE)
        {
            forAll(mesh_.cells(), celli)
            {
                if (layerCells[celli] == -1)
                {
                    const cell& c = mesh_.cells()[celli];
                    label nLayerFaces = 0;
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        const face& f = mesh_.faces()[facei];

                        if (layerFaceType[facei] == 3)
                        {
                            if (f.size() == 3)
                            {
                                nLayerFaces++;
                            }
                            else
                            {
                                nLayerFaces += 2;
                            }
                        }
                    }

                    if (nLayerFaces > 2)
                    {
                        const labelList& cPts = mesh_.cellPoints()[celli];
                        point newCC = newCellCentres[celli];
                        scalar cellLength = edge0Length_
                            / (1<<cellLevel_[celli]);
                        bool warpedCell = false;

                        forAll(cPts, cPtI)
                        {
                            label pointi = cPts[cPtI];
                            if (layerPointType[pointi] >=2)
                            {
                                point pN = pointSurfNormal[pointi];
                                scalarField vProj(cPts.size(), Zero);
                                forAll(cPts, cPtJ)
                                {
                                    label pointj = cPts[cPtJ];
                                    vector n = newPoints[pointj] - newCC;
                                    vProj[cPtJ] = (pN & n);
                                }
                                // Get normal 'span' of cell
                                scalar minVal = min(vProj);
                                scalar maxVal = max(vProj);
                                scalar dW = (maxVal - minVal);
                                if (dW < 0.2*cellLength)
                                {
                                    warpedCell = true;
                                    break;
                                }
                            }
                        }
                        if (warpedCell)
                        {
                            forAll(cPts, cPtI)
                            {
                                resetPts[cPts[cPtI]] = label(0);
                            }
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            resetPts,
            maxEqOp<label>(),
            label(-1)
        );

        forAll(newPoints, pointi)
        {
            if (resetPts[pointi] == 0 && layerPointType[pointi] >= 2)
            {
                newPoints[pointi] = outerPos[pointi];
            }
        }
    }
}


void Foam::layerManipulate::projectOuterLayerFaces
(
    const autoPtr<indirectPrimitivePatch>& outerShellPtr,
    const labelList& layerPointType,
    const labelList& layerFaceType,
    const labelList& ftrPointOrigin,
    const labelList& stationaryPts,
    const boolList& stackEdges,
    const PackedBoolList& isMasterFace,
    const vectorField& newFaceAreas,
    const pointField& newCellCentres,
    const pointField& pointOrigin,
    const vectorField& pointSurfNormal,
    const pointField& tPts,
    pointField& newPoints,
    labelList& resetPts
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& owners = mesh_.faceOwner();

    pointField neiNewCellCentres(mesh_.nFaces()-mesh_.nInternalFaces());
    for
    (
        label faceI = mesh_.nInternalFaces();
        faceI < mesh_.nFaces();
        faceI++
    )
    {
        neiNewCellCentres[faceI-mesh_.nInternalFaces()] =
            newCellCentres[owners[faceI]];
    }
    syncTools::swapBoundaryFaceList(mesh_, neiNewCellCentres);

    scalarField weightSum(mesh_.nPoints(), 0);

    forAll(mesh_.faces(), faceI)
    {
        if
        (
            !isMasterFace[faceI] || layerFaceType[faceI] != 3
        )
        {
            continue;
        }

        label patchI = patches.whichPatch(faceI);

        if (patchI == -1 || patches[patchI].coupled())
        {
            label own = owners[faceI];

            point nCC = vector::zero;
            if (patchI == -1)
            {
                label nei = mesh_.faceNeighbour()[faceI];
                nCC = newCellCentres[nei];
            }
            else
            {
                nCC = neiNewCellCentres[faceI-mesh_.nInternalFaces()];
            }

            point midPoint = 0.5*(newCellCentres[own]+nCC);
            scalar wt = mag(newFaceAreas[faceI])+ SMALL;
            face f = mesh_.faces()[faceI];
            vector pNorm = vector::zero;
            forAll(f, fp)
            {
                label pointi = f[fp];
                pNorm += pointSurfNormal[pointi];
            }
            scalar magPNorm = mag(pNorm);

            if (magPNorm > SMALL)
            {
                pNorm /= magPNorm;
                plane fPlane(midPoint, pNorm);

                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (stationaryPts[pointi] < 0)
                    {
                        scalar fWeight = (1./wt);
                        weightSum[pointi] += fWeight;
                        vector disp = fPlane.nearestPoint(tPts[pointi])
                            - tPts[pointi];
                        newPoints[pointi] += fWeight*disp;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        weightSum,
        plusEqOp<scalar>(),
        scalar(0)
     );

    syncTools::syncPointList
    (
        mesh_,
        newPoints,
        plusEqOp<point>(),
        vector::zero
     );

    forAll(newPoints, pointi)
    {
        if
        (
            layerPointType[pointi] >= 2 && weightSum[pointi] > SMALL
        )
        {
            point updatedPt = tPts[pointi]
                + (newPoints[pointi]/weightSum[pointi]);

            newPoints[pointi] = updatedPt;

            labelList pEdges = mesh_.pointEdges()[pointi];

            forAll(pEdges, pEI)
            {
                label edgeI = pEdges[pEI];

                if (stackEdges[edgeI])
                {
                    edge e = mesh_.edges()[edgeI];
                    label otherPt = (e[0] ==  pointi ? e[1] : e[0]);

                    vector toSurf = updatedPt - pointOrigin[otherPt];
                    scalar toSurfMag = mag(toSurf);
                    toSurf /= (toSurfMag + SMALL);

                    vector sNorm = pointSurfNormal[otherPt];
                    sNorm /= (mag(sNorm) + SMALL);

                    vector toSurfOrig = tPts[pointi]
                        - pointOrigin[otherPt];
                    scalar toSurfOrigMag = mag(toSurfOrig);
                    toSurfOrig /= (toSurfOrigMag  + SMALL);

                    scalar dP1 = (toSurf&sNorm);
                    scalar dP2 = (toSurfOrig&sNorm);

                    if (dP1 < dP2 && dP1 < 0.866)
                    {
                        resetPts[pointi] = label(1);
                        break;
                    }
                }
            }
        }
        else
        {
            newPoints[pointi] = tPts[pointi];
        }
    }

    //Check for inverted edges dues to projection
    if (outerShellPtr.valid())
    {
        const indirectPrimitivePatch& outerShell = outerShellPtr();

        const labelList& outerPts =  outerShell.meshPoints();
        const edgeList& outerEdges = outerShell.edges();
        forAll(outerPts, pti)
        {
            label pointi = outerPts[pti];

            const labelList& pEdges = outerShell.pointEdges()[pti];
            forAll(pEdges, pEI)
            {
                const edge& e = outerEdges[pEdges[pEI]];
                label otherpointi = outerPts[e.otherVertex(pti)];

                if
                (
                    ftrPointOrigin[pointi] != 0
                    && ftrPointOrigin[otherpointi] != 0
                 )
                {
                    continue;
                }

                vector eVec = newPoints[otherpointi]
                    - newPoints[pointi];
                scalar eVecMag  = mag(eVec);

                vector fVec = pointOrigin[otherpointi]
                    - pointOrigin[pointi];
                scalar fVecMag  = mag(fVec);

                if (eVecMag < SMALL || fVecMag < SMALL)
                {
                    resetPts[pointi] = label(0);
                    resetPts[otherpointi] = label(0);
                }
                else
                {
                    eVec /= eVecMag;
                    fVec /= fVecMag;

                    if
                    (
                        eVecMag < 0.25*fVecMag
                        || (fVec & eVec) < 0
                    )
                    {
                        resetPts[pointi] = label(0);
                        resetPts[otherpointi] = label(0);
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            resetPts,
            maxEqOp<label>(),
            label(-1)
        );
    }

    syncTools::syncPointList
    (
        mesh_,
        resetPts,
        maxEqOp<label>(),
        label(-1)
    );

    forAll(newPoints, pointi)
    {
        if (resetPts[pointi] == 0)
        {
            newPoints[pointi] = tPts[pointi];
        }
    }
}


void Foam::layerManipulate::limitLayerHeight
(
    const labelList& ftrPointOrigin,
    const labelList& layerPointType,
    const scalarField& maxLayerThickness,
    const vectorField& pointSurfNormal,
    const pointField& pointOrigin,
    pointField& pts
)
{
    forAll(layerPointType, pointi)
    {
        if (layerPointType[pointi] == 2)
        {
            vector layerVec = pts[pointi]
                - pointOrigin[pointi];
            scalar lheight = mag(layerVec);
            scalar mlh = maxLayerThickness[pointi];

            if (lheight > mlh && mlh > 0)
            {
               if (ftrPointOrigin[pointi] == -1)
                {
                    vector sNorm = pointSurfNormal[pointi];
                    sNorm /= (mag(sNorm) + SMALL);
                    plane fPlane(pointOrigin[pointi], sNorm);
                    point hPt =
                        fPlane.nearestPoint(pts[pointi]);
                    scalar nheight = mag(pts[pointi]-hPt);
                    if (nheight > mlh)
                    {
                        scalar correction(mlh/nheight);
                        pts[pointi] = pointOrigin[pointi]
                            + correction*layerVec;
                    }
                }
                else
                {
                    scalar correction(mlh/lheight);
                    pts[pointi]  = pointOrigin[pointi]
                        + correction*layerVec;
                }
            }
        }
    }
}


void Foam::layerManipulate::fitLayerPointStack()
{
    Info<<"Fitting Layer Stacks"<<endl;

    autoPtr<indirectPrimitivePatch> ppPtr = makeLayerPatch();
    indirectPrimitivePatch& pp = ppPtr();

    //Create a unique tag for the points based on the labels of layerCells
    labelListList pointTagList = createLayerPointTagList();
    pointField state = mesh_.points();

    boolList needForSynchPoints(mesh_.nPoints(), false);
    boolList alreadyAddedToStack(mesh_.nPoints(), false);
    List<layerPointStack> layerPointStackList(mesh_.nPoints());


    labelList totalNuOfLayers(mesh_.nPoints(), -1);
    forAll(pp.meshPoints(), pI)
    {
        const labelList& pointFaces = pp.pointFaces()[pI];
        forAll(pointFaces, fI)
        {
            const label& fL = pp.addressing()[fI];
            const label& patch = mesh_.boundaryMesh().whichPatch(fL);
            totalNuOfLayers[pI] =
                max(totalNuOfLayers[pI], (layerParams_.numLayers()[patch] + 1));
        }
    }

    forAll(pp.meshPoints(), pI)
    {
        const label& pL = pp.meshPoints()[pI];
        layerPointStack lPS = buildLayerStack
        (
            pI,
            pointTagList,
            pL,
            alreadyAddedToStack,
            totalNuOfLayers
        );
        const label& size = lPS.pointStackLabels().size();
        if (size<lPS.layerNu())
        {
            const label& lastLabel = lPS.pointStackLabels()[size-1];
            needForSynchPoints[lastLabel] = true;
        }
        forAll(lPS.pointStackLabels(), i)
        {
            const label& pL = lPS.pointStackLabels()[i];
            layerPointStackList[pL] = lPS;
            layerPointStackList[pL].setAddress(pL);
        }
        const label& baseLabel = lPS.pointStackLabels()[0];
        for (int i=1; i<size;i++)
        {
            totalNuOfLayers[i] = totalNuOfLayers[baseLabel];
        }
    }


    bool fullySynced = false;
    label numberOfSyncs = 0;
    while (!fullySynced)
    {
        label nuOfSyncPointsInit = 0;
        forAll(mesh_.points(),pI)
        {
            if (needForSynchPoints[pI])
            {
                nuOfSyncPointsInit++;
            }
        }
        reduce(nuOfSyncPointsInit, sumOp<label>());

        syncTools::syncPointList
        (
            mesh_,
            needForSynchPoints,
            maxEqOp<bool>(),
            false
        );

        syncTools::syncPointList
        (
            mesh_,
            totalNuOfLayers,
            maxEqOp<label>(),
            label(0)
        );

        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& polyP = mesh_.boundaryMesh()[patchI];
            if (polyP.coupled())
            {
                forAll(polyP.meshPoints(), pI)
                {
                    label pL = polyP.meshPoints()[pI];
                    if (needForSynchPoints[pL] && !alreadyAddedToStack[pL])
                    {
                        layerPointStack lPS = buildLayerStack
                        (
                            pI,
                            pointTagList,
                            pL,
                            alreadyAddedToStack,
                            totalNuOfLayers
                        );

                        forAll(lPS.pointStackLabels(), i)
                        {
                            const label& pLabel = lPS.pointStackLabels()[i];
                            layerPointStackList[pLabel] = lPS;
                            layerPointStackList[pLabel].setAddress(pL);
                        }
                        const label& size = lPS.pointStackLabels().size();
                        if (size<lPS.layerNu())
                        {
                            const label&
                                lastLabel =lPS.pointStackLabels()[size-1];
                            needForSynchPoints[lastLabel] = true;
                        }
                        const label& baseLabel = lPS.pointStackLabels()[0];
                        for (int i=1; i<size;i++)
                        {
                            totalNuOfLayers[i] = totalNuOfLayers[baseLabel];
                        }
                    }
                }
            }
        }

        label nuOfSyncPointsFinal = 0;
        forAll(mesh_.points(),pI)
        {
            if (needForSynchPoints[pI])
            {
                nuOfSyncPointsFinal++;
            }
        }
        reduce(nuOfSyncPointsFinal, sumOp<label>());

        label newSyncPoints = nuOfSyncPointsFinal - nuOfSyncPointsInit;
        reduce(newSyncPoints, sumOp<label>());

        if (newSyncPoints==0)
        {
            fullySynced = true;
        }
        numberOfSyncs++;
    }

    syncTools::syncPointList
    (
        mesh_,
        alreadyAddedToStack,
        maxEqOp<bool>(),
        false
    );

    //Need to create a copy
    label counter = 0;
    forAll(mesh_.points(), pI)
    {
//        if (needForSynchPoints[pI]||alreadyAddedToStack[pI])
            if (true)
        {
            counter++;
        }
    }

    labelList syncPointsLabelList(counter);
    List<List<point>> stackListPoints(counter);

    for (int i=0;i<numberOfSyncs;i++)
    {
        counter = 0;
        forAll(mesh_.points(), pI)
        {
            if (true)
            {
                syncPointsLabelList[counter] = pI;
                const layerPointStack& lPS = layerPointStackList[pI];
                stackListPoints[counter] = lPS.stackPoints();
                counter++;
            }
        }

        syncTools::syncPointList
        (
            mesh_,
    //        syncPointsLabelList,
            stackListPoints,
            mergePointStacks(),
            List<point>(0)
        );

        forAll(stackListPoints, pI)
        {
            const label& pL = syncPointsLabelList[pI];
            layerPointStack& lPS = layerPointStackList[pL];
            lPS.substitute(stackListPoints[pI]);
        }
        forAll(layerPointStackList, stackI)
        {
            layerPointStack& lPS = layerPointStackList[stackI];
            lPS.merge(layerPointStackList);
        }
    }

    labelList stackIndexing(mesh_.nPoints(), -1);

    forAll(layerPointStackList, stackI)
    {
        const layerPointStack& lPS = layerPointStackList[stackI];
        stackIndexing[stackI] = lPS.getIndex(state[stackI]);
    }

    forAll(layerPointStackList, stackI)
    {
        const layerPointStack& lPS = layerPointStackList[stackI];
        if (lPS.size()>0 && stackIndexing[stackI]==-1)
        {
            Pout<<" COULD NOT FIND POINT, CHECK TOLERANCE"<<endl;
            Pout<<" SEARCHING FOR "<<state[stackI]<<endl;
            Pout<<lPS.stackPoints()<<endl;
        }
    }
    forAll(layerPointStackList, stackI)
    {
        if (stackIndexing[stackI]>-1)
        {
            const layerPointStack& lPS = layerPointStackList[stackI];
            label size = lPS.size();
            pointField pF(size, vector::zero);
            forAll(pF, i)
            {
                pF[i] = lPS.stackPoints()[i];
            }
            leastSquaresCurveFit fit(pF);
            pointField fitPoints = fit.getNewPoints();
            state[stackI] = fitPoints[stackIndexing[stackI]];
        }
    }

    mesh_.movePoints(state);
}

Foam::labelListList Foam::layerManipulate::createLayerPointTagList()
{
    const volScalarField& layerCells =
        mesh_.lookupObject<volScalarField>("layerStacks");

    labelListList  pointTagList(mesh_.points().size());

    forAll(mesh_.points(), pI)
    {
        const labelList& pointCells = mesh_.pointCells()[pI];
        DynamicList<label> pointTag;
        forAll(pointCells, cI)
        {
            const label& cL = pointCells[cI];
            if (layerCells[cL] != -1)
            {
                if (pointTag.size() == 0)
                {
                    pointTag.append(layerCells[cL]);
                }
                else
                {
                    label counter = 0;
                    forAll(pointTag, elem)
                    {
                        if (pointTag[elem] != layerCells[cL])
                        {
                            counter++;
                        }
                    }
                    if (counter == pointTag.size())
                    {
                        pointTag.append(layerCells[cL]);
                    }
                }
            }
        }
        SortableList<label> sorted(pointTag);
        pointTagList[pI] = sorted;
    }

    labelList zeroSizeList(0);

    syncTools::syncPointList
    (
        mesh_,
        pointTagList,
        appendUniqueElements(),
        zeroSizeList
    );

    forAll(mesh_.points(), pI)
    {
        SortableList<label> sorted(pointTagList[pI]);
        pointTagList[pI] = sorted;
    }

    return pointTagList;
}

Foam::layerPointStack Foam::layerManipulate::buildLayerStack
(
    const label& basePoint,
    const labelListList& pointTagList,
    const label& ppL,
    boolList& addedToStack,
    labelList& totalNuOfLayers
)
{
    label pL = ppL;
    label prevL = -1;
//    label counter = 1;

    DynamicList<vector> pointStack(1, mesh_.points()[pL]);
    DynamicList<label> pointStackLabels(1, pL);

    addedToStack[basePoint] = true;

    bool finished = false;

    while (!finished)
    {
        finished = true;
        const labelList& pointPoints = mesh_.pointPoints()[pL];
        forAll(pointPoints, pPI)
        {
            const label& pPL = pointPoints[pPI];
            const labelList& pointTagInit = pointTagList[pL];
            const labelList& pointTagNext = pointTagList[pPL];
            bool matchedTags = compareTags(pointTagInit, pointTagNext);
            if (matchedTags && pPL!=prevL)
            {
//                counter++;
                pointStack.append(mesh_.points()[pPL]);
                pointStackLabels.append(pPL);
                addedToStack[pPL] = true;
                prevL = pL;
                pL = pPL;
                finished = false;
                break;
            }
        }
    }
    const label& layerNu = totalNuOfLayers[basePoint];
    layerPointStack lPS(true, pointStack, pointStackLabels, layerNu);
    return lPS;
}

bool Foam::layerManipulate::compareTags
(
    const labelList& pt1,
    const labelList& pt2
)
{
    bool matchedTags = false;
    if (pt1.size() == pt2.size())
    {
        forAll(pt1, tag)
        {
            matchedTags = true;
            if
            (
                pt2.size() > 0
             && pt1[tag]!=pt2[tag]
            )
            {
                matchedTags = false;
                break;
            }
        }
    }
    return matchedTags;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::layerManipulate::layerManipulate
(
    fvMesh& mesh,
    const layerParameters& layerParams,
    const labelList& cellLevel,
    const labelList& pointLevel,
    const scalar edge0Length
)
:
    mesh_(mesh),
    layerParams_(layerParams),
    cellLevel_(cellLevel),
    pointLevel_(pointLevel),
    edge0Length_(edge0Length)
{}


// ************************************************************************* //
