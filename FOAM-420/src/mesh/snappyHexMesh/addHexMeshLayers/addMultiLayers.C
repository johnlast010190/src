/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM® : Professional Open-source CFD
|   o   O   o    |  Version : 4.2.0
|    o     o     |  Copyright © 2016 ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------

License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM® <http://www.openfoam.org/>.

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

\*---------------------------------------------------------------------------*/

#include "addHexMeshLayers/addMultiLayers.H"
#include "addHexMeshLayers/addHexMeshLayer.H"
#include "layerManipulate/layerManipulate.H"
#include "sets/topoSets/pointSet.H"
#include "sets/topoSets/faceSet.H"
#include "motionSmoother/motionSmoother.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveCell.H"
#include "edgeClassification/edgeClassification.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "autoFanExtrude/autoFanExtrude.H"
#include "autoLayerCellsMerge/autoLayerCellsMerge.H"
#include "autoOptimize/autoOptimize.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(addMultiLayers, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::addMultiLayers::addMultiLayers
(
    const meshControl& controller,
    const labelList& globalToMasterPatch,
    const refinementParameters& refineParams,
    const layerParameters& layerParams,
    const dictionary& grownUpGeometryDict,
    const dictionary& grownUpZoneGeometryDict,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    meshRefinement& meshRefiner
)
:
    controller_(controller),
    globalToMasterPatch_(globalToMasterPatch),
    refineParams_(refineParams),
    layerParams_(layerParams),
    grownUpGeometryDict_(grownUpGeometryDict),
    grownUpZoneGeometryDict_(grownUpZoneGeometryDict),
    decomposer_(decomposer),
    distributor_(distributor),
    meshRefiner_(meshRefiner)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::addMultiLayers::checkCollapses
(
    const indirectPrimitivePatch& pp,
    labelList& ppLayers
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& meshPoints = pp.meshPoints();
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    const List<Switch> cornerCollapse = layerParams_.cornerCollapse();
    const List<Switch> baffleCollapse = layerParams_.baffleCollapse();
    scalar baffleCollapseAngle = layerParams_.dualConcaveCollapse();

    bool tryBaffleCollapse = layerParams_.tryBaffleCollapse();
    bool tryCornerCollapse = layerParams_.tryCornerCollapse();

    if (tryBaffleCollapse || tryCornerCollapse)
    {
        Info<<"Checking for edges to collapse layers on"<<endl;

        boolList collapsedPts(mesh.nPoints(), false);

        scalar baffleAngleCos = -0.5;
        if (tryBaffleCollapse)
        {
            baffleAngleCos = -Foam::cos(degToRad(baffleCollapseAngle));
        }

        boolList excludedFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            pp,
            meshEdges,
            excludedFaces,
            0.707,
            0.707,
            baffleAngleCos
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        boolList checkEdges(mesh.nEdges(), false);
        boolList checkPts(meshPoints.size(), false);

        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];
            label patchi = patches.whichPatch(meshFaceI);

            if (patchi != -1)
            {
                if (cornerCollapse[patchi])
                {
                    const face f = pp.localFaces()[i];
                    forAll(f,fp)
                    {
                        checkPts[f[fp]] = true;
                    }
                }
                if (tryBaffleCollapse && baffleCollapse[patchi])
                {
                    const labelList& fEdges = mesh.faceEdges()[meshFaceI];

                    forAll(fEdges, fEI)
                    {
                        checkEdges[fEdges[fEI]] = true;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            checkPts,
            orEqOp<bool>(),
            false // initial value
         );

        syncTools::syncEdgeList
        (
            mesh,
            checkEdges,
            orEqOp<bool>(),
            false // initial value
         );

        forAll(meshEdges, edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            if
            (
                checkEdges[meshEdgeI]
                && eType[edgeI].first() == edgeClassification::BAFFLE
            )
            {
                const edge& e = pp.edges()[edgeI];
                collapsedPts[meshPoints[e[0]]] = true;
                collapsedPts[meshPoints[e[1]]] = true;
            }
        }

        const PackedBoolList isMasterPPEdge
        (
            meshRefinement::getMasterEdges
            (
                mesh,
                meshEdges
             )
        );

        labelList nConvexEdges(pp.meshPoints().size(),0);
        labelList nConcaveEdges(pp.meshPoints().size(),0);
        labelList nBaffleEdges(pp.meshPoints().size(),0);

        forAll(meshPoints, ptI)
        {
            if (!checkPts[ptI])
            {
                continue;
            }
            const labelList& pEdges = pp.pointEdges()[ptI];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if (!isMasterPPEdge[edgei])
                {
                    continue;
                }

                if (eType[edgei].first() == edgeClassification::CONVEX)
                {
                    nConvexEdges[ptI]++;
                }
                else if (eType[edgei].first() == edgeClassification::CONCAVE)
                {
                    nConcaveEdges[ptI]++;
                }
                else if (eType[edgei].first() == edgeClassification::BAFFLE)
                {
                    nBaffleEdges[ptI]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nConvexEdges,
            plusEqOp<label>(),  // combine op
            label(0)            // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nConcaveEdges,
            plusEqOp<label>(),  // combine op
            label(0)            // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nBaffleEdges,
            plusEqOp<label>(),  // combine op
            label(0)            // null value
        );

        forAll(meshPoints, ptI)
        {
            if
            (
                (nConvexEdges[ptI] == 2 && nConcaveEdges[ptI] == 2)
                ||
                (
                    nConvexEdges[ptI] == 1 && nConcaveEdges[ptI] > 0
                    && nBaffleEdges[ptI] == 1
                )
            )
            {
                const labelList& pEdges = pp.pointEdges()[ptI];
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    const edge& e = pp.edges()[edgei];
                    collapsedPts[meshPoints[e[0]]] = true;
                    collapsedPts[meshPoints[e[1]]] = true;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh,
            collapsedPts,
            orEqOp<bool>(),
            false // initial value
         );

        forAll(meshPoints, ptI)
        {
            label meshPointI = meshPoints[ptI];
            if (collapsedPts[meshPointI])
            {
                const labelList& pFaces = pp.pointFaces()[ptI];
                forAll(pFaces, pFI)
                {
                    label facei = pFaces[pFI];
                    if (ppLayers[facei] > 1)
                    {
                        ppLayers[facei] = label(1);
                    }
                }
            }
        }
    }
}


void Foam::addMultiLayers::checkCollapses
(
    const indirectPrimitivePatch& pp,
    const labelList& ppLevel,
    labelList& minEdgeLevel,
    boolList& collapsedPts
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& meshPoints = pp.meshPoints();
    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));
    const labelListList& edgeFaces = pp.edgeFaces();

    const List<Switch> cornerCollapse = layerParams_.cornerCollapse();
    const List<Switch> baffleCollapse = layerParams_.baffleCollapse();
    scalar baffleCollapseAngle = layerParams_.dualConcaveCollapse();

    bool tryBaffleCollapse = layerParams_.tryBaffleCollapse();
    bool tryCornerCollapse = layerParams_.tryCornerCollapse();

    if (tryBaffleCollapse || tryCornerCollapse)
    {
        Info<<"Checking for edges to collapse layers on"<<endl;

        scalar baffleAngleCos = -0.5;
        if (tryBaffleCollapse)
        {
            baffleAngleCos = -Foam::cos(degToRad(baffleCollapseAngle));
        }

        boolList excludedFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            pp,
            meshEdges,
            excludedFaces,
            0.707,
            0.707,
            baffleAngleCos
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        boolList checkEdges(mesh.nEdges(), false);
        boolList checkPts(meshPoints.size(), false);

        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];
            label patchi = patches.whichPatch(meshFaceI);

            if (patchi != -1)
            {
                if (cornerCollapse[patchi])
                {
                    const face f = pp.localFaces()[i];
                    forAll(f,fp)
                    {
                        checkPts[f[fp]] = true;
                    }
                }
                if (tryBaffleCollapse && baffleCollapse[patchi])
                {
                    const labelList& fEdges = mesh.faceEdges()[meshFaceI];

                    forAll(fEdges, fEI)
                    {
                        checkEdges[fEdges[fEI]] = true;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            checkPts,
            orEqOp<bool>(),
            false // initial value
         );

        syncTools::syncEdgeList
        (
            mesh,
            checkEdges,
            orEqOp<bool>(),
            false // initial value
         );

        forAll(meshEdges, edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            const labelList& eFaces = edgeFaces[edgeI];

            if
            (
                checkEdges[meshEdgeI]
                && eType[edgeI].first() == edgeClassification::BAFFLE
            )
            {
                minEdgeLevel[meshEdgeI] = min
                (
                    minEdgeLevel[meshEdgeI],
                    label(1)
                 );
                const edge& e = pp.edges()[edgeI];
                collapsedPts[meshPoints[e[0]]] = true;
                collapsedPts[meshPoints[e[1]]] = true;
            }
            else
            {
                forAll(eFaces, eFI)
                {
                    minEdgeLevel[meshEdgeI] = min
                    (
                        minEdgeLevel[meshEdgeI],
                        ppLevel[eFaces[eFI]]
                    );
                }
            }
        }

        const PackedBoolList isMasterPPEdge
        (
            meshRefinement::getMasterEdges
            (
                mesh,
                meshEdges
             )
        );

        labelList nConvexEdges(pp.meshPoints().size(),0);
        labelList nConcaveEdges(pp.meshPoints().size(),0);
        labelList nBaffleEdges(pp.meshPoints().size(),0);

        forAll(pp.meshPoints(), ptI)
        {
            if (!checkPts[ptI])
            {
                continue;
            }
            const labelList& pEdges = pp.pointEdges()[ptI];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if (!isMasterPPEdge[edgei])
                {
                    continue;
                }

                if (eType[edgei].first() == edgeClassification::CONVEX)
                {
                    nConvexEdges[ptI]++;
                }
                else if (eType[edgei].first() == edgeClassification::CONCAVE)
                {
                    nConcaveEdges[ptI]++;
                }
                else if (eType[edgei].first() == edgeClassification::BAFFLE)
                {
                    nBaffleEdges[ptI]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nConvexEdges,
            plusEqOp<label>(),  // combine op
            label(0)            // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nConcaveEdges,
            plusEqOp<label>(),  // combine op
            label(0)            // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nBaffleEdges,
            plusEqOp<label>(),  // combine op
            label(0)            // null value
        );

        if (debug)
        {
            simpleVTKWriter layerVTK
            (
                pp.localFaces(),
                pp.localPoints()
            );
            layerVTK.addPointData("nConcave", nConcaveEdges);
            layerVTK.addPointData("nConvex", nConvexEdges);
            layerVTK.addPointData("nBaffle", nBaffleEdges);
            layerVTK.write("pointEdgeFeatureTypes.vtk");
        }

        forAll(meshPoints, ptI)
        {
            if
            (
                (nConvexEdges[ptI] == 2 && nConcaveEdges[ptI] == 2)
                ||
                (
                    nConvexEdges[ptI] == 1 && nConcaveEdges[ptI] > 0
                    && nBaffleEdges[ptI] == 1
                )
            )
            {
                const labelList& pEdges = pp.pointEdges()[ptI];
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    label meshEdgeI = meshEdges[edgei];

                    minEdgeLevel[meshEdgeI] = min
                    (
                        minEdgeLevel[meshEdgeI],
                        label(1)
                    );
                    const edge& e = pp.edges()[edgei];
                    collapsedPts[meshPoints[e[0]]] = true;
                    collapsedPts[meshPoints[e[1]]] = true;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh,
            collapsedPts,
            orEqOp<bool>(),
            false // initial value
         );
    }
    else
    {
        forAll(pp.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            const labelList& eFaces = edgeFaces[edgeI];
            forAll(eFaces, eFI)
            {
                minEdgeLevel[meshEdgeI] = min
                (
                    minEdgeLevel[meshEdgeI],
                    ppLevel[eFaces[eFI]]
                );
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        minEdgeLevel,
        minEqOp<label>(),
        labelMax // initial value
    );

    return;
}


Foam::labelList Foam::addMultiLayers::calculateLayers
(
    const labelList& layerPatchIDs,
    volScalarField& targetLayers,
    volScalarField& layerCount,
    volScalarField& actualLayers
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& numLayers = layerParams_.numLayers();
    DynamicList<label> layerFaces(mesh.nFaces()-mesh.nInternalFaces());

    forAll(patches, patchI)
    {
        if
        (
            !patches[patchI].coupled()
            && numLayers[patchI] > 1
        )
        {
            scalarField& targetPatchLayers =
                targetLayers.boundaryFieldRef()[patchI];
            scalarField& layerPatchCounter =
                layerCount.boundaryFieldRef()[patchI];
            scalarField& actualPatchLayers =
                actualLayers.boundaryFieldRef()[patchI];
            const polyPatch& patch = patches[patchI];

            forAll(patch, i)
            {
                if
                (
                    layerPatchCounter[i] <= targetPatchLayers[i]
                    && layerPatchCounter[i] == actualPatchLayers[i]
                )
                {
                    actualPatchLayers[i] += 1;
                }
            }
        }
    }

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            layerPatchIDs
         )
     );
    indirectPrimitivePatch& pp = ppPtr();

    labelList ppLevel(pp.size(), -1);
    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        label patchI = patches.whichPatch(meshFaceI);

        const polyPatch& patch = patches[patchI];
        const scalarField& actualPatchLayers =
            actualLayers.boundaryField()[patchI];
        ppLevel[i] = label(actualPatchLayers[meshFaceI-patch.start()]);
    }

    labelList minEdgeLevel(mesh.nEdges(), labelMax);
    boolList collapsedPts(mesh.nPoints(), false);

    //Check if edges to collapse layers on
    checkCollapses(pp,ppLevel,minEdgeLevel,collapsedPts);

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && numLayers[patchI] > 0)
        {
            const polyPatch& patch = patches[patchI];
            label startIndex = patch.start();
            scalarField& actualPatchLayers =
                actualLayers.boundaryFieldRef()[patchI];

            forAll(patch, i)
            {
                label meshFaceI = startIndex + i;
                const labelList& fEdges = mesh.faceEdges()[meshFaceI];
                forAll(fEdges, fEI)
                {
                    label edgeI = fEdges[fEI];
                    scalar& aLevel = actualPatchLayers[i];

                    if (minEdgeLevel[edgeI] < aLevel)
                    {
                        aLevel -= 1;
                        break;
                    }
                }
            }
        }
    }

    boolList markedFaces(mesh.nFaces(), false);
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && numLayers[patchI] > 0)
        {
            scalarField& targetPatchLayers =
                targetLayers.boundaryFieldRef()[patchI];
            scalarField& layerPatchCounter =
                layerCount.boundaryFieldRef()[patchI];
            scalarField& actualPatchLayers =
                actualLayers.boundaryFieldRef()[patchI];
            const polyPatch& patch = patches[patchI];
            label patchStart = patch.start();
            forAll(patch, i)
            {
                if (layerPatchCounter[i] <= targetPatchLayers[i])
                {
                    layerPatchCounter[i] += 1;
                    if (layerPatchCounter[i] == actualPatchLayers[i])
                    {
                        label meshFaceI = patchStart + i;
                        const face& f = mesh.faces()[meshFaceI];
                        bool collapsePt = false;
                        forAll(f,fp)
                        {
                            if (collapsedPts[f[fp]])
                            {
                                collapsePt = true;
                                break;
                            }
                        }
                        if (!collapsePt)
                        {
                            markedFaces[meshFaceI] = true;
                            layerFaces.append(meshFaceI);
                        }
                    }
                }
            }
        }
    }

    boolList markedPts(mesh.nPoints(), false);
    forAll(layerFaces, lFI)
    {
        label meshFaceI = layerFaces[lFI];
        face f = mesh.faces()[meshFaceI];
        forAll(f, fp)
        {
            markedPts[f[fp]] = true;
        }
    }
    syncTools::syncPointList(mesh, markedPts, orEqOp<bool>(), false);

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            const polyPatch& patch = patches[patchI];
            label patchStart = patches[patchI].start();
            forAll(patch, i)
            {
                label meshFaceI = patchStart + i;
                if (markedFaces[meshFaceI])
                {
                    continue;
                }

                face f = mesh.faces()[meshFaceI];
                bool allCut = true;
                forAll(f,fp)
                {
                    if (!markedPts[f[fp]])
                    {
                        allCut = false;
                        break;
                    }
                }
                if (allCut)
                {
                    layerFaces.append(meshFaceI);
                }
            }
        }
    }

    return layerFaces.shrink();
}


void Foam::addMultiLayers::markCellCuttingEdges
(
    const label& celli,
    boolList& cutEdges
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    const cell& c = mesh.cells()[celli];

    DynamicList<label> triFaces(c.size());
    label cLevel = cellLevel[celli];

    forAll(c, cFI)
    {
        label facei = c[cFI];
        const face& f = mesh.faces()[facei];

        label nAnchors = 0;
        forAll(f,fp)
        {
            if (pointLevel[f[fp]] <= cLevel)
            {
                nAnchors++;
            }
        }
        if (nAnchors == 3)
        {
            triFaces.append(facei);
        }
    }

    label nCuts = triFaces.size()/2;
    if (nCuts > 0 && nCuts < 5)
    {
        const labelList& cEdges = mesh.cellEdges()[celli];
        const labelList& cPts = mesh.cellPoints()[celli];
        labelHashSet cFaces(c);
        labelHashSet triFaceSet(triFaces);

        DynamicList<label> innerTriEdges(cEdges.size());
        DynamicList<label> outerTriEdges(cEdges.size());
        DynamicList<label> outerTriPts(cPts.size());
        DynamicList<label> outerSingleTriPts(cPts.size());

        forAll(cEdges, cEI)
        {
            label edgei = cEdges[cEI];
            const labelList& eFaces = mesh.edgeFaces()[edgei];
            label nEdgeTriFaces = 0;

            forAll(eFaces, eFI)
            {
                label facei = eFaces[eFI];
                if (triFaceSet.found(facei) && cFaces.found(facei))
                {
                    nEdgeTriFaces++;
                }
            }

            if (nEdgeTriFaces == 1)
            {
                outerTriEdges.append(edgei);
            }
            else if (nEdgeTriFaces == 2)
            {
                innerTriEdges.append(edgei);
            }
        }

        forAll(outerTriEdges,oTEI)
        {
            edge e = mesh.edges()[outerTriEdges[oTEI]];

            forAll(e, ei)
            {
                bool addPt = true;
                forAll(outerTriPts, pti)
                {
                    if (outerTriPts[pti] == e[ei])
                    {
                        addPt = false;
                        break;
                    }
                }
                if (addPt)
                {
                    outerTriPts.append(e[ei]);
                }
            }
        }

        forAll(cPts, cPI)
        {
            label pointi = cPts[cPI];
            const labelList& pFaces = mesh.pointFaces()[pointi];
            label nPointTriFaces = 0;

            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                if (triFaceSet.found(facei) && cFaces.found(facei))
                {
                    nPointTriFaces++;
                }
            }

            if (nPointTriFaces == 1  && pointLevel[pointi] <= cLevel)
            {
                outerSingleTriPts.append(pointi);
            }
        }

        DynamicList<label>  exEdges(innerTriEdges.size());
        if (innerTriEdges.size() != nCuts)
        {
            //remove inavlid internal edges
            forAll(outerSingleTriPts, oSTI)
            {
                label pointi = outerSingleTriPts[oSTI];
                const labelList& pFaces = mesh.pointFaces()[pointi];
                forAll(pFaces, pFI)
                {
                    label facei = pFaces[pFI];
                    if
                    (
                        triFaceSet.found(facei) && cFaces.found(facei)
                    )
                    {
                        const labelList& fEdges = mesh.faceEdges()[facei];
                        label intEdge = -1;
                        forAll(fEdges, fEI)
                        {
                            label edgei = fEdges[fEI];
                            if (innerTriEdges.found(edgei))
                            {
                                intEdge = edgei;
                                break;
                            }
                        }
                        if (intEdge != -1)
                        {
                            const labelList& eFaces =
                                mesh.edgeFaces()[intEdge];
                            forAll(eFaces, eFI)
                            {
                                label otherFace = eFaces[eFI];
                                if
                                (
                                    otherFace != facei
                                    && cFaces.found(otherFace)
                                )
                                {
                                    const labelList& nbrEdges =
                                        mesh.faceEdges()[otherFace];
                                    forAll(nbrEdges, nEI)
                                    {
                                        label nbrEdge = nbrEdges[nEI];
                                        if (nbrEdge != intEdge)
                                        {
                                            bool addEdge = true;
                                            forAll(exEdges, excI)
                                            {
                                                if
                                                (
                                                    exEdges[excI] == nbrEdge
                                                )
                                                {
                                                    addEdge = false;
                                                    break;
                                                }
                                            }
                                            if (addEdge)
                                            {
                                                exEdges.append(nbrEdge);
                                            }
                                        }
                                    }
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            labelHashSet cEdgesSet(cEdges);
            labelHashSet outerTriPtSet(outerTriPts);
            forAll(innerTriEdges, iEI)
            {
                label edgei = innerTriEdges[iEI];
                const edge& e = mesh.edges()[edgei];

                if
                (
                    nCuts == 4 && !outerTriPtSet.found(e[0])
                    && !outerTriPtSet.found(e[1])
                )
                {
                    bool addEdge = true;
                    forAll(exEdges, excI)
                    {
                        if (exEdges[excI] == edgei)
                        {
                            addEdge = false;
                            break;
                        }
                    }
                    if (addEdge)
                    {
                        exEdges.append(edgei);
                    }
                }

                forAll(e, eI)
                {
                    const labelList& pEdges = mesh.pointEdges()[e[eI]];
                    label nPtEdges = 0;
                    forAll(pEdges, pEI)
                    {
                        if (cEdgesSet.found(pEdges[pEI]))
                        {
                            nPtEdges++;
                        }
                    }

                    bool hangingPt = true;

                    if
                    (
                        pointLevel[e[eI]] <= cLevel
                    )
                    {
                        hangingPt = false;
                    }

                    if (nPtEdges == 3 || hangingPt)
                    {
                        bool addEdge = true;
                        forAll(exEdges, excI)
                        {
                            if (exEdges[excI] == edgei)
                            {
                                addEdge = false;
                                break;
                            }
                        }
                        if (addEdge)
                        {
                            exEdges.append(edgei);
                        }
                        break;
                    }
                }
            }
        }

        DynamicList<label> selectedEdges(innerTriEdges.size());
        labelHashSet excludeSet(exEdges);

        forAll(innerTriEdges, iTE)
        {
            label edgei = innerTriEdges[iTE];
            if (!excludeSet.found(edgei))
            {
                selectedEdges.append(edgei);
            }
        }

        if (selectedEdges.size() == nCuts)
        {
            forAll(selectedEdges, cEI)
            {
                cutEdges[selectedEdges[cEI]] = true;
            }
        }
    }
    return;
}


void Foam::addMultiLayers::handleWrongOrientedFaces()
{
    if
    (
        controller_.algorithm() != meshControl::EXTRUDE
    )
    {
        return;
    }

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& owners = mesh.faceOwner();
    const labelList& neighbors = mesh.faceNeighbour();

    const volScalarField& layerCells = const_cast<volScalarField&>
        (mesh.lookupObject<volScalarField>("layerStacks"));

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+1);
    mesh.checkFacePyramids(true, -SMALL, &errorFaces);

    scalarField neiLayerCells;
    syncTools::swapBoundaryCellList(mesh, layerCells, neiLayerCells);
    pointField neiCC;
    syncTools::swapBoundaryCellList(mesh, mesh.cellCentres(), neiCC);

    scalarField weightSum(mesh.nPoints(), 0);
    pointField newPoints(mesh.nPoints(), vector::zero);

    label nMoved = 0;
    forAll(mesh.faces(), facei)
    {
        label own = owners[facei];
        scalar ownLayerCell = layerCells[owners[facei]];

        label patchI = -1;
        label checkCell = -1;
        point nCC = vector::zero;
        point oCC = mesh.cellCentres()[own];
        if (mesh.isInternalFace(facei))
        {
            label nei = neighbors[facei];
            label neiLayerCell = layerCells[nei];
            nCC = mesh.cellCentres()[nei];
            if (ownLayerCell > -1 && neiLayerCell < 0)
            {
                checkCell = own;
            }
            else if (ownLayerCell < 0 && neiLayerCell > -1)
            {
                checkCell = nei;
            }
        }
        else
        {
            patchI = patches.whichPatch(facei);
            if (patches[patchI].coupled())
            {
                label bFaceI = facei-mesh.nInternalFaces();
                label neiLayerCell = neiLayerCells[bFaceI];
                if (ownLayerCell > -1 && neiLayerCell < 0)
                {
                    checkCell = own;
                    nCC = neiCC[bFaceI];
                }
            }
        }

        if (checkCell != -1)
        {
            const cell& c = mesh.cells()[checkCell];
            forAll(c, cfi)
            {
                label facej = c[cfi];
                if (facei != facej && errorFaces.found(facej))
                {
                    face fj = mesh.faces()[facej];
                    DynamicList<label> anchorPts(fj.size());
                    {
                        label cLevel = cellLevel[checkCell];
                        forAll(fj,fp)
                        {
                            label pointi = fj[fp];
                            label pLevel = pointLevel[pointi];
                            if (pLevel <= cLevel)
                            {
                                anchorPts.append(pointi);
                            }
                        }
                    }

                    if (anchorPts.size() == 3)
                    {
                        point midPoint = 0.5*(nCC + oCC);
                        point fA = mesh.faceAreas()[facei];
                        scalar area = mag(fA);
                        if (area > SMALL)
                        {
                            nMoved++;
                            fA /= area;
                            plane fPlane(midPoint, fA);
                            face f = mesh.faces()[facei];
                            forAll(f, fp)
                            {
                                label pointi = f[fp];
                                scalar weight = scalar(1)/area;
                                point pt = mesh.points()[pointi];

                                vector disp = fPlane.nearestPoint(pt) - pt;
                                weightSum[pointi] += weight;
                                newPoints[pointi] += weight*disp;
                            }
                        }
                    }
                }
            }
        }
    }

    label totMoved = returnReduce(nMoved, sumOp<label>());
    if (totMoved > 0)
    {
        Info<<"Moving "<<totMoved<<" inverted faces neighboring points"<<endl;

        syncTools::syncPointList
        (
            mesh,
            weightSum,
            plusEqOp<scalar>(),
            scalar(0)
        );

        syncTools::syncPointList
        (
            mesh,
            newPoints,
            plusEqOp<point>(),
            vector::zero
         );

        forAll(newPoints, pointi)
        {
            point pt = mesh.points()[pointi];

            if (weightSum[pointi] > SMALL)
            {
                newPoints[pointi] = pt + (newPoints[pointi] / weightSum[pointi]);
            }
            else
            {
                newPoints[pointi] = pt;
            }
        }
        mesh.movePoints(newPoints);
    }
}


void Foam::addMultiLayers::remergeOuterSplitCells()
{
    if
    (
        !layerParams_.mergeSplitCells()
        || controller_.algorithm() != meshControl::EXTRUDE
    )
    {
        return;
    }

    fvMesh& mesh = meshRefiner_.mesh();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbors = mesh.faceNeighbour();

    volScalarField& layerCells = const_cast<volScalarField&>
        (mesh.lookupObject<volScalarField>("layerStacks"));

    boolList blockedFaces(mesh.nFaces(), true);

    //Identify split cells form outer non-layer cells
    {
        scalarField neiLayerCells;
        syncTools::swapBoundaryCellList(mesh, layerCells, neiLayerCells);

        boolList cutEdges(mesh.nEdges(), false);

        forAll(mesh.cells(), celli)
        {
            const cell& c = mesh.cells()[celli];

            if (layerCells[celli] > -1)
            {
                continue;
            }

            bool interfaceCell = false;
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label neiLayerCell = -1;
                if (mesh.isInternalFace(facei))
                {
                    label nei
                    (
                        owners[facei] == celli ?
                        neighbors[facei] : owners[facei]
                    );
                    neiLayerCell = layerCells[nei];
                }
                else
                {
                    label patchi = patches.whichPatch(facei);
                    if (patches[patchi].coupled())
                    {
                        label bFaceI = facei-mesh.nInternalFaces();
                        neiLayerCell = neiLayerCells[bFaceI];
                    }
                }
                if (neiLayerCell  > -1)
                {
                    interfaceCell = true;
                    break;
                }
            }

            if (interfaceCell)
            {
                markCellCuttingEdges(celli,cutEdges);
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            cutEdges,
            orEqOp<bool>(),
            false               // null value
        );

        boolList boundaryCells(mesh.nCells(), false);
        forAll(mesh.cells(), celli)
        {
            const cell& c = mesh.cells()[celli];
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchi = patches.whichPatch(facei);

                if (patchi != -1 && !patches[patchi].coupled())
                {
                    boundaryCells[celli] = true;
                    break;
                }
            }
        }
        boolList neiBoundaryCells;
        syncTools::swapBoundaryCellList(mesh, boundaryCells,neiBoundaryCells);

        forAll(cutEdges, edgei)
        {
            if (cutEdges[edgei])
            {
                const labelList& eFaces = mesh.edgeFaces()[edgei];

                forAll(eFaces, eFI)
                {
                    label facei = eFaces[eFI];

                    label own = owners[facei];

                    scalar ownLayerCell = layerCells[own];
                    scalar neiLayerCell = -1;

                    bool ownBoundaryCell = boundaryCells[own];
                    bool neiBoundaryCell = false;

                    if (mesh.isInternalFace(facei))
                    {
                        label nei = neighbors[facei];
                        neiLayerCell = layerCells[nei];
                        neiBoundaryCell = boundaryCells[nei];
                    }
                    else
                    {
                        label patchI = patches.whichPatch(facei);
                        if (patches[patchI].coupled())
                        {
                            label bFaceI = facei-mesh.nInternalFaces();
                            neiLayerCell = neiLayerCells[bFaceI];
                            neiBoundaryCell = neiBoundaryCells[bFaceI];
                        }
                    }

                    if
                    (
                        ownLayerCell != -1 && neiLayerCell != -1
                        && !ownBoundaryCell && !neiBoundaryCell
                    )
                    {
                        blockedFaces[facei] = false;
                    }
                }
            }
        }
    }

    if (Pstream::parRun())
    {
        syncTools::syncFaceList
        (
            mesh,
            blockedFaces,
            andEqOp<bool>()     // combine operator
        );

        autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
        (
            blockedFaces,
            scalarField(mesh.nCells(), 1), // dummy weights
            decomposer_,
            distributor_,
            false
        );

        map().distributeFaceData(blockedFaces);
    }

    //Merge cells
    {
        labelList nLayerCells(mesh.nEdges(), 0);
        labelList nNonLayerCells(mesh.nEdges(), 0);

        forAll(mesh.edges(), edgei)
        {
            const labelList& eCells = mesh.edgeCells()[edgei];
            forAll(eCells, eci)
            {
                label celli = eCells[eci];
                if (layerCells()[celli] != -1)
                {
                    nLayerCells[edgei]++;
                }
                else
                {
                    nNonLayerCells[edgei]++;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            nLayerCells,
            plusEqOp<label>(),
            label(0)               // null value
        );

        syncTools::syncEdgeList
        (
            mesh,
            nNonLayerCells,
            plusEqOp<label>(),
            label(0)               // null value
        );

        scalarField neiLayerCells;
        syncTools::swapBoundaryCellList(mesh, layerCells, neiLayerCells);

        polyTopoChange meshMod(mesh);

        boolList checkedFaces(mesh.nFaces(), false);
        labelList keepCell(mesh.nCells(), -1);
        labelList removeCell(mesh.nCells(), -1);

        boolList topEdges(mesh.nEdges(), false);

        forAll(mesh.faces(), facei)
        {
            if (!blockedFaces[facei])
            {
                const labelList& fEdges = mesh.faceEdges()[facei];
                forAll(fEdges, fEI)
                {
                    label meshedgei = fEdges[fEI];
                    if
                    (
                        nLayerCells[meshedgei] == 2
                        && nNonLayerCells[meshedgei] == 1
                    )
                    {
                        topEdges[meshedgei] = true;
                    }
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            topEdges,
            orEqOp<bool>(),
            false  // null value
        );

        label nMerged = 0;
        forAll(topEdges, edgei)
        {
            if (topEdges[edgei])
            {
                const edge& e = mesh.edges()[edgei];
                const labelList& eFaces = mesh.edgeFaces()[edgei];
                DynamicList<label> mergeTopFaces(2);

                forAll(eFaces, eFI)
                {
                    label facei = eFaces[eFI];

                    scalar ownLayerCell = layerCells[owners[facei]];
                    scalar neiLayerCell = -1;

                    label patchI = -1;

                    if (mesh.isInternalFace(facei))
                    {
                        neiLayerCell = layerCells[neighbors[facei]];
                    }
                    else
                    {
                        patchI = patches.whichPatch(facei);
                        if (patches[patchI].coupled())
                        {
                            label bFaceI = facei-mesh.nInternalFaces();
                            neiLayerCell = neiLayerCells[bFaceI];
                        }
                    }

                    if (ownLayerCell != -1 && neiLayerCell != -1)
                    {
                        // Remove face and cell
                        nMerged++;
                        label own = owners[facei];
                        label nei = neighbors[facei];

                        meshMod.setAction(polyRemoveCell(own));
                        meshMod.setAction(polyRemoveFace(facei));
                        checkedFaces[facei] = true;

                        keepCell[nei] = own;
                        removeCell[own] = nei;
                    }
                    else
                    {
                        mergeTopFaces.append(facei);
                    }
                }

                if (mergeTopFaces.size() == 2)
                {
                    //check if first merge faces associated with kept
                    //cell otherwise swap
                    {
                        label facei = mergeTopFaces[0];
                        bool swapAnchor(false);
                        label own = mesh.faceOwner()[facei];

                        if (mesh.isInternalFace(facei))
                        {
                            label nei = neighbors[facei];

                            if (keepCell[own] == -1 && keepCell[nei] == -1)
                            {
                                swapAnchor = true;
                            }
                        }
                        else
                        {
                            if (keepCell[own] == -1)
                            {
                                swapAnchor = true;
                            }
                        }

                        if (swapAnchor)
                        {
                            label tmpFace = mergeTopFaces[0];
                            mergeTopFaces[0] = mergeTopFaces[1];
                            mergeTopFaces[1] = tmpFace;
                        }
                    }

                    face f0 = mesh.faces()[mergeTopFaces[0]];
                    face f1 = mesh.faces()[mergeTopFaces[1]];
                    face newFace(f0.size() + f1.size() -2);

                    label start = findIndex(f0, e[0]);
                    label end = findIndex(f0, e[1]);
                    if (f0[f0.fcIndex(start)] == e[1])
                    {
                        label oldStart = start;
                        start = end;
                        end = oldStart;
                    }

                    label next = start;
                    forAll(f0, fp)
                    {
                        newFace[fp] = f0[next];
                        next = f0.fcIndex(next);
                    }

                    label startPt = f0[start];
                    label endPt = f0[end];

                    start = findIndex(f1, endPt);
                    bool reverse = false;
                    if (f1[f1.fcIndex(start)] == startPt)
                    {
                        reverse = true;
                    }

                    forAll(f1, fp)
                    {
                        if (reverse)
                        {
                            start = f1.rcIndex(start);
                        }
                        else
                        {
                            start = f1.fcIndex(start);
                        }

                        if (f1[start] == startPt)
                        {
                            break;
                        }
                        else
                        {
                            newFace[f0.size()+fp] = f1[start];
                        }
                    }

                    meshMod.setAction(polyRemoveFace(mergeTopFaces[1]));
                    checkedFaces[mergeTopFaces[0]] = true;
                    checkedFaces[mergeTopFaces[1]] = true;
                    label modFace = mergeTopFaces[0];

                    label zoneID = mesh.faceZones().whichZone(modFace);
                    bool zoneFlip = false;
                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = mesh.faceZones()[zoneID];
                        zoneFlip = fZone.flipMap()[fZone.whichFace(modFace)];
                    }

                    label own = owners[modFace];

                    if (removeCell[own] != -1)
                    {
                        WarningInFunction
                            <<"Split edge merge, Owner cell: "<< own
                            <<" has been removed from mesh "
                            <<"kept cell is "<< removeCell[own]<<endl;
                    }

                    if (mesh.isInternalFace(modFace))
                    {
                        label nei = neighbors[modFace];

                        if (removeCell[nei] != -1)
                        {
                            WarningInFunction
                                <<"Split edge merge, Neighbour cell: "<< nei
                                <<" has been removed from mesh "
                                <<" kept cell is "<< removeCell[nei]<<endl;
                        }

                        // Modify the master face.
                        meshMod.setAction
                        (
                            polyModifyFace
                            (
                                newFace,       // original face
                                modFace,    // label of face
                                own,            // owner
                                nei,             // neighbour
                                false,          // face flip
                                -1,         // patch for face
                                false,          // remove from zone
                                zoneID,         // zone for face
                                zoneFlip        // face flip in zone
                             )
                         );
                    }
                    else
                    {
                        label patchI = mesh.boundaryMesh().whichPatch(modFace);

                        // Modify the master face.
                        meshMod.setAction
                        (
                            polyModifyFace
                            (
                                newFace,       // original face
                                modFace,    // label of face
                                own,            // owner
                                -1,             // neighbour
                                false,          // face flip
                                patchI,         // patch for face
                                false,          // remove from zone
                                zoneID,         // zone for face
                                zoneFlip        // face flip in zone
                             )
                         );
                    }
                }
            }
        }

        forAll(keepCell, cellI)
        {
            if (keepCell[cellI] != -1)
            {
                label removedCell = keepCell[cellI];
                const labelList cFaces = mesh.cells()[removedCell];
                forAll(cFaces, cFI)
                {
                    label facei = cFaces[cFI];
                    if (!checkedFaces[facei])
                    {
                        checkedFaces[facei] = true;
                        label zoneID = mesh.faceZones().whichZone(facei);
                        bool zoneFlip = false;
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone =
                                mesh.faceZones()[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                        }

                        if (mesh.isInternalFace(facei))
                        {
                            label own = cellI;
                            label nei = -1;
                            bool flipFace = false;
                            if (mesh.faceOwner()[facei] == removedCell)
                            {
                                nei = neighbors[facei];

                                if (removeCell[nei] != -1)
                                {
                                    nei = removeCell[nei];
                                }

                                if (own > nei)
                                {
                                    own = nei;
                                    nei = cellI;
                                    flipFace = true;
                                }
                            }
                            else
                            {
                                nei = owners[facei];

                                if (removeCell[nei] != -1)
                                {
                                    nei = removeCell[nei];
                                }

                                if (own > nei)
                                {
                                    own = nei;
                                    nei = cellI;
                                }
                                else
                                {
                                    flipFace = true;
                                }
                            }

                            face newFace(mesh.faces()[facei]);
                            if (flipFace)
                            {
                                newFace.flip();
                            }

                            if (removeCell[own] != -1)
                            {
                                WarningInFunction
                                    <<"Deleted cell faces update, Owner cell: "
                                    << own << " has been removed from mesh "
                                    <<" kept cell is "<< removeCell[own]<<endl;
                            }

                            if (removeCell[nei] != -1)
                            {
                                WarningInFunction
                                    <<"Deleted cell faces update, "
                                    <<"Neighbour cell: "<< nei
                                    <<" has been removed from mesh "
                                    <<" kept cell is "<< removeCell[nei]<<endl;
                            }

                            // Modify the master face.
                            meshMod.setAction
                            (
                                polyModifyFace
                                (
                                    newFace,       // original face
                                    facei,    // label of face
                                    own,            // owner
                                    nei,             // neighbour
                                    false,          // face flip
                                    -1,         // patch for face
                                    false,          // remove from zone
                                    zoneID,         // zone for face
                                    zoneFlip        // face flip in zone
                                )
                             );
                        }
                        else
                        {
                            label patchI =
                                mesh.boundaryMesh().whichPatch(facei);

                            // Modify the master face.
                            meshMod.setAction
                            (
                                polyModifyFace
                                (
                                    mesh.faces()[facei], // original face
                                    facei,    // label of face
                                    cellI,            // owner
                                    -1,             // neighbour
                                    false,          // face flip
                                    patchI,         // patch for face
                                    false,          // remove from zone
                                    zoneID,         // zone for face
                                    zoneFlip        // face flip in zone
                                 )
                            );
                        }
                    }
                }
            }
        }

        label totMerged = returnReduce(nMerged, sumOp<label>());
        if (totMerged > 0)
        {
            Info<<"Merging "<<totMerged<<" final layer cells"<<endl;

            // Change the mesh (no inflation)
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

            // Update fields (problem when face added to zero sized patch)
            mesh.updateMesh(map);

            // Move mesh if in inflation mode
            if (map().hasMotionPoints())
            {
                mesh.movePoints(map().preMotionPoints());
            }
            else
            {
                mesh.clearOut();
            }

            meshRefiner_.updateMesh(map,labelList(0));
        }
    }

    return;
}


void Foam::addMultiLayers::reorientateFaces()
{
    fvMesh& mesh = meshRefiner_.mesh();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    scalarField openness;
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const vectorField& areas = mesh.faceAreas();

    vectorField sumClosed(mesh.nCells(), vector::zero);
    vectorField sumMagClosed(mesh.nCells(), vector::zero);

    forAll(own, faceI)
    {
        // Add to owner
        sumClosed[own[faceI]] += areas[faceI];
        sumMagClosed[own[faceI]] += cmptMag(areas[faceI]);
    }

    forAll(nei, faceI)
    {
        // Subtract from neighbour
        sumClosed[nei[faceI]] -= areas[faceI];
        sumMagClosed[nei[faceI]] += cmptMag(areas[faceI]);
    }
    boolList openCells(mesh.nCells(), false);
    label nOpen = 0;

    forAll(sumClosed, cellI)
    {
        scalar maxOpenness = 0;

        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            maxOpenness = max
            (
                maxOpenness,
                mag(sumClosed[cellI][cmpt])
                /(sumMagClosed[cellI][cmpt] + ROOTVSMALL)
            );
        }

        if (maxOpenness > 1e-6)
        {
            openCells[cellI] = true;
            nOpen++;
        }
    }

    reduce(nOpen, sumOp<label>());
    if (nOpen > 1)
    {
        Info<<"Re-orientate wrongly orientated faces"<<endl;
        polyTopoChange meshMod(mesh);
        volScalarField& layerCells = const_cast<volScalarField&>
            (mesh.lookupObject<volScalarField>("layerStacks"));

        boolList neiOpenCells;
        syncTools::swapBoundaryCellList(mesh, openCells, neiOpenCells);
        scalarField neiLayerCells;
        syncTools::swapBoundaryCellList(mesh, layerCells, neiLayerCells);

        forAll(mesh.faces(), faceI)
        {
            bool ownOpen = openCells[own[faceI]];
            scalar ownLayerCell = layerCells[own[faceI]];
            bool neiOpen = false;
            scalar neiLayerCell = -1;

            label patchI = -1;

            if (mesh.isInternalFace(faceI))
            {
                neiOpen = openCells[nei[faceI]];
                neiLayerCell = layerCells[nei[faceI]];
            }
            else
            {
                patchI = patches.whichPatch(faceI);
                if (patches[patchI].coupled())
                {
                    label bFaceI = faceI-mesh.nInternalFaces();
                    neiOpen = neiOpenCells[bFaceI];
                    neiLayerCell = neiLayerCells[bFaceI];
                }
            }

            if
            (
                neiOpen && ownOpen
                && ownLayerCell != -1 && neiLayerCell != -1
                && (ownLayerCell == neiLayerCell)
            )
            {
                label zoneID = mesh.faceZones().whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }

                face modFace = mesh.faces()[faceI];
                modFace.flip();

                meshMod.modifyFace
                (
                    modFace, // modified face
                    faceI,   // label of face being modified
                    own[faceI], // owner
                    (patchI == -1 ? nei[faceI] : -1), // neighbour
                    false,    // face flip
                    patchI,   // new patch for face
                    zoneID,   // zone for face
                    zoneFlip  // face flip in zone
                );
            }
        }


        // Create mesh (no inflation), return map from old to new mesh.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map);

        // Optionally inflate mesh
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh.setInstance(meshRefiner_.timeName());

        // Update intersection info
        meshRefiner_.updateMesh(map, labelList(0));
    }
}


//convert baffle zones into patches
void Foam::addMultiLayers::convertBaffleZones()
{
    fvMesh& mesh = meshRefiner_.mesh();

    List<labelPair> baffles;
    labelList originatingFaceZone;

    List<surfaceZonesInfo::faceZoneType> fzType(1);
    fzType[0] = surfaceZonesInfo::BAFFLE;

    const labelList zonesToBaffle
    (
        meshRefiner_.getZones(fzType)
    );

    if (zonesToBaffle.size() > 0)
    {
        scalarField weights(mesh.nCells(), 1.0);

        autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
        (
            false,
            true,
            false,
            weights,
            decomposer_,
            distributor_,
            false
        );

        meshRefiner_.createZoneBaffles
        (
            zonesToBaffle,
            baffles,
            originatingFaceZone,
            false
        );
    }
}


void Foam::addMultiLayers::updateExtrudeMesh
(
    const labelList& layerPatchIDs,
    const labelList& ppLayers,
    const labelList& ppPtLayers,
    indirectPrimitivePatch& pp,
    labelList& ppLayersGlobal,
    labelList& ppPtLayersGlobal
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& numLayers = layerParams_.numLayers();
    const bool truncateFromWall = layerParams_.truncateFromWall();

    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    label maxNumLayers = gMax(ppPtLayers);

    labelList nPPMeshEdgeLayers(mesh.nEdges(), 0);
    forAll(pp.edges(), edgei)
    {
        label meshedgei = meshEdges[edgei];
        const labelList& eFaces = pp.edgeFaces()[edgei];
        forAll(eFaces, eFI)
        {
            label facei = eFaces[eFI];
            nPPMeshEdgeLayers[meshedgei] = max
            (
                nPPMeshEdgeLayers[meshedgei],
                ppLayers[facei]
            );
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        nPPMeshEdgeLayers,
        maxEqOp<label>(),
        labelMin // initial value
    );

    labelList nPPEdgeLayer(pp.nEdges(), 0);
    forAll(pp.edges(), edgei)
    {
        label meshedgei = meshEdges[edgei];
        nPPEdgeLayer[edgei] = nPPMeshEdgeLayers[meshedgei];
    }

    polyTopoChange meshMod(mesh);

    //Add new stack cells
    List<DynamicList<label>> addedCells
        (pp.size(), DynamicList<label>(maxNumLayers));

    label nCellsSplit = 0;
    forAll(pp, i)
    {
        label nLay = ppLayers[i];
        if (nLay > 1)
        {
            label facei = pp.addressing()[i];
            label own = mesh.faceOwner()[facei];
            label ownZoneI = mesh.cellZones().whichZone(own);
            for (label stackI=0; stackI < nLay-1; ++stackI)
            {
                label newCellI = meshMod.setAction
                (
                    polyAddCell
                    (
                        -1,             // master point
                        -1,             // master edge
                        -1,             // master face
                        own,            //master
                        ownZoneI        // zone for cell
                     )
                );
                addedCells[i].append(newCellI);
                nCellsSplit++;
            }
        }
    }

    //Add stack points
    boolList boundaryPoints(mesh.nPoints(), false);
    forAll(pp.localPoints(), i)
    {
        boundaryPoints[pp.meshPoints()[i]] = true;
    }

    syncTools::syncPointList
    (
        mesh,
        boundaryPoints,
        orEqOp<bool>(),
        false
    );

    boolList stackEdges(mesh.nEdges(), false);
    forAll(mesh.edges(), edgei)
    {
        const edge& e = mesh.edges()[edgei];
        if
        (
            (boundaryPoints[e[0]] && !boundaryPoints[e[1]])
            || (boundaryPoints[e[1]] && !boundaryPoints[e[0]])
        )
        {
            stackEdges[edgei] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        stackEdges,
        orEqOp<bool>(),
        false
    );

    List<DynamicList<label>> stackPts
    (
        pp.localPoints().size(),
        DynamicList<label>(maxNumLayers+1)
    );

    forAll(pp.localPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];
        label zoneI = mesh.pointZones().whichZone(meshPointI);

        const labelList& pEdges = mesh.pointEdges()[meshPointI];
        stackPts[ptI].append(meshPointI);
        forAll(pEdges, pEI)
        {
            label edgei = pEdges[pEI];
            if (stackEdges[edgei])
            {
                edge e = mesh.edges()[edgei];

                label nStackPts = ppPtLayers[ptI];
                label otherPt(e[0] == meshPointI ? e[1] : e[0]);
                if (nStackPts > 1)
                {
                    point innerPt = mesh.points()[meshPointI];
                    point outerPt = mesh.points()[otherPt];
                    vector stackDir = outerPt-innerPt;
                    for (label stackI=0; stackI < nStackPts-1; ++stackI)
                    {
                        scalar ratio = scalar(stackI+1)/scalar(nStackPts);
                        point pt = innerPt + (ratio*stackDir);
                        label newPointI = meshMod.setAction
                        (
                            polyAddPoint
                            (
                                pt,         // point
                                meshPointI,  // master point
                                zoneI,      // zone for point
                                true        // supports a cell
                             )
                        );
                        stackPts[ptI].append(newPointI);
                    }
                }
                stackPts[ptI].append(otherPt);
            }
        }
    }

    //Modify bottom faces
    forAll(pp,i)
    {
        if (ppLayers[i] > 1)
        {
            label facei = pp.addressing()[i];
            const face& f = mesh.faces()[facei];
            label zoneID = mesh.faceZones().whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }
            label patchi = patches.whichPatch(facei);
            label newCellI = addedCells[i][0];

            meshMod.setAction
            (
                polyModifyFace
                (
                    f,                // modified face
                    facei,            // label of face being modified
                    newCellI,         // owner
                    -1,               // neighbour
                    false,            // face flip
                    patchi,           // patch for face
                    false,            // remove from zone
                    zoneID,           // zone for face
                    zoneFlip          // face flip in zone
                 )
            );
        }
    }

    //Add new stack faces
    forAll(pp,i)
    {
        label nLay = ppLayers[i];
        if (nLay > 1)
        {
            for (label stackI=0; stackI < nLay-1; ++stackI)
            {
                const face& f = pp.localFaces()[i];
                const label facei = pp.addressing()[i];
                face newFace(f.size());
                forAll(f,fp)
                {
                    label ptI = f[fp];
                    if (truncateFromWall)
                    {
                        label offset = stackPts[ptI].size()-nLay;
                        newFace[fp] = stackPts[ptI][stackI+offset];
                    }
                    else
                    {
                        newFace[fp] = stackPts[ptI][stackI+1];
                    }
                }

                label own, nei;
                if (stackI == nLay-2)
                {
                    own = mesh.faceOwner()[facei];
                    nei = addedCells[i][stackI];
                }
                else
                {
                    own = addedCells[i][stackI];
                    nei = addedCells[i][stackI+1];
                    newFace.flip();
                }

                label zoneI = -1;
                bool flip = false;
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,   // face
                        own,       // owner
                        nei,       // neighbour
                        -1,        // master point
                        -1,        // master edge
                        facei,     // master face
                        false,     // flux flip
                        -1,        // patch for face
                        zoneI,     // zone for face
                        flip       // face zone flip
                     )
                 );
            }
        }
    }

    // Delete/Add new side faces
    boolList markedFaces(mesh.nFaces(), false);
    forAll(pp, i)
    {
        const face& f = pp.localFaces()[i];
        forAll(f,fp)
        {
            label next = f.fcIndex(fp);
            label startPt = f[fp];
            label nextPt = f[next];
            label edgei =  meshTools::findEdge
            (
                pp.edges(),
                pp.pointEdges()[startPt],
                startPt,
                nextPt
            );
            label nLay = nPPEdgeLayer[edgei];
            if
            (
                nLay > 1
                || stackPts[startPt].size() > 2
                || stackPts[nextPt].size() > 2
            )
            {
                label meshedgei = meshEdges[edgei];
                const labelList& eFaces = mesh.edgeFaces()[meshedgei];
                forAll(eFaces, eFI)
                {
                    label sfacei = eFaces[eFI];
                    if (!markedFaces[sfacei])
                    {
                        label patchi = patches.whichPatch(sfacei);
                        bool sideFace = false;
                        if (patchi == -1)
                        {
                            sideFace = true;
                        }
                        else if (patches[patchi].coupled())
                        {
                            sideFace = true;
                        }
                        else if (numLayers[patchi] < 0)
                        {
                            sideFace = true;
                        }

                        if (sideFace)
                        {
                            meshMod.setAction(polyRemoveFace(sfacei));
                            markedFaces[sfacei] = true;
                            label startOffset = stackPts[startPt].size()-nLay;
                            label nextOffset = stackPts[nextPt].size()-nLay;

                            for (label stackI=0; stackI < nLay; ++stackI)
                            {
                                DynamicList<label> facePts(4);

                                //check for top hanging points
                                label hangingPt = -1;
                                label mshPtStart =
                                    pp.meshPoints()[startPt];
                                label mshPtNext =
                                    pp.meshPoints()[nextPt];

                                const labelList& pEdgesStart =
                                    mesh.pointEdges()[mshPtStart];
                                const labelList& pEdgesNext =
                                    mesh.pointEdges()[mshPtNext];

                                label startEdge = -1;
                                forAll(pEdgesStart, pEI)
                                {
                                    label edgei = pEdgesStart[pEI];
                                    if (stackEdges[edgei])
                                    {
                                        startEdge = edgei;
                                    }
                                }

                                label nextEdge = -1;
                                forAll(pEdgesNext, pEI)
                                {
                                    label edgei = pEdgesNext[pEI];
                                    if (stackEdges[edgei])
                                    {
                                        nextEdge = edgei;
                                    }
                                }

                                edge se = mesh.edges()[startEdge];
                                label otherStart
                                (
                                    se[0] == mshPtStart
                                    ? se[1] : se[0]
                                );

                                edge ne = mesh.edges()[nextEdge];
                                label otherNext
                                (
                                    ne[0] == mshPtNext
                                    ? ne[1] : ne[0]
                                );

                                face sf = mesh.faces()[sfacei];
                                label startSide =
                                    findIndex(sf,mshPtNext);
                                if
                                (
                                    sf[sf.fcIndex(startSide)]
                                    != otherNext
                                )
                                {
                                    sf.flip();
                                }
                                startSide = findIndex(sf,otherNext);
                                forAll(sf, fp)
                                {
                                    label nextSide =
                                        sf.fcIndex(startSide);
                                    if (sf[nextSide] == otherStart)
                                    {
                                        break;
                                    }
                                    else
                                    {
                                        hangingPt = sf[nextSide];
                                    }
                                    startSide = nextSide;
                                }

                                if (truncateFromWall)
                                {
                                    if (stackI == 0)
                                    {
                                        for (label j=0; j < nextOffset+1; ++j)
                                        {
                                            facePts.append(stackPts[nextPt][j]);
                                        }

                                        if (nLay == 1 && hangingPt != -1)
                                        {
                                            facePts.append(hangingPt);
                                        }

                                        for (label j= startOffset; j >= 0; --j)
                                        {
                                            facePts.append(stackPts[startPt][j]);
                                        }
                                    }
                                    else
                                    {
                                        facePts.append
                                        (
                                           stackPts[startPt][startOffset+stackI-1]
                                        );
                                        facePts.append
                                        (
                                            stackPts[nextPt][nextOffset+stackI-1]
                                        );
                                        facePts.append
                                        (
                                            stackPts[nextPt][nextOffset+stackI]
                                        );

                                        if (stackI == nLay-1 && hangingPt != -1)
                                        {
                                            facePts.append(hangingPt);
                                        }

                                        facePts.append
                                        (
                                            stackPts[startPt][startOffset+stackI]
                                        );
                                    }
                                }
                                else
                                {
                                    if (stackI == nLay-1)
                                    {
                                        label sz = stackPts[nextPt].size();
                                        for (label j=stackI; j < sz; ++j)
                                        {
                                            facePts.append(stackPts[nextPt][j]);
                                        }

                                        if (hangingPt != -1)
                                        {
                                            facePts.append(hangingPt);
                                        }

                                        sz = stackPts[startPt].size();
                                        for (label j=sz-1; j >= stackI; --j)
                                        {
                                            facePts.append
                                            (
                                                stackPts[startPt][j]
                                            );
                                        }
                                    }
                                    else
                                    {
                                        facePts.append
                                        (
                                            stackPts[startPt][stackI]
                                        );
                                        facePts.append
                                        (
                                            stackPts[nextPt][stackI]
                                        );
                                        facePts.append
                                        (
                                            stackPts[nextPt][stackI+1]
                                        );
                                        facePts.append
                                        (
                                            stackPts[startPt][stackI+1]
                                        );
                                    }
                                }

                                face newFace(facePts);

                                label own = -1;
                                label nei = -1;

                                if (truncateFromWall)
                                {
                                    if
                                    (
                                        stackI == nLay-1
                                        || addedCells[i].size() == 0
                                    )
                                    {
                                        label oBFace = pp.addressing()[i];
                                        own = mesh.faceOwner()[oBFace];
                                    }
                                    else
                                    {
                                        label ownOffset = nLay
                                            - addedCells[i].size();
                                        if (stackI < ownOffset)
                                        {
                                            own = addedCells[i][0];
                                        }
                                        else
                                        {
                                            own =
                                                addedCells[i][stackI-ownOffset+1];
                                        }
                                    }
                                }
                                else
                                {
                                    if (stackI >= addedCells[i].size())
                                    {
                                        own =
                                            mesh.faceOwner()[pp.addressing()[i]];
                                    }
                                    else
                                    {
                                        own = addedCells[i][stackI];
                                    }
                                }

                                if (patchi == -1)
                                {
                                    const labelList& ledgeFaces =
                                        pp.edgeFaces()[edgei];
                                    label otherFace
                                    (
                                        ledgeFaces[0] == i
                                        ? ledgeFaces[1] : ledgeFaces[0]
                                    );

                                    if (truncateFromWall)
                                    {
                                        if
                                        (
                                            stackI == nLay-1
                                            || addedCells[otherFace].size() == 0
                                        )
                                        {
                                            label nBFace =
                                                pp.addressing()[otherFace];
                                            nei = mesh.faceOwner()[nBFace];
                                        }
                                        else
                                        {
                                            label neiOffset = nLay
                                                - addedCells[otherFace].size();
                                            if (stackI < neiOffset)
                                            {
                                                nei = addedCells[otherFace][0];
                                            }
                                            else
                                            {
                                                label nOff = stackI-neiOffset+1;
                                                nei = addedCells[otherFace][nOff];
                                            }
                                        }
                                    }
                                    else
                                    {
                                        if (stackI >= addedCells[otherFace].size())
                                        {
                                            nei =
                                                mesh.faceOwner()
                                                [pp.addressing()[otherFace]];
                                        }
                                        else
                                        {
                                            nei = addedCells[otherFace][stackI];
                                        }
                                    }

                                    if (nei < own)
                                    {
                                        label tnei = nei;
                                        nei = own;
                                        own = tnei;
                                    }
                                    else
                                    {
                                        newFace.flip();
                                    }
                                }
                                else
                                {
                                    newFace.flip();
                                }

                                label zoneID =
                                    mesh.faceZones().whichZone(sfacei);
                                bool zoneFlip = false;

                                meshMod.setAction
                                (
                                    polyAddFace
                                    (
                                        newFace,   // face
                                        own,       // owner
                                        nei,       // neighbour
                                        -1,        // master point
                                        -1,        // master edge
                                        sfacei,     // master face
                                        false,     // flux flip
                                        patchi,        // patch for face
                                        zoneID,     // zone for face
                                        zoneFlip       // face zone flip
                                     )
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    if (returnReduce(nCellsSplit, sumOp<label>()) != 0)
    {
        Map<label> mapMeshToPatch(pp.size());
        forAll(pp, facei)
        {
            mapMeshToPatch.insert(pp.addressing()[facei], facei);
        }

        Map<label> mapMeshPointToPatch(pp.nPoints());
        forAll(pp.meshPoints(), pointi)
        {
            mapMeshPointToPatch.insert(pp.meshPoints()[pointi], pointi);
        }

        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

        // Update fields (problem when face added to zero sized patch)
        mesh.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            mesh.clearOut();
        }

        meshRefiner_.updateMesh(map,labelList(0));

        //Reset primitive patch addressing
        pp.clearOut();
        labelList newAddressing;
        meshRefinement::calcPatchAddressing(mesh,layerPatchIDs,newAddressing);
        pp.resetAddressing(newAddressing);
        pp.localPoints();

        const labelList& pointMap = map().pointMap();
        const labelList& faceMap = map().faceMap();

        labelList newPPLayersGlobal(pp.size(), -1);
        labelList newPPPtLayersGlobal(pp.localPoints().size(), -1);
        forAll(pp, facei)
        {
            label meshfacei = pp.addressing()[facei];
            label origMeshFacei = faceMap[meshfacei];
            newPPLayersGlobal[facei] =
                ppLayersGlobal[mapMeshToPatch[origMeshFacei]];
        }

        forAll(pp.localPoints(),pointi)
        {
            label meshpointi = pp.meshPoints()[pointi];
            label origMeshPointi = pointMap[meshpointi];
            newPPPtLayersGlobal[pointi] =
                ppPtLayersGlobal[mapMeshPointToPatch[origMeshPointi]];
        }

        ppLayersGlobal = newPPLayersGlobal;
        ppPtLayersGlobal = newPPPtLayersGlobal;

        Info<< "Layers added. Number of cells : "
            << returnReduce(mesh.nCells(), sumOp<label>())
            << endl;
    }

    return;
}


void Foam::addMultiLayers::incrementLayers
(
    const indirectPrimitivePatch& pp,
    labelList& minPointLayers,
    labelList& maxPointLayers,
    labelList& ppLayers
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    //At patch transition from lower to higher number of layers
    // increment number of layers on lower layer count patch
    const bool incrementLower = layerParams_.incrementLower();

    bool firstPass = true;
    while (true)
    {
        minPointLayers = labelMax;
        maxPointLayers = labelMin;
        forAll(pp.meshPoints(), ptI)
        {
            const labelList& pFaces = pp.pointFaces()[ptI];
            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                label nLay =  ppLayers[facei];
                minPointLayers[ptI] = min(minPointLayers[ptI], nLay);
                maxPointLayers[ptI] = max(maxPointLayers[ptI], nLay);
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            minPointLayers,
            minEqOp<label>(),  // combine op
            labelMax            // null value
        );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            maxPointLayers,
            maxEqOp<label>(),  // combine op
            labelMin            // null value
        );

        label nFacesUpdated = 0;
        forAll(pp, facei)
        {
            if (ppLayers[facei] > 1)
            {
                const face& f = pp.localFaces()[facei];
                label minFaceLayers = labelMax;
                label maxFaceLayers = labelMin;
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    minFaceLayers = min(minFaceLayers, minPointLayers[pointi]);
                    maxFaceLayers = max(maxFaceLayers, maxPointLayers[pointi]);
                }

                label layerDif = maxFaceLayers - minFaceLayers;

                if (firstPass)
                {
                    if (incrementLower)
                    {
                        if (layerDif > 1 && ppLayers[facei] < maxFaceLayers-1)
                        {
                            ppLayers[facei]++;
                            nFacesUpdated++;
                        }
                    }
                    else
                    {
                        if (layerDif > 1 &&  ppLayers[facei] > minFaceLayers+1)
                        {
                            ppLayers[facei]--;
                            nFacesUpdated++;
                        }
                    }
                }
                else
                {
                    if (incrementLower)
                    {
                        if (layerDif > 1 &&  ppLayers[facei] > minFaceLayers+1)
                        {
                            ppLayers[facei]--;
                            nFacesUpdated++;
                        }
                    }
                    else
                    {
                        if (layerDif > 1 && ppLayers[facei] < maxFaceLayers-1)
                        {
                            ppLayers[facei]++;
                            nFacesUpdated++;
                        }
                    }
                }
            }
        }

        if (returnReduce(nFacesUpdated, sumOp<label>()) == 0)
        {
            if (firstPass)
            {
                firstPass = false;
            }
            else
            {
                break;
            }
        }
    }

    return;
}

void Foam::addMultiLayers::setRefinementExtrude()
{
    Info<< "Adding final layer cells: "<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& numLayers = layerParams_.numLayers();

    scalarField cellWeights(mesh.nCells(), 1);
    forAll(numLayers, patchi)
    {
        label nLayers = numLayers[patchi];
        if (nLayers > 1)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];
            forAll(pp.faceCells(), i)
            {
                cellWeights[pp.faceCells()[i]] += nLayers;
            }
        }
    }

    meshRefiner_.balance
    (
        false,
        false,
        false,
        cellWeights,
        decomposer_,
        distributor_,
        false
    );

    DynamicList<label> layerPatchIDs(patches.size());
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            if (numLayers[patchI] >= 0)
            {
                layerPatchIDs.append(patchI);
            }
        }
    }
    layerPatchIDs.shrink();

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            layerPatchIDs
         )
     );
    indirectPrimitivePatch& pp = ppPtr();

    if (returnReduce(pp.size(), sumOp<label>()) == 0)
    {
        Info<<"No layers to add." <<endl;
        //convert baffle zones to patches
        convertBaffleZones();
        return;
    }

    labelList ppLayers(pp.size(), 0);
    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        label patchI = patches.whichPatch(meshFaceI);
        label nLayers = numLayers[patchI];
        if (nLayers > -1)
        {
            ppLayers[i] = max(nLayers,label(1));
        }
    }

    //Check points to collapse to a single layer
    checkCollapses(pp, ppLayers);

    labelList minPointLayers(pp.meshPoints().size(), labelMax);
    labelList maxPointLayers(pp.meshPoints().size(), labelMin);

    //Adjust changing layers
    incrementLayers
    (
        pp,
        minPointLayers,
        maxPointLayers,
        ppLayers
    );

    label maxNumLayers = gMax(maxPointLayers);
    const label& nSmoothPerLayer = layerParams_.nSmoothPerDualLayer();
    label nSmooth = max((nSmoothPerLayer*maxNumLayers),label(10));

    //Perform two stage layer addition
    if (layerParams_.twoStageExtrusion())
    {
        bool adjustFinalNLayers = false;
        const List<Tuple2<scalar,scalar>>& targetExpansions =
           layerParams_.targetExpansions();
        const List<wordList>& layerSpec = layerParams_.layerSpec();
        forAll(targetExpansions, patchi)
        {
            const Tuple2<scalar,scalar>& targetExp = targetExpansions[patchi];
            word method = layerSpec[patchi][1];
            if
            (
                (method == "fch" || method == "rfch")
                && (targetExp.first() > -GREAT || targetExp.second() < GREAT)
            )
            {
                Info<< "Adjusting layers based on target stretch" <<endl;
                adjustFinalNLayers = true;
                break;
            }
        }
        labelList ppLayersLocal(pp.size(), -1);
        labelList ppPtLayersLocal(pp.meshPoints().size(), -1);
        forAll(pp, i)
        {
            label nLay = ppLayers[i];
            if (nLay > 2)
            {
                ppLayersLocal[i] = 2;
            }
            else
            {
                ppLayersLocal[i] = nLay;
            }
        }

        forAll(pp.meshPoints(), i)
        {
            label nLay = maxPointLayers[i];
            if (nLay > 2)
            {
                ppPtLayersLocal[i] = 2;
            }
            else
            {
                ppPtLayersLocal[i] = nLay;
            }
        }

        //Add layers
        updateExtrudeMesh
        (
            layerPatchIDs,
            ppLayersLocal,
            ppPtLayersLocal,
            pp,
            ppLayers,
            maxPointLayers
        );

        layerManipulate layerManip
        (
            mesh,
            layerParams_,
            meshRefiner_.meshCutter().cellLevel(),
            meshRefiner_.meshCutter().pointLevel(),
            meshRefiner_.meshCutter().level0EdgeLength()
        );

        ppPtLayersLocal = -1;
        forAll(pp.meshPoints(), i)
        {
            label nLay = maxPointLayers[i];
            if (nLay > 2)
            {
                ppPtLayersLocal[i] = nLay - 2;
            }
            else
            {
                ppPtLayersLocal[i] = 1;
            }
        }

        labelList adjustedNLayers = maxPointLayers;
        layerManip.smoothLayerStack
        (
            controller_,
            nSmooth,
            grownUpGeometryDict_,
            grownUpZoneGeometryDict_,
            ppPtLayersLocal,
            adjustedNLayers,
            adjustFinalNLayers
        );

        if (adjustFinalNLayers)
        {
            forAll(ppLayers, i)
            {
                if(ppLayers[i] > 1)
                {
                    ppLayers[i] = labelMax;
                }
            }
            forAll(pp, i)
            {
                const face& lf =  pp.localFaces()[i];
                forAll(lf, fp)
                {
                    ppLayers[i] = min(adjustedNLayers[lf[fp]],ppLayers[i]);
                }
            }

            labelList adjustedMinLayers(pp.meshPoints().size(), labelMax);
            labelList adjustedMaxLayers(pp.meshPoints().size(), labelMin);

            //Adjust changing layers
            incrementLayers
            (
                pp,
                adjustedMinLayers,
                adjustedMaxLayers,
                ppLayers
            );
            maxPointLayers = adjustedMaxLayers;
        }

        dictionary meshOptimDict = dictionary();
        if (!meshOptimDict.found("cfMeshOptimizeCoeffs"))
        {
            meshOptimDict.add
            (
                "cfMeshOptimizeCoeffs",
                dictionary(),
                true
             );
        }
        dictionary& coeffsDict =
            meshOptimDict.subDict("cfMeshOptimizeCoeffs");
        coeffsDict.add
        (
            "relaxedCheck",
            true,
            true
        );

        autoPtr<autoOptimize> optimMeshPtr
            = autoOptimize::New(mesh, meshOptimDict);
        optimMeshPtr->optimize();

        ppLayersLocal = -1;
        ppPtLayersLocal = -1;
        forAll(pp, i)
        {
            label nLay = ppLayers[i];
            if (nLay > 2)
            {
                ppLayersLocal[i] = nLay - 1;
            }
            else
            {
                ppLayersLocal[i] = 1;
            }
        }

        forAll(pp.meshPoints(), i)
        {
            label nLay = maxPointLayers[i];
            if (nLay > 2)
            {
                ppPtLayersLocal[i] = nLay - 1;
            }
            else
            {
                ppPtLayersLocal[i] = 1;
            }
        }

        //Add final layers
        updateExtrudeMesh
        (
            layerPatchIDs,
            ppLayersLocal,
            ppPtLayersLocal,
            pp,
            ppLayers,
            maxPointLayers
        );
        layerManip.stretchLayers(controller_);
    }
    else
    {
        labelList ppLayersLocal(pp.size(), 0);
        labelList ppPtLayersLocal(pp.meshPoints().size(), 0);

        ppLayersLocal = ppLayers;
        ppPtLayersLocal = maxPointLayers;
        //Add layers
        updateExtrudeMesh
        (
            layerPatchIDs,
            ppLayersLocal,
            ppPtLayersLocal,
            pp,
            ppLayers,
            maxPointLayers
        );

        layerManipulate layerManip
        (
            mesh,
            layerParams_,
            meshRefiner_.meshCutter().cellLevel(),
            meshRefiner_.meshCutter().pointLevel(),
            meshRefiner_.meshCutter().level0EdgeLength()
         );

        labelList layerOffset(pp.meshPoints().size(), 0);
        labelList adjustedNLayers(pp.meshPoints().size(), 0);
        layerManip.smoothLayerStack
        (
            controller_,
            nSmooth,
            grownUpGeometryDict_,
            grownUpZoneGeometryDict_,
            layerOffset,
            adjustedNLayers
        );
    }

    //convert baffle zones to patches
    convertBaffleZones();

    //Check if want to fan layers around convex edges
    autoFanExtrude fExtrude
    (
        meshRefiner_,
        decomposer_,
        distributor_,
        layerParams_
    );
    fExtrude.setRefinement();

    label maxMergeIter = layerParams_.maxMergePreIter();
    if (maxMergeIter > 0)
    {
        autoLayerCellsMerge autoMerge
        (
            meshRefiner_,
            decomposer_,
            distributor_,
            layerParams_
        );
        autoMerge.merge(maxMergeIter);
    }

    return;
}


void Foam::addMultiLayers::setRefinementDual()
{
    Info<< "Adding final layer cells: "<<endl;

    fvMesh& mesh = meshRefiner_.mesh();

    const labelList& numLayers = layerParams_.numLayers();
    const label& nSmoothPerLayer = layerParams_.nSmoothPerDualLayer();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    volScalarField targetLayers
    (
        IOobject
        (
            "targetLayers",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    volScalarField layerCount
    (
        IOobject
        (
            "layerCount",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    volScalarField actualLayers
    (
        IOobject
        (
            "actualLayers",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0)
    );

    DynamicList<label> layerPatchIDs(patches.size());
    label maxLayers = -1;
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            if (numLayers[patchI] > 1)
            {
                scalarField& targetPatchLayers =
                    targetLayers.boundaryFieldRef()[patchI];
                forAll(targetPatchLayers, i)
                {
                    targetPatchLayers[i] = numLayers[patchI]-2;
                }
                maxLayers = max(maxLayers, numLayers[patchI]-1);
            }
            else if (numLayers[patchI] == 1)
            {
                maxLayers = max(maxLayers, numLayers[patchI]-1);
            }

            if (numLayers[patchI] >= 0)
            {
                layerPatchIDs.append(patchI);
            }
        }
    }
    layerPatchIDs.shrink();

    if (maxLayers <= 0)
    {
        if (maxLayers == 0)
        {
            layerManipulate layerManip
            (
                mesh,
                layerParams_,
                meshRefiner_.meshCutter().cellLevel(),
                meshRefiner_.meshCutter().pointLevel(),
                meshRefiner_.meshCutter().level0EdgeLength()
            );

            autoPtr<indirectPrimitivePatch> ppPtr
            (
                meshRefinement::makePatch
                (
                    mesh,
                    layerPatchIDs
                 )
            );
            indirectPrimitivePatch& pp = ppPtr();

            labelList layerOffset(pp.meshPoints().size(), 0);
            labelList adjustedNLayers(pp.meshPoints().size(), 0);
            layerManip.smoothLayerStack
            (
                controller_,
                nSmoothPerLayer,
                grownUpGeometryDict_,
                grownUpZoneGeometryDict_,
                layerOffset,
                adjustedNLayers
            );
        }
        remergeOuterSplitCells();
        return;
    }

    label nLayerIter = 0;
    while (true)
    {
        nLayerIter++;

        labelList layerFaces = calculateLayers
        (
            layerPatchIDs,
            targetLayers,
            layerCount,
            actualLayers
        );

        if (returnReduce(layerFaces.size(), sumOp<label>()) == 0)
        {
            label nSmooth = nSmoothPerLayer * nLayerIter;
            layerManipulate layerManip
            (
                mesh,
                layerParams_,
                meshRefiner_.meshCutter().cellLevel(),
                meshRefiner_.meshCutter().pointLevel(),
                meshRefiner_.meshCutter().level0EdgeLength()
            );

            autoPtr<indirectPrimitivePatch> ppPtr
            (
                meshRefinement::makePatch
                (
                    mesh,
                    layerPatchIDs
                 )
            );
            indirectPrimitivePatch& pp = ppPtr();

            labelList layerOffset(pp.meshPoints().size(), 0);
            labelList adjustedNLayers(pp.meshPoints().size(), 0);
            layerManip.smoothLayerStack
            (
                controller_,
                nSmooth,
                grownUpGeometryDict_,
                grownUpZoneGeometryDict_,
                layerOffset,
                adjustedNLayers
            );

            //convert baffle zones to patches
            convertBaffleZones();

            break;
        }

        autoPtr<indirectPrimitivePatch> ppGrownPtr
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>(mesh.faces(), layerFaces),
                mesh.points()
            )
        );
        indirectPrimitivePatch& ppGrown = ppGrownPtr();

        addHexMeshLayer addRefLayers
        (
            globalToMasterPatch_,
            refineParams_,
            decomposer_,
            distributor_,
            meshRefiner_,
            ppGrown,
            0.8,//cutRatio,
            false,//geometric checks on cuts
            layerPatchIDs,
            true
        );
        addRefLayers.setRefinement();
        reorientateFaces();

        Info<< "Layer Iteration: "<< nLayerIter
            << " Number of cells "
            << returnReduce(mesh.nCells(), sumOp<label>())
            << endl;
    }

    //remergeOuterSplitCells();

    //Move points for final layer wrong oriented prism cells
    //handleWrongOrientedFaces();

    return;
}

// ************************************************************************* //
