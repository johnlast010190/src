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
    (c) 2020 Esi Ltd.

InClass
    autoLayerCellsMerge

\*---------------------------------------------------------------------------*/

#include "autoLayerCellsMerge/autoLayerCellsMerge.H"

#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveCell.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"

#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(autoLayerCellsMerge, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::autoLayerCellsMerge::markAndRedistribute
(
    boolList& markedCells
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbours = mesh.faceNeighbour();

    cellSet errorCells(mesh, "errorCells", mesh.nCells()/100+1);
    faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+1);

    label noFailedChecks(0);
    if (mesh.checkCellVolumes(true, &errorCells)) noFailedChecks++;
    if (mesh.checkFacePyramids(true, -SMALL, &errorFaces)) noFailedChecks++;

    if (noFailedChecks == 0)
    {
        return false;
    }

    const volScalarField& layerCells =
        mesh.lookupObject<volScalarField>("layerStacks");

    markedCells.setSize(mesh.nCells());
    markedCells = false;

    forAllConstIter(labelHashSet, errorFaces, iter)
    {
        label facei = iter.key();
        label own = owners[facei];
        if (layerCells[own] > -1)
        {
            markedCells[own] = true;
        }
        if (mesh.isInternalFace(facei))
        {
            label nei = neighbours[facei];
            if (layerCells[nei] > -1)
            {
                markedCells[nei] = true;
            }
        }
    }

    forAllConstIter(labelHashSet, errorCells, iter)
    {
        label celli = iter.key();
        if (layerCells[celli] > -1)
        {
            markedCells[celli] = true;
        }
    }

    labelList neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());
    labelList neiMarkedCells(mesh.nFaces()-mesh.nInternalFaces());
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        neiLayerCells[facei-mesh.nInternalFaces()] =
            layerCells[owners[facei]];
        neiMarkedCells[facei-mesh.nInternalFaces()] =
            markedCells[owners[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLayerCells);
    syncTools::swapBoundaryFaceList(mesh, neiMarkedCells);

    forAll(mesh.cells(), celli)
    {
        if (!markedCells[celli] && layerCells[celli] > -1)
        {
            const cell& c = mesh.cells()[celli];
            label nMarkedNbr = 0;
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchi = patches.whichPatch(facei);
                if (patchi == -1 || patches[patchi].coupled())
                {
                    label neiLayerID = -1;
                    bool neiMarked = false;
                    if (patchi == -1)
                    {
                        if (owners[facei] == celli)
                        {
                            label nei = neighbours[facei];
                            neiLayerID = layerCells[nei];
                            neiMarked = markedCells[nei];
                        }
                        else
                        {
                            label nei = owners[facei];
                            neiLayerID = layerCells[nei];
                            neiMarked = markedCells[nei];
                        }
                    }
                    else if (patches[patchi].coupled())
                    {
                        neiLayerID =
                            neiLayerCells[facei-mesh.nInternalFaces()];
                        neiMarked =
                            neiMarkedCells[facei-mesh.nInternalFaces()];
                    }

                    if
                    (
                        neiMarked
                        && (layerCells[celli] == neiLayerID)
                     )
                    {
                        nMarkedNbr++;
                    }
                }
            }

            if (nMarkedNbr == 2)
            {
                markedCells[celli] = true;
            }
        }
    }

    boolList marked(mesh.nFaces(), true);
    forAll(markedCells, celli)
    {
        if (markedCells[celli])
        {
            const cell& c = mesh.cells()[celli];
            forAll(c, cFI)
            {
                marked[c[cFI]] = false;
            }
        }
    }
    syncTools::syncFaceList
    (
        mesh,
        marked,
        andEqOp<bool>()
     );

    neiMarkedCells = -1;
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        neiMarkedCells[facei-mesh.nInternalFaces()] =
            markedCells[owners[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiMarkedCells);
    forAll(mesh.cells(), celli)
    {
        if (!markedCells[celli] && layerCells[celli] > -1)
        {
            const cell& c = mesh.cells()[celli];
            bool found = false;
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchi = patches.whichPatch(facei);
                label neiLayerID = -1;
                bool neiMarked = false;
                if (patchi == -1)
                {
                    if (owners[facei] == celli)
                    {
                        label nei = neighbours[facei];
                        neiLayerID = layerCells[nei];
                        neiMarked = markedCells[nei];
                    }
                    else
                    {
                        label nei = owners[facei];
                        neiLayerID = layerCells[nei];
                        neiMarked = markedCells[nei];
                    }
                }
                else if (patches[patchi].coupled())
                {
                    neiLayerID =
                        neiLayerCells[facei-mesh.nInternalFaces()];
                    neiMarked =
                        neiMarkedCells[facei-mesh.nInternalFaces()];
                }

                if
                (
                    neiMarked
                    && (layerCells[celli] == neiLayerID)
                )
                {
                    found = true;
                    break;
                }
            }

            if (found)
            {
                const cell& c = mesh.cells()[celli];
                forAll(c, cFI)
                {
                    marked[c[cFI]] = false;
                }
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh,
        marked,
        andEqOp<bool>()
     );

    autoPtr<mapDistributePolyMesh> distMap = meshRefiner_.balance
    (
        marked,
        scalarField(mesh.nCells(), 1), // dummy weights
        decomposer_,
        distributor_,
        false
     );

    if (distMap.valid())
    {
        distMap().distributeCellData(markedCells);
    }

    return true;
}


label Foam::autoLayerCellsMerge::updateMesh
(
    const boolList& markedCells
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbours = mesh.faceNeighbour();

    const labelList& patchToNLayers = layerParams_.numLayers();

    const volScalarField& layerCells =
        mesh.lookupObject<volScalarField>("layerStacks");

    label nMerged = 0;

    labelList neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
     )
    {
        neiLayerCells[facei-mesh.nInternalFaces()] =
            layerCells[owners[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLayerCells);

    boolList topFaces(mesh.nFaces(), false);
    forAll(mesh.faces(), facei)
    {
        label own = owners[facei];
        scalar ownLayerID = layerCells[own];

        label patchi = patches.whichPatch(facei);
        scalar neiLayerID = -1;
        if (patchi == -1)
        {
            label nei = neighbours[facei];
            neiLayerID = layerCells[nei];
        }
        else if (patches[patchi].coupled())
        {
            neiLayerID =
                neiLayerCells[facei-mesh.nInternalFaces()];
        }
        else
        {
            continue;
        }

        if
        (
            (ownLayerID > -1 && neiLayerID < 0)
            || (ownLayerID < 0 && neiLayerID > -1)
            || (ownLayerID == neiLayerID && ownLayerID > -1)
         )
        {
            topFaces[facei] = true;
        }
    }

    DynamicList<label> allPatches(patches.size());
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        if (!pp.coupled() && patchToNLayers[patchi] > -1)
        {
            allPatches.append(patchi);
        }
    }

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            allPatches
        )
    );
    indirectPrimitivePatch& pp = ppPtr();

    boolList seed(mesh.nFaces(), false);
    forAll(pp, i)
    {
        label facei = pp.addressing()[i];
        topFaces[facei] = true;
        seed[facei] = true;
    }

    boolList visitedCells(mesh.nCells(), false);

    DynamicList<Tuple2<labelList,labelList>> mergeCells(pp.size());
    boolList newSeedFaces(mesh.nFaces(), false);
    DynamicList<label> mCells(50);
    DynamicList<label> mFaces(50);

    while (true)
    {
        newSeedFaces = false;

        label nTransfers = 0;
        forAll(seed, facei)
        {
            if (seed[facei])
            {
                label currFace = facei;

                while (currFace != -1)
                {
                    label own = mesh.faceOwner()[currFace];
                    label currCell = -1;
                    if (mesh.isInternalFace(currFace))
                    {
                        if (!visitedCells[own])
                        {
                            currCell = own;
                        }
                        else
                        {
                            label nei = neighbours[currFace];
                            if (!visitedCells[nei])
                            {
                                currCell = nei;
                            }
                        }
                    }
                    else
                    {
                        if (!visitedCells[own])
                        {
                            currCell = own;
                        }
                    }

                    if (currCell != -1 && layerCells[currCell] > -1)
                    {
                        visitedCells[currCell] = true;
                        if (markedCells[currCell])
                        {
                            mFaces.append(currFace);
                            mCells.append(currCell);
                        }
                        else if (mFaces.size() > 0)
                        {
                            mFaces.append(currFace);
                            mergeCells.append
                            (
                                Tuple2<labelList,labelList>
                                (
                                    mFaces,
                                    mCells
                                 )
                            );
                            mFaces.clear();
                            mCells.clear();
                        }

                        const cell& c = mesh.cells()[currCell];
                        label oppFaceI = -1;

                        forAll(c, cFI)
                        {
                            label otherFace = c[cFI];
                            if
                            (
                                otherFace != currFace
                                && topFaces[otherFace]
                            )
                            {
                                oppFaceI = otherFace;
                            }
                        }

                        if (oppFaceI != -1)
                        {
                            label oppPatch =
                                patches.whichPatch(oppFaceI);
                            if
                            (
                                oppPatch != -1
                                && patches[oppPatch].coupled()
                            )
                            {
                                currFace = -1;
                                newSeedFaces[oppFaceI] = true;
                                nTransfers++;
                            }
                            else
                            {
                                currFace = oppFaceI;
                            }
                        }
                    }
                    else
                    {
                        if (mFaces.size() > 0)
                        {
                            mFaces.append(currFace);
                            mergeCells.append
                            (
                                Tuple2<labelList,labelList>
                                (
                                    mFaces,
                                    mCells
                                 )
                            );
                            mFaces.clear();
                            mCells.clear();
                        }
                        break;
                    }
                }
            }
        }

        if (returnReduce(nTransfers, sumOp<label>()) == 0)
        {
            break;
        }

        syncTools::syncFaceList
        (
            mesh,
            newSeedFaces,
            orEqOp<bool>()
        );
        seed = newSeedFaces;
    }

    forAll(mergeCells, i)
    {
        Tuple2<labelList,labelList>& stack =  mergeCells[i];

        if (stack.second().size() == 1 && stack.first().size() == 2)
        {
            label startFace = stack.first()[0];
            label endFace = stack.first()[1];
            label startCell = stack.second()[0];

            label startPatch = patches.whichPatch(startFace);
            label endPatch = patches.whichPatch(endFace);

            label addedCell = -1;
            label addedFace = -1;

            label matchFace  = -1;
            if (startPatch == -1)
            {
                matchFace = startFace;
            }
            else if (endPatch == -1)
            {
                matchFace = endFace;
            }

            if (matchFace != -1)
            {
                label own = owners[matchFace];
                if (own == startCell)
                {
                    addedCell = neighbours[matchFace];
                }
                else
                {
                    addedCell = own;
                }

                if (layerCells[addedCell] < 0)
                {
                    addedCell = -1;
                }
            }

            if (addedCell != -1)
            {
                const cell& aCell = mesh.cells()[addedCell];
                forAll(aCell, aCI)
                {
                    label facei = aCell[aCI];
                    if (facei != matchFace && topFaces[facei])
                    {
                        addedFace = facei;
                        break;
                    }
                }
            }

            if (addedCell != -1 && addedFace != -1)
            {
                labelList sFaces(3);
                labelList sCells(2);
                if (startPatch == -1)
                {
                    sFaces[0] = addedFace;
                    sFaces[1] = stack.first()[0];
                    sFaces[2] = stack.first()[1];
                    sCells[0] = addedCell;
                    sCells[1] = startCell;
                }
                else
                {
                    sFaces[0] = stack.first()[0];
                    sFaces[1] = stack.first()[1];
                    sFaces[2] = addedFace;
                    sCells[0] = startCell;
                    sCells[1] = addedCell;
                }
                stack.first() = sFaces;
                stack.second() = sCells;
            }
        }
    }

    polyTopoChange meshMod(mesh);

    //layer edge type non-layer (-1), top (0), side (1), removed (2)
    labelList layerEdgeType(mesh.nEdges(), -1);

    forAll(mesh.faces(), facei)
    {
        if (topFaces[facei])
        {
            const labelList& fEdges = mesh.faceEdges()[facei];
            forAll(fEdges, fEI)
            {
                layerEdgeType[fEdges[fEI]] = 0;
            }
        }
    }

    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            const labelList& cEdges = mesh.cellEdges()[celli];
            forAll(cEdges, cEI)
            {
                label edgei = cEdges[cEI];
                if (layerEdgeType[edgei] == -1)
                {
                    layerEdgeType[cEdges[cEI]] = 1;
                }
            }
        }
    }

    boolList removedFaces(mesh.nFaces(), false);
    forAll(mergeCells, i)
    {
        Tuple2<labelList,labelList>& stack =  mergeCells[i];
        if (stack.second().size() > 1)
        {
            const labelList& sFaces = stack.first();
            forAll(sFaces, sFI)
            {
                label facei = sFaces[sFI];
                if (sFI > 0 && sFI < sFaces.size()-1)
                {
                    removedFaces[facei] = true;
                }
            }
        }
    }

    forAll(mesh.edges(), edgei)
    {
        if (layerEdgeType[edgei] == 0)
        {
            const labelList& eFaces = mesh.edgeFaces()[edgei];
            label nEdgeFaces = 0;
            label nRemovedFaces = 0;
            forAll(eFaces, eFI)
            {
                label facei = eFaces[eFI];
                if (removedFaces[facei])
                {
                    nRemovedFaces++;
                }
                else
                {
                    nEdgeFaces++;
                }
            }
            if (nRemovedFaces > 0 && nEdgeFaces == 2)
            {
                layerEdgeType[edgei] = 2;
            }
        }
    }

    //Remove cells and connecting faces
    labelList newCellID(mesh.nCells(), -1);
    forAll(mergeCells, i)
    {
        Tuple2<labelList,labelList>& stack =  mergeCells[i];
        if (stack.second().size() > 1)
        {
            const labelList& sFaces = stack.first();
            const labelList& sCells = stack.second();

            label keptCell = sCells[0];
            forAll(sCells, sCI)
            {
                label celli = sCells[sCI];
                if (sCI > 0)
                {
                    newCellID[celli] = keptCell;
                    meshMod.setAction(polyRemoveCell(celli));
                }
            }
            nMerged += sCells.size();

            forAll(sFaces, sFI)
            {
                label facei = sFaces[sFI];
                if (sFI > 0 && sFI < sFaces.size()-1)
                {
                    meshMod.setAction(polyRemoveFace(facei));
                }
            }
        }
    }

    //Add new side faces
    boolList updatedFaces(mesh.nFaces(), false);
    forAll(mergeCells, i)
    {
        Tuple2<labelList,labelList>& stack =  mergeCells[i];
        if (stack.second().size() > 1)
        {
            const labelList& sCells = stack.second();

            DynamicList<label> mergeSideFaces(sCells.size()*6);
            forAll(sCells, sCI)
            {
                label celli = sCells[sCI];
                const cell& c = mesh.cells()[celli];
                forAll(c, cFI)
                {
                    label facei = c[cFI];
                    if (!topFaces[facei])
                    {
                        mergeSideFaces.append(facei);
                    }
                }
            }

            labelHashSet mergeSideFacesSet(mergeSideFaces);

            const labelList& sFaces = stack.first();
            label startFace = sFaces[0];
            label keptCell = sCells[0];
            bool insidePointing = true;
            if (owners[startFace] == sCells[0])
            {
                insidePointing = false;
            }

            const face& bf = mesh.faces()[startFace];
            forAll(bf, bFI)
            {
                label nextFp = bf.fcIndex(bFI);
                label startPt = bf[bFI];
                label endPt = bf[nextFp];

                label edgei = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[startPt],
                    startPt,
                    endPt
                 );

                edge e(startPt,endPt);

                DynamicList<labelList> layerUp(sCells.size());
                DynamicList<labelList> layerDown(sCells.size());
                DynamicList<label> layerTopEdges(sCells.size());
                DynamicList<label> layerSideFaces(sCells.size());
                label prevFace = -1;
                while (true)
                {
                    const labelList& eFaces = mesh.edgeFaces()[edgei];

                    label sFaceI = -1;
                    forAll(eFaces, eFI)
                    {
                        label facei = eFaces[eFI];
                        if
                        (
                            facei != prevFace
                            && mergeSideFacesSet.found(facei)
                         )
                        {
                            sFaceI = facei;
                            prevFace = facei;
                            break;
                        }
                    }

                    if (sFaceI != -1)
                    {
                        face f = mesh.faces()[sFaceI];
                        label start = findIndex(f, e[0]);
                        label endPt = e[1];
                        if (f[f.fcIndex(start)] == e[1])
                        {
                            f.flip();
                        }
                        start = findIndex(f, e[0]);

                        DynamicList<label> upPts(2);
                        DynamicList<label> downPts(2);

                        forAll(f, fp)
                        {
                            label meshPointI = f[start];
                            upPts.append(meshPointI);

                            start = f.fcIndex(start);
                            label nextPt = f[start];
                            label meshEdgeI = meshTools::findEdge
                            (
                                mesh.edges(),
                                mesh.pointEdges()[meshPointI],
                                meshPointI,
                                nextPt
                             );

                            if (layerEdgeType[meshEdgeI] != 1)
                            {
                                edgei = meshEdgeI;
                                e[0] = meshPointI;
                                e[1] = nextPt;
                                break;
                            }
                        }

                        start = findIndex(f, endPt);
                        forAll(f, fp)
                        {
                            label meshPointI = f[start];
                            downPts.append(meshPointI);
                            if (meshPointI == e[1])
                            {
                                break;
                            }
                            start = f.rcIndex(start);
                        }
                        layerUp.append(upPts);
                        layerDown.append(downPts);
                        layerTopEdges.append(edgei);
                        layerSideFaces.append(sFaceI);
                    }
                    else
                    {
                        break;
                    }
                }
                //add faces
                DynamicList<label> fup;
                DynamicList<label> fdown;
                DynamicList<label> sideFaces;

                forAll(layerUp, layeri)
                {
                    const labelList& lUpPts = layerUp[layeri];
                    const labelList& lDownPts = layerDown[layeri];
                    label ledgei = layerTopEdges[layeri];

                    label startIndex(fup.size() > 0 ? 1 : 0);

                    for (label i = startIndex; i < lUpPts.size(); i++)
                    {
                        fup.append(lUpPts[i]);
                    }
                    for (label i = startIndex; i < lDownPts.size(); i++)
                    {
                        fdown.append(lDownPts[i]);
                    }
                    sideFaces.append(layerSideFaces[layeri]);

                    if (layerEdgeType[ledgei] != 2)
                    {
                        label nUpdated = 0;
                        label masterFace = -1;
                        forAll(sideFaces, sFI)
                        {
                            label facei = sideFaces[sFI];
                            if (updatedFaces[facei])
                            {
                                nUpdated++;
                            }
                            else
                            {
                                meshMod.setAction
                                (
                                    polyRemoveFace(facei)
                                );
                                updatedFaces[facei] = true;
                                masterFace = facei;
                            }
                        }

                        if (nUpdated == 0)
                        {
                            DynamicList<label>
                                newFace(fup.size()+fdown.size());
                            forAll(fup,fp)
                            {
                                newFace.append(fup[fp]);
                            }

                            forAllReverse(fdown,fp)
                            {
                                newFace.append(fdown[fp]);
                            }
                            face f(newFace);

                            label own = owners[masterFace];

                            label newOwn
                            (
                                newCellID[own] != -1 ?
                                newCellID[own] : own
                            );

                            label patchi =
                                patches.whichPatch(masterFace);

                            bool flip = false;
                            label nei = -1;
                            if (patchi == -1)
                            {
                                nei = neighbours[masterFace];
                                label newNei
                                (
                                    newCellID[nei] != -1 ?
                                    newCellID[nei] : nei
                                );
                                if (newOwn == keptCell)
                                {
                                    if (newOwn < newNei)
                                    {
                                        own = newOwn;
                                        nei = newNei;
                                        if (insidePointing)
                                        {
                                            flip = true;
                                        }
                                    }
                                    else
                                    {
                                        own = newNei;
                                        nei = newOwn;
                                        if (!insidePointing)
                                        {
                                            flip = true;
                                        }
                                    }
                                }
                                else
                                {
                                    if (newOwn < newNei)
                                    {
                                        own = newOwn;
                                        nei = newNei;
                                        if (!insidePointing)
                                        {
                                            flip = true;
                                        }
                                    }
                                    else
                                    {
                                        own = newNei;
                                        nei = newOwn;
                                        if (insidePointing)
                                        {
                                            flip = true;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                own = newOwn;
                                if (insidePointing)
                                {
                                    flip = true;
                                }
                            }

                            if (flip)
                            {
                                f.flip();
                            }

                            label zoneI =
                                mesh.faceZones().whichZone(masterFace);
                            bool zoneFlip = false;
                            if (zoneI >= 0)
                            {
                                const faceZone& fZone =
                                    mesh.faceZones()[zoneI];
                                zoneFlip =
                                    fZone.flipMap()[fZone.whichFace(masterFace)];
                            }

                            meshMod.setAction
                            (
                                polyAddFace
                                (
                                    f,   // face
                                    own, // owner
                                    nei,  // neighbour
                                    -1,        // master point
                                    -1,        // master edge
                                    masterFace,     // master face
                                    false,     // flux flip
                                    patchi,        // patch for face
                                    zoneI,     // zone for face
                                    zoneFlip       // face zone flip
                                 )
                             );
                        }
                        fup.clear();
                        fdown.clear();
                        sideFaces.clear();
                    }
                }
            }
        }
    }

    //Modify top faces
    forAll(mergeCells, i)
    {
        Tuple2<labelList,labelList>& stack =  mergeCells[i];
        if (stack.second().size() > 1)
        {
            const labelList& sFaces = stack.first();
            labelList topFaces(2);
            topFaces[0] = sFaces[0];
            topFaces[1] = sFaces[sFaces.size()-1];

            forAll(topFaces, tFI)
            {
                label facei = topFaces[tFI];
                label own = owners[facei];
                label patchi = patches.whichPatch(facei);

                label zoneID = mesh.faceZones().whichZone(facei);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone =
                        mesh.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                }

                if (patchi == -1)
                {
                    label nei = neighbours[facei];
                    label newOwn
                    (
                        newCellID[own] != -1 ?
                        newCellID[own] : own
                     );
                    label newNei
                    (
                        newCellID[nei] != -1 ?
                        newCellID[nei] : nei
                     );

                    face f = mesh.faces()[facei];
                    if (newOwn > newNei)
                    {
                        f.flip();
                        label newNeiTmp = newNei;
                        newNei =  newOwn;
                        newOwn = newNeiTmp;
                    }

                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            f, // original face
                            facei,    // label of face
                            newOwn, // owner
                            newNei,             // neighbour
                            false,          // face flip
                            patchi,         // patch for face
                            false,          // remove from zone
                            zoneID,         // zone for face
                            zoneFlip        // face flip in zone
                         )
                     );
                }
                else
                {
                    label newOwn
                    (
                        newCellID[own] != -1 ?
                        newCellID[own] : own
                    );
                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            mesh.faces()[facei], // original face
                            facei,    // label of face
                            newOwn,         // owner
                            -1,             // neighbour
                            false,          // face flip
                            patchi,         // patch for face
                            false,          // remove from zone
                            zoneID,         // zone for face
                            zoneFlip        // face flip in zone
                         )
                     );
                }
            }
        }
    }

    reduce(nMerged, sumOp<label>());
    Info<<"Merged " << nMerged <<" layer cells"<<endl;

    if (nMerged > 0)
    {
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

    return nMerged;
}


void Foam::autoLayerCellsMerge::removePoints()
{
    fvMesh& mesh = meshRefiner_.mesh();

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));

    labelList nPointEdges(mesh.nPoints(), 0);

    forAll(mesh.points(), pointI)
    {
        const labelList& pEdges = mesh.pointEdges()[pointI];
        forAll(pEdges, pEI)
        {
            if (isMasterEdge[pEdges[pEI]])
            {
                nPointEdges[pointI]++;
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        nPointEdges,
        plusEqOp<label>(), // combine op
        label(0)     // null value
     );

    polyTopoChange meshMod(mesh);
    forAll(mesh.faces(), faceI)
    {
        const face f = mesh.faces()[faceI];
        DynamicList<label> keepPts(f.size());
        forAll(f, fp)
        {
            if (nPointEdges[f[fp]] > 2)
            {
                keepPts.append(f[fp]);
            }
        }

        if (keepPts.size() != f.size())
        {
            if (keepPts.size() < 3)
            {
                meshMod.setAction
                (
                    polyRemoveFace(faceI)
                );
            }
            else
            {
                label patchID = mesh.boundaryMesh().whichPatch(faceI);
                label nei =
                    (patchID == -1 ? mesh.faceNeighbour()[faceI] : -1);
                label zoneID = mesh.faceZones().whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        face(keepPts),           // modified face
                        faceI,                   // label of face
                        mesh.faceOwner()[faceI], // owner
                        nei,                     // neighbour
                        false,                    // face flip
                        patchID,                  // patch for face
                        false,                    // remove from zone
                        zoneID,                   // zone for face
                        zoneFlip                  // face flip in zone
                     )
                 );
            }
        }
    }

    label nPtsRemoved = 0;
    forAll(mesh.points(), pointI)
    {
        if (nPointEdges[pointI] <= 2)
        {
            nPtsRemoved++;
            meshMod.setAction
            (
                polyRemovePoint(pointI)
            );
        }
    }

    reduce(nPtsRemoved, sumOp<label>());

    if (nPtsRemoved > 0)
    {
        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

        // Update fields
        mesh.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            mesh.clearOut();
        }
        meshRefiner_.updateMesh(map, labelList(0));
    }

    return;
}


void Foam::autoLayerCellsMerge::merge
(
    const label nMergeIter
)
{
    label minMergedCells = labelMax;

    for (label i = 0; i < nMergeIter; i++)
    {
        boolList markedCells;
        if (markAndRedistribute(markedCells))
        {
            label nMerged = updateMesh(markedCells);
            if (nMerged > 0 && nMerged < minMergedCells)
            {
                minMergedCells = nMerged;
                removePoints();
            }
            else
            {
                break;
            }
        }
        else
        {
            break;
        }
    }
    return;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::autoLayerCellsMerge::autoLayerCellsMerge
(
    meshRefinement& meshRefiner,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const layerParameters& layerParams
)
:
    meshRefiner_(meshRefiner),
    decomposer_(decomposer),
    distributor_(distributor),
    layerParams_(layerParams)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::autoLayerCellsMerge::~autoLayerCellsMerge()
{}

// ************************************************************************* //
