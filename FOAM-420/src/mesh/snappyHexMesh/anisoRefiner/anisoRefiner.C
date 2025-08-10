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
    (c) 2014-1014 Esi Ltd.
\*---------------------------------------------------------------------------*/

#include "anisoRefiner/anisoRefiner.H"
#include "fvMesh/fvMesh.H"
#include "shellSurfaces/shellSurfaces.H"
#include "meshCut/cellLooper/geomCellLooper.H"
#include "meshCut/cellCuts/cellCuts.H"
#include "meshCut/meshModifiers/meshCutter/meshCutter.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "sets/topoSets/cellSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(anisoRefiner, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Updates refineCell (cells marked for refinement) so across all faces
// there will be 2:1 consistency after refinement.
Foam::label Foam::anisoRefiner::faceConsistentRefinement
(
    const direction cmpt,
    PackedBoolList& isProtectedCell,
    boolList& protectedFaces,
    PackedBoolList& refineCell
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    label nChanged = 0;

    for (int iter = 0; iter < 2; iter++)
    {
        // Internal faces.
        for (label faceI = 0; faceI < mesh.nInternalFaces(); faceI++)
        {
            if (iter == 0 && !protectedFaces[faceI])
            {
                continue;
            }

            label own = mesh.faceOwner()[faceI];
            label ownLevel = cellLevel_[own][cmpt] + refineCell.get(own);

            label nei = mesh.faceNeighbour()[faceI];
            label neiLevel = cellLevel_[nei][cmpt] + refineCell.get(nei);

            if (ownLevel > (neiLevel+1))
            {
                if (isProtectedCell[nei])
                {
                    isProtectedCell.set(own);
                    const cell c = mesh.cells()[own];
                    forAll(c, cFI)
                    {
                        protectedFaces[c[cFI]] = true;
                    }
                    refineCell.unset(own);
                }
                else
                {
                    refineCell.set(nei);
                }
                nChanged++;
            }
            else if (neiLevel > (ownLevel+1))
            {
                if (isProtectedCell[own])
                {
                    isProtectedCell.set(nei);
                    const cell c = mesh.cells()[nei];
                    forAll(c, cFI)
                    {
                        protectedFaces[c[cFI]] = true;
                    }
                    refineCell.unset(nei);
                }
                else
                {
                    refineCell.set(own);
                }
                nChanged++;
            }
        }

        // Coupled faces. Swap owner level to get neighbouring cell level.
        // (only boundary faces of neiLevel used)
        labelList neiLevel(mesh.nFaces()-mesh.nInternalFaces());
        boolList neiProtected(mesh.nFaces()-mesh.nInternalFaces(),false);

        forAll(neiLevel, i)
        {
            label own = mesh.faceOwner()[i+mesh.nInternalFaces()];

            neiLevel[i] = cellLevel_[own][cmpt] + refineCell.get(own);
            if (isProtectedCell[own])
            {
                neiProtected[i] = true;
            }
        }

            // Swap to neighbour
        syncTools::swapBoundaryFaceList(mesh, neiLevel);
        syncTools::swapBoundaryFaceList(mesh, neiProtected);

        // Now we have neighbour value see which cells need refinement
        forAll(neiLevel, i)
        {
            label faceI = i+mesh.nInternalFaces();
            label patchI = mesh.boundaryMesh().whichPatch(faceI);
            if
            (
                (iter == 0 && !protectedFaces[faceI])
                || !isA<processorPolyPatch>(mesh.boundaryMesh()[patchI])
            )
            {
                continue;
            }

            label own = mesh.faceOwner()[faceI];
            label ownLevel = cellLevel_[own][cmpt] + refineCell.get(own);

            if (neiLevel[i] > (ownLevel+1))
            {
                if (!isProtectedCell[own])
                {
                    refineCell.set(own);
                    nChanged++;
                }
            }
            else if (ownLevel > (neiLevel[i]+1))
            {
                if (neiProtected[i])
                {
                    isProtectedCell.set(own);
                    const cell c = mesh.cells()[own];
                    forAll(c, cFI)
                    {
                        protectedFaces[c[cFI]] = true;
                    }
                    refineCell.unset(own);
                    nChanged++;
                }
            }
        }
    }

    return nChanged;
}


Foam::scalar Foam::anisoRefiner::edgeTol(direction cmpt) const
{
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    label maxLevel = labelMin;

    forAll(cellLevel_, cellI)
    {
        label level = cellLevel_[cellI][cmpt];

        maxLevel = max(level,maxLevel);
    }
    reduce(maxLevel, maxOp<label>());

    return 0.001*edge0Len / (1<<maxLevel);
}


void Foam::anisoRefiner::modifyProcessorBoundaries
(
    const PackedBoolList& refineCell,
    const direction cmpt
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const pointField& cellCentres = mesh.cellCentres();

    geomCellLooper cellCutter(mesh);

    //Cut direction
    vector dir = vector::zero;
    dir[cmpt] = 1;

    const globalMeshData& pd = mesh.globalData();
    const labelList& coupledEdges = pd.coupledPatchMeshEdges();
    const labelList& coupledPoints = pd.coupledPatch().meshPoints();
    const edgeVertex ev(mesh);

    boolList vertIsCut(mesh.nPoints(), false);
    boolList edgeIsCut(mesh.nEdges(), false);
    scalarField edgeWeight(mesh.nEdges(), -GREAT);

    // Cut information per cut cell
    DynamicList<label> cutCells(mesh.nCells()/10 + 10);
    DynamicList<labelList> cellLoops(mesh.nCells()/10 + 10);
    DynamicList<scalarField> cellEdgeWeights(mesh.nCells()/10 + 10);

    // At least some cells are cut.
    polyTopoChange meshMod(mesh);

    boolList splitProcessorFaces(mesh.nFaces(), false);
    labelList newEdgePts(mesh.nEdges(), -1);
    labelList newEdgeLevels(mesh.nEdges(), -1);
    labelList cutPoints(mesh.nPoints(), -1);

    PackedBoolList updatedCells(mesh.nCells());

    boolList procEdges(mesh.nEdges(), false);
    forAll(coupledEdges, cEI)
    {
        label edgeI = coupledEdges[cEI];
        procEdges[edgeI] = true;
    }

    boolList procPoints(mesh.nEdges(), false);
    forAll(coupledPoints, cPI)
    {
        label pointI = coupledPoints[cPI];
        procPoints[pointI] = true;
    }

    forAll(coupledEdges, cEI)
    {
        label edgeI = coupledEdges[cEI];
        const labelList& eCells = mesh.edgeCells()[edgeI];
        forAll(eCells, eCI)
        {
            label cellI = eCells[eCI];

            if (refineCell.get(cellI) == 1 && updatedCells.get(cellI) == 0)
            {
                updatedCells.set(cellI);
                // Define plane through cc
                plane cutPlane(cellCentres[cellI], dir);

                labelList loop;
                scalarField loopWeights;

                if
                (
                    cellCutter.cut
                    (
                        cutPlane,
                        cellI,
                        vertIsCut,
                        edgeIsCut,
                        edgeWeight,
                        loop,
                        loopWeights
                     )
                )
                {
                    // Did manage to cut cell. Copy into overall list.
                    cutCells.append(cellI);
                    cellLoops.append(loop);
                    cellEdgeWeights.append(loopWeights);

                    DynamicList<label> cEdges(loop.size());
                    forAll(loop, i)
                    {
                        label cut = loop[i];

                        if (ev.isEdge(cut))
                        {
                            label edgeI = ev.getEdge(cut);

                            if (procEdges[edgeI] && !edgeIsCut[edgeI])
                            {
                                edgeIsCut[edgeI] = true;
                                edge e = mesh.edges()[edgeI];
                                point eC = e.centre(mesh.points());

                                label masterPoint
                                (
                                    pointLevel_[e[0]][cmpt] > pointLevel_[e[1]][cmpt]
                                    ? e[0] : e[1]
                                );

                                newEdgePts[edgeI] = meshMod.setAction
                                (
                                    polyAddPoint
                                    (
                                        eC,         // point
                                        masterPoint,// master point
                                        -1,         // zone for point
                                        true        // supports a cell
                                    )
                                );
                                newEdgeLevels[edgeI] = cellLevel_[cellI][cmpt];
                            }
                        }
                        else
                        {
                            cutPoints[ev.getVertex(cut)] = cellLevel_[cellI][cmpt];
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        cutPoints,
        maxEqOp<label>(),
        label(-1)
    );

    syncTools::syncEdgeList
    (
        mesh,
        edgeIsCut,
        orEqOp<bool>(),
        false              // null value
    );

    syncTools::syncEdgeList
    (
        mesh,
        newEdgeLevels,
        maxEqOp<label>(),
        label(-1)              // null value
    );

    forAll(coupledEdges, cEI)
    {
        label edgeI = coupledEdges[cEI];
        if (newEdgePts[edgeI] == -1 && edgeIsCut[edgeI])
        {
            edge e = mesh.edges()[edgeI];
            point eC = e.centre(mesh.points());

            label masterPoint
            (
                pointLevel_[e[0]][cmpt] > pointLevel_[e[1]][cmpt]
                ? e[0] : e[1]
            );

            newEdgePts[edgeI] = meshMod.setAction
            (
                polyAddPoint
                (
                    eC,         // point
                    masterPoint,       // master point
                    -1,         // zone for point
                    true        // supports a cell
                 )
             );
        }
    }

    forAll(cutCells, i)
    {
        label cellI = cutCells[i];

        const labelList& loop = cellLoops[i];

        DynamicList<label> cEdges(loop.size());
        DynamicList<label> cPoints(loop.size());
        forAll(loop, i)
        {
            label cut = loop[i];

            if (ev.isEdge(cut))
            {
                label edgeI = ev.getEdge(cut);
                if (procEdges[edgeI])
                {
                    cEdges.append(edgeI);
                }
            }
            else
            {
                label pointI = ev.getVertex(cut);
                if (procPoints[pointI])
                {
                    cPoints.append(pointI);
                }
            }
        }
        cEdges.shrink();
        cPoints.shrink();

        labelHashSet cellCutEdges(cEdges);
        labelHashSet cellCutPoints(cPoints);

        const labelList& cFaces = mesh.cells()[cellI];
        forAll(cFaces, cFI)
        {
            label faceI = cFaces[cFI];
            label patchI = patches.whichPatch(faceI);
            if (patchI != -1 && isA<processorPolyPatch>(patches[patchI]))
            {
                const labelList& fEdges = mesh.faceEdges()[faceI];
                label nEdgeCuts = 0;
                forAll(fEdges, fEI)
                {
                    label edgeI = fEdges[fEI];
                    if (cellCutEdges.found(edgeI))
                    {
                        nEdgeCuts++;
                    }
                }
                face f = mesh.faces()[faceI];
                label nPointCuts = 0;
                DynamicList<label> pCuts(2);

                forAll(f, fp)
                {
                    if (cellCutPoints.found(f[fp]))
                    {
                        pCuts.append(f[fp]);
                        nPointCuts++;
                    }
                }

                //split face
                if (nEdgeCuts + nPointCuts == 2)
                {
                    if (nPointCuts == 2)
                    {
                        label edgeI = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[pCuts[0]],
                            pCuts[0],
                            pCuts[1]
                        );
                        if (edgeI != -1)
                        {
                            continue;
                        }
                    }

                    splitFace
                    (
                        faceI,
                        cellI,
                        patchI,
                        cmpt,
                        cutPoints,
                        newEdgePts,
                        newEdgeLevels,
                        ev,
                        cellLevel_[cellI][cmpt],
                        meshMod
                    );
                    splitProcessorFaces[faceI] = true;
                }
            }
        }
    }

    boolList unsyncSplitFaces = splitProcessorFaces;
    syncTools::syncFaceList
    (
        mesh,
        splitProcessorFaces,
        orEqOp<bool>()
    );

    PackedBoolList updatedFaces(mesh.nFaces());
    forAll(coupledEdges, cEI)
    {
        label edgeI = coupledEdges[cEI];

        if (newEdgePts[edgeI] != -1)
        {
            const labelList& eFaces = mesh.edgeFaces()[edgeI];
            forAll(eFaces, eFI)
            {
                label faceI = eFaces[eFI];
                if (!splitProcessorFaces[faceI] && updatedFaces.get(faceI) == 0)
                {
                    updatedFaces.set(faceI);
                    face f = mesh.faces()[faceI];
                    label patchI = patches.whichPatch(faceI);

                    DynamicList<label> f0(f.size()+2);
                    forAll(f, fp)
                    {
                        label nextPt = f.fcIndex(fp);
                        f0.append(f[fp]);

                        label edgeI = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[f[fp]],
                            f[fp],
                            f[nextPt]
                         );

                        if (newEdgePts[edgeI] != -1)
                        {
                            f0.append(newEdgePts[edgeI]);
                        }
                    }

                    label zoneI = mesh.faceZones().whichZone(faceI);
                    bool flip = false;
                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh.faceZones()[zoneI];
                        flip = fz.flipMap()[fz.whichFace(faceI)];
                    }

                    label own = mesh.faceOwner()[faceI];
                    label nei = -1;
                    if (patchI == -1)
                    {
                        nei = mesh.faceNeighbour()[faceI];
                    }

                    face face0(f0.shrink());

                    //Modify face
                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            face0,      // modified face
                            faceI,      // label of face
                            own,      // owner
                            nei,         // neighbour
                            false,      // face flip
                            patchI, // patch for face
                            false,      // remove from zone
                            zoneI,      // zone for face
                            flip        // face flip in zone
                         )
                     );
                }
            }
        }
    }

    labelList nbrCellLevel(mesh.nFaces(), -1);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label meshFaceI = pp.start() + i;
                label own = mesh.faceOwner()[meshFaceI];
                nbrCellLevel[meshFaceI] = cellLevel_[own][cmpt];
            }
        }
    }
    syncTools::syncFaceList
    (
        mesh,
        nbrCellLevel,
        maxEqOp<label>()
    );

    forAll(mesh.faces(), faceI)
    {
        label own = mesh.faceOwner()[faceI];

        if (splitProcessorFaces[faceI] && !unsyncSplitFaces[faceI])
        {
            label patchI = patches.whichPatch(faceI);

            splitFace
            (
                faceI,
                own,
                patchI,
                cmpt,
                cutPoints,
                newEdgePts,
                newEdgeLevels,
                ev,
                nbrCellLevel[faceI],
                meshMod
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

    // Update intersection info
    meshRefiner_.updateMesh(map, labelList(0));

    List<labelVector> newPointLevel(mesh.nPoints());
    label nOldPoints = map().nOldPoints();

    forAll(map().reversePointMap(), pointI)
    {
        if (pointI >= nOldPoints)
        {
            label oldPointI = map().pointMap()[pointI];
            newPointLevel[pointI] = pointLevel_[oldPointI];
            newPointLevel[pointI][cmpt]++;
        }
        else
        {
            newPointLevel[pointI] = pointLevel_[pointI];
        }
    }

    pointLevel_.transfer(newPointLevel);

}


void Foam::anisoRefiner::refine
(
    const labelList& cellsToRefine,
    const direction cmpt
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const pointField& cellCentres = mesh.cellCentres();

    geomCellLooper cellCutter(mesh);

    vector dir = vector::zero;
    dir[cmpt] = 1;

    // Cut information per mesh entity
    boolList vertIsCut(mesh.nPoints(), false);
    boolList edgeIsCut(mesh.nEdges(), false);
    scalarField edgeWeight(mesh.nEdges(), -GREAT);

    // Cut information per cut cell
    DynamicList<label> cutCells(mesh.nCells()/10 + 10);
    DynamicList<labelList> cellLoops(mesh.nCells()/10 + 10);
    DynamicList<scalarField> cellEdgeWeights(mesh.nCells()/10 + 10);

    forAll(cellsToRefine, i)
    {
        label cellI = cellsToRefine[i];
        // Define plane through cc
        plane cutPlane(cellCentres[cellI], dir);

        labelList loop;
        scalarField loopWeights;

        if
        (
            cellCutter.cut
            (
                cutPlane,
                cellI,
                vertIsCut,
                edgeIsCut,
                edgeWeight,
                loop,
                loopWeights
             )
        )
        {
            // Did manage to cut cell. Copy into overall list.
            cutCells.append(cellI);
            cellLoops.append(loop);
            cellEdgeWeights.append(loopWeights);
        }
    }

    // At least some cells are cut.
    polyTopoChange meshMod(mesh);

    cellCuts cuts
    (
        mesh,
        cutCells,       // cells candidate for cutting
        cellLoops,
        cellEdgeWeights
    );

    // Cutting engine
    meshCutter cutter(mesh);

    // Insert mesh refinement into polyTopoChange.
    cutter.setRefinement(cuts, meshMod);

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

    // Update intersection info
    meshRefiner_.updateMesh(map, labelList(0));

    List<labelVector> newCellLevel(mesh.nCells());
    label nOldCells = map().nOldCells();


    forAll(map().reverseCellMap(), cellI)
    {
        if (cellI >= nOldCells)
        {
            label oldCellI = map().cellMap()[cellI];
            newCellLevel[cellI] = cellLevel_[oldCellI];
            //increment split cell level
            newCellLevel[oldCellI][cmpt]++;
            newCellLevel[cellI][cmpt]++;
        }
        else
        {
            newCellLevel[cellI] = cellLevel_[cellI];
        }
    }

    List<labelVector> newPointLevel(mesh.nPoints());
    label nOldPoints = map().nOldPoints();

    forAll(map().reversePointMap(), pointI)
    {
        if (pointI >= nOldPoints)
        {
            label oldPointI = map().pointMap()[pointI];
            newPointLevel[pointI] = pointLevel_[oldPointI];
            newPointLevel[pointI][cmpt]++;
        }
        else
        {
            newPointLevel[pointI] = pointLevel_[pointI];
        }
    }

    pointLevel_.transfer(newPointLevel);
    cellLevel_.transfer(newCellLevel);
}


void Foam::anisoRefiner::splitFace
(
    const label faceI,
    const label own,
    const label patchI,
    const direction cmpt,
    const labelList& cutPoints,
    const labelList& newEdgePts,
    const labelList& newEdgeLevels,
    const edgeVertex& ev,
    const label level,
    polyTopoChange& meshMod
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    face f = mesh.faces()[faceI];

    DynamicList<label> cuts(2);

    forAll(f, fp)
    {
        label nextPt = f.fcIndex(fp);

        label edgeNext = meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[f[fp]],
            f[fp],
            f[nextPt]
        );

        if (newEdgePts[edgeNext] != -1 && newEdgeLevels[edgeNext] == level)
        {
            cuts.append(ev.edgeToEVert(edgeNext));
        }
        else if (cutPoints[f[fp]] == level)
        {
            cuts.append(ev.vertToEVert(f[fp]));
        }
    }
    cuts.shrink();

    if (cuts.size() == 2)
    {
        DynamicList<label> f0(f.size()+2);
        DynamicList<label> f1(f.size()+2);
        label curPt = 0;
        label index = 0;
        bool start = false;
        //Calculate first cut face
        while (true)
        {
            label nextPt = f.fcIndex(curPt);

            label edgeI = meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[f[curPt]],
                f[curPt],
                f[nextPt]
            );

            if (ev.isEdge(cuts[index]) && (edgeI == ev.getEdge(cuts[index])))
            {
                index++;
                if (index == 1)
                {
                    f0.append(newEdgePts[edgeI]);
                    start = true;
                }
                else if (index == 2)
                {
                    f0.append(f[curPt]);
                    if (newEdgePts[edgeI] != -1)
                    {
                        f0.append(newEdgePts[edgeI]);
                    }
                    break;
                }
            }
            else if (!ev.isEdge(cuts[index]) && (f[curPt] == ev.getVertex(cuts[index])))
            {
                index++;
                f0.append(f[curPt]);
                if (index == 1)
                {
                    if (newEdgePts[edgeI] != -1)
                    {
                        f0.append(newEdgePts[edgeI]);
                    }

                    start = true;
                }
                else if (index == 2)
                {
                    break;
                }
            }
            else if (start)
            {
                f0.append(f[curPt]);
                if (newEdgePts[edgeI] != -1)
                {
                    f0.append(newEdgePts[edgeI]);
                }
            }
            curPt = nextPt;
        }

        curPt = 0;
        index = 0;
        start = false;
        label cutTemp = cuts[0];
        cuts[0] = cuts[1];
        cuts[1] = cutTemp;
        //Calculate second cut face
        while (true)
        {
            label nextPt = f.fcIndex(curPt);

            label edgeI = meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[f[curPt]],
                f[curPt],
                f[nextPt]
             );

            if (ev.isEdge(cuts[index]) && (edgeI == ev.getEdge(cuts[index])))
            {
                index++;
                if (index == 1)
                {
                    f1.append(newEdgePts[edgeI]);
                    start = true;
                }
                else if (index == 2)
                {
                    f1.append(f[curPt]);
                    if (newEdgePts[edgeI] != -1)
                    {
                        f1.append(newEdgePts[edgeI]);
                    }
                    break;
                }
            }
            else if (!ev.isEdge(cuts[index]) && (f[curPt] == ev.getVertex(cuts[index])))
            {
                index++;
                f1.append(f[curPt]);
                if (index == 1)
                {
                    if (newEdgePts[edgeI] != -1)
                    {
                        f1.append(newEdgePts[edgeI]);
                    }
                    start = true;
                }
                else if (index == 2)
                {
                    break;
                }
            }
            else if (start)
            {
                f1.append(f[curPt]);
                if (newEdgePts[edgeI] != -1)
                {
                    f1.append(newEdgePts[edgeI]);
                }
            }
            curPt = nextPt;
        }


        face face0(f0);
        face face1(f1);

        label zoneI = mesh.faceZones().whichZone(faceI);
        bool flip = false;
        if (zoneI != -1)
        {
            const faceZone& fz = mesh.faceZones()[zoneI];
            flip = fz.flipMap()[fz.whichFace(faceI)];
        }

        //Modify face
        meshMod.setAction
        (
            polyModifyFace
            (
                face0,      // modified face
                faceI,      // label of face
                own,      // owner
                -1,         // neighbour
                false,      // face flip
                patchI,     // patch for face
                false,      // remove from zone
                zoneI,      // zone for face
                flip        // face flip in zone
            )
        );

        //Addface
        meshMod.setAction
        (
            polyAddFace
            (
                face1,   // face
                own,     // owner
                -1,      // neighbour
                -1,      // master point
                -1,      // master edge
                faceI,   // master face
                false,   // flux flip
                patchI,  // patch for face
                zoneI,   // zone for face
                flip     // face zone flip
             )
         );
    }
}


// * * * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * //


void Foam::anisoRefiner::distribute(const mapDistributePolyMesh& map)
{
    // Update celllevel
    map.distributeCellData(cellLevel_);
    // Update pointlevel
    map.distributePointData(pointLevel_);
}


Foam::label Foam::anisoRefiner::setRefinement
(
    const direction cmpt
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const shellSurfaces& shells = meshRefiner_.shells();

    const vectorField& cellCentres = mesh.cellCentres();

    PackedBoolList isBoundaryCell(mesh.nCells());
    PackedBoolList isProtectedCell(mesh.nCells());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (!pp.coupled())
        {
            forAll(pp, i)
            {
                label meshFaceI = pp.start() + i;
                label own = mesh.faceOwner()[meshFaceI];
                isBoundaryCell.set(own);
                isProtectedCell.set(own);
            }
        }
    }

    // Collect cells to test
    pointField testCc(cellLevel_.size());
    labelList testLevels(cellLevel_.size());
    labelList testMap(cellLevel_.size());
    label testI = 0;

    forAll(mesh.cells(), cellI)
    {
        testCc[testI] = cellCentres[cellI];
        testLevels[testI] = cellLevel_[cellI][cmpt];
        testMap[testI] = cellI;
        testI++;
    }
    testCc.setSize(testI);
    testLevels.setSize(testI);
    testMap.setSize(testI);

    // Do test to see whether cells is inside/outside shell with higher level
    labelList maxLevel;
    labelList minCmptLevel;
    labelList shellCheck;
    shells.findHigherLevel
    (
        mesh,
        testCc,
        testLevels,
        testMap,
        maxLevel,
        minCmptLevel,
        shellCheck,
        shellSurfaces::ANISO,
        cmpt
    );

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.

    PackedBoolList refineCell(mesh.nCells());
    PackedBoolList unmarkedCell(mesh.nCells());
    testI = 0;
    label nBufferLayers = -1;
    forAll(testMap, i)
    {
        label cellI = testMap[i];
        label shellI = shellCheck[i];

        if (maxLevel[i] > testLevels[i])
        {
            if (isBoundaryCell.get(cellI) == 0)
            {
                refineCell.set(cellI);
            }
            else if
            (
                minCmptLevel[i] > testLevels[i]
                || (shellI != -1 && shells.refineSurface()[shellI])
            )
            {
                if (shellI != -1 && shells.refineSurface()[shellI])
                {
                    nBufferLayers = max(nBufferLayers, maxLevel[i]);
                    unmarkedCell.set(cellI);
                }

                if (isProtectedCell[cellI])
                {
                    isProtectedCell.unset(cellI);
                }
                refineCell.set(cellI);
            }
        }
        testI++;
    }

    //Allow cells switched to unprotected to also have neighbours unprotected
    nBufferLayers = returnReduce(nBufferLayers, maxOp<label>());

    boolList unmarkedFaces(mesh.nFaces(), false);
    for (int i = 0; i < nBufferLayers; i++)
    {
        forAll(unmarkedCell, cellI)
        {
            if (unmarkedCell.get(cellI) == 1)
            {
                const cell c = mesh.cells()[cellI];
                forAll(c, cFI)
                {
                    unmarkedFaces[c[cFI]] = true;
                }
            }
        }
        syncTools::syncFaceList
        (
            mesh,
            unmarkedFaces,
            orEqOp<bool>()
        );

        forAll(mesh.cells(), cellI)
        {
            if (isProtectedCell[cellI])
            {
                const cell c = mesh.cells()[cellI];
                forAll(c, cFI)
                {
                    if (unmarkedFaces[c[cFI]])
                    {
                        unmarkedCell.set(cellI);
                        isProtectedCell.unset(cellI);
                        break;
                    }
                }
            }
        }
    }

    boolList protectedFaces(mesh.nFaces(), false);
    forAll(mesh.cells(), cellI)
    {
        if (isProtectedCell[cellI])
        {
            const cell c = mesh.cells()[cellI];
            forAll(c, cFI)
            {
                protectedFaces[c[cFI]] = true;
            }
        }
    }
    syncTools::syncFaceList
    (
        mesh,
        protectedFaces,
        orEqOp<bool>()
    );

    while (true)
    {
        label nChanged = faceConsistentRefinement
        (
            cmpt,
            isProtectedCell,
            protectedFaces,
            refineCell
        );

        reduce(nChanged, sumOp<label>());

        if (nChanged == 0)
        {
            break;
        }
    }

    label nCellsToRefine(0);

    DynamicList<label> cellsToRefine(testI);
    forAll(refineCell, cellI)
    {
        if (refineCell.get(cellI) == 1)
        {
            cellsToRefine.append(cellI);
            nCellsToRefine++;
        }
    }
    cellsToRefine.shrink();

    if (returnReduce(nCellsToRefine, sumOp<label>()) != 0)
    {
        geomCellLooper::setSnapTol(edgeTol(cmpt));

        modifyProcessorBoundaries(refineCell,cmpt);

        refine(cellsToRefine, cmpt);

        //Updated mesh cutter point and cell level
        meshRefiner_.meshCutter().updateLevels
        (
            cmptMax(pointLevel_),
            cmptMax(cellLevel_)
        );
    }

    return nCellsToRefine;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::anisoRefiner::anisoRefiner
(
    meshRefinement& meshRefiner,
    const labelList cellLevel,
    const labelList pointLevel
)
:
    meshRefiner_(meshRefiner)
{
    //construct aniostropic point and cell level lists
    cellLevel_.setSize(cellLevel.size());
    forAll(cellLevel_, cellI)
    {
        cellLevel_[cellI] = labelVector(cellLevel[cellI]);
    }

    pointLevel_.setSize(pointLevel.size());
    forAll(pointLevel_, pointI)
    {
        pointLevel_[pointI] = labelVector(pointLevel[pointI]);
    }
}

// ************************************************************************* //
