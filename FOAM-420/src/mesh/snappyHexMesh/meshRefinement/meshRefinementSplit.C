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
    (c) 2010-2011, Esi Ltd

\*---------------------------------------------------------------------------*/

#include "meshRefinement/meshRefinement.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChanger/polyTopoChanger.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/cellSet.H"
#include "regionSplit/regionSplit.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "meshCut/cellCuts/cellCuts.H"
#include "meshCut/meshModifiers/meshCutter/meshCutter.H"
#include "meshTools/meshTools.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::meshRefinement::diagonalCut
(
    const cellFeatures& cellFeat,
    const Map<label>& pointMap,
    const boolList&  featureEdgePoints,
    const label start,
    const label end,
    DynamicList<label>& loop,
    boolList& cutFaces
)
{
    forAll(cellFeat.faces(), faceI)
    {
        face f1 = cellFeat.faces()[faceI];
        label i1 = findIndex(f1, start);
        if (i1 == -1)
        {
            continue;
        }

        label next = f1.fcIndex(i1);
        next = f1.fcIndex(next);
        if (f1[next] == end)
        {
            //Mark all side faces as being cut
            const labelList& sideFaces = cellFeat.faceMap()[faceI];
            forAll(sideFaces, sFI)
            {
                cutFaces[sideFaces[sFI]] = true;
            }

            if (cellFeat.faceMap()[faceI].size() > 1)
            {
                labelHashSet sFaces(cellFeat.faceMap()[faceI]);
                const labelList& pFaces = mesh_.pointFaces()[start];
                forAll(pFaces, pFI)
                {
                    bool found = false;
                    if (sFaces.found(pFaces[pFI]))
                    {
                        face f2 = mesh_.faces()[pFaces[pFI]];
                        label i2 = findIndex(f2, start);
                        forAll(f2, fp)
                        {
                            i2 = f2.fcIndex(i2);
                            if (!featureEdgePoints[pointMap[f2[i2]]])
                            {
                                loop.append(f2[i2]);
                                found = true;
                                break;
                            }
                        }
                    }
                    if (found)
                    {
                        break;
                    }
                }
            }
            break;
        }
    }

    return;
}


void Foam::meshRefinement::walkEdgeDir
(
    const label celli,
    const vector dir,
    const label start,
    const label end,
    DynamicList<label>& loop
)
{
    //Search only on cell edges
    const labelHashSet cellEdgesSet(mesh_.cellEdges()[celli]);

    label next = start;
    while (true)
    {
        const labelList& pEdges = mesh_.pointEdges()[next];
        bool foundEnd = false;
        forAll(pEdges, pEI)
        {
            label meshedgei = pEdges[pEI];
            if (!cellEdgesSet.found(meshedgei))
            {
                continue;
            }
            edge e = mesh_.edges()[meshedgei];

            label otherPt = (e[0] ==  next ? e[1] : e[0]);

            if (otherPt == end)
            {
                foundEnd = true;
                break;
            }
            else
            {
                vector edir = mesh_.points()[otherPt] - mesh_.points()[next];
                edir /= (mag(edir) + SMALL);
                if ((edir & dir) > 0.707)
                {
                    loop.append(otherPt);
                    next = otherPt;
                    break;
                }
            }
        }

        if (foundEnd)
        {
            break;
        }
    }
}


Foam::labelList Foam::meshRefinement::pyrSplitSideFace
(
    const cellFeatures& cellFeat,
    const label& celli,
    const cell& c,
    const label& pyrStart,
    const label& pyrMid,
    const label& pyrEnd,
    boolList& cutFaces
)
{
    DynamicList<label> loop(8);

    const cell superCell(identity(6));
    labelList superPts = superCell.labels(cellFeat.faces());

    labelList cPts = c.labels(mesh_.faces());
    Map<label> pointMap(cPts.size());
    forAll(cPts, cPtI)
    {
        label meshPointI = cPts[cPtI];
        pointMap.insert(meshPointI,cPtI);
    }

    edgeList cEdges = c.edges(mesh_.faces());
    const labelHashSet& fEdges = cellFeat.featureEdge();
    boolList featureEdgePoints(cPts.size(), false);
    forAll(cEdges, cEI)
    {
        edge e = cEdges[cEI];
        label meshedgei = meshTools::findEdge
        (
            mesh_.edges(),
            mesh_.pointEdges()[e[0]],
            e[0],
            e[1]
         );
        if (fEdges.found(meshedgei))
        {
            featureEdgePoints[pointMap[e[0]]] = true;
            featureEdgePoints[pointMap[e[1]]] = true;
        }
    }

    vector dir = mesh_.points()[superPts[pyrStart]]
        - mesh_.points()[superPts[pyrMid]];
    dir /= (mag(dir) + SMALL);

    loop.append(superPts[pyrStart]);
    loop.append(superPts[pyrEnd]);

    diagonalCut
    (
        cellFeat,
        pointMap,
        featureEdgePoints,
        superPts[pyrEnd],
        superPts[pyrMid],
        loop,
        cutFaces
    );

    loop.append(superPts[pyrMid]);

    walkEdgeDir
    (
        celli,
        dir,
        superPts[pyrMid],
        superPts[pyrStart],
        loop
     );

    return loop.shrink();
}


Foam::labelList Foam::meshRefinement::pyrSplit
(
    const cellFeatures& cellFeat,
    const boolList& cutNodes,
    const label& celli,
    const cell& c,
    const label& pyrStart,
    const label& pyrEnd,
    boolList& cutFaces
)
{
    DynamicList<label> loop(8);

    const cell superCell(identity(6));
    labelList superPts = superCell.labels(cellFeat.faces());
    edgeList superEdges = superCell.edges(cellFeat.faces());

    Map<label> pointMapSuperCell(superPts.size());

    forAll(superPts, cPtI)
    {
        label meshPointI = superPts[cPtI];
        pointMapSuperCell.insert(meshPointI,cPtI);
    }

    labelList cPts = c.labels(mesh_.faces());
    Map<label> pointMap(cPts.size());
    forAll(cPts, cPtI)
    {
        label meshPointI = cPts[cPtI];
        pointMap.insert(meshPointI,cPtI);
    }

    edgeList cEdges = c.edges(mesh_.faces());
    const labelHashSet& fEdges = cellFeat.featureEdge();
    boolList featureEdgePoints(cPts.size(), false);
    forAll(cEdges, cEI)
    {
        edge e = cEdges[cEI];
        label meshedgei = meshTools::findEdge
        (
            mesh_.edges(),
            mesh_.pointEdges()[e[0]],
            e[0],
            e[1]
         );
        if (fEdges.found(meshedgei))
        {
            featureEdgePoints[pointMap[e[0]]] = true;
            featureEdgePoints[pointMap[e[1]]] = true;
        }
    }

    label otherPt1 = -1;
    label otherPt2 = -1;
    label start = superPts[pyrStart];

    bool foundStart = false;

    forAll(superEdges, cEI)
    {
        edge e = superEdges[cEI];
        if
        (
            !foundStart && cutNodes[pointMapSuperCell[e[0]]]
            && cutNodes[pointMapSuperCell[e[1]]]
        )
        {
            otherPt1 = (e[0] ==  start ? e[1] : e[0]);
            foundStart =  true;
        }
        else if
        (
            foundStart && cutNodes[pointMapSuperCell[e[0]]]
            && cutNodes[pointMapSuperCell[e[1]]]
        )
        {
            otherPt2 = (e[0] ==  start ? e[1] : e[0]);
            break;
        }
    }

    vector dir1 = mesh_.points()[otherPt1] - mesh_.points()[start];
    dir1 /= (mag(dir1) + SMALL);

    vector dir2 = mesh_.points()[start] - mesh_.points()[otherPt2];
    dir2 /= (mag(dir2) + SMALL);

    loop.append(start);
    walkEdgeDir
    (
        celli,
        dir1,
        start,
        otherPt1,
        loop
     );

    start = otherPt1;
    loop.append(otherPt1);
    diagonalCut
    (
        cellFeat,
        pointMap,
        featureEdgePoints,
        start,
        superPts[pyrEnd],
        loop,
        cutFaces
    );

    start = superPts[pyrEnd];
    loop.append(start);
    diagonalCut
    (
        cellFeat,
        pointMap,
        featureEdgePoints,
        start,
        otherPt2,
        loop,
        cutFaces
    );

    start = otherPt2;
    loop.append(start);
    walkEdgeDir
    (
        celli,
        dir2,
        start,
        superPts[pyrStart],
        loop
     );

    return loop.shrink();
}


Foam::labelList Foam::meshRefinement::prismSplit
(
    const cellFeatures& cellFeat,
    const boolList& cutNodes,
    const label& celli,
    const cell& c,
    boolList& cutFaces
)
{
    DynamicList<label> loop(8);

    const cell superCell(identity(6));
    labelList superPts = superCell.labels(cellFeat.faces());
    edgeList superEdges = superCell.edges(cellFeat.faces());

    Map<label> pointMapSuperCell(superPts.size());

    forAll(superPts, cPtI)
    {
        label meshPointI = superPts[cPtI];
        pointMapSuperCell.insert(meshPointI,cPtI);
    }

    labelList cPts = c.labels(mesh_.faces());
    Map<label> pointMap(cPts.size());
    forAll(cPts, cPtI)
    {
        label meshPointI = cPts[cPtI];
        pointMap.insert(meshPointI,cPtI);
    }

    edgeList cEdges = c.edges(mesh_.faces());
    const labelHashSet& fEdges = cellFeat.featureEdge();
    boolList featureEdgePoints(cPts.size(), false);
    forAll(cEdges, cEI)
    {
        edge e = cEdges[cEI];
        label meshedgei = meshTools::findEdge
        (
            mesh_.edges(),
            mesh_.pointEdges()[e[0]],
            e[0],
            e[1]
         );
        if (fEdges.found(meshedgei))
        {
            featureEdgePoints[pointMap[e[0]]] = true;
            featureEdgePoints[pointMap[e[1]]] = true;
        }
    }

    edge firstEdge(-1, -1), secondEdge(-1, -1);
    bool foundStart = false;

    forAll(superEdges, cEI)
    {
        edge e = superEdges[cEI];
        if
        (
            !foundStart && cutNodes[pointMapSuperCell[e[0]]]
            && cutNodes[pointMapSuperCell[e[1]]]
        )
        {
            firstEdge = edge(e[0], e[1]);
            foundStart =  true;
        }
        else if
        (
            foundStart && cutNodes[pointMapSuperCell[e[0]]]
            && cutNodes[pointMapSuperCell[e[1]]]
        )
        {
            secondEdge = edge(e[0], e[1]);
            break;
        }
    }

    label start = firstEdge[0];
    label end = firstEdge[1];
    vector dir = mesh_.points()[end] - mesh_.points()[start];
    dir /= (mag(dir) + SMALL);

    //swap second edge direction if required
    vector dir2 = mesh_.points()[secondEdge[1]]
        - mesh_.points()[secondEdge[0]];

    dir2 /= (mag(dir2) + SMALL);
    if ((dir & dir2) > 0)
    {
        label tmp = secondEdge[0];
        secondEdge[0] = secondEdge[1];
        secondEdge[1] = tmp;
        dir2 = -dir2;
    }

    loop.append(start);

    walkEdgeDir
    (
        celli,
        dir,
        start,
        end,
        loop
     );
    loop.append(end);

    start = end;
    diagonalCut
    (
        cellFeat,
        pointMap,
        featureEdgePoints,
        start,
        secondEdge[0],
        loop,
        cutFaces
    );

    start = secondEdge[0];
    loop.append(start);

    walkEdgeDir
    (
        celli,
        dir2,
        start,
        secondEdge[1],
        loop
     );
    loop.append(secondEdge[1]);

    start = secondEdge[1];
    diagonalCut
    (
        cellFeat,
        pointMap,
        featureEdgePoints,
        start,
        firstEdge[0],
        loop,
        cutFaces
    );

    return loop.shrink();
}


Foam::labelList Foam::meshRefinement::tetSplit
(
    const cellFeatures& cellFeat,
    const boolList& cutNodes,
    const cell& c,
    boolList& cutFaces
)
{
    DynamicList<label> loop(6);

    labelList cPts = c.labels(mesh_.faces());
    Map<label> pointMap(cPts.size());
    forAll(cPts, cPtI)
    {
        label meshPointI = cPts[cPtI];
        pointMap.insert(meshPointI,cPtI);
    }

    DynamicList<label> cNodes;
    labelList cFPts = cellFeat.labels();
    forAll(cFPts, cPtI)
    {
        label meshPointI = cFPts[cPtI];
        if (cutNodes[cPtI])
        {
            cNodes.append(meshPointI);
        }
    }
    cNodes.shrink();

    edgeList cEdges = c.edges(mesh_.faces());
    const labelHashSet& fEdges = cellFeat.featureEdge();
    boolList featureEdgePoints(cPts.size(), false);
    forAll(cEdges, cEI)
    {
        edge e = cEdges[cEI];
        label meshedgei = meshTools::findEdge
        (
            mesh_.edges(),
            mesh_.pointEdges()[e[0]],
            e[0],
            e[1]
         );
        if (fEdges.found(meshedgei))
        {
            featureEdgePoints[pointMap[e[0]]] = true;
            featureEdgePoints[pointMap[e[1]]] = true;
        }
    }

    for (int i = 0; i < 3; i++)
    {
        label start = cNodes[i];
        loop.append(start);

        if (c.size() != 6)
        {
            label end = (i ==  2 ? cNodes[0] : cNodes[i+1]);

            diagonalCut
            (
                cellFeat,
                pointMap,
                featureEdgePoints,
                start,
                end,
                loop,
                cutFaces
            );
        }
    }

    return loop.shrink();
}


void Foam::meshRefinement::splitPyramids
(
    const List<labelList>& pyrStartEnd,
    const List<scalarField>& cellWeightsPyr,
    List<labelList>& cellLoopsPyr,
    labelList& cutCellsPyr,
    pointField& movedPoints
)
{
    //split manifold face before splitting cells
    {
        polyTopoChange meshMod(mesh_);
        forAll(cutCellsPyr, i)
        {
            labelList startEnd = pyrStartEnd[i];
            label splitFace = -1;

            label celli = cutCellsPyr[i];
            cell c = mesh_.cells()[celli];

            forAll(c, cI)
            {
                label meshFaceI = c[cI];
                const face& f = mesh_.faces()[meshFaceI];
                labelHashSet fPoints(f);

                if (fPoints.found(startEnd[0]) && fPoints.found(startEnd[1]))
                {
                    splitFace = meshFaceI;
                }
            }

            if (splitFace == -1)
            {
                FatalErrorInFunction
                    << "Cannot find face to"
                    << " split for pyramid splitting of cell: " << celli
                    << " using anchors " <<  pyrStartEnd[i]
                    << abort(FatalError);
            }

            face sFace = mesh_.faces()[splitFace];

            DynamicList<label> newFace(sFace.size());
            DynamicList<label> modFace(sFace.size());

            label next = findIndex(sFace, startEnd[0]);
            newFace.append(sFace[next]);
            forAll(sFace, sFI)
            {
                next = sFace.fcIndex(next);
                newFace.append(sFace[next]);
                if (sFace[next] == startEnd[1])
                {
                    break;
                }
            }
            newFace.shrink();

            next = findIndex(sFace, startEnd[1]);
            modFace.append(sFace[next]);
            forAll(sFace, sFI)
            {
                next = sFace.fcIndex(next);
                modFace.append(sFace[next]);
                if (sFace[next] == startEnd[0])
                {
                    break;
                }
            }
            modFace.shrink();

            label own = mesh_.faceOwner()[splitFace];
            label patchID = mesh_.boundaryMesh().whichPatch(splitFace);
            label nei = (patchID == -1 ? mesh_.faceNeighbour()[splitFace] : -1);
            label zoneID = mesh_.faceZones().whichZone(splitFace);
            bool zoneFlip = false;
            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(splitFace)];
            }

            meshMod.modifyFace
            (
                face(modFace),    //vertices
                splitFace,  // label of face being modified
                own,        // owner
                nei,        // neighbour
                false,      // face flip
                patchID,    // new patch for face
                zoneID,     // zone for face
                zoneFlip    // face flip in zone
            );

            meshMod.addFace
            (
                face(newFace),    // vertices
                own,        // owner,
                nei,        // neighbour,
                -1,         // masterPointID,
                -1,         // masterEdgeID,
                splitFace,  // masterFaceID,
                false,      // flipFaceFlux,
                patchID,    // patchID,
                zoneID,     // zoneID,
                zoneFlip    // zoneFlip
             );
        }

        // Create mesh (no inflation), return map from old to new mesh.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

        // Update fields
        mesh_.updateMesh(map);

        // Optionally inflate mesh
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersection info
        updateMesh(map, getChangedFaces(map, labelList(0)));


        forAll(cellLoopsPyr, i)
        {
            cutCellsPyr[i] = map().reverseCellMap()[cutCellsPyr[i]];
            labelList& cLoop = cellLoopsPyr[i];
            forAll(cLoop, cLI)
            {
                cLoop[cLI] = map().reversePointMap()[cLoop[cLI]];
            }
        }

        pointField tPoints(movedPoints);
        forAll(map().reversePointMap(), pointI)
        {
            tPoints[pointI] = movedPoints[map().pointMap()[pointI]];
        }
        movedPoints = tPoints;
    }

    //Now split the cells
    {
        forAll(cutCellsPyr, i)
        {
            label celli = cutCellsPyr[i];
            cell c = mesh_.cells()[celli];
            const labelList& meshPoints = c.labels(mesh_.faces());
            labelHashSet mPoints(meshPoints);

            const labelList loop = cellLoopsPyr[i];

            label nFound = 0;

            forAll(loop, j)
            {
                if (mPoints.found(loop[j]))
                {
                    nFound++;
                }
            }
            if (nFound != loop.size())
            {
                labelList cellCells = mesh_.cellCells()[celli];
                forAll(cellCells, j)
                {
                    label nei = cellCells[j];
                    cell neiC = mesh_.cells()[nei];
                    const labelList& neiPoints = neiC.labels(mesh_.faces());
                    labelHashSet nPoints(neiPoints);
                    nFound = 0;

                    forAll(loop, j)
                    {
                        if (nPoints.found(loop[j]))
                        {
                            nFound++;
                        }
                    }
                    if (nFound == loop.size())
                    {
                        cutCellsPyr[i] = nei;
                        break;
                    }
                }
            }
        }

        cellCuts cuts(mesh_, cutCellsPyr, cellLoopsPyr, cellWeightsPyr);

        Foam::meshCutter mCut(mesh_);

        polyTopoChange meshMod(mesh_);
        mCut.setRefinement(cuts, meshMod);

        // Create mesh (no inflation), return map from old to new mesh.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

        // Update fields
        mesh_.updateMesh(map);

        // Optionally inflate mesh
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersection info
        updateMesh(map, getChangedFaces(map, labelList(0)));

        pointField tPoints(movedPoints);
        forAll(map().reversePointMap(), pointI)
        {
            tPoints[pointI] = movedPoints[map().pointMap()[pointI]];
        }
        movedPoints = tPoints;
    }
}


void Foam::meshRefinement::selectCutType
(
    const boolList& excludeTetSplit,
    const boolList& excludedCells,
    const boolList& cutNodes,
    const vectorField& cutNormals,
    DynamicList<label>& cutCells,
    DynamicList<label>& cutCellsPyr,
    DynamicList<labelList>& cellLoops,
    DynamicList<labelList>& cellLoopsPyr,
    DynamicList<labelList>& pyrStartEnd
)
{
    label nTets = 0;
    label nPyrs = 0;
    label nPrisms = 0;

    boolList cutFaces(mesh_.nFaces(), false);
    DynamicList<label> cut54;
    DynamicList<label> cut66;
    DynamicList<label> cut67;

    //Record split faces where all nodes are cutNodes
    labelList facesAllNodesCut(mesh_.nFaces(), 0);

    forAll(mesh_.cells(), celli)
    {
        const cell& c = mesh_.cells()[celli];

        cellFeatures cellFeat(mesh_, 0.707, celli);
        if (excludedCells[celli] || cellFeat.faces().size() != 6)
        {
            continue;
        }

        const cell superCell(identity(6));
        labelList superPts = superCell.labels(cellFeat.faces());
        edgeList superEdges = superCell.edges(cellFeat.faces());

        label nCut = 0;
        boolList cutSuperCell(superPts.size(), false);
        Map<label> pointMap(superPts.size());
        forAll(superPts, cPtI)
        {
            label meshPointI = superPts[cPtI];

            pointMap.insert(meshPointI,cPtI);
            if (cutNodes[meshPointI])
            {
                cutSuperCell[cPtI] = true;
                nCut++;
            }
        }

        labelList nConn(superPts.size(), 0);
        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (cutNodes[e[0]] && cutNodes[e[1]])
            {
                nConn[pointMap[e[0]]]++;
                nConn[pointMap[e[1]]]++;
            }
        }

        label nConEdges = 0;
        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (cutNodes[e[0]] && cutNodes[e[1]])
            {
                nConEdges++;
            }
        }

        /*
        labelList cellMeshPts = c.labels(mesh_.faces());
        label nTotalCuts = 0;
        forAll(cellMeshPts, cPtI)
        {
            if (cutNodes[cellMeshPts[cPtI]])
            {
                nTotalCuts++;
            }
        }
        */

        bool splitCell = false;
        if
        (
            (nCut == 3 && nConEdges == 0)
            || (nCut == 4 && nConEdges == 0)
            || (nCut == 4 && nConEdges == 2)
            || (nCut == 5 && nConEdges == 4)
            || (nCut == 6 && nConEdges == 6)
            || (nCut == 6 && nConEdges == 7)
         )
        {
            splitCell = true;
        }

        if (splitCell)
        {
            if (nCut == 3 && nConEdges == 0 && !excludeTetSplit[celli])
            {
                nTets++;
                cutCells.append(celli);
                cellLoops.append(tetSplit(cellFeat, cutSuperCell, c, cutFaces));
            }
            else if (nCut == 4 && nConEdges == 0)
            {
                DynamicList<label> cut4(4);
                forAll(superPts, cPtI)
                {
                    label meshPointI = superPts[cPtI];
                    if (cutNodes[meshPointI])
                    {
                        cut4.append(meshPointI);
                    }
                }
                cut4.shrink();

                List<face> possibleFaces(4);

                if (cut4.size() == 4)
                {
                    forAll(cut4, i)
                    {
                        label np = 0;
                        face f(3);
                        forAll(cut4, j)
                        {
                            if (j != i)
                            {
                                f[np++] = cut4[j];
                            }
                        }
                        possibleFaces[i] = f;
                    }
                    scalar maxDotProd = -GREAT;
                    label maxFace = -1;
                    forAll(possibleFaces, faceI)
                    {
                        face f = possibleFaces[faceI];
                        vector normal = f.unitNormal(mesh_.points());

                        vector aveSurfaceNormal = vector::zero;

                        label nCut = 0;
                        forAll(f, fp)
                        {
                            if (cutNormals[f[fp]] != vector(GREAT, GREAT, GREAT))
                            {
                                nCut++;
                                aveSurfaceNormal += cutNormals[f[fp]];
                            }
                        }
                        if (nCut > 0)
                        {
                            aveSurfaceNormal /= nCut;
                            aveSurfaceNormal /= (mag(aveSurfaceNormal) + SMALL);
                            scalar dotProd = mag(aveSurfaceNormal & normal);
                            if (dotProd > maxDotProd)
                            {
                                maxDotProd = dotProd;
                                maxFace = faceI;
                            }
                        }
                    }
                    if (maxFace != -1)
                    {
                        nTets++;
                        cutCells.append(celli);
                        boolList newCutSuperCell(superPts.size(), false);
                        face f = possibleFaces[maxFace];
                        forAll(f, fp)
                        {
                            newCutSuperCell[pointMap[f[fp]]] = true;
                        }

                        cellLoops.append
                        (
                            tetSplit(cellFeat, newCutSuperCell, c, cutFaces)
                        );
                    }
                }
            }
            else if (nCut == 4 && nConEdges == 2)
            {
                bool prism = true;
                label pyrStart = -1;
                label pyrEnd = -1;
                label pyrMid = -1;

                forAll(superPts, cPtI)
                {
                    if (nConn[cPtI] == 2)
                    {
                        prism = false;
                        pyrStart = cPtI;
                    }
                    if (nConn[cPtI] == 0 && cutSuperCell[cPtI])
                    {
                        pyrEnd = cPtI;
                    }
                }

                forAll(superEdges, cEI)
                {
                    edge e = superEdges[cEI];
                    if (pointMap[e[0]] == pyrStart)
                    {
                        if (!cutSuperCell[pointMap[e[1]]])
                        {
                            pyrMid = pointMap[e[1]];
                            break;
                        }
                    }
                    else if (pointMap[e[1]] == pyrStart)
                    {
                        if (!cutSuperCell[pointMap[e[0]]])
                        {
                            pyrMid = pointMap[e[0]];
                            break;
                        }
                    }
                }

                if (prism)
                {
                    nPrisms++;
                    cutCells.append(celli);
                    cellLoops.append
                    (
                       prismSplit(cellFeat, cutSuperCell, celli, c, cutFaces)
                    );
                }
                else
                {
                    bool pyr = true;
                    vectorField eNormals(2);
                    pointField ePts(2);
                    bool firstPass = true;
                    bool foundEdge = false;

                    boolList newCutSuperCell(superPts.size(), false);
                    newCutSuperCell[pyrEnd] = true;

                    forAll(superPts, cPtI)
                    {
                        if (nConn[cPtI] == 1)
                        {
                            if (firstPass)
                            {
                                eNormals[0] = cutNormals[superPts[cPtI]];
                                ePts[0] = mesh_.points()[superPts[cPtI]];
                                newCutSuperCell[cPtI] = true;
                                firstPass = false;
                            }
                            else
                            {
                                eNormals[1] = cutNormals[superPts[cPtI]];
                                ePts[1] = mesh_.points()[superPts[cPtI]];
                                newCutSuperCell[cPtI] = true;
                                foundEdge = true;
                            }
                        }
                    }

                    if
                    (
                        foundEdge
                        && eNormals[0] != vector(GREAT, GREAT, GREAT)
                        && eNormals[1] != vector(GREAT, GREAT, GREAT)
                    )
                    {
                        vector dir =  ePts[0] - mesh_.cellCentres()[celli];
                        dir +=  ePts[1] - mesh_.cellCentres()[celli];
                        dir /= 2.0;

                        if ((dir & eNormals[0]) > 0.)
                        {
                            eNormals[0] = -eNormals[0];
                        }

                        if ((dir & eNormals[1]) > 0.)
                        {
                            eNormals[1] = -eNormals[1];
                        }
                        vector eDir = ePts[1] - ePts[0];
                        eDir /= (mag(eDir) + SMALL);

                        eNormals /= (mag(eNormals) + SMALL);

                        if
                        (
                            ((eNormals[0] & eDir) - (eNormals[1] & eDir))
                            > -0.1736
                        )
                        {
                            pyr = false;
                        }
                    }

                    if (pyr)//if(!refinedCell && pyr)
                    {
                        nPyrs++;
                        cutCells.append(celli);
                        cellLoops.append
                        (
                            pyrSplit
                            (
                                cellFeat,
                                cutSuperCell,
                                celli,
                                c,
                                pyrStart,
                                pyrEnd,
                                cutFaces
                             )
                         );

                        cutCellsPyr.append(celli);
                        cellLoopsPyr.append
                        (
                            pyrSplitSideFace
                            (
                                cellFeat,
                                celli,
                                c,
                                pyrStart,
                                pyrMid,
                                pyrEnd,
                                cutFaces
                             )
                         );
                        labelList startEnd(2);
                        startEnd[0] = superPts[pyrStart];
                        startEnd[1] = superPts[pyrEnd];
                        pyrStartEnd.append(startEnd);
                    }
                    else if (!excludeTetSplit[celli])
                    {
                        nTets++;
                        cutCells.append(celli);
                        cellLoops.append
                        (
                            tetSplit(cellFeat, newCutSuperCell, c, cutFaces)
                        );
                    }
                }
            }
            else if (nCut == 6 && nConEdges == 6)
            {
                bool prism = false;
                forAll(superPts, cPtI)
                {
                    if (nConn[cPtI] == 1 || nConn[cPtI] == 3)
                    {
                        prism = true;
                    }
                }

                if (prism)
                {
                    forAll(c, cFI)
                    {
                        label faceI = c[cFI];
                        face f = mesh_.faces()[faceI];
                        bool allCuts = true;
                        forAll(f, fp)
                        {
                            if (!cutNodes[f[fp]])
                            {
                                allCuts = false;
                                break;
                            }
                        }
                        if (allCuts)
                        {
                            facesAllNodesCut[faceI]++;
                            break;
                        }
                    }
                    cut66.append(celli);
                }
            }
            else if (nCut == 5 && nConEdges == 4)
            {
                cut54.append(celli);
            }
            else if (nCut == 6 && nConEdges == 7)
            {
                cut67.append(celli);
            }
        }
    }

    syncTools::syncFaceList(mesh_, facesAllNodesCut, plusEqOp<label>());

    cut66.shrink();
    forAll(cut66, i)
    {
        label celli = cut66[i];
        const cell& c = mesh_.cells()[celli];

        bool validCut = true;
        forAll(c, cFI)
        {
            label faceI = c[cFI];
            face f = mesh_.faces()[faceI];
            bool allCuts = true;
            forAll(f, fp)
            {
                if (!cutNodes[f[fp]])
                {
                    allCuts = false;
                    break;
                }
            }
            if (allCuts)
            {
                if (facesAllNodesCut[faceI] == 2)
                {
                    validCut = false;
                }
                break;
            }
        }

        if (validCut)
        {
            cellFeatures cellFeat(mesh_, 0.707, celli);
            const cell superCell(identity(6));
            labelList superPts = superCell.labels(cellFeat.faces());
            edgeList superEdges = superCell.edges(cellFeat.faces());

            Map<label> pointMap(superPts.size());
            forAll(superPts, cPtI)
            {
                label meshPointI = superPts[cPtI];
                pointMap.insert(meshPointI,cPtI);
            }

            labelList nConn(superPts.size(), 0);
            forAll(superEdges, cEI)
            {
                edge e = superEdges[cEI];
                if (cutNodes[e[0]] && cutNodes[e[1]])
                {
                    nConn[pointMap[e[0]]]++;
                    nConn[pointMap[e[1]]]++;
                }
            }

            boolList newCutSuperCell(superPts.size(), false);
            forAll(superPts, cPtI)
            {
                if (nConn[cPtI] == 1 || nConn[cPtI] == 3)
                {
                    newCutSuperCell[cPtI] = true;
                }
            }

            nPrisms++;
            cutCells.append(celli);
            cellLoops.append
            (
               prismSplit(cellFeat, newCutSuperCell, celli, c, cutFaces)
             );
        }
    }

    cut67.shrink();
    pointField start(cut67.size());
    pointField end(cut67.size());
    forAll(cut67, i)
    {
        label celli = cut67[i];
        point cc = mesh_.cellCentres()[celli];

        cellFeatures cellFeat(mesh_, 0.707, celli);
        const cell superCell(identity(6));
        labelList superPts = superCell.labels(cellFeat.faces());
        edgeList superEdges = superCell.edges(cellFeat.faces());

        Map<label> pointMap(superPts.size());
        forAll(superPts, cPtI)
        {
            label meshPointI = superPts[cPtI];
            pointMap.insert(meshPointI,cPtI);
        }

        labelList nConn(superPts.size(), 0);
        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (cutNodes[e[0]] && cutNodes[e[1]])
            {
                nConn[pointMap[e[0]]]++;
                nConn[pointMap[e[1]]]++;
            }
        }

        point edgeCentre(vector::zero);
        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (nConn[pointMap[e[0]]] == 3 && nConn[pointMap[e[1]]] == 3)
            {
                edgeCentre = 0.5*(mesh_.points()[e[0]] + mesh_.points()[e[1]]);
            }
        }

        vector dir = cc - edgeCentre;

        start[i] = edgeCentre + 0.66*dir;
        end[i] = start[i] + 0.66*dir;
    }


    const labelList surfacesToBaffle(identity(surfaces_.surfaces().size()));

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;

    surfaces_.findNearestIntersection
    (
        surfacesToBaffle,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
     );


    vectorField normals67(start.size(), vector(GREAT, GREAT, GREAT));

    forAll(surfacesToBaffle, sI)
    {
        label surfI = surfacesToBaffle[sI];
        DynamicList<pointIndexHit> localHits;
        forAll(surface1, pointI)
        {
            if (surface1[pointI] == surfI)
            {
                localHits.append(hit1[pointI]);
            }
        }
        localHits.shrink();

        pointField localNormals;
        label geomI = surfaces_.surfaces()[surfI];
        surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

        label localI = 0;
        forAll(surface1, pointI)
        {
            if (surface1[pointI] == surfI)
            {
                normals67[pointI] = localNormals[localI];
                normals67[pointI] /= (mag(normals67[pointI]) + SMALL);
                localI++;
            }
        }
    }

    forAll(cut67, i)
    {
        label celli = cut67[i];
        const cell& c = mesh_.cells()[celli];
        point cc = mesh_.cellCentres()[celli];

        cellFeatures cellFeat(mesh_, 0.707, celli);

        const cell superCell(identity(6));
        labelList superPts = superCell.labels(cellFeat.faces());
        edgeList superEdges = superCell.edges(cellFeat.faces());

        boolList cutSuperCell(superPts.size(), false);
        Map<label> pointMap(superPts.size());
        forAll(superPts, cPtI)
        {
            label meshPointI = superPts[cPtI];
            pointMap.insert(meshPointI,cPtI);
            if (cutNodes[meshPointI])
            {
                cutSuperCell[cPtI] = true;
            }
        }

        labelList nConn(superPts.size(), 0);
        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (cutNodes[e[0]] && cutNodes[e[1]])
            {
                nConn[pointMap[e[0]]]++;
                nConn[pointMap[e[1]]]++;
            }
        }

        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (nConn[pointMap[e[0]]] == 3 && nConn[pointMap[e[1]]] == 3)
            {
                cutSuperCell[pointMap[e[0]]] = false;
                cutSuperCell[pointMap[e[1]]] = false;
            }
        }

        vector fNorm = cc - start[i];
        fNorm /= (mag(fNorm) + SMALL);

        if (hit1[i].hit())
        {
            if (mag(fNorm&normals67[i]) > 0.923)
            {
                nPrisms++;
                cutCells.append(celli);
                cellLoops.append
                (
                    prismSplit(cellFeat, cutSuperCell, celli, c, cutFaces)
                );
            }
        }
    }

    syncTools::syncFaceList(mesh_, cutFaces, orEqOp<bool>());
    cut54.shrink();
    while (true)
    {
        DynamicList<label> newCut54;
        label nFound = 0;

        forAll(cut54, i)
        {
            label celli = cut54[i];

            const cell& c = mesh_.cells()[celli];

            cellFeatures cellFeat(mesh_, 0.707, celli);
            const cell superCell(identity(6));
            labelList superPts = superCell.labels(cellFeat.faces());
            edgeList superEdges = superCell.edges(cellFeat.faces());
            faceList superFaces = cellFeat.faces();

//            label nCut = 0;
            boolList cutSuperCell(superPts.size(), false);
            Map<label> pointMap(superPts.size());
            forAll(superPts, cPtI)
            {
                label meshPointI = superPts[cPtI];
                pointMap.insert(meshPointI,cPtI);
                if (cutNodes[meshPointI])
                {
                    cutSuperCell[cPtI] = true;
//                    nCut++;
                }
            }

            labelList nConn(superPts.size(), 0);
            forAll(superEdges, cEI)
            {
                edge e = superEdges[cEI];
                if (cutNodes[e[0]] && cutNodes[e[1]])
                {
                    nConn[pointMap[e[0]]]++;
                    nConn[pointMap[e[1]]]++;
                }
            }

            boolList visited = cutSuperCell;
            label nCutFaces = 0;
            forAll(superFaces, sFI)
            {
                face superFace = superFaces[sFI];
                if (cutFaces[cellFeat.faceMap()[sFI][0]])
                {
                    nCutFaces++;
                    forAll(superFace, fpI)
                    {
                        visited[pointMap[superFace[fpI]]] = false;
                    }
                }
            }

            bool prism = true;
            if (nCutFaces > 1)
            {
                forAll(visited, cPtI)
                {
                    if (visited[cPtI])
                    {
                        prism = false;
                        break;
                    }
                }
            }

            if (prism)
            {
                boolList newCutSuperCell(superPts.size(), false);

                forAll(superEdges, cEI)
                {
                    edge e = superEdges[cEI];
                    if (cutNodes[e[0]] && cutNodes[e[1]])
                    {
                        if
                        (
                            (nConn[pointMap[e[0]]] == 1)
                            || (nConn[pointMap[e[1]]] == 1)
                         )
                        {
                            newCutSuperCell[pointMap[e[0]]] = true;
                            newCutSuperCell[pointMap[e[1]]] = true;
                        }
                    }
                }

                nPrisms++;
                cutCells.append(celli);
                cellLoops.append
                (
                    prismSplit(cellFeat, newCutSuperCell, celli, c, cutFaces)
                );
                nFound++;
            }
            else
            {
                newCut54.append(celli);
            }
        }

        newCut54.shrink();
        cut54.clear();
        cut54 = newCut54;

        syncTools::syncFaceList(mesh_, cutFaces, orEqOp<bool>());
        if (returnReduce(nFound, sumOp<label>()) == 0)
        {
            break;
        }
    }


    forAll(cut54, i)
    {
        label celli = cut54[i];
        const cell& c = mesh_.cells()[celli];

        cellFeatures cellFeat(mesh_, 0.707, celli);
        const cell superCell(identity(6));
        labelList superPts = superCell.labels(cellFeat.faces());
        edgeList superEdges = superCell.edges(cellFeat.faces());
        faceList superFaces = cellFeat.faces();

//        label nCut = 0;
        boolList cutSuperCell(superPts.size(), false);
        Map<label> pointMap(superPts.size());
        forAll(superPts, cPtI)
        {
            label meshPointI = superPts[cPtI];
            pointMap.insert(meshPointI,cPtI);
            if (cutNodes[meshPointI])
            {
                cutSuperCell[cPtI] = true;
//                nCut++;
            }
        }

        labelList nConn(superPts.size(), 0);
        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (cutNodes[e[0]] && cutNodes[e[1]])
            {
                nConn[pointMap[e[0]]]++;
                nConn[pointMap[e[1]]]++;
            }
        }

        boolList visited = cutSuperCell;
        forAll(superFaces, sFI)
        {
            face superFace = superFaces[sFI];
            if (cutFaces[cellFeat.faceMap()[sFI][0]])
            {
                forAll(superFace, fpI)
                {
                    visited[pointMap[superFace[fpI]]] = false;
                }
            }
        }

        bool prism = true;
        forAll(visited, cPtI)
        {
            if (visited[cPtI])
            {
                prism = false;
                break;
            }
        }

        if (prism)
        {
            boolList newCutSuperCell(superPts.size(), false);

            forAll(superEdges, cEI)
            {
                edge e = superEdges[cEI];
                if (cutNodes[e[0]] && cutNodes[e[1]])
                {
                    if
                    (
                        (nConn[pointMap[e[0]]] == 1)
                        || (nConn[pointMap[e[1]]] == 1)
                     )
                    {
                        newCutSuperCell[pointMap[e[0]]] = true;
                        newCutSuperCell[pointMap[e[1]]] = true;
                    }
                }
            }

            nPrisms++;
            cutCells.append(celli);
            cellLoops.append
            (
                prismSplit(cellFeat, newCutSuperCell, celli, c, cutFaces)
            );
        }
        else
        {
            boolList newCutSuperCell(cutSuperCell);

            forAll(superEdges, cEI)
            {
                edge e = superEdges[cEI];
                if (cutNodes[e[0]] && cutNodes[e[1]])
                {
                    if (nConn[pointMap[e[0]]] == 1)
                    {
                        newCutSuperCell[pointMap[e[1]]] = false;
                    }
                    else if (nConn[pointMap[e[1]]] == 1)
                    {
                        newCutSuperCell[pointMap[e[0]]] = false;
                    }
                }
            }
            if (!excludeTetSplit[celli])
            {
                nTets++;
                cutCells.append(celli);

                cellLoops.append
                (
                    tetSplit(cellFeat, newCutSuperCell, c, cutFaces)
                );
            }
        }
    }

    cutCells.shrink();
    cellLoops.shrink();
    cutCellsPyr.shrink();
    cellLoopsPyr.shrink();
    pyrStartEnd.shrink();

    Info<<"Types of cuts: "<<endl;
    Info<<"Number tetrahedrals : "<<returnReduce(nTets, sumOp<label>())<<endl;
    Info<<"Number pyramids : "<<returnReduce(nPyrs, sumOp<label>())<<endl;
    Info<<"Number prisms : "<<returnReduce(nPrisms, sumOp<label>())<<endl;
}


void Foam::meshRefinement::filterCutCells
(
    List<label>& cutCells,
    List<label>& cutCellsPyr,
    List<labelList>& cellLoops,
    List<labelList>& cellLoopsPyr,
    List<labelList>& pyrStartEnd
)
{
    boolList stopCellCut(mesh_.nCells(), false);
    boolList cutEdges(mesh_.nEdges(), false);

    forAll(cutCells, i)
    {
        const labelList& loop = cellLoops[i];
        forAll(loop, i)
        {
            label start = loop[i];
            label end = -1;
            if (i == loop.size() -1)
            {
                end = loop[0];
            }
            else
            {
                end = loop[i+1];
            }

            label meshedgei = meshTools::findEdge
            (
                mesh_.edges(),
                mesh_.pointEdges()[start],
                start,
                end
             );

            if (meshedgei != -1)
            {
                cutEdges[meshedgei] = true;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        cutEdges,
        orEqOp<bool>(),
        false        // null value
     );

    label nPrevented = 0;
    forAll(cutCells, i)
    {
        label celli = cutCells[i];

        const labelList& loop = cellLoops[i];
        label nCutEdges = 0;
        forAll(loop, i)
        {
            label start = loop[i];
            label end = -1;
            if (i == loop.size() -1)
            {
                end = loop[0];
            }
            else
            {
                end = loop[i+1];
            }

            label meshedgei = meshTools::findEdge
            (
                mesh_.edges(),
                mesh_.pointEdges()[start],
                start,
                end
             );

            if (meshedgei != -1)
            {
                nCutEdges++;
            }
        }

        const cell& c = mesh_.cells()[celli];
        edgeList cEdges = c.edges(mesh_.faces());

        label allCutEdges = 0;
        forAll(cEdges, cEI)
        {
            edge e = cEdges[cEI];
            label meshedgei = meshTools::findEdge
            (
                mesh_.edges(),
                mesh_.pointEdges()[e[0]],
                e[0],
                e[1]
             );
            if (cutEdges[meshedgei])
            {
                allCutEdges++;
            }
        }

        if (nCutEdges > 0 && (allCutEdges - nCutEdges) > 1)
        {
            stopCellCut[celli] = true;
            nPrevented++;
        }
    }

    Info<<"Prevented splitting at "<<returnReduce(nPrevented, sumOp<label>())
        <<" cells"<<endl;

    DynamicList<label> newCutCells(cutCells.size());
    DynamicList<label> newCutCellsPyr(cutCellsPyr.size());
    DynamicList<labelList> newCellLoops(cellLoops.size());
    DynamicList<labelList> newCellLoopsPyr(cellLoopsPyr.size());
    DynamicList<labelList> newPyrStartEnd(pyrStartEnd.size());

    forAll(cutCells, i)
    {
        label celli = cutCells[i];

        if (!stopCellCut[celli])
        {
            newCutCells.append(cutCells[i]);
            newCellLoops.append(cellLoops[i]);
        }
    }

    forAll(cutCellsPyr, i)
    {
        label celli = cutCellsPyr[i];

        if (!stopCellCut[celli])
        {
            newCutCellsPyr.append(cutCellsPyr[i]);
            newCellLoopsPyr.append(cellLoopsPyr[i]);
            newPyrStartEnd.append(pyrStartEnd[i]);
        }
    }

    cutCells.transfer(newCutCells.shrink());
    cutCellsPyr.transfer(newCutCellsPyr.shrink());
    cellLoops.transfer(newCellLoops.shrink());
    cellLoopsPyr.transfer(newCellLoopsPyr.shrink());
    pyrStartEnd.transfer(newPyrStartEnd.shrink());
}


Foam::boolList Foam::meshRefinement::setRefinementCells() const
{
    boolList excludeTetSplit(mesh_.nCells(), false);

    const labelList cellLevel = meshCutter_.cellLevel();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Calculate cell level
    scalarField neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        neiLevel[faceI-mesh_.nInternalFaces()] =
           cellLevel[mesh_.faceOwner()[faceI]];
    }
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);

    forAll(mesh_.cells(), celli)
    {
        const cell& c = mesh_.cells()[celli];
        forAll(c, cfI)
        {
            label meshFaceI = c[cfI];
            if (mesh_.isInternalFace(meshFaceI))
            {
                label own = mesh_.faceOwner()[meshFaceI];
                label nbr
                (
                   own == celli ? mesh_.faceNeighbour()[meshFaceI] : own
                );

                if (cellLevel[celli] > cellLevel[nbr])
                {
                    excludeTetSplit[celli] = true;
                    break;
                }
            }
            else
            {
                label patchI = patches.whichPatch(meshFaceI);
                if (patches[patchI].coupled())
                {
                    if
                    (
                       cellLevel[celli] >
                       neiLevel[meshFaceI-mesh_.nInternalFaces()]
                    )
                    {
                        excludeTetSplit[celli] = true;
                        break;
                    }
                }
            }
        }
    }

    return excludeTetSplit;
}


// split processor faces as meshCutter not parallel aware
void Foam::meshRefinement::splitProcessorFaces
(
    List<label>& cutCells,
    List<label>& cutCellsPyr,
    List<labelList>& cellLoops,
    List<labelList>& cellLoopsPyr,
    List<labelList>& pyrStartEnd,
    pointField& movedPoints
)
{
    polyTopoChange meshMod(mesh_);
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Calculate cut face
    pointField startField
    (
        mesh_.nFaces() - mesh_.nInternalFaces(),
        vector(GREAT, GREAT, GREAT)
    );
    pointField endField
    (
        mesh_.nFaces() - mesh_.nInternalFaces(),
        vector(GREAT, GREAT, GREAT)
    );
    boolList processorCutFaces(mesh_.nFaces(), false);

    forAll(cellLoops, cLI)
    {
        const labelList& loop = cellLoops[cLI];

        forAll(loop, i)
        {
            label start = loop[i];
            label end = -1;
            if (i == loop.size() -1)
            {
                end = loop[0];
            }
            else
            {
                end = loop[i+1];
            }

            label meshedgei = meshTools::findEdge
            (
                mesh_.edges(),
                mesh_.pointEdges()[start],
                start,
                end
             );

            if (meshedgei == -1)
            {
                const labelList& pFaces = mesh_.pointFaces()[start];

                forAll(pFaces, pFI)
                {
                    face pf = mesh_.faces()[pFaces[pFI]];
                    labelHashSet fpSet(pf);
                    if (fpSet.found(start) && fpSet.found(end))
                    {
                        label patchI = patches.whichPatch(pFaces[pFI]);
                        if
                        (
                            patchI != -1
                            && isA<processorPolyPatch>(patches[patchI])
                        )
                        {
                            label startFace = pFaces[pFI]-mesh_.nInternalFaces();
                            startField[startFace] = mesh_.points()[start];
                            endField[startFace] = mesh_.points()[end];
                            processorCutFaces[pFaces[pFI]] = true;
                        }
                    }
                }
            }
        }
    }

    forAll(cellLoopsPyr, cLI)
    {
        const labelList& loop = cellLoopsPyr[cLI];

        forAll(loop, i)
        {
            label start = loop[i];
            label end = -1;
            if (i == loop.size() -1)
            {
                end = loop[0];
            }
            else
            {
                end = loop[i+1];
            }

            label meshedgei = meshTools::findEdge
            (
                mesh_.edges(),
                mesh_.pointEdges()[start],
                start,
                end
             );

            if (meshedgei == -1)
            {
                const labelList& pFaces = mesh_.pointFaces()[start];

                forAll(pFaces, pFI)
                {
                    face pf = mesh_.faces()[pFaces[pFI]];
                    labelHashSet fpSet(pf);
                    if (fpSet.found(start) && fpSet.found(end))
                    {
                        label patchI = patches.whichPatch(pFaces[pFI]);
                        if
                        (
                            patchI != -1
                            && isA<processorPolyPatch>(patches[patchI])
                        )
                        {
                            label startFace = pFaces[pFI]
                                -mesh_.nInternalFaces();
                            startField[startFace] = mesh_.points()[start];
                            endField[startFace] = mesh_.points()[end];
                            processorCutFaces[pFaces[pFI]] = true;
                        }
                    }
                }
            }
        }
    }

    pointField swappedStart(startField);
    pointField swappedEnd(endField);

    syncTools::swapBoundaryFacePositions(mesh_, swappedStart);
    syncTools::swapBoundaryFacePositions(mesh_, swappedEnd);
    syncTools::syncFaceList(mesh_, processorCutFaces, orEqOp<bool>());

    forAll(mesh_.faces(), meshFaceI)
    {
        const face& f = mesh_.faces()[meshFaceI];

        if (processorCutFaces[meshFaceI])
        {
            label patchID = mesh_.boundaryMesh().whichPatch(meshFaceI);
            const polyPatch& pp = mesh_.boundaryMesh()[patchID];
            const processorPolyPatch& ppp =
                dynamic_cast<const processorPolyPatch&>(pp);

            scalar minEdgeLength = GREAT;

            for (const edge e : f.walkEdges())
            {
                scalar eLen = e.mag(mesh_.points());
                minEdgeLength = min(eLen, minEdgeLength);
            }

            bool foundSplit =  false;
            point startPt, endPt;

            bool master = (ppp.myProcNo() < ppp.neighbProcNo() ? true : false);

            if (master)
            {
                if
                (
                    startField[meshFaceI - mesh_.nInternalFaces()]
                    != vector(GREAT, GREAT, GREAT)
                 )
                {
                    startPt = startField[meshFaceI - mesh_.nInternalFaces()];
                    endPt = endField[meshFaceI - mesh_.nInternalFaces()];
                    foundSplit = true;
                }
                else if
                (
                    swappedStart[meshFaceI - mesh_.nInternalFaces()]
                    != vector(GREAT, GREAT, GREAT)
                 )
                {
                    startPt = swappedStart[meshFaceI - mesh_.nInternalFaces()];
                    endPt = swappedEnd[meshFaceI - mesh_.nInternalFaces()];
                    foundSplit = true;
                }
            }
            else
            {
                if
                (
                    swappedStart[meshFaceI - mesh_.nInternalFaces()]
                    == vector(GREAT, GREAT, GREAT)
                 )
                {
                    startPt = startField[meshFaceI - mesh_.nInternalFaces()];
                    endPt = endField[meshFaceI - mesh_.nInternalFaces()];
                    foundSplit = true;
                }
                else if
                (
                    swappedStart[meshFaceI - mesh_.nInternalFaces()]
                    != vector(GREAT, GREAT, GREAT)
                 )
                {
                    startPt = swappedStart[meshFaceI - mesh_.nInternalFaces()];
                    endPt = swappedEnd[meshFaceI - mesh_.nInternalFaces()];
                    foundSplit = true;
                }
            }

            label start = -1;
            label end = -1;
            if (foundSplit)
            {
                bool firstPoint = true;

                label next = 0;
                while (true)
                {
                    if (firstPoint)
                    {
                        if
                        (
                            mag(mesh_.points()[f[next]] - startPt)
                           < 0.1*minEdgeLength
                        )
                        {
                            start = f[next];
                            firstPoint = false;
                        }
                    }
                    if (!firstPoint)
                    {
                        if
                        (
                            mag(mesh_.points()[f[next]] - endPt)
                            < 0.1*minEdgeLength
                        )
                        {
                            end = f[next];
                            break;
                        }
                    }
                    next = f.fcIndex(next);
                }
            }
            else
            {
                WarningInFunction
                    << "Could not find valid split for face "
                    << meshFaceI<< " On patch " <<pp.name()
                    << endl;
            }

            if (start != -1 && end != -1)
            {
                DynamicList<label> newFace(f.size());
                DynamicList<label> modFace(f.size());

                label next = findIndex(f, start);
                newFace.append(f[next]);
                forAll(f, fp)
                {
                    next = f.fcIndex(next);
                    newFace.append(f[next]);
                    if (f[next] == end)
                    {
                        break;
                    }
                }
                newFace.shrink();

                next = findIndex(f, end);
                modFace.append(f[next]);
                forAll(f, sFI)
                {
                    next = f.fcIndex(next);
                    modFace.append(f[next]);
                    if (f[next] == start)
                    {
                        break;
                    }
                }
                modFace.shrink();

                label own = mesh_.faceOwner()[meshFaceI];
                label nei =
                (
                    patchID == -1 ? mesh_.faceNeighbour()[meshFaceI] : -1
                );
                label zoneID = mesh_.faceZones().whichZone(meshFaceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh_.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(meshFaceI)];
                }

                meshMod.modifyFace
                (
                    face(modFace),    //vertices
                    meshFaceI,  // label of face being modified
                    own,        // owner
                    nei,        // neighbour
                    false,      // face flip
                    patchID,    // new patch for face
                    zoneID,     // zone for face
                    zoneFlip    // face flip in zone
                 );

                meshMod.addFace
                (
                    face(newFace),    // vertices
                    own,        // owner,
                    nei,        // neighbour,
                    -1,         // masterPointID,
                    -1,         // masterEdgeID,
                    meshFaceI,  // masterFaceID,
                    false,      // flipFaceFlux,
                    patchID,    // patchID,
                    zoneID,     // zoneID,
                    zoneFlip    // zoneFlip
                 );
            }
        }
    }

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    // Update fields
    mesh_.updateMesh(map);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update intersection info
    updateMesh(map, getChangedFaces(map, labelList(0)));

    forAll(cellLoops, i)
    {
        cutCells[i] = map().reverseCellMap()[cutCells[i]];
        labelList& cLoop = cellLoops[i];
        forAll(cLoop, cLI)
        {
            cLoop[cLI] = map().reversePointMap()[cLoop[cLI]];
        }
    }

    forAll(cellLoopsPyr, i)
    {
        cutCellsPyr[i] = map().reverseCellMap()[cutCellsPyr[i]];
        labelList& cLoop = cellLoopsPyr[i];
        forAll(cLoop, cLI)
        {
            cLoop[cLI] = map().reversePointMap()[cLoop[cLI]];
        }
    }

    forAll(pyrStartEnd, i)
    {
        pyrStartEnd[i][0] = map().reversePointMap()[pyrStartEnd[i][0]];
        pyrStartEnd[i][1] = map().reversePointMap()[pyrStartEnd[i][1]];
    }

    pointField tPoints(movedPoints);
    forAll(map().reversePointMap(), pointI)
    {
        tPoints[pointI] = movedPoints[map().pointMap()[pointI]];
    }
    movedPoints = tPoints;
}


Foam::boolList Foam::meshRefinement::calculateCutNodes
(
    vectorField& cutNormals
) const
{
    boolList cutNodes(mesh_.nPoints(), false);

    const labelList surfacesToBaffle(identity(surfaces_.surfaces().size()));

    label nIter = 0;
    label maxIter = 3;//5;

    vectorField dispVec(mesh_.nPoints(), vector::zero);
    vectorField edgeDispVec(mesh_.nEdges(), vector::zero);
    pointField newPoints = mesh_.points();

    labelList hitSurf(mesh_.nPoints(), -1);
    List<pointIndexHit> hitInfo(mesh_.nPoints(), pointIndexHit());

    PackedBoolList isMasterEdges = syncTools::getMasterEdges(mesh_);

    while (nIter < maxIter)
    {
        Info<<"Iteration: "<<nIter<<endl;

        cutNodes = false;
        hitSurf = -1;
        hitInfo = pointIndexHit();

        pointField start(mesh_.nEdges());
        pointField end(mesh_.nEdges());

        edgeDispVec = vector::zero;
        labelList nEdgeCuts(mesh_.nEdges(), 0);

        for (int i = 0; i < 2; i++)
        {
            for (direction j = 0; j < vector::nComponents; j++)
            {
                vector dir = vector(0., 0., 0.);
                forAll(mesh_.edges(), edgeI)
                {
                    edge e = mesh_.edges()[edgeI];
                    scalar edgeLen =  e.mag(mesh_.points());
                    if (i == 0)
                    {
                        dir[j] = 0.05*edgeLen;
                    }
                    else
                    {
                        dir[j] = -0.05*edgeLen;
                    }

                    start[edgeI] = mesh_.points()[e[0]] + dir;
                    end[edgeI] = mesh_.points()[e[1]] + dir;
                }

                // Do test for intersections
                // ~~~~~~~~~~~~~~~~~~~~~~~~~
                labelList surface1;
                List<pointIndexHit> hit1;
                labelList region1;
                labelList surface2;
                List<pointIndexHit> hit2;
                labelList region2;

                surfaces_.findNearestIntersection
                (
                    surfacesToBaffle,
                    start,
                    end,

                    surface1,
                    hit1,
                    region1,

                    surface2,
                    hit2,
                    region2
                 );


                forAll(mesh_.edges(), meshedgei)
                {
                    edge e = mesh_.edges()[meshedgei];
                    if (hit1[meshedgei].hit() && hit2[meshedgei].hit())
                    {
                        point hPoint1 = hit1[meshedgei].hitPoint();
                        point hPoint2 = hit2[meshedgei].hitPoint();

                        vector eVec = e.vec(mesh_.points());
                        scalar eLength = e.mag(mesh_.points());
                        vector unitEdgeVec = eVec / (eLength + SMALL);

                       //project cut points
                        scalar eCut1 =
                        (
                           (hPoint1 - start[meshedgei]) & unitEdgeVec
                        ) /eLength;
                        scalar eCut2 =
                        (
                           (hPoint2 - start[meshedgei]) & unitEdgeVec
                        ) /eLength;

                        if (eCut1 >= 0.5 && eCut2 >= 0.5)
                        {
                            edgeDispVec[meshedgei] +=
                            (hPoint2 - mesh_.points()[e[1]]);
                            nEdgeCuts[meshedgei]++;

                            cutNodes[e[1]] = true;

                            hitInfo[e[1]] = hit2[meshedgei];
                            hitSurf[e[1]] = surface2[meshedgei];
                        }
                        else if (eCut1 < 0.5 && eCut2 < 0.5)
                        {
                            edgeDispVec[meshedgei] +=
                            (hPoint1 - mesh_.points()[e[0]]);
                            nEdgeCuts[meshedgei]++;

                            cutNodes[e[0]] = true;

                            hitInfo[e[0]] = hit1[meshedgei];
                            hitSurf[e[0]] = surface1[meshedgei];
                        }
                        else if (eCut1 <= 0.5 && eCut2 >= 0.5)
                        {
                            cutNodes[e[0]] = true;
                            cutNodes[e[1]] = true;

                            hitInfo[e[0]] = hit1[meshedgei];
                            hitSurf[e[0]] = surface1[meshedgei];
                            hitInfo[e[1]] = hit2[meshedgei];
                            hitSurf[e[1]] = surface2[meshedgei];
                        }
                    }
                }
            }
        }

        forAll(mesh_.edges(), edgeI)
        {
            if (nEdgeCuts[edgeI] > 0)
            {
                edgeDispVec[edgeI] /= nEdgeCuts[edgeI];
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            cutNodes,
            orEqOp<bool>(),
            false        // null value
         );

        syncTools::syncEdgeList
        (
            mesh_,
            edgeDispVec,
            maxMagEqOp(),
            vector::zero        // null value
         );

        labelList nFound(mesh_.nPoints(), 0);
        nFound = 0;
        dispVec = vector::zero;

        forAll(mesh_.points(), pointI)
        {
            const labelList& pEdges = mesh_.pointEdges()[pointI];
            forAll(pEdges, pEI)
            {
                label meshedgei = pEdges[pEI];
                if
                (
                    mag(edgeDispVec[meshedgei]) > SMALL
                    && isMasterEdges.get(meshedgei) == 1
                )
                {
                    dispVec[pointI] += edgeDispVec[meshedgei];
                    nFound[pointI]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            nFound,
            plusEqOp<label>(),
            label(0)      // null value
         );

        syncTools::syncPointList
        (
            mesh_,
            dispVec,
            plusEqOp<vector>(),
            vector::zero        // null value
         );


        forAll(mesh_.points(), pointI)
        {
            if (nFound[pointI] > 1)
            {
                dispVec[pointI] /= nFound[pointI];
            }
        }

        //smooth displacement once
        {
            forAll(mesh_.edges(), edgeI)
            {
                edge e = mesh_.edges()[edgeI];

                edgeDispVec[edgeI] = 0.5*(dispVec[e[0]]+dispVec[e[1]]);
            }

            dispVec = vector::zero;
            nFound = 0;

            forAll(mesh_.points(), pointI)
            {
                const labelList& pEdges = mesh_.pointEdges()[pointI];
                forAll(pEdges, pEI)
                {
                    label meshedgei = pEdges[pEI];
                    if
                    (
                        mag(edgeDispVec[meshedgei]) > SMALL
                        && isMasterEdges.get(meshedgei) == 1
                    )
                    {
                        dispVec[pointI] += edgeDispVec[meshedgei];
                        nFound[pointI]++;
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh_,
                nFound,
                plusEqOp<label>(),
                label(0)      // null value
             );

            syncTools::syncPointList
            (
                mesh_,
                dispVec,
                plusEqOp<vector>(),
                vector::zero        // null value
             );

            forAll(mesh_.points(), pointI)
            {
                if (nFound[pointI] > 1)
                {
                    dispVec[pointI] /= nFound[pointI];
                }
            }
        }

        newPoints = mesh_.points();
        forAll(newPoints, pointI)
        {
            newPoints[pointI] = mesh_.points()[pointI] + 1.0*dispVec[pointI];
        }
        mesh_.movePoints(newPoints);

        nIter++;
    }

    //check intersection to cell vertices
    {
        dispVec = vector::zero;
        scalarField snapDist(mesh_.nPoints(), GREAT);

        const scalar edge0Len = meshCutter_.level0EdgeLength();
        const labelList cellLevel = meshCutter_.cellLevel();

        const labelListList& pointCells = mesh_.pointCells();
        forAll(pointCells, pointI)
        {
            const labelList& pCells = pointCells[pointI];

            forAll(pCells, pCellI)
            {
                label celli = pCells[pCellI];
                label  level = cellLevel[celli];
                scalar len = edge0Len / pow(2., level);

                snapDist[pointI] = min(snapDist[pointI], len);
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            snapDist,
            minEqOp<scalar>(),  // combine op
            GREAT              // null value
         );

        List<pointIndexHit> hitInfoNearest;
        labelList hitSurface;

        surfaces_.findNearest
        (
            surfacesToBaffle,
            0.5*(mesh_.points()+oldPoints()),//mesh_.points(),
            sqr(0.45*snapDist),        // sqr of attract distance
            hitSurface,
            hitInfoNearest
         );

        forAll(mesh_.points(), pointI)
        {
            if (!cutNodes[pointI] && hitInfoNearest[pointI].hit())
            {
                cutNodes[pointI] = true;
                dispVec[pointI] = hitInfoNearest[pointI].hitPoint()
                    - mesh_.points()[pointI];
                hitInfo[pointI] = hitInfoNearest[pointI];
                hitSurf[pointI] = hitSurface[pointI];
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            dispVec,
            minMagSqrEqOp<vector>(),
            vector(GREAT, GREAT, GREAT)
        );

        newPoints = mesh_.points() + dispVec;
        mesh_.movePoints(newPoints);
    }

    syncTools::syncPointList
    (
        mesh_,
        cutNodes,
        orEqOp<bool>(),
        false        // null value
     );

    forAll(surfacesToBaffle, sI)
    {
        label surfI = surfacesToBaffle[sI];
        DynamicList<pointIndexHit> localHits;
        forAll(hitSurf, pointI)
        {
            if (hitSurf[pointI] == surfI)
            {
                localHits.append(hitInfo[pointI]);
            }
        }
        localHits.shrink();

        pointField localNormals;
        label geomI = surfaces_.surfaces()[surfI];
        surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

        label localI = 0;
        forAll(hitSurf, pointI)
        {
            if (hitSurf[pointI] == surfI)
            {
                cutNormals[pointI] = localNormals[localI];
                localI++;
            }
        }
    }

    autoPtr<indirectPrimitivePatch> ppUnmeshedPtr
    (
        meshRefinement::makePatch
        (
            mesh_,
            unmeshedPatches()
         )
     );
    indirectPrimitivePatch& ppUnmeshed = ppUnmeshedPtr();
    pointField pointNormals(PatchTools::pointNormals(mesh_, ppUnmeshed));

    forAll(ppUnmeshed.meshPoints(), ptI)
    {
        label meshPointI = ppUnmeshed.meshPoints()[ptI];
        cutNodes[meshPointI] = true;
        cutNormals[meshPointI] = pointNormals[ptI];
    }

    return cutNodes;
}


Foam::boolList Foam::meshRefinement::updateCutNodes
(
    const vectorField& cutNormals,
    boolList& cutNodes,
    pointField& movedPoints
) const
{
    boolList excludedCells(mesh_.nCells(), false);

    bool hexMesh = true;
    //Exclude all non-hex cells from cutting
    forAll(mesh_.cells(), celli)
    {
        const cell& c = mesh_.cells()[celli];
        cellFeatures cellFeat(mesh_, 0.707, celli);

        if (cellFeat.faces().size() != 6)
        {
            hexMesh = false;
            labelList meshPts = c.labels(mesh_.faces());
            forAll(meshPts, ptI)
            {
                cutNodes[meshPts[ptI]] = false;
            }
        }
    }

    if (!returnReduce(hexMesh, orOp<bool>()))
    {
        Info<<"Detected non-hex base cells. "
            <<"Excluding these cells from cutting"<<endl;
    }

    forAll(mesh_.cells(), celli)
    {
        const cell& c = mesh_.cells()[celli];

        cellFeatures cellFeat(mesh_, 0.707, celli);

        if (cellFeat.faces().size() != 6)
        {
            continue;
        }

        const cell superCell(identity(6));
        labelList superPts = superCell.labels(cellFeat.faces());

        edgeList superEdges = superCell.edges(cellFeat.faces());

        label nCut = 0;
        boolList cutSuperCell(superPts.size(), false);
        Map<label> pointMap(superPts.size());
        forAll(superPts, cPtI)
        {
            label meshPointI = superPts[cPtI];
            pointMap.insert(meshPointI,cPtI);
            if (cutNodes[meshPointI])
            {
                cutSuperCell[cPtI] = true;
                nCut++;
            }
        }

        label nConEdges = 0;
        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (cutNodes[e[0]] && cutNodes[e[1]])
            {
                nConEdges++;
            }
        }

        labelList nConn(superPts.size(), 0);
        forAll(superEdges, cEI)
        {
            edge e = superEdges[cEI];
            if (cutNodes[e[0]] && cutNodes[e[1]])
            {
                nConn[pointMap[e[0]]]++;
                nConn[pointMap[e[1]]]++;
            }
        }

        /*
        labelList cellMeshPts = c.labels(mesh_.faces());
        label nTotalCuts = 0;
        forAll(cellMeshPts, cPtI)
        {
            if (cutNodes[cellMeshPts[cPtI]])
            {
                nTotalCuts++;
            }
        }
        */

        bool refinedCell = true;
        if (c.size() == 6 && c.labels(mesh_.faces()).size() == 8)
        {
            refinedCell = false;
        }

        if (nCut == 4 && nConEdges == 3)
        {
            forAll(superPts, cPtI)
            {
                label meshPointI = superPts[cPtI];
                if (nConn[cPtI] == 3)
                {
                    cutNodes[meshPointI] = false;
                    break;
                }
            }
        }
        else if (nCut == 4 && nConEdges == 2)
        {
            label basePt1 = -1;
            label basePt2 = -1;
            forAll(superPts, cPtI)
            {
                label meshPointI = superPts[cPtI];
                if (nConn[cPtI] == 0 && cutSuperCell[cPtI])
                {
                    basePt1 = meshPointI;
                }
                if (nConn[cPtI] == 2)
                {
                    basePt2 = meshPointI;
                }
            }

            if (basePt1 != -1 && basePt2 != -1)
            {
                if
                (
                    cutNormals[basePt1] != vector(GREAT, GREAT, GREAT)
                    && cutNormals[basePt2] != vector(GREAT, GREAT, GREAT)
                 )
                {
                    vector cNorm1 =
                        cutNormals[basePt1]/(mag(cutNormals[basePt1]) + SMALL);
                    vector cNorm2 =
                        cutNormals[basePt2]/(mag(cutNormals[basePt2]) + SMALL);

                    if (mag((cNorm1&cNorm2)) < 0.1)
                    {
                        forAll(superPts, cPtI)
                        {
                            label meshPointI = superPts[cPtI];
                            if (!cutSuperCell[cPtI])
                            {
                                DynamicList<label> nbrNodes;
                                forAll(superEdges, cEI)
                                {
                                    edge e = superEdges[cEI];
                                    if (e[0] == meshPointI || e[1] == meshPointI)
                                    {
                                        label otherPt =
                                            (e[0] == meshPointI ? e[1] : e[0]);
                                        if (cutNodes[otherPt])
                                        {
                                            nbrNodes.append(otherPt);
                                        }
                                    }
                                }
                                nbrNodes.shrink();

                                if (nbrNodes.size() == 2)
                                {
                                    if
                                    (
                                        cutNormals[nbrNodes[0]] !=
                                        vector(GREAT, GREAT, GREAT)
                                        && cutNormals[nbrNodes[1]] !=
                                        vector(GREAT, GREAT, GREAT)
                                     )
                                    {
                                        vector cNbrNorm1 =
                                            cutNormals[nbrNodes[0]]
                                            /(mag(cutNormals[nbrNodes[0]])
                                              + SMALL);
                                        vector cNbrNorm2 =
                                            cutNormals[nbrNodes[1]]
                                            /(mag(cutNormals[nbrNodes[1]])
                                              + SMALL);
                                        if (mag(cNbrNorm1 & cNbrNorm2) < 0.1)
                                        {
                                            cutNodes[meshPointI] = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                DynamicList<vector> edgeNorm;
                vector eDir0 = vector(GREAT, GREAT, GREAT);
                bool firstEdge = true;
                forAll(superEdges, cEI)
                {
                    edge e = superEdges[cEI];
                    if (cutNodes[e[0]] && cutNodes[e[1]])
                    {
                        if
                        (
                            cutNormals[e[0]] != vector(GREAT, GREAT, GREAT)
                            && cutNormals[e[1]] != vector(GREAT, GREAT, GREAT)
                         )
                        {
                            if (firstEdge)
                            {
                                eDir0 = mesh_.points()[e[1]]
                                    - mesh_.points()[e[0]];

                                firstEdge = false;

                                vector norm = cutNormals[e[0]];
                                norm /= (mag(norm) + SMALL);
                                edgeNorm.append(norm);
                                norm = cutNormals[e[1]];
                                norm /= (mag(norm) + SMALL);
                                edgeNorm.append(norm);
                            }
                            else
                            {
                                vector eDir = mesh_.points()[e[1]]
                                    - mesh_.points()[e[0]];
                                label pt0 = e[0];
                                label pt1 = e[1];
                                if ((eDir&eDir0)<0)
                                {
                                    pt0 = e[1];
                                    pt1 = e[0];
                                }

                                vector norm = cutNormals[pt0];
                                norm /= (mag(norm) + SMALL);
                                edgeNorm.append(norm);
                                norm = cutNormals[pt1];
                                norm /= (mag(norm) + SMALL);
                                edgeNorm.append(norm);
                            }
                        }
                    }
                }

                edgeNorm.shrink();
                if (edgeNorm.size() == 4)
                {
                    if
                    (
                        mag(edgeNorm[0]&edgeNorm[2]) < 0.342
                        &&mag(edgeNorm[1]&edgeNorm[3]) < 0.342
                     )
                    {
                        excludedCells[celli] = true;
                    }
                }
            }
        }
        else if (nCut == 5 && nConEdges == 4)
        {
            DynamicList<vector> edgeNorm;
            vector eDir0 = vector(GREAT, GREAT, GREAT);
            bool firstEdge = true;
            forAll(superEdges, cEI)
            {
                edge e = superEdges[cEI];

                if
                (
                    cutNodes[e[0]] && cutNodes[e[1]]
                    &&  ((nConn[pointMap[e[0]]] == 1)
                        || (nConn[pointMap[e[1]]] == 1))
                 )
                {
                    if
                    (
                        cutNormals[e[0]] != vector(GREAT, GREAT, GREAT)
                        && cutNormals[e[1]] != vector(GREAT, GREAT, GREAT)
                     )
                    {
                        if (firstEdge)
                        {
                            eDir0 = mesh_.points()[e[1]]
                                - mesh_.points()[e[0]];

                            firstEdge = false;

                            vector norm = cutNormals[e[0]];
                            norm /= (mag(norm) + SMALL);
                            edgeNorm.append(norm);
                            norm = cutNormals[e[1]];
                            norm /= (mag(norm) + SMALL);
                            edgeNorm.append(norm);
                        }
                        else
                        {
                            vector eDir = mesh_.points()[e[1]]
                                - mesh_.points()[e[0]];
                            label pt0 = e[0];
                            label pt1 = e[1];
                            if ((eDir&eDir0)<0)
                            {
                                pt0 = e[1];
                                pt1 = e[0];
                            }

                            vector norm = cutNormals[pt0];
                            norm /= (mag(norm) + SMALL);
                            edgeNorm.append(norm);
                            norm = cutNormals[pt1];
                            norm /= (mag(norm) + SMALL);
                            edgeNorm.append(norm);
                        }
                    }
                }
            }

            edgeNorm.shrink();
            if (edgeNorm.size() == 4)
            {
                if
                (
                    mag(edgeNorm[0]&edgeNorm[2]) < 0.342
                    &&mag(edgeNorm[1]&edgeNorm[3]) < 0.342
                 )
                {
                    labelList nConnNbr(superPts.size(), 0);
                    forAll(superEdges, cEI)
                    {
                        edge e = superEdges[cEI];

                        if (cutNodes[e[0]])
                        {
                            nConnNbr[pointMap[e[1]]]++;
                        }

                        if (cutNodes[e[1]])
                        {
                            nConnNbr[pointMap[e[0]]]++;
                        }
                    }

                    forAll(nConnNbr, ptI)
                    {
                        if (nConnNbr[ptI] == 3)
                        {
                            cutNodes[superPts[ptI]] = true;
                        }
                    }
                    excludedCells[celli] = true;
                }
            }
        }
        else if (refinedCell && nCut == 3 && nConEdges == 0)
        {
            DynamicList<label> cNodes(3);
            forAll(superPts, i)
            {
                if (cutNodes[superPts[i]])
                {
                    cNodes.append(superPts[i]);
                }
            }

            cNodes.shrink();
            edgeList cEdges(3);
            cEdges[0] = edge(cNodes[0],cNodes[1]);
            cEdges[1] = edge(cNodes[0],cNodes[2]);
            cEdges[2] = edge(cNodes[1],cNodes[2]);

            bool valid = true;
            forAll(cEdges, i)
            {
                forAll(cellFeat.faces(), faceI)
                {
                    face f1 = cellFeat.faces()[faceI];
                    label i1 = findIndex(f1, cEdges[i][0]);
                    label i2 = findIndex(f1, cEdges[i][1]);

                    if (i1 != -1 && i2 != -1)
                    {
                        const labelList& sideFaces = cellFeat.faceMap()[faceI];
                        if (sideFaces.size() > 1)
                        {
                            forAll(sideFaces, sFI)
                            {
                                face sf1 = mesh_.faces()[sideFaces[sFI]];
                                label isf1 = findIndex(sf1, cEdges[i][0]);
                                if (isf1 != -1)
                                {
                                    label next = sf1.fcIndex(isf1);
                                    next = sf1.fcIndex(next);
                                    if (!cutNodes[sf1[next]])
                                    {
                                        valid = false;
                                    }
                                }
                                if (!valid)
                                {
                                    break;
                                }
                            }
                        }
                    }
                    if (!valid)
                    {
                        break;
                    }
                }
                if (!valid)
                {
                    break;
                }
            }
            if (!valid)
            {
                excludedCells[celli] = true;
            }
        }
        else if (refinedCell && nCut == 3 && nConEdges == 1)
        {
            bool prism = false;

            boolList newCuts = cutSuperCell;
            forAll(superPts, cPtI)
            {
                if (!cutSuperCell[cPtI])
                {
                    newCuts =  cutSuperCell;
                    newCuts[cPtI] = true;

                    label nConEdgesNew = 0;
                    forAll(superEdges, cEI)
                    {
                        edge e = superEdges[cEI];
                        if
                        (
                            newCuts[pointMap[e[0]]]
                            && newCuts[pointMap[e[1]]]
                         )
                        {
                            nConEdgesNew++;
                        }
                    }

                    if (nConEdgesNew == 2)
                    {
                        labelList nConnNew(superPts.size(), 0);
                        forAll(superEdges, cEI)
                        {
                            edge e = superEdges[cEI];
                            if
                            (
                                newCuts[pointMap[e[0]]]
                                && newCuts[pointMap[e[1]]]
                             )
                            {
                                nConnNew[pointMap[e[0]]]++;
                                nConnNew[pointMap[e[1]]]++;
                            }
                        }
                        prism = true;
                        forAll(superPts, cPtI)
                        {
                            if (nConnNew[cPtI] == 2)
                            {
                                prism = false;
                                break;
                            }
                        }
                        if (prism)
                        {
                            break;
                        }
                    }
                    else
                    {
                        continue;
                    }
                }
            }

            if (prism)
            {
                forAll(newCuts, ptI)
                {
                    if (newCuts[ptI])
                    {
                        label meshPointI = superPts[ptI];
                        cutNodes[meshPointI] = true;
                    }
                }
            }
        }
    }

    //Set refined cell non master points to cut on master cut edges
    {
        const labelList& pointLevel = meshCutter_.pointLevel();
        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(mesh_.cells(), celli)
        {
            const cell& c = mesh_.cells()[celli];
            label level = cellLevel[celli];

            bool refinedCell = true;
            if (c.size() == 6 && c.labels(mesh_.faces()).size() == 8)
            {
                refinedCell = false;
            }

            if (refinedCell)
            {
                edgeList cEdges = c.edges(mesh_.faces());
                List<label> cellMeshEdges(cEdges.size());

                forAll(cEdges, cEI)
                {
                    edge e = cEdges[cEI];
                    cellMeshEdges[cEI] = meshTools::findEdge
                    (
                        mesh_.edges(),
                        mesh_.pointEdges()[e[0]],
                        e[0],
                        e[1]
                     );
                }

                labelHashSet cellEdges(cellMeshEdges);
                labelList mPoints = c.labels(mesh_.faces());
                forAll(mPoints, mptI)
                {
                    label meshPointI = mPoints[mptI];

                    if (pointLevel[meshPointI] > level)
                    {
                        labelList pEdges = mesh_.pointEdges()[meshPointI];
                        bool allCuts = true;
                        forAll(pEdges, pEI)
                        {
                            label meshedgei = pEdges[pEI];
                            if (cellEdges.found(meshedgei))
                            {
                                edge e = mesh_.edges()[meshedgei];

                                label otherPt(e[0] ==  meshPointI ? e[1] : e[0]);
                                if (!cutNodes[otherPt])
                                {
                                    allCuts = false;
                                }
                            }
                        }
                        if (allCuts)
                        {
                            cutNodes[meshPointI] = true;
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        cutNodes,
        orEqOp<bool>(),
        false        // null value
     );

    //check refined faces where neighbours might produce inconsistent cuts
    forAll(mesh_.cells(), celli)
    {
        const cell& c = mesh_.cells()[celli];

        bool refinedCell = true;
        if (c.size() == 6 && c.labels(mesh_.faces()).size() == 8)
        {
            refinedCell = false;
        }

        if (refinedCell)
        {
            cellFeatures cellFeat(mesh_, 0.707, celli);
            forAll(cellFeat.faces(), faceI)
            {
                const face f1 = cellFeat.faces()[faceI];

                labelHashSet masterCuts(f1.size());
                label nMasterCuts = 0;
                forAll(f1, fp)
                {
                    if (cutNodes[f1[fp]])
                    {
                        masterCuts.insert(f1[fp]);
                        nMasterCuts++;
                    }
                }

                const labelList& slaveFaces = cellFeat.faceMap()[faceI];

                bool diagCut = true;
                for (const edge e : f1.walkEdges())
                {
                    label pt0 = e[0];
                    label pt1 = e[1];
                    if (cutNodes[pt0] && cutNodes[pt1])
                    {
                        diagCut = false;
                    }
                }

                if
                (
                    (nMasterCuts > 2 ||  (nMasterCuts == 2 && diagCut))
                    &&  slaveFaces.size() >1
                )
                {
                    forAll(slaveFaces, sFI)
                    {
                        face f2 = mesh_.faces()[slaveFaces[sFI]];
                        labelHashSet slaveCuts(f1.size());
                        label nSlaveCuts = 0;
                        label nSlaveMasterCuts = 0;

                        forAll(f2, fp)
                        {
                            if
                            (
                                cutNodes[f2[fp]]
                                && !masterCuts.found(f2[fp])
                                && !slaveCuts.found(f2[fp])
                            )
                            {
                                slaveCuts.insert(f2[fp]);
                                nSlaveCuts++;
                            }
                            else if
                            (
                                cutNodes[f2[fp]]
                                && masterCuts.found(f2[fp])
                            )
                            {
                                nSlaveMasterCuts++;
                            }
                        }

                        bool slaveDiagCut = true;

                        for (const edge e : f2.walkEdges())
                        {
                            label pt0 = e[0];
                            label pt1 = e[1];
                            if
                            (
                                slaveCuts.found(pt0)
                                && slaveCuts.found(pt1)
                            )
                            {
                                slaveDiagCut = false;
                            }
                        }

                        if
                        (
                            (
                                nSlaveCuts > 2
                                || (nSlaveCuts == 2 && slaveDiagCut)
                            )
                            && nSlaveMasterCuts > 0
                        )
                        {
                            excludedCells[celli] = true;
                            break;
                        }
                    }
                }
            }
        }
    }

    const pointField& oldPts = oldPoints();
    forAll(cutNodes, pointI)
    {
        if (!cutNodes[pointI])
        {
            movedPoints[pointI] = oldPts[pointI];
        }
    }

    return excludedCells;
}


void Foam::meshRefinement::resetNonManifoldSplits
(
    const labelList& newFaceOwner,
    const labelList& newFaceNeighbour,
    List<labelList>& cellLoopsPyr,
    List<label>& cutCellsPyr,
    List<labelList>& pyrStartEnd,
    pointField& movedPoints
)
{
    polyTopoChange meshMod(mesh_);

    forAll(newFaceOwner, faceI)
    {
        if (newFaceOwner[faceI] != -1 || newFaceNeighbour[faceI] != -1)
        {
            if (mesh_.isInternalFace(faceI))
            {
                label own =
                (
                    newFaceOwner[faceI] != -1
                    ? newFaceOwner[faceI]
                    : mesh_.faceOwner()[faceI]
                 );

                label nei =
                (
                    newFaceNeighbour[faceI] != -1
                    ? newFaceNeighbour[faceI]
                    : mesh_.faceNeighbour()[faceI]
                 );

                if (own > nei)
                {
                    label town = own;
                    own = nei;
                    nei = town;
                }

                vector fNorm = mesh_.faceAreas()[faceI];
                vector ccVec = mesh_.cellCentres()[own]
                    - mesh_.cellCentres()[nei];

                bool reverseFace = false;
                if ((ccVec & fNorm) < 0)
                {
                    reverseFace = true;
                }

                label zoneID = mesh_.faceZones().whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh_.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }


                face modFace = mesh_.faces()[faceI];
                if (reverseFace)
                {
                    modFace = modFace.reverseFace();
                }
                meshMod.modifyFace
                (
                    modFace,
                    faceI,  // label of face being modified
                    own,        // owner
                    nei,        // neighbour
                    false,      // face flip
                    -1,    // new patch for face
                    zoneID,     // zone for face
                    zoneFlip    // face flip in zone
                 );
            }
            else
            {
                label patchID = mesh_.boundaryMesh().whichPatch(faceI);
                label zoneID = mesh_.faceZones().whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh_.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }

                meshMod.modifyFace
                (
                    mesh_.faces()[faceI],    //vertices
                    faceI,  // label of face being modified
                    newFaceOwner[faceI],        // owner
                    -1,        // neighbour
                    false,      // face flip
                    patchID,    // new patch for face
                    zoneID,     // zone for face
                    zoneFlip    // face flip in zone
                 );
            }
        }
    }

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    // Update fields
    mesh_.updateMesh(map);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update intersection info
    updateMesh(map, getChangedFaces(map, labelList(0)));

    pointField tPoints(movedPoints);
    forAll(map().reversePointMap(), pointI)
    {
        tPoints[pointI] = movedPoints[map().pointMap()[pointI]];
    }
    movedPoints = tPoints;

    forAll(cellLoopsPyr, i)
    {
        cutCellsPyr[i] = map().reverseCellMap()[cutCellsPyr[i]];
        labelList& cLoop = cellLoopsPyr[i];
        forAll(cLoop, cLI)
        {
            cLoop[cLI] = map().reversePointMap()[cLoop[cLI]];
        }
    }
    forAll(pyrStartEnd, i)
    {
        pyrStartEnd[i][0] = map().reversePointMap()[pyrStartEnd[i][0]];
        pyrStartEnd[i][1] = map().reversePointMap()[pyrStartEnd[i][1]];
    }
}


//Perform cutting of cells
void Foam::meshRefinement::splitCells
(
    List<label>& cutCells,
    List<label>& cutCellsPyr,
    List<labelList>& cellLoops,
    List<labelList>& cellLoopsPyr,
    List<labelList>& pyrStartEnd,
    pointField& movedPoints
)
{
    List<scalarField> cellWeights(cutCells.size());
    forAll(cellWeights, cwI)
    {
        cellWeights[cwI].setSize(cellLoops[cwI].size(), -GREAT);
    }

    List<scalarField> cellWeightsPyr(cutCellsPyr.size());
    forAll(cellWeightsPyr, cwI)
    {
        cellWeightsPyr[cwI].setSize(cellLoopsPyr[cwI].size(), -GREAT);
    }

    cellCuts cuts(mesh_, cutCells, cellLoops, cellWeights);

    Foam::meshCutter mCut(mesh_);

    polyTopoChange meshMod(mesh_);
    mCut.setRefinement(cuts, meshMod);

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    // Update fields
    mesh_.updateMesh(map);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update intersection info
    updateMesh(map, getChangedFaces(map, labelList(0)));

    pointField tPoints(movedPoints);
    forAll(map().reversePointMap(), pointI)
    {
       tPoints[pointI] = movedPoints[map().pointMap()[pointI]];
    }
    movedPoints = tPoints;

    forAll(cellLoopsPyr, i)
    {
        cutCellsPyr[i] = map().reverseCellMap()[cutCellsPyr[i]];
        labelList& cLoop = cellLoopsPyr[i];
        forAll(cLoop, cLI)
        {
            cLoop[cLI] = map().reversePointMap()[cLoop[cLI]];
        }
    }
    forAll(pyrStartEnd, i)
    {
        pyrStartEnd[i][0] = map().reversePointMap()[pyrStartEnd[i][0]];
        pyrStartEnd[i][1] = map().reversePointMap()[pyrStartEnd[i][1]];
    }

    //meshCutter does not correctly handle non-manifold splits
    //which are currently generated before splitting into pyramids
    DynamicList<labelPair> updatedFaces(mesh_.nFaces()/100 + 100);
    forAll(mesh_.cells(), celli)
    {
        if (celli < map().nOldCells())
        {
            const cell& c  = mesh_.cells()[celli];
            labelHashSet cfSet(c);
            labelHashSet removedFaces(c.size());
            const edgeList& cEdges = c.edges(mesh_.faces());

            while (true)
            {
                label nSet = 0;
                forAll(cEdges, cEI)
                {
                    edge e = cEdges[cEI];
                    label meshedgei = meshTools::findEdge
                    (
                        mesh_.edges(),
                        mesh_.pointEdges()[e[0]],
                        e[0],
                        e[1]
                     );

                    labelList eFaces = mesh_.edgeFaces()[meshedgei];
                    label nFound = 0;
                    label nonManifoldFace = -1;

                    forAll(eFaces, eFI)
                    {
                        if
                        (
                            cfSet.found(eFaces[eFI])
                           && !removedFaces.found(eFaces[eFI])
                        )
                        {
                            nonManifoldFace = eFaces[eFI];
                            nFound++;
                        }
                    }

                    if (nFound == 1)
                    {
                        label patchID =
                            mesh_.boundaryMesh().whichPatch(nonManifoldFace);
                        if
                        (
                            patchID != -1
                            && isA<processorPolyPatch>
                            (mesh_.boundaryMesh()[patchID])
                         )
                        {
                            removedFaces.insert(nonManifoldFace);
                            updatedFaces.append
                            (
                                labelPair(nonManifoldFace, celli)
                            );
                            nSet++;
                        }
                    }
                }
                if (nSet == 0)
                {
                    break;
                }
            }
        }
    }
    updatedFaces.shrink();

    if (returnReduce(updatedFaces.size(), sumOp<label>()) > 0)
    {
        Info<<"Found incorrectly oriented cells - re-orienting"<<endl;

        labelList newFaceOwner(mesh_.nFaces(), -1);
        labelList newFaceNeighbour(mesh_.nFaces(), -1);

        Map<label> newCells(map().nOldCells());
        forAll(map().cellMap(), celli)
        {
            if (celli >= map().nOldCells())
            {
                label originalCell = map().cellMap()[celli];
                newCells.insert(originalCell,celli);
            }
        }

        forAll(updatedFaces, i)
        {
            label faceI = updatedFaces[i].first();
            label own = mesh_.faceOwner()[faceI];

            if (own == updatedFaces[i].second())
            {
                newFaceOwner[faceI] = newCells[own];
            }
            else
            {
                label nei = mesh_.faceNeighbour()[faceI];
                newFaceNeighbour[faceI] = newCells[nei];
            }
        }

        resetNonManifoldSplits
        (
            newFaceOwner,
            newFaceNeighbour,
            cellLoopsPyr,
            cutCellsPyr,
            pyrStartEnd,
            movedPoints
         );
    }

    splitPyramids
    (
        pyrStartEnd,
        cellWeightsPyr,
        cellLoopsPyr,
        cutCellsPyr,
        movedPoints
     );
}


void Foam::meshRefinement::resetCutDisplacements
(
    const List<labelList>& cellLoops,
    pointField& movedPoints
)
{
   boolList resetDisplacement(mesh_.nPoints(), true);

    forAll(cellLoops, i)
    {
        const labelList& cLoop = cellLoops[i];
        forAll(cLoop, i)
        {
            resetDisplacement[cLoop[i]] = false;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        resetDisplacement,
        andEqOp<bool>(),
        false        // null value
     );

    pointField& oPoints = oldPoints();
    label nCut = 0;
    forAll(resetDisplacement, pointI)
    {
        if (!resetDisplacement[pointI])
        {
            nCut++;
        }
    }

    pointField checkPoints(nCut);

    nCut = 0;
    forAll(resetDisplacement, pointI)
    {
        if (!resetDisplacement[pointI])
        {
            checkPoints[nCut] = movedPoints[pointI];
            nCut++;
        }
    }
    const labelList surfacesToBaffle(identity(surfaces_.surfaces().size()));

    labelList hitSurface;
    List<pointIndexHit> hitInfo;

    surfaces_.findNearest
    (
        surfacesToBaffle,
        checkPoints,
        scalarField(nCut, sqr(GREAT)),        // sqr of attract distance
        hitSurface,
        hitInfo
    );

    nCut = 0;
    forAll(mesh_.points(), pointI)
    {
        if (resetDisplacement[pointI])
        {
            movedPoints[pointI] = oPoints[pointI];
        }
        else
        {
            if (hitInfo[nCut].hit())
            {
                movedPoints[pointI] = hitInfo[nCut].hitPoint();
            }
            nCut++;
        }
    }
}


void Foam::meshRefinement::snapAndCut()
{
    //define which nodes to cut
    vectorField cutNormals(mesh_.nPoints(), vector(GREAT, GREAT, GREAT));

    Info<<"Calculating Cut Nodes"<<endl;

    boolList cutNodes = calculateCutNodes(cutNormals);

    //Move points back to original positions for topological checks.
    //Store moved positions for later intersections checks
    pointField movedPoints = mesh_.points();
    mesh_.movePoints(oldPoints());

    //try and improve the cut nodes
    boolList excludedCells = updateCutNodes(cutNormals,cutNodes,movedPoints);

    boolList excludeTetSplit = setRefinementCells();

    DynamicList<label> cutCells;
    DynamicList<labelList> cellLoops;

    DynamicList<label> cutCellsPyr;
    DynamicList<labelList> cellLoopsPyr;
    DynamicList<labelList> pyrStartEnd;

    Info<<"Selecting new cut cell types"<<endl;

    //select how to split the cells
    selectCutType
    (
        excludeTetSplit,
        excludedCells,
        cutNodes,
        cutNormals,
        cutCells,
        cutCellsPyr,
        cellLoops,
        cellLoopsPyr,
        pyrStartEnd
     );

    filterCutCells
    (
        cutCells,
        cutCellsPyr,
        cellLoops,
        cellLoopsPyr,
        pyrStartEnd
    );

    resetCutDisplacements(cellLoops, movedPoints);

    //split processor faces
    splitProcessorFaces
    (
        cutCells,
        cutCellsPyr,
        cellLoops,
        cellLoopsPyr,
        pyrStartEnd,
        movedPoints
     );

    Info<<"Performing cell cuts"<<endl;

    //split the cells
    splitCells
    (
        cutCells,
        cutCellsPyr,
        cellLoops,
        cellLoopsPyr,
        pyrStartEnd,
        movedPoints
     );

    Info<<"Move mesh to cut positions"<<endl;

    //move points to cut positions
    mesh_.movePoints(movedPoints);

    updateIntersectionsByBlock(identity(mesh_.nFaces()),20);
}
