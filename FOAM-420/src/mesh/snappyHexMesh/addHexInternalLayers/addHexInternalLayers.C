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
    (c) 2010 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "addHexInternalLayers/addHexInternalLayers.H"
#include "meshTools/meshTools.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "meshCut/cellCuts/cellCuts.H"
#include "meshCut/meshModifiers/meshCutter/meshCutter.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(addHexInternalLayers, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::addHexInternalLayers::addHexInternalLayers
(
    meshRefinement& meshRefiner,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
:
    meshRefiner_(meshRefiner),
    decomposer_(decomposer),
    distributor_(distributor),
    baffles_(createBaffles())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


Foam::List<Foam::labelPair> Foam::addHexInternalLayers::createBaffles()
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeighbour = mesh.faceNeighbour();

    //First we need to re-distribute to ensure baffles are not
    // processor boundaries
    if (Pstream::parRun())
    {
        boolList blockedFaces(mesh.nFaces(), true);
        label nBlocked = 0;

        const faceZoneMesh& fZones = mesh.faceZones();
        forAll(fZones, fZI)
        {
            const faceZone& fZone = fZones[fZI];
            forAll(fZone, i)
            {
                label faceI = fZone[i];
                blockedFaces[faceI] = false;
                nBlocked++;
            }
        }

        // Calculate neighbour cell level
        scalarField neiLevel(mesh.nFaces()-mesh.nInternalFaces());
        for
        (
            label faceI = mesh.nInternalFaces();
            faceI < mesh.nFaces();
            faceI++
        )
        {
            neiLevel[faceI-mesh.nInternalFaces()] =
                cellLevel[mesh.faceOwner()[faceI]];
        }
        syncTools::swapBoundaryFaceList(mesh, neiLevel);

        forAll(mesh.faces(), faceI)
        {
            label ownLevel = cellLevel[faceOwner[faceI]];
            if (mesh.isInternalFace(faceI))
            {
                label neiLevel = cellLevel[faceNeighbour[faceI]];

                if (ownLevel != neiLevel)
                {
                    blockedFaces[faceI] = false;
                    nBlocked++;
                }
            }
            else
            {
                label patchI = patches.whichPatch(faceI);
                if (patches[patchI].coupled())
                {
                    if
                    (
                        ownLevel != neiLevel[faceI-mesh.nInternalFaces()]
                    )
                    {
                        blockedFaces[faceI] = false;
                        nBlocked++;
                    }
                }
            }
        }

        if (returnReduce(nBlocked, sumOp<label>()) > 0)
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
        }
    }

    labelList ownPatch(mesh.nFaces(), -1);
    labelList nbrPatch(mesh.nFaces(), -1);

    // Add patch
    dictionary patchInfo;
    patchInfo.set("type", wallPolyPatch::typeName);
    label lowerID = meshRefiner_.addPatch(mesh, "lowerLevel", patchInfo);
    label higherID = meshRefiner_.addPatch(mesh, "higherLevel", patchInfo);

    label nBaffles = 0;
    forAll(mesh.faces(), faceI)
    {
        if (mesh.isInternalFace(faceI))
        {
            label ownLevel = cellLevel[faceOwner[faceI]];
            label neiLevel = cellLevel[faceNeighbour[faceI]];

            if (ownLevel != neiLevel)
            {
                ownPatch[faceI] = (ownLevel > neiLevel ? higherID : lowerID);
                nbrPatch[faceI] = (ownLevel > neiLevel ? lowerID : higherID);
                nBaffles++;
            }
        }
    }

    List<labelPair> baffles(nBaffles);
    if (returnReduce(nBaffles, sumOp<label>()) == 0)
    {
        return baffles;
    }

    // Create baffles. both sides same patch.
    {
        autoPtr<mapPolyMesh> map =
            meshRefiner_.createBaffles(ownPatch, nbrPatch, false);
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        nBaffles = 0;
        forAll(faceMap, faceI)
        {
            label oldFaceI = faceMap[faceI];

            if (ownPatch[oldFaceI] != -1)
            {
                label masterFaceI = reverseFaceMap[oldFaceI];
                if (faceI != masterFaceI)
                {
                    baffles[nBaffles++] = labelPair(masterFaceI, faceI);
                }
            }
        }
    }

    // Duplicate points
    {
        autoPtr<mapPolyMesh> map =
            meshRefiner_.dupNonManifoldPoints();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        forAll(baffles, i)
        {
            labelPair& baffle = baffles[i];
            baffle.first() = reverseFaceMap[baffle.first()];
            baffle.second() = reverseFaceMap[baffle.second()];
        }
    }

    return baffles;
}


void Foam::addHexInternalLayers::addCells
(
    const indirectPrimitivePatch& pp,
    const labelList& pointType,
    const labelList& edgeType,
    const labelList& faceType,
    List<labelPair>& newFaceCells,
    List<labelPair>& newEdgeCells,
    List<labelPair>& newPointCells,
    polyTopoChange& meshMod
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    label nFaceCells = 0;
    label nEdgeCells = 0;
    label nPointCells = 0;

    //Add face cells
    forAll(pp, patchFaceI)
    {
        label faceI = pp.addressing()[patchFaceI];

        const face f =  mesh.faces()[faceI];

        bool addCell = true;
        forAll(f,fp)
        {
            label meshPointI = f[fp];
            label nextPt = f[f.fcIndex(fp)];

            label edgeI = meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[meshPointI],
                meshPointI,
                nextPt
            );

            if (edgeType[edgeI] == 1 || pointType[meshPointI] == 5)
            {
                addCell = false;
                break;
            }
        }

        if (addCell)
        {
            label own = mesh.faceOwner()[faceI];
            label ownZoneI = mesh.cellZones().whichZone(own);

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
            newFaceCells[faceI] = labelPair(newCellI,own);
            nFaceCells++;
        }
    }

    //Add concave edge cells (no parallel treatment needed)
    forAll(mesh.edges(), edgeI)
    {
        if (edgeType[edgeI] == 1)
        {
            const edge e = mesh.edges()[edgeI];
            if
            (
                pointType[e[0]] == 2
                || pointType[e[1]] == 2
            )
            {
                continue;
            }

            label own = mesh.edgeCells()[edgeI][0];
            label ownZoneI = mesh.cellZones().whichZone(own);
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

            newEdgeCells[edgeI] = labelPair(newCellI, own);
            nEdgeCells++;
        }
    }

    //Add convex edge cells (no parallel treatment needed)
    forAll(mesh.edges(), edgeI)
    {
        if (edgeType[edgeI] == 2)
        {
            const edge e = mesh.edges()[edgeI];

            if
            (
                (pointType[e[0]] >= 5 && pointType[e[0]] <= 7)
                || (pointType[e[1]] >= 5 && pointType[e[1]] <= 7)
            )
            {
                continue;
            }

            const labelList& eCells = mesh.edgeCells()[edgeI];
            const labelList& eFaces = mesh.edgeFaces()[edgeI];

            forAll(eCells, eCI)
            {
                label own = eCells[eCI];
                labelHashSet cFaces(mesh.cells()[own]);
                bool validCell = true;

                forAll(eFaces, eFI)
                {
                    label faceI = eFaces[eFI];

                    if (cFaces.found(faceI) && faceType[faceI] != -1)
                    {
                        validCell = false;
                        break;
                    }
                }

                if (validCell)
                {
                    label ownZoneI = mesh.cellZones().whichZone(own);
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

                    newEdgeCells[edgeI] = labelPair(newCellI, own);
                    nEdgeCells++;
                    break;
                }
            }
        }
    }

    //Add point cells
    forAll(mesh.points(), pointI)
    {
        if (pointType[pointI] == 2)
        {
            label own = mesh.pointCells()[pointI][0];
            label ownZoneI = mesh.cellZones().whichZone(own);
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

            newPointCells[pointI] = labelPair(newCellI, own);
            nPointCells++;
        }
        else if
        (
            pointType[pointI] >= 5 && pointType[pointI] <= 7
        )
        {
            label own = -1;
            const labelList& pCells = mesh.pointCells()[pointI];
            const labelList& pEdges = mesh.pointEdges()[pointI];
            forAll(pCells, pCI)
            {
                label cellI = pCells[pCI];
                labelHashSet cEdges(mesh.cellEdges()[cellI]);

                label nConcave = 0;
                label nConvex = 0;
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cEdges.found(edgeI))
                    {
                        if (edgeType[edgeI] == 1)
                        {
                            nConcave++;
                        }
                        else if (edgeType[edgeI] == 2)
                        {
                            nConvex++;
                        }
                    }
                }
                if
                (
                    (pointType[pointI] == 5 && nConvex == 1 && nConcave == 0)
                    || (pointType[pointI] == 6 && nConvex == 2 && nConcave == 0)
                    || (pointType[pointI] == 7 && nConvex == 3 && nConcave == 0)
                )
                {
                    own = cellI;
                    break;
                }
            }

            if (own != -1)
            {
                label ownZoneI = mesh.cellZones().whichZone(own);
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

                newPointCells[pointI] = labelPair(newCellI, own);
                nPointCells++;
            }
        }
        else if (pointType[pointI] == 4)
        {
            const labelList& pCells = mesh.pointCells()[pointI];
            const labelList& pEdges = mesh.pointEdges()[pointI];

            forAll(pCells, pCI)
            {
                label own = pCells[pCI];
                labelHashSet cEdges(mesh.cellEdges()[own]);

                bool validCell = true;
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cEdges.found(edgeI) && edgeType[edgeI] != -1)
                    {
                        validCell = false;
                        break;
                    }
                }

                if (validCell)
                {
                    label ownZoneI = mesh.cellZones().whichZone(own);
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

                    newPointCells[pointI] = labelPair(newCellI, own);
                    nPointCells++;
                    break;
                }
            }
        }
    }

    Info<<"Number of face cells created: "
        <<returnReduce(nFaceCells, sumOp<label>())<<endl;
    Info<<"Number of edge cells created: "
        <<returnReduce(nEdgeCells, sumOp<label>())<<endl;
    Info<<"Number of point cells created: "
        <<returnReduce(nPointCells, sumOp<label>())<<endl;
}


void Foam::addHexInternalLayers::calculateTypes
(
    const indirectPrimitivePatch& pp,
    const vectorField& meshEdgeNormals,
    const labelList& meshEdges,
    labelList& pointType,
    labelList& edgeType,
    labelList& faceType,
    DynamicList<label>& boundaryCornerPoints
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    boolList interfaceEdges(mesh.nEdges(), false);
    forAll(meshEdges, patchEdgeI)
    {
        interfaceEdges[meshEdges[patchEdgeI]] = true;
    }

    syncTools::syncEdgeList
    (
        mesh,
        interfaceEdges,
        orEqOp<bool>(),
        false
     );

    PackedBoolList isMasterEdges = syncTools::getMasterEdges(mesh);

    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        faceType[meshFaceI] = 0;
    }

    forAll(mesh.edges(), edgeI)
    {
        edge e = mesh.edges()[edgeI];
        labelList eFaces = mesh.edgeFaces()[edgeI];

        if (interfaceEdges[edgeI])
        {
            label bFace = -1;

            forAll(eFaces, eFI)
            {
                label faceI = eFaces[eFI];
                if (faceType[faceI] == 0)
                {
                    bFace = faceI;
                    break;
                }
            }

            if (bFace != -1)
            {
                vector eNorm = meshEdgeNormals[edgeI];
                scalar magNorm = mag(eNorm);
                if (magNorm > SMALL)
                {
                    eNorm /= magNorm;

                    vector ePerp = mesh.faceCentres()[bFace] -
                        e.centre(mesh.points());
                    ePerp /= mag(ePerp);

                    scalar dotProd(ePerp&eNorm);

                    if (dotProd > 0.5)
                    {
                        pointType[e[0]] = max(pointType[e[0]],3);
                        pointType[e[1]] = max(pointType[e[1]],3);
                        edgeType[edgeI] = 2;
                    }
                    else if (dotProd < -0.5)
                    {
                        pointType[e[0]] = max(pointType[e[0]],1);
                        pointType[e[1]] = max(pointType[e[1]],1);
                        edgeType[edgeI] = 1;
                    }
                }
            }
        }
    }

    forAll(meshEdges, patchEdgeI)
    {
        label meshEdgeI = meshEdges[patchEdgeI];
        if (edgeType[meshEdgeI] == -1)
        {
            edgeType[meshEdgeI] = 0;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        edgeType,
        maxEqOp<label>(),
        label(-1)
    );

    labelList nConvexEdges(mesh.nPoints(), 0);
    forAll(mesh.points(), pointI)
    {
        const labelList& pEdges = mesh.pointEdges()[pointI];
        forAll(pEdges, pEI)
        {
            label meshEdgeI = pEdges[pEI];
            if (edgeType[meshEdgeI] == 2  && isMasterEdges.get(meshEdgeI) == 1)
            {
                nConvexEdges[pointI]++;
            }
        }
    }

    labelList nConcaveEdges(mesh.nPoints(), 0);
    forAll(mesh.points(), pointI)
    {
        const labelList& pEdges = mesh.pointEdges()[pointI];
        forAll(pEdges, pEI)
        {
            label meshEdgeI = pEdges[pEI];

            if (edgeType[meshEdgeI] == 1  && isMasterEdges.get(meshEdgeI) == 1)
            {
                nConcaveEdges[pointI]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nConvexEdges,
        plusEqOp<label>(),
        label(0)
    );

    syncTools::syncPointList
    (
        mesh,
        nConcaveEdges,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh.points(), pointI)
    {
        if (nConvexEdges[pointI] == 0 && nConcaveEdges[pointI] == 3)
        {
            //Tries to add cells at concave corners
            pointType[pointI] = 2;
        }
        else if (nConvexEdges[pointI] == 3)
        {
            if (nConcaveEdges[pointI] == 3)
            {
                pointType[pointI] = 7;
            }
            else if (nConcaveEdges[pointI] == 0)
            {
                pointType[pointI] = 4;
            }
            else
            {
                FatalErrorInFunction
                    <<"Cannot set pointType: "<<pointI
                    <<" location: "<<mesh.points()[pointI]
                    <<" Number concave edges: "<<nConcaveEdges[pointI]
                    << exit(FatalError);
                pointType[pointI] = -1;
            }
        }
        else if (nConvexEdges[pointI] == 1 && nConcaveEdges[pointI] == 2)
        {
            pointType[pointI] = 5;
        }
        else if (nConvexEdges[pointI] == 2 && nConcaveEdges[pointI] == 2)
        {
            FatalErrorInFunction
                <<"Interface with 2 convex and 2 concave edges found : "
                <<pointI<<" "<<mesh.points()[pointI]
                << exit(FatalError);
            pointType[pointI] = -1;
        }
        else if (nConvexEdges[pointI] == 1 && nConcaveEdges[pointI] == 1)
        {
            FatalErrorInFunction
                <<"Interface with 1 convex and 1 concave edges found : "
                <<pointI<<" "<<mesh.points()[pointI]
                << exit(FatalError);
            pointType[pointI] = -1;
        }
        else if (nConvexEdges[pointI] == 2 && nConcaveEdges[pointI] == 1)
        {
            pointType[pointI] = 6;
        }
    }

    //look for top corner points and mark as convex edges
    //to force creation of corner cell
    labelList nConnections(mesh.nPoints(), 0);
    forAll(pp.localPoints(), i)
    {
        label meshPointI = pp.meshPoints()[i];
        const labelList& pFaces = pp.pointFaces()[i];
        nConnections[meshPointI] += pFaces.size();
    }
    syncTools::syncPointList
    (
        mesh,
        nConnections,
        plusEqOp<label>(),
        label(0)
    );

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label patchLower = mesh.boundaryMesh().findPatchID("lowerLevel");
    label patchHigher = mesh.boundaryMesh().findPatchID("higherLevel");

    DynamicList<label> bPatches(patches.size());
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            if (patchI != patchLower && patchI != patchHigher)
            {
                bPatches.append(patchI);
            }
        }
    }
    bPatches.shrink();

    autoPtr<indirectPrimitivePatch> ppBoundPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            bPatches
        )
    );
    indirectPrimitivePatch& ppBound = ppBoundPtr();

    labelList nBoundConnections(mesh.nPoints(), 0);
    forAll(ppBound.localPoints(), i)
    {
        label meshPointI = ppBound.meshPoints()[i];
        const labelList& pFaces = ppBound.pointFaces()[i];
        nBoundConnections[meshPointI] += pFaces.size();
    }
    syncTools::syncPointList
    (
        mesh,
        nBoundConnections,
        plusEqOp<label>(),
        label(0)
    );

    forAll(ppBound.meshPoints(), i)
    {
        label meshPointI = ppBound.meshPoints()[i];

        if
        (
            pointType[meshPointI] == 3
            && nConnections[meshPointI] == 3 && nBoundConnections[meshPointI] == 3
            && nConvexEdges[meshPointI] == 2 && nConcaveEdges[meshPointI] == 0
        )
        {
            boundaryCornerPoints.append(meshPointI);
            pointType[meshPointI] = 4;
        }
    }
    boundaryCornerPoints.shrink();


    //reset all unset points
    forAll(pp.localPoints(), patchPointI)
    {
        label meshPointI = pp.meshPoints()[patchPointI];
        if (pointType[meshPointI] == -1)
        {
            pointType[meshPointI] = 0;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pointType,
        maxEqOp<label>(),
        label(-1)
    );
}


void Foam::addHexInternalLayers::addPoints
(
    const vectorField& meshNormals,
    const labelList& pointType,
    const labelList& edgeType,
    const labelList& baffleLevel,
    labelList& addedPoints,
    List<DynamicList<labelPair>>& convexEdgePts,
    List<DynamicList<labelPair>>& convexFacePts,
    polyTopoChange& meshMod
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    label nAddedPts = 0;

    forAll(mesh.points(), meshPointI)
    {
        label zoneI = mesh.pointZones().whichZone(meshPointI);

        if
        (
            pointType[meshPointI] == 3
            || pointType[meshPointI] == 4
            || pointType[meshPointI] == 6
        )
        {
            const labelList& pEdges = mesh.pointEdges()[meshPointI];
            forAll(pEdges, pEI)
            {
                label edgeI = pEdges[pEI];

                if (edgeType[edgeI] == -1)
                {
                    point pt = 0.5*
                    (
                        mesh.points()[meshPointI]
                        +mesh.edges()[edgeI].centre(mesh.points())
                    );

                    convexEdgePts[meshPointI].append
                    (
                        labelPair
                        (
                            edgeI,
                            meshMod.setAction
                            (
                                polyAddPoint
                                (
                                    pt,         // point
                                    meshPointI,         // master point
                                    zoneI,      // zone for point
                                    true        // supports a cell
                                 )
                             )
                         )
                    );

                    nAddedPts++;
                }
            }
            convexEdgePts[meshPointI].shrink();

            if
            (
                pointLevel[meshPointI] > baffleLevel[meshPointI]
                && pointType[meshPointI] == 3
            )
            {
                const labelList& pFaces = mesh.pointFaces()[meshPointI];

                forAll(pFaces, pFI)
                {
                    label faceI = pFaces[pFI];
                    if
                    (
                        mesh.isInternalFace(faceI)
                        || patches[patches.whichPatch(faceI)].coupled()
                    )
                    {
                        point pt = 0.5*
                        (
                            mesh.points()[meshPointI]
                            +mesh.faceCentres()[faceI]
                        );

                        convexFacePts[meshPointI].append
                        (
                            labelPair
                            (
                                faceI,
                                meshMod.setAction
                                (
                                    polyAddPoint
                                    (
                                        pt,         // point
                                        meshPointI,
                                        zoneI,      // zone for point
                                        true        // supports a cell
                                    )
                                 )
                             )
                        );
                        nAddedPts++;
                    }
                }
            }
            convexFacePts[meshPointI].shrink();
        }
        else if (pointType[meshPointI] == 0)
        {
            label cellI = mesh.pointCells()[meshPointI][0];

            scalar scalingFac(pointType[meshPointI] == 0 ? 0.5 : 0.5);

            point pt = mesh.points()[meshPointI]
                - meshNormals[meshPointI]
                *(0.5*scalingFac*edge0Len/(1<<cellLevel[cellI]));

            addedPoints[meshPointI] = meshMod.setAction
            (
                polyAddPoint
                (
                    pt,         // point
                    meshPointI, // master point
                    zoneI,      // zone for point
                    true        // supports a cell
                 )
            );
            nAddedPts++;
        }
    }

    Info<<"Number of added points: "
        <<returnReduce(nAddedPts, sumOp<label>())<<endl;
}


void Foam::addHexInternalLayers::mergeRefinementBaffles()
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    polyTopoChange meshMod(mesh);

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            identity(patches.size())
        )
    );
    indirectPrimitivePatch& pp = ppPtr();

    labelList pointMap;
    pointField newPoints;
    mergePoints
    (
        pp.localPoints(),               // points
        1e-10,                        // mergeTol
        false,                        // verbose
        pointMap,
        newPoints
    );

    labelList newMeshPoints(mesh.nPoints(), -1);
    labelList firstOldPoint(newPoints.size(), -1);

    forAll(pp.localPoints(), pointI)
    {
        label meshPointI = pp.meshPoints()[pointI];
        label newPointI = pointMap[pointI];

        if (firstOldPoint[newPointI] == -1)
        {
            firstOldPoint[newPointI] = meshPointI;
            newMeshPoints[meshPointI] = meshPointI;
        }
        else
        {
            newMeshPoints[meshPointI] = firstOldPoint[newPointI];
            // remove duplicate points if already set
            meshMod.setAction
            (
                polyRemovePoint(meshPointI)
            );
        }
    }

    //set all remaining points
    forAll(newMeshPoints, pointI)
    {
        if (newMeshPoints[pointI] == -1)
        {
            newMeshPoints[pointI] = pointI;
        }
    }

    PackedList<1> visitedFace(mesh.nFaces(), 0);
    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();
    const faceZoneMesh& faceZones = mesh.faceZones();
    forAll(baffles_, i)
    {
        label face0 = baffles_[i].first();
        label face1 = baffles_[i].second();

        visitedFace.set(face0, 1);
        visitedFace.set(face1, 1);

        // face1 < 0 signals a coupled face that has been converted to baffle.

        label own0 = faceOwner[face0];
        label own1 = faceOwner[face1];

        if (face1 < 0 || own0 < own1)
        {
            // Use face0 as the new internal face.
            label zoneID = faceZones.whichZone(face0);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
            }

            label nei = (face1 < 0 ? -1 : own1);

            meshMod.setAction(polyRemoveFace(face1));

            face newFace(faces[face0]);
            forAll(newFace, i)
            {
                newFace[i] = newMeshPoints[faces[face0][i]];
            }
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,                // modified face
                    face0,                  // label of face being modified
                    own0,                   // owner
                    nei,                    // neighbour
                    false,                  // face flip
                    -1,                     // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
        else
        {
            // Use face1 as the new internal face.
            label zoneID = faceZones.whichZone(face1);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
            }

            meshMod.setAction(polyRemoveFace(face0));

            face newFace(faces[face1]);
            forAll(newFace, i)
            {
                newFace[i] = newMeshPoints[faces[face1][i]];
            }
            meshMod.setAction
            (
                polyModifyFace
                (
                    newFace,                // modified face
                    face1,                  // label of face being modified
                    own1,                   // owner
                    own0,                   // neighbour
                    false,                  // face flip
                    -1,                     // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }
    }

    forAll(pp.meshPoints(), pointI)
    {
        label meshPointI = pp.meshPoints()[pointI];
        const labelList& pointFaces = mesh.pointFaces()[meshPointI];

        forAll(pointFaces, pI)
        {
            label meshFaceI = pointFaces[pI];
            if (visitedFace.get(meshFaceI) == 0)
            {
                visitedFace.set(meshFaceI, 1);
                label patchID = patches.whichPatch(meshFaceI);
                label own = mesh.faceOwner()[meshFaceI];
                label nei
                (
                    patchID == -1
                    ? mesh.faceNeighbour()[meshFaceI] : -1
                );

                // Use face1 as the new internal face.
                label zoneID = faceZones.whichZone(meshFaceI);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(meshFaceI)];
                }

                face newFace(faces[meshFaceI]);
                forAll(newFace, i)
                {
                    newFace[i] = newMeshPoints[faces[meshFaceI][i]];
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        newFace,                // modified face
                        meshFaceI,              // label of face being modified
                        own,                    // owner
                        nei,                    // neighbour
                        false,                  // face flip
                        patchID,                // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                     )
                );
            }
        }
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

    return;
}


label Foam::addHexInternalLayers::modifyBoundaryFaces
(
    const indirectPrimitivePatch& pp,
    const labelList& pointType,
    const labelList& edgeType,
    const labelList& faceType,
    const List<labelPair>& newFaceCells,
    const List<labelPair>& newEdgeCells,
    const List<labelPair>& newPointCells,
    polyTopoChange& meshMod
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    //Modify boundary faces
    label nModifiedFaces = 0;
    forAll(pp, i)
    {
        label faceI = pp.addressing()[i];
        if (newFaceCells[faceI].first() != -1)
        {
            label newCellI = newFaceCells[faceI].first();
            face f = mesh.faces()[faceI];
            label patchI = patches.whichPatch(faceI);
            label zoneI = mesh.faceZones().whichZone(faceI);
            bool flipZone = false;
            if (zoneI != -1)
            {
                const faceZone& fz = mesh.faceZones()[zoneI];
                flipZone = fz.flipMap()[fz.whichFace(faceI)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    f,          // modified face
                    faceI,      // label of face
                    newCellI,   // owner
                    -1,         // neighbour
                    false,      // face flip
                    patchI, // patch for face
                    false,      // remove from zone
                    zoneI,      // zone for face
                    flipZone        // face flip in zone
                )
            );
            nModifiedFaces++;
        }
    }

    forAll(mesh.edges(), edgeI)
    {
        if (edgeType[edgeI] == 1 && newEdgeCells[edgeI].first() != -1)
        {
            label newCellI = newEdgeCells[edgeI].first();

            const labelList& eFaces = mesh.edgeFaces()[edgeI];

            forAll(eFaces, eFI)
            {
                label faceI = eFaces[eFI];
                face f = mesh.faces()[faceI];
                label patchI = patches.whichPatch(faceI);
                label zoneI = mesh.faceZones().whichZone(faceI);
                bool flipZone = false;
                if (zoneI != -1)
                {
                    const faceZone& fz = mesh.faceZones()[zoneI];
                    flipZone = fz.flipMap()[fz.whichFace(faceI)];
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        f,          // modified face
                        faceI,      // label of face
                        newCellI,   // owner
                        -1,         // neighbour
                        false,      // face flip
                        patchI, // patch for face
                        false,      // remove from zone
                        zoneI,      // zone for face
                        flipZone        // face flip in zone
                     )
                );
                nModifiedFaces++;
            }
        }
    }

    forAll(mesh.points(), pointI)
    {
        if
        (
            (pointType[pointI] == 2 || pointType[pointI] == 5)
            && newPointCells[pointI].first() != -1
        )
        {
            label own = newPointCells[pointI].second();
            label newCellI = newPointCells[pointI].first();

            const labelList& pFaces = mesh.pointFaces()[pointI];
            const labelHashSet cFaceSet(mesh.cells()[own]);

            forAll(pFaces, pFI)
            {
                label faceI = pFaces[pFI];
                if (faceType[faceI] == 0 && cFaceSet.found(faceI))
                {
                    face f = mesh.faces()[faceI];
                    label patchI = patches.whichPatch(faceI);
                    label zoneI = mesh.faceZones().whichZone(faceI);
                    bool flipZone = false;
                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh.faceZones()[zoneI];
                        flipZone = fz.flipMap()[fz.whichFace(faceI)];
                    }

                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            f,          // modified face
                            faceI,      // label of face
                            newCellI,   // owner
                            -1,         // neighbour
                            false,      // face flip
                            patchI, // patch for face
                            false,      // remove from zone
                            zoneI,      // zone for face
                            flipZone        // face flip in zone
                         )
                     );
                    nModifiedFaces++;
                }
            }
        }
    }

    return nModifiedFaces;
}

label Foam::addHexInternalLayers::addTopFaces
(
    const indirectPrimitivePatch& pp,
    const labelList& pointType,
    const labelList& edgeType,
    const labelList& faceType,
    const List<labelPair>& newFaceCells,
    const List<labelPair>& newEdgeCells,
    const List<labelPair>& newPointCells,
    const labelList& addedPoints,
    const List<DynamicList<labelPair>>& convexEdgePts,
    const List<DynamicList<labelPair>>& convexFacePts,
    const labelList& baffleLevel,
    polyTopoChange& meshMod
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const labelListList& cEdgeList = mesh.cellEdges();
    const cellList& cells = mesh.cells();

    label nAddedFaces = 0;

    forAll(pp,i)
    {
        label faceI = pp.addressing()[i];
        if (newFaceCells[faceI].first() != -1)
        {
            label own = newFaceCells[faceI].second();
            label newCellI = newFaceCells[faceI].first();

            face f = mesh.faces()[faceI];
            labelHashSet cEdges(cEdgeList[own]);
            labelHashSet cFaces(cells[own]);
            face newFace(f.size());

            forAll(f, fp)
            {
                label pointI = f[fp];
                if (pointType[pointI] == 0)
                {
                    newFace[fp] = addedPoints[pointI];
                }
                else
                {
                    if (pointLevel[pointI] > baffleLevel[pointI])
                    {
                        forAll(convexFacePts[pointI], cFPI)
                        {
                            label cFaceI = convexFacePts[pointI][cFPI].first();
                            label newPointI = convexFacePts[pointI][cFPI].second();
                            if (cFaces.found(cFaceI))
                            {
                                newFace[fp] = newPointI;
                                break;
                            }
                        }
                    }
                    else
                    {
                        forAll(convexEdgePts[pointI], cCEPI)
                        {
                            label cEdgeI = convexEdgePts[pointI][cCEPI].first();
                            label newPointI = convexEdgePts[pointI][cCEPI].second();
                            if (cEdges.found(cEdgeI))
                            {
                                newFace[fp] = newPointI;
                                break;
                            }
                        }
                    }
                }
            }

            label zoneI = -1;
            bool flip = false;

            //Add top face
            meshMod.setAction
            (
                polyAddFace
                (
                    newFace,   // face
                    own,       // owner
                    newCellI,  // neighbour
                    -1,        // master point
                    -1,        // master edge
                    faceI,     // master face
                    false,     // flux flip
                    -1,        // patch for face
                    zoneI,     // zone for face
                    flip       // face zone flip
                 )
            );
            nAddedFaces++;
        }
    }

    forAll(mesh.edges(), edgeI)
    {
        edge e = mesh.edges()[edgeI];

        if (newEdgeCells[edgeI].first() != -1)
        {
            label own = newEdgeCells[edgeI].second();
            label newCellI = newEdgeCells[edgeI].first();

            face newFace(4);

            if
            (
                edgeType[edgeI] == 1
            )
            {
                labelHashSet cEdges(cEdgeList[own]);
                labelHashSet cFaces(cells[own]);

                labelList grownPts(4);

                label startFaceI = mesh.edgeFaces()[edgeI][0];
                face startFace = mesh.faces()[startFaceI];

                label endPt = e[1];

                label fIndex = findIndex(startFace, e[0]);
                label cIndex = startFace.fcIndex(fIndex);

                if (startFace[cIndex] == e[1])
                {
                    endPt = e[0];
                    cIndex = startFace.fcIndex(cIndex);
                }

                grownPts[0] = startFace[cIndex];
                cIndex = startFace.fcIndex(cIndex);
                grownPts[1] = startFace[cIndex];

                label endFaceI = mesh.edgeFaces()[edgeI][1];
                face endFace = mesh.faces()[endFaceI];

                cIndex = findIndex(endFace, endPt);
                cIndex = endFace.fcIndex(cIndex);
                grownPts[2] = endFace[cIndex];
                cIndex = endFace.fcIndex(cIndex);
                grownPts[3] = endFace[cIndex];

                forAll(grownPts, gPI)
                {
                    label pointI = grownPts[gPI];

                    if (pointType[pointI] == 0)
                    {
                        newFace[gPI] = addedPoints[pointI];
                    }
                    else
                    {
                        if (pointLevel[pointI] > baffleLevel[pointI])
                        {
                            forAll(convexFacePts[pointI], cFPI)
                            {
                                label cFaceI = convexFacePts[pointI][cFPI].first();
                                label newPointI = convexFacePts[pointI][cFPI].second();
                                if (cFaces.found(cFaceI))
                                {
                                    newFace[gPI] = newPointI;
                                    break;
                                }
                            }
                        }
                        else
                        {
                            forAll(convexEdgePts[pointI], cCEPI)
                            {
                                label cEdgeI = convexEdgePts[pointI][cCEPI].first();
                                label newPointI = convexEdgePts[pointI][cCEPI].second();
                                if (cEdges.found(cEdgeI))
                                {
                                    newFace[gPI] = newPointI;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
            else  if (edgeType[edgeI] == 2)
            {
                label startPt = e[0];
                label endPt = e[1];
                if (pointLevel[startPt] <= baffleLevel[startPt])
                {
                    startPt = e[1];
                    endPt = e[0];
                }

                label startFaceI = convexFacePts[startPt][0].first();

                labelHashSet startFaceEdges(mesh.faceEdges()[startFaceI]);
                newFace[0] = convexFacePts[startPt][0].second();
                label endFaceI = convexFacePts[startPt][1].first();
                labelHashSet endFaceEdges(mesh.faceEdges()[endFaceI]);
                newFace[3] = convexFacePts[startPt][1].second();

                forAll(convexEdgePts[endPt], cCEPI)
                {
                    label cEdgeI = convexEdgePts[endPt][cCEPI].first();
                    label newPointI = convexEdgePts[endPt][cCEPI].second();
                    if (startFaceEdges.found(cEdgeI))
                    {
                        newFace[1] = newPointI;
                    }
                    else if (endFaceEdges.found(cEdgeI))
                    {
                        newFace[2] = newPointI;
                    }
                }

                face f = mesh.faces()[startFaceI];

                label cIndex = findIndex(f, startPt);
                label nextPt = f[f.fcIndex(cIndex)];
                if (nextPt == endPt && mesh.faceOwner()[startFaceI] == own)
                {
                    newFace.flip();
                }
                else if (nextPt != endPt && mesh.faceOwner()[startFaceI] != own)
                {
                    newFace.flip();
                }
            }

            label zoneI = -1;
            bool flipZone = false;

            //Add top face
            meshMod.setAction
            (
                polyAddFace
                (
                    newFace,   // face
                    own,       // owner
                    newCellI,  // neighbour
                    -1,        // master point
                    edgeI,     // master edge
                    -1,     // master face
                    false,     // flux flip
                    -1,        // patch for face
                    zoneI,     // zone for face
                    flipZone   // face zone flip
                 )
            );

            nAddedFaces++;
        }
    }


    forAll(mesh.points(), pointI)
    {
        if (newPointCells[pointI].first() != -1)
        {
            label own = newPointCells[pointI].second();
            label newCellI = newPointCells[pointI].first();

            DynamicList<label> newPts(6);
            bool flip(false);

            if (pointType[pointI] == 2)
            {

                const labelList& pFaces = mesh.pointFaces()[pointI];
                face f0 = mesh.faces()[pFaces[0]];
                face f1 = mesh.faces()[pFaces[1]];
                face f2 = mesh.faces()[pFaces[2]];

                label fIndex = findIndex(f0, pointI);
                label cIndex = f0.fcIndex(fIndex);
                cIndex = f0.fcIndex(cIndex);

                newPts.append(addedPoints[f0[cIndex]]);

                cIndex = f0.fcIndex(cIndex);
                labelHashSet fPts(f1);

                if (fPts.found(f0[cIndex]))
                {
                    fIndex = findIndex(f1, pointI);
                    cIndex = f1.fcIndex(fIndex);
                    cIndex = f1.fcIndex(cIndex);
                    newPts.append(addedPoints[f1[cIndex]]);
                    fIndex = findIndex(f2, pointI);
                    cIndex = f2.fcIndex(fIndex);
                    cIndex = f2.fcIndex(cIndex);
                    newPts.append(addedPoints[f2[cIndex]]);
                }
                else
                {
                    fIndex = findIndex(f2, pointI);
                    cIndex = f2.fcIndex(fIndex);
                    cIndex = f2.fcIndex(cIndex);
                    newPts.append(addedPoints[f2[cIndex]]);
                    fIndex = findIndex(f1, pointI);
                    cIndex = f1.fcIndex(fIndex);
                    cIndex = f1.fcIndex(cIndex);
                    newPts.append(addedPoints[f1[cIndex]]);
                }
            }
            else if (pointType[pointI] == 4)
            {
                labelHashSet cFaces(cells[own]);
                const labelList& pFaces = mesh.pointFaces()[pointI];

                label startFaceI = -1;

                forAll(pFaces, pFI)
                {
                    label faceI = pFaces[pFI];

                    if (cFaces.found(faceI))
                    {
                        startFaceI = faceI;
                        break;
                    }
                }
                labelHashSet fEdges(mesh.faceEdges()[startFaceI]);
                label endPt = -1;
                label startEdge = -1;

                forAll(convexEdgePts[pointI], cCEP)
                {
                    label edgeI = convexEdgePts[pointI][cCEP].first();
                    label newPointI = convexEdgePts[pointI][cCEP].second();
                    if (fEdges.found(edgeI))
                    {
                        newPts.append(newPointI);
                        if (startEdge == -1)
                        {
                            startEdge = edgeI;
                        }
                    }
                    else
                    {
                        endPt = newPointI;
                    }
                }
                newPts.append(endPt);

                edge e = mesh.edges()[startEdge];
                label otherPt(e[0] ==  pointI ? e[1] : e[0]);

                const face f = mesh.faces()[startFaceI];

                label cIndex = findIndex(f, pointI);
                label nextPt = f[f.fcIndex(cIndex)];
                if (nextPt != otherPt && mesh.faceOwner()[startFaceI] == own)
                {
                    flip = true;
                }
                else if (nextPt == otherPt && mesh.faceOwner()[startFaceI] != own)
                {
                    flip = true;
                }
            }
            else if (pointType[pointI] == 5)
            {
                labelHashSet cFaces(cells[own]);
                const labelList& pFaces = mesh.pointFaces()[pointI];
                const labelList& pEdges = mesh.pointEdges()[pointI];

                label baseFaceI = -1;

                forAll(pFaces, pFI)
                {
                    label faceI = pFaces[pFI];
                    if (faceType[faceI] == 0 && cFaces.found(faceI))
                    {
                        baseFaceI = faceI;
                    }
                }

                const face baseFace = mesh.faces()[baseFaceI];
                label cIndex = findIndex(baseFace, pointI);
                cIndex = baseFace.fcIndex(cIndex);
                label nextPt = baseFace[cIndex];

                label nextEdgeI = meshTools::findEdge
                (
                    mesh.edges(),
                    pEdges,
                    pointI,
                    nextPt
                );

                const labelList& eFaces = mesh.edgeFaces()[nextEdgeI];
                label startFaceI = -1;
                forAll(eFaces, eFI)
                {
                    label faceI = eFaces[eFI];

                    if (faceI != baseFaceI && cFaces.found(faceI))
                    {
                        startFaceI = faceI;
                        break;
                    }
                }

                label startEdgeI = -1;
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (edgeType[edgeI] == 2)
                    {
                        startEdgeI = edgeI;
                        break;
                    }
                }
                edge e = mesh.edges()[startEdgeI];
                label otherPt(e[0] ==  pointI ? e[1] : e[0]);

                if (convexFacePts[otherPt][0].first() == startFaceI)
                {
                    label newPointI = convexFacePts[otherPt][1].second();
                    newPts.append(newPointI);
                    newPointI = convexFacePts[otherPt][0].second();
                    newPts.append(newPointI);
                }
                else
                {
                    label newPointI = convexFacePts[otherPt][0].second();
                    newPts.append(newPointI);
                    newPointI = convexFacePts[otherPt][1].second();
                    newPts.append(newPointI);
                }

                newPts.append(addedPoints[nextPt]);
                cIndex = baseFace.fcIndex(cIndex);
                nextPt = baseFace[cIndex];
                newPts.append(addedPoints[nextPt]);
                cIndex = baseFace.fcIndex(cIndex);
                nextPt = baseFace[cIndex];
                newPts.append(addedPoints[nextPt]);
            }
            else if (pointType[pointI] == 6)
            {
                label newPointI = convexEdgePts[pointI][0].second();
                label baseEdgeI = convexEdgePts[pointI][0].first();
                newPts.append(newPointI);

                label otherPt = -1;

                const labelList& pEdges = mesh.pointEdges()[pointI];
                labelHashSet cEdges(cEdgeList[own]);
                DynamicList<label> otherEdges(2);

                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (edgeI != baseEdgeI && cEdges.found(edgeI))
                    {
                        otherEdges.append(edgeI);
                    }
                }

                label startFaceI = -1;
                label nextPt = -1;
                forAll(otherEdges, oEI)
                {
                    edge nextEdge = mesh.edges()[otherEdges[oEI]];
                    nextPt = (nextEdge[0] == pointI ? nextEdge[1] : nextEdge[0]);

                    label faceI  = convexFacePts[nextPt][0].first();

                    labelHashSet faceEdges(mesh.faceEdges()[faceI]);

                    bool reverse = false;

                    if (oEI == 0)
                    {
                        startFaceI = faceI;
                        otherPt = nextPt;
                    }

                    if (!faceEdges.found(baseEdgeI) && oEI == 0)
                    {
                        startFaceI = convexFacePts[nextPt][1].first();
                        reverse = true;
                    }
                    else if (faceEdges.found(baseEdgeI) && oEI == 1)
                    {
                        reverse = true;
                    }

                    if (reverse)
                    {
                        newPointI = convexFacePts[nextPt][1].second();
                        newPts.append(newPointI);
                        newPointI = convexFacePts[nextPt][0].second();
                        newPts.append(newPointI);
                    }
                    else
                    {
                        newPointI = convexFacePts[nextPt][0].second();
                        newPts.append(newPointI);
                        newPointI = convexFacePts[nextPt][1].second();
                        newPts.append(newPointI);
                    }
                }

                const face f = mesh.faces()[startFaceI];

                label cIndex = findIndex(f, pointI);
                nextPt = f[f.fcIndex(cIndex)];

                if (nextPt == otherPt && mesh.faceOwner()[startFaceI] == own)
                {
                    flip = true;
                }
                else if (nextPt != otherPt && mesh.faceOwner()[startFaceI] != own)
                {
                    flip = true;
                }
            }
            else if (pointType[pointI] == 7)
            {
                labelHashSet cEdges(cEdgeList[own]);
                const labelList& pEdges = mesh.pointEdges()[pointI];

                label baseEdgeI = -1;

                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cEdges.found(edgeI))
                    {
                        baseEdgeI = edgeI;
                        break;
                    }
                }

                edge e = mesh.edges()[baseEdgeI];
                label otherPt(e[0] ==  pointI ? e[1] : e[0]);

                label newPointI = convexFacePts[otherPt][0].second();
                label startFaceI = convexFacePts[otherPt][0].first();
                newPts.append(newPointI);
                newPointI = convexFacePts[otherPt][1].second();
                newPts.append(newPointI);

                DynamicList<label> otherEdges(2);

                bool switchOrder = false;
                label nextPt = -1;
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (edgeI != baseEdgeI && cEdges.found(edgeI))
                    {
                        otherEdges.append(edgeI);
                        if (otherEdges.size() == 1)
                        {
                            edge e = mesh.edges()[edgeI];

                            nextPt = (e[0] ==  pointI ? e[1] : e[0]);
                            if
                            (
                                convexFacePts[nextPt][0].first() == startFaceI
                                || convexFacePts[nextPt][1].first() == startFaceI
                            )
                            {
                                switchOrder = true;
                            }
                        }
                    }
                }

                if (switchOrder)
                {
                    label oEdgeI =  otherEdges[0];
                    otherEdges[0] = otherEdges[1];
                    otherEdges[1] = oEdgeI;
                }

                forAll(otherEdges, oEI)
                {
                    edge nextEdge = mesh.edges()[otherEdges[oEI]];
                    nextPt = (nextEdge[0] == pointI ? nextEdge[1] : nextEdge[0]);
                    label faceI  = convexFacePts[nextPt][0].first();

                    labelHashSet faceEdges(mesh.faceEdges()[faceI]);

                    bool reverse = false;

                    if (!faceEdges.found(baseEdgeI) && oEI == 0)
                    {
                        reverse = true;
                    }
                    else if (faceEdges.found(baseEdgeI) && oEI == 1)
                    {
                        reverse = true;
                    }

                    if (reverse)
                    {
                        newPointI = convexFacePts[nextPt][1].second();
                        newPts.append(newPointI);
                        newPointI = convexFacePts[nextPt][0].second();
                        newPts.append(newPointI);
                    }
                    else
                    {
                        newPointI = convexFacePts[nextPt][0].second();
                        newPts.append(newPointI);
                        newPointI = convexFacePts[nextPt][1].second();
                        newPts.append(newPointI);
                    }
                }

                const face f = mesh.faces()[startFaceI];

                label cIndex = findIndex(f, pointI);
                nextPt = f[f.fcIndex(cIndex)];
                if (nextPt == otherPt && mesh.faceOwner()[startFaceI] == own)
                {
                    flip = true;
                }
                else if (nextPt != otherPt && mesh.faceOwner()[startFaceI] != own)
                {
                    flip = true;
                }
            }

            face newFace(newPts.shrink());

            if (flip)
            {
                newFace.flip();
            }

            label zoneI = -1;
            bool flipZone = false;

            //Add top face
            meshMod.setAction
            (
                polyAddFace
                (
                    newFace,   // face
                    own,       // owner
                    newCellI,  // neighbour
                    pointI,        // master point
                    -1,        // master edge
                    -1,     // master face
                    false,     // flux flip
                    -1,        // patch for face
                    zoneI,     // zone for face
                    flipZone       // face zone flip
                 )
            );
            nAddedFaces++;
        }
    }

    return nAddedFaces;
}


label Foam::addHexInternalLayers::addSideFaces
(
    const indirectPrimitivePatch& pp,
    const labelList& pointType,
    const labelList& edgeType,
    const labelList& faceType,
    const List<labelPair>& newFaceCells,
    const List<labelPair>& newEdgeCells,
    const List<labelPair>& newPointCells,
    const labelList& addedPoints,
    const List<DynamicList<labelPair>>& convexEdgePts,
    const List<DynamicList<labelPair>>& convexFacePts,
    const labelList& baffleLevel,
    const labelList& meshEdges,
    polyTopoChange& meshMod
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const label patchLower = mesh.boundaryMesh().findPatchID("lowerLevel");
    const label patchHigher = mesh.boundaryMesh().findPatchID("higherLevel");

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const labelListList& cEdgeList = mesh.cellEdges();
    const cellList& cells = mesh.cells();

    label nAddedFaces = 0;

    //Add edge faces
    forAll(pp.edges(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        const edge e = mesh.edges()[meshEdgeI];

        if
        (
            pointType[e[0]] == 1 || pointType[e[1]] == 1
            || pointType[e[0]] == 2 || pointType[e[1]] == 2
            || pointType[e[0]] == 5  || pointType[e[1]] == 5
            || pointType[e[0]] == 7 || pointType[e[1]] == 7
        )
        {
            continue;
        }

        bool validEdge = true;
        DynamicList<label> grownUpFaces(2);
        label internalFaceI = -1;
        label cornerPt = -1;

        if (edgeType[meshEdgeI] == 0 || edgeType[meshEdgeI] == 2)
        {
            const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];
            forAll(eFaces, eFI)
            {
                label faceI = eFaces[eFI];
                label patchI = patches.whichPatch(faceI);
                if (patchI == -1)
                {
                    internalFaceI = faceI;
                    grownUpFaces.append(faceI);
                }
                else if (patchI != -1 && faceType[faceI] == -1)
                {
                    if (pointType[e[0]] == 6 || pointType[e[1]] == 6)
                    {
                        cornerPt = (pointType[e[0]] == 6 ? e[0] : e[1]);
                        face f = mesh.faces()[faceI];
                        label index = findIndex(f, cornerPt);
                        label nextPt = f[f.fcIndex(index)];
                        label prevPt = f[f.rcIndex(index)];

                        label nextEdge = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[cornerPt],
                            cornerPt,
                            nextPt
                        );

                        label prevEdge = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[cornerPt],
                            cornerPt,
                            prevPt
                        );

                        if (edgeType[prevEdge] != 2 || edgeType[nextEdge] != 2)
                        {
                            grownUpFaces.append(faceI);
                        }
                    }
                    else
                    {
                        grownUpFaces.append(faceI);
                    }
                }
            }
            if (edgeType[meshEdgeI] == 2 && internalFaceI != -1)
            {
                validEdge = false;
            }
        }
        else
        {
            validEdge = false;
        }

        if (cornerPt != -1 && validEdge && grownUpFaces.size() == 0)
        {
            validEdge = false;
        }
        else if (validEdge && grownUpFaces.size() == 0)
        {
            grownUpFaces.append(-1);
        }
        grownUpFaces.shrink();

        if (validEdge)
        {
            forAll(grownUpFaces, gFI)
            {
                bool flip = false;
                face newFace(4);
                label own = -1;
                label nei = -1;

                label origFace = grownUpFaces[gFI];

                label patchI = -1;
                if (origFace != -1)
                {
                    patchI = patches.whichPatch(origFace);
                }
                const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];

                DynamicList<label> boundaryFaces(2);

                if (patchI != -1)
                {
                    labelHashSet cFacesOwn(cells[mesh.faceOwner()[origFace]]);

                    forAll(eFaces, eFI)
                    {
                        label faceI = eFaces[eFI];
                        if (faceType[faceI] == 0 && cFacesOwn.found(faceI))
                        {
                           boundaryFaces.append(faceI);
                           break;
                        }
                    }
                }
                else
                {
                    forAll(eFaces, eFI)
                    {
                        label faceI = eFaces[eFI];
                        if (faceType[faceI] == 0)
                        {
                            boundaryFaces.append(faceI);
                        }
                    }
                }

                boundaryFaces.shrink();

                label checkCellI = mesh.faceOwner()[boundaryFaces[0]];

                labelHashSet cEdges(cEdgeList[checkCellI]);
                labelHashSet cFaces(cells[checkCellI]);

                bool swap = false;
                //Set owner and neighbour
                if (boundaryFaces.size() == 1)
                {
                    own = newFaceCells[boundaryFaces[0]].first();
                    if (own == -1)
                    {
                        face f = mesh.faces()[boundaryFaces[0]];
                        label fIndex0 = findIndex(f,e[0]);
                        label fIndex1 = findIndex(f,e[1]);

                        fIndex0 = f.fcIndex(fIndex0);
                        fIndex0 = f.fcIndex(fIndex0);
                        fIndex1 = f.fcIndex(fIndex1);
                        fIndex1 = f.fcIndex(fIndex1);

                        label oppEdge = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[f[fIndex0]],
                            f[fIndex0],
                            f[fIndex1]
                         );

                        own = newEdgeCells[oppEdge].second();

                        if (own == -1)
                        {
                            forAll(f,fp)
                            {
                                label pointI = f[fp];

                                if (newPointCells[pointI].second() != -1)
                                {
                                    own = newPointCells[pointI].second();
                                    break;
                                }
                            }
                        }
                    }
                }
                else
                {
                    label oldCell0 = newFaceCells[boundaryFaces[0]].second();
                    label newCell0 = -1;

                    if (oldCell0 == -1)
                    {
                        face f = mesh.faces()[boundaryFaces[0]];
                        label fIndex0 = findIndex(f,e[0]);
                        label fIndex1 = findIndex(f,e[1]);

                        fIndex0 = f.fcIndex(fIndex0);
                        fIndex0 = f.fcIndex(fIndex0);
                        fIndex1 = f.fcIndex(fIndex1);
                        fIndex1 = f.fcIndex(fIndex1);

                        label oppEdge = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[f[fIndex0]],
                            f[fIndex0],
                            f[fIndex1]
                        );

                        oldCell0 = newEdgeCells[oppEdge].second();

                        if (oldCell0 == -1)
                        {
                            forAll(f,fp)
                            {
                                label pointI = f[fp];

                                if (newPointCells[pointI].second() != -1)
                                {
                                    oldCell0 = newPointCells[pointI].second();
                                    newCell0 = newPointCells[pointI].first();
                                    break;
                                }
                            }
                        }
                        else
                        {
                            newCell0 = newEdgeCells[oppEdge].first();
                        }
                    }
                    else
                    {
                        newCell0 = newFaceCells[boundaryFaces[0]].first();
                    }

                    label oldCell1 = newFaceCells[boundaryFaces[1]].second();
                    label newCell1 = -1;

                    if (oldCell1 == -1)
                    {
                        face f = mesh.faces()[boundaryFaces[1]];
                        label fIndex0 = findIndex(f,e[0]);
                        label fIndex1 = findIndex(f,e[1]);

                        fIndex0 = f.fcIndex(fIndex0);
                        fIndex0 = f.fcIndex(fIndex0);
                        fIndex1 = f.fcIndex(fIndex1);
                        fIndex1 = f.fcIndex(fIndex1);

                        label oppEdge = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[f[fIndex0]],
                            f[fIndex0],
                            f[fIndex1]
                        );

                        oldCell1 = newEdgeCells[oppEdge].second();

                        if (oldCell1 == -1)
                        {
                            forAll(f,fp)
                            {
                                label pointI = f[fp];

                                if (newPointCells[pointI].second() != -1)
                                {
                                    oldCell1 = newPointCells[pointI].second();
                                    newCell1 = newPointCells[pointI].first();
                                    break;
                                }
                            }
                        }
                        else
                        {
                            newCell1 = newEdgeCells[oppEdge].first();
                        }
                    }
                    else
                    {
                        newCell1 = newFaceCells[boundaryFaces[1]].first();
                    }

                    if (newCell1 > newCell0)
                    {
                        swap = false;
                            //swap  = (oldCell1 < oldCell0 ? true : false);
                        own = newCell0;
                        nei = newCell1;
                    }
                    else
                    {
                        swap = true;
                            //swap  = (oldCell0 < oldCell1 ? true : false);
                        own = newCell1;
                        nei = newCell0;
                    }
                }

                const face f = mesh.faces()[boundaryFaces[0]];
                label cIndex = findIndex(f, e[0]);

                if
                (
                    (f[f.fcIndex(cIndex)] == e[1] && swap)
                    || (f[f.rcIndex(cIndex)] == e[1] && !swap)
                )
                {
                    flip = true;
                }

                label nCount = 0;
                newFace[nCount++] = e[1];
                newFace[nCount++] = e[0];

                forAll(e, ePI)
                {
                    label pointI = e[ePI];

                    if (addedPoints[pointI] != -1)
                    {
                        newFace[nCount++] = addedPoints[pointI];
                    }
                    else if (pointLevel[pointI] > baffleLevel[pointI])
                    {
                        forAll(convexFacePts[pointI], cFPI)
                        {
                            label faceI = convexFacePts[pointI][cFPI].first();
                            if (cFaces.found(faceI))
                            {
                                label newPointI =
                                    convexFacePts[pointI][cFPI].second();
                                newFace[nCount++] = newPointI;
                                break;
                            }
                        }
                    }
                    else
                    {
                        forAll(convexEdgePts[pointI], cCEPI)
                        {
                            label nbrEdgeI =
                                convexEdgePts[pointI][cCEPI].first();
                            if (cEdges.found(nbrEdgeI))
                            {
                                label newPointI =
                                    convexEdgePts[pointI][cCEPI].second();
                                newFace[nCount++] = newPointI;
                                break;
                            }
                        }
                    }
                }

                if (flip)
                {
                    newFace.flip();
                }

                label zoneI = -1;
                bool flipZone = false;
                if (origFace != -1)
                {
                    zoneI = mesh.faceZones().whichZone(origFace);

                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh.faceZones()[zoneI];
                        flipZone = fz.flipMap()[fz.whichFace(origFace)];
                    }
                }

                //Add top face
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,   // face
                        own,       // owner
                        nei,  // neighbour
                        -1,        // master point
                        -1,        // master edge
                        origFace,     // master face
                        false,     // flux flip
                        patchI,        // patch for face
                        zoneI,     // zone for face
                        flipZone       // face zone flip
                     )
                 );
                nAddedFaces++;
            }
        }
    }


    //Add convex edge faces
    forAll(mesh.edges(), edgeI)
    {
        edge e = mesh.edges()[edgeI];
        if
        (
            (edgeType[edgeI] == 2 && newEdgeCells[edgeI].first() != -1)
            ||
            (
                edgeType[edgeI] == 2 && pointType[e[0]] == 6
                && newPointCells[e[0]].first() != -1
            )
            ||
            (
                edgeType[edgeI] == 2 && pointType[e[1]] == 6
                && newPointCells[e[1]].first() != -1
            )
        )
        {
            label origOwnI = -1;
            label newEdgeCellI = -1;

            label cornerPt = -1;
            if (newEdgeCells[edgeI].first() != -1)
            {
                origOwnI = newEdgeCells[edgeI].second();
                newEdgeCellI = newEdgeCells[edgeI].first();
            }
            else if (pointType[e[0]] == 6)
            {
                origOwnI = newPointCells[e[0]].second();
                newEdgeCellI = newPointCells[e[0]].first();
                cornerPt = e[0];
            }
            else if (pointType[e[1]] == 6)
            {
                origOwnI = newPointCells[e[1]].second();
                newEdgeCellI = newPointCells[e[1]].first();
                cornerPt = e[1];
             }
            else
            {
                continue;
            }

            labelHashSet cFaces(cells[origOwnI]);
            const labelList& eFaces = mesh.edgeFaces()[edgeI];

            forAll(eFaces, eFI)
            {
                label faceI = eFaces[eFI];
                const face f = mesh.faces()[faceI];

                if (!cFaces.found(faceI))
                {
                    continue;
                }
                else if (cornerPt != -1)
                {
                    label index = findIndex(f, cornerPt);
                    label nextPt = f[f.fcIndex(index)];
                    label prevPt = f[f.rcIndex(index)];

                    label nextEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[cornerPt],
                        cornerPt,
                        nextPt
                    );

                    label prevEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[cornerPt],
                        cornerPt,
                        prevPt
                    );

                    if (edgeType[nextEdge] != -1 && edgeType[prevEdge] != -1)
                    {
                        continue;

                    }
                }

                bool flip = false;
                face newFace(4);
                label own = -1;
                label nei = -1;

                label patchI = patches.whichPatch(faceI);

                bool swap = false;
                 //Set owner and neighbour
                if (patchI != -1)
                {
                    own = newEdgeCellI;
                }
                else
                {
                    own = newEdgeCellI;

                    label origCornerNei
                    (
                        mesh.faceOwner()[faceI] == origOwnI
                        ? mesh.faceNeighbour()[faceI]
                        : mesh.faceOwner()[faceI]
                    );
                    labelHashSet cFacesNei(cells[origCornerNei]);
                    label origNei = -1;

                    forAll(eFaces, eFJ)
                    {
                        label nbrFace = eFaces[eFJ];

                        if (faceType[nbrFace] == 0 && cFacesNei.found(nbrFace))
                        {
                            if (newFaceCells[nbrFace].first() != -1)
                            {
                                nei = newFaceCells[nbrFace].first();
                                origNei = newFaceCells[nbrFace].second();
                            }
                            else
                            {
                                const labelList& fNeiEdges = mesh.faceEdges()[nbrFace];
                                forAll(fNeiEdges, fNEI)
                                {
                                    label nEdgeI = fNeiEdges[fNEI];
                                    if (nEdgeI != edgeI && newEdgeCells[nEdgeI].first() != -1)
                                    {
                                        nei = newEdgeCells[nEdgeI].first();
                                        origNei = newEdgeCells[nEdgeI].second();
                                        break;
                                    }
                                }
                            }

                            break;
                        }
                    }

                    if (own > nei)
                    {
                        label ownTmp = own;
                        own = nei;
                        nei = ownTmp;
                        swap  = (origOwnI > origNei ? true : false);
                    }
                    else
                    {
                        swap  = (origOwnI < origNei ? true : false);
                    }

                }

                label cIndex = findIndex(f, e[0]);

                if (patchI != -1)
                {
                    if (f[f.fcIndex(cIndex)] == e[1])
                    {
                        flip = true;
                    }
                }
                else
                {
                    if
                    (
                        (f[f.fcIndex(cIndex)] == e[1] && swap)
                        || (f[f.rcIndex(cIndex)] == e[1] && !swap)
                     )
                    {
                        flip = true;
                    }
                }

                label nCount = 0;
                newFace[nCount++] = e[1];
                newFace[nCount++] = e[0];

                forAll(e, ePI)
                {
                    label pointI = e[ePI];

                    if (addedPoints[pointI] != -1)
                    {
                        newFace[nCount++] = addedPoints[pointI];
                    }
                    else if (pointLevel[pointI] > baffleLevel[pointI])
                    {
                        forAll(convexFacePts[pointI], cFPI)
                        {
                            label convexFaceI =
                                convexFacePts[pointI][cFPI].first();
                            if (convexFaceI == faceI)
                            {
                                label newPointI =
                                    convexFacePts[pointI][cFPI].second();
                                newFace[nCount++] = newPointI;
                                break;
                            }
                        }
                    }
                    else
                    {
                        labelHashSet fEdges(mesh.faceEdges()[faceI]);

                        forAll(convexEdgePts[pointI], cCEPI)
                        {
                            label nbrEdgeI =
                                convexEdgePts[pointI][cCEPI].first();
                            if (fEdges.found(nbrEdgeI))
                            {
                                label newPointI =
                                    convexEdgePts[pointI][cCEPI].second();
                                newFace[nCount++] = newPointI;
                                break;
                            }
                        }
                    }
                }

                if (flip)
                {
                    newFace.flip();
                }


                label zoneI = -1;
                bool flipZone = false;
                if (faceI != -1)
                {
                    zoneI = mesh.faceZones().whichZone(faceI);

                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh.faceZones()[zoneI];
                        flipZone = fz.flipMap()[fz.whichFace(faceI)];
                    }
                }

                //Add top face
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,   // face
                        own,       // owner
                        nei,  // neighbour
                        -1,        // master point
                        -1,        // master edge
                        faceI,     // master face
                        false,     // flux flip
                        patchI,        // patch for face
                        zoneI,     // zone for face
                        flipZone       // face zone flip
                     )
                );
                nAddedFaces++;
            }
        }
    }

    //Add point faces
    forAll(mesh.points(), pointI)
    {
        bool cornerPt = false;

        DynamicList<label>  grownFaces(3);

        const labelList pFaces = mesh.pointFaces()[pointI];

        if
        (
            pointType[pointI] == 1
            || (pointType[pointI] >= 5 && newPointCells[pointI].first() == -1)
        )
        {
            forAll(pFaces, pFI)
            {
                label faceI = pFaces[pFI];
                if (faceType[faceI] == -1)
                {
                    face f = mesh.faces()[faceI];
                    label index = findIndex(f, pointI);
                    label nextPt = f[f.fcIndex(index)];
                    label prevPt = f[f.rcIndex(index)];

                    label nextEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[pointI],
                        pointI,
                        nextPt
                    );

                    label prevEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[pointI],
                        pointI,
                        prevPt
                    );

                    if (edgeType[nextEdge] != -1 && edgeType[prevEdge] != -1)
                    {
                        grownFaces.append(faceI);
                    }
                }
            }
            if (pointType[pointI] == 1 && grownFaces.size() == 0)
            {
                grownFaces.append(-1);
            }
            else if
            (
                pointType[pointI] >= 5 && newPointCells[pointI].first() == -1
            )
            {
                cornerPt = true;
            }
        }

        if (grownFaces.size() > 0)
        {

            forAll(grownFaces, gFI)
            {
                bool flip = false;
                face newFace(5);
                label origFace = grownFaces[gFI];
                label own = -1;
                label nei = -1;
                label patchI = -1;

                if
                (
                    pointType[pointI] == 1
                    && pointLevel[pointI] > baffleLevel[pointI]
                )
                {
                    label nCount = 0;
                    newFace[nCount++] = pointI;
                    const labelList& pEdges = mesh.pointEdges()[pointI];

                    label ownEdge = -1;
                    label firstEdge = -1;
                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];
                        edge e = mesh.edges()[edgeI];
                        label otherPt(e[0] == pointI ? e[1] : e[0]);

                        if (edgeType[edgeI] != 1)
                        {
                            if (nCount == 1)
                            {
                                firstEdge = edgeI;
                                newFace[nCount++] = otherPt;
                                newFace[nCount++] = addedPoints[otherPt];
                            }
                            else
                            {
                                newFace[nCount++] = addedPoints[otherPt];
                                newFace[nCount++] = otherPt;
                            }
                        }
                        else
                        {
                            label newCellI = newEdgeCells[edgeI].first();
                            if (newCellI == -1)
                            {
                                newCellI = newPointCells[otherPt].first();
                            }

                            if (own == -1)
                            {
                                ownEdge = edgeI;
                                own = newCellI;
                                nei = -1;
                            }
                            else
                            {
                                if (newCellI > own)
                                {
                                    nei = newCellI;
                                }
                                else
                                {
                                    ownEdge = edgeI;
                                    nei = own;
                                    own = newCellI;
                                }
                            }
                        }
                    }

                    const labelList& eFaces = mesh.edgeFaces()[ownEdge];
                    forAll(eFaces, eFI)
                    {
                        label faceI = eFaces[eFI];
                        labelHashSet fEdges(mesh.faceEdges()[faceI]);
                        if (fEdges.found(firstEdge))
                        {
                            const face f = mesh.faces()[faceI];
                            label cIndex = findIndex(f, pointI);
                            if (f[f.fcIndex(cIndex)] == newFace[1])
                            {
                                flip = true;
                            }
                        }
                    }
                }
                else
                {
                    patchI = patches.whichPatch(origFace);
                    label checkCellI = mesh.faceOwner()[origFace];

                    labelHashSet cEdges(cEdgeList[checkCellI]);
                    labelHashSet cFaces(cells[checkCellI]);

                    label nCount = 0;
                    newFace[nCount++] = pointI;
                    const labelList& pEdges = mesh.pointEdges()[pointI];

                    label ownEdge = -1;
                    label firstEdge = -1;

                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];

                        if (cornerPt && !cEdges.found(edgeI))
                        {
                            continue;
                        }

                        if (edgeType[edgeI] != 1)
                        {
                            edge e = mesh.edges()[edgeI];
                            label otherPt(e[0] == pointI ? e[1] : e[0]);

                            label newAddedPt = -1;

                            if (addedPoints[otherPt] != -1)
                            {
                                newAddedPt = addedPoints[otherPt];
                            }
                            else if (pointLevel[otherPt] > baffleLevel[otherPt])
                            {
                                forAll(convexFacePts[otherPt], cFPI)
                                {
                                    label faceI =
                                        convexFacePts[otherPt][cFPI].first();
                                    if (cFaces.found(faceI))
                                    {
                                        newAddedPt =
                                            convexFacePts[otherPt][cFPI].second();
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                forAll(convexEdgePts[otherPt], cCEPI)
                                {
                                    label nbrEdgeI =
                                        convexEdgePts[otherPt][cCEPI].first();
                                    if (cEdges.found(nbrEdgeI))
                                    {
                                        newAddedPt =
                                            convexEdgePts[otherPt][cCEPI].second();
                                        break;
                                    }
                                }
                            }

                            if (nCount == 1)
                            {
                                firstEdge = edgeI;
                                newFace[nCount++] = otherPt;
                                newFace[nCount++] = newAddedPt;
                            }
                            else
                            {
                                newFace[nCount++] = newAddedPt;
                                newFace[nCount++] = otherPt;
                            }
                        }
                        else if (edgeType[edgeI] == 1)
                        {
                            label newCellI = newEdgeCells[edgeI].first();
                            if (own == -1)
                            {
                                ownEdge = edgeI;
                                own = newCellI;
                                nei = -1;
                            }
                            else
                            {
                                if (newCellI > own)
                                {
                                    nei = newCellI;
                                }
                                else
                                {
                                    ownEdge = edgeI;
                                    nei = own;
                                    own = newCellI;
                                }
                            }
                        }
                    }

                    const labelList& eFaces = mesh.edgeFaces()[ownEdge];
                    forAll(eFaces, eFI)
                    {
                        label faceI = eFaces[eFI];
                        labelHashSet fEdges(mesh.faceEdges()[faceI]);
                        if (fEdges.found(firstEdge))
                        {
                            const face f = mesh.faces()[faceI];
                            label cIndex = findIndex(f, pointI);
                            if (f[f.fcIndex(cIndex)] == newFace[1])
                            {
                                flip = true;
                            }
                        }
                    }
                }

                if (flip)
                {
                    newFace.flip();
                }

                label zoneI = -1;
                bool flipZone = false;
                if (origFace != -1)
                {
                    zoneI = mesh.faceZones().whichZone(origFace);

                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh.faceZones()[zoneI];
                        flipZone = fz.flipMap()[fz.whichFace(origFace)];
                    }
                }

                //Add top face
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,   // face
                        own,       // owner
                        nei,  // neighbour
                        -1,        // master point
                        -1,        // master edge
                        origFace,     // master face
                        false,     // flux flip
                        patchI,        // patch for face
                        zoneI,     // zone for face
                        flipZone       // face zone flip
                     )
                 );
                nAddedFaces++;
            }
        }
        else if
        (
            pointType[pointI] == 3
            || (pointType[pointI] == 4 && newPointCells[pointI].first() == -1)
        )
        {

            bool bPt = false;
            const labelList& pEdges = mesh.pointEdges()[pointI];

            if (pointType[pointI] == 3)
            {
                const labelList& pFaces = mesh.pointFaces()[pointI];
                forAll(pFaces, pFI)
                {
                    const label faceI = pFaces[pFI];
                    label patchI = patches.whichPatch(faceI);
                    if (patchI != -1 && !patches[patchI].coupled())
                    {
                        if (patchI != patchLower && patchI != patchHigher)
                        {
                            bPt = true;
                            break;
                        }
                    }
                }
            }

            if (bPt)
            {
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if
                    (
                        edgeType[edgeI] == 2
                        && newEdgeCells[edgeI].second() != -1
                    )
                    {
                        label checkCellI = newEdgeCells[edgeI].second();
                        labelHashSet cFaces(cells[checkCellI]);

                        const labelList& pFaces = mesh.pointFaces()[pointI];
                        forAll(pFaces, pFI)
                        {
                            label faceI = pFaces[pFI];
                            if (cFaces.found(faceI))
                            {
                                face f = mesh.faces()[faceI];
                                label index = findIndex(f, pointI);
                                label nextPt = f[f.fcIndex(index)];
                                label prevPt = f[f.rcIndex(index)];

                                label nextEdge = meshTools::findEdge
                                (
                                    mesh.edges(),
                                    mesh.pointEdges()[pointI],
                                    pointI,
                                    nextPt
                                );

                                label prevEdge = meshTools::findEdge
                                (
                                    mesh.edges(),
                                    mesh.pointEdges()[pointI],
                                    pointI,
                                    prevPt
                                );

                                if
                                (
                                    edgeType[nextEdge] == -1
                                    && edgeType[prevEdge] == -1
                                )
                                {
                                    if (patches.whichPatch(faceI) == -1)
                                    {
                                        bPt = false;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            DynamicList<label> addedEdges(3);

            bool pType4 =  false;
            if (pointType[pointI] == 4 && newPointCells[pointI].first() == -1)
            {
                pType4 = true;
            }

            if (pType4)
            {
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (edgeType[edgeI] == 2 && newEdgeCells[edgeI].second() != -1)
                    {
                        addedEdges.append(edgeI);
                    }
                }
            }
            else
            {
                label cornerEdge = -1;

                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];

                    if
                    (
                        edgeType[edgeI] == 2
                        && newEdgeCells[edgeI].second() != -1
                    )
                    {
                        addedEdges.append(edgeI);
                        if (!bPt)
                        {
                            break;
                        }
                    }
                    else if (edgeType[edgeI] == 2)
                    {
                        if (newEdgeCells[edgeI].second() == -1)
                        {
                            edge e = mesh.edges()[edgeI];
                            label cornerPt(e[0] ==  pointI ? e[1] : e[0]);


                            const labelHashSet
                                pCellSet(mesh.pointCells()[pointI]);

                            if (pCellSet.found(newPointCells[cornerPt].second()))
                            {
                                cornerEdge = edgeI;
                            }
                        }
                    }
                }
                if (addedEdges.size() == 0 && cornerEdge != -1)
                {
                    addedEdges.append(cornerEdge);
                }
            }
            addedEdges.shrink();

            if (addedEdges.size() == 0)
            {
                continue;
            }

            forAll(addedEdges, aEI)
            {
                bool flip = false;
                face newFace(3);
                label origFace = -1;
                label own = -1;
                label nei = -1;
                label patchI = -1;

                pointField normPts(3);
                face norm(identity(3));

                label nCount = 0;
                newFace[nCount] = pointI;
                normPts[nCount] = mesh.points()[pointI];
                nCount++;

                if
                (
                    pointType[pointI] == 3
                    && pointLevel[pointI] > baffleLevel[pointI]
                )
                {
                    forAll(convexFacePts[pointI], cFPI)
                    {
                        label faceI = convexFacePts[pointI][cFPI].first();
                        newFace[nCount] = convexFacePts[pointI][cFPI].second();
                        normPts[nCount] = mesh.faceCentres()[faceI];
                        nCount++;
                    }

                    const labelList& pEdges = mesh.pointEdges()[pointI];

                    label ownEdge = -1;
                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];
                        if (edgeType[edgeI] == 2)
                        {
                            label newCellI = newEdgeCells[edgeI].first();

                            if (newCellI == -1)
                            {
                                edge e = mesh.edges()[edgeI];
                                label cornerPt(e[0] ==  pointI ? e[1] : e[0]);

                                newCellI = newPointCells[cornerPt].first();
                            }

                            if (newCellI == -1)
                            {
                                continue;
                            }

                            if (own == -1)
                            {
                                ownEdge = edgeI;
                                own = newCellI;
                                nei = -1;
                            }
                            else
                            {
                                if (newCellI > own)
                                {
                                    nei = newCellI;
                                }
                                else
                                {
                                    ownEdge = edgeI;
                                    nei = own;
                                    own = newCellI;
                                }
                            }
                        }
                    }

                    edge oEdge = mesh.edges()[ownEdge];
                    label otherPt(oEdge[0] ==  pointI ? oEdge[1] : oEdge[0]);
                    vector eVec = mesh.points()[pointI] - mesh.points()[otherPt];

                    if ((eVec & norm.areaNormal(normPts)) < 0)
                    {
                        flip = true;
                    }
                }
                else
                {
                    label checkCellI = newEdgeCells[addedEdges[aEI]].second();

                    labelHashSet cEdges(cEdgeList[checkCellI]);
                    labelHashSet cFaces(cells[checkCellI]);

                    forAll(convexEdgePts[pointI], cCEPI)
                    {
                        label nbrEdgeI = convexEdgePts[pointI][cCEPI].first();
                        if (cEdges.found(nbrEdgeI))
                        {
                            newFace[nCount] =
                                convexEdgePts[pointI][cCEPI].second();
                            normPts[nCount] =
                                mesh.edges()[nbrEdgeI].centre(mesh.points());
                            nCount++;
                        }
                    }

                    label ownEdge = -1;

                    if (pType4 || bPt)
                    {
                        ownEdge = addedEdges[aEI];
                        own = newEdgeCells[ownEdge].first();
                        nei = -1;
                    }
                    else
                    {
                        forAll(pEdges, pEI)
                        {
                            label edgeI = pEdges[pEI];
                            if (edgeType[edgeI] == 2)
                            {
                                label newCellI = newEdgeCells[edgeI].first();

                                if (newCellI == -1)
                                {
                                    continue;
                                }

                                if (own == -1)
                                {
                                    ownEdge = edgeI;
                                    own = newCellI;
                                    nei = -1;
                                }
                                else
                                {
                                    if (newCellI > own)
                                    {
                                        nei = newCellI;
                                    }
                                    else
                                    {
                                        ownEdge = edgeI;
                                        nei = own;
                                        own = newCellI;
                                    }
                                }
                            }
                        }
                    }

                    if (nei == -1)
                    {
                        const labelList& pFaces = mesh.pointFaces()[pointI];
                        forAll(pFaces, pFI)
                        {
                            label faceI = pFaces[pFI];
                            if (cFaces.found(faceI))
                            {
                                face f = mesh.faces()[faceI];
                                label index = findIndex(f, pointI);
                                label nextPt = f[f.fcIndex(index)];
                                label prevPt = f[f.rcIndex(index)];

                                label nextEdge = meshTools::findEdge
                                (
                                    mesh.edges(),
                                    mesh.pointEdges()[pointI],
                                    pointI,
                                    nextPt
                                );

                                label prevEdge = meshTools::findEdge
                                (
                                    mesh.edges(),
                                    mesh.pointEdges()[pointI],
                                    pointI,
                                    prevPt
                                );

                                if
                                (
                                    edgeType[nextEdge] == -1
                                    && edgeType[prevEdge] == -1
                                )
                                {
                                    patchI = patches.whichPatch(faceI);
                                    origFace = faceI;
                                    break;
                                }
                            }
                        }
                    }

                    edge oEdge = mesh.edges()[ownEdge];

                    label otherPt(oEdge[0] ==  pointI ? oEdge[1] : oEdge[0]);
                    vector eVec = mesh.points()[pointI] - mesh.points()[otherPt];

                    if ((eVec & norm.areaNormal(normPts)) < 0)
                    {
                        flip = true;
                    }
                }

                if (flip)
                {
                    newFace.flip();
                }

                label zoneI = -1;
                bool flipZone = false;
                if (origFace != -1)
                {
                    zoneI = mesh.faceZones().whichZone(origFace);

                    if (zoneI != -1)
                    {
                        const faceZone& fz = mesh.faceZones()[zoneI];
                        flipZone = fz.flipMap()[fz.whichFace(origFace)];
                    }
                }

                //Add top face
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,   // face
                        own,       // owner
                        nei,  // neighbour
                        -1,        // master point
                        -1,        // master edge
                        origFace,     // master face
                        false,     // flux flip
                        patchI,        // patch for face
                        zoneI,     // zone for face
                        flipZone       // face zone flip
                     )
                 );
                nAddedFaces++;
            }
        }
        else if (pointType[pointI] == 4 && newPointCells[pointI].first() != -1)
        {
            label newCellI = newPointCells[pointI].first();
            label origCellI = newPointCells[pointI].second();

            labelHashSet cEdges(cEdgeList[origCellI]);
            labelHashSet cFaces(cells[origCellI]);

            const labelList& pFaces = mesh.pointFaces()[pointI];

            forAll(pFaces, pFI)
            {
                label faceI = pFaces[pFI];
                if (cFaces.found(faceI))
                {
                    face newFace(3);
                    label own = -1;
                    label nei = -1;
                    label patchI = -1;

                    bool flip = false;

                    label nCount = 0;
                    newFace[nCount++] = pointI;
                    label startPt = -1;

                    labelHashSet fEdges(mesh.faceEdges()[faceI]);
                    forAll(convexEdgePts[pointI], cCEPI)
                    {
                        label nbrEdgeI = convexEdgePts[pointI][cCEPI].first();
                        if (fEdges.found(nbrEdgeI))
                        {
                            if (nCount == 1)
                            {
                                edge e = mesh.edges()[nbrEdgeI];

                                startPt = (e[0] ==  pointI ? e[1] : e[0]);
                            }
                            newFace[nCount++] =
                                convexEdgePts[pointI][cCEPI].second();
                        }
                    }

                    bool swap  = false;
                    if (mesh.isInternalFace(faceI))
                    {
                        nei = newCellI;
                        label origNei
                        (
                            mesh.faceOwner()[faceI] == origCellI
                            ? mesh.faceNeighbour()[faceI]
                            : mesh.faceOwner()[faceI]
                        );
                        swap = (origCellI < origNei ? true : false);

                        const labelList& pEdges = mesh.pointEdges()[pointI];
                        forAll(pEdges, pEI)
                        {
                            label edgeI = pEdges[pEI];
                            if
                            (
                                edgeType[edgeI] == 2
                                && newEdgeCells[edgeI].second() == origNei
                            )
                            {
                                own = newEdgeCells[edgeI].first();
                                break;
                            }
                        }
                    }
                    else
                    {
                        patchI = patches.whichPatch(faceI);
                        own = newCellI;
                    }

                    const face f = mesh.faces()[faceI];
                    label cIndex = findIndex(f, pointI);

                    if
                    (
                        (f[f.fcIndex(cIndex)] != startPt && !swap)
                        || (f[f.rcIndex(cIndex)] != startPt && swap)
                    )
                    {
                        flip = true;
                    }

                    if (flip)
                    {
                        newFace.flip();
                    }

                    label zoneI = -1;
                    bool flipZone = false;
                    if (faceI != -1)
                    {
                        zoneI = mesh.faceZones().whichZone(faceI);

                        if (zoneI != -1)
                        {
                            const faceZone& fz = mesh.faceZones()[zoneI];
                            flipZone = fz.flipMap()[fz.whichFace(faceI)];
                        }
                    }

                    //Add top face
                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            newFace,   // face
                            own,       // owner
                            nei,  // neighbour
                            -1,        // master point
                            -1,        // master edge
                            faceI,     // master face
                            false,     // flux flip
                            patchI,        // patch for face
                            zoneI,     // zone for face
                            flipZone       // face zone flip
                         )
                    );
                    nAddedFaces++;
                }
            }
        }
        else if (pointType[pointI] > 4 && newPointCells[pointI].first() != -1)
        {
            label newCellI = newPointCells[pointI].first();
            label origCellI = newPointCells[pointI].second();
            labelHashSet cEdges(cEdgeList[origCellI]);
            labelHashSet cFaces(cells[origCellI]);

            const labelList& pFaces = mesh.pointFaces()[pointI];

            forAll(pFaces, pFI)
            {
                label faceI = pFaces[pFI];

                if (cFaces.found(faceI) && faceType[faceI] == -1)
                {
                    face f = mesh.faces()[faceI];
                    label index = findIndex(f, pointI);
                    label nextPt = f[f.fcIndex(index)];
                    label prevPt = f[f.rcIndex(index)];

                    label nextEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[pointI],
                        pointI,
                        nextPt
                    );

                    label prevEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[pointI],
                        pointI,
                        prevPt
                    );

                    if (edgeType[nextEdge] != -1 && edgeType[prevEdge] != -1)
                    {
                        face newFace(5);
                        label own = -1;
                        label nei = -1;
                        label patchI = -1;

                        bool flip = false;

                        label nCount = 0;
                        newFace[nCount++] = pointI;

                        labelList otherPts(2);
                        otherPts[0] = nextPt;
                        otherPts[1] = prevPt;

                        label startPt = -1;

                        labelHashSet fEdges(mesh.faceEdges()[faceI]);
                        forAll(otherPts, i)
                        {
                            label otherPt = otherPts[i];

                            label newAddedPt = -1;

                            if (addedPoints[otherPt] != -1)
                            {
                                newAddedPt = addedPoints[otherPt];
                            }
                            else if (pointLevel[otherPt] > baffleLevel[otherPt])
                            {
                                forAll(convexFacePts[otherPt], cFPI)
                                {
                                    label nbrFaceI =
                                        convexFacePts[otherPt][cFPI].first();
                                    if (nbrFaceI == faceI)
                                    {
                                        newAddedPt =
                                            convexFacePts[otherPt][cFPI].second();
                                        break;
                                    }
                                }
                            }
                            else
                            {
                                forAll(convexEdgePts[otherPt], cCEPI)
                                {
                                    label nbrEdgeI =
                                        convexEdgePts[otherPt][cCEPI].first();
                                    if (fEdges.found(nbrEdgeI))
                                    {
                                        newAddedPt =
                                            convexEdgePts[otherPt][cCEPI].second();
                                        break;
                                    }
                                }
                            }

                            if (nCount == 1)
                            {
                                startPt = otherPt;
                                newFace[nCount++] = otherPt;
                                newFace[nCount++] = newAddedPt;
                            }
                            else
                            {
                                newFace[nCount++] = newAddedPt;
                                newFace[nCount++] = otherPt;
                            }
                        }

                        bool swap  = false;
                        if (mesh.isInternalFace(faceI))
                        {
                            nei = newCellI;
                            label origNei
                            (
                                mesh.faceOwner()[faceI] == origCellI
                                ? mesh.faceNeighbour()[faceI]
                                : mesh.faceOwner()[faceI]
                            );
                            swap = (origCellI > origNei ? true : false);

                            const labelList& pEdges = mesh.pointEdges()[pointI];
                            forAll(pEdges, pEI)
                            {
                                label edgeI = pEdges[pEI];
                                if
                                (
                                    edgeType[edgeI] == 1
                                    && newEdgeCells[edgeI].second() == origNei
                                )
                                {
                                    own = newEdgeCells[edgeI].first();
                                    break;
                                }
                            }
                        }
                        else
                        {
                            patchI = patches.whichPatch(faceI);
                            own = newCellI;
                        }

                        label cIndex = findIndex(f, pointI);

                        if (patchI != -1)
                        {
                            if (f[f.fcIndex(cIndex)] != startPt)
                            {
                                flip = true;
                            }
                        }
                        else
                        {
                            if
                            (
                                (f[f.rcIndex(cIndex)] != startPt && !swap)
                                || (f[f.fcIndex(cIndex)] != startPt && swap)
                             )
                            {
                                flip = true;
                            }
                        }

                        if (flip)
                        {
                            newFace.flip();
                        }

                        label zoneI = -1;
                        bool flipZone = false;
                        if (faceI != -1)
                        {
                            zoneI = mesh.faceZones().whichZone(faceI);

                            if (zoneI != -1)
                            {
                                const faceZone& fz = mesh.faceZones()[zoneI];
                                flipZone = fz.flipMap()[fz.whichFace(faceI)];
                            }
                        }

                        meshMod.setAction
                        (
                            polyAddFace
                            (
                                newFace,   // face
                                own,       // owner
                                nei,  // neighbour
                                -1,        // master point
                                -1,        // master edge
                                faceI,     // master face
                                false,     // flux flip
                                patchI,        // patch for face
                                zoneI,     // zone for face
                                flipZone       // face zone flip
                             )
                         );
                         nAddedFaces++;
                    }
                }
            }
        }
    }

    return nAddedFaces;
}


label Foam::addHexInternalLayers::modifyOtherFaces
(
    const labelList& pointType,
    const labelList& edgeType,
    const labelList& faceType,
    const labelList& addedPoints,
    const List<DynamicList<labelPair>>& convexEdgePts,
    const List<DynamicList<labelPair>>& convexFacePts,
    const labelList& baffleLevel,
    polyTopoChange& meshMod
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    label nModifiedFaces = 0;

    forAll(mesh.faces(), faceI)
    {
        const face f = mesh.faces()[faceI];

        if (faceType[faceI] != -1)
        {
            continue;
        }

        bool boundPt = false;
        forAll(f, fp)
        {
            if (pointType[f[fp]] != -1)
            {
                boundPt = true;
                break;
            }
        }

        if (!boundPt)
        {
            continue;
        }

        DynamicList<label> newPts(2*f.size());
        labelHashSet fEdges(mesh.faceEdges()[faceI]);

        forAll(f, fp)
        {
            label pointI = f[fp];
            if (pointType[pointI] == -1)
            {
                newPts.append(pointI);
            }
            else
            {
                label prevIndex = f.rcIndex(fp);
                label nextIndex = f.fcIndex(fp);

                label prevPt = f[prevIndex];
                label nextPt = f[nextIndex];

                label nextEdge = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[pointI],
                    pointI,
                    nextPt
                );

                label prevEdge = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[pointI],
                    pointI,
                    prevPt
                );

                if (addedPoints[pointI] != -1)
                {
                   newPts.append(addedPoints[pointI]);
                }
                else if (pointLevel[pointI] > baffleLevel[pointI])
                {
                    forAll(convexFacePts[pointI], cFPI)
                    {
                        if (convexFacePts[pointI][cFPI].first() == faceI)
                        {
                            newPts.append(convexFacePts[pointI][cFPI].second());
                        }
                    }
                }
                else
                {
                    label prevEdgePt = -1;
                    label nextEdgePt = -1;

                    forAll(convexEdgePts[pointI], cCEPI)
                    {
                        label nbrEdgeI = convexEdgePts[pointI][cCEPI].first();
                        if (nbrEdgeI == prevEdge)
                        {
                            prevEdgePt = convexEdgePts[pointI][cCEPI].second();
                        }
                        if (nbrEdgeI == nextEdge)
                        {
                            nextEdgePt = convexEdgePts[pointI][cCEPI].second();
                        }
                    }

                    if (prevEdgePt != -1)
                    {
                        newPts.append(prevEdgePt);
                    }
                    if (nextEdgePt != -1)
                    {
                        newPts.append(nextEdgePt);
                    }
                }
            }
        }

        label patchI = patches.whichPatch(faceI);
        label own = mesh.faceOwner()[faceI];
        label nei(patchI == -1 ? mesh.faceNeighbour()[faceI] : -1);

        label zoneI = mesh.faceZones().whichZone(faceI);
        bool flip = false;
        if (zoneI != -1)
        {
            const faceZone& fz = mesh.faceZones()[zoneI];
            flip = fz.flipMap()[fz.whichFace(faceI)];
        }

        face newFace(newPts.shrink());


        meshMod.setAction
        (
            polyModifyFace
            (
                newFace,          // modified face
                faceI,      // label of face
                own,   // owner
                nei,         // neighbour
                false,      // face flip
                patchI, // patch for face
                false,      // remove from zone
                zoneI,      // zone for face
                flip        // face flip in zone
            )
        );
        nModifiedFaces++;
    }

    return nModifiedFaces;
}


void Foam::addHexInternalLayers::setRefinement()
{
    Info<< nl << "Adding layers to refinement  interfaces"<<  nl <<endl;

    label nBaffles = returnReduce(baffles_.size(), sumOp<label>());

    if (nBaffles == 0)
    {
        Info<< nl << "No refinement interfaces to add layers too"<<endl;

        return;
    }
    else
    {
        Info<< nl << "Adding layers to "<<nBaffles <<" faces"<<endl;
    }

    fvMesh& mesh = meshRefiner_.mesh();
    polyTopoChange meshMod(mesh);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    label patchLower = mesh.boundaryMesh().findPatchID("lowerLevel");

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            labelList(1, patchLower)
        )
    );
    indirectPrimitivePatch& pp = ppPtr();

    // Determine edge normals on original patch
    labelList patchEdges;
    labelList coupledEdges;
    PackedBoolList sameEdgeOrientation;
    PatchTools::matchEdges
    (
        pp,
        mesh.globalData().coupledPatch(),
        patchEdges,
        coupledEdges,
        sameEdgeOrientation
    );

    pointField ppEdgeNormals
    (
        PatchTools::edgeNormals
        (
            mesh,
            pp,
            patchEdges,
            coupledEdges
        )
    );

    labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    vectorField meshEdgeNormals(mesh.nEdges(), vector::zero);
    forAll(ppEdgeNormals, edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        meshEdgeNormals[meshEdgeI] = ppEdgeNormals[edgeI];
    }
    syncTools::syncEdgeList
    (
        mesh,
        meshEdgeNormals,
        maxMagSqrEqOp<point>(),
        vector::zero // null value (note: cannot use VGREAT)
    );

    pointField ppNormals(PatchTools::pointNormals(mesh, pp));
    vectorField meshNormals(mesh.nPoints(), vector::zero);
    forAll(ppNormals, patchPointI)
    {
        label meshPointI = pp.meshPoints()[patchPointI];
        meshNormals[meshPointI] = ppNormals[patchPointI];
    }
    syncTools::syncPointList
    (
        mesh,
        meshNormals,
        maxMagSqrEqOp<point>(),
        vector::zero // null value (note: cannot use VGREAT)
    );

    if (debug)
    {
        simpleVTKWriter
        (
            pp.localFaces(),
            pp.localPoints()
         ).write("ppRefinementInterface.vtk");
    }

    // Point status:
    //  -1 : is internal
    //  0 : is non feature boundary point
    //  1 : is concave point
    //  2 : is concave corner point
    //  3 : is convex point
    //  4 : is convex corner point (3 convex edges)
    //  5 : is concave/convex corner point (2 concave edges and 1 convex edge)
    //  6 : is concave/convex corner point (1 concave edges and 2 convex edge)
    //  7 : is concave/convex corner point (3 concave edges and 3 convex edge)
    labelList pointType(mesh.nPoints(), -1);

    // Edge status:
    //  -1 : is internal
    //  0 : is non feature boundary edge
    //  1 : is concave edge
    //  2 : is convex edge
    labelList edgeType(mesh.nEdges(), -1);

    // Face status:
    //  -1 : no added layers
    //  0 : is boundary face
    labelList faceType(mesh.nFaces(), -1);

    DynamicList<label> boundaryCornerPoints(mesh.nPoints()/100);
    calculateTypes
    (
        pp,
        meshEdgeNormals,
        meshEdges,
        pointType,
        edgeType,
        faceType,
        boundaryCornerPoints
    );

    //Add patch face cells
    List<labelPair> newFaceCells(mesh.nFaces(), labelPair(-1, -1));

    //Add convex feature edge cells (first: new cell ID second: original cell ID)
    List<labelPair> newEdgeCells(mesh.nEdges(), labelPair(-1, -1));

    //Add convex and concave corner cells
    List<labelPair> newPointCells(mesh.nPoints(), labelPair(-1, -1));

    addCells
    (
        pp,
        pointType,
        edgeType,
        faceType,
        newFaceCells,
        newEdgeCells,
        newPointCells,
        meshMod
    );

    labelList addedPoints(mesh.nPoints(),-1);
    List<DynamicList<labelPair>>
        convexEdgePts(mesh.nPoints(), DynamicList<labelPair>(3));
    List<DynamicList<labelPair>>
        convexFacePts(mesh.nPoints(), DynamicList<labelPair>(3));

    labelList baffleLevel(mesh.nPoints(),-1);
    forAll(pp.meshPoints(), i)
    {
        label pointI = pp.meshPoints()[i];
        const labelList& pCells = mesh.pointCells()[pointI];

        baffleLevel[pointI] = max(baffleLevel[pointI], cellLevel[pCells[0]]);
    }

    syncTools::syncPointList
    (
        mesh,
        baffleLevel,
        maxEqOp<label>(),
        labelMin
    );

    boolList procEdges(mesh.nEdges(), false);
    forAll(meshEdges, edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];
        forAll(eFaces, eFI)
        {
            label faceI = eFaces[eFI];
            label patchI = patches.whichPatch(faceI);
            if (patchI != -1 && patches[patchI].coupled())
            {
                procEdges[meshEdgeI] = true;
                break;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh,
        procEdges,
        orEqOp<bool>(),
        false
    );

    addPoints
    (
        meshNormals,
        pointType,
        edgeType,
        baffleLevel,
        addedPoints,
        convexEdgePts,
        convexFacePts,
        meshMod
    );

    label nModifiedFaces = modifyBoundaryFaces
    (
        pp,
        pointType,
        edgeType,
        faceType,
        newFaceCells,
        newEdgeCells,
        newPointCells,
        meshMod
    );

    if (nModifiedFaces != pp.size())
    {
        FatalErrorInFunction
            <<"Failed to update all boundary faces. "
            <<" Modified : " << nModifiedFaces
            <<" Boundary faces: "<<pp.size()
            << exit(FatalError);
    }

    label nAddedFaces = addTopFaces
    (
        pp,
        pointType,
        edgeType,
        faceType,
        newFaceCells,
        newEdgeCells,
        newPointCells,
        addedPoints,
        convexEdgePts,
        convexFacePts,
        baffleLevel,
        meshMod
    );

    nAddedFaces += addSideFaces
    (
        pp,
        pointType,
        edgeType,
        faceType,
        newFaceCells,
        newEdgeCells,
        newPointCells,
        addedPoints,
        convexEdgePts,
        convexFacePts,
        baffleLevel,
        meshEdges,
        meshMod
    );

    nModifiedFaces += modifyOtherFaces
    (
        pointType,
        edgeType,
        faceType,
        addedPoints,
        convexEdgePts,
        convexFacePts,
        baffleLevel,
        meshMod
    );

    Info<<"Number of new faces: "
        <<returnReduce(nAddedFaces, sumOp<label>())<<endl;
    Info<<"Number of modified faces: "
        <<returnReduce(nModifiedFaces, sumOp<label>())<<endl;

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

    // Update and merge baffles
    forAll(baffles_, i)
    {
        label& face0 = baffles_[i].first();
        label& face1 = baffles_[i].second();
        face0 = map().reverseFaceMap()[face0];
        face1 = map().reverseFaceMap()[face1];
    }

    mergeRefinementBaffles();

    Info<<"Finished refining interface cells"<<  nl <<endl;

    return;
}



// ************************************************************************* //
