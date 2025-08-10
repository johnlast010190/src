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

#include "addHexMeshLayers/addHexMeshLayer.H"
#include "meshTools/meshTools.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "meshCut/cellCuts/cellCuts.H"
#include "meshCut/meshModifiers/meshCutter/meshCutter.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "sets/topoSets/cellSet.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "meshes/polyMesh/syncTools/dummyTransform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(addHexMeshLayer, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::addHexMeshLayer::addHexMeshLayer
(
    const labelList& globalToMasterPatch,
    const refinementParameters& refineParams,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    meshRefinement& meshRefiner,
    indirectPrimitivePatch& pp,
    const scalar cutRatio,
    const bool geometricChecks,
    const labelList grownPatchIDs,
    const bool trackLayerCells
)
:
    globalToMasterPatch_(globalToMasterPatch),
    refineParams_(refineParams),
    decomposer_(decomposer),
    distributor_(distributor),
    meshRefiner_(meshRefiner),
    pp_(pp),
    cutRatio_(cutRatio),
    geometricChecks_(geometricChecks),
    grownPatches_(grownPatchIDs),
    trackLayerCells_(trackLayerCells)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::addHexMeshLayer::calculateBoundary
(
    boolList& boundaryPts,
    boolList& boundaryEdges,
    boolList& boundaryFaces
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    boundaryPts.setSize(mesh.nPoints(), false);
    boundaryEdges.setSize(mesh.nEdges(), false);
    boundaryFaces.setSize(mesh.nFaces(), false);

    forAll(pp_.meshPoints(), ptI)
    {
        label meshPointI = pp_.meshPoints()[ptI];

        boundaryPts[meshPointI] = true;
    }

    syncTools::syncPointList
    (
        mesh,
        boundaryPts,
        orEqOp<bool>(),
        false
    );

    labelList meshEdges(pp_.meshEdges(mesh.edges(), mesh.pointEdges()));
    forAll(meshEdges, edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        boundaryEdges[meshEdgeI] = true;
    }

    if (grownPatches_.size())
    {
        autoPtr<indirectPrimitivePatch> grownPPPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                grownPatches_
             )
        );
        indirectPrimitivePatch& grownPP = grownPPPtr();

        labelList grownMeshEdges
        (
            grownPP.meshEdges(mesh.edges(), mesh.pointEdges())
        );
        forAll(grownMeshEdges, edgeI)
        {
            label meshEdgeI = grownMeshEdges[edgeI];
            boundaryEdges[meshEdgeI] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    forAll(pp_, faceI)
    {
        label meshFaceI = pp_.addressing()[faceI];
        boundaryFaces[meshFaceI] = true;
    }

    syncTools::syncFaceList
    (
        mesh,
        boundaryFaces,
        orEqOp<bool>()
    );
}


bool Foam::addHexMeshLayer::checkLoop
(
    const fvMesh& mesh,
    const List<labelList>& cutLoops,
    const labelList& edgeLoop,
    const boolList& boundaryPts,
    const label cellI
)
{
    const cell c = mesh.cells()[cellI];

    if (edgeLoop.size() <= 2)
    {
        return false;
    }

    forAll(c, cFI)
    {
        label faceI = c[cFI];
        labelHashSet fPts(mesh.faces()[faceI]);

        label nFound = 0;

        forAll(edgeLoop, eLI)
        {
            if (fPts.found(edgeLoop[eLI]))
            {
                nFound++;
            }
        }

        if (nFound > 2)
        {
            return false;
        }
    }
/*
    forAll(cutLoops, loopI)
    {
        const labelList& cLoop = cutLoops[loopI];
        labelHashSet cLoopSet(cLoop);
        bool foundAll = true;
        label nShared = 0;

        forAll(edgeLoop, eLI)
        {
            if (!cLoopSet.found(edgeLoop[eLI]))
            {
                foundAll = false;
            }
            else
            {
                nShared++;
            }
        }


        if (foundAll || nShared > 2)
        {
            return false;
        }
    }
*/
    //Check if duplicate loop exists
    forAll(cutLoops, loopI)
    {
        const labelList& cLoop = cutLoops[loopI];
        labelHashSet cLoopSet(cLoop);

        forAll(edgeLoop, eLI)
        {
            if (cLoopSet.found(edgeLoop[eLI]))
            {
                return false;
            }
        }
    }

    if (geometricChecks_)
    {
        labelHashSet cellEdges(mesh.cellEdges()[cellI]);
        label nPerp = 0;

        forAll(edgeLoop, eLI)
        {
            label start = -1;
            label end = -1;
            label mid = edgeLoop[eLI];
            if (eLI == 0)
            {
                start = edgeLoop[edgeLoop.size()-1];
                end = edgeLoop[eLI+1];
            }
            else if (eLI == edgeLoop.size()-1)
            {
                start = edgeLoop[eLI-1];
                end = edgeLoop[0];
            }
            else
            {
                start = edgeLoop[eLI-1];
                end = edgeLoop[eLI+1];
            }

            //check neighbouring cut edges
            vector eDirStart = vector::zero;
            vector eDirNext = vector::zero;

            const labelList& pEdgesStart = mesh.pointEdges()[start];

            forAll(pEdgesStart, pEI)
            {
                label edgeI = pEdgesStart[pEI];
                edge e = mesh.edges()[edgeI];
                label otherPt = (e[0] ==  start ? e[1] : e[0]);
                if (boundaryPts[otherPt])
                {
                    eDirStart = mesh.points()[otherPt] - mesh.points()[start];
                    eDirStart /= (mag(eDirStart) + SMALL);
                    break;
                }
            }

            const labelList& pEdgesNext = mesh.pointEdges()[mid];

            forAll(pEdgesNext, pEI)
            {
                label edgeI = pEdgesNext[pEI];
                edge e = mesh.edges()[edgeI];
                label otherPt = (e[0] ==  mid ? e[1] : e[0]);
                if (boundaryPts[otherPt])
                {
                    eDirNext = mesh.points()[otherPt] - mesh.points()[mid];
                    eDirNext /= (mag(eDirNext) + SMALL);
                    break;
                }
            }

            if ((eDirStart & eDirNext) < -0.99)
            {
                return false;
            }

            vector e1 = mesh.points()[end] - mesh.points()[mid];
            vector e2 = mesh.points()[start] - mesh.points()[mid];

            vector e3 = (e2^e1);
            scalar magE3 = mag(e3);

            if (magE3 > SMALL)
            {
                e3 /= magE3;

                const labelList& pEdges = mesh.pointEdges()[mid];

                vector pVec = vector::zero;

                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cellEdges.found(edgeI))
                    {
                        pVec = mesh.edges()[edgeI].vec(mesh.points());
                        pVec /= (mag(pVec) + SMALL);
                        break;
                    }
                }
                if (mag(pVec & e3)<0.1)
                {
                    nPerp++;
                }
            }
        }

        //geometric checks on cut face
        if (nPerp > 1)
        {
            face newFace(edgeLoop);
            point newFC = newFace.centre(mesh.points());

            scalar minDisp = GREAT;
            scalar maxDisp = -GREAT;

            forAll(c, cFI)
            {
                label faceI = c[cFI];
                plane fPlane(mesh.faceCentres()[faceI], mesh.faceAreas()[faceI]);

                scalar disp = mag(fPlane.nearestPoint(newFC)
                                  - newFC);

                if (disp < minDisp)
                {
                    minDisp =  disp;
                }

                if (disp > maxDisp)
                {
                    maxDisp =  disp;
                }
            }

            scalar ratio = minDisp/maxDisp;

            if (ratio < 0.085)
            {
                return false;
            }
        }
    }

    return true;
}


void Foam::addHexMeshLayer::balance()
{
    if (Pstream::nProcs() > 1)
    {
        fvMesh& mesh = meshRefiner_.mesh();

        scalar nIdealCells =
            mesh.globalData().nTotalCells()
          / Pstream::nProcs();

        scalar unbalance = returnReduce
        (
            scalar(mag(1.0-mesh.nCells()/nIdealCells)),
            maxOp<scalar>()
        );

        if (unbalance >= refineParams_.maxLoadUnbalance())
        {
            Info<<"Performing Re-balancing of mesh "<<endl;

            boolList ppFaceAddressing(mesh.nFaces(), false);
            forAll(pp_, i)
            {
                ppFaceAddressing[pp_.addressing()[i]] = true;
            }

            scalarField weights(mesh.nCells(), 1.0);
            boolList boundaryPts,boundaryEdges,boundaryFaces;
            calculateBoundary(boundaryPts, boundaryEdges, boundaryFaces);

            forAll(mesh.cells(), cellI)
            {
                const labelList& cPts = mesh.cellPoints()[cellI];
                forAll(cPts, cPtI)
                {
                    if (boundaryPts[cPts[cPtI]])
                    {
                        weights[cellI] += 1;
                        break;
                    }
                }
            }

            autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
            (
                false,
                false,
                false,
                weights,
                decomposer_,
                distributor_,
                false
            );

            //re-distrubute pp_
            map().distributeFaceData(ppFaceAddressing);

            DynamicList<label> newAddressing(ppFaceAddressing.size()/10);
            forAll(ppFaceAddressing, faceI)
            {
                if (ppFaceAddressing[faceI])
                {
                    newAddressing.append(faceI);
                }
            }

            pp_.clearOut();
            pp_.resetAddressing(newAddressing);
            pp_.localPoints();
        }
    }
}


Foam::labelList Foam::addHexMeshLayer::createCutLoop
(
    const fvMesh& mesh,
    const boolList& cutPts,
    const Map<label>& cellPts,
    const label startPt,
    const label startEdge,
    const label startFace,
    const labelHashSet& cellFaces,
    boolList& anchors
)
{
    DynamicList<label> edgeLoop(anchors.size());

    edge e = mesh.edges()[startEdge];
    label otherPt(e[0] == startPt ? e[1] : e[0]);
    label endPt = otherPt;
    edgeLoop.append(otherPt);

    label nextPt = startPt;
    label nextEdge = startEdge;
    label nextFace = startFace;

    while (true)
    {
        const labelList& eFaces = mesh.edgeFaces()[nextEdge];
        forAll(eFaces, eFI)
        {
            label faceI = eFaces[eFI];
            if (faceI != nextFace && cellFaces.found(faceI))
            {
                nextFace = faceI;
                break;
            }
        }

        face f = mesh.faces()[nextFace];
        label start = findIndex(f, otherPt);
        label fPt = f[f.fcIndex(start)];

        if (nextPt != fPt)
        {
            f.flip();
        }

        label fp = findIndex(f, nextPt);
        label fn = f.fcIndex(fp);

        forAll(f, fI)
        {
            label mEdge = meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[f[fp]],
                f[fp],
                f[fn]
             );

            anchors[cellPts[f[fp]]] = true;
            anchors[cellPts[f[fn]]] = true;

            if (cutPts[cellPts[f[fn]]])
            {
                nextEdge = mEdge;
                otherPt = f[fn];
                nextPt = f[fp];
                break;
            }

            fp = fn;
            fn = f.fcIndex(fp);
        }

        if (otherPt != endPt)
        {
            edgeLoop.append(otherPt);
        }
        else
        {
            break;
        }
    }

    return edgeLoop.shrink();
}


Foam::labelList Foam::addHexMeshLayer::calculateCutEdges()
{
    fvMesh& mesh = meshRefiner_.mesh();

    boolList boundaryPts,boundaryEdges,boundaryFaces;
    calculateBoundary(boundaryPts, boundaryEdges, boundaryFaces);

    // Edge status:
    //  >0 : label of number of cuts (1 or 2)
    //  -1 : no cuts
    labelList cutEdges(mesh.nEdges(), 0);

    forAll(mesh.edges(), edgeI)
    {
        if (!boundaryEdges[edgeI])
        {
            const edge& e = mesh.edges()[edgeI];
            if (boundaryPts[e[0]] || boundaryPts[e[1]])
            {
                if (boundaryPts[e[0]])
                {
                    cutEdges[edgeI]++;
                }

                if (boundaryPts[e[1]])
                {
                    cutEdges[edgeI]++;
                }
            }
            else
            {
                cutEdges[edgeI] = -1;
            }
        }
        else
        {
            cutEdges[edgeI] = -1;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        cutEdges,
        maxEqOp<label>(),
        labelMin
    );

    return cutEdges;
}



void Foam::addHexMeshLayer::removeProblemCells()
{
    fvMesh& mesh = meshRefiner_.mesh();
    const labelListList& cEdgeList = mesh.cellEdges();

    labelList cutEdges = calculateCutEdges();

    boolList boundaryPts(mesh.nPoints(), false);

    forAll(pp_.meshPoints(), ptI)
    {
        label meshPointI = pp_.meshPoints()[ptI];

        boundaryPts[meshPointI] = true;
    }

    syncTools::syncPointList
    (
        mesh,
        boundaryPts,
        orEqOp<bool>(),
        false
    );

    label nChecks = 0;
    forAll(mesh.cells(), cellI)
    {
        const labelList& cellEdges = cEdgeList[cellI];

        forAll(cellEdges, cEI)
        {
            label edgeI = cellEdges[cEI];
            if (cutEdges[edgeI] != -1)
            {
                if (cutEdges[edgeI] == 1)
                {
                    nChecks++;
                }
                else
                {
                    nChecks += 2;
                }
            }
        }
    }

    labelList checkedCells(nChecks);
    pointField start(nChecks);
    pointField end(nChecks);

    nChecks = 0;
    forAll(mesh.cells(), cellI)
    {
        const labelList& cellEdges = cEdgeList[cellI];
        const point cc = mesh.cellCentres()[cellI];

        forAll(cellEdges, cEI)
        {
            label edgeI = cellEdges[cEI];
            if (cutEdges[edgeI] != -1)
            {
                edge e = mesh.edges()[edgeI];
                vector eVec = e.vec(mesh.points());

                if (cutEdges[edgeI] == 1)
                {
                    checkedCells[nChecks] = cellI;
                    start[nChecks] = cc;

                    vector midPoint = e[0];

                    if (boundaryPts[e[0]])
                    {
                        midPoint = mesh.points()[e[0]]
                            + cutRatio_*eVec;
                    }
                    else
                    {
                        midPoint = mesh.points()[e[0]]
                            + (1.-cutRatio_)*eVec;
                    }

                    end[nChecks] = midPoint;
                    nChecks++;
                }
                else
                {
                    checkedCells[nChecks] = cellI;
                    start[nChecks] = cc;
                    vector cutPoint = mesh.points()[e[0]]
                        + 0.33*eVec;
                    end[nChecks] = cutPoint;
                    nChecks++;

                    checkedCells[nChecks] = cellI;
                    start[nChecks] = cc;
                    cutPoint = mesh.points()[e[0]]
                        + 0.66*eVec;
                    end[nChecks] = cutPoint;
                    nChecks++;
                }
            }
        }
    }

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;

    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    // Surfaces that need to be baffled
    const labelList surfacesToBaffle
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfaces.surfZones())
    );

    surfaces.findNearestIntersection
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

    // Get cells to remove
    label foundCell = -1;

    labelList ownPatch(mesh.nFaces(), 0);
    boolList blockedFace(mesh.nFaces(),false);
    boolList problemCells(mesh.nCells(), false);
    forAll(start, i)
    {
        label cellI = checkedCells[i];
        if (hit1[i].hit() && cellI != foundCell)
        {
            label patchI = globalToMasterPatch_
            [
                surfaces.globalRegion(surface1[i], region1[i])
            ];

            const cell& c = mesh.cells()[cellI];
            forAll(c, cFI)
            {
                label faceI = c[cFI];
                ownPatch[faceI] = patchI;
                blockedFace[faceI] = true;
            }

            problemCells[cellI] = true;
            foundCell = cellI;
        }
    }

    // Set region per cell based on walking
    regionSplit cellRegion(mesh, blockedFace);

    DynamicList<label> keepRegions(cellRegion.nRegions());
    const pointField& locationsInMesh =
        refineParams_.locationsInMesh();

    forAll(locationsInMesh, locationI)
    {
        point locationInMesh = locationsInMesh[locationI];

        // Find the region containing the keepPoint
        label keepRegionI = -1;

        label cellI = meshRefinement::findCell
        (
            locationInMesh,
            mesh,
            meshRefiner_.meshCutter()
        );

        if (cellI != -1)
        {
            keepRegionI = cellRegion[cellI];
        }
        reduce(keepRegionI, maxOp<label>());


        // Find the region containing the keepPoint
        keepRegions.append(keepRegionI);
    }
    labelHashSet keepRegionsSet(keepRegions.shrink());

    DynamicList<label> cellsToRemove(mesh.nCells()/100);
    forAll(mesh.cells(), cellI)
    {
        if (problemCells[cellI] || !keepRegionsSet.found(cellRegion[cellI]))
        {
            cellsToRemove.append(cellI);
        }
    }
    cellsToRemove.shrink();

    // Remove cells
    removeCells cellRemover(mesh);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatchIDs(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label faceI = exposedFaces[i];

        exposedPatchIDs[i] = ownPatch[faceI];
    }

    autoPtr<mapPolyMesh> map = meshRefiner_.doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        cellRemover
    );

    updateMesh(map);
}


void Foam::addHexMeshLayer::addSplitEdgePoints()
{
    fvMesh& mesh = meshRefiner_.mesh();

    boolList boundaryPts,boundaryEdges,boundaryFaces;
    calculateBoundary(boundaryPts, boundaryEdges, boundaryFaces);

    // Edge status:
    //  >0 : label of number of cuts (1 or 2)
    //  -1 : no cuts
    labelList cutEdges(mesh.nEdges(), 0);

    forAll(mesh.edges(), edgeI)
    {
        if (!boundaryEdges[edgeI])
        {
            const edge& e = mesh.edges()[edgeI];
            if (boundaryPts[e[0]] || boundaryPts[e[1]])
            {
                if (boundaryPts[e[0]])
                {
                    cutEdges[edgeI]++;
                }

                if (boundaryPts[e[1]])
                {
                    cutEdges[edgeI]++;
                }
            }
            else
            {
                cutEdges[edgeI] = -1;
            }
        }
        else
        {
            cutEdges[edgeI] = -1;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        cutEdges,
        maxEqOp<label>(),
        labelMin
    );

    polyTopoChange meshMod(mesh);

    List<labelPair> newEdgePts(mesh.nEdges(), labelPair());

    boolList affectedFaces(mesh.nFaces(), false);
    forAll(mesh.edges(), edgeI)
    {
        if (cutEdges[edgeI] != -1)
        {
            const edge e = mesh.edges()[edgeI];

            const labelList& eFaces = mesh.edgeFaces()[edgeI];
            forAll(eFaces, eFI)
            {
                affectedFaces[eFaces[eFI]] = true;
            }

            if (cutEdges[edgeI] == 1)
            {
                vector eVec = e.vec(mesh.points());
                vector cutPt = e[0];

                if (boundaryPts[e[0]])
                {
                    cutPt = mesh.points()[e[0]]
                            + cutRatio_*eVec;
                }
                else
                {
                    cutPt = mesh.points()[e[0]]
                        + (1.-cutRatio_)*eVec;
                }

                newEdgePts[edgeI][0] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        cutPt,         // point
                        e[0],         // master point
                        -1,      // zone for point
                        true        // supports a cell
                     )
                 );
            }
            else
            {
                vector eVec = e.vec(mesh.points());
                point cutPt = mesh.points()[e[0]]+(1./3.)*eVec;

                newEdgePts[edgeI][0] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        cutPt,         // point
                        e[0],         // master point
                        -1,      // zone for point
                        true        // supports a cell
                     )
                );
                cutPt = mesh.points()[e[0]]+(2./3.)*eVec;
                newEdgePts[edgeI][1] = meshMod.setAction
                (
                    polyAddPoint
                    (
                        cutPt,         // point
                        e[1],         // master point
                        -1,      // zone for point
                        true        // supports a cell
                     )
                 );
            }
        }
    }

    forAll(mesh.faces(), faceI)
    {
        if (affectedFaces[faceI])
        {
            face f = mesh.faces()[faceI];

            DynamicList<label> newFace(f.size()+4);
            label prevFp = f[0];
            forAll(f, fp)
            {
                newFace.append(prevFp);

                label nextFp = f[f.fcIndex(fp)];
                label edgeI = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[prevFp],
                    prevFp,
                    nextFp
                );

                if (cutEdges[edgeI] != -1)
                {
                    if (cutEdges[edgeI] == 1)
                    {
                        newFace.append(newEdgePts[edgeI][0]);
                    }
                    else
                    {
                        edge e = mesh.edges()[edgeI];

                        if (prevFp == e[0])
                        {
                            newFace.append(newEdgePts[edgeI][0]);
                            newFace.append(newEdgePts[edgeI][1]);
                        }
                        else
                        {
                            newFace.append(newEdgePts[edgeI][1]);
                            newFace.append(newEdgePts[edgeI][0]);
                        }
                    }
                }
                prevFp = nextFp;
            }
            newFace.shrink();

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
                    face(newFace),              // modified face
                    faceI,                      // label of face
                    mesh.faceOwner()[faceI],    // owner
                    nei,                        // neighbour
                    false,                      // face flip
                    patchID,                    // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                 )
             );
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
    updateMesh(map);

    boolList isAddedPt(mesh.nPoints(), false);

    forAll(map().pointMap(), pointI)
    {
        if (map().pointMap()[pointI] == -1)
        {
            isAddedPt[pointI] = true;
        }
    }
}


void Foam::addHexMeshLayer::calculateSplitLoops
(
    DynamicList<DynamicList<label>>& allCutCells,
    DynamicList<DynamicList<labelList>>& allCellLoops
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    boolList boundaryPts,boundaryEdges,boundaryFaces;
    calculateBoundary(boundaryPts, boundaryEdges, boundaryFaces);

    // Edge status:
    //  >0 : label of number of cuts (1 or 2)
    //  -1 : no cuts
    labelList cutEdges(mesh.nEdges(), 0);

    forAll(mesh.edges(), edgeI)
    {
        if (!boundaryEdges[edgeI])
        {
            const edge& e = mesh.edges()[edgeI];
            if (boundaryPts[e[0]] || boundaryPts[e[1]])
            {
                if (boundaryPts[e[0]])
                {
                    cutEdges[edgeI]++;
                }

                if (boundaryPts[e[1]])
                {
                    cutEdges[edgeI]++;
                }
            }
            else
            {
                cutEdges[edgeI] = -1;
            }
        }
        else
        {
            cutEdges[edgeI] = -1;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        cutEdges,
        maxEqOp<label>(),
        labelMin
    );

    const labelListList& cEdgeList = mesh.cellEdges();
    const labelListList& cPointList = mesh.cellPoints();

    label maxCuts = labelMin;

    List<DynamicList<labelList>>
        cutLoops(mesh.nCells(), DynamicList<labelList>());

    forAll(mesh.cells(), cellI)
    {
        //check faces first
        const cell& cFaces = mesh.cells()[cellI];

        const labelList& cEdges = cEdgeList[cellI];
        const labelList& cPoints = cPointList[cellI];

        labelHashSet cellFaces(cFaces);
        labelHashSet cellEdges(cEdges);
        labelHashSet cellPoints(cPoints);
        label nCuts = 0;

        Map<label> cellPts(cPoints.size());

        DynamicList<label> alreadyCut(cPoints.size());
        forAll(cPoints, cPtI)
        {
            label meshPointI = cPoints[cPtI];
            cellPts.insert(meshPointI,cPtI);
        }

        labelList nConnections(cPoints.size(), -1);
//        label nBoundaryPts = 0;

        forAll(cPoints, cPtI)
        {
            label meshPointI = cPoints[cPtI];
            if (boundaryPts[meshPointI])
            {
//                nBoundaryPts++;
                nConnections[cPtI] = 0;
                const labelList& pEdges = mesh.pointEdges()[meshPointI];
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cellEdges.found(edgeI))
                    {
                        if (boundaryEdges[edgeI])
                        {
                            nConnections[cPtI]++;
                        }
                    }
                }
            }
        }

        boolList cutPts(cPoints.size(), false);
        forAll(cPoints, cPtI)
        {
            label meshPointI = cPoints[cPtI];
            if (boundaryPts[meshPointI])
            {
                const labelList& pEdges = mesh.pointEdges()[meshPointI];
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cellEdges.found(edgeI) && !boundaryEdges[edgeI])
                    {
                        edge e = mesh.edges()[edgeI];
                        label otherPt(e[0] == meshPointI ? e[1] : e[0]);
                        cutPts[cellPts[otherPt]] = true;
                    }
                }
            }
        }

        boolList anchors(cPoints.size(),false);
        //corner points
        forAll(cPoints, cPtI)
        {
            if (anchors[cPtI])
            {
                continue;
            }
            label startPt = cPoints[cPtI];
            label startEdge = -1;
            label startFace = -1;

            if
            (
                boundaryPts[startPt] && nConnections[cPtI] == 0
            )
            {
                const labelList& pEdges = mesh.pointEdges()[startPt];
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cellEdges.found(edgeI) && cutEdges[edgeI] != -1)
                    {
                        startEdge = edgeI;
                        break;
                    }
                }
                const labelList& eFaces = mesh.edgeFaces()[startEdge];
                forAll(eFaces, eFI)
                {
                    label eFaceI = eFaces[eFI];

                    if (cellFaces.found(eFaceI))
                    {
                        startFace = eFaceI;
                        break;
                    }
                }
            }
            if (startEdge != -1)
            {
                boolList anchorsOrig = anchors;

                labelList edgeLoop = createCutLoop
                (
                    mesh,
                    cutPts,
                    cellPts,
                    startPt,
                    startEdge,
                    startFace,
                    cellFaces,
                    anchors
                );

                if
                (
                    checkLoop
                    (
                        mesh,
                        cutLoops[cellI],
                        edgeLoop,
                        boundaryPts,
                        cellI
                    )
                )
                {
                    cutLoops[cellI].append(edgeLoop);
                    nCuts++;
                }
                else
                {
                    anchors = anchorsOrig;
                }
            }
        }

        //boundary edges
        forAll(cPoints, cPtI)
        {
            if (anchors[cPtI])
            {
                continue;
            }
            label startPt = cPoints[cPtI];
            label startEdge = -1;
            label startFace = -1;
            bool valid = true;

            if (nConnections[cPtI] == 1)
            {
                const labelList& pEdges = mesh.pointEdges()[startPt];
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cellEdges.found(edgeI))
                    {
                        if (cutEdges[edgeI] != -1)
                        {
                            startEdge = edgeI;
                        }
                        if (boundaryEdges[edgeI])
                        {
                            edge e = mesh.edges()[edgeI];
                            label otherPt(e[0] == startPt ? e[1] : e[0]);
                            if (nConnections[cellPts[otherPt]] != 1)
                            {
                                valid = false;
                            }
                        }
                    }
                }
                const labelList& eFaces = mesh.edgeFaces()[startEdge];
                forAll(eFaces, eFI)
                {
                    label eFaceI = eFaces[eFI];

                    if (cellFaces.found(eFaceI))
                    {
                        startFace = eFaceI;
                        break;
                    }
                }
            }

            if (startEdge != -1 && valid)
            {
                boolList anchorsOrig = anchors;

                labelList edgeLoop = createCutLoop
                (
                    mesh,
                    cutPts,
                    cellPts,
                    startPt,
                    startEdge,
                    startFace,
                    cellFaces,
                    anchors
                );

                if
                (
                    checkLoop
                    (
                        mesh,
                        cutLoops[cellI],
                        edgeLoop,
                        boundaryPts,
                        cellI
                    )
                )
                {
                    cutLoops[cellI].append(edgeLoop);
                    nCuts++;
                }
                else
                {
                    anchors = anchorsOrig;
                }
            }
        }

        forAll(cPoints, cPtI)
        {
            if (anchors[cPtI])
            {
                continue;
            }

            label startPt = cPoints[cPtI];
            label startEdge = -1;
            label startFace = -1;

            if
            (
                boundaryPts[startPt] && nConnections[cellPts[startPt]] == 2
            )
            {
                const labelList& pEdges = mesh.pointEdges()[startPt];
                label nSingleEdges = 0;

                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cellEdges.found(edgeI))
                    {
                        edge e = mesh.edges()[edgeI];
                        label otherPt(e[0] == startPt ? e[1] : e[0]);
                        if (boundaryEdges[edgeI])
                        {
                            if (nConnections[cellPts[otherPt]] == 1)
                            {
                                nSingleEdges++;
                            }
                        }
                        else if (cutEdges[edgeI] != -1)
                        {
                            startEdge = edgeI;
                        }
                    }
                }

                if (startEdge != -1 && nSingleEdges == 2)
                {
                    const labelList& eFaces = mesh.edgeFaces()[startEdge];
                    forAll(eFaces, eFI)
                    {
                        label eFaceI = eFaces[eFI];

                        if (cellFaces.found(eFaceI))
                        {
                            startFace = eFaceI;
                            break;
                        }
                    }

                    boolList anchorsOrig = anchors;

                    labelList edgeLoop = createCutLoop
                    (
                        mesh,
                        cutPts,
                        cellPts,
                        startPt,
                        startEdge,
                        startFace,
                        cellFaces,
                        anchors
                     );

                    if
                    (
                        checkLoop
                        (
                            mesh,
                            cutLoops[cellI],
                            edgeLoop,
                            boundaryPts,
                            cellI
                         )
                    )
                    {
                        cutLoops[cellI].append(edgeLoop);
                        nCuts++;
                    }
                    else
                    {
                        anchors = anchorsOrig;
                    }
                }
            }
        }

        forAll(cPoints, cPtI)
        {
            if (anchors[cPtI])
            {
                continue;
            }

            label startPt = cPoints[cPtI];
            label startEdge = -1;
            label startFace = -1;

            if
            (
                boundaryPts[startPt] && nConnections[cellPts[startPt]] == 3
            )
            {
                label nSingleEdges = 0;
                label otherPtI = -1;

                {
                    const labelList& pEdges = mesh.pointEdges()[startPt];

                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];
                        if (cellEdges.found(edgeI))
                        {
                            edge e = mesh.edges()[edgeI];
                            label otherPt(e[0] == startPt ? e[1] : e[0]);
                            if (boundaryEdges[edgeI])
                            {
                                if (nConnections[cellPts[otherPt]] == 1)
                                {
                                    nSingleEdges++;
                                    otherPtI = otherPt;
                                }
                            }
                        }
                    }
                }

                if (nSingleEdges == 3)
                {
                    const labelList& pEdges = mesh.pointEdges()[otherPtI];
                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];
                        if (cellEdges.found(edgeI) && cutEdges[edgeI] != -1)
                        {
                            startPt = otherPtI;
                            startEdge = edgeI;
                            break;
                        }
                    }
                }

                if (startEdge != -1)
                {
                    const labelList& eFaces = mesh.edgeFaces()[startEdge];
                    forAll(eFaces, eFI)
                    {
                        label eFaceI = eFaces[eFI];

                        if (cellFaces.found(eFaceI))
                        {
                            startFace = eFaceI;
                            break;
                        }
                    }


                    boolList anchorsOrig = anchors;

                    labelList edgeLoop = createCutLoop
                    (
                        mesh,
                        cutPts,
                        cellPts,
                        startPt,
                        startEdge,
                        startFace,
                        cellFaces,
                        anchors
                     );

                    if
                    (
                        checkLoop
                        (
                            mesh,
                            cutLoops[cellI],
                            edgeLoop,
                            boundaryPts,
                            cellI
                        )
                    )
                    {
                        cutLoops[cellI].append(edgeLoop);
                        nCuts++;
                    }
                    else
                    {
                        anchors = anchorsOrig;
                    }
                }
            }
        }

        forAll(cPoints, cPtI)
        {
            if (anchors[cPtI])
            {
                continue;
            }

            label startPt = cPoints[cPtI];
            label startEdge = -1;
            label startFace = -1;

            if (boundaryPts[startPt])
            {
                const labelList& pFaces = mesh.pointFaces()[startPt];

                forAll(pFaces, pFI)
                {
                    label faceI = pFaces[pFI];

                    if (cellFaces.found(faceI) && boundaryFaces[faceI])
                    {
                        const labelList& pEdges = mesh.pointEdges()[startPt];
                        forAll(pEdges, pEI)
                        {
                            label edgeI = pEdges[pEI];
                            if (cellEdges.found(edgeI))
                            {
                                if (cutEdges[edgeI] != -1)
                                {
                                    startEdge = edgeI;
                                    break;
                                }
                            }
                        }
                    }

                    if (startEdge != -1)
                    {
                        break;
                    }
                }

                if (startEdge != -1)
                {
                    const labelList& eFaces = mesh.edgeFaces()[startEdge];
                    forAll(eFaces, eFI)
                    {
                        label eFaceI = eFaces[eFI];

                        if (cellFaces.found(eFaceI))
                        {
                            startFace = eFaceI;
                            break;
                        }
                    }

                    boolList anchorsOrig = anchors;

                    labelList edgeLoop = createCutLoop
                    (
                        mesh,
                        cutPts,
                        cellPts,
                        startPt,
                        startEdge,
                        startFace,
                        cellFaces,
                        anchors
                     );

                    if
                    (
                        checkLoop
                        (
                            mesh,
                            cutLoops[cellI],
                            edgeLoop,
                            boundaryPts,
                            cellI
                        )
                    )
                    {
                        cutLoops[cellI].append(edgeLoop);
                        nCuts++;
                    }
                    else
                    {
                        anchors = anchorsOrig;
                    }
                }
            }
        }

        bool foundAnchors = true;
        forAll(cPoints, cPtI)
        {
            if (!anchors[cPtI] && boundaryPts[cPoints[cPtI]])
            {
                foundAnchors = false;
                break;
            }
        }

        if (foundAnchors)
        {
            continue;
        }

        bool foundCut = false;

        // walk from non boundary points
        forAll(cPoints, cPtI)
        {
            if (anchors[cPtI])
            {
                continue;
            }

            label startPt = cPoints[cPtI];
            label startEdge = -1;
            label startFace = -1;

            if (!boundaryPts[startPt] && !cutPts[cellPts[startPt]])
            {
                const labelList& pEdges = mesh.pointEdges()[startPt];
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (cellEdges.found(edgeI))
                    {
                        edge e = mesh.edges()[edgeI];
                        label otherPt(e[0] == startPt ? e[1] : e[0]);

                        if (cutPts[cellPts[otherPt]])
                        {
                            startEdge = edgeI;
                            break;
                        }
                    }
                }

                if (startEdge != -1)
                {
                    const labelList& eFaces = mesh.edgeFaces()[startEdge];
                    forAll(eFaces, eFI)
                    {
                        label eFaceI = eFaces[eFI];

                        if (cellFaces.found(eFaceI))
                        {
                            startFace = eFaceI;
                            break;
                        }
                    }
                    boolList anchorsOrig = anchors;

                    labelList edgeLoop = createCutLoop
                    (
                        mesh,
                        cutPts,
                        cellPts,
                        startPt,
                        startEdge,
                        startFace,
                        cellFaces,
                        anchors
                     );

                    if
                    (
                        checkLoop
                        (
                            mesh,
                            cutLoops[cellI],
                            edgeLoop,
                            boundaryPts,
                            cellI
                        )
                    )
                    {
                        foundCut = true;
                        cutLoops[cellI].append(edgeLoop);
                        nCuts++;
                    }
                    else
                    {
                        anchors = anchorsOrig;
                    }
                }
            }
        }

        if (!foundCut)
        {
            // walk from non boundary points
            forAll(cPoints, cPtI)
            {
                if (anchors[cPtI])
                {
                    continue;
                }

                label startPt = cPoints[cPtI];
                label startEdge = -1;
                label startFace = -1;

                if (boundaryPts[startPt])
                {
                    const labelList& pEdges = mesh.pointEdges()[startPt];
                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];
                        if (cellEdges.found(edgeI))
                        {
                            edge e = mesh.edges()[edgeI];
                            label otherPt(e[0] == startPt ? e[1] : e[0]);

                            if (cutPts[cellPts[otherPt]])
                            {
                                startEdge = edgeI;
                                break;
                            }
                        }
                    }

                    if (startEdge != -1)
                    {
                        const labelList& eFaces = mesh.edgeFaces()[startEdge];
                        forAll(eFaces, eFI)
                        {
                            label eFaceI = eFaces[eFI];

                            if (cellFaces.found(eFaceI))
                            {
                                startFace = eFaceI;
                                break;
                            }
                        }
                        boolList anchorsOrig = anchors;

                        labelList edgeLoop = createCutLoop
                        (
                            mesh,
                            cutPts,
                            cellPts,
                            startPt,
                            startEdge,
                            startFace,
                            cellFaces,
                            anchors
                         );

                        if
                        (
                            checkLoop
                            (
                                mesh,
                                cutLoops[cellI],
                                edgeLoop, boundaryPts,
                                cellI
                             )
                        )
                        {
                            foundCut = true;
                            cutLoops[cellI].append(edgeLoop);
                            nCuts++;
                        }
                        else
                        {
                            anchors = anchorsOrig;
                        }
                    }
                }
            }
        }

        if (nCuts > maxCuts)
        {
            maxCuts = nCuts;
        }
    }

    labelList nCuts(mesh.nCells(), 0);
    {
        DynamicList<label> cCells(mesh.nCells()/100);
        DynamicList<labelList> cCellLoops(mesh.nCells()/100);

        forAll(mesh.cells(), cellI)
        {
            cutLoops[cellI].shrink();
            if (cutLoops[cellI].size() == 1)
            {
                nCuts[cellI]++;
                cCells.append(cellI);
                cCellLoops.append(cutLoops[cellI][0]);
            }
        }
        allCutCells.append(cCells);
        allCellLoops.append(cCellLoops);
    }


    {
        boolList noSplit(mesh.nCells(), false);
//        label nLoops = 0;
        while (true)
        {
            DynamicList<label> cCells(mesh.nCells()/100);
            DynamicList<labelList> cCellLoops(mesh.nCells()/100);

//            nLoops++;
            label nSplit = 0;
            forAll(mesh.cells(), cellI)
            {
                if
                (
                    cutLoops[cellI].size() > 1
                    &&  nCuts[cellI] < cutLoops[cellI].size()
                    && !noSplit[cellI]
                 )
                {
                    nSplit++;
                    cCells.append(cellI);
                    cCellLoops.append(cutLoops[cellI][nCuts[cellI]]);
                    nCuts[cellI]++;

                    const labelList& cFaces = mesh.cells()[cellI];
                    forAll(cFaces, cFI)
                    {
                        label faceI = cFaces[cFI];
                        if (mesh.isInternalFace(faceI))
                        {
                            label own = mesh.faceOwner()[faceI];
                            label nbrCell
                            (
                                own == cellI ?
                                mesh.faceNeighbour()[faceI] :
                                own
                            );
                            noSplit[nbrCell] = true;
                        }
                    }
                }
            }
            reduce(nSplit, sumOp<label>());

            if (nSplit == 0)
            {
                break;
            }

            allCutCells.append(cCells);
            allCellLoops.append(cCellLoops);

            noSplit = false;
        }
    }

    allCutCells.shrink();

    forAll(allCutCells, loopI)
    {
        allCutCells[loopI].shrink();
        allCellLoops[loopI].shrink();
    }
}


void Foam::addHexMeshLayer::splitProcessorFaces
(
    const DynamicList<label>& cutCells,
    const DynamicList<labelList>& cellLoops,
    const List<scalarField>& cellWeights
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    polyTopoChange meshMod(mesh);

    DynamicList<label> procCutCells(cutCells.size()/10);
    DynamicList<labelList> procCellLoops(cutCells.size()/10);
    DynamicList<scalarField> procCellWeights(cutCells.size()/10);
    forAll(cutCells, i)
    {
        label cellI = cutCells[i];
        const labelList& cFaces = mesh.cells()[cellI];
        forAll(cFaces, cFI)
        {
            label faceI = cFaces[cFI];
            label patchI = patches.whichPatch(faceI);
            if (patchI != -1 && patches[patchI].coupled())
            {
                procCutCells.append(cellI);
                procCellLoops.append(cellLoops[i]);
                procCellWeights.append(cellWeights[i]);
                break;
            }
        }
    }

    cellCuts cuts(mesh, procCutCells, procCellLoops, procCellWeights);

    const labelListList& meshCellLoops = cuts.cellLoops();

    List<labelPair> masterCuts
    (
        mesh.nFaces()-mesh.nInternalFaces(),
        labelPair(-1,-1)
    );
    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);
            label startI = procPatch.start();

            forAll(procPatch, i)
            {
                label faceI = startI + i;
                label bFaceI = faceI-mesh.nInternalFaces();

                label cellI = mesh.faceOwner()[faceI];
                const labelList& cLoop = meshCellLoops[cellI];

                if (cLoop.size() == 0)
                {
                    continue;
                }

                labelHashSet cLoopSet(cLoop);
                face f = mesh.faces()[faceI];

                label nCutPts = 0;
                labelPair cut(-1, -1);
                forAll(f, fp)
                {
                    label pointI = f[fp];
                    if (cLoopSet.found(pointI))
                    {
                        cut[nCutPts] = pointI;
                        nCutPts++;
                    }
                    if (nCutPts == 2)
                    {
                        break;
                    }
                }
                if (nCutPts == 2)
                {
                    label meshEdgeI = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[cut[0]],
                        cut[0],
                        cut[1]
                    );
                    if (meshEdgeI == -1)
                    {
                        masterCuts[bFaceI] = cut;
                    }
                }
            }
        }
    }

    List<label> neiCuts0(mesh.nFaces()-mesh.nInternalFaces(),-1);
    List<label> neiCuts1(mesh.nFaces()-mesh.nInternalFaces(),-1);

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);
            const labelList& nbrPts = procPatch.nbrPoints();
            const labelList& meshPts = procPatch.meshPoints();

            label startI = procPatch.start();

            forAll(procPatch, i)
            {
                label faceI = startI + i;
                label bFaceI = faceI-mesh.nInternalFaces();
                if (masterCuts[bFaceI][0] != -1)
                {
                    face f = procPatch.localFaces()[i];
                    forAll(f, fp)
                    {
                        label pointI = f[fp];
                        label meshPointI = meshPts[pointI];

                        if (meshPointI == masterCuts[bFaceI][0])
                        {
                            neiCuts0[bFaceI] =
                                nbrPts[pointI];
                        }
                        else if (meshPointI == masterCuts[bFaceI][1])
                        {
                            neiCuts1[bFaceI] =
                                nbrPts[pointI];
                        }
                    }
                }
            }
        }
    }

    syncTools::swapBoundaryFaceList
    (
        mesh,
        neiCuts0
    );

    syncTools::swapBoundaryFaceList
    (
        mesh,
        neiCuts1
    );

    forAll(patches, patchI)
    {
        if (patches[patchI].coupled())
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);
            const labelList& meshPts = procPatch.meshPoints();

            label startI = procPatch.start();

            forAll(procPatch, i)
            {
                label faceI = startI + i;
                label bFaceI = faceI-mesh.nInternalFaces();

                DynamicList<labelPair> cuts(2);

                if (masterCuts[bFaceI][0] != -1)
                {
                    if (neiCuts0[bFaceI] != -1)
                    {
                        label meshPt0 =  meshPts[neiCuts0[bFaceI]];
                        label meshPt1 =  meshPts[neiCuts1[bFaceI]];
                        label mCut0 = masterCuts[bFaceI][0];
                        label mCut1 = masterCuts[bFaceI][1];

                        if
                        (
                            (mCut0 != meshPt0 && mCut0 != meshPt1)
                            || (mCut1 != meshPt0 && mCut1 != meshPt1)
                        )
                        {
                            cuts.append(masterCuts[bFaceI]);
                            cuts.append(labelPair(meshPt0,meshPt1));
                        }
                    }
                    else
                    {
                        cuts.append(masterCuts[bFaceI]);
                    }
                }
                else if (neiCuts0[bFaceI] != -1)
                {
                    label meshPt0 =  meshPts[neiCuts0[bFaceI]];
                    label meshPt1 =  meshPts[neiCuts1[bFaceI]];
                    cuts.append(labelPair(meshPt0,meshPt1));
                }

                if (cuts.size() == 0)
                {
                    continue;
                }

                label own = mesh.faceOwner()[faceI];

                label zoneID = mesh.faceZones().whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }

                face f = mesh.faces()[faceI];

                if (cuts.size() == 2)
                {
                    DynamicList<label> modFace(f.size());
                    DynamicList<label> newFace(f.size());

                    label meshPt0 = cuts[1][0];
                    label meshPt1 = cuts[1][1];

                    modFace.append(meshPt0);
                    label nextPt = findIndex(f, meshPt0);
                    forAll(f, fp)
                    {
                        nextPt = f.fcIndex(nextPt);
                        if (f[nextPt] == meshPt1)
                        {
                            modFace.append(meshPt1);
                            break;
                        }
                        else
                        {
                            modFace.append(f[nextPt]);
                        }
                    }

                    newFace.append(meshPt1);
                    nextPt = findIndex(f, meshPt1);
                    forAll(f, fp)
                    {
                        nextPt = f.fcIndex(nextPt);
                        if (f[nextPt] == meshPt0)
                        {
                            newFace.append(meshPt0);
                            break;
                        }
                        else
                        {
                            newFace.append(f[nextPt]);
                        }
                    }

                    label foundMod = 0;
                    forAll(modFace, mFI)
                    {
                        label pointI = modFace[mFI];
                        if (pointI == cuts[0][0] || pointI == cuts[0][1])
                        {
                            foundMod++;
                        }
                    }

                    if (foundMod == 2)
                    {
                        f = face(modFace);
                        meshMod.setAction
                        (
                            polyAddFace
                            (
                                face(newFace), // vertices
                                own,            // owner,
                                -1,             // neighbour,
                                -1,             // masterPointID,
                                -1,             // masterEdgeID,
                                faceI,    // masterFaceID,
                                false,          // flipFaceFlux,
                                patchI,         // patchID,
                                zoneID,         // zoneID,
                                zoneFlip        // zoneFlip
                             )
                         );
                    }
                    else
                    {
                        f = face(newFace);

                        meshMod.setAction
                        (
                            polyAddFace
                            (
                                face(modFace), // vertices
                                own,            // owner,
                                -1,             // neighbour,
                                -1,             // masterPointID,
                                -1,             // masterEdgeID,
                                faceI,    // masterFaceID,
                                false,          // flipFaceFlux,
                                patchI,         // patchID,
                                zoneID,         // zoneID,
                                zoneFlip        // zoneFlip
                            )
                        );
                    }
                }

                //Split face
                DynamicList<label> modFace(f.size());
                DynamicList<label> newFace(f.size());

                label meshPt0 = cuts[0][0];
                label meshPt1 = cuts[0][1];

                modFace.append(meshPt0);
                label nextPt = findIndex(f, meshPt0);
                forAll(f, fp)
                {
                    nextPt = f.fcIndex(nextPt);
                    if (f[nextPt] == meshPt1)
                    {
                        modFace.append(meshPt1);
                        break;
                    }
                    else
                    {
                        modFace.append(f[nextPt]);
                    }
                }
                newFace.append(meshPt1);
                nextPt = findIndex(f, meshPt1);
                forAll(f, fp)
                {
                    nextPt = f.fcIndex(nextPt);
                    if (f[nextPt] == meshPt0)
                    {
                        newFace.append(meshPt0);
                        break;
                    }
                    else
                    {
                        newFace.append(f[nextPt]);
                    }
                }

                // Modify the master face.
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        face(modFace),  // original face
                        faceI,          // label of face
                        own,            // owner
                        -1,             // neighbour
                        false,          // face flip
                        patchI,         // patch for face
                        false,          // remove from zone
                        zoneID,         // zone for face
                        zoneFlip        // face flip in zone
                     )
                );

                meshMod.setAction
                (
                    polyAddFace
                    (
                        face(newFace), // vertices
                        own,           // owner,
                        -1,            // neighbour,
                        -1,            // masterPointID,
                        -1,            // masterEdgeID,
                        faceI,         // masterFaceID,
                        false,         // flipFaceFlux,
                        patchI,        // patchID,
                        zoneID,        // zoneID,
                        zoneFlip       // zoneFlip
                     )
                 );
            }
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
    updateMesh(map);
}


void Foam::addHexMeshLayer::createSplits
(
    DynamicList<DynamicList<label>>& allCutCells,
    DynamicList<DynamicList<labelList>>& allCellLoops
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    forAll(allCutCells, loopI)
    {
        const DynamicList<label>& cutCells = allCutCells[loopI];
        const DynamicList<labelList>& cellLoops = allCellLoops[loopI];

        List<scalarField> cellWeights(cutCells.size());
        forAll(cellWeights, cwI)
        {
            cellWeights[cwI].setSize(cellLoops[cwI].size(), -GREAT);
        }

        //Pre split processor faces
        splitProcessorFaces(cutCells,cellLoops,cellWeights);

        //Cut cells
        {
            cellCuts cuts(mesh, cutCells, cellLoops, cellWeights);
            meshCutter mCut(mesh);

            polyTopoChange meshMod(mesh);
            mCut.setRefinement(cuts, meshMod);

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
            updateMesh(map);

            Map<label> newCells(map().nOldCells());
            forAll(map().cellMap(), cellI)
            {
                if (cellI >= map().nOldCells())
                {
                    label originalCell = map().cellMap()[cellI];
                    newCells.insert(originalCell,cellI);
                }
            }

            //update layer field
            if (trackLayerCells_)
            {
                volScalarField& layerCells = const_cast<volScalarField&>
                    (mesh.lookupObject<volScalarField>("layerStacks"));

                forAll(map().cellMap(), cellI)
                {
                    if (cellI >= map().nOldCells())
                    {
                        label originalCell = map().cellMap()[cellI];

                        scalar layerID = layerCells[originalCell];

                        if (layerID != -1)
                        {
                            layerCells[cellI] = layerID;
                        }
                    }
                }
            }

            for (int i = loopI+1; i < allCutCells.size(); i++)
            {
                labelList& updateCutCells = allCutCells[i];
                List<labelList>& updateCellLoops = allCellLoops[i];

                forAll(updateCellLoops, cLI)
                {
                    labelList& loop = updateCellLoops[cLI];
                    forAll(loop, j)
                    {
                        loop[j] = map().reversePointMap()[loop[j]];
                    }
                }

                forAll(updateCutCells, cLI)
                {
                    const labelList& loop = updateCellLoops[cLI];

                    label oldCellI = updateCutCells[cLI];
                    label newCellI = map().reverseCellMap()[oldCellI];

                    const labelList& newCellPts = mesh.cellPoints()[newCellI];
                    labelHashSet newPts(newCellPts);

                    bool allPts = true;
                    forAll(loop, j)
                    {
                        if (!newPts.found(loop[j]))
                        {
                            allPts = false;
                            break;
                        }
                    }

                    if (allPts)
                    {
                        updateCutCells[cLI] = newCellI;
                    }
                    else
                    {
                        Map<label>::const_iterator iter =
                            newCells.find(oldCellI);
                        if (iter != newCells.end())
                        {
                            updateCutCells[cLI] = iter();
                        }
                    }
                }
            }
        }
    }
}


void Foam::addHexMeshLayer::removeUnusedPoints()
{
    fvMesh& mesh = meshRefiner_.mesh();
    polyTopoChange meshMod(mesh);

    label nPointsRemoved = 0;
    boolList removedPts(mesh.nPoints(), false);

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    labelList nPointEdges(mesh.nPoints(), 0);

    forAll(mesh.points(), pointI)
    {
        const labelList pEdges = mesh.pointEdges()[pointI];
        forAll(pEdges, pEI)
        {
            label edgeI = pEdges[pEI];
            if (isMasterEdge[edgeI])
            {
                nPointEdges[pointI]++;
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        nPointEdges,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh.faces(), faceI)
    {
        face f = mesh.faces()[faceI];

        DynamicList<label> newFacePts(f.size());

        forAll(f, fp)
        {
            if (nPointEdges[f[fp]] == 2)
            {
                if (!removedPts[f[fp]])
                {
                    meshMod.setAction(polyRemovePoint(f[fp]));
                    removedPts[f[fp]] = true;
                    nPointsRemoved++;
                }
            }
            else
            {
                newFacePts.append(f[fp]);
            }
        }
        newFacePts.shrink();

        face newFace(newFacePts);
        if (newFace.size() != f.size())
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
                    face(newFace),              // modified face
                    faceI,                      // label of face
                    mesh.faceOwner()[faceI],    // owner
                    nei,                        // neighbour
                    false,                      // face flip
                    patchID,                    // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                 )
             );
        }
    }

    if (returnReduce(nPointsRemoved, sumOp<label>()))
    {
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
        updateMesh(map);
    }
}


void Foam::addHexMeshLayer::splitProblemCells()
{
    fvMesh& mesh = meshRefiner_.mesh();

    polyTopoChange meshMod(mesh);

    DynamicList<label>  cutCells(100);
    DynamicList<labelList>  cellLoops(100);

    labelList cutEdge = calculateCutEdges();

    boolList boundaryPts,boundaryEdges,boundaryFaces;
    calculateBoundary(boundaryPts, boundaryEdges, boundaryFaces);

    const labelListList& cEdgeList = mesh.cellEdges();

    labelList candidateCells;
    cellSet candidateCellSet(mesh, "candidateCells", mesh.nCells()/1000);
    labelList newEdgePts(mesh.nEdges(), -1);
    PackedList<1> changedFace(mesh.nFaces(), 0);

    forAll(mesh.cells(), cellI)
    {
        const labelList& cEdges = cEdgeList[cellI];
        label nCutEdges = 0;
        label startEdge = -1;

        forAll(cEdges, cEI)
        {
            label edgeI = cEdges[cEI];
            if (cutEdge[edgeI] != -1)
            {
                nCutEdges++;
                startEdge = edgeI;
            }
        }

        if (nCutEdges > 0 && nCutEdges < 3)
        {
            labelHashSet cellFaces(mesh.cells()[cellI]);
            edge cutEdge = mesh.edges()[startEdge];
            const labelList& eFaces = mesh.edgeFaces()[startEdge];

            DynamicList<label> eFIC(eFaces.size());
            forAll(eFaces, eFI)
            {
                label faceI = eFaces[eFI];
                if (cellFaces.found(faceI))
                {
                    eFIC.append(faceI);
                }
            }
            eFIC.shrink();
            labelHashSet edgeFacesInCell(eFIC);

            label boundPt(boundaryPts[cutEdge[0]] ? cutEdge[0] : cutEdge[1]);

            const labelList& pFaces = mesh.pointFaces()[boundPt];

            label grownFace = -1;
            forAll(pFaces, pFI)
            {
                label faceI = pFaces[pFI];
                if (cellFaces.found(faceI) && !edgeFacesInCell.found(faceI))
                {
                    grownFace = faceI;
                    break;
                }
            }

            if (grownFace != -1)
            {
                cutCells.append(cellI);

                labelHashSet cellEdges(cEdges);
                labelHashSet faceEdges(mesh.faceEdges()[grownFace]);

                face f = mesh.faces()[grownFace];
                label start = findIndex(f, boundPt);
                labelList cLoop(f.size());

                forAll(f, fp)
                {
                    label curPt = f[start];

                    const labelList& pEdges = mesh.pointEdges()[curPt];

                    label splitEdge = -1;
                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];
                        if (cellEdges.found(edgeI) && !faceEdges.found(edgeI))
                        {
                            splitEdge = edgeI;
                            break;
                        }
                    }

                    if (splitEdge != -1)
                    {
                        edge e = mesh.edges()[splitEdge];

                        if (newEdgePts[splitEdge] == -1)
                        {

                            vector eVec = e.vec(mesh.points());
                            vector cutPt = e[0];

                            if (boundaryPts[e[0]])
                            {
                                cutPt = mesh.points()[e[0]]
                                    + cutRatio_*eVec;
                            }
                            else
                            {
                                cutPt = mesh.points()[e[0]]
                                    + (1.-cutRatio_)*eVec;
                            }

                            newEdgePts[splitEdge] = meshMod.setAction
                            (
                                polyAddPoint
                                (
                                    cutPt,   // point
                                    e[0],       // master point
                                    -1,         // zone for point
                                    true        // supports a cell
                                 )
                            );
                            const labelList& eFaces =
                                mesh.edgeFaces()[splitEdge];
                            forAll(eFaces, eFI)
                            {
                                label faceI = eFaces[eFI];
                                changedFace.set(faceI, 1);
                            }
                        }
                        cLoop[fp] = newEdgePts[splitEdge];
                    }

                    start = f.fcIndex(start);
                }
                cellLoops.append(cLoop);
            }
        }

    }

    forAll(mesh.faces(), faceI)
    {
        if (changedFace.get(faceI) == 1)
        {
            const face f = mesh.faces()[faceI];
            DynamicList<label> newFace(2*f.size());

            forAll(f,fp)
            {
                label curPt = f[fp];
                newFace.append(curPt);

                label nextPt = f[f.fcIndex(fp)];
                label edgeI = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[curPt],
                    curPt,
                    nextPt
                 );

                if (newEdgePts[edgeI] != -1)
                {
                    newFace.append(newEdgePts[edgeI]);
                }
            }
            newFace.shrink();

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
                    face(newFace),              // modified face
                    faceI,                      // label of face
                    mesh.faceOwner()[faceI],    // owner
                    nei,                        // neighbour
                    false,                      // face flip
                    patchID,                    // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                 )
             );
        }
    }

    //Update mesh with new cut nodes
    {
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
        updateMesh(map);
    }

    //Split cells
    {
        List<scalarField> cellWeights(cutCells.size());
        forAll(cellWeights, cwI)
        {
            cellWeights[cwI].setSize(cellLoops[cwI].size(), -GREAT);
        }

        cellCuts cuts(mesh, cutCells, cellLoops, cellWeights);

        meshCutter mCut(mesh);

        polyTopoChange meshMod(mesh);
        mCut.setRefinement(cuts, meshMod);

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
        updateMesh(map);
    }
}


void Foam::addHexMeshLayer::checkForGapSideWalls()
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    if (mesh.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh.lookupObject<volScalarField>("gapCells");
        labelList gapEdges(mesh.nEdges(), 0);
        labelList gapFaces(mesh.nFaces(), 0);

        forAll(gapCells, celli)
        {
            if (gapCells[celli])
            {
                const cell& cFaces = mesh.cells()[celli];
                forAll(cFaces, cFI)
                {
                    gapFaces[cFaces[cFI]] = 1;
                }
                const labelList& cEdges = mesh.cellEdges()[celli];
                forAll(cEdges, cEI)
                {
                    gapEdges[cEdges[cEI]] = 1;
                }

            }
        }

        forAll(gapCells, celli)
        {
            if (gapCells[celli])
            {
                const cell& cFaces = mesh.cells()[celli];
                forAll(cFaces, cFI)
                {
                    label facei = cFaces[cFI];
                    label patchi = patches.whichPatch(facei);
                    if (patchi != -1 && !patches[patchi].coupled())
                    {
                        gapFaces[facei] = 2;
                        const labelList& fEdges = mesh.faceEdges()[facei];
                        forAll(fEdges,fEI)
                        {
                            gapEdges[fEdges[fEI]] = 2;
                        }
                    }
                }
            }
        }

        forAll(gapEdges, edgei)
        {
            if (gapEdges[edgei] == 1)
            {
                const edge e = mesh.edges()[edgei];
                forAll(e, eI)
                {
                    label pointi = e[eI];
                    const labelList pFaces = mesh.pointFaces()[pointi];
                    forAll(pFaces, pFI)
                    {
                        label facei = pFaces[pFI];
                        if (gapFaces[facei] > 0)
                        {
                            gapFaces[facei] = 3;
                        }
                    }
                }
            }
        }

        DynamicList<label> newAddressing(pp_.size());
        forAll(pp_, i)
        {
            label facei = pp_.addressing()[i];
            if (gapFaces[facei] != 2)
            {
                newAddressing.append(facei);
            }
        }

        pp_.clearOut();
        pp_.resetAddressing(newAddressing);
        pp_.localPoints();
    }
}


void Foam::addHexMeshLayer::updateMesh(const mapPolyMesh& map)
{
//    labelList newAddressing(pp_.size());
    DynamicList<label> newAddressing(pp_.size());

    forAll(pp_, i)
    {
        label oldFaceI = pp_.addressing()[i];

        label newFaceI = map.reverseFaceMap()[oldFaceI];
        if (newFaceI != -1)
        {
            newAddressing.append(newFaceI);
        }
    }
    newAddressing.shrink();

    pp_.clearOut();
    pp_.resetAddressing(newAddressing);
    pp_.localPoints();
}


void Foam::addHexMeshLayer::setRefinement()
{
    Info<< nl << "Adding layers to boundary touching cells"<<endl;

    //Check for gap side walls to allow them to be split
    checkForGapSideWalls();

    //Perform redistribution of mesh if required
    balance();

    // Add split edge points
    addSplitEdgePoints();

    DynamicList<DynamicList<label>> allCutCells(100);
    DynamicList<DynamicList<labelList>> allCellLoops(100);

    calculateSplitLoops(allCutCells,allCellLoops);

    createSplits(allCutCells,allCellLoops);

    removeUnusedPoints();

    return;
}


// ************************************************************************* //
