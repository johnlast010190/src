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
    (c) 2016-2017 Esi Ltd.

InClass
    autoHexDualiser

\*---------------------------------------------------------------------------*/

#include "autoHexDualiser/autoHexDualiser.H"
#include "meshTools/meshTools.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveCell.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "meshes/polyMesh/syncTools/dummyTransform.H"
#include "regionSplit/regionSplit.H"
#include "refinementSurfaces/refinementSurfaces.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(autoHexDualiser, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::autoHexDualiser::calcFaceLoop
(
    const extendedEdgeCells& eEC,
    face& f,
    pointField& pts
)
{
    const pointField& ePts = eEC.edgeCentres();
    const labelList& eCells = eEC.edgeCells();
    const List<labelPair>& eFOwnNei = eEC.edgeFacesOwnNei();
    f.setSize(eCells.size());
    pts.setSize(eCells.size());

    label start = eFOwnNei[0][0];
    label next = eFOwnNei[0][1];
    label index = 0;
    f[index++] = start;
    f[index++] = next;

    bool finished = false;
    boolList marked(eFOwnNei.size(), false);
    marked[0] = true;

    while (!finished)
    {
        forAll(eFOwnNei, i)
        {
            if (!marked[i])
            {
                const labelPair ownNei = eFOwnNei[i];
                if (ownNei[0] == next)
                {
                    marked[i] = true;
                    if (ownNei[1] == start)
                    {
                        finished = true;
                    }
                    else
                    {
                        next = ownNei[1];
                        f[index++] = next;
                    }
                    break;
                }
                else if (ownNei[1] == next)
                {
                    marked[i] = true;
                    if (ownNei[0] == start)
                    {
                        finished = true;
                    }
                    else
                    {
                        next = ownNei[0];
                        f[index++] = next;
                    }
                    break;
                }
            }
        }
    }

    Map<label> pointMap(f.size());
    forAll(eCells, eCI)
    {
        pointMap.insert(eCells[eCI],eCI);
    }
    forAll(f, fp)
    {
        pts[fp] = ePts[pointMap[f[fp]]];
    }
}


void Foam::autoHexDualiser::calcBoundaryZones
(
    labelList& pointCellZones
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const cellZoneMesh& cellZones = mesh.cellZones();

    pointCellZones.setSize(mesh.nPoints(), -1);

    //Check for named surfaces
    const PtrList<surfaceZonesInfo>& surfZones =
        meshRefiner_.surfaces().surfZones();

    DynamicList<label> boundaryNamedZones(cellZones.size());
    forAll(surfZones, surfI)
    {
        const word& faceZoneName = surfZones[surfI].faceZoneName();

        if (faceZoneName.size())
        {
            const surfaceZonesInfo::faceZoneType& faceType =
                surfZones[surfI].faceType();

            if (faceType == surfaceZonesInfo::BOUNDARY)
            {
                const word& cellZoneName = surfZones[surfI].cellZoneName();
                if (cellZoneName.size())
                {
                    label zoneI = cellZones.findZoneID(cellZoneName);
                    if (zoneI != -1)
                    {
                        boundaryNamedZones.append(zoneI);
                    }
                }
            }
        }
    }
    labelHashSet boundaryNamedZonesSet(boundaryNamedZones);

    forAll(cellZones, cellZoneI)
    {
        if (boundaryNamedZonesSet.found(cellZoneI))
        {
            const cellZone& cZone = cellZones[cellZoneI];
            forAll(cZone, czi)
            {
                label cellI = cZone[czi];
                const labelList& cPts = mesh.cellPoints()[cellI];
                forAll(cPts, cptI)
                {
                    label pointI = cPts[cptI];
                    pointCellZones[pointI] = cellZoneI;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pointCellZones,
        maxEqOp<label>(),  // combine op
        label(-1) // null value
    );
}


void Foam::autoHexDualiser::calcDual()
{
    Info<<"Calculating dualised mesh"<<nl<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    polyTopoChange meshMod(mesh);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    //update gap cells field if present
    boolList allGapPts(mesh.nPoints(), false);
    bool updateGapField = false;

    if (mesh.foundObject<volScalarField>("gapCells"))
    {
        updateGapField = true;
        const volScalarField& gapCells = const_cast<volScalarField&>
            (mesh.lookupObject<volScalarField>("gapCells"));

        forAll(gapCells, celli)
        {
            if (gapCells[celli] > -1)
            {
                const labelList& cPts = mesh.cellPoints()[celli];
                forAll(cPts, cPI)
                {
                    allGapPts[cPts[cPI]] = true;
                }
            }
        }

        boolList internalPts(mesh.nPoints(),true);
        forAll(patches, patchI)
        {
            if (!patches[patchI].coupled())
            {
                const polyPatch& pp = patches[patchI];
                forAll(pp.meshPoints(), ptI)
                {
                    label meshPointI = pp.meshPoints()[ptI];
                    const labelList& pCells = mesh.pointCells()[meshPointI];
                    forAll(pCells, pCI)
                    {
                        const labelList& cPts = mesh.cellPoints()[pCells[pCI]];
                        forAll(cPts, cPtI)
                        {
                            internalPts[cPts[cPtI]] = false;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            internalPts,
            andEqOp<bool>(),  // combine op
            false // null value
        );

        forAll(internalPts, pointi)
        {
            if (allGapPts[pointi])
            {
                const labelList& pEdges = mesh.pointEdges()[pointi];
                forAll(pEdges, pEI)
                {
                    const edge e = mesh.edges()[pEdges[pEI]];
                    label otherPt(e[0] == pointi ? e[1] : e[0]);
                    if (internalPts[otherPt])
                    {
                        allGapPts[pointi] = false;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            allGapPts,
            orEqOp<bool>(),  // combine op
            false // null value
        );
    }

    //Calculate the mesh regions and mark point regions
    const labelIOList& cellRegionID = meshRefiner_.cellRegionID();

    labelList pointRegion(mesh.nPoints(), -1);
    forAll(mesh.cells(), cellI)
    {
        label cRegion = cellRegionID[cellI];
        const labelList& cPts = mesh.cellPoints()[cellI];
        forAll(cPts, cPtI)
        {
            pointRegion[cPts[cPtI]] = cRegion;
        }
    }

    labelList pointCellZones;
    calcBoundaryZones(pointCellZones);

    const PackedBoolList isMasterPoint(syncTools::getMasterPoints(mesh));
    const PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    DynamicList<label> nonCoupledPatches(patches.size());
    DynamicList<label> coupledPatches(patches.size());
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            nonCoupledPatches.append(patchI);
        }
        else
        {
            coupledPatches.append(patchI);
        }
    }
    nonCoupledPatches.shrink();
    coupledPatches.shrink();

    autoPtr<indirectPrimitivePatch> ppBoundPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            nonCoupledPatches
        )
    );
    indirectPrimitivePatch& ppBound = ppBoundPtr();

    labelList meshEdges(ppBound.meshEdges(mesh.edges(), mesh.pointEdges()));

    // Edge status:
    //  0 : is boundary edge.
    //  -1 : is internal edge.
    labelList edgeToDualPoint(mesh.nEdges(), -1);
    forAll(meshEdges, patchEdgeI)
    {
        label edgeI = meshEdges[patchEdgeI];

        edgeToDualPoint[edgeI] = 0;
    }

    syncTools::syncEdgeList
    (
        mesh,
        edgeToDualPoint,
        maxEqOp<label>(),
        label(-1)
    );

    // Point status:
    //  >0 : dual point patchID
    //  -1 : is coupled point
    //  -2 : is internal point
    labelList pointToDualPoint(mesh.nPoints(), -2);

    forAll(coupledPatches, i)
    {
        label patchI = coupledPatches[i];
        const polyPatch& pp = patches[patchI];
        forAll(pp.meshPoints(), ptI)
        {
            label meshPointI = pp.meshPoints()[ptI];
            pointToDualPoint[meshPointI] = -1;
        }
    }

    forAll(nonCoupledPatches, i)
    {
        label patchI = nonCoupledPatches[i];
        const polyPatch& pp = patches[patchI];
        forAll(pp.meshPoints(), ptI)
        {
            pointToDualPoint[pp.meshPoints()[ptI]] = patchI;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pointToDualPoint,
        maxEqOp<label>(),  // combine op
        labelMin // null value
    );


    //Remove oiginal mesh points, faces, cells
    forAll(mesh.points(), pointI)
    {
        meshMod.setAction(polyRemovePoint(pointI));
    }

    forAll(mesh.faces(), faceI)
    {
        meshMod.setAction(polyRemoveFace(faceI));
    }

    forAll(mesh.cells(), cellI)
    {
        meshMod.setAction(polyRemoveCell(cellI));
    }

    // Assign cells
    // ~~~~~~~~~~~~

    labelList dualCells(mesh.nPoints(),-1);
    DynamicList<label> dualCellLevel(mesh.nPoints());
    DynamicList<label> dualCellRegion(mesh.nPoints());
    DynamicList<scalar> dualGapCells(mesh.nPoints());

    labelList pointCellMaxLevel(mesh.nPoints(), labelMin);

    forAll(mesh.points(), pointI)
    {
        const labelList& pCells = mesh.pointCells()[pointI];
        forAll(pCells, pCI)
        {
            label cellI = pCells[pCI];
            pointCellMaxLevel[pointI] = max
            (
                pointCellMaxLevel[pointI],
                cellLevel[cellI]
            );
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pointCellMaxLevel,
        maxEqOp<label>(),  // combine op
        labelMin // null value
    );

    const pointField& pts = mesh.points();
    labelList newCellCount(mesh.nPoints(),-1);
    labelList reverseNewCellCount(mesh.nPoints(),-1);

    label newCellCountI = 0;
    forAll(pts, pointI)
    {
        if (pointToDualPoint[pointI] < 0 && isMasterPoint.get(pointI) == 1)
        {
            dualCellLevel.append(pointCellMaxLevel[pointI]);

            if (updateGapField)
            {
                if (allGapPts[pointI])
                {
                    dualGapCells.append(scalar(1));
                }
                else
                {
                    dualGapCells.append(scalar(-1));
                }
            }

            dualCellRegion.append(pointRegion[pointI]);
            label cellZoneI = pointCellZones[pointI];
            dualCells[pointI] = meshMod.setAction
            (
                polyAddCell
                (
                    -1,             // master point
                    -1,             // master edge
                    -1,             // master face
                    -1,             // master cell
                    cellZoneI       // zone for cell
                )
            );

            newCellCount[pointI] = newCellCountI;
            reverseNewCellCount[newCellCountI] = pointI;
            newCellCountI++;
        }
    }
    dualCellLevel.shrink();
    dualCellRegion.shrink();

    globalIndex globalNewCells(newCellCountI);

    labelList globalNewCellIndex(mesh.nPoints(), -1);

    forAll(pts, pointI)
    {
        label newCellI = newCellCount[pointI];
        if (newCellI != -1)
        {
            globalNewCellIndex[pointI] = globalNewCells.toGlobal(newCellI);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        globalNewCellIndex,
        maxEqOp<label>(),  // combine op
        labelMin // null value
    );

    // Assign faces
    // ~~~~~~~~~~~~
    globalIndex globalCells(mesh.nCells());
    globalIndex globalEdges(mesh.nEdges());

    labelList globalNei(mesh.nFaces()-mesh.nInternalFaces());
    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        label own = mesh.faceOwner()[faceI];
        globalNei[faceI-mesh.nInternalFaces()] = globalCells.toGlobal(own);
    }
    syncTools::swapBoundaryFaceList(mesh, globalNei);

    List<face> newFaces(mesh.nEdges(), face(0));
    List<pointField> newFacePts(mesh.nEdges(),pointField(0));

    // Calculate face loops
    {
        EdgeMap<extendedEdgeCells> extendedEdges;

        forAll(mesh.edges(), edgeI)
        {
            const edge e =  mesh.edges()[edgeI];
            if (pointToDualPoint[e[0]] != -1 && pointToDualPoint[e[1]] != -1)
            {
                continue;
            }

            if (edgeToDualPoint[edgeI] == -1)
            {
                const labelList& eCells = mesh.edgeCells()[edgeI];
                labelList eCellsGlobal(eCells.size());
                pointField eCC(eCells.size());

                forAll(eCells, eCI)
                {
                    label cellI = eCells[eCI];
                    eCellsGlobal[eCI] = globalCells.toGlobal(cellI);
                    eCC[eCI] = mesh.cellCentres()[cellI];
                }

                const labelList& eFaces = mesh.edgeFaces()[edgeI];
                DynamicList<labelPair> eFacesOwnNeiGlobal(eFaces.size());

                forAll(eFaces, eFI)
                {
                    label faceI = eFaces[eFI];
                    if (isMasterFace.get(faceI) == 1)
                    {
                        label gOwn =
                            globalCells.toGlobal(mesh.faceOwner()[faceI]);
                        label gNei = -1;
                        if (mesh.isInternalFace(faceI))
                        {
                            gNei = globalCells.toGlobal
                            (mesh.faceNeighbour()[faceI]);
                        }
                        else
                        {
                            gNei = globalNei[faceI-mesh.nInternalFaces()];
                        }
                        eFacesOwnNeiGlobal.append(labelPair(gOwn,gNei));
                    }
                }
                eFacesOwnNeiGlobal.shrink();

                label globalEdgeID = globalEdges.toGlobal(edgeI);
                extendedEdgeCells eEC
                (
                    eCellsGlobal,
                    eCC,
                    eFacesOwnNeiGlobal,
                    globalEdgeID
                );

                extendedEdges.insert
                (
                    mesh.edges()[edgeI],
                    eEC
                );
            }
        }

        syncTools::syncEdgeMap
        (
            mesh,
            extendedEdges,
            extendedEdgeCells::sumEqOp(),
            dummyTransform()
        );

        forAll(mesh.edges(), edgeI)
        {
            if (edgeToDualPoint[edgeI] == -1 && isMasterEdge.get(edgeI) == 1)
            {
                const edge e =  mesh.edges()[edgeI];
                if (pointToDualPoint[e[0]] != -1 && pointToDualPoint[e[1]] != -1)
                {
                    continue;
                }
                const extendedEdgeCells& eEC = extendedEdges[e];

                pointField facePts(0);
                face f(0);
                calcFaceLoop(eEC,f,facePts);
                newFaces[edgeI] = f;
                newFacePts[edgeI] = facePts;
            }
        }
    }

    DynamicList<label> dualPointLevel(mesh.nCells());

    // Calculate new points and faces
    {
        Map<extendedPointEdgeFaces> extendedPoints;
        forAll(pts, pointI)
        {
            if (pointToDualPoint[pointI] == -1)
            {
                const labelList& pEdges = mesh.pointEdges()[pointI];

                DynamicList<face> pEdgeFaces(pEdges.size());
                DynamicList<pointField> pEdgeFacePts(pEdges.size());
                DynamicList<word> edgePatchNames(pEdges.size());
                DynamicList<label> neighbours(pEdges.size());

                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if
                    (
                        isMasterEdge.get(edgeI) == 1
                        && newFaces[edgeI].size() > 0
                    )
                    {
                        const edge e = mesh.edges()[edgeI];
                        face f = newFaces[edgeI];
                        pointField facePts = newFacePts[edgeI];
                        label otherPt = (e[0] ==  pointI ? e[1] : e[0]);

                        label nei = -1;

                        label otherGlobal = globalNewCellIndex[otherPt];
                        label masterGlobal = globalNewCellIndex[pointI];
                        vector eVec = mesh.points()[otherPt]
                            - mesh.points()[pointI];
                        vector fNorm =
                            face(identity(f.size())).areaNormal(facePts);

                        word patchName("internal");
                        bool flip = false;
                        if (otherGlobal == -1)
                        {
                            label patchI = pointToDualPoint[otherPt];
                            patchName = patches[patchI].name();

                            if ((eVec & fNorm) < 0)
                            {
                                flip = true;
                            }
                        }
                        else
                        {
                            label otherProcID =
                                globalNewCells.whichProcID(otherGlobal);
                            label masterProcID =
                                globalNewCells.whichProcID(masterGlobal);

                            if (otherProcID == masterProcID) //internal face
                            {
                                if (masterGlobal > otherGlobal)
                                {
                                    continue;
                                }
                                else
                                {
                                    nei = otherGlobal;
                                    if ((eVec & fNorm) < 0)
                                    {
                                        flip = true;
                                    }
                                }
                            }
                            else
                            {
                                patchName =
                                    "procBoundary"
                                    + Foam::name(masterProcID)
                                    + "to"
                                    + Foam::name(otherProcID);

                                if ((eVec & fNorm) < 0)
                                {
                                    flip = true;
                                }
                            }
                        }
                        if (flip)
                        {
                            label n = f.size();
                            for (label i=1; i < (n+1)/2; ++i)
                            {
                                Swap(f[i], f[n-i]);
                                Swap(facePts[i],facePts[n-i]);
                            }
                        }
                        pEdgeFaces.append(f);
                        pEdgeFacePts.append(facePts);
                        edgePatchNames.append(patchName);
                        neighbours.append(nei);
                    }
                }
                pEdgeFaces.shrink();
                edgePatchNames.shrink();
                neighbours.shrink();
                DynamicList<label> pPtsAddr(pEdgeFaces.size()*10);
                DynamicList<point> pPts(pEdgeFaces.size()*10);
                if (pEdgeFaces.size() > 0)
                {
                    labelHashSet gPts;
                    forAll(pEdgeFaces, pEFI)
                    {
                        face f = pEdgeFaces[pEFI];
                        const pointField& facePts = pEdgeFacePts[pEFI];
                        forAll(f,fp)
                        {
                            label gPointI = f[fp];
                            if (gPts.insert(gPointI))
                            {
                                pPtsAddr.append(gPointI);
                                pPts.append(facePts[fp]);
                            }
                        }
                    }
                    extendedPointEdgeFaces ePEF
                    (
                        pEdgeFaces,
                        pPts,
                        pPtsAddr,
                        edgePatchNames,
                        neighbours
                     );

                    extendedPoints.insert
                    (
                        pointI,
                        ePEF
                     );
                }
            }
        }

        syncTools::syncPointMap
        (
            mesh,
            extendedPoints,
            extendedPointEdgeFaces::sumEqOp(),
            dummyTransform()
        );

        // Create new processor patches and add points and faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Map<label> mapPointId(2*mesh.nCells());
        List<List<labelPair>>
            pointCellsLevel(mesh.nPoints(),List<labelPair>(0));

        forAll(mesh.points(), pointI)
        {
            if (pointToDualPoint[pointI] == -1)
            {
                const labelList& pCells = mesh.pointCells()[pointI];
                List<labelPair> pCellLevel(pCells.size());
                forAll(pCells, pCellI)
                {
                    label cellI = pCells[pCellI];
                    pCellLevel[pCellI]= labelPair
                    (
                        globalCells.toGlobal(cellI),
                        cellLevel[cellI]
                    );
                }
                pointCellsLevel[pointI] = pCellLevel;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointCellsLevel,
            unionEqOp(),  // combine op
            List<labelPair>(0), // null value
            dummyTransform()
        );

        //Add points
        forAll(pts, pointI)
        {
            if (dualCells[pointI] != -1)
            {
                if
                (
                    pointToDualPoint[pointI] == -1
                    && extendedPoints.found(pointI)
                )
                {
                    const extendedPointEdgeFaces& ePEF = extendedPoints[pointI];
                    const labelList& addressing = ePEF.addressing();
                    const List<point>& pPts = ePEF.points();
                    forAll(addressing, ptI)
                    {
                        label globalPtIndex = addressing[ptI];
                        if (!mapPointId.found(globalPtIndex))
                        {
                            if (globalCells.isLocal(globalPtIndex))
                            {
                                label cellI =
                                    globalCells.toLocal
                                    (
                                        Pstream::myProcNo(),
                                        globalPtIndex
                                    );
                                dualPointLevel.append(cellLevel[cellI]);
                            }
                            else
                            {
                                const List<labelPair>& pCL =
                                    pointCellsLevel[pointI];
                                forAll(pCL, pCLI)
                                {
                                    if (pCL[pCLI].first() == globalPtIndex)
                                    {
                                        dualPointLevel.append
                                        (pCL[pCLI].second());
                                        break;
                                    }
                                }
                            }
                            label newPointI = meshMod.setAction
                            (
                                polyAddPoint
                                (
                                    pPts[ptI], // point
                                    -1,         // master point
                                    -1,      // zone for point
                                    true        // supports a cell
                                )
                            );
                            mapPointId.insert(globalPtIndex, newPointI);
                        }
                    }
                }
            }
        }

        //Add all other untouched cells
        forAll(mesh.cells(), cellI)
        {
            label globalCellI = globalCells.toGlobal(cellI);
            if (!mapPointId.found(globalCellI))
            {
                const labelList& cPts = mesh.cellPoints()[cellI];
                label nAdded = 0;
                forAll(cPts, cPtI)
                {
                    if (dualCells[cPts[cPtI]] == -1)
                    {
                        nAdded++;
                    }
                    else
                    {
                        break;
                    }
                }

                if (nAdded == cPts.size())
                {
                    continue;
                }

                dualPointLevel.append(cellLevel[cellI]);
                label newPointI = meshMod.setAction
                (
                    polyAddPoint
                    (
                        mesh.cellCentres()[cellI], // point
                        -1,         // master point
                        -1,      // zone for point
                        true        // supports a cell
                    )
                 );
                mapPointId.insert(globalCellI, newPointI);
            }
        }

        //Add any new processor patches that may be generated by dualisation
        forAll(pts, pointI)
        {
            if (dualCells[pointI] != -1)
            {
                if
                (
                    pointToDualPoint[pointI] == -1
                    && extendedPoints.found(pointI)
                )
                {
                    const extendedPointEdgeFaces& ePEF = extendedPoints[pointI];
                    const faceList& edgeFaces = ePEF.edgeFaces();
                    const wordList& edgePatchNames = ePEF.edgePatchNames();

                    forAll(edgeFaces, eFI)
                    {
                        word patchName = edgePatchNames[eFI];
                        if (patchName != "internal")
                        {
                            label patchI = patches.findPatchID(patchName);
                            if (patchI  == -1)
                            {
                                label nbrProcI = readLabel
                                (
                                    IStringStream
                                    (
                                        patchName.substr
                                        (patchName.find("to") + 2)
                                    )()
                                );

                                dictionary patchDict;
                                patchDict.add
                                (
                                    "type", processorPolyPatch::typeName
                                );
                                patchDict.add("myProcNo", Pstream::myProcNo());
                                patchDict.add("neighbProcNo", nbrProcI);
                                patchDict.add("nFaces", 0);
                                patchDict.add("startFace", mesh.nFaces());

                                meshRefiner_.appendPatch
                                (
                                    mesh,
                                    patchName,
                                    patchDict
                                );
                            }
                        }
                    }
                }
            }
        }
        //Update number of patches in meshMod
        meshMod.setNumPatches(mesh.boundaryMesh().size());

        //Add faces
        forAll(pts, pointI)
        {
            if (dualCells[pointI] != -1)
            {
                if
                (
                    pointToDualPoint[pointI] == -1
                    && extendedPoints.found(pointI)
                )
                {
                    const extendedPointEdgeFaces& ePEF = extendedPoints[pointI];
                    const faceList& edgeFaces = ePEF.edgeFaces();
                    const wordList& edgePatchNames = ePEF.edgePatchNames();
                    const labelList& neighbours = ePEF.neighbours();

                    forAll(edgeFaces, eFI)
                    {
                        face newFace = edgeFaces[eFI];
                        forAll(newFace, fp)
                        {
                            newFace[fp] = mapPointId[newFace[fp]];
                        }

                        word patchName = edgePatchNames[eFI];
                        label patchI = -1;
                        label own = dualCells[pointI];
                        label nei = -1;

                        if (patchName != "internal")
                        {
                            patchI = patches.findPatchID(patchName);
                        }
                        else
                        {
                            label otherPt =
                                globalNewCells.toLocal(neighbours[eFI]);
                            nei = dualCells[reverseNewCellCount[otherPt]];
                        }

                        label zoneI = -1;
                        bool flipZone = false;

                        meshMod.setAction
                        (
                            polyAddFace
                            (
                                newFace,          // face
                                own,              // owner
                                nei,              // neighbour
                                -1,               // master point
                                -1,               // master edge
                                -1,               // master face
                                false,            // flux flip
                                patchI,           // patch for face
                                zoneI,            // zone for face
                                flipZone          // face zone flip
                             )
                         );
                    }
                }
            }
        }

        //Add all other faces
        forAll(pts, pointI)
        {
            if (dualCells[pointI] != -1)
            {
                if (pointToDualPoint[pointI] == -2)
                {
                    const labelList& pEdges = mesh.pointEdges()[pointI];

                    forAll(pEdges, pEI)
                    {
                        label edgeI = pEdges[pEI];

                        if
                        (
                            isMasterEdge.get(edgeI) == 1
                            && edgeToDualPoint[edgeI] == -1
                        )
                        {
                            const edge e = mesh.edges()[edgeI];
                            label otherPt = (e[0] ==  pointI ? e[1] : e[0]);
                            label otherGlobal = globalNewCellIndex[otherPt];
                            label masterGlobal = globalNewCellIndex[pointI];

                            if
                            (
                                otherGlobal != -1
                                && (masterGlobal > otherGlobal)
                                && dualCells[otherPt] != -1
                            )
                            {
                                continue;
                            }

                            // Find dual owner, neighbour
                            label owner = -1;
                            label neighbour = -1;
                            label patchI = -1;

                            if (pointToDualPoint[otherPt] > -1)
                            {
                                patchI = pointToDualPoint[otherPt];
                            }
                            else if
                            (
                                pointToDualPoint[otherPt] == -1
                                && dualCells[otherPt] == -1
                            )
                            {
                                label otherProcID =
                                    globalNewCells.whichProcID(otherGlobal);
                                label masterProcID =
                                    globalNewCells.whichProcID(masterGlobal);

                                word patchName =
                                    "procBoundary"
                                    + Foam::name(masterProcID)
                                    + "to"
                                    + Foam::name(otherProcID);

                                patchI = patches.findPatchID(patchName);

                                if (patchI == -1)
                                {
                                    FatalErrorInFunction
                                        << "Requires processor patch but could"
                                        << " not find processor patch : "
                                        << patchName
                                        << " check dual faces constructed "
                                        << " at point : " << pointI
                                        << abort(FatalError);
                                }
                            }

                            if (e[0] < e[1])
                            {
                                owner = e[0];
                                neighbour = e[1];
                            }
                            else
                            {
                                owner = e[1];
                                neighbour = e[0];
                            }

                            // Get a starting cell
                            const labelList& eCells = mesh.edgeCells()[edgeI];

                            label cellI = eCells[0];

                            // Get the two faces on the cell and edge.
                            label face0, face1;
                            meshTools::getEdgeFaces
                                (mesh, cellI, edgeI, face0, face1);

                            // Find the starting face by looking at the order
                            // in which the face uses the owner, neighbour
                            const face& f0 = mesh.faces()[face0];

                            label index = findIndex(f0, neighbour);

                            bool f0OrderOk = (f0.nextLabel(index) == owner);

                            label startFaceI = -1;
                            if (f0OrderOk == (mesh.faceOwner()[face0] == cellI))
                            {
                                startFaceI = face0;
                            }
                            else
                            {
                                startFaceI = face1;
                            }

                            // Walk face-cell-face until starting face reached.
                            DynamicList<label> dualFace
                            (
                                mesh.edgeCells()[edgeI].size()
                            );

                            label faceI = startFaceI;
                            while (true)
                            {
                                label globalCellI = globalCells.toGlobal(cellI);
                                // Store dual vertex from cellI.
                                dualFace.append(mapPointId[globalCellI]);

                                // Cross cell to other face on edge.
                                label f0, f1;
                                meshTools::getEdgeFaces
                                    (mesh, cellI, edgeI, f0, f1);
                                if (f0 == faceI)
                                {
                                    faceI = f1;
                                }
                                else
                                {
                                    faceI = f0;
                                }

                                // Cross face to other cell.
                                if (faceI == startFaceI)
                                {
                                    break;
                                }

                                if (mesh.faceOwner()[faceI] == cellI)
                                {
                                    cellI = mesh.faceNeighbour()[faceI];
                                }
                                else
                                {
                                    cellI = mesh.faceOwner()[faceI];
                                }
                            }

                            label zoneI = -1;
                            bool flipZone = false;
                            bool flipFace = false;
                            face newFace(dualFace.shrink());
                            label dualOwner = dualCells[owner];
                            label dualNeighbour = -1;

                            if (patchI == -1)
                            {
                                dualNeighbour = dualCells[neighbour];
                            }
                            else
                            {
                                label pown = (e[0] != pointI ? e[1]  : e[0]);
                                if (owner != pown)
                                {
                                    flipFace = true;
                                }
                                dualOwner = dualCells[pown];
                            }

                            if (flipFace)
                            {
                                newFace.flip();
                            }

                            meshMod.setAction
                            (
                                polyAddFace
                                (
                                    newFace,         // face
                                    dualOwner,       // owner
                                    dualNeighbour,   // neighbour
                                    -1,              // master point
                                    -1,              // master edge
                                    -1,              // master face
                                    false,           // flux flip
                                    patchI,          // patch for face
                                    zoneI,           // zone for face
                                    flipZone         // face zone flip
                                 )
                             );
                        }
                    }
                }
            }
        }
    }

    dualPointLevel.shrink();

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

    meshRefiner_.meshCutter().updateLevels(dualPointLevel,dualCellLevel);

    if (updateGapField)
    {
        volScalarField& gapCells = const_cast<volScalarField&>
            (mesh.lookupObject<volScalarField>("gapCells"));

        forAll(dualGapCells, celli)
        {
            gapCells[celli] = dualGapCells[celli];
        }
    }

    meshRefiner_.updateRegionID(dualCellRegion);
}


void Foam::autoHexDualiser::splitBoundaryFaces()
{
    fvMesh& mesh = meshRefiner_.mesh();
    polyTopoChange meshMod(mesh);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (!pp.coupled())
        {
            label startFaceI = pp.start();

            forAll(pp, i)
            {
                label faceI = startFaceI+i;
                label own = mesh.faceOwner()[faceI];
                const cell& ownFaces = mesh.cells()[own];
                if (ownFaces.size() != 4)
                {
                    continue;
                }

                labelHashSet cFaces(ownFaces);

                face f = mesh.faces()[faceI];
                label degPt = -1;

                forAll(f, fp)
                {
                    const labelList& pFaces = mesh.pointFaces()[f[fp]];
                    label nPointCellFaces = 0;
                    forAll(pFaces, pFI)
                    {
                        if (cFaces.found(pFaces[pFI]))
                        {
                            nPointCellFaces++;
                        }
                    }
                    if (nPointCellFaces == 2)
                    {
                        degPt = f[fp];
                        break;
                    }
                }

                if (degPt != -1)
                {
                    label start = findIndex(f, degPt);
                    face newFace(3);
                    newFace[0] = f[start];
                    label next = f.fcIndex(start);
                    label nAdded = 0;

                    label zoneID = mesh.faceZones().whichZone(faceI);

                    while (true)
                    {
                        newFace[1] = f[next];
                        next = f.fcIndex(next);
                        newFace[2] = f[next];

                        bool zoneFlip = false;
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = mesh.faceZones()[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                        }


                        if (nAdded == 0)
                        {
                            meshMod.setAction
                            (
                                polyModifyFace
                                (
                                    newFace,  // modified face
                                    faceI,    // label of face
                                    own,      // owner
                                    -1,       // neighbour
                                    false,    // face flip
                                    patchI,   // patch for face
                                    false,    // remove from zone
                                    zoneID,   // zone for face
                                    zoneFlip  // face flip in zone
                                 )
                             );
                        }
                        else
                        {
                            meshMod.setAction
                            (
                                polyAddFace
                                (
                                    newFace,  // vertices
                                    own,      // owner,
                                    -1,       // neighbour,
                                    -1,       // masterPointID,
                                    -1,       // masterEdgeID,
                                    faceI,    // masterFaceID,
                                    false,    // flipFaceFlux,
                                    patchI,   // patchID,
                                    zoneID,   // zoneID,
                                    zoneFlip  // zoneFlip
                                 )
                             );
                        }


                        if (f.fcIndex(next) == start)
                        {
                            break;
                        }
                        nAdded++;
                    }
                }
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
}


void Foam::autoHexDualiser::updateGapCells()
{
    const fvMesh& mesh = meshRefiner_.mesh();

    //update gap cells field if present
    if (mesh.foundObject<volScalarField>("gapCells"))
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        volScalarField& gapCells = const_cast<volScalarField&>
            (mesh.lookupObject<volScalarField>("gapCells"));

        boolList boundaryPts(mesh.nPoints(), false);
        forAll(mesh.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);

            if (patchi != -1 && !patches[patchi].coupled())
            {
                const labelList& f = mesh.faces()[facei];
                forAll(f,fp)
                {
                    boundaryPts[f[fp]] =  true;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            boundaryPts,
            orEqOp<bool>(),
            false           // null value
        );

        //look for purely internal faces gap faces
        scalarField neiGapCells;
        syncTools::swapBoundaryCellList(mesh, gapCells, neiGapCells);

        const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));
        boolList internalPts(mesh.nPoints(), false);

        forAll(mesh.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);

            if (patchi != -1 && !patches[patchi].coupled())
            {
                continue;
            }

            label own = mesh.faceOwner()[facei];
            scalar gapOwn = gapCells[own];
            scalar gapNei = 0;

            if (patchi == -1)
            {
                label nei = mesh.faceNeighbour()[facei];
                gapNei = gapCells[nei];
            }
            else if (patches[patchi].coupled())
            {
                gapNei = neiGapCells[facei-mesh.nInternalFaces()];
            }

            if (gapOwn > -1 && gapNei > -1)
            {
                const face f = mesh.faces()[facei];
                bool foundInternal = true;
                forAll(f,fp)
                {
                    if (boundaryPts[f[fp]])
                    {
                        foundInternal = false;
                        break;
                    }
                }
                if (foundInternal)
                {
                    const labelList& ownPts = mesh.cellPoints()[own];
                    forAll(ownPts, oPI)
                    {
                        internalPts[ownPts[oPI]] = true;
                    }

                    if (patchi == -1)
                    {
                        label nei = mesh.faceNeighbour()[facei];
                        const labelList& neiPts = mesh.cellPoints()[nei];
                        forAll(neiPts, nPI)
                        {
                            internalPts[neiPts[nPI]] = true;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            internalPts,
            orEqOp<bool>(),
            false           // null value
        );

        DynamicList<label> internalGapFaces(mesh.nFaces()/10);
        forAll(mesh.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            if
            (
                !isMasterFace[facei]
                || (patchi != -1 && !patches[patchi].coupled())
            )
            {
                continue;
            }

            label own = mesh.faceOwner()[facei];
            scalar gapOwn = gapCells[own];
            scalar gapNei = 0;

            if (patchi == -1)
            {
                label nei = mesh.faceNeighbour()[facei];
                gapNei = gapCells[nei];
            }
            else if (patches[patchi].coupled())
            {
                gapNei = neiGapCells[facei-mesh.nInternalFaces()];
            }

            if (gapOwn > -1 && gapNei > -1)
            {
                const face f = mesh.faces()[facei];
                bool foundInternal = true;
                forAll(f,fp)
                {
                    label pointi = f[fp];

                    if (internalPts[pointi] && boundaryPts[pointi])
                    {
                        foundInternal = false;
                        break;
                    }
                }
                if (foundInternal)
                {
                    internalGapFaces.append(facei);
                }
            }
        }

        labelList gapCellGlobalIndex(mesh.nFaces(), -1);
        globalIndex globalGapCells(internalGapFaces.size());
        label offset = globalGapCells.offset(Pstream::myProcNo());

        forAll(internalGapFaces, i)
        {
            label facei = internalGapFaces[i];
            gapCellGlobalIndex[facei] = i + offset;
        }

        syncTools::syncFaceList
        (
            mesh,
            gapCellGlobalIndex,
            maxEqOp<label>()
        );

        forAll(mesh.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            label gIndex = gapCellGlobalIndex[facei];
            if
            (
                (patchi != -1 && !patches[patchi].coupled()) || gIndex == -1
            )
            {
                continue;
            }

            label own = mesh.faceOwner()[facei];
            gapCells[own] = gIndex;

            if (patchi == -1)
            {
                label nei = mesh.faceNeighbour()[facei];
                gapCells[nei] = gIndex;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::autoHexDualiser::autoHexDualiser(meshRefinement& meshRefiner)
:
    meshRefiner_(meshRefiner)
{
    fvMesh& mesh = meshRefiner_.mesh();

    //Clear out on-demand addressing
    mesh.clearOut();

    //Create dual mesh
    calcDual();

    //Update gap cell field if present
    updateGapCells();

    splitBoundaryFaces();

    //Need to re-order points after dualisation
    polyTopoChange meshMod(mesh);
    // Change the mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh
    (
        mesh,
        false,      // inflate
        true,       // parallel sync
        false,       // cell ordering
        true       // point ordering
    );

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh.clearOut();
    }
    meshRefiner.updateMesh(map, labelList(0));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::autoHexDualiser::~autoHexDualiser()
{}


// ************************************************************************* //
