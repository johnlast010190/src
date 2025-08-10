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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "autoFanExtrude/autoFanExtrude.H"
#include "edgeClassification/edgeClassification.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addPatchCellLayer.H"
#include "regionSplit/localPointRegion.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(autoFanExtrude, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::autoFanExtrude::autoFanExtrude
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::autoFanExtrude::distribute
(
    boolList& marked
)
{
    fvMesh& mesh = meshRefiner_.mesh();

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
        distMap().distributeFaceData(marked);
    }
}


bool Foam::autoFanExtrude::filterAndSplit
(
    boolList& marked
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));
    labelList nPointBaffle(mesh.nPoints(), 0);

    forAll(mesh.points(), pointi)
    {
        const labelList& pFaces = mesh.pointFaces()[pointi];
        forAll(pFaces, pFI)
        {
            label facei = pFaces[pFI];
            if (isMasterFace.get(facei) && !marked[facei])
            {
                nPointBaffle[pointi]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nPointBaffle,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh.faces(), facei)
    {
        if (!marked[facei])
        {
            const face& f = mesh.faces()[facei];
            label nCorner = 0;
            label maxConnections = 0;

            forAll(f,fp)
            {
                label pointi = f[fp];
                maxConnections = max(maxConnections,nPointBaffle[pointi]);
                if (nPointBaffle[pointi] == 1)
                {
                    nCorner++;
                }
            }

            if (nCorner > 1 || maxConnections < 4)
            {
                marked[facei] = true;
            }
        }
    }

    label nMarked = 0;
    forAll(mesh.faces(), facei)
    {
        if (!marked[facei])
        {
            nMarked++;
        }
    }

    if (returnReduce(nMarked, sumOp<label>()) == 0)
    {
        return false;
    }

    if (layerParams_.fanTetSplit())
    {
        polyTopoChange meshMod(mesh);

        forAll(mesh.faces(), facei)
        {
            if (!marked[facei])
            {
                const face& f = mesh.faces()[facei];
                label nCorner = 0;
                label cornerPt = -1;

                forAll(f,fp)
                {
                    label pointi = f[fp];
                    if (nPointBaffle[pointi] == 1)
                    {
                        nCorner++;
                        cornerPt = pointi;
                    }
                }

                if (nCorner == 1 && f.size() == 4)
                {
                    label zoneID = mesh.faceZones().whichZone(facei);
                    bool zoneFlip = false;

                    label own = mesh.faceOwner()[facei];
                    label nei = mesh.faceNeighbour()[facei];

                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = mesh.faceZones()[zoneID];
                        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                    }

                    label nextFp = findIndex(f, cornerPt);

                    face f0(3);
                    face f1(3);

                    forAll(f0, fp)
                    {
                        if (fp > 0)
                        {
                            nextFp = f.fcIndex(nextFp);
                        }
                        f0[fp] = f[nextFp];
                    }

                    forAll(f1, fp)
                    {
                        f1[fp] = f[nextFp];
                        nextFp = f.fcIndex(nextFp);
                    }

                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            f0,        // modified face
                            facei,     // label of face
                            own,       // owner
                            nei,       // neighbour
                            false,     // face flip
                            -1,        // patch for face
                            false,     // remove from zone
                            zoneID,    // zone for face
                            zoneFlip   // face flip in zone
                        )
                    );

                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            f1,        // modified face
                            own,       // owner
                            nei,       // neighbour
                            -1,        // masterPointID
                            -1,        // masterEdgeID
                            facei,     // masterFaceID,
                            false,     // face flip
                            -1,        // patch for face
                            zoneID,    // zone for face
                            zoneFlip       // face flip in zone
                         )
                    );
                }
            }
        }

        // Change the mesh (no inflation, parallel sync)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

        // Update fields
        mesh.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh.clearOut();
        }

        meshRefiner_.updateMesh(map,labelList(0));

        meshRefinement::updateList
        (
            map().faceMap(),
            true,
            marked
        );
    }

    return true;
}


void Foam::autoFanExtrude::duplicate
(
    labelList& bafflePtType,
    scalarField& minExternalLen,
    List<labelPair>& baffles
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    boolList marked(mesh.nFaces(), false);
    forAll(baffles, bI)
    {
        marked[baffles[bI].first()] = true;
        marked[baffles[bI].second()] = true;
    }

    labelList nEdgeBaffleFaces(mesh.nEdges(), 0);
    forAll(mesh.edges(), edgei)
    {
        const labelList& eFaces = mesh.edgeFaces()[edgei];
        forAll(eFaces, eFI)
        {
            label facei = eFaces[eFI];
            if (marked[facei])
            {
                nEdgeBaffleFaces[edgei]++;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh,
        nEdgeBaffleFaces,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh.edges(), edgei)
    {
        const edge e = mesh.edges()[edgei];

        if
        (
            nEdgeBaffleFaces[edgei] == 0
        )
        {
            scalar eLen = e.mag(mesh.points());
            minExternalLen[e[0]] = max(minExternalLen[e[0]], eLen);
            minExternalLen[e[1]] = max(minExternalLen[e[1]], eLen);
            bafflePtType[e[0]] = max(bafflePtType[e[0]], label(-1));
            bafflePtType[e[1]] = max(bafflePtType[e[1]], label(-1));
        }
        else if (nEdgeBaffleFaces[edgei] == 4)
        {
            bafflePtType[e[0]] = max(bafflePtType[e[0]], label(0));
            bafflePtType[e[1]] = max(bafflePtType[e[1]], label(0));
        }
        else
        {
            bafflePtType[e[0]] = max(bafflePtType[e[0]], label(1));
            bafflePtType[e[1]] = max(bafflePtType[e[1]], label(1));
        }
    }
    syncTools::syncPointList
    (
        mesh,
        minExternalLen,
        maxEqOp<scalar>(),
        scalar(0)
    );
    syncTools::syncPointList
    (
        mesh,
        bafflePtType,
        maxEqOp<label>(),
        label(-1)
    );

    autoPtr<mapPolyMesh> mapPts = meshRefiner_.dupNonManifoldPoints();

    if (mapPts.valid())
    {
        // Update baffle numbering
        {
            const labelList& reverseFaceMap = mapPts().reverseFaceMap();
            forAll(baffles, i)
            {
                label f0 = reverseFaceMap[baffles[i].first()];
                label f1 = reverseFaceMap[baffles[i].second()];
                baffles[i] = labelPair(f0, f1);
            }
        }

        meshRefinement::updateList
        (
            mapPts().pointMap(),
            label(-1),
            bafflePtType
        );

        meshRefinement::updateList
        (
            mapPts().pointMap(),
            scalar(),
            minExternalLen
        );
    }
}


void Foam::autoFanExtrude::createBaffles
(
    const boolList& marked,
    List<labelPair>& baffles
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    DynamicList<label> extFaces(mesh.nFaces()/10);
    forAll(mesh.faces(), facei)
    {
        if (!marked[facei])
        {
            extFaces.append(facei);
        }
    }

    autoPtr<indirectPrimitivePatch> externalFaces
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), extFaces),
            mesh.points()
         )
    );

    if (debug)
    {
        simpleVTKWriter extFacesVTK
        (
            externalFaces().localFaces(),
            externalFaces().localPoints()
         );
        extFacesVTK.write("internalFaceLayers.vtk");
    }

    baffles.setSize(externalFaces().size());

    polyTopoChange meshMod(mesh);
    forAll(externalFaces(), i)
    {
        label facei = externalFaces().addressing()[i];

        label ownPatch = 0;
        label nbrPatch = 0;


        const face& f = mesh.faces()[facei];
        label zoneID = mesh.faceZones().whichZone(facei);
        bool zoneFlip = false;

        if (zoneID >= 0)
        {
            const faceZone& fZone = mesh.faceZones()[zoneID];
            zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
        }

        meshMod.setAction
        (
            polyModifyFace
            (
                f,                          // modified face
                facei,                      // label of face
                mesh.faceOwner()[facei],   // owner
                -1,                         // neighbour
                false,                      // face flip
                ownPatch,                   // patch for face
                false,                      // remove from zone
                zoneID,                     // zone for face
                zoneFlip                    // face flip in zone
             )
         );

        if (mesh.isInternalFace(facei))
        {
            if (nbrPatch == -1)
            {
                FatalErrorInFunction
                    << "No neighbour patch for internal face " << facei
                    << " fc:" << mesh.faceCentres()[facei]
                    << " ownPatch:" << ownPatch << abort(FatalError);
            }

            bool reverseFlip = false;
            if (zoneID >= 0)
            {
                reverseFlip = !zoneFlip;
            }

            meshMod.setAction
            (
                polyAddFace
                (
                    f.reverseFace(),            // modified face
                    mesh.faceNeighbour()[facei],// owner
                    -1,                         // neighbour
                    -1,                         // masterPointID
                    -1,                         // masterEdgeID
                    facei,                      // masterFaceID,
                    true,                       // face flip
                    nbrPatch,                   // patch for face
                    zoneID,                     // zone for face
                    reverseFlip                 // face flip in zone
                 )
             );
        }
    }
    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

    // Update fields
    mesh.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh.clearOut();
    }

    meshRefiner_.updateMesh(map,labelList(0));

    const labelList& faceMap = map().faceMap();
    const labelList& reverseFaceMap = map().reverseFaceMap();
    label baffleI = 0;
    forAll(mesh.faces(), facei)
    {
        label oldFaceI = faceMap[facei];
        label masterFaceI = reverseFaceMap[oldFaceI];

        if (masterFaceI != facei && !marked[oldFaceI])
        {
            baffles[baffleI] = labelPair(masterFaceI, facei);
            baffleI++;
        }
    }
}


bool Foam::autoFanExtrude::facesToExtrude
(
    const scalar& fanAngleCos,
    boolList& marked
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbours = mesh.faceNeighbour();

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    volScalarField& layerCells = const_cast<volScalarField&>
        (mesh.lookupObject<volScalarField>("layerStacks"));

    labelList neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());
    for
    (
        label faceI = mesh.nInternalFaces();
        faceI < mesh.nFaces();
        faceI++
     )
    {
        neiLayerCells[faceI-mesh.nInternalFaces()] =
            layerCells[owners[faceI]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLayerCells);

    labelList layerFaceType(mesh.nFaces(), -1);

    forAll(layerCells, celli)
    {
        if (layerCells[celli] > -1)
        {
            label cLevel = cellLevel[celli];
            const labelList& c = mesh.cells()[celli];
            forAll(c, i)
            {
                label facei = c[i];

                const face& f = mesh.faces()[facei];
                label nAnchors = 0;
                forAll(f,fp)
                {
                    if (pointLevel[f[fp]] <= cLevel)
                    {
                        nAnchors++;
                    }
                }

                label patchI = patches.whichPatch(facei);

                if (nAnchors == 3)
                {
                    layerFaceType[facei] = 0;
                }
                else if (patchI == -1 || patches[patchI].coupled())
                {
                    label nLayerCell = -1;
                    if (patchI == -1)
                    {
                        label nei = neighbours[facei];
                        if (nei == celli)
                        {
                            nei = owners[facei];
                        }
                        nLayerCell = layerCells[nei];
                    }
                    else
                    {
                        nLayerCell =
                            neiLayerCells[facei-mesh.nInternalFaces()];
                    }
                    if
                    (
                        layerCells[celli] == nLayerCell
                        || nLayerCell < 0
                    )
                    {
                        layerFaceType[facei] = 0;
                    }
                    else
                    {
                        layerFaceType[facei] = 1;
                    }
                }
                else
                {
                    layerFaceType[facei] = 0;
                }
            }
        }
    }

    const labelList& patchToNLayers = layerParams_.numLayers();
    DynamicList<label> adaptPatches(patches.size());

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && patchToNLayers[patchI] != -1)
        {
            adaptPatches.append(patchI);
        }
    }

    indirectPrimitivePatch grownPatch
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatches
        )
    );
    const labelList grownMeshEdges
    (
        grownPatch.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    boolList excludedFaces(grownPatch.size(), false);
    edgeClassification eClass
    (
        mesh,
        mesh.points(),
        grownPatch,
        grownMeshEdges,
        excludedFaces,
        0.93969,
        0.93969,
        fanAngleCos
    );
    const List<Tuple2<edgeClassification::edgeType,scalar>>&
        eType = eClass.edgeTypes();

    boolList markedEdge(mesh.nEdges(), false);
    forAll(grownMeshEdges, edgei)
    {
        label meshEdgeI = grownMeshEdges[edgei];
        if (eType[edgei].first() == edgeClassification::BAFFLE)
        {
            markedEdge[meshEdgeI] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        markedEdge,
        orEqOp<bool>(),
        false
    );

    PackedBoolList isMasterEdges = syncTools::getMasterEdges(mesh);

    //Check for variable number of layers and unmark grown edges
    {
        boolList failedEdges(mesh.nEdges(), false);
        boolList checkFaces(mesh.nFaces(), true);
        boolList marchedEdges = markedEdge;

        while (true)
        {
            label nChanged = 0;
            forAll(marchedEdges, edgei)
            {
                if (marchedEdges[edgei])
                {
                    labelHashSet edgePts(mesh.edges()[edgei]);

                    const labelList& eFaces = mesh.edgeFaces()[edgei];
                    forAll(eFaces, eFI)
                    {
                        label facei = eFaces[eFI];
                        if (layerFaceType[facei] == 1)
                        {
                            if (checkFaces[facei])
                            {
                                const labelList& fEdges =
                                    mesh.faceEdges()[facei];
                                if (fEdges.size() != 4)
                                {
                                    failedEdges[edgei] = true;
                                }
                                else
                                {
                                    forAll(fEdges, fEI)
                                    {
                                        label oppEdge = fEdges[fEI];
                                        if (oppEdge != edgei)
                                        {
                                            if (marchedEdges[oppEdge])
                                            {
                                                continue;
                                            }
                                            const edge& oe =
                                                mesh.edges()[oppEdge];
                                            label nShared = 0;

                                            forAll(oe, oei)
                                            {
                                                if (edgePts.found(oe[oei]))
                                                {
                                                    nShared++;
                                                }
                                            }
                                            if (nShared == 0)
                                            {
                                                nChanged++;
                                                checkFaces[facei] = false;
                                                marchedEdges[oppEdge] = true;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (returnReduce(nChanged, sumOp<label>()) == 0)
            {
                break;
            }
            syncTools::syncEdgeList
            (
                mesh,
                failedEdges,
                orEqOp<bool>(),
                false
             );
        }


        checkFaces = true;
        while (true)
        {
            label nChanged = 0;
            forAll(failedEdges, edgei)
            {
                if (failedEdges[edgei])
                {
                    labelHashSet edgePts(mesh.edges()[edgei]);

                    const labelList& eFaces = mesh.edgeFaces()[edgei];
                    forAll(eFaces, eFI)
                    {
                        label facei = eFaces[eFI];
                        if (layerFaceType[facei] == 1)
                        {
                            if (checkFaces[facei])
                            {
                                const labelList& fEdges =
                                    mesh.faceEdges()[facei];
                                if (fEdges.size() != 4)
                                {
                                    continue;
                                }
                                else
                                {
                                    forAll(fEdges, fEI)
                                    {
                                        label oppEdge = fEdges[fEI];
                                        if (oppEdge != edgei)
                                        {
                                            if (failedEdges[oppEdge])
                                            {
                                                continue;
                                            }
                                            const edge& oe =
                                                mesh.edges()[oppEdge];
                                            label nShared = 0;

                                            forAll(oe, oei)
                                            {
                                                if (edgePts.found(oe[oei]))
                                                {
                                                    nShared++;
                                                }
                                            }
                                            if (nShared == 0)
                                            {
                                                nChanged++;
                                                checkFaces[facei] = false;
                                                failedEdges[oppEdge] = true;
                                                if (markedEdge[oppEdge])
                                                {
                                                    markedEdge[oppEdge] = false;
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (returnReduce(nChanged, sumOp<label>()) == 0)
            {
                break;
            }

            syncTools::syncEdgeList
            (
                mesh,
                failedEdges,
                orEqOp<bool>(),
                false
             );

            syncTools::syncEdgeList
            (
                mesh,
                markedEdge,
                andEqOp<bool>(),
                false
             );
        }
    }

    //Filter out small feature edges
    while (true)
    {
        label nReset = 0;

        syncTools::syncEdgeList
        (
            mesh,
            markedEdge,
            orEqOp<bool>(),
            false
        );

        labelList nPtFtrEdges(mesh.nPoints(),0);

        forAll(nPtFtrEdges, pointi)
        {
            const labelList& pEdges = mesh.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if (isMasterEdges.get(edgei) == 1 &&  markedEdge[edgei])
                {
                    nPtFtrEdges[pointi]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nPtFtrEdges,
            plusEqOp<label>(),
            label(0)
        );

        labelList nSingleConnections(mesh.nPoints(), 0);
        forAll(nSingleConnections, pointi)
        {
            if (nPtFtrEdges[pointi] == 0)
            {
                continue;
            }
            const labelList& pEdges = mesh.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if (isMasterEdges.get(edgei) == 1 &&  markedEdge[edgei])
                {
                    const edge e = mesh.edges()[edgei];
                    label otherPt(e[0] == pointi ? e[1] : e[0]);
                    if (nPtFtrEdges[otherPt] == 1)
                    {
                        nSingleConnections[pointi]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nSingleConnections,
            plusEqOp<label>(),
            label(0)
         );

        forAll(mesh.edges(), edgei)
        {
            if (markedEdge[edgei])
            {
                const edge e = mesh.edges()[edgei];
                bool excludeEdge = false;
                if
                (
                    nSingleConnections[e[0]] == 1
                    && nSingleConnections[e[1]] == 1
                )
                {
                    excludeEdge = true;
                }
                else if
                (
                    nSingleConnections[e[0]] > 1
                    ||  nSingleConnections[e[1]] > 1
                )
                {
                    excludeEdge = true;
                }
                if (excludeEdge)
                {
                    markedEdge[edgei] = false;
                    nReset++;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            markedEdge,
            orEqOp<bool>(),
            false
        );

        forAll(mesh.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            if (patchi == -1 || patches[patchi].coupled())
            {
                continue;
            }

            const face& f = mesh.faces()[facei];
            forAll(f,fp)
            {
                label pointi = f[fp];
                if (nPtFtrEdges[pointi] == 2)
                {
                    label nextFp = f.fcIndex(fp);
                    label prevFp = f.rcIndex(fp);
                    label nextEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[pointi],
                        pointi,
                        f[nextFp]
                    );

                    label prevEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[pointi],
                        pointi,
                        f[prevFp]
                    );
                    if (markedEdge[nextEdge] && markedEdge[prevEdge])
                    {
                        point pt = mesh.points()[pointi];
                        vector nextVec = mesh.points()[f[nextFp]] - pt;
                        vector prevVec = pt - mesh.points()[f[prevFp]];
                        nextVec /= Foam::mag(nextVec) + VSMALL;
                        prevVec /= Foam::mag(prevVec) + VSMALL;
                        if ((nextVec&prevVec) < 0)
                        {
                            markedEdge[nextEdge] = false;
                            nReset++;
                            markedEdge[prevEdge] = false;
                            nReset++;
                        }
                    }
                }
            }
        }

        if (returnReduce(nReset, sumOp<label>()) == 0)
        {
            break;
        }

        syncTools::syncEdgeList
        (
            mesh,
            markedEdge,
            andEqOp<bool>(),
            false
        );
    }

    bool extrudedFaces = false;
    while (true)
    {
        syncTools::syncEdgeList
        (
            mesh,
            markedEdge,
            orEqOp<bool>(),
            false
        );
        syncTools::syncFaceList
        (
            mesh,
            marked,
            andEqOp<bool>()
        );

        label nUpdatedEdges = 0;
        forAll(markedEdge, edgei)
        {
            if (markedEdge[edgei])
            {
                labelHashSet edgePts(mesh.edges()[edgei]);

                const labelList& eFaces = mesh.edgeFaces()[edgei];
                forAll(eFaces, eFI)
                {
                    label facei = eFaces[eFI];
                    if (layerFaceType[facei] == 1)
                    {
                        if (marked[facei])
                        {
                            const labelList& fEdges = mesh.faceEdges()[facei];
                            forAll(fEdges, fEI)
                            {
                                label oppEdge = fEdges[fEI];
                                if (oppEdge != edgei)
                                {
                                    if (markedEdge[oppEdge])
                                    {
                                        continue;
                                    }
                                    const edge& oe = mesh.edges()[oppEdge];
                                    label nShared = 0;

                                    forAll(oe, oei)
                                    {
                                        if (edgePts.found(oe[oei]))
                                        {
                                            nShared++;
                                        }
                                    }
                                    if (nShared == 0)
                                    {
                                        marked[facei] = false;
                                        markedEdge[oppEdge] = true;
                                        nUpdatedEdges++;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (returnReduce(nUpdatedEdges, sumOp<label>()) == 0)
        {
            break;
        }
        else
        {
            extrudedFaces = true;
        }
    }

    return extrudedFaces;
}


void Foam::autoFanExtrude::move
(
    const indirectPrimitivePatch& bafflePP,
    const labelList& bafflePtType,
    const scalarField& minExternalLen,
    vectorField& disp
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    const labelList baffleMeshEdges
    (
        bafflePP.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    //Calculate baffle normals for projection
    boolList excludedFaces(bafflePP.size(), false);
    edgeClassification eClass
    (
        mesh,
        mesh.points(),
        bafflePP,
        baffleMeshEdges,
        excludedFaces,
        0.766,
        0.766
    );
    const pointField unsmoothedPointNormals =
        eClass.calculatePointNormals(excludedFaces, 0, true);

    pointField sNormals(mesh.nPoints(), vector::zero);
    forAll(bafflePP.meshPoints(), pti)
    {
        label meshPointI = bafflePP.meshPoints()[pti];
        sNormals[meshPointI] = unsmoothedPointNormals[pti];
    }
    syncTools::syncPointList
    (
        mesh,
        sNormals,
        maxMagSqrEqOp<point>(),
        vector::zero
    );

    scalarField minEdgeLength(mesh.nPoints(), GREAT);
    forAll(baffleMeshEdges, edgei)
    {
        label meshEdgeI = baffleMeshEdges[edgei];
        const edge& e = mesh.edges()[meshEdgeI];
        scalar eLen = e.mag(mesh.points());

        minEdgeLength[e[0]] = min(minEdgeLength[e[0]],eLen);
        minEdgeLength[e[1]] = min(minEdgeLength[e[1]],eLen);
    }
    syncTools::syncPointList
    (
        mesh,
        minEdgeLength,
        minEqOp<scalar>(),
        GREAT
    );


    forAll(mesh.points(), pointi)
    {
        if (bafflePtType[pointi] == 0)
        {
            scalar nLen = 0.2*minExternalLen[pointi];
            scalar projDist = min(nLen,0.5*minEdgeLength[pointi]);
            disp[pointi] = -projDist*sNormals[pointi];
        }
    }

    syncTools::syncPointList
    (
        mesh,
        disp,
        maxMagSqrEqOp<point>(),
        vector::zero
    );

    pointField newPoints = mesh.points();
    newPoints += disp;
    mesh.movePoints(newPoints);

}


Foam::autoPtr<Foam::mapPolyMesh> Foam::autoFanExtrude::extrudeLayer
(
    const indirectPrimitivePatch& bafflePP,
    const labelList& bafflePtType,
    const vectorField& disp
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches)
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            bafflePP
        )
    );

    labelList edgePatchID;
    labelList edgeZoneID;
    boolList edgeFlip;
    labelList inflateFaceID;
    label nPatches;
    Map<label> nbrProcToPatch;
    Map<label> patchToNbrProc;
    addPatchCellLayer::calcExtrudeInfo
    (
        true,   // for internal edges get zone info from any face

        mesh,
        globalFaces,
        edgeGlobalFaces,
        bafflePP,
        labelList(0), //list of grown-up patches

        edgePatchID,
        nPatches,
        nbrProcToPatch,
        patchToNbrProc,

        edgeZoneID,
        edgeFlip,
        inflateFaceID
    );

    label nAdded = nPatches - mesh.boundaryMesh().size();
    reduce(nAdded, sumOp<label>());

    Info<< "Adding overall " << nAdded << " processor patches." << endl;

    if (nAdded > 0)
    {
        DynamicList<polyPatch*> newPatches(nPatches);
        forAll(mesh.boundaryMesh(), patchi)
        {
            newPatches.append
            (
                mesh.boundaryMesh()[patchi].clone
                (
                    mesh.boundaryMesh()
                ).ptr()
            );
        }

        for
        (
            label patchi = mesh.boundaryMesh().size();
            patchi < nPatches;
            patchi++
        )
        {
            label nbrProci = patchToNbrProc[patchi];
            word name
            (
                processorPolyPatch::newName(Pstream::myProcNo(), nbrProci)
            );

            dictionary patchDict;
            patchDict.add("type", processorPolyPatch::typeName);
            patchDict.add("myProcNo", Pstream::myProcNo());
            patchDict.add("neighbProcNo", nbrProci);
            patchDict.add("nFaces", 0);
            patchDict.add("startFace", mesh.nFaces());

            meshRefiner_.appendPatch
            (
                mesh,
                name,
                patchDict
            );
        }
        mesh.clearOut();
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh()).updateMesh();
    }

    label nLayers = 1;
    // Layers per face
    labelList nFaceLayers(bafflePP.size(), nLayers);

    // Layers per point
    labelList nPointLayers(bafflePP.nPoints(), 0);
    scalarField ratio(bafflePP.nPoints(), 1.0);
    pointField ppDisp(bafflePP.nPoints(),vector::zero);

    forAll(bafflePP.meshPoints(), pointi)
    {
        label meshPointI = bafflePP.meshPoints()[pointi];
        if (bafflePtType[meshPointI] == 0)
        {
            nPointLayers[pointi] = nLayers;
            ppDisp[pointi] = -disp[meshPointI];
        }
    }

    // Topo change container.
    autoPtr<polyTopoChange> extrudeMod
    (
        new polyTopoChange(mesh)
    );
    addPatchCellLayer layerExtrude(mesh, true);

    labelList exposedPatchID(0);

    layerExtrude.setRefinement
    (
        globalFaces,
        edgeGlobalFaces,

        ratio,              // expansion ratio
        bafflePP,       // patch faces to extrude

        edgePatchID,        // if boundary edge: patch for extruded face
        edgeZoneID,         // optional zone for extruded face
        edgeFlip,
        inflateFaceID,      // mesh face that zone/patch info is from

        exposedPatchID,     // if new mesh: patches for exposed faces
        nFaceLayers,
        nPointLayers,
        ppDisp,
        true,
        extrudeMod()
    );

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> extrudeMap = extrudeMod().changeMesh(mesh, false);

    // Update fields
    mesh.updateMesh(extrudeMap);

    // Optionally inflate mesh
    if (extrudeMap().hasMotionPoints())
    {
        mesh.movePoints(extrudeMap().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh.setInstance(meshRefiner_.timeName());

    // Update intersection info
    meshRefiner_.updateMesh(extrudeMap, labelList(0));

    return extrudeMap;
}


void Foam::autoFanExtrude::mergeBaffles
(
    const mapPolyMesh& extrudeMap,
    List<labelPair>& baffles
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Collect candidate if
    //  - point on boundary face originating from baffle
    //  - and point originating from duplicate

    // Estimate number of points-to-be-merged
    DynamicList<label> candidates(baffles.size()*4);
    // Mark whether old face was on baffle
    PackedBoolList oldBaffleFace(extrudeMap.nOldFaces());
    forAll(baffles, i)
    {
        const labelPair& baffle = baffles[i];
        oldBaffleFace[baffle[0]] = true;
        oldBaffleFace[baffle[1]] = true;
    }

    boolList checkPts(mesh.nPoints(), false);
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        label oldFacei = extrudeMap.faceMap()[facei];
        if (oldFacei != -1 && oldBaffleFace[oldFacei])
        {
            const face& f = mesh.faces()[facei];
            forAll(f, fp)
            {
                label pointi = f[fp];
                checkPts[pointi] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        checkPts,
        orEqOp<bool>(),
        false
    );

    forAll(checkPts, pointi)
    {
        if (checkPts[pointi])
        {
            candidates.append(pointi);
        }
    }

    scalar minEdgeLength = GREAT;
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        label oldFacei = extrudeMap.faceMap()[facei];
        if (oldFacei != -1 && oldBaffleFace[oldFacei])
        {
            const labelList& fEdges = mesh.faceEdges()[facei];
            forAll(fEdges, fei)
            {
                label edgei = fEdges[fei];
                const edge& e = mesh.edges()[edgei];
                scalar eLen = e.mag(mesh.points());
                minEdgeLength = min(eLen, minEdgeLength);
            }
        }
    }
    reduce(minEdgeLength, minOp<scalar>());

    scalar mergeTol =
        min(0.5*minEdgeLength, meshRefiner_.mergeDistance());
    labelList oldToNew;
    label nNew = mergePoints
    (
        pointField(mesh.points(), candidates),
        mergeTol,
        false,
        oldToNew
    );

    // Extract points to be merged (i.e. multiple points originating
    // from a single one)

    labelListList newToOld(invertOneToMany(nNew, oldToNew));


    // Extract points with more than one old one
    labelList pointToMaster(mesh.nPoints(),-1);
    forAll(newToOld, newi)
    {
        const labelList& oldPoints = newToOld[newi];
        if (oldPoints.size() > 1)
        {
            labelList meshPoints
            (
                UIndirectList<label>(candidates, oldPoints)
            );
            label masteri = min(meshPoints);
            forAll(meshPoints, i)
            {
                pointToMaster[meshPoints[i]] = masteri;
            }
        }
    }

    // Count duplicate points
    label nPointPairs = 0;
    forAll(pointToMaster, pointi)
    {
        label otherPointi = pointToMaster[pointi];
        if (otherPointi != -1)
        {
            nPointPairs++;
        }
    }
    reduce(nPointPairs, sumOp<label>());
    if (nPointPairs > 0)
    {
        // Merge any duplicated points
        Info<< "Merging " << nPointPairs << " duplicated points ..." << endl;

        autoPtr<mapPolyMesh> map = meshRefiner_.mergePoints(pointToMaster);
    }

    autoPtr<mapPolyMesh> mapMergePtr = meshRefiner_.mergeBaffles
    (
        localPointRegion::findDuplicateFacePairs(mesh),
        Map<label>(0),
        false, //perform faceZone checks for new faces
        false
    );
}



void Foam::autoFanExtrude::setRefinement()
{
    const scalar& fanAngle = layerParams_.fanAngle();
    bool tryBaffleCollapse = layerParams_.tryBaffleCollapse();
    bool tryCornerCollapse = layerParams_.tryCornerCollapse();

    if
    (
        fanAngle < 270 || fanAngle > 360
        || tryBaffleCollapse || tryCornerCollapse
    )
    {
        return;
    }
    scalar fanAngleCos = -Foam::cos(degToRad(fanAngle));

    Info<< nl << "Extruding fan layer at internal faces"<<  nl <<endl;

    fvMesh& mesh = meshRefiner_.mesh();

    boolList marked(mesh.nFaces(), true);
    //Mark faces to extrude layer mesh
    if (!facesToExtrude(fanAngleCos,marked))
    {
        Info<<"No faces to extrude fan layers on. "<<endl;
        return;
    }

    //Ensure marked faces sit on same processor
    distribute(marked);

    //Filter out faces that cannot be extruded
    if (!filterAndSplit(marked))
    {
        Info<<"After filtering faces, no extrude fan mesh to generate. "<<endl;
        return;
    }

    List<labelPair> baffles;
    //Create baffles
    createBaffles(marked,baffles);

    //Duplicate baffle points
    labelList bafflePtType(mesh.nPoints(), -1);
    scalarField minExternalLen(mesh.nPoints(), GREAT);
    duplicate(bafflePtType,minExternalLen,baffles);

    List<label> baffleFaces(2*baffles.size());
    label baffleI = 0;
    forAll(baffles, i)
    {
        baffleFaces[baffleI++] = baffles[i].first();
        baffleFaces[baffleI++] = baffles[i].second();
    }

    autoPtr<indirectPrimitivePatch> bafflePP
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), baffleFaces),
            mesh.points()
         )
    );

    vectorField disp(mesh.nPoints(), vector::zero);

    //Move surface points to create void for layer mesh
    move(bafflePP(), bafflePtType, minExternalLen, disp);

    //Create layer in void
    autoPtr<mapPolyMesh> extrudeMap = extrudeLayer
    (
        bafflePP(),
        bafflePtType,
        disp
    );

    //Final merge of baffle faces
    mergeBaffles(extrudeMap(), baffles);

    return;
}


// ************************************************************************* //
