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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

Description
    All to do with adding cell layers

\*----------------------------------------------------------------------------*/
#include "snappyHexMeshDriver/snappyLayerDriver.H"
#include "snappyHexMeshDriver/snappySnapDriver.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "meshRefinement/meshRefinement.H"
#include "polyTopoChange/removePoints/removePoints.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "motionSmoother/motionSmoother.H"
#include "global/unitConversion/unitConversion.H"
#include "sets/topoSets/pointSet.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/cellSet.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "polyTopoChange/polyTopoChange/addPatchCellLayer.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "surfaceFormats/obj/OBJstream.H"
#include "snappyHexMeshDriver/layerParameters/layerParameters.H"
#include "polyTopoChange/combineFaces/combineFaces.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "fields/Fields/DynamicField/DynamicField.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "meshes/meshTools/mergePoints.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "fields/pointPatchFields/derived/slip/slipPointPatchFields.H"
#include "fields/pointPatchFields/derived/fixedNormalSlip/fixedNormalSlipPointPatchField.H"
#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchFields.H"
#include "externalDisplacementMeshMover/zeroFixedValue/zeroFixedValuePointPatchFields.H"
#include "fields/pointPatchFields/basic/calculated/calculatedPointPatchFields.H"
#include "fields/pointPatchFields/constraint/cyclicSlip/cyclicSlipPointPatchFields.H"
#include "fields/fvPatchFields/basic/fixedValue/fixedValueFvPatchFields.H"
#include "regionSplit/localPointRegion.H"
#include "externalDisplacementMeshMover/externalDisplacementMeshMover.H"
#include "fields/Fields/scalarField/scalarIOField.H"
#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(snappyLayerDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::label Foam::snappyLayerDriver::stopNonConsecutiveExtrusion
(
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,
    const scalarField& minThickness,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    label nChanged = 0;
    while (true)
    {
        syncPatchDisplacement
        (
            pp,
            minThickness,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Make sure that a face doesn't have two non-consecutive areas not
        // extruded (e.g. quad where vertex 0 and 2 are not extruded but
        // 1 and 3 are) since this gives topological  errors.

        label nPinched = 0;
        const faceList& localFaces = pp.localFaces();
        forAll(localFaces, i)
        {
            const face& localF = localFaces[i];

            // Count number of transitions from unsnapped to snapped.
            label nTrans = 0;

            extrudeMode prevMode = extrudeStatus[localF.prevLabel(0)];

            forAll(localF, fp)
            {
                extrudeMode fpMode = extrudeStatus[localF[fp]];

                if (prevMode == NOEXTRUDE && fpMode != NOEXTRUDE)
                {
                    nTrans++;
                }
                prevMode = fpMode;
            }

            if (nTrans > 1)
            {
                // Multiple pinches. Reset whole face as unextruded.
                if
                (
                    unmarkExtrusion
                    (
                        localF,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nPinched++;
                    nChanged++;
                }
            }
        }

        reduce(nPinched, sumOp<label>());

        Info<< "stopNonConsecutiveExtrusion : Unextruded " << nPinched
            << " faces due to non-consecutive vertices being extruded." << endl;

        // Butterfly
        // ~~~~~~~~~

        // Make sure that a string of edges becomes a single face so
        // not a butterfly. Occassionally an 'edge' will have a single dangling
        // vertex due to face combining. These get extruded as a single face
        // (with a dangling vertex) so make sure this extrusion forms a single
        // shape.
        //  - continuous i.e. no butterfly:
        //      +     +
        //      |\   /|
        //      | \ / |
        //      +--+--+
        //  - extrudes from all but the endpoints i.e. no partial
        //    extrude
        //            +
        //           /|
        //          / |
        //      +--+--+
        // The common error topology is a pinch somewhere in the middle
        label nButterFly = 0;
        {
            DynamicList<label> stringedVerts;
            forAll(pp.edges(), edgeI)
            {
                const labelList& globFaces = edgeGlobalFaces[edgeI];

                if (globFaces.size() == 2)
                {
                    label myFaceI = pp.edgeFaces()[edgeI][0];
                    label myGlobalFaceI = globalFaces.toGlobal
                    (
                        pp.addressing()[myFaceI]
                    );
                    label nbrGlobalFaceI =
                    (
                        globFaces[0] != myGlobalFaceI
                      ? globFaces[0]
                      : globFaces[1]
                    );
                    getVertexString
                    (
                        pp,
                        edgeGlobalFaces,
                        myFaceI,
                        edgeI,
                        myGlobalFaceI,
                        nbrGlobalFaceI,
                        stringedVerts
                    );

                    if
                    (
                        extrudeStatus[stringedVerts[0]] != NOEXTRUDE
                     || extrudeStatus[stringedVerts.last()] != NOEXTRUDE
                    )
                    {
                        // Any pinch in the middle
                        bool pinch = false;
                        for (label i = 1; i < stringedVerts.size()-1; i++)
                        {
                            if
                            (
                                extrudeStatus[stringedVerts[i]] == NOEXTRUDE
                            )
                            {
                                pinch = true;
                                break;
                            }
                        }
                        if (pinch)
                        {
                            forAll(stringedVerts, i)
                            {
                                if
                                (
                                    unmarkExtrusion
                                    (
                                        stringedVerts[i],
                                        patchDisp,
                                        patchNLayers,
                                        extrudeStatus
                                    )
                                )
                                {
                                    nButterFly++;
                                    nChanged++;
                                }
                            }
                        }
                    }
                }
            }
        }

        reduce(nButterFly, sumOp<label>());

        Info<< "stopNonConsecutiveExtrusion : Unextruded " << nButterFly
            << " faces due to stringed edges with inconsistent extrusion."
            << endl;

        // Make sure that a face has consistent number of layers for all
        // its vertices.

        if (nPinched + nButterFly == 0)
        {
            break;
        }
    }
    return nChanged;
}


Foam::List<Foam::face> Foam::snappyLayerDriver::calcSplitFace
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const label& patchFaceI,
    vectorField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus,
    PackedList<1>& isTriSplit,
    label& nChanged
)
{
    scalar minFaceArea = 0.15;

    const face& f = pp.localFaces()[patchFaceI];
    label masterFaceI = pp.addressing()[patchFaceI];

    DynamicList<face> updatedFaces(3);
    // splitting zone faces not currently supported
    if (f.size() > 3)
    {
        label nExtruded = 0;
        label fpStart = -1;
        label fpEnd = -1;

        extrudeMode prevMode = extrudeStatus[f.prevLabel(0)];

        forAll(f, fp)
        {
            extrudeMode fpMode = extrudeStatus[f[fp]];

            if (fpMode != NOEXTRUDE)
            {
                nExtruded++;
            }

            if (prevMode == NOEXTRUDE && fpMode != NOEXTRUDE)
            {
                fpStart = fp;
            }

            if (prevMode != NOEXTRUDE && fpMode == NOEXTRUDE)
            {
                fpEnd = f.rcIndex(fp);
            }

            prevMode = fpMode;
        }

        if (f.size() - nExtruded >= 3 && nExtruded > 0)
        {
            label iNext = f.fcIndex(fpEnd);
            label iPrev = f.rcIndex(fpStart);

            face baseFace(nExtruded + 2);
            label j = iPrev;
            forAll(baseFace, fp)
            {
                baseFace[fp] = pp.meshPoints()[f[j]];
                j = f.fcIndex(j);
            }
            scalar faceArea = mag(baseFace.areaNormal(mesh.points()));

            scalar areaRatio = faceArea /
                ( mesh.magFaceAreas()[masterFaceI] + SMALL);

            if (areaRatio >= minFaceArea && areaRatio < (1.-minFaceArea))
            {
                face terminatedFace(f.size() - nExtruded);

                j = iNext;
                forAll(terminatedFace, fp)
                {
                    terminatedFace[fp] = pp.meshPoints()[f[j]];
                    j = f.fcIndex(j);
                }

                updatedFaces.append(baseFace);
                updatedFaces.append(terminatedFace);
            }
            else if (areaRatio < minFaceArea)
            {
                // try splitting in the opposite direction
                face tetFace1(3);
                j = fpEnd;
                forAll(tetFace1, fp)
                {
                    tetFace1[fp] = pp.meshPoints()[f[j]];
                    j = f.fcIndex(j);
                }
                label jEnd = f.rcIndex(j);

                scalar faceArea1 = mag(tetFace1.areaNormal(mesh.points()));
                scalar areaRatio1 = faceArea1 /
                    ( mesh.magFaceAreas()[masterFaceI] + SMALL);

                face tetFace2(3);
                j = fpStart;
                for (label iter = 0; iter < 2; iter++)
                {
                    j = f.rcIndex(j);
                }
                label jStart = j;

                forAll(tetFace2, fp)
                {
                    tetFace2[fp] = pp.meshPoints()[f[j]];
                    j = f.fcIndex(j);
                }

                scalar faceArea2 = mag(tetFace2.areaNormal(mesh.points()));
                scalar areaRatio2 = faceArea2 /
                    ( mesh.magFaceAreas()[masterFaceI] + SMALL);

                scalar areaRatio3(1.0);

                face tetFace3(f.size()-2);

                if (nExtruded > 1)
                {
                    j = fpStart;
                    label i = 0;
                    for (label nExt = 0; nExt < nExtruded; nExt++)
                    {
                        tetFace3[i++] = pp.meshPoints()[f[j]];
                        j = f.fcIndex(j);
                    }
                    j = jEnd;
                    tetFace3[i++] = pp.meshPoints()[f[j]];

                    if (jStart != jEnd)
                    {
                        for
                        (
                            label iter = 1;
                            iter < tetFace3.size() - nExtruded;
                            iter++
                         )
                        {
                            j = f.fcIndex(j);
                            tetFace3[i++] = pp.meshPoints()[f[j]];
                        }
                    }
                    areaRatio3 = 1. - areaRatio2 - areaRatio1;
                }
                else
                {
                    if (jStart != jEnd)
                    {
                        j = jEnd;
                        tetFace3[0] = pp.meshPoints()[f[fpStart]];
                        for (label iter = 1; iter < tetFace3.size(); iter++)
                        {
                            tetFace3[iter] = pp.meshPoints()[f[j]];
                            j = f.fcIndex(j);
                        }

                        areaRatio3 = 1. - areaRatio2 - areaRatio1;
                    }
                }

                if
                (
                    areaRatio1 > minFaceArea
                    && areaRatio2 > minFaceArea
                    && areaRatio3 > minFaceArea
                 )
                {
                    isTriSplit.set(masterFaceI, 1);
                    updatedFaces.append(tetFace1);
                    updatedFaces.append(tetFace2);

                    if (jStart != jEnd || nExtruded > 1)
                    {
                        updatedFaces.append(tetFace3);
                    }
                }
                else
                {
                    if
                    (
                        unmarkExtrusion
                        (
                            pp.localFaces()[patchFaceI],
                            patchDisp,
                            patchNLayers,
                            extrudeStatus
                         )
                     )
                    {
                        nChanged++;
                    }
                }
            }
        }
    }

    return updatedFaces.shrink();
}


void Foam::snappyLayerDriver::splitLayerTerminationFaces
(
    const layerParameters& layerParams,
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,
    const scalarField& minThickness,
    vectorField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus,
    polyTopoChange& meshMod,
    PackedList<1>& isTriSplit
)
{
    stopNonConsecutiveExtrusion
    (
        globalFaces,
        edgeGlobalFaces,
        pp,
        minThickness,
        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    const fvMesh& mesh = meshRefiner_.mesh();

    boolList zoneFace(pp.size(), false);
    labelHashSet zonePatches(layerZonePatchIDs());
    Map<label> baffleMap(pp.size());
    forAll(pp, patchFaceI)
    {
        label meshFaceI = pp.addressing()[patchFaceI];
        label patchi = mesh.boundaryMesh().whichPatch(meshFaceI);
        if (zonePatches.found(patchi))
        {
            zoneFace[patchFaceI] =  true;
            baffleMap.insert(meshFaceI,patchFaceI);
        }
    }

    label nChanged = 0;
    List<List<face>> updatedFaces(pp.size(), List<face>(0));

    forAll(pp, patchFaceI)
    {
        updatedFaces[patchFaceI] = calcSplitFace
        (
            mesh,
            pp,
            patchFaceI,
            patchDisp,
            patchNLayers,
            extrudeStatus,
            isTriSplit,
            nChanged
         );
    }

    //Check split faces to check splitting doesn't form a non-manifold edge
    const labelList& meshPoints = pp.meshPoints();
    const labelList meshEdges
    (
            pp.meshEdges(mesh.edges(), mesh.pointEdges())
    );
    const PackedBoolList isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh,
            meshEdges
        )
    );
    labelList nPPPointEdges(meshPoints.size(), 0);
    labelList nPPPointFacesSplit(meshPoints.size(), 0);

    forAll(meshPoints, ptI)
    {
        const labelList& pFaces = pp.pointFaces()[ptI];
        forAll(pFaces, pFI)
        {
            label facei = pFaces[pFI];
            if (updatedFaces[facei].size() > 0)
            {
                nPPPointFacesSplit[ptI]++;
            }
        }
        const labelList& pEdges = pp.pointEdges()[ptI];
        forAll(pEdges, pEI)
        {
            label edgei = pEdges[pEI];

            if (isPatchMasterEdge[edgei])
            {
                nPPPointEdges[ptI]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        nPPPointFacesSplit,
        plusEqOp<label>(),
        label(0)      // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        nPPPointEdges,
        plusEqOp<label>(),
        label(0)      // null value
    );

    forAll(meshPoints, ptI)
    {
        if
        (
            patchNLayers[ptI] == 0
            && nPPPointFacesSplit[ptI] == 2
            && nPPPointEdges[ptI] == 2
        )
        {
            const labelList& pFaces = pp.pointFaces()[ptI];
            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                if (updatedFaces[facei].size() > 0)
                {
                    updatedFaces[facei] = List<face>(0);
                    const face& lf = pp.localFaces()[facei];
                    unmarkExtrusion
                    (
                        lf,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                     );
                }
            }
        }
    }

    reduce(nChanged, sumOp<label>());
    Info<< "Prevented extrusion on " << nChanged
        << " points due to failure to split face." << nl << endl;

    forAll(pp, patchFaceI)
    {
        if (!zoneFace[patchFaceI] && updatedFaces[patchFaceI].size() > 0)
        {
            label masterFaceI = pp.addressing()[patchFaceI];

            updateCutFaces(updatedFaces[patchFaceI], masterFaceI, meshMod);
        }
    }

    stopNonConsecutiveExtrusion
    (
        globalFaces,
        edgeGlobalFaces,
        pp,
        minThickness,
        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    return;
}


void Foam::snappyLayerDriver::updateCutFaces
(
    const List<face>& updatedFaces,
    const label masterFaceI,
    polyTopoChange& meshMod
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    label own = mesh.faceOwner()[masterFaceI];

    label zoneID = mesh.faceZones().whichZone(masterFaceI);
    bool zoneFlip = false;
    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh.faceZones()[zoneID];
        zoneFlip = fZone.flipMap()[fZone.whichFace(masterFaceI)];
    }

    label patchi = mesh.boundaryMesh().whichPatch(masterFaceI);

    labelList newFaces(updatedFaces.size()-1, -1);

    // Modify the master face.
    meshMod.setAction
    (
        polyModifyFace
        (
            updatedFaces[0],       // original face
            masterFaceI,    // label of face
            own,            // owner
            -1,             // neighbour
            false,          // face flip
            patchi,         // patch for face
            false,          // remove from zone
            zoneID,         // zone for face
            zoneFlip        // face flip in zone
         )
     );

    for (label j=1; j < updatedFaces.size(); j++)
    {
        meshMod.setAction
        (
            polyAddFace
            (
                updatedFaces[j], // vertices
                own,            // owner,
                -1,             // neighbour,
                -1,             // masterPointID,
                -1,             // masterEdgeID,
                masterFaceI,    // masterFaceID,
                false,          // flipFaceFlux,
                patchi,         // patchID,
                zoneID,         // zoneID,
                zoneFlip        // zoneFlip
             )
         );
    }

    return;
}


Foam::List<Foam::edge> Foam::snappyLayerDriver::calculateFaceCuts
(
    const List<face>& updatedFaces
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    DynamicList<edge> cuts(updatedFaces.size()-1);

    forAll(updatedFaces, i)
    {
        face f = updatedFaces[i];
        label start = 0;
        forAll(f,fp)
        {
            label v0 = f[start];
            start = f.fcIndex(start);
            label v1 = f[start];
            label meshEdgeI = meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[v0],
                v0,
                v1
             );
            if (meshEdgeI == -1)
            {
                edge eCut(v0,v1);

                bool foundEdge = false;
                forAll(cuts, mcI)
                {
                    edge e = cuts[mcI];
                    if
                    (
                        (e[0] == eCut[0] && e[1] == eCut[1])
                        || (e[0] == eCut[1] && e[1] == eCut[0])
                     )
                    {
                        foundEdge = true;
                        break;
                    }
                }
                if (!foundEdge)
                {
                    cuts.append(eCut);
                }
            }
        }
    }
    return cuts.shrink();
}

// For debugging: Dump displacement to .obj files
void Foam::snappyLayerDriver::dumpDisplacement
(
    const fileName& prefix,
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp,
    const List<extrudeMode>& extrudeStatus
)
{
    OBJstream dispStr(prefix + "_disp.obj");
    Info<< "Writing all displacements to " << dispStr.name() << endl;

    forAll(patchDisp, patchPointi)
    {
        const point& pt = pp.localPoints()[patchPointi];
        dispStr.write(linePointRef(pt, pt + patchDisp[patchPointi]));
    }


    OBJstream illStr(prefix + "_illegal.obj");
    Info<< "Writing invalid displacements to " << illStr.name() << endl;

    forAll(patchDisp, patchPointi)
    {
        if (extrudeStatus[patchPointi] != EXTRUDE)
        {
            const point& pt = pp.localPoints()[patchPointi];
            illStr.write(linePointRef(pt, pt + patchDisp[patchPointi]));
        }
    }
}


Foam::tmp<Foam::scalarField> Foam::snappyLayerDriver::avgPointData
(
    const indirectPrimitivePatch& pp,
    const scalarField& pointFld
)
{
    tmp<scalarField> tfaceFld(new scalarField(pp.size(), 0.0));
    scalarField& faceFld = tfaceFld.ref();

    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];
        if (f.size())
        {
            forAll(f, fp)
            {
                faceFld[facei] += pointFld[f[fp]];
            }
            faceFld[facei] /= f.size();
        }
    }
    return tfaceFld;
}


// Check that primitivePatch is not multiply connected. Collect non-manifold
// points in pointSet.
void Foam::snappyLayerDriver::checkManifold
(
    const indirectPrimitivePatch& fp,
    pointSet& nonManifoldPoints
)
{
    // Check for non-manifold points (surface pinched at point)
    fp.checkPointManifold(false, &nonManifoldPoints);

    // Check for edge-faces (surface pinched at edge)
    const labelListList& edgeFaces = fp.edgeFaces();

    forAll(edgeFaces, edgei)
    {
        const labelList& eFaces = edgeFaces[edgei];

        if (eFaces.size() > 2)
        {
            const edge& e = fp.edges()[edgei];

            nonManifoldPoints.insert(fp.meshPoints()[e[0]]);
            nonManifoldPoints.insert(fp.meshPoints()[e[1]]);
        }
    }
}


void Foam::snappyLayerDriver::checkMeshManifold() const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Checking mesh manifoldness ..." << endl;

    // Get all outside faces
    labelList outsideFaces(mesh.nFaces() - mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
         outsideFaces[facei - mesh.nInternalFaces()] = facei;
    }

    pointSet nonManifoldPoints
    (
        mesh,
        "nonManifoldPoints",
        mesh.nPoints() / 100
    );

    // Build primitivePatch out of faces and check it for problems.
    checkManifold
    (
        indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), outsideFaces),
            mesh.points()
        ),
        nonManifoldPoints
    );

    label nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    if (nNonManif > 0)
    {
        Info<< "Outside of mesh is multiply connected across edges or"
            << " points." << nl
            << "This is not a fatal error but might cause some unexpected"
            << " behaviour." << nl
            //<< "Writing " << nNonManif
            //<< " points where this happens to pointSet "
            //<< nonManifoldPoints.name()
            << endl;
        if (debug)
        {
            nonManifoldPoints.instance() = meshRefiner_.timeName();
            nonManifoldPoints.write();
        }
    }
    Info<< endl;
}



// Unset extrusion on point. Returns true if anything unset.
bool Foam::snappyLayerDriver::unmarkExtrusion
(
    const label patchPointi,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    if (extrudeStatus[patchPointi] == EXTRUDE)
    {
        extrudeStatus[patchPointi] = NOEXTRUDE;
        patchNLayers[patchPointi] = 0;
        patchDisp[patchPointi] = Zero;
        return true;
    }
    else if (extrudeStatus[patchPointi] == EXTRUDEREMOVE)
    {
        extrudeStatus[patchPointi] = NOEXTRUDE;
        patchNLayers[patchPointi] = 0;
        patchDisp[patchPointi] = Zero;
        return true;
    }
    else
    {
        return false;
    }
}


// Unset extrusion on face. Returns true if anything unset.
bool Foam::snappyLayerDriver::unmarkExtrusion
(
    const face& localFace,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
)
{
    bool unextruded = false;

    forAll(localFace, fp)
    {
        if
        (
            unmarkExtrusion
            (
                localFace[fp],
                patchDisp,
                patchNLayers,
                extrudeStatus
            )
        )
        {
            unextruded = true;
        }
    }
    return unextruded;
}


Foam::label Foam::snappyLayerDriver::constrainFp(const label sz, const label fp)
{
    if (fp >= sz)
    {
        return 0;
    }
    else if (fp < 0)
    {
        return sz-1;
    }
    else
    {
        return fp;
    }
}


void Foam::snappyLayerDriver::countCommonPoints
(
    const indirectPrimitivePatch& pp,
    const label facei,

    Map<label>& nCommonPoints
) const
{
    const faceList& localFaces = pp.localFaces();
    const labelListList& pointFaces = pp.pointFaces();

    const face& f = localFaces[facei];

    nCommonPoints.clear();

    forAll(f, fp)
    {
        label pointi = f[fp];
        const labelList& pFaces = pointFaces[pointi];

        forAll(pFaces, pFacei)
        {
            label nbFacei = pFaces[pFacei];

            if (facei < nbFacei)
            {
                // Only check once for each combination of two faces.

                Map<label>::iterator fnd = nCommonPoints.find(nbFacei);

                if (fnd == nCommonPoints.end())
                {
                    // First common vertex found.
                    nCommonPoints.insert(nbFacei, 1);
                }
                else
                {
                    fnd()++;
                }
            }
        }
    }
}


bool Foam::snappyLayerDriver::checkCommonOrder
(
    const label nCommon,
    const face& curFace,
    const face& nbFace
) const
{
    forAll(curFace, fp)
    {
        // Get the index in the neighbouring face shared with curFace
        const label nb = findIndex(nbFace, curFace[fp]);

        if (nb != -1)
        {

            // Check the whole face from nb onwards for shared vertices
            // with neighbouring face. Rule is that any shared vertices
            // should be consecutive on both faces i.e. if they are
            // vertices fp,fp+1,fp+2 on one face they should be
            // vertices nb, nb+1, nb+2 (or nb+2, nb+1, nb) on the
            // other face.


            // Vertices before and after on curFace
            label fpPlus1 = curFace.fcIndex(fp);
            label fpMin1  = curFace.rcIndex(fp);

            // Vertices before and after on nbFace
            label nbPlus1 = nbFace.fcIndex(nb);
            label nbMin1  = nbFace.rcIndex(nb);

            // Find order of walking by comparing next points on both
            // faces.
            label curInc = labelMax;
            label nbInc = labelMax;

            if (nbFace[nbPlus1] == curFace[fpPlus1])
            {
                curInc = 1;
                nbInc = 1;
            }
            else if (nbFace[nbPlus1] == curFace[fpMin1])
            {
                curInc = -1;
                nbInc = 1;
            }
            else if (nbFace[nbMin1] == curFace[fpMin1])
            {
                curInc = -1;
                nbInc = -1;
            }
            else
            {
                curInc = 1;
                nbInc = -1;
            }


            // Pass1: loop until start of common vertices found.
            label curNb = nb;
            label curFp = fp;

            do
            {
                curFp = constrainFp(curFace.size(), curFp+curInc);
                curNb = constrainFp(nbFace.size(), curNb+nbInc);
            } while (curFace[curFp] == nbFace[curNb]);

            // Pass2: check equality walking from curFp, curNb
            // in opposite order.

            curInc = -curInc;
            nbInc = -nbInc;

            for (label commonI = 0; commonI < nCommon; commonI++)
            {
                curFp = constrainFp(curFace.size(), curFp+curInc);
                curNb = constrainFp(nbFace.size(), curNb+nbInc);

                if (curFace[curFp] != nbFace[curNb])
                {
                    // Error: gap in string of connected vertices
                    return false;
                }
            }

            // Done the curFace - nbFace combination.
            break;
        }
    }

    return true;
}


void Foam::snappyLayerDriver::checkCommonOrder
(
    const indirectPrimitivePatch& pp,
    const label facei,
    const Map<label>& nCommonPoints,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    forAllConstIter(Map<label>, nCommonPoints, iter)
    {
        label nbFacei = iter.key();
        label nCommon = iter();

        const face& curFace = pp[facei];
        const face& nbFace = pp[nbFacei];

        if
        (
            nCommon >= 2
         && nCommon != nbFace.size()
         && nCommon != curFace.size()
        )
        {
            bool stringOk = checkCommonOrder(nCommon, curFace, nbFace);

            if (!stringOk)
            {
                // Note: unmark whole face or just the common points?
                // For now unmark the whole face
                unmarkExtrusion
                (
                    pp.localFaces()[facei],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
                unmarkExtrusion
                (
                    pp.localFaces()[nbFacei],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
            }
        }
    }
}


void Foam::snappyLayerDriver::handleNonStringConnected
(
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    // Detect faces which are connected on non-consecutive vertices.
    // This is the "<<Number of faces with non-consecutive shared points"
    // warning from checkMesh. These faces cannot be extruded so
    // there is no need to even attempt it.

    List<extrudeMode> oldExtrudeStatus;
    autoPtr<OBJstream> str;
    if (debug&meshRefinement::LAYERINFO)
    {
        oldExtrudeStatus = extrudeStatus;
        str.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
               /"nonStringConnected.obj"
            )
        );
        Pout<< "Dumping string edges to " << str().name();
    }


    // 1) Local
    Map<label> nCommonPoints(100);

    forAll(pp, facei)
    {
        countCommonPoints(pp, facei, nCommonPoints);

        // Faces share pointi. Find any more shared points
        // and if not in single string unmark all. See
        // primitiveMesh::checkCommonOrder
        checkCommonOrder
        (
            pp,
            facei,
            nCommonPoints,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );
    }

    // 2) TDB. Other face remote



    if (debug&meshRefinement::LAYERINFO)
    {
        forAll(extrudeStatus, pointi)
        {
            if (extrudeStatus[pointi] != oldExtrudeStatus[pointi])
            {
                str().write
                (
                    meshRefiner_.mesh().points()[pp.meshPoints()[pointi]]
                );
            }
        }
    }
}


// No extrusion at non-manifold points.
void Foam::snappyLayerDriver::handleNonManifolds
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const labelListList& edgeGlobalFaces,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling non-manifold points ..." << endl;

    // Detect non-manifold points
    Info<< nl << "Checking patch manifoldness ..." << endl;

    pointSet nonManifoldPoints(mesh, "nonManifoldPoints", pp.nPoints());

    // 1. Local check
    checkManifold(pp, nonManifoldPoints);

    // 2. Remote check for boundary edges on coupled boundaries
    forAll(edgeGlobalFaces, edgei)
    {
        if
        (
            pp.edgeFaces()[edgei].size() == 1
         && edgeGlobalFaces[edgei].size() > 2
        )
        {
            // So boundary edges that are connected to more than 2 processors
            // i.e. a non-manifold edge which is exactly on a processor
            // boundary.
            const edge& e = pp.edges()[edgei];
            nonManifoldPoints.insert(pp.meshPoints()[e[0]]);
            nonManifoldPoints.insert(pp.meshPoints()[e[1]]);
        }
    }

    // 3. Remote check for end of layer across coupled boundaries
    {
        PackedBoolList isCoupledEdge(mesh.nEdges());

        const labelList& cpEdges = mesh.globalData().coupledPatchMeshEdges();
        forAll(cpEdges, i)
        {
            isCoupledEdge[cpEdges[i]] = true;
        }
        syncTools::syncEdgeList
        (
            mesh,
            isCoupledEdge,
            orEqOp<unsigned int>(),
            0
        );

        forAll(edgeGlobalFaces, edgei)
        {
            label meshEdgei = meshEdges[edgei];

            if
            (
                pp.edgeFaces()[edgei].size() == 1
             && edgeGlobalFaces[edgei].size() == 1
             && isCoupledEdge[meshEdgei]
            )
            {
                // Edge of patch but no continuation across processor.
                const edge& e = pp.edges()[edgei];
                //Pout<< "** Stopping extrusion on edge "
                //    << pp.localPoints()[e[0]]
                //    << pp.localPoints()[e[1]] << endl;
                nonManifoldPoints.insert(pp.meshPoints()[e[0]]);
                nonManifoldPoints.insert(pp.meshPoints()[e[1]]);
            }
        }
    }



    label nNonManif = returnReduce(nonManifoldPoints.size(), sumOp<label>());

    Info<< "Outside of local patch is multiply connected across edges or"
        << " points at " << nNonManif << " points." << endl;

    if (nNonManif > 0)
    {
        const labelList& meshPoints = pp.meshPoints();

        forAll(meshPoints, patchPointi)
        {
            if (nonManifoldPoints.found(meshPoints[patchPointi]))
            {
                unmarkExtrusion
                (
                    patchPointi,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
            }
        }
    }

    Info<< "Set displacement to zero for all " << nNonManif
        << " non-manifold points" << endl;
}


// No extrusion at grown up edges where angle is too large
void Foam::snappyLayerDriver::handleGrownUpEdges
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const labelList& grownUpIDs,
    const scalar grownUpAngleTerminateCos,

    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const labelList& meshPoints = pp.meshPoints();
    const vectorField& faceNormals = pp.faceNormals();

    pointField pointNormals(pp.nPoints(), vector::zero);
    {
        scalarField nPointFaces(pp.nPoints(), 0.);

        forAll(faceNormals, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            forAll(f, fp)
            {
                scalar fA = mesh.magFaceAreas()[pp.addressing()[faceI]];
                pointNormals[f[fp]] += faceNormals[faceI]*fA;
                nPointFaces[f[fp]] += fA;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            pointNormals,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nPointFaces,
            plusEqOp<scalar>(),
            scalar(0.)          // null value
        );

        forAll(pointNormals, i)
        {
            pointNormals[i] /= nPointFaces[i];
        }
    }
    pointNormals /= (mag(pointNormals) + SMALL);

    // Do not smooth on points marked as external
    boolList grownUpPoints(mesh.nPoints(),false);
    boolList isExternalEdge(mesh.nEdges(),false);

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelHashSet grownUpPatches(grownUpIDs);
    label nNonGrownUp = 0;
    {
        labelList nExternalEdge(mesh.nEdges(), 0);

        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(pp.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            const labelList& eFaces = edgeFaces[edgeI];
            nExternalEdge[meshEdgeI] = eFaces.size();
        }

        syncTools::syncEdgeList
        (
            mesh,
            nExternalEdge,
            plusEqOp<label>(),
            label(0)              // null value
        );

        vectorField grownUpNormals(mesh.nPoints(), vector::zero);
        labelList grownUpEdgePatchID(mesh.nEdges(), -1);

        forAll(mesh.edges(), edgeI)
        {
            if (nExternalEdge[edgeI] == 1)
            {
                isExternalEdge[edgeI] = true;
                const labelList& meshEdgeFaces = mesh.edgeFaces()[edgeI];

                forAll(meshEdgeFaces, k)
                {
                    label faceI = meshEdgeFaces[k];
                    if (!mesh.isInternalFace(faceI))
                    {
                        label patchi = patches.whichPatch(faceI);

                        if (grownUpPatches.found(patchi))
                        {
                            grownUpEdgePatchID[edgeI] = patchi;

                            const face& f =  mesh.faces()[faceI];
                            const edge& e = mesh.edges()[edgeI];
                            vector unitEdgeVec = e.unitVec(mesh.points());

                            label startPt = findIndex(f, e[0]);

                            if (f[f.fcIndex(startPt)] != e[1])
                            {
                                unitEdgeVec = -unitEdgeVec;
                            }

                            point norm  = mesh.faceAreas()[faceI];
                            norm /= (mag(norm) + SMALL);

                            vector gDir = (unitEdgeVec ^ norm);
                            gDir /= (mag(gDir) + SMALL);

                            grownUpNormals[e[0]] += gDir;
                            grownUpNormals[e[1]] += gDir;

                            grownUpPoints[e[0]] = true;
                            grownUpPoints[e[1]] = true;
                            break;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            grownUpPoints,
            orEqOp<bool>(),
            false              // null value
        );

        syncTools::syncPointList
        (
            mesh,
            grownUpNormals,
            plusEqOp<vector>(),
            vector::zero
        );

        syncTools::syncEdgeList
        (
            mesh,
            grownUpEdgePatchID,
            maxEqOp<label>(),
            label(-1)            // initial value
        );

       forAll(pp.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];

            if
            (
                nExternalEdge[meshEdgeI] == 1
                && grownUpEdgePatchID[meshEdgeI] == -1
            )
            {
                edge e = pp.edges()[edgeI];
                if
                (
                    unmarkExtrusion
                    (
                        e[0],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                     )
                )
                {
                    nNonGrownUp++;
                }

                if
                (
                    unmarkExtrusion
                    (
                        e[1],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                     )
                )
                {
                    nNonGrownUp++;
                }
            }
        }

        forAll(pointNormals, pointI)
        {
            label meshPointI = pp.meshPoints()[pointI];

            if (grownUpPoints[meshPointI])
            {
                vector gUpDir = grownUpNormals[meshPointI];
                scalar gUpDirMag = mag(gUpDir);

                if (gUpDirMag > SMALL)
                {
                    gUpDir /=  gUpDirMag;
                }

                if ((gUpDir & pointNormals[pointI]) <= grownUpAngleTerminateCos)
                {
                    if
                    (
                        unmarkExtrusion
                        (
                            pointI,
                            patchDisp,
                            patchNLayers,
                            extrudeStatus
                         )
                     )
                    {
                        nNonGrownUp++;
                    }
                }
            }
        }
    }

    Info<< "Set displacement to zero on "
        << returnReduce(nNonGrownUp, sumOp<label>())
        << " problem grown up points" << endl;
}


//Prevent layers being generated on both sides of single sided patches
void Foam::snappyLayerDriver::handleSingleSidedPatches
(
    const indirectPrimitivePatch& pp,
    const layerParameters& layerParams,

    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const labelList& singleSidedIDs = layerParams.singleSidedIDs();

    if (singleSidedIDs.size())
    {
        const fvMesh& mesh = meshRefiner_.mesh();
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const refinementSurfaces& surfaces = meshRefiner_.surfaces();

        boolList markedBoundary(mesh.boundaryMesh().size(), false);
        forAll(singleSidedIDs, ssID)
        {
            markedBoundary[singleSidedIDs[ssID]] = true;
        }

        labelList checkedFaces(pp.size());
        pointField checkedFaceCentres(pp.size());
        scalarField searchDistance(pp.size());
        label sz = 0;

        const scalar edge0Len =
            meshRefiner_.meshCutter().level0EdgeLength();
        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];
            label patchI = patches.whichPatch(meshFaceI);
            if (markedBoundary[patchI])
            {
                checkedFaces[sz] = i;
                checkedFaceCentres[sz] = mesh.faceCentres()[meshFaceI];
                label own = mesh.faceOwner()[meshFaceI];

                const label cellLevel =
                    meshRefiner_.meshCutter().cellLevel()[own];

                const scalar edgeLen = edge0Len/(1<<cellLevel);
                searchDistance[sz] = sqr(5*edgeLen);
                sz++;
            }
        }
        checkedFaces.setSize(sz);
        checkedFaceCentres.setSize(sz);
        searchDistance.setSize(sz);

        const labelList
            allSurfaces(identity(surfaces.surfaces().size()));

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces.findNearest
        (
            allSurfaces,
            checkedFaceCentres,
            searchDistance, // sqr of attract distance
            hitSurface,
            hitInfo
        );

        label nExcludedFaces = 0;
        labelList pRegions(hitInfo.size());
        vectorField sNormals(hitInfo.size());
        forAll(allSurfaces, sI)
        {
            DynamicList<pointIndexHit> localHits;
            label surfI = allSurfaces[sI];

            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    localHits.append(hitInfo[i]);
                }
            }
            labelList localRegions;
            vectorField localNormals;
            label geomI = surfaces.surfaces()[surfI];
            surfaces.geometry()[geomI].getRegion(localHits, localRegions);
            surfaces.geometry()[geomI].getNormal(localHits, localNormals);

            label localI = 0;
            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    pRegions[i] = localRegions[localI];
                    sNormals[i] = localNormals[localI];
                    localI++;

                    label global = surfaces.globalRegion(surfI, pRegions[i]);

                    if (markedBoundary[globalToSlavePatch_[global]])
                    {
                        label ppI = checkedFaces[i];
                        label meshFaceI = pp.addressing()[ppI];
                        vector fArea = mesh.faceAreas()[meshFaceI];
                        if ((fArea & sNormals[i]) > 0.)
                        {
                            if
                            (
                                unmarkExtrusion
                                (
                                    pp.localFaces()[ppI],
                                    patchDisp,
                                    patchNLayers,
                                    extrudeStatus
                                 )
                             )
                            {
                                nExcludedFaces++;
                            }
                        }
                    }
                }
            }
        }

        Info<< "Set displacement to zero on "
            << returnReduce(nExcludedFaces, sumOp<label>())
            << " faces set as single sided patches " << endl;
    }
}


//Prevent layers being generated in excluded regions
void Foam::snappyLayerDriver::handleExcludedRegions
(
    const indirectPrimitivePatch& pp,
    const layerParameters& layerParams,

    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const PtrList<entry>& excludeRegions = layerParams.excludedRegions();

    if (excludeRegions.size())
    {
        Info<< nl << "Handling Excluded regions ..." << endl;

        const fvMesh& mesh = meshRefiner_.mesh();

        label nExcludedFaces(0);

        boolList excludedCells(mesh.nCells(), false);

        forAll(excludeRegions, regionI)
        {
            const entry& region = excludeRegions[regionI];

            autoPtr<topoSetSource> cellSelector =
                topoSetSource::New(region.keyword(), mesh, region.dict());

            cellSet selectedCellSet
            (
                mesh,
                "cellSet",
                mesh.nCells()/10+1  // Reasonable size estimate.
             );

            cellSelector->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
             );
            forAllConstIter(labelHashSet, selectedCellSet, iter)
            {
                excludedCells[iter.key()] = true;
            }
        }

        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];
            label own = mesh.faceOwner()[meshFaceI];

            if (excludedCells[own])
            {
                if
                (
                    unmarkExtrusion
                    (
                        pp.localFaces()[i],
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                     )
                 )
                {
                    nExcludedFaces++;
                }
            }
        }

        Info<< "Set displacement to zero on "
            << returnReduce(nExcludedFaces, sumOp<label>())
            << " faces inside excluded regions " << endl;
    }
}


// Parallel feature edge detection. Assumes non-manifold edges already handled.
void Foam::snappyLayerDriver::handleFeatureAngle
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const scalar minCos,
    const bool growConvexEdge,
    const bool growConcaveEdge,

    pointField& patchDisp,
    vectorField& patchEdgeNormals,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus,
    PackedList<1>& isConvexEdgePoint,
    PackedList<1>& isConcaveEdgePoint
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling feature edges ..." << endl;

    // Normal component of normals of connected faces.
    vectorField edgeNormal(mesh.nEdges(), vector(GREAT, GREAT, GREAT));

    const labelListList& edgeFaces = pp.edgeFaces();

    forAll(edgeFaces, edgei)
    {
         const labelList& eFaces = pp.edgeFaces()[edgei];

         label meshEdgei = meshEdges[edgei];

         forAll(eFaces, i)
         {
             nomalsCombine()
             (
                 edgeNormal[meshEdgei],
                 pp.faceNormals()[eFaces[i]]
             );
         }
    }

    syncTools::syncEdgeList
    (
        mesh,
        edgeNormal,
        nomalsCombine(),
        vector(GREAT, GREAT, GREAT)          // null value
     );

    autoPtr<OBJstream> str;

    forAll(edgeFaces, edgeI)
    {
        label meshEdgei = meshEdges[edgeI];
        patchEdgeNormals[edgeI] = edgeNormal[meshEdgei];
    }

    if (debug&meshRefinement::MESH)
    {
        str.reset
        (
            new OBJstream
            (
                mesh.time().path()
              / "featureEdges_"
              + meshRefiner_.timeName()
              + ".obj"
            )
        );
        Info<< "Writing feature edges to " << str().name() << endl;
    }

    label nFeats = 0;

    if (minCos < 1-SMALL)
    {
        // Now on coupled edges the edgeNormal will have been truncated and
        // only be still be the old value where two faces have the same normal
        forAll(edgeFaces, edgei)
        {
            const labelList& eFaces = pp.edgeFaces()[edgei];

            label meshEdgei = meshEdges[edgei];

            const vector n = patchEdgeNormals[edgei];

            if (n != vector(GREAT, GREAT, GREAT))
            {
                scalar cos = n & pp.faceNormals()[eFaces[0]];
                if (cos < minCos)
                {
                    const edge& e = pp.edges()[edgei];

                    const point& edgeCentre =
                        mesh.edges()[meshEdgei].centre(mesh.points());
                    const point& faceCentre =
                        mesh.faceCentres()[pp.addressing()[eFaces[0]]];

                    if (((faceCentre - edgeCentre) & n) > 0.)
                    {
                        if (growConvexEdge)
                        {
                            isConvexEdgePoint.set(pp.meshPoints()[e[0]], 1);
                            isConvexEdgePoint.set(pp.meshPoints()[e[1]], 1);
                        }
                        else
                        {
                            unmarkExtrusion
                            (
                                e[0],
                                patchDisp,
                                patchNLayers,
                                extrudeStatus
                             );

                            unmarkExtrusion
                            (
                                e[1],
                                patchDisp,
                                patchNLayers,
                                extrudeStatus
                             );

                            nFeats++;

                            if (str.valid())
                            {
                                const point& p0 = pp.localPoints()[e[0]];
                                const point& p1 = pp.localPoints()[e[1]];
                                str().write(linePointRef(p0, p1));
                            }
                        }
                    }
                    else
                    {
                        if (growConcaveEdge)
                        {
                            isConcaveEdgePoint.set(pp.meshPoints()[e[0]], 1);
                            isConcaveEdgePoint.set(pp.meshPoints()[e[1]], 1);
                        }
                        else
                        {
                            unmarkExtrusion
                            (
                                e[0],
                                patchDisp,
                                patchNLayers,
                                extrudeStatus
                             );

                            unmarkExtrusion
                            (
                                e[1],
                                patchDisp,
                                patchNLayers,
                                extrudeStatus
                             );

                            nFeats++;
                        }
                    }
                }
            }
        }
        Info<< "Set displacement to zero for points on "
            << returnReduce(nFeats, sumOp<label>())
            << " feature edges" << endl;
    }

}


//stop extrusion at non-manifold edges
void Foam::snappyLayerDriver::handleNonManifoldEdges
(
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Precalculate meshEdge per pp edge
    labelList meshEdges(pp.nEdges());

    forAll(meshEdges, patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];

        label v0 = pp.meshPoints()[e[0]];
        label v1 = pp.meshPoints()[e[1]];
        meshEdges[patchEdgeI] = meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[v0],
            v0,
            v1
        );
    }

    labelList nPPEdgeFaces(mesh.nEdges(), 0);
    forAll(pp.edges(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        nPPEdgeFaces[meshEdgeI] = pp.edgeFaces()[edgeI].size();
    }

    syncTools::syncEdgeList
    (
        mesh,
        nPPEdgeFaces,
        plusEqOp<label>(),
        label(0)               // null value
    );

    label nNonManifoldEdges = 0;

    forAll(pp.edges(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        if (nPPEdgeFaces[meshEdgeI] > 2)
        {
            edge e = pp.edges()[edgeI];

            if
            (
                unmarkExtrusion
                (
                    e[0],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                nNonManifoldEdges++;
            }

            if
            (
                unmarkExtrusion
                (
                    e[1],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                nNonManifoldEdges++;
            }
        }
    }

    Info<< "Set displacement to zero on "
        << returnReduce(nNonManifoldEdges, sumOp<label>())
        << " non manifold edges." << endl;
}


// No extrusion on cells with warped faces. Calculates the thickness of the
// layer and compares it to the space the warped face takes up. Disables
// extrusion if layer thickness is more than faceRatio of the thickness of
// the face.
void Foam::snappyLayerDriver::handleWarpedFaces
(
    const indirectPrimitivePatch& pp,
    const scalar faceRatio,
    const scalar edge0Len,
    const labelList& cellLevel,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    if (faceRatio > SMALL)
    {
        const fvMesh& mesh = meshRefiner_.mesh();

        Info<< nl << "Handling cells with warped patch faces ..." << nl;

        const pointField& points = mesh.points();

        label nWarpedFaces = 0;

        forAll(pp, i)
        {
            const face& f = pp[i];

            if (f.size() > 3)
            {
                label facei = pp.addressing()[i];

                label ownLevel = cellLevel[mesh.faceOwner()[facei]];
                scalar edgeLen = edge0Len/(1<<ownLevel);

                // Normal distance to face centre plane
                const point& fc = mesh.faceCentres()[facei];
                const vector& fn = pp.faceNormals()[i];
                scalarField vProj(f.size());

                forAll(f, fp)
                {
                    vector n = points[f[fp]] - fc;
                    vProj[fp] = (n & fn);
                }

                // Get normal 'span' of face
                scalar minVal = min(vProj);
                scalar maxVal = max(vProj);

                if ((maxVal - minVal) > faceRatio * edgeLen)
                {
                    if
                    (
                        unmarkExtrusion
                        (
                            pp.localFaces()[i],
                            patchDisp,
                            patchNLayers,
                            extrudeStatus
                        )
                    )
                    {
                        nWarpedFaces++;
                    }
                }
            }
        }

        Info<< "Set displacement to zero on "
            << returnReduce(nWarpedFaces, sumOp<label>())
            << " warped faces since layer would be > " << faceRatio
            << " of the size of the bounding box." << endl;
    }
}


//stop extrusion at selected edges
void Foam::snappyLayerDriver::handleMiscEdges
(
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Precalculate meshEdge per pp edge
    labelList meshEdges(pp.nEdges());

    forAll(meshEdges, patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];

        label v0 = pp.meshPoints()[e[0]];
        label v1 = pp.meshPoints()[e[1]];
        meshEdges[patchEdgeI] = meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[v0],
            v0,
            v1
        );
    }

    label nMisc = 0;

    //check for zero edges and disable extrusion
    forAll(meshEdges, patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];
        scalar eLen = e.mag(pp.localPoints());
        if (eLen < SMALL)
        {
            if
            (
                unmarkExtrusion
                (
                    e[0],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                nMisc++;
            }

            if
            (
                unmarkExtrusion
                (
                    e[1],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                nMisc++;
            }
        }
    }

    // Check for single face baffles
    // Normal component of normals of connected faces.
    vectorField edgeDisp(mesh.nEdges(), vector(GREAT, GREAT, GREAT));

    const labelListList& edgeFaces = pp.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
         const labelList& eFaces = pp.edgeFaces()[edgeI];
         label meshEdgeI = meshEdges[edgeI];
         point eC = mesh.edges()[meshEdgeI].centre(mesh.points());

         forAll(eFaces, i)
         {
             label meshFaceI = pp.addressing()[eFaces[i]];
             point fC = mesh.faceCentres()[meshFaceI];

             vector eCfC = fC-eC;
             minusCombine()
             (
                 edgeDisp[meshEdgeI],
                 eCfC
             );
         }
    }

    syncTools::syncEdgeList
    (
        mesh,
        edgeDisp,
        minusCombine(),
        vector(GREAT, GREAT, GREAT)          // null value
     );

    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        face f = mesh.faces()[meshFaceI];

        bool baffle = true;
        forAll(f,fp)
        {
            label v0 = f[fp];
            label next = f.fcIndex(fp);
            label v1 = f[next];

            label meshEdgeI = meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[v0],
                v0,
                v1
             );

            if (edgeDisp[meshEdgeI] != vector(GREAT, GREAT, GREAT))
            {
                if (mag(edgeDisp[meshEdgeI])> SMALL)
                {
                    baffle = false;
                    break;
                }
            }
            else
            {
                baffle = false;
                break;
            }
        }


        if (baffle)
        {
            if
            (
                unmarkExtrusion
                (
                    pp.localFaces()[i],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                 )
             )
            {
                nMisc++;
            }
        }
    }

    Info<< "Set displacement to zero on "
        << returnReduce(nMisc, sumOp<label>())
        << " faces and edges with single baffle extrusions and "
        << "zero sized edges" << endl;

}


//// No extrusion on cells with multiple patch faces. There ususally is a reason
//// why combinePatchFaces hasn't succeeded.
//void Foam::snappyLayerDriver::handleMultiplePatchFaces
//(
//    const indirectPrimitivePatch& pp,
//    pointField& patchDisp,
//    labelList& patchNLayers,
//    List<extrudeMode>& extrudeStatus
//) const
//{
//    const fvMesh& mesh = meshRefiner_.mesh();
//
//    Info<< nl << "Handling cells with multiple patch faces ..." << nl;
//
//    const labelListList& pointFaces = pp.pointFaces();
//
//    // Cells that should not get an extrusion layer
//    cellSet multiPatchCells(mesh, "multiPatchCells", pp.size());
//
//    // Detect points that use multiple faces on same cell.
//    forAll(pointFaces, patchPointi)
//    {
//        const labelList& pFaces = pointFaces[patchPointi];
//
//        labelHashSet pointCells(pFaces.size());
//
//        forAll(pFaces, i)
//        {
//            label celli = mesh.faceOwner()[pp.addressing()[pFaces[i]]];
//
//            if (!pointCells.insert(celli))
//            {
//                // Second or more occurrence of cell so cell has two or more
//                // pp faces connected to this point.
//                multiPatchCells.insert(celli);
//            }
//        }
//    }
//
//    label nMultiPatchCells = returnReduce
//    (
//        multiPatchCells.size(),
//        sumOp<label>()
//    );
//
//    Info<< "Detected " << nMultiPatchCells
//        << " cells with multiple (connected) patch faces." << endl;
//
//    label nChanged = 0;
//
//    if (nMultiPatchCells > 0)
//    {
//        multiPatchCells.instance() = meshRefiner_.timeName();
//        Info<< "Writing " << nMultiPatchCells
//            << " cells with multiple (connected) patch faces to cellSet "
//            << multiPatchCells.objectPath() << endl;
//        multiPatchCells.write();
//
//
//        // Go through all points and remove extrusion on any cell in
//        // multiPatchCells
//        // (has to be done in separate loop since having one point on
//        // multipatches has to reset extrusion on all points of cell)
//
//        forAll(pointFaces, patchPointi)
//        {
//            if (extrudeStatus[patchPointi] != NOEXTRUDE)
//            {
//                const labelList& pFaces = pointFaces[patchPointi];
//
//                forAll(pFaces, i)
//                {
//                    label celli =
//                        mesh.faceOwner()[pp.addressing()[pFaces[i]]];
//
//                    if (multiPatchCells.found(celli))
//                    {
//                        if
//                        (
//                            unmarkExtrusion
//                            (
//                                patchPointi,
//                                patchDisp,
//                                patchNLayers,
//                                extrudeStatus
//                            )
//                        )
//                        {
//                            nChanged++;
//                        }
//                    }
//                }
//            }
//        }
//
//        reduce(nChanged, sumOp<label>());
//    }
//
//    Info<< "Prevented extrusion on " << nChanged
//        << " points due to multiple patch faces." << nl << endl;
//}


void Foam::snappyLayerDriver::setNumLayers
(
    const labelList& patchToNLayers,
    const labelList& patchIDs,
    const indirectPrimitivePatch& pp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus,
    label& nAddedCells
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl << "Handling points with inconsistent layer specification ..."
        << endl;

    // Get for every point (really only necessary on patch external points)
    // the max and min of any patch faces using it.
    labelList maxLayers(patchNLayers.size(), labelMin);
    labelList minLayers(patchNLayers.size(), labelMax);

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];

        const labelList& meshPoints = mesh.boundaryMesh()[patchi].meshPoints();

        label wantedLayers = patchToNLayers[patchi];

        forAll(meshPoints, patchPointi)
        {
            label ppPointi = pp.meshPointMap()[meshPoints[patchPointi]];

            maxLayers[ppPointi] = max(wantedLayers, maxLayers[ppPointi]);
            minLayers[ppPointi] = min(wantedLayers, minLayers[ppPointi]);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxLayers,
        maxEqOp<label>(),
        labelMin            // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minLayers,
        minEqOp<label>(),
        labelMax            // null value
    );

    // Unmark any point with different min and max
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //label nConflicts = 0;

    forAll(maxLayers, i)
    {
        if (maxLayers[i] == labelMin || minLayers[i] == labelMax)
        {
            FatalErrorInFunction
                << "Patchpoint:" << i << " coord:" << pp.localPoints()[i]
                << " maxLayers:" << maxLayers[i]
                << " minLayers:" << minLayers[i]
                << abort(FatalError);
        }
        else
        {
            // Ok setting.
            patchNLayers[i] = maxLayers[i];
        }
    }


    // Calculate number of cells to create
    nAddedCells = 0;
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        // Get max of extrusion per point
        label nCells = 0;
        forAll(f, fp)
        {
            nCells = max(nCells, patchNLayers[f[fp]]);
        }

        nAddedCells += nCells;
    }
    reduce(nAddedCells, sumOp<label>());

    //reduce(nConflicts, sumOp<label>());
    //
    //Info<< "Set displacement to zero for " << nConflicts
    //    << " points due to points being on multiple regions"
    //    << " with inconsistent nLayers specification." << endl;
}


// Construct pointVectorField with correct boundary conditions for adding
// layers
Foam::tmp<Foam::pointVectorField>
Foam::snappyLayerDriver::makeLayerDisplacementField
(
    const pointMesh& pMesh,
    const layerParameters& layerParams
)
{
    // Construct displacement field.
    const labelList& numLayers = layerParams.numLayers();
    const labelList& grownUpIDs = layerParams.grownUpIDs();
    const bool useFixedNormalSlip = layerParams.fixedNormalSlip();

    const pointBoundaryMesh& pointPatches = pMesh.boundary();
    const polyMesh& mesh = pMesh();

    labelHashSet grownUpPatches(grownUpIDs);

    wordList patchFieldTypes
    (
        pointPatches.size(),
        slipPointPatchField<vector>::typeName
    );

    if (useFixedNormalSlip)
    {
        forAll(patchFieldTypes, patchi)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (!pp.coupled())
            {
                vector pNorm0(vector::zero);
                if (pp.size() > 0)
                {
                    vector fA = pp.faceAreas()[0];
                    scalar area = mag(fA);
                    if (area > SMALL)
                    {
                        fA /= area;
                        pNorm0 = fA;
                    }
                }
                reduce(pNorm0,maxMagSqrOp<point>());

                bool valid = true;
                forAll(pp, i)
                {
                    vector fA = pp.faceAreas()[i];
                    scalar area = mag(fA);
                    if (area > SMALL)
                    {
                        fA /= area;
                        if ((pNorm0&fA)<0.984)
                        {
                            valid = false;
                            break;
                        }
                    }
                }
                reduce(valid,andOp<bool>());

                if (valid)
                {
                    patchFieldTypes[patchi] =
                        fixedNormalSlipPointPatchField<vector>::typeName;
                }
            }
        }
    }

    wordList actualPatchTypes(patchFieldTypes.size());
    forAll(pointPatches, patchi)
    {
        actualPatchTypes[patchi] = pointPatches[patchi].type();
    }

    forAll(numLayers, patchi)
    {
        if (!grownUpPatches.found(patchi))
        {
            //  0 layers: do not allow slip so fixedValue 0
            // >0 layers: fixedValue which gets adapted
            if (numLayers[patchi] == 0)
            {
                patchFieldTypes[patchi] =
                    zeroFixedValuePointPatchVectorField::typeName;
            }
            else if (numLayers[patchi] > 0)
            {
                patchFieldTypes[patchi] =
                    fixedValuePointPatchVectorField::typeName;
            }
        }
    }

    forAll(pointPatches, patchi)
    {
        if (isA<processorPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = calculatedPointPatchVectorField::typeName;
        }
        else if (isA<cyclicPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = cyclicSlipPointPatchVectorField::typeName;
        }
    }


    // Note: time().timeName() instead of meshRefinement::timeName() since
    // postprocessable field.
    tmp<pointVectorField> tfld
    (
        new pointVectorField
        (
            IOobject
            (
                "pointDisplacement",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("displacement", dimLength, Zero),
            patchFieldTypes,
            actualPatchTypes
        )
    );

    if (useFixedNormalSlip)
    {
        const pointVectorField::Boundary& dispBf =
            tfld().boundaryField();
        forAll(dispBf, patchi)
        {
            if
            (
                isA<fixedNormalSlipPointPatchField<vector>>(dispBf[patchi])
            )
            {
                fixedNormalSlipPointPatchField<vector>& patchFNS =
                    const_cast<fixedNormalSlipPointPatchField<vector>&>
                    (
                        refCast<const fixedNormalSlipPointPatchField<vector>>
                        (
                            dispBf[patchi]
                         )
                     );


                const polyPatch& pp = mesh.boundaryMesh()[patchi];
                point pNorm(vector::zero);
                scalar pAreas(scalar(0));

                forAll(pp, facei)
                {
                    pNorm += pp.faceAreas()[facei];
                    pAreas += pp.magFaceAreas()[facei];
                }
                reduce(
                    std::tie(pNorm, pAreas),
                    ParallelOp<sumOp<point>, sumOp<scalar>>{}
                );

                if (pAreas > SMALL)
                {
                    pNorm /= pAreas;
                    pNorm /= (mag(pNorm) + SMALL);
                    patchFNS.setNormal(pNorm);
                }
            }
        }
    }

    return tfld;
}


void Foam::snappyLayerDriver::growNoExtrusion
(
    const indirectPrimitivePatch& pp,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    Info<< nl << "Growing non-extrusion points by one layer ..." << endl;

    List<extrudeMode> grownExtrudeStatus(extrudeStatus);

    const faceList& localFaces = pp.localFaces();

    label nGrown = 0;

    forAll(localFaces, facei)
    {
        const face& f = localFaces[facei];

        bool hasSqueeze = false;
        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] == NOEXTRUDE)
            {
                hasSqueeze = true;
                break;
            }
        }

        if (hasSqueeze)
        {
            // Squeeze all points of face
            forAll(f, fp)
            {
                if
                (
                    extrudeStatus[f[fp]] == EXTRUDE
                 && grownExtrudeStatus[f[fp]] != NOEXTRUDE
                )
                {
                    grownExtrudeStatus[f[fp]] = NOEXTRUDE;
                    nGrown++;
                }
            }
        }
    }

    extrudeStatus.transfer(grownExtrudeStatus);


    // Synchronise since might get called multiple times.
    // Use the fact that NOEXTRUDE is the minimum value.
    {
        labelList status(extrudeStatus.size());
        forAll(status, i)
        {
            status[i] = extrudeStatus[i];
        }
        syncTools::syncPointList
        (
            meshRefiner_.mesh(),
            pp.meshPoints(),
            status,
            minEqOp<label>(),
            labelMax            // null value
        );
        forAll(status, i)
        {
            extrudeStatus[i] = extrudeMode(status[i]);
        }
    }


    forAll(extrudeStatus, patchPointi)
    {
        if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            patchDisp[patchPointi] = Zero;
            patchNLayers[patchPointi] = 0;
        }
    }

    reduce(nGrown, sumOp<label>());

    Info<< "Set displacement to zero for an additional " << nGrown
        << " points." << endl;
}


void Foam::snappyLayerDriver::determineSidePatches
(
    fvMesh& mesh,

    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,
    const labelList grownUpPatchIDs,
    labelList& edgePatchID,
    labelList& edgeZoneID,
    boolList& edgeFlip,
    labelList& inflateFaceID
)
{
    // Sometimes edges-to-be-extruded are on more than 2 processors.
    // Work out which 2 hold the faces to be extruded and thus which procpatch
    // the edge-face should be in. As an additional complication this might
    // mean that 2 procesors that were only edge-connected now suddenly need
    // to become face-connected i.e. have a processor patch between them.

    // Determine edgePatchID. Any additional processor boundary gets added to
    // patchToNbrProc,nbrProcToPatch and nPatches gets set to the new number
    // of patches.
    label nPatches;
    Map<label> nbrProcToPatch;
    Map<label> patchToNbrProc;
    addPatchCellLayer::calcExtrudeInfo
    (
        true,           // zoneFromAnyFace

        mesh,
        globalFaces,
        edgeGlobalFaces,
        pp,
        grownUpPatchIDs,
        edgePatchID,
        nPatches,
        nbrProcToPatch,
        patchToNbrProc,
        edgeZoneID,
        edgeFlip,
        inflateFaceID
    );

    label nOldPatches = mesh.boundaryMesh().size();
    label nAdded = returnReduce(nPatches-nOldPatches, sumOp<label>());
    Info<< nl << "Adding in total " << nAdded/2 << " inter-processor patches to"
        << " handle extrusion of non-manifold processor boundaries."
        << endl;

    if (nAdded > 0)
    {
        // We might not add patches in same order as in patchToNbrProc
        // so prepare to renumber edgePatchID
        Map<label> wantedToAddedPatch;

        for (label patchi = nOldPatches; patchi < nPatches; patchi++)
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

            //Pout<< "Adding patch " << patchi
            //    << " name:" << name
            //    << " between " << Pstream::myProcNo()
            //    << " and " << nbrProci << endl;

            label procPatchi = meshRefiner_.appendPatch
            (
                mesh,
                name,
                patchDict
            );
            wantedToAddedPatch.insert(patchi, procPatchi);
        }

        // Renumber edgePatchID
        forAll(edgePatchID, i)
        {
            label patchi = edgePatchID[i];
            Map<label>::const_iterator fnd = wantedToAddedPatch.find(patchi);
            if (fnd != wantedToAddedPatch.end())
            {
                edgePatchID[i] = fnd();
            }
        }

        mesh.clearOut();
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh()).updateMesh();
    }
}


void Foam::snappyLayerDriver::calculateLayerThickness
(
    const indirectPrimitivePatch& pp,
    const labelList& patchIDs,
    const layerParameters& layerParams,
    const labelList& cellLevel,
    const scalar edge0Len,

    labelList& patchNLayers,
    scalarField& thickness,
    scalarField& minThickness,
    scalarField& targetExpansion,
    scalarField& targetfch,
    boolList& fixedfch
) const
{
    const List<wordList>& layerSpec = layerParams.layerSpec();
    const labelList& patchToNLayers = layerParams.numLayers();
    const labelList& patchToMinLayers = layerParams.minLayers();
    const scalarField& patchFCH = layerParams.fch();
    const scalarField& patchMaxLayerThickness = layerParams.maxLayerThickness();
    const scalarField& patchExpansionRatio = layerParams.expansionRatio();
    const scalarField& patchFinalLayerThickness =
        layerParams.finalLayerThickness();
    const scalarField& patchMinThickness = layerParams.minThickness();

    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Rework patch-wise layer parameters into minimum per point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Rework relative thickness into absolute
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // by multiplying with the internal cell size.

    // Reuse input fields
    fixedfch.setSize(pp.nPoints());
    fixedfch = false;
    targetfch.setSize(pp.nPoints());
    targetfch = GREAT;
    targetExpansion.setSize(pp.nPoints());
    targetExpansion = GREAT;
    thickness.setSize(pp.nPoints());
    thickness = GREAT;
    minThickness.setSize(pp.nPoints());
    minThickness = GREAT;
    patchNLayers.setSize(pp.nPoints());
    patchNLayers = labelMax;

    labelList maxPointLevel(pp.nPoints(), labelMin);
    labelList minPointLevel(pp.nPoints(), labelMax);

    forAll(pp, i)
    {
        label ownLevel = cellLevel[mesh.faceOwner()[pp.addressing()[i]]];

        const face& f = pp.localFaces()[i];

        forAll(f, fp)
        {
            maxPointLevel[f[fp]] = max(maxPointLevel[f[fp]], ownLevel);
            minPointLevel[f[fp]] = min(minPointLevel[f[fp]], ownLevel);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxPointLevel,
        maxEqOp<label>(),
        labelMin           // null value
     );

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minPointLevel,
        minEqOp<label>(),
        labelMax           // null value
     );

    forAll(patchIDs, i)
    {
        label patchi = patchIDs[i];
        const word& patchName = patches[patchi].name();

        HashSet<word> order(layerSpec[patchi]);

        const bool patchFixFCH = layerParams.fixedFCH()[patchi];

        const labelList& meshPoints = patches[patchi].meshPoints();

        forAll(meshPoints, pointi)
        {
            label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
            scalar edgeLen = 0.5*
                (edge0Len/(1<<maxPointLevel[ppPointi])
                 +edge0Len/(1<<minPointLevel[ppPointi]));

            minThickness[ppPointi] = min
            (
                minThickness[ppPointi],
                edgeLen * patchMinThickness[patchi]
             );

            fixedfch[ppPointi] = patchFixFCH;
        }

        if (order.found("numLayers") && patchToNLayers[patchi])
        {
            if (order.found("expansionRatio"))
            {
                if (order.found("finalLayerThickness"))
                {
                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));

                        label wantedLayers = patchToNLayers[patchi];
                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            edgeLen * patchFinalLayerThickness[patchi]
                        );
                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            patchExpansionRatio[patchi]
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                        );

                        scalar pfch;
                        if (patchNLayers[ppPointi] <= 1)
                        {
                            pfch = thickness[ppPointi];
                        }
                        else
                        {
                            scalar strfac = pow
                            (
                                targetExpansion[ppPointi],patchNLayers[ppPointi]-1
                            );
                            pfch = thickness[ppPointi] / strfac;
                        }

                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            pfch
                        );
                    }
                }
                else if (order.found("fch") || order.found("rfch"))
                {
                    label wantedLayers = patchToNLayers[patchi];
                    scalar flh = patchFCH[patchi]
                        * pow(patchExpansionRatio[patchi], wantedLayers -1);

                    bool relative(order.found("rfch") ? true : false);
                    scalar relativeSclaing = 1;

                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        if (relative)
                        {
                            relativeSclaing = 0.5*
                                (edge0Len/(1<<maxPointLevel[ppPointi])
                                 +edge0Len/(1<<minPointLevel[ppPointi]));
                        }

                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            relativeSclaing*flh
                        );
                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            patchExpansionRatio[patchi]
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                        );

                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            patchFCH[patchi] * relativeSclaing
                        );
                    }
                }
                else if (order.found("maxLayerThickness"))
                {
                    label wantedLayers = patchToNLayers[patchi];
                    scalar flt;
                    if (patchExpansionRatio[patchi] == 1)
                    {
                        flt = patchMaxLayerThickness[patchi] / wantedLayers;
                    }
                    else
                    {
                        flt = patchMaxLayerThickness[patchi]
                            * (1 - patchExpansionRatio[patchi])
                            / (
                                  1 - pow
                                  (
                                      patchExpansionRatio[patchi],
                                      wantedLayers
                                  )
                              );
                        flt *= pow(patchExpansionRatio[patchi], wantedLayers-1);
                    }

                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));

                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            flt*edgeLen
                        );

                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            patchExpansionRatio[patchi]
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                        );

                        scalar pfch;
                        if (patchNLayers[ppPointi] <= 1)
                        {
                            pfch = thickness[ppPointi];
                        }
                        else
                        {
                            scalar strfac = pow
                            (
                                targetExpansion[ppPointi],
                                patchNLayers[ppPointi]-1
                            );
                            pfch = thickness[ppPointi] / strfac;
                        }

                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            pfch
                        );
                    }
                }
            }
            else if (order.found("finalLayerThickness"))
            {
                if (order.found("fch") || order.found("rfch"))
                {
                    label wantedLayers = patchToNLayers[patchi];

                    bool relative(order.found("rfch") ? true : false);
                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));

                        scalar ratio = 1.0;
                        if (relative)
                        {
                            ratio = patchFinalLayerThickness[patchi]
                                / patchFCH[patchi];
                        }
                        else
                        {
                            ratio = patchFinalLayerThickness[patchi] * edgeLen
                                / (patchFCH[patchi]);
                        }

                        scalar er = 1.0;
                        if (wantedLayers > 1)
                        {
                            er = pow(ratio, 1./(wantedLayers-1));
                        }

                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            er
                        );

                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                        );
                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            edgeLen * patchFinalLayerThickness[patchi]
                        );

                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            (
                                relative ?
                                patchFCH[patchi] * edgeLen
                                : patchFCH[patchi]
                             )
                        );
                    }
                }
                else if (order.found("maxLayerThickness"))
                {
                    scalar xn = 0.0;
                    scalar x = -1.0;
                    label wantedLayers = patchToNLayers[patchi];

                    scalar a1 = patchFinalLayerThickness[patchi];
                    scalar a2 = patchMaxLayerThickness[patchi];
                    label repeat = 0;

                    scalar tol = 0.001;

                    while (!(repeat++ >= 100 || mag(x - xn) < tol))
                    {
                        x = xn;
                        scalar f1 = a1 * pow(x, wantedLayers) + a2*(1 - x) - a1;
                        scalar f2 = a1 * wantedLayers
                                  * pow(x, wantedLayers -1) - a2;
                        xn = x - f1/(f2+SMALL);
                    }

                    if (xn < 1.0 + tol && xn > 1.0 - tol)
                    {
                        xn = 10.0;
                        x = -1.0;
                        repeat = 0;
                        while (!(repeat++ >= 100 || mag(x - xn) < tol))
                        {
                            x = xn;
                            scalar f1 = a1 * pow(x, wantedLayers)
                                      + a2*(1 - x) - a1;
                            scalar f2 = a1 * wantedLayers
                                * pow(x, wantedLayers -1) - a2;
                            xn = x - f1/(f2+SMALL);
                        }
                    }

                    //now calculate inverse
                    scalar er = 1. / (xn+SMALL);

                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));

                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            er
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                         );
                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            edgeLen * patchFinalLayerThickness[patchi]
                        );

                        scalar pfch;
                        if (patchNLayers[ppPointi] <= 1)
                        {
                            pfch = thickness[ppPointi];
                        }
                        else
                        {
                            scalar strfac = pow
                            (
                                targetExpansion[ppPointi],
                                patchNLayers[ppPointi]-1
                            );
                            pfch = thickness[ppPointi] / strfac;
                        }

                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            pfch
                        );
                    }
                }
            }
            else if (order.found("fch") || order.found("rfch"))
            {
                bool relative(order.found("rfch") ? true : false);
                scalar relativeScaling = 1.0;
                if (order.found("maxLayerThickness"))
                {
                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));

                        scalar xn = 10.0;
                        scalar x = -1.0;
                        label wantedLayers = patchToNLayers[patchi];

                        if (relative)
                        {
                            relativeScaling = edgeLen;
                        }

                        scalar a1 = patchFCH[patchi] * relativeScaling;
                        scalar a2 = patchMaxLayerThickness[patchi] * edgeLen;
                        label repeat = 0;

                        scalar tol = 0.001;

                        while (!(repeat++ >= 100 || mag(x - xn) < tol))
                        {
                            x = xn;
                            scalar f1 = a1 * pow(x, wantedLayers)
                                + a2*(1 - x) - a1;
                            scalar f2 = a1 * wantedLayers
                                * pow(x, wantedLayers -1) - a2;
                            xn = x - f1/(f2+SMALL);
                        }
                        if (xn < 1.0 + tol && xn > 1.0 - tol)
                        {
                            xn = 0.0;
                            x = -1.0;
                            repeat = 0;
                            while (!(repeat++ >= 100 || mag(x - xn) < tol))
                            {
                                x = xn;
                                scalar f1 = a1 * pow(x, wantedLayers)
                                          + a2*(1 - x) - a1;
                                scalar f2 = a1 * wantedLayers
                                          * pow(x, wantedLayers -1) - a2;
                                xn = x - f1/(f2+SMALL);
                            }
                        }

                        scalar er = xn;
                        scalar flt = a1 * pow(xn, wantedLayers - 1);

                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            er
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                         );
                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            flt
                        );
                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            patchFCH[patchi] * relativeScaling
                        );
                    }
                }
            }
        }
        else if (order.found("expansionRatio"))
        {
            if (order.found("finalLayerThickness"))
            {
                if (order.found("fch"))
                {
                    bool relative(order.found("rfch") ? true : false);
                    scalar relativeScaling = 1.0;
                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));

                        if (relative)
                        {
                            relativeScaling = edgeLen;
                        }

                        scalar flt = patchFinalLayerThickness[patchi] * edgeLen;
                        scalar ratio = flt
                            / (patchFCH[patchi] * relativeScaling);
                        label wantedLayers =
                            (log(ratio) / log(patchExpansionRatio[patchi])) + 1;
                        label minLayers = patchToMinLayers[patchi];
                        wantedLayers = max(minLayers,wantedLayers);

                        if (wantedLayers <= 0)
                        {
                            wantedLayers = 1;
                        }

                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            patchExpansionRatio[patchi]
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                        );
                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            flt
                        );
                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            patchFCH[patchi] * relativeScaling
                        );
                    }
                }
                else if (order.found("maxLayerThickness"))
                {
                    scalar invStr = 1. / patchExpansionRatio[patchi];
                    scalar flt = patchFinalLayerThickness[patchi];
                    scalar mlt = patchMaxLayerThickness[patchi];
                    scalar ratio = 1.  - (mlt * (1. - invStr) / flt);

                    label wantedLayers = 0;
                    if (ratio < SMALL)
                    {
                        WarningInFunction
                            << "No solution for growing patch layers with "
                            << "current settings for patch: "
                            << patchName
                            << " Switching off layer growth for this patch"
                            << endl;
                    }
                    else
                    {
                        wantedLayers = log(ratio)
                            / log(invStr);

                        label minLayers = patchToMinLayers[patchi];
                        wantedLayers = max(minLayers,wantedLayers);

                        if (wantedLayers <= 0)
                        {
                            wantedLayers = 1;
                        }
                    }

                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));

                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            patchExpansionRatio[patchi]
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                        );
                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            flt*edgeLen
                        );

                        scalar pfch;
                        if (patchNLayers[ppPointi] <= 1)
                        {
                            pfch = thickness[ppPointi];
                        }
                        else
                        {
                            scalar strfac = pow
                            (
                                targetExpansion[ppPointi],
                                patchNLayers[ppPointi]-1
                            );
                            pfch = thickness[ppPointi] / strfac;
                        }

                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            pfch
                        );
                    }
                }
            }
            if (order.found("fch") || order.found("rfch"))
            {
                bool relative(order.found("rfch") ? true : false);
                scalar relativeScaling = 1.0;
                if (order.found("maxLayerThickness"))
                {
                    scalar str = patchExpansionRatio[patchi];
                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));
                        if (relative)
                        {
                            relativeScaling = edgeLen;
                        }

                        scalar fch = patchFCH[patchi] * relativeScaling;
                        scalar mlt = patchMaxLayerThickness[patchi] * edgeLen;
                        scalar ratio = 1. - (mlt * (1. - str) / fch);

                        label wantedLayers = 0;
                        if (ratio < SMALL)
                        {
                            WarningInFunction
                                << "No solution for growing patch layers with "
                                << "current settings for patch: "
                                << patchName
                                << " Switching off layer growth for this patch"
                                << endl;
                        }
                        else
                        {
                            wantedLayers = log(ratio) / log(str);

                            label minLayers = patchToMinLayers[patchi];
                            wantedLayers = max(minLayers,wantedLayers);

                            if (wantedLayers <= 0)
                            {
                                wantedLayers = 1;
                            }
                        }
                        scalar flt = fch * pow(str, wantedLayers -1);

                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            patchExpansionRatio[patchi]
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                        );
                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            flt
                        );
                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            (patchFCH[patchi] * relativeScaling)
                        );
                    }
                }
            }
        }
        else if (order.found("finalLayerThickness"))
        {
            if (order.found("fch") || order.found("rfch"))
            {
                bool relative(order.found("rfch") ? true : false);
                scalar relativeScaling = 1.0;
                if (order.found("maxLayerThickness"))
                {
                    forAll(meshPoints, pointi)
                    {
                        label ppPointi = pp.meshPointMap()[meshPoints[pointi]];
                        scalar edgeLen = 0.5*
                            (edge0Len/(1<<maxPointLevel[ppPointi])
                             +edge0Len/(1<<minPointLevel[ppPointi]));

                        if (relative)
                        {
                            relativeScaling = edgeLen;
                        }
                        scalar flt = patchFinalLayerThickness[patchi] * edgeLen;
                        scalar mlt = patchMaxLayerThickness[patchi] * edgeLen;
                        scalar fch = patchFCH[patchi] * relativeScaling;
                        scalar str = (fch - mlt)/(flt - mlt);

                        label wantedLayers = 0;
                        if (str < SMALL)
                        {
                            WarningInFunction
                                << "No solution for growing patch layers with "
                                << "current settings for patch: "
                                << patchName
                                << " Switching off layer growth for this patch"
                                << endl;
                        }
                        else
                        {
                            scalar ratio = flt * str / fch;
                            wantedLayers = log(ratio) / log(str);
                            label minLayers = patchToMinLayers[patchi];
                            wantedLayers = max(minLayers,wantedLayers);
                        }

                        targetExpansion[ppPointi] = min
                        (
                            targetExpansion[ppPointi],
                            str
                        );
                        patchNLayers[ppPointi] = min
                        (
                            patchNLayers[ppPointi],
                            wantedLayers
                        );
                        thickness[ppPointi] = min
                        (
                            thickness[ppPointi],
                            flt
                        );
                        targetfch[ppPointi] = min
                        (
                            targetfch[ppPointi],
                            patchFCH[patchi] * relativeScaling
                        );
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        fixedfch,
        orEqOp<bool>(),
        false           // null value
     );

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        targetExpansion,
        minEqOp<scalar>(),
        GREAT               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        thickness,
        minEqOp<scalar>(),
        GREAT               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minThickness,
        minEqOp<scalar>(),
        GREAT               // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchNLayers,
        minEqOp<label>(),
        labelMax            // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        targetfch,
        minEqOp<scalar>(),
        GREAT               // null value
    );

    // Rework thickness (of final layer) into overall thickness of all layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(thickness, pointi)
    {
        // Calculate layer thickness based on expansion ratio
        // and final layer height
        if (targetExpansion[pointi] == 1)
        {
            thickness[pointi] *= patchNLayers[pointi];
        }
        else
        {
            scalar invExpansion = 1.0 / targetExpansion[pointi];
            label nLay = patchNLayers[pointi];
            thickness[pointi] *=
                (1.0 - pow(invExpansion, nLay))
              / (1.0 - invExpansion);
        }
    }

    //Info<< "calculateLayerThickness : min:" << gMin(thickness)
    //    << " max:" << gMax(thickness) << endl;
}


// Synchronize displacement among coupled patches.
void Foam::snappyLayerDriver::syncPatchDisplacement
(
    const indirectPrimitivePatch& pp,
    const scalarField& minThickness,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& meshPoints = pp.meshPoints();

//    label nChangedTotal = 0;

    while (true)
    {
        label nChanged = 0;

        // Sync displacement (by taking min)
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            patchDisp,
            minMagSqrEqOp<vector>(),
            point::rootMax      // null value
        );

        // Unmark if displacement too small
        forAll(patchDisp, i)
        {
            if (mag(patchDisp[i]) < minThickness[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }

        labelList syncPatchNLayers(patchNLayers);

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            syncPatchNLayers,
            minEqOp<label>(),
            labelMax            // null value
        );

        // Reset if differs
        // 1. take max
        forAll(syncPatchNLayers, i)
        {
            if (syncPatchNLayers[i] != patchNLayers[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            syncPatchNLayers,
            maxEqOp<label>(),
            labelMin            // null value
        );

        // Reset if differs
        // 2. take min
        forAll(syncPatchNLayers, i)
        {
            if (syncPatchNLayers[i] != patchNLayers[i])
            {
                if
                (
                    unmarkExtrusion
                    (
                        i,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    )
                )
                {
                    nChanged++;
                }
            }
        }
//        nChangedTotal += nChanged;

        if (!returnReduce(nChanged, sumOp<label>()))
        {
            break;
        }
    }

    //Info<< "Prevented extrusion on "
    //    << returnReduce(nChangedTotal, sumOp<label>())
    //    << " coupled patch points during syncPatchDisplacement." << endl;
}


// Calculate displacement vector for all patch points. Uses pointNormal.
// Checks that displaced patch point would be visible from all centres
// of the faces using it.
// extrudeStatus is both input and output and gives the status of each
// patch point.
void Foam::snappyLayerDriver::getPatchDisplacement
(
    const indirectPrimitivePatch& pp,
    const scalarField& thickness,
    const scalarField& minThickness,
    const labelList& grownUpIDs,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    Info<< nl << "Determining displacement for added points"
        << " according to pointNormal ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelListList& edgeFaces = pp.edgeFaces();
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();
    const pointField& localPoints = pp.localPoints();

    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~

    const labelList& meshPoints = pp.meshPoints();
    pointField pointNormals(pp.nPoints(), vector::zero);
    boolList isExternalPt(pp.meshPoints().size(),false);

    {
        labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));
        labelList nExternalEdge(mesh.nEdges(), 0);
        forAll(pp.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            const labelList& eFaces = edgeFaces[edgeI];
            nExternalEdge[meshEdgeI] = eFaces.size();
        }

        syncTools::syncEdgeList
        (
            mesh,
            nExternalEdge,
            plusEqOp<label>(),
            label(0)              // null value
        );

        forAll(pp.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            if (nExternalEdge[meshEdgeI] == 1)
            {
                edge e = pp.edges()[edgeI];
                isExternalPt[e[0]] = true;
                isExternalPt[e[1]] = true;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            isExternalPt,
            orEqOp<bool>(),
            false               // null value
        );

        pointField pointNormalsMesh(mesh.nPoints(), vector::zero);
        scalarField nPointFacesMesh(mesh.nPoints(), 0.);

        forAll(faceNormals, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            forAll(f, fp)
            {
                label pointI = f[fp];
                if (!isExternalPt[pointI])
                {
                    label meshPointI = meshPoints[pointI];

                    scalar fA = mesh.magFaceAreas()[pp.addressing()[faceI]];
                    pointNormalsMesh[meshPointI] += faceNormals[faceI]*fA;
                    nPointFacesMesh[meshPointI] += fA;
                }
            }
        }

        labelHashSet grownPatches(grownUpIDs);

        forAll(mesh.edges(), meshEdgeI)
        {
            if (nExternalEdge[meshEdgeI] == 1)
            {
                const labelList& meshEdgeFaces = mesh.edgeFaces()[meshEdgeI];

                forAll(meshEdgeFaces, k)
                {
                    label faceI = meshEdgeFaces[k];
                    if (!mesh.isInternalFace(faceI))
                    {
                        label patchi = patches.whichPatch(faceI);

                        if (grownPatches.found(patchi))
                        {
                            edge e = mesh.edges()[meshEdgeI];
                            scalar eLength = e.mag(mesh.points());

                            vector eVec = e.vec(mesh.points());
                            eVec /= (eLength + SMALL);
                            point norm  = -mesh.faceAreas()[faceI];
                            norm /= (mag(norm) + SMALL);

                            face f = mesh.faces()[faceI];

                            label i1 = findIndex(f, e[0]);

                            vector gDir(norm^eVec);
                            gDir /= (mag(gDir) + SMALL);

                            if (f[f.fcIndex(i1)] != e[1])
                            {
                                gDir = -gDir;
                            }

                            pointNormalsMesh[e[0]] += gDir*eLength;
                            nPointFacesMesh[e[0]] += eLength;
                            pointNormalsMesh[e[1]] += gDir*eLength;
                            nPointFacesMesh[e[1]] += eLength;
                        }
                    }
                }
            }
        }


        syncTools::syncPointList
        (
            mesh,
            pointNormalsMesh,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            nPointFacesMesh,
            plusEqOp<scalar>(),
            scalar(0.)          // null value
        );

        forAll(pointNormals, i)
        {
            label meshPointI = meshPoints[i];
            if (nPointFacesMesh[meshPointI] > SMALL)
            {
                pointNormals[i] = pointNormalsMesh[meshPointI]
                    /nPointFacesMesh[meshPointI];
            }
        }
    }
    pointNormals /= (mag(pointNormals) + SMALL);

    // Determine local length scale on patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Start off from same thickness everywhere (except where no extrusion)
    patchDisp = thickness*pointNormals;

    label nNoVisNormal = 0;
    label nExtrudeRemove = 0;


    // Check if no extrude possible.
    forAll(pointNormals, patchPointi)
    {
        label meshPointi = pp.meshPoints()[patchPointi];

        if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            // Do not use unmarkExtrusion; forcibly set to zero extrusion.
            patchNLayers[patchPointi] = 0;
            patchDisp[patchPointi] = Zero;
        }
        else if (!isExternalPt[patchPointi])
        {
            // Get normal
            const vector& n = pointNormals[patchPointi];

            if (!meshTools::visNormal(n, faceNormals, pointFaces[patchPointi]))
            {
                if (debug&meshRefinement::ATTRACTION)
                {
                    Pout<< "No valid normal for point " << meshPointi
                        << ' ' << pp.points()[meshPointi]
                        << "; setting displacement to "
                        << patchDisp[patchPointi]
                        << endl;
                }

                extrudeStatus[patchPointi] = EXTRUDEREMOVE;
                nNoVisNormal++;
            }
        }
    }

    // At illegal points make displacement average of new neighbour positions
    forAll(extrudeStatus, patchPointi)
    {
        if (extrudeStatus[patchPointi] == EXTRUDEREMOVE)
        {
            point avg(Zero);
            label nPoints = 0;

            const labelList& pEdges = pp.pointEdges()[patchPointi];

            forAll(pEdges, i)
            {
                label edgei = pEdges[i];

                label otherPointi = pp.edges()[edgei].otherVertex(patchPointi);

                if (extrudeStatus[otherPointi] != NOEXTRUDE)
                {
                    avg += localPoints[otherPointi] + patchDisp[otherPointi];
                    nPoints++;
                }
            }

            if (nPoints > 0)
            {
                if (debug&meshRefinement::ATTRACTION)
                {
                    Pout<< "Displacement at illegal point "
                        << localPoints[patchPointi]
                        << " set to "
                        << (avg / nPoints - localPoints[patchPointi])
                        << endl;
                }

                patchDisp[patchPointi] =
                    avg / nPoints
                  - localPoints[patchPointi];

                nExtrudeRemove++;
            }
            else
            {
                // All surrounding points are not extruded. Leave patchDisp
                // intact.
            }
        }
    }

    Info<< "Detected " << returnReduce(nNoVisNormal, sumOp<label>())
        << " points with point normal pointing through faces." << nl
        << "Reset displacement at "
        << returnReduce(nExtrudeRemove, sumOp<label>())
        << " points to average of surrounding points." << endl;

    // Make sure displacement is equal on both sides of coupled patches.
    syncPatchDisplacement
    (
        pp,
        minThickness,
        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    Info<< endl;
}


bool Foam::snappyLayerDriver::sameEdgeNeighbour
(
    const labelListList& globalEdgeFaces,
    const label myGlobalFacei,
    const label nbrGlobFacei,
    const label edgei
) const
{
    const labelList& eFaces = globalEdgeFaces[edgei];
    if (eFaces.size() == 2)
    {
        return edge(myGlobalFacei, nbrGlobFacei) == edge(eFaces[0], eFaces[1]);
    }
    else
    {
        return false;
    }
}


void Foam::snappyLayerDriver::getVertexString
(
    const indirectPrimitivePatch& pp,
    const labelListList& globalEdgeFaces,
    const label facei,
    const label edgei,
    const label myGlobFacei,
    const label nbrGlobFacei,
    DynamicList<label>& vertices
) const
{
    const labelList& fEdges = pp.faceEdges()[facei];
    label fp = findIndex(fEdges, edgei);

    if (fp == -1)
    {
        FatalErrorInFunction
            << "problem." << abort(FatalError);
    }

    // Search back
    label startFp = fp;

    forAll(fEdges, i)
    {
        label prevFp = fEdges.rcIndex(startFp);
        if
        (
           !sameEdgeNeighbour
            (
                globalEdgeFaces,
                myGlobFacei,
                nbrGlobFacei,
                fEdges[prevFp]
            )
        )
        {
            break;
        }
        startFp = prevFp;
    }

    label endFp = fp;
    forAll(fEdges, i)
    {
        label nextFp = fEdges.fcIndex(endFp);
        if
        (
           !sameEdgeNeighbour
            (
                globalEdgeFaces,
                myGlobFacei,
                nbrGlobFacei,
                fEdges[nextFp]
            )
        )
        {
            break;
        }
        endFp = nextFp;
    }

    const face& f = pp.localFaces()[facei];
    vertices.clear();
    fp = startFp;
    while (fp != endFp)
    {
        vertices.append(f[fp]);
        fp = f.fcIndex(fp);
    }
    vertices.append(f[fp]);
    fp = f.fcIndex(fp);
    vertices.append(f[fp]);
}


// Truncates displacement
// - for all patchFaces in the faceset displacement gets set to zero
// - all displacement < minThickness gets set to zero
Foam::label Foam::snappyLayerDriver::truncateDisplacement
(
    const globalIndex& globalFaces,
    const labelListList& edgeGlobalFaces,
    const indirectPrimitivePatch& pp,
    const scalarField& minThickness,
    const faceSet& illegalPatchFaces,
    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    label nChanged = 0;

    const Map<label>& meshPointMap = pp.meshPointMap();

    forAllConstIter(faceSet, illegalPatchFaces, iter)
    {
        label facei = iter.key();

        if (mesh.isInternalFace(facei))
        {
            FatalErrorInFunction
                << "Faceset " << illegalPatchFaces.name()
                << " contains internal face " << facei << nl
                << "It should only contain patch faces" << abort(FatalError);
        }

        const face& f = mesh.faces()[facei];


        forAll(f, fp)
        {
            if (meshPointMap.found(f[fp]))
            {
                label patchPointi = meshPointMap[f[fp]];

                if (extrudeStatus[patchPointi] != NOEXTRUDE)
                {
                    unmarkExtrusion
                    (
                        patchPointi,
                        patchDisp,
                        patchNLayers,
                        extrudeStatus
                    );
                    nChanged++;
                }
            }
        }
    }

    forAll(patchDisp, patchPointi)
    {
        if (mag(patchDisp[patchPointi]) < minThickness[patchPointi])
        {
            if
            (
                unmarkExtrusion
                (
                    patchPointi,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                nChanged++;
            }
        }
        else if (extrudeStatus[patchPointi] == NOEXTRUDE)
        {
            // Make sure displacement is 0. Should already be so but ...
            patchDisp[patchPointi] = Zero;
            patchNLayers[patchPointi] = 0;
        }
    }

    nChanged += stopNonConsecutiveExtrusion
    (
        globalFaces,
        edgeGlobalFaces,
        pp,
        minThickness,
        patchDisp,
        patchNLayers,
        extrudeStatus
    );

    return returnReduce(nChanged, sumOp<label>());
}


// Setup layer information (at points and faces) to modify mesh topology in
// regions where layer mesh terminates.
void Foam::snappyLayerDriver::setupLayerInfoTruncation
(
    const indirectPrimitivePatch& pp,
    const labelList& patchNLayers,
    const List<extrudeMode>& extrudeStatus,
    const label& nBufferCellsNoExtrude,
    const label& terminationStrategy,
    const label& layerRecovery,
    const PackedList<1>& isConvexEdgePoint,
    const PackedList<1>& isConcaveEdgePoint,
    labelList& nPatchPointLayers,
    labelList& nPatchFaceLayers
) const
{
    Info<< nl << "Setting up information for layer truncation ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (nBufferCellsNoExtrude < 0)
    {
        Info<< nl << "Performing no layer truncation."
            << " nBufferCellsNoExtrude set to less than 0  ..." << endl;

        // Face layers if any point gets extruded
        forAll(pp.localFaces(), patchFacei)
        {
            const face& f = pp.localFaces()[patchFacei];

            forAll(f, fp)
            {
                if (patchNLayers[f[fp]] > 0)
                {
                    nPatchFaceLayers[patchFacei] = patchNLayers[f[fp]];
                    break;
                }
            }
        }
        nPatchPointLayers = patchNLayers;

        // Set any unset patch face layers
        forAll(nPatchFaceLayers, patchFacei)
        {
            if (nPatchFaceLayers[patchFacei] == -1)
            {
                nPatchFaceLayers[patchFacei] = 0;
            }
        }
    }
    else
    {
        // Determine max point layers per face.
        labelList maxLevel(pp.size(), 0);

        forAll(pp.localFaces(), patchFacei)
        {
            const face& f = pp.localFaces()[patchFacei];

            // find patch faces where layer terminates (i.e contains extrude
            // and noextrude points).

            bool noExtrude = false;
            label mLevel = 0;

            forAll(f, fp)
            {
                if
                (
                    extrudeStatus[f[fp]] == NOEXTRUDE
                    //  || isConcaveEdgePoint.get(pp.meshPoints()[f[fp]]) == 1
                    //  || isConvexEdgePoint.get(pp.meshPoints()[f[fp]]) == 1
                )
                {
                    noExtrude = true;
                }
                mLevel = max(mLevel, patchNLayers[f[fp]]);
            }

            if (mLevel > 0)
            {
                // So one of the points is extruded. Check if all are extruded
                // or is a mix.

                if (noExtrude)
                {
                    nPatchFaceLayers[patchFacei] = 1;
                    maxLevel[patchFacei] = mLevel;
                }
                else
                {
                    maxLevel[patchFacei] = mLevel;
                }
            }
        }

        // We have the seed faces (faces with nPatchFaceLayers != maxLevel)
        // Now do a meshwave across the patch where we pick up neighbours
        // of seed faces.
        // Note: quite inefficient. Could probably be coded better.

        const labelListList& pointFaces = pp.pointFaces();

        label nLevels = gMax(patchNLayers);

        // flag neighbouring patch faces with number of layers to grow
        for (label ilevel = 1; ilevel < nLevels; ilevel+=layerRecovery)
        {
            label nBuffer;

            if (ilevel == 1  && terminationStrategy == 0)
            {
                nBuffer = nBufferCellsNoExtrude - 1;
            }
            else
            {
                nBuffer = nBufferCellsNoExtrude;
            }

            for (label ibuffer = 0; ibuffer < nBuffer + 1; ibuffer++)
            {
                labelList tempCounter(nPatchFaceLayers);

                boolList foundNeighbour(pp.nPoints(), false);

                forAll(pp.meshPoints(), patchPointi)
                {
                    forAll(pointFaces[patchPointi], pointFacei)
                    {
                        label facei = pointFaces[patchPointi][pointFacei];

                        if
                        (
                            nPatchFaceLayers[facei] != -1
                         && maxLevel[facei] > 0
                        )
                        {
                            foundNeighbour[patchPointi] = true;
                            break;
                        }
                    }
                }

                syncTools::syncPointList
                (
                    mesh,
                    pp.meshPoints(),
                    foundNeighbour,
                    orEqOp<bool>(),
                    false               // null value
                );

                forAll(pp.meshPoints(), patchPointi)
                {
                    if (foundNeighbour[patchPointi])
                    {
                        forAll(pointFaces[patchPointi], pointFacei)
                        {
                            label facei = pointFaces[patchPointi][pointFacei];
                            if
                            (
                                nPatchFaceLayers[facei] == -1
                             && maxLevel[facei] > 0
                             && ilevel < maxLevel[facei]
                            )
                            {
                                tempCounter[facei] = ilevel;
                            }
                        }
                    }
                }
                nPatchFaceLayers = tempCounter;
            }
        }

        boolList marked(pp.size(), false);
        forAll(pp.localFaces(), patchFacei)
        {
            if (nPatchFaceLayers[patchFacei] == -1)
            {
                nPatchFaceLayers[patchFacei] = maxLevel[patchFacei];
                marked[patchFacei] = true;
            }
        }

        forAll(pp.meshPoints(), patchPointi)
        {
            if (extrudeStatus[patchPointi] != NOEXTRUDE)
            {
                forAll(pointFaces[patchPointi], pointFacei)
                {
                    label face = pointFaces[patchPointi][pointFacei];
                    nPatchPointLayers[patchPointi] = max
                    (
                        nPatchPointLayers[patchPointi],
                        nPatchFaceLayers[face]
                    );
                }
            }
            else
            {
                nPatchPointLayers[patchPointi] = 0;
            }
        }
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPatchPointLayers,
            maxEqOp<label>(),
            label(0)        // null value
         );

        while (true)
        {
            label nChanged = 0;

            labelList minPatchPointLayers(pp.meshPoints().size(), labelMax);

            forAll(pp.meshPoints(), patchPointI)
            {
                if (extrudeStatus[patchPointI] != NOEXTRUDE)
                {
                    forAll(pointFaces[patchPointI], pointFaceI)
                    {
                        label face = pointFaces[patchPointI][pointFaceI];
                        minPatchPointLayers[patchPointI] = min
                        (
                            minPatchPointLayers[patchPointI],
                            nPatchFaceLayers[face]
                         );
                    }
                }
                else
                {
                    minPatchPointLayers[patchPointI] = 0;
                }
            }
            syncTools::syncPointList
            (
                mesh,
                pp.meshPoints(),
                minPatchPointLayers,
                minEqOp<label>(),
                labelMax        // null value
             );


            forAll(pp.localFaces(), patchFaceI)
            {
                if (marked[patchFaceI])
                {
                    const face& f = pp.localFaces()[patchFaceI];
                    label minLayers = labelMax;
                    label maxLayers = labelMin;
                    forAll(f, fp)
                    {
                        minLayers = min(minPatchPointLayers[f[fp]], minLayers);
                        maxLayers = max(minPatchPointLayers[f[fp]], maxLayers);
                    }

                    if (maxLayers - minLayers > 1)
                    {
                        if
                        (
                            nPatchFaceLayers[patchFaceI] != maxLayers -1
                            && nPatchFaceLayers[patchFaceI] >  minLayers
                         )
                        {
                            nPatchFaceLayers[patchFaceI] = maxLayers -1;
                            nChanged++;
                        }
                    }
                }
            }

            if (!returnReduce(nChanged, sumOp<label>()))
            {
                break;
            }

            forAll(pp.meshPoints(), patchPointI)
            {
                if (extrudeStatus[patchPointI] != NOEXTRUDE)
                {
                    nPatchPointLayers[patchPointI] = labelMin;
                    forAll(pointFaces[patchPointI], pointFaceI)
                    {
                        label face = pointFaces[patchPointI][pointFaceI];
                        nPatchPointLayers[patchPointI] = max
                        (
                            nPatchPointLayers[patchPointI],
                            nPatchFaceLayers[face]
                        );
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                pp.meshPoints(),
                nPatchPointLayers,
                maxEqOp<label>(),
                label(0)                  // null value
            );
        }
    }
}


// Does any of the cells use a face from faces?
bool Foam::snappyLayerDriver::cellsUseFace
(
    const polyMesh& mesh,
    const labelList& cellLabels,
    const labelList& faceLabels,
    const labelHashSet& faces
)
{
    forAll(cellLabels, i)
    {
        const cell& cFaces = mesh.cells()[cellLabels[i]];

        forAll(cFaces, cFacei)
        {
            if (faces.found(cFaces[cFacei]))
            {
                return true;
            }
        }
    }

    // Look at first cell above layer
    if (faceLabels.size())
    {
        const label& own = mesh.faceOwner()[faceLabels[0]];
        const cell& cFaces = mesh. cells()[own];
        forAll(cFaces, cFaceI)
        {
            if (faces.found(cFaces[cFaceI]))
            {
                return true;
            }
        }
    }

    return false;
}


// Checks the newly added cells and locally unmarks points so they
// will not get extruded next time round. Returns global number of unmarked
// points (0 if all was fine)
Foam::label Foam::snappyLayerDriver::checkAndUnmark
(
    const addPatchCellLayer& addLayer,
    const dictionary& meshQualityDict,
    const bool additionalReporting,
    const List<labelPair>& baffles,
    const indirectPrimitivePatch& pp,
    const fvMesh& newMesh,
    const labelList& faceMap,

    pointField& patchDisp,
    labelList& patchNLayers,
    List<extrudeMode>& extrudeStatus,

    boolList& nonConvergedAfterNLayerIter,
    const label nLayerIter
)
{
    // Check the resulting mesh for errors
    Info<< nl << "Checking mesh with layer ..." << endl;
    faceSet wrongFaces(newMesh, "wrongFaces", newMesh.nFaces()/1000);
    label illegalFaceCount = motionSmoother::checkMesh
    (
        false,
        newMesh,
        meshQualityDict,
        identity(newMesh.nFaces()),
        baffles,
        wrongFaces
    );
    Info<< "Detected " << illegalFaceCount
        << " illegal faces"
        << " (concave, zero area or negative cell pyramid volume)"
        << endl;

    // Undo local extrusion if
    // - any of the added cells in error

    label nChanged = 0;
    nonConvergedAfterNLayerIter= false;

    // Check if any of the faces in error uses any face of an added cell
    // - if additionalReporting print the few remaining areas for ease of
    //   finding out where the problems are.

    const label nReportMax = 10;
    DynamicField<point> disabledFaceCentres(nReportMax);

    forAll(addLayer.layerCells(), oldPatchFacei)
    {
        // Get the cells (in newMesh labels) per old patch face (in mesh
        // labels)
        const labelList& fCells = addLayer.layerCells()[oldPatchFacei];
        const labelList& fFaces = addLayer.layerFaces()[oldPatchFacei];

        if (cellsUseFace(newMesh, fCells, fFaces, wrongFaces))
        {
            label origPatchFacei = faceMap[oldPatchFacei];
            // Unmark points on old mesh
            if
            (
                unmarkExtrusion
                (
                    pp.localFaces()[origPatchFacei],
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                )
            )
            {
                if (additionalReporting && (nChanged < nReportMax))
                {
                    disabledFaceCentres.append
                    (
                        pp.faceCentres()[oldPatchFacei]
                    );
                }

                nChanged++;
                if (nLayerIter > 10)
                {
                    const face f = pp.localFaces()[origPatchFacei];
                    forAll(f, fp)
                    {
                        nonConvergedAfterNLayerIter[f[fp]] = true;
                    }
                }
            }
        }
    }

    label nChangedTotal = returnReduce(nChanged, sumOp<label>());

    if (additionalReporting)
    {
        // Limit the number of points to be printed so that
        // not too many points are reported when running in parallel
        // Not accurate, i.e. not always nReportMax points are written,
        // but this estimation avoid some communication here.
        // The important thing, however, is that when only a few faces
        // are disabled, their coordinates are printed, and this should be
        // the case
        label nReportLocal = nChanged;
        if (nChangedTotal > nReportMax)
        {
            nReportLocal = min
            (
                max(nChangedTotal / Pstream::nProcs(), 1),
                min
                (
                    nChanged,
                    max(nReportMax / Pstream::nProcs(), 1)
                )
            );
        }

        if (nReportLocal)
        {
            Pout<< "Checked mesh with layers. Disabled extrusion at " << endl;
            for (label i=0; i < nReportLocal; i++)
            {
                Pout<< "    " << disabledFaceCentres[i] << endl;
            }
        }

        label nReportTotal = returnReduce(nReportLocal, sumOp<label>());

        if (nReportTotal < nChangedTotal)
        {
            Info<< "Suppressed disabled extrusion message for other "
                << nChangedTotal - nReportTotal << " faces." << endl;
        }
    }

    return nChangedTotal;
}


//- Count global number of extruded faces
Foam::label Foam::snappyLayerDriver::countExtrusion
(
    const indirectPrimitivePatch& pp,
    const List<extrudeMode>& extrudeStatus
)
{
    // Count number of extruded patch faces
    label nExtruded = 0;
    {
        const faceList& localFaces = pp.localFaces();

        forAll(localFaces, i)
        {
            const face& localFace = localFaces[i];

            forAll(localFace, fp)
            {
                if (extrudeStatus[localFace[fp]] != NOEXTRUDE)
                {
                    nExtruded++;
                    break;
                }
            }
        }
    }

    return returnReduce(nExtruded, sumOp<label>());
}


Foam::List<Foam::labelPair> Foam::snappyLayerDriver::getBafflesOnAddedMesh
(
    const polyMesh& mesh,
    const labelList& newToOldFaces,
    const List<labelPair>& baffles
)
{
    // The problem is that the baffle faces are now inside the
    // mesh (addPatchCellLayer modifies original boundary faces and
    // adds new ones. So 2 pass:
    // - find the boundary face for all faces originating from baffle
    // - use the boundary face for the new baffles

    Map<label> baffleSet(4*baffles.size());
    forAll(baffles, bafflei)
    {
        baffleSet.insert(baffles[bafflei][0], bafflei);
        baffleSet.insert(baffles[bafflei][1], bafflei);
    }


    List<labelPair> newBaffles(baffles.size(), labelPair(-1, -1));
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        label oldFacei = newToOldFaces[facei];

        Map<label>::const_iterator faceFnd = baffleSet.find(oldFacei);
        if (faceFnd != baffleSet.end())
        {
            label bafflei = faceFnd();
            labelPair& p = newBaffles[bafflei];
            if (p[0] == -1)
            {
                p[0] = facei;
            }
            else if (p[1] == -1)
            {
                p[1] = facei;
            }
            else
            {
                FatalErrorInFunction
                    << "Problem:" << facei << " at:"
                    << mesh.faceCentres()[facei]
                    << " is on same baffle as " << p[0]
                    << " at:" << mesh.faceCentres()[p[0]]
                    << " and " << p[1]
                    << " at:" << mesh.faceCentres()[p[1]]
                    << exit(FatalError);
            }
        }
    }
    return newBaffles;
}


// Collect layer faces and layer cells into mesh fields for ease of handling
void Foam::snappyLayerDriver::getLayerCellsFaces
(
    const polyMesh& mesh,
    const addPatchCellLayer& addLayer,
    const scalarField& oldRealThickness,

    labelList& cellNLayers,
    scalarField& faceRealThickness
)
{
    cellNLayers.setSize(mesh.nCells());
    cellNLayers = 0;
    faceRealThickness.setSize(mesh.nFaces());
    faceRealThickness = 0;

    // Mark all faces in the layer
    const labelListList& layerFaces = addLayer.layerFaces();

    // Mark all cells in the layer.
    labelListList addedCells(addPatchCellLayer::addedCells(mesh, layerFaces));

    forAll(addedCells, oldPatchFacei)
    {
        const labelList& added = addedCells[oldPatchFacei];

        const labelList& layer = layerFaces[oldPatchFacei];

        if (layer.size())
        {
            // Leave out original internal face
            forAll(added, i)
            {
                cellNLayers[added[i]] = layer.size()-1;
            }
        }
    }

    forAll(layerFaces, oldPatchFacei)
    {
        const labelList& layer = layerFaces[oldPatchFacei];
        const scalar realThickness = oldRealThickness[oldPatchFacei];

        if (layer.size())
        {
            // Layer contains both original boundary face and new boundary
            // face so is nLayers+1. Leave out old internal face.
            for (label i = 1; i < layer.size(); i++)
            {
                faceRealThickness[layer[i]] = realThickness;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::snappyLayerDriver::snappyLayerDriver
(
    meshRefinement& meshRefiner,
    autoPtr<searchableSurfaces>& allGeometryPtr,
    hexReport& stats,
    const meshControl& controller,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
)
:
    meshRefiner_(meshRefiner),
    allGeometryPtr_(allGeometryPtr),
    stats_(stats),
    controller_(controller),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::snappyLayerDriver::addLayers
(
    layerParameters& layerParams,
    const dictionary& motionDict,
    const labelList& patchIDs,
    const label multiLayerIter
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // faceZones of type internal or baffle (for merging points across)
    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = meshRefiner_.getZones(fzTypes);
    }

    // faceZones of type internal (for checking mesh quality across and
    // merging baffles)
    const labelList internalFaceZones
    (
        meshRefiner_.getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::INTERNAL
            )
        )
    );

    // Create baffles (pairs of faces that share the same points)
    // Baffles stored as owner and neighbour face that have been created.
    List<labelPair> baffles;
    {
        labelList originatingFaceZone;
        meshRefiner_.createZoneBaffles
        (
            identity(mesh.faceZones().size()),
            baffles,
            originatingFaceZone,
            false
        );

        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing baffled mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner_.timeName(),
                false
            );
        }
    }

    // Check initial mesh
    Info<< "Checking initial mesh ..." << endl;
    labelHashSet wrongFaces(mesh.nFaces()/100);
    label nAllowableErrors = motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);

    Info<< "Detected " << nAllowableErrors << " illegal faces"
        << " (concave, zero area or negative cell pyramid volume)"
        << endl;

    // Duplicate points on faceZones of type boundary. Should normally already
    // be done by snapping phase
    {
        autoPtr<mapPolyMesh> map = meshRefiner_.dupNonManifoldBoundaryPoints();
        if (map.valid())
        {
            const labelList& reverseFaceMap = map().reverseFaceMap();
            forAll(baffles, i)
            {
                label f0 = reverseFaceMap[baffles[i].first()];
                label f1 = reverseFaceMap[baffles[i].second()];
                baffles[i] = labelPair(f0, f1);
            }
        }
    }

    //setup which patches to try growing up
    if (layerParams.growUpPatches())
    {
        layerParams.setGrownUpIDs(mesh.boundaryMesh());
    }

    const labelList& grownUpIDs = layerParams.grownUpIDs();
    labelHashSet grownPatches(grownUpIDs);

    if (layerParams.projectGrownUp() > SMALL)
    {
        //Calculate grown up planes for subsequent projections onto these planes
        List<plane> grownUpPlanes = calculateGrownUpPlanes(grownUpIDs);

        projectToGrownUpPatches(layerParams, grownUpIDs, grownUpPlanes);
    }


    // Now we have all patches determine the number of layer per patch
    // This will be layerParams.numLayers() except for faceZones where one
    // side does get layers and the other not in which case we want to
    // suppress movement by explicitly setting numLayers 0
    labelList numLayers(layerParams.numLayers());
    {
        labelHashSet layerIDs(patchIDs);
        forAll(mesh.faceZones(), zonei)
        {
            label mpi, spi;
            surfaceZonesInfo::faceZoneType fzType;
            bool hasInfo = meshRefiner_.getFaceZoneInfo
            (
                mesh.faceZones()[zonei].name(),
                mpi,
                spi,
                fzType
            );
            if (hasInfo)
            {
                const polyBoundaryMesh& pbm = mesh.boundaryMesh();
                if (layerIDs.found(mpi) && !layerIDs.found(spi))
                {
                    // Only layers on master side. Fix points on slave side
                    Info<< "On faceZone " << mesh.faceZones()[zonei].name()
                        << " adding layers to master patch " << pbm[mpi].name()
                        << " only. Freezing points on slave patch "
                        << pbm[spi].name() << endl;
                    numLayers[spi] = 0;
                }
                else if (!layerIDs.found(mpi) && layerIDs.found(spi))
                {
                    // Only layers on slave side. Fix points on master side
                    Info<< "On faceZone " << mesh.faceZones()[zonei].name()
                        << " adding layers to slave patch " << pbm[spi].name()
                        << " only. Freezing points on master patch "
                        << pbm[mpi].name() << endl;
                    numLayers[mpi] = 0;
                }
            }
        }
    }



    // Duplicate points on faceZones that layers are added to
    labelList pointToMaster;

    {
        labelList candidatePoints;

        // Do full analysis to see if we need to extrude points
        // so have to duplicate them
        autoPtr<indirectPrimitivePatch> pp
        (
            meshRefinement::makePatch
            (
                mesh,
                patchIDs
            )
        );

        List<extrudeMode> extrudeStatus(pp().nPoints(), EXTRUDE);

        // Do duplication only if all patch points decide to extrude. Ignore
        // contribution from non-patch points. Note that we need to
        // apply this to all mesh points
        labelList minPatchState(mesh.nPoints(), labelMax);
        forAll(extrudeStatus, patchPointi)
        {
            label pointi = pp().meshPoints()[patchPointi];
            minPatchState[pointi] = extrudeStatus[patchPointi];
        }

        syncTools::syncPointList
        (
            mesh,
            minPatchState,
            minEqOp<label>(),   // combine op
            labelMax            // null value
        );

        // So now minPatchState:
        // - labelMax on non-patch points
        // - NOEXTRUDE if any patch point was not extruded
        // - EXTRUDE or EXTRUDEREMOVE if all patch points are extruded/
        //   extrudeRemove.

        label n = 0;
        forAll(minPatchState, pointi)
        {
            label state = minPatchState[pointi];
            if (state == EXTRUDE || state == EXTRUDEREMOVE)
            {
                n++;
            }
        }
        candidatePoints.setSize(n);
        n = 0;
        forAll(minPatchState, pointi)
        {
            label state = minPatchState[pointi];
            if (state == EXTRUDE || state == EXTRUDEREMOVE)
            {
                candidatePoints[n++] = pointi;
            }
        }

        // Not duplicate points on either side of baffles that don't get any
        // layers
        labelPairList nonDupBaffles;

        {
            // faceZones that are not being duplicated
            DynamicList<label> nonDupZones(mesh.faceZones().size());

            labelHashSet layerIDs(patchIDs);
            forAll(mesh.faceZones(), zonei)
            {
                label mpi, spi;
                surfaceZonesInfo::faceZoneType fzType;
                bool hasInfo = meshRefiner_.getFaceZoneInfo
                (
                    mesh.faceZones()[zonei].name(),
                    mpi,
                    spi,
                    fzType
                );
                if (hasInfo && !layerIDs.found(mpi) && !layerIDs.found(spi))
                {
                    nonDupZones.append(zonei);
                }
            }
            nonDupBaffles = meshRefinement::subsetBaffles
            (
                mesh,
                nonDupZones,
                localPointRegion::findDuplicateFacePairs(mesh)
            );
        }

        const localPointRegion regionSide(mesh, nonDupBaffles, candidatePoints);

        autoPtr<mapPolyMesh> map = meshRefiner_.dupNonManifoldPoints
        (
            regionSide
        );

        if (map.valid())
        {
            // Store point duplication
            pointToMaster.setSize(mesh.nPoints(), -1);

            const labelList& pointMap = map().pointMap();
            const labelList& reversePointMap = map().reversePointMap();

            forAll(pointMap, pointi)
            {
                label oldPointi = pointMap[pointi];
                label newMasterPointi = reversePointMap[oldPointi];

                if (newMasterPointi != pointi)
                {
                    // Found slave. Mark both master and slave
                    pointToMaster[pointi] = newMasterPointi;
                    pointToMaster[newMasterPointi] = newMasterPointi;
                }
            }

            // Update baffle numbering
            {
                const labelList& reverseFaceMap = map().reverseFaceMap();
                forAll(baffles, i)
                {
                    label f0 = reverseFaceMap[baffles[i].first()];
                    label f1 = reverseFaceMap[baffles[i].second()];
                    baffles[i] = labelPair(f0, f1);
                }
            }


            if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing point-duplicate mesh to time "
                    << meshRefiner_.timeName() << endl;

                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    mesh.time().path()/meshRefiner_.timeName(),
                    false
                );

                OBJstream str
                (
                    mesh.time().path()
                  / "duplicatePoints_"
                  + meshRefiner_.timeName()
                  + ".obj"
                );
                Info<< "Writing point-duplicates to " << str.name() << endl;
                const pointField& p = mesh.points();
                forAll(pointMap, pointi)
                {
                    label newMasteri = reversePointMap[pointMap[pointi]];

                    if (newMasteri != pointi)
                    {
                        str.write(linePointRef(p[pointi], p[newMasteri]));
                    }
                }
            }
        }
    }


    // Add layers to patches
    // ~~~~~~~~~~~~~~~~~~~~~

    // Now we have
    // - mesh with optional baffles and duplicated points for faceZones
    //   where layers are to be added
    // - pointToMaster    : correspondence for duplicated points
    // - baffles          : list of pairs of faces

    autoPtr<indirectPrimitivePatch> pp
    (
        meshRefinement::makePatch
        (
            mesh,
            patchIDs
        )
    );

    //If protected cells field exist make sure all layer patch cells set
    if (mesh.foundObject<volScalarField>("protectedCells"))
    {
        volScalarField& protectedCells =
           const_cast<volScalarField&>
           (
                mesh.lookupObject<volScalarField>("protectedCells")
           );

        const labelList& meshPoints = pp().meshPoints();
        forAll(meshPoints, ptI)
        {
            label meshPointI = meshPoints[ptI];
            const labelList& pCells = mesh.pointCells()[meshPointI];
            forAll(pCells, pCellI)
            {
                label cellI = pCells[pCellI];
                if (protectedCells[cellI] == -1)
                {
                    protectedCells[cellI] = scalar(1);
                }
            }
        }
    }

    // Global face indices engine
    const globalIndex globalFaces(mesh.nFaces());

    // Determine extrudePatch.edgeFaces in global numbering (so across
    // coupled patches). This is used only to string up edges between coupled
    // faces (all edges between same (global)face indices get extruded).
    labelListList edgeGlobalFaces
    (
        addPatchCellLayer::globalEdgeFaces
        (
            mesh,
            globalFaces,
            pp
        )
    );

    // Point-wise extrusion data
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    // Displacement for all pp.localPoints. Set to large value so does
    // not trigger the minThickness truncation (see syncPatchDisplacement below)
    vectorField patchDisp(pp().nPoints(), vector(GREAT, GREAT, GREAT));

    // Number of layers for all pp.localPoints. Note: only valid if
    // extrudeStatus = EXTRUDE.
    labelList patchNLayers(pp().nPoints(), 0);

    // Whether to add edge for all pp.localPoints.
    List<extrudeMode> extrudeStatus(pp().nPoints(), EXTRUDE);

    // Disable extrusion on feature angles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Points on a convex edge.
    PackedList<1> isConvexEdgePoint(mesh.nPoints(), 0);
    // Points on a concave edge.
    PackedList<1> isConcaveEdgePoint(mesh.nPoints(), 0);
    // Edge normal for all pp. edges.
    vectorField patchEdgeNormals(pp().nEdges(), vector::zero);

    // Undistorted edge length
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    // Determine (wanted) point-wise layer thickness and expansion ratio
    scalarField thickness(pp().nPoints());
    scalarIOField minThickness
    (
        IOobject
        (
            "minThickness",
            meshRefiner_.timeName(),
            mesh,
            IOobject::NO_READ
        ),
        pp().nPoints()
    );
    scalarField targetExpansion(pp().nPoints());
    scalarField targetfch(pp().nPoints());
    boolList fixedfch(pp().nPoints());

    calculateLayerThickness
    (
        pp,
        layerPatchIDs(),
        layerParams,
        cellLevel,
        edge0Len,

        patchNLayers,
        thickness,
        minThickness,
        targetExpansion,
        targetfch,
        fixedfch
    );

    {
        // Precalculate mesh edge labels for patch edges
        labelList meshEdges(pp().meshEdges(mesh.edges(), mesh.pointEdges()));

        //Collapse any single sided patch faces
        if (!allGeometryPtr_.empty())
        {
            handleSingleSidedPatches
            (
                pp,
                layerParams,
                patchDisp,
                patchNLayers,
                extrudeStatus
            );

            //Now clear up surface geometry  as no longer needed
            allGeometryPtr_.clear();
        }

        // Disable extrusion on split strings of common points
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        handleNonStringConnected
        (
            pp,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        handleExcludedRegions
        (
            pp,
            layerParams,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        handleFeatureAngle
        (
            pp,
            meshEdges,
            layerParams.cosAngleTermination(),
            layerParams.growConvexEdge(),
            layerParams.growConcaveEdge(),

            patchDisp,
            patchEdgeNormals,
            patchNLayers,
            extrudeStatus,
            isConvexEdgePoint,
            isConcaveEdgePoint
        );

        // Disable extrusion on warped faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        handleWarpedFaces
        (
            pp,
            layerParams.maxFaceThicknessRatio(),
            edge0Len,
            cellLevel,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        #if !defined( WIN32 ) && !defined( WIN64 )
        addProfiling(grow, "snappyHexMesh::layers::grow");
        #endif

        // Disable extrusion on non manifold edges
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        handleNonManifoldEdges
        (
            pp,
            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        // Disable extrusion on non manifold edges
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        handleGrownUpEdges
        (
            pp,
            meshEdges,
            layerParams.grownUpIDs(), // patches to grow layers up
            layerParams.grownUpAngleTerminateCos(),

            patchDisp,
            patchNLayers,
            extrudeStatus
        );

        handleMiscEdges
        (
            pp,

            patchDisp,
            patchNLayers,
            extrudeStatus
        );
    }

    // Print a bit
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Find maximum length of a patch name, for a nicer output
        label maxPatchNameLen = 0;
        forAll(layerPatchIDs(), i)
        {
            label patchi = layerPatchIDs()[i];
            word patchName = patches[patchi].name();
            maxPatchNameLen = max(maxPatchNameLen, label(patchName.size()));
        }

        Info<< nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "patch"
            << setw(0) << " faces    layers avg thickness[m]" << nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << " "
            << setw(0) << "                 near-wall overall" << nl
            << setf(ios_base::left) << setw(maxPatchNameLen) << "-----"
            << setw(0) << " -----    ------ --------- -------" << endl;

        forAll(layerPatchIDs(), i)
        {
            label patchi = layerPatchIDs()[i];
            const labelList& meshPoints = patches[patchi].meshPoints();

            scalar sumThickness = 0;
            scalar sumNearWallThickness = 0;
            label sumNumLayers = 0;

            forAll(meshPoints, patchPointI)
            {
                label ppPointI = pp().meshPointMap()[meshPoints[patchPointI]];

                sumThickness += thickness[ppPointI];

                label nLay = patchNLayers[ppPointI];
                sumNumLayers += nLay;
                 if (nLay > 0)
                {
                    if (targetExpansion[ppPointI] == 1)
                    {
                        sumNearWallThickness += thickness[ppPointI]/nLay;
                    }
                    else
                    {
                        scalar s =
                            (1.0-pow(targetExpansion[ppPointI], nLay))
                          / (1.0-targetExpansion[ppPointI]);
                        sumNearWallThickness += thickness[ppPointI]/s;
                    }
                }
            }

            label totNPoints = returnReduce(meshPoints.size(), sumOp<label>());

            // For empty patches, totNPoints is 0.
            scalar avgThickness = 0;
            scalar avgNearWallThickness = 0;
            label  avgNumLayers = 0;

            if (totNPoints > 0)
            {
                std::tuple<scalar, scalar, label> sums{sumThickness, sumNearWallThickness, sumNumLayers};
                reduce(sums, ParallelOp<sumOp<scalar>, sumOp<scalar>, sumOp<label>>{});
                auto [summedThickness, summedNearWallThickness, summedNumLayers] = sums;

                avgThickness = summedThickness / totNPoints;
                avgNearWallThickness = summedNearWallThickness / totNPoints;
                avgNumLayers = summedNumLayers / totNPoints;
            }

            Info<< setf(ios_base::left) << setw(maxPatchNameLen)
                << patches[patchi].name() << setprecision(3)
                << " " << setw(8)
                << returnReduce(patches[patchi].size(), sumOp<label>())
                << " " << setw(6) << avgNumLayers
                << " " << setw(8) << avgNearWallThickness
                << "  " << setw(8) << avgThickness
                << endl;
        }
        Info<< endl;
    }

    // Current set of topology changes. (changing mesh clears out
    // polyTopoChange)
    polyTopoChange savedMeshMod(mesh.boundaryMesh().size());
    // First set of topology changes. (pp face splitting)
    polyTopoChange savedSplitMod(mesh.boundaryMesh().size());

    label nRefinedCells = mesh.nCells();
    reduce(nRefinedCells, sumOp<label>());

    {
        // Overall displacement field
        pointVectorField displacement
        (
            makeLayerDisplacementField
            (
                pointMesh::New(mesh),
                layerParams
            )
        );

        // Allocate run-time selectable mesh mover
        autoPtr<externalDisplacementMeshMover> medialAxisMoverPtr;
        {
            // Set up controls for meshMover
            dictionary combinedDict(layerParams.dict());
            // Add mesh quality constraints
            combinedDict.merge(motionDict);
            // Where to get minThickness from
            combinedDict.add("minThicknessName", minThickness.name());
            // Access to level 0 edge length
            combinedDict.add("edge0Len", edge0Len);

            const List<labelPair> internalBaffles
            (
                meshRefinement::subsetBaffles
                (
                    mesh,
                    internalFaceZones,
                    localPointRegion::findDuplicateFacePairs(mesh)
                )
            );

            // Take over patchDisp as boundary conditions on displacement
            // pointVectorField
            medialAxisMoverPtr = externalDisplacementMeshMover::New
            (
                layerParams.meshShrinker(),
                combinedDict,
                internalBaffles,
                displacement
            );

            if (controller_.mode() == meshControl::DRYRUN)
            {
                string errorMsg(FatalError.message());
                string IOerrorMsg(FatalIOError.message());

                if (errorMsg.size() || IOerrorMsg.size())
                {
                    //errorMsg = "[dryRun] " + errorMsg;
                    //errorMsg.replaceAll("\n", "\n[dryRun] ");
                    //IOerrorMsg = "[dryRun] " + IOerrorMsg;
                    //IOerrorMsg.replaceAll("\n", "\n[dryRun] ");

                    IOWarningInFunction(combinedDict)
                        << nl
                        << "Missing/incorrect required dictionary entries:"
                        << nl << nl
                        << IOerrorMsg.c_str() << nl
                        << errorMsg.c_str() << nl << nl
                        << "Exiting dry-run" << nl << endl;

                    if (Pstream::parRun())
                    {
                        Perr<< "\nFOAM parallel run exiting\n" << endl;
                        Pstream::exit(0);
                    }
                    else
                    {
                        Perr<< "\nFOAM exiting\n" << endl;
                        std::exit(0);
                    }
                }
            }
        }

        // Saved old points
        const pointField oldPoints(mesh.points());

        label nLayerIter = 1;
        label nPrevErrors = labelMax;


        boolList nonConvergedAfterNLayerIter(pp().localPoints().size(), false);

        label checkConvIter
        (
            controller_.mode() == meshControl::QUALITY ? 15 : 8
        );

        label minErr = labelMax;
        label minCounter = 0;

        while (true)
        {
            Info<< nl
                << "Layer addition iteration " << nLayerIter << nl
                << "--------------------------" << endl;
            scalar startTime = mesh.time().elapsedCpuTime();

            // Unset the extrusion at the pp.
            const dictionary& meshQualityDict =
            (
                nLayerIter < layerParams.nRelaxedIter()
                ? motionDict
                : motionDict.subDict("relaxed")
             );

            if (nLayerIter >= layerParams.nRelaxedIter())
            {
                Info<< "Switched to relaxed meshQuality constraints." << endl;
            }

            // Displacement acc. to pointnormals
            getPatchDisplacement
            (
                pp,
                thickness,
                minThickness,
                layerParams.grownUpIDs(), // patches to grow layers up
                patchDisp,
                patchNLayers,
                extrudeStatus
            );

            // Shrink mesh by displacement value first.
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            {
                pointField oldPatchPos(pp().localPoints());

                // We have patchDisp which is the outwards pointing
                // extrusion distance. Convert into an inwards pointing
                // shrink distance
                patchDisp = -patchDisp;

                // Take over patchDisp into pointDisplacement field and
                // adjust both for multi-patch constraints
                motionSmootherAlgo::setDisplacement
                (
                    patchIDs,
                    pp,
                    patchDisp,
                    displacement
                );

                // Move mesh
                // ~~~~~~~~~

                // Set up controls for meshMover
                dictionary combinedDict(layerParams.dict());
                // Add standard quality constraints
                combinedDict.merge(motionDict);
                // Add relaxed constraints (overrides standard ones)
                combinedDict.merge(meshQualityDict);
                // Where to get minThickness from
                combinedDict.add("minThicknessName", minThickness.name());
                // Access to level 0 edge length
                combinedDict.add("edge0Len", edge0Len);

                labelList checkFaces(identity(mesh.nFaces()));
                medialAxisMoverPtr().move
                (
                    combinedDict,
                    nAllowableErrors,
                    isConvexEdgePoint,
                    isConcaveEdgePoint,
                    checkFaces
                 );

                pp().clearGeom();

                    // Update patchDisp (since not all might have been honoured)
                patchDisp = oldPatchPos - pp().localPoints();
            }

            // Truncate displacements that are too small (this will do internal
            // ones, coupled ones have already been truncated by
            // syncPatchDisplacement)
            faceSet dummySet(mesh, "wrongPatchFaces", 0);

            label nTrunc = truncateDisplacement
            (
                globalFaces,
                edgeGlobalFaces,
                pp,
                minThickness,
                dummySet,
                patchDisp,
                patchNLayers,
                extrudeStatus
             );

            // Dump to .obj file for debugging.
            if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
            {
                dumpDisplacement
                (
                    mesh.time().path()/"layer_" + meshRefiner_.timeName(),
                    pp(),
                    patchDisp,
                    extrudeStatus
                 );

                const_cast<Time&>(mesh.time())++;
                Info<< "Writing shrunk mesh to time "
                    << meshRefiner_.timeName() << endl;

                // See comment in snappySnapDriver why we should not remove
                // using mesh.clearOut().

                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                        | meshRefinement::WRITEMESH
                     ),
                    mesh.time().path()/meshRefiner_.timeName(),
                    false
                 );
            }

            // Mesh topo change engine
            polyTopoChange splitMod(mesh);

            PackedList<1> isTriSplit(mesh.nFaces(), 0);

            if (controller_.topoChanges())
            {
                //Split faces if possible at layer terminations to
                //improve layer cell quality
                splitLayerTerminationFaces
                (
                    layerParams,
                    globalFaces,
                    edgeGlobalFaces,
                    pp,
                    minThickness,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus,
                    splitMod,
                    isTriSplit
                 );
            }

            // Determine per point/per face number of layers to extrude. Also
            // handles the slow termination of layers when going switching layers
            labelList nPatchPointLayers(pp().nPoints(), -1);
            labelList nPatchFaceLayers(pp().size(), -1);
            setupLayerInfoTruncation
            (
                pp,
                patchNLayers,
                extrudeStatus,
                layerParams.nBufferCellsNoExtrude(),
                layerParams.terminationStrategy(),
                layerParams.layerRecovery(),
                isConvexEdgePoint,
                isConcaveEdgePoint,
                nPatchPointLayers,
                nPatchFaceLayers
            );

            // Calculate displacement for final layer for addPatchLayer.
            // (layer of cells next to the original mesh)
            vectorField finalDisp(patchNLayers.size(), Zero);

            forAll(nPatchPointLayers, i)
            {
                if (nPatchPointLayers[i] > 0)
                {
                    if (targetExpansion[i] == 1.0)
                    {
                        finalDisp[i] = patchDisp[i] /nPatchPointLayers[i];
                    }
                    else
                    {
                        label nLay = nPatchPointLayers[i];

                        //if wanting fixed fch recalculate expansion ratio
                        // to achieve this
                        if (fixedfch[i])
                        {
                            scalar initStr = targetExpansion[i];
                            targetExpansion[i] = calculateStretching
                            (
                                initStr,
                                targetfch[i],
                                mag(patchDisp[i]),
                                nLay
                            );
                        }

                        if (targetExpansion[i] == 1.0)
                        {
                            finalDisp[i] = patchDisp[i] /nPatchPointLayers[i];
                        }
                        else
                        {
                            scalar h =
                                pow(targetExpansion[i], nLay - 1)
                                * (1.0 - targetExpansion[i])
                                / (1.0 - pow(targetExpansion[i], nLay));
                            finalDisp[i] = h*patchDisp[i];
                        }
                    }
                }
            }

            // Store mesh changes for if mesh is correct.
            savedSplitMod = splitMod;

            // clear primitive mesh demans driven data
            static_cast<primitiveMesh&>(mesh).clearOut();

            autoPtr<fvMesh> newMeshPtr;
            autoPtr<mapPolyMesh> splitMap = splitMod.makeMesh
            (
                newMeshPtr,
                IOobject
                (
                    //mesh.name()+"_split",
                    mesh.name(),
                    static_cast<polyMesh&>(mesh).instance(),
                    mesh.time(),  // register with runTime
                    IOobject::NO_READ,
                    static_cast<polyMesh&>(mesh).writeOpt()
                ),          // io params from original mesh but new name
                mesh,       // original mesh
                true        // parallel sync
             );
            fvMesh& newMesh = newMeshPtr();

            #if !defined( WIN32 ) && !defined( WIN64 )
            addProfiling(grow, "snappyHexMesh::layers::updateMesh");
            #endif

            List<labelPair> updatedBaffles(baffles.size());
            if (splitMap.valid())
            {
                const labelList& reverseFaceMap = splitMap().reverseFaceMap();
                forAll(baffles, i)
                {
                    label f0 = reverseFaceMap[baffles[i].first()];
                    label f1 = reverseFaceMap[baffles[i].second()];
                    updatedBaffles[i] = labelPair(f0, f1);
                }
            }

            autoPtr<indirectPrimitivePatch> ppSplitPtr
            (
                meshRefinement::makePatch
                (
                    newMesh,
                    layerPatchIDs()
                 )
             );
            indirectPrimitivePatch& ppSplit = ppSplitPtr();

            // Calculate displacement for first layer for addPatchLayer
            labelList nPatchFaceSplit(ppSplit.size(),-1);

            Map<label> mapMeshToPatch(pp().size());
            forAll(pp(), faceI)
            {
                mapMeshToPatch.insert(pp().addressing()[faceI], faceI);
            }

            Map<label> mapMeshPointToPatch(pp().nPoints());
            forAll(pp().localPoints(), pointI)
            {
                mapMeshPointToPatch.insert(pp().meshPoints()[pointI], pointI);
            }

            labelList splitFaceMap(ppSplit.size(), -1);
            forAll(ppSplit, faceI)
            {
                label meshFaceI = ppSplit.addressing()[faceI];
                label origMeshFaceI = splitMap().faceMap()[meshFaceI];

                if
                (
                    splitMap().reverseFaceMap()[origMeshFaceI] == meshFaceI
                    || isTriSplit.get(origMeshFaceI) == 1
                 )
                {
                    nPatchFaceSplit[faceI] =
                        nPatchFaceLayers[mapMeshToPatch[origMeshFaceI]];
                }
                else
                {
                    nPatchFaceSplit[faceI] = 0;
                }
                splitFaceMap[faceI] = mapMeshToPatch[origMeshFaceI];
            }

            scalarField invExpansionRatioSplit(ppSplit.nPoints());
            vectorField finalDispSplit(ppSplit.nPoints(), vector::zero);
            vectorField patchDispSplit(ppSplit.nPoints(), vector::zero);
            scalarField fchSplit(ppSplit.nPoints(), 0);

            labelList nPatchPointLayersSplit(ppSplit.nPoints(),-1);

            labelList patchSplitLevel(ppSplit.localFaces().size(), 0);

            forAll(ppSplit.localFaces(), faceI)
            {
                label meshFaceI = ppSplit.addressing()[faceI];
                label own = newMesh.faceOwner()[meshFaceI];
                patchSplitLevel[faceI] = cellLevel[own];
            }

            forAll(ppSplit.localPoints(), pointI)
            {
                label meshPointI = ppSplit.meshPoints()[pointI];
                label origMeshPointI = splitMap().pointMap()[meshPointI];
                label origPatchPointI = mapMeshPointToPatch[origMeshPointI];
                nPatchPointLayersSplit[pointI] =
                    nPatchPointLayers[origPatchPointI];

                finalDispSplit[pointI] = finalDisp[origPatchPointI];
                patchDispSplit[pointI] = patchDisp[origPatchPointI];
                invExpansionRatioSplit[pointI] =
                    1.0 / targetExpansion[origPatchPointI];
                scalar fld = mag(finalDispSplit[pointI]);
                fchSplit[pointI] = fld*
                    pow
                    (
                        invExpansionRatioSplit[pointI],
                        nPatchPointLayersSplit[pointI]-1
                     );
            }

            // Global face indices engine (split mesh)
            const globalIndex globalFacesSplit(newMesh.nFaces());

            // Determine extrudePatch.edgeFaces in global numbering (so across
            // coupled patches). This is used only to string up edges between
            // coupled faces (all edges between same (global)face indices get
            // extruded).
            labelListList edgeGlobalFacesSplit
            (
                addPatchCellLayer::globalEdgeFaces
                (
                    newMesh,
                    globalFacesSplit,
                    ppSplit
                 )
            );

            // Determine patches for extruded boundary edges. Calculates if any
            // additional processor patches need to be constructed.
            labelList edgePatchID;
            labelList edgeZoneID;
            boolList edgeFlip;
            labelList inflateFaceID;
            determineSidePatches
            (
                newMesh,

                globalFacesSplit,
                edgeGlobalFacesSplit,
                ppSplit,
                layerParams.grownUpIDs(),

                edgePatchID,
                edgeZoneID,
                edgeFlip,
                inflateFaceID
            );

            // Mesh topo change engine
            polyTopoChange meshMod(newMesh);

            // Grow layer of cells on to patch. Handles zero sized displacement.
            addPatchCellLayer addLayer(newMesh);
            newMesh.clearOut();

            // Add topo regardless of whether extrudeStatus is extruderemove.
            // Not add layer if patchDisp is zero.
            addLayer.setRefinement
            (
                globalFacesSplit,
                edgeGlobalFacesSplit,

                invExpansionRatioSplit,
                ppSplit,

                edgePatchID,    // boundary patch for extruded boundary edges
                edgeZoneID,     // zone for extruded edges
                edgeFlip,
                inflateFaceID,

                labelList(0),  // exposed patchIDs, not used for adding layers
                nPatchFaceSplit,          // layers per face
                nPatchPointLayersSplit,   // layers per point
                finalDispSplit,           // thickness of first layer
                layerParams.truncateFromWall(),
                meshMod
            );

            // Store mesh changes for if mesh is correct.
            savedMeshMod = meshMod;

            if (debug)
            {
                const_cast<Time&>(mesh.time())++;
            }

            // Apply the stored topo changes to the current mesh.
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(newMesh, false);
            // Update numbering on addLayer:
            // - cell/point labels to be newMesh.
            // - patchFaces to remain in oldMesh order.

            addLayer.updateMesh
            (
                map,
                identity(ppSplit.size()),
                identity(ppSplit.nPoints()),
                identity(ppSplit.size())
             );

            if (debug)
            {
                newMesh.write();
                faceSet wrongFaces(newMesh,"wrongFaces",newMesh.nFaces()/1000);
                motionSmoother::checkMesh
                (
                    false,
                    newMesh,
                    meshQualityDict,
                    wrongFaces
                );
                wrongFaces.write();
            }

            //- Get baffles in newMesh numbering. Note that we cannot detect
            //  baffles here since the points are duplicated
            List<labelPair> internalBaffles;
            {
                // From old mesh face to corresponding newMesh boundary face
                labelList meshToNewMesh(newMesh.nFaces(), -1);
                for
                (
                    label facei = newMesh.nInternalFaces();
                    facei < newMesh.nFaces();
                    facei++
                )
                {
                    label newMeshFacei = map().faceMap()[facei];

                    if (newMeshFacei != -1)
                    {
                        meshToNewMesh[newMeshFacei] = facei;
                    }
                }

                List<labelPair> newMeshBaffles(updatedBaffles.size());
                label newi = 0;
                forAll(updatedBaffles, i)
                {
                    const labelPair& p = updatedBaffles[i];
                    labelPair newMeshBaffle
                    (
                        meshToNewMesh[p[0]],
                        meshToNewMesh[p[1]]
                    );
                    if (newMeshBaffle[0] != -1 && newMeshBaffle[1] != -1)
                    {
                        newMeshBaffles[newi++] = newMeshBaffle;
                    }
                }
                newMeshBaffles.setSize(newi);

                internalBaffles = meshRefinement::subsetBaffles
                (
                    newMesh,
                    internalFaceZones,
                    newMeshBaffles
                );

                Info<< "Detected "
                    << returnReduce(internalBaffles.size(), sumOp<label>())
                    << " baffles across faceZones of type internal" << nl
                    << endl;
            }

            // Unset the extrusion at the pp.
            label nTotChanged = checkAndUnmark
            (
                addLayer,
                meshQualityDict,
                layerParams.additionalReporting(),
                internalBaffles,
                pp,
                newMesh,
                splitFaceMap,

                patchDisp,
                patchNLayers,
                extrudeStatus,

                nonConvergedAfterNLayerIter,
                nLayerIter
             );

            label nExtruded = countExtrusion(pp, extrudeStatus);
            label nTotFaces = returnReduce(pp().size(), sumOp<label>());
            Info<< "Extruding " << nExtruded
                << " out of " << nTotFaces
                << " faces (" << 100.0*nExtruded/nTotFaces << "%)."
                << " Removed extrusion at " << nTotChanged << " faces."
                << endl;


            //expand stencil of non converged points
            {
                syncTools::syncPointList
                (
                    mesh,
                    pp().meshPoints(),
                    nonConvergedAfterNLayerIter,
                    orEqOp<bool>(),
                    false               // null value
                 );

                boolList nonConvergedAfterNLayerIterOld =
                    nonConvergedAfterNLayerIter;

                forAll(pp().localPoints(), pointI)
                {
                    if (nonConvergedAfterNLayerIterOld[pointI])
                    {
                        labelList pFaces = pp().pointFaces()[pointI];
                        forAll(pFaces, pfI)
                        {
                            const face f = pp().localFaces()[pFaces[pfI]];
                            forAll(f, fp)
                            {
                                nonConvergedAfterNLayerIter[f[fp]] = true;
                            }
                        }
                    }
                }

                syncTools::syncPointList
                (
                    mesh,
                    pp().meshPoints(),
                    nonConvergedAfterNLayerIter,
                    orEqOp<bool>(),
                    false               // null value
                 );
            }

            if (controller_.mode() == meshControl::FAST)
            {
                if (nTotChanged < minErr)
                {
                    minErr = nTotChanged;
                    minCounter = 0;
                }
                else
                {
                   minCounter++;
                }
            }

            if
            (
                (
                    controller_.mode() == meshControl::DRYRUN
                    && nLayerIter == 1
                )
                ||
                (
                    nLayerIter == layerParams.maxLayerIter()
                    && !layerParams.noErrors()
                 )
                || (nTotChanged == 0 && nTrunc == 0)
                ||
                (
                    nLayerIter > checkConvIter &&  nTotChanged >= nPrevErrors
                    && !layerParams.noErrors()
                 )
                || minCounter > 2
            )
            {
                    //create patch layer info
                label maxBoundarySize = mesh.boundaryMesh().size();
                reduce(maxBoundarySize, maxOp<label>());
                scalarField fch(maxBoundarySize, 0);
                scalarField flt(maxBoundarySize, 0);
                scalarField expansion(maxBoundarySize, 0);
                scalarField tlt(maxBoundarySize, 0);
                scalarField nl(maxBoundarySize, 0);
                scalarField coverage(maxBoundarySize, 0);
                labelList nf(maxBoundarySize, 0);

                forAll(pp(), faceI)
                {
                    label meshFaceI = pp().addressing()[faceI];
                    face f = pp().localFaces()[faceI];

                    scalar aveFCH = 0;
                    scalar aveFLT = 0;
                    scalar aveEXP = 0;
                    scalar aveTLT = 0;
                    scalar aveNL = 0;
                    bool grown = false;


                    forAll(f, fp)
                    {
                        label ptI = f[fp];
                        scalar fl = mag(finalDisp[ptI]);

                        if (targetExpansion[ptI] > SMALL)
                        {
                            aveFCH += fl*pow
                            (
                                (1./targetExpansion[ptI]),
                                nPatchPointLayers[ptI]-1
                            );
                        }
                        aveFLT += fl;
                        aveEXP += targetExpansion[ptI];
                        aveTLT += mag(patchDisp[ptI]);
                        aveNL += nPatchPointLayers[ptI];

                        if (extrudeStatus[ptI] != NOEXTRUDE)
                        {
                            grown = true;
                        }
                    }

                    label patchID = mesh.boundaryMesh().whichPatch(meshFaceI);
                    if (f.size())
                    {
                        fch[patchID] += aveFCH/f.size();
                        flt[patchID] += aveFLT/f.size();
                        expansion[patchID] += aveEXP/f.size();
                        tlt[patchID] += aveTLT/f.size();
                        nl[patchID] += aveNL/f.size();
                        if (grown)
                        {
                            coverage[patchID]++;
                        }
                        nf[patchID]++;
                    }
                }

                Pstream::listCombineGather(fch, plusEqOp<scalar>());
                Pstream::listCombineScatter(fch);
                Pstream::listCombineGather(flt, plusEqOp<scalar>());
                Pstream::listCombineScatter(flt);
                Pstream::listCombineGather(expansion, plusEqOp<scalar>());
                Pstream::listCombineScatter(expansion);
                Pstream::listCombineGather(tlt, plusEqOp<scalar>());
                Pstream::listCombineScatter(tlt);
                Pstream::listCombineGather(nl, plusEqOp<scalar>());
                Pstream::listCombineScatter(nl);
                Pstream::listCombineGather(coverage, plusEqOp<scalar>());
                Pstream::listCombineScatter(coverage);
                Pstream::listCombineGather(nf, plusEqOp<label>());
                Pstream::listCombineScatter(nf);

                DynamicList<Tuple2<word, scalarField>>
                    patchLayerInfo(maxBoundarySize);
                forAll(nf, patchi)
                {
                    if (nf[patchi] > 0)
                    {
                        word name = mesh.boundaryMesh()[patchi].name();

                        scalarField linfo(7);
                        linfo[0] = nf[patchi];
                        linfo[1] = fch[patchi]/nf[patchi];
                        linfo[2] = flt[patchi]/nf[patchi];
                        linfo[3] = tlt[patchi]/nf[patchi];
                        linfo[4] = expansion[patchi]/nf[patchi];
                        linfo[5] = nl[patchi]/nf[patchi];
                        linfo[6] = 100.0*coverage[patchi]/nf[patchi];

                        patchLayerInfo.append
                        (
                            Tuple2<word, scalarField>(name,linfo)
                        );
                    }
                }
                patchLayerInfo.shrink();

                stats_.setLayerPatchInfo(patchLayerInfo);
                stats_.setLayerCoverage(100.0*nExtruded/nTotFaces);

                if
                (
                    layerParams.writeVTK()
                    && controller_.mode() != meshControl::DRYRUN
                )
                {
                    scalarField fchFaceSplit(nPatchFaceSplit.size(), 0);
                    forAll(fchFaceSplit, faceI)
                    {
                        const face& f = ppSplit.localFaces()[faceI];
                        forAll(f, fp)
                        {
                            fchFaceSplit[faceI] += fchSplit[f[fp]];
                        }
                        fchFaceSplit[faceI] /= f.size();
                    }

                    simpleVTKWriter layerVTK
                    (
                        ppSplit.localFaces(),
                        ppSplit.localPoints()+patchDispSplit
                     );

                    layerVTK.addFaceData("numLayers", nPatchFaceSplit);
                    layerVTK.addFaceData("fch", fchFaceSplit);
                    layerVTK.addFaceData("level", patchSplitLevel);

                    fileName linfoPath("layerInfo"/mesh.name());
                    if (Pstream::master())
                    {
                        if (!isDir(linfoPath))
                        {
                            mkDir(linfoPath);
                        }
                    }
                    if (multiLayerIter > 0)
                    {
                        layerVTK.write
                        (
                            linfoPath/"layerInfo-"
                            + Foam::name(multiLayerIter) + ".vtk"
                        );
                    }
                    else
                    {
                        layerVTK.write(linfoPath/"layerInfo.vtk");
                    }
                }

                // Add new processors patches ready for population by mesh
                // modifier.
                const polyBoundaryMesh& layerPolyPatches =
                    newMesh.boundaryMesh();
                polyBoundaryMesh& meshPolyPatches =
                    const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());

                forAll(layerPolyPatches, patchi)
                {
                    if
                    (
                        isA<processorPolyPatch>(layerPolyPatches[patchi])
                     )
                    {
                        const processorPolyPatch& procPatch =
                            refCast<const processorPolyPatch>
                            (layerPolyPatches[patchi]);
                        label indexI =
                            meshPolyPatches.findPatchID(procPatch.name());
                        if (indexI == -1)
                        {
                            dictionary patchDict;
                            patchDict.add("type", processorPolyPatch::typeName);
                            patchDict.add("myProcNo", Pstream::myProcNo());
                            patchDict.add
                            (
                                "neighbProcNo",
                                procPatch.neighbProcNo()
                            );
                            patchDict.add("nFaces", 0);
                            patchDict.add("startFace", mesh.nFaces());

                            meshRefiner_.appendPatch
                            (
                                mesh,
                                procPatch.name(),
                                patchDict
                            );
                        }
                    }
                }
                savedSplitMod.setNumPatches(mesh.boundaryMesh().size());
                savedMeshMod.setNumPatches(mesh.boundaryMesh().size());
                break;
            }
            nPrevErrors = nTotChanged;
            nLayerIter++;
            // Reset mesh points and start again
            mesh.movePoints(oldPoints);
            pp().clearGeom();
            medialAxisMoverPtr().movePoints(mesh.points());

            // Grow out region of non-extrusion
            for (label i = 0; i < layerParams.nGrow(); i++)
            {
                growNoExtrusion
                (
                    pp,
                    patchDisp,
                    patchNLayers,
                    extrudeStatus
                );
            }

            scalar endTime = mesh.time().elapsedCpuTime();
            Info<< "Layer iteration performed in "
                 << endTime - startTime << "s" << endl << endl;
        }
    }

    // At this point we have a (shrunk) mesh and a set of topology changes
    // which will make a valid mesh with layer. Apply these changes to the
    // current mesh.

    // Added clearOut because of problem with meshPhi
    mesh.clearOut();

    {
        // Apply the stored topo changes (split faces) to the current mesh.
        autoPtr<mapPolyMesh> mapSplit = savedSplitMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(mapSplit);

        // Move mesh (since morphing does not do this)
        if (mapSplit().hasMotionPoints())
        {
            mesh.movePoints(mapSplit().preMotionPoints());
        }

        meshRefiner_.updateMesh(mapSplit, labelList(0));

        if (mapSplit.valid())
        {
            const labelList& reverseFaceMap = mapSplit().reverseFaceMap();
            forAll(baffles, i)
            {
                label f0 = reverseFaceMap[baffles[i].first()];
                label f1 = reverseFaceMap[baffles[i].second()];
                baffles[i] = labelPair(f0, f1);
            }
        }
    }

    {
        // Apply the stored topo changes (layers) to the current mesh.
        autoPtr<mapPolyMesh> map = savedMeshMod.changeMesh
        (
            mesh,
            false,
            true,
            true,
            false
        );
        // Update fields
        mesh.updateMesh(map);

        // Move mesh (since morphing does not do this)
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

        meshRefiner_.updateMesh(map, labelList(0));

        // Use geometric detection of points-to-be-merged
        //  - detect any boundary face created from a duplicated face (=baffle)
        //  - on these mark any point created from a duplicated point
        if (returnReduce(pointToMaster.size(), sumOp<label>()))
        {
            // Estimate number of points-to-be-merged
            DynamicList<label> candidates(baffles.size()*4);
            // Mark whether old face was on baffle
            PackedBoolList oldBaffleFace(map().nOldFaces());
            forAll(baffles, i)
            {
                const labelPair& baffle = baffles[i];
                oldBaffleFace[baffle[0]] = true;
                oldBaffleFace[baffle[1]] = true;
            }

            boolList markedBafflePts(mesh.nPoints(), false);

            // Collect candidate if
            //  - point on boundary face originating from baffle
            //  - and point originating from duplicate
            for
            (
                label facei = mesh.nInternalFaces();
                facei < mesh.nFaces();
                facei++
            )
            {
                label oldFacei = map().faceMap()[facei];
                if (oldFacei != -1 && oldBaffleFace[oldFacei])
                {
                    const face& f = mesh.faces()[facei];
                    forAll(f, fp)
                    {
                        label pointi = f[fp];
                        label oldPointi = map().pointMap()[pointi];

                        if (pointToMaster[oldPointi] != -1)
                        {
                            markedBafflePts[pointi]= true;
                        }
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                markedBafflePts,
                orEqOp<bool>(),
                false              // null value
            );

            forAll(markedBafflePts, pointi)
            {
                if (markedBafflePts[pointi])
                {
                    candidates.append(pointi);
                }
            }

            // Do geometric merge. Ideally we'd like to use a topological
            // merge but we've thrown away all layer-wise addressing when
            // throwing away addPatchCellLayer engine. Also the addressing
            // is extremely complicated. There is no problem with merging
            // too many points; the problem would be if merging baffles.
            // Trust mergeZoneBaffles to do sufficient checks.

            scalar minEdgeLength = GREAT;
            for
            (
                label facei = mesh.nInternalFaces();
                facei < mesh.nFaces();
                facei++
            )
            {
                label oldFacei = map().faceMap()[facei];
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
            pointToMaster.setSize(mesh.nPoints());
            pointToMaster = -1;

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

        //Re sync of points to be merged
        labelList nAve(mesh.nPoints(), 0);
        pointField avePt(mesh.nPoints(), vector::zero);

        forAll(pointToMaster, pointi)
        {
            label otherPointI = pointToMaster[pointi];
            if (otherPointI != -1 && pointi != otherPointI)
            {
                nAve[pointi] += 2;
                nAve[otherPointI] += 2;
                point sPt = mesh.points()[pointi]+mesh.points()[otherPointI];
                avePt[pointi] += sPt;
                avePt[otherPointI] += sPt;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nAve,
            plusEqOp<label>(),
            label(0)      // null value
        );

        syncTools::syncPointList
        (
            mesh,
            avePt,
            plusEqOp<point>(),
            vector::zero      // null value
        );

        forAll(mesh.points(), pointi)
        {
            if (nAve[pointi]>0)
            {
                avePt[pointi] /= nAve[pointi];
            }
            else
            {
                avePt[pointi] = mesh.points()[pointi];
            }
        }
        mesh.movePoints(avePt);

        if (debug&meshRefinement::MESH || debug&meshRefinement::LAYERINFO)
        {
            OBJstream str
            (
                mesh.time().path()
              / "mergePoints_"
              + meshRefiner_.timeName()
              + ".obj"
            );
            Info<< "Points to be merged to " << str.name() << endl;
            forAll(pointToMaster, pointi)
            {
                label otherPointi = pointToMaster[pointi];
                if (otherPointi != -1)
                {
                    const point& pt = mesh.points()[pointi];
                    const point& otherPt = mesh.points()[otherPointi];
                    str.write(linePointRef(pt, otherPt));
                }
            }
        }

        autoPtr<mapPolyMesh> map = meshRefiner_.mergePoints(pointToMaster);
        if (map.valid())
        {
            Info<< "Merged points in = "
                << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
        }
    }

    if (mesh.faceZones().size() > 0)
    {
        // Merge any baffles
        Info<< "Converting baffles back into zoned faces ..."
            << endl;

        autoPtr<mapPolyMesh> map = meshRefiner_.mergeZoneBaffles
        (
            true,   // internal zones
            false,  // baffle zones
            false   // whether to update surface intersections
        );

        Info<< "Converted baffles in = "
            << meshRefiner_.mesh().time().cpuTimeIncrement()
            << " s\n" << nl << endl;
    }

    label totalNCells = mesh.nCells();
    reduce(totalNCells, sumOp<label>());
    label numLayerCells = totalNCells - nRefinedCells;

    Info<< "Total number of layer cells: " << numLayerCells << endl;

    stats_.setNumLayerCells(numLayerCells);
}


Foam::List<Foam::plane> Foam::snappyLayerDriver::calculateGrownUpPlanes
(
    const labelList& grownUpIDs
) const
{
    fvMesh& mesh = meshRefiner_.mesh();

    List<plane> grownUpPlanes(grownUpIDs.size());

    forAll(grownUpIDs, i)
    {
        label patchi = grownUpIDs[i];
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        point patchCentre(vector::zero);
        point patchNormal(vector::zero);
        scalar patchArea(0.);

        forAll(pp, ppI)
        {
            label faceI = pp.start() + ppI;
            point faceCentre = mesh.faceCentres()[faceI];
            vector faceNormal = mesh.faceAreas()[faceI];
            scalar faceArea = mag(faceNormal);

            patchNormal += faceNormal;
            patchCentre += faceCentre*faceArea;
            patchArea += faceArea;
        }

        reduce(
            std::tie(patchCentre, patchNormal, patchArea),
            ParallelOp<sumOp<point>, sumOp<point>, sumOp<scalar>>{}
        );

        if (patchArea > SMALL)
        {
            vector basePoint = patchCentre/patchArea;
            vector normalVector = patchNormal/patchArea;
            grownUpPlanes[i] = plane(basePoint, normalVector);
        }
    }

    return grownUpPlanes;
}


void Foam::snappyLayerDriver::projectToGrownUpPatches
(
    const layerParameters& layerParams,
    const labelList& grownUpIDs,
    const List<plane>& grownUpPlanes
) const
{
    Info<< "Projecting grown up patch points onto plane ..." << endl;

    fvMesh& mesh = meshRefiner_.mesh();

    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    scalar displRatio = layerParams.projectGrownUp();

    pointField displVec(mesh.nPoints(), vector::zero);

    forAll(grownUpIDs, i)
    {
        label patchi = grownUpIDs[i];
        const polyPatch& pp = mesh.boundaryMesh()[patchi];

        forAll(pp.meshPoints(), pI)
        {
            label meshPointI = pp.meshPoints()[pI];
            point pt = mesh.points()[meshPointI];

            vector ptVec = grownUpPlanes[i].nearestPoint(pt)-pt;

            scalar maxDispl = displRatio*edge0Len/(1<<pointLevel[meshPointI]);
            if (mag(ptVec) < maxDispl)
            {
                displVec[meshPointI] = ptVec;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        displVec,
        maxMagEqOp(),
        vector(GREAT, GREAT, GREAT)
     );

    mesh.movePoints(mesh.points()+displVec);
}


Foam::scalar Foam::snappyLayerDriver::calculateStretching
(
    const scalar initStr,
    const scalar a1,
    const scalar a2,
    const label nLay
) const
{
    if (nLay == 1)
    {
        if (a1 > a2)
        {
            return scalar(1);
        }
        else if (a1 < (a2 + SMALL) && a1 > (a2 - SMALL))
        {
            return scalar(1);
        }
    }
    else if (nLay > 1)
    {
        if (a1 >= a2)
        {
            return scalar(1);
        }
    }

    bool invStr(initStr <= 1.0 ? true : false);
    scalar xn = initStr;
    scalar x = -1.0;

    label repeat = 0;
    scalar tol = 0.001;

    while (!(repeat++ >= 100 || mag(x - xn) < tol))
    {
        x = xn;
        scalar f1 = a1 * pow(x, nLay) + a2*(1 - x) - a1;
        scalar f2 = a1 * nLay
            * pow(x, nLay -1) - a2;
        xn = x - f1/(f2+SMALL);
    }

    //first re-try with changed initial value if needed
    if ((xn < 1.0 + tol && xn > 1.0 - tol) || xn < 0)
    {
        xn = (invStr ? 0.0 : 10.0);
        x = -1.0;
        repeat = 0;
        while (!(repeat++ >= 100 || mag(x - xn) < tol))
        {
            x = xn;
            scalar f1 = a1 * pow(x, nLay)
                + a2*(1 - x) - a1;
            scalar f2 = a1 * nLay
                * pow(x, nLay -1) - a2;
            xn = x - f1/(f2+SMALL);
        }
    }

    //finally invert stretching if needed
    if ((xn < 1.0 + tol && xn > 1.0 - tol) || xn < 0)
    {
        xn = (invStr ? 10.0 : 0.0);
        x = -1.0;
        repeat = 0;
        while (!(repeat++ >= 100 || mag(x - xn) < tol))
        {
            x = xn;
            scalar f1 = a1 * pow(x, nLay)
                + a2*(1 - x) - a1;
            scalar f2 = a1 * nLay
                * pow(x, nLay -1) - a2;
            xn = x - f1/(f2+SMALL);
        }
    }

    if (xn < SMALL)
    {
        if (initStr < SMALL)
        {
            return scalar(1);
        }
        else
        {
            return initStr;
        }
    }
    else
    {
        return xn;
    }
}


bool Foam::snappyLayerDriver::doLayers
(
    const dictionary& shrinkDict,
    const dictionary& motionDict,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const label multiLayerIter
)
{
    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(layers, "snappyHexMesh::layers");
    #endif
    const fvMesh& mesh = meshRefiner_.mesh();

    Info<< nl
        << "Shrinking and layer addition phase" << nl
        << "----------------------------------" << nl
        << endl;

    Info<< "Using mesh parameters " << motionDict << nl << endl;

    //Read in layer parameters
    labelList numLayers;
    bool preBalance;
    {
        // Layer addition parameters
        layerParameters layerParams(shrinkDict, mesh.boundaryMesh());
        numLayers = layerParams.numLayers();
        preBalance = layerParams.preBalance();
    }

    // Per patch the number of layers (0 if no layer)

    // Patches that need to get a layer
    DynamicList<label> patchIDs(numLayers.size());
    DynamicList<label> zonePatchIDs(numLayers.size());
    label nFacesWithLayers = 0;

    forAll(numLayers, patchi)
    {
        if (numLayers[patchi] > 0)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            if (!pp.coupled())
            {
                patchIDs.append(patchi);
                nFacesWithLayers += mesh.boundaryMesh()[patchi].size();
            }
            else
            {
                WarningInFunction
                    << "Ignoring layers on coupled patch " << pp.name()
                    << endl;
            }
        }
    }

    // Add contributions from faceZones that get layers
    const faceZoneMesh& fZones = mesh.faceZones();

    forAll(fZones, zonei)
    {
        label mpi, spi;
        surfaceZonesInfo::faceZoneType fzType;
        meshRefiner_.getFaceZoneInfo(fZones[zonei].name(), mpi, spi, fzType);

        if (numLayers[mpi] > 0)
        {
            zonePatchIDs.append(mpi);
            nFacesWithLayers += fZones[zonei].size();
        }
        if (numLayers[spi] > 0)
        {
            zonePatchIDs.append(spi);
            nFacesWithLayers += fZones[zonei].size();
        }
    }
    zonePatchIDs.shrink();
    patchIDs.shrink();

    setZoneLayerPatchIDs(zonePatchIDs);
    setLayerPatchIDs(patchIDs);
    if (returnReduce(nFacesWithLayers, sumOp<label>()) == 0)
    {
        Info<< nl << "No layers to generate ..." << endl;

        return false;
    }
    else
    {
        bool faceZoneOnCoupledFace = false;

        if (!preBalance)
        {
            // Check if there are faceZones on processor boundaries. This
            // requires balancing to move them off the processor boundaries.

            // Is face on a faceZone
            PackedBoolList isExtrudedZoneFace(mesh.nFaces());
            {
                // Add contributions from faceZones that get layers
                const faceZoneMesh& fZones = mesh.faceZones();
                forAll(fZones, zonei)
                {
                    const faceZone& fZone = fZones[zonei];
                    const word& fzName = fZone.name();

                    label mpi, spi;
                    surfaceZonesInfo::faceZoneType fzType;
                    meshRefiner_.getFaceZoneInfo(fzName, mpi, spi, fzType);

                    if (numLayers[mpi] > 0 || numLayers[spi])
                    {
                        forAll(fZone, i)
                        {
                            isExtrudedZoneFace[fZone[i]] = true;
                        }
                    }
                }
            }

            PackedBoolList intOrCoupled
            (
                syncTools::getInternalOrCoupledFaces(mesh)
            );

            for
            (
                label facei = mesh.nInternalFaces();
                facei < mesh.nFaces();
                facei++
            )
            {
                if (intOrCoupled[facei] && isExtrudedZoneFace[facei])
                {
                    faceZoneOnCoupledFace = true;
                    break;
                }
            }

            reduce(faceZoneOnCoupledFace, orOp<bool>());
        }

        // Balance
        if (Pstream::parRun() && (preBalance || faceZoneOnCoupledFace))
        {
            Info<< nl
                << "Doing initial balancing" << nl
                << "-----------------------" << nl
                << endl;

            scalarField cellWeights(mesh.nCells(), 1);
            forAll(numLayers, patchi)
            {
                if (numLayers[patchi] > 0)
                {
                    const polyPatch& pp = mesh.boundaryMesh()[patchi];
                    forAll(pp.faceCells(), i)
                    {
                        cellWeights[pp.faceCells()[i]] += numLayers[patchi];
                    }
                }
            }

            // Add contributions from faceZones that get layers
            const faceZoneMesh& fZones = mesh.faceZones();
            forAll(fZones, zonei)
            {
                const faceZone& fZone = fZones[zonei];
                const word& fzName = fZone.name();

                label mpi, spi;
                surfaceZonesInfo::faceZoneType fzType;
                meshRefiner_.getFaceZoneInfo(fzName, mpi, spi, fzType);

                if (numLayers[mpi] > 0)
                {
                    // Get the owner side for unflipped faces, neighbour side
                    // for flipped ones
                    const labelList& cellIDs = fZone.slaveCells();
                    forAll(cellIDs, i)
                    {
                        if (cellIDs[i] >= 0)
                        {
                            cellWeights[cellIDs[i]] += numLayers[mpi];
                        }
                    }
                }
                if (numLayers[spi] > 0)
                {
                    const labelList& cellIDs = fZone.masterCells();
                    forAll(cellIDs, i)
                    {
                        if (cellIDs[i] >= 0)
                        {
                            cellWeights[cellIDs[i]] += numLayers[mpi];
                        }
                    }
                }
            }
            autoPtr<mapDistributePolyMesh> map = meshRefiner_.balance
            (
                false,
                true,           // keepZoneFaces
                false,
                cellWeights,
                decomposer,
                distributor,
                false
            );
        }

        // Check that outside of mesh is not multiply connected.
        checkMeshManifold();

        // Layer addition parameters
        layerParameters layerParams(shrinkDict, mesh.boundaryMesh());

        // Do all topo changes
        addLayers
        (
            layerParams,
            motionDict,
            patchIDs,
            multiLayerIter
        );

        return true;
    }
}


Foam::PackedBoolList Foam::snappyLayerDriver::getMasterPPEdges
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp
)
{
    PackedBoolList isMasterEdge(mesh.nEdges(),false);

    labelList procID(mesh.nEdges(), -1);

    labelList meshEdges(pp.nEdges());

    forAll(pp.edges(), patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];
        label v0 = pp.meshPoints()[e[0]];
        label v1 = pp.meshPoints()[e[1]];
        meshEdges[patchEdgeI] = meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[v0],
            v0,
            v1
        );
        procID[meshEdges[patchEdgeI]] = Pstream::myProcNo();
    }

    syncTools::syncEdgeList
    (
        mesh,
        procID,
        maxEqOp<label>(),
        label(-1)            // initial value
     );

    forAll(pp.edges(), patchEdgeI)
    {
        label meshEdgeI = meshEdges[patchEdgeI];

        if (procID[meshEdgeI] == Pstream::myProcNo())
        {
            isMasterEdge[meshEdgeI] = true;
        }
    }

    return isMasterEdge;
}


Foam::PackedBoolList Foam::snappyLayerDriver::getMasterPPPoints
(
    const indirectPrimitivePatch& pp
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    PackedBoolList isMasterPt(mesh.nPoints(),false);

    labelList procID(mesh.nPoints(), -1);

    labelList meshEdges(pp.nPoints());

    forAll(pp.localPoints(), patchPointI)
    {
        procID[pp.meshPoints()[patchPointI]] = Pstream::myProcNo();
    }

    syncTools::syncPointList
    (
        mesh,
        procID,
        maxEqOp<label>(),
        label(-1)            // initial value
     );

    forAll(pp.localPoints(), patchPointI)
    {
        label meshPointI = pp.meshPoints()[patchPointI];

        if (procID[meshPointI] == Pstream::myProcNo())
        {
            isMasterPt[meshPointI] = true;
        }
    }

    return isMasterPt;
}


// ************************************************************************* //
