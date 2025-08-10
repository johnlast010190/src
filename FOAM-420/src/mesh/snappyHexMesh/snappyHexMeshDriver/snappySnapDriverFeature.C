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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "snappyHexMeshDriver/snappySnapDriver.H"
#include "snappyHexMeshDriver/snappyLayerDriver.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "fvMesh/fvMesh.H"
#include "surfaceFormats/obj/OBJstream.H"
#include "motionSmoother/motionSmoother.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "refinementFeatures/refinementFeatures.H"
#include "global/unitConversion/unitConversion.H"
#include "meshes/primitiveShapes/plane/plane.H"
#include "edgeMesh/featureEdgeMesh/featureEdgeMesh.H"
#include "indexedOctree/treeDataPoint.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "snappyHexMeshDriver/snapParameters/snapParameters.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "meshes/meshShapes/cell/pyramidPointFaceRef.H"
#include "regionSplit/localPointRegion.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "autoOptimize/autoOptimize.H"
#include "meshes/polyMesh/syncTools/dummyTransform.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"

#ifdef FOAM_USE_TBB
  #include "include/TBBTimer.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<class T>
    class listPlusEqOp
    {
    public:

        void operator()(List<T>& x, const List<T>& y) const
        {
            label sz = x.size();
            x.setSize(sz+y.size());
            forAll(y, i)
            {
                x[sz++] = y[i];
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::snappySnapDriver::isFeaturePoint
(
    const scalar featureCos,
    const indirectPrimitivePatch& pp,
    const PackedBoolList& isFeatureEdge,
    const label pointi
) const
{
    const pointField& points = pp.localPoints();
    const edgeList& edges = pp.edges();
    const labelList& pEdges = pp.pointEdges()[pointi];

    label nFeatEdges = 0;

    forAll(pEdges, i)
    {
        if (isFeatureEdge[pEdges[i]])
        {
            nFeatEdges++;

            for (label j = i+1; j < pEdges.size(); j++)
            {
                if (isFeatureEdge[pEdges[j]])
                {
                    const edge& ei = edges[pEdges[i]];
                    const edge& ej = edges[pEdges[j]];

                    const point& p = points[pointi];
                    const point& pi = points[ei.otherVertex(pointi)];
                    const point& pj = points[ej.otherVertex(pointi)];

                    vector vi = p-pi;
                    scalar viMag = mag(vi);

                    vector vj = pj-p;
                    scalar vjMag = mag(vj);

                    if
                    (
                        viMag > SMALL
                     && vjMag > SMALL
                     && ((vi/viMag & vj/vjMag) < featureCos)
                    )
                    {
                        return true;
                    }
                }
            }
        }
    }

    if (nFeatEdges == 1)
    {
        // End of feature-edge string
        return true;
    }

    return false;
}


void Foam::snappySnapDriver::smoothAndConstrain
(
    const PackedBoolList& isPatchMasterEdge,
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const List<pointConstraint>& constraints,
    vectorField& disp
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    for (label avgIter = 0; avgIter < 20; avgIter++)
    {
        // Calculate average displacement of neighbours
        // - unconstrained (i.e. surface) points use average of all
        //   neighbouring points
        // - from testing it has been observed that it is not beneficial
        //   to have edge constrained points use average of all edge or point
        //   constrained neighbours since they're already attracted to
        //   the nearest point on the feature.
        //   Having them attract to point-constrained neighbours does not
        //   make sense either since there is usually just one of them so
        //   it severely distorts it.
        // - same for feature points. They are already attracted to the
        //   nearest feature point.

        vectorField dispSum(pp.nPoints(), Zero);
        labelList dispCount(pp.nPoints(), 0);

        const labelListList& pointEdges = pp.pointEdges();
        const edgeList& edges = pp.edges();

        forAll(pointEdges, pointi)
        {
            const labelList& pEdges = pointEdges[pointi];

            label nConstraints = constraints[pointi].first();

            if (nConstraints <= 1)
            {
                forAll(pEdges, i)
                {
                    label edgei = pEdges[i];

                    if (isPatchMasterEdge[edgei])
                    {
                        label nbrPointi = edges[edgei].otherVertex(pointi);
                        if (constraints[nbrPointi].first() >= nConstraints)
                        {
                            dispSum[pointi] += disp[nbrPointi];
                            dispCount[pointi]++;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispSum,
            plusEqOp<point>(),
            vector::zero,
            mapDistribute::transform()
        );
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            dispCount,
            plusEqOp<label>(),
            label(0),
            mapDistribute::transform()
        );

        // Constraints
        forAll(constraints, pointi)
        {
            if (dispCount[pointi] > 0)
            {
                // Mix my displacement with neighbours' displacement
                disp[pointi] =
                    0.5
                   *(disp[pointi] + dispSum[pointi]/dispCount[pointi]);
            }
        }
    }
}


void Foam::snappySnapDriver::calcNearestFace
(
    const label iter,
    const indirectPrimitivePatch& pp,
    const scalarField& faceSnapDist,
    vectorField& faceDisp,
    vectorField& faceSurfaceNormal,
    labelList& faceSurfaceGlobalRegion
    //vectorField& faceRotation
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    // Displacement and orientation per pp face.
    faceDisp.setSize(pp.size());
    faceDisp = Zero;
    faceSurfaceNormal.setSize(pp.size());
    faceSurfaceNormal = Zero;
    faceSurfaceGlobalRegion.setSize(pp.size());
    faceSurfaceGlobalRegion = -1;

    // Divide surfaces into zoned and unzoned
    const labelList zonedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones());
    const labelList unzonedSurfaces =
        surfaceZonesInfo::getUnnamedSurfaces(surfaces.surfZones());

    // Per pp face the current surface snapped to
    labelList snapSurf(pp.size(), -1);


    // Do zoned surfaces
    // ~~~~~~~~~~~~~~~~~
    // Zoned faces only attract to corresponding surface

    // Extract faces per zone
    const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

    forAll(zonedSurfaces, i)
    {
        label zoneSurfi = zonedSurfaces[i];

        const word& faceZoneName = surfZones[zoneSurfi].faceZoneName();

        // Get indices of faces on pp that are also in zone
        label zonei = mesh.faceZones().findZoneID(faceZoneName);
        if (zonei == -1)
        {
            FatalErrorInFunction
                << "Problem. Cannot find zone " << faceZoneName
                << exit(FatalError);
        }
        const faceZone& fZone = mesh.faceZones()[zonei];
        PackedBoolList isZonedFace(mesh.nFaces());
        forAll(fZone, i)
        {
            isZonedFace[fZone[i]] = 1;
        }

        DynamicList<label> ppFaces(fZone.size());
        DynamicList<label> meshFaces(fZone.size());
        forAll(pp.addressing(), i)
        {
            if (isZonedFace[pp.addressing()[i]])
            {
                snapSurf[i] = zoneSurfi;
                ppFaces.append(i);
                meshFaces.append(pp.addressing()[i]);
            }
        }

        //Pout<< "For faceZone " << fZone.name()
        //    << " found " << ppFaces.size() << " out of " << pp.size()
        //    << endl;

        pointField fc
        (
            indirectPrimitivePatch
            (
                IndirectList<face>(mesh.faces(), meshFaces),
                mesh.points()
            ).faceCentres()
        );

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        labelList hitRegion;
        vectorField hitNormal;
        surfaces.findNearestRegion
        (
            labelList(1, zoneSurfi),
            fc,
            sqr(faceSnapDist),// sqr of attract dist
            hitSurface,
            hitInfo,
            hitRegion,
            hitNormal
        );

        forAll(hitInfo, hiti)
        {
            if (hitInfo[hiti].hit())
            {
                label facei = ppFaces[hiti];
                faceDisp[facei] = hitInfo[hiti].hitPoint() - fc[hiti];
                faceSurfaceNormal[facei] = hitNormal[hiti];
                faceSurfaceGlobalRegion[facei] = surfaces.globalRegion
                (
                    hitSurface[hiti],
                    hitRegion[hiti]
                );
            }
            else
            {
                WarningInFunction
                    << "Did not find surface near face centre " << fc[hiti]
                    << endl;
            }
        }
    }


    // Do unzoned surfaces
    // ~~~~~~~~~~~~~~~~~~~
    // Unzoned faces attract to any unzoned surface

    DynamicList<label> ppFaces(pp.size());
    DynamicList<label> meshFaces(pp.size());
    forAll(pp.addressing(), i)
    {
        if (snapSurf[i] == -1)
        {
            ppFaces.append(i);
            meshFaces.append(pp.addressing()[i]);
        }
    }
    //Pout<< "Found " << ppFaces.size() << " unzoned faces out of "
    //   << pp.size() << endl;

    pointField fc
    (
        indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), meshFaces),
            mesh.points()
        ).faceCentres()
    );

    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    labelList hitRegion;
    vectorField hitNormal;
    surfaces.findNearestRegion
    (
        unzonedSurfaces,
        fc,
        sqr(faceSnapDist),// sqr of attract dist
        hitSurface,
        hitInfo,
        hitRegion,
        hitNormal
    );

    forAll(hitInfo, hiti)
    {
        if (hitInfo[hiti].hit())
        {
            label facei = ppFaces[hiti];
            faceDisp[facei] = hitInfo[hiti].hitPoint() - fc[hiti];
            faceSurfaceNormal[facei] = hitNormal[hiti];
            faceSurfaceGlobalRegion[facei] = surfaces.globalRegion
            (
                hitSurface[hiti],
                hitRegion[hiti]
            );
        }
        else
        {
            WarningInFunction
                << "Did not find surface near face centre " << fc[hiti]
                << endl;
        }
    }


    //// Determine rotation
    //// ~~~~~~~~~~~~~~~~~~
    //
    //// Determine rotation axis
    //faceRotation.setSize(pp.size());
    //faceRotation = Zero;
    //
    //forAll(faceRotation, facei)
    //{
    //    // Note: extend to >180 degrees checking
    //    faceRotation[facei] =
    //        pp.faceNormals()[facei]
    //      ^ faceSurfaceNormal[facei];
    //}
    //
    //if (debug&meshRefinement::ATTRACTION)
    //{
    //    dumpMove
    //    (
    //        mesh.time().path()
    //      / "faceDisp_" + name(iter) + ".obj",
    //        pp.faceCentres(),
    //        pp.faceCentres() + faceDisp
    //    );
    //    dumpMove
    //    (
    //        mesh.time().path()
    //      / "faceRotation_" + name(iter) + ".obj",
    //        pp.faceCentres(),
    //        pp.faceCentres() + faceRotation
    //    );
    //}
}


// Collect (possibly remote) per point data of all surrounding faces
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// - faceSurfaceNormal
// - faceDisp
// - faceCentres&faceNormal
void Foam::snappySnapDriver::calcNearestFacePointProperties
(
    const label iter,
    const indirectPrimitivePatch& pp,

    const vectorField& faceDisp,
    const vectorField& faceSurfaceNormal,
    const labelList& faceSurfaceGlobalRegion,

    List<List<point>>& pointFaceSurfNormals,
    List<List<point>>& pointFaceDisp,
    List<List<point>>& pointFaceCentres,
    List<labelList>&    pointFacePatchID
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));


    // For now just get all surrounding face data. Expensive - should just
    // store and sync data on coupled points only
    // (see e.g PatchToolsNormals.C)

    pointFaceSurfNormals.setSize(pp.nPoints());
    pointFaceDisp.setSize(pp.nPoints());
    pointFaceCentres.setSize(pp.nPoints());
    pointFacePatchID.setSize(pp.nPoints());

    // Fill local data
    forAll(pp.pointFaces(), pointi)
    {
        const labelList& pFaces = pp.pointFaces()[pointi];

        // Count valid face normals
        label nFaces = 0;
        forAll(pFaces, i)
        {
            label facei = pFaces[i];
            if (isMasterFace[facei] && faceSurfaceGlobalRegion[facei] != -1)
            {
                nFaces++;
            }
        }


        List<point>& pNormals = pointFaceSurfNormals[pointi];
        pNormals.setSize(nFaces);
        List<point>& pDisp = pointFaceDisp[pointi];
        pDisp.setSize(nFaces);
        List<point>& pFc = pointFaceCentres[pointi];
        pFc.setSize(nFaces);
        labelList& pFid = pointFacePatchID[pointi];
        pFid.setSize(nFaces);

        nFaces = 0;
        forAll(pFaces, i)
        {
            label facei = pFaces[i];
            label globalRegioni = faceSurfaceGlobalRegion[facei];

            if (isMasterFace[facei] && globalRegioni != -1)
            {
                pNormals[nFaces] = faceSurfaceNormal[facei];
                pDisp[nFaces] = faceDisp[facei];
                pFc[nFaces] = pp.faceCentres()[facei];
                pFid[nFaces] = globalToMasterPatch_[globalRegioni];
                nFaces++;
            }
        }
    }


    // Collect additionally 'normal' boundary faces for boundaryPoints of pp
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // points on the boundary of pp should pick up non-pp normals
    // as well for the feature-reconstruction to behave correctly.
    // (the movement is already constrained outside correctly so it
    //  is only that the unconstrained attraction vector is calculated
    //  correctly)
    {
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        labelList patchID(pbm.patchID());

        // Unmark all non-coupled boundary faces
        forAll(pbm, patchi)
        {
            const polyPatch& pp = pbm[patchi];

            if (pp.coupled() || isA<emptyPolyPatch>(pp))
            {
                forAll(pp, i)
                {
                    label meshFacei = pp.start()+i;
                    patchID[meshFacei-mesh.nInternalFaces()] = -1;
                }
            }
        }

        // Remove any meshed faces
        forAll(pp.addressing(), i)
        {
            label meshFacei = pp.addressing()[i];
            patchID[meshFacei-mesh.nInternalFaces()] = -1;
        }



        // See if edge of pp uses any non-meshed boundary faces. If so add the
        // boundary face as additional constraint. Note that we account for
        // both 'real' boundary edges and boundary edge of baffles

        const labelList bafflePair
        (
            localPointRegion::findDuplicateFaces(mesh, pp.addressing())
        );


        // Mark all points on 'boundary' edges
        PackedBoolList isBoundaryPoint(pp.nPoints());

        const labelListList& edgeFaces = pp.edgeFaces();
        const edgeList& edges = pp.edges();

        forAll(edgeFaces, edgei)
        {
            const edge& e = edges[edgei];
            const labelList& eFaces = edgeFaces[edgei];

            if (eFaces.size() == 1)
            {
                // 'real' boundary edge
                isBoundaryPoint[e[0]] = true;
                isBoundaryPoint[e[1]] = true;
            }
            else if (eFaces.size() == 2 && bafflePair[eFaces[0]] == eFaces[1])
            {
                // 'baffle' boundary edge
                isBoundaryPoint[e[0]] = true;
                isBoundaryPoint[e[1]] = true;
            }
        }


        // Construct labelList equivalent of meshPointMap
        labelList meshToPatchPoint(mesh.nPoints(), -1);
        forAll(pp.meshPoints(), pointi)
        {
            meshToPatchPoint[pp.meshPoints()[pointi]] = pointi;
        }

        forAll(patchID, bFacei)
        {
            label patchi = patchID[bFacei];

            if (patchi != -1)
            {
                label facei = mesh.nInternalFaces()+bFacei;
                const face& f = mesh.faces()[facei];

                forAll(f, fp)
                {
                    label pointi = meshToPatchPoint[f[fp]];

                    if (pointi != -1 && isBoundaryPoint[pointi])
                    {
                        List<point>& pNormals = pointFaceSurfNormals[pointi];
                        List<point>& pDisp = pointFaceDisp[pointi];
                        List<point>& pFc = pointFaceCentres[pointi];
                        labelList& pFid = pointFacePatchID[pointi];

                        const point& pt = mesh.points()[f[fp]];
                        vector fn = mesh.faceAreas()[facei];

                        pNormals.append(fn/mag(fn));
                        pDisp.append(mesh.faceCentres()[facei]-pt);
                        pFc.append(mesh.faceCentres()[facei]);
                        pFid.append(patchi);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceSurfNormals,
        listPlusEqOp<point>(),
        List<point>(),
        mapDistribute::transform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceDisp,
        listPlusEqOp<point>(),
        List<point>(),
        mapDistribute::transform()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFaceCentres,
        listPlusEqOp<point>(),
        List<point>(),
        mapDistribute::transformPosition()
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointFacePatchID,
        listPlusEqOp<label>(),
        List<label>()
    );


    // Sort the data according to the face centres. This is only so we get
    // consistent behaviour serial and parallel.
    labelList visitOrder;
    forAll(pointFaceDisp, pointi)
    {
        List<point>& pNormals = pointFaceSurfNormals[pointi];
        List<point>& pDisp = pointFaceDisp[pointi];
        List<point>& pFc = pointFaceCentres[pointi];
        labelList& pFid = pointFacePatchID[pointi];

        sortedOrder(mag(pFc)(), visitOrder);

        pNormals = List<point>(pNormals, visitOrder);
        pDisp = List<point>(pDisp, visitOrder);
        pFc = List<point>(pFc, visitOrder);
        pFid = UIndirectList<label>(pFid, visitOrder)();
    }
}


// Gets passed in offset to nearest point on feature edge. Calculates
// if the point has a different number of faces on either side of the feature
// and if so attracts the point to that non-dominant plane.
void Foam::snappySnapDriver::correctAttraction
(
    const DynamicList<point>& surfacePoints,
    const DynamicList<label>& surfaceCounts,
    const point& edgePt,
    const vector& edgeNormal,       // normalised normal
    const point& pt,

    vector& edgeOffset              // offset from pt to point on edge
) const
{
    // Tangential component along edge
    scalar tang = ((pt-edgePt)&edgeNormal);

    labelList order;
    Foam::sortedOrder(surfaceCounts, order);

    if (order[0] < order[1])
    {
        // There is a non-dominant plane. Use the point on the plane to
        // attract to.
        vector attractD = surfacePoints[order[0]]-edgePt;
        // Tangential component along edge
        scalar tang2 = (attractD&edgeNormal);
        // Normal component
        attractD -= tang2*edgeNormal;
        // Calculate fraction of normal distances
        scalar magAttractD = mag(attractD);
        scalar fraction = magAttractD/(magAttractD+mag(edgeOffset));

        point linePt =
            edgePt
          + ((1.0-fraction)*tang2 + fraction*tang)*edgeNormal;
        edgeOffset = linePt-pt;
    }
}


Foam::pointIndexHit Foam::snappySnapDriver::findMultiPatchPoint
(
    const point& pt,
    const labelList& patchIDs,
    const List<point>& faceCentres
) const
{
    // Determine if multiple patchIDs
    if (patchIDs.size())
    {
        label patch0 = patchIDs[0];

        for (label i = 1; i < patchIDs.size(); i++)
        {
            if (patchIDs[i] != patch0)
            {
                return pointIndexHit(true, pt, labelMax);
            }
        }
    }
    return pointIndexHit(false, Zero, labelMax);
}


Foam::label Foam::snappySnapDriver::findNormal
(
    const scalar featureCos,
    const vector& n,
    const DynamicList<vector>& surfaceNormals
) const
{
    label index = -1;

    forAll(surfaceNormals, j)
    {
        scalar cosAngle = (n&surfaceNormals[j]);

        if
        (
            (cosAngle >= featureCos)
         || (cosAngle < (-1+0.001)) // triangle baffles
        )
        {
            index = j;
            break;
        }
    }
    return index;
}


// Detect multiple patches. Returns pointIndexHit:
// - false, index=-1 : single patch
// - true , index=0  : multiple patches but on different normals planes
//                     (so geometric feature edge is also a region edge)
// - true , index=1  : multiple patches on same normals plane i.e. flat region
//                     edge
Foam::pointIndexHit Foam::snappySnapDriver::findMultiPatchPoint
(
    const point& pt,
    const labelList& patchIDs,
    const DynamicList<vector>& surfaceNormals,
    const labelList& faceToNormalBin
) const
{
    if (patchIDs.empty())
    {
        return pointIndexHit(false, pt, -1);
    }

    // Detect single patch situation (to avoid allocation)
    label patch0 = patchIDs[0];

    for (label i = 1; i < patchIDs.size(); i++)
    {
        if (patchIDs[i] != patch0)
        {
            patch0 = -1;
            break;
        }
    }

    if (patch0 >= 0)
    {
        // Single patch
        return pointIndexHit(false, pt, -1);
    }
    else
    {
        if (surfaceNormals.size() == 1)
        {
            // Same normals plane, flat region edge.
            return pointIndexHit(true, pt, 1);
        }
        else
        {
            // Detect per normals bin
            labelList normalToPatch(surfaceNormals.size(), -1);
            forAll(faceToNormalBin, i)
            {
                if (faceToNormalBin[i] != -1)
                {
                    label& patch = normalToPatch[faceToNormalBin[i]];
                    if (patch == -1)
                    {
                        // First occurence
                        patch = patchIDs[i];
                    }
                    else if (patch == -2)
                    {
                        // Already marked as being on multiple patches
                    }
                    else if (patch != patchIDs[i])
                    {
                        // Mark as being on multiple patches
                        patch = -2;
                    }
                }
            }

            forAll(normalToPatch, normali)
            {
                if (normalToPatch[normali] == -2)
                {
                    // Multiple patches on same normals plane, flat region
                    // edge
                    return pointIndexHit(true, pt, 1);
                }
            }

            // All patches on either side of geometric feature anyway
            return pointIndexHit(true, pt, 0);
        }
    }
}


void Foam::snappySnapDriver::writeStats
(
    const indirectPrimitivePatch& pp,
    const PackedBoolList& isPatchMasterPoint,
    const List<pointConstraint>& patchConstraints
) const
{
    label nMasterPoints = 0;
    label nPlanar = 0;
    label nEdge = 0;
    label nPoint = 0;

    forAll(patchConstraints, pointi)
    {
        if (isPatchMasterPoint[pointi])
        {
            nMasterPoints++;

            if (patchConstraints[pointi].first() == 1)
            {
                nPlanar++;
            }
            else if (patchConstraints[pointi].first() == 2)
            {
                nEdge++;
            }
            else if (patchConstraints[pointi].first() == 3)
            {
                nPoint++;
            }
        }
    }

    reduce(
        std::tie(nMasterPoints, nPlanar, nEdge, nPoint),
        UniformParallelOp<sumOp<label>, 4>{}
    );

    Info<< "total master points :" << nMasterPoints
        << " of which attracted to :" << nl
        << "    feature point   : " << nPoint << nl
        << "    feature edge    : " << nEdge << nl
        << "    nearest surface : " << nPlanar << nl
        << "    rest            : " << nMasterPoints-nPoint-nEdge-nPlanar
        << nl
        << endl;
}


void Foam::snappySnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    const label pointi,

    const List<List<point>>& pointFaceSurfNormals,
    const List<List<point>>& pointFaceDisp,
    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    DynamicList<point>& surfacePoints,
    DynamicList<vector>& surfaceNormals,
    labelList& faceToNormalBin,

    vector& patchAttraction,
    pointConstraint& patchConstraint
) const
{
    patchAttraction = Zero;
    patchConstraint = pointConstraint();

    const List<point>& pfSurfNormals = pointFaceSurfNormals[pointi];
    const List<point>& pfDisp = pointFaceDisp[pointi];
    const List<point>& pfCentres = pointFaceCentres[pointi];

    // Bin according to surface normal
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Bins of differing normals:
    //  - one normal   : flat(tish) surface
    //  - two normals  : geometric feature edge
    //  - three normals: geometric feature point
    //  - four normals : too complex a feature
    surfacePoints.clear();
    surfaceNormals.clear();

    //- From face to above normals bin
    faceToNormalBin.setSize(pfDisp.size());
    faceToNormalBin = -1;

    forAll(pfSurfNormals, i)
    {
        const point& fc = pfCentres[i];
        const vector& fSNormal = pfSurfNormals[i];
        const vector& fDisp = pfDisp[i];

        // What to do with very far attraction? For now just ignore the face
        if (magSqr(fDisp) < sqr(snapDist[pointi]) && mag(fSNormal) > VSMALL)
        {
            const point pt = fc + fDisp;

            // Do we already have surface normal?
            faceToNormalBin[i] = findNormal
            (
                featureCos,
                fSNormal,
                surfaceNormals
            );

            if (faceToNormalBin[i] != -1)
            {
                // Same normal
            }
            else
            {
                // Now check if the planes go through the same edge or point

                if (surfacePoints.size() <= 1)
                {
                    surfacePoints.append(pt);
                    faceToNormalBin[i] = surfaceNormals.size();
                    surfaceNormals.append(fSNormal);
                }
                else if (surfacePoints.size() == 2)
                {
                    plane pl0(surfacePoints[0], surfaceNormals[0]);
                    plane pl1(surfacePoints[1], surfaceNormals[1]);
                    plane::ray r(pl0.planeIntersect(pl1));
                    vector featureNormal = r.dir() / mag(r.dir());

                    if (mag(fSNormal&featureNormal) >= 0.001)
                    {
                        // Definitely makes a feature point
                        surfacePoints.append(pt);
                        faceToNormalBin[i] = surfaceNormals.size();
                        surfaceNormals.append(fSNormal);
                    }
                }
                else if (surfacePoints.size() == 3)
                {
                    // Have already feature point. See if this new plane is
                    // the same point or not.
                    plane pl0(surfacePoints[0], surfaceNormals[0]);
                    plane pl1(surfacePoints[1], surfaceNormals[1]);
                    plane pl2(surfacePoints[2], surfaceNormals[2]);
                    point p012(pl0.planePlaneIntersect(pl1, pl2));

                    plane::ray r(pl0.planeIntersect(pl1));
                    vector featureNormal = r.dir() / mag(r.dir());
                    if (mag(fSNormal&featureNormal) >= 0.001)
                    {
                        plane pl3(pt, fSNormal);
                        point p013(pl0.planePlaneIntersect(pl1, pl3));

                        if (mag(p012-p013) > snapDist[pointi])
                        {
                            // Different feature point
                            surfacePoints.append(pt);
                            faceToNormalBin[i] = surfaceNormals.size();
                            surfaceNormals.append(fSNormal);
                        }
                    }
                }
            }
        }
    }


    const point& pt = pp.localPoints()[pointi];

    // Check the number of directions
    if (surfaceNormals.size() == 1)
    {
        // Normal distance to plane
        vector d =
            ((surfacePoints[0]-pt) & surfaceNormals[0])
           *surfaceNormals[0];

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointi]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointi])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
    }
    else if (surfaceNormals.size() == 2)
    {
        plane pl0(surfacePoints[0], surfaceNormals[0]);
        plane pl1(surfacePoints[1], surfaceNormals[1]);
        plane::ray r(pl0.planeIntersect(pl1));
        vector n = r.dir() / mag(r.dir());

        // Get nearest point on infinite ray
        vector d = r.refPoint()-pt;
        d -= (d&n)*n;

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointi]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointi])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
        patchConstraint.applyConstraint(surfaceNormals[1]);
    }
    else if (surfaceNormals.size() == 3)
    {
        // Calculate point from the faces.
        plane pl0(surfacePoints[0], surfaceNormals[0]);
        plane pl1(surfacePoints[1], surfaceNormals[1]);
        plane pl2(surfacePoints[2], surfaceNormals[2]);
        point cornerPt(pl0.planePlaneIntersect(pl1, pl2));
        vector d = cornerPt - pt;

        // Trim to snap distance
        if (magSqr(d) > sqr(snapDist[pointi]))
        {
            d *= Foam::sqrt(sqr(snapDist[pointi])/magSqr(d));
        }

        patchAttraction = d;

        // Store constraints
        patchConstraint.applyConstraint(surfaceNormals[0]);
        patchConstraint.applyConstraint(surfaceNormals[1]);
        patchConstraint.applyConstraint(surfaceNormals[2]);
    }
}


// Special version that calculates attraction in one go
void Foam::snappySnapDriver::featureAttractionUsingReconstruction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestDisp,

    const List<List<point>>& pointFaceSurfNormals,
    const List<List<point>>& pointFaceDisp,
    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> feStr;
    autoPtr<OBJstream> fpStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        feStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "implicitFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping implicit feature-edge direction to "
            << feStr().name() << endl;

        fpStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "implicitFeaturePoint_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping implicit feature-point direction to "
            << fpStr().name() << endl;
    }


    DynamicList<point> surfacePoints(4);
    DynamicList<vector> surfaceNormals(4);
    labelList faceToNormalBin;

    forAll(pp.localPoints(), pointi)
    {
        vector attraction = Zero;
        pointConstraint constraint;

        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            snapDist,
            nearestDisp,

            pointi,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            surfacePoints,
            surfaceNormals,
            faceToNormalBin,

            attraction,
            constraint
        );

        if
        (
            (constraint.first() > patchConstraints[pointi].first())
         || (
                (constraint.first() == patchConstraints[pointi].first())
             && (magSqr(attraction) < magSqr(patchAttraction[pointi]))
            )
        )
        {
            patchAttraction[pointi] = attraction;
            patchConstraints[pointi] = constraint;

            const point& pt = pp.localPoints()[pointi];

            if (patchConstraints[pointi].first() == 2 && feStr.valid())
            {
                feStr().write(linePointRef(pt, pt+patchAttraction[pointi]));
            }
            else if (patchConstraints[pointi].first() == 3 && fpStr.valid())
            {
                fpStr().write(linePointRef(pt, pt+patchAttraction[pointi]));
            }
        }
    }
}


void Foam::snappySnapDriver::stringFeatureEdges
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    // Snap edges to feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Walk existing edges and snap remaining ones (that are marked as
    // feature edges in rawPatchConstraints)

    // What this does is fill in any faces where not all points
    // on the face are being attracted:
    /*
           +
          / \
         /   \
      ---+    +---
         \   /
          \ /
           +
    */
    // so the top and bottom will never get attracted since the nearest
    // back from the feature edge will always be one of the left or right
    // points since the face is diamond like. So here we walk the feature edges
    // and add any non-attracted points.


    while (true)
    {
        label nChanged = 0;

        const labelListList& pointEdges = pp.pointEdges();
        forAll(pointEdges, pointi)
        {
            if (patchConstraints[pointi].first() == 2)
            {
                const point& pt = pp.localPoints()[pointi];
                const labelList& pEdges = pointEdges[pointi];
                const vector& featVec = patchConstraints[pointi].second();

                // Detect whether there are edges in both directions.
                // (direction along the feature edge that is)
                bool hasPos = false;
                bool hasNeg = false;

                forAll(pEdges, pEdgei)
                {
                    const edge& e = pp.edges()[pEdges[pEdgei]];
                    label nbrPointi = e.otherVertex(pointi);

                    if (patchConstraints[nbrPointi].first() > 1)
                    {
                        const point& nbrPt = pp.localPoints()[nbrPointi];
                        const point featPt =
                            nbrPt + patchAttraction[nbrPointi];
                        const scalar cosAngle = (featVec & (featPt-pt));

                        if (cosAngle > 0)
                        {
                            hasPos = true;
                        }
                        else
                        {
                            hasNeg = true;
                        }
                    }
                }

                if (!hasPos || !hasNeg)
                {
                    //Pout<< "**Detected feature string end at  "
                    //    << pp.localPoints()[pointi] << endl;

                    // No string. Assign best choice on either side
                    label bestPosPointi = -1;
                    scalar minPosDistSqr = GREAT;
                    label bestNegPointi = -1;
                    scalar minNegDistSqr = GREAT;

                    forAll(pEdges, pEdgei)
                    {
                        const edge& e = pp.edges()[pEdges[pEdgei]];
                        label nbrPointi = e.otherVertex(pointi);

                        if
                        (
                            patchConstraints[nbrPointi].first() <= 1
                         && rawPatchConstraints[nbrPointi].first() > 1
                        )
                        {
                            const vector& nbrFeatVec =
                                rawPatchConstraints[pointi].second();

                            if (mag(featVec&nbrFeatVec) > featureCos)
                            {
                                // nbrPointi attracted to sameish feature
                                // Note: also check on position.

                                scalar d2 = magSqr
                                (
                                    rawPatchAttraction[nbrPointi]
                                );

                                const point featPt =
                                    pp.localPoints()[nbrPointi]
                                  + rawPatchAttraction[nbrPointi];
                                const scalar cosAngle =
                                    (featVec & (featPt-pt));

                                if (cosAngle > 0)
                                {
                                    if (!hasPos && d2 < minPosDistSqr)
                                    {
                                        minPosDistSqr = d2;
                                        bestPosPointi = nbrPointi;
                                    }
                                }
                                else
                                {
                                    if (!hasNeg && d2 < minNegDistSqr)
                                    {
                                        minNegDistSqr = d2;
                                        bestNegPointi = nbrPointi;
                                    }
                                }
                            }
                        }
                    }

                    if (bestPosPointi != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt =
                        //    pp.localPoints()[bestPosPointi];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << rawPatchAttraction[bestPosPointi]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestPosPointi] =
                            0.5*rawPatchAttraction[bestPosPointi];
                        patchConstraints[bestPosPointi] =
                            rawPatchConstraints[bestPosPointi];

                        nChanged++;
                    }
                    if (bestNegPointi != -1)
                    {
                        // Use reconstructed-feature attraction. Use only
                        // part of it since not sure...
                        //const point& bestPt =
                        //    pp.localPoints()[bestNegPointi];
                        //Pout<< "**Overriding point " << bestPt
                        //    << " on reconstructed feature edge at "
                        //    << rawPatchAttraction[bestNegPointi]+bestPt
                        //    << " to attracted-to-feature-edge." << endl;
                        patchAttraction[bestNegPointi] =
                            0.5*rawPatchAttraction[bestNegPointi];
                        patchConstraints[bestNegPointi] =
                            rawPatchConstraints[bestNegPointi];

                        nChanged++;
                    }
                }
            }
        }


        reduce(nChanged, sumOp<label>());
        Info<< "Stringing feature edges : changed " << nChanged << " points"
            << endl;
        if (nChanged == 0)
        {
            break;
        }
    }
}


void Foam::snappySnapDriver::releasePointsNextToMultiPatch
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> multiPatchStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        multiPatchStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "multiPatch_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping removed constraints due to same-face"
            << " multi-patch points to "
            << multiPatchStr().name() << endl;
    }


    // 1. Mark points on multiple patches
    PackedBoolList isMultiPatchPoint(pp.size());

    forAll(pointFacePatchID, pointi)
    {
        pointIndexHit multiPatchPt = findMultiPatchPoint
        (
            pp.localPoints()[pointi],
            pointFacePatchID[pointi],
            pointFaceCentres[pointi]
        );
        isMultiPatchPoint[pointi] = multiPatchPt.hit();
    }

    // 2. Make sure multi-patch points are also attracted
    forAll(isMultiPatchPoint, pointi)
    {
        if (isMultiPatchPoint[pointi])
        {
            if
            (
                patchConstraints[pointi].first() <= 1
             && rawPatchConstraints[pointi].first() > 1
            )
            {
                patchAttraction[pointi] = rawPatchAttraction[pointi];
                patchConstraints[pointi] = rawPatchConstraints[pointi];

                //if (multiPatchStr.valid())
                //{
                //    Pout<< "Adding constraint on multiPatchPoint:"
                //        << pp.localPoints()[pointi]
                //        << " constraint:" << patchConstraints[pointi]
                //        << " attraction:" << patchAttraction[pointi]
                //        << endl;
                //}
            }
        }
    }

    // Up to here it is all parallel ok.


    // 3. Knock out any attraction on faces with multi-patch points
    label nChanged = 0;
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        label nMultiPatchPoints = 0;
        forAll(f, fp)
        {
            label pointi = f[fp];
            if
            (
                isMultiPatchPoint[pointi]
             && patchConstraints[pointi].first() > 1
            )
            {
                nMultiPatchPoints++;
            }
        }

        if (nMultiPatchPoints > 0)
        {
            forAll(f, fp)
            {
                label pointi = f[fp];
                if
                (
                   !isMultiPatchPoint[pointi]
                 && patchConstraints[pointi].first() > 1
                )
                {
                    //Pout<< "Knocking out constraint"
                    //    << " on non-multiPatchPoint:"
                    //    << pp.localPoints()[pointi] << endl;
                    patchAttraction[pointi] = Zero;
                    patchConstraints[pointi] = pointConstraint();
                    nChanged++;

                    if (multiPatchStr.valid())
                    {
                        multiPatchStr().write(pp.localPoints()[pointi]);
                    }
                }
            }
        }
    }

    reduce(nChanged, sumOp<label>());
    Info<< "Removing constraints near multi-patch points : changed "
        << nChanged << " points" << endl;
}


Foam::labelPair Foam::snappySnapDriver::findDiagonalAttraction
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchAttraction,
    const List<pointConstraint>& patchConstraints,
    const label facei
) const
{
    const face& f = pp.localFaces()[facei];
    // For now just detect any attraction. Improve this to look at
    // actual attraction position and orientation

    labelPair attractIndices(-1, -1);

    if (f.size() >= 4)
    {
        for (label startFp = 0; startFp < f.size()-2; startFp++)
        {
            label minFp = f.rcIndex(startFp);

            for
            (
                label endFp = f.fcIndex(f.fcIndex(startFp));
                endFp < f.size() && endFp != minFp;
                endFp++
            )
            {
                if
                (
                    patchConstraints[f[startFp]].first() >= 2
                 && patchConstraints[f[endFp]].first() >= 2
                )
                {
                    attractIndices = labelPair(startFp, endFp);
                    break;
                }
            }
        }
    }
    return attractIndices;
}


bool Foam::snappySnapDriver::isSplitAlignedWithFeature
(
    const scalar featureCos,
    const point& p0,
    const pointConstraint& pc0,
    const point& p1,
    const pointConstraint& pc1
) const
{
    vector d(p1-p0);
    scalar magD = mag(d);
    if (magD < VSMALL)
    {
        // Two diagonal points already colocated?
        return false;
    }
    else
    {
        d /= magD;

        // Is diagonal d aligned with at least one of the feature
        // edges?

        if (pc0.first() == 2 && mag(d & pc0.second()) > featureCos)
        {
            return true;
        }
        else if (pc1.first() == 2 && mag(d & pc1.second()) > featureCos)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


// Is situation very concave
bool Foam::snappySnapDriver::isConcave
(
    const point& c0,
    const vector& area0,
    const point& c1,
    const vector& area1,
    const scalar concaveCos
) const
{
    vector n0 = area0;
    scalar magN0 = mag(n0);
    if (magN0 < VSMALL)
    {
        // Zero area face. What to return? For now disable splitting.
        return true;
    }
    n0 /= magN0;

    // Distance from c1 to plane of face0
    scalar d = (c1-c0)&n0;

    if (d <= 0)
    {
        // Convex (face1 centre on 'inside' of face0)
        return false;
    }
    else
    {
        // Is a bit or very concave?
        vector n1 = area1;
        scalar magN1 = mag(n1);
        if (magN1 < VSMALL)
        {
            // Zero area face. See above
            return true;
        }
        n1 /= magN1;

        if ((n0&n1) < concaveCos)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}


Foam::labelPair Foam::snappySnapDriver::findDiagonalAttraction
(
    const scalar featureCos,
    const scalar concaveCos,
    const scalar minAreaRatio,
    const indirectPrimitivePatch& pp,
    const vectorField& patchAttr,
    const List<pointConstraint>& patchConstraints,
    const vectorField& nearestAttr,
    const vectorField& nearestNormal,
    const label facei,

    DynamicField<point>& points0,
    DynamicField<point>& points1
) const
{
    const face& localF = pp.localFaces()[facei];

    labelPair attractIndices(-1, -1);

    if (localF.size() >= 4)
    {
        const pointField& localPts = pp.localPoints();

        //// Estimate cell centre taking patchAttraction into account
        //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        //// (is this necessary?)
        //const polyMesh& mesh = meshRefiner_.mesh();
        //label meshFacei = pp.addressing()[facei];
        //const face& meshF = mesh.faces()[meshFacei];
        //label celli = mesh.faceOwner()[meshFacei];
        //const labelList& cPoints = mesh.cellPoints(celli);
        //
        //point cc(mesh.points()[meshF[0]]);
        //for (label i = 1; i < meshF.size(); i++)
        //{
        //    cc += mesh.points()[meshF[i]]+patchAttr[localF[i]];
        //}
        //forAll(cPoints, i)
        //{
        //    label pointi = cPoints[i];
        //    if (findIndex(meshF, pointi) == -1)
        //    {
        //        cc += mesh.points()[pointi];
        //    }
        //}
        //cc /= cPoints.size();
        ////const point& cc = mesh.cellCentres()[celli];
        //
        //const scalar vol = pyrVol(pp, patchAttr, localF, cc);
        //const scalar area = localF.mag(localPts);



        // Try all diagonal cuts
        // ~~~~~~~~~~~~~~~~~~~~~

        face f0(3);
        face f1(3);

        for (label startFp = 0; startFp < localF.size()-2; startFp++)
        {
            label minFp = localF.rcIndex(startFp);

            for
            (
                label endFp = localF.fcIndex(localF.fcIndex(startFp));
                endFp < localF.size() && endFp != minFp;
                endFp++
            )
            {
                label startPti = localF[startFp];
                label endPti = localF[endFp];

                const pointConstraint& startPc = patchConstraints[startPti];
                const pointConstraint& endPc = patchConstraints[endPti];

                if (startPc.first() >= 2 && endPc.first() >= 2)
                {
                    if (startPc.first() == 2 || endPc.first() == 2)
                    {
                        // Check if
                        // - sameish feature edge normal
                        // - diagonal aligned with feature edge normal
                        point start = localPts[startPti]+patchAttr[startPti];
                        point end = localPts[endPti]+patchAttr[endPti];

                        if
                        (
                           !isSplitAlignedWithFeature
                            (
                                featureCos,
                                start,
                                startPc,
                                end,
                                endPc
                            )
                        )
                        {
                            // Attract to different features. No need to
                            // introduce split
                            continue;
                        }
                    }



                    // Form two faces
                    // ~~~~~~~~~~~~~~
                    // Predict position of faces. End points of the faces
                    // attract to the feature
                    // and all the other points just attract to the nearest

                    // face0

                    f0.setSize(endFp-startFp+1);
                    label i0 = 0;
                    for (label fp = startFp; fp <= endFp; fp++)
                    {
                        f0[i0++] = localF[fp];
                    }

                    // Get compact face and points
                    const face compact0(identity(f0.size()));
                    points0.clear();
                    points0.append(localPts[f0[0]] + patchAttr[f0[0]]);
                    for (label fp=1; fp < f0.size()-1; fp++)
                    {
                        label pi = f0[fp];
                        points0.append(localPts[pi] + nearestAttr[pi]);
                    }
                    points0.append
                    (
                        localPts[f0.last()] + patchAttr[f0.last()]
                    );


                    // face1

                    f1.setSize(localF.size()+2-f0.size());
                    label i1 = 0;
                    for
                    (
                        label fp = endFp;
                        fp != startFp;
                        fp = localF.fcIndex(fp)
                    )
                    {
                        f1[i1++] = localF[fp];
                    }
                    f1[i1++] = localF[startFp];


                    // Get compact face and points
                    const face compact1(identity(f1.size()));
                    points1.clear();
                    points1.append(localPts[f1[0]] + patchAttr[f1[0]]);
                    for (label fp=1; fp < f1.size()-1; fp++)
                    {
                        label pi = f1[fp];
                        points1.append(localPts[pi] + nearestAttr[pi]);
                    }
                    points1.append
                    (
                        localPts[f1.last()] + patchAttr[f1.last()]
                    );



                    // Avoid splitting concave faces
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    if
                    (
                        isConcave
                        (
                            compact0.centre(points0),
                            compact0.areaNormal(points0),
                            compact1.centre(points1),
                            compact1.areaNormal(points1),
                            concaveCos
                        )
                    )
                    {
                        //// Dump to obj
                        //Pout<< "# f0:" << f0 << endl;
                        //forAll(p, i)
                        //{
                        //    meshTools::writeOBJ(Pout, points0[i]);
                        //}
                        //Pout<< "# f1:" << f1 << endl;
                        //forAll(p, i)
                        //{
                        //    meshTools::writeOBJ(Pout, points1[i]);
                        //}
                    }
                    else
                    {
                        // Existing areas
                        const scalar area0 = f0.mag(localPts);
                        const scalar area1 = f1.mag(localPts);

                        if
                        (
                            area0/area1 >= minAreaRatio
                         && area1/area0 >= minAreaRatio
                        )
                        {
                            attractIndices = labelPair(startFp, endFp);
                        }
                    }
                }
            }
        }
    }
    return attractIndices;
}


void Foam::snappySnapDriver::splitDiagonals
(
    const scalar featureCos,
    const scalar concaveCos,
    const scalar minAreaRatio,

    const indirectPrimitivePatch& pp,
    const vectorField& nearestAttraction,
    const vectorField& nearestNormal,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints,
    DynamicList<label>& splitFaces,
    DynamicList<labelPair>& splits
) const
{
    const labelList& bFaces = pp.addressing();

    splitFaces.clear();
    splitFaces.setCapacity(bFaces.size());
    splits.clear();
    splits.setCapacity(bFaces.size());


    // Work arrays for storing points of face
    DynamicField<point> facePoints0;
    DynamicField<point> facePoints1;

    forAll(bFaces, facei)
    {
        const labelPair split
        (
            findDiagonalAttraction
            (
                featureCos,
                concaveCos,
                minAreaRatio,

                pp,
                patchAttraction,
                patchConstraints,

                nearestAttraction,
                nearestNormal,
                facei,

                facePoints0,
                facePoints1
            )
        );

        if (split != labelPair(-1, -1))
        {
            splitFaces.append(bFaces[facei]);
            splits.append(split);

            const face& f = pp.localFaces()[facei];

            // Knock out other attractions on face
            forAll(f, fp)
            {
                // Knock out any other constraints
                if
                (
                    fp != split[0]
                 && fp != split[1]
                 && patchConstraints[f[fp]].first() >= 2
                )
                {
                    patchConstraints[f[fp]] = pointConstraint
                    (
                        Tuple2<label, vector>
                        (
                            1,
                            nearestNormal[f[fp]]
                        )
                    );
                    patchAttraction[f[fp]] = nearestAttraction[f[fp]];
                }
            }
        }
    }
}


void Foam::snappySnapDriver::avoidDiagonalAttraction
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        labelPair diag = findDiagonalAttraction
        (
            pp,
            patchAttraction,
            patchConstraints,
            facei
        );

        if (diag[0] != -1 && diag[1] != -1)
        {
            // Found two diagonal points that being attracted.
            // For now just attract my one to the average of those.
            const label i0 = f[diag[0]];
            const point pt0 =
                pp.localPoints()[i0]+patchAttraction[i0];
            const label i1 = f[diag[1]];
            const point pt1 =
                pp.localPoints()[i1]+patchAttraction[i1];
            const point mid = 0.5*(pt0+pt1);

            const scalar cosAngle = mag
            (
                patchConstraints[i0].second()
              & patchConstraints[i1].second()
            );

            //Pout<< "Found diagonal attraction at indices:"
            //    << diag[0]
            //    << " and " << diag[1]
            //    << " with cosAngle:" << cosAngle
            //    << " mid:" << mid << endl;

            if (cosAngle > featureCos)
            {
                // The two feature edges are roughly in the same direction.
                // Add the nearest of the other points in the face as
                // attractor
                label minFp = -1;
                scalar minDistSqr = GREAT;
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    if (patchConstraints[pointi].first() <= 1)
                    {
                        const point& pt = pp.localPoints()[pointi];
                        scalar distSqr = magSqr(mid-pt);
                        if (distSqr < minDistSqr)
                        {
                            distSqr = minDistSqr;
                            minFp = fp;
                        }
                    }
                }
                if (minFp != -1)
                {
                    label minPointi = f[minFp];
                    patchAttraction[minPointi] =
                        mid-pp.localPoints()[minPointi];
                    patchConstraints[minPointi] = patchConstraints[f[diag[0]]];
                }
            }
            else
            {
                //Pout<< "Diagonal attractors at" << nl
                //    << "    pt0:" << pt0
                //    << "    constraint:"
                //    << patchConstraints[i0].second() << nl
                //    << "    pt1:" << pt1
                //    << "    constraint:"
                //    << patchConstraints[i1].second() << nl
                //    << "    make too large an angle:"
                //    <<  mag
                //        (
                //            patchConstraints[i0].second()
                //          & patchConstraints[i1].second()
                //        )
                //    << endl;
            }
        }
    }
}


Foam::Tuple2<Foam::label, Foam::pointIndexHit>
Foam::snappySnapDriver::findNearFeatureEdge
(
    const bool isRegionEdge,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const label pointi,
    const point& estimatedPt,

    List<List<DynamicList<point>>>& edgeAttractors,
    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    labelList nearEdgeFeat;
    List<pointIndexHit> nearEdgeInfo;
    vectorField nearNormal;

    if (isRegionEdge)
    {
        features.findNearestRegionEdge
        (
            pointField(1, estimatedPt),
            scalarField(1, sqr(snapDist[pointi])),
            nearEdgeFeat,
            nearEdgeInfo,
            nearNormal
        );
    }
    else
    {
        features.findNearestEdge
        (
            pointField(1, estimatedPt),
            scalarField(1, sqr(snapDist[pointi])),
            nearEdgeFeat,
            nearEdgeInfo,
            nearNormal
        );
    }

    const pointIndexHit& nearInfo = nearEdgeInfo[0];
    label feati = nearEdgeFeat[0];

    if (nearInfo.hit())
    {
        // So we have a point on the feature edge. Use this
        // instead of our estimate from planes.
        edgeAttractors[feati][nearInfo.index()].append
        (
            nearInfo.hitPoint()
        );
        pointConstraint c(Tuple2<label, vector>(2, nearNormal[0]));
        edgeConstraints[feati][nearInfo.index()].append(c);

        // Store for later use
        patchAttraction[pointi] = nearInfo.hitPoint()-pp.localPoints()[pointi];
        patchConstraints[pointi] = c;
    }
    return Tuple2<label, pointIndexHit>(feati, nearInfo);
}


Foam::Tuple2<Foam::label, Foam::pointIndexHit>
Foam::snappySnapDriver::findNearFeaturePoint
(
    const bool isRegionPoint,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const label pointi,
    const point& estimatedPt,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint>>& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point>>>& edgeAttractors,
    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    // Search for for featurePoints only! This ignores any non-feature points.

    labelList nearFeat;
    List<pointIndexHit> nearInfo;
    features.findNearestPoint
    (
        pointField(1, estimatedPt),
        scalarField(1, sqr(snapDist[pointi])),
        nearFeat,
        nearInfo
    );

    label feati = nearFeat[0];

    if (feati != -1)
    {
        const point& pt = pp.localPoints()[pointi];

        label featPointi = nearInfo[0].index();
        const point& featPt = nearInfo[0].hitPoint();
        scalar distSqr = magSqr(featPt-pt);

        // Check if already attracted
        label oldPointi = pointAttractor[feati][featPointi];

        if (oldPointi != -1)
        {
            // Check distance
            if (distSqr >= magSqr(featPt-pp.localPoints()[oldPointi]))
            {
                // oldPointi nearest. Keep.
                feati = -1;
                featPointi = -1;
            }
            else
            {
                // Current pointi nearer.
                pointAttractor[feati][featPointi] = pointi;
                pointConstraints[feati][featPointi].first() = 3;
                pointConstraints[feati][featPointi].second() = Zero;

                // Store for later use
                patchAttraction[pointi] = featPt-pt;
                patchConstraints[pointi] = pointConstraints[feati][featPointi];

                // Reset oldPointi to nearest on feature edge
                patchAttraction[oldPointi] = Zero;
                patchConstraints[oldPointi] = pointConstraint();

                findNearFeatureEdge
                (
                    isRegionPoint,      // search region edges only

                    pp,
                    snapDist,
                    oldPointi,
                    pp.localPoints()[oldPointi],

                    edgeAttractors,
                    edgeConstraints,
                    patchAttraction,
                    patchConstraints
                );
            }
        }
        else
        {
            // Current pointi nearer.
            pointAttractor[feati][featPointi] = pointi;
            pointConstraints[feati][featPointi].first() = 3;
            pointConstraints[feati][featPointi].second() = Zero;

            // Store for later use
            patchAttraction[pointi] = featPt-pt;
            patchConstraints[pointi] = pointConstraints[feati][featPointi];
        }
    }

    return Tuple2<label, pointIndexHit>(feati, nearInfo[0]);
}


// Determines for every pp point - that is on multiple faces that form
// a feature - the nearest feature edge/point.
void Foam::snappySnapDriver::determineFeatures
(
    const label iter,
    const scalar featureCos,
    const bool multiRegionFeatureSnap,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestDisp,

    const List<List<point>>& pointFaceSurfNormals,
    const List<List<point>>& pointFaceDisp,
    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint>>& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point>>>& edgeAttractors,
    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> featureEdgeStr;
    autoPtr<OBJstream> missedEdgeStr;
    autoPtr<OBJstream> featurePointStr;
    autoPtr<OBJstream> missedMP0Str;
    autoPtr<OBJstream> missedMP1Str;

    if (debug&meshRefinement::ATTRACTION)
    {
        featureEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "featureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-edge sampling to "
            << featureEdgeStr().name() << endl;

        missedEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-edges that are too far away to "
            << missedEdgeStr().name() << endl;

        featurePointStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "featurePoint_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping feature-point sampling to "
            << featurePointStr().name() << endl;

        missedMP0Str.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdgeFromMPEdge_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping region-edges that are too far away to "
            << missedMP0Str().name() << endl;

        missedMP1Str.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "missedFeatureEdgeFromMPPoint_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping region-points that are too far away to "
            << missedMP1Str().name() << endl;
    }


    DynamicList<point> surfacePoints(4);
    DynamicList<vector> surfaceNormals(4);
    labelList faceToNormalBin;

    forAll(pp.localPoints(), pointi)
    {
        const point& pt = pp.localPoints()[pointi];


        // Determine the geometric planes the point is (approximately) on.
        // This is returned as a
        // - attraction vector
        // - and a constraint
        //   (1: attract to surface, constraint is normal of plane
        //    2: attract to feature line, constraint is feature line direction
        //    3: attract to feature point, constraint is zero)

        vector attraction = Zero;
        pointConstraint constraint;

        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            snapDist,
            nearestDisp,

            pointi,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            surfacePoints,
            surfaceNormals,
            faceToNormalBin,

            attraction,
            constraint
        );

        // Now combine the reconstruction with the current state of the
        // point. The logic is quite complicated:
        // - the new constraint (from reconstruction) will only win if
        //   - the constraint is higher (feature-point wins from feature-edge
        //     etc.)
        //   - or the constraint is the same but the attraction distance is less
        //
        // - then this will be combined with explicit searching on the
        //   features and optionally the analysis of the patches using the
        //   point. This analysis can do three thing:
        //      - the point is not on multiple patches
        //      - the point is on multiple patches but these are also
        //        different planes (so the region feature is also a geometric
        //        feature)
        //      - the point is on multiple patches some of which are on
        //        the same plane. This is the problem one - do we assume it is
        //        an additional constraint (feat edge upgraded to region point,
        //        see below)?
        //
        //      Reconstruction  MultiRegionFeatureSnap          Attraction
        //      -------         ----------------------          -----------
        //      surface         false                           surface
        //      surface         true                            region edge
        //      feat edge       false                           feat edge
        //      feat edge       true and no planar regions      feat edge
        //      feat edge       true and yes planar regions     region point
        //      feat point      false                           feat point
        //      feat point      true                            region point


        if
        (
            (constraint.first() > patchConstraints[pointi].first())
         || (
                (constraint.first() == patchConstraints[pointi].first())
             && (magSqr(attraction) < magSqr(patchAttraction[pointi]))
            )
        )
        {
            patchAttraction[pointi] = attraction;
            patchConstraints[pointi] = constraint;

            // Check the number of directions
            if (patchConstraints[pointi].first() == 1)
            {
                // Flat surface. Check for different patchIDs
                if (multiRegionFeatureSnap)
                {
                    const point estimatedPt(pt + nearestDisp[pointi]);
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointi],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );

                    if (multiPatchPt.hit())
                    {
                        // Behave like when having two surface normals so
                        // attract to nearest feature edge (with a guess for
                        // the multipatch point as starting point)
                        Tuple2<label, pointIndexHit> nearInfo =
                        findNearFeatureEdge
                        (
                            true,                       // isRegionEdge
                            pp,
                            snapDist,
                            pointi,
                            multiPatchPt.hitPoint(),    // estimatedPt

                            edgeAttractors,
                            edgeConstraints,

                            patchAttraction,
                            patchConstraints
                        );

                        const pointIndexHit& info = nearInfo.second();
                        if (info.hit())
                        {
                            // Dump
                            if (featureEdgeStr.valid())
                            {
                                featureEdgeStr().write
                                (
                                    linePointRef(pt, info.hitPoint())
                                );
                            }
                        }
                        else
                        {
                            if (missedEdgeStr.valid())
                            {
                                missedEdgeStr().write
                                (
                                    linePointRef(pt, multiPatchPt.hitPoint())
                                );
                            }
                        }
                    }
                }
            }
            else if (patchConstraints[pointi].first() == 2)
            {
                // Mark point on the nearest feature edge. Note that we
                // only search within the surrounding since the plane
                // reconstruction might find a feature where there isn't one.
                const point estimatedPt(pt + patchAttraction[pointi]);

                Tuple2<label, pointIndexHit> nearInfo(-1, pointIndexHit());

                // Geometric feature edge. Check for different patchIDs
                bool hasSnapped = false;
                if (multiRegionFeatureSnap)
                {
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointi],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );
                    if (multiPatchPt.hit())
                    {
                        if (multiPatchPt.index() == 0)
                        {
                            // Region edge is also a geometric feature edge
                            nearInfo = findNearFeatureEdge
                            (
                                true,               // isRegionEdge
                                pp,
                                snapDist,
                                pointi,
                                estimatedPt,

                                edgeAttractors,
                                edgeConstraints,

                                patchAttraction,
                                patchConstraints
                            );
                            hasSnapped = true;

                            // Debug: dump missed feature point
                            if
                            (
                                missedMP0Str.valid()
                            && !nearInfo.second().hit()
                            )
                            {
                                missedMP0Str().write
                                (
                                    linePointRef(pt, estimatedPt)
                                );
                            }
                        }
                        else
                        {
                            // One of planes of feature contains multiple
                            // regions. We assume (contentious!) that the
                            // separation between
                            // the regions is not aligned with the geometric
                            // feature so is an additional constraint on the
                            // point -> is region-feature-point.
                            nearInfo = findNearFeaturePoint
                            (
                                true,           // isRegionPoint
                                pp,
                                snapDist,
                                pointi,
                                estimatedPt,

                                // Feature-point to pp point
                                pointAttractor,
                                pointConstraints,
                                // Feature-edge to pp point
                                edgeAttractors,
                                edgeConstraints,
                                // pp point to nearest feature
                                patchAttraction,
                                patchConstraints
                            );
                            hasSnapped = true;

                            // More contentious: if we don't find
                            // a near feature point we will never find the
                            // attraction to a feature edge either since
                            // the edgeAttractors/edgeConstraints do not get
                            // filled and we're using reverse attraction
                            // Note that we're in multiRegionFeatureSnap which
                            // where findMultiPatchPoint can decide the
                            // wrong thing. So: if failed finding a near
                            // feature point try for a feature edge
                            if (!nearInfo.second().hit())
                            {
                                nearInfo = findNearFeatureEdge
                                (
                                    true,           // isRegionEdge
                                    pp,
                                    snapDist,
                                    pointi,
                                    estimatedPt,

                                    // Feature-edge to pp point
                                    edgeAttractors,
                                    edgeConstraints,
                                    // pp point to nearest feature
                                    patchAttraction,
                                    patchConstraints
                                );
                            }

                            // Debug: dump missed feature point
                            if
                            (
                                missedMP1Str.valid()
                            && !nearInfo.second().hit()
                            )
                            {
                                missedMP1Str().write
                                (
                                    linePointRef(pt, estimatedPt)
                                );
                            }
                        }
                    }
                }

                if (!hasSnapped)
                {
                    // Determine nearest point on feature edge. Store
                    // constraint
                    // (calculated from feature edge, alternative would be to
                    //  use constraint calculated from both surfaceNormals)
                    nearInfo = findNearFeatureEdge
                    (
                        false,      // isRegionPoint
                        pp,
                        snapDist,
                        pointi,
                        estimatedPt,

                        edgeAttractors,
                        edgeConstraints,

                        patchAttraction,
                        patchConstraints
                    );
                    hasSnapped = true;
                }

                // Dump to obj
                const pointIndexHit& info = nearInfo.second();
                if (info.hit())
                {
                    if
                    (
                        patchConstraints[pointi].first() == 3
                     && featurePointStr.valid()
                    )
                    {
                        featurePointStr().write
                        (
                            linePointRef(pt, info.hitPoint())
                        );
                    }
                    else if
                    (
                        patchConstraints[pointi].first() == 2
                     && featureEdgeStr.valid()
                    )
                    {
                        featureEdgeStr().write
                        (
                            linePointRef(pt, info.hitPoint())
                        );
                    }
                }
                else
                {
                    if (missedEdgeStr.valid())
                    {
                        missedEdgeStr().write
                        (
                            linePointRef(pt, estimatedPt)
                        );
                    }
                }
            }
            else if (patchConstraints[pointi].first() == 3)
            {
                // Mark point on the nearest feature point.
                const point estimatedPt(pt + patchAttraction[pointi]);

                Tuple2<label, pointIndexHit> nearInfo(-1, pointIndexHit());

                if (multiRegionFeatureSnap)
                {
                    pointIndexHit multiPatchPt
                    (
                        findMultiPatchPoint
                        (
                            estimatedPt,
                            pointFacePatchID[pointi],
                            surfaceNormals,
                            faceToNormalBin
                        )
                    );
                    if (multiPatchPt.hit())
                    {
                        // Multiple regions
                        nearInfo = findNearFeaturePoint
                        (
                            true,           // isRegionPoint
                            pp,
                            snapDist,
                            pointi,
                            estimatedPt,

                            // Feature-point to pp point
                            pointAttractor,
                            pointConstraints,
                            // Feature-edge to pp point
                            edgeAttractors,
                            edgeConstraints,
                            // pp point to nearest feature
                            patchAttraction,
                            patchConstraints
                        );
                    }
                    else
                    {
                        nearInfo = findNearFeaturePoint
                        (
                            false,              // isRegionPoint
                            pp,
                            snapDist,
                            pointi,
                            estimatedPt,

                            // Feature-point to pp point
                            pointAttractor,
                            pointConstraints,
                            // Feature-edge to pp point
                            edgeAttractors,
                            edgeConstraints,
                            // pp point to nearest feature
                            patchAttraction,
                            patchConstraints
                        );
                    }
                }
                else
                {
                    // No multi-patch snapping
                    nearInfo = findNearFeaturePoint
                    (
                        false,              // isRegionPoint
                        pp,
                        snapDist,
                        pointi,
                        estimatedPt,

                        // Feature-point to pp point
                        pointAttractor,
                        pointConstraints,
                        // Feature-edge to pp point
                        edgeAttractors,
                        edgeConstraints,
                        // pp point to nearest feature
                        patchAttraction,
                        patchConstraints
                    );
                }

                const pointIndexHit& info = nearInfo.second();
                if (info.hit() && featurePointStr.valid())
                {
                    featurePointStr().write
                    (
                        linePointRef(pt, info.hitPoint())
                    );
                }
            }
        }
    }
}


// Baffle handling
// ~~~~~~~~~~~~~~~
// Override pointAttractor, edgeAttractor, patchAttration etc. to
// implement 'baffle' handling.
// Baffle: the mesh pp point originates from a loose standing
// baffle.
// Sampling the surface with the surrounding face-centres only picks up
// a single triangle normal so above determineFeatures will not have
// detected anything. So explicitly pick up feature edges on the pp
// (after duplicating points & smoothing so will already have been
// expanded) and match these to the features.
void Foam::snappySnapDriver::determineBaffleFeatures
(
    const label iter,
    const bool baffleFeaturePoints,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    // Feature-point to pp point
    List<labelList>& pointAttractor,
    List<List<pointConstraint>>& pointConstraints,
    // Feature-edge to pp point
    List<List<DynamicList<point>>>& edgeAttractors,
    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementFeatures& features = meshRefiner_.features();

    // Calculate edge-faces
    List<List<point>> edgeFaceNormals(pp.nEdges());

    // Fill local data
    forAll(pp.edgeFaces(), edgei)
    {
        const labelList& eFaces = pp.edgeFaces()[edgei];
        List<point>& eFc = edgeFaceNormals[edgei];
        eFc.setSize(eFaces.size());
        forAll(eFaces, i)
        {
            label facei = eFaces[i];
            eFc[i] = pp.faceNormals()[facei];
        }
    }

    {
        // Precalculate mesh edges for pp.edges.
        const labelList meshEdges
        (
            pp.meshEdges(mesh.edges(), mesh.pointEdges())
        );
        syncTools::syncEdgeList
        (
            mesh,
            meshEdges,
            edgeFaceNormals,
            listPlusEqOp<point>(),
            List<point>(),
            mapDistribute::transform()
        );
    }

    // Detect baffle edges. Assume initial mesh will have 0,90 or 180
    // (baffle) degree angles so smoothing should make 0,90
    // to be less than 90. Choose reasonable value
    const scalar baffleFeatureCos = Foam::cos(degToRad(110));


    autoPtr<OBJstream> baffleEdgeStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        baffleEdgeStr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "baffleEdge_" + name(iter) + ".obj"
            )
        );
        Info<< nl << "Dumping baffle-edges to "
            << baffleEdgeStr().name() << endl;
    }


    // Is edge on baffle
    PackedBoolList isBaffleEdge(pp.nEdges());
    label nBaffleEdges = 0;
    // Is point on
    //  0 : baffle-edge (0)
    //  1 : baffle-feature-point (1)
    // -1 : rest
    labelList pointStatus(pp.nPoints(), -1);

    forAll(edgeFaceNormals, edgei)
    {
        const List<point>& efn = edgeFaceNormals[edgei];

        if (efn.size() == 2 && (efn[0]&efn[1]) < baffleFeatureCos)
        {
            isBaffleEdge[edgei] = true;
            nBaffleEdges++;
            const edge& e = pp.edges()[edgei];
            pointStatus[e[0]] = 0;
            pointStatus[e[1]] = 0;

            if (baffleEdgeStr.valid())
            {
                const point& p0 = pp.localPoints()[e[0]];
                const point& p1 = pp.localPoints()[e[1]];
                baffleEdgeStr().write(linePointRef(p0, p1));
            }
        }
    }

    reduce(nBaffleEdges, sumOp<label>());

    Info<< "Detected " << nBaffleEdges
        << " baffle edges out of "
        << returnReduce(pp.nEdges(), sumOp<label>())
        << " edges." << endl;


    //- Baffle edges will be too ragged to sensibly determine feature points
    //forAll(pp.pointEdges(), pointi)
    //{
    //    if
    //    (
    //        isFeaturePoint
    //        (
    //            featureCos,
    //            pp,
    //            isBaffleEdge,
    //            pointi
    //        )
    //    )
    //    {
    //        //Pout<< "Detected feature point:" << pp.localPoints()[pointi]
    //        //    << endl;
    //        //-TEMPORARILY DISABLED:
    //        //pointStatus[pointi] = 1;
    //    }
    //}


    label nBafflePoints = 0;
    forAll(pointStatus, pointi)
    {
        if (pointStatus[pointi] != -1)
        {
            nBafflePoints++;
        }
    }
    reduce(nBafflePoints, sumOp<label>());


    label nPointAttract = 0;
    label nEdgeAttract = 0;

    forAll(pointStatus, pointi)
    {
        const point& pt = pp.localPoints()[pointi];

        if (pointStatus[pointi] == 0)   // baffle edge
        {
            // 1: attract to near feature edge first

            Tuple2<label, pointIndexHit> nearInfo = findNearFeatureEdge
            (
                false,          // isRegionPoint?
                pp,
                snapDist,
                pointi,
                pt,

                edgeAttractors,
                edgeConstraints,
                patchAttraction,
                patchConstraints
            );


            //- MEJ:
            // 2: optionally override with nearest feature point.
            //    On baffles we don't have enough normals to construct a feature
            //    point so assume all feature edges are close to feature points
            if (nearInfo.second().hit())
            {
                nEdgeAttract++;

                if (baffleFeaturePoints)
                {
                    nearInfo = findNearFeaturePoint
                    (
                        false,          // isRegionPoint,

                        pp,
                        snapDist,
                        pointi,
                        pt,             // estimatedPt,

                        // Feature-point to pp point
                        pointAttractor,
                        pointConstraints,
                        // Feature-edge to pp point
                        edgeAttractors,
                        edgeConstraints,
                        // pp point to nearest feature
                        patchAttraction,
                        patchConstraints
                    );

                    if (nearInfo.first() != -1)
                    {
                        nEdgeAttract--;
                        nPointAttract++;
                    }
                }
            }
        }
        else if (pointStatus[pointi] == 1)   // baffle point
        {
            labelList nearFeat;
            List<pointIndexHit> nearInfo;
            features.findNearestPoint
            (
                pointField(1, pt),
                scalarField(1, sqr(snapDist[pointi])),
                nearFeat,
                nearInfo
            );

            label feati = nearFeat[0];

            if (feati != -1)
            {
                nPointAttract++;

                label featPointi = nearInfo[0].index();
                const point& featPt = nearInfo[0].hitPoint();
                scalar distSqr = magSqr(featPt-pt);

                // Check if already attracted
                label oldPointi = pointAttractor[feati][featPointi];

                if
                (
                    oldPointi == -1
                 || (
                        distSqr
                      < magSqr(featPt-pp.localPoints()[oldPointi])
                    )
                )
                {
                    pointAttractor[feati][featPointi] = pointi;
                    pointConstraints[feati][featPointi].first() = 3;
                    pointConstraints[feati][featPointi].second() = Zero;

                    // Store for later use
                    patchAttraction[pointi] = featPt-pt;
                    patchConstraints[pointi] =
                        pointConstraints[feati][featPointi];

                    if (oldPointi != -1)
                    {
                        // The current point is closer so wins. Reset
                        // the old point to attract to nearest edge
                        // instead.
                        findNearFeatureEdge
                        (
                            false,              // isRegionPoint
                            pp,
                            snapDist,
                            oldPointi,
                            pp.localPoints()[oldPointi],

                            edgeAttractors,
                            edgeConstraints,
                            patchAttraction,
                            patchConstraints
                        );
                    }
                }
                else
                {
                    // Make it fall through to check below
                    feati = -1;
                }
            }

            // Not found a feature point or another point is already
            // closer to that feature
            if (feati == -1)
            {
                //Pout<< "*** Falling back to finding nearest feature"
                //    << " edge"
                //    << " for baffle-feature-point " << pt
                //    << endl;

                Tuple2<label, pointIndexHit> nearInfo = findNearFeatureEdge
                (
                    false,                  // isRegionPoint
                    pp,
                    snapDist,
                    pointi,
                    pt,                     // starting point

                    edgeAttractors,
                    edgeConstraints,
                    patchAttraction,
                    patchConstraints
                );

                if (nearInfo.first() != -1)
                {
                    nEdgeAttract++;
                }
            }
        }
    }

    reduce(
        std::tie(nPointAttract, nEdgeAttract),
        UniformParallelOp<sumOp<label>, 2>{}
    );

    Info<< "Baffle points     : " << nBafflePoints
        << " of which attracted to :" << nl
        << "    feature point : " << nPointAttract << nl
        << "    feature edge  : " << nEdgeAttract << nl
        << "    rest          : " << nBafflePoints-nPointAttract-nEdgeAttract
        << nl
        << endl;
}


void Foam::snappySnapDriver::reverseAttractMeshPoints
(
    const label iter,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,

    // Feature-point to pp point
    const List<labelList>& pointAttractor,
    const List<List<pointConstraint>>& pointConstraints,
    // Feature-edge to pp point
    const List<List<DynamicList<point>>>& edgeAttractors,
    const List<List<DynamicList<pointConstraint>>>& edgeConstraints,

    const vectorField& rawPatchAttraction,
    const List<pointConstraint>& rawPatchConstraints,

    // pp point to nearest feature
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();

    // Find nearest mesh point to feature edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Reverse lookup : go through all edgeAttractors and find the
    // nearest point on pp

    // Get search domain and extend it a bit
    treeBoundBox bb(pp.localPoints());
    {
        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        bb = bb.extend(rndGen, 1e-4);
        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    }

    // Collect candidate points for attraction
    DynamicList<label> attractPoints(pp.nPoints());
    {
        const fvMesh& mesh = meshRefiner_.mesh();

        boolList isFeatureEdgeOrPoint(pp.nPoints(), false);
        label nFeats = 0;
        forAll(rawPatchConstraints, pointi)
        {
            if (rawPatchConstraints[pointi].first() >= 2)
            {
                isFeatureEdgeOrPoint[pointi] = true;
                nFeats++;
            }
        }

        Info<< "Initially selected " << returnReduce(nFeats, sumOp<label>())
            << " mesh points out of "
            << returnReduce(pp.nPoints(), sumOp<label>())
            << " for reverse attraction." << endl;

        // Make sure is synchronised (note: check if constraint is already
        // synced in which case this is not needed here)
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            isFeatureEdgeOrPoint,
            orEqOp<bool>(),         // combine op
            false
        );

        for (label nGrow = 0; nGrow < 1; nGrow++)
        {
            boolList newIsFeatureEdgeOrPoint(isFeatureEdgeOrPoint);

            forAll(pp.localFaces(), facei)
            {
                const face& f = pp.localFaces()[facei];

                forAll(f, fp)
                {
                    if (isFeatureEdgeOrPoint[f[fp]])
                    {
                        // Mark all points on face
                        forAll(f, fp)
                        {
                            newIsFeatureEdgeOrPoint[f[fp]] = true;
                        }
                        break;
                    }
                }
            }

            isFeatureEdgeOrPoint = newIsFeatureEdgeOrPoint;

            syncTools::syncPointList
            (
                mesh,
                pp.meshPoints(),
                isFeatureEdgeOrPoint,
                orEqOp<bool>(),         // combine op
                false
            );
        }


        // Collect attractPoints
        forAll(isFeatureEdgeOrPoint, pointi)
        {
            if (isFeatureEdgeOrPoint[pointi])
            {
                attractPoints.append(pointi);
            }
        }

        Info<< "Selected "
            << returnReduce(attractPoints.size(), sumOp<label>())
            << " mesh points out of "
            << returnReduce(pp.nPoints(), sumOp<label>())
            << " for reverse attraction." << endl;
    }


    indexedOctree<treeDataPoint> ppTree
    (
        treeDataPoint(pp.localPoints(), attractPoints),
        bb,                             // overall search domain
        8,                              // maxLevel
        10,                             // leafsize
        3.0                             // duplicity
    );

    // Per mesh point the point on nearest feature edge.
    patchAttraction.setSize(pp.nPoints());
    patchAttraction = Zero;
    patchConstraints.setSize(pp.nPoints());
    patchConstraints = pointConstraint();

    forAll(edgeAttractors, feati)
    {
        const List<DynamicList<point>>& edgeAttr = edgeAttractors[feati];
        const List<DynamicList<pointConstraint>>& edgeConstr =
            edgeConstraints[feati];

        forAll(edgeAttr, featEdgei)
        {
            const DynamicList<point>& attr = edgeAttr[featEdgei];
            forAll(attr, i)
            {
                // Find nearest pp point
                const point& featPt = attr[i];
                pointIndexHit nearInfo = ppTree.findNearest
                (
                    featPt,
                    sqr(GREAT)
                );

                if (nearInfo.hit())
                {
                    label pointi =
                        ppTree.shapes().pointLabels()[nearInfo.index()];
                    const point attraction = featPt-pp.localPoints()[pointi];

                    // Check if this point is already being attracted. If so
                    // override it only if nearer.
                    if
                    (
                        patchConstraints[pointi].first() <= 1
                     || magSqr(attraction) < magSqr(patchAttraction[pointi])
                    )
                    {
                        patchAttraction[pointi] = attraction;
                        patchConstraints[pointi] = edgeConstr[featEdgei][i];
                    }
                }
                else
                {
                    WarningInFunction
                        << "Did not find pp point near " << featPt
                        << endl;
                }
            }
        }
    }


    // Different procs might have different patchAttraction,patchConstraints
    // however these only contain geometric information, no topology
    // so as long as we synchronise after overriding with feature points
    // there is no problem, just possibly a small error.


    // Find nearest mesh point to feature point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (overrides attraction to feature edge)
    forAll(pointAttractor, feati)
    {
        const labelList& pointAttr = pointAttractor[feati];
        const List<pointConstraint>& pointConstr = pointConstraints[feati];

        forAll(pointAttr, featPointi)
        {
            if (pointAttr[featPointi] != -1)
            {
                const point& featPt = features[feati].points()
                [
                    featPointi
                ];

                // Find nearest pp point
                pointIndexHit nearInfo = ppTree.findNearest
                (
                    featPt,
                    sqr(GREAT)
                );

                if (nearInfo.hit())
                {
                    label pointi =
                        ppTree.shapes().pointLabels()[nearInfo.index()];

                    const point& pt = pp.localPoints()[pointi];
                    const point attraction = featPt-pt;

                    // - already attracted to feature edge : point always wins
                    // - already attracted to feature point: nearest wins

                    if (patchConstraints[pointi].first() <= 1)
                    {
                        patchAttraction[pointi] = attraction;
                        patchConstraints[pointi] = pointConstr[featPointi];
                    }
                    else if (patchConstraints[pointi].first() == 2)
                    {
                        patchAttraction[pointi] = attraction;
                        patchConstraints[pointi] = pointConstr[featPointi];
                    }
                    else if (patchConstraints[pointi].first() == 3)
                    {
                        // Only if nearer
                        if
                        (
                            magSqr(attraction)
                          < magSqr(patchAttraction[pointi])
                        )
                        {
                            patchAttraction[pointi] = attraction;
                            patchConstraints[pointi] =
                                pointConstr[featPointi];
                        }
                    }
                }
            }
        }
    }
}


void Foam::snappySnapDriver::featureAttractionUsingFeatureEdges
(
    const label iter,
    const bool multiRegionFeatureSnap,

    const bool detectBaffles,
    const bool baffleFeaturePoints,

    const bool releasePoints,
    const bool stringFeatures,
    const bool avoidDiagonal,

    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    const vectorField& nearestNormal,

    const List<List<point>>& pointFaceSurfNormals,
    const List<List<point>>& pointFaceDisp,
    const List<List<point>>& pointFaceCentres,
    const labelListList& pointFacePatchID,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    const refinementFeatures& features = meshRefiner_.features();
    const fvMesh& mesh = meshRefiner_.mesh();

    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh,
            pp.meshPoints()
        )
    );


    // Collect ordered attractions on feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Per feature, per feature-edge a list of attraction points and their
    // originating vertex.
    List<List<DynamicList<point>>> edgeAttractors(features.size());
    List<List<DynamicList<pointConstraint>>> edgeConstraints
    (
        features.size()
    );
    forAll(features, feati)
    {
        label nFeatEdges = features[feati].edges().size();
        edgeAttractors[feati].setSize(nFeatEdges);
        edgeConstraints[feati].setSize(nFeatEdges);
    }

    // Per feature, per feature-point the pp point that is attracted to it.
    // This list is only used to subset the feature-points that are actually
    // used.
    List<labelList> pointAttractor(features.size());
    List<List<pointConstraint>> pointConstraints(features.size());
    forAll(features, feati)
    {
        label nFeatPoints = features[feati].points().size();
        pointAttractor[feati].setSize(nFeatPoints, -1);
        pointConstraints[feati].setSize(nFeatPoints);
    }

    // Reverse: from pp point to nearest feature
    vectorField rawPatchAttraction(pp.nPoints(), Zero);
    List<pointConstraint> rawPatchConstraints(pp.nPoints());

    determineFeatures
    (
        iter,
        featureCos,
        multiRegionFeatureSnap,

        pp,
        snapDist,               // per point max distance and nearest surface
        nearestDisp,

        pointFaceSurfNormals,   // per face nearest surface
        pointFaceDisp,
        pointFaceCentres,
        pointFacePatchID,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,
        // pp point to nearest feature
        rawPatchAttraction,
        rawPatchConstraints
    );

    // Print a bit about the attraction from patch point to feature
    if (debug)
    {
        Info<< "Raw geometric feature analysis : ";
        writeStats(pp, isPatchMasterPoint, rawPatchConstraints);
    }

    // Baffle handling
    // ~~~~~~~~~~~~~~~
    // Override pointAttractor, edgeAttractor, rawPatchAttration etc. to
    // implement 'baffle' handling.
    // Baffle: the mesh pp point originates from a loose standing
    // baffle.
    // Sampling the surface with the surrounding face-centres only picks up
    // a single triangle normal so above determineFeatures will not have
    // detected anything. So explicitly pick up feature edges on the pp
    // (after duplicating points & smoothing so will already have been
    // expanded) and match these to the features.
    if (detectBaffles)
    {
        determineBaffleFeatures
        (
            iter,
            baffleFeaturePoints,
            featureCos,

            pp,
            snapDist,

            // Feature-point to pp point
            pointAttractor,
            pointConstraints,
            // Feature-edge to pp point
            edgeAttractors,
            edgeConstraints,
            // pp point to nearest feature
            rawPatchAttraction,
            rawPatchConstraints
        );
    }

    // Print a bit about the attraction from patch point to feature
    if (debug)
    {
        Info<< "After baffle feature analysis : ";
        writeStats(pp, isPatchMasterPoint, rawPatchConstraints);
    }


    // Reverse lookup: Find nearest mesh point to feature edge
    // ~~~~~~~~~~~~~~~~----------------~~~~~~~~~~~~~~~~~~~~~~~
    // go through all edgeAttractors and find the nearest point on pp

    reverseAttractMeshPoints
    (
        iter,

        pp,
        snapDist,

        // Feature-point to pp point
        pointAttractor,
        pointConstraints,
        // Feature-edge to pp point
        edgeAttractors,
        edgeConstraints,

        // Estimated feature point
        rawPatchAttraction,
        rawPatchConstraints,

        // pp point to nearest feature
        patchAttraction,
        patchConstraints
    );

    // Print a bit about the attraction from patch point to feature
    if (debug)
    {
        Info<< "Reverse attract feature analysis : ";
        writeStats(pp, isPatchMasterPoint, patchConstraints);
    }

    // Dump
    if (debug&meshRefinement::ATTRACTION)
    {
        OBJstream featureEdgeStr
        (
            meshRefiner_.mesh().time().path()
          / "edgeAttractors_" + name(iter) + ".obj"
        );
        Info<< "Dumping feature-edge attraction to "
            << featureEdgeStr.name() << endl;

        OBJstream featurePointStr
        (
            meshRefiner_.mesh().time().path()
          / "pointAttractors_" + name(iter) + ".obj"
        );
        Info<< "Dumping feature-point attraction to "
            << featurePointStr.name() << endl;

        forAll(patchConstraints, pointi)
        {
            const point& pt = pp.localPoints()[pointi];
            const vector& attr = patchAttraction[pointi];

            if (patchConstraints[pointi].first() == 2)
            {
                featureEdgeStr.write(linePointRef(pt, pt+attr));
            }
            else if (patchConstraints[pointi].first() == 3)
            {
                featurePointStr.write(linePointRef(pt, pt+attr));
            }
        }
    }


    //MEJ: any faces that have multi-patch points only keep the
    //     multi-patch
    //     points. The other points on the face will be dragged along
    //     (hopefully)
    if (releasePoints)
    {
        releasePointsNextToMultiPatch
        (
            iter,
            featureCos,

            pp,
            snapDist,

            pointFaceCentres,
            pointFacePatchID,

            rawPatchAttraction,
            rawPatchConstraints,

            patchAttraction,
            patchConstraints
        );
    }


    // Snap edges to feature edges
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Walk existing edges and snap remaining ones (that are marked as
    // feature edges in rawPatchConstraints)
    if (stringFeatures)
    {
        stringFeatureEdges
        (
            iter,
            featureCos,

            pp,
            snapDist,

            rawPatchAttraction,
            rawPatchConstraints,

            patchAttraction,
            patchConstraints
        );
    }


    // Avoid diagonal attraction
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Attract one of the non-diagonal points.
    if (avoidDiagonal)
    {
        avoidDiagonalAttraction
        (
            iter,
            featureCos,
            pp,
            patchAttraction,
            patchConstraints
        );
    }


    if (debug&meshRefinement::ATTRACTION)
    {
        dumpMove
        (
            meshRefiner_.mesh().time().path()
          / "patchAttraction_" + name(iter) + ".obj",
            pp.localPoints(),
            pp.localPoints() + patchAttraction
        );
    }
}


// Correct for squeezing of face
void Foam::snappySnapDriver::preventFaceSqueeze
(
    const label iter,
    const scalar featureCos,

    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const vectorField& nearestAttraction,

    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints
) const
{
    autoPtr<OBJstream> strPtr;
    if (debug&meshRefinement::ATTRACTION)
    {
        strPtr.reset
        (
            new OBJstream
            (
                meshRefiner_.mesh().time().path()
              / "faceSqueeze_" + name(iter) + ".obj"
            )
        );
        Info<< "Dumping faceSqueeze corrections to "
            << strPtr().name() << endl;
    }

    pointField points;
    face singleF;
    forAll(pp.localFaces(), facei)
    {
        const face& f = pp.localFaces()[facei];

        if (f.size() != points.size())
        {
            points.setSize(f.size());
            singleF.setSize(f.size());
            for (label i = 0; i < f.size(); i++)
            {
                singleF[i] = i;
            }
        }
        label nConstraints = 0;
        forAll(f, fp)
        {
            label pointi = f[fp];
            const point& pt = pp.localPoints()[pointi];

            if (patchConstraints[pointi].first() > 1)
            {
                points[fp] = pt + patchAttraction[pointi];
                nConstraints++;
            }
            else
            {
                points[fp] = pt;
            }
        }

        if (nConstraints == f.size())
        {
            if (f.size() == 3)
            {
                // Triangle: knock out attraction altogether

                // For now keep the points on the longest edge
                label maxFp = -1;
                scalar maxS = -1;
                forAll(f, fp)
                {
                    const point& pt = pp.localPoints()[f[fp]];
                    const point& nextPt = pp.localPoints()[f.nextLabel(fp)];

                    scalar s = magSqr(pt-nextPt);
                    if (s > maxS)
                    {
                        maxS = s;
                        maxFp = fp;
                    }
                }
                if (maxFp != -1)
                {
                    label pointi = f.prevLabel(maxFp);

                    // Reset attraction on pointi to nearest

                    const point& pt = pp.localPoints()[pointi];

                    //Pout<< "** on triangle " << pp.faceCentres()[facei]
                    //    << " knocking out attraction to " << pointi
                    //    << " at:" << pt
                    //    << endl;

                    patchAttraction[pointi] = nearestAttraction[pointi];

                    if (strPtr.valid())
                    {
                        strPtr().write
                        (
                            linePointRef(pt, pt+patchAttraction[pointi])
                        );
                    }
                }
            }
            else
            {
                scalar oldArea = f.mag(pp.localPoints());
                scalar newArea = singleF.mag(points);
                if (newArea < 0.1*oldArea)
                {
                    // For now remove the point with largest distance
                    label maxFp = -1;
                    scalar maxS = -1;
                    forAll(f, fp)
                    {
                        scalar s = magSqr(patchAttraction[f[fp]]);
                        if (s > maxS)
                        {
                            maxS = s;
                            maxFp = fp;
                        }
                    }
                    if (maxFp != -1)
                    {
                        label pointi = f[maxFp];
                        // Lower attraction on pointi
                        patchAttraction[pointi] *= 0.5;
                    }
                }
            }
        }
    }
}


Foam::vectorField Foam::snappySnapDriver::calcNearestSurfaceFeature
(
    const snapParameters& snapParams,
    const bool alignMeshEdges,
    const label iter,
    const scalar featureCos,
    const scalar featureAttract,
    const scalarField& snapDist,
    const vectorField& nearestDisp,
    const vectorField& nearestNormal,
    motionSmoother& meshMover,
    vectorField& patchAttraction,
    List<pointConstraint>& patchConstraints,

    DynamicList<label>& splitFaces,
    DynamicList<labelPair>& splits

) const
{
    const Switch implicitFeatureAttraction = snapParams.implicitFeatureSnap();
    const Switch explicitFeatureAttraction = snapParams.explicitFeatureSnap();
    const Switch multiRegionFeatureSnap = snapParams.multiRegionFeatureSnap();

    Info<< "Overriding displacement on features :" << nl
        << "   implicit features    : " << implicitFeatureAttraction << nl
        << "   explicit features    : " << explicitFeatureAttraction << nl
        << "   multi-patch features : " << multiRegionFeatureSnap << nl
        << endl;


    const indirectPrimitivePatch& pp = meshMover.patch();
    const pointField& localPoints = pp.localPoints();
    const fvMesh& mesh = meshRefiner_.mesh();


    //const PackedBoolList isMasterPoint(syncTools::getMasterPoints(mesh));
    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh,
            pp.meshPoints()
        )
    );

    // Per point, per surrounding face:
    // - faceSurfaceNormal
    // - faceDisp
    // - faceCentres
    List<List<point>> pointFaceSurfNormals;
    List<List<point>> pointFaceDisp;
    List<List<point>> pointFaceCentres;
    List<labelList> pointFacePatchID;

    {
        // Calculate attraction distance per face (from the attraction distance
        // per point)
        scalarField faceSnapDist(pp.size(), -GREAT);
        forAll(pp.localFaces(), facei)
        {
            const face& f = pp.localFaces()[facei];
            forAll(f, fp)
            {
                faceSnapDist[facei] = max(faceSnapDist[facei], snapDist[f[fp]]);
            }
        }


        // Displacement and orientation per pp face
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // vector from point on surface back to face centre
        vectorField faceDisp(pp.size(), Zero);
        // normal of surface at point on surface
        vectorField faceSurfaceNormal(pp.size(), Zero);
        labelList faceSurfaceGlobalRegion(pp.size(), -1);
        //vectorField faceRotation(pp.size(), Zero);

        calcNearestFace
        (
            iter,
            pp,
            faceSnapDist,
            faceDisp,
            faceSurfaceNormal,
            faceSurfaceGlobalRegion
            //faceRotation
        );


        // Collect (possibly remote) per point data of all surrounding faces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // - faceSurfaceNormal
        // - faceDisp
        // - faceCentres
        calcNearestFacePointProperties
        (
            iter,
            pp,

            faceDisp,
            faceSurfaceNormal,
            faceSurfaceGlobalRegion,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID
        );
    }


    // Start off with nearest point on surface
    vectorField patchDisp = nearestDisp;


    // Main calculation
    // ~~~~~~~~~~~~~~~~
    // This is the main intelligence which calculates per point the vector to
    // attract it to the nearest surface. There are lots of possibilities
    // here.

    // Nearest feature
    patchAttraction.setSize(localPoints.size());
    patchAttraction = Zero;
    // Constraints at feature
    patchConstraints.setSize(localPoints.size());
    patchConstraints = pointConstraint();

    if (implicitFeatureAttraction)
    {
        // Sample faces around each point and see if nearest surface normal
        // differs. Reconstruct a feature edge/point if possible and snap to
        // it.
        featureAttractionUsingReconstruction
        (
            iter,
            featureCos,

            pp,
            snapDist,
            nearestDisp,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            patchAttraction,
            patchConstraints
        );
    }

    if (explicitFeatureAttraction)
    {
        // Only do fancy stuff if alignMeshEdges
        bool releasePoints = false;
        bool stringFeatures = false;
        bool avoidDiagonal = false;
        if (alignMeshEdges)
        {
            releasePoints = snapParams.releasePoints();
            stringFeatures = snapParams.stringFeatures();
            avoidDiagonal = snapParams.avoidDiagonal();
        }


        // Sample faces around each point and see if nearest surface normal
        // differs. For those find the nearest real feature edge/point and
        // store the correspondence. Then loop over feature edge/point
        // and attract those nearest mesh point. (the first phase just is
        // a subsetting of candidate points, the second makes sure that only
        // one mesh point gets attracted per feature)
        featureAttractionUsingFeatureEdges
        (
            iter,
            multiRegionFeatureSnap,

            snapParams.detectBaffles(),
            snapParams.baffleFeaturePoints(),   // all points on baffle edges
                                                // are attracted to feature pts

            releasePoints,
            stringFeatures,
            avoidDiagonal,

            featureCos,

            pp,
            snapDist,
            nearestDisp,
            nearestNormal,

            pointFaceSurfNormals,
            pointFaceDisp,
            pointFaceCentres,
            pointFacePatchID,

            patchAttraction,
            patchConstraints
        );
    }

    if (!alignMeshEdges)
    {
        const scalar concaveCos = Foam::cos
        (
            degToRad(snapParams.concaveAngle())
        );
        const scalar minAreaRatio = snapParams.minAreaRatio();

        Info<< "Experimental: introducing face splits to avoid rotating"
            << " mesh edges. Splitting faces when" << nl
            << indent << "- angle not concave by more than "
            << snapParams.concaveAngle() << " degrees" << nl
            << indent << "- resulting triangles of similar area "
            << " (ratio within " << minAreaRatio << ")" << nl
            << endl;

        splitDiagonals
        (
            featureCos,
            concaveCos,
            minAreaRatio,
            pp,

            nearestDisp,
            nearestNormal,

            patchAttraction,
            patchConstraints,
            splitFaces,
            splits
        );

        if (debug)
        {
            Info<< "Diagonal attraction feature correction : ";
            writeStats(pp, isPatchMasterPoint, patchConstraints);
        }
    }


    preventFaceSqueeze
    (
        iter,
        featureCos,

        pp,
        snapDist,
        nearestDisp,

        patchAttraction,
        patchConstraints
    );

    {
        vector avgPatchDisp = meshRefinement::gAverage
        (
            isPatchMasterPoint,
            patchDisp
        );
        vector avgPatchAttr = meshRefinement::gAverage
        (
            isPatchMasterPoint,
            patchAttraction
        );

        Info<< "Attraction:" << endl
            << "    linear   : max:" << gMaxMagSqr(patchDisp)
            << " avg:" << avgPatchDisp << endl
            << "    feature  : max:" << gMaxMagSqr(patchAttraction)
            << " avg:" << avgPatchAttr << endl;
    }

    // So now we have:
    // - patchDisp          : point movement to go to nearest point on surface
    //                       (either direct or through interpolation of
    //                        face nearest)
    // - patchAttraction    : direct attraction to features
    // - patchConstraints   : type of features

    // Use any combination of patchDisp and direct feature attraction.


    // Mix with direct feature attraction
    forAll(patchConstraints, pointi)
    {
        if (patchConstraints[pointi].first() > 1)
        {
            patchDisp[pointi] =
                (1.0-featureAttract)*patchDisp[pointi]
              + featureAttract*patchAttraction[pointi];
        }
    }



    // Count
    {
        Info<< "Feature analysis : ";
        writeStats(pp, isPatchMasterPoint, patchConstraints);
    }


    // Now we have the displacement per patch point to move onto the surface
    // Split into tangential and normal direction.
    // - start off with all non-constrained points following the constrained
    //   ones since point normals not relevant.
    // - finish with only tangential component smoothed.
    // Note: tangential is most
    // likely to come purely from face-centre snapping, not face rotation.
    // Note: could use the constraints here (constraintTransformation())
    //       but this is not necessarily accurate and we're smoothing to
    //       get out of problems.

    if (featureAttract < 1-0.001)
    {
        //const PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
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

        const vectorField pointNormals
        (
            PatchTools::pointNormals
            (
                mesh,
                pp
            )
        );

        // 1. Smoothed all displacement
        vectorField smoothedPatchDisp = patchDisp;
        smoothAndConstrain
        (
            isPatchMasterEdge,
            pp,
            meshEdges,
            patchConstraints,
            smoothedPatchDisp
        );


        // 2. Smoothed tangential component
        vectorField tangPatchDisp = patchDisp;
        tangPatchDisp -= (pointNormals & patchDisp) * pointNormals;
        smoothAndConstrain
        (
            isPatchMasterEdge,
            pp,
            meshEdges,
            patchConstraints,
            tangPatchDisp
        );

        // Re-add normal component
        tangPatchDisp += (pointNormals & patchDisp) * pointNormals;

        if (debug&meshRefinement::ATTRACTION)
        {
            dumpMove
            (
                mesh.time().path()
              / "tangPatchDispConstrained_" + name(iter) + ".obj",
                pp.localPoints(),
                pp.localPoints() + tangPatchDisp
            );
        }

        patchDisp =
             (1.0-featureAttract)*smoothedPatchDisp
           + featureAttract*tangPatchDisp;
    }


    const scalar relax = featureAttract;
    patchDisp *= relax;


    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value (note: cant use VGREAT)
    );

    return patchDisp;
}


void Foam::snappySnapDriver::findSurfaceNormals
(
    const labelList& checkSurfaces,
    const indirectPrimitivePatch& pp,
    const scalarField& snapDist,
    const pointField& newPoints,
    const pointField& ppFaceCentres,
    const pointField& ppCellCentres,
    const vectorField& ppFaceAreas,
    const boolList& stationaryFaces,
    const bool& fastFeatureSnap,
    List<pointIndexHit>& finalHitInfo,
    labelList& pointIndex,
    vectorField& pNormals
) const
{
//#ifdef FOAM_USE_TBB
//    Timer timer("snappySnapDriver::findSurfaceNormals");
//#endif

    const fvMesh& mesh = meshRefiner_.mesh();

    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    // Check for alignement of face with surface
    scalar alignmentAngle = Foam::cos(degToRad(45.));

    pointField interPoints(pp.size(),vector::zero);

    vectorField pointCellCentres(pp.nPoints(), vector::zero);
    vectorField pointFaceNormals(pp.nPoints(), vector::zero);
    labelList nPointCounter(pp.nPoints(), 0);

    vectorField surfHitPatchPointNormals(pp.nPoints(), vector::zero);
    vectorField surfHitPatchCentreNormals(pp.size(), vector::zero);

    // Displacement data from interrogation of surface
    pointField patchPointDisp(pp.nPoints(), vector::zero);
    pointField patchFaceDisp(pp.size(), vector::zero);

    vectorField alignedSurfaceNormal(pp.size(), vector::zero);

    vectorField aveFaceNormals(pp.size(), vector::zero);
    vectorField avePointNormals(pp.nPoints(), vector::zero);

    // 1. All points to non-interface surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (!fastFeatureSnap)
    {
        pointCellCentres = vector::zero;
        pointFaceNormals = vector::zero;
        nPointCounter = 0;

        forAll(pp.pointFaces(), pointI)
        {
            const labelList& pFaces = pp.pointFaces()[pointI];

            forAll(pFaces, faceI)
            {
                pointCellCentres[pointI] +=
                    (ppCellCentres[pFaces[faceI]]-ppFaceCentres[pFaces[faceI]]);
                pointFaceNormals[pointI] += ppFaceAreas[pFaces[faceI]];
            }
            nPointCounter[pointI] += pFaces.size();
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            pointCellCentres,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            pointFaceNormals,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPointCounter,
            plusEqOp<label>(),
            label(0) // null value
        );

        forAll(pointCellCentres, pointI)
        {
            pointCellCentres[pointI] /= nPointCounter[pointI];
            pointCellCentres[pointI] += newPoints[pp.meshPoints()[pointI]];
            pointFaceNormals[pointI] /= nPointCounter[pointI];
        }

        // Use first intersected point if available,
        // otherwise find nearest surface
        vectorField endPoints(pp.nPoints(), vector::zero);
        forAll(endPoints, pointI)
        {
            endPoints[pointI] = newPoints[pp.meshPoints()[pointI]];
        }

        labelList surface1;
        List<pointIndexHit> hitPoint1;
        labelList pointRegion1;
        labelList surface2;
        List<pointIndexHit> hitPoint2;
        labelList pointRegion2;

        surfaces.findNearestIntersection
        (
            checkSurfaces,
            pointCellCentres,
            endPoints,

            surface1,
            hitPoint1,
            pointRegion1,

            surface2,
            hitPoint2,
            pointRegion2
        );

        forAll(hitPoint1, pointI)
        {
            if (hitPoint1[pointI].hit())
            {
                endPoints[pointI] = hitPoint1[pointI].hitPoint();
            }
            else
            {
                endPoints[pointI] = newPoints[pp.meshPoints()[pointI]];
            }
        }

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces.findNearest
        (
            checkSurfaces,
            endPoints,
            sqr(snapDist),        // sqr of attract distance
            hitSurface,
            hitInfo
         );

        vectorField pHitNormals(hitInfo.size());
        forAll(checkSurfaces, sI)
        {
            DynamicList<pointIndexHit> localHits;
            label surfI = checkSurfaces[sI];

            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    localHits.append(hitInfo[i]);
                }
            }
            pointField localNormals;
            label geomI = surfaces.surfaces()[surfI];
            surfaces.geometry()[geomI].getNormal(localHits, localNormals);

            label localI = 0;
            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    pHitNormals[i] = localNormals[localI];
                    localI++;
                }
            }
        }

        forAll(hitInfo, pointI)
        {
            if (hitInfo[pointI].hit())
            {
                surfHitPatchPointNormals[pointI] = pHitNormals[pointI];

                //make the hit point normal have opposite orientation
                //to the mean surrounding face normal
                if ((surfHitPatchPointNormals[pointI] & pointFaceNormals[pointI]) > 0)
                {
                    surfHitPatchPointNormals[pointI] =
                        -surfHitPatchPointNormals[pointI];
                }

                patchPointDisp[pointI] = hitInfo[pointI].hitPoint()
                    - newPoints[pp.meshPoints()[pointI]];
            }
            else
            {
                WarningInFunction
                    << "For point:" << pointI
                    << " coordinate:" <<  newPoints[pp.meshPoints()[pointI]]
                    << " did not find any surface within:"
                    << sqr(snapDist[pointI])
                    << " meter." << endl;
            }
        }
    }

    // 2. All face centres to non-interface surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        scalarField maxSnapDist(pp.size(), -GREAT);

        forAll(pp, i)
        {
            const face& f = pp.localFaces()[i];
            forAll(f, fp)
            {
                maxSnapDist[i] = max(snapDist[f[fp]], maxSnapDist[i]);
            }
        }

        // Use first intersected point if available,
        // otherwise find nearest surface
        pointField endPoints(pp.size(), vector::zero);
        forAll(endPoints, pointI)
        {
            endPoints[pointI] = ppFaceCentres[pointI];
        }

        labelList surface1;
        List<pointIndexHit> hitPoint1;
        labelList pointRegion1;
        labelList surface2;
        List<pointIndexHit> hitPoint2;
        labelList pointRegion2;

        surfaces.findNearestIntersection
        (
            checkSurfaces,
            ppCellCentres,
            ppFaceCentres,

            surface1,
            hitPoint1,
            pointRegion1,

            surface2,
            hitPoint2,
            pointRegion2
         );

        forAll(hitPoint1, pointI)
        {
            if (hitPoint1[pointI].hit())
            {
                endPoints[pointI] = hitPoint1[pointI].hitPoint();
            }
            else
            {
                endPoints[pointI] = ppFaceCentres[pointI];
            }
        }

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces.findNearest
        (
            checkSurfaces,
            endPoints,
            sqr(maxSnapDist),        // sqr of attract distance
            hitSurface,
            hitInfo
         );

        labelList pointCentreIndex(pp.size());
        pointField startPoints(pp.size());

        label testI = 0;
        forAll(pp, i)
        {
            if (hitInfo[i].hit())
            {
                vector areaNorm = ppFaceAreas[i]
                    / (mag(ppFaceAreas[i]) + SMALL);

                vector hitDir = hitInfo[i].hitPoint() - ppFaceCentres[i];

                if ((hitDir & areaNorm) > 0.)
                {
                    scalar distToSurf = mag(hitDir);
                    vector searchDir = 2.0 * distToSurf * areaNorm;
                    startPoints[testI] = ppFaceCentres[i];
                    endPoints[testI] = ppFaceCentres[i] + searchDir;
                    pointCentreIndex[testI] = i;
                    testI++;
                }
            }
        }
        pointCentreIndex.setSize(testI);
        startPoints.setSize(testI);
        endPoints.setSize(testI);

        labelList surf1;
        List<pointIndexHit> hPoint1;
        labelList pRegion1;
        labelList surf2;
        List<pointIndexHit> hPoint2;
        labelList pRegion2;

        surfaces.findNearestIntersection
        (
            checkSurfaces,
            startPoints,
            endPoints,

            surf1,
            hPoint1,
            pRegion1,

            surf2,
            hPoint2,
            pRegion2
         );

        vectorField nearestNormal(startPoints.size());
        vectorField origNormal(startPoints.size());

        forAll(checkSurfaces, sI)
        {
            DynamicList<pointIndexHit> localHits;
            DynamicList<pointIndexHit> prevHits;
            label surfI = checkSurfaces[sI];

            forAll(surf1, i)
            {
                if (surf1[i] == surfI)
                {
                    localHits.append(hPoint1[i]);
                }
                if (hitSurface[pointCentreIndex[i]] == surfI)
                {
                    prevHits.append(hitInfo[pointCentreIndex[i]]);
                }
            }
            label geomI = surfaces.surfaces()[surfI];
            pointField localNormals;
            surfaces.geometry()[geomI].getNormal(localHits, localNormals);
            pointField prevNormals;
            surfaces.geometry()[geomI].getNormal(prevHits, prevNormals);

            label localI = 0;
            label prevI = 0;
            forAll(surf1, i)
            {
                if (surf1[i] == surfI)
                {
                    nearestNormal[i] = localNormals[localI];
                    localI++;
                }
                if (hitSurface[pointCentreIndex[i]] == surfI)
                {
                    origNormal[i] = prevNormals[prevI];
                    prevI++;
                }
            }
        }

        forAll(surf1, i)
        {
            if (hPoint1[i].hit())
            {
                vector areaNorm = ppFaceAreas[pointCentreIndex[i]]
                    / (mag(ppFaceAreas[pointCentreIndex[i]]) + SMALL);
                vector vec1 = hPoint1[i].hitPoint()
                    - ppCellCentres[pointCentreIndex[i]];
                vector projSurfNorm = nearestNormal[i];
                if ((projSurfNorm & vec1) > 0)
                {
                    projSurfNorm  = -projSurfNorm;
                }

                vector surfNorm = origNormal[i];
                vec1 = hitInfo[pointCentreIndex[i]].hitPoint()
                    - ppCellCentres[pointCentreIndex[i]];
                if ((surfNorm & vec1) > 0)
                {
                    surfNorm  = -surfNorm;
                }
                scalar dot1 = -(projSurfNorm
                                / (mag(projSurfNorm) + SMALL)) & areaNorm;
                scalar dot2 = -(surfNorm
                                / (mag(surfNorm) + SMALL)) & areaNorm;

                if (dot1 > 1.1*dot2)
                {
                    hitSurface[pointCentreIndex[i]] = surf1[i];
                    hitInfo[pointCentreIndex[i]] = hPoint1[i];
                }
            }
        }

        forAll(checkSurfaces, sI)
        {
            DynamicList<pointIndexHit> localHits;
            label surfI = checkSurfaces[sI];

            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    localHits.append(hitInfo[i]);
                }
            }
            pointField localNormals;
            label geomI = surfaces.surfaces()[surfI];
            surfaces.geometry()[geomI].getNormal(localHits, localNormals);

            label localI = 0;
            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    surfHitPatchCentreNormals[i] = localNormals[localI];
                    vector vec1 =
                        ppFaceAreas[i] / (mag(ppFaceAreas[i]) + SMALL);
                    if ((surfHitPatchCentreNormals[i] & vec1) > 0)
                    {
                        surfHitPatchCentreNormals[i] =
                            -surfHitPatchCentreNormals[i];
                    }

                    if (hitSurface[i] != -1)
                    {
                        patchFaceDisp[i] =
                            hitInfo[i].hitPoint() - ppFaceCentres[i];
                    }
                    else
                    {
                        WarningInFunction
                            << "For point:" << i
                            << " coordinate:" << ppFaceCentres[i]
                            << " did not find any surface within:"
                            << sqr(maxSnapDist[i])
                            << " meter." << endl;
                    }
                    localI++;
                }
            }
        }
    }


    vector fNorm(vector::zero);
    vector sNorm(vector::zero);
    vector alignedNormal(vector::zero);
    scalar dotProd;

    if (!fastFeatureSnap)
    {
        forAll(pp, i)
        {
            const face& f = pp.localFaces()[i];

            SortableList<scalar> maxList(f.size()+1);

            fNorm = ppFaceAreas[i];
            fNorm = fNorm / (mag(fNorm) + SMALL);
            // check face points
            forAll(f, fp)
            {
                sNorm = surfHitPatchPointNormals[f[fp]];
                sNorm = sNorm / (mag(sNorm) + SMALL);
                dotProd = -(fNorm & sNorm);
                maxList[fp] = dotProd;
            }

            // check face centres
            sNorm = surfHitPatchCentreNormals[i];
            sNorm = sNorm / (mag(sNorm) + SMALL);

            dotProd = -(fNorm & sNorm);
            maxList[f.size()] = dotProd;
            maxList.sort();

            label maxI = maxList.indices()[f.size()];

            if (maxI == f.size())
            {
                alignedSurfaceNormal[i] = surfHitPatchCentreNormals[i];
            }
            else
            {
                alignedSurfaceNormal[i] = surfHitPatchPointNormals[f[maxI]];
            }
        }

        avePointNormals = vector::zero;
        nPointCounter = 0;

        forAll(pp.pointFaces(), pointI)
        {
            const labelList& pFaces = pp.pointFaces()[pointI];

            forAll(pFaces, faceI)
            {
                fNorm = ppFaceAreas[pFaces[faceI]];
                fNorm = fNorm / (mag(fNorm) + SMALL);
                sNorm = alignedSurfaceNormal[pFaces[faceI]];
                dotProd = -fNorm & sNorm;

                if (dotProd > alignmentAngle)
                {
                    avePointNormals[pointI] +=
                        alignedSurfaceNormal[pFaces[faceI]];
                    nPointCounter[pointI]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            avePointNormals,
            plusEqOp<vector>(),
            vector::zero       // null value
         );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPointCounter,
            plusEqOp<label>(),
            label(0) // null value
         );

        // find average of neighbouring faces aligned with alignedSurfaceNormal
        forAll(pp, i)
        {
            const face& f = pp.localFaces()[i];
            label nCount = 0;
            forAll(f, fp)
            {
                aveFaceNormals[i] += avePointNormals[f[fp]];
                nCount += nPointCounter[f[fp]];
            }

            if (nCount > 0 && mag(aveFaceNormals[i]) > 0)
            {
                aveFaceNormals[i] /= nCount;
            }
            else
            {
                aveFaceNormals[i] = alignedSurfaceNormal[i];
            }
        }
    }

    scalarField maxSnap(pp.size());
    label testI = 0;

    forAll(pp, i)
    {
        if (!stationaryFaces[i])
        {
            const face& f = pp.localFaces()[i];
            const point& fc = ppFaceCentres[i];

            point avePoint = vector::zero;
            if (fastFeatureSnap)
            {
                avePoint = (patchFaceDisp[i] + fc);
            }
            else
            {
                fNorm = ppFaceAreas[i];
                fNorm = fNorm / (mag(fNorm) + SMALL);
                sNorm = alignedSurfaceNormal[i];
                dotProd = -fNorm & sNorm;
                // if not well aligned use average of neighbouring face
                //normals to find closest alignment

                if (dotProd < alignmentAngle)
                {
                    SortableList<scalar> maxList(f.size()+1);
                    fNorm = aveFaceNormals[i];
                    fNorm = fNorm / (mag(fNorm) + SMALL);
                    // check face points
                    forAll(f, fp)
                    {
                        sNorm = surfHitPatchPointNormals[f[fp]];
                        sNorm = sNorm / (mag(sNorm) + SMALL);
                        dotProd = (fNorm & sNorm);
                        maxList[fp] = dotProd;
                    }
                    // check face centres
                    sNorm = surfHitPatchCentreNormals[i];
                    sNorm = sNorm / (mag(sNorm) + SMALL);

                    dotProd = (fNorm & sNorm);

                    maxList[f.size()] = dotProd;
                    maxList.sort();
                    label maxI = maxList.indices()[f.size()];
                    if (maxI == f.size())
                    {
                        alignedNormal = surfHitPatchCentreNormals[i];
                    }
                    else
                    {
                        alignedNormal = surfHitPatchPointNormals[f[maxI]];
                    }
                }
                else
                {
                    alignedNormal = alignedSurfaceNormal[i];
                }
                alignedNormal /= (mag(alignedNormal) + SMALL);

                label nPointsFound = 0;
                //use face points that are aligned to define where to move face
                forAll(f, fp)
                {
                    sNorm = surfHitPatchPointNormals[f[fp]];
                    sNorm = sNorm / (mag(sNorm) + SMALL);
                    dotProd = alignedNormal & sNorm;
                    if (dotProd > alignmentAngle)
                    {
                        avePoint += newPoints[pp.meshPoints()[f[fp]]]
                            + patchPointDisp[f[fp]];
                        nPointsFound++;
                    }
                }
                sNorm = surfHitPatchCentreNormals[i];
                sNorm = sNorm / (mag(sNorm) + SMALL);
                dotProd = alignedNormal & sNorm;

                if (dotProd > alignmentAngle || nPointsFound == 0)
                {
                    avePoint = (patchFaceDisp[i] + fc);
                }
                else
                {
                    avePoint /= nPointsFound;
                }
            }

            scalar maxSnapDist = -GREAT;
            forAll(f, fp)
            {
                maxSnapDist = max(snapDist[f[fp]], maxSnapDist);
            }
            interPoints[testI] = avePoint;
            maxSnap[testI] = sqr(maxSnapDist);
            pointIndex[testI] = i;

            testI++;
        }
    }

    interPoints.setSize(testI);
    maxSnap.setSize(testI);
    pointIndex.setSize(testI);

    labelList hitSurface;
    // find closest point on surface using average point
    surfaces.findNearest
    (
        checkSurfaces,
        interPoints,
        maxSnap,        // sqr of attract distance
        hitSurface,
        finalHitInfo
     );

    pNormals.setSize(interPoints.size(),vector::zero);
    forAll(checkSurfaces, sI)
    {
        DynamicList<pointIndexHit> localHits;
        label surfI = checkSurfaces[sI];

        forAll(hitSurface, i)
        {
            if (hitSurface[i] == surfI)
            {
                localHits.append(finalHitInfo[i]);
            }
        }

        pointField localNormals;
        label geomI = surfaces.surfaces()[surfI];
        surfaces.geometry()[geomI].getNormal(localHits, localNormals);

        label localI = 0;
        forAll(hitSurface, i)
        {
            if (hitSurface[i] == surfI)
            {
                pNormals[i] = localNormals[localI];
                localI++;
            }
        }
    }
}


void Foam::snappySnapDriver::featureSnapInit
(
    const bool initialSnap,
    pointField& newPoints,
    vectorField& pointDisp,
    scalarField& sumAreasMesh,
    vectorField& sumDispMesh,
    boolList& gapBoundaryPts
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    if (mesh.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh.lookupObject<volScalarField>("gapCells");
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(mesh.cells(), celli)
        {
            if (gapCells[celli] > -1)
            {
                const labelList& cFaces = mesh.cells()[celli];
                forAll(cFaces, cFI)
                {
                    label facei = cFaces[cFI];

                    label patchi = patches.whichPatch(facei);

                    if (patchi != -1 && !patches[patchi].coupled())
                    {
                        const labelList& f = mesh.faces()[facei];
                        forAll(f,fp)
                        {
                            gapBoundaryPts[f[fp]] =  true;
                        }
                    }
                }
            }
        }

        boolList internalFaces(mesh.nFaces(),false);
        scalarField neiGapCells;
        syncTools::swapBoundaryCellList(mesh, gapCells, neiGapCells);

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

            if ((gapOwn > -1) && (gapOwn == gapNei))
            {
                internalFaces[facei] = true;
            }
        }

        syncTools::syncFaceList(mesh, internalFaces, orEqOp<bool>());
        forAll(mesh.faces(), facei)
        {
            if (!internalFaces[facei])
            {
                continue;
            }

            const face f = mesh.faces()[facei];
            forAll(f,fp)
            {
                label pointi = f[fp];
                if (gapBoundaryPts[pointi])
                {
                    gapBoundaryPts[pointi] = false;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            gapBoundaryPts,
            orEqOp<bool>(),
            false           // null value
         );
    } // if foundGapCells

    if (initialSnap)
    {
        const refinementSurfaces& surfaces = meshRefiner_.surfaces();

        // Undistorted edge length
        const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
        const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

        //First snap zone then non zoned surfaces
        for (label iter = 0; iter < 2; iter++)
        {
            pointDisp = vector::zero;
            labelList checkSurfaces;
            labelList initSnapPatches;

            labelList zonePatches = zonePatchIDs();

            if (iter == 0)
            {
                checkSurfaces = surfaceZonesInfo::getNonBoundaryNamedSurfaces
                (
                    surfaces.surfZones()
                );
                initSnapPatches = zonePatchSlaveIDs();
            }
            else
            {
                checkSurfaces =
                    surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
                    (
                        surfaces.surfZones()
                    );

                labelList adaptPatchIDs
                (
                    meshRefiner_.meshedPatches()
                );
                labelHashSet excludePatches(zonePatches);

                initSnapPatches.setSize(adaptPatchIDs.size());
                label nInitPatches = 0;
                forAll(adaptPatchIDs, i)
                {
                    if (!excludePatches.found(adaptPatchIDs[i]))
                    {
                        initSnapPatches[nInitPatches++] = adaptPatchIDs[i];
                    }
                }
                initSnapPatches.setSize(nInitPatches);
            }

            if (initSnapPatches.size() == 0)
            {
                continue;
            }

            autoPtr<indirectPrimitivePatch> ppInitPtr
            (
                meshRefinement::makePatch
                (
                    mesh,
                    initSnapPatches
                )
            );
            indirectPrimitivePatch& ppInit = ppInitPtr();
            const labelList& meshPointsInit = ppInit.meshPoints();
            const labelList meshEdgesInit
            (
                ppInit.meshEdges(mesh.edges(), mesh.pointEdges())
             );
            boolList exclFaces(ppInit.size(), false);
            edgeClassification eClass
            (
                mesh,
                mesh.points(),
                ppInit,
                meshEdgesInit,
                exclFaces,
                -0.5,
                -0.5
             );
            const List<Tuple2<edgeClassification::edgeType,scalar>>&
                eType = eClass.edgeTypes();

            pointField pointNormals =
                eClass.calculatePointNormals(exclFaces, 8, false);

            scalarField maxProjection(meshPointsInit.size(), GREAT);
            boolList nonManifold(meshPointsInit.size(), false);
            forAll(ppInit.edges(), patchEdgeI)
            {
                const edge& e = ppInit.edges()[patchEdgeI];
                label v0 = e[0];
                label v1 = e[1];

                if
                (
                    eType[patchEdgeI].first()
                    == edgeClassification::NONMANIFOLD
                )
                {
                    nonManifold[v0] = true;
                    nonManifold[v1] = true;
                }

                vector n0 = pointNormals[v0];
                vector n1 = pointNormals[v1];

                if ((n0&n1) < 0.99)
                {
                    vector eVec = e.vec(ppInit.localPoints());
                    scalar eMag = mag(eVec) + SMALL;
                    eVec /= eMag;
                    scalar dProd0 = (n0&eVec);
                    scalar dProd1 = (n1&eVec);
                    scalar theta0 =
                        Foam::acos(min(max(dProd0, -scalar(1)), scalar(1)));
                    scalar theta1 =
                        Foam::acos(min(max(dProd1, -scalar(1)), scalar(1)));

                    theta1 = Foam::constant::mathematical::pi - theta1;

                    scalar theta2 = Foam::constant::mathematical::pi
                        - theta0 - theta1;

                    if (theta2 < SMALL)
                    {
                        continue;
                    }

                    scalar l1 = eMag/Foam::sin(theta2);
                    scalar inter0 = Foam::sin(theta1)*l1;
                    scalar inter1 = Foam::sin(theta0)*l1;

                    maxProjection[v0] = min(maxProjection[v0], inter0);
                    maxProjection[v1] = min(maxProjection[v1], inter1);
                }
            }

            syncTools::syncPointList
            (
                mesh,
                meshPointsInit,
                maxProjection,
                minEqOp<scalar>(),
                GREAT        // null value
            );

            pointField start(ppInit.size());
            pointField end(ppInit.size());
            forAll(ppInit, i)
            {
                face f = ppInit.localFaces()[i];
                vector aveNormal = vector::zero;
                forAll(f, fp)
                {
                    aveNormal += pointNormals[f[fp]];
                }
                aveNormal /= f.size();
                aveNormal /= (mag(aveNormal) + SMALL);

                label meshFaceI = ppInit.addressing()[i];
                label own = mesh.faceOwner()[meshFaceI];

                start[i] = mesh.cellCentres()[own];

                label  level = cellLevel[own];
                scalar len = edge0Len / pow(2., level);

                scalar backwardLength(iter == 0 ? 2*len : 0.5*len);
                scalar forwardLength(iter == 0 ? 2*len : 5*len);

                end[i] = start[i] + forwardLength*aveNormal;
                start[i] = start[i] - backwardLength*aveNormal;
            }

            labelList surface1;
            List<pointIndexHit> hitPoint1;
            labelList pointRegion1;
            labelList surface2;
            List<pointIndexHit> hitPoint2;
            labelList pointRegion2;

            //GGG
            surfaces.findNearestIntersection
            (
                checkSurfaces,
                start,
                end,
                surface1,
                hitPoint1,
                pointRegion1,

                surface2,
                hitPoint2,
                pointRegion2
             );

            forAll(start, i)
            {
                if (hitPoint1[i].hit())
                {
                    face f = ppInit.localFaces()[i];
                    scalar maxProjectionDistance = GREAT;
                    forAll(f,fp)
                    {
                        maxProjectionDistance =
                            min(maxProjection[f[fp]],maxProjectionDistance);
                    }

                    vector hitVec = hitPoint1[i].hitPoint()
                        - ppInit.faceCentres()[i];
                    scalar hitDist = mag(hitVec);

                    if (hitDist > SMALL && hitDist > maxProjectionDistance)
                    {
                        hitVec /= hitDist;
                        hitVec*= maxProjectionDistance;
                    }

                    label meshFaceI = ppInit.addressing()[i];

                    forAll(f,fp)
                    {
                        scalar weight = 1.
                            / (mesh.magFaceAreas()[meshFaceI] + SMALL);
                        if (nonManifold[f[fp]])
                        {
                            continue;
                        }
                        label meshPointI = ppInit.meshPoints()[f[fp]];
                        sumDispMesh[meshPointI] += weight*hitVec;
                        sumAreasMesh[meshPointI] += weight;
                    }
                }
            }
        } // for iter (initialSnap)

        syncTools::syncPointList
        (
            mesh,
            sumDispMesh,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            sumAreasMesh,
            plusEqOp<scalar>(),
            scalar(0.0)         // null value
        );

        forAll(mesh.points(), meshPointI)
        {
            if (sumAreasMesh[meshPointI] > SMALL && !gapBoundaryPts[meshPointI])
            {
                pointDisp[meshPointI] = sumDispMesh[meshPointI]
                    / sumAreasMesh[meshPointI];
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointDisp,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        forAll(mesh.points(), meshPointI)
        {
            newPoints[meshPointI] = newPoints[meshPointI]
                + pointDisp[meshPointI];
        }
    }

    return;
}


void Foam::snappySnapDriver::featureOptimise
(
    const label iter,
    const optStages& ostages,
    const snapParameters& snapParams,
    const dictionary& motionDict,
    const indirectPrimitivePatch& pp,
    const labelList& snapSurfaces,
    const labelList& meshEdges,
    const boolList& boundaryPoints,
    const boolList& gapBoundaryPts,
    pointField& newPoints,
    motionSmoother& meshMover,
    boolList& stationaryPoints
) const
{
    const direction opt = ostages.find(iter);

    fvMesh& mesh = meshRefiner_.mesh();
    const labelList& meshPoints = pp.meshPoints();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    if
    (
        (
            controller_.algorithm() == meshControl::EXTRUDE
            || controller_.algorithm() == meshControl::DUAL
         )
        && controller_.mode() != meshControl::DRYRUN
    )
    {
        //final iteration optimisation
        if
        (
            !(mesh.foundObject<volScalarField>("gapCells"))
            &&
            (
                controller_.algorithm() == meshControl::EXTRUDE
                || controller_.mode() == meshControl::QUALITY
                || controller_.mode() == meshControl::BALANCED
             )
            && opt != optStages::NONE
        )
        {
            if
            (
                opt & optStages::EXTRUDE
                && controller_.algorithm() == meshControl::EXTRUDE
            )
            {
                extrudeCorrect
                (
                    pp,
                    motionDict,
                    boundaryPoints,
                    meshMover,
                    newPoints
                );

                cornerCorrect
                (
                    pp,
                    meshEdges,
                    boundaryPoints,
                    newPoints
                );
            }

            //foam optimisation step
            if (opt & optStages::SMOOTH)
            {
                dictionary dummyOptimDict;
                dummyOptimDict.add("type","foamOptimize",true);
                dummyOptimDict.add
                    ("foamOptimizeCoeffs",dictionary(), true);
                dummyOptimDict.subDict("foamOptimizeCoeffs").add
                (
                    "optimizationMethod",
                    "fullMeshOptimization",
                    true
                );
                dummyOptimDict.subDict("foamOptimizeCoeffs").add
                (
                    "anisotropy",
                    "isotropy",
                    true
                );
                dummyOptimDict.subDict("foamOptimizeCoeffs").add
                (
                    "iterations",
                    "5",
                    true
                );
                autoPtr<autoOptimize> optimMeshPtr
                    = autoOptimize::New(mesh, dummyOptimDict);
                optimMeshPtr->movePoints(newPoints);
            }

            //optimisation step of cells in error
            if
            (
                (opt & optStages::ERROR)
                || (opt & optStages::RELAXEDERROR)
                || (opt & optStages::BOUNDARYRELAXEDERROR)
            )
            {
                //Optimise the mesh
                dictionary meshOptimDict
                (
                    meshDict_.found("meshOptimization") ?
                    meshDict_.subDict("meshOptimization") :
                    dictionary()
                );

                word optimType =
                    meshOptimDict.lookupOrDefault<word>
                    (
                        "type",
                        "cfMeshOptimize"
                    );

                if (optimType == "foamOptimize")
                {
                    meshOptimDict.subDict("foamOptimizeCoeffs").add
                    (
                        "optimizationMethod",
                        "fullMeshOptimization",
                        true
                    );
                }
                else if
                (
                    optimType == "cfMeshOptimize"
                )
                {
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
                    if (opt & optStages::RELAXEDERROR)
                    {
                        coeffsDict.add
                        (
                            "relaxedCheck",
                            true,
                            true
                        );
                    }
                    else if (opt & optStages::BOUNDARYRELAXEDERROR)
                    {
                        coeffsDict.add
                        (
                            "relaxedBoundaryCheck",
                            true,
                            true
                        );
                    }
                }
                autoPtr<autoOptimize> optimMeshPtr
                    = autoOptimize::New(mesh, meshOptimDict);
                optimMeshPtr->movePoints(newPoints);
            }

            //Update internal field values
            pointVectorField& disp = meshMover.displacement();
            disp.primitiveFieldRef() = newPoints - meshMover.oldPoints();
        }
        else if
        (
            mesh.foundObject<volScalarField>("gapCells") && (iter % 10 == 0)
        )
        {
            const volScalarField& gapCells =
                mesh.lookupObject<volScalarField>("gapCells");

            boolList internalFaces(mesh.nFaces(), false);
            label nInternalFaces = 0;
            forAll(gapCells, celli)
            {
                if (gapCells[celli] > -1)
                {
                    const cell& cFaces = mesh.cells()[celli];
                    forAll(cFaces, cFI)
                    {
                        label facei = cFaces[cFI];
                        const face& f = mesh.faces()[facei];
                        bool foundBoundaryPt = false;
                        forAll(f, fp)
                        {
                            if (gapBoundaryPts[f[fp]])
                            {
                                foundBoundaryPt = true;
                                break;
                            }
                        }
                        if (!foundBoundaryPt)
                        {
                            internalFaces[facei] = true;
                            nInternalFaces++;
                        }
                    }
                }
            }

            syncTools::syncFaceList(mesh, internalFaces, orEqOp<bool>());
            DynamicList<label> testFaces(nInternalFaces);
            forAll(internalFaces, facei)
            {
                if (internalFaces[facei])
                {
                    testFaces.append(facei);
                }
            }

            const indirectPrimitivePatch internalPatch
            (
                IndirectList<face>
                (
                    mesh.faces(),
                    testFaces
                ),
                mesh.points()
             );

            const labelList& internalPts = internalPatch.meshPoints();
            pointField midPoint(internalPts.size(),vector::zero);
            scalarField midPointWeight(internalPts.size(),Zero);
            forAll(internalPts, ptI)
            {
                label meshPointI = internalPts[ptI];
                const labelList& pEdges = mesh.pointEdges()[meshPointI];
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    edge e = mesh.edges()[edgei];
                    label otherPt(meshPointI == e[0] ? e[1] : e[0]);
                    if (gapBoundaryPts[otherPt])
                    {
                        midPoint[ptI] += newPoints[otherPt];
                        midPointWeight[ptI] += scalar(1.0);
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                internalPts,
                midPoint,
                plusEqOp<point>(),
                vector::zero         // null value
            );

            syncTools::syncPointList
            (
                mesh,
                internalPts,
                midPointWeight,
                plusEqOp<scalar>(),
                scalar(0)      // null value
            );

            forAll(internalPts, ptI)
            {
                label meshPointI = internalPts[ptI];
                if (midPointWeight[ptI] > SMALL)
                {
                    newPoints[meshPointI] = midPoint[ptI]
                        / midPointWeight[ptI];
                }
            }

            boolList exclFaces(pp.size(), false);
            edgeClassification eClass
            (
                mesh,
                newPoints,
                pp,
                meshEdges,
                exclFaces,
                0.8191,
                0.8191
            );
            const List<Tuple2<edgeClassification::edgeType,scalar>>&
                eType = eClass.edgeTypes();

            pointField pointNormals =
                eClass.calculatePointNormals(exclFaces, 8, false);

            boolList featurePts(meshPoints.size(), false);
            forAll(pp.edges(), patchEdgeI)
            {
                const edge& e = pp.edges()[patchEdgeI];
                label v0 = e[0];
                label v1 = e[1];

                if
                (
                    eType[patchEdgeI].first() == edgeClassification::CONCAVE
                    || eType[patchEdgeI].first() == edgeClassification::CONVEX
                 )
                {
                    featurePts[v0] = true;
                    featurePts[v1] = true;
                }
            }

            forAll(meshPoints, ptI)
            {
                if (featurePts[ptI])
                {
                    label meshPointI = meshPoints[ptI];
                    const labelList& pEdges =
                        mesh.pointEdges()[meshPointI];
                    forAll(pEdges, pEI)
                    {
                        const edge e = mesh.edges()[pEdges[pEI]];
                        label otherPt(e[0] == meshPointI ? e[1] : e[0]);
                        if (!boundaryPoints[otherPt])
                        {
                            const scalar eLen = e.mag(newPoints);
                            newPoints[otherPt] = newPoints[meshPointI]
                                - eLen*pointNormals[ptI];
                        }
                    }
                }
            }

            DynamicList<label> gapSnapPts(meshPoints.size());
            DynamicList<point> gapSnapLoc(meshPoints.size());
            forAll(meshPoints, ptI)
            {
                label meshPointI = meshPoints[ptI];
                if (!featurePts[ptI] && gapBoundaryPts[meshPointI])
                {
                    const labelList& pEdges =
                        mesh.pointEdges()[meshPointI];
                    forAll(pEdges, pEI)
                    {
                        const edge e = mesh.edges()[pEdges[pEI]];
                        label otherPt(e[0] == meshPointI ? e[1] : e[0]);
                        if (!boundaryPoints[otherPt])
                        {
                            plane pl
                            (
                                newPoints[meshPointI],
                                pointNormals[ptI]
                            );
                            point snPt =
                                pl.nearestPoint(newPoints[otherPt]);
                            gapSnapPts.append(meshPointI);
                            gapSnapLoc.append(snPt);
                            break;
                        }
                    }
                }
            }

            pointField appPts(gapSnapLoc.xfer());
            List<pointIndexHit> gapHitInfo;
            labelList gapHitSurface;

            surfaces.findNearest
            (
                snapSurfaces,
                appPts,
                scalarField(gapSnapPts.size(), GREAT),
                gapHitSurface,
                gapHitInfo
             );

            forAll(gapSnapPts, ptI)
            {
                if (gapHitInfo[ptI].hit())
                {
                    newPoints[gapSnapPts[ptI]] =
                        gapHitInfo[ptI].hitPoint();
                }
            }

            //first global optimisation step using sphericity
            {
                dictionary dummyOptimDict;
                dummyOptimDict.add("type","foamOptimize",true);
                dummyOptimDict.add
                    ("foamOptimizeCoeffs",dictionary(), true);
                dummyOptimDict.subDict("foamOptimizeCoeffs").add
                (
                    "optimizationMethod",
                    "fullMeshOptimization",
                    true
                 );
                dummyOptimDict.subDict("foamOptimizeCoeffs").add
                (
                    "anisotropy",
                    "isotropy",
                    true
                );
                dummyOptimDict.subDict("foamOptimizeCoeffs").add
                (
                    "iterations",
                    "5",
                    true
                );
                autoPtr<autoOptimize> optimMeshPtr
                    = autoOptimize::New(mesh, dummyOptimDict);
                optimMeshPtr->movePoints(newPoints);
            }

            //second optimisation step of cells in error
            {
                //Optimise the mesh
                dictionary meshOptimDict
                (
                    meshDict_.found("meshOptimization") ?
                    meshDict_.subDict("meshOptimization") :
                    dictionary()
                );

                word optimType =
                    meshOptimDict.lookupOrDefault<word>("type","cfMesh");
                if (optimType == "foamOptimize")
                {
                    meshOptimDict.subDict("foamOptimizeCoeffs").add
                    (
                        "optimizationMethod",
                        "fullMeshOptimization",
                        true
                     );
                }
                autoPtr<autoOptimize> optimMeshPtr
                    = autoOptimize::New(mesh, meshOptimDict);
                optimMeshPtr->movePoints(newPoints);
            }
        }

        //Update internal field values
        pointVectorField& disp = meshMover.displacement();
        disp.primitiveFieldRef() = newPoints - meshMover.oldPoints();
    }
    else if
    (
        controller_.mode() == meshControl::QUALITY
        && snapParams.featureSnapChecks() && iter > 0 && iter % 10 == 0
    )
    {
        pointField newDisp( newPoints - meshMover.oldPoints() );
        syncTools::syncPointList
        (
            mesh,
            newDisp,
            minMagSqrEqOp<point>(),
            vector(GREAT, GREAT, GREAT)
        );

        newDisp += meshMover.oldPoints();
        mesh.movePoints(newDisp);

        dictionary updatedMotionDict(motionDict);
        updatedMotionDict.add("minTwist",-1,true);

        labelHashSet wrongFaces(mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh,updatedMotionDict,wrongFaces);
        forAll(meshPoints, i)
        {
            label meshPointI = meshPoints[i];
            const labelList& pointFaces = mesh.pointFaces()[meshPointI];

            forAll(pointFaces, faceI)
            {
                if (wrongFaces.found(pointFaces[faceI]))
                {
                    if (!stationaryPoints[meshPointI])
                    {
                        stationaryPoints[meshPointI] = true;
                    }
                }
            }
        }

        mesh.movePoints(meshMover.oldPoints());
    }

    syncTools::syncPointList
    (
        mesh,
        stationaryPoints,
        orEqOp<bool>(),
        false               // null value
    );

    return;
}



void Foam::snappySnapDriver::tightAngleCorrect
(
    const bool weakSnap,
    const scalar weakAngleMin,
    const scalar weakAngleMax,
    const scalar tightAngle,
    const scalar minTightAngle,
    const pointField& newPoints,
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const List<Tuple2<edgeClassification::edgeType,scalar>>& eType,
    const List<pointIndexHit>& finalHitInfo,
    const labelList& pointIndex,
    const vectorField& ppFaceAreas,
    const vectorField& pNormals,
    scalarField& sumAreasPP,
    vectorField& sumDispPP
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    // Undistorted edge length
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    labelList reversePtMap(pp.size(), -1);
    forAll(finalHitInfo, fI)
    {
        label i = pointIndex[fI];
        if (finalHitInfo[fI].hit())
        {
            reversePtMap[i]= fI;
        }
    }

    typedef std::tuple<pointField, pointField, pointField> eData;
    eData eDataNullVal = eData(pointField(0),pointField(0),pointField(0));

    List<eData> edgeData(pp.edges().size(),eDataNullVal);
    forAll(pp.edges(), edgei)
    {
        const labelList& eFaces = pp.edgeFaces()[edgei];
        pointField bVecs(eFaces.size());
        pointField sVecs(eFaces.size());
        pointField sPts(eFaces.size());
        label nPts = 0;
        forAll(eFaces, efi)
        {
            label facei = eFaces[efi];
            label indexi = reversePtMap[facei];
            if (indexi != -1)
            {
                bVecs[nPts] = ppFaceAreas[facei];
                sVecs[nPts] = pNormals[indexi];
                sPts[nPts] = finalHitInfo[indexi].hitPoint();
                nPts++;
            }
        }
        bVecs.setSize(nPts);
        sVecs.setSize(nPts);
        sPts.setSize(nPts);

        edgeData[edgei] = std::make_tuple
        (
            std::move(bVecs),
            std::move(sVecs),
            std::move(sPts)
        );
    }

    syncTools::syncEdgeList
    (
        mesh,
        meshEdges,
        edgeData,
        packedEltwiseAggOp<typename edgeClassification::sumEqOp>{},
        eDataNullVal,
        dummyTransform()
    );

    forAll(pp.edges(), edgei)
    {
        const auto& [eN, sN, sP] = edgeData[edgei];
        if (eN.size() == 2)
        {
            const vector fArea0 = eN[0];
            const vector fArea1 = eN[1];

            const vector fNorm0 = fArea0/(mag(fArea0)+ SMALL);
            const vector fNorm1 = fArea1/(mag(fArea1)+ SMALL);
            scalar dProd = (fNorm0 & fNorm1);

            vector sNorm0 = sN[0];
            vector sNorm1 = sN[1];

            scalar magSNorm0 = mag(sNorm0);
            scalar magSNorm1 = mag(sNorm1);

            if (magSNorm0 < SMALL || magSNorm1 < SMALL)
            {
                continue;
            }
            sNorm0 /= magSNorm0;
            sNorm1 /= magSNorm1;

            scalar sProd  = mag(sNorm0&sNorm1);
            bool baffleEdge = false;
            if
            (
                eType[edgei].first() == edgeClassification::BAFFLE
                && eType[edgei].second() < minTightAngle
            )
            {
                baffleEdge = true;
            }

            bool weakSnapCheck = false;
            bool tightAngleCheck = false;
            if
            (
                weakSnap
                && dProd < weakAngleMin && dProd > weakAngleMax
                && sProd < weakAngleMin && sProd > weakAngleMax
            )
            {
                weakSnapCheck = true;
            }

            if
            (
                dProd < tightAngle
                && (dProd > minTightAngle || baffleEdge)
            )
            {
                tightAngleCheck = true;
            }

            if (weakSnapCheck || tightAngleCheck)
            {
                if
                (
                    tightAngleCheck &&
                    (
                        (!baffleEdge && sProd > -minTightAngle)
                        || (baffleEdge && sProd > 0.996)
                    )
                )
                {
                    continue;
                }

                plane pl0(sP[0],sNorm0);
                plane pl1(sP[1],sNorm1);

                plane::ray r(pl0.planeIntersect(pl1));
                vector n = r.dir() / mag(r.dir());

                edge e = pp.edges()[edgei];
                forAll(e, ei)
                {
                    label meshpointi = pp.meshPoints()[e[ei]];
                    point pt = newPoints[meshpointi];
                    // Get nearest point on infinite ray
                    vector d = r.refPoint()-pt;
                    d -= (d&n)*n;

                    scalar edgeLen = edge0Len /
                        (1<<pointLevel[meshpointi]);
                    scalar magD = mag(d);

                    scalar limitDisp  =  0.5*edgeLen;
                    scalar maxDisp  =  4*edgeLen;
                    if (magD < maxDisp)
                    {
                        if (magD > limitDisp)
                        {
                            d *= (limitDisp/magD);
                        }
                        scalar weight =
                            0.5*(mag(fArea0)+mag(fArea1));
                        sumDispPP[e[ei]] += d / (weight + SMALL);
                        sumAreasPP[e[ei]] += 1.0 / (weight + SMALL);
                    }
                }
            }
        }
    }

    return;
}


void Foam::snappySnapDriver::recalculatePointNormals
(
    const snapParameters& snapParams,
    const indirectPrimitivePatch& pp,
    const labelList& snapSurfaces,
    const scalarField& snapDist,
    const pointField& newPoints,
    const List<pointIndexHit>& finalHitInfo,
    const labelList& pointIndex,
    const pointField& ppFaceCentres,
    const vectorField& ppFaceAreas,
    vectorField& pNormals
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const labelList& meshPoints = pp.meshPoints();

    label sz = 0;
    forAll(finalHitInfo, fI)
    {
        if (!finalHitInfo[fI].hit())
        {
            continue;
        }
        label i = pointIndex[fI];
        const face& f = pp.localFaces()[i];
        label meshfacei = pp.addressing()[i];
        label own = mesh.faceOwner()[meshfacei];
        label cLevel = cellLevel[own];
        label nAnchors = 0;
        forAll(f,fp)
        {
            label pLevel = pointLevel[meshPoints[f[fp]]];
            if (pLevel <= cLevel)
            {
                nAnchors++;
            }
        }
        sz += nAnchors;
        if (snapParams.enlargeStencil())
        {
            sz += nAnchors;
        }
    }

    pointField nbrPts(sz);
    scalarField nbrSnapDist(sz);

    sz = 0;
    forAll(finalHitInfo, fI)
    {
        if (!finalHitInfo[fI].hit())
        {
            continue;
        }
        label i = pointIndex[fI];
        const face& f = pp.localFaces()[i];
        const point& fc = ppFaceCentres[i];
        label meshfacei = pp.addressing()[i];
        label own = mesh.faceOwner()[meshfacei];
        label cLevel = cellLevel[own];
        DynamicList<label> anchorPts(f.size());

        forAll(f, fp)
        {
            label pointi = f[fp];
            label meshpointi = meshPoints[pointi];
            label pLevel = pointLevel[meshpointi];
            if (pLevel <= cLevel)
            {
                anchorPts.append(pointi);
                point mid =
                    fc + 0.5*(newPoints[meshpointi] -fc);
                nbrPts[sz] = mid;
                nbrSnapDist[sz] = snapDist[pointi];
                sz++;
            }
        }

        if (snapParams.enlargeStencil())
        {
            label start = 0;
            face fAnchors(anchorPts);
            forAll(fAnchors, fp)
            {
                label next = fAnchors.fcIndex(start);
                point edgeMid = 0.5*
                (
                    newPoints[meshPoints[fAnchors[start]]]
                    +newPoints[meshPoints[fAnchors[next]]]
                 );

                point mid = fc + 0.5*(edgeMid -fc);

                nbrPts[sz] = mid;
                nbrSnapDist[sz] = snapDist[fAnchors[fp]];
                start = fAnchors.fcIndex(start);
                sz++;
            }
        }
    }

    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    surfaces.findNearest
    (
        snapSurfaces,
        nbrPts,
        sqr(nbrSnapDist),        // sqr of attract distance
        hitSurface,
        hitInfo
    );

    pointField pNormalsNbr(nbrPts.size(), vector(GREAT, GREAT, GREAT));
    forAll(snapSurfaces, sI)
    {
        DynamicList<pointIndexHit> localHits;
        label surfI = snapSurfaces[sI];
        forAll(hitSurface, i)
        {
            if (hitSurface[i] == surfI)
            {
                localHits.append(hitInfo[i]);
            }
        }
        localHits.shrink();

        pointField localNormals;
        label geomI = surfaces.surfaces()[surfI];
        surfaces.geometry()[geomI].getNormal
        (
            localHits,
            localNormals
        );
        label localI = 0;
        forAll(hitSurface, i)
        {
            if (hitSurface[i] == surfI)
            {
                pNormalsNbr[i] = localNormals[localI];
                localI++;
            }
        }
    }
    sz = 0;
    forAll(finalHitInfo, fI)
    {
        if (!finalHitInfo[fI].hit())
        {
            continue;
        }
        label i = pointIndex[fI];
        const face& f = pp.localFaces()[i];

        label meshfacei = pp.addressing()[i];
        label own = mesh.faceOwner()[meshfacei];
        label cLevel = cellLevel[own];
        DynamicList<label> anchorPts(f.size());
        forAll(f, fp)
        {
            label pointi = f[fp];
            label meshpointi = meshPoints[pointi];
            label pLevel = pointLevel[meshpointi];
            if (pLevel <= cLevel)
            {
                anchorPts.append(pointi);
            }
        }
        face fAnchors(anchorPts);

        point aveNorm = pNormals[fI] / (mag(pNormals[fI]) + SMALL);
        label nHit = 0;
        vector fArea = ppFaceAreas[i] / (mag(ppFaceAreas[i]) + SMALL);

        scalar dotProd = mag(fArea & aveNorm);
        vector maxNorm = pNormals[fI] / (mag(pNormals[fI]) + SMALL);
        scalar aveDot = dotProd;
        scalar maxDot = dotProd;

        forAll(fAnchors, fp)
        {
            if (hitInfo[sz].hit())
            {
                vector nNorm = pNormalsNbr[sz]
                    / (mag(pNormalsNbr[sz]) + SMALL);

                if ((nNorm & pNormals[fI]) > 0.)
                {
                    aveNorm +=  nNorm;
                }
                else
                {
                    aveNorm -=  nNorm;
                }
                nHit++;

                dotProd = mag(fArea & nNorm);
                if (dotProd > maxDot)
                {
                    maxDot = dotProd;
                    maxNorm = nNorm;
                }
                aveDot += dotProd;
            }
            sz++;
        }
        if (snapParams.enlargeStencil())
        {
            forAll(fAnchors, fp)
            {
                if (hitInfo[sz].hit())
                {
                    vector nNorm = pNormalsNbr[sz]
                        / (mag(pNormalsNbr[sz]) + SMALL);

                    if ((nNorm & pNormals[fI]) > 0.)
                    {
                        aveNorm +=  nNorm;
                    }
                    else
                    {
                        aveNorm -=  nNorm;
                    }
                    nHit++;
                    dotProd = mag(fArea & nNorm);
                    if (dotProd > maxDot)
                    {
                        maxDot = dotProd;
                        maxNorm = nNorm;
                    }
                    aveDot += dotProd;
                }
                sz++;
            }
        }

        aveDot /= (nHit + 1);
        if (aveDot > 0.9848)
        {
            pNormals[fI] = aveNorm / (nHit + 1) ;
        }
        else
        {
            pNormals[fI] = maxNorm;
        }
    }

    return;
}


void Foam::snappySnapDriver::featureSnapSurfaceSmooth
(
    const label iter,
    const label nSurfSmooth,
    const indirectPrimitivePatch& pp,
    const List<Tuple2<edgeClassification::edgeType,scalar>>& eType,
    const PackedBoolList& isMasterEdge,
    const PackedBoolList& isMasterFace,
    const PackedBoolList& isMasterPPEdge,
    const scalarField& snapDist,
    const boolList& boundaryPoints,
    const scalarField& smoothingWeights,
    const boolList& excludeFaces,
    const boolList& zonedEdges,
    pointField& newPoints,
    vectorField& pointDisp,
    pointField& ppFaceCentres,
    vectorField& ppFaceAreas,
    boolList& stationaryPoints,
    scalarField& sumAreasMesh,
    vectorField& sumDispMesh,
    scalarField& sumAreasPP,
    vectorField& sumDispPP
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();

    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));
    const labelList& meshPoints = pp.meshPoints();

    // Update cell centre, face centre and face areas
    forAll(pp, i)
    {
        const label& faceI = pp.addressing()[i];
        ppFaceCentres[i] =
            mesh.faces()[faceI].centre(newPoints);
        ppFaceAreas[i] =
            mesh.faces()[faceI].areaNormal(newPoints);
    }

    boolList featureEdge(mesh.nEdges(), false);
    labelList featurePoint(mesh.nPoints(), -2);
    labelList nFeatureEdges(mesh.nPoints(), 0);

    stationaryPoints = false;

    forAll(pp.edgeFaces(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        const edge e = mesh.edges()[meshEdgeI];
        if
        (
            eType[edgeI].first() == edgeClassification::CONCAVE
        )
        {
            featureEdge[meshEdgeI] = true;
            featurePoint[e[0]] = max(label(1),featurePoint[e[0]]);
            featurePoint[e[1]] = max(label(1),featurePoint[e[1]]);
        }
        else if
        (
            eType[edgeI].first() == edgeClassification::CONVEX
        )
        {
            featureEdge[meshEdgeI] = true;
            featurePoint[e[0]] = max(label(0),featurePoint[e[0]]);
            featurePoint[e[1]] = max(label(0),featurePoint[e[1]]);
        }
        else if
        (
            eType[edgeI].first() == edgeClassification::BAFFLE
        )
        {
            if (eType[edgeI].second() < -0.9848)
            {
                featureEdge[meshEdgeI] = true;
                featurePoint[e[0]] = max(label(2),featurePoint[e[0]]);
                featurePoint[e[1]] = max(label(2),featurePoint[e[1]]);
            }
            else
            {
                featureEdge[meshEdgeI] = true;
                featurePoint[e[0]] = max(label(0),featurePoint[e[0]]);
                featurePoint[e[1]] = max(label(0),featurePoint[e[1]]);
            }
        }
        else if
        (
            eType[edgeI].first() == edgeClassification::NONMANIFOLD
        )
        {
            featureEdge[meshEdgeI] = true;
            featurePoint[e[0]] = max(label(-1),featurePoint[e[0]]);
            featurePoint[e[1]] = max(label(-1),featurePoint[e[1]]);
        }
        else if
        (
            eType[edgeI].first() == edgeClassification::BOUNDARY
            && zonedEdges[edgeI]
        )
        {
            featureEdge[meshEdgeI] = true;
            featurePoint[e[0]] = max(label(2),featurePoint[e[0]]);
            featurePoint[e[1]] = max(label(2),featurePoint[e[1]]);
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        featureEdge,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh,
        featurePoint,
        maxEqOp<label>(),
        label(-1)
    );

    forAll(pp.edgeFaces(), edgeI)
    {
        if
        (
            eType[edgeI].first() == edgeClassification::NONMANIFOLD
            ||
            (
               eType[edgeI].first() == edgeClassification::BOUNDARY
               && !zonedEdges[edgeI]
            )
        )
        {
            label meshEdgeI = meshEdges[edgeI];
            const edge e = mesh.edges()[meshEdgeI];
            stationaryPoints[e[0]] = true;
            stationaryPoints[e[1]] = true;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        stationaryPoints,
        orEqOp<bool>(),
        false               // null value
    );

    nFeatureEdges = 0;
    forAll(mesh.points(), pointI)
    {
        labelList pEdges = mesh.pointEdges()[pointI];
        forAll(pEdges, pEI)
        {
            if (isMasterEdge[pEdges[pEI]] && featureEdge[pEdges[pEI]])
            {
                nFeatureEdges[pointI]++;
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        nFeatureEdges,
        plusEqOp<label>(),
        label(0)
    );

    //Detect feature points to be deactivated for smoothing
    forAll(pp, i)
    {
        const label& faceI = pp.addressing()[i];

        face f = mesh.faces()[faceI];
        label fsz = f.size();

        label nFtrPts = 0;
        label nBafflePts = 0;
        label nNonMan = 0;
        label nCornerPts = 0;
        label nSingleFtrPts = 0;
        forAll(f,fp)
        {
            label pointi = f[fp];
            label ftrPtType = featurePoint[f[fp]];
            if (ftrPtType > -1)
            {
                if (ftrPtType == 2)
                {
                    nBafflePts++;
                }
                nFtrPts++;
            }
            else if (ftrPtType == -1)
            {
                nNonMan++;
            }

            if (nFeatureEdges[pointi] > 2)
            {
                nCornerPts++;
            }
            else if (nFeatureEdges[pointi] == 1)
            {
                nSingleFtrPts++;
            }
        }

        const labelList& fEdges = mesh.faceEdges()[faceI];
        label nFtrEdges = 0;
        forAll(fEdges, fEI)
        {
            label edgei = fEdges[fEI];
            if (featureEdge[edgei])
            {
                nFtrEdges++;
            }
        }
        bool resetFtrFacePts = false;
        bool featureTermination = false;
        if (nFtrPts > 0 && nNonMan > 0)
        {
            resetFtrFacePts = true;
        }
        else if (nFtrPts == fsz)
        {
            if (nFtrEdges == fEdges.size() && nCornerPts != fsz)
            {
                resetFtrFacePts = true;
            }
            else if (nSingleFtrPts > 0 && nFtrEdges == fEdges.size()-1)
            {
                featureTermination = true;
            }
        }
        else if (nFtrPts == 2)
        {
            forAll(fEdges, fEI)
            {
                label edgeI = fEdges[fEI];
                if (featureEdge[edgeI])
                {
                    edge e = mesh.edges()[edgeI];
                    if
                    (
                        nFeatureEdges[e[0]] == 1
                        && nFeatureEdges[e[1]] == 1
                    )
                    {
                        resetFtrFacePts = true;
                        break;
                    }
                }
            }
        }

        if (featureTermination)
        {
            forAll(f,fp)
            {
                label pointi = f[fp];
                label ftrPtType = featurePoint[pointi];
                if (ftrPtType > -1  && nFeatureEdges[pointi] == 1)
                {
                    featurePoint[pointi] = -2;
                }
            }
        }
        else if (resetFtrFacePts)
        {
            bool allBaffle(fsz == nBafflePts ? true : false);
            forAll(f,fp)
            {
                label pointi = f[fp];
                label ftrPtType = featurePoint[pointi];
                if (ftrPtType > -1)
                {
                    if (allBaffle)
                    {
                        featurePoint[pointi] = -2;
                    }
                    else if (nFeatureEdges[pointi] == 1)
                    {
                        featurePoint[pointi] = -2;
                    }
                    else if (ftrPtType != 2 && nNonMan == 0)
                    {
                        featurePoint[pointi] = -2;
                    }
                }
            }
        }
        else if (f.size() == 3 && nFtrEdges == 2)
        {
            forAll(f, fp)
            {
                label pointi = f[fp];
                if (nFeatureEdges[pointi] == 2)
                {
                    point currPt = newPoints[pointi];
                    label nextFp = f.fcIndex(fp);
                    label prevFp = f.rcIndex(fp);

                    vector nextEdge = newPoints[f[nextFp]]
                       - currPt;
                    nextEdge /= (mag(nextEdge) + SMALL);
                    vector prevEdge = currPt
                       - newPoints[f[prevFp]];
                    prevEdge /= (mag(prevEdge) + SMALL);

                    if ((nextEdge & prevEdge) > 0.939)
                    {
                        featurePoint[pointi] = -2;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        featurePoint,
        minEqOp<label>(),
        label(-1)
    );

    sumDispPP = vector::zero;
    sumAreasPP = 0.0;

    boolList switchPts(mesh.nPoints(), false);
    if (controller_.algorithm() == meshControl::EXTRUDE && iter > 20)
    {
        pointField pEdgeCentre(meshPoints.size(), vector::zero);
        labelList nPEdges(meshPoints.size(), 0);

        forAll(meshPoints, pointi)
        {
            const labelList& pEdges = pp.pointEdges()[pointi];
            forAll(pEdges, pei)
            {
                label edgei = pEdges[pei];
                label meshedgei = meshEdges[edgei];
                if (!featureEdge[meshedgei])
                {
                    pEdgeCentre[pointi] +=
                        mesh.edges()[meshedgei].centre(newPoints);
                    nPEdges[pointi]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPEdges,
            plusEqOp<label>(),
            label(0) // null value
        );
        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            pEdgeCentre,
            plusEqOp<vector>(),
            vector::zero // null value
        );

        labelList nPPPointFaces(meshPoints.size(), 0);
        forAll(meshPoints, pointi)
        {
            nPPPointFaces[pointi] += pp.pointFaces()[pointi].size();
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPPPointFaces,
            plusEqOp<label>(),
            label(0) // null value
        );

        labelList nEdgeFlipFaces(mesh.nEdges(), 0);
        boolList markedForFlip(pp.size(), false);
        forAll(pp, i)
        {
            const face& lf = pp.localFaces()[i];

            if (lf.size() != 4)
            {
                continue;
            }

            forAll(lf, lfi)
            {
                label pointi = lf[lfi];
                label meshpointi = meshPoints[pointi];
                if
                (
                    featurePoint[meshpointi] == 1
                    && nFeatureEdges[meshpointi] == 2
                    && nPPPointFaces[pointi] == 4
                )
                {
                    label nextFp = lf.fcIndex(lfi);
                    label oppFp = lf.fcIndex(nextFp);
                    label prevFp = lf.rcIndex(lfi);
                    label oppmeshpointi = meshPoints[lf[oppFp]];
                    if
                    (
                        featurePoint[oppmeshpointi] == -2
                        && nPPPointFaces[lf[oppFp]] == 3
                     )
                    {
                        label nextPt = lf[nextFp];
                        label prevPt = lf[prevFp];
                        label meshpointnext = meshPoints[nextPt];
                        label meshpointprev = meshPoints[prevPt];
                        vector nextEdge = newPoints[meshpointnext]
                            - newPoints[meshpointi];
                        nextEdge /= (mag(nextEdge) + SMALL);

                        vector prevEdge = newPoints[meshpointi]
                            - newPoints[meshpointprev];
                        prevEdge /= (mag(prevEdge) + SMALL);

                        if ((nextEdge & prevEdge) > 0.939)
                        {
                            markedForFlip[i] = true;
                        }
                    }
                }
            }
        }

        forAll(pp.edges(), edgei)
        {
            label meshedgei = meshEdges[edgei];
            const labelList& eFaces = pp.edgeFaces()[edgei];
            forAll(eFaces, eFI)
            {
                if (markedForFlip[eFaces[eFI]])
                {
                    nEdgeFlipFaces[meshedgei]++;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            nEdgeFlipFaces,
            plusEqOp<label>(),
            label(0)
        );

        forAll(pp, i)
        {
            const face& lf = pp.localFaces()[i];
            if (lf.size() != 4 || !markedForFlip[i])
            {
                continue;
            }
            forAll(lf, lfi)
            {
                label pointi = lf[lfi];
                label meshpointi = meshPoints[pointi];
                if
                (
                    featurePoint[meshpointi] == 1
                    && nFeatureEdges[meshpointi] == 2
                    && nPPPointFaces[pointi] == 4
                )
                {
                    label nextFp = lf.fcIndex(lfi);
                    label oppFp = lf.fcIndex(nextFp);
                    label prevFp = lf.rcIndex(lfi);
                    label oppmeshpointi = meshPoints[lf[oppFp]];
                    if
                    (
                        featurePoint[oppmeshpointi] == -2
                        && nPPPointFaces[lf[oppFp]] == 3
                    )
                    {
                        label nextPt = lf[nextFp];
                        label prevPt = lf[prevFp];
                        label nextEdge = meshTools::findEdge
                        (
                            pp.edges(),
                            pp.pointEdges()[pointi],
                            pointi,
                            nextPt
                        );

                        label prevEdge = meshTools::findEdge
                        (
                            pp.edges(),
                            pp.pointEdges()[pointi],
                            pointi,
                            prevPt
                        );

                        if
                        (
                            nEdgeFlipFaces[meshEdges[nextEdge]] > 1
                            || nEdgeFlipFaces[meshEdges[prevEdge]] > 1
                        )
                        {
                            continue;
                        }

                        sumDispPP[lf[oppFp]] += newPoints[meshpointi];
                        sumAreasPP[lf[oppFp]] += scalar(1);

                        if (nPEdges[pointi] > 0)
                        {
                            sumDispPP[pointi] +=
                                pEdgeCentre[pointi]/nPEdges[pointi];
                            sumAreasPP[pointi] += scalar(1);
                        }
                        switchPts[meshpointi] = true;
                        switchPts[oppmeshpointi] = true;
                        break;
                    }
                }
            }
        }

        forAll(pp, i)
        {
            const face& lf = pp.localFaces()[i];
            if (lf.size() != 4 || markedForFlip[i])
            {
                continue;
            }

            forAll(lf, lfi)
            {
                label pointi = lf[lfi];
                label meshpointi = meshPoints[pointi];

                if
                (
                    nFeatureEdges[meshpointi] == 1
                    || nFeatureEdges[meshpointi] >= 3
                )
                {
                    label nextFp = lf.fcIndex(lfi);
                    label prevFp = lf.rcIndex(lfi);

                    label nextMeshPt = meshPoints[lf[nextFp]];
                    label prevMeshPt = meshPoints[lf[prevFp]];

                    label nextEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[meshpointi],
                        meshpointi,
                        nextMeshPt
                    );
                    label prevEdge = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[meshpointi],
                        meshpointi,
                        prevMeshPt
                    );

                    if (prevEdge != -1  && nextEdge != -1)
                    {
                        bool nextFE = featureEdge[nextEdge];
                        bool prevFE = featureEdge[prevEdge];
                        bool valid = false;
                        if
                        (
                            nFeatureEdges[meshpointi] == 1
                            && !nextFE && !prevFE
                        )
                        {
                            valid = true;
                        }
                        else if
                        (
                            nFeatureEdges[meshpointi] == 3
                            && nextFE && prevFE
                        )
                        {
                            valid = true;
                        }

                        if (valid)
                        {
                            label oppFp = lf.fcIndex(nextFp);
                            label oppmeshpointi =
                                meshPoints[lf[oppFp]];
                            if
                            (
                                nFeatureEdges[oppmeshpointi] == 1
                                || nFeatureEdges[oppmeshpointi] >= 3
                            )
                            {
                                nextEdge = meshTools::findEdge
                                (
                                    mesh.edges(),
                                    mesh.pointEdges()[oppmeshpointi],
                                    oppmeshpointi,
                                    nextMeshPt
                                );

                                prevEdge = meshTools::findEdge
                                (
                                    mesh.edges(),
                                    mesh.pointEdges()[oppmeshpointi],
                                    oppmeshpointi,
                                    prevMeshPt
                                );

                                if (prevEdge != -1  && nextEdge != -1)
                                {
                                    nextFE = featureEdge[nextEdge];
                                    prevFE = featureEdge[prevEdge];
                                    valid = false;
                                    if
                                    (
                                        nFeatureEdges[oppmeshpointi] == 1
                                        && !nextFE && !prevFE
                                    )
                                    {
                                        valid = true;
                                    }
                                    else if
                                    (
                                        nFeatureEdges[oppmeshpointi] == 3
                                        && nextFE && prevFE
                                    )
                                    {
                                        valid = true;
                                    }

                                    if (valid)
                                    {
                                        sumDispPP[lf[nextFp]] += 0.5*
                                        (
                                            newPoints[oppmeshpointi]
                                            + newPoints[meshpointi]
                                        );
                                        sumAreasPP[lf[nextFp]] +=
                                            scalar(1);
                                        switchPts[meshPoints[lf[nextFp]]]
                                            = true;
                                        break;
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
            mesh,
            switchPts,
            orEqOp<bool>(),
            false
        );

        markedForFlip = false;
        forAll(pp, i)
        {
            const face& lf = pp.localFaces()[i];
            if (lf.size() != 4)
            {
                continue;
            }

            bool foundSwitch = false;
            bool flipFace = false;
            forAll(lf, lfi)
            {
                label pointi = lf[lfi];
                label meshpointi = meshPoints[pointi];
                if (switchPts[meshpointi])
                {
                    foundSwitch = true;
                    break;
                }

                if
                (
                    featurePoint[meshpointi] == 1
                    && nFeatureEdges[meshpointi] == 2
                    && nPPPointFaces[pointi] == 4
                 )
                {
                    label nextFp = lf.fcIndex(lfi);
                    label oppFp = lf.fcIndex(nextFp);
                    label prevFp = lf.rcIndex(lfi);
                    label oppmeshpointi = meshPoints[lf[oppFp]];
                    if
                    (
                        featurePoint[oppmeshpointi] == -2
                        && nPPPointFaces[lf[oppFp]] == 4
                    )
                    {
                        label meshpointnext = meshPoints[lf[nextFp]];
                        label meshpointprev = meshPoints[lf[prevFp]];
                        label nedgei = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[meshpointi],
                            meshpointi,
                            meshpointnext
                        );
                        label pedgei = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[meshpointi],
                            meshpointi,
                            meshpointprev
                        );

                        if
                        (
                            featureEdge[nedgei] && featureEdge[pedgei]
                        )
                        {
                            vector nextEdge = newPoints[meshpointnext]
                                - newPoints[meshpointi];
                            nextEdge /= (mag(nextEdge) + SMALL);

                            vector prevEdge = newPoints[meshpointi]
                                - newPoints[meshpointprev];
                            prevEdge /= (mag(prevEdge) + SMALL);
                            if ((nextEdge & prevEdge) > 0.939)
                            {
                                flipFace = true;
                            }
                        }
                    }
                }
            }
            if (flipFace && !foundSwitch)
            {
                markedForFlip[i] = true;
            }
        }

        //Flip two feature faces connected by internal edge
        nEdgeFlipFaces = 0;
        forAll(pp.edges(), edgei)
        {
            label meshedgei = meshEdges[edgei];

            const labelList& eFaces = pp.edgeFaces()[edgei];
            forAll(eFaces, eFI)
            {
                if (markedForFlip[eFaces[eFI]])
                {
                    nEdgeFlipFaces[meshedgei]++;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            nEdgeFlipFaces,
            plusEqOp<label>(),
            label(0)
        );

        forAll(pp.edges(), edgei)
        {
            label meshedgei = meshEdges[edgei];
            if
            (
                nEdgeFlipFaces[meshedgei] == 2 && !featureEdge[meshedgei]
            )
            {
                const labelList& eFaces = pp.edgeFaces()[edgei];
                forAll(eFaces, eFI)
                {
                    label facei = eFaces[eFI];
                    if (markedForFlip[eFaces[eFI]])
                    {
                        bool flipFace = false;
                        face lf = pp.localFaces()[facei];
                        label startPt = -1;
                        forAll(lf,fp)
                        {
                            label pointi = lf[fp];
                            label meshpointi = meshPoints[pointi];
                            if (featurePoint[meshpointi] == -2)
                            {
                                startPt = pointi;
                                label nextFp = lf.fcIndex(fp);
                                label nextEdge = meshTools::findEdge
                                (
                                    pp.edges(),
                                    pp.pointEdges()[pointi],
                                    pointi,
                                    lf[nextFp]
                                );
                                if (nextEdge != edgei)
                                {
                                    flipFace = true;
                                }
                            }
                        }
                        if (flipFace)
                        {
                            lf.flip();
                        }

                        if (startPt != -1)
                        {
                            label fIndex = findIndex(lf, startPt);
                            fIndex = lf.fcIndex(fIndex);
                            label nextPt1 = lf[fIndex];
                            fIndex = lf.fcIndex(fIndex);
                            label nextPt2 = lf[fIndex];

                            if (nPEdges[nextPt1] > 0 && nPEdges[nextPt2] > 0)
                            {
                                sumDispPP[startPt] +=
                                    newPoints[meshPoints[nextPt1]];
                                sumAreasPP[startPt] += scalar(1);
                                switchPts[meshPoints[startPt]] = true;

                                sumDispPP[nextPt1] +=
                                    pEdgeCentre[nextPt1]/nPEdges[nextPt1];
                                sumAreasPP[nextPt1] += scalar(1);
                                switchPts[meshPoints[nextPt1]] = true;
                                sumDispPP[nextPt2] +=
                                    pEdgeCentre[nextPt2]/nPEdges[nextPt2];
                                sumAreasPP[nextPt2] += scalar(1);
                                switchPts[meshPoints[nextPt2]] = true;
                            }
                        }
                    }
                }
            }
        }

        //Fix convex faces
        forAll(pp, i)
        {
            const face& lf = pp.localFaces()[i];
            vector areaNorm = ppFaceAreas[i]
                / (mag(ppFaceAreas[i]) + SMALL);

            forAll(lf, fp)
            {
                label nextFp = lf.fcIndex(fp);
                label prevFp = lf.rcIndex(fp);
                label pointi = lf[fp];
                label meshpointi = meshPoints[lf[fp]];
                label meshpointnext = meshPoints[lf[nextFp]];
                label meshpointprev = meshPoints[lf[prevFp]];
                vector nextEdge = newPoints[meshpointnext]
                    - newPoints[meshpointi];
                nextEdge /= (mag(nextEdge) + SMALL);

                vector prevEdge = newPoints[meshpointprev]
                    - newPoints[meshpointi];
                prevEdge /= (mag(prevEdge) + SMALL);
                if
                (
                    mag(nextEdge & prevEdge) < 0.707
                    && !switchPts[meshpointi]
                 )
                {
                    vector pNorm = (nextEdge ^ prevEdge);
                    pNorm /= (mag(pNorm) + SMALL);
                    if ((pNorm & areaNorm) < 0)
                    {
                        switchPts[meshpointi] = true;
                        sumDispPP[pointi] += 0.5*
                        (
                            newPoints[meshpointprev]
                            +newPoints[meshpointnext]
                         );
                        sumAreasPP[pointi] += scalar(1);
                    }
                }
            }
        }
    } //if meshControl::EXTRUDE && iter > 20

    bool baffleProject((nSurfSmooth % 4) == 0 ? true : false);
    forAll(meshPoints, pointI)
    {
        const label meshPointI = meshPoints[pointI];
        if
        (
            !switchPts[meshPointI]
            && !stationaryPoints[meshPointI]
            && nFeatureEdges[meshPointI] < 3
        )
        {
            if (featurePoint[meshPointI] > -1)
            {
                if (baffleProject && featurePoint[meshPointI] == 2)
                {
                    const refinementFeatures& features =
                        meshRefiner_.features();

                    labelList nearEdgeFeat;
                    List<pointIndexHit> nearEdgeInfo;
                    vectorField nearNormal;

                    features.findNearestEdge
                    (
                        pointField(1, newPoints[meshPointI]),
                        scalarField(1, sqr(snapDist[pointI])),
                        nearEdgeFeat,
                        nearEdgeInfo,
                        nearNormal
                     );
                    const pointIndexHit& nearInfo = nearEdgeInfo[0];

                    if (nearInfo.hit())
                    {
                        sumDispPP[pointI] +=  nearInfo.hitPoint();
                        sumAreasPP[pointI] += scalar(1);
                    }
                }
                else
                {
                    const labelList& pEdges = pp.pointEdges()[pointI];
                    forAll(pEdges, eI)
                    {
                        label meshEdgeI = meshEdges[pEdges[eI]];
                        if
                        (
                            featureEdge[meshEdgeI]
                            && isMasterPPEdge[meshEdgeI]
                            && nFeatureEdges[meshPointI] == 2
                        )
                        {
                            edge e = mesh.edges()[meshEdgeI];
                            point eC = e.centre(newPoints);
                            scalar weight = e.mag(newPoints);
                            sumDispPP[pointI] +=  eC * weight;
                            sumAreasPP[pointI] += weight;
                        }
                    }
                }
            }
            else
            {
                const labelList& pFaces = pp.pointFaces()[pointI];
                forAll(pFaces, faceI)
                {
                    const label pFaceI = pFaces[faceI];
                    const label meshFaceI = pp.addressing()[pFaceI];
                    const label own = mesh.faceOwner()[meshFaceI];
                    label  level = cellLevel[own];

                    scalar len = edge0Len / pow(2., level);
                    scalar weight =
                        sqrt(mag(ppFaceAreas[pFaceI]))/ len;

                    sumDispPP[pointI] += ppFaceCentres[pFaceI]
                        * weight;
                    sumAreasPP[pointI] += weight;
                }
            }
        }
    }
    pointDisp = vector::zero;

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        sumDispPP,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        sumAreasPP,
        plusEqOp<scalar>(),
        scalar(0.0)         // null value
     );

    forAll(meshPoints, i)
    {
        label meshPointI = meshPoints[i];
        if (!stationaryPoints[meshPointI] && sumAreasPP[i] > SMALL)
        {
            pointDisp[meshPointI] = sumDispPP[i] / (sumAreasPP[i] + SMALL)
                - newPoints[meshPointI];
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pointDisp,
        maxMagSqrEqOp<point>(),
        vector::zero // null value (note: cannot use VGREAT)
    );

    forAll(mesh.points(), meshPointI)
    {
        if (switchPts[meshPointI])
        {
            newPoints[meshPointI] = newPoints[meshPointI]
                + pointDisp[meshPointI];
        }
        else if (!stationaryPoints[meshPointI])
        {
            newPoints[meshPointI] = newPoints[meshPointI]
                + smoothingWeights[meshPointI]*pointDisp[meshPointI];
        }
    }

    //try flattening internal faces
    if (iter > 10)
    {
        sumDispMesh = vector::zero;
        sumAreasMesh = 0.0;

        forAll(mesh.points(), meshPointI)
        {
            if (boundaryPoints[meshPointI] && !stationaryPoints[meshPointI])
            {
                labelList pFaces = mesh.pointFaces()[meshPointI];
                forAll(pFaces, pFI)
                {
                    if
                    (
                        !excludeFaces[pFaces[pFI]]
                        && isMasterFace[pFaces[pFI]]
                    )
                    {
                        // flatten the face
                        face f = mesh.faces()[pFaces[pFI]];
                        point fc = f.centre(newPoints);

                        point fn =  f.areaNormal(newPoints);
                        scalar fa = mag(fn);
                        fn /= (fa + SMALL);
                        vector pos =  fc - newPoints[meshPointI];

                        if (fa > SMALL)
                        {
                            sumDispMesh[meshPointI] += (1/fa)*(pos & fn)*fn;
                            sumAreasMesh[meshPointI] += (1/fa);
                        }
                    }
                }
            }
        }
        pointDisp = vector::zero;

        syncTools::syncPointList
        (
            mesh,
            sumDispMesh,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            sumAreasMesh,
            plusEqOp<scalar>(),
            scalar(0.0)         // null value
         );

        forAll(mesh.points(), meshPointI)
        {
            if (boundaryPoints[meshPointI])
            {
                if (sumAreasMesh[meshPointI] > SMALL)
                {
                    pointDisp[meshPointI] = 0.2*(sumDispMesh[meshPointI]
                                                 / sumAreasMesh[meshPointI]);
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointDisp,
            maxMagSqrEqOp<point>(),
            vector::zero // null value (note: cannot use VGREAT)
         );

        forAll(mesh.points(), meshPointI)
        {
            if (!stationaryPoints[meshPointI])
            {
                newPoints[meshPointI] += pointDisp[meshPointI];
            }
        }
    } //try flattening internal faces

    stationaryPoints = false;

    return;
}


// Snap face to surface and preserve face flatness
void Foam::snappySnapDriver::featureSnap
(
    const label maxIter,
    const label smoothStartIter,
    const labelList& snapSurfaces,
    const indirectPrimitivePatch& pp,
    const snapParameters& snapParams,
    const dictionary& motionDict,
    motionSmoother& meshMover,
    bool initialSnap,
    const optStages& ostages
) const
{
    Info<< "snappySnapDriver : Snap to Preserve Surface Features"<< endl;

    //Calculate the snap distance
    scalarField snapDist =
        calcSnapDistance(meshRefiner_, snapParams, pp);

    fvMesh& mesh = meshRefiner_.mesh();
    scalar startTime = mesh.time().elapsedCpuTime();

    // Undistorted edge length
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));
    const labelList& meshPoints = pp.meshPoints();

    // Storage used when moving points but not updating the mesh
    pointField newPoints = mesh.points();

    //Additional implicit snapping to capture acute and reflex angles
    scalar tightAngle = -1;
    scalar minTightAngle = -1;
    {
        scalar acuteAngle = snapParams.acuteReflexSnap();
        scalar maxAcuteAngle = 45.;
        if (acuteAngle > 0 && acuteAngle < maxAcuteAngle)
        {
            tightAngle = -Foam::cos(degToRad(acuteAngle));
            scalar minAngle = snapParams.minAcuteReflexSnap();
            if (minAngle >= 0 && minAngle < acuteAngle)
            {
                minTightAngle = -Foam::cos(degToRad(minAngle));
            }
            else
            {
                minTightAngle = tightAngle;
            }
        }
        else if (acuteAngle > maxAcuteAngle)
        {
            WarningInFunction
                << "Disabling acute/reflex angle snap"
                << " acuteReflexSnapAngle needs to be set > 0 "
                << " and < "<<maxAcuteAngle<<endl;
        }
    }

    scalar weakAngleMin = -1;
    scalar weakAngleMax = -1;
    bool weakSnap = false;

    const Tuple2<scalar,scalar> weakFeatureSnap = snapParams.weakFeatureSnap();

    if
    (
        weakFeatureSnap.first() >= 15 && weakFeatureSnap.first() <= 60
        && weakFeatureSnap.second() >= 15 && weakFeatureSnap.second() <= 60
    )
    {
        if (weakFeatureSnap.first() < weakFeatureSnap.second())
        {
            weakAngleMin = Foam::cos(degToRad(weakFeatureSnap.first()));
            weakAngleMax = Foam::cos(degToRad(weakFeatureSnap.second()));
            weakSnap = true;
        }
        else
        {
            WarningInFunction
                << "Minimum weak feature angle " << weakFeatureSnap.first()
                << "Greater than maximum weak feature angle "
                << weakFeatureSnap.second() << endl;
        }
    }

    vectorField pointDisp(mesh.nPoints(), vector::zero);

    // Updated pp face centre, areas and cell centres
    pointField ppFaceCentres(pp.size(), vector::zero);
    pointField ppCellCentres(pp.size(), vector::zero);
    vectorField ppFaceAreas(pp.size(), vector::zero);
    vectorField ppOrigFaceAreas(pp.size(), vector::zero);

    // Faces an points flagged as no longer requiring moving
    boolList stationaryFaces(pp.size(), false);
    boolList stationaryPoints(mesh.nPoints(), false);

    // For calculating feature edge points and edges
    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));

    // Get labels of faces to count (master of coupled faces and baffle pairs)
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    PackedBoolList isMasterPPEdge(snappyLayerDriver::getMasterPPEdges(mesh,pp));

    boolList excludeFaces(mesh.nFaces(), false);
    boolList boundaryPoints(mesh.nPoints(), false);

    scalar baffleAngle = -0.9848;
    if (minTightAngle > baffleAngle)
    {
        baffleAngle = minTightAngle;
    }

    List<Tuple2<edgeClassification::edgeType,scalar>> eType
    (
        meshEdges.size(),
        Tuple2<edgeClassification::edgeType,scalar>
        (
            edgeClassification::MANIFOLD,
            -GREAT
        )
    );

    //Check primitive patch edge types
    {
        boolList exclFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh,
            newPoints,
            pp,
            meshEdges,
            exclFaces,
            0.93969,//convex
            0.93969,//concave
            baffleAngle //baffle
         );
        eType = eClass.edgeTypes();
    }

    forAll(pp.localPoints(), i)
    {
        label meshPointI = meshPoints[i];
        boundaryPoints[meshPointI] = true;
    }

    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        excludeFaces[meshFaceI] = true;
    }

    syncTools::syncPointList
    (
        mesh,
        boundaryPoints,
        orEqOp<bool>(),
        false
    );

    forAll(mesh.points(), meshPointI)
    {
        if (boundaryPoints[meshPointI])
        {
            labelList pFaces = mesh.pointFaces()[meshPointI];

            forAll(pFaces, pFI)
            {
                if (!excludeFaces[pFaces[pFI]])
                {
                    face f = mesh.faces()[pFaces[pFI]];
                    label start = findIndex(f, meshPointI);
                    label next = f.fcIndex(start);
                    label prev = f.rcIndex(start);
                    if (boundaryPoints[f[next]] && boundaryPoints[f[prev]])
                    {
                        excludeFaces[pFaces[pFI]] = true;
                    }
                    else if
                    (
                        !boundaryPoints[f[next]] && !boundaryPoints[f[prev]]
                    )
                    {
                        excludeFaces[pFaces[pFI]] = true;
                    }
                }
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh,
        excludeFaces,
        orEqOp<bool>()
    );

    // Mark pp zone edges
    boolList zonedEdges(meshEdges.size(), false);
    {
        PackedList<1> isZonedFace(setZonedFaces(meshRefiner_));
        forAll(pp, i)
        {
            label meshfacei = pp.addressing()[i];
            if (isZonedFace.get(meshfacei))
            {
                const labelList& fEdges = pp.faceEdges()[i];
                forAll(fEdges, fEI)
                {
                    zonedEdges[fEdges[fEI]] = true;
                }
            }
        }
    }


    scalarField smoothingWeights(mesh.nPoints(), 0.2);
    forAll(pp.edges(), edgei)
    {
        if (eType[edgei].first() == edgeClassification::NONMANIFOLD)
        {
            const edge e = pp.edges()[edgei];
            forAll(e, eI)
            {
                const labelList pEdges = pp.pointEdges()[e[eI]];
                forAll(pEdges, pEI)
                {
                    label oEdge = pEdges[pEI];
                    if
                    (
                        eType[oEdge].first() != edgeClassification::NONMANIFOLD
                        && eType[oEdge].first() != edgeClassification::BOUNDARY
                    )
                    {
                        edge oe = pp.edges()[oEdge];
                        label otherPt(oe[0] == e[eI] ? oe[1] : oe[0]);
                        label otherMeshPt = meshPoints[otherPt];
                        smoothingWeights[otherMeshPt] = 1.0;
                    }
                }
            }
        }
        else if
        (
            eType[edgei].first() == edgeClassification::BAFFLE
            ||
            (
                eType[edgei].first() == edgeClassification::BOUNDARY
                &&  zonedEdges[edgei]
            )
        )
        {
            const edge e = pp.edges()[edgei];
            forAll(e, ei)
            {
                smoothingWeights[meshPoints[e[ei]]] = 0.5;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        smoothingWeights,
        maxEqOp<scalar>(),
        scalar(0)
    );

    // After fastIter feature snapping iterations perform
    // fast surface interrogation
    label fastIter = 3;
    bool fastFeatureSnap = false;

    scalarField sumAreasPP(pp.nPoints(), 0.0);
    vectorField sumDispPP(pp.nPoints(), vector::zero);

    scalarField sumAreasMesh(mesh.nPoints(), 0.0);
    vectorField sumDispMesh(mesh.nPoints(), vector::zero);

    bool smoothSurface = snapParams.smoothSnappedSurface();

    boolList gapBoundaryPts(mesh.nPoints(), false);

    featureSnapInit
    (
        initialSnap,
        newPoints,
        pointDisp,
        sumAreasMesh,
        sumDispMesh,
        gapBoundaryPts
    );

    sumAreasMesh = 0.0;
    sumDispMesh = vector::zero;

    label nSurfSmooth = 0;
    for (label iter = 0; iter < maxIter; iter++)
    {
        //Relax snap distance when starting to lock onto the surface
        if (iter == 20)
        {
            snapDist *= 0.25;
        }

        if (snapParams.writeSnapVTK())
        {
            simpleVTKWriter
            (
                pp(),
                newPoints
            ).write("snappyFtrIter"+Foam::name(iter)+".vtk");
        }

        if (iter % 10 == 0)
        {
            Info<<"Snapping iteration: "<<iter<<endl;
        }

        // Speed up feature snapping iteration by only considering
        // face centres. Initially  perform snapping by using face
        // centres and face points to find best alignment of face with surface
        if (iter >= fastIter)
        {
            fastFeatureSnap = true;
        }

        // Smooth surface points (feature edge aware)
        if (smoothSurface && (iter % 2 == 0) && iter > smoothStartIter)
        {
            nSurfSmooth++;
            boolList exclFaces(pp.size(), false);
            edgeClassification eClass
            (
                mesh,
                newPoints,
                pp,
                meshEdges,
                exclFaces,
                0.93969,//convex
                0.93969,//concave
                baffleAngle //baffle
            );
            eType = eClass.edgeTypes();

            featureSnapSurfaceSmooth
            (
                iter,
                nSurfSmooth,
                pp,
                eType,
                isMasterEdge,
                isMasterFace,
                isMasterPPEdge,
                snapDist,
                boundaryPoints,
                smoothingWeights,
                excludeFaces,
                zonedEdges,
                newPoints,
                pointDisp,
                ppFaceCentres,
                ppFaceAreas,
                stationaryPoints,
                sumAreasMesh,
                sumDispMesh,
                sumAreasPP,
                sumDispPP
            );
        } //smooth surface point

        // Update cell centre, face centre and face areas
        // after each feature snapping iteration.
        forAll(pp, i)
        {
            if (iter == 0 && !initialSnap)
            {
                const label& faceI = pp.addressing()[i];
                const label& own = mesh.faceOwner()[faceI];
                ppCellCentres[i] =
                    mesh.cellCentres()[own];
                ppFaceCentres[i] =
                    mesh.faceCentres()[faceI];
                ppFaceAreas[i] =
                    mesh.faceAreas()[faceI];
            }
            else
            {
                const label& faceI = pp.addressing()[i];
                const label& own = mesh.faceOwner()[faceI];

                ppCellCentres[i] =
                    mesh.cells()[own].centre(newPoints, mesh.faces());
                ppFaceCentres[i] =
                    mesh.faces()[faceI].centre(newPoints);
                ppFaceAreas[i] =
                    mesh.faces()[faceI].areaNormal(newPoints);
            }
        }
        if (iter == 0)
        {
            ppOrigFaceAreas = ppFaceAreas;
        }

        // Stop feature snapping faces where face area is getting too
        // large or small
        if (iter > 5)
        {
            forAll(pp, i)
            {
                const face& f = pp.localFaces()[i];
                point fc = ppFaceCentres[i];
                vectorField flatFace(f.size(), vector::zero);
                Pair<vector> n1n2;
                n1n2[0] =  ppFaceAreas[i];
                n1n2[0] /= (mag(n1n2[0]) + SMALL);

                forAll(f, fp)
                {
                    flatFace[fp] = newPoints[meshPoints[f[fp]]];
                    vector pos = newPoints[meshPoints[f[fp]]] - fc;
                    flatFace[fp] -= (pos & n1n2[0])*n1n2[0];
                }

                point newFc = face(identity(f.size())).centre(flatFace);

                scalar minLength = GREAT;
                scalar maxLength = -GREAT;
                forAll(f, fp)
                {
                    point curPt = flatFace[fp];
                    point nextPt = flatFace[f.fcIndex(fp)];
                    scalar projDist = linePointRef(curPt, nextPt)
                        .nearestDist(newFc).distance();
                    minLength = min(minLength, projDist);
                    maxLength = max(maxLength, projDist);
                }

                scalar aRatio = minLength/(maxLength + SMALL);

                scalar areaRatio = mag(ppFaceAreas[i])
                    /(mag(ppOrigFaceAreas[i]) + SMALL);
                if (areaRatio < 0.05 || areaRatio > 20. || aRatio < 0.05)
                {
                    stationaryFaces[i] = true;
                    forAll(f, fp)
                    {
                        stationaryPoints[meshPoints[f[fp]]] = true;
                    }
                }
            }
        } // Stop feature snapping faces where face area is getting too large or small

        syncTools::syncPointList
        (
            mesh,
            stationaryPoints,
            orEqOp<bool>(),
            false               // null value
        );

        sumAreasPP = 0.0;
        sumDispPP = vector::zero;
        List<pointIndexHit> finalHitInfo;
        labelList pointIndex(pp.size(),label(-1));
        vectorField pNormals(pp.size(),vector::zero);

        findSurfaceNormals
        (
            snapSurfaces,
            pp,
            snapDist,
            newPoints,
            ppFaceCentres,
            ppCellCentres,
            ppFaceAreas,
            stationaryFaces,
            fastFeatureSnap,
            finalHitInfo,
            pointIndex,
            pNormals
        );

        //Calculate average normal
        if
        (
            iter > 10 && snapParams.averageSurfaceNormal()
            && controller_.mode() == meshControl::QUALITY
        )
        {
            recalculatePointNormals
            (
                snapParams,
                pp,
                snapSurfaces,
                snapDist,
                newPoints,
                finalHitInfo,
                pointIndex,
                ppFaceCentres,
                ppFaceAreas,
                pNormals
            );
        }

        List<DynamicList<point>> allHits(pp.meshPoints().size());
        forAll(pointIndex, fI)
        {
           label i = pointIndex[fI];
           const face& f = pp.localFaces()[i];
           const point& fc = ppFaceCentres[i];
           const  vector fArea = ppFaceAreas[i] / (mag(ppFaceAreas[i]) + SMALL);

           vector surfNormToSnap = pNormals[fI];

           if (!finalHitInfo[fI].hit())
           {
               continue;
           }
           else if
           (
               (
                   controller_.algorithm() == meshControl::DUAL
                   || controller_.algorithm() == meshControl::EXTRUDE
               )
               && fastFeatureSnap
           )
           {
               //check for snap moving wrong side of cell-centre
               // after initial snapping phase has completed
               vector snapDir = finalHitInfo[fI].hitPoint() - fc;

               if ((snapDir& fArea) < 0)
               {
                   vector ccDir = ppCellCentres[i] - fc;
                   scalar ccNormDist = mag(ccDir & fArea);
                   scalar snapNormDist = mag(snapDir & fArea);
                   if (snapNormDist > 1.5*ccNormDist)
                   {
                       continue;
                   }
               }
           }

           point surfPointToSnap = finalHitInfo[fI].hitPoint();

           // flatten the face before moving
           vectorField flatFace(f.size(), vector::zero);
           Pair<vector> n1n2;
           n1n2[0] =  fArea;

           forAll(f, fp)
           {
              flatFace[fp] = newPoints[meshPoints[f[fp]]];
              vector pos = newPoints[meshPoints[f[fp]]] - fc;
              flatFace[fp] -= (pos & n1n2[0])*n1n2[0];
           }

           label nPoints = f.size();
           // calculate new face centre and area of flattened face
           vector newFaceCentre, newFaceArea;

           // If the face is a triangle, do a direct calculation for
           // efficiency and to avoid round-off error-related problems
           if (nPoints == 3)
           {
               newFaceCentre = (1.0/3.0)*(flatFace[0] + flatFace[1]
                             + flatFace[2]);
               newFaceArea = 0.5*((flatFace[1] - flatFace[0])^(flatFace[2]
                           - flatFace[0]));
           }
           else
           {
               vector sumN = vector::zero;
               scalar sumA = 0.0;
               vector sumAc = vector::zero;
               point fCentre = flatFace[0];
               for (label pi = 1; pi < nPoints; pi++)
               {
                   fCentre += flatFace[pi];
               }
               fCentre /= nPoints;

               for (label pi = 0; pi < nPoints; pi++)
               {
                   const point& nextPoint = flatFace[(pi + 1) % nPoints];
                   vector c = flatFace[pi] + nextPoint + fCentre;
                   vector n = (nextPoint - flatFace[pi]) ^
                      (fCentre - flatFace[pi]);
                   scalar a = mag(n);
                   sumN += n;
                   sumA += a;
                   sumAc += a*c;
               }

               newFaceCentre = (1.0/3.0)*sumAc/(sumA + VSMALL);
               newFaceArea = 0.5*sumN;
           }

           n1n2[0] = newFaceArea;
           n1n2[1] = surfNormToSnap;
           n1n2[0] /= (mag(n1n2[0]) + VSMALL);
           n1n2[1] /= (mag(n1n2[1]) + VSMALL);

           scalar dProd = (n1n2[1] & n1n2[0]);
           if (dProd > 0)
           {
              n1n2[1] = -n1n2[1];
           }

           scalar fWeight = max(mag(dProd),scalar(0.25));

           tensor T = rotationTensor(n1n2[0], -n1n2[1]);

           forAll(f, fp)
           {
               // rotate points
               point p = flatFace[fp];

               p -= newFaceCentre;
               p = transform(T, p);
               p += newFaceCentre;

               // translate to flatten face
               p += (flatFace[fp] - newPoints[meshPoints[f[fp]]]);

               // translate points
               p += (surfPointToSnap - newFaceCentre);

               vector disp = p - newPoints[meshPoints[f[fp]]];

               allHits[f[fp]].append(p);
               scalar weight = mag(newFaceArea);

               sumDispPP[f[fp]] += fWeight * disp / (weight + SMALL);
               sumAreasPP[f[fp]] += fWeight / (weight + SMALL);
           }
        }

        if
        (
            tightAngle > -1 && smoothSurface
            && iter > smoothStartIter
            && ((iter % 4 == 0) || (iter == maxIter -1))
        )
        {
            tightAngleCorrect
            (
                weakSnap,
                weakAngleMin,
                weakAngleMax,
                tightAngle,
                minTightAngle,
                newPoints,
                pp,
                meshEdges,
                eType,
                finalHitInfo,
                pointIndex,
                ppFaceAreas,
                pNormals,
                sumAreasPP,
                sumDispPP
            );
        }

        forAll(pp, i)
        {
            if (stationaryFaces[i])
            {
                const face& f = pp.localFaces()[i];
                forAll(f, fp)
                {
                    sumAreasPP[f[fp]] += 1. /(mag(ppFaceAreas[i]) + SMALL);
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            sumDispPP,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            sumAreasPP,
            plusEqOp<scalar>(),
            scalar(0.0)         // null value
        );

        pointDisp = vector::zero;
        forAll(meshPoints, i)
        {
            if (!stationaryPoints[meshPoints[i]])
            {
                label meshPointI = meshPoints[i];
                pointDisp[meshPointI] = sumDispPP[i] / (sumAreasPP[i] + SMALL);
            }
            else
            {
                pointDisp[meshPoints[i]] = vector::zero;
            }

        }

        scalarField allHitsDisp(meshPoints.size(), 0);
        labelList nHits(meshPoints.size(), 0);
        forAll(meshPoints, i)
        {
            label meshPointI = meshPoints[i];
            if (!stationaryPoints[meshPointI] && allHits[i].size()>0)
            {
                point projPt = newPoints[meshPointI] +  pointDisp[meshPointI];

                forAll(allHits[i], hitI)
                {
                    allHitsDisp[i] += mag(allHits[i][hitI]-projPt);
                    nHits[i]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            allHitsDisp,
            plusEqOp<scalar>(),
            scalar(0.0)         // null value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nHits,
            plusEqOp<label>(),
            label(0)         // null value
        );

        forAll(meshPoints, i)
        {
            label meshPointI = meshPoints[i];

            // The last condition keeps the denominator of edgeLevelLen from
            // being undefined behaviour. If we presume
            // that left shift by a negative number is intended to be right
            // shift, then it'll be infinite. Which is then divided by to
            // calculate relativeDisp, which gives zero, which is <= 0.25.
            if
            (
                !stationaryPoints[meshPointI]
                && nHits[i] >0 && pointLevel[meshPointI] >= 0
            )
            {
                scalar edgeLevelLen = edge0Len /
                    (1<<pointLevel[meshPointI]);
                scalar aveDisp = allHitsDisp[i] / nHits[i];

                scalar relativeDisp =  aveDisp/edgeLevelLen;
                if (relativeDisp > 0.25)
                {
                    pointDisp[meshPointI] *= 0.5;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointDisp,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        forAll(mesh.points(), meshPointI)
        {
            if (!stationaryPoints[meshPointI])
            {
                newPoints[meshPointI] = newPoints[meshPointI]
                    + pointDisp[meshPointI];
            }
        }

        //Perform intermediate optimisation stages
        featureOptimise
        (
            iter,
            ostages,
            snapParams,
            motionDict,
            pp,
            snapSurfaces,
            meshEdges,
            boundaryPoints,
            gapBoundaryPts,
            newPoints,
            meshMover,
            stationaryPoints
        );

        // stop generation of small faces (causes sync issues)
        boolList smallPointFaces(mesh.nPoints(), false);
        forAll(meshPoints, i)
        {
            label meshPointI = meshPoints[i];
            const labelList& pointFaces = mesh.pointFaces()[meshPointI];

            forAll(pointFaces, faceI)
            {
                scalar faceArea =
                    mag(mesh.faces()[pointFaces[faceI]].areaNormal(newPoints));

                if (faceArea < SMALL)
                {
                    if (!stationaryPoints[meshPointI])
                    {
                        stationaryPoints[meshPointI] = true;
                        smallPointFaces[meshPointI] = true;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            smallPointFaces,
            orEqOp<bool>(),
            false               // null value
        );

        forAll(smallPointFaces, pointI)
        {
            if (smallPointFaces[pointI])
            {
                newPoints[pointI] = newPoints[pointI]
                            - pointDisp[pointI];
            }
        }

        stationaryFaces = false;
    }

    // Displacement at every patch point
    pointField patchDisp(meshMover.patch().nPoints(), vector::zero);

    pointDisp = vector::zero;

    forAll(pp.localPoints(), i)
    {
        pointDisp[meshPoints[i]] =
            newPoints[meshPoints[i]] - pp.localPoints()[i];
    }

    forAll(patchDisp, i)
    {
        patchDisp[i] = pointDisp[meshMover.patch().meshPoints()[i]];
    }

    // The current mesh is the starting mesh to smooth from.
    meshMover.setDisplacement(patchDisp);

    scalar endTime = mesh.time().elapsedCpuTime();
    Info<< "Snap to preserve surface features in " << endTime - startTime
         << "s" << endl << endl;
}


void Foam::snappySnapDriver::extrudeCorrect
(
    const indirectPrimitivePatch& pp,
    const dictionary& motionDict,
    const boolList& boundaryPoints,
    motionSmoother& meshMover,
    pointField& newPoints
) const
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));
    const labelList& meshPoints = pp.meshPoints();

    boolList layerEdges(mesh.nEdges(), false);
    boolList staticPts(mesh.nPoints(), false);

    forAll(mesh.edges(), edgei)
    {
        const edge& e = mesh.edges()[edgei];
        if
            (
                (boundaryPoints[e[0]] && !boundaryPoints[e[1]])
                || (boundaryPoints[e[1]] && !boundaryPoints[e[0]])
             )
        {
            layerEdges[edgei] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        layerEdges,
        orEqOp<bool>(),
        false
    );

    boolList exclFaces(pp.size(), false);
    edgeClassification eClass
    (
        mesh,
        newPoints,
        pp,
        meshEdges,
        exclFaces,
        0.5,
        -0.5
    );
    const List<Tuple2<edgeClassification::edgeType,scalar>>&
        eType = eClass.edgeTypes();

    boolList featurePts(mesh.nPoints(), false);
    forAll(eType, edgei)
    {
        if
        (
            eType[edgei].first() == edgeClassification::CONCAVE
            || eType[edgei].first() == edgeClassification::CONVEX
        )
        {
            label meshEdgeI = meshEdges[edgei];
            edge e = mesh.edges()[meshEdgeI];
            featurePts[e[0]] = true;
            featurePts[e[1]] = true;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        featurePts,
        orEqOp<bool>(),
        false
    );

    const volScalarField& layerCells =
        mesh.lookupObject<volScalarField>("layerStacks");
    scalarField layerPts(mesh.nPoints(), -1);
    forAll(layerCells, celli)
    {
        scalar layerCellID = layerCells[celli];

        if (layerCellID  > -1)
        {
            const labelList& cPts = mesh.cellPoints()[celli];
            forAll(cPts, cPtI)
            {
                label pointi = cPts[cPtI];
                layerPts[pointi] =
                    max(layerCellID,layerPts[pointi]);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        layerPts,
        maxEqOp<scalar>(),
        scalar(-1) // null value
    );

    boolList gapPts(mesh.nPoints(), false);
    forAll(layerCells, celli)
    {
        scalar layerCellID = layerCells[celli];

        if (layerCellID  < 0)
        {
            const labelList& cPts = mesh.cellPoints()[celli];
            label nLayerPts = 0;
            forAll(cPts, cPtI)
            {
                label pointi = cPts[cPtI];
                if (layerPts[pointi] > -1)
                {
                    nLayerPts++;
                }
            }
            if (nLayerPts == cPts.size())
            {
                forAll(cPts, cPtI)
                {
                    label pointi = cPts[cPtI];
                    gapPts[pointi] = true;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        gapPts,
        orEqOp<bool>(),
        false // null value
    );

    pointField newDisp(mesh.nPoints(), vector::zero);
    labelList nMoved(mesh.nPoints(), 0);

    edgeClassification eClassNormals
    (
        mesh,
        newPoints,
        pp,
        meshEdges,
        exclFaces,
        0.5,
        0.5
    );
    pointField surfNormals = eClassNormals.calculatePointNormals
    (
        exclFaces,
        0,
        true
    );

    forAll(surfNormals, pti)
    {
        label meshpointi = meshPoints[pti];
        const labelList& pEdges = mesh.pointEdges()[meshpointi];
        forAll(pEdges, pEI)
        {
            label edgei = pEdges[pEI];
            if (layerEdges[edgei])
            {
                edge e = mesh.edges()[edgei];
                label otherPt
                (
                    e[0] == meshpointi ? e[1] : e[0]
                );

                if (!featurePts[meshpointi] && !gapPts[otherPt])
                {
                    continue;
                }

                point startPt = newPoints[meshpointi];
                point outerPt = newPoints[otherPt];
                if (gapPts[otherPt])
                {
                    vector sN = -surfNormals[pti];
                    scalar magSn = mag(sN);
                    if (magSn < SMALL)
                    {
                       continue;
                    }
                    sN /= magSn;

                    const labelList& pCells = mesh.pointCells()[otherPt];
                    forAll(pCells, pCI)
                    {
                        label celli = pCells[pCI];
                        scalar layerCellID = layerCells[celli];

                        if (layerCellID  < 0)
                        {
                            point cc = mesh.cells()[celli].centre
                            (
                                newPoints,
                                mesh.faces()
                            );

                            plane fPlane(startPt, sN);
                            scalar hitDist = mag(fPlane.nearestPoint(cc)-cc);
                            point endPt = startPt + 0.5*hitDist*sN;
                            staticPts[otherPt] = true;
                            nMoved[otherPt]++;
                            newDisp[otherPt] += (endPt-outerPt);
                        }
                    }
                }
                else
                {
                    vector layerVec = outerPt - startPt;
                    scalar layerHeight = mag(layerVec);
                    if (layerHeight > SMALL)
                    {
                        layerVec /= layerHeight;
                        vector sN = -surfNormals[pti];
                        sN /= mag(sN);
                        scalar dProd = (layerVec & sN);
                        if (dProd >SMALL)
                        {
                            point endPt = startPt
                                + dProd*layerHeight*sN;
                            staticPts[otherPt] = true;
                            nMoved[otherPt]++;
                            newDisp[otherPt] += (endPt-outerPt);
                        }
                    }
                }
             }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        newDisp,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nMoved,
        plusEqOp<label>(),
        label(0)     // null value
    );

    forAll(nMoved, pointi)
    {
        if (nMoved[pointi] > 0)
        {
            newPoints[pointi] += newDisp[pointi]/nMoved[pointi];
        }
    }

    syncTools::syncPointList
    (
        mesh,
        staticPts,
        orEqOp<bool>(),
        false     // null value
    );

    //re-project non-feature pt on faces with two feature edges
    newDisp = vector::zero;
    nMoved = 0;
    DynamicList<label> ftrCells(mesh.nCells()/100);

    forAll(pp, i)
    {
        const face& f = pp.localFaces()[i];
        label meshfacei = pp.addressing()[i];

        if (f.size() != 4)
        {
            continue;
        }
        label midPt = -1;
        label prevPt = -1;
        label nextPt = -1;
        label otherPt = -1;
        bool valid = false;
        label nFtrEdges = 0;

        forAll(f, fp)
        {
            label pointi = f[fp];
            label meshpointi = meshPoints[pointi];
            point pt = newPoints[meshpointi];
            label nextFp = f.fcIndex(fp);
            label prevFp = f.rcIndex(fp);
            label nextEdge = meshTools::findEdge
            (
                pp.edges(),
                pp.pointEdges()[pointi],
                pointi,
                f[nextFp]
            );

            label prevEdge = meshTools::findEdge
            (
                pp.edges(),
                pp.pointEdges()[pointi],
                pointi,
                f[prevFp]
            );

            bool nextFE
            (
                eType[nextEdge].first() == edgeClassification::CONCAVE
                || eType[nextEdge].first() == edgeClassification::CONVEX
            );

            bool prevFE
            (
                eType[prevEdge].first() == edgeClassification::CONCAVE
                || eType[prevEdge].first() == edgeClassification::CONVEX
            );

            if (nextFE)
            {
                nFtrEdges++;
            }

            if (!nextFE && !prevFE)
            {
                otherPt = meshpointi;
            }
            else if (nextFE && prevFE)
            {
                label mpn = meshPoints[f[nextFp]];
                label mpp = meshPoints[f[prevFp]];
                vector nextVec = newPoints[mpn] - pt;
                vector prevVec = pt - newPoints[mpp];
                nextVec /= Foam::mag(nextVec) + VSMALL;
                prevVec /= Foam::mag(prevVec) + VSMALL;
                if ((nextVec & prevVec) > 0.939)
                {
                    prevPt = mpp;
                    nextPt = mpn;
                    midPt = meshpointi;
                    valid = true;
                }
            }
        }
        if (nFtrEdges == 2 && valid)
        {
            vector edgeVec = newPoints[nextPt] - newPoints[prevPt];
            scalar edgeLen = mag(edgeVec);
            if (edgeLen > SMALL)
            {
                edgeVec /= edgeLen;
                vector otherDir = newPoints[otherPt] - newPoints[prevPt];
                scalar tanLen = mag(otherDir & edgeVec);
                scalar hypLen = mag(otherDir);
                scalar normLen = sqrt(hypLen*hypLen-tanLen*tanLen);
                face mf = mesh.faces()[meshfacei];

                vector fA = mf.areaNormal(newPoints);
                scalar area = mag(fA);
                if (area > SMALL)
                {
                    vector unitNorm = fA / area;
                    vector edgeDir = (unitNorm ^ edgeVec);
                    point projPt = newPoints[midPt] + normLen*edgeDir;
                    newDisp[otherPt] += (projPt-newPoints[otherPt]);
                    nMoved[otherPt]++;
                    //add owner cell for smoothing later
                    ftrCells.append(mesh.faceOwner()[meshfacei]);
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        newDisp,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh,
        nMoved,
        plusEqOp<label>(),
        label(0)     // null value
    );

    forAll(nMoved, pointi)
    {
        if (nMoved[pointi] > 0)
        {
            newPoints[pointi] += newDisp[pointi]/nMoved[pointi];
        }
    }

    pointField avePt(mesh.nPoints(), vector::zero);
    scalarField aveWt(mesh.nPoints(), 0);

    newDisp = newPoints - meshMover.oldPoints();
    syncTools::syncPointList
    (
        mesh,
        newDisp,
        minMagSqrEqOp<point>(),
        vector(GREAT, GREAT, GREAT)
    );
    newDisp += meshMover.oldPoints();
    mesh.movePoints(newDisp);



    faceSet wrongFaces
    (
        mesh,
        "wrongFaces",
        mesh.nFaces()/100+100
    );

    //Smooth cells with volume errors
    dictionary updatedMotionDict;
    if (motionDict.found("minFaceWeight"))
    {
        updatedMotionDict.add
        (
            "minFaceWeight",
            motionDict.lookup("minFaceWeight")
            ,true
        );
    }
    if (motionDict.found("minVolRatio"))
    {
        updatedMotionDict.add
        (
            "minVolRatio",
            motionDict.lookup("minVolRatio")
            ,true
        );
    }
    if (motionDict.found("maxConcave"))
    {
        updatedMotionDict.add
        (
            "maxConcave",
            motionDict.lookup("maxConcave")
            ,true
        );
    }

    motionSmoother::checkMesh
    (
        false,
        mesh,
        updatedMotionDict,
        wrongFaces
    );

    labelHashSet usedPoints(mesh.nPoints()/100);

    forAllConstIter(labelHashSet, wrongFaces, iter)
    {
        usedPoints.insert(mesh.faces()[iter.key()]);
    }

    forAll(ftrCells, i)
    {
        const cell& c = mesh.cells()[ftrCells[i]];
        forAll(c, cfi)
        {
            usedPoints.insert(mesh.faces()[c[cfi]]);
        }
    }

    pointSet wrongPoints
    (
        mesh,
        "wrongPoints",
        usedPoints
    );
    wrongPoints.sync(mesh);

    forAllConstIter(pointSet, wrongPoints, iter)
    {
        label pointi = iter.key();
        if (!boundaryPoints[pointi] && !staticPts[pointi])
        {
            const labelList& pCells =
                mesh.pointCells()[pointi];
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                scalar vol = mesh.cellVolumes()[celli];
                point cc = mesh.cellCentres()[celli];
                avePt[pointi] += cc*vol;
                aveWt[pointi] += vol;
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        avePt,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh,
        aveWt,
        plusEqOp<scalar>(),
        scalar(0)     // null value
    );

    forAll(mesh.points(), pointi)
    {
        if (aveWt[pointi] > SMALL)
        {
            newPoints[pointi] = 0.5*
            (
                newPoints[pointi]
                +(avePt[pointi]/aveWt[pointi])
            );
        }
    }

    // Calculate coupled layerID
    scalarField
        neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());

    avePt = vector::zero;
    aveWt = scalar(0);

    for
        (
            label facei = mesh.nInternalFaces();
            facei < mesh.nFaces();
            facei++
         )
    {
        neiLayerCells[facei-mesh.nInternalFaces()] =
            layerCells[mesh.faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLayerCells);

    pointField bNorm(mesh.nFaces(), vector::zero);
    pointField bCentre(mesh.nFaces(), vector::zero);

    forAll(mesh.cells(), celli)
    {
        scalar layerCellID = layerCells[celli];
        const cell& c = mesh.cells()[celli];
        if (layerCellID  > -1 && c.size() == 6)
        {
            point fC = vector::zero;
            point fN = vector::zero;
            label outerFace = -1;

            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchI = patches.whichPatch(facei);
                label own = mesh.faceOwner()[facei];
                scalar neiLayer = GREAT;

                if (patchI == -1)
                {
                    neiLayer =
                    (
                        own == celli ?
                        layerCells[mesh.faceNeighbour()[facei]]
                        : layerCells[own]
                    );
                }
                else if (patches[patchI].coupled())
                {
                    label bFaceI = facei-mesh.nInternalFaces();
                    neiLayer = neiLayerCells[bFaceI];
                }
                else
                {
                    const face f = mesh.faces()[facei];
                    fC = f.centre(newPoints);
                    fN = f.unitNormal(newPoints);
                }

                if (neiLayer < 0)
                {
                    outerFace = facei;
                }
            }

            if (outerFace != -1)
            {
                bCentre[outerFace] = fC;
                bNorm[outerFace] = fN;
            }
        }
    }
    syncTools::syncFaceList
    (
        mesh,
        bCentre,
        maxMagSqrEqOp<point>()
     );

    syncTools::syncFaceList
    (
        mesh,
        bNorm,
        maxMagSqrEqOp<point>()
    );

    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] < 0)
        {
            const labelList& cPts = mesh.cellPoints()[celli];
            bool allLayers = true;

            forAll(cPts, cPtI)
            {
                if (layerPts[cPts[cPtI]] < 0)
                {
                    allLayers = false;
                    break;
                }
            }
            if (allLayers)
            {
                const cell& c = mesh.cells()[celli];
                DynamicList<point> bN(c.size());
                DynamicList<point> bC(c.size());
                forAll(c, cFI)
                {
                    label facei = c[cFI];
                    if (mag(bNorm[facei]) > SMALL)
                    {
                        bN.append(bNorm[facei]);
                        bC.append(bCentre[facei]);
                    }
                }

                if (bN.size() == 2)
                {
                    point midPt = (bC[0] + bC[1]) / 2.0;
                    point midNorm = (bN[0] - bN[1]) / 2.0;

                    DynamicList<label> oppFaces(2);
                    DynamicList<label> oppEdges(8);
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        if (mag(bNorm[facei]) > SMALL)
                        {
                            oppFaces.append(facei);
                            const labelList& fEdges =
                                mesh.faceEdges()[facei];
                            forAll(fEdges, fEI)
                            {
                                oppEdges.append(fEdges[fEI]);
                            }
                        }
                    }

                    labelHashSet cEdgesSet
                    (
                        mesh.cellEdges()[celli]
                    );
                    labelHashSet oppEdgesSet(oppEdges);

                    forAll(oppFaces, oFI)
                    {
                        label facei = oppFaces[oFI];
                        const face f = mesh.faces()[facei];
                        point plPt = (2.0/3.0)*midPt
                            + (1.0/3.0)*bCentre[facei];

                        plane pl
                        (
                            plPt,
                            midNorm
                        );
                        forAll(f,fp)
                        {
                            label pointi =f[fp];
                            const labelList& pEdges =
                                mesh.pointEdges()[pointi];

                            label cEdge = -1;
                            forAll(pEdges, pEI)
                            {
                                label edgei = pEdges[pEI];
                                if
                                (
                                    cEdgesSet.found(edgei)
                                    && !oppEdgesSet.found(edgei)
                                )
                                {
                                    cEdge = edgei;
                                    break;
                                }
                            }

                            if (cEdge != -1)
                            {
                                const edge e =
                                    mesh.edges()[cEdge];
                                point aveEdgePts = 0.5*
                                (
                                    newPoints[e[0]]
                                    +newPoints[e[1]]
                                );

                                point projPt =
                                    pl.nearestPoint(aveEdgePts);
                                avePt[pointi] +=
                                    (projPt-newPoints[pointi]);
                                aveWt[pointi] += scalar(1);
                            }
                        }
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        avePt,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh,
        aveWt,
        plusEqOp<scalar>(),
        scalar(0)     // null value
    );

    forAll(mesh.points(), pointi)
    {
        if (aveWt[pointi] > SMALL)
        {
            newPoints[pointi] =
                (
                    newPoints[pointi]
                    +(avePt[pointi]/aveWt[pointi])
                 );
        }
    }

    mesh.movePoints(meshMover.oldPoints());

    return;
}


void Foam::snappySnapDriver::cornerCorrect
(
    const indirectPrimitivePatch& pp,
    const labelList& meshEdges,
    const boolList& boundaryPoints,
    pointField& newPoints
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& meshPoints = pp.meshPoints();

    const volScalarField& layerCells =
        mesh.lookupObject<volScalarField>("layerStacks");
    // Calculate coupled layerID
    scalarField neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());

    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        neiLayerCells[facei-mesh.nInternalFaces()] =
            layerCells[mesh.faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLayerCells);

    pointField avePt(mesh.nPoints(), vector::zero);
    scalarField aveWt(mesh.nPoints(), 0);

    labelHashSet usedPoints(mesh.nPoints()/100);
    boolList markedPts(mesh.nPoints(),false);

    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] != -1)
        {
            const labelList& cPts  = mesh.cellPoints()[celli];
            forAll(cPts, cPI)
            {
                markedPts[cPts[cPI]] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        markedPts,
        orEqOp<bool>(),
        false               // null value
    );

    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] == -1)
        {
            const labelList& cPts  = mesh.cellPoints()[celli];
            label nLayerPts = 0;

            forAll(cPts, cPI)
            {
                if (markedPts[cPts[cPI]])
                {
                    nLayerPts++;
                }
            }

            if (nLayerPts == cPts .size())
            {
                const cell& c = mesh.cells()[celli];
                label nLayerStacks = 0;
                forAll(c, cfi)
                {
                    label facei = c[cfi];
                    label patchi = patches.whichPatch(facei);
                    if (patchi == -1)
                    {
                        label nei =  mesh.faceNeighbour()[facei];
                        if (nei == celli)
                        {
                            nei = mesh.faceOwner()[facei];
                        }
                        if (layerCells[nei] > -1)
                        {
                            nLayerStacks++;
                        }
                    }
                    else if (patches[patchi].coupled())
                    {
                        label bFaceI = facei-mesh.nInternalFaces();
                        if (neiLayerCells[bFaceI] > -1)
                        {
                            nLayerStacks++;
                        }
                    }
                }

                if (nLayerStacks > 2)
                {
                    forAll(c, cfi)
                    {
                        usedPoints.insert(mesh.faces()[c[cfi]]);
                    }
                }
            }
        }
    }

    pointSet cornerPoints
    (
        mesh,
        "cornerPoints",
        usedPoints
    );
    cornerPoints.sync(mesh);


    scalarField cVolumes(mesh.nCells(), Zero);
    pointField cCentres(mesh.nCells(), Zero);

    boolList updatedCells(mesh.nCells(),false);

    //Reset marked pts from layer points to corner cell points
    markedPts = false;

    forAllConstIter(pointSet, cornerPoints, iter)
    {
        label pointi = iter.key();
        markedPts[pointi] = true;
        const labelList& pCells =
            mesh.pointCells()[pointi];
        forAll(pCells, pCI)
        {
            label celli = pCells[pCI];
            updatedCells[celli] = true;
        }
    }

    boolList layerEdges(mesh.nEdges(), false);
    forAll(mesh.edges(), edgei)
    {
        const edge& e = mesh.edges()[edgei];
        if
        (
            (boundaryPoints[e[0]] && !boundaryPoints[e[1]])
            || (boundaryPoints[e[1]] && !boundaryPoints[e[0]])
        )
        {
            layerEdges[edgei] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        layerEdges,
        orEqOp<bool>(),
        false
    );

    boolList exclFaces(pp.size(), false);
    edgeClassification eClass
    (
        mesh,
        newPoints,
        pp,
        meshEdges,
        exclFaces,
        0.707,
        0.707
    );
    const List<Tuple2<edgeClassification::edgeType,scalar>>&
        eType = eClass.edgeTypes();

    pointField surfNormals = eClass.calculatePointNormals
    (
        exclFaces,
        0,
        true
    );
    boolList featurePts(meshPoints.size(), false);
    forAll(pp.edges(), patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];
        label v0 = e[0];
        label v1 = e[1];

        if
        (
            eType[patchEdgeI].first() != edgeClassification::MANIFOLD
        )
        {
            featurePts[v0] = true;
            featurePts[v1] = true;
        }
    }


    for (label iter = 0; iter < 5; iter++)
    {
        avePt = vector::zero;
        aveWt= Zero;

        forAll(updatedCells, celli)
        {
            if (updatedCells[celli])
            {
                cCentres[celli] =
                    mesh.cells()[celli].centre(newPoints, mesh.faces());
                cVolumes[celli] =
                    mesh.cells()[celli].mag(newPoints, mesh.faces());
            }
        }
        forAllConstIter(pointSet, cornerPoints, iter)
        {
            label pointi = iter.key();
            const labelList& pCells =
                mesh.pointCells()[pointi];
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                scalar vol = mag(cVolumes[celli]);
                point cc = cCentres[celli];
                avePt[pointi] += cc*vol;
                aveWt[pointi] += vol;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            avePt,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            aveWt,
            plusEqOp<scalar>(),
            scalar(0)     // null value
        );

        forAll(mesh.points(), pointi)
        {
            if (aveWt[pointi] > SMALL)
            {
                newPoints[pointi] = 0.5*
                (
                    newPoints[pointi]
                    +(avePt[pointi]/aveWt[pointi])
                );
            }
        }
    }
    avePt = vector::zero;
    aveWt = 0;

    //Move non-feature surface points at corner cells
    forAll(featurePts, pti)
    {
        if (!featurePts[pti])
        {
            label meshpointi = meshPoints[pti];
            const labelList& pEdges = mesh.pointEdges()[meshpointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if (layerEdges[edgei])
                {
                    plane fPlane(newPoints[meshpointi], surfNormals[pti]);

                    const edge e = mesh.edges()[edgei];
                    label otherPt(e[0] == meshpointi ? e[1] : e[0]);
                    if (markedPts[otherPt])
                    {
                        point hPt = fPlane.nearestPoint(newPoints[otherPt]);

                        avePt[meshpointi] += hPt;
                        aveWt[meshpointi] += scalar(1);

                        avePt[otherPt] += 0.5*(hPt+newPoints[otherPt]);
                        aveWt[otherPt] += scalar(1);
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        avePt,
        plusEqOp<vector>(),
        vector::zero  // null value
    );

    syncTools::syncPointList
    (
        mesh,
        aveWt,
        plusEqOp<scalar>(),
        scalar(0)     // null value
    );

    forAll(mesh.points(), pointi)
    {
        if (aveWt[pointi] > SMALL)
        {
            newPoints[pointi] = 0.5*
            (
                newPoints[pointi] +(avePt[pointi]/aveWt[pointi])
            );
        }
    }

    return;
}



//************************************************************************* //
