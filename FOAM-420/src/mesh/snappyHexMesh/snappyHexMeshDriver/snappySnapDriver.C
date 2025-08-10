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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

Description
    All to do with snapping to the surface

\*----------------------------------------------------------------------------*/

#include "snappyHexMeshDriver/snappySnapDriver.H"
#include "snappyHexMeshDriver/snappyLayerDriver.H"
#include "motionSmoother/motionSmoother.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "surfaceFormats/obj/OBJstream.H"
#include "algorithms/PointEdgeWave/PointEdgeWave.H"
#include "snappyHexMeshDriver/snapParameters/snapParameters.H"
#include "regionSplit/localPointRegion.H"
#include "boundaryLayerRefinement/boundaryLayerRefinement.H"

#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif

#include "featureLineSnapping/featureLinePrep.H"
#include "featureLineSnapping/featureLineSnapper.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "snappyHexMeshDriver/keepData/keepData.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "algorithms/FaceCellWave/FaceCellWave.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "meshes/primitiveShapes/tetrahedron/tetrahedron.H"
#include "externalDisplacementMeshMover/medialAxisMeshMover.H"
#include "layerManipulate/layerManipulate.H"
#include "autoOptimize/autoOptimize.H"
#include "autoSplitCells/autoSplitCells.H"
#include "edgeClassification/edgeClassification.H"
#include <list>
#include "motionSmoother/polyMeshGeometry/polyMeshGeometry.H"
#include "sets/topoSets/faceSet.H"

#ifdef FOAM_USE_TBB
  #include "include/TBBTimer.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(snappySnapDriver, 0);

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Calculate geometrically collocated points, Requires PackedList to be
// sized and initalised!
Foam::label Foam::snappySnapDriver::getCollocatedPoints
(
    const scalar tol,
    const pointField& points,
    PackedBoolList& isCollocatedPoint
)
{
    labelList pointMap;
    label nUnique = mergePoints
    (
        points,                         // points
        tol,                            // mergeTol
        false,                          // verbose
        pointMap
    );
    bool hasMerged = (nUnique < points.size());

    if (!returnReduce(hasMerged, orOp<bool>()))
    {
        return 0;
    }

    // Determine which merged points are referenced more than once
    label nCollocated = 0;

    // Per old point the newPoint. Or -1 (not set yet) or -2 (already seen
    // twice)
    labelList firstOldPoint(nUnique, -1);
    forAll(pointMap, oldPointi)
    {
        label newPointi = pointMap[oldPointi];

        if (firstOldPoint[newPointi] == -1)
        {
            // First use of oldPointi. Store.
            firstOldPoint[newPointi] = oldPointi;
        }
        else if (firstOldPoint[newPointi] == -2)
        {
            // Third or more reference of oldPointi -> non-manifold
            isCollocatedPoint.set(oldPointi, 1u);
            nCollocated++;
        }
        else
        {
            // Second reference of oldPointi -> non-manifold
            isCollocatedPoint.set(firstOldPoint[newPointi], 1u);
            nCollocated++;

            isCollocatedPoint.set(oldPointi, 1u);
            nCollocated++;

            // Mark with special value to save checking next time round
            firstOldPoint[newPointi] = -2;
        }
    }
    return returnReduce(nCollocated, sumOp<label>());
}


Foam::tmp<Foam::pointField> Foam::snappySnapDriver::smoothInternalDisplacement
(
    const meshRefinement& meshRefiner,
    const motionSmoother& meshMover
)
{
    const indirectPrimitivePatch& pp = meshMover.patch();
    const polyMesh& mesh = meshMover.mesh();

    // Get neighbour refinement
    const hexRef8& cutter = meshRefiner.meshCutter();
    const labelList& cellLevel = cutter.cellLevel();


    // Get the faces on the boundary
    PackedBoolList isFront(mesh.nFaces());
    forAll(pp.addressing(), i)
    {
        isFront[pp.addressing()[i]] = true;
    }

    // Walk out from the surface a bit. Poor man's FaceCellWave.
    // Commented out for now - not sure if needed and if so how much
    //for (label iter = 0; iter < 2; iter++)
    //{
    //    PackedBoolList newIsFront(mesh.nFaces());
    //
    //    forAll(isFront, facei)
    //    {
    //        if (isFront[facei])
    //        {
    //            label own = mesh.faceOwner()[facei];
    //            const cell& ownFaces = mesh.cells()[own];
    //            forAll(ownFaces, i)
    //            {
    //                newIsFront[ownFaces[i]] = true;
    //            }
    //
    //            if (mesh.isInternalFace(facei))
    //            {
    //                label nei = mesh.faceNeighbour()[facei];
    //                const cell& neiFaces = mesh.cells()[nei];
    //                forAll(neiFaces, i)
    //                {
    //                    newIsFront[neiFaces[i]] = true;
    //                }
    //            }
    //        }
    //    }
    //
    //    syncTools::syncFaceList
    //    (
    //        mesh,
    //        newIsFront,
    //        orEqOp<unsigned int>()
    //    );
    //
    //    isFront = newIsFront;
    //}

    // Mark all points on faces
    //  - not on the boundary
    //  - inbetween differing refinement levels
    PackedBoolList isMovingPoint(mesh.nPoints());

    label nInterface = 0;

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        label ownLevel = cellLevel[mesh.faceOwner()[facei]];
        label neiLevel = cellLevel[mesh.faceNeighbour()[facei]];

        if (!isFront[facei] && ownLevel != neiLevel)
        {
            const face& f = mesh.faces()[facei];
            forAll(f, fp)
            {
                isMovingPoint[f[fp]] = true;
            }

            nInterface++;
        }
    }

    labelList neiCellLevel;
    syncTools::swapBoundaryCellList(mesh, cellLevel, neiCellLevel);

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        label ownLevel = cellLevel[mesh.faceOwner()[facei]];
        label neiLevel = neiCellLevel[facei-mesh.nInternalFaces()];

        if (!isFront[facei] && ownLevel != neiLevel)
        {
            const face& f = mesh.faces()[facei];
            forAll(f, fp)
            {
                isMovingPoint[f[fp]] = true;
            }

            nInterface++;
        }
    }

    if (debug)
    {
        reduce(nInterface, sumOp<label>());
        Info<< "Found " << nInterface << " faces out of "
            << mesh.globalData().nTotalFaces()
            << " inbetween refinement regions." << endl;
    }

    // Make sure that points that are coupled to a moving point are marked
    // as well
    syncTools::syncPointList(mesh, isMovingPoint, maxEqOp<unsigned int>(), 0);

    // Unmark any point on the boundary. If we're doing zero iterations of
    // face-cell wave we might have coupled points not being unmarked.
    forAll(pp.meshPoints(), pointi)
    {
        isMovingPoint[pp.meshPoints()[pointi]] = false;
    }

    // Make sure that points that are coupled to meshPoints but not on a patch
    // are unmarked as well
    syncTools::syncPointList(mesh, isMovingPoint, minEqOp<unsigned int>(), 1);


    // Calculate average of connected cells
    labelList nCells(mesh.nPoints(), 0);
    pointField sumLocation(mesh.nPoints(), Zero);

    forAll(isMovingPoint, pointi)
    {
        if (isMovingPoint[pointi])
        {
            const labelList& pCells = mesh.pointCells(pointi);

            forAll(pCells, i)
            {
                sumLocation[pointi] += mesh.cellCentres()[pCells[i]];
                nCells[pointi]++;
            }
        }
    }

    // Sum
    syncTools::syncPointList(mesh, nCells, plusEqOp<label>(), label(0));
    syncTools::syncPointList
    (
        mesh,
        sumLocation,
        plusEqOp<point>(),
        vector::zero
    );

    tmp<pointField> tdisplacement(new pointField(mesh.nPoints(), Zero));
    pointField& displacement = tdisplacement.ref();

    label nAdapted = 0;

    forAll(displacement, pointi)
    {
        if (nCells[pointi] > 0)
        {
            displacement[pointi] =
                sumLocation[pointi]/nCells[pointi]-mesh.points()[pointi];
            nAdapted++;
        }
    }

    reduce(nAdapted, sumOp<label>());
    Info<< "Smoothing " << nAdapted << " points inbetween refinement regions."
        << endl;

    return tdisplacement;
}


// Calculate displacement as average of patch points.
Foam::tmp<Foam::pointField> Foam::snappySnapDriver::smoothPatchDisplacement
(
    const motionSmoother& meshMover,
    const List<labelPair>& baffles,
    const bool preSmoothBaffles
)
{
    const indirectPrimitivePatch& pp = meshMover.patch();

    // Calculate geometrically non-manifold points on the patch to be moved.
    PackedBoolList nonManifoldPoint(pp.nPoints());
    label nNonManifoldPoints = getCollocatedPoints
    (
        SMALL,
        pp.localPoints(),
        nonManifoldPoint
    );
    Info<< "Found " << nNonManifoldPoints << " non-manifold point(s)."
        << endl;


    // Average points
    // ~~~~~~~~~~~~~~

    // We determine three points:
    // - average of (centres of) connected patch faces
    // - average of (centres of) connected internal mesh faces
    // - as fallback: centre of any connected cell
    // so we can do something moderately sensible for non/manifold points.

    // Note: the averages are calculated properly parallel. This is
    // necessary to get the points shared by processors correct.


    const labelListList& pointFaces = pp.pointFaces();
    const labelList& meshPoints = pp.meshPoints();
    const pointField& points = pp.points();
    const polyMesh& mesh = meshMover.mesh();

    // Get labels of faces to count (master of coupled faces and baffle pairs)
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    {
        forAll(baffles, i)
        {
            label f0 = baffles[i].first();
            label f1 = baffles[i].second();

            if (isMasterFace.get(f0))
            {
                // Make f1 a slave
                isMasterFace.unset(f1);
            }
            else if (isMasterFace.get(f1))
            {
                isMasterFace.unset(f0);
            }
            else
            {
                FatalErrorInFunction
                    << "Both sides of baffle consisting of faces " << f0
                    << " and " << f1 << " are already slave faces."
                    << abort(FatalError);
            }
       }
   }


    // Get average position of boundary face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vectorField avgBoundary(pointFaces.size(), Zero);
    labelList nBoundary(pointFaces.size(), 0);

    forAll(pointFaces, patchPointi)
    {
        const labelList& pFaces = pointFaces[patchPointi];

        forAll(pFaces, pfi)
        {
            label facei = pFaces[pfi];

            if (isMasterFace.get(pp.addressing()[facei]))
            {
                avgBoundary[patchPointi] += pp[facei].centre(points);
                nBoundary[patchPointi]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        avgBoundary,
        plusEqOp<point>(),  // combine op
        vector::zero        // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nBoundary,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    forAll(avgBoundary, i)
    {
        avgBoundary[i] /= nBoundary[i];
    }


    // Get average position of internal face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    vectorField avgInternal;
    labelList nInternal;
    {
        vectorField globalSum(mesh.nPoints(), Zero);
        labelList globalNum(mesh.nPoints(), 0);

        // Note: no use of pointFaces
        const faceList& faces = mesh.faces();

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            const face& f = faces[facei];
            const point& fc = mesh.faceCentres()[facei];

            forAll(f, fp)
            {
                globalSum[f[fp]] += fc;
                globalNum[f[fp]]++;
            }
        }

        // Count coupled faces as internal ones (but only once)
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(patches, patchi)
        {
            if
            (
                patches[patchi].coupled()
             && refCast<const coupledPolyPatch>(patches[patchi]).owner()
            )
            {
                const coupledPolyPatch& pp =
                    refCast<const coupledPolyPatch>(patches[patchi]);

                const vectorField::subField faceCentres = pp.faceCentres();

                forAll(pp, i)
                {
                    const face& f = pp[i];
                    const point& fc = faceCentres[i];

                    forAll(f, fp)
                    {
                        globalSum[f[fp]] += fc;
                        globalNum[f[fp]]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            globalSum,
            plusEqOp<vector>(), // combine op
            vector::zero        // null value
        );
        syncTools::syncPointList
        (
            mesh,
            globalNum,
            plusEqOp<label>(),  // combine op
            label(0)            // null value
        );

        avgInternal.setSize(meshPoints.size());
        nInternal.setSize(meshPoints.size());

        forAll(avgInternal, patchPointi)
        {
            label meshPointi = meshPoints[patchPointi];

            nInternal[patchPointi] = globalNum[meshPointi];

            if (nInternal[patchPointi] == 0)
            {
                avgInternal[patchPointi] = globalSum[meshPointi];
            }
            else
            {
                avgInternal[patchPointi] =
                    globalSum[meshPointi]
                  / nInternal[patchPointi];
            }
        }
    }


    // Precalculate any cell using mesh point (replacement of pointCells()[])
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList anyCell(mesh.nPoints(), -1);
    forAll(mesh.faceOwner(), facei)
    {
        label own = mesh.faceOwner()[facei];
        const face& f = mesh.faces()[facei];

        forAll(f, fp)
        {
            anyCell[f[fp]] = own;
        }
    }


    // Displacement to calculate.
    tmp<pointField> tpatchDisp(new pointField(meshPoints.size(), Zero));
    pointField& patchDisp = tpatchDisp.ref();

    forAll(pointFaces, i)
    {
        label meshPointi = meshPoints[i];
        const point& currentPos = pp.points()[meshPointi];

        // Now we have the two average points: avgBoundary and avgInternal
        // and how many boundary/internal faces connect to the point
        // (nBoundary, nInternal)
        // Do some blending between the two.
        // Note: the following section has some reasoning behind it but the
        // blending factors can be experimented with.

        point newPos;

        if (!nonManifoldPoint.get(i))
        {
            // Points that are manifold. Weight the internal and boundary
            // by their number of faces and blend with
            scalar internalBlend = 0.1;
            scalar blend = 0.75;

            point avgPos =
                (
                   internalBlend*nInternal[i]*avgInternal[i]
                  +(1-internalBlend)*nBoundary[i]*avgBoundary[i]
                )
              / (internalBlend*nInternal[i]+(1-internalBlend)*nBoundary[i]);

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else if (nInternal[i] == 0)
        {
            // Non-manifold without internal faces. Use any connected cell
            // as internal point instead. Use precalculated any cell to avoid
            // e.g. pointCells()[meshPointi][0]

            const point& cc = mesh.cellCentres()[anyCell[meshPointi]];

            scalar cellCBlend = 0.1;
            scalar blend = 0.75;

            point avgPos = (1-cellCBlend)*avgBoundary[i] + cellCBlend*cc;

            newPos = (1-blend)*avgPos + blend*currentPos;
        }
        else
        {
            // Non-manifold point with internal faces connected to them
            scalar internalBlend = 0.9;
            scalar blend = 0.75;

            point avgPos =
                internalBlend*avgInternal[i]
              + (1-internalBlend)*avgBoundary[i];

            newPos = (1-blend)*avgPos + blend*currentPos;
        }

        patchDisp[i] = newPos - currentPos;
    }

    //Detect baffle edges and smooth adjacent points
    if (preSmoothBaffles)
    {
        labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));
        //Set angle for detecting baffle edges
        scalar baffleCos = -Foam::cos(degToRad(315));
        boolList excludedFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            pp,
            meshEdges,
            excludedFaces,
            baffleCos,
            baffleCos
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        labelList pType(meshPoints.size(), -1);
        forAll(pp.edges(), edgeI)
        {
            if (eType[edgeI].first() == edgeClassification::CONVEX)
            {
                const edge& e = pp.edges()[edgeI];
                pType[e[0]] = 0;
                pType[e[1]] = 0;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            pType,
            maxEqOp<label>(),  // combine op
            label(-1)        // null value
        );

        forAll(pp.edges(), edgeI)
        {
            const edge& e = pp.edges()[edgeI];
            label n0 = e[0];
            label n1 = e[1];
            if (pType[n0] == 0 && pType[n1] == -1)
            {
                pType[n1] = 1;
            }
            else if (pType[n1] == 0 && pType[n0] == -1)
            {
                pType[n0] = 1;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            pType,
            maxEqOp<label>(),  // combine op
            label(-1)        // null value
        );

        pointField avePoint(meshPoints.size(),vector::zero);
        pointField aveNorm(meshPoints.size(),vector::zero);
        labelList wtSum(meshPoints.size(),0);

        forAll(pointFaces, i)
        {
            if (pType[i] == 1)
            {
                const labelList& pFaces = pp.pointFaces()[i];

                forAll(pFaces, pFI)
                {
                    aveNorm[i] += pp.faceAreas()[pFaces[pFI]];
                    avePoint[i] += pp.faceCentres()[pFaces[pFI]];
                    wtSum[i]++;
                }
            }
        }


        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            avePoint,
            plusEqOp<vector>(), // combine op
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            aveNorm,
            plusEqOp<vector>(), // combine op
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            wtSum,
            plusEqOp<label>(), // combine op
            label(0)        // null value
        );

        forAll(pointFaces, i)
        {
            if (wtSum[i] > 0 && mag(aveNorm[i]) > SMALL)
            {
                label meshPointi = meshPoints[i];
                point currentPos = pp.points()[meshPointi];

                point pPt = avePoint[i]/wtSum[i];
                point pNorm = aveNorm[i]/wtSum[i];

                plane pl(pPt, pNorm);
                patchDisp[i] = pl.nearestPoint(currentPos) - currentPos;
            }
        }
    }

    return tpatchDisp;
}
//XXXXXXX
//Foam::tmp<Foam::pointField> Foam::snappySnapDriver::avg
//(
//    const indirectPrimitivePatch& pp,
//    const pointField& localPoints
//)
//{
//    const labelListList& pointEdges = pp.pointEdges();
//    const edgeList& edges = pp.edges();
//
//    tmp<pointField> tavg(new pointField(pointEdges.size(), Zero));
//    pointField& avg = tavg();
//
//    forAll(pointEdges, verti)
//    {
//        vector& avgPos = avg[verti];
//
//        const labelList& pEdges = pointEdges[verti];
//
//        forAll(pEdges, myEdgei)
//        {
//            const edge& e = edges[pEdges[myEdgei]];
//
//            label otherVerti = e.otherVertex(verti);
//
//            avgPos += localPoints[otherVerti];
//        }
//
//        avgPos /= pEdges.size();
//    }
//    return tavg;
//}
//Foam::tmp<Foam::pointField>
//Foam::snappySnapDriver::smoothLambdaMuPatchDisplacement
//(
//    const motionSmoother& meshMover,
//    const List<labelPair>& baffles
//)
//{
//    const indirectPrimitivePatch& pp = meshMover.patch();
//    pointField newLocalPoints(pp.localPoints());
//
//    const label iters = 90;
//    const scalar lambda = 0.33;
//    const scalar mu = 0.34;
//
//    for (label iter = 0; iter < iters; iter++)
//    {
//        // Lambda
//        newLocalPoints =
//            (1 - lambda)*newLocalPoints
//          + lambda*avg(pp, newLocalPoints);
//
//        // Mu
//        newLocalPoints =
//            (1 + mu)*newLocalPoints
//          - mu*avg(pp, newLocalPoints);
//    }
//    return newLocalPoints-pp.localPoints();
//}
//XXXXXXX


Foam::tmp<Foam::scalarField> Foam::snappySnapDriver::edgePatchDist
(
    const pointMesh& pMesh,
    const indirectPrimitivePatch& pp
)
{
    const polyMesh& mesh = pMesh();

    // Set initial changed points to all the patch points
    List<pointEdgePoint> wallInfo(pp.nPoints());

    forAll(pp.localPoints(), ppi)
    {
        wallInfo[ppi] = pointEdgePoint(pp.localPoints()[ppi], 0.0);
    }

    // Current info on points
    List<pointEdgePoint> allPointInfo(mesh.nPoints());

    // Current info on edges
    List<pointEdgePoint> allEdgeInfo(mesh.nEdges());

    PointEdgeWave<pointEdgePoint> wallCalc
    (
        mesh,
        pp.meshPoints(),
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh.globalData().nTotalPoints()  // max iterations
    );

    // Copy edge values into scalarField
    tmp<scalarField> tedgeDist(new scalarField(mesh.nEdges()));
    scalarField& edgeDist = tedgeDist.ref();

    forAll(allEdgeInfo, edgei)
    {
        edgeDist[edgei] = Foam::sqrt(allEdgeInfo[edgei].distSqr());
    }

    return tedgeDist;
}


void Foam::snappySnapDriver::dumpMove
(
    const fileName& fName,
    const pointField& meshPts,
    const pointField& surfPts
)
{
    // Dump direction of growth into file
    Info<< "Dumping move direction to " << fName << endl;

    OFstream nearestStream(fName);

    label verti = 0;

    forAll(meshPts, pti)
    {
        meshTools::writeOBJ(nearestStream, meshPts[pti]);
        verti++;

        meshTools::writeOBJ(nearestStream, surfPts[pti]);
        verti++;

        nearestStream<< "l " << verti-1 << ' ' << verti << nl;
    }
}


// Check whether all displacement vectors point outwards of patch. Return true
// if so.
bool Foam::snappySnapDriver::outwardsDisplacement
(
    const indirectPrimitivePatch& pp,
    const vectorField& patchDisp
)
{
    const vectorField& faceNormals = pp.faceNormals();
    const labelListList& pointFaces = pp.pointFaces();

    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];

        vector disp(patchDisp[pointi]);

        scalar magDisp = mag(disp);

        if (magDisp > SMALL)
        {
            disp /= magDisp;

            bool outwards = meshTools::visNormal(disp, faceNormals, pFaces);

            if (!outwards)
            {
                Warning<< "Displacement " << patchDisp[pointi]
                    << " at mesh point " << pp.meshPoints()[pointi]
                    << " coord " << pp.points()[pp.meshPoints()[pointi]]
                    << " points through the surrounding patch faces" << endl;
                return false;
            }
        }
        else
        {
            //? Displacement small but in wrong direction. Would probably be ok.
        }
    }
    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::snappySnapDriver::snappySnapDriver
(
    meshRefinement& meshRefiner,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const dictionary& meshDict,
    const meshControl& controller,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
)
:
    meshRefiner_(meshRefiner),
    decomposer_(decomposer),
    distributor_(distributor),
    meshDict_(meshDict),
    controller_(controller),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::snappySnapDriver::calcSnapDistance
(
    const meshRefinement& meshRefiner,
    const snapParameters& snapParams,
    const indirectPrimitivePatch& pp
)
{
//#ifdef FOAM_USE_TBB
//    Timer timer("snappySnapDriver::calcSnapDistance");
//#endif

    const fvMesh& mesh = meshRefiner.mesh();

    // Undistorted edge length
    const scalar edge0Len = meshRefiner.meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshRefiner.meshCutter().cellLevel();

    const pointField& localPoints = pp.localPoints();
    scalarField maxEdgeLen(localPoints.size(), -GREAT);
    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        label own = mesh.faceOwner()[meshFaceI];
        label  level = cellLevel[own];
        scalar len = edge0Len / pow(2., level);
        const face& f = pp.localFaces()[i];
        forAll(f,fp)
        {
            label pointi = f[fp];
            maxEdgeLen[pointi] = max(maxEdgeLen[pointi], len);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxEdgeLen,
        maxEqOp<scalar>(),  // combine op
        -GREAT              // null value
    );

    return scalarField(snapParams.snapTol()*maxEdgeLen);
}


void Foam::snappySnapDriver::preSmoothPatch
(
    const meshRefinement& meshRefiner,
    const snapParameters& snapParams,
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(smooth, "snappyHexMesh::snap::smoothing");
    #endif

    const fvMesh& mesh = meshRefiner.mesh();

    labelList checkFaces(meshRefiner.selectInnerFaces());

    if (snapParams.nSmoothInternal() > 0)
    {
        Info<< "Smoothing patch and internal points ..." << endl;
    }
    else
    {
        Info<< "Smoothing patch points ..." << endl;
    }

    bool preSmoothBaffles = snapParams.preSmoothBaffles();

    vectorField& pointDisp = meshMover.pointDisplacement().primitiveFieldRef();

    for
    (
        label smoothIter = 0;
        smoothIter < snapParams.nSmoothPatch();
        smoothIter++
    )
    {
        Info<< "Smoothing iteration " << smoothIter << endl;

        // If enabled smooth the internal points
        if (snapParams.nSmoothInternal() > smoothIter)
        {
            // Override values on internal points on refinement interfaces
            pointDisp = smoothInternalDisplacement(meshRefiner, meshMover);
        }

        // Smooth the patch points
        pointField patchDisp
        (
            smoothPatchDisplacement
            (
                meshMover,
                baffles,
                preSmoothBaffles
            )
        );
        //pointField patchDisp
        //(
        //  smoothLambdaMuPatchDisplacement(meshMover, baffles)
        //);

        // Take over patch displacement as boundary condition on
        // pointDisplacement
        meshMover.setDisplacement(patchDisp);
        // Start off from current mesh.points()
        meshMover.correct();

        scalar oldErrorReduction = -1;

        for (label snapIter = 0; snapIter < 2*snapParams.nSnap(); snapIter++)
        {
            Info<< nl << "Scaling iteration " << snapIter << endl;

            if (snapIter == snapParams.nSnap())
            {
                Info<< "Displacement scaling for error reduction set to 0."
                    << endl;
                oldErrorReduction = meshMover.setErrorReduction(0.0);
            }

            // Try to adapt mesh to obtain displacement by smoothly
            // decreasing displacement at error locations.
            if (meshMover.scaleMesh(checkFaces, baffles, true, nInitErrors))
            {
                Info<< "Successfully moved mesh" << endl;
                break;
            }
        }

        if (oldErrorReduction >= 0)
        {
            meshMover.setErrorReduction(oldErrorReduction);
        }
        Info<< endl;
    }


    // The current mesh is the starting mesh to smooth from.
    meshMover.correct();

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing patch smoothed mesh to time "
            << meshRefiner.timeName() << '.' << endl;
        meshRefiner.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner.timeName()
        );
        Info<< "Dumped mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }

    Info<< "Patch points smoothed in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
}


// Get (pp-local) indices of points that are both on zone and on patched surface
Foam::labelList Foam::snappySnapDriver::getZoneSurfacePoints
(
    const fvMesh& mesh,
    const indirectPrimitivePatch& pp,
    const word& zoneName
)
{
    label zonei = mesh.faceZones().findZoneID(zoneName);

    if (zonei == -1)
    {
        FatalErrorInFunction
            << "Cannot find zone " << zoneName
            << exit(FatalError);
    }

    const faceZone& fZone = mesh.faceZones()[zonei];


    // Could use PrimitivePatch & localFaces to extract points but might just
    // as well do it ourselves.

    boolList pointOnZone(pp.nPoints(), false);

    forAll(fZone, i)
    {
        const face& f = mesh.faces()[fZone[i]];

        forAll(f, fp)
        {
            label meshPointi = f[fp];

            Map<label>::const_iterator iter =
                pp.meshPointMap().find(meshPointi);

            if (iter != pp.meshPointMap().end())
            {
                label pointi = iter();
                pointOnZone[pointi] = true;
            }
        }
    }

    return findIndices(pointOnZone, true);
}


Foam::tmp<Foam::pointField> Foam::snappySnapDriver::avgCellCentres
(
    const fvMesh& mesh,
    const indirectPrimitivePatch& pp
)
{
    const labelListList& pointFaces = pp.pointFaces();


    tmp<pointField> tavgBoundary
    (
        new pointField(pointFaces.size(), Zero)
    );
    pointField& avgBoundary = tavgBoundary.ref();
    labelList nBoundary(pointFaces.size(), 0);

    forAll(pointFaces, pointi)
    {
        const labelList& pFaces = pointFaces[pointi];

        forAll(pFaces, pfi)
        {
            label facei = pFaces[pfi];
            label meshFacei = pp.addressing()[facei];

            label own = mesh.faceOwner()[meshFacei];
            avgBoundary[pointi] += mesh.cellCentres()[own];
            nBoundary[pointi]++;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        avgBoundary,
        plusEqOp<point>(),  // combine op
        vector::zero        // null value
    );
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nBoundary,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    forAll(avgBoundary, i)
    {
        avgBoundary[i] /= nBoundary[i];
    }
    return tavgBoundary;
}


void Foam::snappySnapDriver::detectNearSurfaces
(
    const scalar planarCos,
    const indirectPrimitivePatch& pp,
    const pointField& nearestPoint,
    const vectorField& nearestNormal,

    vectorField& disp
) const
{
    Info<< "Detecting near surfaces ..." << endl;

    const pointField& localPoints = pp.localPoints();
    const labelList& meshPoints = pp.meshPoints();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const fvMesh& mesh = meshRefiner_.mesh();

    //// Get local edge length based on refinement level
    //const scalarField edgeLen(calcEdgeLen(pp));
    //
    //// Generate rays for every surface point
    //// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //
    //{
    //    const scalar cos45 = Foam::cos(degToRad(45));
    //    vector n(cos45, cos45, cos45);
    //    n /= mag(n);
    //
    //    pointField start(14*pp.nPoints());
    //    pointField end(start.size());
    //
    //    label rayi = 0;
    //    forAll(localPoints, pointi)
    //    {
    //        const point& pt = localPoints[pointi];
    //
    //        // Along coordinate axes
    //
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.y() -= edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.y() += edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.z() -= edgeLen[pointi];
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.z() += edgeLen[pointi];
    //        }
    //
    //        // At 45 degrees
    //
    //        const vector vec(edgeLen[pointi]*n);
    //
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() += vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() += vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() += vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //        {
    //            start[rayi] = pt;
    //            point& endPt = end[rayi++];
    //            endPt = pt;
    //            endPt.x() -= vec.x();
    //            endPt.y() -= vec.y();
    //            endPt.z() -= vec.z();
    //        }
    //    }
    //
    //    labelList surface1;
    //    List<pointIndexHit> hit1;
    //    labelList region1;
    //    vectorField normal1;
    //
    //    labelList surface2;
    //    List<pointIndexHit> hit2;
    //    labelList region2;
    //    vectorField normal2;
    //    surfaces.findNearestIntersection
    //    (
    //        unzonedSurfaces,    // surfacesToTest,
    //        start,
    //        end,
    //
    //        surface1,
    //        hit1,
    //        region1,
    //        normal1,
    //
    //        surface2,
    //        hit2,
    //        region2,
    //        normal2
    //    );
    //
    //    // All intersections
    //    {
    //        OBJstream str
    //        (
    //            mesh.time().path()
    //          / "surfaceHits_" + meshRefiner_.timeName() + ".obj"
    //        );
    //
    //        Info<< "Dumping intersections with rays to " << str.name()
    //            << endl;
    //
    //        forAll(hit1, i)
    //        {
    //            if (hit1[i].hit())
    //            {
    //                str.write(linePointRef(start[i], hit1[i].hitPoint()));
    //            }
    //            if (hit2[i].hit())
    //            {
    //                str.write(linePointRef(start[i], hit2[i].hitPoint()));
    //            }
    //        }
    //    }
    //
    //    // Co-planar intersections
    //    {
    //        OBJstream str
    //        (
    //            mesh.time().path()
    //          / "coplanarHits_" + meshRefiner_.timeName() + ".obj"
    //        );
    //
    //        Info<< "Dumping intersections with co-planar surfaces to "
    //            << str.name() << endl;
    //
    //        forAll(localPoints, pointi)
    //        {
    //            bool hasNormal = false;
    //            point surfPointA;
    //            vector surfNormalA;
    //            point surfPointB;
    //            vector surfNormalB;
    //
    //            bool isCoplanar = false;
    //
    //            label rayi = 14*pointi;
    //            for (label i = 0; i < 14; i++)
    //            {
    //                if (hit1[rayi].hit())
    //                {
    //                    const point& pt = hit1[rayi].hitPoint();
    //                    const vector& n = normal1[rayi];
    //
    //                    if (!hasNormal)
    //                    {
    //                        hasNormal = true;
    //                        surfPointA = pt;
    //                        surfNormalA = n;
    //                    }
    //                    else
    //                    {
    //                        if
    //                        (
    //                            meshRefiner_.isGap
    //                            (
    //                                planarCos,
    //                                surfPointA,
    //                                surfNormalA,
    //                                pt,
    //                                n
    //                            )
    //                        )
    //                        {
    //                            isCoplanar = true;
    //                            surfPointB = pt;
    //                            surfNormalB = n;
    //                            break;
    //                        }
    //                    }
    //                }
    //                if (hit2[rayi].hit())
    //                {
    //                    const point& pt = hit2[rayi].hitPoint();
    //                    const vector& n = normal2[rayi];
    //
    //                    if (!hasNormal)
    //                    {
    //                        hasNormal = true;
    //                        surfPointA = pt;
    //                        surfNormalA = n;
    //                    }
    //                    else
    //                    {
    //                        if
    //                        (
    //                            meshRefiner_.isGap
    //                            (
    //                                planarCos,
    //                                surfPointA,
    //                                surfNormalA,
    //                                pt,
    //                                n
    //                            )
    //                        )
    //                        {
    //                            isCoplanar = true;
    //                            surfPointB = pt;
    //                            surfNormalB = n;
    //                            break;
    //                        }
    //                    }
    //                }
    //
    //                rayi++;
    //            }
    //
    //            if (isCoplanar)
    //            {
    //                str.write(linePointRef(surfPointA, surfPointB));
    //            }
    //        }
    //    }
    //}


    const pointField avgCc(avgCellCentres(mesh, pp));

    // Construct rays through localPoints to beyond cell centre
    pointField start(pp.nPoints());
    pointField end(pp.nPoints());
    forAll(localPoints, pointi)
    {
        const point& pt = localPoints[pointi];
        const vector d = 2*(avgCc[pointi]-pt);
        start[pointi] = pt - d;
        end[pointi] = pt + d;
    }


    autoPtr<OBJstream> gapStr;
    if (debug&meshRefinement::ATTRACTION)
    {
        gapStr.reset
        (
            new OBJstream
            (
                mesh.time().path()
              / "detectNearSurfaces_" + meshRefiner_.timeName() + ".obj"
            )
        );
    }


    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh,
            meshPoints
        )
    );

    label nOverride = 0;

    // 1. All points to non-interface surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    {
        const labelList unzonedSurfaces =
            surfaceZonesInfo::getUnnamedSurfaces
            (
                meshRefiner_.surfaces().surfZones()
            );

        // Do intersection test
        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;
        vectorField normal2;
        surfaces.findNearestIntersection
        (
            unzonedSurfaces,
            start,
            end,

            surface1,
            hit1,
            region1,
            normal1,

            surface2,
            hit2,
            region2,
            normal2
        );


        forAll(localPoints, pointi)
        {
            // Current location
            const point& pt = localPoints[pointi];

            bool override = false;

            //if (hit1[pointi].hit())
            //{
            //    if
            //    (
            //        meshRefiner_.isGap
            //        (
            //            planarCos,
            //            nearestPoint[pointi],
            //            nearestNormal[pointi],
            //            hit1[pointi].hitPoint(),
            //            normal1[pointi]
            //        )
            //    )
            //    {
            //        disp[pointi] = hit1[pointi].hitPoint()-pt;
            //        override = true;
            //    }
            //}
            //if (hit2[pointi].hit())
            //{
            //    if
            //    (
            //        meshRefiner_.isGap
            //        (
            //            planarCos,
            //            nearestPoint[pointi],
            //            nearestNormal[pointi],
            //            hit2[pointi].hitPoint(),
            //            normal2[pointi]
            //        )
            //    )
            //    {
            //        disp[pointi] = hit2[pointi].hitPoint()-pt;
            //        override = true;
            //    }
            //}

            if (hit1[pointi].hit() && hit2[pointi].hit())
            {
                if
                (
                    meshRefiner_.isGap
                    (
                        planarCos,
                        hit1[pointi].hitPoint(),
                        normal1[pointi],
                        hit2[pointi].hitPoint(),
                        normal2[pointi]
                    )
                )
                {
                    // TBD: check if the attraction (to nearest) would attract
                    // good enough and not override attraction

                    if (gapStr.valid())
                    {
                        const point& intPt = hit2[pointi].hitPoint();
                        gapStr().write(linePointRef(pt, intPt));
                    }

                    // Choose hit2 : nearest to end point (so inside the domain)
                    disp[pointi] = hit2[pointi].hitPoint()-pt;
                    override = true;
                }
            }

            if (override && isPatchMasterPoint[pointi])
            {
                nOverride++;
            }
        }
    }


    // 2. All points on zones to their respective surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Surfaces with zone information
        const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

        const labelList zonedSurfaces = surfaceZonesInfo::getNamedSurfaces
        (
            surfZones
        );

        forAll(zonedSurfaces, i)
        {
            label zoneSurfi = zonedSurfaces[i];

            const word& faceZoneName = surfZones[zoneSurfi].faceZoneName();

            const labelList surfacesToTest(1, zoneSurfi);

            // Get indices of points both on faceZone and on pp.
            labelList zonePointIndices
            (
                getZoneSurfacePoints
                (
                    mesh,
                    pp,
                    faceZoneName
                )
            );

            // Do intersection test
            labelList surface1;
            List<pointIndexHit> hit1;
            labelList region1;
            vectorField normal1;

            labelList surface2;
            List<pointIndexHit> hit2;
            labelList region2;
            vectorField normal2;
            surfaces.findNearestIntersection
            (
                surfacesToTest,
                pointField(start, zonePointIndices),
                pointField(end, zonePointIndices),

                surface1,
                hit1,
                region1,
                normal1,

                surface2,
                hit2,
                region2,
                normal2
            );


            forAll(hit1, i)
            {
                label pointi = zonePointIndices[i];

                // Current location
                const point& pt = localPoints[pointi];

                bool override = false;

                //if (hit1[i].hit())
                //{
                //    if
                //    (
                //        meshRefiner_.isGap
                //        (
                //            planarCos,
                //            nearestPoint[pointi],
                //            nearestNormal[pointi],
                //            hit1[i].hitPoint(),
                //            normal1[i]
                //        )
                //    )
                //    {
                //        disp[pointi] = hit1[i].hitPoint()-pt;
                //        override = true;
                //    }
                //}
                //if (hit2[i].hit())
                //{
                //    if
                //    (
                //        meshRefiner_.isGap
                //        (
                //            planarCos,
                //            nearestPoint[pointi],
                //            nearestNormal[pointi],
                //            hit2[i].hitPoint(),
                //            normal2[i]
                //        )
                //    )
                //    {
                //        disp[pointi] = hit2[i].hitPoint()-pt;
                //        override = true;
                //    }
                //}

                if (hit1[i].hit() && hit2[i].hit())
                {
                    if
                    (
                        meshRefiner_.isGap
                        (
                            planarCos,
                            hit1[i].hitPoint(),
                            normal1[i],
                            hit2[i].hitPoint(),
                            normal2[i]
                        )
                    )
                    {
                        if (gapStr.valid())
                        {
                            const point& intPt = hit2[i].hitPoint();
                            gapStr().write(linePointRef(pt, intPt));
                        }

                        disp[pointi] = hit2[i].hitPoint()-pt;
                        override = true;
                    }
                }

                if (override && isPatchMasterPoint[pointi])
                {
                    nOverride++;
                }
            }
        }
    }

    Info<< "Overriding nearest with intersection of close gaps at "
        << returnReduce(nOverride, sumOp<label>())
        << " out of " << returnReduce(pp.nPoints(), sumOp<label>())
        << " points." << endl;
} // detectNearSurfaces


void Foam::snappySnapDriver::calcNearestSurface
(
    const refinementSurfaces& surfaces,

    const labelList& surfacesToTest,
    const labelListList& regionsToTest,

    const pointField& localPoints,
    const labelList& zonePointIndices,

    scalarField& minSnapDist,
    labelList& snapSurf,
    vectorField& patchDisp,

    // Optional: nearest point, normal
    pointField& nearestPoint,
    vectorField& nearestNormal
)
{
//#ifdef FOAM_USE_TBB
//    Timer timer("snappySnapDriver::calcNearestSurface 1");
//#endif
    // Find nearest for points both on faceZone and pp.
    List<pointIndexHit> hitInfo;
    labelList hitSurface;

    if (nearestNormal.size() == localPoints.size())
    {
        labelList hitRegion;
        vectorField hitNormal;
        surfaces.findNearestRegion
        (
            surfacesToTest,
            regionsToTest,

            pointField(localPoints, zonePointIndices),
            sqr(scalarField(minSnapDist, zonePointIndices)),

            hitSurface,
            hitInfo,
            hitRegion,
            hitNormal
        );

        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                label pointi = zonePointIndices[i];
                nearestPoint[pointi] = hitInfo[i].hitPoint();
                nearestNormal[pointi] = hitNormal[i];
            }
        }
    }
    else
    {
        surfaces.findNearest
        (
            surfacesToTest,
            regionsToTest,

            pointField(localPoints, zonePointIndices),
            sqr(scalarField(minSnapDist, zonePointIndices)),

            hitSurface,
            hitInfo
        );
    }

    forAll(hitInfo, i)
    {
        if (hitInfo[i].hit())
        {
            label pointi = zonePointIndices[i];

            patchDisp[pointi] = hitInfo[i].hitPoint() - localPoints[pointi];
            minSnapDist[pointi] = mag(patchDisp[pointi]);
            snapSurf[pointi] = hitSurface[i];
        }
    }
} // void calcNearestSurface


Foam::vectorField Foam::snappySnapDriver::calcNearestSurface
(
    const scalarField& snapDist,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
//#ifdef FOAM_USE_TBB
//    Timer timer("snappySnapDriver::calcNearestSurface 2");
//#endif
    Info<< "Calculating patchDisplacement as distance to nearest surface"
        << " point ..." << endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const indirectPrimitivePatch& pp = meshMover.patch();
    const pointField& localPoints = pp.localPoints();

    // Divide surfaces into zoned and unzoned
    const labelList
        allSurfaces(identity(surfaces.surfaces().size()));

    // Displacement per patch point
    vectorField patchDisp(localPoints.size(), vector::zero);
    vectorField vecToSurface(pp.nPoints(), vector(GREAT, GREAT, GREAT));
    vectorField pointCellCentres(pp.nPoints(), vector::zero);
    labelList nPointFaces(pp.nPoints(), 0);
    scalar tol = 1e-10;

    forAll(localPoints, pointI)
    {
        const labelList& pFaces = pp.pointFaces()[pointI];

        forAll(pFaces, faceI)
        {
            const label& own = mesh.faceOwner()[pp.addressing()[pFaces[faceI]]];
            pointCellCentres[pointI] += mesh.cellCentres()[own];
            nPointFaces[pointI]++;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        pointCellCentres,
        plusEqOp<vector>(),
        vector::zero       // null value
    );

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        nPointFaces,
        plusEqOp<label>(),
        label(0)       // null value
    );

    forAll(pointCellCentres, i)
    {
        pointCellCentres[i] /= nPointFaces[i];
    }

    const vectorField& faceNormals = pp.faceNormals();

    pointField pointNormals(pp.nPoints(), vector::zero);
    {
        labelList nPointFaces(pp.nPoints(), 0);

        forAll(faceNormals, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            forAll(f, fp)
            {
                pointNormals[f[fp]] += faceNormals[faceI];
                nPointFaces[f[fp]] ++;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            pointNormals,
            plusEqOp<vector>(),
            vector::zero       // null value
        );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            nPointFaces,
            plusEqOp<label>(),
            label(0)            // null value
        );
        forAll(pointNormals, i)
        {
            pointNormals[i] /= nPointFaces[i];
        }
    }

    pointNormals /= (mag(pointNormals) + SMALL);

    List<pointIndexHit> hitInfo;
    List<pointIndexHit> hit1;

    {
       labelList surface1;
       labelList region1;
       labelList surface2;
       List<pointIndexHit> hit2;
       labelList region2;

       surfaces.findNearestIntersection
       (
          allSurfaces,
          pointCellCentres,
          localPoints,

          surface1,
          hit1,
          region1,

          surface2,
          hit2,
          region2
       );

       labelList hitSurface;
       surfaces.findNearest
       (
          allSurfaces,
          localPoints,
          sqr(snapDist),        // sqr of attract distance
          hitSurface,
          hitInfo
       );

       forAll(localPoints, i)
       {
           if (hitInfo[i].hit())
           {
               patchDisp[i] = hitInfo[i].hitPoint() - localPoints[i];
           }
       }
    }

    {
       pointField endPoints(localPoints.size());
       pointField endPointsReverse(localPoints.size());
       pointField startPoints(localPoints.size());
       labelList pointIndex(localPoints.size());

       label testI = 0;
       forAll(localPoints, i)
       {
           if (hitInfo[i].hit())
           {
               scalar hitDist = 2.*mag(hitInfo[i].hitPoint() - localPoints[i]);

               if
               (
                   hitDist > tol && !hit1[i].hit()
               )
               {
                   endPoints[testI] =
                       localPoints[i] + (hitDist*pointNormals[i]);
                   endPointsReverse[testI] =
                       localPoints[i] - (hitDist*pointNormals[i]);
                   startPoints[testI] = localPoints[i];
                   pointIndex[testI] = i;
                   testI++;
               }
           }
       }

       startPoints.setSize(testI);
       endPoints.setSize(testI);
       endPointsReverse.setSize(testI);
       pointIndex.setSize(testI);

       labelList surface1;
       List<pointIndexHit> hitPoint1;
       labelList pointRegion1;
       labelList surface2;
       List<pointIndexHit> hitPoint2;
       labelList pointRegion2;

       surfaces.findNearestIntersection
       (
          allSurfaces,
          startPoints,
          endPoints,

          surface1,
          hitPoint1,
          pointRegion1,

          surface2,
          hitPoint2,
          pointRegion2
       );

       labelList surface1r;
       List<pointIndexHit> hitPoint1r;
       labelList pointRegion1r;
       labelList surface2r;
       List<pointIndexHit> hitPoint2r;
       labelList pointRegion2r;

       surfaces.findNearestIntersection
       (
          allSurfaces,
          startPoints,
          endPointsReverse,

          surface1r,
          hitPoint1r,
          pointRegion1r,

          surface2r,
          hitPoint2r,
          pointRegion2r
       );

       forAll(hitPoint1, ip)
       {
          const label i = pointIndex[ip];
          if (hitPoint1[ip].hit())
          {
              patchDisp[i] = hitPoint1[ip].hitPoint()
                   - localPoints[i];
          }
          else if (hitPoint1r[ip].hit())
          {
              patchDisp[i] = hitPoint1r[ip].hitPoint()
                  - localPoints[i];
          }
       }
    }

    {
        label displacement = returnReduce(patchDisp.size(), sumOp<label>());
        if (displacement)
        {
            scalarField magDisp(mag(patchDisp));
            auto[mSum, mMin, mMax] = gFuncs<GSum, GMin, GMax>(magDisp);

            Info<< "Wanted displacement : average:"
                << mSum / displacement
                << " min:" << mMin
                << " max:" << mMax << endl;
        }
    }

    Info<< "Calculated surface displacement in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;


    // Limit amount of movement.
    forAll(patchDisp, patchPointI)
    {
        scalar magDisp = mag(patchDisp[patchPointI]);

        if (magDisp > snapDist[patchPointI])
        {
            patchDisp[patchPointI] *= snapDist[patchPointI] / magDisp;

            Pout<< "Limiting displacement for " << patchPointI
                << " from " << magDisp << " to " << snapDist[patchPointI]
                << " at location: "<< localPoints[patchPointI]
                << endl;
        }
    }
    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value
    );

    // Check for displacement being outwards.
    outwardsDisplacement(pp, patchDisp);

    // Set initial distribution of displacement field (on patches) from
    // patchDisp and make displacement consistent with b.c. on displacement
    // pointVectorField.
    meshMover.setDisplacement(patchDisp);

    if (debug)
    {
        dumpMove
        (
            mesh.time().path()/"patchDisplacement.obj",
            pp.localPoints(),
            pp.localPoints() + patchDisp
        );
    }

    return patchDisp;
} // vectorField calcNearestSurface


Foam::vectorField Foam::snappySnapDriver::calcNearestSurface
(
    const bool strictRegionSnap,
    const meshRefinement& meshRefiner,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const scalarField& snapDist,
    const indirectPrimitivePatch& pp,
    pointField& nearestPoint,
    vectorField& nearestNormal
)
{
//#ifdef FOAM_USE_TBB
//    Timer timer("snappySnapDriver::calcNearestSurface 3");
//#endif
    Info<< "Calculating patchDisplacement as distance to nearest surface"
        << " point ..." << endl;
    if (strictRegionSnap)
    {
        Info<< "    non-zone points : attract to local region on surface only"
            << nl
            << "    zone points     : attract to local region on surface only"
            << nl
            << endl;
    }
    else
    {
        Info<< "    non-zone points :"
            << " attract to nearest of all non-zone surfaces"
            << nl
            << "    zone points     : attract to zone surface only" << nl
            << endl;
    }


    const pointField& localPoints = pp.localPoints();
    const refinementSurfaces& surfaces = meshRefiner.surfaces();
    const fvMesh& mesh = meshRefiner.mesh();

    // Displacement per patch point
    vectorField patchDisp(localPoints.size(), Zero);

    if (returnReduce(localPoints.size(), sumOp<label>()) > 0)
    {
        // Current surface snapped to. Used to check whether points have been
        // snapped at all
        labelList snapSurf(localPoints.size(), -1);

        // Current best snap distance (since point might be on multiple
        // regions)
        scalarField minSnapDist(snapDist);


        if (strictRegionSnap)
        {
            // Attract patch points to same region only

            forAll(surfaces.surfaces(), surfi)
            {
                label geomi = surfaces.surfaces()[surfi];
                label nRegions = surfaces.geometry()[geomi].regions().size();

                const labelList surfacesToTest(1, surfi);

                for (label regioni = 0; regioni < nRegions; regioni++)
                {
                    label globali = surfaces.globalRegion(surfi, regioni);
                    label masterPatchi = globalToMasterPatch[globali];

                    // Get indices of points both on patch and on pp
                    labelList zonePointIndices
                    (
                        getFacePoints
                        (
                            pp,
                            mesh.boundaryMesh()[masterPatchi]
                        )
                    );

                    calcNearestSurface
                    (
                        surfaces,

                        surfacesToTest,
                        labelListList(1, labelList(1, regioni)), //regionsToTest

                        localPoints,
                        zonePointIndices,

                        minSnapDist,
                        snapSurf,
                        patchDisp,

                        // Optional: nearest point, normal
                        nearestPoint,
                        nearestNormal
                    );

                    if (globalToSlavePatch[globali] != masterPatchi)
                    {
                        label slavePatchi = globalToSlavePatch[globali];

                        // Get indices of points both on patch and on pp
                        labelList zonePointIndices
                        (
                            getFacePoints
                            (
                                pp,
                                mesh.boundaryMesh()[slavePatchi]
                            )
                        );

                        calcNearestSurface
                        (
                            surfaces,

                            surfacesToTest,
                            labelListList(1, labelList(1, regioni)),

                            localPoints,
                            zonePointIndices,

                            minSnapDist,
                            snapSurf,
                            patchDisp,

                            // Optional: nearest point, normal
                            nearestPoint,
                            nearestNormal
                        );
                    }
                }
            }
        }
        else
        {
            // Divide surfaces into zoned and unzoned
            const labelList unzonedSurfaces =
                surfaceZonesInfo::getUnnamedSurfaces
                (
                    meshRefiner.surfaces().surfZones()
                );


            // 1. All points to non-interface surfaces
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            List<pointIndexHit> hitInfo;
            labelList hitSurface;

            if (nearestNormal.size() == localPoints.size())
            {
                labelList hitRegion;
                vectorField hitNormal;
                surfaces.findNearestRegion
                (
                    unzonedSurfaces,
                    localPoints,
                    sqr(snapDist),
                    hitSurface,
                    hitInfo,
                    hitRegion,
                    hitNormal
                );

                forAll(hitInfo, pointi)
                {
                    if (hitInfo[pointi].hit())
                    {
                        nearestPoint[pointi] = hitInfo[pointi].hitPoint();
                        nearestNormal[pointi] = hitNormal[pointi];
                    }
                }
            }
            else
            {
                surfaces.findNearest
                (
                    unzonedSurfaces,
                    localPoints,
                    sqr(snapDist),        // sqr of attract distance
                    hitSurface,
                    hitInfo
                );
            }

            forAll(hitInfo, pointi)
            {
                if (hitInfo[pointi].hit())
                {
                    patchDisp[pointi] =
                        hitInfo[pointi].hitPoint()
                      - localPoints[pointi];

                    snapSurf[pointi] = hitSurface[pointi];
                }
            }


            const labelList zonedSurfaces =
                surfaceZonesInfo::getNamedSurfaces
                (
                    meshRefiner.surfaces().surfZones()
                );


            // 2. All points on zones to their respective surface
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            // Surfaces with zone information
            const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

            forAll(zonedSurfaces, i)
            {
                label surfi = zonedSurfaces[i];

                const word& faceZoneName = surfZones[surfi].faceZoneName();

                const labelList surfacesToTest(1, surfi);

                label geomi = surfaces.surfaces()[surfi];
                label nRegions = surfaces.geometry()[geomi].regions().size();


                // Get indices of points both on faceZone and on pp.
                labelList zonePointIndices
                (
                    getZoneSurfacePoints
                    (
                        mesh,
                        pp,
                        faceZoneName
                    )
                );


                calcNearestSurface
                (
                    surfaces,

                    surfacesToTest,
                    labelListList(1, identity(nRegions)),

                    localPoints,
                    zonePointIndices,

                    minSnapDist,
                    snapSurf,
                    patchDisp,

                    // Optional: nearest point, normal
                    nearestPoint,
                    nearestNormal
                );
            }
        }


        // Check if all points are being snapped
        forAll(snapSurf, pointi)
        {
            if (snapSurf[pointi] == -1)
            {
                WarningInFunction
                    << "For point:" << pointi
                    << " coordinate:" << localPoints[pointi]
                    << " did not find any surface within:"
                    << minSnapDist[pointi]
                    << " metre." << endl;
            }
        }

        {
            const PackedBoolList isPatchMasterPoint
            (
                meshRefinement::getMasterPoints
                (
                    mesh,
                    pp.meshPoints()
                )
            );

            scalarField magDisp(mag(patchDisp));

            auto[mMin, mMax] = gFuncs<GMin, GMax>(magDisp);
            Info<< "Wanted displacement : average:"
                <<  meshRefinement::gAverage(isPatchMasterPoint, magDisp)
                << " min:" << mMin
                << " max:" << mMax << endl;
        }
    }

    Info<< "Calculated surface displacement in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;


    // Limit amount of movement. Can not happen for triSurfaceMesh but
    // can happen for some analytical shapes?
    forAll(patchDisp, patchPointi)
    {
        scalar magDisp = mag(patchDisp[patchPointi]);

        if (magDisp > snapDist[patchPointi])
        {
            patchDisp[patchPointi] *= snapDist[patchPointi] / magDisp;

            Pout<< "Limiting displacement for " << patchPointi
                << " from " << magDisp << " to " << snapDist[patchPointi]
                << endl;
        }
    }

    // Points on zones in one domain but only present as point on other
    // will not do condition 2 on all. Sync explicitly.
    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        patchDisp,
        minMagSqrEqOp<point>(),         // combine op
        vector(GREAT, GREAT, GREAT)     // null value
    );

    return patchDisp;
} // vectorField calcNearestSurface


void Foam::snappySnapDriver::smoothDisplacement
(
    const snapParameters& snapParams,
    motionSmoother& meshMover
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const indirectPrimitivePatch& pp = meshMover.patch();

    Info<< "Smoothing displacement ..." << endl;

    // Set edge diffusivity as inverse of distance to patch
    scalarField edgeGamma(1.0/(edgePatchDist(meshMover.pMesh(), pp) + SMALL));
    //scalarField edgeGamma(mesh.nEdges(), 1.0);
    //scalarField edgeGamma(wallGamma(mesh, pp, 10, 1));

    // Get displacement field
    pointVectorField& disp = meshMover.displacement();

    for (label iter = 0; iter < snapParams.nSmoothDispl(); iter++)
    {
        pointVectorField oldDisp(disp);
        meshMover.smooth(oldDisp, edgeGamma, disp);
    }
    Info<< "Displacement smoothed in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;
        Info<< "Writing smoothed mesh to time " << meshRefiner_.timeName()
            << endl;

        // Moving mesh creates meshPhi. Can be cleared out by a mesh.clearOut
        // but this will also delete all pointMesh but not pointFields which
        // gives an illegal situation.

        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner_.timeName()
        );
        Info<< "Writing displacement field ..." << endl;
        disp.write();
        tmp<pointScalarField> magDisp(mag(disp));
        magDisp().write();

        Info<< "Writing actual patch displacement ..." << endl;
        vectorField actualPatchDisp(disp, pp.meshPoints());
        dumpMove
        (
            mesh.time().path()
          / "actualPatchDisplacement_" + meshRefiner_.timeName() + ".obj",
            pp.localPoints(),
            pp.localPoints() + actualPatchDisp
        );
    }
}


/**
 * Looks for quads in the feature patch that are being forced non-flat,
 * and splits them. The idea is this should help snapping where faces
 * straddle corners in the mesh.
 */
void Foam::snappySnapDriver::splitNonplanarQuads
(
    const indirectPrimitivePatch &pp,
    const pointField &newPoints,
    const vectorField &surfHitPatchPointNormals
) const
{
  //1. Detect which faces need to be split

  //Find faces which would snap much better if they were split -
  //i.e. if the face normals would be much closer to being psrallel to
  //the direction the respective adjacent nodes want to go in.

    const scalar minFlatness = 0.8; //somewhat arbitrary magic number
    const scalar goodAlignment = 0.8; //ditto

    DynamicList<labelledTri> splitFaces;
    // successive elements are split face pairs
    // label splitFaceRegionNumber = 99;
    // labelledTris require a "region number".  Might as well make one up.

    // loop over patch faces
    forAll(pp,iface) {

    const face &f = pp[iface];

        if (f.size() == 4) {

        //Try the two possible splits
        //0 1 2 / 2 3 0

        labelList tri1Nodes(3);
        tri1Nodes[0] = f[0];
        tri1Nodes[1] = f[1];
        tri1Nodes[2] = f[2];
        vector normal1 = face(tri1Nodes).unitNormal(newPoints);

        labelList tri2Nodes(3);
        tri2Nodes[0] = f[2];
        tri2Nodes[1] = f[3];
        tri2Nodes[2] = f[0];
        vector normal2 = face(tri2Nodes).unitNormal(newPoints);

        const scalar flatness12 = normal1 & normal2;

        // 1 2 3 / 3 0 1

        labelList tri3Nodes(3);
        tri3Nodes[0] = f[1];
        tri3Nodes[1] = f[2];
        tri3Nodes[2] = f[3];
        vector normal3 = face(tri3Nodes).unitNormal(newPoints);

        labelList tri4Nodes(3);
        tri4Nodes[0] = f[3];
        tri4Nodes[1] = f[0];
        tri4Nodes[2] = f[1];
        vector normal4 = face(tri4Nodes).unitNormal(newPoints);

        const scalar flatness34 = normal3 & normal4;


        if (flatness12 < minFlatness || flatness34 < minFlatness)
        {
            //split whichever way has less flatness *if* doing so makes
            //the new faces' normals close to the opposite nodes' hit
            //point normals
            if (flatness12 < flatness34)
            {
                scalar alignment1 =
                    -normal1 &
                    surfHitPatchPointNormals[pp.meshPointMap()[f[1]]];
                scalar alignment2 =
                    -normal2
                    & surfHitPatchPointNormals[pp.meshPointMap()[f[3]]];

                if
                (
                    alignment1 > goodAlignment
                    && alignment2 > goodAlignment
                 )
                {
                    Info<< "Alignment after splitting: " << alignment1
                         << "," << alignment2 << endl;
                    splitFaces.append
                        (
                            labelledTri
                            (
                                tri1Nodes[0],
                                tri1Nodes[1],
                                tri1Nodes[2],
                                int(flatness12)*100000
                            )
                        );
                    splitFaces.append
                        (
                            labelledTri
                            (
                                tri2Nodes[0],
                                tri2Nodes[1],
                                tri2Nodes[2],
                                int(flatness12)*100000
                            )
                        );
                }
            }
            else
            {
                scalar alignment3 =
                    -normal3
                    & surfHitPatchPointNormals[pp.meshPointMap()[f[2]]];
                scalar alignment4 =
                    -normal4
                    & surfHitPatchPointNormals[pp.meshPointMap()[f[0]]];

                if
                (
                    alignment3 > goodAlignment
                    && alignment4 > goodAlignment
                 )
                {
                    Info<< "Alignment after splitting: " << alignment3
                         << "," << alignment4 << endl;
                    splitFaces.append
                        (
                            labelledTri
                            (
                                tri3Nodes[0],
                                tri3Nodes[1],
                                tri3Nodes[2],
                                int(flatness34)*100000
                            )
                        );
                    splitFaces.append
                        (
                            labelledTri
                            (
                                tri4Nodes[0],
                                tri4Nodes[1],
                                tri4Nodes[2],
                                int(flatness34)*100000
                            )
                        );
                }
            }

        } // if not flat enough

        } //if a quad
    } // loop over faces in patch


    //convert dynamic list to list and use to instantiate a triSurface object
    List<labelledTri> faces(splitFaces.size());
    forAll(splitFaces,faceI) {
    faces[faceI] = splitFaces[faceI];
    }

    triSurface visSurface(faces,newPoints);

    visSurface.write("splitFaces.vtk",false);



  //2a. Split them in the main mesh

  //2b. Regenerate the patch.  (Question: will this mess up any arrays
  //in use which are based on edges or faces in the mesh or patch?

}


bool Foam::snappySnapDriver::scaleMesh
(
    const snapParameters& snapParams,
    const label nInitErrors,
    const List<labelPair>& baffles,
    motionSmoother& meshMover
)
{
    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(scale, "snappyHexMesh::snap::scale");
    #endif

    const fvMesh& mesh = meshRefiner_.mesh();

    // Relax displacement until correct mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList checkFaces(meshRefiner_.selectInnerFaces());

    scalar oldErrorReduction = -1;

    bool meshOk = false;

    Info<< "Moving mesh ..." << endl;
    for (label iter = 0; iter < 3*snapParams.nSnap(); iter++)
    {
        Info<< nl << "Iteration " << iter << endl;

        if (iter == snapParams.nSnap())
        {
            Info<< "Displacement scaling for error reduction set to 0." << endl;
            oldErrorReduction = meshMover.setErrorReduction(0.0);
        }

        meshOk = meshMover.scaleMesh(checkFaces, baffles, /*smoothInternalMesh*/true, nInitErrors);

        if (meshOk)
        {
            Info<< "Successfully moved mesh" << endl;

            break;
        }
        if (debug&meshRefinement::MESH)
        {
            const_cast<Time&>(mesh.time())++;
            Info<< "Writing scaled mesh to time " << meshRefiner_.timeName()
                << endl;
            mesh.write();

            Info<< "Writing displacement field ..." << endl;
            meshMover.displacement().write();
            tmp<pointScalarField> magDisp(mag(meshMover.displacement()));
            magDisp().write();
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover.setErrorReduction(oldErrorReduction);
    }
    Info<< "Moved mesh in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;

    return meshOk;
}


void Foam::snappySnapDriver::updateBaffles
(
    const mapPolyMesh &map,
    List<labelPair>& baffles
)
{
    DynamicList<labelPair> newBaffles(baffles.size());
    forAll(baffles, i)
    {
        label& face0 = baffles[i].first();
        label& face1 = baffles[i].second();
        face0 = map.reverseFaceMap()[face0];
        face1 = map.reverseFaceMap()[face1];
        //interior faces can be converted into boundary
        //faces when using cell removal
        if (face0 >= 0 && face1 >= 0)
        {
            newBaffles.append(baffles[i]);
        }
    }
    newBaffles.shrink();
    baffles.transfer(newBaffles);
}


Foam::autoPtr<mapPolyMesh> Foam::snappySnapDriver::surfaceToPatch
(
    const snapParameters& snapParams,
    const refinementParameters& refineParams
)
{

    Info<< nl
        << "Checking faces are correctly patched" << nl
        << "------------------------------------" << nl
        << endl;

    fvMesh& mesh = meshRefiner_.mesh();

    // Get the labels of added patches.
    labelList adaptPatchIDs
    (
       meshRefiner_.meshedPatches()
    );
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatchIDs
        )
    );
    indirectPrimitivePatch& pp = ppPtr();

    const scalarField snapDist
    (
        calcSnapDistance(meshRefiner_, snapParams, pp)
    );
    // Topology changes container
    polyTopoChange meshMod(mesh);

    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

    const labelList unzonedSurfaces =
        surfaceZonesInfo::getUnnamedSurfaces(surfZones);
    const labelList zonedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(surfZones);

    // Points that do not change.
    // Faces that do not move
    boolList boundaryZone(surfaces.surfaces().size(), false);
    bool foundBoundaryZone = false;
    PackedBoolList isZonedFace
    (
        setZonedFacesExclBoundaryNamed
        (
            meshRefiner_,
            refineParams
        )
    );

    forAll(surfZones, surfI)
    {
        const word& faceZoneName = surfZones[surfI].faceZoneName();

        if (faceZoneName.size())
        {
            const surfaceZonesInfo::faceZoneType& faceType =
                surfZones[surfI].faceType();

            if (faceType == surfaceZonesInfo::BOUNDARY)
            {
                boundaryZone[surfI] = true;
                foundBoundaryZone = true;
            }
        }
    }

     // 1. All faces to non-interface surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    pointField ppFaceCentres(pp.size());
    scalarField maxSnapDist(pp.size(), -GREAT);

    forAll(pp, i)
    {
        label faceI = pp.addressing()[i];
        const face& f = pp.localFaces()[i];

        ppFaceCentres[i] = mesh.faceCentres()[faceI];

        forAll(f, fp)
        {
            maxSnapDist[i] = max(snapDist[f[fp]], maxSnapDist[i]);
        }
    }

    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    surfaces.findNearest
    (
        unzonedSurfaces,
        ppFaceCentres,
        sqr(maxSnapDist), // sqr of attract distance
        hitSurface,
        hitInfo
    );

    labelList pRegions(hitInfo.size());
    vectorField sNormals(hitInfo.size());
    forAll(unzonedSurfaces, sI)
    {
       DynamicList<pointIndexHit> localHits;
       label surfI = unzonedSurfaces[sI];

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
          }
       }
    }

    //use surface normal to define side of boundary zone
    labelList twoSidedPatchID(pp.size(), -1);

    List<pointIndexHit> zoneHitInfo(pp.size(), pointIndexHit());
    vectorField sZoneNormals(pp.size(),vector::zero);

    if (foundBoundaryZone)
    {
        labelList zoneHitSurface;
        surfaces.findNearest
        (
            zonedSurfaces,
            ppFaceCentres,
            sqr(maxSnapDist), // sqr of attract distance
            zoneHitSurface,
            zoneHitInfo
         );

        labelList pZoneRegions(zoneHitInfo.size());
        forAll(zonedSurfaces, sI)
        {
            DynamicList<pointIndexHit> localHits;
            label surfI = zonedSurfaces[sI];

            forAll(zoneHitSurface, i)
            {
                if (zoneHitSurface[i] == surfI)
                {
                    localHits.append(zoneHitInfo[i]);
                }
            }
            labelList localRegions;
            vectorField localNormals;
            label geomI = surfaces.surfaces()[surfI];
            surfaces.geometry()[geomI].getRegion(localHits, localRegions);
            surfaces.geometry()[geomI].getNormal(localHits, localNormals);

            label localI = 0;
            forAll(zoneHitSurface, i)
            {
                if (zoneHitSurface[i] == surfI)
                {
                    pZoneRegions[i] = localRegions[localI];
                    sZoneNormals[i] = localNormals[localI];
                    localI++;
                }
            }
        }

        forAll(zoneHitInfo, i)
        {
            if (zoneHitInfo[i].hit())
            {
                label faceI = pp.addressing()[i];
                if (boundaryZone[zoneHitSurface[i]])
                {
                    vector fArea = mesh.faceAreas()[faceI];

                    label global =
                        surfaces.globalRegion
                        (zoneHitSurface[i], pZoneRegions[i]);

                    if ((fArea & sZoneNormals[i]) > 0.)
                    {
                        twoSidedPatchID[i] = globalToSlavePatch_[global];
                    }
                    else
                    {
                        twoSidedPatchID[i] = globalToMasterPatch_[global];
                    }
                }
            }
        }
    }

    // 2. Find if there is a fc-cc intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    pointField ppCellCentres(pp.size());

    forAll(pp, i)
    {
        label faceI = pp.addressing()[i];
        label cellI = mesh.faceOwner()[faceI];
        ppCellCentres[i] = mesh.cellCentres()[cellI];
    }

    labelList surf1;
    List<pointIndexHit> hPoint1;
    labelList pRegion1;
    labelList surf2;
    List<pointIndexHit> hPoint2;
    labelList pRegion2;

    surfaces.findNearestIntersection
    (
        unzonedSurfaces,
        ppCellCentres,
        ppFaceCentres,

        surf1,
        hPoint1,
        pRegion1,

        surf2,
        hPoint2,
        pRegion2
    );

    vectorField pNormals(hPoint1.size());
    forAll(unzonedSurfaces, sI)
    {
       DynamicList<pointIndexHit> localHits;
       label surfI = unzonedSurfaces[sI];

       forAll(surf1, i)
       {
          if (surf1[i] == surfI)
          {
             localHits.append(hPoint1[i]);
          }
       }
       pointField localNormals;
       label geomI = surfaces.surfaces()[surfI];
       surfaces.geometry()[geomI].getNormal(localHits, localNormals);

       label localI = 0;
       forAll(surf1, i)
       {
          if (surf1[i] == surfI)
          {
             pNormals[i] = localNormals[localI];
             localI++;
          }
       }
    }


    // Undistorted edge length
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    labelList closestPatch(pp.size(), -1);
    pointField closestPoint(pp.size(), vector::min);
    forAll(pp, i)
    {
        label faceI = pp.addressing()[i];
        if (!isZonedFace.get(faceI))
        {
            bool twoSided = false;
            if (twoSidedPatchID[i] != -1)
            {
                if (zoneHitInfo[i].hit())
                {
                    if (hitInfo[i].hit())
                    {
                        scalar hitZoneDist = mag
                        (
                            zoneHitInfo[i].hitPoint() - ppFaceCentres[i]
                        );
                        scalar hitDist = mag
                        (
                            hitInfo[i].hitPoint() - ppFaceCentres[i]
                        );

                        label own = mesh.faceOwner()[faceI];
                        scalar len = edge0Len / pow(2., cellLevel[own]);

                        vector fArea = mesh.faceAreas()[faceI];
                        fArea /= (mag(fArea)+SMALL);

                        if
                        (
                            (
                                hitZoneDist <= 1.01*hitDist
                             || hitZoneDist < 0.1*len
                            )
                            && mag(sZoneNormals[i]&fArea) > 0.707
                        )
                        {
                            twoSided = true;
                        }
                    }
                    else
                    {
                        twoSided = true;
                    }
                }
            }

            if (twoSided)
            {
                closestPatch[i] = twoSidedPatchID[i];
            }
            else if (hPoint1[i].hit() && hitInfo[i].hit())
            {
               vector surfNorm = pNormals[i];
               vector vec1 = hitInfo[i].hitPoint() - ppCellCentres[i];
               if ((surfNorm & vec1) < 0)
               {
                   surfNorm  = -surfNorm;
               }
               surfNorm /= (mag(surfNorm) + SMALL);

               vector faceArea = mesh.faceAreas()[faceI];
               faceArea /= (mag(faceArea) + SMALL);
               if ((surfNorm & faceArea) > 0.707)
               {
                   closestPatch[i] = globalToMasterPatch_
                   [
                       surfaces.globalRegion(surf1[i], pRegion1[i])
                   ];
                   closestPoint[i] = hPoint1[i].hitPoint();
               }
               else
               {
                   label global =
                       surfaces.globalRegion(hitSurface[i], pRegions[i]);
                   closestPatch[i] = globalToMasterPatch_[global];
                   closestPoint[i] = hitInfo[i].hitPoint();
               }
            }
            else if (hitInfo[i].hit())
            {
                label global =
                    surfaces.globalRegion(hitSurface[i], pRegions[i]);
                closestPatch[i] = globalToMasterPatch_[global];
                closestPoint[i] = hitInfo[i].hitPoint();
            }
        }
    }

    //Check for overlapping surfaces and set to minimum patch ID
    if (snapParams.repatchOverlapping())
    {
        scalar overlappingTol = scalar(1e-8);
        pointField overlapStart(pp.size(),vector::zero);
        pointField overlapEnd(pp.size(),vector::zero);
        labelList mapToPP(pp.size(),-1);
        label nChecks = 0;
        forAll(pp, i)
        {
            if (closestPoint[i] != vector::min)
            {
                label faceI = pp.addressing()[i];
                vector faceArea = mesh.faceAreas()[faceI];
                faceArea /= (mag(faceArea) + SMALL);
                vector projVec = overlappingTol*faceArea;
                overlapStart[nChecks] = closestPoint[i] - projVec;
                overlapEnd[nChecks] = closestPoint[i] + projVec;
                mapToPP[nChecks] = i;
                nChecks++;
            }
        }
        overlapStart.setSize(nChecks);
        overlapEnd.setSize(nChecks);
        mapToPP.setSize(nChecks);

        forAll(unzonedSurfaces, i)
        {
            label surfI = unzonedSurfaces[i];
            label geomI = surfaces.surfaces()[surfI];

            List<List<pointIndexHit>> hitInfo;
            surfaces.geometry()[geomI].findLineAll
            (
                overlapStart,
                overlapEnd,
                hitInfo
            );

            label n = 0;
            forAll(hitInfo, pointI)
            {
                n += hitInfo[pointI].size();
            }

            List<pointIndexHit> surfInfo(n);
            labelList pointMap(n);
            n = 0;

            forAll(hitInfo, pointI)
            {
                const List<pointIndexHit>& pHits = hitInfo[pointI];

                forAll(pHits, i)
                {
                    surfInfo[n] = pHits[i];
                    pointMap[n] = pointI;
                    n++;
                }
            }

            labelList surfRegion(n);
            surfaces.geometry()[geomI].getRegion(surfInfo,surfRegion);

            forAll(surfRegion, i)
            {
                label region = surfaces.globalRegion(surfI, surfRegion[i]);
                label pointI = mapToPP[pointMap[i]];
                if (region >= 0)
                {
                    closestPatch[pointI] = min
                    (
                        globalToMasterPatch_[region],
                        closestPatch[pointI]
                    );
                }
            }
        }
    }

    forAll(pp, i)
    {
        label closestPatchID = closestPatch[i];
        if (closestPatchID != -1)
        {
            label faceI = pp.addressing()[i];
            label patchI = mesh.boundaryMesh().whichPatch(faceI);

            if (patchI != closestPatchID)
            {
                label own = mesh.faceOwner()[faceI];
                label zoneID = mesh.faceZones().whichZone(faceI);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                }

                meshMod.modifyFace
                (
                    mesh.faces()[faceI], // modified face
                    faceI,          // label of face being modified
                    own,            // owner
                    -1,             // neighbour
                    false,          // face flip
                    closestPatchID,   // new patch for face
                    zoneID,         // zone for face
                    zoneFlip        // face flip in zone
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

    return map;
}


void Foam::snappySnapDriver::mergePatchFacesUndo
(
    const dictionary& motionDict,
    const List<labelList> &featureShadows,
    const bool maintainPatches
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    scalar minCos =
        Foam::cos(degToRad(20));

    scalar concaveCos =
        Foam::cos(degToRad(80));

    Info<< nl
        << "Merging all faces of a cell" << nl
        << "---------------------------" << nl
        << "    - which are on the same patch" << nl
        << "    - which make an angle < 45"
        << " degrees"
        << nl
        << "      (cos:" << minCos << ')' << nl
        << "    - as long as the resulting face doesn't become concave"
        << " by more than "
        << " 80 degrees" << nl
        << "      (0=straight, 180=fully concave)" << nl
        << endl;

    PackedList<1> isZonedFace(setZonedFaces(meshRefiner_));
    vectorField surfNormals(mesh.faceAreas());

    forAll(mesh.faces(), faceI)
    {
        if (isZonedFace.get(faceI))
        {
            surfNormals[faceI] = point::max;
        }
    }

    meshRefiner_.mergePatchFacesUndo
    (
         minCos,
         concaveCos,
         0.075, //max area ratio (-1 to disable)
         meshRefiner_.meshedPatches(),

         motionDict,
         false,
         maintainPatches,
         featureShadows,
         &surfNormals
    );

    meshRefiner_.mergeEdgesUndo(-1., motionDict, false);
}


void Foam::snappySnapDriver::mergePatchFacesUndo
(
    const dictionary& motionDict,
    const scalarField& snapDist,
    indirectPrimitivePatch &pp,
    List<labelPair>& baffles,
    bool calcNormals
)
{
//#ifdef FOAM_USE_TBB
//    Timer timer("snappySnapDriver::mergePatchFacesUndo");
//#endif

    fvMesh& mesh = meshRefiner_.mesh();
    PackedList<1> isZonedFace(setZonedFaces(meshRefiner_));

    boolList stationaryFaces(pp.size(), false);
    // Feature snapping currently not performed on zoned faces
    forAll(pp, i)
    {
        label faceI = pp.addressing()[i];
        if (isZonedFace.get(faceI))
        {
            stationaryFaces[i] = true;
        }
    }

    vectorField surfNormals(mesh.nFaces(), point::max);
    if (calcNormals)
    {
        pointField ppFaceCentres(pp.size(), vector::zero);
        pointField ppCellCentres(pp.size(), vector::zero);
        vectorField ppFaceAreas(pp.size(), vector::zero);
        forAll(pp, i)
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

        // Divide surfaces into zoned and unzoned
        const labelList
            allSurfaces(identity(meshRefiner_.surfaces().surfaces().size()));

        labelList pointIndex(pp.size(),label(-1));
        vectorField pNormals(pp.size(),vector::zero);
        List<pointIndexHit> finalHitInfo;

        findSurfaceNormals
        (
            allSurfaces,
            pp,
            snapDist,
            mesh.points(),
            ppFaceCentres,
            ppCellCentres,
            ppFaceAreas,
            stationaryFaces,
            false,
            finalHitInfo,
            pointIndex,
            pNormals
         );

        forAll(pointIndex, fI)
        {
            if (!finalHitInfo[fI].hit())
            {
                continue;
            }
            label i = pointIndex[fI];
            label faceI = pp.addressing()[i];

            point fA = mesh.faceAreas()[faceI];
            fA /= (mag(fA) + SMALL);

            //filter out problem surface normals
            if (mag(fA & pNormals[fI])<0.342)
            {
               continue;
            }

            surfNormals[faceI] = pNormals[fI];
        }
    }

    scalar minCos = Foam::cos(degToRad(45));

    scalar concaveCos = Foam::cos(degToRad(80));

    meshRefiner_.mergePatchFacesUndo
    (
        minCos,
        concaveCos,
        0.075, //max area ratio (-1 to disable)
        meshRefiner_.meshedPatches(),
        motionDict,
        true, //update intersections
        false,
        List<labelList>(),
        &surfNormals,
        &baffles
    );

    meshRefiner_.mergeEdgesUndo
    (
        -1.,
        motionDict,
        true, //update intersections
        &baffles
    );
}


void Foam::snappySnapDriver::splitCells
(
    const List<labelList>& shadowChains,
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const labelList& adaptPatchIDs,
    const PackedList<1>& isZonedFace,
    indirectPrimitivePatch &pp
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    DynamicList<Tuple2<label,labelPair>> cutFaces(mesh.nBoundaryFaces());

    if
    (
        snapParams.regionSnappedIDs().size() > 0
        && returnReduce(shadowChains.size(), sumOp<label>()) > 0
    )
    {
        Info<<"Splitting region interface boundary face cells"<<endl;
        forAll(shadowChains,scI)
        {
            const labelList &chain = shadowChains[scI];
            if (chain.size() < 2)
            {
                continue;
            }

            for (label pI=0; pI < chain.size()-1; ++pI)
            {
                //chain is stored as local points
                label pt1 = chain[pI];
                label pt2 = chain[pI+1];

                label mpt1 = pp.meshPoints()[pt1];
                label mpt2 = pp.meshPoints()[pt2];
                label meshedgei = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[mpt1],
                    mpt1,
                    mpt2
                );

                if (meshedgei == -1)
                {
                    labelList faces1 = pp.pointFaces()[pt1];
                    labelList faces2 = pp.pointFaces()[pt2];
                    std::list<label> intx;
                    std::set_intersection
                    (
                        faces1.begin(),faces1.end(),
                        faces2.begin(),faces2.end(),
                        std::back_inserter(intx)
                    );
                    label numCommon = intx.size();

                    if (numCommon == 1)
                    {
                        label meshfacei = pp.addressing()[intx.front()];
                        label patchi = patches.whichPatch(meshfacei);
                        if
                        (
                           patchi != -1 && !patches[patchi].coupled()
                           && !isZonedFace[meshfacei]
                        )
                        {
                            cutFaces.append
                            (
                                Tuple2<label,labelPair>
                                (
                                    meshfacei,
                                    labelPair(mpt1,mpt2)
                                )
                            );
                        }
                    }
                }
            }
        }
    }

    label nSplits = returnReduce(cutFaces.size(), sumOp<label>());

    if (nSplits > 0)
    {
        Info<<"Splitting " << nSplits <<" region boundary faces"<<endl;

        //Split boundary feature cells
        autoSplitCells splitRegionCells
        (
            meshRefiner_,
            labelList(0),
            (controller_.algorithm() == meshControl::EXTRUDE),
            snapParams.preMergeExtrude()
        );

        splitRegionCells.splitPreCalculated(cutFaces);
        resetPrimitivePatchAddressing(adaptPatchIDs, pp);
        //Repatch surface
        autoPtr<mapPolyMesh> map = surfaceToPatch
        (
            snapParams,
            refineParams
        );
        resetPrimitivePatchAddressing(adaptPatchIDs, pp);
    }
}


void Foam::snappySnapDriver::smoothSliverFaces
(
    const dictionary& motionDict,
    const snapParameters& snapParams,
    const label nInitErrors,
    const labelList& adaptPatchIDs,
    const List<labelPair>& baffles,
    indirectPrimitivePatch &pp
)
{
    Info<< nl
        << "Trying to smooth sliver faces" << nl
        << "-----------------------------" << nl
        << endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const pointMesh& pMesh = pointMesh::New(mesh);

    pointField newPoints = mesh.points();
    vectorField pointDisp(mesh.nPoints(), vector::zero);

    scalarField sumAreas(pp.nPoints(), 0.0);
    vectorField sumDisp(pp.nPoints(), vector::zero);

    boolList stationaryPoints(pp.meshPoints().size(), true);

    scalarField minArea(pp.meshPoints().size(), GREAT);
    scalarField maxArea(pp.meshPoints().size(), -GREAT);

    forAll(pp.meshPoints(), pointI)
    {
        const labelList& pFaces = pp.pointFaces()[pointI];
        forAll(pFaces, pFI)
        {
            label meshFaceI = pp.addressing()[pFaces[pFI]];
            scalar fArea = mesh.magFaceAreas()[meshFaceI];

            minArea[pointI] = min(minArea[pointI],fArea);
            maxArea[pointI] = max(maxArea[pointI],fArea);
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        minArea,
        minEqOp<scalar>(),
        GREAT            // null value
    );

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        maxArea,
        maxEqOp<scalar>(),
        -GREAT            // null value
    );

    forAll(pp.meshPoints(), pointI)
    {
        scalar ratio = minArea[pointI] / maxArea[pointI];
        if (ratio < 0.05)
        {
            stationaryPoints[pointI] = false;
        }
    }

    // Precalculate meshEdge per pp edge
    labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    labelList nExternalEdge(mesh.nEdges(), 0);
    forAll(pp.edges(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        const labelList& eFaces = pp.edgeFaces()[edgeI];
        nExternalEdge[meshEdgeI] = eFaces.size();
    }

    syncTools::syncEdgeList
    (
        mesh,
        nExternalEdge,
        plusEqOp<label>(),
        label(0)               // null value
    );

    forAll(pp.edges(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        if (nExternalEdge[meshEdgeI] == 1)
        {
            const edge& e = pp.edges()[edgeI];
            stationaryPoints[e[0]] = true;
            stationaryPoints[e[1]] = true;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        stationaryPoints,
        orEqOp<bool>(),
        false
    );

    forAll(pp.meshPoints(), pointI)
    {
        if (!stationaryPoints[pointI])
        {
            const labelList& pFaces = pp.pointFaces()[pointI];
            label meshPointI = pp.meshPoints()[pointI];

            bool foundCandidate = true;
            forAll(pFaces, pFI)
            {
                label meshFaceI = pp.addressing()[pFaces[pFI]];
                face f = mesh.faces()[meshFaceI];
                label start = 0;
                forAll(f,fp)
                {
                    if (f[start] == meshPointI)
                    {
                        label next = f.fcIndex(start);
                        label prev = f.rcIndex(start);

                        vector edgeVec1 = mesh.points()[f[next]] -
                            mesh.points()[f[start]];
                        vector edgeVec2 = mesh.points()[f[start]] -
                            mesh.points()[f[prev]];

                        edgeVec1 /= (mag(edgeVec1) + SMALL);
                        edgeVec2 /= (mag(edgeVec2) + SMALL);

                        if ((edgeVec1 & edgeVec2) < -0.944)
                        {
                            foundCandidate = false;
                            break;
                        }
                    }
                    start = f.fcIndex(start);
                }
            }
            if (!foundCandidate)
            {
                stationaryPoints[pointI] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        stationaryPoints,
        orEqOp<bool>(),
        false
    );

    // Extract faces per zone
    PackedList<1> isZonedFace(setZonedFaces(meshRefiner_));

    for (label i = 0; i < snapParams.nSliverSmooths(); i++)
    {
        pointField ppFaceCentres(pp.size());
        pointField ppFaceAreas(pp.size());

        forAll(pp, i)
        {
            const label& faceI = pp.addressing()[i];

            ppFaceCentres[i] =
                mesh.faces()[faceI].centre(newPoints);
            ppFaceAreas[i] =
                mesh.faces()[faceI].areaNormal(newPoints);
        }

        const scalar edge0Len =
            meshRefiner_.meshCutter().level0EdgeLength();
        const labelList cellLevel =
            meshRefiner_.meshCutter().cellLevel();
        sumDisp = vector::zero;
        sumAreas = 0.0;

        forAll(pp.meshPoints(), pointI)
        {
            if (!stationaryPoints[pointI])
            {
                const labelList& pFaces = pp.pointFaces()[pointI];
                forAll(pFaces, faceI)
                {
                    const label pFaceI = pFaces[faceI];
                    const label meshFaceI = pp.addressing()[pFaceI];

                    if (isZonedFace[meshFaceI])
                    {
                        continue;
                    }

                    const label own = mesh.faceOwner()[meshFaceI];
                    label  level = cellLevel[own];
                    scalar len = edge0Len / pow(2., level);
                    scalar weight = sqrt(mag(ppFaceAreas[pFaceI]))
                        / len;

                    sumDisp[pointI] +=  ppFaceCentres[pFaceI]
                        * weight;
                    sumAreas[pointI] += weight;
                }
            }
        }
        pointDisp = vector::zero;

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            sumDisp,
            plusEqOp<vector>(),
            vector::zero        // null value
         );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            sumAreas,
            plusEqOp<scalar>(),
            scalar(0.0)         // null value
         );

        forAll(pp.meshPoints(), i)
        {
            label meshPointI = pp.meshPoints()[i];

            if (!stationaryPoints[i] && sumAreas[i] > SMALL)
            {
                pointDisp[meshPointI] = sumDisp[i] / (sumAreas[i] + SMALL)
                    - newPoints[meshPointI];
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointDisp,
            maxMagSqrEqOp<point>(),
            vector::zero
         );

        forAll(pp.meshPoints(), i)
        {
            if (!stationaryPoints[i])
            {
                label meshPointI = pp.meshPoints()[i];
                newPoints[meshPointI] = newPoints[meshPointI]
                    + 0.1*pointDisp[meshPointI];
            }
        }
    }


    // The current mesh is the starting mesh to smooth from.
    motionSmoother meshMover
    (
        mesh,
        pp,
        adaptPatchIDs,
        meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs),
        motionDict
     );

    // Displacement at every patch point
    pointField patchDisp(meshMover.patch().nPoints(), vector::zero);

    pointDisp = vector::zero;

    forAll(pp.localPoints(), i)
    {
        pointDisp[pp.meshPoints()[i]] =
            newPoints[pp.meshPoints()[i]] - pp.localPoints()[i];
    }

    forAll(patchDisp, i)
    {
        patchDisp[i] = pointDisp[meshMover.patch().meshPoints()[i]];
    }

    // The current mesh is the starting mesh to smooth from.
    meshMover.setDisplacement(patchDisp);

    scaleMesh(snapParams, nInitErrors, baffles, meshMover);
    meshMover.correct();
}


void Foam::snappySnapDriver::smoothNonConformalFaces
(
    const dictionary& motionDict,
    const snapParameters& snapParams,
    const label nInitErrors,
    const labelList& adaptPatchIDs,
    const List<labelPair>& baffles,
    indirectPrimitivePatch &pp
)
{
    Info<< nl
        << "Trying to smooth non-conformal faces" << nl
        << "------------------------------------" << nl
        << endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const pointMesh& pMesh = pointMesh::New(mesh);
    const scalar edge0Len =
        meshRefiner_.meshCutter().level0EdgeLength();
    const labelList cellLevel =
        meshRefiner_.meshCutter().cellLevel();
    const labelList pointLevel =
        meshRefiner_.meshCutter().pointLevel();

    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    labelList snapSurfaces =
        surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
        (
            surfaces.surfZones()
        );

    const labelList& meshPoints = pp.meshPoints();

    pointField newPoints = mesh.points();
    vectorField pointDisp(mesh.nPoints(), vector::zero);

    scalarField sumAreas(pp.nPoints(), 0.0);
    vectorField sumDisp(pp.nPoints(), vector::zero);

    boolList stationaryPoints(pp.meshPoints().size(), true);

    //Find face centre nearest point normal
    vectorField fNormals(pp.size(), vector::zero);
    labelList faceHitSurface;
    {
        List<pointIndexHit> faceHitInfo;
        scalarField maxSnapDist(pp.size());
        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];
            label own = mesh.faceOwner()[meshFaceI];
            label  level = cellLevel[own];
            scalar len = edge0Len / pow(2., level);
            maxSnapDist[i] = 5*len*len;
        }

        surfaces.findNearest
        (
            snapSurfaces,
            pp.faceCentres(),
            maxSnapDist,
            faceHitSurface,
            faceHitInfo
        );

        forAll(snapSurfaces, sI)
        {
            DynamicList<pointIndexHit> localHits;
            label surfI = snapSurfaces[sI];

            forAll(faceHitSurface, i)
            {
                if (faceHitSurface[i] == surfI)
                {
                    localHits.append(faceHitInfo[i]);
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
            forAll(faceHitSurface, i)
            {
                if (faceHitSurface[i] == surfI)
                {
                    fNormals[i] = localNormals[localI];
                    localI++;
                }
            }
        }
    }

    label nFacePts = 0;
    forAll(pp, i)
    {
        nFacePts += pp.localFaces()[i].size();
    }

    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    scalarField sDist(nFacePts, 1);
    pointField sPts(nFacePts);
    vectorField pNormals(nFacePts, vector::zero);
    {
        nFacePts = 0;
        scalarField maxSnapDist(sPts.size());
        forAll(pp, i)
        {
            const face& f = pp.localFaces()[i];
            const point& fC = pp.faceCentres()[i];
            forAll(f,fp)
            {
                label localpointi = f[fp];
                label meshPointI = meshPoints[localpointi];
                label  level = pointLevel[meshPointI];
                scalar len = edge0Len / pow(2., level);
                maxSnapDist[nFacePts] = 5*len*len;
                sPts[nFacePts] = 0.5*(fC+newPoints[meshPointI]);
                nFacePts++;
            }
        }

        surfaces.findNearest
        (
            snapSurfaces,
            sPts,
            scalarField(sPts.size(), GREAT),
            hitSurface,
            hitInfo
        );

        nFacePts = 0;
        forAll(pp, i)
        {
            const face& f = pp.localFaces()[i];
            forAll(f,fp)
            {
                if (hitInfo[nFacePts].hit())
                {
                   scalar sd = mag
                   (
                      sPts[nFacePts]
                      -hitInfo[nFacePts].hitPoint()
                   );
                   label localpointi = f[fp];
                   label meshPointI = meshPoints[localpointi];
                   scalar len = edge0Len /
                      (1<<pointLevel[meshPointI]);
                   sDist[nFacePts] = min(scalar(1), (sd/len));
                }
                nFacePts++;
            }
        }

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
                    pNormals[i] = localNormals[localI];
                    localI++;
                }
            }
        }
    }

    scalar conformityTol = snapParams.conformityTol();
    nFacePts = 0;
    forAll(pp, i)
    {
        face f = pp.localFaces()[i];
        vector fArea = pp.faceAreas()[i]
            / (pp.magFaceAreas()[i] + SMALL);

        label nNonAligned = 0;
        label nAligned = 0;
        scalar fAve = 0;
        forAll(f,fp)
        {
            fAve += sDist[nFacePts];
            vector sNorm = pNormals[nFacePts];

            if (mag(fArea & sNorm) < 0.766)
            {
                nNonAligned++;
            }
            else
            {
                nAligned++;
            }
            nFacePts++;
        }

        if (mag(fArea & fNormals[i]) < 0.766)
        {
            nNonAligned++;
        }
        else
        {
            nAligned++;
        }

        scalar fConform = (fAve / f.size());

        bool movePts = false;
        if (faceHitSurface[i] == -1)
        {
            movePts = true;
        }
        else if (fConform > conformityTol)
        {
            movePts = true;
        }
        else if (nNonAligned > nAligned)
        {
            movePts = true;
        }

        if (movePts)
        {
            const face& f = pp.localFaces()[i];
            forAll(f,fp)
            {
                stationaryPoints[f[fp]] = false;
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        stationaryPoints,
        andEqOp<bool>(),
        false           // null value
    );

    // Precalculate meshEdge per pp edge
    labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));
    labelList nExternalEdge(mesh.nEdges(), 0);
    forAll(pp.edges(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        const labelList& eFaces = pp.edgeFaces()[edgeI];
        nExternalEdge[meshEdgeI] = eFaces.size();
    }

    syncTools::syncEdgeList
    (
        mesh,
        nExternalEdge,
        plusEqOp<label>(),
        label(0)               // null value
    );

    forAll(pp.edges(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        if (nExternalEdge[meshEdgeI] == 1)
        {
            const edge& e = pp.edges()[edgeI];
            stationaryPoints[e[0]] = true;
            stationaryPoints[e[1]] = true;
        }
    }

    syncTools::syncPointList
    (
        mesh,
        pp.meshPoints(),
        stationaryPoints,
        orEqOp<bool>(),
        false
    );

    // Extract faces per zone
    PackedList<1> isZonedFace(setZonedFaces(meshRefiner_));
    scalar weight = 0.2;

    for (label i = 0; i < snapParams.nSliverSmooths(); i++)
    {
        pointField ppFaceCentres(pp.size());
        pointField ppFaceAreas(pp.size());

        forAll(pp, i)
        {
            const label& faceI = pp.addressing()[i];

            ppFaceCentres[i] =
                mesh.faces()[faceI].centre(newPoints);
            ppFaceAreas[i] =
                mesh.faces()[faceI].areaNormal(newPoints);
        }

        sumDisp = vector::zero;
        sumAreas = 0.0;

        forAll(pp.meshPoints(), pointI)
        {
            if (!stationaryPoints[pointI])
            {
                const labelList& pFaces = pp.pointFaces()[pointI];
                forAll(pFaces, faceI)
                {
                    const label pFaceI = pFaces[faceI];
                    const label meshFaceI = pp.addressing()[pFaceI];

                    if (isZonedFace[meshFaceI])
                    {
                        continue;
                    }

                    const label own = mesh.faceOwner()[meshFaceI];
                    label  level = cellLevel[own];
                    scalar len = edge0Len / pow(2., level);
                    scalar weight = sqrt(mag(ppFaceAreas[pFaceI]))
                        / len;

                    sumDisp[pointI] +=  ppFaceCentres[pFaceI]
                        * weight;
                    sumAreas[pointI] += weight;
                }
            }
        }
        pointDisp = vector::zero;

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            sumDisp,
            plusEqOp<vector>(),
            vector::zero        // null value
         );

        syncTools::syncPointList
        (
            mesh,
            pp.meshPoints(),
            sumAreas,
            plusEqOp<scalar>(),
            scalar(0.0)         // null value
         );

        forAll(pp.meshPoints(), i)
        {
            label meshPointI = pp.meshPoints()[i];

            if (!stationaryPoints[i] && sumAreas[i] > SMALL)
            {
                pointDisp[meshPointI] = sumDisp[i] / (sumAreas[i] + SMALL)
                    - newPoints[meshPointI];
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointDisp,
            maxMagSqrEqOp<point>(),
            vector::zero
         );

        forAll(pp.meshPoints(), i)
        {
            if (!stationaryPoints[i])
            {
                label meshPointI = pp.meshPoints()[i];
                newPoints[meshPointI] = newPoints[meshPointI]
                    + weight*pointDisp[meshPointI];
            }
        }
    }

    // The current mesh is the starting mesh to smooth from.
    motionSmoother meshMover
    (
        mesh,
        pp,
        adaptPatchIDs,
        meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs),
        motionDict
     );

    // Displacement at every patch point
    pointField patchDisp(meshMover.patch().nPoints(), vector::zero);

    pointDisp = vector::zero;

    forAll(pp.localPoints(), i)
    {
        pointDisp[pp.meshPoints()[i]] =
            newPoints[pp.meshPoints()[i]] - pp.localPoints()[i];
    }

    forAll(patchDisp, i)
    {
        patchDisp[i] = pointDisp[meshMover.patch().meshPoints()[i]];
    }

    // The current mesh is the starting mesh to smooth from.
    meshMover.setDisplacement(patchDisp);

    scaleMesh(snapParams, nInitErrors, baffles, meshMover);
    meshMover.correct();
}


void Foam::snappySnapDriver::collapseSmallEdges
(
    const dictionary& motionDict,
    const labelList& adaptPatchIDs,
    const snapParameters& snapParams,
    indirectPrimitivePatch &pp
)
{
//#ifdef FOAM_USE_TBB
//    Timer t("snappySnapDriver::collapseSmallEdges");
//#endif

    fvMesh& mesh = meshRefiner_.mesh();

    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    label nCollapseIter = 0;
    boolList notMerged(pp.localPoints().size(),false);

    //Try and collapse small edges
    while (nCollapseIter < 3)
    {
        Info<<"Edge collapsing Iter: "<< nCollapseIter<< endl;

        labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

        // Normal component of normals of connected faces.
        vectorField edgeNormal(mesh.nEdges(), vector(GREAT, GREAT, GREAT));

        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            forAll(eFaces, i)
            {
                nomalsCombine()
                (
                    edgeNormal[meshEdgeI],
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

        boolList featureEdge(mesh.nEdges(), false);
        boolList featurePoint(mesh.nPoints(), false);
        labelList nPPEdgeFaces(mesh.nEdges(), 0);

        forAll(pp.edgeFaces(), edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];
            label meshEdgeI = meshEdges[edgeI];

            nPPEdgeFaces[meshEdgeI] += eFaces.size();

            const vector n = edgeNormal[meshEdgeI];

            if (n != vector(GREAT, GREAT, GREAT))
            {
                scalar cos = n & pp.faceNormals()[eFaces[0]];
                if (cos < 0.939)
                {
                    featureEdge[meshEdgeI] = true;
                    const edge e = mesh.edges()[meshEdgeI];
                    featurePoint[e[0]] = true;
                    featurePoint[e[1]] = true;
                }
            }
        }


        syncTools::syncEdgeList
        (
            mesh,
            nPPEdgeFaces,
            plusEqOp<label>(),
            label(0)
        );

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
            orEqOp<bool>(),
            false
        );

        syncTools::syncEdgeList
        (
            mesh,
            nPPEdgeFaces,
            plusEqOp<label>(),
            label(0)
        );

        labelList edgeLevel(mesh.nEdges(), -1);
        const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
        forAll(cellLevel, cellI)
        {
            const label cLevel = cellLevel[cellI];
            const labelList& cEdges = mesh.cellEdges(cellI);

            forAll(cEdges, i)
            {
                label edgeI = cEdges[i];

                if (edgeLevel[edgeI] == -1)
                {
                    edgeLevel[edgeI] = cLevel;
                }
                else if (edgeLevel[edgeI] < cLevel)
                {
                    edgeLevel[edgeI] = cLevel;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            edgeLevel,
            maxEqOp<label>(),
            label(-1)
        );
        label nCollapses = 0;
        boolList collapseEdges(mesh.nEdges(), false);

        const labelList& fOwn = mesh.faceOwner();
        const labelList& fNei = mesh.faceNeighbour();
        const labelListList& fEdges = mesh.faceEdges();
        const labelListList& eFaces = mesh.edgeFaces();
        const labelListList& pFaces = mesh.pointFaces();

        forAll(pp.edges(), edgeI)
        {
            const edge& e = pp.edges()[edgeI];
            scalar edgeLength = e.mag(pp.localPoints());

            label v0 = pp.meshPoints()[e[0]];
            label v1 = pp.meshPoints()[e[1]];

            label meshEdgeI =  meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[v0],
                v0,
                v1
            );

            scalar edgeLengthStart = edge0Len /
                (1<<edgeLevel[meshEdgeI]);

            bool previousCollapse = false;
            if (notMerged[e[0]] && notMerged[e[1]])
            {
                previousCollapse = true;
            }

            //Check for non-orthogonality and try merging these edges too
            point eC = e.centre(pp.localPoints());
            const labelList& eFaces = pp.edgeFaces()[edgeI];
            vector eV = e.vec(pp.localPoints())/(e.mag(pp.localPoints())+SMALL);

            bool nonOrtho = false;
            forAll(eFaces, eFI)
            {
                label meshFaceI = pp.addressing()[eFaces[eFI]];
                point fC = mesh.faceCentres()[meshFaceI];
                vector fcToEc = fC-eC;
                fcToEc /= (mag(fcToEc) + SMALL);

                if (mag(fcToEc & eV) > 0.766)
                {
                    nonOrtho = true;
                    break;
                }
            }

            if
            (
                (
                    (nonOrtho && edgeLength < 0.75*edgeLengthStart)
                    || (edgeLength < snapParams.collapseTol()*edgeLengthStart)
                )
                && !previousCollapse
                && nPPEdgeFaces[meshEdgeI] < 3
            )
            {
                collapseEdges[meshEdgeI] = true;
                nCollapses++;
            }
        }

        //collapse non-manifold edges
        forAll(pp, i)
        {
            const label meshFaceI = pp.addressing()[i];
            const label own = fOwn[meshFaceI];
            face f = pp.localFaces()[i];
            labelHashSet cFaces(mesh.cells()[own]);

            label foundEdge = -1;
            scalar minEdgeLength = GREAT;

            forAll(f, fp)
            {
                label meshPointI = pp.meshPoints()[f[fp]];

                label prev = f.rcIndex(fp);
                label next = f.fcIndex(fp);

                label ppEdge1 =  meshTools::findEdge
                (
                    pp.edges(),
                    pp.pointEdges()[f[next]],
                    f[next],
                    f[fp]
                );

                label meshEdge1 =  meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[pp.meshPoints()[f[next]]],
                    pp.meshPoints()[f[next]],
                    pp.meshPoints()[f[fp]]
                );


                vector evec1 = pp.localPoints()[f[next]]
                    -pp.localPoints()[f[fp]];
                scalar eLength1 = mag(evec1);
                evec1 /= (eLength1 + SMALL);

                label ppEdge2 =  meshTools::findEdge
                (
                    pp.edges(),
                    pp.pointEdges()[f[prev]],
                    f[prev],
                    f[fp]
                );

                label meshEdge2 =  meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[pp.meshPoints()[f[prev]]],
                    pp.meshPoints()[f[prev]],
                    pp.meshPoints()[f[fp]]
                );

                vector evec2 = pp.localPoints()[f[fp]]
                    -pp.localPoints()[f[prev]];
                scalar eLength2 = mag(evec2);
                evec2 /= (eLength2 + SMALL);

                label npf = 0;
                forAll(pFaces[meshPointI], pFI)
                {
                    if (cFaces.found(pFaces[meshPointI][pFI]))
                    {
                        npf++;
                    }
                }

                if (npf == 2)
                {
                    //look for edge to collapse. weight feature edge higher
                    if ((evec1&evec2)<0.647)
                    {
                        if (featureEdge[meshEdge1])
                        {
                            if (0.01*eLength1 <  minEdgeLength)
                            {
                                minEdgeLength = 0.01*eLength1;
                                foundEdge = ppEdge1;
                            }
                        }
                        else if (eLength1 <  minEdgeLength)
                        {
                            minEdgeLength = eLength1;
                            foundEdge = ppEdge1;
                        }

                        if (featureEdge[meshEdge2])
                        {
                            if (0.01*eLength2 <  minEdgeLength)
                            {
                                minEdgeLength = 0.01*eLength2;
                                foundEdge = ppEdge2;
                            }
                        }
                        else if (eLength2 <  minEdgeLength)
                        {
                            minEdgeLength = eLength2;
                            foundEdge = ppEdge2;
                        }
                    }
                }
            }

            if (foundEdge != -1)
            {
                const edge e = pp.edges()[foundEdge];

                if (!(notMerged[e[0]] && notMerged[e[1]]))
                {
                    label meshEdgeI =  meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[pp.meshPoints()[e[0]]],
                        pp.meshPoints()[e[0]],
                        pp.meshPoints()[e[1]]
                    );

                    if (nPPEdgeFaces[meshEdgeI] < 3)
                    {
                        collapseEdges[meshEdgeI] = true;
                        nCollapses++;
                    }
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            collapseEdges,
            orEqOp<bool>(),
            false
        );

        //Check for degenerate faceFaces formed from collapsing
        labelList cellOwnerMarker(fOwn.size(), -1);
        labelList cellNeighMarker(fOwn.size(), -1);

//        label nDegenerateFaceCells = 0;

        labelList checkFaces(meshRefiner_.selectInnerFaces());

        forAll(checkFaces, fI)
        {
            label faceI = checkFaces[fI];

            label nAddOwnFaces = 0;

            label fownCellI = fOwn[faceI];

            forAll(fEdges[faceI], feI)
            {
                label cgEdge(fEdges[faceI][feI]);

                if (collapseEdges[cgEdge])
                {
                    continue;
                }

                forAll(eFaces[cgEdge], efI)
                {
                    label cgFace(eFaces[cgEdge][efI]);

                    if (cgFace != faceI && cellOwnerMarker[cgFace] != faceI)
                    {
                        if (fOwn[cgFace] == fownCellI)
                        {
                            nAddOwnFaces++;
                            cellOwnerMarker[cgFace] = faceI;
                        }
                        else if (cgFace < mesh.nInternalFaces())
                        {
                            if (fNei[cgFace] == fownCellI)
                            {
                                nAddOwnFaces++;
                                cellOwnerMarker[cgFace] = faceI;
                            }
                        }
                    }
                }
            }

            if (nAddOwnFaces < 3)
            {
                forAll(fEdges[faceI], feI)
                {
                    label cgEdge(fEdges[faceI][feI]);

                    if (collapseEdges[cgEdge])
                    {
                        collapseEdges[cgEdge] = false;
                        nCollapses--;
                    }
                }
//                nDegenerateFaceCells++;
            }


            if (faceI < mesh.nInternalFaces())
            {
                label nAddNeiFaces = 0;

                label fneiCellI = fNei[faceI];

                forAll(fEdges[faceI], feI)
                {
                    label cgEdge(fEdges[faceI][feI]);

                    if (collapseEdges[cgEdge])
                    {
                        continue;
                    }

                    forAll(eFaces[cgEdge], efI)
                    {
                        label cgFace(eFaces[cgEdge][efI]);

                        if (cgFace != faceI && cellNeighMarker[cgFace] != faceI)
                        {
                            if (fOwn[cgFace] == fneiCellI)
                            {
                                nAddNeiFaces++;
                                cellNeighMarker[cgFace] = faceI;
                            }
                            else if (cgFace < mesh.nInternalFaces())
                            {
                                if (fNei[cgFace] == fneiCellI)
                                {
                                    nAddNeiFaces++;
                                    cellNeighMarker[cgFace] = faceI;
                                }
                            }
                        }
                    }
                }

                if (nAddNeiFaces < 3)
                {
                    forAll(fEdges[faceI], feI)
                    {
                        label cgEdge(fEdges[faceI][feI]);

                        if (collapseEdges[cgEdge])
                        {
                            collapseEdges[cgEdge] = false;
                            nCollapses--;
                        }
                    }
//                    nDegenerateFaceCells++;
                }
            }
        }

        if (!returnReduce(nCollapses, sumOp<label>()))
        {
            break;
        }

        syncTools::syncEdgeList
        (
            mesh,
            collapseEdges,
            andEqOp<bool>(),
            false
        );

/*
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        forAll(collapseEdges, edgeI)
        {
            if (collapseEdges[edgeI])
            {
                const labelList& eFaces = mesh.edgeFaces()[edgeI];
                forAll(eFaces, efI)
                {
                    label eFace = eFaces[efI];
                    label patchI = patches.whichPatch(eFace);

                    if
                    (
                        patchI != -1
                        && isA<processorPolyPatch>(patches[patchI])
                    )
                    {
                        if (mesh.faces()[eFace].size() == 3)
                        {
                            collapseEdges[edgeI] = false;
                        }
                    }
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            collapseEdges,
            andEqOp<bool>(),
            false
        );
*/
        boolList collapsedPoint(mesh.nPoints(), false);
        boolList marked(mesh.nEdges(), false);

        //sort edges so edges collapsed consistently on multi-processor cases
        label nCollapsed = 0;
        forAll(mesh.edges(), edgeI)
        {
            if (collapseEdges[edgeI])
            {
                nCollapsed++;
            }
        }
        List<label> cEdges(nCollapsed);
        SortableList<scalar> maxEdgeLengthList(nCollapsed);

        nCollapsed = 0;
        forAll(mesh.edges(), edgeI)
        {
            if (collapseEdges[edgeI])
            {
                cEdges[nCollapsed] = edgeI;
                maxEdgeLengthList[nCollapsed] =
                    mesh.edges()[edgeI].mag(mesh.points());

                nCollapsed++;
            }
        }
        maxEdgeLengthList.reverseSort();

        const labelList& sortedEdgeIndices = maxEdgeLengthList.indices();

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (Pstream::myProcNo() == procI)
            {
                forAll(sortedEdgeIndices, indexI)
                {
                    label edgeI = cEdges[sortedEdgeIndices[indexI]];
                    const edge& e = mesh.edges()[edgeI];
                    if
                    (
                        collapseEdges[edgeI]
                        && !collapsedPoint[e[0]] && !collapsedPoint[e[1]]
                     )
                    {
                        marked[edgeI] = true;
                        collapsedPoint[e[0]] = true;
                        collapsedPoint[e[1]] = true;
                    }
                }
            }

            syncTools::syncEdgeList
            (
                mesh,
                marked,
                orEqOp<bool>(),
                false           // null value
            );

            syncTools::syncPointList
            (
                mesh,
                collapsedPoint,
                orEqOp<bool>(),
                false           // null value
            );
        }

        const pointMesh& pMeshCollapse = pointMesh::New(mesh);

        dictionary updatedMotionDict(motionDict);
        //Don't perform edge length check
        updatedMotionDict.add("minEdgeLength",-1,true);

        motionSmoother collapseMover
        (
            mesh,
            pp,
            adaptPatchIDs,
            meshRefinement::makeDisplacementField(pMeshCollapse, adaptPatchIDs),
            updatedMotionDict
        );

        pointVectorField& displacement = collapseMover.displacement();

        boolList midPointCollapse(mesh.nPoints(), true);
        forAll(mesh.edges(), edgeI)
        {
            if (collapseEdges[edgeI] && marked[edgeI])
            {
                const edge& e = mesh.edges()[edgeI];
                vector edgeCentre = e.centre(mesh.points());

                if (featureEdge[edgeI])
                {
                    displacement[e[0]] = edgeCentre - mesh.points()[e[0]];
                    displacement[e[1]] = edgeCentre - mesh.points()[e[1]];
                }
                else
                {
                    if (featurePoint[e[0]] && featurePoint[e[1]])
                    {
                        collapseEdges[edgeI] = false;
                        nCollapses--;
                        continue;
                    }
                    else if (featurePoint[e[0]] || featurePoint[e[1]])
                    {
                        midPointCollapse[e[0]] = false;
                        midPointCollapse[e[1]] = false;
                        displacement[e[0]] =
                        (
                            featurePoint[e[0]] ?
                            vector::zero : mesh.points()[e[1]] - mesh.points()[e[0]]
                        );

                        displacement[e[1]] =
                        (
                            featurePoint[e[0]] ?
                            mesh.points()[e[0]] - mesh.points()[e[1]] : vector::zero
                        );
                    }
                    else
                    {
                        displacement[e[0]] = edgeCentre - mesh.points()[e[0]];
                        displacement[e[1]] = edgeCentre - mesh.points()[e[1]];
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            midPointCollapse,
            andEqOp<bool>(),
            false
        );

        syncTools::syncEdgeList
        (
            mesh,
            collapseEdges,
            andEqOp<bool>(),
            false
        );

        syncTools::syncPointList
        (
            mesh,
            displacement,
            maxMagEqOp(),
            vector::zero
        );

        labelHashSet wrongFaces(mesh.nFaces()/100);
        //check whole mesh
        label nInitErrors = motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);

        checkFaces = labelList(meshRefiner_.selectInnerFaces());
        scalar oldErrorReduction = collapseMover.setErrorReduction(0.0);
        label oldSmoothScale = collapseMover.setSmoothScale(0);

        collapseMover.setErrorReduction(0.0);

        for (label snapIter = 0; snapIter < 2*snapParams.nSnap(); snapIter++)
        {
            Info<< nl << "Scaling iteration " << snapIter << endl;

            if
            (
                collapseMover.scaleMesh
                (
                    checkFaces,
                    List<labelPair>(0),
                    true,
                    nInitErrors
                 )
            )
            {
                Info<< "Successfully moved mesh" << endl;
                break;
            }
            else
            {
                forAll(mesh.edges(), edgeI)
                {
                    if (collapseEdges[edgeI] && marked[edgeI])
                    {
                        const edge& e = mesh.edges()[edgeI];
                        scalar mag0 = mag(displacement[e[0]]);
                        scalar mag1 = mag(displacement[e[1]]);
                        if (mag0 > SMALL && (mag1 < SMALL && midPointCollapse[e[1]]))
                        {
                            displacement[e[0]] = vector::zero;
                        }
                        if (mag1 > SMALL && (mag0 < SMALL && midPointCollapse[e[0]]))
                        {
                            displacement[e[1]] = vector::zero;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            displacement,
            maxMagEqOp(),
            vector::zero
        );

        collapseMover.setErrorReduction(oldErrorReduction);
        collapseMover.setSmoothScale(oldSmoothScale);
        collapseMover.correct();

        DynamicList<point> possibleMerge(pp.localPoints().size());
        DynamicList<label> meshPointMap(pp.localPoints().size());
        forAll(mesh.points(), pointI)
        {
            if (mag(displacement[pointI]) > SMALL || !midPointCollapse[pointI])
            {
                possibleMerge.append(mesh.points()[pointI]);
                meshPointMap.append(pointI);
            }
        }
        possibleMerge.shrink();
        meshPointMap.shrink();

        labelList pointMap;
        pointField newPoints;
        bool hasMerged = mergePoints
        (
            possibleMerge,                  // points
            SMALL,                          // mergeTol
            false,                          // verbose
            pointMap,
            newPoints
        );

        if (!returnReduce(hasMerged, orOp<bool>()))
        {
            break;
        }

        labelList firstOldPoint(newPoints.size(), -1);
        polyTopoChange collapseMod(mesh);
        Map<label> pointToMaster(pointMap.size() - newPoints.size());

        const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

        forAll(pointMap, i)
        {
            label oldPointI = meshPointMap[i];
            label newPointI = pointMap[i];

            if (firstOldPoint[newPointI] == -1)
            {
                // First use of oldPointI. Store.
                firstOldPoint[newPointI] = oldPointI;
            }
            else
            {
                label meshPointCurrent = oldPointI;
                label meshPointPrev = firstOldPoint[newPointI];
                label pointToRemove = -1;
                label pointToKeep = -1;

                if
                (
                    pointLevel[meshPointPrev]
                    > pointLevel[meshPointCurrent]
                )
                {
                    pointToRemove = meshPointCurrent;
                    pointToKeep = meshPointPrev;
                }
                else
                {
                    pointToRemove = meshPointPrev;
                    pointToKeep = meshPointCurrent;
                }

                collapseMod.setAction
                (
                    polyRemovePoint(pointToRemove)
                );

                pointToMaster.insert
                (
                    meshPointCurrent,
                    pointToKeep
                );

                pointToMaster.insert
                (
                    meshPointPrev,
                    pointToKeep
                );
            }
        }

        const faceList& faces = mesh.faces();

        checkFaces = labelList(meshRefiner_.selectInnerFaces());

        forAll(checkFaces, fI)
        {
            const label faceI = checkFaces[fI];

            const face& f = faces[faceI];

            bool hasMerged = false;

            forAll(f, fp)
            {
                label pointI = f[fp];
                Map<label>::const_iterator iter = pointToMaster.find(pointI);

                if (iter != pointToMaster.end())
                {
                    hasMerged = true;
                }
            }

            if (hasMerged)
            {
                face newF(f);

                forAll(f, fp)
                {
                    label pointI = f[fp];
                    Map<label>::const_iterator iter =
                        pointToMaster.find(pointI);

                    if (iter != pointToMaster.end())
                    {
                        newF[fp] = iter();
                    }
                }
                newF.collapse();

                if (newF.size() > 2)
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

                    collapseMod.setAction
                    (
                        polyModifyFace
                        (
                            newF,                       // modified face
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
                else
                {
                    collapseMod.setAction
                    (
                        polyRemoveFace(faceI)
                     );
                }
            }
        }
        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> map = collapseMod.changeMesh(mesh, false, true);

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


        resetPrimitivePatchAddressing(adaptPatchIDs,pp);

        meshRefiner_.updateMesh(map,labelList(0));

        notMerged.resize(pp.localPoints().size());
        notMerged = false;
        // Map data
        const labelList& origPointMap = map().pointMap();

        forAll(pp.meshPoints(), pointI)
        {
            label meshPointI = pp.meshPoints()[pointI];
            label origPointI = origPointMap[meshPointI];

            if (origPointI != -1)
            {
                scalar magDisp = mag(displacement[meshPointI]);
                if
                (
                    collapsedPoint[origPointI]
                    && (magDisp < SMALL && midPointCollapse[origPointI])
                )
                {
                    notMerged[pointI] = true;
                }
            }
        }
        nCollapseIter++;
    }

    meshRefiner_.meshCutter().syncLevel();
}


void Foam::snappySnapDriver::correctConcaveFaces
(
    const dictionary& motionDict,
    const labelList& adaptPatchIDs,
    const snapParameters& snapParams,
    indirectPrimitivePatch &pp
)
{
//#ifdef FOAM_USE_TBB
//    Timer t("snappySnapDriver::correctConcaveFaces");
//#endif

    Info<<"Correcting concave faces "<<endl;

    fvMesh& mesh = meshRefiner_.mesh();

    vectorField pointDisp(mesh.nPoints(), vector::zero);
    const vectorField& faceAreas = mesh.faceAreas();
    const pointField& p = mesh.points();

    scalar concaveTol = snapParams.concaveTol();

    forAll(pp, i)
    {
        label faceI = pp.addressing()[i];

        vector nf = faceAreas[faceI];
        scalar magNf = mag(nf);
        nf /= magNf + VSMALL;

        const face& f = mesh.faces()[faceI];

        forAll(f, fp)
        {
            vector prevN = triPointRef
            (
                p[f[fp]],
                p[f.nextLabel(fp)],
                p[f.prevLabel(fp)]
             ).areaNormal();

            scalar triArea = mag(prevN);

            if (mag(triArea) > SMALL)
            {
                prevN /= triArea;

                if
                (
                    (prevN&nf) < 0
                    && (triArea/mag(faceAreas[faceI])) > concaveTol
                )
                {
                    point midPt = 0.5*
                        (p[f.nextLabel(fp)]+ p[f.prevLabel(fp)]);
                    pointDisp[f[fp]] = midPt - p[f[fp]];
                }
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

    const pointMesh& pMesh= pointMesh::New(mesh);

    // The current mesh is the starting mesh to smooth from.
    motionSmoother meshMover
    (
        mesh,
        pp,
        adaptPatchIDs,
        meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs),
        motionDict
     );

    // Displacement at every patch point
    pointField patchDisp(meshMover.patch().nPoints(), vector::zero);

    forAll(patchDisp, i)
    {
        patchDisp[i] = pointDisp[meshMover.patch().meshPoints()[i]];
    }

    // The current mesh is the starting mesh to smooth from.
    meshMover.setDisplacement(patchDisp);

    labelHashSet wrongFaces(mesh.nFaces()/100);
    label nInitErrors = motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);

    scaleMesh
    (
        snapParams,
        nInitErrors,
        List<labelPair>(0),
        meshMover
    );

    meshMover.correct();

    return;
}


void Foam::snappySnapDriver::addTetrahedralCellsToSplitMesh
(
    const labelList& adaptPatchIDs
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    polyTopoChange meshMod(mesh);
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatchIDs
         )
     );
    indirectPrimitivePatch& pp = ppPtr();

    PackedBoolList isMasterPPEdge(snappyLayerDriver::getMasterPPEdges(mesh,pp));

    // Precalculate meshEdge per pp edge
    labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));

    labelList nPPEdgeFaces(mesh.nEdges(), 0);
    forAll(pp.edgeFaces(), patchEdgeI)
    {
        label meshEdgeI = meshEdges[patchEdgeI];
        if (isMasterPPEdge[meshEdgeI])
        {
            nPPEdgeFaces[meshEdgeI] += pp.edgeFaces()[patchEdgeI].size();
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        nPPEdgeFaces,
        plusEqOp<label>(),
        label(0)              // null value
     );

    pointField pointNormals(PatchTools::pointNormals(mesh, pp));

    const labelList
        allSurfaces(identity(meshRefiner_.surfaces().surfaces().size()));

    List<pointIndexHit> hitInfo;
    labelList hitSurf;
    meshRefiner_.surfaces().findNearest
    (
        allSurfaces,
        pp.localPoints(),
        scalarField(pp.localPoints().size(), GREAT),
        hitSurf,
        hitInfo
     );

    vectorField fNormals(pp.localPoints().size(), vector(GREAT, GREAT, GREAT));

    forAll(allSurfaces, sI)
    {
        label surfI = allSurfaces[sI];
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
        label geomI = meshRefiner_.surfaces().surfaces()[surfI];
        meshRefiner_.surfaces().geometry()[geomI].getNormal
        (
            localHits,
            localNormals
        );

        label localI = 0;
        forAll(hitSurf, pointI)
        {
            if (hitSurf[pointI] == surfI)
            {
                fNormals[pointI] = localNormals[localI];
                localI++;
            }
        }
    }

    boolList changedFaces(mesh.nFaces(), false);

    label nChanged = 0;
    forAll(pp.edges(), edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        if (nPPEdgeFaces[meshEdgeI] == 2 && pp.edgeFaces()[edgeI].size() == 2)
        {
            edge ppEdge = pp.edges()[edgeI];

            label face0 = pp.addressing()[pp.edgeFaces()[edgeI][0]];
            label face1 = pp.addressing()[pp.edgeFaces()[edgeI][1]];

            face fpp0 = pp.localFaces()[pp.edgeFaces()[edgeI][0]];
            face fpp1 = pp.localFaces()[pp.edgeFaces()[edgeI][1]];

            if
            (
                changedFaces[face0] || changedFaces[face1]
                || fNormals[ppEdge[0]] == vector(GREAT, GREAT, GREAT)
                || fNormals[ppEdge[1]] == vector(GREAT, GREAT, GREAT)
            )
            {
                continue;
            }

            face f0 = mesh.faces()[face0];
            face f1 = mesh.faces()[face1];

            if (f0.size() == 3 && f1.size() == 3)
            {
                vector norm0 = mesh.faceAreas()[face0];
                norm0 /= (mag(norm0) + SMALL);

                vector norm1 = mesh.faceAreas()[face1];
                norm1 /= (mag(norm1) + SMALL);

                vector edgeNorm = (norm0 + norm1) / 2.;

                vector sn0 = fNormals[ppEdge[0]];
                vector sn1 = fNormals[ppEdge[1]];
                sn0 /= (mag(sn0) + SMALL);
                sn1 /= (mag(sn1) + SMALL);

                if ((sn0 & pointNormals[ppEdge[0]]) < 0.)
                {
                    sn0 = -sn0;
                }

                if ((sn1 & pointNormals[ppEdge[1]]) < 0.)
                {
                    sn1 = -sn1;
                }

                vector eVec = pp.localPoints()[ppEdge[1]]
                    - pp.localPoints()[ppEdge[0]];

                scalar vSum = (sn1 & eVec) - (sn0 & eVec);

                point eC = mesh.edges()[meshEdgeI].centre(mesh.points());
                vector ecfc0 = mesh.faceCentres()[face0] - eC;

                scalar pSum = ecfc0 & edgeNorm;

                if
                (
                    mag(norm0 & norm1) < 0.1 && pSum > 0.
                    && mag(sn0 & sn1) < 0.707 && vSum > 0.
                )
                {
                    changedFaces[face0] = true;
                    changedFaces[face1] = true;

                    label meshPoint0 = pp.meshPoints()[ppEdge[0]];
                    label meshPoint1 = pp.meshPoints()[ppEdge[1]];

                    label otherPoint0 = -1;
                    label otherPointPP0 = -1;
                    forAll(f0, fp)
                    {
                        if (f0[fp] != meshPoint0 && f0[fp] != meshPoint1)
                        {
                            otherPoint0 = f0[fp];
                            otherPointPP0 = fpp0[fp];
                        }
                    }

                    label otherPoint1 = -1;
                    label otherPointPP1 = -1;
                    forAll(f1, fp)
                    {
                        if (f1[fp] != meshPoint0 && f1[fp] != meshPoint1)
                        {
                            otherPoint1 = f1[fp];
                            otherPointPP1 = fpp1[fp];
                        }
                    }

                    point a0 = mesh.points()[meshPoint0];
                    point a1 = mesh.points()[meshPoint1];
                    point a2 = mesh.points()[otherPoint0];
                    point a3 = mesh.points()[otherPoint1];

                    tetrahedron<point, point> newTet(a0,a1,a2,a3);
                    scalar tetVol = newTet.mag();
                    point tetCentre = newTet.centre();

                    vector eCtC = eC - tetCentre;

                    if
                    (
                        hitInfo[ppEdge[0]].hit()
                        && hitInfo[ppEdge[1]].hit()
                        && hitInfo[otherPointPP0].hit()
                        && hitInfo[otherPointPP1].hit()
                    )
                    {
                        a0 = hitInfo[ppEdge[0]].hitPoint();
                        a1 = hitInfo[ppEdge[1]].hitPoint();
                        a2 = hitInfo[otherPointPP0].hitPoint();
                        a3 = hitInfo[otherPointPP1].hitPoint();

                        tetrahedron<point, point>  deformedTet(a0,a1,a2,a3);
                        scalar deformedVol = deformedTet.mag();

                        if (deformedVol/tetVol < 0.3)
                        {
                            continue;
                        }

                        //check whether tet centre sits on wrong side of surface
                        {
                            pointField projPts(3);
                            projPts[0] = a2;
                            projPts[1] = a3;
                            projPts[2] = a0;

                            triFace projFace(identity(3));
                            scalar projArea = projFace.mag(projPts);
                            if (projArea > SMALL)
                            {
                                pointHit hit =
                                    projFace.nearestPoint(tetCentre, projPts);
                                if (hit.hit())
                                {
                                    vector hitVec = hit.hitPoint() - tetCentre;
                                    if ((hitVec & eCtC) > 0)
                                    {
                                        continue;
                                    }
                                }
                            }
                        }
                        {
                            pointField projPts(3);
                            projPts[0] = a2;
                            projPts[1] = a3;
                            projPts[2] = a1;

                            triFace projFace(identity(3));
                            scalar projArea = projFace.mag(projPts);
                            if (projArea > SMALL)
                            {
                                pointHit hit =
                                    projFace.nearestPoint(tetCentre, projPts);
                                if (hit.hit())
                                {
                                    vector hitVec = hit.hitPoint() - tetCentre;
                                    if ((hitVec & eCtC) > 0)
                                    {
                                        continue;
                                    }
                                }
                            }
                        }
                    }

                    nChanged++;

                    face newFace0(3);
                    {
                        label start = findIndex(f0, otherPoint0);
                        newFace0[0] = otherPoint0;
                        newFace0[1] = f0[f0.fcIndex(start)];
                        newFace0[2] = otherPoint1;
                    }

                    face newFace1(3);
                    {
                        label start = findIndex(f1, otherPoint1);
                        newFace1[0] = otherPoint1;
                        newFace1[1] = f1[f1.fcIndex(start)];
                        newFace1[2] = otherPoint0;
                    }

                    label ownZoneI = mesh.cellZones().whichZone
                    (
                        mesh.faceOwner()[face0]
                     );

                    label newCellI = meshMod.setAction
                    (
                        polyAddCell
                        (
                            -1,             // master point
                            -1,             // master edge
                            -1,             // master face
                            mesh.faceOwner()[face0],//master
                            ownZoneI        // zone for cell
                         )
                     );

                    //convert face0 to internal face
                    {
                        label zoneID = mesh.faceZones().whichZone(face0);
                        bool zoneFlip = false;
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = mesh.faceZones()[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                        }

                        meshMod.setAction
                        (
                            polyModifyFace
                            (
                                f0,                       // modified face
                                face0,                    // label of face
                                mesh.faceOwner()[face0], // owner
                                newCellI,                 // neighbour
                                false,                    // face flip
                                -1,                       // patch for face
                                false,                    // remove from zone
                                zoneID,                   // zone for face
                                zoneFlip                  // face flip in zone
                             )
                         );
                    }

                    //convert face1 to internal face
                    {
                        label zoneID = mesh.faceZones().whichZone(face1);
                        bool zoneFlip = false;
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = mesh.faceZones()[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                        }

                        meshMod.setAction
                        (
                            polyModifyFace
                            (
                                f1,                       // modified face
                                face1,                    // label of face
                                mesh.faceOwner()[face1], // owner
                                newCellI,                 // neighbour
                                false,                    // face flip
                                -1,                       // patch for face
                                false,                    // remove from zone
                                zoneID,                   // zone for face
                                zoneFlip                  // face flip in zone
                             )
                         );
                    }

                    //Add new face0
                    {
                        label zoneID = mesh.faceZones().whichZone(face0);
                        bool zoneFlip = false;
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = mesh.faceZones()[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                        }
                        label patchI = mesh.boundaryMesh().whichPatch(face0);

                        meshMod.setAction
                        (
                            polyAddFace
                            (
                                newFace0, // vertices
                                newCellI,      // owner,
                                -1,             // neighbour,
                                -1,             // masterPointID,
                                -1,             // masterEdgeID,
                                face0,    // masterFaceID,
                                false,          // flipFaceFlux,
                                patchI,         // patchID,
                                zoneID,         // zoneID,
                                zoneFlip        // zoneFlip
                             )
                         );
                    }

                    //Add new face1
                    {
                        label zoneID = mesh.faceZones().whichZone(face1);
                        bool zoneFlip = false;
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = mesh.faceZones()[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                        }
                        label patchI = mesh.boundaryMesh().whichPatch(face1);

                        meshMod.setAction
                        (
                            polyAddFace
                            (
                                newFace1, // vertices
                                newCellI,      // owner,
                                -1,             // neighbour,
                                -1,             // masterPointID,
                                -1,             // masterEdgeID,
                                face1,    // masterFaceID,
                                false,          // flipFaceFlux,
                                patchI,         // patchID,
                                zoneID,         // zoneID,
                                zoneFlip        // zoneFlip
                             )
                         );
                    }
                }
            }
        }
    }

    nChanged = returnReduce(nChanged, sumOp<label>());

    if (nChanged != 0)
    {
        Info<<"Constructing additional : "<<nChanged<<" tetrahedral cells"<<endl;

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


Foam::labelList Foam::snappySnapDriver::zonePatchIDs() const
{
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const labelList& surfaceGeometry = surfaces.surfaces();

    const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();
    DynamicList<label> zonePatches(globalToMasterPatch_.size());

    forAll(surfZones, surfI)
    {
        const word& faceZoneName = surfZones[surfI].faceZoneName();

        if (faceZoneName.size())
        {
            const surfaceZonesInfo::faceZoneType& faceType =
                surfZones[surfI].faceType();

            if (faceType != surfaceZonesInfo::BOUNDARY)
            {
                label geomI = surfaceGeometry[surfI];
                const wordList& regNames =
                    surfaces.geometry().regionNames()[geomI];

                forAll(regNames, i)
                {
                    label globalRegionI = surfaces.globalRegion(surfI, i);
                    zonePatches.append(globalToMasterPatch_[globalRegionI]);
                    zonePatches.append(globalToSlavePatch_[globalRegionI]);
                }
            }
        }
    }
    return zonePatches.shrink();
}


Foam::labelList Foam::snappySnapDriver::zonePatchSlaveIDs() const
{
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const labelList& surfaceGeometry = surfaces.surfaces();

    const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();
    DynamicList<label> zonePatches(globalToMasterPatch_.size());

    forAll(surfZones, surfI)
    {
        const word& faceZoneName = surfZones[surfI].faceZoneName();

        if (faceZoneName.size())
        {
            const surfaceZonesInfo::faceZoneType& faceType =
                surfZones[surfI].faceType();

            if (faceType != surfaceZonesInfo::BOUNDARY)
            {
                label geomI = surfaceGeometry[surfI];
                const wordList& regNames =
                    surfaces.geometry().regionNames()[geomI];

                forAll(regNames, i)
                {
                    label globalRegionI = surfaces.globalRegion(surfI, i);
                    zonePatches.append(globalToSlavePatch_[globalRegionI]);
                }
            }
        }
    }
    return zonePatches.shrink();
}


void Foam::snappySnapDriver::detectWarpedFaces
(
    const scalar featureCos,
    const indirectPrimitivePatch& pp,

    DynamicList<label>& splitFaces,
    DynamicList<labelPair>& splits
) const
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const faceList& localFaces = pp.localFaces();
    const pointField& localPoints = pp.localPoints();
    const labelList& bFaces = pp.addressing();

    splitFaces.clear();
    splitFaces.setCapacity(bFaces.size());
    splits.clear();
    splits.setCapacity(bFaces.size());

    // Determine parallel consistent normals on points
    const vectorField pointNormals(PatchTools::pointNormals(mesh, pp));

    face f0(4);
    face f1(4);

    forAll(localFaces, facei)
    {
        const face& f = localFaces[facei];

        if (f.size() >= 4)
        {
            // See if splitting face across diagonal would make two faces with
            // biggish normal angle

            labelPair minDiag(-1, -1);
            scalar minCos(GREAT);

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
                    // Form two faces
                    f0.setSize(endFp-startFp+1);
                    label i0 = 0;
                    for (label fp = startFp; fp <= endFp; fp++)
                    {
                        f0[i0++] = f[fp];
                    }
                    f1.setSize(f.size()+2-f0.size());
                    label i1 = 0;
                    for (label fp = endFp; fp != startFp; fp = f.fcIndex(fp))
                    {
                        f1[i1++] = f[fp];
                    }
                    f1[i1++] = f[startFp];

                    //Info<< "Splitting face:" << f << " into f0:" << f0
                    //    << " f1:" << f1 << endl;

                    vector n0 = f0.areaNormal(localPoints);
                    scalar n0Mag = mag(n0);
                    vector n1 = f1.areaNormal(localPoints);
                    scalar n1Mag = mag(n1);

                    if (n0Mag > ROOTVSMALL && n1Mag > ROOTVSMALL)
                    {
                        scalar cosAngle = (n0/n0Mag) & (n1/n1Mag);
                        if (cosAngle < minCos)
                        {
                            minCos = cosAngle;
                            minDiag = labelPair(startFp, endFp);
                        }
                    }
                }
            }


            if (minCos < featureCos)
            {
                splitFaces.append(bFaces[facei]);
                splits.append(minDiag);
            }
        }
    }
}


Foam::labelList Foam::snappySnapDriver::getInternalOrBaffleDuplicateFace() const
{
    const fvMesh& mesh = meshRefiner_.mesh();

    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = meshRefiner_.getZones(fzTypes);
    }

    List<labelPair> baffles
    (
        meshRefiner_.subsetBaffles
        (
            mesh,
            internalOrBaffleFaceZones,
            localPointRegion::findDuplicateFacePairs(mesh)
        )
    );

    labelList faceToDuplicate(mesh.nFaces(), -1);
    forAll(baffles, i)
    {
        const labelPair& p = baffles[i];
        faceToDuplicate[p[0]] = p[1];
        faceToDuplicate[p[1]] = p[0];
    }

    return faceToDuplicate;
}


void Foam::snappySnapDriver::removeLayersExtrude
(
    const dictionary& motionDict,
    const snapParameters& snapParams,
    const label nInitErrors,
    const dictionary& layerDict,
    const labelList& adaptPatchIDs,
    indirectPrimitivePatch &pp
)
{
    Info<< "Removing layers at grown up patches "<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    layerParameters layerParams(layerDict, mesh.boundaryMesh(), true);
    const labelList& numLayers = layerParams.numLayers();

    DynamicList<label> noLayerPatchIDs(patches.size());
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && numLayers[patchI] == -1)
        {
            noLayerPatchIDs.append(patchI);
        }
    }
    noLayerPatchIDs.shrink();

    labelList pointFtrType(mesh.nPoints(), -1);
    pointField pointFtrDir(mesh.nPoints(), vector::zero);
    pointField pointFtrOrigin(mesh.nPoints(), vector::zero);
    boolList visited(mesh.nPoints(), false);

    //remove layers already generated on grown up layer patches
    {
        autoPtr<indirectPrimitivePatch> ppCellRemovalPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                noLayerPatchIDs
            )
        );

        const indirectPrimitivePatch& ppCellRemoval = ppCellRemovalPtr();

        if (returnReduce(ppCellRemoval.size(), sumOp<label>()) == 0)
        {
            //Make sure we try an optmize mesh at least once
            dictionary meshOptimDict
            (
                meshDict_.found("meshOptimization") ?
                meshDict_.subDict("meshOptimization") :
                dictionary()
            );

            if (meshOptimDict.found("foamOptimizeCoeffs"))
            {
                dictionary& coeffsDict =
                    meshOptimDict.subDict("foamOptimizeCoeffs");
                if (!coeffsDict.found("meshQualityControls"))
                {
                    coeffsDict.add("meshQualityControls",motionDict,true);
                }
            }

            autoPtr<autoOptimize> optimMeshPtr
                = autoOptimize::New(mesh, meshOptimDict);
            optimMeshPtr->optimize();

            return;
        }

        labelList boundaryEdges(mesh.nEdges(), -1);
        //TODO loop only through innerGrid faces //GGG
        forAll(mesh.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            if (patchi != -1 && !patches[patchi].coupled())
            {
                const labelList& fEdges = mesh.faceEdges()[facei];
                forAll(fEdges, fe)
                {
                    boundaryEdges[fEdges[fe]] = patchi;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            boundaryEdges,
            maxEqOp<label>(),
            label(-1)               // null value
        );

        const labelList& meshPoints = ppCellRemoval.meshPoints();
        const labelList meshEdges
        (
            ppCellRemoval.meshEdges(mesh.edges(), mesh.pointEdges())
        );

        boolList excludedFaces(ppCellRemoval.size(), false);
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            ppCellRemoval,
            meshEdges,
            excludedFaces,
            0.93969,
            0.93969
         );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        pointField pointNormals =
            eClass.calculatePointNormals(excludedFaces, 0, false);

        PackedBoolList isMasterPPEdge
        (
            snappyLayerDriver::getMasterPPEdges(mesh,ppCellRemoval)
        );

        labelList nFtrEdges(meshPoints.size(), 0);
        List<pointField> pECC(meshPoints.size(), pointField(0));


        List<surfaceZonesInfo::faceZoneType> fzType(2);
        fzType[0] = surfaceZonesInfo::INTERNAL;
        fzType[1] = surfaceZonesInfo::BAFFLE;
        const labelList zonesToBaffle
        (
            meshRefiner_.getZones(fzType)
        );

        boolList zoneEdges(mesh.nEdges(), false);
        forAll(zonesToBaffle, j)
        {
            label zonei = zonesToBaffle[j];
            const faceZone& fz = mesh.faceZones()[zonei];
            forAll(fz, i)
            {
                label facei = fz[i];
                const labelList& fEdges = mesh.faceEdges()[facei];
                forAll(fEdges, fei)
                {
                    zoneEdges[fEdges[fei]] = true;
                }
            }
        }
        syncTools::syncEdgeList
        (
            mesh,
            zoneEdges,
            orEqOp<bool>(),
            false               // null value
        );

        //Set boundary points
        forAll(ppCellRemoval.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            const edge e = mesh.edges()[meshEdgeI];

            if (eType[edgeI].first() == edgeClassification::BOUNDARY)
            {
                pointFtrType[e[0]] = label(0);
                pointFtrType[e[1]] = label(0);
            }
            else if (zoneEdges[meshEdgeI])
            {
                pointFtrType[e[0]] = label(4);
                pointFtrType[e[1]] = label(4);
            }
        }
        syncTools::syncPointList
        (
            mesh,
            pointFtrType,
            maxEqOp<label>(),
            label(-1)
        );

        forAll(ppCellRemoval.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            const edge e = mesh.edges()[meshEdgeI];
            const edge le = ppCellRemoval.edges()[edgeI];

            label pType0 = -1;
            label pType1 = -1;

            if
            (
                eType[edgeI].first() == edgeClassification::CONCAVE
                || eType[edgeI].first() == edgeClassification::CONVEX
            )
            {
                pType0 = label(2);
                pType1 = label(2);
                if (isMasterPPEdge[meshEdgeI])
                {
                    nFtrEdges[le[0]]++;
                    nFtrEdges[le[1]]++;
                }
            }
            else
            {
                pType0 = label(1);
                pType1 = label(1);
            }

            if (pointFtrType[e[0]] != 0)
            {
                pointFtrType[e[0]] = max(pointFtrType[e[0]],pType0);
            }

            if (pointFtrType[e[1]] != 0)
            {
                pointFtrType[e[1]] = max(pointFtrType[e[1]],pType1);
            }
        }

        forAll(meshPoints, pti)
        {
            const labelList& pEdges = ppCellRemoval.pointEdges()[pti];
            pointField checkPts(pEdges.size());
            label nPts = 0;
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                label meshEdgeI = meshEdges[edgei];

                if
                (
                    isMasterPPEdge[meshEdgeI] &&
                    (
                        eType[edgei].first() == edgeClassification::CONCAVE
                        || eType[edgei].first() == edgeClassification::CONVEX
                    )
                )
                {
                    const edge& e = mesh.edges()[meshEdgeI];
                    point eC = e.centre(mesh.points());
                    checkPts[nPts++] = eC;
                }
            }
            checkPts.setSize(nPts);
            pECC[pti] = checkPts;
        }

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            pECC,
            meshRefinement::pointFieldCombine(),
            pointField() // initial value
        );

        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nFtrEdges,
            plusEqOp<label>(),
            label(0) // null value
        );

        forAll(meshPoints, pti)
        {
            if (nFtrEdges[pti] > 2)
            {
                label pointi = meshPoints[pti];
                pointFtrType[pointi] = label(3);
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointFtrType,
            maxEqOp<label>(),
            label(-1)
        );

        forAll(meshPoints, pti)
        {
            label pointi = meshPoints[pti];
            visited[pointi] = true;
            pointFtrOrigin[pointi] = mesh.points()[pointi];
            const pointField& eCC = pECC[pti];
            if (pointFtrType[pointi] == 2 && eCC.size() == 2)
            {
                pointFtrDir[pointi] = eCC[1]-eCC[0];
            }
            else if (pointFtrType[pointi] == 0 && eCC.size() == 1)
            {
                pointFtrDir[pointi] = pointFtrOrigin[pointi] - eCC[0];
            }
            else
            {
                pointFtrDir[pointi] = pointNormals[pti];
            }
        }

        forAll(meshPoints, pti)
        {
            label pointi = meshPoints[pti];
            const labelList& pEdges = mesh.pointEdges()[pointi];
            forAll(pEdges, pei)
            {
                label edgei = pEdges[pei];
                if (boundaryEdges[edgei] == -1)
                {
                    const edge& e = mesh.edges()[edgei];
                    label otherPt(pointi == e[0] ? e[1] : e[0]);
                    label pType = pointFtrType[pointi];
                    if (pType == 0)
                    {
                        if (pECC[pti].size() > 0)
                        {
                            pType = 2;
                        }
                        else
                        {
                            pType = 1;
                        }
                    }

                    point pDir = pointFtrDir[pointi];
                    point pOrig = pointFtrOrigin[pointi];
                    pointFtrType[otherPt] = pType;
                    pointFtrDir[otherPt] = pDir;
                    pointFtrOrigin[otherPt] = pOrig;
                    visited[otherPt] = true;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            visited,
            orEqOp<bool>(),
            false
        );

        syncTools::syncPointList
        (
            mesh,
            pointFtrDir,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            pointFtrOrigin,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh,
            pointFtrType,
            maxEqOp<label>(),
            label(-1)
        );

        DynamicList<label> cellsToRemove(ppCellRemoval.size());
        labelList exposedPatchID(mesh.nFaces(), -1);

        forAll(ppCellRemoval, i)
        {
            label facei = ppCellRemoval.addressing()[i];
            label patchI = patches.whichPatch(facei);
            label own = mesh.faceOwner()[facei];
            cellsToRemove.append(own);
            cell c = mesh.cells()[own];
            forAll(c, cFI)
            {
                exposedPatchID[c[cFI]] = patchI;
            }
        }
        cellsToRemove.shrink();

        syncTools::syncFaceList
        (
            mesh,
            exposedPatchID,
            maxEqOp<label>()
        );

        removeCells cellRemover(mesh);
        labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
        labelList exposedPatches(exposedFaces.size());

        forAll(exposedFaces, i)
        {
            label faceI = exposedFaces[i];
            exposedPatches[i] = exposedPatchID[faceI];
        }

        polyTopoChange meshMod(mesh);

        cellRemover.setRefinement
        (
            cellsToRemove,
            exposedFaces,
            exposedPatches,
            meshMod
        );

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

        meshRefinement::updateList
        (
            map().pointMap(),
            label(-1),
            pointFtrType
        );
        meshRefinement::updateList
        (
            map().pointMap(),
            vector::zero,
            pointFtrDir
        );
        meshRefinement::updateList
        (
            map().pointMap(),
            vector::zero,
            pointFtrOrigin
        );
        meshRefinement::updateList
        (
            map().pointMap(),
            false,
            visited
        );
        resetPrimitivePatchAddressing(adaptPatchIDs,pp);
    }

    //Project and smooth
    {
        autoPtr<indirectPrimitivePatch> ppCellRemovalPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                noLayerPatchIDs
            )
        );
        const indirectPrimitivePatch& ppCellRemoval = ppCellRemovalPtr();

        PackedBoolList isMasterPPEdge
        (
            snappyLayerDriver::getMasterPPEdges(mesh,ppCellRemoval)
        );

        const labelList& meshPoints = ppCellRemoval.meshPoints();
        const labelList meshEdges
        (
            ppCellRemoval.meshEdges(mesh.edges(), mesh.pointEdges())
        );

        labelList nReset(meshPoints.size(), 0);
        pointField resetOrigin(meshPoints.size(), vector::zero);
        vectorField resetDir(meshPoints.size(), vector::zero);
        labelList resetType(meshPoints.size(), -1);

        boolList boundaryPts(mesh.nPoints(), false);

        forAll(meshPoints, pti)
        {
            label pointi = meshPoints[pti];
            boundaryPts[pointi] = true;
            if (!visited[pointi])
            {
                const labelList& pEdges = ppCellRemoval.pointEdges()[pti];
                forAll(pEdges, pEI)
                {
                    label edgei = meshEdges[pEdges[pEI]];
                    edge e = mesh.edges()[edgei];
                    label otherPt(e[0] == pointi ? e[1] : e[0]);
                    if (visited[otherPt])
                    {
                        nReset[pti]++;
                        resetOrigin[pti] += pointFtrOrigin[otherPt];
                        resetDir[pti] += pointFtrDir[otherPt];
                        resetType[pti] = max
                        (
                            pointFtrType[otherPt],
                            resetType[pti]
                        );
                    }
                }
            }
        }
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            nReset,
            plusEqOp<label>(),
            label(0) // null value
        );
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            resetOrigin,
            plusEqOp<point>(),
            vector::zero // null value
        );
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            resetDir,
            plusEqOp<point>(),
            vector::zero // null value
        );
        syncTools::syncPointList
        (
            mesh,
            meshPoints,
            resetType,
            maxEqOp<label>(),
            label(-1) // null value
        );
        syncTools::syncPointList
        (
            mesh,
            boundaryPts,
            orEqOp<bool>(),
            false // null value
        );

        bool foundAllBdy = false;
        forAll(ppCellRemoval, i)
        {
            label facei = ppCellRemoval.addressing()[i];
            label own = mesh.faceOwner()[facei];
            const labelList& cellPts = mesh.cellPoints()[own];
            label nCellBdy = 0;
            forAll(cellPts, ptI)
            {
                if (boundaryPts[cellPts[ptI]])
                {
                    nCellBdy++;
                }
            }
            if (nCellBdy == cellPts.size())
            {
                foundAllBdy = true;
                break;
            }
        }

        if (returnReduce(foundAllBdy, orOp<bool>()))
        {
            WarningInFunction
               << nl
               << "Cells have been detected where all points are boundary."
               << " This can generate squashed cells at grown up interface."
               << nl
               << " Try re-running and setting castellatedMeshControls keyword"
               << " removeBoundaryAnyGrownUp to true"
               << nl << endl;
        }

        forAll(meshPoints, pti)
        {
            label pointi = meshPoints[pti];
            if (nReset[pti] > 0)
            {
                pointFtrOrigin[pointi] = resetOrigin[pti]/nReset[pti];
                pointFtrDir[pointi] = resetDir[pti]/nReset[pti];
                pointFtrType[pointi] = resetType[pti];
            }
        }

        pointField newPoints =  mesh.points();
        vectorField pointDisp(mesh.nPoints(), vector::zero);
        labelList nSet(mesh.nPoints(), 0);

        forAll(meshPoints, pti)
        {
            label pointi = meshPoints[pti];

            point pOrig = pointFtrOrigin[pointi];
            point pDir = pointFtrDir[pointi];

            if (pointFtrType[pointi] == 1)
            {
                plane fPlane(pOrig, pDir);
                vector disp = fPlane.nearestPoint(newPoints[pointi])
                    - newPoints[pointi];
                pointDisp[pointi] += disp;
                nSet[pointi]++;
            }
            else if (pointFtrType[pointi] == 2)
            {
                point startPt = pOrig + 100*pDir;
                point endPt = pOrig - 100*pDir;
                pointHit lHit = linePointRef(startPt, endPt).nearestDist
                    (newPoints[pointi]);
                if (lHit.hit())
                {
                    vector disp = lHit.hitPoint()
                    - newPoints[pointi];
                    pointDisp[pointi] += disp;
                    nSet[pointi]++;
                }
            }
            else if (pointFtrType[pointi] == 3 || pointFtrType[pointi] == 4)
            {
                vector disp = pointFtrOrigin[pointi]
                    - newPoints[pointi];
                pointDisp[pointi] += disp;
                nSet[pointi]++;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pointDisp,
            plusEqOp<point>(),
            vector::zero         // null value
        );

        syncTools::syncPointList
        (
            mesh,
            nSet,
            plusEqOp<label>(),
            label(0)         // null value
        );

        forAll(pointDisp, pointi)
        {
            if (nSet[pointi] > 0)
            {
                newPoints[pointi] += (pointDisp[pointi] / nSet[pointi]);
            }
        }

        pointField ppFaceCentres(ppCellRemoval.size(), vector::zero);
        vectorField ppFaceAreas(ppCellRemoval.size(), vector::zero);
        scalarField sumAreas(ppCellRemoval.nPoints(), 0.0);
        vectorField sumDisp(ppCellRemoval.nPoints(), vector::zero);

        const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
        const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

        for (label nSmooth=0; nSmooth < 4; ++nSmooth)
        {
            sumAreas = 0;
            sumDisp =  vector::zero;
            forAll(ppCellRemoval, i)
            {
                const label& faceI = ppCellRemoval.addressing()[i];
                ppFaceCentres[i] =
                    mesh.faces()[faceI].centre(newPoints);
                ppFaceAreas[i] =
                    mesh.faces()[faceI].areaNormal(newPoints);
            }

            forAll(meshPoints, pti)
            {
                label pointi = meshPoints[pti];

                if (pointFtrType[pointi] == 1)
                {
                    const labelList pFaces = ppCellRemoval.pointFaces()[pti];
                    forAll(pFaces, pfi)
                    {
                        label pFaceI = pFaces[pfi];
                        const label facei =
                            ppCellRemoval.addressing()[pFaceI];
                        const label own = mesh.faceOwner()[facei];
                        label  level = cellLevel[own];

                        scalar len = edge0Len / pow(2., level);
                        scalar weight =
                            sqrt(mag(ppFaceAreas[pFaceI]))/ len;

                        sumDisp[pti] += ppFaceCentres[pFaceI]
                            * weight;
                        sumAreas[pti] += weight;
                    }
                }
                else if (pointFtrType[pointi] == 2)
                {
                    const labelList& pEdges =
                        ppCellRemoval.pointEdges()[pti];
                    forAll(pEdges, eI)
                    {
                        label edgei = pEdges[eI];
                        label meshEdgeI = meshEdges[edgei];
                        const edge& e = mesh.edges()[meshEdgeI];
                        label otherPt(pointi == e[0] ? e[1] : e[0]);
                         if
                        (
                            (
                                pointFtrType[otherPt] == 0
                                || pointFtrType[otherPt] >= 2
                            )
                            && isMasterPPEdge[meshEdgeI]
                        )
                        {
                            point eC = e.centre(newPoints);
                            scalar weight = e.mag(newPoints);
                            sumDisp[pti] +=  eC * weight;
                            sumAreas[pti] += weight;
                        }
                    }
                }
            }
            syncTools::syncPointList
            (
                mesh,
                meshPoints,
                sumDisp,
                plusEqOp<vector>(),
                vector::zero        // null value
             );

            syncTools::syncPointList
            (
                mesh,
                meshPoints,
                sumAreas,
                plusEqOp<scalar>(),
                scalar(0.0)         // null value
             );

            pointDisp = vector::zero;
            forAll(meshPoints, i)
            {
                label meshPointI = meshPoints[i];

                if (sumAreas[i] > SMALL)
                {
                    pointDisp[meshPointI] = (sumDisp[i] / sumAreas[i])
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

            newPoints += pointDisp;
        }

        mesh.movePoints(newPoints);

        //Final re-snap
        {
            const pointMesh& pMeshSnap = pointMesh::New(mesh);

            autoPtr<indirectPrimitivePatch> ppCellRemovalPtr
            (
                meshRefinement::makePatch
                (
                    mesh,
                    noLayerPatchIDs
                 )
            );
            indirectPrimitivePatch& ppCellRemoval = ppCellRemovalPtr();

            motionSmoother meshMoverSnap
            (
                mesh,
                ppCellRemoval,
                noLayerPatchIDs,
                meshRefinement::makeDisplacementField
                (
                    pMeshSnap,
                    noLayerPatchIDs
                ),
                motionDict
            );

            labelList snapSurfaces =
                surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
            (
                surfaces.surfZones()
            );

            label nFeatureIter = 1;
            optStages ostages;

            direction ot = optStages::RELAXEDERROR;
            ostages.add(nFeatureIter-1,ot);

            // Snap faces to lock on to feature edges
            featureSnap
            (
                nFeatureIter,
                1,
                snapSurfaces,
                ppCellRemoval,
                snapParams,
                motionDict,
                meshMoverSnap,
                false, //initial snap
                ostages
            );

            scaleMesh
            (
                snapParams,
                nInitErrors,
                List<labelPair>(0),
                meshMoverSnap
            );
            meshMoverSnap.correct();
        }
    }
    resetPrimitivePatchAddressing(adaptPatchIDs,pp);
}


void Foam::snappySnapDriver::removeLayersDual
(
    const dictionary& motionDict,
    const snapParameters& snapParams,
    const label nInitErrors,
    const dictionary& layerDict,
    const labelList& adaptPatchIDs,
    const labelList& featureIDs,
    indirectPrimitivePatch &pp
)
{
    Info<< "Removing layers at grown up patches "<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    layerParameters layerParams(layerDict, mesh.boundaryMesh(), true);
    const labelList& numLayers = layerParams.numLayers();

    DynamicList<label> noLayerPatchIDs(patches.size());
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && numLayers[patchI] == -1)
        {
            noLayerPatchIDs.append(patchI);
        }
    }
    noLayerPatchIDs.shrink();

    label nDualLayers = 1;
    if (snapParams.nAdditionalDualLayers() > 0)
    {
        nDualLayers += snapParams.nAdditionalDualLayers();
    }

    labelList origPatchNewMesh(0);
    pointField origNormalNewMesh(0);
    pointField origCentreNewMesh(0);

    vector greatPoint = vector(GREAT, GREAT, GREAT);

    labelList originatingPatch(mesh.nFaces(), -1);
    pointField origCentre(mesh.nFaces(), greatPoint);
    pointField origNormal(mesh.nFaces(), greatPoint);

    //remove layers already generated on zero layer patches
    {
        autoPtr<indirectPrimitivePatch> ppCellRemovalPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                noLayerPatchIDs
            )
        );

        const indirectPrimitivePatch& ppCellRemoval = ppCellRemovalPtr();

        if (returnReduce(ppCellRemoval.size(), sumOp<label>()) == 0)
        {
            //Make sure we try an optmize mesh at least once
            dictionary meshOptimDict
            (
                meshDict_.found("meshOptimization") ?
                meshDict_.subDict("meshOptimization") :
                dictionary()
            );

            if (meshOptimDict.found("foamOptimizeCoeffs"))
            {
                dictionary& coeffsDict =
                    meshOptimDict.subDict("foamOptimizeCoeffs");
                if (!coeffsDict.found("meshQualityControls"))
                {
                    coeffsDict.add("meshQualityControls",motionDict,true);
                }
            }

            autoPtr<autoOptimize> optimMeshPtr
                = autoOptimize::New(mesh, meshOptimDict);
            optimMeshPtr->optimize();

            return;
        }

        //Calculate patchID for mesh points
        labelList pointPatchID(mesh.nPoints(), -1);
        boolList markedPatches(patches.size(), false);
        forAll(patches, patchI)
        {
            if (!patches[patchI].coupled() && numLayers[patchI] != -1)
            {
                const polyPatch& gP = patches[patchI];
                label start = gP.start();
                forAll(gP, gPI)
                {
                    label meshFaceI = start + gPI;
                    const face f = mesh.faces()[meshFaceI];
                    forAll(f,fp)
                    {
                        pointPatchID[f[fp]] = patchI;
                    }
                }
            }
        }
        syncTools::syncPointList
        (
            mesh,
            pointPatchID,
            maxEqOp<label>(), // combine op
            label(-1)     // null value
        );

        forAll(patches, patchI)
        {
            if (!patches[patchI].coupled() && numLayers[patchI] == -1)
            {
                markedPatches[patchI] = true;
                const polyPatch& guP = patches[patchI];
                label start = guP.start();
                forAll(guP, gPI)
                {
                    label meshFaceI = start + gPI;
                    const face f = mesh.faces()[meshFaceI];
                    forAll(f,fp)
                    {
                        label meshPointI = f[fp];
                        if (pointPatchID[meshPointI] == -1)
                        {
                            pointPatchID[meshPointI] = patchI;
                        }
                    }
                }
            }
        }
        syncTools::syncPointList
        (
            mesh,
            pointPatchID,
            maxEqOp<label>(), // combine op
            label(-1)     // null value
        );

        //Calculate grown up patch normals
        const vectorField& faceNormals = ppCellRemoval.faceNormals();
        pointField pointNormals(ppCellRemoval.nPoints(), vector::zero);
        scalarField nPointFaces(ppCellRemoval.nPoints(), 0.);

        forAll(faceNormals, faceI)
        {
            const face& f = ppCellRemoval.localFaces()[faceI];

            forAll(f, fp)
            {
                scalar fA = mag
                (
                    mesh.faceAreas()[ppCellRemoval.addressing()[faceI]]
                );
                pointNormals[f[fp]] += faceNormals[faceI]*fA;
                nPointFaces[f[fp]] += fA;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            ppCellRemoval.meshPoints(),
            pointNormals,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh,
            ppCellRemoval.meshPoints(),
            nPointFaces,
            plusEqOp<scalar>(),
            scalar(0.)          // null value
        );

        forAll(pointNormals, i)
        {
            if (nPointFaces[i] > 0)
            {
                pointNormals[i] /= nPointFaces[i];
            }
        }
        pointNormals /= (mag(pointNormals) + SMALL);

        labelList cellRemovalID(mesh.nCells(), -1);

        const volScalarField& layerCells =
            mesh.lookupObject<volScalarField>("layerStacks");

        boolList markedCells(mesh.nCells(), false);

        forAll(ppCellRemoval, i)
        {
            label meshFaceI = ppCellRemoval.addressing()[i];
            label patchI = mesh.boundaryMesh().whichPatch(meshFaceI);
            label own = mesh.faceOwner()[meshFaceI];
            cellRemovalID[own] = patchI;
            markedCells[own]  = true;

            point fC = ppCellRemoval.faceCentres()[i];
            vector fN = ppCellRemoval.faceNormals()[i];

            const labelList& cFaces = mesh.cells()[own];
            forAll(cFaces, cFI)
            {
                label cFaceI = cFaces[cFI];
                originatingPatch[cFaceI] = patchI;
                origCentre[cFaceI] = fC;
                origNormal[cFaceI] = fN;
            }
        }

        //Remove cells not connected by a face
        forAll(ppCellRemoval.meshPoints(), ptI)
        {
            label meshPointI = ppCellRemoval.meshPoints()[ptI];
            label patchI = pointPatchID[meshPointI];
            if (patchI != -1 && markedPatches[patchI])
            {
                const labelList& pCells = mesh.pointCells()[meshPointI];
                forAll(pCells, pCI)
                {
                    label cellI = pCells[pCI];
                    if (cellRemovalID[cellI] != -1)
                    {
                        continue;
                    }
                    cellRemovalID[cellI] = patchI;
                    markedCells[cellI]  = true;

                    const labelList& cFaces = mesh.cells()[cellI];
                    forAll(cFaces, cFI)
                    {
                        label cFaceI = cFaces[cFI];
                        originatingPatch[cFaceI] = patchI;
                        origCentre[cFaceI] = mesh.points()[meshPointI];
                        origNormal[cFaceI] = pointNormals[ptI];
                    }
                }
            }
        }

        const labelList& owners = mesh.faceOwner();

        // Calculate coupled layerID
        scalarField neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());

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

        scalarField neiOriginatingPatch(mesh.nFaces()-mesh.nInternalFaces());
        pointField neiOrigCentre(mesh.nFaces()-mesh.nInternalFaces());
        pointField neiOrigNormal(mesh.nFaces()-mesh.nInternalFaces());

        while (true)
        {
            label nMarked = 0;

            for
            (
                label faceI = mesh.nInternalFaces();
                faceI < mesh.nFaces();
                faceI++
            )
            {
                label bFaceI = faceI-mesh.nInternalFaces();

                neiOriginatingPatch[bFaceI] = originatingPatch[faceI];
                neiOrigCentre[bFaceI] = origCentre[owners[faceI]];
                neiOrigNormal[bFaceI] = origNormal[owners[faceI]];
            }
            syncTools::swapBoundaryFaceList(mesh, neiOriginatingPatch);
            syncTools::swapBoundaryFacePositions(mesh, neiOrigCentre);
            syncTools::swapBoundaryFacePositions(mesh, neiOrigNormal);

            //TODO loop only through innerGrid faces //GGG
            forAll(mesh.faces(), faceI)
            {
                label own = owners[faceI];
                label nei = -1;

                scalar layerOwnID = layerCells[own];
                scalar layerNeiID = -1;

                if (mesh.isInternalFace(faceI))
                {
                    nei = mesh.faceNeighbour()[faceI];
                    layerNeiID = layerCells[nei];
                }
                else if (patches[patches.whichPatch(faceI)].coupled())
                {
                    layerNeiID  = neiLayerCells[faceI-mesh.nInternalFaces()];
                }

                if
                (
                    layerOwnID != -1 && layerNeiID != -1
                    && (layerOwnID == layerNeiID)
                )
                {
                    if (mesh.isInternalFace(faceI))
                    {
                        if (originatingPatch[faceI] != -1)
                        {
                            if (!markedCells[own])
                            {
                                markedCells[own] = true;
                                const labelList& cFaces = mesh.cells()[own];
                                forAll(cFaces, cFI)
                                {
                                    originatingPatch[cFaces[cFI]] =
                                        originatingPatch[faceI];
                                    origNormal[cFaces[cFI]] = origNormal[faceI];
                                    origCentre[cFaces[cFI]] = origCentre[faceI];
                                }
                                nMarked++;
                            }
                            else if (!markedCells[nei])
                            {
                                markedCells[nei] = true;
                                const labelList& cFaces = mesh.cells()[nei];
                                forAll(cFaces, cFI)
                                {
                                    originatingPatch[cFaces[cFI]] =
                                        originatingPatch[faceI];
                                    origNormal[cFaces[cFI]] = origNormal[faceI];
                                    origCentre[cFaces[cFI]] = origCentre[faceI];
                                }
                                nMarked++;
                            }
                        }
                    }
                    else if (patches[patches.whichPatch(faceI)].coupled())
                    {
                        label neiOrigPatch =
                            neiOriginatingPatch[faceI-mesh.nInternalFaces()];

                        if (neiOrigPatch != -1)
                        {
                            if (!markedCells[own])
                            {
                                markedCells[own] = true;
                                const labelList& cFaces = mesh.cells()[own];
                                forAll(cFaces, cFI)
                                {
                                    originatingPatch[cFaces[cFI]] =
                                        neiOrigPatch;
                                    origNormal[cFaces[cFI]] =
                                        neiOrigNormal[faceI];
                                    origCentre[cFaces[cFI]] =
                                        neiOrigCentre[faceI];
                                }
                                nMarked++;
                            }
                        }
                    }
                }
            }

            if (returnReduce(nMarked, sumOp<label>()) == 0)
            {
                break;
            }
        }


        syncTools::syncFaceList
        (
            mesh,
            originatingPatch,
            maxEqOp<label>()
        );

        syncTools::syncFaceList
        (
            mesh,
            origCentre,
            minMagSqrEqOp<point>()
        );

        syncTools::syncFaceList
        (
            mesh,
            origNormal,
            minMagSqrEqOp<point>()
        );

        //Remove already added surface layers if not required
        for (int nIter = 0; nIter < nDualLayers; nIter++)
        {
            labelList neiCellRemovalID(mesh.nFaces()-mesh.nInternalFaces(),-1);
            labelList neiLayerCells(mesh.nFaces()-mesh.nInternalFaces(),-1);
            forAll(patches, patchI)
            {
                const polyPatch& patch = patches[patchI];

                const labelUList& faceCells = patch.faceCells();

                label bFaceI = patch.start()-mesh.nInternalFaces();

                if (patch.coupled())
                {
                    forAll(faceCells, i)
                    {
                        neiCellRemovalID[bFaceI] = cellRemovalID[faceCells[i]];
                        neiLayerCells[bFaceI] = layerCells[faceCells[i]];
                        bFaceI++;
                    }
                }
            }
            syncTools::swapBoundaryFaceList(mesh, neiCellRemovalID);
            syncTools::swapBoundaryFaceList(mesh, neiLayerCells);

            label nChanged =  0;

            //TODO loop only through innerGrid faces //GGG
            forAll(mesh.faces(), faceI)
            {
                label own = owners[faceI];
                label ownID = cellRemovalID[own];
                label ownStackID = layerCells[own];
                label neiID = -1;
                label neiStackID = -1;

                if (mesh.isInternalFace(faceI))
                {
                    label nei = mesh.faceNeighbour()[faceI];
                    neiID = cellRemovalID[nei];
                    neiStackID = layerCells[nei];

                    if (ownID == -1 && neiID != -1)
                    {
                        if (ownStackID == neiStackID)
                        {
                            cellRemovalID[own] = neiID;
                            nChanged++;
                        }
                    }
                    else if (ownID != -1 && neiID == -1)
                    {
                        if (ownStackID == neiStackID)
                        {
                            cellRemovalID[nei] = ownID;
                            nChanged++;
                        }
                    }
                }
                else
                {
                    neiID = neiCellRemovalID[faceI-mesh.nInternalFaces()];
                    neiStackID = neiLayerCells[faceI-mesh.nInternalFaces()];
                    if (ownID == -1 && neiID != -1)
                    {
                        if (ownStackID == neiStackID)
                        {
                            cellRemovalID[own] = neiID;
                            nChanged++;
                        }
                    }
                }
            }

            if (returnReduce(nChanged, sumOp<label>()) == 0)
            {
                break;
            }
        }

        DynamicList<label> cellsToRemove(nDualLayers*ppCellRemoval.size());
        labelList exposedPatchID(mesh.nFaces(), -1);

        //TODO loop only through innerGrid cells //GGG
        forAll(mesh.cells(), cellI)
        {
            label patchI = cellRemovalID[cellI];
            if (patchI != -1)
            {
                cellsToRemove.append(cellI);
                cell c = mesh.cells()[cellI];
                forAll(c, cFI)
                {
                    exposedPatchID[c[cFI]] = patchI;
                }
            }
        }
        cellsToRemove.shrink();

        syncTools::syncFaceList
        (
            mesh,
            exposedPatchID,
            maxEqOp<label>()
        );

        removeCells cellRemover(mesh);
        labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
        labelList exposedPatches(exposedFaces.size());

        forAll(exposedFaces, i)
        {
            label faceI = exposedFaces[i];
            exposedPatches[i] = exposedPatchID[faceI];
        }

        polyTopoChange meshMod(mesh);

        cellRemover.setRefinement
        (
            cellsToRemove,
            exposedFaces,
            exposedPatches,
            meshMod
        );

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

        resetPrimitivePatchAddressing(adaptPatchIDs,pp);

        origPatchNewMesh.setSize(mesh.nFaces(), -1);
        origCentreNewMesh.setSize(mesh.nFaces(), greatPoint);
        origNormalNewMesh.setSize(mesh.nFaces(), greatPoint);

        //TODO loop only through innerGrid faces //GGG
        forAll(mesh.faces(), faceI)
        {
            label oldFaceI = map().faceMap()[faceI];
            label origPatch = originatingPatch[oldFaceI];

            if (origPatch != -1)
            {
                origPatchNewMesh[faceI] = origPatch;
                origCentreNewMesh[faceI] = origCentre[oldFaceI];
                origNormalNewMesh[faceI] = origNormal[oldFaceI];
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh,
        origPatchNewMesh,
        maxEqOp<label>()
    );

    syncTools::syncFaceList
    (
        mesh,
        origCentreNewMesh,
        minMagSqrEqOp<point>()
    );

    syncTools::syncFaceList
    (
        mesh,
        origNormalNewMesh,
        minMagSqrEqOp<point>()
    );

    {
        const pointMesh& pMeshSnap = pointMesh::New(mesh);

        autoPtr<indirectPrimitivePatch> ppCellRemovalPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                noLayerPatchIDs
            )
        );
        indirectPrimitivePatch& ppCellRemoval = ppCellRemovalPtr();

        dictionary updatedMotionDict(motionDict);
        updatedMotionDict.add("minTwist",-1,true);

        motionSmoother meshMoverSnap
        (
            mesh,
            ppCellRemoval,
            noLayerPatchIDs,
            meshRefinement::makeDisplacementField(pMeshSnap, noLayerPatchIDs),
            updatedMotionDict
        );

        const pointField& localPoints = ppCellRemoval.localPoints();
        vectorField patchDisp(localPoints.size(), vector::zero);

        pointField curPts(localPoints);
        pointField curFCs(ppCellRemoval.size(),vector::zero);

        const labelList meshEdges
        (
            ppCellRemoval.meshEdges(mesh.edges(), mesh.pointEdges())
        );
        const labelListList& edgeFaces = ppCellRemoval.edgeFaces();

        labelList nExternalEdge(mesh.nEdges(), 0);
        forAll(ppCellRemoval.edges(), edgeI)
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

        boolList stationaryPoints(ppCellRemoval.localPoints().size(), false);

        forAll(ppCellRemoval.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            if (nExternalEdge[meshEdgeI] == 1)
            {
                const edge& e = ppCellRemoval.edges()[edgeI];
                stationaryPoints[e[0]] = true;
                stationaryPoints[e[1]] = true;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            ppCellRemoval.meshPoints(),
            stationaryPoints,
            orEqOp<bool>(), // combine op
            false     // null value
        );

        pointField prevPts(curPts.size(), vector::zero);
        labelList weights(curPts.size(), 0);

        label maxIter = 10;//6;
        for (label iter = 0; iter < maxIter; iter++)
        {
            prevPts = curPts;
            curPts = vector::zero;
            weights = 0;

            forAll(localPoints, i)
            {
                if (stationaryPoints[i])
                {
                    continue;
                }

                const labelList& pFaces = ppCellRemoval.pointFaces()[i];

                forAll(pFaces, pFI)
                {
                    label meshFaceI = ppCellRemoval.addressing()[pFaces[pFI]];

                    if (origCentreNewMesh[meshFaceI] != greatPoint)
                    {
                        point fC = origCentreNewMesh[meshFaceI];
                        vector fN = origNormalNewMesh[meshFaceI];
                        if (mag(fN) < VSMALL)
                        {
                            continue;
                        }
                        plane snapPlane(fC,fN);
                        curPts[i] += snapPlane.nearestPoint(prevPts[i]);
                        weights[i]++;
                     }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                ppCellRemoval.meshPoints(),
                curPts,
                plusEqOp<point>(), // combine op
                vector::zero     // null value
            );

            syncTools::syncPointList
            (
                mesh,
                ppCellRemoval.meshPoints(),
                weights,
                plusEqOp<label>(), // combine op
                label(0)     // null value
            );

            forAll(localPoints, i)
            {
                if (weights[i] > 0)
                {
                    curPts[i] /= weights[i];
                }
                else
                {
                   curPts[i] = localPoints[i];
                }
            }

            if (iter == maxIter -1)
            {
                break;
            }

            forAll(ppCellRemoval, i)
            {
                curFCs[i] = ppCellRemoval.localFaces()[i].centre(curPts);
            }

            curPts = vector::zero;
            weights = 0;

            forAll(localPoints, i)
            {
                if (stationaryPoints[i])
                {
                    curPts[i] = localPoints[i];
                    continue;
                }

                const labelList& pFaces = ppCellRemoval.pointFaces()[i];

                forAll(pFaces, pFI)
                {
                    label faceI = pFaces[pFI];
                    curPts[i] += curFCs[faceI];
                    weights[i]++;
                }
            }

            syncTools::syncPointList
            (
                mesh,
                ppCellRemoval.meshPoints(),
                curPts,
                plusEqOp<point>(), // combine op
                vector::zero     // null value
            );

            syncTools::syncPointList
            (
                mesh,
                ppCellRemoval.meshPoints(),
                weights,
                plusEqOp<label>(), // combine op
                label(0)     // null value
            );

            forAll(localPoints, i)
            {
                if (weights[i] > 0)
                {
                    curPts[i] /= weights[i];
                }
            }
        }

        patchDisp = curPts - localPoints;

        syncTools::syncPointList
        (
            mesh,
            ppCellRemoval.meshPoints(),
            patchDisp,
            minMagSqrEqOp<point>(),         // combine op
            greatPoint     // null value
        );

        dictionary meshOptimDict
        (
            meshDict_.found("meshOptimization") ?
            meshDict_.subDict("meshOptimization") :
            dictionary()
        );

        if (meshOptimDict.found("foamOptimizeCoeffs"))
        {
            dictionary& coeffsDict =
                meshOptimDict.subDict("foamOptimizeCoeffs");
            if (!coeffsDict.found("meshQualityControls"))
            {
                coeffsDict.add("meshQualityControls",motionDict,true);
            }
        }
        else if (meshOptimDict.found("cfMeshOptimizeCoeffs"))
        {
            dictionary& coeffsDict =
                meshOptimDict.subDict("cfMeshOptimizeCoeffs");

            scalar meshDictOrtho = motionDict.lookupOrDefault<scalar>
            (
                "maxNonOrtho",
                scalar(70),
                true
            );

            if (!coeffsDict.found("maxNonOrtho"))
            {
                coeffsDict.add("maxNonOrtho",meshDictOrtho,true);
            }
            else
            {
                scalar origOrth = readScalar
                (
                    coeffsDict.lookup("maxNonOrtho")
                );
                if (origOrth > 180-SMALL)
                {
                    coeffsDict.add("maxNonOrtho",meshDictOrtho,true);
                }
            }
        }

        vectorField pointDisp(mesh.nPoints(), vector::zero);
        forAll(patchDisp, ptI)
        {
            label meshPointI = ppCellRemoval.meshPoints()[ptI];
            pointDisp[meshPointI] = patchDisp[ptI];
        }

        syncTools::syncPointList
        (
            mesh,
            pointDisp,
            maxMagSqrEqOp<point>(),
            vector::zero
        );

        pointField newPoints = meshMoverSnap.oldPoints();
        newPoints += pointDisp;

        autoPtr<autoOptimize> optimMeshPtr
            = autoOptimize::New(mesh, meshOptimDict);
        optimMeshPtr->movePoints(newPoints);

        pointVectorField& disp = meshMoverSnap.displacement();
        disp.primitiveFieldRef() = newPoints - meshMoverSnap.oldPoints();

        forAll(patchDisp, ptI)
        {
            label meshPointI = ppCellRemoval.meshPoints()[ptI];
            patchDisp[ptI] = newPoints[meshPointI] - localPoints[ptI];
        }

        syncTools::syncPointList
        (
            mesh,
            ppCellRemoval.meshPoints(),
            patchDisp,
            minMagSqrEqOp<point>(), // combine op
            greatPoint     // null value
        );

        // Set initial distribution of displacement field (on patches) from
        // patchDisp and make displacement consistent with b.c. on displacement
        // pointVectorField.
        meshMoverSnap.setDisplacement(patchDisp);

        // Apply internal displacement to mesh.
        scaleMesh(snapParams, nInitErrors, List<labelPair>(0), meshMoverSnap);
    }

    {
        const pointMesh& pMeshSnap = pointMesh::New(mesh);

        autoPtr<indirectPrimitivePatch> ppFtrPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                noLayerPatchIDs//featureIDs
             )
        );
        indirectPrimitivePatch& ppFtr = ppFtrPtr();

        motionSmoother meshMoverSnap
        (
            mesh,
            pp,
            adaptPatchIDs,
            meshRefinement::makeDisplacementField(pMeshSnap, adaptPatchIDs),
            motionDict
        );

        labelList snapSurfaces =
            surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
            (
                surfaces.surfZones()
            );

        label nFeatureIter = 10;
        optStages ostages;
        if (controller_.algorithm() == meshControl::EXTRUDE)
        {
            if (controller_.mode() == meshControl::FAST)
            {
                nFeatureIter = 5;
                direction ot = optStages::ERROR;
                ostages.add(nFeatureIter-1,ot);
            }
            else
            {
                direction ot = optStages::ERROR | optStages::SMOOTH;
                ostages.add(nFeatureIter-1,ot);
            }
        }

        // Snap faces to lock on to feature edges
        featureSnap
        (
            nFeatureIter,
            5,
            snapSurfaces,
            ppFtr,
            snapParams,
            motionDict,
            meshMoverSnap,
            false, //initial snap
            ostages
         );

        //perform volumetric smoothing
        if (controller_.algorithm() != meshControl::EXTRUDE)
        {
            medialAxisMeshMover::smoothDisplacement
            (
                motionDict,
                labelList(0),
                meshRefiner_.meshCutter().level0EdgeLength(),

                mesh,
                meshMoverSnap.curPoints()(),
                meshMoverSnap.displacement()
            );
        }

        scaleMesh(snapParams, nInitErrors, List<labelPair>(0), meshMoverSnap);
        meshMoverSnap.correct();
    }

    //optimize the mesh
    {
        dictionary meshOptimDict
        (
            meshDict_.found("meshOptimization") ?
            meshDict_.subDict("meshOptimization") :
            dictionary()
        );

        if (meshOptimDict.found("foamOptimizeCoeffs"))
        {
            dictionary& coeffsDict =
                meshOptimDict.subDict("foamOptimizeCoeffs");
            if (!coeffsDict.found("meshQualityControls"))
            {
                coeffsDict.add("meshQualityControls",motionDict,true);
            }
        }

        autoPtr<autoOptimize> optimMeshPtr
            = autoOptimize::New(mesh, meshOptimDict);
        optimMeshPtr->optimize();
    }

    {
        //For final snap re-introduce zone baffles
        List<labelPair> baffles;
        {
            labelList originatingFaceZone;

            List<surfaceZonesInfo::faceZoneType> fzType(2);
            fzType[0] = surfaceZonesInfo::INTERNAL;
            fzType[1] = surfaceZonesInfo::BAFFLE;

            const labelList zonesToBaffle
            (
                meshRefiner_.getZones(fzType)
            );

            meshRefiner_.createZoneBaffles
            (
                zonesToBaffle,//identity(mesh.faceZones().size()),
                baffles,
                originatingFaceZone
            );
        }

        const pointMesh& pMeshSnap = pointMesh::New(mesh);

        autoPtr<indirectPrimitivePatch> ppFtrPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                adaptPatchIDs//noLayerPatchIDs
             )
         );
        indirectPrimitivePatch& ppFtr = ppFtrPtr();

        motionSmoother meshMoverSnap
        (
            mesh,
            pp,
            adaptPatchIDs,
            meshRefinement::makeDisplacementField(pMeshSnap, adaptPatchIDs),
            motionDict
        );

        labelList snapSurfaces = identity(surfaces.surfaces().size());

        label nFeatureIter = 10;
        optStages ostages;
        if (controller_.algorithm() == meshControl::EXTRUDE)
        {
            if (controller_.mode() == meshControl::FAST)
            {
                nFeatureIter = 5;
            }
            direction ot = optStages::ERROR;
            ostages.add(nFeatureIter-1,ot);
        }

        // Snap faces to lock on to feature edges
        featureSnap
        (
            nFeatureIter,
            5,
            snapSurfaces,
            ppFtr,
            snapParams,
            motionDict,
            meshMoverSnap,
            false, //initial snap
            ostages
         );

        //perform volumetric smoothing
        if (controller_.algorithm() != meshControl::EXTRUDE)
        {
            medialAxisMeshMover::smoothDisplacement
            (
                motionDict,
                labelList(0),
                meshRefiner_.meshCutter().level0EdgeLength(),

                mesh,
                meshMoverSnap.curPoints()(),
                meshMoverSnap.displacement()
            );
        }

        scaleMesh(snapParams, nInitErrors, baffles, meshMoverSnap);
        meshMoverSnap.correct();

        //Re merge zone baffles
        {
            autoPtr<mapPolyMesh> map = meshRefiner_.mergeZoneBaffles
            (
                true,   // internal zones
                true   // baffle zones
            );
        }
    }
}


void Foam::snappySnapDriver::snapAndMerge
(
    snapParameters& snapParams
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict_.subDict("meshQualityControls");

    // Get the labels of added patches.
    labelList adaptPatchIDs(meshRefiner_.meshedPatches());
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatchIDs
         )
    );
    indirectPrimitivePatch& pp = ppPtr();

    const pointMesh& pMeshSnapNearest = pointMesh::New(mesh);
    motionSmoother meshMoverSnapNearest
    (
        mesh,
        pp,
        adaptPatchIDs,
        meshRefinement::makeDisplacementField
        (
            pMeshSnapNearest,
            adaptPatchIDs
        ),
        motionDict
    );

    // Pre-smooth patch vertices (so before determining nearest)
    preSmoothPatch
    (
        meshRefiner_,
        snapParams,
        label(0),
        List<labelPair>(0),
        meshMoverSnapNearest
    );

    meshMoverSnapNearest.correct();

    // Calculate displacement at every patch point. Insert into
    // meshMover.
    calcNearestSurface
    (
        calcSnapDistance(meshRefiner_, snapParams, pp),
        List<labelPair>(0),
        meshMoverSnapNearest
    );

    // Apply internal displacement to mesh.
    scaleMesh
    (
        snapParams,
        label(0),
        List<labelPair>(0),
        meshMoverSnapNearest
    );

    //update pp addressing
    resetPrimitivePatchAddressing(adaptPatchIDs,pp);

    mergePatchFacesUndo
    (
        motionDict,
        List<labelList>(0),
        true
     );

    return;
}


void Foam::snappySnapDriver::testSnapPerformance
(
    const indirectPrimitivePatch& pp,
    const snapParameters& snapParams
)
{
    Info<<"Checking snap performance of individual surfaces" << nl << endl;

    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    labelList snapSurfaces = identity(surfaces.surfaces().size());
    List<pointIndexHit> hitInfo;
    labelList hitSurface;
    scalarField maxSnapDist =
        calcSnapDistance(meshRefiner_, snapParams, pp);

    scalar minPPsz = returnReduce(pp.size(), minOp<label>());
    scalar maxPPSz = returnReduce(pp.size(), maxOp<label>());
    Info<<"Snap patch size processor (min) : "
        <<  minPPsz << " processor (max) : "<< maxPPSz << nl << endl;

    scalar totalSnapTime(0);
    forAll(snapSurfaces, surfi)
    {
        scalar snapStartTime = mesh.time().elapsedCpuTime();
        labelList snapSurf(1, snapSurfaces[surfi]);
        surfaces.findNearest
        (
            snapSurf,
            pp.faceCentres(),
            sqr(maxSnapDist), // sqr of attract
            hitSurface,
            hitInfo
        );
        scalar snapTime = mesh.time().elapsedCpuTime()
            - snapStartTime;

        totalSnapTime += snapTime;
        scalar minSnapTime = returnReduce(snapTime, minOp<scalar>());
        scalar maxSnapTime = returnReduce(snapTime, maxOp<scalar>());

        Info<< "Snap surface : " << surfaces.names()[surfi]
            << " : time (processsor min) : " << minSnapTime
            << " : time (processor max) : " << maxSnapTime
            << " : processor balance : " << minSnapTime/(maxSnapTime+SMALL)
            <<endl;
    }

    scalar minTotalSnapTime = returnReduce(totalSnapTime, minOp<scalar>());
    scalar maxTotalSnapTime = returnReduce(totalSnapTime, maxOp<scalar>());

    Info<<nl<<"Overall test snapping"
        << " : time (processor min) " << minTotalSnapTime
        << " : time (processor max) " << maxTotalSnapTime << nl << endl;

    return;
}


void Foam::snappySnapDriver::doSnap
(
    snapParameters& snapParams,
    const refinementParameters& refineParams
)
{
    bool dryRun = false;
    if (controller_.mode() == meshControl::DRYRUN)
    {
        dryRun = true;
    }
    else if
    (
        controller_.mode() == meshControl::FAST
        || controller_.algorithm() == meshControl::EXTRUDE
    )
    {
        scalar& snapTol = snapParams.snapTol();
        if (snapTol > 2)
        {
            snapTol = scalar(2);
        }
    }

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict_.subDict("meshQualityControls");

    // snap-to-surface parameters
    const dictionary& snapDict = meshDict_.subDict("snapControls");

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(snap, "snappyHexMesh::snap");
    #endif

    fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    scalar featureCos = surfaces.minFeatureRefineAngle();
    scalar planarAngle = refineParams.planarAngle();

    meshRefiner_.features()().setCheckRefinementOnly(true);

    bool writeSnapVTK(dryRun ? false : bool(snapParams.writeSnapVTK()));

    Info<< nl
        << "Morphing phase" << nl
        << "--------------" << nl
        << endl;

    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    // faceZone handling
    // ~~~~~~~~~~~~~~~~~
    //
    // We convert all faceZones into baffles during snapping so we can use
    // a standard mesh motion (except for the mesh checking which for baffles
    // created from internal faces should check across the baffles). The state
    // is stored in two variables:
    //      baffles : pairs of boundary faces
    //      duplicateFace : from mesh face to its baffle colleague (or -1 for
    //                      normal faces)
    // There are three types of faceZones according to the faceType property:
    //
    // internal
    // --------
    // - baffles: need to be checked across
    // - duplicateFace: from face to duplicate face. Contains
    //   all faces on faceZone to prevents merging patch faces.
    //
    // baffle
    // ------
    // - baffles: no need to be checked across
    // - duplicateFace: contains all faces on faceZone to prevent
    //   merging patch faces.
    //

    // faceZones of type internal
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

    bool internalOrBaffleZone = false;

    // Create baffles (pairs of faces that share the same points)
    // Baffles stored as owner and neighbour face that have been created.
    {
        List<labelPair> baffles;
        labelList originatingFaceZone;

        List<surfaceZonesInfo::faceZoneType> fzType(2);

        fzType[0] = surfaceZonesInfo::INTERNAL;
        fzType[1] = surfaceZonesInfo::BAFFLE;

        const labelList zonesToBaffle
        (
            meshRefiner_.getZones(fzType)
        );

        meshRefiner_.createZoneBaffles
        (
            zonesToBaffle,//identity(mesh.faceZones().size()),
            baffles,
            originatingFaceZone
        );

        if (returnReduce(baffles.size(), sumOp<label>()) > 0)
        {
            internalOrBaffleZone = true;
        }
    }

    // Get the labels of added patches.
    labelList adaptPatchIDs(meshRefiner_.meshedPatches());

    // Duplicate points on faceZones of type boundary
    meshRefiner_.dupNonManifoldBoundaryPoints();

    if (refineParams.splitCells() && snapParams.addTetsToSplitMesh())
    {
        labelList zonePatches = zonePatchIDs();
        labelHashSet excludePatches(zonePatches);

        DynamicList<label> tetSplitPatches(adaptPatchIDs.size());

        forAll(adaptPatchIDs, i)
        {
            if (!excludePatches.found(adaptPatchIDs[i]))
            {
                tetSplitPatches.append(adaptPatchIDs[i]);
            }
        }
        tetSplitPatches.shrink();

        addTetrahedralCellsToSplitMesh(tetSplitPatches);
    }

    //Create initial volume field if checking cell distortion during snapping
    autoPtr<volScalarField> preSnapVol;
    scalar minSnapRelativeVolume =
        motionDict.lookupOrDefault<scalar>("minSnapRelativeVolume", -1, true);
    scalar minSnapRelativeTetVolume =
        motionDict.lookupOrDefault<scalar>("minSnapRelativeTetVolume", -1, true);
    scalar maxGaussGreenCentroid =
        motionDict.lookupOrDefault<scalar>("maxGaussGreenCentroid", -1, true);


    if
    (
        minSnapRelativeVolume > SMALL || minSnapRelativeTetVolume > SMALL
     || maxGaussGreenCentroid > SMALL
    )
    {
        preSnapVol.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "preSnapVolume",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                mesh,
                dimensionedScalar("scalar", dimless, 0),
                zeroGradientFvPatchScalarField::typeName
            )
        );
        preSnapVol().primitiveFieldRef() = mesh.V();
    }

    // Check initial mesh
    Info<< "Checking initial mesh ..." << endl;
    labelHashSet wrongFaces(mesh.nFaces()/100);
    label errorsFound = motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);

    label nInitErrors = 0;
    if
    (
        controller_.algorithm() == meshControl::STANDARD
     || controller_.algorithm() == meshControl::SHELL
    )
    {
        nInitErrors = errorsFound;

        Info<< "Detected " << nInitErrors << " illegal faces"
            << " (concave, zero area or negative cell pyramid volume)"
            << endl;
    }

    Info<< "Checked initial mesh in = "
        << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;


    bool doFeatures = false;
    label nPreFeatIter = 0;
    if (snapParams.nPreFeatureIter() > 0)
    {
        doFeatures = true;
        nPreFeatIter = snapParams.nPreFeatureIter();

        Info<< "Snapping to features in " << nPreFeatIter
            << " iterations ..." << endl;
    }

    bool  meshOk = false;

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            adaptPatchIDs
         )
     );
    indirectPrimitivePatch& pp = ppPtr();

    // Per patch whether to feature snap
    const boolList& featureEdgePatches = snapParams.featureEdgePatches();

    // Patches that need to feature snap
    DynamicList<label> featureIDs(featureEdgePatches.size());
    forAll(featureEdgePatches, patchI)
    {
        if (featureEdgePatches[patchI] == true)
        {
            featureIDs.append(patchI);
        }
    }
    featureIDs.shrink();

    // Extract baffles across internal faceZones (for checking mesh quality
    // across
    labelPairList internalBaffles
    (
        meshRefiner_.subsetBaffles
        (
            mesh,
            internalFaceZones,
            localPointRegion::findDuplicateFacePairs(mesh)
        )
    );

    //Pre smooth and merge patch faces
    if (controller_.algorithm() != meshControl::DUAL)
    {
        autoPtr<motionSmoother> meshMoverPtr
        (
            new motionSmoother
            (
                mesh,
                pp,
                adaptPatchIDs,
                meshRefinement::makeDisplacementField
                (
                    pointMesh::New(mesh),
                    adaptPatchIDs
                ),
                motionDict
            )
        );

        // Pre-smooth patch vertices (so before determining nearest)
        preSmoothPatch
        (
           meshRefiner_,
           snapParams,
           nInitErrors,
           internalBaffles,
           meshMoverPtr()
        );

        if
        (
            (
                controller_.algorithm() == meshControl::STANDARD
             || controller_.algorithm() == meshControl::SHELL
            )
         && controller_.topoChanges() && snapParams.mergeBoundaryFaces()
        )
        {
            autoPtr<indirectPrimitivePatch> ppFtrPtr
            (
                meshRefinement::makePatch
                (
                    mesh,
                    featureIDs
                )
            );
            indirectPrimitivePatch& ppFtr = ppFtrPtr();

            // Distance to attract to nearest feature on surface
            const scalarField snapDist
            (
                calcSnapDistance(meshRefiner_, snapParams, ppFtr)
            );
            mergePatchFacesUndo
            (
                motionDict,
                snapDist,
                ppFtr,
                internalBaffles
            );

            resetPrimitivePatchAddressing(adaptPatchIDs,pp);
        }
    }

    if (writeSnapVTK)
    {
        simpleVTKWriter
        (
            pp.localFaces(),
            pp.localPoints()
         ).write("snappyFeatureFirstMerge.vtk");
    }

    boolList regionSnap(globalToMasterPatch_.size(), false);
    labelHashSet regionSnappedPatches
    (
        UList<label>(snapParams.regionSnappedIDs())
    );

    forAll(regionSnap, i)
    {
        label patchI = globalToMasterPatch_[i];

        if (regionSnappedPatches.found(patchI))
        {
            regionSnap[i] = true;
        }
    }

    label nFeatIter = snapParams.nFeatureIter();

    for (label i = 0; i < snapParams.nOuterIter(); i++)
    {
//#ifdef FOAM_USE_TBB
//        Timer timer("for nOuterIter()");
//#endif
        autoPtr<indirectPrimitivePatch> ppFtrPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                featureIDs
            )
        );

        indirectPrimitivePatch& ppFtr = ppFtrPtr();
        const pointMesh& pMeshSnap = pointMesh::New(mesh);

        if (i == 0 && snapParams.testPerformance())
        {
            testSnapPerformance(ppFtr,snapParams);
        }

        //Do an optional initial snap of non-zoned patches
        if (i == 0 && snapParams.preZoneSnap() && !dryRun)
        {
//#ifdef FOAM_USE_TBB
//            Timer timer("optional initial snap of non-zoned patches");
//#endif
            labelList zonePatches = zonePatchIDs();
            if (zonePatches.size())
            {
                //First zoned surfaces
                for (label iter = 0; iter < 2; iter++)
                {
                    labelList initSnapSurfaces;
                    labelList initSnapPatches;
                    if (iter == 0)
                    {
                        initSnapSurfaces =
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
                                initSnapPatches[nInitPatches++] =
                                    adaptPatchIDs[i];
                            }
                        }
                        initSnapPatches.setSize(nInitPatches);
                    }
                    else
                    {
                        initSnapSurfaces =
                            surfaceZonesInfo::getNonBoundaryNamedSurfaces
                            (
                                surfaces.surfZones()
                            );
                        initSnapPatches = zonePatchIDs();
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

                    motionSmoother meshMoverInitSnap
                    (
                        mesh,
                        ppInit,
                        initSnapPatches,
                        meshRefinement::makeDisplacementField
                        (
                            pMeshSnap,
                            initSnapPatches
                        ),
                        motionDict
                    );
                    label maxPreZoneSnap = 4;
                    optStages ostages;
                     // Snap faces to lock on to feature edges
                    featureSnap
                    (
                        maxPreZoneSnap,
                        (iter == 1 ? maxPreZoneSnap : 0),
                        initSnapSurfaces,
                        ppInit,
                        snapParams,
                        motionDict,
                        meshMoverInitSnap,
                        false, //initial snap
                        ostages
                    );

                    if
                    (
                        controller_.algorithm() == meshControl::STANDARD
                     || controller_.algorithm() == meshControl::SHELL
                    )
                    {
                        // Get smoothly varying internal displacement field.
                        medialAxisMeshMover::smoothDisplacement
                        (
                            motionDict,
                            labelList(0),
                            meshRefiner_.meshCutter().level0EdgeLength(),

                            mesh,
                            meshMoverInitSnap.curPoints()(),
                            meshMoverInitSnap.displacement()
                        );
                    }

                    scaleMesh
                    (
                        snapParams,
                        nInitErrors,
                        internalBaffles,
                        meshMoverInitSnap
                    );
                    meshMoverInitSnap.correct();
                }
            }
        } //optional initial snap of non-zoned patches

        if (nFeatIter && returnReduce(ppFtr.size(), sumOp<label>()))
        {
//#ifdef FOAM_USE_TBB
//            Timer timer("snap to features");
//#endif

            motionSmoother meshMoverSnap
            (
                mesh,
                pp,
                adaptPatchIDs,
                meshRefinement::makeDisplacementField(pMeshSnap, adaptPatchIDs),
                motionDict
             );

            //label smoothStartIter = snapParams.nFeatureIter() - 20;
            label smoothStartIter = 0;

            labelList snapSurfaces = identity(surfaces.surfaces().size());

            optStages ostages;
            if (controller_.algorithm() == meshControl::EXTRUDE)
            {
                direction ot = optStages::ERROR | optStages::SMOOTH;
                ostages.add(10, ot);
                ot = optStages::BOUNDARYRELAXEDERROR | optStages::EXTRUDE;
                ostages.add(snapParams.nFeatureIter()-1, ot);
            }
            else if (controller_.algorithm() == meshControl::DUAL)
            {
                direction ot = optStages::ERROR;
                ostages.add(10, ot);
            }

            // Snap faces to lock on to feature edges
            featureSnap
            (
                (dryRun ? 1 : snapParams.nFeatureIter()),
                smoothStartIter,
                snapSurfaces,
                ppFtr,
                snapParams,
                motionDict,
                meshMoverSnap,
                (controller_.algorithm() == meshControl::DUAL),
                ostages
            );

            // Get smoothly varying internal displacement field.
            if (controller_.algorithm() != meshControl::EXTRUDE)
            {
                medialAxisMeshMover::smoothDisplacement
                (
                    motionDict,
                    meshRefiner_.unmeshedPatches(),
                    meshRefiner_.meshCutter().level0EdgeLength(),

                    mesh,
                    meshMoverSnap.curPoints()(),
                    meshMoverSnap.displacement()
                );
            }

            //Disable mesh checks across baffles where checks have been disabled
            PackedList<1> isZonedDisabledFace
            (
                setZonedBaffleCheckFaces(meshRefiner_,refineParams)
            );

            DynamicList<labelPair> checkedBaffles(internalBaffles.size());
            forAll(internalBaffles, i)
            {
                const labelPair& p = internalBaffles[i];
                if (isZonedDisabledFace.get(p[0]) == 0)
                {
                    checkedBaffles.append(p);
                }
            }
            checkedBaffles.shrink();
            scaleMesh(snapParams, nInitErrors, checkedBaffles, meshMoverSnap);

            meshMoverSnap.correct();
        } // snap to features

        const pointMesh& pMeshSnapNearest = pointMesh::New(mesh);
        motionSmoother meshMoverSnapNearest
        (
            mesh,
            pp,
            adaptPatchIDs,
            meshRefinement::makeDisplacementField
            (
                pMeshSnapNearest,
                adaptPatchIDs
            ),
            motionDict
         );

        // Calculate displacement at every patch point. Insert into
        // meshMover.
        calcNearestSurface
        (
            calcSnapDistance(meshRefiner_, snapParams, pp),
            internalBaffles,
            meshMoverSnapNearest
        );

            //if(true)controller_.algorithm() == meshControl::STANDARD
        if
        (
            (
                controller_.algorithm() == meshControl::STANDARD
                || controller_.algorithm() == meshControl::SHELL
            )
            && !dryRun
        )
        {
//#ifdef FOAM_USE_TBB
//            Timer timer("snappySnapDriver::medialAxisMeshMover::smoothDisplacement");
//#endif
            // Get smoothly varying internal displacement field.
            medialAxisMeshMover::smoothDisplacement
            (
                motionDict,
                labelList(0),
                meshRefiner_.meshCutter().level0EdgeLength(),

                mesh,
                meshMoverSnapNearest.curPoints()(),
                meshMoverSnapNearest.displacement()
             );
        }

        // Apply internal displacement to mesh.
        scaleMesh(snapParams, nInitErrors, internalBaffles, meshMoverSnapNearest);

        //update pp addressing
        resetPrimitivePatchAddressing(adaptPatchIDs,pp);

        if
        (
            (
                controller_.algorithm() == meshControl::STANDARD
             || controller_.algorithm() == meshControl::SHELL
            )
         && i==0 && snapParams.nOuterIter() > 1
         && controller_.topoChanges() && snapParams.mergeBoundaryFaces()
         && !dryRun
        )
        {
            // Distance to attract to nearest feature on surface
            const scalarField snapDist
            (
                calcSnapDistance(meshRefiner_, snapParams, pp)
            );

            mergePatchFacesUndo
            (
                motionDict,
                snapDist,
                pp,
                internalBaffles
            );

            //update pp addressing
            resetPrimitivePatchAddressing(adaptPatchIDs,pp);
        }
    } // for snapParams.nOuterIter()

    if (writeSnapVTK)
    {
        simpleVTKWriter
        (
            pp.localFaces(),
            pp.localPoints()
         ).write("snappyFeatureSnap.vtk");
    }

    if (!dryRun)
    {
        autoPtr<motionSmoother> meshMoverPtr
        (
            new motionSmoother
            (
                mesh,
                ppPtr(),
                adaptPatchIDs,
                meshRefinement::makeDisplacementField
                (
                    pointMesh::New(mesh),
                    adaptPatchIDs
                ),
                motionDict
            )
        );

        // Distance to attract to nearest feature on surface
        scalarField snapDist
        (
            calcSnapDistance(meshRefiner_, snapParams, ppPtr())
        );

        //- Only if in feature attraction mode:
        // Nearest feature
        vectorField patchAttraction;
        // Constraints at feature
        List<pointConstraint> patchConstraints;


        //- Any faces to split
        DynamicList<label> splitFaces;
        //- Indices in face to split across
        DynamicList<labelPair> splits;

        for (label iter = 0; iter < nPreFeatIter; iter++)
        {
            Info<< nl
                << "Morph iteration " << iter << nl
                << "-----------------" << endl;

            // Splitting iteration?
            bool doSplit = false;
            if
            (
                doFeatures
             && snapParams.nFaceSplitInterval() > 0
             && (
                    (iter == nFeatIter-1)
                 || (iter > 0 && (iter%snapParams.nFaceSplitInterval()) == 0)
                )
            )
            {
                doSplit = true;
            }

            indirectPrimitivePatch& pp = ppPtr();
            motionSmoother& meshMover = meshMoverPtr();


            // Calculate displacement at every patch point. Insert into
            // meshMover.
            // Calculate displacement at every patch point
            pointField nearestPoint;
            vectorField nearestNormal;

            if (snapParams.detectNearSurfacesSnap())
            {
                nearestPoint.setSize(pp.nPoints(), vector::max);
                nearestNormal.setSize(pp.nPoints(), Zero);
            }

            vectorField disp = calcNearestSurface
            (
                snapParams.strictRegionSnap(),  // attract points to region only
                meshRefiner_,
                globalToMasterPatch_,           // for if strictRegionSnap
                globalToSlavePatch_,            // for if strictRegionSnap
                snapDist,
                pp,

                nearestPoint,
                nearestNormal
            );


            // Override displacement at thin gaps
            if (snapParams.detectNearSurfacesSnap())
            {
                detectNearSurfaces
                (
                    Foam::cos(degToRad(planarAngle)),// planar cos for gaps
                    pp,
                    nearestPoint,   // surfacepoint from nearest test
                    nearestNormal,  // surfacenormal from nearest test

                    disp
                );
            }

            // Override displacement with feature edge attempt
            if (doFeatures)
            {
                splitFaces.clear();
                splits.clear();
                disp = calcNearestSurfaceFeature
                (
                    snapParams,
                    !doSplit,       // alignMeshEdges
                    iter,
                    featureCos,
                    scalar(iter+1)/nFeatIter,

                    snapDist,
                    disp,
                    nearestNormal,
                    meshMover,

                    patchAttraction,
                    patchConstraints,

                    splitFaces,
                    splits
                );
            }

            // Check for displacement being outwards.
            outwardsDisplacement(pp, disp);

            // Set initial distribution of displacement field (on patches)
            // from patchDisp and make displacement consistent with b.c.
            // on displacement pointVectorField.
            meshMover.setDisplacement(disp);


            if (debug&meshRefinement::ATTRACTION)
            {
                dumpMove
                (
                    mesh.time().path()
                  / "patchDisplacement_" + name(iter) + ".obj",
                    pp.localPoints(),
                    pp.localPoints() + disp
                );
            }

            // Get smoothly varying internal displacement field.
            smoothDisplacement(snapParams, meshMover);

            // Apply internal displacement to mesh.
            meshOk = scaleMesh
            (
                snapParams,
                nInitErrors,
                internalBaffles,
                meshMover
            );

            if (!meshOk)
            {
                WarningInFunction
                    << "Did not succesfully snap mesh."
                    << " Continuing to snap to resolve easy" << nl
                    << "    surfaces but the"
                    << " resulting mesh will not satisfy your quality"
                    << " constraints" << nl << endl;
            }

            if (debug&meshRefinement::MESH)
            {
                const_cast<Time&>(mesh.time())++;
                Info<< "Writing scaled mesh to time "
                    << meshRefiner_.timeName() << endl;
                meshRefiner_.write
                (
                    meshRefinement::debugType(debug),
                    meshRefinement::writeType
                    (
                        meshRefinement::writeLevel()
                      | meshRefinement::WRITEMESH
                    ),
                    mesh.time().path()/meshRefiner_.timeName()
                );
                Info<< "Writing displacement field ..." << endl;
                meshMover.displacement().write();
                tmp<pointScalarField> magDisp
                (
                    mag(meshMover.displacement())
                );
                magDisp().write();
            }

            // Use current mesh as base mesh
            meshMover.correct();



            // See if any faces need splitting
            label nTotalSplit = returnReduce(splitFaces.size(), sumOp<label>());
            if (nTotalSplit && doSplit)
            {
                // Filter out baffle faces from faceZones of type
                // internal/baffle

                labelList duplicateFace(getInternalOrBaffleDuplicateFace());

                {
                    labelList oldSplitFaces(splitFaces.xfer());
                    List<labelPair> oldSplits(splits.xfer());
                    forAll(oldSplitFaces, i)
                    {
                        if (duplicateFace[oldSplitFaces[i]] == -1)
                        {
                            splitFaces.append(oldSplitFaces[i]);
                            splits.append(oldSplits[i]);
                        }
                    }
                    nTotalSplit = returnReduce
                    (
                        splitFaces.size(),
                        sumOp<label>()
                    );
                }

                // Update mesh
                meshRefiner_.splitFacesUndo
                (
                    splitFaces,
                    splits,
                    motionDict,

                    duplicateFace,
                    internalBaffles
                );

                // Redo meshMover
                meshMoverPtr.clear();
                ppPtr.clear();

                // Update mesh mover
                ppPtr = meshRefinement::makePatch(mesh, adaptPatchIDs);
                meshMoverPtr.reset
                (
                    new motionSmoother
                    (
                        mesh,
                        ppPtr(),
                        adaptPatchIDs,
                        meshRefinement::makeDisplacementField
                        (
                            pointMesh::New(mesh),
                            adaptPatchIDs
                        ),
                        motionDict
                    )
                );

                // Update snapping distance
                snapDist = calcSnapDistance(meshRefiner_, snapParams, ppPtr());

                if (debug&meshRefinement::MESH)
                {
                    const_cast<Time&>(mesh.time())++;
                    Info<< "Writing split-faces mesh to time "
                        << meshRefiner_.timeName() << endl;
                    meshRefiner_.write
                    (
                        meshRefinement::debugType(debug),
                        meshRefinement::writeType
                        (
                            meshRefinement::writeLevel()
                          | meshRefinement::WRITEMESH
                        ),
                        mesh.time().path()/meshRefiner_.timeName()
                    );
                }
            }


            if (debug&meshRefinement::MESH)
            {
                forAll(internalBaffles, i)
                {
                    const labelPair& p = internalBaffles[i];
                    const point& fc0 = mesh.faceCentres()[p[0]];
                    const point& fc1 = mesh.faceCentres()[p[1]];

                    if (mag(fc0-fc1) > meshRefiner_.mergeDistance())
                    {
                        FatalErrorInFunction
                            << "Separated baffles : f0:" << p[0]
                            << " centre:" << fc0
                            << " f1:" << p[1] << " centre:" << fc1
                            << " distance:" << mag(fc0-fc1)
                            << exit(FatalError);
                    }
                }
            }
        } // for nPreFeatIter
    }

    if (writeSnapVTK)
    {
        simpleVTKWriter
        (
            pp.localFaces(),
            pp.localPoints()
         ).write("snappyFeatureEmesh.vtk");
    }

    if (!dryRun)
    {
        // Check that faces are correctly patched
        autoPtr<mapPolyMesh> map = surfaceToPatch(snapParams,refineParams);

        updateBaffles(map, internalBaffles);
        //update pp addressing
        resetPrimitivePatchAddressing(adaptPatchIDs,pp);
    }

    if
    (
        (
            controller_.algorithm() == meshControl::STANDARD
            || controller_.algorithm() == meshControl::SHELL
        )
        && controller_.topoChanges() && snapParams.mergeBoundaryFaces()
        && !dryRun
    )
    {
        autoPtr<indirectPrimitivePatch> ppFtrPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                featureIDs
            )
        );
        indirectPrimitivePatch& ppFtr = ppFtrPtr();

        mergePatchFacesUndo
        (
            motionDict,
            calcSnapDistance(meshRefiner_, snapParams, ppFtr),
            ppFtr,
            internalBaffles,
            false
        );

        //update pp addressing
        resetPrimitivePatchAddressing(adaptPatchIDs,pp);
    }

    featureLinePrep flPrep(meshRefiner_,pp,snapDict,regionSnap);
    featureLineSnapper flSnap
    (
        flPrep,
        snapDict,
        motionDict,
        adaptPatchIDs,
        pp,
        meshRefiner_,
        internalBaffles,
        (
            controller_.algorithm() == meshControl::EXTRUDE
         || controller_.algorithm() == meshControl::DUAL
        ),
        controller_.topoChanges()
    );

    // Snap onto feature lines directly before snapping anything else
    // (note: this might not be the best place for this!)
    if (snapParams.directFeatureSnapping() && !dryRun)
    {
        scalar startTime = mesh.time().elapsedCpuTime();

        refinementFeatures& rFeatures = meshRefiner_.features()();
        rFeatures.trim(mesh);

        manifoldFeatures& mFeatures = rFeatures.manFeatures();
        PtrList<Tuple2<labelList, label>>& featureMeshes = mFeatures.features
        (
            labelMax, //filter by size of feature
            refineParams.minFeatureLength() //filter by minimum feature length
        );

        flPrep.prepareFeatureLines(rFeatures, featureMeshes);

        scalar endTime = mesh.time().elapsedCpuTime();
        Info<< "Feature lines prepared in " << endTime - startTime
             << "s" << endl << endl;
        startTime = endTime;

        PackedList<1> isZonedFace(setZonedFaces(meshRefiner_));
        flSnap.constructShadowChains(isZonedFace);

        const pointMesh& pMesh = pointMesh::New(mesh);

        // The current mesh is the starting mesh to smooth from.
        motionSmoother meshMover
        (
            mesh,
            pp,
            adaptPatchIDs,
            meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs),
            motionDict
        );

        meshMover.setDisplacement(flSnap.pointDisp());
        scaleMesh(snapParams, nInitErrors, internalBaffles, meshMover);
        meshMover.correct();

        if
        (
            controller_.algorithm() == meshControl::DUAL
         || controller_.algorithm() == meshControl::EXTRUDE
        )
        {
            //Try splitting cells based on region boundaries
            List<labelList> shadowChains = flSnap.getShadowChains();
            splitCells
            (
                shadowChains,
                refineParams,
                snapParams,
                adaptPatchIDs,
                isZonedFace,
                pp
            );
        }

        //Clear up manifold features
        mFeatures.clear();

        endTime = mesh.time().elapsedCpuTime();
        Info<< "Feature shadow chains constructed in "
             << endTime - startTime << "s" << endl << endl;
    }

    if (writeSnapVTK)
    {
        simpleVTKWriter
        (
            pp.localFaces(),
            pp.localPoints()
         ).write("snappyFeatureDirect.vtk");
    }

    label nInternalBaffles = internalBaffles.size();
    reduce(nInternalBaffles, sumOp<label>());

    if
    (
        (
            controller_.algorithm() == meshControl::EXTRUDE
         || controller_.algorithm() == meshControl::DUAL
        )
        && !dryRun
    )
    {
        if (internalOrBaffleZone)
        {
            autoPtr<mapPolyMesh> map = meshRefiner_.mergeZoneBaffles
            (
                true,  // internal zones
                true   // baffle zones
            );
            updateAfterMeshMods(map(),mesh,adaptPatchIDs,pp,flSnap);
            internalBaffles.clear();
        }

        // layer addition parameters
        const dictionary& layerDict = meshDict_.subDict("addLayersControls");
        if
        (
            controller_.algorithm() == meshControl::EXTRUDE
            && !snapParams.dualLayerRemoval()
        )
        {
            removeLayersExtrude
            (
                motionDict,
                snapParams,
                nInitErrors,
                layerDict,
                adaptPatchIDs,
                pp
            );
        }
        else
        {
            removeLayersDual
            (
                motionDict,
                snapParams,
                nInitErrors,
                layerDict,
                adaptPatchIDs,
                featureIDs,
                pp
            );
        }
    }

    if
    (
        (
            controller_.algorithm() == meshControl::EXTRUDE
            ||
            (
                controller_.algorithm() == meshControl::STANDARD
                && controller_.mode() != meshControl::FAST
            )
            || controller_.algorithm() == meshControl::SHELL
        )
        && snapParams.nSliverSmooths() > 0
    )
    {
        smoothSliverFaces
        (
            motionDict,
            snapParams,
            nInitErrors,
            adaptPatchIDs,
            internalBaffles,
            pp
         );

        if (writeSnapVTK)
        {
            simpleVTKWriter
            (
                pp.localFaces(),
                pp.localPoints()
             ).write("snappyFeatureSliver.vtk");
        }
    }

    if
    (
        controller_.algorithm() == meshControl::EXTRUDE
     && snapParams.nConformitySmooths() > 0
    )
    {
        smoothNonConformalFaces
        (
            motionDict,
            snapParams,
            nInitErrors,
            adaptPatchIDs,
            internalBaffles,
            pp
        );

        if (writeSnapVTK)
        {
            simpleVTKWriter
            (
                pp.localFaces(),
                pp.localPoints()
             ).write("snappyNonConformalSmooth.vtk");
        }
    }

    // Merge any introduced internalBaffles.
    if
    (
        (
            controller_.algorithm() == meshControl::STANDARD
         || controller_.algorithm() == meshControl::SHELL
        )
        && nInternalBaffles > 0
    )
    {
        autoPtr<mapPolyMesh> map = meshRefiner_.mergeZoneBaffles
        (
            true,   // internal zones
            false   // baffle zones
        );
        updateAfterMeshMods(map(),mesh,adaptPatchIDs,pp,flSnap);
    }

    // Merge coplanar boundary faces
    if
    (
        (
            controller_.algorithm() == meshControl::STANDARD
         || controller_.algorithm() == meshControl::SHELL
        )
     && controller_.topoChanges()
     && snapParams.mergeBoundaryFaces()
    )
    {
        mergePatchFacesUndo
        (
            motionDict,
            (
                snapParams.mergeAcrossPatches()
                ? List<labelList>(0)
                : flSnap.getMeshNumberedShadowChains()
            ),
            true
         );

        //update pp addressing
        resetPrimitivePatchAddressing(adaptPatchIDs,pp);
    }

    if (writeSnapVTK)
    {
        simpleVTKWriter
        (
            pp.localFaces(),
            pp.localPoints()
         ).write("snappyFeatureSecondMerge.vtk");
    }

    //Look for degenerate elements and try face splitting (no reversal yet)
    if
    (
        (
            controller_.algorithm() == meshControl::STANDARD
         || controller_.algorithm() == meshControl::SHELL
        )
     && snapParams.splitDegenerateCells()
     && controller_.topoChanges() && controller_.mode() != meshControl::FAST
     && !dryRun
    )
    {
        splitDegenerateCells(adaptPatchIDs, motionDict, pp);
    }

    if (writeSnapVTK)
    {
        simpleVTKWriter
        (
            pp.localFaces(),
            pp.localPoints()
         ).write("snappyFeatureSplitDegenerate.vtk");
    }

    if
    (
        (
            controller_.algorithm() == meshControl::STANDARD
            || controller_.algorithm() == meshControl::SHELL
        )
        && snapParams.collapseTol() > SMALL
        && controller_.topoChanges() && controller_.mode() != meshControl::FAST
        && !dryRun
    )
    {
        collapseSmallEdges
        (
            motionDict,
            adaptPatchIDs,
            snapParams,
            pp
        );
    }

    if
    (
        (
            controller_.algorithm() == meshControl::STANDARD
            || controller_.algorithm() == meshControl::SHELL
        )
        && snapParams.concaveTol() >= 0
        && controller_.mode() != meshControl::FAST
        && !dryRun
    )
    {
        correctConcaveFaces
        (
            motionDict,
            adaptPatchIDs,
            snapParams,
            pp
        );
    }

    if (writeSnapVTK)
    {
        simpleVTKWriter
        (
            pp.localFaces(),
            pp.localPoints()
         ).write("snappyFeatureCollapseEdges.vtk");
    }

    // Final check that faces are correctly patched
    if
    (
        controller_.algorithm() == meshControl::STANDARD
        || controller_.algorithm() == meshControl::SHELL
    )
    {
        autoPtr<mapPolyMesh> map = surfaceToPatch(snapParams,refineParams);
        //update pp addressing
        resetPrimitivePatchAddressing(adaptPatchIDs, pp);
    }

    if (snapParams.checkConformity())
    {
        checkSurfaceConformity(snapParams, pp);
    }

    //Merge faces created by degenerate cell splitting
    if
    (
        (
            controller_.algorithm() == meshControl::STANDARD
            || controller_.algorithm() == meshControl::SHELL
        )
        && controller_.topoChanges()
        && snapParams.mergeBoundaryFaces()
        && !dryRun
    )
    {
        DynamicList<labelList> edgesToKeep(pp.edges().size());
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        forAll(pp.edges(), edgeI)
        {
            labelList eFaces = pp.edgeFaces()[edgeI];
            label meshFaceI = pp.addressing()[eFaces[0]];
            label patch0 = patches.whichPatch(meshFaceI);
            label cell0 = mesh.faceOwner()[meshFaceI];
            bool noMerge = false;
            if (patch0 != -1)
            {
                if (eFaces.size() == 1)
                {
                    noMerge =  true;
                }
                else
                {
                    forAll(eFaces, eFI)
                    {
                        label patchI =
                            patches.whichPatch(pp.addressing()[eFaces[eFI]]);
                        label cellI =
                            mesh.faceOwner()[pp.addressing()[eFaces[eFI]]];
                        if (patchI != patch0 && cellI == cell0)
                        {
                            noMerge =  true;
                        }
                    }
                }
            }

            if (noMerge)
            {
                labelList e(2);
                e[0] = pp.meshPoints()[pp.edges()[edgeI][0]];
                e[1] = pp.meshPoints()[pp.edges()[edgeI][1]];
                edgesToKeep.append(e);
            }
        }
        edgesToKeep.shrink();

        // Merge coplanar boundary faces
        mergePatchFacesUndo
        (
            motionDict,
            edgesToKeep,
            true
        );
    }

    if (controller_.algorithm() == meshControl::EXTRUDE)
    {
        correctConvexLayerCells();
    }

    if
    (
        (
            controller_.algorithm() == meshControl::STANDARD
            || controller_.algorithm() == meshControl::SHELL
        )
        && controller_.topoChanges()
        && controller_.mode() != meshControl::FAST
        && !dryRun
    )
    {
        cleanupMesh(motionDict);
    }

    if
    (
        controller_.algorithm() == meshControl::EXTRUDE
        && meshDict_.found("addLayerRefinements")
    )
    {
        refineLayerCells();
    }

} // doSnap


void Foam::snappySnapDriver::refineLayerCells()
{
    Info<<"Refining layer cells "<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    const dictionary& layerRefDict = meshDict_.subDict("addLayerRefinements");
    const dictionary& refPatchDict = layerRefDict.subDict("patches");
    const dictionary& layerDict = meshDict_.subDict("addLayersControls");
    layerParameters layerParams(layerDict, patches, true);
    const labelList& numLayers = layerParams.numLayers();
    DynamicList<label> layerPatchIDs(patches.size());

    List<scalar> patchConcaveAngle(patches.size(), scalar(-1));
    List<scalar> patchConvexAngle(patches.size(), scalar(-1));
    List<label> patchLevel(patches.size(), label(-1));
    forAll(patches, patchi)
    {
        const word& patchName = patches[patchi].name();
        if
        (
            !patches[patchi].coupled()
            && numLayers[patchi] > -1
            && refPatchDict.found(patchName)
        )
        {
            const dictionary& refDict = refPatchDict.subDict(patchName);
            scalar concaveAngle = -Foam::cos
            (
                degToRad
                (
                    refDict.lookupOrDefault<scalar>
                    (
                        "concaveAngle", scalar(0.)
                    )
                )
            );
            scalar convexAngle = -Foam::cos
            (
                degToRad
                (
                    refDict.lookupOrDefault<scalar>
                    (
                        "convexAngle", scalar(0.)
                    )
                )
            );
            label level = refDict.lookupOrDefault<label>("level", label(-1));
            patchConcaveAngle[patchi] = concaveAngle;
            patchConvexAngle[patchi] = convexAngle;
            patchLevel[patchi] = level;
            layerPatchIDs.append(patchi);
        }
    }
    layerPatchIDs.shrink();

    while (true)
    {
        autoPtr<indirectPrimitivePatch> ppLayPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                layerPatchIDs
             )
         );
        indirectPrimitivePatch& ppLay = ppLayPtr();
        const labelList meshEdges
        (
            ppLay.meshEdges(mesh.edges(), mesh.pointEdges())
        );

        List<scalar> ppEdgeConcaveAngle(meshEdges.size(), scalar(-1));
        List<scalar> ppEdgeConvexAngle(meshEdges.size(), scalar(-1));
        List<label> ppEdgeMaxLevel(meshEdges.size(), -1);
        forAll(ppLay, i)
        {
            label facei = ppLay.addressing()[i];
            label patchi = patches.whichPatch(facei);
            const labelList& fEdges = ppLay.faceEdges()[i];
            scalar pConcave = patchConcaveAngle[patchi];
            scalar pConvex = patchConvexAngle[patchi];
            label pLevel = patchLevel[patchi];

            forAll(fEdges, fei)
            {
                label edgei = fEdges[fei];
                ppEdgeConcaveAngle[edgei] =
                    max(ppEdgeConcaveAngle[edgei],pConcave);
                ppEdgeConvexAngle[edgei] =
                    max(ppEdgeConcaveAngle[edgei],pConvex);
                ppEdgeMaxLevel[edgei] = max(ppEdgeConcaveAngle[edgei],pLevel);
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            meshEdges,
            ppEdgeConcaveAngle,
            maxEqOp<scalar>(),
            scalar(-1)
         );

        syncTools::syncEdgeList
        (
            mesh,
            meshEdges,
            ppEdgeConvexAngle,
            maxEqOp<scalar>(),
            scalar(-1)
         );

        syncTools::syncEdgeList
        (
            mesh,
            meshEdges,
            ppEdgeMaxLevel,
            maxEqOp<label>(),
            label(-1)
         );

        boolList exclFaces(ppLay.size(), false);
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            ppLay,
            meshEdges,
            exclFaces,
            0.99939,//convex
            0.99939//concave
        );
        List<Tuple2<edgeClassification::edgeType,scalar>>
            eType = eClass.edgeTypes();
        DynamicList<label> ftrFaces(ppLay.size());
        const volScalarField& layerCells =
            mesh.lookupObject<volScalarField>("layerStacks");
        boolList refineCell(mesh.nCells(), false);

        forAll(ppLay, i)
        {
            const labelList& fEdges = ppLay.faceEdges()[i];
            forAll(fEdges, fEI)
            {
                label edgei = fEdges[fEI];
                scalar eAngle = eType[edgei].second();
                if
                (
                    (
                        eType[edgei].first() == edgeClassification::CONCAVE
                        && eAngle < ppEdgeConcaveAngle[edgei]
                     )
                    ||
                    (
                        eType[edgei].first() == edgeClassification::CONVEX
                        && eAngle < ppEdgeConvexAngle[edgei]
                     )
                 )
                {
                    label maxEdgeLevel = ppEdgeMaxLevel[edgei];
                    label meshedgei = meshEdges[edgei];
                    const edge& e = mesh.edges()[meshedgei];
                    forAll(e,ei)
                    {
                        const labelList& pCells = mesh.pointCells()[e[ei]];
                        forAll(pCells, pCI)
                        {
                            label celli = pCells[pCI];
                            label cLevel = cellLevel[celli];
                            if (layerCells[celli] > -1 && cLevel < maxEdgeLevel)
                            {
                                refineCell[celli] = true;
                            }
                        }
                    }
                }
            }
        }

        labelList neiCellLevel(mesh.nFaces()-mesh.nInternalFaces(),-1);
        boolList neiRefineCell(mesh.nFaces()-mesh.nInternalFaces(),false);
        while (true)
        {
            neiCellLevel = -1;
            neiRefineCell = false;
            forAll(patches, patchI)
            {
                const polyPatch& patch = patches[patchI];
                const labelUList& faceCells = patch.faceCells();
                label bFaceI = patch.start()-mesh.nInternalFaces();

                if (patch.coupled())
                {
                    forAll(faceCells, i)
                    {
                        neiCellLevel[bFaceI] = cellLevel[faceCells[i]];
                        neiRefineCell[bFaceI] = refineCell[faceCells[i]];
                        bFaceI++;
                    }
                }
            }
            syncTools::swapBoundaryFaceList(mesh, neiCellLevel);
            syncTools::swapBoundaryFaceList(mesh, neiRefineCell);

            label nRefined = 0;
            forAll(mesh.cells(), celli)
            {
                if (!refineCell[celli] && layerCells[celli] > -1)
                {
                    const label cLevel = cellLevel[celli];
                    const cell& c = mesh.cells()[celli];

                    forAll(c,cFI)
                    {
                        label facei = c[cFI];
                        label neiLevel = -1;
                        bool neiRefined = false;
                        label patchi = patches.whichPatch(facei);
                        if (patchi == -1)
                        {
                            label own = mesh.faceOwner()[facei];
                            label nei
                            (
                                own == celli
                                ? mesh.faceNeighbour()[facei]
                                : own
                             );
                            neiLevel = cellLevel[nei];
                            neiRefined = refineCell[nei];
                        }
                        else  if (patches[patchi].coupled())
                        {
                            label bfacei = facei-mesh.nInternalFaces();
                            neiLevel = neiCellLevel[bfacei];
                            neiRefined = neiRefineCell[bfacei];
                        }

                        if (neiRefined && !refineCell[celli])
                        {
                            if (neiLevel+1-cLevel > 1)
                            {
                                refineCell[celli] = true;
                                nRefined++;
                                break;
                            }
                        }
                    }
                }
            }
            if (returnReduce(nRefined, sumOp<label>()) == 0)
            {
                break;
            }
        }

        label nRefinedCells = 0;
        DynamicList<label> refineFaces(ppLay.size());
        forAll(ppLay, i)
        {
            label facei = ppLay.addressing()[i];
            label own = mesh.faceOwner()[facei];

            if (refineCell[own])
            {
                nRefinedCells++;
                refineFaces.append(facei);
            }
        }

        reduce(nRefinedCells, sumOp<label>());
        if (nRefinedCells > 0)
        {
            Info<<"Refining " << nRefinedCells
                <<" boundary cells"<<endl;
            autoPtr<indirectPrimitivePatch> ppRefPtr
            (
                new indirectPrimitivePatch
                (
                    IndirectList<face>(mesh.faces(), refineFaces),
                    mesh.points()
                 )
             );
            indirectPrimitivePatch& ppRef = ppRefPtr();
            boundaryLayerRefinement bRef(ppRef,meshRefiner_);
            bRef.setRefinement();
        }
        else
        {
            break;
        }
    }

    //Smooth surfaces
    Switch smoothSurf = layerRefDict.lookupOrDefault<Switch>
    (
       "smoothSurface",
       true
    );
    if (smoothSurf)
    {
        const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
        const labelList pointLevel =
            meshRefiner_.meshCutter().pointLevel();

        autoPtr<indirectPrimitivePatch> ppLayPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                layerPatchIDs
             )
         );
        indirectPrimitivePatch& ppLay = ppLayPtr();
        const labelList meshEdges
        (
            ppLay.meshEdges(mesh.edges(), mesh.pointEdges())
        );

        boolList exclFaces(ppLay.size(), false);
        edgeClassification eClass
        (
            mesh,
            mesh.points(),
            ppLay,
            meshEdges,
            exclFaces,
            0.93969,//convex
            0.93969//concave
        );
        List<Tuple2<edgeClassification::edgeType,scalar>>
            eType = eClass.edgeTypes();

        boolList stationaryPoints(ppLay.localPoints().size(), false);
        forAll(eType, edgei)
        {
            if
            (
                eType[edgei].first() != edgeClassification::MANIFOLD
            )
            {
                const edge& e = ppLay.edges()[edgei];
                stationaryPoints[e[0]] = true;
                stationaryPoints[e[1]] = true;
            }
        }

        syncTools::syncPointList
        (
            mesh,
            ppLay.meshPoints(),
            stationaryPoints,
            orEqOp<bool>(),
            false
        );

        scalar weight = 0.2;
        pointField newPoints =  mesh.points();
        vectorField pointDisp(mesh.nPoints(), vector::zero);
        scalarField sumAreas(ppLay.nPoints(), 0.0);
        vectorField sumCentres(ppLay.nPoints(), vector::zero);

        const refinementSurfaces& surfaces = meshRefiner_.surfaces();
        labelList snapSurfaces =
            surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
            (
                surfaces.surfZones()
            );

        for (label i = 0; i < 5; i++)
        {
            pointField ppFaceCentres(ppLay.size());
            pointField ppFaceAreas(ppLay.size());

            forAll(ppLay, i)
            {
                const label& faceI = ppLay.addressing()[i];
                ppFaceCentres[i] =
                    mesh.faces()[faceI].centre(newPoints);
                ppFaceAreas[i] =
                    mesh.faces()[faceI].areaNormal(newPoints);
            }

            sumCentres = vector::zero;
            sumAreas = 0.0;

            forAll(ppLay.meshPoints(), pointI)
            {
                if (!stationaryPoints[pointI])
                {
                    const labelList& pFaces = ppLay.pointFaces()[pointI];
                    forAll(pFaces, faceI)
                    {
                        const label pFaceI = pFaces[faceI];
                        const label meshFaceI = ppLay.addressing()[pFaceI];
                        const label own = mesh.faceOwner()[meshFaceI];
                        label  level = cellLevel[own];
                        scalar len = edge0Len / pow(2., level);
                        scalar weight = sqrt(mag(ppFaceAreas[pFaceI]))
                            / len;
                        sumCentres[pointI] +=  ppFaceCentres[pFaceI]
                            * weight;
                        sumAreas[pointI] += weight;
                    }
                }
            }
            pointDisp = vector::zero;

            syncTools::syncPointList
            (
                mesh,
                ppLay.meshPoints(),
                sumCentres,
                plusEqOp<vector>(),
                vector::zero        // null value
             );

            syncTools::syncPointList
            (
                mesh,
                ppLay.meshPoints(),
                sumAreas,
                plusEqOp<scalar>(),
                scalar(0.0)         // null value
             );

            pointField snapPts(ppLay.meshPoints().size());
            scalarField snapDist(ppLay.meshPoints().size());
            label nSet = 0;
            forAll(ppLay.meshPoints(), i)
            {
                label meshPointI = ppLay.meshPoints()[i];
                if (!stationaryPoints[i] && sumAreas[i] > SMALL)
                {
                    snapPts[nSet] = sumCentres[i] / (sumAreas[i] + SMALL);
                    snapDist[nSet] =
                       sqr(edge0Len / (1<<pointLevel[meshPointI]));
                    nSet++;
                }
            }
            snapPts.setSize(nSet);
            snapDist.setSize(nSet);

            List<pointIndexHit> hitInfo;
            labelList hitSurface;
            surfaces.findNearest
            (
                snapSurfaces,
                snapPts,
                snapDist,        // sqr of attract distance
                hitSurface,
                hitInfo
            );

            nSet = 0;
            forAll(ppLay.meshPoints(), i)
            {
                label meshPointI = ppLay.meshPoints()[i];
                if (!stationaryPoints[i] && sumAreas[i] > SMALL)
                {
                    if (hitInfo[nSet].hit())
                    {
                        pointDisp[meshPointI] = hitInfo[nSet].hitPoint()
                           - newPoints[meshPointI];
                    }
                    else
                    {
                        pointDisp[meshPointI] =
                           sumCentres[i]/(sumAreas[i] + SMALL)
                           - newPoints[meshPointI];
                    }
                    nSet++;
                }
            }

            syncTools::syncPointList
            (
                mesh,
                pointDisp,
                maxMagSqrEqOp<point>(),
                vector::zero
             );

            forAll(ppLay.meshPoints(), i)
            {
                if (!stationaryPoints[i])
                {
                    label meshPointI = ppLay.meshPoints()[i];
                    newPoints[meshPointI] = newPoints[meshPointI]
                        + weight*pointDisp[meshPointI];
                }
            }
        }

        mesh.movePoints(newPoints);
    }

    return;
}


void Foam::snappySnapDriver::checkSurfaceConformity
(
    const snapParameters& snapParams,
    const indirectPrimitivePatch& pp
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const refinementSurfaces& surfaces = meshRefiner_.surfaces();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();

    List<pointIndexHit> hitInfo;
    labelList hitSurface;

    scalarField maxSnapDist =
        calcSnapDistance(meshRefiner_, snapParams, pp);
    labelList snapSurfaces = identity(surfaces.surfaces().size());
    surfaces.findNearest
    (
        snapSurfaces,
        pp.faceCentres(),
        sqr(maxSnapDist), // sqr of attract
        hitSurface,
        hitInfo
    );

    scalarField sDist(pp.size(), GREAT);
    scalarField rDist(pp.size(), GREAT);

    label nHits = 0;
    scalar sumDisp = scalar(0);
    scalarField maxDist(Pstream::nProcs(),-GREAT);
    vectorField maxLoc(Pstream::nProcs(),vector::zero);
    const label proci = Pstream::myProcNo();
    forAll(pp, i)
    {
        if (hitInfo[i].hit())
        {
            nHits++;
            point fC = pp.faceCentres()[i];
            scalar sd = mag(fC-hitInfo[i].hitPoint());
            if (sd > maxDist[proci])
            {
                maxDist[proci] = sd;
                maxLoc[proci] = fC;
            }
            sumDisp += sd;

            const label meshFaceI = pp.addressing()[i];
            const label own = mesh.faceOwner()[meshFaceI];
            label  level = cellLevel[own];
            scalar len = edge0Len / pow(2., level);
            scalar rd = (sd/len);

            sDist[i] = sd;
            rDist[i] = rd;
        }
    }

    Pstream::allGatherList(maxDist);
    Pstream::allGatherList(maxLoc);
    label minProci = findMax(maxDist);

    reduce(
        std::tie(nHits, sumDisp),
        ParallelOp<sumOp<label>, sumOp<scalar>>{}
    );

    if (nHits > 0)
    {
        sumDisp /= nHits;
        Info<<"Max snap error : "
            << maxDist[minProci] << " at location : "<< maxLoc[minProci] <<endl;
        Info<<"Average snap error : "
            << sumDisp << nl <<endl;

        simpleVTKWriter conformVTK
        (
            pp.localFaces(),
            pp.localPoints()
        );
        conformVTK.addFaceData("relative", rDist);
        conformVTK.addFaceData("absolute", sDist);
        conformVTK.write("snapConformity.vtk");
    }
}


void Foam::snappySnapDriver::correctConvexLayerCells()
{
    fvMesh& mesh = meshRefiner_.mesh();

    const volScalarField& layerCells =
        mesh.lookupObject<volScalarField>("layerStacks");
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    boolList boundaryPoints(mesh.nPoints(), false);
    boolList boundaryEdges(mesh.nEdges(), false);
    boolList boundaryFaces(mesh.nFaces(), false);
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (!pp.coupled())
        {
            label startFaceI = pp.start();
            forAll(pp, i)
            {
                label faceI = startFaceI+i;

                boundaryFaces[faceI] = true;

                const face& f = mesh.faces()[faceI];
                forAll(f,fp)
                {
                    boundaryPoints[f[fp]] = true;
                }
                const labelList& fEdges = mesh.faceEdges()[faceI];
                forAll(fEdges, fEI)
                {
                    boundaryEdges[fEdges[fEI]] = true;
                }
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        boundaryPoints,
        orEqOp<bool>(),
        false
    );
    syncTools::syncEdgeList
    (
        mesh,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    vectorField newPoints(mesh.nPoints(), vector::zero);
    scalarField dispWeight(mesh.nPoints(), 0);
    label nMoved = 0;
    forAll(mesh.cells(), celli)
    {
        if (layerCells[celli] > -1)
        {
            const cell& c = mesh.cells()[celli];
            const labelList& cEdges = mesh.cellEdges()[celli];
            labelHashSet cFaces(c);
            forAll(cEdges, cEI)
            {
                label edgei = cEdges[cEI];
                if (boundaryEdges[edgei])
                {
                    const edge& e = mesh.edges()[edgei];
                    const labelList& eFaces = mesh.edgeFaces()[edgei];
                    DynamicList<label> nbrFaces(2);

                    forAll(eFaces,eFI)
                    {
                        label facei = eFaces[eFI];
                        if (cFaces.found(facei))
                        {
                            nbrFaces.append(facei);
                        }
                    }
                    label boundFace = -1;
                    label intFace = -1;

                    if (nbrFaces.size() == 2)
                    {
                        if
                        (
                            boundaryFaces[nbrFaces[0]]
                            && !boundaryFaces[nbrFaces[1]]
                        )
                        {
                            boundFace = nbrFaces[0];
                            intFace  = nbrFaces[1];
                        }
                        else if
                        (
                            boundaryFaces[nbrFaces[1]]
                            && !boundaryFaces[nbrFaces[0]]
                        )
                        {
                            boundFace = nbrFaces[1];
                            intFace  = nbrFaces[0];
                        }
                    }

                    if (boundFace != -1 && intFace != -1)
                    {
                        vector bN = mesh.faceAreas()[boundFace];
                        bN /= (mag(bN) + SMALL);
                        vector iN = mesh.faceAreas()[intFace];
                        iN /= (mag(iN) + SMALL);
                        label intPatch = patches.whichPatch(intFace);

                        if
                        (
                            intPatch == -1
                            && mesh.faceOwner()[intFace] != celli
                        )
                        {
                            iN  = -iN;
                        }
                        point fcToEc = e.centre(mesh.points()) -
                                mesh.faceCentres()[intFace];
                        point aveNorm = 0.5*(bN+iN);
                        if ((fcToEc&aveNorm) < 0)
                        {
                            point bFC = mesh.faceCentres()[boundFace];
                            if (mag(bN) < VSMALL)
                            {
                                continue;
                            }
                            plane bPlane(bFC, bN);
                            face f = mesh.faces()[intFace];
                            forAll(f,fp)
                            {
                                label pointi = f[fp];
                                if (!boundaryPoints[pointi])
                                {
                                    point pt = mesh.points()[pointi];
                                    vector dV = bPlane.nearestPoint(pt) -pt;
                                    newPoints[pointi] += dV;
                                    dispWeight[pointi] += scalar(1);
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
        dispWeight,
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
        if (dispWeight[pointi] > 0)
        {
            newPoints[pointi] /= dispWeight[pointi];
        }
    }

    newPoints += mesh.points();
    pointField newCellCentres(mesh.nCells(), vector::zero);
    scalarField newCellVolumes(mesh.nCells(), 0);
    pointField newFaceCentres(mesh.nFaces(), vector::zero);
    vectorField newFaceAreas(mesh.nFaces(), vector::zero);
    scalarField newMagFaceAreas(mesh.nFaces(), Zero);

    mesh.makeFaceCentresAndAreas
    (
        newPoints,
        newFaceCentres,
        newFaceAreas,
        newMagFaceAreas
    );
    mesh.makeCellCentresAndVols
    (
        newFaceCentres,
        newFaceAreas,
        newCellCentres,
        newCellVolumes
    );

    faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+1);
    labelList checkFaces(meshRefiner_.selectInnerFaces());

    Foam::polyMeshGeometry::checkFacePyramids
    (
        false,
        -SMALL,
        mesh,
        newFaceCentres,
        newCellCentres,
        newFaceAreas,
        newPoints,
        checkFaces,
        List<labelPair>(0),
        &errorFaces
    );

    Foam::polyMeshGeometry::checkFaceDotProduct
    (
        false,
        scalar(70),
        scalar(180),
        mesh,
        newFaceCentres,
        newCellCentres,
        newFaceAreas,
        checkFaces,
        List<labelPair>(0),
        &errorFaces
    );

    boolList errorPts(mesh.nPoints(), false);
    forAllConstIter(labelHashSet, errorFaces, iter)
    {
        label facei = iter.key();
        label own = mesh.faceOwner()[facei];
        const labelList& ownPts = mesh.cellPoints()[own];
        forAll(ownPts, oPI)
        {
            errorPts[ownPts[oPI]] = true;
        }

        if (mesh.isInternalFace(facei))
        {
            label nei = mesh.faceNeighbour()[facei];
            const labelList& neiPts = mesh.cellPoints()[nei];
            forAll(neiPts, nPI)
            {
                errorPts[neiPts[nPI]] = true;
            }
        }
    }
    syncTools::syncPointList
    (
        mesh,
        errorPts,
        orEqOp<bool>(),
        false
    );

    forAll(errorPts, pointi)
    {
        if (errorPts[pointi])
        {
            newPoints[pointi] = mesh.points()[pointi];
        }
        else if (dispWeight[pointi] > 0)
        {
            nMoved++;
        }
    }

    Info<<"Moving convex layer cell internal points : "
        << returnReduce(nMoved, sumOp<label>()) << nl <<endl;

    mesh.movePoints(newPoints);
}


// Clean up final mesh
void Foam::snappySnapDriver::cleanupMesh(const dictionary& motionDict)
{
    Info<< nl
        << "Cleaning up final snapped mesh" << nl
        << "------------------------------" << nl
        << endl;

    fvMesh& mesh = meshRefiner_.mesh();

    //Remove hanging nodes (connected to only two edges)
    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));

    labelList nPointEdges(mesh.nPoints(), 0);

    //TODO loop only through innerGrid points //GGG
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

    pointField oldPoints = mesh.points();
    pointField newPoints(mesh.nPoints(), vector::zero);
    //TODO loop only through innerGrid points //GGG
    forAll(mesh.points(), pointi)
    {
        if (nPointEdges[pointi] == 2)
        {
            const labelList& pEdges = mesh.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if (isMasterEdge[edgei])
                {
                    edge e = mesh.edges()[edgei];
                    label otherPt(e[0] == pointi ? e[1] : e[0]);
                    newPoints[pointi] += mesh.points()[otherPt];
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        newPoints,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    forAll(mesh.points(), pointi)
    {
        if (nPointEdges[pointi] == 2)
        {
            newPoints[pointi] /= scalar(2);
        }
        else
        {
            newPoints[pointi] = mesh.points()[pointi];
        }
    }

    //Move points and check if errors introduced
    mesh.movePoints(newPoints);
    labelHashSet wrongFaces(mesh.nFaces()/100);
    //TODO check only innerGrid faces //GGG
    motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces);

    //TODO loop only through innerGrid points //GGG
    forAll(mesh.points(), pointi)
    {
        if (nPointEdges[pointi] == 2)
        {
            const labelList& pointFaces = mesh.pointFaces()[pointi];
            bool reset = false;
            forAll(pointFaces, facei)
            {
                if (wrongFaces.found(pointFaces[facei]))
                {
                    //Mark as valid to prevent removal
                    reset = true;
                    break;
                }
            }
            if (!reset)
            {
                const labelList& pointCells = mesh.pointCells()[pointi];
                forAll(pointCells, i)
                {
                    label celli = pointCells[i];
                    const cell& c = mesh.cells()[celli];
                    forAll(c, cfi)
                    {
                        if (wrongFaces.found(c[cfi]))
                        {
                            reset = true;
                            break;
                        }
                    }
                    if (reset)
                    {
                        break;
                    }
                }
            }
            if (reset)
            {
                nPointEdges[pointi] = 3;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nPointEdges,
        maxEqOp<label>(), // combine op
        label(0)     // null value
    );

    //move points back to original position
    mesh.movePoints(oldPoints);

    //TODO loop only through innerGrid faces //GGG
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

        bool preventRemoval = false;
        if (keepPts.size() != f.size())
        {
            if (keepPts.size() >= 3)
            {
                face nf(keepPts);
                scalar fA = nf.mag(mesh.points());
                if (fA < SMALL)
                {
                    preventRemoval=true;
                }
            }

            if (preventRemoval)
            {
                forAll(f, fp)
                {
                    if (nPointEdges[f[fp]] <= 2)
                    {
                        //Mark as valid to prevent removal
                        nPointEdges[f[fp]] = 3;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh,
        nPointEdges,
        maxEqOp<label>(), // combine op
        label(0)     // null value
    );

    polyTopoChange meshMod(mesh);

    //TODO loop only through innerGrid faces //GGG
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
    //TODO loop only through innerGrid points //GGG
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

    Info<<"Removing "<<nPtsRemoved<<" hanging nodes"<<endl;

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
}


void Foam::snappySnapDriver::splitDegenerateCells
(
    const labelList& adaptPatchIDs,
    const dictionary& motionDict,
    indirectPrimitivePatch &pp
)
{
    Info<< nl
        << "Splitting degenerate cell faces" << nl
        << "-------------------------------" << nl
        << endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    //decompose to keep degenerate faces on single core
    if (Pstream::parRun())
    {
        boolList blockedFaces(mesh.nFaces(), true);
        label nProcBlocked = 0;

        const faceZoneMesh& fZones = mesh.faceZones();
        forAll(fZones, fZI)
        {
            const faceZone& fZone = fZones[fZI];
            forAll(fZone, i)
            {
                label faceI = fZone[i];
                blockedFaces[faceI] = false;
            }
        }

        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];
            label ownerI = mesh.faceOwner()[meshFaceI];

            const face& f = mesh.faces()[meshFaceI];

            labelHashSet edgeFaces(pp[i].size());

            label prevFp = f[0];
            forAll(f, fp)
            {
                label nextFp = f[f.fcIndex(fp)];
                label meshEdgeI = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[prevFp],
                    prevFp,
                    nextFp
                );

                const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];

                forAll(eFaces, efI)
                {
                    label eFace = eFaces[efI];
                    label patchI = patches.whichPatch(eFace);

                    if
                    (
                        mesh.isInternalFace(eFace)
                        || (patchI != -1
                            && isA<processorPolyPatch>(patches[patchI]))
                    )
                    {
                        label own = mesh.faceOwner()[eFace];
                        label nei =  -1;

                        if (mesh.isInternalFace(eFace))
                        {
                            nei = mesh.faceNeighbour()[eFace];
                        }

                        if
                        (
                            (eFace != meshFaceI)
                            && ((own == ownerI) || (nei == ownerI))
                         )
                        {
                            if (!edgeFaces.insert(eFace))
                            {
                                blockedFaces[eFace] = false;
                                if
                                (
                                    patchI != -1
                                    && isA<processorPolyPatch>(patches[patchI])
                                )
                                {
                                    nProcBlocked++;
                                }
                            }
                        }
                    }
                    prevFp = nextFp;
                }
            }
        }

        if (returnReduce(nProcBlocked, sumOp<label>()) > 0)
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

            //update pp addressing
            resetPrimitivePatchAddressing(adaptPatchIDs,pp);
        }
    }

    PackedList<1> isZonedFace(setZonedFaces(meshRefiner_));

    boolList setFaces(mesh.nFaces(), false);
    bool firstPass =  true;
    label finalCheck = 0;

    while (true)
    {
        polyTopoChange meshMod(mesh);
        labelHashSet splitFaces(pp.size());
        label nChanged = 0;
        Map<label> changeMap(pp.size());
        DynamicList<List<Tuple2<label,face>>> originalFaces(pp.size());

        forAll(pp, i)
        {
            label meshFaceI = pp.addressing()[i];

            if (isZonedFace.get(meshFaceI))
            {
                continue;
            }
            label ownerI = mesh.faceOwner()[meshFaceI];

            const face& f = mesh.faces()[meshFaceI];

            labelHashSet edgeFaces(pp[i].size());
            label prevFp = f[0];


            label startFp = -1;
            label endFp = -1;
            label nDegenerateEdges = 0;
            label degenerateFace = -1;
            label neighbourI = -1;
            boolList dEdges(f.size(), false);

            forAll(f, fp)
            {
                label nextFp = f[f.fcIndex(fp)];
                label meshEdgeI = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[prevFp],
                    prevFp,
                    nextFp
                 );

                const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];

                if (eFaces.size() >= 3)
                {
                    forAll(eFaces, efI)
                    {
                        label eFace = eFaces[efI];

                        if (mesh.isInternalFace(eFace))
                        {
                            label own = mesh.faceOwner()[eFace];
                            label nei = mesh.faceNeighbour()[eFace];
                            if
                            (
                                (eFace != meshFaceI)
                                && ((own == ownerI) || (nei == ownerI))
                             )
                            {
                                if (!edgeFaces.insert(eFace) && !setFaces[eFace])
                                {
                                    if (nDegenerateEdges == 0)
                                    {
                                        degenerateFace = eFace;
                                        neighbourI
                                            = (own == ownerI ? nei : own);
                                    }
                                    if (eFace == degenerateFace)
                                    {
                                        nDegenerateEdges++;
                                    }
                                }
                            }
                        }
                    }
                }
                prevFp = nextFp;
            }

            if (nDegenerateEdges > 0)
            {
                nDegenerateEdges = 0;
                prevFp = f[0];
                forAll(f, fp)
                {
                    label nextFp = f[f.fcIndex(fp)];
                    label meshEdgeI = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[prevFp],
                        prevFp,
                        nextFp
                     );

                    const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];

                    forAll(eFaces, efI)
                    {
                        label eFace = eFaces[efI];

                        if (eFace == degenerateFace)
                        {
                            dEdges[fp] = true;
                            nDegenerateEdges++;
                            break;
                        }
                    }
                    prevFp = nextFp;
                }
                bool found = false;
                label fp = 0;
                while (true)
                {
                    if (!dEdges[fp])
                    {
                        found = true;
                    }
                    if (found && dEdges[fp])
                    {
                        if (startFp == -1)
                        {
                            startFp = f[fp];
                        }
                    }
                    else if (found && startFp != -1 && !dEdges[fp])
                    {
                        endFp = f[fp];
                        break;
                    }
                    fp = f.fcIndex(fp);
                }
            }

            if (nDegenerateEdges > 0)
            {
                setFaces[degenerateFace] = true;

                DynamicList<label> modPts(f.size());
                label fp = 0;
                bool found = false;
                while (true)
                {
                    label pointI = f[fp];
                    if (pointI == endFp || found)
                    {
                        found = true;
                        modPts.append(pointI);
                        if (pointI == startFp)
                        {
                            break;
                        }
                    }
                    fp = f.fcIndex(fp);
                }
                face modF(modPts.shrink());

                DynamicList<label> newPts(f.size());
                fp = 0;
                found = false;
                while (true)
                {
                    label pointI = f[fp];
                    if (pointI == startFp || found)
                    {
                        found = true;
                        newPts.append(pointI);
                        if (pointI == endFp)
                        {
                            break;
                        }
                    }
                    fp = f.fcIndex(fp);
                }
                face newF(newPts.shrink());

                face dFace = mesh.faces()[degenerateFace];
                DynamicList<label> intPts(dFace.size());
                fp = 0;
                found = false;
                bool reverse = false;
                while (true)
                {
                    label pointI = dFace[fp];
                    if (pointI == startFp || found)
                    {
                        found = true;
                        if
                        (
                            pointI == startFp
                            && f[f.fcIndex(f.which(pointI))] ==
                            dFace[dFace.fcIndex(fp)]
                         )
                        {
                            reverse = true;
                        }

                        intPts.append(pointI);
                        if (pointI == endFp)
                        {
                            break;
                        }
                    }
                    if (reverse)
                    {
                        fp = dFace.rcIndex(fp);
                    }
                    else
                    {
                        fp = dFace.fcIndex(fp);
                    }
                }
                face intF(intPts.shrink());

                if (reverse)
                {
                    intF = intF.reverseFace();
                }

                //look for convex face
                bool concave = true;
                vector newFNormal = newF.areaNormal(mesh.points());
                vector modFNormal = modF.areaNormal(mesh.points());

                if ((modFNormal & newFNormal) < 0)
                {
                    concave = false;
                }

                //look at face area change
                scalar ac1 =  modF.mag(mesh.points())
                    / (f.mag(mesh.points())+ SMALL);
                scalar ac2 =  intF.mag(mesh.points())
                    / (mesh.faces()[degenerateFace].mag(mesh.points())+ SMALL);

                if
                (
                    modF.size() > 2 && newF.size() > 2 && intF.size() > 2
                    && (ac1 < 0.80 && ac1 > 0.20 ) && (ac2 < 0.80 && ac2 > 0.20)
                    && splitFaces.insert(degenerateFace)
                    && concave
                 )
                {
                    label patchID = mesh.boundaryMesh().whichPatch(meshFaceI);
                    label zoneID = mesh.faceZones().whichZone(meshFaceI);
                    bool zoneFlip = false;

                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = mesh.faceZones()[zoneID];
                        zoneFlip = fZone.flipMap()[fZone.whichFace(meshFaceI)];
                    }

                    //modify boundary face
                    meshMod.modifyFace
                    (
                        modF, // modified face
                        meshFaceI,      // label of face being modified
                        ownerI,         // owner
                        -1,             // neighbour
                        false,          // face flip
                        patchID,        // patch for face
                        zoneID,         // zone for face
                        zoneFlip        // face flip in zone
                    );

                    //Add new boundary face
                    label newFaceI = meshMod.addFace
                    (
                        newF,           // face
                        neighbourI,     // owner
                        -1,             // neighbour
                        -1,             // master point
                        -1,             // master edge
                        meshFaceI,      // master face
                        false,          // flux flip
                        patchID,        // patch for face
                        zoneID,         // zone for face
                        zoneFlip        // face zone flip
                    );

                    zoneID = mesh.faceZones().whichZone(degenerateFace);
                    zoneFlip = false;

                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = mesh.faceZones()[zoneID];
                        zoneFlip = fZone.flipMap()
                            [fZone.whichFace(degenerateFace)];
                    }

                    //modify interior face
                    meshMod.modifyFace
                    (
                        intF,   // modified face
                        degenerateFace,        // label of face being modified
                        mesh.faceOwner()[degenerateFace],     // owner
                        mesh.faceNeighbour()[degenerateFace], // neighbour
                        false,                 // face flip
                        -1,                    // new patch for face
                        zoneID,                // zone for face
                        zoneFlip               // face flip in zone
                    );
                    changeMap.insert(meshFaceI, nChanged);
                    changeMap.insert(newFaceI, nChanged);
                    changeMap.insert(degenerateFace, nChanged);

                    List<Tuple2<label,face>> oldFaces(3);
                    oldFaces[0] = Tuple2<label, face>(meshFaceI, f);
                    oldFaces[1] = Tuple2<label, face>
                    (degenerateFace, mesh.faces()[degenerateFace]);
                    oldFaces[2] = Tuple2<label, face>(newFaceI, f);
                    originalFaces.append(oldFaces);

                    nChanged++;
                }
            }
        }
        originalFaces.shrink();

        if (returnReduce(nChanged, sumOp<label>()) == 0)
        {
            if (firstPass)
            {
                setFaces = false;
                firstPass = false;
            }
            else if (finalCheck == 0)
            {
                finalCheck++;
            }
            else
            {
                break;
            }
        }
        else
        {
            Info<<"Updating "<<returnReduce(nChanged, sumOp<label>())
                <<" Degenerate faces"<<endl;

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

            //update pp addressing
            resetPrimitivePatchAddressing(adaptPatchIDs,pp);

            boolList updateSetFaces(mesh.nFaces(), false);

            forAll(setFaces, faceI)
            {
                if (setFaces[faceI])
                {
                    updateSetFaces[map().reverseFaceMap()[faceI]] = true;
                }
            }
            setFaces = updateSetFaces;

            labelHashSet wrongFaces(mesh.nFaces()/100);
            labelList checkFaces(meshRefiner_.selectInnerFaces());
            label nErrors = motionSmoother::checkMesh(false, mesh, motionDict, checkFaces, wrongFaces);

            if (nErrors > 0)
            {
                polyTopoChange meshReset(mesh);
                boolList reset(nChanged, false);

                forAll(map().reverseFaceMap(), oldFaceI)
                {
                    label newFaceI = map().reverseFaceMap()[oldFaceI];

                    bool foundError =  false;
                    label own = mesh.faceOwner()[newFaceI];
                    const cell& ownFaces = mesh.cells()[own];
                    forAll(ownFaces, cfI)
                    {
                        if (wrongFaces.found(ownFaces[cfI]))
                        {
                            foundError = true;
                        }
                    }

                    if (newFaceI < mesh.nInternalFaces())
                    {
                        label nei = mesh.faceNeighbour()[newFaceI];
                        const cell& neiFaces = mesh.cells()[nei];
                        forAll(neiFaces, cfI)
                        {
                            if (wrongFaces.found(neiFaces[cfI]))
                            {
                                foundError = true;
                            }
                        }
                    }

                    if (foundError && changeMap.found(oldFaceI))
                    {
                        label changeID = changeMap[oldFaceI];
                        if (!reset[changeID])
                        {
                            const List<Tuple2<label,face>> originalFace =
                                originalFaces[changeID];
                            const label i1 = originalFace[0].first();
                            const face face1 = originalFace[0].second();
                            const label i2 = originalFace[1].first();
                            const face face2 = originalFace[1].second();
                            const label i3 = originalFace[2].first();

                            label ni1 = map().reverseFaceMap()[i1];
                            label ni2 = map().reverseFaceMap()[i2];
                            label ni3 = map().reverseFaceMap()[i3];

                            meshReset.removeFace(ni3, -1);

                            face nface1(face1.size());
                            face nface2(face2.size());

                            forAll(nface1, i)
                            {
                                nface1[i] = map().reversePointMap()[face1[i]];
                            }

                            forAll(nface2, i)
                            {
                                nface2[i] = map().reversePointMap()[face2[i]];
                            }

                            label patchID = mesh.boundaryMesh().whichPatch(ni1);
                            label zoneID = mesh.faceZones().whichZone(ni1);
                            bool zoneFlip = false;

                            if (zoneID >= 0)
                            {
                                const faceZone& fZone =
                                    mesh.faceZones()[zoneID];
                                zoneFlip =
                                    fZone.flipMap()[fZone.whichFace(ni1)];
                            }

                            //modify boundary face
                            meshReset.modifyFace
                            (
                                nface1,   // modified face
                                ni1,        // label of face being modified
                                mesh.faceOwner()[ni1],     // owner
                                -1, // neighbour
                                false,                 // face flip
                                patchID,               // new patch for face
                                zoneID,                // zone for face
                                zoneFlip               // face flip in zone
                             );

                            patchID = mesh.boundaryMesh().whichPatch(ni2);
                            zoneID = mesh.faceZones().whichZone(ni2);
                            zoneFlip = false;

                            if (zoneID >= 0)
                            {
                                const faceZone& fZone =
                                    mesh.faceZones()[zoneID];
                                zoneFlip =
                                    fZone.flipMap()[fZone.whichFace(ni2)];
                            }

                            //modify interior face
                            meshReset.modifyFace
                            (
                                nface2,   // modified face
                                ni2,        // label of face being modified
                                mesh.faceOwner()[ni2],     // owner
                                mesh.faceNeighbour()[ni2],     // owner
                                false,                 // face flip
                                patchID,               // new patch for face
                                zoneID,                // zone for face
                                zoneFlip               // face flip in zone
                             );

                            reset[changeID] = true;
                        }
                    }
                }
                // Change the mesh (no inflation)
                autoPtr<mapPolyMesh> resetMap =
                    meshReset.changeMesh(mesh, false, true);

                // Update fields (problem when face added to zero sized patch)
                mesh.updateMesh(resetMap);

                // Move mesh if in inflation mode
                if (resetMap().hasMotionPoints())
                {
                    mesh.movePoints(resetMap().preMotionPoints());
                }
                else
                {
                    mesh.clearOut();
                }
                meshRefiner_.updateMesh(resetMap,labelList(0));

                resetPrimitivePatchAddressing(adaptPatchIDs,pp);

                boolList updateSetFaces(mesh.nFaces(), false);

                forAll(setFaces, faceI)
                {
                    if (setFaces[faceI])
                    {
                        updateSetFaces[resetMap().reverseFaceMap()[faceI]]
                            = true;
                    }
                }
                setFaces = updateSetFaces;
            }
        }
    }

    return;
}


//updates relevant objects to take account of face splits in the patch surface
//(which may have altered the set of nodes that lie in the patch)
void Foam::snappySnapDriver::updateAfterMeshMods
(
    const mapPolyMesh &map,
    const fvMesh &mesh,
    const labelList &adaptPatchIDs,
    indirectPrimitivePatch &pp,
    featureLineSnapper &flSnap
)
{
    //modify the surface patch so it sees the new faces and edges

    // save the old mesh points map
    labelList ppPointMap(pp.meshPoints().begin(),pp.meshPoints().end());

    //update pp addressing
    resetPrimitivePatchAddressing(adaptPatchIDs,pp);

    //create mapping from old mesh numbers to new mesh numbers

    //now renumber the mesh point map so it becomes an old-to-new local
    //point mapping
    forAll(ppPointMap,pI)
    {
        label newMeshPointI = map.reversePointMap()[ppPointMap[pI]];

        Map<label>::const_iterator findEntry =
            pp.meshPointMap().find(newMeshPointI);

        if (findEntry != pp.meshPointMap().end())
        {
            ppPointMap[pI] = *findEntry;
        }
        else
        {
            ppPointMap[pI] = -1; //node has been removed from the patch
        }
    }

    //reset the shadow chains according to the new local numbering
    flSnap.renumber(ppPointMap);
}


Foam::PackedList<1> Foam::snappySnapDriver::setZonedFaces
(
    const meshRefinement& meshRefiner
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    PackedList<1> isZonedFace(mesh.nFaces(), 0);

    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = meshRefiner.getZones(fzTypes);
    }

    const faceZoneMesh& fZones = mesh.faceZones();

    forAll(internalOrBaffleFaceZones, i)
    {
        label zoneI = internalOrBaffleFaceZones[i];
        const faceZone& fZone = fZones[zoneI];
        forAll(fZone, i)
        {
            label faceI = fZone[i];
            isZonedFace.set(faceI, 1);
        }
    }

    return isZonedFace;
}

Foam::PackedList<1> Foam::snappySnapDriver::setZonedFacesExclBoundaryNamed
(
    const meshRefinement& meshRefiner,
    const refinementParameters& refineParams
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    PackedList<1> isZonedFace(mesh.nFaces(), 0);

    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = meshRefiner.getZones(fzTypes);
    }

    const faceZoneMesh& fZones = mesh.faceZones();

    forAll(internalOrBaffleFaceZones, i)
    {
        label zoneI = internalOrBaffleFaceZones[i];
        const faceZone& fZone = fZones[zoneI];

        forAll(fZone, i)
        {
            label faceI = fZone[i];
            isZonedFace.set(faceI, 1);
        }
    }

    forAll(fZones, zoneI)
    {
        const faceZone& fZone = fZones[zoneI];
        if (refineParams.isFaceZoneControlsBoundary(fZone.name()))
        {
            forAll(fZone, i)
            {
                label faceI = fZone[i];
                isZonedFace.set(faceI, 1);
            }
        }
    }

    return isZonedFace;
}


Foam::PackedList<1> Foam::snappySnapDriver::setZonedBaffleCheckFaces
(
    const meshRefinement& meshRefiner,
    const refinementParameters& refineParams
)
{
    const fvMesh& mesh = meshRefiner.mesh();

    PackedList<1> isZonedFace(mesh.nFaces(), 0);

    const faceZoneMesh& fZones = mesh.faceZones();
    const refinementSurfaces& surfaces = meshRefiner.surfaces();

    const PtrList<surfaceZonesInfo>& surfZones = surfaces.surfZones();

    forAll(surfZones, surfI)
    {
        if (!surfZones[surfI].baffleCheck())
        {
            const word& faceZoneName = surfZones[surfI].faceZoneName();

            if (faceZoneName.size())
            {
                // Filter out all faces for this zone.
                label zoneI = fZones.findZoneID(faceZoneName);
                const faceZone& fZone = fZones[zoneI];
                forAll(fZone, i)
                {
                    label faceI = fZone[i];
                    isZonedFace.set(faceI, 1);
                }
            }
        }
    }

    forAll(fZones, zoneI)
    {
        const word& faceZoneName = fZones[zoneI].name();

        if (!refineParams.isFaceZoneBaffleChecks(faceZoneName))
        {
            const faceZone& fZone = fZones[zoneI];
            forAll(fZone, i)
            {
                label faceI = fZone[i];
                isZonedFace.set(faceI, 1);
            }
        }
    }

    return isZonedFace;
}


void Foam::snappySnapDriver::resetPrimitivePatchAddressing
(
    const labelList& adaptPatchIDs,
    indirectPrimitivePatch &pp
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    pp.clearOut();
    labelList newAddressing;
    meshRefinement::calcPatchAddressing(mesh,adaptPatchIDs,newAddressing);
    pp.resetAddressing(newAddressing);
    //need to setup new addressing just in case required during a subsequent map
    pp.localPoints();
}
