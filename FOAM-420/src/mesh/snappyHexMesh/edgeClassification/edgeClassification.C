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
    (c) 2017-1017 Esi Ltd.
\*---------------------------------------------------------------------------*/

#include "edgeClassification/edgeClassification.H"
#include "meshes/polyMesh/syncTools/dummyTransform.H"
#include "meshTools/meshTools.H"
#include "meshes/meshShapes/edge/EdgeMap.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "externalDisplacementMeshMover/fieldSmoother/fieldSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(edgeClassification, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::vectorField Foam::edgeClassification::calculatePointNormals
(
    const boolList& excludedFaces,
    const label nSmoothIter,
    const bool correctBoundaryNormals
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& meshPoints = pp_.meshPoints();
    vectorField faceNormals(pp_.size(),vector::zero);
    forAll(faceNormals, facei)
    {
        faceNormals[facei] =
            mesh_.faces()[pp_.addressing()[facei]].unitNormal(pts_);
    }

    vectorField pointNormals(pp_.nPoints(), vector::zero);
    scalarField nPointFaces(pp_.nPoints(), 0.);

    boolList externalPoints(mesh_.nPoints(), false);
    boolList externalEdges(mesh_.nEdges(), false);

    labelList nConvexEdges(pp_.nPoints(),0);
    labelList nConcaveEdges(pp_.nPoints(),0);

    const PackedBoolList isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh_,
            ppMeshEdges_
         )
    );

    forAll(pp_.meshPoints(), ptI)
    {
        const labelList& pEdges = pp_.pointEdges()[ptI];
        forAll(pEdges, pEI)
        {
            label edgei = pEdges[pEI];
            if (!isPatchMasterEdge[edgei])
            {
                continue;
            }

            if
            (
                eType_[edgei].first() == edgeClassification::CONVEX
                || eType_[edgei].first() == edgeClassification::BAFFLE
            )
            {
                nConvexEdges[ptI]++;
            }
            else if (eType_[edgei].first() == edgeClassification::CONCAVE)
            {
                nConcaveEdges[ptI]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        nConvexEdges,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        nConcaveEdges,
        plusEqOp<label>(),  // combine op
        label(0)            // null value
    );

    forAll(meshPoints, pointi)
    {
        label meshPointI = meshPoints[pointi];
        if (nConcaveEdges[pointi] > 3 && nConvexEdges[pointi] == 0)
        {
            const labelList& pFaces = pp_.pointFaces()[pointi];
            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                label meshFaceI = pp_.addressing()[facei];
                vector fN = mesh_.faces()[meshFaceI].unitNormal(pts_);
                externalPoints[meshPointI] = true;
                pointNormals[pointi] += fN;
                nPointFaces[pointi] += scalar(1.);
            }
            const labelList& pEdges = pp_.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                if
                (
                    eType_[edgei].first() == edgeClassification::CONCAVE
                )
                {
                    label meshedgei = ppMeshEdges_[edgei];
                    externalEdges[meshedgei] = true;
                }
            }
        }
        else if ((nConcaveEdges[pointi]+nConvexEdges[pointi]) > 0)
        {
            const labelList& pEdges = pp_.pointEdges()[pointi];
            forAll(pEdges, pEI)
            {
                label edgei = pEdges[pEI];
                bool concaveEdge = false;
                if (eType_[edgei].first() == edgeClassification::CONCAVE)
                {
                    concaveEdge = true;
                }

                bool convexEdge = false;
                if
                (
                    eType_[edgei].first() == edgeClassification::CONVEX
                    || eType_[edgei].first() == edgeClassification::BAFFLE
                )
                {
                    convexEdge = true;
                }

                if (concaveEdge || convexEdge)
                {
                    label meshedgei = ppMeshEdges_[edgei];
                    edge e = mesh_.edges()[meshedgei];
                    const labelList& edgeFaces = pp_.edgeFaces()[edgei];
                    forAll(edgeFaces, ppEFI)
                    {
                        label meshFaceI = pp_.addressing()[edgeFaces[ppEFI]];
                        face f = mesh_.faces()[meshFaceI];
                        vector eVec = e.vec(pts_);
                        label fp = findIndex(f, e[0]);

                        vector eNorm =
                            (eVec ^ mesh_.faces()[meshFaceI].unitNormal(pts_));

                        if (concaveEdge && f[f.fcIndex(fp)] != e[1])
                        {
                            eNorm = -eNorm;
                        }
                        else if (convexEdge && f[f.fcIndex(fp)] == e[1])
                        {
                            eNorm = -eNorm;
                        }

                        eNorm /= (mag(eNorm) + SMALL);

                        pointNormals[pointi] += eNorm;
                        nPointFaces[pointi] += scalar(1.);
                    }

                    externalEdges[meshedgei] = true;
                    externalPoints[meshPointI] = true;
                }
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        externalEdges,
        orEqOp<bool>(),
        false
    );

    syncTools::syncPointList
    (
        mesh_,
        externalPoints,
        orEqOp<bool>(),
        false
    );

    forAll(faceNormals, facei)
    {
        if (excludedFaces[facei])
        {
            continue;
        }
        const face& f = pp_.localFaces()[facei];

        forAll(f, fp)
        {
            label meshPointI = meshPoints[f[fp]];

            if (!externalPoints[meshPointI])
            {
                scalar fA = mesh_.faces()[pp_.addressing()[facei]].mag(pts_);
                pointNormals[f[fp]] += faceNormals[facei]*fA;
                nPointFaces[f[fp]] += fA;
            }
        }
    }

    forAll(eType_, edgei)
    {
        if
        (
            eType_[edgei].first() == edgeClassification::BOUNDARY
            || eType_[edgei].first() == edgeClassification::NONMANIFOLD
        )
        {
             //Set stationary external points
            label meshedgei = ppMeshEdges_[edgei];
            edge e = mesh_.edges()[meshedgei];
            externalEdges[meshedgei] = true;
            externalPoints[e[0]] = true;
            externalPoints[e[1]] = true;
        }
    }


    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        pointNormals,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        meshPoints,
        nPointFaces,
        plusEqOp<scalar>(),
        scalar(0.)          // null value
    );

    forAll(pointNormals, i)
    {
        if (nPointFaces[i] > SMALL)
        {
            pointNormals[i] /= nPointFaces[i];
        }
    }
    pointNormals /= (mag(pointNormals) + SMALL);

    if (correctBoundaryNormals)
    {
        boolList markedFaces(mesh_.nFaces(), false);
        boolList bdyPts(pp_.meshPoints().size(), false);

        forAll(pp_, i)
        {
            markedFaces[pp_.addressing()[i]] = true;
        }

        forAll(eType_, edgei)
        {
            if (eType_[edgei].first() == edgeClassification::BOUNDARY)
            {
                const edge e = pp_.edges()[edgei];
                bdyPts[e[0]] = true;
                bdyPts[e[1]] = true;
            }
        }

        labelList nNbrFaces(pp_.meshPoints().size(), 0);
        vectorField aveNbrNorm(pp_.meshPoints().size(), vector::zero);

        forAll(bdyPts, pti)
        {
            if (bdyPts[pti])
            {
                label meshpointi = meshPoints[pti];
                const labelList& pFaces = mesh_.pointFaces()[meshpointi];
                forAll(pFaces, pfi)
                {
                    label facei = pFaces[pfi];
                    if (!markedFaces[facei])
                    {
                        label patchi = patches.whichPatch(facei);
                        if (patchi != -1 && !patches[patchi].coupled())
                        {
                            nNbrFaces[pti]++;
                            aveNbrNorm[pti] += mesh_.faceAreas()[facei];
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            meshPoints,
            aveNbrNorm,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh_,
            meshPoints,
            nNbrFaces,
            plusEqOp<label>(),
            label(0)        // null value
        );

        forAll(aveNbrNorm, i)
        {
            if (nNbrFaces[i] > 0)
            {
                aveNbrNorm[i] /= nNbrFaces[i];
                aveNbrNorm[i] /= (mag(aveNbrNorm[i]) + SMALL);
            }
        }

        forAll(pointNormals, pti)
        {
            if (bdyPts[pti])
            {
                if (nNbrFaces[pti] == 2)
                {
                    label meshpointi = meshPoints[pti];
                    const labelList& pFaces = mesh_.pointFaces()[meshpointi];
                    forAll(pFaces, pfi)
                    {
                        label facei = pFaces[pfi];
                        if (!markedFaces[facei])
                        {
                            label patchi = patches.whichPatch(facei);
                            if (patchi != -1 && !patches[patchi].coupled())
                            {
                                vector fA = mesh_.faceAreas()[facei];
                                fA /= (mag(fA) + SMALL);
                                if ((fA & aveNbrNorm[pti]) > 0.984)
                                {
                                    point pt = pts_[meshpointi];
                                    plane pl(pt, aveNbrNorm[pti]);
                                    point outerPt = pt + pointNormals[pti];
                                    point nearestPlanePt =
                                        pl.nearestPoint(outerPt);
                                    vector gUp = nearestPlanePt - pt;
                                    gUp /= (mag(gUp) + SMALL);
                                    pointNormals[pti] = gUp;
                                    break;
                                }
                                else
                                {
                                    const labelList& pEdges =
                                        mesh_.pointEdges()[meshpointi];
                                    label markedEdge = -1;
                                    forAll(pEdges, pei)
                                    {
                                        label edgei = pEdges[pei];
                                        const labelList& eFaces =
                                            mesh_.edgeFaces()[edgei];
                                        bool foundMarkedFace = false;
                                        bool foundBoundaryFace = false;
                                        forAll(eFaces, efi)
                                        {
                                            label efacei = eFaces[efi];
                                            patchi = patches.whichPatch(facei);
                                            if
                                            (
                                                !markedFaces[efacei]
                                                && patchi != -1
                                                && !patches[patchi].coupled()
                                            )
                                            {
                                                foundBoundaryFace = true;
                                            }

                                            if (markedFaces[efacei])
                                            {
                                                foundMarkedFace = true;
                                            }
                                        }
                                        if
                                        (
                                            !foundMarkedFace
                                            && foundBoundaryFace
                                        )
                                        {
                                            markedEdge = edgei;
                                            break;
                                        }
                                    }
                                    if (markedEdge != -1)
                                    {
                                        edge e = mesh_.edges()[markedEdge];
                                        label otherPt
                                        (
                                            e[0] == meshpointi ? e[1] : e[0]
                                        );
                                        vector eDir = pts_[meshpointi]
                                            - pts_[otherPt];
                                        eDir /= (mag(eDir) + SMALL);
                                        pointNormals[pti] = eDir;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (debug)
    {
        pointField sPts(meshPoints.size());
        forAll(meshPoints, pti)
        {
            sPts[pti] = pts_[meshPoints[pti]];
        }
        simpleVTKWriter layerVTK
        (
            pp_.localFaces(),
            sPts
        );
        layerVTK.addPointData("pNormals", pointNormals);
        layerVTK.write("pNormUnSmooth.vtk");
    }

    if (nSmoothIter > 0)
    {
        fieldSmoother fSmoother(mesh_);

        // Precalulate (patch) master point/edge
        const PackedBoolList isPatchMasterPoint
        (
            meshRefinement::getMasterPoints
            (
                mesh_,
                meshPoints
             )
        );

        boolList cornerPoints(mesh_.nPoints(), false);
        // Smooth patch normal vectors
        fSmoother.smoothPatchNormals
        (
            nSmoothIter,
            isPatchMasterPoint,
            isPatchMasterEdge,
            externalEdges,
            externalPoints,
            cornerPoints,
            ppMeshEdges_,
            pp_,
            pointNormals
        );
    }

    if (debug)
    {
        pointField sPts(meshPoints.size());
        forAll(meshPoints, pti)
        {
            sPts[pti] = pts_[meshPoints[pti]];
        }
        simpleVTKWriter layerVTK
        (
            pp_.localFaces(),
            sPts
        );

        layerVTK.addPointData("pNormals", pointNormals);
        layerVTK.write("pNormSmoothed.vtk");
    }

    return pointNormals;
}


Foam::List<Foam::Tuple2<Foam::edgeClassification::edgeType,Foam::scalar>>
Foam::edgeClassification::setEdgeType
(
    const boolList& excludedFaces,
    const scalar convexFeatureAngle,
    const scalar concaveFeatureAngle,
    const scalar baffleFeatureAngle,
    const bool zoneFlipFaces
)
{
    scalar featureAngle = max(convexFeatureAngle,concaveFeatureAngle);

    label nPPEdges = pp_.edges().size();
    List<pointField> edgeNormals(nPPEdges,pointField(0));

    labelList neiCellZone;
    if (zoneFlipFaces)
    {
        labelList cellToZone(mesh_.nCells(), -1);
        forAll(mesh_.cells(), celli)
        {
           label zonei = mesh_.cellZones().whichZone(celli);
           cellToZone[celli] = zonei;
        }
        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
    }

    forAll(pp_.edges(), edgei)
    {
        const labelList& eFaces = pp_.edgeFaces()[edgei];
        pointField bVecs(eFaces.size());
        label nPts = 0;
        forAll(eFaces, eFI)
        {
            label facei = pp_.addressing()[eFaces[eFI]];
            bool flip = false;
            if (zoneFlipFaces)
            {
                label ownZoneI = mesh_.cellZones().whichZone
                (
                    mesh_.faceOwner()[facei]
                );

                label  neiZone = -1;
                if (mesh_.isInternalFace(facei))
                {
                    neiZone = mesh_.cellZones().whichZone
                    (
                        mesh_.faceNeighbour()[facei]
                    );
                }
                else
                {
                    neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                }
                if (ownZoneI < neiZone)
                {
                    flip = true;
                }
            }
            point fA = mesh_.faces()[facei].areaNormal(pts_);

            vector fN
            (
                flip ? -fA : fA
            );
            fN /= (mag(fN) + SMALL);
            bVecs[nPts++] = fN;
        }
        bVecs.setSize(nPts);
        edgeNormals[edgei]  = std::move(bVecs);
    }

    syncTools::syncEdgeList
    (
        mesh_,
        ppMeshEdges_,
        edgeNormals,
        edgeClassification::sumEqOp(),
        pointField(0)          // null value
    );

    scalarField eDotProd(nPPEdges, GREAT);
    forAll(pp_.edges(), edgei)
    {
        const pointField& eN =  edgeNormals[edgei];
        for (int i = 0; i < eN.size()-1; i++)
        {
            for (int j = i+1; j < eN.size(); j++)
            {
                eDotProd[edgei] = min(eDotProd[edgei],(eN[i] & eN[j]));
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        ppMeshEdges_,
        eDotProd,
        minEqOp<scalar>(),
        GREAT          // null value
    );

    labelList eT(nPPEdges, -1);
    forAll(pp_, i)
    {
        label facei = pp_.addressing()[i];
        const labelList& fEdges = pp_.faceEdges()[i];
        forAll(fEdges, fEI)
        {
            label edgei = fEdges[fEI];
            label meshedgei = ppMeshEdges_[edgei];
            scalar eAngle = eDotProd[edgei];

            if (eAngle < featureAngle)
            {
                if (eAngle < -0.9998)
                {
                    if (eAngle < baffleFeatureAngle)
                    {
                        eT[edgei] = 3;
                    }
                    else if (eAngle < convexFeatureAngle)
                    {
                        eT[edgei] = 2;
                    }
                    else
                    {
                        eT[edgei] = 0;
                    }
                }
                else
                {
                    edge e = mesh_.edges()[meshedgei];
                    point fC = mesh_.faces()[facei].centre(pts_);

                    point lp0 = pts_[e[0]];
                    point lp1 = pts_[e[1]];
                    vector offset = (lp1-lp0);
                    lp0 -= offset;
                    lp1 += offset;
                    pointHit lHit = linePointRef(lp0, lp1).nearestDist(fC);

                    vector eCtofC = vector::zero;
                    if (lHit.hit())
                    {
                        eCtofC = lHit.hitPoint() - fC ;
                        eCtofC /= (mag(eCtofC) + SMALL);
                    }
                    else
                    {
                        point eC = e.centre(pts_);
                        eCtofC = eC - fC;
                        eCtofC /= (mag(eCtofC) + SMALL);
                    }
                    const pointField& eN =  edgeNormals[edgei];
                    scalar maxDotProd = -GREAT;
                    scalar minDotProd = GREAT;
                    forAll(eN, eNI)
                    {
                        scalar dProd = (eN[eNI] & eCtofC);

                        if (mag(dProd) > maxDotProd)
                        {
                            maxDotProd = mag(dProd);
                            minDotProd = dProd;
                        }
                    }

                    if (minDotProd > 0)
                    {
                        if (eAngle < concaveFeatureAngle)
                        {
                            eT[edgei] = 1;
                        }
                        else
                        {
                            eT[edgei] = 0;
                        }
                    }
                    else
                    {
                        if (eAngle < baffleFeatureAngle)
                        {
                            eT[edgei] = 3;
                        }
                        else if (eAngle < convexFeatureAngle)
                        {
                            eT[edgei] = 2;
                        }
                        else
                        {
                            eT[edgei] = 0;
                        }
                    }
                }
            }
            else
            {
                eT[edgei] = 0;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        ppMeshEdges_,
        eT,
        maxEqOp<label>(),
        label(-1)
    );

    labelList nExternalEdge(nPPEdges, 0);
    const labelListList& edgeFaces = pp_.edgeFaces();
    forAll(pp_.edges(), edgei)
    {
        const labelList& eFaces = edgeFaces[edgei];
        nExternalEdge[edgei] = eFaces.size();
    }

    syncTools::syncEdgeList
    (
        mesh_,
        ppMeshEdges_,
        nExternalEdge,
        plusEqOp<label>(),
        label(0)              // null value
    );

    List<Tuple2<edgeType,scalar>> eType
    (
        pp_.edges().size(),
        Tuple2<edgeType,scalar>(MANIFOLD,-GREAT)
    );

    forAll(ppMeshEdges_, edgei)
    {
        scalar eAngle = eDotProd[edgei];
        if (nExternalEdge[edgei] == 1)
        {
            eType[edgei] = Tuple2<edgeType,scalar>(BOUNDARY,-GREAT);
        }
        else if (nExternalEdge[edgei] > 2)
        {
            eType[edgei] = Tuple2<edgeType,scalar>(NONMANIFOLD,-GREAT);
        }
        else
        {
            if (eT[edgei] == 0)
            {
                eType[edgei] = Tuple2<edgeType,scalar>(MANIFOLD,eAngle);
            }
            if (eT[edgei] == 1)
            {
                eType[edgei] = Tuple2<edgeType,scalar>(CONCAVE,eAngle);
            }
            else if (eT[edgei] == 2)
            {
                eType[edgei] = Tuple2<edgeType,scalar>(CONVEX,eAngle);
            }
            else if (eT[edgei] == 3)
            {
                eType[edgei] = Tuple2<edgeType,scalar>(BAFFLE,eAngle);
            }
            else if (eT[edgei] != 0)
            {
                label meshedgei = ppMeshEdges_[edgei];
                edge e = mesh_.edges()[meshedgei];
                WarningInFunction
                    << "Cannot characterise feature type for boundary edge : "
                    << meshedgei << " at location : " << e.centre(pts_)
                    << endl;
            }
        }
    }

    return eType;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeClassification::edgeClassification
(
    const polyMesh& mesh,
    const pointField& pts,
    const indirectPrimitivePatch& pp,
    const labelList& ppMeshEdges,
    const boolList& excludedFaces,
    const scalar convexFeatureAngle,
    const scalar concaveFeatureAngle,
    const scalar baffleFeatureAngle,
    const bool zoneFlipFaces
)
:
    mesh_(mesh),
    pts_(pts),
    pp_(pp),
    ppMeshEdges_(ppMeshEdges),
    eType_
    (
        setEdgeType
        (
            excludedFaces,
            convexFeatureAngle,
            concaveFeatureAngle,
            baffleFeatureAngle,
            zoneFlipFaces
        )
    )
{}

// ************************************************************************* //
