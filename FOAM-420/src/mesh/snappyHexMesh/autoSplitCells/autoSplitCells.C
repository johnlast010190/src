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
    (c) 2020 Esi Ltd.

InClass
    autoSplitCells

\*---------------------------------------------------------------------------*/

#include "autoSplitCells/autoSplitCells.H"

#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveCell.H"

#include "edgeClassification/edgeClassification.H"

#include "snappyHexMeshDriver/snappySnapDriver.H"

#include "refinementSurfaces/refinementSurfaces.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"

#include "meshCut/cellCuts/cellCuts.H"
#include "meshCut/meshModifiers/meshCutter/meshCutter.H"

#include "meshes/polyMesh/syncTools/dummyTransform.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"
#include <list>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(autoSplitCells, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::autoSplitCells::isValidPatch(const label patchi)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const polyPatch& pp = patches[patchi];

    if (pp.coupled() || excludePatchSet_.found(patchi))
    {
        return false;
    }
    else
    {
        return true;
    }
}


Foam::label Foam::autoSplitCells::checkEdges
(
    const pointField& origPoints,
    labelHashSet& setPtr
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();
    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();

    boolList markedPoints(mesh.nPoints(), false);
    forAll(mesh.edges(), edgei)
    {
        const edge& e =  mesh.edges()[edgei];
        scalar newEdgeLength = e.mag(mesh.points());
        scalar oldEdgeLength = e.mag(origPoints);

        label pLev = max(pointLevel[e[0]],pointLevel[e[1]]);
        scalar len = edge0Len / pow(2., pLev);
        if (newEdgeLength < 0.2*len && newEdgeLength < oldEdgeLength)
        {
            markedPoints[e[0]] = true;
            markedPoints[e[1]] = true;
        }
    }
    syncTools::syncPointList(mesh, markedPoints, orEqOp<bool>(), false);

    label nReset = 0;
    forAll(markedPoints, pointi)
    {
       if (markedPoints[pointi])
       {
           setPtr.insert(pointi);
           nReset++;
       }
    }
    reduce(nReset, sumOp<label>());

    return nReset;
}


void Foam::autoSplitCells::correctPoints
(
    vectorField& disp
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    pointField origPoints = mesh.points();
    pointField newPoints( mesh.points() + disp );

    mesh.movePoints(newPoints);

    //Undo point motion if errors generated
    disp = vector::zero;

    cellSet errorCells(mesh, "errorCells", mesh.nCells()/100+1);
    faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+1);
    pointSet errorPoints(mesh, "errorPoints", mesh.nPoints()/100+1);

    label noFailedChecks(0);
    if (mesh.checkClosedCells(true, &errorCells)) noFailedChecks++;
    if (revertErrors_)
    {
        if (mesh.checkCellVolumes(true, &errorCells)) noFailedChecks++;
        if (mesh.checkFaceAreas(true, &errorFaces)) noFailedChecks++;
        if (mesh.checkFaceOrthogonality(true, &errorFaces)) noFailedChecks++;
        if (mesh.checkFacePyramids(true, -SMALL, &errorFaces)) noFailedChecks++;
        if (checkEdges(origPoints, errorPoints) > 0) noFailedChecks++;
    }

    if (noFailedChecks > 0)
    {
        forAllConstIter(labelHashSet, errorPoints, iter)
        {
            label pointI = iter.key();
            disp[pointI] = origPoints[pointI] - newPoints[pointI];
        }

        forAllConstIter(labelHashSet, errorFaces, iter)
        {
            label own = mesh.faceOwner()[iter.key()];
            const labelList& ownPts = mesh.cellPoints()[own];
            forAll(ownPts, ptI)
            {
                label pointI = ownPts[ptI];
                disp[pointI] = origPoints[pointI] - newPoints[pointI];
            }

            if (iter.key() < mesh.nInternalFaces())
            {
                label nei = mesh.faceNeighbour()[iter.key()];
                const labelList& neiPts = mesh.cellPoints()[nei];
                forAll(neiPts, ptI)
                {
                    label pointI = neiPts[ptI];
                    disp[pointI] = origPoints[pointI] - newPoints[pointI];
                }
            }
        }

        forAllConstIter(labelHashSet, errorCells, iter)
        {
            const labelList& cellPts = mesh.cellPoints()[iter.key()];
            forAll(cellPts, ptI)
            {
                label pointI = cellPts[ptI];
                disp[pointI] = origPoints[pointI] - newPoints[pointI];
            }
        }

        //move points
        syncTools::syncPointList
        (
            mesh,
            disp,
            maxMagSqrEqOp<point>(),
            vector::zero
        );
        newPoints = mesh.points() + disp;
        mesh.movePoints(newPoints);
    }
}

void Foam::autoSplitCells::moveFeaturePoints
(
    const List<Tuple2<label,labelPair>>& cutFaces,
    const scalar featureAngle,
    const bool moveFeaturePts,
    const bool moveConvexPts
)
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    boolList moved(mesh.nPoints(), false);
    vectorField disp(mesh.nPoints(), vector::zero);

    if (moveFeaturePts)
    {
        forAll(cutFaces, i)
        {
            label facei = cutFaces[i].first();
            face f = mesh.faces()[facei];

            DynamicList<label> anchorPts(f.size());
            {
                label celli = mesh.faceOwner()[facei];
                label cLevel = cellLevel[celli];

                forAll(f,fp)
                {
                    label pointi = f[fp];
                    label pLevel = pointLevel[pointi];
                    if (extruded_)
                    {
                        if (pLevel <= cLevel)
                        {
                            anchorPts.append(pointi);
                        }
                    }
                    else
                    {
                        anchorPts.append(pointi);
                    }
                }
            }

            face af(anchorPts);
            if (af.size() != 4)
            {
                continue;
            }

            label mPt = -1;
            vector midPoint;
            scalar maxDProd = -GREAT;
            forAll(af,fp)
            {
                label pointi = af[fp];
                point pt = mesh.points()[pointi];

                label nextFp = af.fcIndex(fp);
                label prevFp = af.rcIndex(fp);

                vector nextVec = mesh.points()[af[nextFp]] - pt;
                vector prevVec = pt - mesh.points()[af[prevFp]];
                nextVec /= Foam::mag(nextVec) + VSMALL;
                prevVec /= Foam::mag(prevVec) + VSMALL;
                scalar dProd = (nextVec & prevVec);

                if (dProd > featureAngle)
                {
                    if (dProd > maxDProd)
                    {
                        maxDProd = dProd;
                        mPt = pointi;

                        midPoint = 0.5*
                        (
                            mesh.points()[af[prevFp]] +
                            mesh.points()[af[nextFp]]
                        );
                    }
                }
            }

            if (mPt != -1)
            {
                disp[mPt] = midPoint - mesh.points()[mPt];
                moved[mPt] = true;
            }
        }
        //sync moved points
        syncTools::syncPointList
        (
            mesh,
            moved,
            orEqOp<bool>(),
            false
        );
    }

    if (moveConvexPts)
    {
        //Move convex face points to mid point
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (isValidPatch(patchI))
            {
                label startFaceI = pp.start();

                forAll(pp, i)
                {
                    label faceI = startFaceI+i;
                    face f = mesh.faces()[faceI];

                    vectorField fEdgeVec(f.size());
                    forAll(f, fp)
                    {
                        label nextFp = f.fcIndex(fp);

                        point thisPt = mesh.points()[f[fp]];
                        point nextPt = mesh.points()[f[nextFp]];
                        fEdgeVec[fp] = (nextPt - thisPt);
                    }

                    vector fA = mesh.faceAreas()[faceI];
                    forAll(f, fp)
                    {
                        label pointi = f[fp];
                        if (moved[pointi])
                        {
                            continue;
                        }
                        label prevFp = f.rcIndex(fp);
                        vector thisVec = fEdgeVec[fp];
                        vector prevVec = fEdgeVec[prevFp];

                        vector pN(thisVec ^ prevVec);
                        if ((pN&fA) > SMALL)
                        {
                            label nextFp = f.fcIndex(fp);
                            pointHit lHit = linePointRef
                            (
                                mesh.points()[f[prevFp]],
                                mesh.points()[f[nextFp]]
                            ).nearestDist(mesh.points()[pointi]);

                            if (lHit.hit())
                            {
                                scalar eLength1 = mag
                                (
                                    lHit.hitPoint()-mesh.points()[f[prevFp]]
                                );
                                scalar eLength2 = mag
                                (
                                    lHit.hitPoint()-mesh.points()[f[nextFp]]
                                );

                                scalar ratio = min(eLength1, eLength2)
                                    / max(eLength1, eLength2);

                                if (ratio > 0.1)
                                {
                                    point midPoint = 0.5*
                                    (
                                        mesh.points()[f[prevFp]] +
                                        mesh.points()[f[nextFp]]
                                    );
                                    disp[pointi] = midPoint
                                    - mesh.points()[pointi];
                                    moved[pointi] = true;
                                }
                                else
                                {
                                    if (eLength1 > SMALL && eLength2 > SMALL)
                                    {
                                        disp[pointi] = lHit.hitPoint()
                                            - mesh.points()[pointi];
                                        moved[pointi] = true;
                                    }
                                    else
                                    {
                                        disp[pointi] =  vector::zero;
                                    }
                                }
                            }
                            else
                            {
                                disp[pointi] =  vector::zero;
                            }
                        }
                    }
                }
            }
        }
        //sync moved points
        syncTools::syncPointList
        (
            mesh,
            moved,
            orEqOp<bool>(),
            false
        );
    }

    //move points
    syncTools::syncPointList
    (
        mesh,
        disp,
        maxMagSqrEqOp<point>(),
        vector::zero
    );

    correctPoints(disp);
}


void Foam::autoSplitCells::moveSplitFaceHangingPts()
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    vectorField disp(mesh.nPoints(), vector::zero);

    //Move hanging boundary face points to mid point
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isValidPatch(patchI))
        {
            label startFacei = pp.start();

            forAll(pp, i)
            {
                label facei = startFacei+i;
                label own = mesh.faceOwner()[facei];
                label cLevel = cellLevel[own];

                face f = mesh.faces()[facei];

                DynamicList<label> anchors(f.size());
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    label pLevel = pointLevel[pointi];
                    if (pLevel <= cLevel)
                    {
                        anchors.append(pointi);
                    }
                }
                label nAnchors = anchors.size();
                if (nAnchors == 3 && f.size() != nAnchors)
                {
                    face af(anchors);
                    point fc = af.centre(mesh.points());
                    point fn = af.areaNormal(mesh.points());

                    if (magSqr(fn) < SMALL)
                    {
                        continue;
                    }

                    plane pl(fc, fn);
                    forAll(f, fp)
                    {
                        label pointi = f[fp];
                        label pLevel = pointLevel[pointi];
                        if (pLevel > cLevel)
                        {
                            point pt = mesh.points()[pointi];
                            disp[pointi] = pl.nearestPoint(pt) - pt;
                        }
                    }
                }
            }
        }
    }

    //sync displacement
    syncTools::syncPointList
    (
        mesh,
        disp,
        maxMagSqrEqOp<point>(),
        vector::zero
    );

    //Move and correct point displacement
    correctPoints(disp);
}


void Foam::autoSplitCells::flattenRefinementFaces()
{
    fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    vectorField disp(mesh.nPoints(), vector::zero);

    //Move hanging boundary face points to mid point
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isValidPatch(patchI))
        {
            label startFacei = pp.start();

            forAll(pp, i)
            {
                label facei = startFacei+i;
                label own = mesh.faceOwner()[facei];
                label cLevel = cellLevel[own];

                face f = mesh.faces()[facei];

                bool hangingFace = false;
                forAll(f, fp)
                {
                    label pointi = f[fp];
                    label pLevel = pointLevel[pointi];
                    if (pLevel > cLevel)
                    {
                        hangingFace = true;
                        break;
                    }
                }

                if (hangingFace)
                {
                    plane pl
                    (
                        mesh.faceCentres()[facei],
                        mesh.faceAreas()[facei]
                    );
                    forAll(f, fp)
                    {
                        label pointi = f[fp];
                        point pt = mesh.points()[pointi];

                        disp[pointi] = pl.nearestPoint(pt) - pt;
                    }
                }
            }
        }
    }

    //Sync point displacement
    syncTools::syncPointList
    (
        mesh,
        disp,
        maxMagSqrEqOp<point>(),
        vector::zero
    );

    //Move and correct point displacement
    correctPoints(disp);
}


Foam::List<Foam::labelPair> Foam::autoSplitCells::split
(
    const List<Tuple2<label,labelPair>>& cutFaces,
    const bool checkCutAreas
)
{
    if (returnReduce(cutFaces.size(), sumOp<label>()) == 0)
    {
        return List<labelPair>(0);
    }

    fvMesh& mesh = meshRefiner_.mesh();

    polyTopoChange meshMod(mesh);
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    DynamicList<label> cutCells(cutFaces.size());
    DynamicList<labelList> cellLoops(cutFaces.size());
    DynamicList<labelList> cellAnchorPts(cutFaces.size());
    boolList cutTopFace(mesh.nFaces(), false);
    List<DynamicList<labelPair>> cutTopFacePts(mesh.nFaces());

    forAll(cutFaces, i)
    {
        labelList loop(4);
        label facei = cutFaces[i].first();
        const face& f =  mesh.faces()[facei];
        edge cutEdge = cutFaces[i].second();
        label own = mesh.faceOwner()[facei];

        if (checkCutAreas)
        {
            label nextFp = findIndex(f, cutEdge[1]);
            DynamicList<label> cutPts0(f.size());
            forAll(f,fp)
            {
                label pointi = f[nextFp];
                cutPts0.append(pointi);
                if (pointi == cutEdge[0])
                {
                    break;
                }
                nextFp = f.fcIndex(nextFp);
            }
            face cutFace0(cutPts0);

            nextFp = findIndex(f, cutEdge[0]);
            DynamicList<label> cutPts1(f.size());
            forAll(f,fp)
            {
                label pointi = f[nextFp];
                cutPts1.append(pointi);
                if (pointi == cutEdge[1])
                {
                    break;
                }
                nextFp = f.fcIndex(nextFp);
            }
            face cutFace1(cutPts1);

            scalar areaA = cutFace0.mag(mesh.points())+SMALL;
            scalar areaB = cutFace1.mag(mesh.points())+SMALL;
            scalar splitRatio = min(areaA,areaB) / max(areaA,areaB);

            if (splitRatio < 0.25)
            {
                continue;
            }
        }

        label nPts = 0;
        loop[nPts++] = cutEdge[1];
        loop[nPts++] = cutEdge[0];

        forAll(cutEdge, cPt)
        {
            label pointi = cutEdge[cPt];
            const labelList& pedges = mesh.pointEdges()[pointi];
            labelHashSet cEdgeSet(mesh.cellEdges()[own]);
            labelHashSet fEdgeSet(mesh.faceEdges()[facei]);

            forAll(pedges, pei)
            {
                label meshedgei = pedges[pei];
                if (cEdgeSet.found(meshedgei) && !fEdgeSet.found(meshedgei))
                {
                    const edge& e = mesh.edges()[meshedgei];
                    label otherPt(e[0] == pointi ? e[1] : e[0]);
                    loop[nPts++] = otherPt;
                    break;
                }
            }
        }

        labelList faces1 = mesh.pointFaces()[loop[2]];
        labelList faces2 = mesh.pointFaces()[loop[3]];
        std::list<label> intx;
        std::set_intersection
        (
            faces1.begin(),
            faces1.end(),
            faces2.begin(),
            faces2.end(),
            std::back_inserter(intx)
        );

        if (intx.size() == 1)
        {
            //Calculate anchor points for cut cell
            label nextFp = findIndex(f, cutEdge[0]);
            DynamicList<label> sidedPts(f.size());
            forAll(f,fp)
            {
                label pointi = f[nextFp];
                if (pointi == cutEdge[1])
                {
                    break;
                }
                if (pointi != cutEdge[0])
                {
                    sidedPts.append(pointi);
                }
                nextFp = f.fcIndex(nextFp);
            }
            labelHashSet cPtSet(mesh.cellPoints()[own]);
            labelHashSet bFaceSet(f);
            labelHashSet loopSet(loop);
            DynamicList<label> anchorPts(sidedPts.size()*2);
            forAll(sidedPts, pti)
            {
                label pointi = sidedPts[pti];
                anchorPts.append(pointi);
                const labelList& pEdges = mesh.pointEdges()[pointi];
                forAll(pEdges, pEI)
                {
                    label meshedgei = pEdges[pEI];
                    const edge& e = mesh.edges()[meshedgei];
                    label otherPt(e[0] == pointi ? e[1] : e[0]);
                    if
                    (
                        cPtSet.found(otherPt)
                        && !bFaceSet.found(otherPt)
                        && !loopSet.found(otherPt)
                    )
                    {
                        anchorPts.append(otherPt);
                    }
                }
            }

            label topFace = intx.front();
            cutCells.append(own);
            cellLoops.append(loop);
            cellAnchorPts.append(anchorPts);
            cutTopFace[topFace] = true;
            cutTopFacePts[topFace].append
            (
                labelPair(loop[2], loop[3])
            );
        }
    }
    syncTools::syncFaceList(mesh, cutTopFace, orEqOp<bool>());

    DynamicList<labelPair> splitCells(cutCells.size());

    label nCellsToSplit = returnReduce(cutCells.size(), sumOp<label>());
    if (nCellsToSplit > 0)
    {
        labelList neiCuts0(mesh.nFaces()-mesh.nInternalFaces(),-1);
        labelList neiCuts1(mesh.nFaces()-mesh.nInternalFaces(),-1);
        List<labelPair> masterCuts
            (mesh.nFaces()-mesh.nInternalFaces(),labelPair(-1,-1));
        forAll(patches, patchI)
        {
            if (patches[patchI].coupled())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(patches[patchI]);

                label startI = procPatch.start();

                const labelList& nbrPts = procPatch.nbrPoints();
                const labelList& meshPts = procPatch.meshPoints();

                forAll(procPatch, i)
                {
                    label faceI = startI + i;

                    if (cutTopFacePts[faceI].size() == 1)
                    {
                        labelPair cutsPts = cutTopFacePts[faceI][0];

                        face f = procPatch.localFaces()[i];
                        DynamicList<label> localCutsMaster(2);
                        DynamicList<label> localCuts(2);

                        forAll(f, fp)
                        {
                            label meshPointI = meshPts[f[fp]];

                            if
                            (
                                meshPointI == cutsPts[0]
                                || meshPointI == cutsPts[1]
                            )
                            {
                                localCutsMaster.append(f[fp]);
                                localCuts.append(nbrPts[f[fp]]);
                            }
                        }
                        neiCuts0[faceI-mesh.nInternalFaces()] = localCuts[0];
                        neiCuts1[faceI-mesh.nInternalFaces()] = localCuts[1];
                        masterCuts[faceI-mesh.nInternalFaces()] =
                            labelPair(localCutsMaster[0], localCutsMaster[1]);
                    }
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh, neiCuts0);
        syncTools::swapBoundaryFaceList(mesh, neiCuts1);

        List<labelPair> preCutProcFaces
            (mesh.nFaces()-mesh.nInternalFaces(),labelPair(-1,-1));

        DynamicList<label> newCutCells(cutCells.size());
        DynamicList<labelList> newCellLoops(cellLoops.size());
        DynamicList<labelList> newAnchorPts(cellLoops.size());

        boolList disable(mesh.nFaces(), false);

        forAll(mesh.faces(), faceI)
        {
            if (!cutTopFace[faceI])
            {
                continue;
            }

            label patchI = patches.whichPatch(faceI);

            if (patchI == -1 || patches[patchI].coupled())
            {
                if (patchI == -1)
                {
                    const List<labelPair>& cutFaces = cutTopFacePts[faceI];
                    if (cutFaces.size() == 2)
                    {
                        if
                        (
                            cutFaces[0][0] != cutFaces[1][0]
                            && cutFaces[0][0] != cutFaces[1][1]
                        )
                        {
                            disable[faceI] = true;
                        }
                    }
                }
                else
                {
                    const processorPolyPatch& procPatch =
                        refCast<const processorPolyPatch>(patches[patchI]);
                    const labelList& meshPts = procPatch.meshPoints();
                    label bFaceI = faceI-mesh.nInternalFaces();

                    labelPair mCut = masterCuts[bFaceI];
                    label sCut0 = neiCuts0[bFaceI];
                    label sCut1 = neiCuts1[bFaceI];

                    if (mCut[0] != -1 && sCut0 != -1)
                    {
                        if
                        (
                            (mCut[0] != sCut0 && mCut[0] != sCut1)
                        )
                        {
                            disable[faceI] = true;
                        }
                    }
                    else if
                    (
                        (mCut[0] == -1 && sCut0 != -1)
                        || (sCut1 == -1 && mCut[0] != -1)
                    )
                    {
                        if (sCut0 != -1)
                        {
                            preCutProcFaces[bFaceI][0] = meshPts[sCut0];
                            preCutProcFaces[bFaceI][1] = meshPts[sCut1];
                        }
                        else
                        {
                            preCutProcFaces[bFaceI][0] = meshPts[mCut[0]];
                            preCutProcFaces[bFaceI][1] = meshPts[mCut[1]];
                        }
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh, disable, orEqOp<bool>());

        forAll(cutCells, i)
        {
            label cellI = cutCells[i];

            const labelList& cFaces = mesh.cells()[cellI];
            bool removeCut = false;

            forAll(cFaces, cFI)
            {
                label faceI = cFaces[cFI];
                if (disable[faceI])
                {
                    removeCut = true;
                    break;
                }
            }

            if (!removeCut)
            {
                newCutCells.append(cutCells[i]);
                newCellLoops.append(cellLoops[i]);
                newAnchorPts.append(cellAnchorPts[i]);
            }
        }

        cutCells = newCutCells;
        cellLoops = newCellLoops;
        cellAnchorPts = newAnchorPts;

        nCellsToSplit = returnReduce(cutCells.size(), sumOp<label>());

        if (nCellsToSplit == 0)
        {
            return List<labelPair>(splitCells, true);
        }

        Info<<"Selected "<< nCellsToSplit <<" cells for splitting"<<endl;

        //Split processor faces
        {
            polyTopoChange meshMod(mesh);
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
                        labelPair cFaces =
                            preCutProcFaces[faceI-mesh.nInternalFaces()];
                        if (cFaces[0] != -1)
                        {
                            face f = mesh.faces()[faceI];

                            DynamicList<label> origFacePts(f.size());
                            label nextFp = findIndex(f, cFaces[0]);
                            forAll(f,fp)
                            {
                                label pointi = f[nextFp];
                                origFacePts.append(pointi);
                                if (pointi == cFaces[1])
                                {
                                    break;
                                }
                                nextFp = f.fcIndex(nextFp);
                            }
                            face origFace(origFacePts);

                            DynamicList<label> newFacePts(f.size());
                            nextFp = findIndex(f, cFaces[1]);
                            forAll(f,fp)
                            {
                                label pointi = f[nextFp];
                                newFacePts.append(pointi);
                                if (pointi == cFaces[0])
                                {
                                    break;
                                }
                                nextFp = f.fcIndex(nextFp);
                            }
                            face newFace(newFacePts);

                            label own = mesh.faceOwner()[faceI];

                            label zoneID = mesh.faceZones().whichZone(faceI);
                            bool zoneFlip = false;
                            if (zoneID >= 0)
                            {
                                const faceZone& fZone =
                                    mesh.faceZones()[zoneID];
                                zoneFlip =
                                    fZone.flipMap()[fZone.whichFace(faceI)];
                            }

                            // Modify the master face.
                            meshMod.setAction
                            (
                                polyModifyFace
                                (
                                    origFace,       // original face
                                    faceI,    // label of face
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
                                    newFace, // vertices
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

        {
            List<scalarField> cellWeights(cutCells.size());
            forAll(cellWeights, cwI)
            {
                cellWeights[cwI].setSize(cellLoops[cwI].size(), -GREAT);
            }

            cellCuts cuts(mesh, cutCells, cellLoops, cellAnchorPts, cellWeights);

            meshCutter mCut(mesh);

            polyTopoChange meshMod(mesh);
            mCut.setRefinement(cuts, meshMod);

            // Change the mesh (no inflation)
            autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

            // Update fields (problem when face added to zero sized patch)
            mesh.updateMesh(map);

            //Update layer cell field
            volScalarField& layerCells = const_cast<volScalarField&>
                (mesh.lookupObject<volScalarField>("layerStacks"));

            scalar maxLayerCount(gMax(layerCells));
            label nNewCells = map().cellMap().size() - map().nOldCells();

            globalIndex globalNewCells(nNewCells);
            label start  = maxLayerCount
                + globalNewCells.offset(Pstream::myProcNo());
            label count = 1;
            forAll(map().cellMap(), cellI)
            {
                if (cellI >= map().nOldCells())
                {
                    label oldCellI = map().cellMap()[cellI];

                    splitCells.append(labelPair(oldCellI, cellI));
                    if (oldCellI != -1 && layerCells[oldCellI] == -1)
                    {
                        continue;
                    }

                    layerCells[cellI] = start + count;
                    count++;
                }
            }

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

    return List<labelPair>(splitCells, true);
}


void Foam::autoSplitCells::remerge
(
    const List<labelPair>& splitCells
)
{
    Info<<"Checking if boundary cells need merging"<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    cellSet errorCells(mesh, "errorCells", mesh.nCells()/100+1);
    faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+1);

    label noFailedChecks(0);
    if (mesh.checkClosedCells(true, &errorCells)) noFailedChecks++;
    if (revertErrors_)
    {
        if (mesh.checkCellVolumes(true, &errorCells)) noFailedChecks++;
        if (mesh.checkFaceAreas(true, &errorFaces)) noFailedChecks++;
        if (mesh.checkFaceOrthogonality(true, &errorFaces)) noFailedChecks++;
        if (mesh.checkFacePyramids(true, -SMALL, &errorFaces)) noFailedChecks++;
    }

    if (noFailedChecks > 0)
    {
        polyTopoChange meshMod(mesh);
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        boolList errors(mesh.nCells(), false);
        boolList splitEdges(mesh.nEdges(), false);
        boolList checkedFaces(mesh.nFaces(), false);

        labelList keepCell(mesh.nCells(), -1);
        labelList removeCell(mesh.nCells(), -1);

        boolList boundaryEdges(mesh.nEdges(), false);
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (isValidPatch(patchI))
            {
                label startFaceI = pp.start();

                forAll(pp, i)
                {
                    label faceI = startFaceI+i;
                    const labelList& fEdges = mesh.faceEdges()[faceI];
                    forAll(fEdges, fEI)
                    {
                        boundaryEdges[fEdges[fEI]] = true;
                    }
                }
            }
        }
        syncTools::syncEdgeList
        (
            mesh,
            boundaryEdges,
            orEqOp<bool>(),
            false              // null value
        );

        forAllConstIter(labelHashSet, errorFaces, iter)
        {
            label own = mesh.faceOwner()[iter.key()];
            errors[own] = true;
            if (iter.key() < mesh.nInternalFaces())
            {
                label nei = mesh.faceNeighbour()[iter.key()];
                errors[nei] = true;
            }
        }

        forAllConstIter(labelHashSet, errorCells, iter)
        {
            errors[iter.key()] = true;
        }

        label nMerged = 0;
        forAll(splitCells, i)
        {
            label own = splitCells[i].first();
            label nei = splitCells[i].second();

            if (errors[own] || errors[nei])
            {
                //rejoin split cells
                nMerged++;
                meshMod.setAction(polyRemoveCell(own));

                keepCell[nei] = own;
                removeCell[own] = nei;
                const labelList& cFaces = mesh.cells()[own];

                label sFace = -1;
                forAll(cFaces, cFI)
                {
                    label faceI = cFaces[cFI];
                    if (mesh.isInternalFace(faceI))
                    {
                        if
                        (
                            mesh.faceOwner()[faceI] == nei
                            || mesh.faceNeighbour()[faceI] == nei
                        )
                        {
                            sFace = faceI;
                            break;
                        }
                    }
                }

                //remove shared faces
                meshMod.setAction(polyRemoveFace(sFace));
                checkedFaces[sFace] = true;

                //mark split edges
                const face f = mesh.faces()[sFace];
                forAll(f, fp)
                {
                    label meshPointI = f[fp];
                    label next = f.fcIndex(fp);

                    label meshEdgeI =  meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[meshPointI],
                        meshPointI,
                        f[next]
                    );

                    if (boundaryEdges[meshEdgeI])
                    {
                        labelHashSet cfOwn(mesh.cells()[own]);
                        labelHashSet cfNei(mesh.cells()[nei]);

                        const labelList& eFaces = mesh.edgeFaces()[meshEdgeI];
                        label nSplit = 0;
                        label nOther = 0;
                        forAll(eFaces, eFI)
                        {
                            label facei = eFaces[eFI];
                            if (facei == sFace)
                            {
                                continue;
                            }
                            const face& f = mesh.faces()[facei];
                            if (extruded_)
                            {
                                label cLevel = -1;

                                if (cfOwn.found(facei))
                                {
                                    cLevel = cellLevel[own];
                                }
                                else
                                {
                                    cLevel = cellLevel[nei];
                                }
                                if (cLevel != -1)
                                {
                                    label nAnchors = 0;
                                    forAll(f,fp)
                                    {
                                        label pointi = f[fp];
                                        label pLevel = pointLevel[pointi];
                                        if (pLevel <= cLevel)
                                        {
                                            nAnchors++;
                                        }
                                    }
                                    if (nAnchors == 3)
                                    {
                                        nSplit++;
                                    }
                                    else if (f.size() > 2)
                                    {
                                        nOther++;
                                    }
                                }
                            }
                            else
                            {
                                if
                                (
                                    (cfOwn.found(facei) || cfNei.found(facei))
                                    && f.size() == 3
                                 )
                                {
                                    nSplit++;
                                }
                            }
                        }

                        if
                        (
                            nSplit == 2
                            || (preMerged_ && (nSplit + nOther == 2) )
                        )
                        {

                            splitEdges[meshEdgeI] = true;
                            label oppositePt0 = f.fcIndex(next);
                            label oppositePt1 = f.fcIndex(oppositePt0);

                            meshEdgeI =  meshTools::findEdge
                            (
                                mesh.edges(),
                                mesh.pointEdges()[f[oppositePt0]],
                                f[oppositePt0],
                                f[oppositePt1]
                            );
                            splitEdges[meshEdgeI] = true;
                            break;
                        }
                    }
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            splitEdges,
            orEqOp<bool>(),
            false              // null value
        );

        boolList cuttingFace = checkedFaces;
        forAll(splitEdges, edgeI)
        {
            if (splitEdges[edgeI])
            {
                const edge e = mesh.edges()[edgeI];
                const labelList& eFaces = mesh.edgeFaces()[edgeI];
                DynamicList<label> mergeFaces(2);
                forAll(eFaces,eFI)
                {
                    label facei = eFaces[eFI];
                    if (!cuttingFace[facei])
                    {
                        mergeFaces.append(facei);
                    }
                }

                //check if first merge faces associated with kept
                //cell otherwise swap
                {
                    label faceI = mergeFaces[0];
                    bool swapAnchor(false);
                    label own = mesh.faceOwner()[faceI];

                    if (mesh.isInternalFace(faceI))
                    {
                        label nei = mesh.faceNeighbour()[faceI];

                        if (keepCell[own] == -1 && keepCell[nei] == -1)
                        {
                            swapAnchor = true;
                        }
                    }
                    else
                    {
                        if (keepCell[own] == -1)
                        {
                            swapAnchor = true;
                        }
                    }

                    if (swapAnchor)
                    {
                        label tmpFace = mergeFaces[0];
                        mergeFaces[0] = mergeFaces[1];
                        mergeFaces[1] = tmpFace;
                    }
                }

                face f0 = mesh.faces()[mergeFaces[0]];
                face f1 = mesh.faces()[mergeFaces[1]];
                face newFace(f0.size() + f1.size() -2);

                label start = findIndex(f0, e[0]);
                label end = findIndex(f0, e[1]);
                if (f0[f0.fcIndex(start)] == e[1])
                {
                    label oldStart = start;
                    start = end;
                    end = oldStart;
                }

                label next = start;
                forAll(f0, fp)
                {
                    newFace[fp] = f0[next];
                    next = f0.fcIndex(next);
                }

                label startPt = f0[start];
                label endPt = f0[end];

                start = findIndex(f1, endPt);
                bool reverse = false;
                if (f1[f1.fcIndex(start)] == startPt)
                {
                    reverse = true;
                }

                forAll(f1, fp)
                {
                    if (reverse)
                    {
                        start = f1.rcIndex(start);
                    }
                    else
                    {
                        start = f1.fcIndex(start);
                    }

                    if (f1[start] == startPt)
                    {
                        break;
                    }
                    else
                    {
                        newFace[f0.size()+fp] = f1[start];
                    }
                }

                meshMod.setAction(polyRemoveFace(mergeFaces[1]));
                checkedFaces[mergeFaces[0]] = true;
                checkedFaces[mergeFaces[1]] = true;
                label modFace = mergeFaces[0];


                label zoneID = mesh.faceZones().whichZone(modFace);
                bool zoneFlip = false;
                if (zoneID >= 0)
                {
                    const faceZone& fZone = mesh.faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(modFace)];
                }

                label own = mesh.faceOwner()[modFace];

                if (removeCell[own] != -1)
                {
                    WarningInFunction
                        <<"Split edge merge, Owner cell: "<< own
                        <<" has been removed from mesh "
                        <<"kept cell is "<< removeCell[own]<<endl;
                }

                if (mesh.isInternalFace(modFace))
                {
                    label nei = mesh.faceNeighbour()[modFace];

                    if (removeCell[nei] != -1)
                    {
                        WarningInFunction
                            <<"Split edge merge, Neighbour cell: "<< nei
                            <<" has been removed from mesh "
                            <<" kept cell is "<< removeCell[nei]<<endl;
                    }

                    // Modify the master face.
                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            newFace,       // original face
                            modFace,    // label of face
                            own,            // owner
                            nei,             // neighbour
                            false,          // face flip
                            -1,         // patch for face
                            false,          // remove from zone
                            zoneID,         // zone for face
                            zoneFlip        // face flip in zone
                         )
                    );
                }
                else
                {
                    label patchI = mesh.boundaryMesh().whichPatch(modFace);

                    // Modify the master face.
                    meshMod.setAction
                    (
                        polyModifyFace
                        (
                            newFace,       // original face
                            modFace,    // label of face
                            own,            // owner
                            -1,             // neighbour
                            false,          // face flip
                            patchI,         // patch for face
                            false,          // remove from zone
                            zoneID,         // zone for face
                            zoneFlip        // face flip in zone
                         )
                    );
                }
            }
        }

        forAll(keepCell, cellI)
        {
            if (keepCell[cellI] != -1)
            {
                label removedCell = keepCell[cellI];
                const labelList cFaces = mesh.cells()[removedCell];
                forAll(cFaces, cFI)
                {
                    label faceI = cFaces[cFI];
                    if (!checkedFaces[faceI])
                    {
                        checkedFaces[faceI] = true;
                        label zoneID = mesh.faceZones().whichZone(faceI);
                        bool zoneFlip = false;
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone =
                                mesh.faceZones()[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
                        }

                        if (mesh.isInternalFace(faceI))
                        {
                            label own = cellI;
                            label nei = -1;
                            bool flipFace = false;
                            if (mesh.faceOwner()[faceI] == removedCell)
                            {
                                nei = mesh.faceNeighbour()[faceI];

                                if (removeCell[nei] != -1)
                                {
                                    nei = removeCell[nei];
                                }

                                if (own > nei)
                                {
                                    own = nei;
                                    nei = cellI;
                                    flipFace = true;
                                }
                            }
                            else
                            {
                                nei = mesh.faceOwner()[faceI];

                                if (removeCell[nei] != -1)
                                {
                                    nei = removeCell[nei];
                                }

                                if (own > nei)
                                {
                                    own = nei;
                                    nei = cellI;
                                }
                                else
                                {
                                    flipFace = true;
                                }
                            }

                            face newFace(mesh.faces()[faceI]);
                            if (flipFace)
                            {
                                newFace.flip();
                            }

                            if (removeCell[own] != -1)
                            {
                                WarningInFunction
                                    <<"Deleted cell faces update, Owner cell: "
                                    << own << " has been removed from mesh "
                                    <<" kept cell is "<< removeCell[own]<<endl;
                            }

                            if (removeCell[nei] != -1)
                            {
                                WarningInFunction
                                    <<"Deleted cell faces update, "
                                    <<"Neighbour cell: "<< nei
                                    <<" has been removed from mesh "
                                    <<" kept cell is "<< removeCell[nei]<<endl;
                            }

                            // Modify the master face.
                            meshMod.setAction
                            (
                                polyModifyFace
                                (
                                    newFace,       // original face
                                    faceI,    // label of face
                                    own,            // owner
                                    nei,             // neighbour
                                    false,          // face flip
                                    -1,         // patch for face
                                    false,          // remove from zone
                                    zoneID,         // zone for face
                                    zoneFlip        // face flip in zone
                                )
                             );
                        }
                        else
                        {
                            label patchI =
                                mesh.boundaryMesh().whichPatch(faceI);

                            // Modify the master face.
                            meshMod.setAction
                            (
                                polyModifyFace
                                (
                                    mesh.faces()[faceI], // original face
                                    faceI,    // label of face
                                    cellI,            // owner
                                    -1,             // neighbour
                                    false,          // face flip
                                    patchI,         // patch for face
                                    false,          // remove from zone
                                    zoneID,         // zone for face
                                    zoneFlip        // face flip in zone
                                 )
                            );
                        }
                    }
                }
            }
        }


        if (returnReduce(nMerged, sumOp<label>()) > 0)
        {
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
}


Foam::List<Foam::Tuple2<Foam::label,Foam::labelPair>>
Foam::autoSplitCells::warpedFacesToSplit
(
    const scalar warpage,
    const scalar fchwarpage,
    const List<wordList>& layerSpec,
    const labelList& patchNumLayers,
    const scalarField& patchFCH,
    const scalarField& patchExpansionRatio,
    const bool curvatureSplit,
    const bool mergedFaces,
    const bool fourAnchors
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const scalar edge0Len = meshRefiner_.meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    DynamicList<Tuple2<label,labelPair>> cutFaces(mesh.nBoundaryFaces());

    DynamicList<label> allPatches(patches.size());
    forAll(patches, patchi)
    {
        if (isValidPatch(patchi))
        {
            allPatches.append(patchi);
        }
    }

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            allPatches
        )
    );
    indirectPrimitivePatch& pp = ppPtr();

    vectorField fNormals
    (
        pp.localPoints().size(),
        vector(GREAT, GREAT, GREAT)
    );
    List<pointIndexHit> hitInfo;

    if (curvatureSplit)
    {
        const labelList allSurfaces
        (
            identity(meshRefiner_.surfaces().surfaces().size())
        );

        labelList hitSurf;
        meshRefiner_.surfaces().findNearest
        (
            allSurfaces,
            pp.localPoints(),
            scalarField(pp.localPoints().size(), GREAT),
            hitSurf,
            hitInfo
        );

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
    }
    scalar minSplitRatio = 0.15;
    scalar alignRatio = 0.1736;
    const labelList& meshPoints = pp.meshPoints();
    forAll(pp, i)
    {
        label facei = pp.addressing()[i];
        label own = mesh.faceOwner()[facei];
        label cLevel = cellLevel[own];
        const face& lf = pp.localFaces()[i];
        DynamicList<label> anchorPts(pp[i].size());
        {
            forAll(lf,fp)
            {
                label lfp = lf[fp];
                label meshpointi = meshPoints[lfp];
                label pLevel = pointLevel[meshpointi];

                if (pLevel <= cLevel)
                {
                    anchorPts.append(lfp);
                }
            }

            if (!fourAnchors)
            {
                if (anchorPts.size() == 1 && lf.size() == 4)
                {
                    anchorPts.clear();
                    forAll(lf,fp)
                    {
                        anchorPts.append(lf[fp]);
                    }
                }
                else if (anchorPts.size() == 2 && lf.size() == 5)
                {
                    label anchor0 = anchorPts[0];
                    label anchor1 = anchorPts[1];

                    label startAnchor = -1;
                    forAll(lf,fp)
                    {
                        label lfp = lf[fp];
                        label nextFp = lf.fcIndex(fp);
                        label nextPt = lf[nextFp];
                        if (anchor0 == lfp && anchor1 == nextPt)
                        {
                            startAnchor = fp;
                        }
                        if (anchor1 == lfp && anchor0 == nextPt)
                        {
                            startAnchor = fp;
                        }
                    }

                    if (startAnchor != -1)
                    {
                        anchorPts.clear();
                        label fp = lf.rcIndex(startAnchor);
                        anchorPts.append(lf[fp]);
                        fp = lf.fcIndex(fp);
                        anchorPts.append(lf[fp]);
                        fp = lf.fcIndex(fp);
                        anchorPts.append(lf[fp]);
                        fp = lf.fcIndex(fp);
                        anchorPts.append(lf[fp]);
                    }
                }
            }
        }

        // Normal distance to face centre plane
        const point& fc = mesh.faceCentres()[facei];
        const vector& fn = pp.faceNormals()[i];

        scalarField vProj(lf.size());
        forAll(lf,fp)
        {
            label lfp = lf[fp];
            label meshpointi = meshPoints[lfp];
            vector n = mesh.points()[meshpointi] - fc;
            vProj[fp] = (n & fn);
        }

        // Get normal 'span' of face
        scalar minVal = min(vProj);
        scalar maxVal = max(vProj);

        bool warpedFace = false;

        scalar dW = (maxVal - minVal);

        if (warpage > 0 && dW > warpage)
        {
            warpedFace = true;
        }
        else if (fchwarpage > 0)
        {
            label patchi = patches.whichPatch(facei);
            word method = layerSpec[patchi][1];
            if (method == "fch")
            {
                if (dW > fchwarpage*patchFCH[patchi])
                {
                    warpedFace = true;
                }
            }
            else if (method == "expansionRatio")
            {
                label nLayers = patchNumLayers[patchi];
                if (nLayers > 0)
                {
                    scalar eR =  patchExpansionRatio[patchi];
                    scalar len = edge0Len / pow(2., cLevel);
                    scalar estFCH = len / pow(eR, nLayers);
                    if (dW > fchwarpage*estFCH)
                    {
                        warpedFace = true;
                    }
                }
            }
        }

        if (warpedFace)
        {
            if (anchorPts.size() != 4 && lf.size() > 3)
            {
                scalar maxSplitRatio = -1;
                edge cutEdge(-1, -1);
                for (int fpi = 0; fpi < lf.size()-2; fpi++)
                {
                    label oppPt = lf.fcIndex(fpi);
                    if (oppPt == 0)
                    {
                        break;
                    }

                    forAll(lf, fp)
                    {
                        oppPt = lf.fcIndex(oppPt);
                        if (oppPt == 0)
                        {
                            break;
                        }
                        DynamicList<label> fpts0;
                        DynamicList<label> fpts1;
                        label startPt = fpi;
                        forAll(lf,fp)
                        {
                            fpts0.append(meshPoints[lf[startPt]]);
                            if (startPt == oppPt)
                            {
                                break;
                            }
                            startPt = lf.fcIndex(startPt);
                        }
                        startPt = oppPt;
                        forAll(lf,fp)
                        {
                            fpts1.append(meshPoints[lf[startPt]]);
                            if (startPt == fpi)
                            {
                                break;
                            }
                            startPt = lf.fcIndex(startPt);
                        }

                        face cf0(fpts0);
                        face cf1(fpts1);

                        vector fN0 = cf0.areaNormal(mesh.points());
                        vector fN1 = cf1.areaNormal(mesh.points());
                        scalar fA0 = mag(fN0);
                        scalar fA1 = mag(fN1);

                        if (fA0 > SMALL && fA1 > SMALL)
                        {
                            scalar splitRatio = min(fA0,fA1) / max(fA0,fA1);
                            if (splitRatio > 0.25)
                            {
                                edge candidateEdge
                                (
                                    meshPoints[lf[fpi]],
                                    meshPoints[lf[oppPt]]
                                );
                                bool splitOK = true;
                                if (mergedFaces)
                                {
                                    DynamicList<label> topEdge(2);
                                    forAll(candidateEdge, cPt)
                                    {
                                        label pointi = candidateEdge[cPt];
                                        const labelList& pedges =
                                            mesh.pointEdges()[pointi];
                                        labelHashSet cEdgeSet
                                        (
                                            mesh.cellEdges()[own]
                                        );
                                        labelHashSet fEdgeSet
                                        (
                                            mesh.faceEdges()[facei]
                                        );

                                        forAll(pedges, pei)
                                        {
                                            label meshedgei = pedges[pei];
                                            if
                                            (
                                                cEdgeSet.found(meshedgei)
                                                && !fEdgeSet.found(meshedgei)
                                            )
                                            {
                                                const edge& e =
                                                    mesh.edges()[meshedgei];
                                                label otherPt
                                                (
                                                    e[0] == pointi ?
                                                    e[1] : e[0]
                                                );
                                                topEdge.append(otherPt);
                                                break;
                                            }
                                        }
                                    }
                                    std::list<label> intx;
                                    if (topEdge.size() == 2)
                                    {
                                        labelList faces1 =
                                            mesh.pointFaces()[topEdge[0]];
                                        labelList faces2 =
                                            mesh.pointFaces()[topEdge[1]];
                                        std::set_intersection
                                        (
                                            faces1.begin(),
                                            faces1.end(),
                                            faces2.begin(),
                                            faces2.end(),
                                            std::back_inserter(intx)
                                         );
                                    }
                                    if (intx.size() != 1)
                                    {
                                        splitOK = false;
                                    }
                                }

                                if (splitOK && splitRatio > maxSplitRatio)
                                {
                                    maxSplitRatio = splitRatio;
                                    cutEdge = candidateEdge;
                                }
                            }
                        }
                    }
                }

                if (cutEdge[0] != -1 && cutEdge[1] != -1)
                {
                    cutFaces.append
                    (
                        Tuple2<label,labelPair>
                        (
                            facei,
                            cutEdge
                         )
                     );
                }
            }
            else if (anchorPts.size() == 4)
            {
                face f(anchorPts);
                if (!curvatureSplit)
                {
                    for (label ci = 0; ci < 2; ci++)
                    {
                        edge cutEdge(meshPoints[f[ci]],meshPoints[f[ci+2]]);

                        label start0 = ci;
                        label start1 = ci+2;

                        face cf0(3);
                        face cf1(3);

                        forAll(cf0, j)
                        {
                            cf0[j] = meshPoints[f[start0]];
                            cf1[j] =meshPoints[f[start1]];
                            start0 = f.fcIndex(start0);
                            start1 = f.fcIndex(start1);
                        }

                        vector fN0 = cf0.areaNormal(mesh.points());
                        vector fN1 = cf1.areaNormal(mesh.points());

                        scalar fA0 = mag(fN0);
                        scalar fA1 = mag(fN1);

                        vector fC0 = cf0.centre(mesh.points());

                        if (fA0 > SMALL && fA1 > SMALL)
                        {
                            fN0 /= fA0;
                            fN1 /= fA1;

                            vector eV = 0.5*(fN0+fN1);
                            point cutEdgeCentre = cutEdge.centre(mesh.points());
                            vector ecTofc = cutEdgeCentre - fC0;
                            if ((ecTofc & eV) > 0)
                            {
                                scalar splitRatio = min(fA0,fA1) / max(fA0,fA1);
                                if (splitRatio > 0.25 && (fN0&fN1) > alignRatio)
                                {
                                    cutFaces.append
                                    (
                                        Tuple2<label,labelPair>
                                        (
                                            facei,
                                            cutEdge
                                        )
                                    );
                                    break;
                                }
                            }
                        }
                    }
                }
                else
                {
                    label nHits = 0;
                    forAll(f, fp)
                    {
                        if (hitInfo[f[fp]].hit())
                        {
                            nHits++;
                        }
                    }
                    if (nHits == f.size())
                    {
                        vector pN0
                        (
                            (fn & fNormals[f[0]]) > 0
                            ? fNormals[f[0]]
                            : -fNormals[f[0]]
                         );

                        vector pN1
                        (
                            (fn & fNormals[f[1]]) > 0
                            ? fNormals[f[1]]
                            : -fNormals[f[1]]
                        );

                        vector pN2
                        (
                            (fn & fNormals[f[2]]) > 0
                            ? fNormals[f[2]]
                            : -fNormals[f[2]]
                        );

                        vector pN3
                        (
                            (fn & fNormals[f[3]]) > 0
                            ? fNormals[f[3]]
                            : -fNormals[f[3]]
                        );

                        scalar dP0 = (pN0 & pN2);
                        scalar dP1 = (pN1 & pN3);

                        scalar splitRatio02;
                        scalar splitRatio13;
                        edge cutEdge02;
                        edge cutEdge13;
                        bool validAlign02 = true;
                        bool validAlign13 = true;

                        {
                            face fA(3);
                            face fB(3);
                            label nextA = 0;
                            label nextB = 2;
                            forAll(fA,fp)
                            {
                                fA[fp] = meshPoints[f[nextA]];
                                fB[fp] = meshPoints[f[nextB]];
                                nextA = f.fcIndex(nextA);
                                nextB = f.fcIndex(nextB);
                            }
                            vector fAA  = fA.areaNormal(mesh.points());
                            vector fAB  = fB.areaNormal(mesh.points());
                            scalar areaA = mag(fAA)+SMALL;
                            scalar areaB = mag(fAB)+SMALL;

                            fAA /= areaA;
                            fAB /= areaB;

                            cutEdge02[0] = meshPoints[f[0]];
                            cutEdge02[1] = meshPoints[f[2]];
                            point cutEdgeCentre = cutEdge02.centre(mesh.points());
                            vector fCA = fA.centre(mesh.points());
                            vector ecTofc = cutEdgeCentre - fCA;
                            vector eV = 0.5*(fAA+fAB);
                            if ((ecTofc & eV) > 0 && (fAA&&fAB) < alignRatio)
                            {
                                validAlign02 = false;
                            }

                            splitRatio02 = min(areaA,areaB) / max(areaA,areaB);

                            nextA = 1;
                            nextB = 3;
                            forAll(fA,fp)
                            {
                                fA[fp] = meshPoints[f[nextA]];
                                fB[fp] = meshPoints[f[nextB]];
                                nextA = f.fcIndex(nextA);
                                nextB = f.fcIndex(nextB);
                            }

                            fAA  = fA.areaNormal(mesh.points());
                            fAB  = fB.areaNormal(mesh.points());
                            areaA = mag(fAA)+SMALL;
                            areaB = mag(fAB)+SMALL;

                            fAA /= areaA;
                            fAB /= areaB;

                            cutEdge13[0] = meshPoints[f[1]];
                            cutEdge13[1] = meshPoints[f[3]];
                            cutEdgeCentre = cutEdge13.centre(mesh.points());
                            fCA = fA.centre(mesh.points());
                            ecTofc = cutEdgeCentre - fCA;
                            eV = 0.5*(fAA+fAB);

                            if ((ecTofc & eV) > 0 && (fAA&&fAB) < alignRatio)
                            {
                                validAlign13 = false;
                            }
                            splitRatio13 = min(areaA,areaB) / max(areaA,areaB);
                        }

                        scalar maxSplitRatio = max(splitRatio13,splitRatio02);
                        if
                        (
                            maxSplitRatio <= minSplitRatio
                            || (!validAlign02 && !validAlign13)
                        )
                        {
                            continue;
                        }

                        edge cutEdge;
                        if (dP0 > dP1)
                        {
                            if
                            (
                                splitRatio02  > 0.75*maxSplitRatio
                                && splitRatio02 > minSplitRatio
                                && validAlign02
                            )
                            {
                                cutEdge = cutEdge02;
                            }
                            else
                            {
                                cutEdge = cutEdge13;
                            }
                        }
                        else
                        {
                            if
                            (
                                splitRatio13  > 0.75*maxSplitRatio
                                && splitRatio13 > minSplitRatio
                                && validAlign13
                            )
                            {
                                cutEdge = cutEdge13;
                            }
                            else
                            {
                                cutEdge = cutEdge02;
                            }
                        }

                        cutFaces.append
                        (
                            Tuple2<label,labelPair>
                            (
                                facei,
                                cutEdge
                             )
                         );
                    }
                }
            }
        }
    }

    return List<Tuple2<label,labelPair>>(cutFaces, true);
}


Foam::List<Foam::Tuple2<Foam::label,Foam::labelPair>>
Foam::autoSplitCells::featureFacesToSplit
(
    const scalar featureAngle
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    DynamicList<label> allPatches(patches.size());
    forAll(patches, patchi)
    {
        if (isValidPatch(patchi))
        {
            allPatches.append(patchi);
        }
    }

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            allPatches
        )
    );

    const labelList ppMeshEdges
    (
        ppPtr().meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    boolList excludedFaces(ppPtr().size(), false);
    edgeClassification eClass
    (
        mesh,
        mesh.points(),
        ppPtr(),
        ppMeshEdges,
        excludedFaces,
        0.8191,
        0.8191
    );
    const List<Tuple2<edgeClassification::edgeType,scalar>>&
        eType = eClass.edgeTypes();

    boolList featureEdges(mesh.nEdges(), false);
    forAll(ppMeshEdges, edgeI)
    {
        if
        (
            eType[edgeI].first() == edgeClassification::CONCAVE
            || eType[edgeI].first() == edgeClassification::CONVEX
        )
        {
            label meshEdgeI = ppMeshEdges[edgeI];
            featureEdges[meshEdgeI] = true;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        featureEdges,
        orEqOp<bool>(),
        false
    );

    DynamicList<Tuple2<label,labelPair>> cutFaces(mesh.nBoundaryFaces());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isValidPatch(patchI))
        {
            label startfacei = pp.start();

            forAll(pp, i)
            {
                label facei = startfacei+i;
                DynamicList<label> anchorPts(mesh.faces()[facei].size());
                {
                    face f = mesh.faces()[facei];
                    label celli = mesh.faceOwner()[facei];

                    label cLevel = cellLevel[celli];

                    forAll(f,fp)
                    {
                        label pointi = f[fp];
                        if (extruded_)
                        {
                            label pLevel = pointLevel[pointi];
                            if (pLevel <= cLevel)
                            {
                                anchorPts.append(pointi);
                            }
                        }
                        else
                        {
                            anchorPts.append(pointi);
                        }
                    }
                }

                face f(anchorPts);
                if (f.size() != 4)
                {
                    continue;
                }

                face mf = mesh.faces()[facei];
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    point pt = mesh.points()[pointi];

                    label nextFp = f.fcIndex(fp);
                    label prevFp = f.rcIndex(fp);

                    vector nextVec = mesh.points()[f[nextFp]] - pt;
                    vector prevVec = pt - mesh.points()[f[prevFp]];
                    nextVec /= Foam::mag(nextVec) + VSMALL;
                    prevVec /= Foam::mag(prevVec) + VSMALL;

                    label currentPt = findIndex(mf, pointi);
                    bool nextFE = true;
                    while (true)
                    {
                        label origNextFp = mf.fcIndex(currentPt);
                        label nextEdge = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[mf[currentPt]],
                            mf[currentPt],
                            mf[origNextFp]
                        );
                        if (!featureEdges[nextEdge])
                        {
                            nextFE = false;
                            break;
                        }
                        else if (mf[origNextFp] == f[nextFp])
                        {
                            break;
                        }
                        currentPt = origNextFp;
                    }

                    currentPt = findIndex(mf, pointi);
                    bool prevFE = true;
                    while (true)
                    {
                        label origPrevFp = mf.rcIndex(currentPt);
                        label prevEdge = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[mf[currentPt]],
                            mf[currentPt],
                            mf[origPrevFp]
                        );
                        if (!featureEdges[prevEdge])
                        {
                            prevFE = false;
                            break;
                        }
                        else if (mf[origPrevFp] == f[prevFp])
                        {
                            break;
                        }
                        currentPt = origPrevFp;
                    }

                    if
                    (
                        prevFE && nextFE
                        && (nextVec & prevVec) > featureAngle
                    )
                    {
                        label oppositePt = f[f.fcIndex(nextFp)];
                        cutFaces.append
                        (
                            Tuple2<label,labelPair>
                            (
                                facei,
                                labelPair(pointi,oppositePt)
                            )
                        );
                        break;
                    }
                }
            }
        }
    }

    return List<Tuple2<label,labelPair>>(cutFaces, true);
}


void Foam::autoSplitCells::splitWarpedCells
(
    const scalar warpage,
    const scalar fchwarpage,
    const List<wordList>& layerSpec,
    const labelList& patchNumLayers,
    const scalarField& patchFCH,
    const scalarField& patchExpansionRatio,
    const bool curvatureSplit,
    const bool convexSplit,
    const bool mergedFaces
)
{
    Info<< "Splitting warped cells" << endl;

    if (extruded_)
    {
        //First flatten faces at refinement interfaces (hanging nodes)
        flattenRefinementFaces();
    }

    //first try cutting any warped faces
    {
        List<Tuple2<label,labelPair>> cutFaces = warpedFacesToSplit
        (
            warpage,
            fchwarpage,
            layerSpec,
            patchNumLayers,
            patchFCH,
            patchExpansionRatio,
            curvatureSplit,
            mergedFaces,
            false
        );
        //Check face cuts are valid
        filterCutFaces(cutFaces);
        update(cutFaces,false);
    }

    //Move hanging points agin and try re-cutting
    if (extruded_)
    {
        moveSplitFaceHangingPts();
        List<Tuple2<label,labelPair>> cutFaces = warpedFacesToSplit
        (
            warpage,
            fchwarpage,
            layerSpec,
            patchNumLayers,
            patchFCH,
            patchExpansionRatio,
            curvatureSplit,
            mergedFaces,
            false
        );
        //Check face cuts are valid
        filterCutFaces(cutFaces);
        update(cutFaces,false);
    }

    if (curvatureSplit && convexSplit)
    {
        List<Tuple2<label,labelPair>> cutFaces = warpedFacesToSplit
        (
            warpage,
            fchwarpage,
            layerSpec,
            patchNumLayers,
            patchFCH,
            patchExpansionRatio,
            false,
            mergedFaces,
            true
        );
        //Check face cuts are valid
        filterCutFaces(cutFaces);
        update(cutFaces,false);
    }
}


void Foam::autoSplitCells::splitFeatureCells
(
    const scalar featureAngle,
    const bool moveFeaturePts,
    const bool moveConvexPts
)
{
    Info<< "Splitting feature cells" << endl;
    List<Tuple2<label,labelPair>> cutFaces = featureFacesToSplit(featureAngle);

    //Check face cuts are valid
    filterCutFaces(cutFaces);

    //Move cells to help cutting
    moveFeaturePoints
    (
        cutFaces,
        featureAngle,
        moveFeaturePts,
        moveConvexPts
    );

    update(cutFaces,false);
}


void Foam::autoSplitCells::splitPreCalculated
(
    List<Tuple2<label,labelPair>>& cutFaces
)
{
    Info<< "Splitting cells using pre-defined cuts" << endl;

    //Check face cuts are valid
    filterCutFaces
    (
        cutFaces,
        true //prevent cutting cells attached to face zones
    );

    update(cutFaces,true);
}

//Check for invalid cuts and remove
void Foam::autoSplitCells::filterCutFaces
(
    List<Tuple2<label,labelPair>>& cutFaces,
    const bool preventZoneCuts
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    boolList blockCut(mesh.nCells(), false);
    if (preventZoneCuts)
    {
        PackedList<1> isZonedFace
        (
            snappySnapDriver::setZonedFaces(meshRefiner_)
        );
        forAll(mesh.cells(), celli)
        {
            const cell& c = mesh.cells()[celli];
            forAll(c,cFI)
            {
                label facei = c[cFI];
                if (isZonedFace.get(facei))
                {
                    blockCut[celli] = true;
                    break;
                }
            }
        }
    }

    label nCorrect = 0;
    forAll(cutFaces, i)
    {
        const label facei = cutFaces[i].first();
        const edge cutEdge = cutFaces[i].second();
        const label own = mesh.faceOwner()[facei];

        if (blockCut[own])
        {
            continue;
        }
        else
        {
            blockCut[own] = true;
        }

        const label cLevel = cellLevel[own];
        const cell& c = mesh.cells()[own];
        const face& f = mesh.faces()[facei];
        if (extruded_)
        {
            if (f.size() < 4)
            {
                continue;
            }

            if ((2*f.size()) <= mesh.cellPoints()[own].size())
            {
                bool validNonAnchor = false;
                if (preMerged_)
                {
                    label nFaceAnchors = 0;
                    forAll(f, fp)
                    {
                        label pointi = f[fp];
                        if (pointLevel[pointi] <= cLevel)
                        {
                            nFaceAnchors++;
                        }
                    }

                    if (nFaceAnchors == 2 && f.size() == 5)
                    {
                        validNonAnchor = true;
                    }
                    else if (nFaceAnchors == 1 && f.size() == 4)
                    {
                        validNonAnchor = true;
                    }
                    else if (nFaceAnchors != 4)
                    {
                        validNonAnchor = true;
                    }
                }

                if (validNonAnchor)
                {
                    label meshEdgei = meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[cutEdge[0]],
                        cutEdge[0],
                        cutEdge[1]
                     );
                    if (meshEdgei == -1)
                    {
                        cutFaces[nCorrect] = cutFaces[i];
                        nCorrect++;
                    }
                }
                else if
                (
                    pointLevel[cutEdge[0]] <= cLevel
                    && pointLevel[cutEdge[1]] <= cLevel
                )
                {
                    label currentPt = findIndex(f, cutEdge[0]);
                    label nAnchors = 0;
                    forAll(f, fp)
                    {
                        label nextPt = f.fcIndex(currentPt);
                        label pointi = f[nextPt];

                        if (pointLevel[pointi] <= cLevel)
                        {
                            nAnchors++;
                            if (pointi == cutEdge[1])
                            {
                                break;
                            }
                        }
                        currentPt = nextPt;
                    }

                    if (nAnchors == 2)
                    {
                        cutFaces[nCorrect] = cutFaces[i];
                        nCorrect++;
                    }
                }
            }
        }
        else if
        (
            f.size() == 4 && c.size() == 6
            && mesh.cellPoints()[own].size() == 8
        )
        {
            label meshEdgei = meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[cutEdge[0]],
                cutEdge[0],
                cutEdge[1]
            );

            if (meshEdgei == -1)
            {
                cutFaces[nCorrect] = cutFaces[i];
                nCorrect++;
            }
        }
    }
    cutFaces.setSize(nCorrect);
}


void Foam::autoSplitCells::update
(
    List<Tuple2<label,labelPair>>& cutFaces,
    bool checkCutAreas
)
{
    //Split the cells based on list of face cuts
    List<labelPair> splitCells = split(cutFaces,checkCutAreas);

    if (returnReduce(splitCells.size(), sumOp<label>()) > 0)
    {
        //Remerge split cells if errors introduced
        remerge(splitCells);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::autoSplitCells::autoSplitCells
(
    meshRefinement& meshRefiner,
    labelList excludePatches,
    bool extruded,
    bool preMerged,
    bool revertErrors
)
:
    meshRefiner_(meshRefiner),
    excludePatchSet_(excludePatches),
    extruded_(extruded),
    preMerged_(preMerged),
    revertErrors_(revertErrors)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::autoSplitCells::~autoSplitCells()
{}

// ************************************************************************* //
