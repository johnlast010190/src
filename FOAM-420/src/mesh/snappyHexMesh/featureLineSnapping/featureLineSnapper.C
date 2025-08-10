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
    (c) 1991-2008 OpenCFD Ltd.

\*----------------------------------------------------------------------------*/

//include this first because of VDB include conflicts down the line
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "featureLineSnapping/featureLineSnapper.H"
#include "snappyHexMeshDriver/snappySnapDriver.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "fvMesh/fvMesh.H"
#include "containers/Lists/List/List.H"
#include "containers/Lists/DynamicList/DynamicList.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "containers/HashTables/Map/Map.H"
#include "meshTools/meshTools.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "indexedOctree/treeDataFace.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "primitives/ops/ops.H"
#include "meshRefinement/meshRefinement.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "primitives/Pair/labelPair.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "db/dictionary/dictionary.H"
#include "edgeMesh/featureEdgeMesh/featureEdgeMesh.H"
#include "sets/topoSets/faceSet.H"
#include "motionSmoother/motionSmoother.H"
#include <list>
#include <algorithm>


Foam::featureLineSnapper::featureLineSnapper
(
    const featureLinePrep &flPrep,
    const dictionary &snapDict,
    const dictionary &motionDict,
    const labelList &adaptPatchIDs,
    indirectPrimitivePatch &pp,
    meshRefinement &meshRefiner,
    List<labelPair>& baffles,
    bool dualMesh,
    bool allowQuadSplits
):
    flPrep_(flPrep),
    motionDict_(motionDict),
    adaptPatchIDs_(adaptPatchIDs),
    pp_(pp),
    meshRefiner_(meshRefiner),
    baffles_(baffles),
    dualMesh_(dualMesh),
    allowQuadSplits_(allowQuadSplits),
    faceIsSplit(pp.size(),0),
    maxShrinkFactor_(0.01),
    minFaceAngle_(0.17453),
    minDihedralForSplits_(-1.5),
    verbose_(false),
    dump_(false),
    nNoValidEndPts_(0),
    nFailedConstr_(0)
{
    if (snapDict.found("verboseFeatureLines"))
    {
        verbose_ = readBool(snapDict.lookup("verboseFeatureLines"));
    }
    if (snapDict.found("dumpFeatureLineData"))
    {
        dump_ = readBool(snapDict.lookup("dumpFeatureLineData"));
    }
}

/**
 * Controls building and clean-up of shadow chains from the supplied
 * set of feature lines
 */
void Foam::featureLineSnapper::constructShadowChains
(
    const PackedList<1>& isZonedFace
)
{
    Info<< "Finding optimal points for direct feature snapping" << endl;
    Info<< "--------------------------------------------------" << endl
         << endl;

    if (dump_)
    {
        simpleVTKWriter(pp_.localFaces(),pp_.localPoints()).write("preChain.vtk");
    }

    if (dump_)
    {
        for (label flI = 0; flI < flPrep_.getNumFL(); ++flI)
        {
            simpleVTKWriter writer;
            const featureLine &featurel = flPrep_.getFeatureLine(flI);
            writer.addMorePoints(featurel,featurel.points());
            writer.write("allFeatureLines"+Foam::name(flI)+".vtk");
        }
    }

    pointDisp_.setSize(pp_.localPoints().size(),vector::zero);

    Info<< " Building node chains" << endl;

    buildShadowChains();

    if (dump_)
    {
        simpleVTKWriter writer;
        forAll(shadowChains_,scI)
        {
            writer.addMorePoints(shadowChains_[scI],pp_.localPoints());
        }
        writer.write("preFilteringChains.vtk");
    }
    Info<< " Filtering out duplicates" << endl;

    filterChains();

    if (dump_)
    {
        simpleVTKWriter writer;
        forAll(shadowChains_,scI)
        {
            writer.addMorePoints(shadowChains_[scI],pp_.localPoints());
        }
        writer.write("postFilteringChains.vtk");
    }

    chainDispToPatchDisp();

    Info<< " Cleaning up results" << endl;

    cleanUpChains();

    syncTools::syncPointList
    (
        meshRefiner_.mesh(),
        pp_.meshPoints(),
        pointDisp_,
        maxMagEqOp(),
        vector::zero
    );

    filterDisplacementField();

    syncTools::syncPointList
    (
        meshRefiner_.mesh(),
        pp_.meshPoints(),
        pointDisp_,
        maxMagEqOp(),
        vector::zero
    );

    if (!dualMesh_)
    {
        //If using standard mesher split boundary faces if region
        // splitting has been selected
        addChainDiagonalEdges(isZonedFace);
    }

    return;
}


/**
 * Return list of neighbouring points
*/
Foam::labelList Foam::featureLineSnapper::neigbouringPoints
(
    const label curPt,
    const bool allowSplits,
    const bool includeCurrentPt
) const
{
    labelList  nodeMfd;
    {
        if (allowSplits)
        {
            //Add all nodes excl/incl curPt lying in faces containing curPt.
            labelList faceMfd = pp_.pointFaces()[curPt];
            DynamicList<label> rawNodes;
            if (includeCurrentPt)
            {
                rawNodes.append(curPt);

            }
            forAll(faceMfd,fI)
            {
                const face &face = pp_.localFaces()[faceMfd[fI]];
                forAll(face,pI)
                {
                    if (face[pI] != curPt)
                    {
                        rawNodes.append(face[pI]);
                    }
                }
            }
            //sort-and-unique the list
            nodeMfd.transfer(rawNodes);
            sort(nodeMfd);
            nodeMfd.setSize(label(std::unique(nodeMfd.begin(),nodeMfd.end())
                                  - nodeMfd.begin()));
        }
        else
        {
            if (includeCurrentPt)
            {
                nodeMfd.setSize(pp_.pointEdges()[curPt].size()+1);
                nodeMfd[0] = curPt;
            }
            else
            {
                nodeMfd.setSize(pp_.pointEdges()[curPt].size());
            }

            label offset = includeCurrentPt ? 1 : 0;

            forAll(nodeMfd,j)
            {
                nodeMfd[j+offset] = pp_.edges()[pp_.pointEdges()[curPt][j]]
                    .otherVertex(curPt);
            }
        }
    }
    return nodeMfd;
}


/**
 * Loops over the feature lines constructing closest patch points to fl
 * end points
*/
void Foam::featureLineSnapper::findEndPts
(
    const treeBoundBox& bb,
    const featureLine& featurel,
    const List<scalar>& flParams,
    const indexedOctree<treeDataEdge>& tree,
    const List<scalar>& lengthAve,
    label& minStartPt,
    label& minEndPt
)
{
    //Find the closest patch points to the fl end points
    //--------------------------------------------------
    const bool lineIsLoop = featurel[0] ==
        featurel[featurel.size()-1];

    scalar scaleFactor = 2.0; //magic number, anyone?
    //Exclude points further from the end point than 200% of the feature
    //line's intrinsic length.
    scalar minStartDist = scaleFactor*flParams[flParams.size()-1];

    minStartPt = -1;
    minEndPt = -1;
    scalar minEndDist = minStartDist;

    point flStart = featurel.point(0);
    point flEnd = featurel.point(featurel.size()-1);
    for (int i = 0; i < featurel.size(); i++)
    {
        point pt = featurel.point(i);
        if (bb.contains(pt))
        {
            flStart = pt;
            break;
        }
    }

    if (!lineIsLoop)
    {
        for (int i = featurel.size()-1; i >= 0; i--)
        {
            point pt = featurel.point(i);
            if (bb.contains(pt))
            {
                flEnd = pt;
                break;
            }
        }
    }

    forAll(pp_.localPoints(),pI)
    {
        scalar distStart = mag(pp_.localPoints()[pI] - flStart);
        if (distStart < minStartDist)
        {
            minStartPt = pI;
            minStartDist = distStart;
        }

        if (!lineIsLoop)
        {
            scalar distEnd = mag
                             (
                                 pp_.localPoints()[pI]
                                 - flEnd
                             );
            if (distEnd < minEndDist)
            {
                minEndPt = pI;
                minEndDist = distEnd;
            }
        }
    }

    if (lineIsLoop)
    {
        minEndPt = minStartPt;
    }

    // try and find a better starting point
    if (minStartPt != -1 && minEndPt != -1)
    {
        bool findNewMin =  true;
        scalar distMag = mag(pp_.localPoints()[minStartPt]
                             - flStart);

        if (distMag < lengthAve[minStartPt])
        {
            findNewMin = false;
        }

        bool findNewMax =  true;
        distMag = mag(pp_.localPoints()[minEndPt]
                      - flEnd);
        if (distMag < lengthAve[minEndPt])
        {
            findNewMax = false;
        }

        label curPt = minStartPt;
        scalar globalMin = GREAT;
        bool onLine = false;

        if (!findNewMin)
        {
            onLine = true;
        }

        while (findNewMin)
        {
            labelList nodeMfd = neigbouringPoints(curPt, true, true);

            scalar minDist =  GREAT;
            label minLoc = 0;
            forAll(nodeMfd,i)
            {
                pointIndexHit hit = tree.findNearest
                    (pp_.localPoints()[nodeMfd[i]],GREAT);
                if (hit.hit())
                {
                    scalar distSqr =
                        (pp_.localPoints()[nodeMfd[i]] - hit.hitPoint())
                        & (pp_.localPoints()[nodeMfd[i]] - hit.hitPoint()) ;
                    if (distSqr < minDist)
                    {
                        minLoc = i;
                        minDist = distSqr;
                    }
                }
            }
            scalar dist = sqrt(minDist);

            if (dist < lengthAve[nodeMfd[minLoc]])
            {
                minStartDist = mag
                (
                    pp_.localPoints()[nodeMfd[minLoc]] -
                    flStart
                 );
                minStartPt = nodeMfd[minLoc];
                onLine = true;
                break;
            }
            else if (dist >= globalMin)
            {
                minStartDist =
                    mag(pp_.localPoints()[curPt] - flStart);
                minStartPt = curPt;
                break;
            }
            else
            {
                globalMin = dist;
                curPt = nodeMfd[minLoc];
            }
        }

        if (onLine)
        {
            curPt = minStartPt;
            bool found = true;
            pointIndexHit currentHit = tree.findNearest
                (pp_.localPoints()[curPt],GREAT);

            while (found)
            {
                found = false;

                if (currentHit.hit())
                {
                    scalar currentFtrDist = mag(currentHit.hitPoint()
                                          - featurel.point(currentHit.index()))
                        + flParams[currentHit.index()];

                    labelList nodeMfd = neigbouringPoints(curPt, true, false);

                    List<nodeMfdProjResult> projResults(nodeMfd.size());
                    forAll(nodeMfd,i)
                    {
                        pointIndexHit hit = tree.findNearest
                            (pp_.localPoints()[nodeMfd[i]],GREAT);

                        scalar distSqr = (pp_.localPoints()[nodeMfd[i]]
                                       - hit.hitPoint())
                            & (pp_.localPoints()[nodeMfd[i]] - hit.hitPoint()) ;

                        projResults[i] = nodeMfdProjResult(i,distSqr,hit);
                    }
                    //sort the results so that the closest comes first
                    std::sort(projResults.begin(),projResults.end(),nMPRCmp());

                    forAll(projResults,pI)
                    {
                        const nodeMfdProjResult &nearest = projResults[pI];
                        label nearestInd = nodeMfd[nearest.index];
                        vector ftrHitPoint = nearest.hit.hitPoint();
                        label nearestFtrEdge = nearest.hit.index();
                        scalar minDistSqr = nearest.distSqr;

                        scalar ftrDist = mag(ftrHitPoint
                                             - featurel.point(nearestFtrEdge))
                                       + flParams[nearestFtrEdge];

                        scalar ftrLength = mag(ftrDist - currentFtrDist);
                        scalar edgeLength = mag(pp_.localPoints()[nearestInd]
                                                - pp_.localPoints()[curPt]);
                        scalar ratio = ftrLength / (edgeLength + SMALL);

                        if
                        (
                            ftrDist < currentFtrDist
                            && ftrDist < 0.5*flParams[flParams.size()-1]
                            && ratio > 0.1
                            && minDistSqr
                            < 0.4 * lengthAve[curPt] * lengthAve[curPt]
                        )
                        {
                            found = true;
                            curPt = nearestInd;
                            break;
                        }
                    }
                }
                currentHit = tree.findNearest
                    (pp_.localPoints()[curPt],GREAT);
            }


            minStartDist =
                mag(pp_.localPoints()[curPt] - flStart);
            minStartPt = curPt;
        }

        if (!lineIsLoop)
        {
            curPt = minEndPt;
            globalMin = GREAT;

            onLine = false;

            if (!findNewMax)
            {
                onLine = true;
            }

            while (true && findNewMax)
            {
                labelList nodeMfd = neigbouringPoints(curPt, true, true);

                scalar minDist =  GREAT;
                label minLoc = 0;
                forAll(nodeMfd,i)
                {
                    pointIndexHit hit = tree.findNearest
                        (pp_.localPoints()[nodeMfd[i]],GREAT);
                    if (hit.hit())
                    {
                        scalar distSqr =
                            (pp_.localPoints()[nodeMfd[i]] - hit.hitPoint())
                            & (pp_.localPoints()[nodeMfd[i]] - hit.hitPoint()) ;
                        if (distSqr < minDist)
                        {
                            minLoc = i;
                            minDist = distSqr;
                        }
                        }
                }
                scalar dist = sqrt(minDist);

                if (dist < lengthAve[nodeMfd[minLoc]])
                {
                    minEndDist =
                        mag(pp_.localPoints()[nodeMfd[minLoc]]
                            - flEnd);
                    minEndPt = nodeMfd[minLoc];
                    onLine = true;
                    break;
                }
                else if (dist >= globalMin)
                {
                    minEndDist = mag(pp_.localPoints()[curPt]
                                     - flEnd);
                    minEndPt = curPt;
                    break;
                }
                else
                {
                    globalMin = dist;
                    curPt = nodeMfd[minLoc];
                }
            }

            if (onLine)
            {
                curPt = minEndPt;
                bool found = true;

                pointIndexHit currentHit = tree.findNearest
                    (pp_.localPoints()[curPt],GREAT);

                while (found)
                {
                    found = false;

                    if (currentHit.hit())
                    {
                        scalar currentFtrDist = mag(currentHit.hitPoint()
                            - featurel.point(currentHit.index()))
                            + flParams[currentHit.index()];

                        labelList nodeMfd = neigbouringPoints
                                            (
                                                curPt,
                                                true,
                                                false
                                            );

                        List<nodeMfdProjResult> projResults(nodeMfd.size());
                        forAll(nodeMfd,i)
                        {
                            pointIndexHit hit = tree.findNearest
                                (pp_.localPoints()[nodeMfd[i]],GREAT);

                            scalar distSqr = (pp_.localPoints()[nodeMfd[i]]
                                              - hit.hitPoint())
                                & (pp_.localPoints()[nodeMfd[i]]
                                   - hit.hitPoint()) ;

                            projResults[i] = nodeMfdProjResult(i,distSqr,hit);
                        }
                        //sort the results so that the closest comes first
                        std::sort
                        (
                            projResults.begin(),
                            projResults.end(),
                            nMPRCmp()
                        );

                        forAll(projResults,pI)
                        {
                            const nodeMfdProjResult &nearest = projResults[pI];
                            label nearestInd = nodeMfd[nearest.index];
                            vector ftrHitPoint = nearest.hit.hitPoint();
                            label nearestFtrEdge = nearest.hit.index();
                            scalar minDistSqr = nearest.distSqr;

                            scalar ftrDist =
                                mag(ftrHitPoint
                                    - featurel.point(nearestFtrEdge))
                                + flParams[nearestFtrEdge];

                            scalar ftrLength = mag(ftrDist - currentFtrDist);
                            scalar edgeLength =
                                mag(pp_.localPoints()[nearestInd] -
                                                    pp_.localPoints()[curPt]);
                            scalar ratio = ftrLength / (edgeLength + SMALL);

                            if
                            (
                                ftrDist > currentFtrDist
                                && ftrDist > 0.5*flParams[flParams.size()-1]
                                && ratio > 0.1
                                && minDistSqr
                                < 0.4 * lengthAve[curPt] * lengthAve[curPt]
                             )
                            {
                                found = true;
                                curPt = nearestInd;
                                break;
                            }
                        }
                    }
                    currentHit = tree.findNearest
                        (pp_.localPoints()[curPt],GREAT);
                }

                minEndDist =
                    mag(pp_.localPoints()[curPt]
                        - flEnd);
                minEndPt = curPt;
            }

        }

    }

    if (lineIsLoop)
    {
        minEndPt = minStartPt;
    }

}


/**
 * Loops over the feature lines constructing shadow chains for each
 * -# Finds closest patch points to fl end points and checks they're actually
 * on the nearest surface
 * -# Loops along calling featureLineSnapper::findNext between start and end
 */
void Foam::featureLineSnapper::buildShadowChains()
{
    simpleVTKWriter failedLine;
    simpleVTKWriter failedProj;
    label trialCount = 0;
    label disableCount = 0;
    label countAllowed = 0;

    const fvMesh& mesh = meshRefiner_.mesh();
    treeBoundBox bb(mesh.points());

    //Checking whether a processor point
    boolList procPoint(pp_.localPoints().size(), false);
    if (Pstream::parRun())
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        const label nInternalFaces = mesh.nInternalFaces();

        forAll(pp_.meshPoints(), i)
        {
            label meshPointI = pp_.meshPoints()[i];
            const labelList pFaces =
                meshRefiner_.mesh().pointFaces()[meshPointI];

            forAll(pFaces, j)
            {
                if (pFaces[j] >= nInternalFaces)
                {
                    label patchI = patches.whichPatch(pFaces[j]);
                    if (patches[patchI].coupled())
                    {
                        procPoint[i] = true;
                        break;
                    }
                }
            }
        }
    }

    //Note: using *local* point indices in the patch (pp_) (but not in the
    //feature lines (surface mesh))

    List<scalar> lengthAve(pp_.localPoints().size(),-GREAT);

    //loop over the displacement field
    forAll(pp_.localPoints(),pI)
    {
        //what is the maximum edge length at this node?  (note: we don't
        //bother with parallel syncing, as this is just an indicative
        //lengthscale and doesn't need to be accurate)
        forAll(pp_.pointEdges()[pI],j)
        {
            lengthAve[pI] = max
            (
                lengthAve[pI],
                mag
                (
                    pp_.localPoints()[pI]
                    - pp_.localPoints()
                      [
                          pp_.edges()[pp_.pointEdges()[pI][j]].otherVertex(pI)
                      ]
                 )
             );
        }
        lengthAve[pI] *= 0.707;
    }


    DynamicList<DynamicList<label>> rawShadows;
    DynamicList<DynamicList<vector>> rawDisp;
    DynamicList<bool> rawSurfBound;
    vectorField allDisp(pp_.localPoints().size(), vector::zero);

    for (label flI = 0; flI < flPrep_.getNumFL(); ++flI)
    {
        const featureLine &featurel = flPrep_.getFeatureLine(flI);
        const List<scalar> &flParams = flPrep_.getFLParams(flI);

        const bool boundFtr = featurel.minDihedral() < -1;

        if (allowQuadSplits_) countAllowed++;
        const indexedOctree<treeDataEdge> tree = flPrep_.getFLSearchTree(flI);

        //Find primitive patch start and end points for feature line
        //----------------------------------------------------------
        label minStartPt = -1;
        label minEndPt = -1;

        findEndPts
        (
            bb,
            featurel,
            flParams,
            tree,
            lengthAve,
            minStartPt,
            minEndPt
        );

        if (verbose_)
        {
            if (minStartPt != -1 && minEndPt != -1)
            {
                Pout<<"Feature line: "<< flI
                    <<" Start point: "<< pp_.localPoints()[minStartPt]
                    <<" End point: "<< pp_.localPoints()[minEndPt]
                    <<endl;
            }
        }

        //Now go ahead and build the chain
        //--------------------------------

        rawShadows.append(DynamicList<label>());
        rawDisp.append(DynamicList<vector>());
        rawSurfBound.append(boundFtr);

        // If we didn't find start and end points close enough, we just
        //leave the shadow chain empty
        if (minStartPt != -1 && minEndPt != -1)
        {
            DynamicList<label> &chain = rawShadows[rawShadows.size()-1];
            DynamicList<vector> &chainDisp = rawDisp[rawDisp.size()-1];

            label curPt = minStartPt;
            if (verbose_)
            {
                Pout<< "Adding start point " << curPt << " to shadow of "
                     << "feature line " << flI << "(which now contains "
                     << chain.size()+1 << " entries)\n\n";
            }

            {
                point cPoint = pp_.localPoints()[curPt];

                if (mag(allDisp[curPt]) > SMALL)
                {
                    cPoint += allDisp[curPt];
                }
                point flp;
                pointIndexHit hitf = tree.findNearest
                    (cPoint,GREAT);
                if (hitf.hit())
                {
                    flp = hitf.hitPoint();
                    if (mag(flp - pp_.localPoints()[curPt]) < lengthAve[curPt])
                    {
                        allDisp[curPt] = flp - pp_.localPoints()[curPt];
                    }
                }
                else
                {
                    flp = featurel.point(0);
                }
                chain.append(curPt);
                chainDisp.append(flp - pp_.localPoints()[curPt]);
            }

            scalar curFtrDist = 0.0;
            scalar minFtrDist = GREAT;
            scalar maxFtrDist = -GREAT;
            vector curTargetPt = vector::zero;

            label foundNext = 0;
            while (foundNext == 0)
            {
                foundNext = findNext
                (
                    curPt,
                    curFtrDist,
                    curTargetPt,
                    tree,
                    featurel,
                    flParams,
                    procPoint[minEndPt],
                    minEndPt,
                    lengthAve,
                    allDisp,
                    false
                );

                minFtrDist = min(minFtrDist, curFtrDist);
                maxFtrDist = max(maxFtrDist, curFtrDist);

                if (foundNext == 0)
                {
                    if (verbose_)
                    {
                        Pout<< "Adding point " << curPt
                             << " to shadow of feature line " << flI
                             << "(which now contains " << chain.size()+1
                             << " entries)\n\n";
                    }
                    chain.append(curPt);
                    chainDisp.append(curTargetPt - pp_.localPoints()[curPt]);
                }
            }
            if (foundNext == 1)
            {
                if (verbose_)
                {
                    Pout<< "Adding end point " << minEndPt
                         << " to shadow of feature line " << flI
                         << "(which now contains " << chain.size()+1
                         << " entries)\n\n";
                }

                {
                    point cPoint = pp_.localPoints()[minEndPt];
                    if (mag(allDisp[minEndPt]) > SMALL)
                    {
                        cPoint += allDisp[minEndPt];
                    }
                    point flp;
                    pointIndexHit hitf = tree.findNearest
                        (cPoint,GREAT);
                    if (hitf.hit())
                    {
                        flp = hitf.hitPoint();
                        if
                        (
                            mag(flp - pp_.localPoints()[minEndPt])
                            < lengthAve[minEndPt]
                        )
                        {
                            allDisp[minEndPt] = flp - pp_.localPoints()[minEndPt];
                        }
                    }
                    else
                    {
                        flp = featurel.point(featurel.size()-1);
                    }
                    chain.append(minEndPt);
                    chainDisp.append(flp - pp_.localPoints()[minEndPt]);
                }
            }
            else if (foundNext == 2)
            {
                if (verbose_)
                {
                    Pout<< "Trying to build more of chain " << flI
                         << " by starting from the end point" << endl;
                }

                //we didn't make it to the end, so try building the chain
                //backwards from the end point
                scalar fwdFtrDist = curFtrDist;

                foundNext = 0;
                DynamicList<label> reverseChain;
                DynamicList<vector> reverseChainDisp;
                curPt = minEndPt;
                curFtrDist = mag(flParams[flParams.size()-1]);
                if (verbose_)
                {
                    Pout<< "Adding 'start' point " << curPt
                         << " to reverse shadow of feature line " << flI
                         << "(which now contains " << reverseChain.size()+1
                         << " entries)\n\n";
                }

                {
                    point cPoint = pp_.localPoints()[curPt];
                    if (mag(allDisp[curPt]) > SMALL)
                    {
                        cPoint += allDisp[curPt];
                    }
                    point flp;
                    pointIndexHit hitf = tree.findNearest
                        (pp_.localPoints()[curPt],GREAT);
                    if (hitf.hit())
                    {
                        flp = hitf.hitPoint();
                        if (mag(flp - pp_.localPoints()[curPt]) < lengthAve[curPt])
                        {
                            allDisp[curPt] = flp - pp_.localPoints()[curPt];
                        }
                    }
                    else
                    {
                        flp = featurel.point(featurel.size()-1);
                    }
                    reverseChain.append(curPt);
                    reverseChainDisp.append(flp - pp_.localPoints()[curPt]);
                }

                while (foundNext == 0)
                {
                    foundNext = findNext
                    (
                        curPt,
                        curFtrDist,
                        curTargetPt,
                        tree,
                        featurel,
                        flParams,
                        procPoint[chain[chain.size()-1]],
                        chain[chain.size()-1],
                        lengthAve,
                        allDisp,
                        true
                    );

                    //if we've gone back beyond where the forward chain got to,
                    //then stop(would be more elegant to build this into
                    //findNext but I don't have time)
                    if
                    (
                        curFtrDist < fwdFtrDist
                        && curFtrDist > minFtrDist
                        && (maxFtrDist - minFtrDist)
                        > 0.02* mag(flParams[flParams.size()-1])
                    )
                    {
                        foundNext = 1;
                    }

                    if (foundNext == 0)
                    {
                        if (verbose_)
                        {
                            Pout<< "Adding point " << curPt
                                 << " to shadow of feature line " << flI
                                 << "(which now contains " << chain.size()+1
                                 << " entries)\n\n";
                        }
                        reverseChain.append(curPt);
                        reverseChainDisp.append
                            (curTargetPt - pp_.localPoints()[curPt]);
                    }
                }

                //if the forward chain only contains the start point,
                //throw it out
                if (chain.size() == 1)
                {
                    chain.clear();
                }
                //if the reverse chain has more than just the end point in it,
                //and it doesn't overlap the forward part of the chain, add it
                //to the back of the main chain
                if (reverseChain.size() > 1)
                {
                    if (verbose_)
                    {
                        Pout<< " Managed to add " << reverseChain.size()
                             << " nodes to the chain by growing back from "
                             << "the end point" << endl << endl;
                    }
                    label totalSize = chain.size() + reverseChain.size();
                    label origChainSize = chain.size();
                    chain.setSize(totalSize);
                    chainDisp.setSize(totalSize);
                    for (label j = 0; j < reverseChain.size(); ++j)
                    {
                        chain[origChainSize + reverseChain.size()-j-1] =
                            reverseChain[j];
                        chainDisp[origChainSize + reverseChain.size()-j-1] =
                            reverseChainDisp[j];
                    }
                }
            }

            chain.shrink();
            chainDisp.shrink();

            if (!chainIsValid(chain, chainDisp, lengthAve))
            {
                chain.clear();
                chainDisp.clear();
            }
        }
        else
        {
            nNoValidEndPts_++;
        }
    }

    if (verbose_)
    {
        Pout<< countAllowed << " of " << flPrep_.getNumFL()
             << " feature lines were flat enough for quad splitting" << endl;
    }

    if (dump_)
    {
        failedLine.write("failedLines.vtk");
        failedProj.write("failedProjs.vtk");
    }

    if (verbose_) {
        Pout<< "Disabled " << disableCount << " of " << trialCount
             << " otherwise valid feature lines" << endl;
    }

    forAll(rawShadows, chainI)
    {
       const labelList& nchain = rawShadows[chainI];
       forAll(nchain, i)
       {
          if (mag(allDisp[nchain[i]]) > SMALL)
          {
             rawDisp[chainI][i] = allDisp[nchain[i]];
          }
       }
    }

    //compact the shadow chains
    shadowChains_.setSize(rawShadows.size());
    shadowChainDisp_.setSize(rawShadows.size());
    shadowChainSurfBound_.setSize(rawShadows.size());

    forAll(shadowChains_,scI)
    {
        shadowChains_[scI].transfer(rawShadows[scI]);
        shadowChainDisp_[scI].transfer(rawDisp[scI]);
        shadowChainSurfBound_[scI] = rawSurfBound[scI];
    }

}


/**
 * The guts of the shadow chain construction.  Given the current shadow chain
 * point and the feature line arc length to its hit point, finds the next one
 * by testing all the points in the current point's extended node manifold to
 * find the closest valid point to the feature line.  Jumps out early if the
 * previously-found end point is in the node manifold.  Can fail to find any
 * valid point under some circumstances.
 */
Foam::label Foam::featureLineSnapper::findNext
(
    label &curPt,
    scalar &curFtrDist,
    vector &curTargetPt,
    const indexedOctree<treeDataEdge> &ftrLineTree,
    const featureLine &featurel,
    const List<scalar> &flParams,
    const bool procPt,
    const label endPt,
    const List<scalar> &lengthAve,
    vectorField& allDisp,
    bool reverse
) const
{
    scalar ftrDirSign = reverse ? -1.0 : 1.0; //flip feature distance checks
                          //if searching backwards along
                          //the fl

    labelList nodeMfd = neigbouringPoints(curPt, allowQuadSplits_, false);

    if (verbose_)
    {
        Pout<< " Found " << nodeMfd.size() << " surrounding points\n";
    }

    //If the node manifold contains the end point then we're at the end if the
    //distance from here to the end of the line is roughly correct.
    forAll(nodeMfd,i)
    {
        if (nodeMfd[i] == endPt)
        {
            scalar ftrDist = !reverse ? mag(flParams[flParams.size()-1]) : 0.0;
            scalar length = flParams[flParams.size()-1] + SMALL;
            scalar ratio = ftrDirSign*(ftrDist - curFtrDist) / length;

            if
            (
                (ratio < 0.5 && procPt) ||
                (ftrDirSign*(ftrDist - curFtrDist)
                 < 2.0*mag(pp_.localPoints()[endPt]
                           - pp_.localPoints()[curPt]))
            )
            {
                curPt = endPt;
                curFtrDist = ftrDist;
                //don't need to set curTargetPt because it won't be used
                //(end point's location already known)

                if (verbose_)
                {
                    Pout<< " End point is amongst surrounding points: "
                         << "algorithm is complete\n";
                }

                return 1;
            }
            //if the end point is there but it's a red herring, we need to
            //stop it being a temptation...
            else
            {
                if (verbose_)
                {
                    Pout<< " Found end point amongst surrounding points but "
                         << "it is too far away on the line\n";
                }
                labelList temp = nodeMfd;
                nodeMfd.setSize(temp.size()-1);
                forAll(nodeMfd,ii)
                {
                    nodeMfd[ii] = ii < i ? temp[ii] : temp[ii+1];
                }
                break;
            }
        }
    }

    if (verbose_)
    {
        Pout<< " Finding nearest point in manifold\n";
    }

    List<nodeMfdProjResult> projResults(nodeMfd.size());

    //Loop over all the nodes that are candidates for being the next point in
    //the shadow chain finding for each one the closest point on the feature
    //line
    forAll(nodeMfd,i)
    {
        point cPoint = pp_.localPoints()[nodeMfd[i]];
        if (mag(allDisp[nodeMfd[i]]) > SMALL)
        {
            cPoint += allDisp[nodeMfd[i]];
        }
        pointIndexHit hit = ftrLineTree.findNearest
            (cPoint,GREAT);

        scalar distSqr = (cPoint - hit.hitPoint())
            & (cPoint - hit.hitPoint()) ;

        if (verbose_)
        {
            Pout<< "  surrounding node " << i << "(" <<  nodeMfd[i]
                 << ") is a " << (hit.hit() ? "hit" : "miss")
                 << " at location: "<< pp_.localPoints()[nodeMfd[i]]
                 << ". Distance is " << sqrt(distSqr) << "\n";
        }
        projResults[i] = nodeMfdProjResult(i,distSqr,hit);
    }

    //sort the results so that the closest comes first
    std::sort(projResults.begin(),projResults.end(),nMPRCmp());

    //now choose the point closest to the line that is also a positive
    //parametric distance from the current point, and not at so great a
    //parametric distance that it is likely we've short-circuited the feature
    //line

    for (label repeat = 0; repeat < 2; ++repeat)
    {
        forAll(projResults,pI)
        {
            const nodeMfdProjResult &nearest = projResults[pI];
            label nearestInd = nearest.index;
            vector ftrHitPoint = nearest.hit.hitPoint();
            label nearestFtrEdge = nearest.hit.index();
            scalar minDistSqr = nearest.distSqr;

            if (verbose_)
            {
                Pout<< " Nearest point is surrounding node " << nearest.index
                     << "(" <<  nodeMfd[nearestInd] << ") at a distance "
                     << sqrt(minDistSqr) << "\n";
                Pout<< " Nearest feature edge is edge " << nearestFtrEdge
                     << "\n";
                Pout<< " Feature hit point is "
                     << mag(ftrHitPoint - featurel.point(nearestFtrEdge))
                     << " from previous feature line node\n";
            }

            //now need to work out the distance along the line from curFtrPt.
            //The tree tells us which edge we hit.  By construction edge n in
            //the feature line tree spans nodes n and n+1, so we add the
            //distance to node n to the normal parameter at node n to get our
            //new feature line distance, and then simply compare it with the
            //current one
            scalar ftrDist = mag(ftrHitPoint - featurel.point(nearestFtrEdge))
                + flParams[nearestFtrEdge];

            point nPoint = pp_.localPoints()[nodeMfd[nearestInd]];
            point cPoint = pp_.localPoints()[curPt];
            scalar geometricDist = mag(nPoint - cPoint);
            scalar parametricDist = ftrDirSign*(ftrDist - curFtrDist);

            if
            (
                parametricDist > 0.0 &&
                (
                    repeat > 0
                 || parametricDist < 2.0*geometricDist
                )
            )
            {
                if
                (
                    (!reverse && ftrDist + SMALL < flParams[flParams.size()-1])
                    || (reverse && ftrDist - SMALL > flParams[0])
                )
                {
                    //we have our new point!
                    if (verbose_)
                    {
                        Pout<< "SUCCESS: Found valid point " << curPt
                             << " at geometric distance "
                             << geometricDist << " and parametric distance "
                             << parametricDist << " from previous point"
                             << endl;
                        Pout<< endl;
                    }

                    if (isEdge(curPt,nodeMfd[nearestInd]))
                    {
                        curPt = nodeMfd[nearestInd];
                        curFtrDist = ftrDist;
                        curTargetPt = ftrHitPoint;
                        if
                        (
                            mag(ftrHitPoint - pp_.localPoints()[curPt])
                            < lengthAve[curPt]
                        )
                        {
                            allDisp[curPt] =  ftrHitPoint
                                - pp_.localPoints()[curPt];
                        }
                        return 0;
                    }
                    else
                    {
                        label origLocalID = -1;
                        face origNew;
                        face newNew;
                        if
                        (
                            getFaceSplit
                            (
                                curPt,
                                nodeMfd[nearestInd],
                                origLocalID,

                                origNew,
                                newNew
                             )
                         )
                        {

                            face origNewLocal(origNew);
                            face newNewLocal(newNew);

                            //map face node numbers to mesh numbering:
                            forAll(origNew, i)
                            {
                                origNew[i] = pp_.meshPoints()[origNew[i]];
                            }

                            forAll(newNew, i)
                            {
                                newNew[i] = pp_.meshPoints()[newNew[i]];
                            }

                            label origID = pp_.addressing()[origLocalID];

                            vector fNorm =
                               origNew.areaNormal(meshRefiner_.mesh().points());
                            vector nNorm =
                                newNew.areaNormal(meshRefiner_.mesh().points());

                            scalar faceArea = mag(fNorm);
                            scalar newFaceArea = mag(nNorm);

                            scalar areaRatio = faceArea /
                                ( meshRefiner_.mesh().magFaceAreas()[origID]
                                  + SMALL);

                            if (areaRatio >= 0.25 && areaRatio < 0.75)
                            {
                                fNorm /= (faceArea + SMALL);
                                nNorm /= (newFaceArea + SMALL);

                                if ((fNorm & nNorm) > 0.965)
                                {
                                    curPt = nodeMfd[nearestInd];
                                    curFtrDist = ftrDist;
                                    curTargetPt = ftrHitPoint;
                                    if
                                    (
                                        mag(ftrHitPoint - pp_.localPoints()[curPt])
                                       < lengthAve[curPt]
                                    )
                                    {
                                        allDisp[curPt] =  ftrHitPoint
                                        - pp_.localPoints()[curPt];
                                    }
                                    return 0;
                                }
                            }
                        }
                    }
                }
                else
                {
                    //we've got to the end but this point would get
                    //snapped to the same position as the end pt
                    //so let's leave it at that.
                    return 1;
                }
            }
            else
            {
                //This is not the point we are looking for...
                //loop around to try the next one
                if (verbose_)
                {
                    Pout<< " Rejecting point " << nearestInd << "("
                         << nodeMfd[nearestInd] << ") because it is "
                         << geometricDist << "\n";
                    Pout<< "  from the previous point but " << parametricDist
                         << " away along the line\n";
                }
            }
        }
        if (verbose_ && repeat == 0)
        {
            Pout<< " Failed to find a point on first pass.  Going round "
                 << "again but lifting upper limit on parametric distance"
                 << endl;
        }
    }
    if (verbose_)
    {
        Pout<< "ERROR: featureLineSnapper::findNext: failed to find a "
             << "valid next node!\n";
    }
    nFailedConstr_++;
    return 2;
}



/**
 * Checks a couple of chain properties to filter out obviously invalid chains
 */
bool Foam::featureLineSnapper::chainIsValid
(
    DynamicList<label> &chain,
    DynamicList<vector> &chainDisp,
    const List<scalar> &lengthAve
) const
{
    //Throw away chains that form impossible loops
    if
    (
        chain.size() < 4
        && chain.size() > 0 && chain[0] == chain[chain.size()-1]
    )
    {
        return false;
    }

    DynamicList<label> validChain(chain.size());
    DynamicList<vector> validDisp(chain.size());

    forAll(chain, i)
    {
        const label pointI = chain[i];

        if (mag(chainDisp[i]) < lengthAve[pointI])
        {
            validChain.append(chain[i]);
            validDisp.append(chainDisp[i]);
        }
    }

    chain = validChain;
    chainDisp = validDisp;

    if
    (
        chain.size() < 4
        && chain.size() > 0 && chain[0] == chain[chain.size()-1]
    )
    {
        return false;
    }
    else
    {
        return true;
    }
}


/**
 * Looks for where shadow chains completely or partially duplicate each other
 * and eliminates the duplication.  Note: this function reconstructs the
 * shadow chain list and as a result removes the correspondence between
 * feature lines and shadow chains via their indices.
 */
Foam::label Foam::featureLineSnapper::filterChains()
{
    //build up chain appearance info for all the nodes in the patch in a
    //ragged array
    labelList chainAppCount(pp_.meshPoints().size(),0);
    label totalNumChainNodes = 0;
    {
        forAll(shadowChains_,scI)
        {
            forAll(shadowChains_[scI],j)
            {
                //don't count end points
                if (j != 0 && j != shadowChains_[scI].size()-1)
                {
                    chainAppCount[shadowChains_[scI][j]]++;
                    totalNumChainNodes++;
                }
            }
        }
    }

    labelList chainAppStartInd(chainAppCount.size()+1,0);
    {
        forAll(chainAppCount,j)
        {
            chainAppStartInd[j+1] = chainAppStartInd[j] + chainAppCount[j];
        }
    }

    labelList chainApps(totalNumChainNodes);
    {
        labelList currentStartOffset(pp_.meshPoints().size(),0);
        forAll(shadowChains_,scI)
        {
            forAll(shadowChains_[scI],j)
            {
                //don't count the end points, because they coincide
                //"naturally" with other feature lines
                if (j != 0 && j != shadowChains_[scI].size()-1)
                {
                    label start = chainAppStartInd[shadowChains_[scI][j]];
                    label offset = currentStartOffset[shadowChains_[scI][j]]++;
                    chainApps[start+offset] = scI;
                }
            }
        }
    }

    //loop over chains.

    List<boolList> deletionMask(shadowChains_.size());

    forAll(shadowChains_,scI)
    {
        const labelList &chain = shadowChains_[scI];
        deletionMask[scI].setSize(chain.size(),false);

        labelList countIndex(shadowChains_.size(),-1);
        DynamicList<label> counts;

        //walk along counting how many chains each one coincides with
        forAll(chain,j)
        {

            if (j != chain.size()-1 || chain[0] != chain[j])
            {

                label node = chain[j];
                //loop through all of this node's chain appearances,
                //adding them up
                for
                (
                    label k = chainAppStartInd[node];
                    k < chainAppStartInd[node+1];
                    ++k
                )
                {
                    //if we encounter chain chainApps[k] for the first time,
                    //we create an entry for it in the index
                    if (countIndex[chainApps[k]] == -1)
                    {
                        countIndex[chainApps[k]] = counts.size();
                        counts.append(1);
                    }
                    else
                    {
                        counts[countIndex[chainApps[k]]]++;
                    }
                }
            }
        }

        labelList reverseIndex(counts.size());
        label jj = 0;
        forAll(countIndex,j)
        {
            if (countIndex[j] != -1)
            {
                reverseIndex[jj++] = j;
            }
        }

        label chainLength = chain.size();

        if (chainLength)
        {

            if (chain[0] == chain[chain.size()-1])
            {
                chainLength--;
            }

            //loop over all the chains
            forAll(counts,j)
            {
                    //the chain this "counts" entry relates to
                label scJ = reverseIndex[j];

                //if this chain is the shorter of the two flag the
                //coincident nodes to be thrown away
                if
                (
                    chain.size() < shadowChains_[scJ].size()
                  || (chain.size() == shadowChains_[scJ].size() && scI < scJ)
                )
                {
                    if (verbose_)
                    {
                        Pout<< "Chain " << scJ << " coincides with chain "
                             << scI << " for " << (100.0*counts[j])/chainLength
                             << "% of its length (" << counts[j] << " nodes)"
                             << endl;
                    }

                    forAll(chain,k)
                    {
                        //loop over all the appearances of each of the chain's
                        //nodes in (other) chains and flag for deletion if
                        //chain scJ is one of them
                        if (k != 0 && k != chain.size()-1)
                        {
                            labelList::const_iterator begin =
                                chainApps.begin()+chainAppStartInd[chain[k]];
                            labelList::const_iterator end =
                                chainApps.begin()+chainAppStartInd[chain[k]+1];
                            if (std::find(begin,end,scJ) != end)
                            {
                                deletionMask[scI][k] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    //Now go along deleting the entries that have been flagged.  We do this by
    //creating an entirely new shadowChain list (applying identical changes to
    //the displacement arrays of course).  Crucially, we will only retain
    //undeleted parts of partially deleted chains if they are more than 10
    //nodes in length

    label origNumShadowChains = shadowChains_.size();

    DynamicList<labelList> newChains;
    DynamicList<List<vector>> newDisps;
    DynamicList<bool> newSurfBound;

    forAll(deletionMask,scI)
    {
        label startRetain = -1;
        label countRetentions = 0;
        forAll(deletionMask[scI],j)
        {
            //if this is a retained node...
            if (!deletionMask[scI][j])
            {
                //...and we aren't currently amongst a set of retained nodes,
                //set the mark...
                if (startRetain == -1)
                {
                    startRetain = j;
                }
                //...but if this is the last node in the chain and the
                //retained section is big enough...
                else if
                (
                    j == deletionMask[scI].size()-1
                    && (j-startRetain > 2
                    || (startRetain == 0 && j > 4))
                )
                {
                    //...retain from startRetain to j (including j)
                    newChains.append(labelList(j-startRetain+1));
                    newDisps.append(List<vector>(j-startRetain+1));
                    newSurfBound.append(shadowChainSurfBound_[scI]);
                    for (label k = startRetain; k <= j; ++k)
                    {
                        newChains[newChains.size()-1][k-startRetain] =
                            shadowChains_[scI][k];
                        newDisps[newDisps.size()-1][k-startRetain] =
                            shadowChainDisp_[scI][k];
                    }
                    countRetentions++;
                }
            }
            //if this isn't a retained node...
            else
            {
                //...but the mark was previously set and the retained section
                //is big enough...
                if (startRetain != -1 && j - startRetain > 3)
                {
                    //...retain from startRetain to j (excluding j)
                    newChains.append(labelList(j-startRetain));
                    newDisps.append(List<vector>(j-startRetain));
                    newSurfBound.append(shadowChainSurfBound_[scI]);
                    for (label k = startRetain; k < j; ++k)
                    {
                        newChains[newChains.size()-1][k-startRetain] =
                            shadowChains_[scI][k];
                        newDisps[newDisps.size()-1][k-startRetain] =
                            shadowChainDisp_[scI][k];
                    }
                    countRetentions++;
                }
                startRetain = -1;
            }
        }

        if (verbose_ && countRetentions > 1)
        {
            Pout<< "Retained " << countRetentions
                 << " section(s) of shadow chain " << scI << endl;
        }
    }

    shadowChains_.transfer(newChains);
    shadowChainDisp_.transfer(newDisps);
    shadowChainSurfBound_.transfer(newSurfBound);

    Info<< "  Before removing duplications: "
         << returnReduce(origNumShadowChains,sumOp<label>()) << " chains"
         << endl;
    Info<< "  After removing duplications : "
         << returnReduce(shadowChains_.size(),sumOp<label>()) << " chains"
         << endl;

    return origNumShadowChains - shadowChains_.size();
}


/**
 * Converts the displacements specified for the individual chains into a patch
 * displacement field.  Implements per-point displacement filter;
 * filterDisplacementField implements edge and face-based checks.
 */
void Foam::featureLineSnapper::chainDispToPatchDisp()
{

    forAll(shadowChains_,scI)
    {
        //copy the chain displacements into the point field
        forAll(shadowChains_[scI],j)
        {
            pointDisp_[shadowChains_[scI][j]] = shadowChainDisp_[scI][j];
        }
    }
}



/**
 * Checks the computed shadow chains for unwanted topologies.
 *  - checks for chain traversing n-1 sides of an n-gon - most likely two sides
 *    of a triangle
 *  - checks for chain traversing a quad that has actually been flagged to be
 *    split the other way because of a two-point singularity (not implemented
 *    yet)
 */
void Foam::featureLineSnapper::cleanUpChains()
{
    //loop through each chain finding where it traverses a face diagonal
    forAll(shadowChains_,scI)
    {
        if (shadowChains_[scI].size() > 0)
        {
            const labelList &chain = shadowChains_[scI];

            // Must avoid silly situation where the chain is a loop that's so
            //short we could "creep up behind" ourselves
            bool chainIsLoop = chain[0] == chain[chain.size()-1];
            if (!chainIsLoop || chain.size() > 5)
            {
                boolList removeNode(chain.size(),false);

                for (label pI=0; pI < chain.size()-2; ++pI)
                {
                    if (isEdge(chain[pI],chain[pI+2]))
                    {
                        //chain traverses two sides of a triangle - or
                        //will do once splits introduced.  Remove node pI+1
                        //removeNode[pI+1] = true;
                    }
                    else if
                    (
                        pI < chain.size()-3
                        && isEdge(chain[pI],chain[(pI+3)%chain.size()])
                    )
                    {
                        //chain traverses three sides of a quad.
                        //Remove nodes pI+1 and pI+2
                        removeNode[pI+1] = true;
                        removeNode[pI+2] = true;
                    }
                }

                //if it's a loop, we need to check across the "join" (remember
                //node n-1 is basically just a dummy that's equal to node 0 so it's
                //skipped in these checks).
                if (chainIsLoop)
                {
                    if (isEdge(chain[chain.size()-2],chain[1]))
                    {
                        removeNode[0] = true;
                    }
                    else if (isEdge(chain[chain.size()-3],chain[1]))
                    {
                        removeNode[chain.size()-2] = true;
                        removeNode[0] = true;
                    }
                    else if (isEdge(chain[chain.size()-2],chain[2]))
                    {
                        removeNode[0] = true;
                        removeNode[1] = true;
                    }
                }

                //Now loop through and reset the chain
                int nToRemove = std::count
                (
                    removeNode.begin(),
                    removeNode.end(),
                    true
                 );

                labelList newChain(chain.size()-nToRemove);

                label newPI = 0;
                forAll(chain,pI)
                {
                    if (!removeNode[pI])
                    {
                        newChain[newPI++] = chain[pI];
                    }
                }
                ///just check...
                if (newPI != newChain.size())
                {
                    Pout<< "WARNING: featureLineSnapper: cleanUpChains: "
                         <<"haven't constructed the new chain properly:\n";
                    Pout<< " one past final assigned index = " << newPI
                         << "; size of newChain = " << newChain.size() << "\n";
                }

                //join loop back up if start node was deleted
                if (chainIsLoop && removeNode[0])
                {
                    newChain[newChain.size()-1] = newChain[0];
                }

                    //overwrite the old chain
                shadowChains_[scI].setSize(newChain.size());
                shadowChains_[scI].transfer(newChain);
            }
        }
    }

    cleanDisplacementField();
}

/**
 *   Zeros displacement at non-chain points - may be required after
 *   topological clean-up of chains
 */
void Foam::featureLineSnapper::cleanDisplacementField()
{
    pointField cleanDisp(pointDisp_.size(),vector::zero);
    forAll(shadowChains_,scI)
    {
        forAll(shadowChains_[scI],j)
        {
            cleanDisp[shadowChains_[scI][j]] =
                pointDisp_[shadowChains_[scI][j]];
        }
    }
    pointDisp_.transfer(cleanDisp);
}


/**
 * Applies local limiter to displacements to avoid excessive distortion of the
 * mesh.  Checks non-pointwise problems, in contrast to chainDispToPatchDisp.
 *
 * Note: though the aim of this filter is to improve mesh quality, this
 * is not a quality-based limiter.  It would be possible to implement that if
 * required, but the way featureLineSnapper is currently used it doesn't make
 * sense because we should wait for the full snapping algorithm to complete
 * before checking quality properly.
 */
void Foam::featureLineSnapper::filterDisplacementField()
{
    // loop over shadow chains and set displacement to zero for surf boundaries
    // which are not mesh region boundaries
    {
        const fvMesh& mesh = meshRefiner_.mesh();

        labelList meshEdges(pp_.nEdges());
        forAll(meshEdges, patchEdgeI)
        {
           const edge& e = pp_.edges()[patchEdgeI];

           label v0 = pp_.meshPoints()[e[0]];
           label v1 = pp_.meshPoints()[e[1]];
           meshEdges[patchEdgeI] = meshTools::findEdge
           (
              mesh.edges(),
              mesh.pointEdges()[v0],
              v0,
              v1
           );
        }

        const labelListList& edgeFaces = pp_.edgeFaces();

        // unique edge centre to face centre vector.
        vectorField edgeCentreFaceCentre
            (mesh.nEdges(), vector(GREAT, GREAT, GREAT));

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = edgeFaces[edgeI];
            label meshEdgeI = meshEdges[edgeI];

            forAll(eFaces, i)
            {
                const point& edgeCentre =
                    mesh.edges()[meshEdgeI].centre(mesh.points());
                const point& faceCentre =
                    mesh.faceCentres()[pp_.addressing()[eFaces[i]]];

                vector dir = edgeCentre - faceCentre;

                uniqueEdgeDir()
                (
                    edgeCentreFaceCentre[meshEdgeI],
                    dir
                 );
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            edgeCentreFaceCentre,
            uniqueEdgeDir(),
            vector(GREAT, GREAT, GREAT)  // null value
         );

        boolList convexEdge(mesh.nEdges(), false);

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = edgeFaces[edgeI];
            label meshEdgeI = meshEdges[edgeI];
            convexEdge[meshEdgeI] = true;

            if (edgeCentreFaceCentre[meshEdgeI] != vector(GREAT, GREAT, GREAT))
            {
                forAll(eFaces, i)
                {
                    const point& edgeCentre =
                        mesh.edges()[meshEdgeI].centre(mesh.points());
                    const point& faceCentre =
                        mesh.faceCentres()[pp_.addressing()[eFaces[i]]];

                    vector dir = edgeCentre - faceCentre;
                    dir /= (mag(dir) + SMALL);
                    vector uniqueEdgeDir = edgeCentreFaceCentre[meshEdgeI];
                    uniqueEdgeDir /= (mag(uniqueEdgeDir) + SMALL);

                    if ((uniqueEdgeDir &  dir) < -0.1736)
                    {
                        convexEdge[meshEdgeI] = false;
                        break;
                    }
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            convexEdge,
            andEqOp<bool>(),
            false           // null value
         );


        List<scalar> lengthAve(pp_.localPoints().size(),-GREAT);

        //loop over the displacement field
        forAll(pp_.localPoints(),pI)
        {
                //what is the maximum edge length at this node?  (note: we don't
                //bother with parallel syncing, as this is just an indicative
                //lengthscale and doesn't need to be accurate)
            forAll(pp_.pointEdges()[pI],j)
            {
                lengthAve[pI] = max
                (
                    lengthAve[pI],
                    mag
                    (
                        pp_.localPoints()[pI]
                        - pp_.localPoints()
                        [
                            pp_.edges()[pp_.pointEdges()[pI][j]].otherVertex(pI)
                         ]
                     )
                 );
            }
        }

        boolList foundPt(pp_.meshPoints().size(), false);
        labelListList edgeRegions;
        List<DynamicList<label>> meshEdgeRegions(mesh.nEdges());

        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        //Look at edge rotation and limit
        forAll(shadowChains_, scI)
        {
           const labelList &chain = shadowChains_[scI];

           for (label pI=0; pI < chain.size()-1; ++pI)
           {
              const label v0 = chain[pI];
              const label v1 = chain[pI+1];

              label edgeI = meshTools::findEdge
              (
                  pp_.edges(),
                  pp_.pointEdges()[v0],
                  v0,
                  v1
              );

              if
              (
                 edgeI != -1
              )
              {
                  label meshEdgeI = meshEdges[edgeI];

                  const labelList& eFaces = edgeFaces[edgeI];
                  const point& edgeCentre =
                      mesh.edges()[meshEdgeI].centre(mesh.points());
                  const point& faceCentre =
                      mesh.faceCentres()[pp_.addressing()[eFaces[0]]];

                  vector eCfC = edgeCentre - faceCentre;
                  eCfC /= (mag(eCfC) + SMALL);

                  point pt0 = pp_.localPoints()[v0];
                  point pt1 = pp_.localPoints()[v1];

                  vector v0v1 = pt1 - pt0;
                  scalar magv0v1 = mag(v0v1) + SMALL;
                  v0v1 /= magv0v1;

                  point pt0m = pp_.localPoints()[v0] + pointDisp_[v0];
                  point pt1m = pp_.localPoints()[v1] + pointDisp_[v1];

                  vector v0v1m = pt1m - pt0m;
                  scalar magv0v1m = mag(v0v1m) + SMALL;
                  v0v1m /= magv0v1m;

                  if ((v0v1 & v0v1m) < 0.9848)
                  {
                      if (mag(eCfC & v0v1m) > mag(eCfC & v0v1))
                      {
                          foundPt[v0] = true;
                          foundPt[v1] = true;
                      }
                  }

              }
           }
        }

        forAll(pp_.edges(), edgeI)
        {
           label meshEdgeI = meshEdges[edgeI];
           const labelList& meshEdgeFaces = pp_.edgeFaces()[edgeI];
           forAll(meshEdgeFaces, i)
           {
              label faceI = pp_.addressing()[meshEdgeFaces[i]];

              if (!mesh.isInternalFace(faceI))
              {
                 label patchI = bMesh.whichPatch(faceI);
                 if (!isA<processorPolyPatch>(bMesh[patchI]))
                 {
                    if (findIndex(meshEdgeRegions[meshEdgeI], patchI) == -1)
                    {
                       meshEdgeRegions[meshEdgeI].append(patchI);
                    }
                 }
              }
           }
           meshEdgeRegions[meshEdgeI].shrink();
        }

        edgeRegions.setSize(meshEdgeRegions.size());
        forAll(edgeRegions, i)
        {
           edgeRegions[i].transfer(meshEdgeRegions[i]);
        }

        // Synchronise across coupled edges.
        syncTools::syncEdgeList
        (
           mesh,
           edgeRegions,
           uniqueEqOp(),
           labelList()    // null value
        );

        forAll(shadowChains_, scI)
        {
           const labelList &chain = shadowChains_[scI];

           for (label pI=0; pI < chain.size()-1; ++pI)
           {
              const label v0 = pp_.meshPoints()[chain[pI]];
              const label v1 = pp_.meshPoints()[chain[pI+1]];

              label meshEdgeI = meshTools::findEdge
              (
                  mesh.edges(),
                  mesh.pointEdges()[v0],
                  v0,
                  v1
              );

              if
              (
                 meshEdgeI != -1
                 && shadowChainSurfBound_[scI]
                 && edgeRegions[meshEdgeI].size() == 1
                 && !convexEdge[meshEdgeI]
              )
              {
                  foundPt[chain[pI]] = true;
                  foundPt[chain[pI+1]] = true;
              }
              else if
              (
                  meshEdgeI == -1
               && shadowChainSurfBound_[scI]
              )
              {
                  //find the face they have in common
                  labelList faces1 = pp_.pointFaces()[chain[pI]];
                  labelList faces2 = pp_.pointFaces()[chain[pI+1]];
                  std::list<label> intx;
                  std::set_intersection
                  (
                      faces1.begin(),
                      faces1.end(),
                      faces2.begin(),
                      faces2.end(),
                      std::back_inserter(intx)
                  );
                  label numCommon = intx.size();

                  if (numCommon == 1)
                  {
                      label origFaceID = intx.front();
                      const labelList& fEdges = pp_.faceEdges()[origFaceID];
                      bool foundMultiRegion = false;

                      forAll(fEdges, fp)
                      {
                          label edgeI = fEdges[fp];
                          const edge& e = pp_.edges()[edgeI];
                          const label mp0 = pp_.meshPoints()[e[0]];
                          const label mp1 = pp_.meshPoints()[e[1]];

                          label faceMeshEdgeI = meshTools::findEdge
                          (
                              mesh.edges(),
                              mesh.pointEdges()[mp0],
                              mp0,
                              mp1
                          );

                          if
                          (
                              faceMeshEdgeI != -1
                              &&
                              (
                                  edgeRegions[faceMeshEdgeI].size() > 1
                                  || convexEdge[faceMeshEdgeI]
                              )
                          )
                          {
                              foundMultiRegion = true;
                              break;
                          }
                      }

                      if (!foundMultiRegion)
                      {
                          foundPt[chain[pI]] = true;
                          foundPt[chain[pI+1]] = true;
                      }
                  }
                  else
                  {
                      foundPt[chain[pI]] = true;
                      foundPt[chain[pI+1]] = true;
                  }
              }
           }
        }

        forAll(shadowChains_, scI)
        {
           const labelList &chain = shadowChains_[scI];

           for (label pI=0; pI < chain.size()-1; ++pI)
           {
              const label v0 = pp_.meshPoints()[chain[pI]];
              const label v1 = pp_.meshPoints()[chain[pI+1]];

              label meshEdgeI = meshTools::findEdge
              (
                  mesh.edges(),
                  mesh.pointEdges()[v0],
                  v0,
                  v1
              );

              if
              (
                 meshEdgeI != -1
                 &&
                 (
                     edgeRegions[meshEdgeI].size() > 1
                     || convexEdge[meshEdgeI]
                 )
              )
              {
                      foundPt[chain[pI]] = false;
                      foundPt[chain[pI+1]] = false;
              }
              else if
              (
                  meshEdgeI == -1
              )
              {
                  //find the face they have in common
                  labelList faces1 = pp_.pointFaces()[chain[pI]];
                  labelList faces2 = pp_.pointFaces()[chain[pI+1]];
                  std::list<label> intx;
                  std::set_intersection
                  (
                      faces1.begin(),
                      faces1.end(),
                      faces2.begin(),
                      faces2.end(),
                      std::back_inserter(intx)
                  );
                  label numCommon = intx.size();

                  if (numCommon == 1)
                  {
                      label origFaceID = intx.front();
                      const labelList& fEdges = pp_.faceEdges()[origFaceID];
                      bool foundMultiRegion = false;

                      forAll(fEdges, fp)
                      {
                          label edgeI = fEdges[fp];
                          const edge& e = pp_.edges()[edgeI];
                          const label mp0 = pp_.meshPoints()[e[0]];
                          const label mp1 = pp_.meshPoints()[e[1]];

                          label faceMeshEdgeI = meshTools::findEdge
                          (
                              mesh.edges(),
                              mesh.pointEdges()[mp0],
                              mp0,
                              mp1
                          );

                          if
                          (
                              faceMeshEdgeI != -1
                              &&
                              (
                                  edgeRegions[faceMeshEdgeI].size() > 1
                                  || convexEdge[faceMeshEdgeI]
                              )
                          )
                          {
                              foundMultiRegion = true;
                              break;
                          }
                      }

                      if (foundMultiRegion)
                      {
                              foundPt[chain[pI]] = false;
                              foundPt[chain[pI+1]] = false;
                      }
                  }
              }
           }
        }
        const scalar validScale(allowQuadSplits_ ? 0.6 : 0.75);
        const scalar validConvexScale = 1.5;//1.5;

        //Unset if edges shrinking too much
        forAll(pp_.edges(), edgeI)
        {
            const edge e = pp_.edges()[edgeI];
            const label v0 = e[0];
            const label v1 = e[1];

            point pt0 = pp_.localPoints()[v0];
            point pt1 = pp_.localPoints()[v1];

            vector v0v1 = pt1 - pt0;
            scalar magv0v1 = mag(v0v1) + SMALL;

            point pt0m = pp_.localPoints()[v0] + pointDisp_[v0];
            point pt1m = pp_.localPoints()[v1] + pointDisp_[v1];

            vector v0v1m = pt1m - pt0m;
            scalar magv0v1m = mag(v0v1m) + SMALL;


            scalar eLengthRatio = magv0v1m / magv0v1;
            if (eLengthRatio < 0.25)
            {
                foundPt[v0] = true;
                foundPt[v1] = true;
            }
        }

        //go along the chain checking how the displacement compares to the local
        //scale, as embodied in the average manifold edge length.

        boolList convexPts(pp_.meshPoints().size(), false);

        forAll(pp_.edges(), edgeI)
        {
            const edge e = pp_.edges()[edgeI];

            const label v0 = pp_.meshPoints()[e[0]];
            const label v1 = pp_.meshPoints()[e[1]];

            label meshEdgeI = meshTools::findEdge
            (
                mesh.edges(),
                mesh.pointEdges()[v0],
                v0,
                v1
            );

            scalar vScale = validConvexScale;
            if (convexEdge[meshEdgeI])
            {
                convexPts[e[0]] = true;
                convexPts[e[1]] = true;
                bool validDisp = mag(pointDisp_[e[0]]) <
                    vScale*lengthAve[e[0]];

                if (!validDisp)
                {
                    foundPt[e[0]] = true;
                }
                validDisp = mag(pointDisp_[e[1]]) <
                    vScale*lengthAve[e[1]];

                if (!validDisp)
                {
                    foundPt[e[1]] = true;
                }
            }
        }

        forAll(pp_.meshPoints(), pointI)
        {
            scalar vScale = validScale;
            if (!convexPts[pointI])
            {
                bool validDisp = mag(pointDisp_[pointI]) <
                    vScale*lengthAve[pointI];

                if (!validDisp)
                {
                    foundPt[pointI] = true;
                }
            }
        }

        //Reset displacement if face area ratio change is too large
        forAll(pp_, facei)
        {
            const face& f = pp_.localFaces()[facei];

            scalar fArea0 = pp_.magFaceAreas()[facei];
            face fNew(identity(f.size()));
            pointField newFacePts(f.size());
            forAll(f,fp)
            {
                label pointi = f[fp];
                newFacePts[fp] =  pp_.localPoints()[pointi]
                    + pointDisp_[pointi];
            }
            scalar fArea1 = fNew.mag(newFacePts);
            scalar aRatio = fArea1 / (fArea0 + SMALL);

            if (aRatio < 0.1)
            {
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    foundPt[pointi] = true;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            pp_.meshPoints(),
            foundPt,
            orEqOp<bool>(),
            false              // null value
        );

        forAll(foundPt, i)
        {
           if (foundPt[i])
           {
               pointDisp_[i] = vector::zero;
           }
        }
    }

    //loop over points looking for any points that have moved through a surface
   pointField start(pp_.localPoints().size());
   pointField end(pp_.localPoints().size());
   labelList pointIndex(pp_.localPoints().size());
   label testI = 0;

   forAll(pp_.localPoints(),pI)
   {
      if (mag(pointDisp_[pI]) > SMALL)
      {
         const labelList& pFaces = pp_.pointFaces()[pI];
         // choose one cc in the domain
         const label owningCell =
               meshRefiner_.mesh().faceOwner()[pp_.addressing()[pFaces[0]]];
         point pointCC = meshRefiner_.mesh().cellCentres()[owningCell];

         start[testI] = pointCC;
         end[testI] = pp_.localPoints()[pI] + pointDisp_[pI];
         pointIndex[testI] = pI;
         testI++;
      }
   }

   start.setSize(testI);
   end.setSize(testI);
   pointIndex.setSize(testI);

   labelList surface1;
   List<pointIndexHit> hitPoint1;
   labelList pointRegion1;
   labelList surface2;
   List<pointIndexHit> hitPoint2;
   labelList pointRegion2;

   const labelList unzonedSurfaces = surfaceZonesInfo::getUnnamedSurfaces
   (
       meshRefiner_.surfaces().surfZones()
   );

   meshRefiner_.surfaces().findNearestIntersection
   (
      unzonedSurfaces,
      start,
      end,
      surface1,
      hitPoint1,
      pointRegion1,

      surface2,
      hitPoint2,
      pointRegion2
   );

   forAll(start, pI)
   {
      label i = pointIndex[pI];
      if (surface1[pI] != -1)
      {
          vector vec1 = end[pI] - start[pI];
          vector vec2 = hitPoint1[pI].hitPoint() - start[pI];

          scalar d1 = mag(vec1);
          scalar d2 = mag(vec2);

          if ((d2 / (d1 + SMALL)) < 0.98)
          {
              pointDisp_[i] = vector::zero;//hitPoint1[pI].hitPoint()-pp_.localPoints()[i];
          }
      }
   }

    //loop over edges looking for any that are too short as a result of the
    //displacement
    label nShortEdges = 0;
    forAll(pp_.edges(),eI)
    {
        vector pt0 = pp_.localPoints()[pp_.edges()[eI][0]];
        vector pt1 = pp_.localPoints()[pp_.edges()[eI][1]];
        vector npt0 = pt0 + pointDisp_[pp_.edges()[eI][0]];
        vector npt1 = pt1 + pointDisp_[pp_.edges()[eI][1]];

        //if the edge length has shrunk by more than a couple of
        //orders of magnitude, undo the snapping
        if (mag(npt1-npt0) < maxShrinkFactor_*mag(pt1-pt0))
        {
            nShortEdges++;
            pointDisp_[pp_.edges()[eI][0]] = vector::zero;
            pointDisp_[pp_.edges()[eI][1]] = vector::zero;
        }
    }

    if (nShortEdges && verbose_)
    {
        Pout<< " Reset " << nShortEdges << " edges due to excessive shrinking"
             << endl;
    }

    //loop over faces looking for any that are too skewed as a result of the displacement
    label nSkewedFaces = 0;
/*
    forAll(pp_.localFaces(),fI)
    {
        List<vector> points(pp_.localFaces()[fI].size());
        List<vector> npoints(points.size());
        forAll(pp_.localFaces()[fI],j)
        {
            points[j] = pp_.localPoints()[pp_.localFaces()[fI][j]];
            npoints[j] = points[j] + pointDisp_[pp_.localFaces()[fI][j]];
        }

        if
        (
            getMinFaceAngle(npoints) < maxShrinkFactor_*getMinFaceAngle(points)
         )
        {
            nSkewedFaces++;
            forAll(pp_.localFaces()[fI],j)
            {
                pointDisp_[pp_.localFaces()[fI][j]] = vector::zero;
            }
        }
    }
*/
    if (nSkewedFaces && verbose_)
    {
        Pout<< " Reset " << nSkewedFaces << " faces due to excessive  "
             << "shrinking of the minimum internal angle" << endl;
    }

}

Foam::scalar Foam::featureLineSnapper::getMinFaceAngle
(
    const List<vector> &points
)
{
    List<vector> unitEdges(points.size());

    forAll(points,j)
    {
        unitEdges[j] = points[(j+1)%points.size()] - points[j];
        unitEdges[j] /= (mag(unitEdges[j]) + SMALL);
    }
    //find smallest dp
    scalar minDp = 1.0;
    forAll(unitEdges,j)
    {
        scalar dp = unitEdges[j] & unitEdges[(j+1)%unitEdges.size()];
        if (dp < minDp)
        {
            minDp = max(dp, -1.0);
        }
    }
    return std::acos(-minDp);
}

/**
 * Experimental: apply a smidgen of Laplacian smoothing to non-chain points in
 * the displacement field just to see what happens. (Note - propagation across
 * parallel boundaries won't be quite right because there are no halos)
 */
void Foam::featureLineSnapper::smoothDisplacementField()
{
    scalar relax = 1.0;
    label numIters = 1;

    const polyBoundaryMesh& patches = meshRefiner_.mesh().boundaryMesh();

    forAll(patches, patchI)
    {
        if
        (
            isA<processorPolyPatch>(patches[patchI])
            && patches[patchI].nPoints() > 0
        )
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(patches[patchI]);

            //dump out the patch to file
            OStringStream name;
            name << "processorPatch." << procPatch.neighbProcNo() << ".vtk";
            simpleVTKWriter
            (
                procPatch.localFaces(),
                procPatch.localPoints()
             ).write(name.str());
        }
    }

    autoPtr<PackedList<1>> isChainNode = flagChainNodes();

    for (label dummy = 0; dummy < numIters; ++dummy)
    {
        List<vector> newPointDisp(pointDisp_.size(),vector::zero);
        //we'll use a Jacobi-type scheme so we're order-independent

        forAll(pp_.pointEdges(),pI)
        {
        //average neighbouring displacements if not on a chain
            if (!isChainNode()[pI])
            {
                forAll(pp_.pointEdges()[pI],j)
                {
                    newPointDisp[pI] +=
                        pointDisp_
                        [
                            pp_.edges()[pp_.pointEdges()[pI][j]].otherVertex(pI)
                        ];
                }
                newPointDisp[pI] /= pp_.pointEdges()[pI].size();
            }
            //fix displacement for chain nodes
            else
            {
                newPointDisp[pI] = pointDisp_[pI];
            }
        }

        //under-relax contributions
        forAll(pointDisp_,j)
        {
            pointDisp_[j] = relax*newPointDisp[j] + (1.0-relax)*pointDisp_[j];
        }

        //transmit points across partition boundaries
        syncTools::syncPointList
        (
            meshRefiner_.mesh(),
            pp_.meshPoints(),
            pointDisp_,
            maxMagEqOp(),
            vector::zero
        );
    }
}


/**
 * Sets up the topological changes to the mesh required to introduce edges
 * where the shadow chains cross faces diagonally.
 */
void Foam::featureLineSnapper::addChainDiagonalEdges
(
    const PackedList<1>& isZonedFace
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    DynamicList<Tuple2<label, face>> origFaces;
    polyTopoChange firstMeshMod(mesh);

    //loop through each chain finding where it traverses a face diagonal
    forAll(shadowChains_,scI)
    {
        const labelList &chain = shadowChains_[scI];

        for (label pI=0; pI < chain.size()-1; ++pI)
        {
            //chain is stored as local points
            label lPt1 = chain[pI];
            label lPt2 = chain[pI+1];

            //Eliminate cases where we're just traversing an edge
            if (allowQuadSplits_ && !isEdge(lPt1,lPt2))
            {
                //For required splits, work out what face it's in and return
                //the two new faces
                label origID = -1;
                label origLocalID = -1;
                face origNew;
                face newNew;

                bool displacedEdge = true;

                if
                (
                    mag(pointDisp_[lPt1]) < SMALL
                    || mag(pointDisp_[lPt2]) < SMALL
                )
                {
                    displacedEdge = false;
                }

                if
                (
                    displacedEdge &&
                    getFaceSplit
                    (
                        lPt1,
                        lPt2,
                        origLocalID,

                        origNew,
                        newNew
                    )
                )
                {
                    if
                    (
                        !faceIsSplit[origLocalID]
                        && !isZonedFace.get(pp_.addressing()[origLocalID])
                    )
                    {
                        faceIsSplit[origLocalID] = 1;

                        //now convert everything to mesh numbering

                        face origNewLocal(origNew);
                        face newNewLocal(newNew);

                        //map face node numbers to mesh numbering:
                        forAll(origNew, i)
                        {
                            origNew[i] = pp_.meshPoints()[origNew[i]];
                        }

                        forAll(newNew, i)
                        {
                            newNew[i] = pp_.meshPoints()[newNew[i]];
                        }

                        origID = pp_.addressing()[origLocalID];

                        scalar faceArea =
                            mag(origNew.areaNormal(mesh.points()));
                        scalar areaRatio = faceArea /
                            ( mesh.magFaceAreas()[origID] + SMALL);

                        if (areaRatio >= 0.25 && areaRatio < 0.75)
                        {
                            //now go ahead and set up the mesh change
                            label owner = mesh.faceOwner()[origID];
                            label patchI =
                                mesh.boundaryMesh().whichPatch(origID);
                            label zoneID = mesh.faceZones().whichZone(origID);

                            const face f =  mesh.faces()[origID];
                            origFaces.append(Tuple2<label, face>(origID, f));

                            firstMeshMod.modifyFace
                            (
                                origNew,
                                origID,
                                owner,
                                -1,
                                false,
                                patchI,
                                zoneID,
                                false
                            );
                            firstMeshMod.addFace
                            (
                                newNew,
                                owner,
                                -1,
                                -1,
                                -1,
                                origID,
                                false,
                                patchI,
                                zoneID,
                                false
                            );
                        }
                    }
                }
            }
        }
    }

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> mapSplit = firstMeshMod.changeMesh
    (
        mesh,
        false,
        true
    );
    updateMesh(mapSplit());

    faceSet errorFaces
    (
        mesh,
        "errorFaces",
        mesh.nFaces()-mesh.nInternalFaces()
    );

    bool hasErrors = motionSmoother::checkMesh
    (
        false,  // report
        mesh,
        motionDict_,
        errorFaces
    );

    polyTopoChange meshMod(mesh);

    if (hasErrors)
    {
        forAll(origFaces, faceI)
        {
            label masterFaceI = origFaces[faceI].first();
            label newFaceI = mapSplit().reverseFaceMap()[masterFaceI];
            label owner = mesh.faceOwner()[newFaceI];
            label patchI = mesh.boundaryMesh().whichPatch(newFaceI);
            label zoneID = mesh.faceZones().whichZone(newFaceI);

            const cell& cFaces = mesh.cells()[owner];

            bool errorCell = false;

            forAll(cFaces, i)
            {
                if (errorFaces.found(cFaces[i]))
                {
                    errorCell = true;
                    break;
                }
            }

            if (errorCell)
            {
                face newFace = origFaces[faceI].second();
                forAll(newFace, i)
                {
                    newFace[i] = mapSplit().reversePointMap()[newFace[i]];
                }

                meshMod.modifyFace
                (
                    newFace,
                    newFaceI,
                    owner,
                    -1,
                    false,
                    patchI,
                    zoneID,
                    false
                 );

                forAll(cFaces, i)
                {
                    if
                    (
                        cFaces[i] != newFaceI
                        && mapSplit().faceMap()[cFaces[i]] == masterFaceI
                     )
                    {
                        meshMod.removeFace(cFaces[i], -1);
                        break;
                    }
                }
            }
        }
    }

    autoPtr<mapPolyMesh> mapReverse = meshMod.changeMesh
    (
        meshRefiner_.mesh(),
        false,
        true
    );

    updateMesh(mapReverse());
}

bool Foam::featureLineSnapper::isEdge(label pt1, label pt2) const
{
    //do the two points form an edge?
    labelList edges1 = pp_.pointEdges()[pt1];
    forAll(edges1,eI)
    {
        if (pp_.edges()[edges1[eI]].otherVertex(pt1) == pt2)
        {
            return true;
        }
    }
    return false;
}

/**
 * Given two points on a face and the face itself, works out the two triangles
 * that will be produced by introducing an edge from pt1 to pt2.  Silently
 * ignores point pairs that don't actually lie in the same face.  Complains if
 * points aren't diagonally opposite.  Checks that the two new faces won't be
 * highly distorted under the intended displacement.  If they would be, zeros
 * the displacement (this mimics what is done for the patch globally in
 * filterDisplacement).
 */
bool Foam::featureLineSnapper::getFaceSplit
(
    label pt1,
    label pt2,
    label &origFaceID,
    face &origNew,
    face &newNew
) const
{
    //find the face they have in common
    labelList faces1 = pp_.pointFaces()[pt1];
    labelList faces2 = pp_.pointFaces()[pt2];
    std::list<label> intx;
    std::set_intersection
    (
        faces1.begin(),
        faces1.end(),
        faces2.begin(),
        faces2.end(),
        std::back_inserter(intx)
    );
    label numCommon = intx.size();

    if (numCommon == 1)
    {
        origFaceID = intx.front();

        const face &origFace = pp_.localFaces()[origFaceID];

        if (origFace.size() > 3)
        {
            //search through the face for the two points
            label ind1 = label(std::find(origFace.begin(),origFace.end(),pt1)
                               - origFace.begin());
            label ind2 = label(std::find(origFace.begin(),origFace.end(),pt2)
                               - origFace.begin());


            label find1 = origFace.fcIndex(ind1);
            label rind1 = origFace.rcIndex(ind1);

            if (find1 == ind2 || rind1 == ind2)
            {
                Pout<< "WARNING: featureLineSnapper::getFaceSplit: Supplied "
                     << "nodes don't lie diagonally opposite (or face is "
                     << "<not ordered).\n";
                Pout<< " indices are: " << ind1 << "," << ind2 << "\n";
                return false;
            }

            DynamicList<label> origNewTemp(origFace.size());
            label count = ind1;
            origNewTemp.append(count % origFace.size());

            while ((count % origFace.size()) != (ind2 % origFace.size()))
            {
               count = origFace.fcIndex(count);
               origNewTemp.append(count % origFace.size());
            }

            DynamicList<label> newNewTemp(origFace.size());
            count = ind2;
            newNewTemp.append(count % origFace.size());

            while ((count % origFace.size()) != (ind1 % origFace.size()))
            {
               count = origFace.fcIndex(count);
               newNewTemp.append(count % origFace.size());
            }

            origNewTemp.shrink();
            newNewTemp.shrink();

            origNew.setSize(origNewTemp.size());
            forAll(origNew, i)
            {
                origNew[i] = origFace[origNewTemp[i]];
            }
            newNew.setSize(newNewTemp.size());
            forAll(newNew, i)
            {
                newNew[i] = origFace[newNewTemp[i]];
            }

            //check neither face is going to get too thin under the displacement
            List<vector> points(origNew.size());
            List<vector> npoints(points.size());
            forAll(points, i)
            {
                points[i] = pp_.localPoints()[origNew[i]];
                npoints[i] = points[i] + pointDisp_[origNew[i]];
            }

            if
            (
                getMinFaceAngle(points) < minFaceAngle_
                || (getMinFaceAngle(npoints) <
                    maxShrinkFactor_*getMinFaceAngle(points))
            )
            {
                return false;
            }

            List<vector> pointsNew(newNew.size());
            List<vector> npointsNew(pointsNew.size());
            forAll(pointsNew, i)
            {
                pointsNew[i] = pp_.localPoints()[newNew[i]];
                npointsNew[i] = pointsNew[i] + pointDisp_[newNew[i]];
            }

            if
            (
                getMinFaceAngle(pointsNew) < minFaceAngle_
                || (getMinFaceAngle(npointsNew) <
                    maxShrinkFactor_*getMinFaceAngle(pointsNew))
            )
            {
                return false;
            }
        }
        else
        {
            if (verbose_)
            {
                //Issue a warning
                Pout<< "WARNING: featureLineSnapper::getFaceSplit: Face "
                     << "to be split is not a quad!\n";
            }
            return false;
        }

        return true;
    }
    else if (numCommon > 1)
    {
        if (verbose_)
        {
            //Issue a warning - split ambiguous.
            Pout<< "WARNING: featureLineSnapper::getFaceSplit: Supplied "
                 << "nodes have more than one common face.\n";
        }
    }

    return false;
}

/**
 * Loops through all the points in shadow chains, moving them to their target
 * locations and flagging them as stationary (for the purposes of further
 * feature snapping) if their displacement is non-zero.  This function is
 * designed to be called from within autoSnapDriver::featureSnap.
 */
void Foam::featureLineSnapper::snapShadowChainsToLines
(
    pointField &newPoints,
    boolList &stationaryPoints
) const
{
    forAll(shadowChains_,scI)
    {
        const labelList &chain = shadowChains_[scI];

        forAll(chain,pI)
        {
            newPoints[pp_.meshPoints()[chain[pI]]] =
                pp_.localPoints()[chain[pI]] + pointDisp_[chain[pI]];
            //only fix points if they were actually snapped - some chain
            //points' displacements are zeroed because they are spurious

            if (mag(pointDisp_[chain[pI]]) > SMALL)
            {
                stationaryPoints[pp_.meshPoints()[chain[pI]]] = true;
            }
        }
    }
}

void Foam::featureLineSnapper::updateMesh
(
    const mapPolyMesh& map
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    // Update fields
    mesh.updateMesh(map);

    snappySnapDriver::updateBaffles(map, baffles_);

    meshRefiner_.updateMesh(map,labelList(0));

    //modify the surface patch so it sees the new faces and edges

    // save the old mesh points map
    labelList ppPointMap(pp_.meshPoints().begin(),pp_.meshPoints().end());
    pp_.clearOut();
    labelList newAddressing;
    meshRefinement::calcPatchAddressing(mesh,adaptPatchIDs_,newAddressing);
    pp_.resetAddressing(newAddressing);

    //create mapping from old mesh numbers to new mesh numbers

    //now renumber the mesh point map so it becomes an old-to-new local
    //point mapping
    forAll(ppPointMap,pI)
    {
        Map<label>::const_iterator findEntry =
            pp_.meshPointMap().find(ppPointMap[pI]);
        if (findEntry != pp_.meshPointMap().end())
        {
            ppPointMap[pI] = *findEntry;
        }
        else
        {
            ppPointMap[pI] = -1; //node has been removed from the patch
        }
    }

    //reset the shadow chains according to the new local numbering
    renumber(ppPointMap);
}


/**
 * Reorders the shadow chains and the displacement field according to the
 * supplied change in the patch node numbering
 */
void Foam::featureLineSnapper::renumber(const labelList &localPointMap)
{
    //check we don't have an empty patch!  (It's quite possible in parallel)
    if (localPointMap.size() > 0)
    {
        label maxNewPtNum = -1;
        forAll(localPointMap,j)
        {
            if (localPointMap[j] > maxNewPtNum)
                maxNewPtNum = localPointMap[j];
        }
        pointField newDisp(maxNewPtNum+1,vector::zero);

        forAll(shadowChains_,scI)
        {
            const labelList &chain = shadowChains_[scI];

            //check whether any nodes have been deleted from the chain
            label count = 0;
            forAll(chain,pI)
            {
                if (localPointMap[chain[pI]] >= 0)
                {
                    count++;
                }
            }
            labelList newChain(count);

            count = 0;
            forAll(chain,pI)
            {
                if (localPointMap[chain[pI]] >= 0)
                {
                    newChain[count++] = localPointMap[chain[pI]];
                    newDisp[localPointMap[chain[pI]]] = pointDisp_[chain[pI]];
                }
            }
            shadowChains_[scI].transfer(newChain);
        }
        pointDisp_.transfer(newDisp);
    }
}

/**
 * Returns a local-node-based list of which nodes are in feature lines
 */
Foam::autoPtr<Foam::PackedList<1>>
Foam::featureLineSnapper::flagChainNodes() const
{
    autoPtr<PackedList<1>> flagList
    (
        new PackedList<1>(pp_.localPoints().size(),0)
    );

    forAll(shadowChains_,scI)
    {
        forAll(shadowChains_[scI],pI)
        {
            flagList()[shadowChains_[scI][pI]] = 1;
        }
    }

    return flagList;
}

/**
 * Returns the shadow chain node lists in mesh
 * (as opposed to their native patch) numbering
 */
Foam::List<labelList> Foam::featureLineSnapper::getMeshNumberedShadowChains()
const
{
    List<labelList> sc;
    DynamicList<labelList> meshChains;
    forAll(shadowChains_,scI)
    {
        if (shadowChains_[scI].size() > 0)
        {
            meshChains.append(labelList());
            meshChains[meshChains.size()-1].setSize(shadowChains_[scI].size());
            forAll(shadowChains_[scI],pI)
            {
                meshChains[meshChains.size()-1][pI] =
                    pp_.meshPoints()[shadowChains_[scI][pI]];
            }
        }
    }

    sc.transfer(meshChains);
    return sc;
}
