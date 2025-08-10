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

#include "featureLineSnapping/featureLinePrep.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "indexedOctree/treeDataFace.H"
#include "meshRefinement/meshRefinement.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "refinementFeatures/refinementFeatures.H"
#include "meshes/meshTools/simpleVTKWriter.H"

/**
 * Builds (if necessary) and returns an edge list for the line by walking
 * along the point list.
 */
const Foam::edgeList &Foam::featureLinePrep::featureLine::edges() const
{
    if (!edgeList_.valid())
    {
        edgeList_.reset(new List<edge>(this->size()-1));
        forAll(edgeList_(),j)
        {
            edgeList_()[j][0] = (*this)[j];
            edgeList_()[j][1] = (*this)[j+1];
        }
    }
    return edgeList_();
}



Foam::featureLinePrep::featureLinePrep
(
    const meshRefinement& meshRefiner,
    const indirectPrimitivePatch& pp,
    const dictionary& dict,
    const boolList& regionSnap
):
    meshRefiner_(meshRefiner),
    pp_(pp),
    regionSnap_(regionSnap),
    regionFLs_(false),
    geometryFLs_(true),
    verbose_(false),
    dump_(false)
{
    //Check if any region snap patches are active
    forAll(regionSnap, regionI)
    {
        if (regionSnap[regionI])
        {
            regionFLs_ = true;
            break;
        }
    }

    if (dict.found("geometryFeatureLines"))
    {
       geometryFLs_ = readBool(dict.lookup("geometryFeatureLines"));
    }
    if (dict.found("verboseFeatureLines"))
    {
       verbose_ = readBool(dict.lookup("verboseFeatureLines"));
    }
    if (dict.found("dumpFeatureLineData"))
    {
       dump_ = readBool(dict.lookup("dumpFeatureLineData"));
    }
}

/**
 * Deletes any trees and arc length distributions that may have been built
 */
Foam::featureLinePrep::~featureLinePrep()
{
    forAll(flParams_,flI)
    {
        delete flParams_[flI];
    }
    forAll(flTrees_,flI)
    {
        delete flTrees_[flI];
    }
}

/**
 * Adds feature lines to the central list by extracting lines marking the
 * boundaries between regions in the surface patch
 */
void Foam::featureLinePrep::extractRegionBdries()
{
    const refinementSurfaces& refineSurfaces = meshRefiner_.surfaces();
    const labelList& surf = refineSurfaces.surfaces();

    label sz = 0;
    forAll(surf, surfI)
    {
        const searchableSurface& geom = refineSurfaces.geometry()[surf[surfI]];
        if (isA<triSurfaceMesh>(geom))
        {
            sz++;
        }
    }

    regionFeatures_.setSize(sz);

    sz = 0;

    forAll(surf, surfI)
    {
        const searchableSurface& geom = refineSurfaces.geometry()[surf[surfI]];

        if (isA<triSurfaceMesh>(geom))
        {
            const triSurfaceMesh& triMesh = refCast<const triSurfaceMesh>(geom);
            const triSurface& s = static_cast<const triSurface&>(triMesh);

            DynamicList<edge> featureEdges(s.nPoints());
            // Construct pointFaces. Let's hope surface has compact point
            // numbering ...
            labelListList pointFaces;
            invertManyToMany
            (
                s.localPoints().size(),
                s.localFaces(),
                pointFaces
            );

            forAll(pointFaces, pointI)
            {
                const labelList& pFaces = pointFaces[pointI];
                boolList visited(pFaces.size(), false);

                forAll(pFaces, i)
                {
                    const labelledTri& f = s.localFaces()[pFaces[i]];
                    const label& region = s.localFaces()[pFaces[i]].region();
                    label globalRegion =
                        refineSurfaces.globalRegion(surfI, region);

                    // Forward edge
                    label fp = findIndex(f, pointI);
                    visited[i] = true;

                    label nextPointI = f[f.fcIndex(fp)];

                    if (nextPointI > pointI)
                    {
                        bool foundNbr = false;
                        forAll(pFaces, j)
                        {
                            const labelledTri& fnbr =
                                s.localFaces()[pFaces[j]];

                            label fpnbr = findIndex(fnbr, pointI);
                            label nPointI = fnbr[fnbr.fcIndex(fpnbr)];
                            label pPointI = fnbr[fnbr.rcIndex(fpnbr)];

                            if
                            (
                                nPointI == nextPointI
                                || pPointI == nextPointI
                            )
                            {
                                if (pFaces[i]!= pFaces[j])
                                {
                                    foundNbr = true;
                                }

                                if (!visited[j])
                                {
                                    const label& nRegion =
                                        s.localFaces()[pFaces[j]].region();
                                    label globalNbrRegion =
                                        refineSurfaces.globalRegion
                                        (
                                            surfI,
                                            nRegion
                                        );

                                    bool foundRegion = false;
                                    if
                                    (
                                        regionSnap_[globalNbrRegion]
                                        || regionSnap_[globalRegion]
                                    )
                                    {
                                        foundRegion = true;
                                    }

                                    if (region != nRegion && foundRegion)
                                    {
                                        const edge e(pointI, nextPointI);
                                        featureEdges.append(e);
                                    }
                                    break;
                                }
                            }
                        }
                        if
                        (
                            !foundNbr
                            && regionSnap_[globalRegion]
                        )
                        {
                            const edge e(pointI, nextPointI);
                            featureEdges.append(e);
                        }
                    }

                    // Reverse edge
                    label prevPointI = f[f.rcIndex(fp)];

                    if (prevPointI > pointI)
                    {
                        bool foundNbr = false;
                        forAll(pFaces, j)
                        {
                            const labelledTri& fnbr =
                                s.localFaces()[pFaces[j]];
                            label fpnbr = findIndex(fnbr, pointI);
                            label nPointI = fnbr[fnbr.fcIndex(fpnbr)];
                            label pPointI = fnbr[fnbr.rcIndex(fpnbr)];

                            if
                            (
                                nPointI == prevPointI
                                || pPointI == prevPointI
                             )
                            {
                                if (pFaces[i]!= pFaces[j])
                                {
                                    foundNbr = true;
                                }

                                if (!visited[j])
                                {
                                    const label& nRegion =
                                        s.localFaces()[pFaces[j]].region();
                                    label globalNbrRegion =
                                        refineSurfaces.globalRegion
                                        (
                                            surfI,
                                            nRegion
                                        );

                                    bool foundRegion = false;
                                    if
                                    (
                                        regionSnap_[globalNbrRegion]
                                        || regionSnap_[globalRegion]
                                    )
                                    {
                                        foundRegion = true;
                                    }

                                    if (region != nRegion && foundRegion)
                                    {
                                        const edge e(pointI, prevPointI);
                                        featureEdges.append(e);
                                    }
                                    break;
                                }
                            }
                        }
                        if
                        (
                            !foundNbr
                            && regionSnap_[globalRegion]
                        )
                        {
                            const edge e(pointI, prevPointI);
                            featureEdges.append(e);
                        }
                    }
                }
            }
            featureEdges.shrink();
            label nFeatureEdges = featureEdges.size();

            Map<label> pointMap(2*nFeatureEdges);
            pointField pts(2*nFeatureEdges);
            edgeList edg(nFeatureEdges);
            label nPoints = 0;
            label nEdges = 0;
            forAll(featureEdges, i)
            {
                edge e = featureEdges[i];
                label v0 = e[0];
                label v1 = e[1];
                edg[nEdges] = e;
                nEdges++;
                if (pointMap.insert(v0,nPoints))
                {
                    pts[nPoints] = s.localPoints()[v0];
                    nPoints++;
                }
                if (pointMap.insert(v1,nPoints))
                {
                    pts[nPoints] = s.localPoints()[v1];
                    nPoints++;
                }
            }
            pts.setSize(nPoints);
            edg.setSize(nEdges);

            forAll(edg, edgeI)
            {
                edg[edgeI][0] = pointMap[edg[edgeI][0]];
                edg[edgeI][1] = pointMap[edg[edgeI][1]];
            }

            regionFeatures_.set
            (
                sz,
                new edgeMesh
                (
                    edgeMesh(pts, edg)
                )
             );
            sz++;
        }
    }
    regionFeatures_.setSize(sz);

    if (Pstream::parRun())
    {
        Random rndGen(653213);

        // Determine mesh bounding boxes:
        List<List<treeBoundBox>> meshBb(Pstream::nProcs());
        {
            meshBb[Pstream::myProcNo()] = List<treeBoundBox>
            (
                1,
                treeBoundBox
                (
                    boundBox(meshRefiner_.mesh().points(), false)
                 ).extend(rndGen, 1E-3)
             );
            Pstream::allGatherList(meshBb);
        }

        forAll(regionFeatures_, fI)
        {
            regionFeatures_[fI].trim(meshBb[Pstream::myProcNo()]);
        }
    }

    forAll(regionFeatures_, featI)
    {
        const edgeMesh& rMesh = regionFeatures_[featI];
        const edgeList& rEdges = rMesh.edges();
        const pointField& rPoints = rMesh.points();

        const labelListList& pointEdges = rMesh.pointEdges();

        // Whether edge has been visited.
        PackedList<1> edgeVisited(rEdges.size(), 0u);
        PackedList<1> visited(rEdges.size(), 0u);

        forAll(pointEdges, pointI)
        {
            if (pointEdges[pointI].size() != 2)
            {
                forAll(pointEdges[pointI], i)
                {
                    label edgeI = pointEdges[pointI][i];

                    if (edgeVisited.set(edgeI, 1u))
                    {
                        PackedList<1> pointVisited(rPoints.size(), 0u);
                        DynamicList<label> meshPoints(rPoints.size());

                        // Unvisited edge. Make the particle go to
                        // the other point on the edge.
                        const edge& e = rEdges[edgeI];
                        label otherPointI = e.otherVertex(pointI);
                        meshPoints.append(pointI);
                        meshPoints.append(otherPointI);

                        pointVisited.set(pointI, 1u);
                        pointVisited.set(otherPointI, 1u);

                        manifoldFeatures::manifoldEdgeWave
                        (
                            pointI,
                            otherPointI,
                            pointEdges,
                            rMesh,
                            edgeVisited,
                            pointVisited,
                            meshPoints,
                            labelMax
                        );
                        meshPoints.shrink();

                        featureLines_.append
                        (
                            featureLine
                            (
                                rPoints,
                                -2.0
                             )
                         );
                        //Now construct a single feature line
                        featureLine &baseLine =
                            featureLines_[featureLines_.size()-1];
                        baseLine.setSize(meshPoints.size());
                        forAll(meshPoints,j)
                        {
                            baseLine[j] = meshPoints[j];
                        }
                    }
                }
            }
        }

        //start from a feature point
        forAll(pointEdges, pointI)
        {
            forAll(pointEdges[pointI], i)
            {
                label edgeI = pointEdges[pointI][i];

                bool foundFeat = false;
                if (edgeVisited.get(edgeI) == 0u)
                {
                    const edge& e = rEdges[edgeI];
                    label otherPointI = e.otherVertex(pointI);

                    forAll(pointEdges[pointI], j)
                    {
                        label edgeJ = pointEdges[pointI][j];
                        if (edgeJ != edgeI)
                        {
                            const edge& ej = rEdges[edgeJ];
                            label otherPointJ = ej.otherVertex(pointI);

                            vector vecA = rPoints[otherPointI]-rPoints[pointI];
                            vecA /= (mag(vecA) + SMALL);

                            vector vecB = rPoints[pointI]-rPoints[otherPointJ];
                            vecB /= (mag(vecB) + SMALL);
                            if ((vecA & vecB) < 0.9659)
                            {
                                foundFeat = true;
                                break;
                            }
                        }
                    }
                }

                if (foundFeat && edgeVisited.set(edgeI, 1u))
                {
                    PackedList<1> pointVisited(rPoints.size(), 0u);
                    DynamicList<label> meshPoints(rPoints.size());
                    // Unvisited edge. Make the particle go to
                    // the other point on the edge.
                    const edge& e = rEdges[edgeI];
                    label otherPointI = e.otherVertex(pointI);
                    meshPoints.append(pointI);
                    meshPoints.append(otherPointI);

                    pointVisited.set(pointI, 1u);
                    pointVisited.set(otherPointI, 1u);

                    manifoldFeatures::manifoldEdgeWave
                    (
                        pointI,
                        otherPointI,
                        pointEdges,
                        rMesh,
                        edgeVisited,
                        pointVisited,
                        meshPoints,
                        labelMax
                    );
                    meshPoints.shrink();

                    featureLines_.append
                    (
                        featureLine
                        (
                            rPoints,
                            -2.0
                         )
                     );
                    //Now construct a single feature line
                    featureLine &baseLine =
                        featureLines_[featureLines_.size()-1];
                    baseLine.setSize(meshPoints.size());
                    forAll(meshPoints,j)
                    {
                        baseLine[j] = meshPoints[j];
                    }
                }
            }
        }

        forAll(pointEdges, pointI)
        {
            forAll(pointEdges[pointI], i)
            {
                label edgeI = pointEdges[pointI][i];

                if (edgeVisited.set(edgeI, 1u))
                {
                    PackedList<1> pointVisited(rPoints.size(), 0u);
                    DynamicList<label> meshPoints(rPoints.size());
                    // Unvisited edge. Make the particle go to
                    // the other point on the edge.
                    const edge& e = rEdges[edgeI];
                    label otherPointI = e.otherVertex(pointI);
                    meshPoints.append(pointI);
                    meshPoints.append(otherPointI);

                    pointVisited.set(pointI, 1u);
                    pointVisited.set(otherPointI, 1u);

                    manifoldFeatures::manifoldEdgeWave
                    (
                        pointI,
                        otherPointI,
                        pointEdges,
                        rMesh,
                        edgeVisited,
                        pointVisited,
                        meshPoints,
                        labelMax
                    );
                    meshPoints.shrink();

                    featureLines_.append
                    (
                        featureLine
                        (
                            rPoints,
                            -2.0
                         )
                     );
                    //Now construct a single feature line
                    featureLine &baseLine =
                        featureLines_[featureLines_.size()-1];
                    baseLine.setSize(meshPoints.size());
                    forAll(meshPoints,j)
                    {
                        baseLine[j] = meshPoints[j];
                    }
                }
            }
        }
    }

    if (dump_)
    {
        simpleVTKWriter writer;
        forAll(featureLines_,flI)
        {
           writer.addMorePoints(featureLines_[flI],featureLines_[flI].points());
        }
        writer.write("regionFeatureLines.vtk");
    }
}


/**
 * Takes feature lines constructed during feature refinement.
 */
void Foam::featureLinePrep::mapFeatureMeshes
(
    const refinementFeatures& rFeatures,
    PtrList<Tuple2<labelList, label>>& featureMeshes
)
{
    const label minValidLength = 2;
    forAll(featureMeshes,i)
    {
        const labelList &meshPoints = featureMeshes[i].first();
        const label origFeatI = featureMeshes[i].second();

        //only consider substantial lines
        if (meshPoints.size() >= minValidLength)
        {
            featureLines_.append
            (
                featureLine
                (
                    rFeatures[origFeatI].points(),
                    -1.0
                 )
            );
            //Now construct a single feature line
            featureLine &baseLine = featureLines_[featureLines_.size()-1];
            baseLine.setSize(meshPoints.size());
            forAll(meshPoints,j)
            {
                baseLine[j] = meshPoints[j];
            }
        }
    }

    if (dump_)
    {
        simpleVTKWriter writer;
        forAll(featureLines_,flI)
        {
            writer.addMorePoints
            (
                featureLines_[flI],
                featureLines_[flI].points()
            );
        }
        writer.write("allFeatureLines.vtk");
    }
}


/**
 * Controlling algorithm for feature line construction and distribution.
 * Currently does the following:
 * -# Extracts region boundaries
 * -# Distributes the global feature lines it produces
 * -# Extracts feature lines from featureEdgeMeshes (these come out distributed)
 * -# Trims lines if they extend outside the global bounding box (e.g. if mesh
 *     is built as half-model on full-model geometry)
 */
void Foam::featureLinePrep::prepareFeatureLines
(
    const refinementFeatures& rFeatures,
    PtrList<Tuple2<labelList, label>>& featureMeshes
)
{
    Info<< "Preparing feature lines for direct snapping" << endl;
    Info<< "-------------------------------------------" << endl << endl;

    if (regionFLs_)
    {
       Info<< " Extracting global region boundaries" << endl;
       extractRegionBdries();
    }

    if (geometryFLs_)
    {
       Info<< " Converting refinement feature meshes to feature lines" << endl;
       mapFeatureMeshes(rFeatures, featureMeshes);
    }

    flParams_.setSize(featureLines_.size(),0);
    flTrees_.setSize(featureLines_.size(),0);

    Info<< endl;
}


/**
 * Builds patch's face search tree if required
 */
const Foam::indexedOctree<Foam::treeDataFace>
&Foam::featureLinePrep::getPatchTree() const
{
    if (!patchTree_.valid())
    {
        treeBoundBox bb(pp_.localPoints());
        patchTree_.reset
        (
            new indexedOctree<treeDataFace>
            (
                treeDataFace
                (
                    false,
                    meshRefiner_.mesh(),
                    pp_.addressing()
                ),
                bb,
                10,      // maxLevel
                1,     // leafsize
                3.0   // duplicity - no idea
            )
        );
    }
    return patchTree_();
}

/**
 * Constructs (if necessary) and returns the arc length distribution
 * of the feature line
 */
const Foam::List<Foam::scalar> &
Foam::featureLinePrep::getFLParams(label flI) const
{
    if (flParams_[flI] == 0)
    {
        flParams_[flI] = new List<scalar>(featureLines_[flI].size());

        //walk along the line calculating its cumulative length at each node
        (*flParams_[flI])[0] = 0.0;
        for (int pI = 0; pI < int(featureLines_[flI].size())-1; ++pI)
        {
            (*flParams_[flI])[pI+1] =
                (
                    *flParams_[flI])[pI]
                      + mag(featureLines_[flI].point(pI+1)
                      -  featureLines_[flI].point(pI)
                );
        }
    }
    return *flParams_[flI];
}

/**
 * Constructs (if necessary) and returns the edge search tree for the
 * feature line
 */
const Foam::indexedOctree<Foam::treeDataEdge> &
Foam::featureLinePrep::getFLSearchTree(label flI) const
{
    if (flTrees_[flI] == 0)
    {
        flTrees_[flI] = createFeatureLineSearchTree(flI);
    }
    return *flTrees_[flI];
}

/**
 * Builds a feature line edge search tree
 */
Foam::indexedOctree<Foam::treeDataEdge> *
Foam::featureLinePrep::createFeatureLineSearchTree(label flI) const
{
    const pointField &points = featureLines_[flI].points();

    // need to compute bounding box for feature line
    treeBoundBox bb(points,featureLines_[flI]);

    //ensure there are no zero dimensions because the tree construction
    //can't handle this
    scalar mySmall = bb.avgDim() < 1.0 ? SMALL : bb.avgDim()*SMALL;
    mySmall *= 1e6;
    for (label i = 0; i < 3; ++i)
    {
        if (bb.min()[i] + mySmall >= bb.max()[i])
        {
            bb.min()[i] -= mySmall;
            bb.max()[i] += mySmall;
        }
    }

    // need to find the edges comprising the feature line
    labelList flEdges(identity(featureLines_[flI].size()-1));

    //build search tree for feature line
    return new indexedOctree<treeDataEdge>
    (
        treeDataEdge
        (
            false,                     // cachebb
            featureLines_[flI].edges(),// edges
            points,                    // points
            flEdges                    // selected edges
        ),
        bb,
        10,    // maxLevel
        1,     // leafsize
        3.0    // duplicity - no idea, this is just the number I
           // copied from elsewhere!!
    );
}
