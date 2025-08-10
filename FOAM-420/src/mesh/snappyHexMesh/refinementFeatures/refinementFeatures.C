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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "refinementFeatures/refinementFeatures.H"
#include "db/Time/Time.H"
#include "primitives/Tuple2/Tuple2.H"
#include "fields/Fields/DynamicField/DynamicField.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "distributedTriSurfaceMesh/distributedTriSurfaceMesh.H"
#include "edgeMesh/featureEdgeMesh/featureEdgeMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::refinementFeatures::read
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts
)
{
   //Check for global geometry dict
    IOobject globalGeomHeader
    (
        "geometryDict",
        io.time().system(),
        io.time(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    dictionary globalGeomDict;
    if (globalGeomHeader.typeHeaderOk<IOdictionary>(true))
    {
        globalGeomDict = IOdictionary(globalGeomHeader);
    }
    else
    {
        globalGeomDict = dictionary();
    }

    forAll(featDicts, featI)
    {
        dictionary dict = featDicts[featI];
        fileName featFileName(dict.lookup("file"));

        // Try reading extendedEdgeMesh first

        IOobject extFeatObj
        (
            featFileName,                       // name
            io.time().constant(),               // instance
            "extendedFeatureEdgeMesh",          // local
            io.time(),                          // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        const fileName fName(typeFilePath<extendedFeatureEdgeMesh>(extFeatObj));

        if (!fName.empty() && extendedEdgeMesh::canRead(fName))
        {
            autoPtr<extendedEdgeMesh> eMeshPtr = extendedEdgeMesh::New
            (
                fName
            );

            Info<< "Read extendedFeatureEdgeMesh " << extFeatObj.name()
                << nl << incrIndent;
            eMeshPtr().writeStats(Info);
            Info<< decrIndent << endl;

            set(featI, new extendedFeatureEdgeMesh(extFeatObj, eMeshPtr()));
        }
        else
        {
            // Try reading edgeMesh

            IOobject featObj
            (
                featFileName,                       // name
                io.time().constant(),               // instance
                "triSurface",                       // local
                io.time(),                          // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            const fileName fName(typeFilePath<featureEdgeMesh>(featObj));

            if (fName.empty())
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Could not open " << featObj.objectPath()
                    << exit(FatalIOError);
            }


            // Read as edgeMesh
            autoPtr<edgeMesh> eMeshPtr = edgeMesh::New(fName);
            const edgeMesh& eMesh = eMeshPtr();

            Info<< "Read edgeMesh " << featObj.name() << nl
                << incrIndent;
            eMesh.writeStats(Info);
            Info<< decrIndent << endl;


            // Analyse for feature points. These are all classified as mixed
            // points for lack of anything better
            const labelListList& pointEdges = eMesh.pointEdges();

            labelList oldToNew(eMesh.points().size(), -1);
            DynamicField<point> newPoints(eMesh.points().size());
            forAll(pointEdges, pointi)
            {
                if (pointEdges[pointi].size() > 2)
                {
                    oldToNew[pointi] = newPoints.size();
                    newPoints.append(eMesh.points()[pointi]);
                }
                //else if (pointEdges[pointi].size() == 2)
                //MEJ: do something based on a feature angle?
            }
            label nFeatures = newPoints.size();
            forAll(oldToNew, pointi)
            {
                if (oldToNew[pointi] == -1)
                {
                    oldToNew[pointi] = newPoints.size();
                    newPoints.append(eMesh.points()[pointi]);
                }
            }


            const edgeList& edges = eMesh.edges();
            edgeList newEdges(edges.size());
            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];
                newEdges[edgeI] = edge
                (
                    oldToNew[e[0]],
                    oldToNew[e[1]]
                );
            }

            // Construct an extendedEdgeMesh with
            // - all points on more than 2 edges : mixed feature points
            // - all edges : external edges

            extendedEdgeMesh eeMesh
            (
                newPoints,          // pts
                newEdges,           // eds
                0,                  // (point) concaveStart
                0,                  // (point) mixedStart
                nFeatures,          // (point) nonFeatureStart
                edges.size(),       // (edge) internalStart
                edges.size(),       // (edge) flatStart
                edges.size(),       // (edge) openStart
                edges.size(),       // (edge) multipleStart
                vectorField(0),     // normals
                List<extendedEdgeMesh::sideVolumeType>(0),// normalVolumeTypes
                vectorField(0),     // edgeDirections
                labelListList(0),   // normalDirections
                labelListList(0),   // edgeNormals
                labelListList(0),   // featurePointNormals
                labelListList(0),   // featurePointEdges
                identity(newEdges.size())   // regionEdges
            );

            //Info<< "Constructed extendedFeatureEdgeMesh " << featObj.name()
            //    << nl << incrIndent;
            //eeMesh.writeStats(Info);
            //Info<< decrIndent << endl;

            set(featI, new extendedFeatureEdgeMesh(featObj, eeMesh));
        }

        extendedEdgeMesh& eMesh = operator[](featI);

        //Add global transforms
        if (globalGeomDict.found("lines"))
        {
            const dictionary& ldict = globalGeomDict.subDict("lines");
            if (ldict.found(featFileName))
            {
                const dictionary& tdict = ldict.subDict(featFileName);
                if (tdict.found("transforms"))
                {
                    List<dictionary> globalTrans(tdict.lookup("transforms"));
                    if (dict.found("transforms"))
                    {
                        List<dictionary> localTrans = dict.lookup("transforms");
                        label globalSz = globalTrans.size();
                        label newSz = globalSz + localTrans.size();
                        globalTrans.setSize(newSz);
                        forAll(localTrans, transi)
                        {
                            globalTrans[globalSz++] = localTrans[transi];
                        }
                    }
                    dict.add("transforms",globalTrans,true);
                }
            }
        }

        //Perform optional transforms
        eMesh.doTransforms(dict);

        if (dict.found("levels"))
        {
            List<Tuple2<scalar, label>> distLevels(dict["levels"]);

            if (dict.size() < 1)
            {
                FatalErrorInFunction
                    << " : levels should be at least size 1" << endl
                    << "levels : "  << dict["levels"]
                    << exit(FatalError);
            }

            distances_[featI].setSize(distLevels.size());
            levels_[featI].setSize(distLevels.size());

            forAll(distLevels, j)
            {
                distances_[featI][j] = distLevels[j].first();
                levels_[featI][j] = distLevels[j].second();

                // Check in incremental order
                if (j > 0)
                {
                    if
                    (
                        (distances_[featI][j] <= distances_[featI][j-1])
                     || (levels_[featI][j] > levels_[featI][j-1])
                    )
                    {
                        FatalErrorInFunction
                            << " : Refinement should be specified in order"
                            << " of increasing distance"
                            << " (and decreasing refinement level)." << endl
                            << "Distance:" << distances_[featI][j]
                            << " refinementLevel:" << levels_[featI][j]
                            << exit(FatalError);
                    }
                }
            }
        }
        else
        {
            // Look up 'level' for single level
            levels_[featI] = labelList(1, readLabel(dict.lookup("level")));
            distances_[featI] = scalarField(1, 0.0);
        }

        Info<< "Refinement level according to distance to "
            << featFileName << " (" << eMesh.points().size() << " points, "
            << eMesh.edges().size() << " edges)." << endl;
        forAll(levels_[featI], j)
        {
            Info<< "    level " << levels_[featI][j]
                << " for all cells within " << distances_[featI][j]
                << " metre." << endl;
        }
    }
}


// Find maximum level of a shell.
void Foam::refinementFeatures::findHigherLevel
(
    const pointField& pt,
    const label featI,
    labelList& maxLevel
)
{
    const labelList& levels = levels_[featI];

    const scalarField& distances = distances_[featI];

    // Collect all those points that have a current maxLevel less than
    // (any of) the shell. Also collect the furthest distance allowable
    // to any shell with a higher level.

    pointField candidates(pt.size());
    labelList candidateMap(pt.size());
    scalarField candidateDistSqr(pt.size());
    label candidateI = 0;

    forAll(maxLevel, pointi)
    {
        forAllReverse(levels, levelI)
        {
            if (levels[levelI] > maxLevel[pointi])
            {
                candidates[candidateI] = pt[pointi];
                candidateMap[candidateI] = pointi;
                candidateDistSqr[candidateI] = sqr(distances[levelI]);
                candidateI++;
                break;
            }
        }
    }
    candidates.setSize(candidateI);
    candidateMap.setSize(candidateI);
    candidateDistSqr.setSize(candidateI);

    // Do the expensive nearest test only for the candidate points.
    const indexedOctree<treeDataEdge>& tree = edgeTrees()[featI];

    List<pointIndexHit> nearInfo(candidates.size());
    forAll(candidates, candidateI)
    {
        nearInfo[candidateI] = tree.findNearest
        (
            candidates[candidateI],
            candidateDistSqr[candidateI]
        );
    }

    // Update maxLevel
    forAll(nearInfo, candidateI)
    {
        if (nearInfo[candidateI].hit())
        {
            // Check which level it actually is in.
            label minDistI = findLower
            (
                distances,
                mag(nearInfo[candidateI].hitPoint()-candidates[candidateI])
            );

            label pointi = candidateMap[candidateI];

            // pt is inbetween shell[minDistI] and shell[minDistI+1]
            maxLevel[pointi] = levels[minDistI+1];
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const Foam::PtrList<Foam::indexedOctree<Foam::treeDataEdge>>&
Foam::refinementFeatures::regionEdgeTrees() const
{
    if (!regionEdgeTrees_.valid())
    {
        regionEdgeTrees_.reset
        (
            new PtrList<indexedOctree<treeDataEdge>>(size())
        );
        PtrList<indexedOctree<treeDataEdge>>& trees = regionEdgeTrees_();

        forAll(*this, featI)
        {
            const extendedEdgeMesh& eMesh = operator[](featI);
            const pointField& points = eMesh.points();
            const edgeList& edges = eMesh.edges();

            // Calculate bb of all points
            treeBoundBox bb(points);

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

            trees.set
            (
                featI,
                new indexedOctree<treeDataEdge>
                (
                    treeDataEdge
                    (
                        false,                  // do not cache bb
                        edges,
                        points,
                        eMesh.regionEdges()
                    ),
                    bb,     // overall search domain
                    8,      // maxLevel
                    10,     // leafsize
                    3.0     // duplicity
                )
            );
        }
    }
    return regionEdgeTrees_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementFeatures::refinementFeatures
(
    const fvMesh& mesh,
    const PtrList<dictionary>& featDicts,
    const refinementSurfaces& refineSurfaces,
    const refinementParameters& refineParams,
    const scalar& mergeDist
)
:
    PtrList<extendedFeatureEdgeMesh>(),
    checkRefinementOnly_(false),
    refinementOnly_(),
    mFeatures_(*this,refinementOnly_,checkRefinementOnly_),
    distances_(),
    levels_(),
    edgeTrees_(0),
    pointTrees_(0),
    minCos_(-1)
{
    const labelList& surf = refineSurfaces.surfaces();

    label sz = 0;
    forAll(surf, surfI)
    {
        const searchableSurface& geom = refineSurfaces.geometry()[surf[surfI]];
        if (isA<triSurfaceMesh>(geom))
        {
            const triSurfaceMesh& triMesh = refCast<const triSurfaceMesh>(geom);
            const triSurface& s = static_cast<const triSurface&>(triMesh);
            forAll(s.patches(), patchI)
            {
                sz++;
            }
        }
    }
    sz += featDicts.size();

    this->setSize(sz);
    refinementOnly_.setSize(sz),
    distances_.setSize(sz);
    levels_.setSize(sz);
    boolList distributed(sz,false);

    scalar edge0Len = -GREAT;
    forAll(mesh.edges(), edgeI)
    {
        scalar eLen = mesh.edges()[edgeI].mag(mesh.points());
        edge0Len = max(eLen, edge0Len);
    }
    edge0Len = returnReduce(edge0Len, maxOp<scalar>());

    sz = 0;
    forAll(surf, surfI)
    {
        const bool refineBoundary = refineSurfaces.refineBoundary()[surfI];
        const searchableSurface& geom = refineSurfaces.geometry()[surf[surfI]];

        if (isA<triSurfaceMesh>(geom))
        {
            const triSurfaceMesh& triMesh = refCast<const triSurfaceMesh>(geom);

            triSurface& s = const_cast<triSurface&>
            (
                refCast<const triSurface>(triMesh)
             );

            const pointField& faceNormals = s.faceNormals();
            const pointField& faceCentres = s.faceCentres();

            EdgeMap<label> featureEdges(s.nPoints());
            label nFeatureEdges = 0;

            // Construct pointFaces. Let's hope surface has compact point
            // numbering ...
            labelListList pointFaces;
            invertManyToMany(s.localPoints().size(), s.localFaces(), pointFaces);

            forAll(pointFaces, pointi)
            {
                const labelList& pFaces = pointFaces[pointi];

                EdgeMap<label> pointFeatureEdges(100);

                boolList visited(pFaces.size(), false);

                forAll(pFaces, i)
                {
                    const labelledTri& f = s.localFaces()[pFaces[i]];
                    const label& region = s[pFaces[i]].region();
                    const label maxLevel = refineSurfaces.maxLevel
                    (
                        surfI,
                        region
                    );
                    const scalar featureAngle =
                        refineSurfaces.featureRefineAngle
                        (
                            surfI,
                            region
                        );

                    point fNorm = faceNormals[pFaces[i]];
                    fNorm /= mag(fNorm) + VSMALL;
                    // Forward edge
                    label fp = findIndex(f, pointi);
                    visited[i] = true;

                    label nextPointI = f[f.fcIndex(fp)];

                    if (nextPointI > pointi)
                    {
                        bool foundNbr = false;
                        bool manifold = true;
                        label nNext = 0;

                        forAll(pFaces, j)
                        {
                            const labelledTri& fnbr =
                                s.localFaces()[pFaces[j]];
                            label fpnbr = findIndex(fnbr, pointi);
                            label nPointI = fnbr[fnbr.fcIndex(fpnbr)];
                            label pPointI = fnbr[fnbr.rcIndex(fpnbr)];
                            if
                            (
                                nPointI == nextPointI
                                || pPointI == nextPointI
                             )
                            {
                                nNext++;
                            }
                        }
                        if (nNext > 2)
                        {
                            manifold = false;
                        }
                        if (manifold)
                        {
                            forAll(pFaces, j)
                            {
                                const labelledTri& fnbr =
                                    s.localFaces()[pFaces[j]];

                                label fpnbr = findIndex(fnbr, pointi);
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
                                        point nNorm = faceNormals[pFaces[j]];
                                        nNorm /= mag(nNorm) + VSMALL;

                                        point edgeMidPoint = 0.5*
                                            (s.localPoints()[pointi]
                                             + s.localPoints()[nextPointI]);
                                        vector eVecOwn = faceCentres[pFaces[i]]
                                            - edgeMidPoint;
                                        vector eVecNei = faceCentres[pFaces[j]]
                                            - edgeMidPoint;

                                        bool wrongOrient =  false;
                                        scalar dotProd = (fNorm & nNorm);
                                        //check for wrongly oriented surface
                                        if
                                        (
                                            dotProd < -0.9659
                                            && (eVecOwn & eVecNei) < 0
                                        )
                                        {
                                            wrongOrient = true;
                                        }

                                        const label& nRegion =
                                            s[pFaces[j]].region();
                                        const scalar nFeatureAngle =
                                                refineSurfaces.featureRefineAngle
                                                (
                                                    surfI,
                                                    nRegion
                                                 );

                                        const label nMaxLevel =
                                            refineSurfaces.maxLevel
                                            (
                                                surfI,
                                                nRegion
                                             );

                                        if (!wrongOrient)
                                        {
                                            const edge e(pointi, nextPointI);

                                            if
                                            (
                                                dotProd < featureAngle
                                                || dotProd < nFeatureAngle
                                             )
                                            {
                                                if (nMaxLevel > maxLevel)
                                                {
                                                    pointFeatureEdges.insert
                                                    (
                                                        e,
                                                        nRegion
                                                     );

                                                }
                                                else
                                                {
                                                    pointFeatureEdges.insert
                                                    (
                                                        e,
                                                        region
                                                     );
                                                }
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }

                        if (!foundNbr && refineBoundary && manifold)
                        {
                            const edge e(pointi, nextPointI);
                            pointFeatureEdges.insert
                            (
                                e,
                                region
                             );
                        }
                    }

                    // Reverse edge
                    label prevPointI = f[f.rcIndex(fp)];

                    if (prevPointI > pointi)
                    {
                        bool foundNbr = false;
                        bool manifold = true;
                        label nNext = 0;

                        forAll(pFaces, j)
                        {
                            const labelledTri& fnbr =
                                s.localFaces()[pFaces[j]];
                            label fpnbr = findIndex(fnbr, pointi);
                            label nPointI = fnbr[fnbr.fcIndex(fpnbr)];
                            label pPointI = fnbr[fnbr.rcIndex(fpnbr)];
                            if
                            (
                                nPointI == prevPointI
                                || pPointI == prevPointI
                             )
                            {
                                nNext++;
                            }
                        }
                        if (nNext > 2)
                        {
                            manifold = false;
                        }
                        if (manifold)
                        {
                            forAll(pFaces, j)
                            {
                                const labelledTri& fnbr =
                                    s.localFaces()[pFaces[j]];
                                label fpnbr = findIndex(fnbr, pointi);
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
                                        point nNorm = faceNormals[pFaces[j]];
                                        nNorm /= mag(nNorm) + VSMALL;

                                        point edgeMidPoint = 0.5*
                                            (s.localPoints()[pointi]
                                             + s.localPoints()[prevPointI]);
                                        vector eVecOwn = faceCentres[pFaces[i]]
                                            - edgeMidPoint;
                                        vector eVecNei = faceCentres[pFaces[j]]
                                            - edgeMidPoint;

                                        bool wrongOrient =  false;
                                        scalar dotProd = (fNorm & nNorm);
                                        //check for wrongly oriented surface
                                        if
                                        (
                                            dotProd < -0.9659
                                            && (eVecOwn & eVecNei) < 0
                                        )
                                        {
                                            wrongOrient = true;
                                        }


                                        const label& nRegion =
                                                      s[pFaces[j]].region();
                                        const scalar nFeatureAngle =
                                                refineSurfaces.featureRefineAngle
                                                (
                                                    surfI,
                                                    nRegion
                                                );

                                        const label nMaxLevel =
                                            refineSurfaces.maxLevel
                                            (
                                                surfI,
                                                nRegion
                                            );

                                        if (!wrongOrient)
                                        {
                                            const edge e(pointi, prevPointI);
                                            if
                                            (
                                                dotProd < featureAngle
                                                || dotProd < nFeatureAngle
                                             )
                                            {
                                                if (nMaxLevel > maxLevel)
                                                {
                                                    pointFeatureEdges.insert
                                                    (
                                                        e,
                                                        nRegion
                                                    );
                                                }
                                                else
                                                {
                                                    pointFeatureEdges.insert
                                                    (
                                                        e,
                                                        region
                                                    );
                                                }
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                        }

                        if (!foundNbr && refineBoundary && manifold)
                        {
                            const edge e(pointi, prevPointI);
                            pointFeatureEdges.insert
                            (
                                e,
                                region
                             );
                        }
                    }

                }
                forAllConstIter
                (
                    EdgeMap<label>,
                    pointFeatureEdges,
                    iter
                )
                {
                    featureEdges.insert(iter.key(), iter());
                    nFeatureEdges++;
                }
            }

            forAll(s.patches(), patchI)
            {
                Map<label> pointMap(2*nFeatureEdges);
                pointField pts(2*nFeatureEdges);
                edgeList edg(nFeatureEdges);
                label nPoints = 0;
                label nEdges = 0;
                forAllConstIter
                (
                    EdgeMap<label>,
                    featureEdges,
                    iter
                )
                {
                    label edgeRegion = iter();
                    if (edgeRegion == patchI)
                    {
                        label v0 = iter.key()[0];
                        label v1 = iter.key()[1];
                        edg[nEdges] = iter.key();
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
                }
                pts.setSize(nPoints);
                edg.setSize(nEdges);

                forAll(edg, edgeI)
                {
                    edg[edgeI][0] = pointMap[edg[edgeI][0]];
                    edg[edgeI][1] = pointMap[edg[edgeI][1]];
                }

                if (returnReduce(pts.size(), sumOp<label>()))
                {
                    fileName fName("allFeature" + Foam::name(sz) + ".ftr");

                    IOobject featObj
                    (
                        fName,                       // name
                        mesh.time().constant(),      // instance
                        "triSurface",                // local
                        mesh.time(),                 // registry
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    );


                    edgeMesh eMesh(pts,edg);

                    if (isA<distributedTriSurfaceMesh>(geom))
                    {
                        eMesh.reconstruct(mergeDist);
                    }

                    set
                    (
                        sz,
                        new extendedFeatureEdgeMesh
                        (
                            featObj,
                            extendedEdgeMesh
                            (
                                eMesh.points(),
                                eMesh.edges(),
                                0,
                                0,
                                eMesh.points().size(), // (point) nonFeatureStart
                                eMesh.edges().size(),  // (edge) internalStart
                                eMesh.edges().size(),  // (edge) flatStart
                                eMesh.edges().size(),  // (edge) openStart
                                eMesh.edges().size(),  // (edge) multipleStart
                                vectorField(0), // normals
                                List<extendedEdgeMesh::sideVolumeType>(0),
                                vectorField(0), // edgeDirections
                                labelListList(0),   // normalDirections
                                labelListList(0),   // edgeNormals
                                labelListList(0),   // featurePointNormals
                                labelListList(0),   // featurePointEdges
                                labelList(0)        // regionEdges
                             )
                        )
                    );

                    // Look up 'level' for single level
                    // distance based feature refinement for geometric
                    // entities not supported
                    levels_[sz] =
                        labelList(1, refineSurfaces.maxLevel(surfI,patchI));
                    distances_[sz] = scalarField(1, 0.0);
                    refinementOnly_[sz] =
                        refineSurfaces.featureRefineOnly(surfI,patchI);
                    sz++;
                }
            }
            s.clearOut();
        }
    }

    Info<< "Reading external feature lines." << endl;

    //Check for global geometry dict
    IOobject globalGeomHeader
    (
        "geometryDict",
        mesh.time().system(),
        mesh.time(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    dictionary globalGeomDict;
    if (globalGeomHeader.typeHeaderOk<IOdictionary>(true))
    {
        globalGeomDict = IOdictionary(globalGeomHeader);
    }
    else
    {
        globalGeomDict = dictionary();
    }

    forAll(featDicts, featI)
    {
        dictionary dict = featDicts[featI];
        fileName featFileName(dict.lookup("file"));

        // Try reading extendedEdgeMesh first

        IOobject extFeatObj
        (
            featFileName,                       // name
            mesh.time().constant(),               // instance
            "extendedFeatureEdgeMesh",          // local
            mesh.time(),                          // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        const fileName fName(typeFilePath<extendedFeatureEdgeMesh>(extFeatObj));

        if (!fName.empty() && extendedEdgeMesh::canRead(fName))
        {
            autoPtr<extendedEdgeMesh> eMeshPtr = extendedEdgeMesh::New
            (
                fName
            );

            Info<< "Read extendedFeatureEdgeMesh " << extFeatObj.name()
                << nl << incrIndent;
            eMeshPtr().writeStats(Info);
            Info<< decrIndent << endl;

            set(sz, new extendedFeatureEdgeMesh(extFeatObj, eMeshPtr()));
        }
        else
        {
            // Try reading edgeMesh

            IOobject featObj
            (
                featFileName,               // name
                mesh.time().constant(),     // instance
                "triSurface",               // local
                mesh.time(),                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            const fileName fName(typeFilePath<featureEdgeMesh>(featObj));

            if (fName.empty())
            {
                FatalIOErrorInFunction
                (
                    dict
                )   << "Could not open " << featObj.objectPath()
                    << exit(FatalIOError);
            }


            // Read as edgeMesh
            autoPtr<edgeMesh> eMeshPtr = edgeMesh::New(fName);
            const edgeMesh& eMesh = eMeshPtr();

            Info<< "Read edgeMesh " << featObj.name() << nl
                << incrIndent;
            eMesh.writeStats(Info);
            Info<< decrIndent << endl;


            // Analyse for feature points. These are all classified as mixed
            // points for lack of anything better
            const labelListList& pointEdges = eMesh.pointEdges();

            labelList oldToNew(eMesh.points().size(), -1);
            DynamicField<point> newPoints(eMesh.points().size());
            forAll(pointEdges, pointi)
            {
                if (pointEdges[pointi].size() > 2)
                {
                    oldToNew[pointi] = newPoints.size();
                    newPoints.append(eMesh.points()[pointi]);
                }
                //else if (pointEdges[pointi].size() == 2)
                //MEJ: do something based on a feature angle?
            }
            label nFeatures = newPoints.size();
            forAll(oldToNew, pointi)
            {
                if (oldToNew[pointi] == -1)
                {
                    oldToNew[pointi] = newPoints.size();
                    newPoints.append(eMesh.points()[pointi]);
                }
            }


            const edgeList& edges = eMesh.edges();
            edgeList newEdges(edges.size());
            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];
                newEdges[edgeI] = edge
                (
                    oldToNew[e[0]],
                    oldToNew[e[1]]
                );
            }

            // Construct an extendedEdgeMesh with
            // - all points on more than 2 edges : mixed feature points
            // - all edges : external edges

            extendedEdgeMesh eeMesh
            (
                newPoints,          // pts
                newEdges,           // eds
                0,                  // (point) concaveStart
                0,                  // (point) mixedStart
                nFeatures,          // (point) nonFeatureStart
                edges.size(),       // (edge) internalStart
                edges.size(),       // (edge) flatStart
                edges.size(),       // (edge) openStart
                edges.size(),       // (edge) multipleStart
                vectorField(0),     // normals
                List<extendedEdgeMesh::sideVolumeType>(0),// normalVolumeTypes
                vectorField(0),     // edgeDirections
                labelListList(0),   // normalDirections
                labelListList(0),   // edgeNormals
                labelListList(0),   // featurePointNormals
                labelListList(0),   // featurePointEdges
                identity(newEdges.size())   // regionEdges
            );

            //Info<< "Constructed extendedFeatureEdgeMesh " << featObj.name()
            //    << nl << incrIndent;
            //eeMesh.writeStats(Info);
            //Info<< decrIndent << endl;

            set(sz, new extendedFeatureEdgeMesh(featObj, eeMesh));
        }

        extendedEdgeMesh& eMesh = operator[](sz);

        //Add global transforms
        if (globalGeomDict.found("lines"))
        {
            const dictionary& ldict = globalGeomDict.subDict("lines");
            if (ldict.found(featFileName))
            {
                const dictionary& tdict = ldict.subDict(featFileName);
                if (tdict.found("transforms"))
                {
                    List<dictionary> globalTrans(tdict.lookup("transforms"));
                    if (dict.found("transforms"))
                    {
                        List<dictionary> localTrans = dict.lookup("transforms");
                        label globalSz = globalTrans.size();
                        label newSz = globalSz + localTrans.size();
                        globalTrans.setSize(newSz);
                        forAll(localTrans, transi)
                        {
                            globalTrans[globalSz++] = localTrans[transi];
                        }
                    }
                    dict.add("transforms",globalTrans,true);
                }
            }
        }

        //Perform optional transforms
        eMesh.doTransforms(dict);

        refinementOnly_[sz] =
            dict.lookupOrDefault("refineFeatureEdgesOnly", false);

        if (dict.found("levels"))
        {
            List<Tuple2<scalar, label>> distLevels(dict["levels"]);

            if (dict.size() < 1)
            {
                FatalErrorInFunction
                    << " : levels should be at least size 1" << endl
                    << "levels : "  << dict["levels"]
                    << exit(FatalError);
            }

            distances_[sz].setSize(distLevels.size());
            levels_[sz].setSize(distLevels.size());

            forAll(distLevels, j)
            {
                distances_[sz][j] = distLevels[j].first();
                levels_[sz][j] = distLevels[j].second();

                // Check in incremental order
                if (j > 0)
                {
                    if
                    (
                        (distances_[sz][j] <= distances_[sz][j-1])
                     || (levels_[sz][j] > levels_[sz][j-1])
                    )
                    {
                        FatalErrorInFunction
                            << " : Refinement should be specified in order"
                            << " of increasing distance"
                            << " (and decreasing refinement level)." << endl
                            << "Distance:" << distances_[sz][j]
                            << " refinementLevel:" << levels_[sz][j]
                            << exit(FatalError);
                    }
                }
            }
        }
        else
        {
            // Look up 'level' for single level
            levels_[sz] = labelList(1, readLabel(dict.lookup("level")));
            distances_[sz] = scalarField(1, 0.0);
        }

        Info<< "Refinement level according to distance to "
            << featFileName << " (" << eMesh.points().size() << " points, "
            << eMesh.edges().size() << " edges)." << endl;
        forAll(levels_[sz], j)
        {
            Info<< "    level " << levels_[sz][j]
                << " for all cells within " << distances_[sz][j]
                << " meter." << endl;
        }
        sz++;
    }

    this->resize(sz);
    levels_.resize(sz);
    distances_.resize(sz);
    refinementOnly_.resize(sz);
}


Foam::refinementFeatures::refinementFeatures
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts
)
:
    PtrList<extendedFeatureEdgeMesh>(featDicts.size()),
    checkRefinementOnly_(false),
    refinementOnly_(featDicts.size()),
    mFeatures_(*this,refinementOnly_,checkRefinementOnly_),
    distances_(featDicts.size()),
    levels_(featDicts.size()),
    edgeTrees_(0),
    pointTrees_(0),
    minCos_(-1)
{
    // Read features
    read(io, featDicts);
}


Foam::refinementFeatures::refinementFeatures
(
    const objectRegistry& io,
    const PtrList<dictionary>& featDicts,
    const scalar minCos
)
:
    PtrList<extendedFeatureEdgeMesh>(featDicts.size()),
    checkRefinementOnly_(false),
    refinementOnly_(featDicts.size()),
    mFeatures_(*this,refinementOnly_,checkRefinementOnly_),
    distances_(featDicts.size()),
    levels_(featDicts.size()),
    edgeTrees_(0),
    pointTrees_(0),
    minCos_(minCos)
{
    // Read features
    read(io, featDicts);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::refinementFeatures::checkSizes
(
    const scalar maxRatio,
    const boundBox& meshBb,
    const bool report,
    Ostream& os
) const
{
    if (report)
    {
        os<< "Checking for size." << endl;
    }

    bool hasError = false;

    forAll(*this, i)
    {
        const extendedFeatureEdgeMesh& em = operator[](i);
        if (returnReduce(em.points().size(), sumOp<label>()) == 0)
        {
            continue;
        }

        const boundBox bb(em.points(), true);

        for (label j = i+1; j < size(); j++)
        {
            const extendedFeatureEdgeMesh& em2 = operator[](j);

            if (returnReduce(em2.points().size(), sumOp<label>()) == 0)
            {
                continue;
            }

            const boundBox bb2(em2.points(), true);

            scalar ratio = bb.mag()/bb2.mag();

            if (ratio > maxRatio || ratio < 1.0/maxRatio)
            {
                hasError = true;

                if (report)
                {
                    os  << "    " << em.name()
                        << " bounds differ from " << em2.name()
                        << " by more than a factor 100:" << nl
                        << "        bounding box : " << bb << nl
                        << "        bounding box : " << bb2
                        << endl;
                }
            }
        }
    }

    forAll(*this, i)
    {
        const extendedFeatureEdgeMesh& em = operator[](i);
        const boundBox bb(em.points(), true);
        if (!meshBb.contains(bb))
        {
            if (report)
            {
                os  << "    " << em.name()
                    << " bounds not fully contained in mesh"<< nl
                    << "        bounding box      : " << bb << nl
                    << "        mesh bounding box : " << meshBb
                    << endl;
            }
        }
    }

    if (report)
    {
        os<< endl;
    }

    return returnReduce(hasError, orOp<bool>());
}


const Foam::PtrList<Foam::indexedOctree<Foam::treeDataEdge>>&
Foam::refinementFeatures::edgeTrees() const
{
    if (edgeTrees_.empty())
    {
        edgeTrees_.setSize(this->size());

        forAll(*this, featI)
        {
            const extendedEdgeMesh& eMesh = operator[](featI);
            const pointField& points = eMesh.points();
            const edgeList& edges = eMesh.edges();

            // Calculate bb of all points
            treeBoundBox bb(points);

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

            edgeTrees_.set
            (
                featI,
                new indexedOctree<treeDataEdge>
                (
                    treeDataEdge
                    (
                        false,                  // do not cache bb
                        edges,
                        points,
                        identity(edges.size())
                     ),
                    bb,     // overall search domain
                    8,      // maxLevel
                    10,     // leafsize
                    3.0     // duplicity
                 )
             );
        }
    }

    return edgeTrees_;
}


const Foam::PtrList<Foam::indexedOctree<Foam::treeDataPoint>>&
Foam::refinementFeatures::pointTrees() const
{
    if (pointTrees_.empty())
    {
        pointTrees_.setSize(this->size());

        forAll(*this, featI)
        {
            const extendedEdgeMesh& eMesh = operator[](featI);
            const pointField& points = eMesh.points();
            const edgeList& edges = eMesh.edges();
            const labelListList& pointEdges = eMesh.pointEdges();

            DynamicList<label> featurePoints;
            forAll(pointEdges, pointi)
            {
                const labelList& pEdges = pointEdges[pointi];
                if (pEdges.size() > 2)
                {
                    featurePoints.append(pointi);
                }
                else if (pEdges.size() == 2 && minCos_ > 0.0)
                {
                    // Check the angle
                    const edge& e0 = edges[pEdges[0]];
                    const edge& e1 = edges[pEdges[1]];

                    const point& p = points[pointi];
                    const point& p0 = points[e0.otherVertex(pointi)];
                    const point& p1 = points[e1.otherVertex(pointi)];

                    vector v0 = p-p0;
                    scalar v0Mag = mag(v0);

                    vector v1 = p1-p;
                    scalar v1Mag = mag(v1);

                    if
                    (
                        v0Mag > SMALL
                        && v1Mag > SMALL
                        && ((v0/v0Mag & v1/v1Mag) < minCos_)
                     )
                    {
                        featurePoints.append(pointi);
                    }
                }
            }

            // Calculate bb of all points
            treeBoundBox bb(points);

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

            pointTrees_.set
            (
                featI,
                new indexedOctree<treeDataPoint>
                (
                    treeDataPoint(points, featurePoints),
                    bb,     // overall search domain
                    8,      // maxLevel
                    10,     // leafsize
                    3.0     // duplicity
                 )
             );
        }
    }

    return pointTrees_;
}


void Foam::refinementFeatures::clearTrees()
{
    if (!pointTrees_.empty())
    {
        pointTrees_.clear();
    }
    if (!edgeTrees_.empty())
    {
        edgeTrees_.clear();
    }
    if (!regionEdgeTrees_.empty())
    {
        regionEdgeTrees_.clear();
    }
}


void Foam::refinementFeatures::findNearestEdge
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    List<pointIndexHit>& nearInfo,
    vectorField& nearNormal
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearInfo.setSize(samples.size());
    nearInfo = pointIndexHit();
    nearNormal.setSize(samples.size());
    nearNormal = Zero;

    forAll(edgeTrees(), featI)
    {
        if (checkRefinementOnly() && refinementOnly()[featI])
        {
            continue;
        }

        const indexedOctree<treeDataEdge>& tree = edgeTrees()[featI];

        if (tree.shapes().size() > 0)
        {
            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                scalar distSqr;
                if (nearInfo[sampleI].hit())
                {
                    distSqr = magSqr(nearInfo[sampleI].hitPoint()-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[sampleI];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearFeature[sampleI] = featI;
                    nearInfo[sampleI] = pointIndexHit
                    (
                        info.hit(),
                        info.hitPoint(),
                        tree.shapes().edgeLabels()[info.index()]
                    );

                    const treeDataEdge& td = tree.shapes();
                    const edge& e = td.edges()[nearInfo[sampleI].index()];
                    nearNormal[sampleI] =  e.vec(td.points());
                    nearNormal[sampleI] /= mag(nearNormal[sampleI])+VSMALL;
                }
            }
        }
    }
}


void Foam::refinementFeatures::findNearestRegionEdge
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    List<pointIndexHit>& nearInfo,
    vectorField& nearNormal
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearInfo.setSize(samples.size());
    nearInfo = pointIndexHit();
    nearNormal.setSize(samples.size());
    nearNormal = Zero;


    const PtrList<indexedOctree<treeDataEdge>>& regionTrees =
        regionEdgeTrees();

    forAll(regionTrees, featI)
    {
        const indexedOctree<treeDataEdge>& regionTree = regionTrees[featI];

        forAll(samples, sampleI)
        {
            const point& sample = samples[sampleI];

            scalar distSqr;
            if (nearInfo[sampleI].hit())
            {
                distSqr = magSqr(nearInfo[sampleI].hitPoint()-sample);
            }
            else
            {
                distSqr = nearestDistSqr[sampleI];
            }

            // Find anything closer than current best
            pointIndexHit info = regionTree.findNearest(sample, distSqr);

            if (info.hit())
            {
                const treeDataEdge& td = regionTree.shapes();

                nearFeature[sampleI] = featI;
                nearInfo[sampleI] = pointIndexHit
                (
                    info.hit(),
                    info.hitPoint(),
                    regionTree.shapes().edgeLabels()[info.index()]
                );

                const edge& e = td.edges()[nearInfo[sampleI].index()];
                nearNormal[sampleI] =  e.vec(td.points());
                nearNormal[sampleI] /= mag(nearNormal[sampleI])+VSMALL;
            }
        }
    }
}


void Foam::refinementFeatures::findNearestPoint
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& nearFeature,
    List<pointIndexHit>& nearInfo
) const
{
    nearFeature.setSize(samples.size());
    nearFeature = -1;
    nearInfo.setSize(samples.size());
    nearInfo = pointIndexHit();

    forAll(pointTrees(), featI)
    {
        if (checkRefinementOnly() && refinementOnly()[featI])
        {
            continue;
        }

        const indexedOctree<treeDataPoint>& tree = pointTrees()[featI];

        if (tree.shapes().pointLabels().size() > 0)
        {
            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                scalar distSqr;
                if (nearFeature[sampleI] != -1)
                {
                    distSqr = magSqr(nearInfo[sampleI].hitPoint()-sample);
                }
                else
                {
                    distSqr = nearestDistSqr[sampleI];
                }

                pointIndexHit info = tree.findNearest(sample, distSqr);

                if (info.hit())
                {
                    nearFeature[sampleI] = featI;
                    nearInfo[sampleI] = pointIndexHit
                    (
                        info.hit(),
                        info.hitPoint(),
                        tree.shapes().pointLabels()[info.index()]
                    );
                }
            }
        }
    }
}


void Foam::refinementFeatures::findHigherLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& maxLevel
)
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;

    forAll(*this, featI)
    {
        findHigherLevel(pt, featI, maxLevel);
    }
}


Foam::scalar Foam::refinementFeatures::maxDistance() const
{
    scalar overallMax = -GREAT;
    forAll(distances_, featI)
    {
        overallMax = max(overallMax, max(distances_[featI]));
    }
    return overallMax;
}


void Foam::refinementFeatures::trim(const polyMesh& mesh)
{
    Random rndGen(653213);

    // Determine mesh bounding boxes:
    List<treeBoundBox> meshBb(1);
    {
        meshBb[0] =
        treeBoundBox
        (
            boundBox(mesh.points(), false)
        ).extend(rndGen, 1E-3);
    }

    forAll(*this, i)
    {
        operator[](i).edgeMesh::trim(meshBb);
    }

    //clear feature line search trees
    clearTrees();
}

// ************************************************************************* //
