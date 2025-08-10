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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include "triSurface/triSurface.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "algorithms/indexedOctree/volumeType.H"

#include "db/IOstreams/IOstreams/IOmanip.H"

#ifdef FOAM_USE_TBB
  //#include <tbb/tick_count.h>
  #include <tbb/parallel_for.h>
  #include "include/TBBTimer.H"
#endif

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::triSurfaceSearch::checkUniqueHit
(
    const pointIndexHit& currHit,
    const DynamicList<pointIndexHit, 1, 1>& hits,
    const vector& lineVec
) const
{
    // Classify the hit
    label nearType = -1;
    label nearLabel = -1;

    const labelledTri& f = surface()[currHit.index()];

    f.nearestPointClassify
    (
        currHit.hitPoint(),
        surface().points(),
        nearType,
        nearLabel
    );

    if (nearType == triPointRef::POINT)
    {
        // near point

        const label nearPointi = f[nearLabel];

        const labelList& pointFaces =
            surface().pointFaces()[surface().meshPointMap()[nearPointi]];

        forAll(pointFaces, pI)
        {
            const label pointFacei = pointFaces[pI];

            if (pointFacei != currHit.index())
            {
                forAll(hits, hI)
                {
                    const pointIndexHit& hit = hits[hI];

                    if (hit.index() == pointFacei)
                    {
                        return false;
                    }
                }
            }
        }
    }
    else if (nearType == triPointRef::EDGE)
    {
        // near edge
        // check if the other face of the edge is already hit

        const labelList& fEdges = surface().faceEdges()[currHit.index()];

        const label edgeI = fEdges[nearLabel];

        const labelList& edgeFaces = surface().edgeFaces()[edgeI];

        forAll(edgeFaces, fI)
        {
            const label edgeFacei = edgeFaces[fI];

            if (edgeFacei != currHit.index())
            {
                forAll(hits, hI)
                {
                    const pointIndexHit& hit = hits[hI];

                    if (hit.index() == edgeFacei)
                    {
                        // Check normals
                        const vector currHitNormal =
                            surface().faceNormals()[currHit.index()];

                        const vector existingHitNormal =
                            surface().faceNormals()[edgeFacei];

                        const label signCurrHit =
                            pos0(currHitNormal & lineVec);

                        const label signExistingHit =
                            pos0(existingHitNormal & lineVec);

                        if (signCurrHit == signExistingHit)
                        {
                            return false;
                        }
                    }
                }
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceSearch::triSurfaceSearch(const triSurface& surface)
:
    surface_(surface),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    maxTreeDepth_(10),
    treePtr_(nullptr)
{}


Foam::triSurfaceSearch::triSurfaceSearch
(
    const triSurface& surface,
    const dictionary& dict
)
:
    surface_(surface),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
    maxTreeDepth_(10),
    treePtr_(nullptr)
{
    //Info<< "triSurfaceSearch: Octree search constructor" << endl;

    // Have optional non-standard search tolerance for gappy surfaces.
    if (dict.readIfPresent("tolerance", tolerance_) && tolerance_ > 0)
    {
        Info<< "    using intersection tolerance " << tolerance_ << endl;
    }

    // Have optional non-standard tree-depth to limit storage.
    if (dict.readIfPresent("maxTreeDepth", maxTreeDepth_) && maxTreeDepth_ > 0)
    {
        Info<< "    using maximum tree depth " << maxTreeDepth_ << endl;
    }
}


Foam::triSurfaceSearch::triSurfaceSearch
(
    const triSurface& surface,
    const scalar tolerance,
    const label maxTreeDepth
)
:
    surface_(surface),
    tolerance_(tolerance),
    maxTreeDepth_(maxTreeDepth),
    treePtr_(nullptr)
{
    if (tolerance_ < 0)
    {
        tolerance_ = indexedOctree<treeDataTriSurface>::perturbTol();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceSearch::~triSurfaceSearch()
{
    clearOut();
}


void Foam::triSurfaceSearch::clearOut()
{
    treePtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataTriSurface>&
Foam::triSurfaceSearch::tree() const
{
    if (treePtr_.empty())
    {
        //Timer t("triSurfaceSearch::tree() Octree");
        // Calculate bb without constructing local point numbering.
        treeBoundBox bb(Zero, Zero);

        if (surface().size())
        {
        //Timer t("    triSurfaceSearch::tree() calcBounds");
            label nPoints;
            PatchTools::calcBounds(surface(), bb, nPoints);

            if (nPoints != surface().points().size())
            {
                WarningInFunction
                    << "Surface does not have compact point numbering."
                    << " Of " << surface().points().size()
                    << " only " << nPoints
                    << " are used. This might give problems in some routines."
                    << endl;
            }

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        }

        const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
        indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

        treePtr_.reset
        (
            new indexedOctree<treeDataTriSurface>
            (
                treeDataTriSurface(false, surface_, tolerance_),
                bb,
                maxTreeDepth_,  // maxLevel
                10,             // leafsize
                3.0             // duplicity
            )
        );

        indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
    }

    return treePtr_();
}


// Determine inside/outside for samples
Foam::boolList Foam::triSurfaceSearch::calcInside
(
    const pointField& samples
) const
{
    boolList inside(samples.size());

    forAll(samples, sampleI)
    {
        const point& sample = samples[sampleI];

        if (!tree().bb().contains(sample))
        {
            inside[sampleI] = false;
        }
        else if (tree().getVolumeType(sample) == volumeType::INSIDE)
        {
            inside[sampleI] = true;
        }
        else
        {
            inside[sampleI] = false;
        }
    }
    return inside;
}


void Foam::triSurfaceSearch::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info,
    const bool threaded, /*= false*/
    const scalar isoValue /*= 0*/
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

//#ifdef FOAM_USE_TBB
//    Timer t("        triSurfaceSearch::findNearest() Octree");
//#endif

    if (samples.size() && tree().nodes().size())
    {
        info.setSize(samples.size());

        const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
        indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

        const treeDataTriSurface::findNearestOp fOp(octree);

        if (threaded && samples.size() > 2048)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    samples.size() / tbb::this_task_arena::max_concurrency(),
                    1024
                );

            tbb::parallel_for
            (
                tbb::blocked_range<size_t>(0, samples.size(), grainSize),
                [&](const tbb::blocked_range<size_t>& r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        scalar maxDist = 1.732*sqrt(mag(nearestDistSqr[i]));
                        point minSamplePt = samples[i];
                        point maxSamplePt = samples[i];
                        minSamplePt -= maxDist*vector::one;
                        maxSamplePt += maxDist*vector::one;
                        boundBox sampleBox(minSamplePt,maxSamplePt);
                        if (octree.bb().overlaps(sampleBox))
                        {
                            info[i] = octree.findNearest
                            (
                                samples[i],
                                nearestDistSqr[i],
                                fOp
                            );
                        }
                        else
                        {
                            info[i].setMiss();
                        }
                    }
                },
                tbb::simple_partitioner()
            );
#endif
        }
        else
        {
            forAll(samples, i)
            {
                scalar maxDist = 1.732*sqrt(mag(nearestDistSqr[i]));
                point minSamplePt = samples[i];
                point maxSamplePt = samples[i];
                minSamplePt -= maxDist*vector::one;
                maxSamplePt += maxDist*vector::one;
                boundBox sampleBox(minSamplePt,maxSamplePt);
                if (octree.bb().overlaps(sampleBox))
                {
                    info[i] = octree.findNearest
                    (
                        samples[i],
                        nearestDistSqr[i],
                        fOp
                    );
                }
                else
                {
                    info[i].setMiss();
                }
            }
        }

        indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
    }
    else
    {
        if (info.size())
        {
            //set to miss
            info = pointIndexHit(false,vector::zero,-1);
        }
    }
}


Foam::pointIndexHit Foam::triSurfaceSearch::nearest
(
    const point& pt,
    const vector& span
)
const
{
    //Timer t("triSurfaceSearch::nearest");
    const scalar nearestDistSqr = 0.25*magSqr(span);

    return tree().findNearest(pt, nearestDistSqr);
}


void Foam::triSurfaceSearch::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    //Timer t("triSurfaceSearch::findLine Octree");
    const indexedOctree<treeDataTriSurface>& octree = tree();

    const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    info.setSize(start.size());

    forAll(start, i)
    {
        info[i] = octree.findLine(start[i], end[i]);
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
} //findLine


void Foam::triSurfaceSearch::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info,
    const bool threaded /*= false*/
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

//#ifdef FOAM_USE_TBB
//    Timer t("triSurfaceSearch::findLineAny Octree");
//#endif

    info.setSize(start.size());

    const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    if (threaded && start.size() > 2048)
    {
#ifndef FOAM_USE_TBB
        WarningInFunction
            << "FOAM not linked against TBB! Cannot run multithreaded."
            << endl;
#else
        const auto grainSize =
            std::max<size_t>
            (
                start.size() / tbb::this_task_arena::max_concurrency(),
                1024
            );

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, start.size(), grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    info[i] = octree.findLineAny(start[i], end[i]);
                }
            },
            tbb::simple_partitioner()
        );
#endif
    }
    else
    {
        forAll(start, i)
        {
            info[i] = octree.findLineAny(start[i], end[i]);
        }
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
} // findLineAny


void Foam::triSurfaceSearch::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();

    //Timer t("triSurfaceSearch::findLineAll Octree");

    info.setSize(start.size());

    const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
    indexedOctree<treeDataTriSurface>::perturbTol() = tolerance();

    // Work array
    DynamicList<pointIndexHit, 1, 1> hits;

    DynamicList<label> shapeMask;

    treeDataTriSurface::findAllIntersectOp allIntersectOp(octree, shapeMask);

    forAll(start, pointi)
    {
        hits.clear();
        shapeMask.clear();

        while (true)
        {
            // See if any intersection between pt and end
            pointIndexHit inter = octree.findLine
            (
                start[pointi],
                end[pointi],
                allIntersectOp
            );

            if (inter.hit())
            {
                vector lineVec = end[pointi] - start[pointi];
                lineVec.normalise();

                if
                (
                    checkUniqueHit
                    (
                        inter,
                        hits,
                        lineVec
                    )
                )
                {
                    hits.append(inter);
                }

                shapeMask.append(inter.index());
            }
            else
            {
                break;
            }
        }

        info[pointi].transfer(hits);
    }

    indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceSearch::findSphere
(
    const pointField& centres,
    const scalarField& radiusSqr,
    List<labelList>& hitIndexes,
    labelList& nOverlapChecks
) const
{
    const indexedOctree<treeDataTriSurface>& octree = tree();
    forAll(centres, centrei)
    {
        hitIndexes[centrei] = octree.findSphere
        (
            centres[centrei],
            radiusSqr[centrei],
            nOverlapChecks[centrei]
        );
    }
}

// ************************************************************************* //
