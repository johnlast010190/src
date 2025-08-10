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
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "searchableSurfaces/searchableBox/searchableBox.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/Lists/SortableList/SortableList.H"

#ifdef FOAM_USE_TBB
  #include <tbb/parallel_for.h>
#endif
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(searchableBox, 0);
    addToRunTimeSelectionTable(searchableSurface, searchableBox, dict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::searchableBox::projectOntoCoordPlane
(
    const direction dir,
    const point& planePt,
    pointIndexHit& info
) const
{
    // Set point
    info.rawPoint()[dir] = planePt[dir];
    // Set face
    if (planePt[dir] == min()[dir])
    {
        info.setIndex(dir*2);
    }
    else if (planePt[dir] == max()[dir])
    {
        info.setIndex(dir*2+1);
    }
    else
    {
        FatalErrorInFunction
            << "Point on plane " << planePt
            << " is not on coordinate " << min()[dir]
            << " nor " << max()[dir] << abort(FatalError);
    }
}


// Returns miss or hit with face (0..5) and region(always 0)
Foam::pointIndexHit Foam::searchableBox::findNearest
(
    const point& bbMid,
    const point& sample,
    const scalar nearestDistSqr
) const
{
    // Point can be inside or outside. For every component direction can be
    // left of min, right of max or inbetween.
    // - outside points: project first one x plane (either min().x()
    // or max().x()), then onto y plane and finally z. You should be left
    // with intersection point
    // - inside point: find nearest side (compare to mid point). Project onto
    //   that.

    // The face is set to the last projected face.


    // Outside point projected onto cube. Assume faces 0..5.
    pointIndexHit info(true, sample, -1);
    bool outside = false;

    // (for internal points) per direction what nearest cube side is
    point near;

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (info.rawPoint()[dir] < min()[dir])
        {
            projectOntoCoordPlane(dir, min(), info);
            outside = true;
        }
        else if (info.rawPoint()[dir] > max()[dir])
        {
            projectOntoCoordPlane(dir, max(), info);
            outside = true;
        }
        else if (info.rawPoint()[dir] > bbMid[dir])
        {
            near[dir] = max()[dir];
        }
        else
        {
            near[dir] = min()[dir];
        }
    }


    // For outside points the info will be correct now. Handle inside points
    // using the three near distances. Project onto the nearest plane.
    if (!outside)
    {
        vector dist(cmptMag(info.rawPoint() - near));

        if (dist.x() < dist.y())
        {
            if (dist.x() < dist.z())
            {
                // Project onto x plane
                projectOntoCoordPlane(vector::X, near, info);
            }
            else
            {
                projectOntoCoordPlane(vector::Z, near, info);
            }
        }
        else
        {
            if (dist.y() < dist.z())
            {
                projectOntoCoordPlane(vector::Y, near, info);
            }
            else
            {
                projectOntoCoordPlane(vector::Z, near, info);
            }
        }
    }


    // Check if outside. Optimisation: could do some checks on distance already
    // on components above
    if (magSqr(info.rawPoint() - sample) > nearestDistSqr)
    {
        info.setMiss();
        info.setIndex(-1);
    }

    return info;
}


Foam::treeBoundBox Foam::searchableBox::calcTransforms
(
    const dictionary& dict
)
{
    point minLoc = dict.lookup("min");
    point maxLoc = dict.lookup("max");

    if (dict.found("transforms"))
    {
        PtrList<dictionary> transforms(dict.lookup("transforms"));

        forAll(transforms, dictI)
        {
            const dictionary& transformDict = transforms[dictI];

            const word type(transformDict.lookup("type"));
            Info<< searchableSurface::name() << " : applying transform " << type
                << endl;

            if (type == "scale")
            {
                const vector scale(transformDict.lookup("scaleVec"));
                const point pt
                (
                    transformDict.lookupOrDefault("aboutPoint", vector(0, 0, 0))
                );
                minLoc -= pt;
                maxLoc -= pt;
                minLoc = cmptMultiply(scale,minLoc);
                maxLoc = cmptMultiply(scale,maxLoc);
                minLoc += pt;
                maxLoc += pt;
            }
            else if (type == "translate")
            {
                const vector v(transformDict.lookup("translateVec"));
                minLoc += v;
                maxLoc += v;
            }
            else
            {
                WarningInFunction
                    << "Surface transform: " << type
                    << " not available for searchableSurface: "
                    << searchableSurface::name()<<endl;
            }
        }
    }

    return treeBoundBox(minLoc, maxLoc);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableBox::searchableBox
(
    const IOobject& io,
    const treeBoundBox& bb
)
:
    searchableSurface(io),
    treeBoundBox(bb)
{
    if (!contains(midpoint()))
    {
        FatalErrorInFunction
            << "Illegal bounding box specification : "
            << static_cast<const treeBoundBox>(*this) << exit(FatalError);
    }

    bounds() = static_cast<boundBox>(*this);
}


Foam::searchableBox::searchableBox
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    treeBoundBox(calcTransforms(dict))
{
    if (!contains(midpoint()))
    {
        FatalErrorInFunction
            << "Illegal bounding box specification : "
            << static_cast<const treeBoundBox>(*this) << exit(FatalError);
    }

    bounds() = static_cast<boundBox>(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableBox::~searchableBox()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableBox::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


Foam::tmp<Foam::pointField> Foam::searchableBox::coordinates() const
{
    tmp<pointField> tCtrs = tmp<pointField>(new pointField(6));
    pointField& ctrs = tCtrs.ref();

    const pointField pts(treeBoundBox::points());
    const faceList& fcs = treeBoundBox::faces;

    forAll(fcs, i)
    {
        ctrs[i] = fcs[i].centre(pts);
    }

    return tCtrs;
}


void Foam::searchableBox::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(size());
    radiusSqr.setSize(size());
    radiusSqr = 0.0;

    const pointField pts(treeBoundBox::points());
    const faceList& fcs = treeBoundBox::faces;

    forAll(fcs, i)
    {
        const face& f = fcs[i];

        centres[i] = f.centre(pts);
        forAll(f, fp)
        {
            const point& pt = pts[f[fp]];

            radiusSqr[i] = Foam::max
            (
                radiusSqr[i],
                Foam::magSqr(pt-centres[i])
            );
        }
    }

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}


Foam::tmp<Foam::pointField> Foam::searchableBox::points() const
{
    return treeBoundBox::points();
}


Foam::pointIndexHit Foam::searchableBox::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    return findNearest(midpoint(), sample, nearestDistSqr);
}


Foam::pointIndexHit Foam::searchableBox::findNearestOnEdge
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    const point bbMid(midpoint());

    // Outside point projected onto cube. Assume faces 0..5.
    pointIndexHit info(true, sample, -1);
    bool outside = false;

    // (for internal points) per direction what nearest cube side is
    point near;

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (info.rawPoint()[dir] < min()[dir])
        {
            projectOntoCoordPlane(dir, min(), info);
            outside = true;
        }
        else if (info.rawPoint()[dir] > max()[dir])
        {
            projectOntoCoordPlane(dir, max(), info);
            outside = true;
        }
        else if (info.rawPoint()[dir] > bbMid[dir])
        {
            near[dir] = max()[dir];
        }
        else
        {
            near[dir] = min()[dir];
        }
    }


    // For outside points the info will be correct now. Handle inside points
    // using the three near distances. Project onto the nearest two planes.
    if (!outside)
    {
        // Get the per-component distance to nearest wall
        vector dist(cmptMag(info.rawPoint() - near));

        SortableList<scalar> sortedDist(3);
        sortedDist[0] = dist[0];
        sortedDist[1] = dist[1];
        sortedDist[2] = dist[2];
        sortedDist.sort();

        // Project onto nearest
        projectOntoCoordPlane(sortedDist.indices()[0], near, info);
        // Project onto second nearest
        projectOntoCoordPlane(sortedDist.indices()[1], near, info);
    }


    // Check if outside. Optimisation: could do some checks on distance already
    // on components above
    if (magSqr(info.rawPoint() - sample) > nearestDistSqr)
    {
        info.setMiss();
        info.setIndex(-1);
    }

    return info;
}


Foam::pointIndexHit Foam::searchableBox::findNearest
(
    const linePointRef& ln,
    treeBoundBox& tightest,
    point& linePoint
) const
{
    NotImplemented;
}


Foam::pointIndexHit Foam::searchableBox::findLine
(
    const point& start,
    const point& end
) const
{
    pointIndexHit info(false, start, -1);

    bool foundInter;

    if (posBits(start) == 0)
    {
        if (posBits(end) == 0)
        {
            // Both start and end inside.
            foundInter = false;
        }
        else
        {
            // end is outside. Clip to bounding box.
            foundInter = intersects(end, start, info.rawPoint());
        }
    }
    else
    {
        // start is outside. Clip to bounding box.
        foundInter = intersects(start, end, info.rawPoint());
    }


    // Classify point
    if (foundInter)
    {
        info.setHit();

        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            if (info.rawPoint()[dir] == min()[dir])
            {
                info.setIndex(2*dir);
                break;
            }
            else if (info.rawPoint()[dir] == max()[dir])
            {
                info.setIndex(2*dir+1);
                break;
            }
        }

        if (info.index() == -1)
        {
            FatalErrorInFunction
                << "point " << info.rawPoint()
                << " on segment " << start << end
                << " should be on face of " << *this
                << " but it isn't." << abort(FatalError);
        }
    }

    return info;
}


Foam::pointIndexHit Foam::searchableBox::findLineAny
(
    const point& start,
    const point& end
) const
{
    return findLine(start, end);
}


void Foam::searchableBox::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info,
    const bool threaded /*= false*/
) const
{
    info.setSize(samples.size());

    const point bbMid(midpoint());

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
                    info[i] = findNearest(bbMid, samples[i], nearestDistSqr[i]);
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
            info[i] = findNearest(bbMid, samples[i], nearestDistSqr[i]);
        }
    }
}


void Foam::searchableBox::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        info[i] = findLine(start[i], end[i]);
    }
}


void Foam::searchableBox::findLineAny
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info,
    const bool threaded /*= false*/
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        info[i] = findLineAny(start[i], end[i]);
    }
}


void Foam::searchableBox::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    // Work array
    DynamicList<pointIndexHit, 1, 1> hits;

    // Tolerances:
    // To find all intersections we add a small vector to the last intersection
    // This is chosen such that
    // - it is significant (SMALL is smallest representative relative tolerance;
    //   we need something bigger since we're doing calculations)
    // - if the start-end vector is zero we still progress
    const vectorField dirVec(end-start);
    const scalarField magSqrDirVec(magSqr(dirVec));
    const vectorField smallVec
    (
        ROOTSMALL*dirVec
      + vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL)
    );

    forAll(start, pointi)
    {
        // See if any intersection between pt and end
        pointIndexHit inter = findLine(start[pointi], end[pointi]);

        if (inter.hit())
        {
            hits.clear();
            hits.append(inter);

            point pt = inter.hitPoint() + smallVec[pointi];

            while (((pt-start[pointi])&dirVec[pointi]) <= magSqrDirVec[pointi])
            {
                // See if any intersection between pt and end
                pointIndexHit inter = findLine(pt, end[pointi]);

                // Check for not hit or hit same face as before (can happen
                // if vector along surface of face)
                if
                (
                    !inter.hit()
                 || (inter.index() == hits.last().index())
                )
                {
                    break;
                }
                hits.append(inter);

                pt = inter.hitPoint() + smallVec[pointi];
            }

            info[pointi].transfer(hits);
        }
        else
        {
            info[pointi].clear();
        }
    }
}


void Foam::searchableBox::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region,
    const bool threaded /*= false*/
) const
{
    region.setSize(info.size());
    region = 0;
}


void Foam::searchableBox::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = Zero;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            normal[i] = treeBoundBox::faceNormals[info[i].index()];
        }
        else
        {
            // Set to what?
        }
    }
}


void Foam::searchableBox::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType,
    const bool threaded /*= false*/
) const
{
    volType.setSize(points.size());
    volType = volumeType::INSIDE;

    if (threaded && points.size() > 2048)
    {
#ifndef FOAM_USE_TBB
        WarningInFunction
            << "FOAM not linked against TBB! Cannot run multithreaded."
            << endl;
#else
        const auto grainSize =
            std::max<size_t>
            (
                points.size() / tbb::this_task_arena::max_concurrency(),
                1024
            );

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, points.size(), grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t pointi = r.begin(); pointi < r.end(); ++pointi)
                {
                    const point& pt = points[pointi];

                    for (direction dir = 0; dir < vector::nComponents; dir++)
                    {
                        if (pt[dir] < min()[dir] || pt[dir] > max()[dir])
                        {
                            volType[pointi] = volumeType::OUTSIDE;
                            break;
                        }
                    }
                }
            },
            tbb::simple_partitioner()
        );
#endif
    }
    else
    {
        forAll(points, pointi)
        {
            const point& pt = points[pointi];

            for (direction dir = 0; dir < vector::nComponents; dir++)
            {
                if (pt[dir] < min()[dir] || pt[dir] > max()[dir])
                {
                    volType[pointi] = volumeType::OUTSIDE;
                    break;
                }
            }
        }
    }
}


// ************************************************************************* //
