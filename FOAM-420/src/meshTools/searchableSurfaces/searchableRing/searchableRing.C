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
    (c) 2011 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "searchableSurfaces/searchableRing/searchableRing.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableRing, 0);
addToRunTimeSelectionTable(searchableSurface, searchableRing, dict);

}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


Foam::tmp<Foam::pointField> Foam::searchableRing::coordinates() const
{
    tmp<pointField> tCtrs(new pointField(1, 0.5*(point1_ + point2_)));

    return tCtrs;
}

void Foam::searchableRing::boundingSpheres
(
    pointField& centres,
    scalarField& radiusSqr
) const
{
    centres.setSize(1);
    centres[0] = 0.5*(point1_ + point2_);

    radiusSqr.setSize(1);
    radiusSqr[0] = Foam::magSqr(point1_-centres[0]) + Foam::sqr(outerRadius_);

    // Add a bit to make sure all points are tested inside
    radiusSqr += Foam::sqr(SMALL);
}

Foam::tmp<Foam::pointField> Foam::searchableRing::points() const
{
    tmp<pointField> tPts(new pointField(2));
    pointField& pts = tPts.ref();

    pts[0] = point1_;
    pts[1] = point2_;

    return tPts;
}

Foam::pointIndexHit Foam::searchableRing::findNearest
(
    const point& sample,
    const scalar nearestDistSqr
) const
{
    pointIndexHit info(false, sample, -1);

    pointIndexHit infoInner =
        innerCylinder_.findNearest(sample, nearestDistSqr);

    pointIndexHit infoOuter =
        outerCylinder_.findNearest(sample, nearestDistSqr);

    if (infoInner.hit() && infoOuter.hit())
    {
        vector v(infoOuter.hitPoint() - point1_);
        // Decompose sample-point1 into normal and parallel component
        scalar parallel = (v & unitDir_);

        // Remove the parallel component
        v -= parallel*unitDir_;
        scalar magV = mag(v);

        if (magV >= (innerRadius_ - SMALL))
        {
            if
            (
                magSqr(infoInner.hitPoint() - sample)
                < magSqr(infoOuter.hitPoint() - sample)
            )
            {
                info.setPoint(infoInner.hitPoint());
                info.setHit();
                info.setIndex(0);
            }
            else
            {
                info.setPoint(infoOuter.hitPoint());
                info.setHit();
                info.setIndex(0);
            }
        }
        else
        {
            info.setPoint(infoInner.hitPoint());
            info.setHit();
            info.setIndex(0);
        }
    }
    else if (infoInner.hit())
    {
        info.setPoint(infoInner.hitPoint());
        info.setHit();
        info.setIndex(0);
    }
    else if (infoOuter.hit())
    {
        vector v(infoOuter.hitPoint() - point1_);
        // Decompose sample-point1 into normal and parallel component
        scalar parallel = (v & unitDir_);

        // Remove the parallel component
        v -= parallel*unitDir_;
        scalar magV = mag(v);

        if (magV >= (innerRadius_ - SMALL))
        {
            info.setPoint(infoOuter.hitPoint());
            info.setHit();
            info.setIndex(0);
        }
    }

    return info;
}


Foam::scalar Foam::searchableRing::radius2(const point& pt) const
{
    const vector x = (pt-point1_) ^ unitDir_;
    return x&x;
}


// From http://www.gamedev.net/community/forums/topic.asp?topic_id=467789 -
// intersection of ring with ray
void Foam::searchableRing::findLineAll
(
    const point& start,
    const point& end,
    pointIndexHit& near,
    pointIndexHit& far
) const
{
    near.setMiss();
    far.setMiss();

    pointIndexHit nearInner, farInner;
    nearInner.setMiss();
    farInner.setMiss();
    innerCylinder_.findLineAll(start, end, nearInner, farInner);

    pointIndexHit nearOuter, farOuter;
    nearOuter.setMiss();
    farOuter.setMiss();
    outerCylinder_.findLineAll(start, end, nearOuter, farOuter);

    if (nearOuter.hit())
    {
        vector v(nearOuter.hitPoint() - point1_);
        // Decompose sample-point1 into normal and parallel component
        scalar parallel = (v & unitDir_);

        // Remove the parallel component
        v -= parallel*unitDir_;
        scalar magV = mag(v);

        if (magV >= (innerRadius_ - SMALL))
        {
            near.setPoint(nearOuter.hitPoint());
            near.setHit();
            near.setIndex(0);
        }
        else
        {
            if (nearInner.hit())
            {
                near.setPoint(nearInner.hitPoint());
                near.setHit();
                near.setIndex(0);
            }
        }
    }
    else  if (nearInner.hit())
    {
        near.setPoint(nearInner.hitPoint());
        near.setHit();
        near.setIndex(0);
    }

    if (farOuter.hit())
    {
        vector v(farOuter.hitPoint() - point1_);
        // Decompose sample-point1 into normal and parallel component
        scalar parallel = (v & unitDir_);

        // Remove the parallel component
        v -= parallel*unitDir_;
        scalar magV = mag(v);

        if (magV >= (innerRadius_ - SMALL))
        {
            far.setPoint(farOuter.hitPoint());
            far.setHit();
            far.setIndex(0);
        }
        else
        {
            if (farInner.hit())
            {
                far.setPoint(farInner.hitPoint());
                far.setHit();
                far.setIndex(0);
            }
        }
    }
    else if (farInner.hit())
    {
        far.setPoint(farInner.hitPoint());
        far.setHit();
        far.setIndex(0);
    }
}


Foam::boundBox Foam::searchableRing::calcBounds() const
{

    // Adapted from
    // http://www.gamedev.net/community/forums
    //       /topic.asp?topic_id=338522&forum_id=20&gforum_id=0

    // Let ring have end points A,B and radius r,

    // Bounds in direction X (same for Y and Z) can be found as:
    // Let A.X<B.X (otherwise swap points)
    // Good approximate lowest bound is A.X-r and highest is B.X+r (precise for
    // capsule). At worst, in one direction it can be larger than needed by 2*r.

    // Accurate bounds for ring is
    // A.X-kx*r, B.X+kx*r
    // where
    // kx=sqrt(((A.Y-B.Y)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))

    // similar thing for Y and Z
    // (i.e.
    // ky=sqrt(((A.X-B.X)^2+(A.Z-B.Z)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
    // kz=sqrt(((A.X-B.X)^2+(A.Y-B.Y)^2)/((A.X-B.X)^2+(A.Y-B.Y)^2+(A.Z-B.Z)^2))
    // )

    // How derived: geometric reasoning. Bounds of ring is same as for 2
    // circles centered on A and B. This sqrt thingy gives sine of angle between
    // axis and direction, used to find projection of radius.

    vector kr
    (
        sqrt(sqr(unitDir_.y()) + sqr(unitDir_.z())),
        sqrt(sqr(unitDir_.x()) + sqr(unitDir_.z())),
        sqrt(sqr(unitDir_.x()) + sqr(unitDir_.y()))
    );

    kr *= outerRadius_;

    point min = point1_ - kr;
    point max = point1_ + kr;

    min = ::Foam::min(min, point2_ - kr);
    max = ::Foam::max(max, point2_ + kr);

    return boundBox(min, max);
}


Foam::point Foam::searchableRing::endPoint
(
    const dictionary& dict
) const
{
    if (dict.found("point2"))
    {
        return dict.lookup("point2");
    }
    else
    {
        point pt1 = dict.lookup("point1");
        vector dir = dict.lookup("direction");
        scalar length = readScalar(dict.lookup("length"));
        vector unitDir = dir/mag(dir);

        return pt1 + (length*unitDir);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableRing::searchableRing
(
    const IOobject& io,
    const point& point1,
    point& point2,
    const scalar innerRadius,
    const scalar outerRadius,
    const dictionary& dict
)
:
    searchableSurface(io),
    innerCylinder_(io, point1, point2, innerRadius, dict, false),
    outerCylinder_(io, point1, point2, outerRadius, dict, true)
{
    //Need to reset to take account of transforms
    point1_ = innerCylinder_.start();
    point2_ = innerCylinder_.end();
    innerRadius_ = innerCylinder_.radius();
    outerRadius_ = outerCylinder_.radius();

    magDir_  = mag(point2_-point1_);
    unitDir_ = (point2_-point1_)/magDir_;

    bounds() = calcBounds();
}


Foam::searchableRing::searchableRing
(
    const IOobject& io,
    const dictionary& dict
)
:
    searchableSurface(io),
    point1_(dict.lookup("point1")),
    point2_(endPoint(dict)),
    innerCylinder_(io, point1_, point2_, readScalar(dict.lookup("innerRadius")), dict, false),
    outerCylinder_(io, point1_, point2_, readScalar(dict.lookup("outerRadius")), dict, true)
{
    //Need to reset to take account of transforms
    point1_ = innerCylinder_.start();
    point2_ = innerCylinder_.end();
    innerRadius_ = innerCylinder_.radius();
    outerRadius_ = outerCylinder_.radius();

    magDir_ = mag(point2_-point1_);
    unitDir_ = (point2_-point1_)/magDir_;

    bounds() = calcBounds();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableRing::~searchableRing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::wordList& Foam::searchableRing::regions() const
{
    if (regions_.empty())
    {
        regions_.setSize(1);
        regions_[0] = "region0";
    }
    return regions_;
}


void Foam::searchableRing::findNearest
(
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info,
    const bool threaded /*= false*/
) const
{
    info.setSize(samples.size());

    //TODO threaded
    forAll(samples, i)
    {
        info[i] = findNearest(samples[i], nearestDistSqr[i]);
    }
}


void Foam::searchableRing::findLine
(
    const pointField& start,
    const pointField& end,
    List<pointIndexHit>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        // Pick nearest intersection. If none intersected take second one.
        pointIndexHit b;
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableRing::findLineAny
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
        // Discard far intersection
        pointIndexHit b;
        findLineAll(start[i], end[i], info[i], b);
        if (!info[i].hit() && b.hit())
        {
            info[i] = b;
        }
    }
}


void Foam::searchableRing::findLineAll
(
    const pointField& start,
    const pointField& end,
    List<List<pointIndexHit>>& info
) const
{
    info.setSize(start.size());

    forAll(start, i)
    {
        pointIndexHit near, far;
        findLineAll(start[i], end[i], near, far);

        if (near.hit())
        {
            if (far.hit())
            {
                info[i].setSize(2);
                info[i][0] = near;
                info[i][1] = far;
            }
            else
            {
                info[i].setSize(1);
                info[i][0] = near;
            }
        }
        else
        {
            if (far.hit())
            {
                info[i].setSize(1);
                info[i][0] = far;
            }
            else
            {
                info[i].clear();
            }
        }
    }
}


void Foam::searchableRing::getRegion
(
    const List<pointIndexHit>& info,
    labelList& region,
    const bool threaded /*= false*/
) const
{
    region.setSize(info.size());
    region = 0;
}

void Foam::searchableRing::getCurvature
(
    const List<pointIndexHit>& info,
    scalarField& curvature
) const
{
    // not currently implemented
    curvature.setSize(info.size());
    curvature = 0.;
}

void Foam::searchableRing::getNormal
(
    const List<pointIndexHit>& info,
    vectorField& normal
) const
{
    normal.setSize(info.size());
    normal = vector::zero;

    scalar midRadius = (innerRadius_ + outerRadius_) / 2.;

    forAll(info, i)
    {
        if (info[i].hit())
        {
            vector v(info[i].hitPoint() - point1_);
            // Decompose sample-point1 into normal and parallel component
            scalar parallel = (v & unitDir_);

            // Remove the parallel component
            v -= parallel*unitDir_;
            scalar magV = mag(v);

            if (parallel <= SMALL)
            {
                normal[i] = -unitDir_;
            }
            else if (parallel >= magDir_ - SMALL)
            {
                normal[i] = unitDir_;
            }
            else
            {
                if (magV > midRadius)
                {
                    normal[i] = v/magV;
                }
                else
                {
                    normal[i] = -v/magV;
                }
            }
        }
    }
}


void Foam::searchableRing::getVolumeType
(
    const pointField& points,
    List<volumeType>& volType,
    const bool threaded /*= false*/
) const
{
    volType.setSize(points.size());
    volType = volumeType::INSIDE;

    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        vector v(pt - point1_);

        // Decompose sample-point1 into normal and parallel component
        scalar parallel = v & unitDir_;

        if (parallel < 0)
        {
            // left of point1 endcap
            volType[pointI] = volumeType::OUTSIDE;
        }
        else if (parallel > magDir_)
        {
            // right of point2 endcap
            volType[pointI] = volumeType::OUTSIDE;
        }
        else
        {
            // Remove the parallel component
            v -= parallel*unitDir_;

            if (mag(v) > outerRadius_ || mag(v) < innerRadius_)
            {
                volType[pointI] = volumeType::OUTSIDE;
            }
            else
            {
                volType[pointI] = volumeType::INSIDE;
            }
        }
    }
}


// ************************************************************************* //
