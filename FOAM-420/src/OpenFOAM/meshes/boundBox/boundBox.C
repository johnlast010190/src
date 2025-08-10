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
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/boundBox/boundBox.H"
#include "db/IOstreams/Pstreams/PstreamReduceOps.H"
#include "memory/tmp/tmp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::boundBox Foam::boundBox::greatBox
(
    point::uniform(-ROOTVGREAT),
    point::uniform(ROOTVGREAT)
);

const Foam::boundBox Foam::boundBox::invertedBox
(
    point::uniform(ROOTVGREAT),
    point::uniform(-ROOTVGREAT)
);

const Foam::faceList Foam::boundBox::faces
({
    // Point and face order as per hex cellmodel
    face{0, 4, 7, 3}, // x-min
    face{1, 2, 6, 5}, // x-max
    face{0, 1, 5, 4}, // y-min
    face{3, 7, 6, 2}, // y-max
    face{0, 3, 2, 1}, // z-min
    face{4, 5, 6, 7}  // z-max
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundBox::boundBox(const UList<point>& points, bool doReduce)
:
    min_(invertedBox.min()),
    max_(invertedBox.max())
{
    add(points);

    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBox::boundBox(const tmp<pointField>& tpoints, bool doReduce)
:
    min_(invertedBox.min()),
    max_(invertedBox.max())
{
    add(tpoints);

    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBox::boundBox
(
    const UList<point>& points,
    const labelUList& indices,
    bool doReduce
)
:
    min_(invertedBox.min()),
    max_(invertedBox.max())
{
    add(points, indices);

    if (doReduce)
    {
        reduce();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::boundBox::points() const
{
    tmp<pointField> tpoints = tmp<pointField>(new pointField(8));
    pointField& pt = tpoints.ref();

    pt[0] = min_;                                   // min-x, min-y, min-z
    pt[1] = point(max_.x(), min_.y(), min_.z());    // max-x, min-y, min-z
    pt[2] = point(max_.x(), max_.y(), min_.z());    // max-x, max-y, min-z
    pt[3] = point(min_.x(), max_.y(), min_.z());    // min-x, max-y, min-z
    pt[4] = point(min_.x(), min_.y(), max_.z());    // min-x, min-y, max-z
    pt[5] = point(max_.x(), min_.y(), max_.z());    // max-x, min-y, max-z
    pt[6] = max_;                                   // max-x, max-y, max-z
    pt[7] = point(min_.x(), max_.y(), max_.z());    // min-x, max-y, max-z

    return tpoints;
}


void Foam::boundBox::inflate(const scalar s)
{
    const vector ext = vector::one*s*mag();

    min_ -= ext;
    max_ += ext;
}


void Foam::boundBox::inflateAbs(const scalar s)
{
    const vector ext = vector::one*s;

    min_ -= ext;
    max_ += ext;
}


void Foam::boundBox::reduce()
{
    Foam::reduce(
        std::tie(min_, max_),
        ParallelOp<minOp<point>, maxOp<point>>{}
    );
}


bool Foam::boundBox::intersect(const boundBox& bb)
{
    min_ = ::Foam::max(min_, bb.min_);
    max_ = ::Foam::min(max_, bb.max_);

    return !empty();
}


bool Foam::boundBox::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    forAll(points, i)
    {
        if (!contains(points[i]))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBox::contains
(
    const UList<point>& points,
    const labelUList& indices
) const
{
    if (points.empty() || indices.empty())
    {
        return true;
    }

    forAll(indices, i)
    {
        if (!contains(points[indices[i]]))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBox::containsAny(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    forAll(points, i)
    {
        if (contains(points[i]))
        {
            return true;
        }
    }

    return false;
}


bool Foam::boundBox::containsAny
(
    const UList<point>& points,
    const labelUList& indices
) const
{
    if (points.empty() || indices.empty())
    {
        return true;
    }

    forAll(indices, i)
    {
        if (contains(points[indices[i]]))
        {
            return true;
        }
    }

    return false;
}


Foam::point Foam::boundBox::nearest(const point& pt) const
{
    // Clip the point to the range of the bounding box
    const scalar surfPtx = Foam::max(Foam::min(pt.x(), max_.x()), min_.x());
    const scalar surfPty = Foam::max(Foam::min(pt.y(), max_.y()), min_.y());
    const scalar surfPtz = Foam::max(Foam::min(pt.z(), max_.z()), min_.z());

    point closest(surfPtx, surfPty, surfPtz);

    //Now move internal points
    scalar minDisp = GREAT;
    label minDir = label(-1);

    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        scalar upper(max_[dir]-closest[dir]);
        scalar lower(closest[dir]-min_[dir]);

        if (upper < 0 || lower < 0)
        {
            continue;
        }

        if
        (
            upper <= lower
            && upper < Foam::mag(minDisp)
        )
        {
            minDisp = upper;
            minDir = dir;
        }
        else if
        (
            lower < Foam::mag(minDisp)
        )
        {
            minDisp = -lower;
            minDir = dir;
        }
    }

    if (minDir != -1)
    {
        closest[minDir] += minDisp;
    }

    return closest;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBox& bb)
{
    if (os.format() == IOstream::ASCII)
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundBox& bb)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> bb.min_ >> bb.max_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
