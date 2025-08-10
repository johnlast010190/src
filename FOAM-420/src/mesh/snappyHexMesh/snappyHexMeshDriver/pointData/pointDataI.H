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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"
#include "primitives/transform/transform.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::pointData::pointData()
:
    pointEdgePoint(),
    s_(GREAT),
    v_(point::max)
{}


// Construct from origin, distance
inline Foam::pointData::pointData
(
    const point& origin,
    const scalar distSqr,
    const scalar s,
    const vector& v
)
:
    pointEdgePoint(origin, distSqr),
    s_(s),
    v_(v)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::pointData::s() const
{
    return s_;
}


inline const Foam::vector& Foam::pointData::v() const
{
    return v_;
}


template<class TrackingData>
inline void Foam::pointData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{
    pointEdgePoint::transform(patch, patchFacei, transform, td);
    v_ = transform.transform(v_);
}


// Update this with information from connected edge
template<class TrackingData>
inline bool Foam::pointData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgeI,
    const pointData& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if
    (
        pointEdgePoint::updatePoint
        (
            mesh,
            pointi,
            edgeI,
            edgeInfo,
            tol,
            td
        )
    )
    {
        s_ = edgeInfo.s_;
        v_ = edgeInfo.v_;
        return true;
    }
    else
    {
        return false;
    }
}

// Update this with new information on same point
template<class TrackingData>
inline bool Foam::pointData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointData& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if
    (
        pointEdgePoint::updatePoint
        (
            mesh,
            pointi,
            newPointInfo,
            tol,
            td
        )
    )
    {
        s_ = newPointInfo.s_;
        v_ = newPointInfo.v_;
        return true;
    }
    else
    {
        return false;
    }
}


// Update this with new information on same point. No extra information.
template<class TrackingData>
inline bool Foam::pointData::updatePoint
(
    const pointData& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (pointEdgePoint::updatePoint(newPointInfo, tol, td))
    {
        s_ = newPointInfo.s_;
        v_ = newPointInfo.v_;
        return true;
    }
    else
    {
        return false;
    }
}


// Update this with information from connected point
template<class TrackingData>
inline bool Foam::pointData::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointi,
    const pointData& pointInfo,
    const scalar tol,
    TrackingData& td

)
{
    if
    (
        pointEdgePoint::updateEdge
        (
            mesh,
            edgeI,
            pointi,
            pointInfo,
            tol,
            td
        )
    )
    {
        s_ = pointInfo.s_;
        v_ = pointInfo.v_;
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointData::operator==(const Foam::pointData& rhs)
const
{
    return
        pointEdgePoint::operator==(rhs)
     && (s() == rhs.s())
     && (v() == rhs.v());
}


inline bool Foam::pointData::operator!=(const Foam::pointData& rhs)
const
{
    return !(*this == rhs);
}


// ************************************************************************* //
