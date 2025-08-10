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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2019 OpenCFD Ltd
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::pointTopoDistanceData::pointTopoDistanceData()
:
    data_(-1),
    distance_(-1)
{}


inline Foam::pointTopoDistanceData::pointTopoDistanceData
(
    const label data,
    const label distance
)
:
    data_(data),
    distance_(distance)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::pointTopoDistanceData::valid(TrackingData& td) const
{
    return distance_ != -1;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::pointTopoDistanceData::sameGeometry
(
    const pointTopoDistanceData&,
    const scalar tol,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::pointTopoDistanceData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


// Update this with information from connected edge
template<class TrackingData>
inline bool Foam::pointTopoDistanceData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const label edgeI,
    const pointTopoDistanceData& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        data_ = edgeInfo.data_;
        distance_ = edgeInfo.distance_ + 1;
        return true;
    }

    return false;
}


// Update this with new information on same point
template<class TrackingData>
inline bool Foam::pointTopoDistanceData::updatePoint
(
    const polyMesh& mesh,
    const label pointi,
    const pointTopoDistanceData& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        operator=(newPointInfo);
        return true;
    }

    return false;
}


// Update this with new information on same point. No extra information.
template<class TrackingData>
inline bool Foam::pointTopoDistanceData::updatePoint
(
    const pointTopoDistanceData& newPointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        operator=(newPointInfo);
        return true;
    }

    return false;
}


// Update this with information from connected point
template<class TrackingData>
inline bool Foam::pointTopoDistanceData::updateEdge
(
    const polyMesh& mesh,
    const label edgeI,
    const label pointi,
    const pointTopoDistanceData& pointInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (distance_ == -1)
    {
        operator=(pointInfo);
        return true;
    }

    return false;
}


template<class TrackingData>
inline bool Foam::pointTopoDistanceData::equal
(
    const pointTopoDistanceData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::pointTopoDistanceData::operator==
(
    const Foam::pointTopoDistanceData& rhs
) const
{
    return data() == rhs.data() && distance() == rhs.distance();
}


inline bool Foam::pointTopoDistanceData::operator!=
(
    const Foam::pointTopoDistanceData& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
