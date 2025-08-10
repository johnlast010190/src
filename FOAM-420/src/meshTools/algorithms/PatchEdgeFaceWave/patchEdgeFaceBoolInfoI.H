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
    (c) 2018 Esi Ltd.
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"
#include "primitives/transform/transform.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::patchEdgeFaceBoolInfo::patchEdgeFaceBoolInfo()
:
    origin_(false)
{}


// Construct from origin, distance
inline Foam::patchEdgeFaceBoolInfo::patchEdgeFaceBoolInfo
(
    const bool& origin
)
:
    origin_(origin)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const bool& Foam::patchEdgeFaceBoolInfo::origin() const
{
    return origin_;
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceBoolInfo::valid(TrackingData& td) const
{
    return origin_;
}


template<class TrackingData>
inline void Foam::patchEdgeFaceBoolInfo::transform
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const tensor& rotTensor,
    const scalar tol,
    TrackingData& td
)
{
    origin_ = true;
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceBoolInfo::updateEdge
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const label edgeI,
    const label faceI,
    const patchEdgeFaceBoolInfo& faceInfo,
    const scalar tol,
    TrackingData& td
)
{
    origin_ = faceInfo.origin_;
    return true;
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceBoolInfo::updateEdge
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const patchEdgeFaceBoolInfo& edgeInfo,
    const bool sameOrientation,
    const scalar tol,
    TrackingData& td
)
{
    origin_ = edgeInfo.origin_;
    return true;
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceBoolInfo::updateFace
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& patch,
    const label faceI,
    const label edgeI,
    const patchEdgeFaceBoolInfo& edgeInfo,
    const scalar tol,
    TrackingData& td
)
{
    origin_ = edgeInfo.origin_;
    return true;
}


template<class TrackingData>
inline bool Foam::patchEdgeFaceBoolInfo::equal
(
    const patchEdgeFaceBoolInfo& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::patchEdgeFaceBoolInfo::operator==
(
    const Foam::patchEdgeFaceBoolInfo& rhs
) const
{
    return origin() == rhs.origin();
}


inline bool Foam::patchEdgeFaceBoolInfo::operator!=
(
    const Foam::patchEdgeFaceBoolInfo& rhs
) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
