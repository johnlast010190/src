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
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Update this with w2 (information on same point)
template<class TrackingData>
inline bool Foam::fieldData::update
(
    const fieldData& w2,
    TrackingData& td
)
{
    if (!valid(td))
    {
        origin_ = w2.origin();
        s_ = w2.s();

        return true;
    }

    scalar ratio = pressure_ / (w2.origin() + SMALL);

    if (ratio < 0.5 || s_)
    {
        // stop as at low pressure
        return false;
    }
    else
    {
        // update with new values
        origin_ = w2.origin();
        s_ = w2.s();

        return true;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::fieldData::fieldData()
:
    pressure_(-GREAT),
    origin_(-GREAT),
    s_(false)
{}


// Construct from origin and pressure
inline Foam::fieldData::fieldData
(
    const scalar pressure,
    const scalar origin,
    const bool s
)
:
    pressure_(pressure),
    origin_(origin),
    s_(s)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::fieldData::pressure() const
{
    return pressure_;
}


inline Foam::scalar Foam::fieldData::origin() const
{
    return origin_;
}


inline bool Foam::fieldData::s() const
{
    return s_;
}

template<class TrackingData>
inline bool Foam::fieldData::valid(TrackingData& td) const
{
    return origin_ != -GREAT;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::fieldData::sameGeometry
(
    const polyMesh&,
    const fieldData&,
    const scalar,
    TrackingData& td
) const
{
    return true;
}


// No geometric data.
template<class TrackingData>
inline void Foam::fieldData::translate
(
    const polyMesh&,
    const vector& transVec,
    TrackingData& td
)
{}


template<class TrackingData>
inline void Foam::fieldData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


// Update cell with neighbouring face information
template<class TrackingData>
inline bool Foam::fieldData::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const fieldData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
  return update(neighbourInfo, td);
}


// Update face with neighbouring cell information
template<class TrackingData>
inline bool Foam::fieldData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const fieldData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
  return update(neighbourInfo, td);
}

// Update face with coupled face information
template<class TrackingData>
inline bool Foam::fieldData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const fieldData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
  return update(neighbourInfo, td);
}


template <class TrackingData>
inline bool Foam::fieldData::equal
(
    const fieldData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::fieldData::operator==(const fieldData& rhs) const
{
    return origin() == rhs.origin();
}


inline bool Foam::fieldData::operator!=(const fieldData& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
