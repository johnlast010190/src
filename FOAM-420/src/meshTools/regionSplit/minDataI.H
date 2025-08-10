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
    (c) 2015 OpenCFD Ltd.
    (c) 2014-2016 OpenFOAM Foundation
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::minData::minData()
:
    data_(labelMax)
{}


inline Foam::minData::minData(const label data)
:
    data_(data)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::minData::valid(TrackingData& td) const
{
    return data_ != labelMax;
}


template<class TrackingData>
inline bool Foam::minData::sameGeometry
(
    const polyMesh&,
    const minData&,
    const scalar,
    TrackingData&
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::minData::translate
(
    const polyMesh&,
    const vector& transVec,
    TrackingData& td
)
{}


template<class TrackingData>
inline void Foam::minData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


template<class TrackingData>
inline bool Foam::minData::updateCell
(
    const polyMesh&,
    const label thisCelli,
    const label neighbourFacei,
    const minData& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    if (neighbourInfo.data_ < data_)
    {
        operator=(neighbourInfo);
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::minData::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const minData& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    // From cell to its faces.

    if (neighbourInfo.data_ < data_)
    {
        operator=(neighbourInfo);
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::minData::updateFace
(
    const polyMesh&,
    const label thisFacei,
    const minData& neighbourInfo,
    const scalar tol,
    TrackingData&
)
{
    // From face to face (e.g. coupled faces)
    if (neighbourInfo.data_ < data_)
    {
        operator=(neighbourInfo);
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::minData::equal
(
    const minData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::minData::operator==
(
    const Foam::minData& rhs
) const
{
    return data() == rhs.data();
}


inline bool Foam::minData::operator!=
(
    const Foam::minData& rhs
) const
{
    return !(*this == rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::minData& wDist
)
{
    return os << wDist.data_;
}


inline Foam::Istream& Foam::operator>>
(
    Foam::Istream& is,
    Foam::minData& wDist
)
{
    return is >> wDist.data_;
}


// ************************************************************************* //
