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
    (c) 1991-2007 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from label
inline Foam::keepData::keepData
(
    const label& found
)
:
    found_(found)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::label& Foam::keepData::found() const
{
    return found_;
}

template<class TrackingData>
inline bool Foam::keepData::valid(TrackingData& td) const
{
    return (found_ != 0);
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::keepData::sameGeometry
(
    const polyMesh&,
    const keepData&,
    const scalar,
    TrackingData& td
) const
{
    return true;
}


// No geometric data.
template<class TrackingData>
inline void Foam::keepData::translate
(
    const polyMesh&,
    const vector& transVec,
    TrackingData& td
)
{}


template<class TrackingData>
inline void Foam::keepData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


// Update cell with neighbouring face information
template<class TrackingData>
inline bool Foam::keepData::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const keepData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
  return update(neighbourInfo, td);
}


// Update face with neighbouring cell information
template<class TrackingData>
inline bool Foam::keepData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const keepData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
  return update(neighbourInfo, td);
}


// Update face with coupled face information
template<class TrackingData>
inline bool Foam::keepData::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const keepData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
  return update(neighbourInfo, td);
}

// Update face and cell information
template<class TrackingData>
inline bool Foam::keepData::update
(
    const keepData& neighbourInfo,
    TrackingData& td
)
{
    if (!valid(td))
    {
        if (!neighbourInfo.valid(td))
        {
             FatalErrorInFunction
                << "problem" << abort(FatalError);
        }
        if (neighbourInfo.found() > 1)
        {
            return false;
        }
        else
        {
            found_ = neighbourInfo.found();
            return true;
        }
    }
    else
    {
        return false;
    }
}

template <class TrackingData>
inline bool Foam::keepData::equal
(
    const keepData& rhs,
    TrackingData& td
) const
{
    return operator==(rhs);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::keepData::operator==(const keepData& rhs) const
{
    return (found() == rhs.found());
}


inline bool Foam::keepData::operator!=(const keepData& rhs) const
{
    return !(*this == rhs);
}


// ************************************************************************* //
