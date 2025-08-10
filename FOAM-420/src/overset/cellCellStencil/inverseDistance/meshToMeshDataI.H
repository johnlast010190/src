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
    (c) 2017 OpenCFD Ltd.
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline Foam::meshToMeshData::meshToMeshData()
:
    tgtCelli_(-2)
{}


// Construct from components
inline Foam::meshToMeshData::meshToMeshData(const label tgtCelli)
:
    tgtCelli_(tgtCelli)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool Foam::meshToMeshData::valid(TrackingData& td) const
{
    return tgtCelli_ != -2;
}


// No geometric data so never any problem on cyclics
template<class TrackingData>
inline bool Foam::meshToMeshData::sameGeometry
(
    const polyMesh&,
    const meshToMeshData&,
    const scalar,
    TrackingData& td
) const
{
    return true;
}


template<class TrackingData>
inline void Foam::meshToMeshData::translate
(
    const polyMesh&,
    const vector& rotTensor,
    TrackingData& td
)
{}


template<class TrackingData>
inline void Foam::meshToMeshData::transform
(
    const polyPatch& patch,
    const label patchFacei,
    const transformer& transform,
    TrackingData& td
)
{}


// Update cell with neighbouring face information
template<class TrackingData>
inline bool Foam::meshToMeshData::updateCell
(
    const polyMesh& mesh,
    const label thisCelli,
    const label neighbourFacei,
    const meshToMeshData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        const point& cc = mesh.cellCentres()[thisCelli];

        // See if we are in the neighbour
        label tgti = neighbourInfo.tgtCelli_;

        // Try match of neighbour
        if (td.tgtMesh_.pointInCell(cc, tgti, polyMesh::CELL_TETS))
        {
            tgtCelli_ = tgti;
            return true;
        }

        // Try match of neighbour's neighbours
        const labelList& tgtNbrs = td.tgtMesh_.cellCells(tgti);
        forAll(tgtNbrs, i)
        {
            label tgti = tgtNbrs[i];

            // Try match of neighbour
            if (td.tgtMesh_.pointInCell(cc, tgti, polyMesh::CELL_TETS))
            {
                tgtCelli_ = tgti;
                return true;
            }
        }

        // Leave unmatched cell for visit from different face
    }
    return false;
}


// Update face with neighbouring cell information
template<class TrackingData>
inline bool Foam::meshToMeshData::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const label neighbourCelli,
    const meshToMeshData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        tgtCelli_ = neighbourInfo.tgtCelli_;
        return true;
    }
    else
    {
        return false;
    }
}


// Update face with coupled face information
template<class TrackingData>
inline bool Foam::meshToMeshData::updateFace
(
    const polyMesh& mesh,
    const label thisFacei,
    const meshToMeshData& neighbourInfo,
    const scalar tol,
    TrackingData& td
)
{
    if (!valid(td))
    {
        tgtCelli_ = neighbourInfo.tgtCelli_;
        return true;
    }
    else
    {
        return false;
    }
}


template<class TrackingData>
inline bool Foam::meshToMeshData::equal
(
    const meshToMeshData& rhs,
    TrackingData& td
) const
{
    if (!valid(td))
    {
        if (!rhs.valid(td))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return operator==(rhs);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::meshToMeshData::operator==
(
    const Foam::meshToMeshData& rhs
)
 const
{
    return tgtCelli_ == rhs.tgtCelli_;
}


inline bool Foam::meshToMeshData::operator!=
(
    const Foam::meshToMeshData& rhs
)
 const
{
    return !(*this == rhs);
}


// ************************************************************************* //
