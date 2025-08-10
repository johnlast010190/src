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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class TrackingData>
inline bool fvWallPointYPlus::update
(
    const point& pt,
    const fvWallPointYPlus& w2,
    const scalar tol,
    TrackingData& td
)
{
    scalar dist2 = magSqr(pt - w2.origin());

    if (valid(td))
    {
        scalar diff = distSqr() - dist2;

        if (diff < 0)
        {
            // Already nearer to point
            return false;
        }

        if ((diff < SMALL) || ((distSqr() > SMALL) && (diff/distSqr() < tol)))
        {
            // Don't propagate small changes
            return false;
        }
    }

    // Either *this is not yet valid or w2 is closer
    {
        // Only propagate if interesting (i.e. y+ < 100)
        scalar yPlus = Foam::sqrt(dist2)/w2.data();

        if (yPlus < yPlusCutOff)
        {
            // Update with new values
            distSqr() = dist2;
            origin() = w2.origin();
            data() = w2.data();

            return true;
        }
        else
        {
            return false;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline fvWallPointYPlus::fvWallPointYPlus()
:
    fvWallPointData<scalar>()
{
    // Important: The value of yStar where the meshWave does not come.
    data() = 1.0;
}


inline fvWallPointYPlus::fvWallPointYPlus
(
    const point& origin,
    const scalar yStar,
    const scalar distSqr
)
:
    fvWallPointData<scalar>(origin, yStar, distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackingData>
inline bool fvWallPointYPlus::updateCell
(
    const fvMesh& mesh,
    const label thisCelli,
    const labelPair& neighbourPatchAndFacei,
    const fvWallPointYPlus& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    return update
    (
        mesh.C()[thisCelli],
        neighbourWallInfo,
        tol,
        td
    );
}


template<class TrackingData>
inline bool fvWallPointYPlus::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const label neighbourCelli,
    const fvWallPointYPlus& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    return update
    (
        thisPatchi == -1
      ? mesh.Cf()[thisFacei]
      : mesh.Cf().boundaryField()[thisPatchi][thisFacei],
        neighbourWallInfo,
        tol,
        td
    );
}


template<class TrackingData>
inline bool fvWallPointYPlus::updateFace
(
    const fvMesh& mesh,
    const labelPair& thisPatchAndFacei,
    const fvWallPointYPlus& neighbourWallInfo,
    const scalar tol,
    TrackingData& td
)
{
    const label thisPatchi = thisPatchAndFacei.first();
    const label thisFacei = thisPatchAndFacei.second();

    return update
    (
        thisPatchi == -1
      ? mesh.Cf()[thisFacei]
      : mesh.Cf().boundaryField()[thisPatchi][thisFacei],
        neighbourWallInfo,
        tol,
        td
    );
}


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
