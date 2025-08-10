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

Class
    Foam::fvWallPointYPlus

Description
    Holds information (coordinate and yStar) regarding nearest wall point.

    Used in VanDriest wall damping where the interest is in y+ but only
    needs to be calculated up to e.g. y+ < 200. In all other cells/faces
    the damping function becomes 1, since y gets initialized to GREAT and
    yStar to 1.

    Note: should feed the additional argument (yPlusCutoff) through as a
    template argument into FvFaceCellWave.

SourceFiles
    fvWallPointYPlusI.H
    fvWallPointYPlus.C

\*---------------------------------------------------------------------------*/

#ifndef fvWallPointYPlus_H
#define fvWallPointYPlus_H

#include "fvMesh/wallDist/fvWallPointData/fvWallPointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fvWallPointYPlus Declaration
\*---------------------------------------------------------------------------*/

class fvWallPointYPlus
:
    public fvWallPointData<scalar>
{
    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer point. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point& pt,
            const fvWallPointYPlus& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    // Static data members

        //- Cut-off value for y+
        static scalar yPlusCutOff;


    // Constructors

        //- Construct null
        inline fvWallPointYPlus();

        //- Construct from origin, yStar, distance
        inline fvWallPointYPlus
        (
            const point& origin,
            const scalar yStar,
            const scalar distSqr
        );


    // Member Functions

        // Needed by FvFaceCellWave

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of celli.
            template<class TrackingData>
            inline bool updateCell
            (
                const fvMesh& mesh,
                const label thisCelli,
                const labelPair& neighbourPatchAndFacei,
                const fvWallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            //  Calls update(...) with faceCentre of facei.
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const label neighbourCelli,
                const fvWallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            //  Merge new and old info.
            //  Calls update(...) with faceCentre of facei.
            template<class TrackingData>
            inline bool updateFace
            (
                const fvMesh& mesh,
                const labelPair& thisPatchAndFacei,
                const fvWallPointYPlus& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );
};


//- Data associated with pointEdgePoint type as contiguous as underlying type
template<>
constexpr bool contiguous<fvWallPointYPlus>()
{
    return contiguous<fvWallPointData<scalar>>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh/wallDist/fvWallPointYPlus/fvWallPointYPlusI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
