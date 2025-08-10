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
    Foam::fvWallPointData

Description
    Holds information (coordinate and normal) regarding nearest wall point.

    Is like fvWallPoint but transfer extra (passive) data.
    Used e.g. in wall distance calculation with wall reflection vectors.

SourceFiles
    fvWallPointDataI.H
    fvWallPointData.C

\*---------------------------------------------------------------------------*/

#ifndef fvWallPointData_H
#define fvWallPointData_H

#include "fvMesh/wallDist/fvWallPoint/fvWallPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class fvWallPointData;

// Forward declaration of friend functions and operators

template<class Type>
Istream& operator>>(Istream&, fvWallPointData<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const fvWallPointData<Type>&);


/*---------------------------------------------------------------------------*\
                       Class fvWallPointData Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fvWallPointData
:
    public fvWallPoint
{
    // Private data

        //- Data at nearest wall center
        Type data_;


    // Private Member Functions

        //- Evaluate distance to point. Update distSqr, origin from whomever
        //  is nearer point. Return true if w2 is closer to point,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const point&,
            const fvWallPointData<Type>& w2,
            const scalar tol,
            TrackingData& td
        );


public:

    typedef Type dataType;


    // Constructors

        //- Construct null
        inline fvWallPointData();

        //- Construct from origin, normal and square distance
        inline fvWallPointData
        (
            const point& origin,
            const Type& data,
            const scalar distSqr
        );


    // Member Functions

        // Access

            inline const Type& data() const;

            inline Type& data();

            template<class TrackingData>
            inline const Type& data(TrackingData& td) const;


        // Needed by meshWave

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const fvPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            //  Calls update(...) with cellCentre of celli.
            template<class TrackingData>
            inline bool updateCell
            (
                const fvMesh& mesh,
                const label thisCelli,
                const labelPair& neighbourPatchAndFacei,
                const fvWallPointData<Type>& neighbourWallInfo,
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
                const fvWallPointData<Type>& neighbourWallInfo,
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
                const fvWallPointData<Type>& neighbourWallInfo,
                const scalar tol,
                TrackingData& td
            );


    // IOstream Operators

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const fvWallPointData<Type>&
        );

        friend Istream& operator>> <Type>
        (
            Istream&,
            fvWallPointData<Type>&
        );
};


//- Data associated with fvWallPointData type are contiguous. List the usual
//  ones.

template<>
constexpr bool contiguous<fvWallPointData<bool>>()
{
    return contiguous<fvWallPoint>();
}
template<>
constexpr bool contiguous<fvWallPointData<label>>()
{
    return contiguous<fvWallPoint>();
}
template<>
constexpr bool contiguous<fvWallPointData<scalar>>()
{
    return contiguous<fvWallPoint>();
}
template<>
constexpr bool contiguous<fvWallPointData<vector>>()
{
    return contiguous<fvWallPoint>();
}
template<>
constexpr bool contiguous<fvWallPointData<sphericalTensor>>()
{
    return contiguous<fvWallPoint>();
}
template<>
constexpr bool contiguous<fvWallPointData<symmTensor>>()
{
    return contiguous<fvWallPoint>();
}
template<>
constexpr bool contiguous<fvWallPointData<tensor>>()
{
    return contiguous<fvWallPoint>();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMesh/wallDist/fvWallPointData/fvWallPointData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fvMesh/wallDist/fvWallPointData/fvWallPointDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
