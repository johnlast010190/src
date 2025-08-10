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
    (c) 2017-2023 Esi Ltd.

Class
    Foam::cellInfo

Description
    Holds information regarding type of cell. Used in inside/outside
    determination in cellClassification.

SourceFiles
    cellInfoI.H
    cellInfo.C

\*---------------------------------------------------------------------------*/

#ifndef cellInfo_H
#define cellInfo_H

#include "meshes/primitiveShapes/point/point.H"
#include "primitives/ints/label/label.H"
#include "primitives/Tensor/tensor/tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;


// Forward declaration of friend functions and operators

class cellInfo;

Istream& operator>>(Istream&, cellInfo&);
Ostream& operator<<(Ostream&, const cellInfo&);


/*---------------------------------------------------------------------------*\
                           Class cellInfo Declaration
\*---------------------------------------------------------------------------*/

class cellInfo
{
    // Private data

        label type_;

    // Private Member Functions

        //- Update current cell/face type with neighbouring
        //  type. Return true if information needs to propagate,
        //  false otherwise.
        template<class TrackingData>
        inline bool update
        (
            const cellInfo& w2,
            const label thisFacei,
            const label thisCelli,
            const label neighbourFacei,
            const label neighbourCelli,
            TrackingData& td
        );

public:

    // Constructors

        //- Construct null
        inline cellInfo();

        //- Construct from cType
        inline cellInfo(const label);


    // Member Functions

        // Access

            inline label type() const
            {
                return type_;
            }


        // Needed by meshWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const cellInfo&,
                const scalar,
                TrackingData& td
            ) const;

            //- Apply translation to any coordinates
            template<class TrackingData>
            inline void translate
            (
                const polyMesh&,
                const vector& transVec,
                TrackingData& td
            );

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCelli,
                const label neighbourFacei,
                const cellInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const label neighbourCelli,
                const cellInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const cellInfo& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const cellInfo&, TrackingData& td) const;

    // Member Operators

        //Note: Used to determine whether to call update.
        inline bool operator==(const cellInfo&) const;

        inline bool operator!=(const cellInfo&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const cellInfo&);
        friend Istream& operator>>(Istream&, cellInfo&);
};


//- Data associated with cellInfo type are contiguous
template<>
constexpr bool contiguous<cellInfo>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cellClassification/cellInfoI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
