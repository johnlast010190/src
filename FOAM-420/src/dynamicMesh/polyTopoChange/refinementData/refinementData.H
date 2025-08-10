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
    Foam::refinementData

Description
    Transfers refinement levels such that slow transition between levels is
    maintained. Used in FaceCellWave.

SourceFiles
    refinementDataI.H
    refinementData.C

\*---------------------------------------------------------------------------*/

#ifndef refinementData_H
#define refinementData_H

#include "meshes/primitiveShapes/point/point.H"
#include "primitives/Tensor/tensor/tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class polyMesh;
class transformer;


// Forward declaration of friend functions and operators

class refinementData;

Istream& operator>>(Istream&, refinementData&);
Ostream& operator<<(Ostream&, const refinementData&);


/*---------------------------------------------------------------------------*\
                           Class refinementData Declaration
\*---------------------------------------------------------------------------*/

class refinementData
{
    // Private data

        //- Count which triggers refinement
        label refinementCount_;

        //- Refinement level
        label count_;

public:

    // Constructors

        //- Construct null
        inline refinementData();

        //- Construct from count
        inline refinementData(const label refinementCount, const label count);


    // Member Functions

        // Access

            inline label refinementCount() const
            {
                return refinementCount_;
            }

            inline label& refinementCount()
            {
                return refinementCount_;
            }

            inline label count() const
            {
                return count_;
            }

            inline label& count()
            {
                return count_;
            }

            inline bool isRefined() const
            {
                return count_ >= refinementCount_;
            }



        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const polyMesh&,
                const refinementData&,
                const scalar,
                TrackingData& td
            ) const;

            //- Apply translation to any coordinates
            template<class TrackingData>
            inline void translate
            (
                const polyMesh&,
                const vector&,
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
                const refinementData& neighbourInfo,
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
                const refinementData& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFacei,
                const refinementData& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const refinementData&, TrackingData& td) const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const refinementData&) const;

        inline bool operator!=(const refinementData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const refinementData&);
        friend Istream& operator>>(Istream&, refinementData&);
};


//- Data associated with refinementData type are contiguous
template<>
constexpr bool contiguous<refinementData>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyTopoChange/refinementData/refinementDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
