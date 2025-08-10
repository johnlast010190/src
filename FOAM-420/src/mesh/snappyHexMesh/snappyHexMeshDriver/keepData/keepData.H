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

Class
    keepData

Description
    Holds information on visited cells.

SourceFiles
    keepDataI.H
    keepData.C

\*---------------------------------------------------------------------------*/

#ifndef keepData_H
#define keepData_H

#include "primitives/ints/label/label.H"
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

class keepData;

Ostream& operator<<(Ostream&, const keepData&);
Istream& operator>>(Istream&, keepData&);

/*---------------------------------------------------------------------------*\
                           Class keepData Declaration
\*---------------------------------------------------------------------------*/

class keepData
{
    // Private data

        //- position of nearest wall center
        label found_;

public:

    // Static data members

    // Constructors

        //- Construct null
        inline keepData() = default;

        //- Construct from label
        inline keepData
        (
            const label& found
        );


    // Member Functions

        // Access

            inline const label& found() const;

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
                const keepData&,
                const scalar,
                TrackingData& td
            ) const;

            //- Apply rotation matrix to any coordinates
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
                const label thisCellI,
                const label neighbourFaceI,
                const keepData& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring cell.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const label neighbourCellI,
                const keepData& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same face.
            template<class TrackingData>
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const keepData& neighbourInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of neighbouring face.
            template<class TrackingData>
            inline bool update
            (
                const keepData& neighbourInfo,
                TrackingData& td
            );

            //- Same (like operator==)
           template<class TrackingData>
           inline bool equal(const keepData&, TrackingData& td) const;

    // Member Operators

        //Note: Used to determine whether to call update.
        inline bool operator==(const keepData&) const;

        inline bool operator!=(const keepData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const keepData&);
        friend Istream& operator>>(Istream&, keepData&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "snappyHexMeshDriver/keepData/keepDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
