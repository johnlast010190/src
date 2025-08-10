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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

Class
    Foam::pointTopoDistanceData

Description
    For use with PointEdgeWave. Determines topological distance to
    starting points

SourceFiles
    pointTopoDistanceDataI.H
    pointTopoDistanceData.C

\*---------------------------------------------------------------------------*/

#ifndef pointTopoDistanceData_H
#define pointTopoDistanceData_H

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

class pointTopoDistanceData;

Istream& operator>>(Istream&, pointTopoDistanceData&);
Ostream& operator<<(Ostream&, const pointTopoDistanceData&);


/*---------------------------------------------------------------------------*\
                       Class pointTopoDistanceData Declaration
\*---------------------------------------------------------------------------*/

class pointTopoDistanceData
{
    // Private data

        //- Starting data
        label data_;

        //- Distance
        label distance_;


public:

    // Constructors

        //- Construct null
        inline pointTopoDistanceData();

        //- Construct from count
        inline pointTopoDistanceData
        (
            const label data,
            const label distance
        );


    // Member Functions

        // Access


            inline label data() const
            {
                return data_;
            }
            inline label distance() const
            {
                return distance_;
            }


        // Needed by PointEdgeWave


            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            template<class TrackingData>
            inline bool valid(TrackingData& td) const;

            //- Check for identical geometrical data. Used for cyclics checking.
            template<class TrackingData>
            inline bool sameGeometry
            (
                const pointTopoDistanceData&,
                const scalar tol,
                TrackingData& td
            ) const;

            //- Transform across an interface
            template<class TrackingData>
            inline void transform
            (
                const polyPatch& patch,
                const label patchFacei,
                const transformer& transform,
                TrackingData& td
            );

            //- Influence of edge on point
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const label edgeI,
                const pointTopoDistanceData& edgeInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  Merge new and old info.
            template<class TrackingData>
            inline bool updatePoint
            (
                const polyMesh& mesh,
                const label pointi,
                const pointTopoDistanceData& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of different value on same point.
            //  No information about current position whatsoever.
            template<class TrackingData>
            inline bool updatePoint
            (
                const pointTopoDistanceData& newPointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Influence of point on edge.
            template<class TrackingData>
            inline bool updateEdge
            (
                const polyMesh& mesh,
                const label edgeI,
                const label pointi,
                const pointTopoDistanceData& pointInfo,
                const scalar tol,
                TrackingData& td
            );

            //- Same (like operator==)
            template<class TrackingData>
            inline bool equal(const pointTopoDistanceData&, TrackingData&)
            const;


    // Member Operators

        // Needed for List IO
        inline bool operator==(const pointTopoDistanceData&) const;
        inline bool operator!=(const pointTopoDistanceData&) const;

    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const pointTopoDistanceData&);
        friend Istream& operator>>(Istream&, pointTopoDistanceData&);
};


//- Data associated with pointTopoDistanceData type are contiguous
template<>
constexpr bool contiguous<pointTopoDistanceData>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshStructure/pointTopoDistanceDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
