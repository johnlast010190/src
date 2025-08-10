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
    (c) 2012-2016 OpenFOAM Foundation

Class
    Foam::pointToPointPlanarInterpolation

Description
    Interpolates between two sets of unstructured points using 2D Delaunay
    triangulation. Used in e.g. timeVaryingMapped bcs.

SourceFiles
    pointToPointPlanarInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef pointToPointPlanarInterpolation_H
#define pointToPointPlanarInterpolation_H

#include "containers/Lists/FixedList/FixedList.H"
#include "coordinate/systems/coordinateSystem.H"
#include "db/Time/instant/instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class pointToPointPlanarInterpolation Declaration
\*---------------------------------------------------------------------------*/

class pointToPointPlanarInterpolation
{
    // Private data

        //- Perturbation factor
        const scalar perturb_;

        //- Whether to use delabella delaunay triangulation
        const bool delabella_;

        //- Whether to use nearest point only (avoids triangulation, projection)
        const bool nearestOnly_;

        //- Coordinate system
        coordinateSystem referenceCS_;

        //- Number of source points (for checking)
        label nPoints_;

        //- Current interpolation addressing to face centres of underlying
        //  patch
        List<FixedList<label, 3>> nearestVertex_;

        //- Current interpolation factors to face centres of underlying
        //  patch
        List<FixedList<scalar, 3>> nearestVertexWeight_;


    // Private Member Functions

        //- Calculate a local coordinate system from set of points
        coordinateSystem calcCoordinateSystem(const pointField&) const;

        //- Calculate addressing and weights
        void calcWeights
        (
            const pointField& sourcePoints,
            const pointField& destPoints
        );


public:

    // Declare name of the class and its debug switch
    ClassName("pointToPointPlanarInterpolation");


    // Constructors

        //- Construct from 3D locations. Determines local coordinate system
        //  from sourcePoints and maps onto that. If nearestOnly skips any
        //  local coordinate system and triangulation and uses nearest vertex
        //  only
        pointToPointPlanarInterpolation
        (
            const pointField& sourcePoints,
            const pointField& destPoints,
            const scalar perturb,
            const bool delabella,
            const bool nearestOnly = false
        );

        //- Construct from coordinate system and locations.
        pointToPointPlanarInterpolation
        (
            const coordinateSystem& referenceCS,
            const pointField& sourcePoints,
            const pointField& destPoints,
            const scalar perturb,
            const bool delabella
        );

        //- Construct from components
        pointToPointPlanarInterpolation
        (
            const scalar perturb,
            const bool delabella,
            const bool nearestOnly,
            const coordinateSystem& referenceCS,
            const label sourceSize,
            const List<FixedList<label, 3>>& nearestVertex,
            const List<FixedList<scalar, 3>>& nearestVertexWeight
        );

        //- Construct and return a clone
        autoPtr<pointToPointPlanarInterpolation> clone() const;


    // Member Functions

        //- Perturbation factor (for triangulation)
        scalar perturb() const
        {
            return perturb_;
        }

        //- Whether to use nearest point only (avoids triangulation, projection)
        bool nearestOnly() const
        {
            return nearestOnly_;
        }

        //- Return the coordinateSystem
        const coordinateSystem& referenceCS() const
        {
            return referenceCS_;
        }

        //- Number of source points
        label sourceSize() const
        {
            return nPoints_;
        }

        //  patch
        const List<FixedList<label, 3>>& nearestVertex() const
        {
            return nearestVertex_;
        }

        //- Current interpolation factors to face centres of underlying
        //  patch
        const List<FixedList<scalar, 3>>& nearestVertexWeight() const
        {
            return nearestVertexWeight_;
        }

        //- Helper: extract words of times
        static wordList timeNames(const instantList&);

        //- Helper: find time. Return true if succesful.
        static bool findTime
        (
            const instantList& times,
            const label startSampleTime,
            const scalar timeVal,
            label& lo,
            label& hi
        );

        //- Interpolate from field on source points to dest points
        template<class Type>
        tmp<Field<Type>> interpolate(const Field<Type>& sourceFld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "triSurface/triSurfaceTools/pointToPointPlanarInterpolationTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
