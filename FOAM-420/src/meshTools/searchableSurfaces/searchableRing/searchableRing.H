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
    (c) 2010-2012 Esi Ltd.
    (c) 2011 OpenFOAM Foundation

Class
    Foam::searchableRing

Description
    Searching on ring

SourceFiles
    searchableRing.C

\*---------------------------------------------------------------------------*/

#ifndef searchableRing_H
#define searchableRing_H

#include "meshes/treeBoundBox/treeBoundBox.H"
#include "searchableSurfaces/searchableSurface/searchableSurface.H"
#include "searchableSurfaces/searchableCylinder/searchableCylinder.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class searchableRing Declaration
\*---------------------------------------------------------------------------*/

class searchableRing
:
    public searchableSurface
{
private:

    // Private Member Data

        //- 'left' point
        point point1_;

        //- 'right' point
        point point2_;

        //- Radius squared
        scalar innerRadius_;

        //- Radius squared
        scalar outerRadius_;

        //- length of vector point2-point1
        scalar magDir_;

        //- normalised vector point2-point1
        vector unitDir_;

        //- Names of regions
        mutable wordList regions_;

        //- Define inner cylinder with out endCaps
        const searchableCylinder innerCylinder_;

        //- Define outer cylinder with endCaps
        const searchableCylinder outerCylinder_;

    // Private Member Functions

        //- Find nearest point on ring.
        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

        scalar radius2(const point& pt) const;

        //- Find intersection with ring
        void findLineAll
        (
            const point& start,
            const point& end,
            pointIndexHit& near,
            pointIndexHit& far
        ) const;

        //- Return the boundBox of the ring
        boundBox calcBounds() const;

        //- Disallow default bitwise copy construct
        searchableRing(const searchableRing&);

        //- Disallow default bitwise assignment
        void operator=(const searchableRing&);

        //- Calculate cylinder end point
        point endPoint(const dictionary& dict) const;

public:

    //- Runtime type information
    TypeName("searchableRing");


    // Constructors

        //- Construct from components
        searchableRing
        (
            const IOobject& io,
            const point&,
            point&,
            const scalar innerRadius,
            const scalar outerRadius,
            const dictionary& dict
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableRing
        (
            const IOobject& io,
            const dictionary& dict
        );

    //- Destructor
    virtual ~searchableRing();


    // Member Functions

        virtual const wordList& regions() const;

        //- Whether supports volume type below
        virtual bool hasVolumeType() const
        {
            return true;
        }

        //- Range of local indices that can be returned.
        virtual label size() const
        {
            return 1;
        }

        //- Get representative set of element coordinates
        //  Usually the element centres (should be of length size()).
        virtual tmp<pointField> coordinates() const;

        //- Get bounding spheres (centre and radius squared), one per element.
        //  Any point on element is guaranteed to be inside.
        virtual void boundingSpheres
        (
            pointField& centres,
            scalarField& radiusSqr
        ) const;

        //- Get the points that define the surface.
        virtual tmp<pointField> points() const;

        //- Does any part of the surface overlap the supplied bound box?
        virtual bool overlaps(const boundBox& bb) const
        {
            NotImplemented;
        }

        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&,
                const bool threaded = false
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&,
                const bool threaded = false
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit>>&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region,
                const bool threaded = false
            ) const;

            //- From a set of points and indices get the normal
            virtual void getCurvature
            (
                const List<pointIndexHit>&,
                scalarField& curvature
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&,
                const bool threaded = false
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                NotImplemented;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
