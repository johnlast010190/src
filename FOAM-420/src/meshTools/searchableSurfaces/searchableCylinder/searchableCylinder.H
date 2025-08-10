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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::searchableCylinder

Description
    Searching on cylinder

SourceFiles
    searchableCylinder.C

\*---------------------------------------------------------------------------*/

#ifndef searchableCylinder_H
#define searchableCylinder_H

#include "meshes/treeBoundBox/treeBoundBox.H"
#include "searchableSurfaces/searchableSurface/searchableSurface.H"
#include "primitives/bools/Switch/Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class searchableCylinder Declaration
\*---------------------------------------------------------------------------*/

class searchableCylinder
:
    public searchableSurface
{
private:

    // Private Member Data

        //- 'left' point
        point point1_;

        //- Radius squared
        scalar radius_;

        //- Whether cylinder has end caps
        const Switch endCaps_;

        //- 'right' point
        point point2_;

        //- Length of vector point2-point1
        scalar magDir_;

        //- Normalised vector point2-point1
        vector unitDir_;

        //- Names of regions
        mutable wordList regions_;


    // Private Member Functions

        //- Inherit findNearest from searchableSurface
        using searchableSurface::findNearest;

        scalar radius2(const point& pt) const;

        //- Return the boundBox of the cylinder
        boundBox calcBounds() const;

        //- Disallow default bitwise copy construct
        searchableCylinder(const searchableCylinder&);

        //- Disallow default bitwise assignment
        void operator=(const searchableCylinder&);


public:

    //- Runtime type information
    TypeName("searchableCylinder");


    // Constructors

        //- Construct from components
        searchableCylinder
        (
            const IOobject& io,
            const point&,
            const point&,
            const scalar radius,
            const dictionary& dict,
            const Switch endCaps = true
        );

        //- Construct from dictionary (used by searchableSurface)
        searchableCylinder
        (
            const IOobject& io,
            const dictionary& dict
        );

    //- Destructor
    virtual ~searchableCylinder();


    // Member Functions

        //- Find intersection with cylinder
        void findLineAll
        (
            const point& start,
            const point& end,
            pointIndexHit& near,
            pointIndexHit& far
        ) const;

        //- Find nearest point on cylinder.
        pointIndexHit findNearest
        (
            const point& sample,
            const scalar nearestDistSqr
        ) const;

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

        //- Start point of cylinder
        const point start() const
        {
            return point1_;
        }

        //- End point of cylinder
        const point end() const
        {
            return point2_;
        }

        //- Radius of cylinder
        scalar radius() const
        {
            return radius_;
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
