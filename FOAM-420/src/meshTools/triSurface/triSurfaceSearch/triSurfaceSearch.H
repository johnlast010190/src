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
    (c) 2017 OpenCFD Ltd.
    (c) 2021 Esi Ltd.

Class
    Foam::triSurfaceSearch

Description
    Helper class to search on triSurface.

SourceFiles
    triSurfaceSearch.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceSearch_H
#define triSurfaceSearch_H

#include "meshes/primitiveShapes/point/pointField.H"
#include "primitives/bools/lists/boolList.H"
#include "meshes/primitiveShapes/objectHit/pointIndexHit.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataTriSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurface;

/*---------------------------------------------------------------------------*\
                      Class triSurfaceSearch Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceSearch
{
    // Private data

        //- Reference to surface to work on
        const triSurface& surface_;

        //- Optional tolerance to use in searches
        scalar tolerance_;

        //- Optional max tree depth of octree
        label maxTreeDepth_;

        //- Octree for searches
        mutable autoPtr<indexedOctree<treeDataTriSurface>> treePtr_;

        //- Check whether the current hit on the surface which lies on lineVec
        //  is unique.
        //  point : return 1 unique hit
        //  edge  : return 1 hit if in the cone of the edge faces
        //        : return 2 hits if outside or on the cone.
        bool checkUniqueHit
        (
            const pointIndexHit& currHit,
            const DynamicList<pointIndexHit, 1, 1>& hits,
            const vector& lineVec
        ) const;

        //- Disallow default bitwise copy construct
        triSurfaceSearch(const triSurfaceSearch&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceSearch&) = delete;


public:

    // Constructors

        //- Construct from surface. Holds reference to surface!
        explicit triSurfaceSearch(const triSurface& surface);

        //- Construct from surface and dictionary.
        triSurfaceSearch(const triSurface& surface, const dictionary& dict);

        //- Construct from components.
        //  A invalid (negative) tolerance uses the default tolerance.
        triSurfaceSearch
        (
            const triSurface& surface,
            const scalar tolerance,
            const label maxTreeDepth = 10
        );


    //- Destructor
    ~triSurfaceSearch();

        //- Clear storage
        void clearOut();


    // Member Functions

        //- Demand driven construction of the octree
        const indexedOctree<treeDataTriSurface>& tree() const;

        //- Return reference to the surface.
        const triSurface& surface() const
        {
            return surface_;
        }

        //- Return tolerance to use in searches
        scalar tolerance() const
        {
            return tolerance_;
        }

        //- Return max tree depth of octree
        label maxTreeDepth() const
        {
            return maxTreeDepth_;
        }

        //- Calculate for each searchPoint inside/outside status.
        boolList calcInside(const pointField& searchPoints) const;

        void findNearest
        (
            const pointField& samples,
            const scalarField& nearestDistSqr,
            List<pointIndexHit>& info,
            const bool threaded = false,
            const scalar isoValue = 0
        ) const;

        //- Calculate nearest point on surface for single searchPoint. Returns
        //  in pointIndexHit:
        //  - hit()      : whether nearest point found within bounding box
        //  - hitPoint() : coordinate of nearest point
        //  - index()    : surface triangle label
        pointIndexHit nearest(const point& pt, const vector& span) const;

        void findLine
        (
            const pointField& start,
            const pointField& end,
            List<pointIndexHit>& info
        ) const;

        void findLineAny
        (
            const pointField& start,
            const pointField& end,
            List<pointIndexHit>& info,
            const bool threaded = false
        ) const;

        //- Calculate all intersections from start to end
        void findLineAll
        (
            const pointField& start,
            const pointField& end,
            List<List<pointIndexHit>>& info
        ) const;

        //- Find all elements intersecting sphere
        void findSphere
        (
            const pointField& centres,
            const scalarField& radiusSqr,
            List<labelList>& hitIndexes,
            labelList& nOverlapChecks
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
