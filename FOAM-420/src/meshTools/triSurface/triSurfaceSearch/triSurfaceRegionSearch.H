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
    (c) 2022 Esi Ltd.

Class
    Foam::triSurfaceRegionSearch

Description
    Helper class to search on triSurface. Creates an octree for each region of
    the surface and only searches on the specified regions.

SourceFiles
    triSurfaceRegionSearch.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceRegionSearch_H
#define triSurfaceRegionSearch_H

#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/primitiveShapes/objectHit/pointIndexHit.H"
#include "meshes/meshShapes/labelledTri/labelledTri.H"
#include "containers/Lists/IndirectList/IndirectList.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "algorithms/indexedOctree/indexedOctree.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class triSurfaceRegionSearch Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceRegionSearch
:
    public triSurfaceSearch
{
    // Private typedefs

        typedef PrimitivePatch<IndirectList<labelledTri>, const pointField&>
            indirectTriSurface;

        typedef treeDataPrimitivePatch<indirectTriSurface>
            treeDataIndirectTriSurface;

        typedef indexedOctree<treeDataIndirectTriSurface> treeType;


    // Private data

        //- Surface is split into patches by region
        mutable PtrList<indirectTriSurface> indirectRegionPatches_;

        //- Search tree for each region
        mutable PtrList<treeType> treeByRegion_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        triSurfaceRegionSearch(const triSurfaceRegionSearch&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceRegionSearch&);


public:

    // Constructors

        //- Construct from surface. Holds reference to surface!
        explicit triSurfaceRegionSearch(const triSurface&);

        //- Construct from surface and dictionary. Holds reference to surface!
        triSurfaceRegionSearch(const triSurface&, const dictionary& dict);


    //- Destructor
    ~triSurfaceRegionSearch();

        //- Clear storage
        void clearOut();


    // Member Functions

        // Access

            //- Demand driven construction of octree for each region.
            //  Currently creates a tree for each region; could optimise
            //  by only constructing trees when they are in regionIndices
            const PtrList<treeType>& treeByRegion() const;

        // Query

            //- Find the nearest point on the surface out of the regions
            //  supplied in the list regionIndices. Ignores regions that are
            //  not specified
            void findNearest
            (
                const pointField& samples,
                const scalarField& nearestDistSqr,
                const labelList& regionIndices,
                List<pointIndexHit>& info
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
