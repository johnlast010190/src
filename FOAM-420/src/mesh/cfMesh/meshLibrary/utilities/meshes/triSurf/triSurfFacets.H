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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Class
    triSurfFacets

Description
    Facets for the triangulated surface

SourceFiles
    triSurfFacets.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfFacets_H
#define triSurfFacets_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/meshes/triSurf/triSurfAddressing.H"
#include "utilities/meshes/primitives/meshSubsets/meshSubset.H"
#include "triSurface/patches/geometricSurfacePatchList.H"
#include "utilities/containers/LongList/LongList.H"
#include "meshes/meshShapes/labelledTri/labelledTri.H"
#include "utilities/containers/DynList/DynList.H"
#include "containers/HashTables/Map/Map.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfFacets
{
protected:

    // Protected data
        //- list of triangles
        LongList<labelledTri> triangles_;

        //- list of boundary patches and their properties
        geometricSurfacePatchList patches_;

        //- map of point subsets
        Map<meshSubset> facetSubsets_;

    // Disallow bitwise assignment
        void operator=(const triSurfFacets&);

        triSurfFacets(const triSurfFacets&);

public:

    // Constructors
        //- Null constructor
        triSurfFacets();

        //- Construct from components without the boundary
        triSurfFacets(const LongList<labelledTri>& triangles);

        //- Construct from components
        triSurfFacets
        (
            const LongList<labelledTri>& triangles,
            const geometricSurfacePatchList& patches
        );

    // Destructor
        ~triSurfFacets();

    // Member functions
        //- return the number of triangles
        inline label size() const;

        //- access to facets
        inline const LongList<labelledTri>& facets() const;

        //- access to patches
        inline const geometricSurfacePatchList& patches() const;

        //- return list of patches in the boundary
        wordList patchNames() const;

        //- return a list of patch indices corresponding to the given
        // name, expanding regular expressions
        labelList findPatches(const word& patchName) const;

        //- append a triangle to the end of the list
        inline void appendTriangle(const labelledTri& tria);

        //- point subsets
        label addFacetSubset(const word&);
        void removeFacetSubset(const label);
        word facetSubsetName(const label) const;
        label facetSubsetIndex(const word&) const;
        inline void addFacetToSubset(const label, const label);
        inline void removeFacetFromSubset(const label, const label);
        inline void facetInSubsets(const label, DynList<label>&) const;
        inline void facetSubsetIndices(DynList<label>&) const;
        template<class ListType>
        inline void facetsInSubset(const label, ListType&) const;
        template<class ListType>
        inline void updateFacetsSubsets(const ListType&);

    // Operators

        //- access to a triangle
        inline const labelledTri& operator[](const label) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/meshes/triSurf/triSurfFacetsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
