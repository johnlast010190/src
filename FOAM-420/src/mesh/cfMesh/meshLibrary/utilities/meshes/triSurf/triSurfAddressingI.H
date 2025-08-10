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

Description


\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "utilities/meshes/triSurf/triSurfAddressing.H"

# ifdef USE_OMP
#include <omp.h>
# endif

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const VRWGraph& triSurfAddressing::pointFacets() const
{
    if (!pointFacetsPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate pointFacets" << abort(FatalError);
        # endif

        calculatePointFacets();
    }

    return *pointFacetsPtr_;
}

inline const edgeLongList& triSurfAddressing::edges() const
{
    if (!edgesPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate edges" << abort(FatalError);
        # endif

        calculateEdges();
    }

    return *edgesPtr_;
}

inline const VRWGraph& triSurfAddressing::facetEdges() const
{
    if (!facetEdgesPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate facetEdges" << abort(FatalError);
        # endif

        calculateFacetEdges();
    }

    return *facetEdgesPtr_;
}

inline const VRWGraph& triSurfAddressing::edgeFacets() const
{
    if (!edgeFacetsPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate edgeFacets" << abort(FatalError);
        # endif

        calculateEdgeFacets();
    }

    return *edgeFacetsPtr_;
}

inline const VRWGraph& triSurfAddressing::pointEdges() const
{
    if (!pointEdgesPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate pointEdges" << abort(FatalError);
        # endif

        calculatePointEdges();
    }

    return *pointEdgesPtr_;
}

inline const VRWGraph& triSurfAddressing::facetFacets() const
{
    if (!facetFacetsEdgesPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate facetFacets" << abort(FatalError);
        # endif

        calculateFacetFacetsEdges();
    }

    return *facetFacetsEdgesPtr_;
}

inline const vectorField& triSurfAddressing::pointNormals() const
{
    if (!pointNormalsPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate pointNormals" << abort(FatalError);
        # endif

        calculatePointNormals();
    }

    return *pointNormalsPtr_;
}

inline const vectorField& triSurfAddressing::facetNormals() const
{
    if (!facetNormalsPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate facetNormals" << abort(FatalError);
        # endif

        calculateFacetNormals();
    }

    return *facetNormalsPtr_;
}

inline const vectorField& triSurfAddressing::facetCentres() const
{
    if (!facetCentresPtr_)
    {
        # ifdef USE_OMP
        if (omp_in_parallel())
            FatalErrorInFunction
                << "Cannot calculate facetCentres" << abort(FatalError);
        # endif

        calculateFacetCentres();
    }

    return *facetCentresPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
