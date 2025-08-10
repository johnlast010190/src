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
    meshSurfaceEdgeExtractor

Description
    Extracts edges and corners to cornform the mesh to the geometry

SourceFiles
    meshSurfaceEdgeExtractor.C

\*---------------------------------------------------------------------------*/

#ifndef meshSurfaceEdgeExtractor_H
#define meshSurfaceEdgeExtractor_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "primitives/bools/lists/boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctree;
class meshSurfaceEngine;

/*---------------------------------------------------------------------------*\
                Class meshSurfaceEdgeExtractor Declaration
\*---------------------------------------------------------------------------*/

class meshSurfaceEdgeExtractor
{
    // Private data
        //- mesh
        polyMeshGen& mesh_;

        label nPoints_;

        boolList boundaryCell_;

        List<direction> nFacesInCell_;

        //- reference to the octree
        const meshOctree& meshOctree_;

        //- regions for boundary vertices
        VRWGraph pointRegions_;

    // Private member functions
        //- remove old boundary faces
        void removeOldBoundaryFaces();

        //- create vertices on surface edges
        void createEdgeVertices();

        //- create new boundary faces
        void createBoundaryFaces();

        //- Disallow default bitwise copy construct
        meshSurfaceEdgeExtractor(const meshSurfaceEdgeExtractor&);

        //- Disallow default bitwise assignment
        void operator=(const meshSurfaceEdgeExtractor&);

public:

    // Constructors

        //- Construct from mesh data
        meshSurfaceEdgeExtractor
        (
            polyMeshGen& mesh,
            const meshOctree& octree,
            const labelList& pointRegion
        );

    // Destructor

        ~meshSurfaceEdgeExtractor();

    // Member Functions
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
