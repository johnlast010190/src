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
    tetMeshExtractorOctree

Description
    A class which extracts tet mesh out of an octree structure

SourceFiles
    tetMeshExtractorOctree.C

\*---------------------------------------------------------------------------*/

#ifndef tetMeshExtractorOctree_H
#define tetMeshExtractorOctree_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/meshes/primitives/partTet/partTet.H"
#include "utilities/tetrahedra/tetCreatorOctree/tetCreatorOctree.H"
#include "utilities/containers/VRWGraph/VRWGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class tetMeshExtractorOctree Declaration
\*---------------------------------------------------------------------------*/

class tetMeshExtractorOctree
{
    // Private data
    //- create tets
    tetCreatorOctree tetCreator_;

    //- reference to the mesh
    polyMeshGen& mesh_;

    // Private member functions
    //- copy tetPoints_ into mesh
    void createPoints();

    //- create mesh data
    void createPolyMesh();

    // Private copy constructor
    //- Disallow default bitwise copy construct
    tetMeshExtractorOctree ( const tetMeshExtractorOctree& );

    //- Disallow default bitwise assignment
    void operator= ( const tetMeshExtractorOctree& );

public:

    // Constructors

    //- Construct from octree and mesh data
    tetMeshExtractorOctree
    (
        const meshOctree& octree,
        const IOdictionary& meshDict,
        polyMeshGen& mesh
    );

    // Destructor

    ~tetMeshExtractorOctree();


    // Member Functions
    void createMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
