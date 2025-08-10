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
    meshOctreeCreator

Description
    Creates octree for mesh generation

SourceFiles
    meshOctreeCreator.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeCreator_H
#define meshOctreeCreator_H

#include "primitives/bools/lists/boolList.H"
#include "utilities/containers/DynList/DynList.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "utilities/octrees/meshOctree/meshOctreeModifier/meshOctreeModifier.H"
#include "utilities/octrees/meshOctree/refinementControls/patchRefinement/patchRefinementList.H"

#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class meshOctreeCube;

/*---------------------------------------------------------------------------*\
                    Class meshOctreeCreator Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeCreator
{
protected:

    // Protected data
        //- Reference to meshOctree
        meshOctree& octree_;

        //- Scaling factor
        scalar scalingFactor_;

        //- Dictionary containing information necessary to perform refinement
        const IOdictionary* meshDictPtr_;

        //- hex refinement flag
        bool hexRefinement_;

private:
    // Private member functions

    // Octree refinement
        //- refine boxes contained inside the objects for refinement
        void refineBoxesContainedInObjects();

        //- refine boxes intersected by surface meshes
        //- used as refinement sources
        void refineBoxesIntersectingSurfaces();

        //- refine boxes by edge meshes
        //- used as refinement sources
        void refineBoxesIntersectingEdgeMeshes();

        //- refine boxes near DATA boxes to get a nice smooth surface
        void refineBoxesNearDataBoxes(const label nLayers = 1);

        //- refine boxes of the given flag to the given size
        void refineBoxes(const direction refLevel, const direction cubeType);

        // Private octree creation methods
        void refineBoundary();
        void createInsideOutsideInformation();
        void refineInsideAndUnknownBoxes();

    // Private member functions for parallel runs
        //- move octree cubes from one processor to another until
        //- each processor contains the same number of leaves
        //- of the required type
        void loadDistribution(const bool distributeUsed = false);

        // information about octree refinement
        //- ref level to achieve max cell size
        direction globalRefLevel_;

        //- this list contains ref level for each surface triangle
        List<DynList<std::pair<direction, scalar>>> surfRefLevel_;

        //- set the boundBox such that maxCellSize is achieved
        void setRootCubeSizeAndRefParameters();

    // Private copy constructor
        //- Disallow default bitwise copy construct
        meshOctreeCreator(const meshOctreeCreator&);

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeCreator&);

public:

    // Constructors

        //- construct from meshOctree
        meshOctreeCreator(meshOctree&);

        //- Construct from meshOctree and dictionary
        meshOctreeCreator(meshOctree& mo, const IOdictionary& dict);

    // Destructor

        ~meshOctreeCreator();


    // Member Functions
        //- set the scaling factor
        void setScalingFactor(const scalar);

        //- activate octree refinement which marks all sons of an octree box
        //- for refinement in case a single son is marked for refinement
        //- this type of refinement is necessary for creating hex meshes
        void activateHexRefinement();

        //- create octree boxes
        void createOctreeBoxes();

        //- refine boxes containing surface elements
        //- this is used for proximity searches
        void createOctreeWithRefinedBoundary
        (
            const direction maxLevel,
            const label nTrianglesInLeaf = 15
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
