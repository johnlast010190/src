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
    (c) 2010 Esi Ltd.

Application
    Foam::addHexMeshLayer

Description
    Add layers of cells

\*---------------------------------------------------------------------------*/

#ifndef addHexMeshLayer_H
#define addHexMeshLayer_H

#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "containers/Lists/PackedList/PackedBoolList.H"
#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyTopoChange;
class mapPolyMesh;
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                           Class addHexMeshLayer Declaration
\*---------------------------------------------------------------------------*/

class addHexMeshLayer
{
        //- From surface region to patch
        const labelList& globalToMasterPatch_;

        //- refinement parameters
        const refinementParameters& refineParams_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        fvMeshDistribute& distributor_;

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Primitive patch containg faces to grow layers
        indirectPrimitivePatch& pp_;

        //- How to cut edges
        const scalar cutRatio_;

        //- whether to perform geometric tests on cuts
        const bool geometricChecks_;

        //-List of all grown patches
        const labelList grownPatches_;

        //- Keep track of layer cells and stacks
        const bool trackLayerCells_;

    // Private Member Functions

public:

    //- Runtime type information
    ClassName("addHexMeshLayer");

    // Constructors

        //- Construct from mesh.
        addHexMeshLayer
        (
            const labelList& globalToMasterPatch,
            const refinementParameters& refineParams,
            decompositionMethod& decomposer,
            fvMeshDistribute& distributor,
            meshRefinement& meshRefiner,
            indirectPrimitivePatch& pp,
            const scalar cutRatio,
            const bool geometricTest,
            const labelList grownPatchIDs = labelList(0),
            const bool trackLayerCells = false
        );

    // Member Functions

        // Helper

           void calculateBoundary
           (
               boolList& boundaryPts,
               boolList& boundaryEdges,
               boolList& boundaryFaces
           );

           bool checkLoop
           (
               const fvMesh& mesh,
               const List<labelList>& cutLoops,
               const labelList& edgeLoop,
               const boolList& boundaryPts,
               const label cellI
           );

           labelList createCutLoop
           (
               const fvMesh& mesh,
               const boolList& cutPts,
               const Map<label>& cellPts,
               const label startPt,
               const label startEdge,
               const label startFace,
               const labelHashSet& cellFaces,
               boolList& visited
           );

           labelList calculateCutEdges();


        // Access

            void splitProblemCells();

            void removeProblemCells();

            void addSplitEdgePoints();

            void removeUnusedPoints();

            void calculateSplitLoops
            (
                DynamicList<DynamicList<label>>& allCutCells,
                DynamicList<DynamicList<labelList>>& allCellLoops
            );

            void splitProcessorFaces
            (
                const DynamicList<label>& cutCells,
                const DynamicList<labelList>& cellLoops,
                const List<scalarField>& cellWeights
            );

            void createSplits
            (
                DynamicList<DynamicList<label>>& allCutCells,
                DynamicList<DynamicList<labelList>>& allCellLoops
            );

            void balance();

            void setRefinement();

            void checkForGapSideWalls();

        // Update after topo changes

            void updateMesh
            (
                const mapPolyMesh&
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
