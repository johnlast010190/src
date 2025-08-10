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
    (c) 2015 ICON CFD
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::snappyRefineDriver

Description

SourceFiles
    snappyRefineDriver.C

\*---------------------------------------------------------------------------*/

#ifndef snappyRefineDriver_H
#define snappyRefineDriver_H

#ifdef FOAM_USE_OPENVDB
#include "foamVDB.H"
#endif

#include "snappyHexMeshDriver/refinementParameters/wordPairHashTable.H"
#include "primitives/ints/lists/labelList.H"
#include "meshes/meshShapes/edge/EdgeMap.H"
#include "primitives/Pair/labelPair.H"
#include "edgeMesh/featureEdgeMesh/featureEdgeMesh.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "meshControl/meshControl.H"
#include "sampledSetWriters/writer.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class refinementParameters;
class snapParameters;

class meshRefinement;
class decompositionMethod;
class fvMeshDistribute;

/*---------------------------------------------------------------------------*\
                           Class snappyRefineDriver Declaration
\*---------------------------------------------------------------------------*/

class snappyRefineDriver
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        fvMeshDistribute& distributor_;

        //- From surface region to patch
        const labelList globalToMasterPatch_;

        //- From surface region to patch
        const labelList globalToSlavePatch_;

        //- Reference to mesh process controller
        const meshControl& controller_;

        //- How to write lines. Used e.g. when writing leak-paths
        const writer<scalar>& setFormatter_;

    // Private Member Functions

        void addEdgeToRegion
        (
            const edge& e,
            const label region,
            const label maxLevel,
            EdgeMap<labelPair>& levelPerEdge
        ) const;

        labelList consistentRefinement
        (
            const label nBuffer,
            const labelList& candidateCells
        ) const;

        //- Refine all cells where refinement surface bounding box contained
        // in cell is smaller than the cell size
        label initialRefine
        (
            const refinementParameters& refineParams,
            const label maxIter,
            const label minRefine
        );

        //- Refine all cells pierced by explicit feature edges
        label featureEdgeRefine
        (
            const refinementParameters& refineParams,
            const label maxIter,
            const label minRefine
        );

        //- Refine gaps where only one cell spans the gap
        void refineGaps
        (
            const refinementParameters& refineParams
        );

        //- Refinecells at removal interface
        void removalInterfaceRefine
        (
            const refinementParameters& refineParams
        );

        //- Refine baffle edges in cell level change around edge
        void refineAtBaffleEdges
        (
            const refinementParameters& refineParams
        );

        //- Refine all cells containing small surface features
        label smallFeatureRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Refine all cells interacting with the surface
        label surfaceOnlyRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Refine based on proximity of surfaces
        label proximityRefine
        (
            const refinementParameters& refineParams,
            const label maxIter,
            const label minRefine,
            const bool oneMachine = false
        );

         //- Refine all cells in small gaps
        label gapOnlyRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );


        //- Refine all cells in large gaps
        label bigGapOnlyRefine
        (
            const refinementParameters& refineParams,
            const bool spreadGapSize,
            const label maxIter
        );

        //- Refine cells with almost all sides refined
        label danglingCellRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Maintain single level jump across edges
        label singleLevelEdgeRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Refine cells with opposite faces with differing refinement level
        label refinementInterfaceRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Remove problem cells when mesh dualised
        void removeDualProblemCells
        (
            const refinementParameters& refineParams
        );

        //- Refine interface cells that might have problems adding layers to
        label interfaceProblemCellRefine
        (
            const refinementParameters& refineParams,
            const dictionary& layerDict,
            const label maxIter
        );

        //- Remove all cells within intersected region (including wrapping)
        void removeWrappedCells
        (
            const refinementParameters& refineParams,
            const dictionary& dict,
            const bool handleSnapProblems,
            const dictionary& motionDict
        );

        //- Remove all cells within intersected region
        void removeInsideCells
        (
            const refinementParameters& refineParams,
            const label nBufferLayers
        );

        //- Remove all cells within intersected region
        void removeOutsideCells
        (
            const refinementParameters& refineParams
        );

        //- Remove all cells inside/outside shell
        label shellRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        label refineBoundaryCells
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Add baffles and remove unreachable cells
        void baffleAndSplitMesh
        (
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool handleSnapProblems,
            const dictionary& motionDict,
            const bool threaded = false
        );

        void removeBoundaryZoneProblemCells
        (
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const dictionary& motionDict,
            List<labelPair>& baffles
        );

        //- Add zones
        void zonify
        (
            const refinementParameters& refineParams,
            wordPairHashTable& zonesToFaceZone
        );

        void splitAndMergeBaffles
        (
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool handleSnapProblems,
            const dictionary& motionDict
        );

        //- Merge refined boundary faces (from exposing coarser cell)
        void mergePatchFaces
        (
            const bool geometricMerge,
            const dictionary& motionDict
        );

        //- Volumetric anisotropic refinement
        label anisoRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- reset boundary faces alloacted in a faceZone
        void resetBoundaryFaceZones();

        //- Disallow default bitwise copy construct
        snappyRefineDriver(const snappyRefineDriver&);

        //- Disallow default bitwise assignment
        void operator=(const snappyRefineDriver&);


public:

    //- Runtime type information
    ClassName("snappyRefineDriver");


    // Constructors

        //- Construct from components
        snappyRefineDriver
        (
            meshRefinement& meshRefiner,
            decompositionMethod& decomposer,
            fvMeshDistribute& distributor,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch,
            const meshControl& controller,
            const writer<scalar>& setFormatter
        );


    // Member Functions

        //- Do all the refinement
        void doRefine
        (
            const dictionary& meshDict,
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool prepareForSnapping
        );

        //- Do all the refinement (for shell meshing)
        void doThinShellRefine
        (
            const dictionary& meshDict,
            const refinementParameters& refineParams
        );

#ifdef FOAM_USE_OPENVDB
        //- Do all the refinement using OpenVDB
        void doVDBRefine
        (
            const dictionary& meshDict,
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool prepareForSnapping
        );

        //- create zones
        void doVDBzones
        (
            const dictionary& meshDict,
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool prepareForSnapping
        );
#endif

        //- Helper: add faceZones and patches
        static void addFaceZones
        (
            meshRefinement& meshRefiner,
            const refinementParameters& refineParams,
            const HashTable<Pair<word>>& faceZoneToPatches
        );

        //- Add cell and face zones to dual mesh
        void zonifyDual
        (
            const refinementParameters& refineParams,
            const bool updateIntersections
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
