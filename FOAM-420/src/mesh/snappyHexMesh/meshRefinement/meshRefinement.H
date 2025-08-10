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
    (c) 2010-2021 Esi Ltd.
    (c) 2010-2012 ICON CFD
    (c) 2015-2023 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

Class
    Foam::meshRefinement

Description
    Helper class which maintains intersections of (changing) mesh with
    (static) surfaces.

    Maintains
    - per face any intersections of the cc-cc segment with any of the surfaces

SourceFiles
    meshRefinement.C
    meshRefinementBaffles.C
    meshRefinementGapRefine.C
    meshRefinementMerge.C
    meshRefinementProblemCells.C
    meshRefinementRefine.C
    meshRefinementWrapper.C
    meshRefinementSplit.C

\*---------------------------------------------------------------------------*/

#ifndef meshRefinement_H
#define meshRefinement_H

#include "polyTopoChange/hexRef8/hexRef8.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "memory/autoPtr/autoPtr.H"
#include "primitives/Pair/labelPairHashes.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "fields/GeometricFields/pointFields/pointFieldsFwd.H"
#include "primitives/Tuple2/Tuple2.H"
#include "meshes/primitiveShapes/objectHit/pointIndexHit.H"
#include "snappyHexMeshDriver/refinementParameters/wordPairHashTable.H"
#include "refinementSurfaces/surfaceZonesInfo.H"
#include "algorithms/indexedOctree/volumeType.H"
#include "fields/Fields/DynamicField/DynamicField.H"
#include "edgeMesh/featureEdgeMesh/featureEdgeMesh.H"
#include "cfdTools/general/include/fvCFD.H"
#include "snappyHexMeshDriver/refinementParameters/refinementParameters.H"
#include "containers/Lists/CompactListList/CompactListList.H"
#include "cellFeatures/cellFeatures.H"
#include "meshControl/meshControl.H"
#include "sampledSetWriters/writer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class fvMesh;
class mapDistributePolyMesh;
class decompositionMethod;
class refinementSurfaces;
class refinementFeatures;
class shellSurfaces;
class removeCells;
class fvMeshDistribute;
class removePoints;
class localPointRegion;
class snapParameters;
class regionSplit;

/*---------------------------------------------------------------------------*\
                       Class meshRefinement Declaration
\*---------------------------------------------------------------------------*/

class meshRefinement
{

public:

    // Public data types
        //- A HashSet with labelPair keys.
        typedef HashSet<labelPair, Hash<labelPair>> labelPairHashSet;

        //- Enumeration for what to debug
        enum IOdebugType
        {
            IOMESH,
            IOOBJINTERSECTIONS,
            IOFEATURESEEDS,
            IOATTRACTION,
            IOLAYERINFO
        };

        static const NamedEnum<IOdebugType, 5> IOdebugTypeNames;
        enum debugType
        {
            MESH = 1<<IOMESH,
            OBJINTERSECTIONS = 1<<IOOBJINTERSECTIONS,
            FEATURESEEDS = 1<<IOFEATURESEEDS,
            ATTRACTION = 1<< IOATTRACTION,
            LAYERINFO = 1<<IOLAYERINFO
        };

        //- Enumeration for what to output
        enum IOoutputType
        {
            IOOUTPUTLAYERINFO
        };

        static const NamedEnum<IOoutputType, 1> IOoutputTypeNames;
        enum outputType
        {
            OUTPUTLAYERINFO = 1<<IOOUTPUTLAYERINFO
        };

        //- Enumeration for what to write
        enum IOwriteType
        {
            IOWRITEMESH,
            IONOWRITEREFINEMENT,
            IOWRITELEVELS
        };

        static const NamedEnum<IOwriteType, 3> IOwriteTypeNames;
        enum writeType
        {
            WRITEMESH = 1<<IOWRITEMESH,
            NOWRITEREFINEMENT = 1<<IONOWRITEREFINEMENT,
            WRITELEVELS = 1<<IOWRITELEVELS
        };

        //- Enumeration for how the userdata is to be mapped upon refinement.
        enum mapType
        {
            MASTERONLY = 1, //!< maintain master only
            KEEPALL = 2,    //!< have slaves (upon refinement) from master
            REMOVE = 4      //!< set value to -1 any face that was refined
        };

        class maxMagEqOp
        {

        public:

            void operator()(vector& x, const vector& y) const
            {
                for (direction i = 0; i < vector::nComponents; i++)
                {
                    scalar magX = mag(x[i]);
                    scalar magY = mag(y[i]);

                    if (magX < magY)
                    {
                        x[i] = y[i];
                    }
                    else if (magX == magY)
                    {
                        if (y[i] > x[i])
                        {
                            x[i] = y[i];
                        }
                    }
                }
            }
        };

        //- Combine operator class to combine normal with other normal.
        class pointFieldCombine
        {
        public:

            void operator()(pointField& x, const pointField& y) const
            {
                if (y.size())
                {
                    if (x.size() == 0)
                    {
                        x = y;
                    }
                    else
                    {
                        x.append(y);
                    }
                 }
            }
        };

private:

    // Static data members

        //- Control of writing level
        static writeType writeLevel_;

        //- Control of output/log level
        static outputType outputLevel_;

    // Private class

        //- Combine operator class for equalizing normals.
        class minMagEqOp
        {
            public:

            void operator()(vector& x, const vector& y) const
            {
                if (magSqr(y) < magSqr(x))
                {
                    x = y;
                }
            }
        };

        class xorEqOp
        {
            public:

            void operator()(bool& x, const bool& y) const
            {
                if (!x && y)
                {
                    x = true;
                }
                else if (x && !y)
                {
                    x = true;
                }
                else
                {
                    x = false;
                }
            }
        };

    // Private data

        //- Reference to mesh
        fvMesh& mesh_;

        //- cell region ID (primarily used for dual mesher)
        labelIOList cellRegionID_;

        //- Tolerance used for sorting coordinates (used in 'less' routine)
        const scalar mergeDistance_;

        //- Overwrite the mesh?
        const bool overwrite_;

        //- Control of the overall mesh generation process
        const meshControl& controller_;

        //- perform crack detection?
        const  Tuple2<bool, scalar> crackDetection_;

        //- Use additional rays in crack test?
        bool addedRays_;

        //- Instance of mesh upon construction. Used when in overwrite_ mode.
        const word oldInstance_;

        //- All surface-intersection interaction
        const refinementSurfaces& surfaces_;

        //- All shell-refinement interaction
        const shellSurfaces& shells_;

        //- All limit-refinement interaction
        const shellSurfaces& limitShells_;

        //- All feature-edge interaction
        autoPtr<refinementFeatures>& featuresPtr_;

        //- Refinement engine
        hexRef8 meshCutter_;

        //- Per cc-cc vector the index of the surface hit
        labelIOList surfaceIndex_;

        //- Per face if face needs to be wrapped
        labelIOList wrapIndex_;

        // For baffle merging

            //- Original patch for baffle faces that used to be on
            //  coupled patches
            Map<label> faceToCoupledPatch_;


        //- User supplied face based data.
        List<Tuple2<mapType, labelList>> userFaceData_;

        //- Meshed patches - are treated differently. Stored as wordList since
        //  order changes.
        wordList meshedPatches_;

        //- Estimate of crack tolerance
        scalar crackTol_;

         //- store oldPoints field (during cell splitting)
        pointField* oldPointsPtr_;

        //- Control to activate cell level wrapping
        bool wrapActive_;

        //- Control to activate gap field creation
        bool createGapField_;

        //- FaceZone to master patch name
        HashTable<word, word> faceZoneToMasterPatch_;

        //- FaceZone to slave patch name
        HashTable<word, word> faceZoneToSlavePatch_;

        //- FaceZone to method to handle faces
        HashTable<surfaceZonesInfo::faceZoneType, word> faceZoneToType_;


    // Private Member Functions

        //- Set whether to update oldPoints field
        void setOldPoints(const pointField& oldPoints)
        {
            if (!oldPointsPtr_)
            {
                oldPointsPtr_ = new pointField(oldPoints);
            }
        };

        //- Set crack tolerance based on minimum edge length
        void setCrackTol();

        //- Add patchfield of given type to all fields on mesh
        template<class GeoField>
        static void addPatchFields(fvMesh&, const word& patchFieldType);

        //- Reorder patchfields of all fields on mesh
        template<class GeoField>
        static void reorderPatchFields(fvMesh&, const labelList& oldToNew);

        //- Find out which faces have changed given cells (old mesh labels)
        //  that were marked for refinement.
        static labelList getChangedFaces
        (
            const mapPolyMesh&,
            const labelList& oldCellsToRefine
        );

        //- Calculate coupled boundary end vector and refinement level
        void calcNeighbourData
        (
            const pointField& cellCentres,
            labelList& neiLevel,
            pointField& neiCc,
            const bool threaded = false
        ) const;

        //- Calculate coupled boundary end vector and refinement level
        // with added perturbation vector
        void calcPerturbedNeighbourData
        (
            const pointField& cellCentres,
            const pointField& perturbCellCentres,
            labelList& neiLevel,
            pointField& neiCc
        ) const;

        //- Calculate rays from cell-centre to cell-centre and corresponding
        //  min cell refinement level
        void calcCellCellRays
        (
            const pointField& cellCentres,
            const pointField& neiCc,
            const labelList& neiLevel,
            const labelList& testFaces,
            pointField& start,
            pointField& end,
            labelList& minLevel,
            const bool threaded = false
        ) const;

        //- Find any intersection of surface (using blocks).
        // Store in surfaceIndex_.
        void updateIntersectionsByBlock
        (
            const labelList& changedFaces,
            const label numberBlocks
        );

        //- Get cells which are inside any closed surface. Note that
        //  all closed surfaces
        //  will have already been oriented to have keepPoint outside.
        labelList getInsideCells(const word&) const;

        //- Do all to remove inside cells
        autoPtr<mapPolyMesh> removeInsideCells
        (
            const string& msg,
            const label exposedPatchi
        );

        // Refinement candidate selection

            //- Mark cell for refinement (if not already marked). Return false
            // if refinelimit hit. Keeps running count (in nRefine) of cells
            // marked for refinement
            static bool markForRefine
            (
                const label markValue,
                const label nAllowRefine,
                label& cellValue,
                label& nRefine
            );

            //- Calculate list of cells to refine based on surface
            // bounding box size with respect to local cell size.
            label markInitialRefinement
            (
                const label nAllowRefine,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark every cell with level of feature passing through it
            //  (or -1 if not passed through). Uses tracking.
            void markFeatureCellLevel
            (
                const refinementParameters& refineParams,
                labelList& maxFeatureLevel
            ) const;

            //- Calculate list of cells to refine based on intersection of
            //  features.
            label markFeatureRefinement
            (
                const refinementParameters& refineParams,
                const label nAllowRefine,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for distance-to-feature based refinement.
            label markInternalDistanceToFeatureRefinement
            (
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for refinement-shells based refinement.
            label markInternalRefinement
            (
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Unmark cells for refinement based on limit-shells. Return number
            //  of limited cells.
            label unmarkInternalRefinement
            (
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Collect faces that are intersected and whose neighbours aren't
            //  yet marked  for refinement.
            labelList getRefineCandidateFaces
            (
                const labelList& refineCell
            ) const;

            //- Mark cells for surface intersection based refinement.
            label markSurfaceRefinement
            (
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells intersected by the surface if they are inside
            //  close gaps
            label markSurfaceGapRefinement
            (
                const scalar planarCos,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Generate single ray from nearPoint in direction of nearNormal
            label generateRays
            (
                const point& nearPoint,
                const vector& nearNormal,
                const FixedList<label, 3>& gapInfo,
                const volumeType& mode,

                const label cLevel,

                DynamicField<point>& start,
                DynamicField<point>& end
            ) const;

            //- Generate pairs of rays through cell centre
            //  Each ray pair has start, end, and expected gap size
            label generateRays
            (
                const bool useSurfaceNormal,

                const point& nearPoint,
                const vector& nearNormal,
                const FixedList<label, 3>& gapInfo,
                const volumeType& mode,

                const point& cc,
                const label cLevel,

                DynamicField<point>& start,
                DynamicField<point>& end,
                DynamicField<scalar>& gapSize,

                DynamicField<point>& start2,
                DynamicField<point>& end2,
                DynamicField<scalar>& gapSize2
            ) const;

            //- Select candidate cells (cells inside a shell with gapLevel
            //  specified)
            void selectGapCandidates
            (
                const labelList& refineCell,
                const label nRefine,

                labelList& cellMap,
                List<FixedList<label, 3>>& shellGapInfo,
                List<volumeType>& shellGapMode
            ) const;

            //- Merge gap information coming from shell and from surface
            //  (surface wins)
            void mergeGapInfo
            (
                const FixedList<label, 3>& shellGapInfo,
                const volumeType shellGapMode,
                const FixedList<label, 3>& surfGapInfo,
                const volumeType surfGapMode,
                FixedList<label, 3>& gapInfo,
                volumeType& gapMode
            ) const;

            //- Mark cells for non-surface intersection based gap refinement
            label markInternalGapRefinement
            (
                const scalar planarCos,
                const bool spreadGapSize,
                const label nAllowRefine,
                labelList& refineCell,
                label& nRefine,
                labelList& numGapCells,
                scalarField& gapSize
            ) const;

            //- Refine cells containing small gaps
            label markSmallFeatureRefinement
            (
                const scalar planarCos,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Helper: count number of normals1 that are in normals2
            label countMatches
            (
                const List<point>& normals1,
                const List<point>& normals2,
                const scalar tol = 1e-6
            ) const;

            //- Mark cells for surface curvature based refinement. Marks if
            //  local curvature > curvature or if on different regions
            //  (markDifferingRegions)
            label markSurfaceCurvatureRefinement
            (
                const scalar curvature,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,
                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cell if local a gap topology or
            bool checkProximity
            (
                const scalar planarCos,
                const label nAllowRefine,

                const label surfaceLevel,
                const vector& surfaceLocation,
                const vector& surfaceNormal,

                const label celli,

                label& cellMaxLevel,
                vector& cellMaxLocation,
                vector& cellMaxNormal,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for surface proximity based refinement.
            label markProximityRefinement
            (
                const scalar curvature,
                const label nAllowRefine,
                const labelList& neiLevel,
                const pointField& neiCc,

                labelList& refineCell,
                label& nRefine
            ) const;

            //- Mark cells for proximity intersection based refinement.
            label markProximityIncrementRefinement
            (
                const refinementParameters& refineParams,
                const label nAllowRefine,

                label iter,
                labelList& refineCell,
                label& nRefine
            );

            void proximityIntersectionChecks
            (
                const refinementParameters& refineParams,
                const boolList& checkCells,
                const label nRays,
                const label maxGapCells,
                boolList& possibles,
                boolList& refine
             );


        // Baffle handling

            //- Keep boundary point cells if point inside geometry
            void keepBoundaryPointCells
            (
                const labelList& testFaces,
                labelList& ownPatch,
                labelList& nbrPatch
            );

            //- Get faces to repatch. Returns map from face to patch.
            Map<labelPair> getZoneBafflePatches
            (
                const bool allowBoundary,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch
            ) const;

            //- Calculate intersections. Return per face -1 or the global
            //  surface region
            void getIntersections
            (
                const pointField& cellCentres,
                const labelList& surfacesToTest,
                const labelList& testFaces,

                labelList& globalRegion1,
                labelList& globalRegion2,
                const bool checkCracks
            ) const;

            //- Calculate intersections on zoned faces. Return per face -1
            // or the index of the surface and the orientation w.r.t. surface
            void getIntersections
            (
                const pointField& cellCentres,
                const labelList& surfacesToTest,
                const labelList& testFaces,

                labelList& namedSurfaceIndex,
                PackedBoolList& posOrientation
            ) const;

            //- Determine patches for baffles
            void getBafflePatches
            (
                const refinementParameters& refineParams,
                const labelList& globalToMasterPatch,
                labelList& ownPatch,
                labelList& nbrPatch
            ) const;

            //- Determine patch for baffle using some heuristic (and not
            //  surface)
            label getBafflePatch
            (
                const labelList& mPatches,
                const labelList& facePatch,
                const label faceI
            ) const;

            //- Repatches external face or creates baffle for internal face
            //  with user specified patches (might be different for both sides).
            //  Returns label of added face.
            label createBaffle
            (
                const label facei,
                const label ownPatch,
                const label nbrPatch,
                polyTopoChange& meshMod
            ) const;

        // Problem cell handling

            //- Helper function to mark face as being on 'boundary'. Used by
            //  markFacesOnProblemCells
            bool markBoundaryFace
            (
                const label facei,
                boolList& isBoundaryFace,
                boolList& isBoundaryEdge,
                boolList& isBoundaryPoint
            ) const;

            void findNearest
            (
                const labelList& meshFaces,
                List<pointIndexHit>& nearestInfo,
                labelList& nearestSurface,
                labelList& nearestRegion,
                vectorField& nearestNormal
            ) const;

            Map<label> findEdgeConnectedProblemCells
            (
                const scalarField& perpendicularAngle,
                const labelList&
            ) const;

            bool isCollapsedFace
            (
                const pointField&,
                const pointField& neiCc,
                const scalar minFaceArea,
                const scalar maxNonOrtho,
                const label facei
            ) const;

            bool isCollapsedCell
            (
                const pointField&,
                const scalar volFraction,
                const label celli
            ) const;

            void markConnected
            (
                const label cellI,
                const boolList& marked,
                const boolList& isBoundaryFace,
                const boolList& toBaffle,
                boolList& visited,
                DynamicList<label>& walkedCells,
                bool& ownerSide
            ) const;

            //- Look for situation that will cause zone snapping problems
            void rezoneProblemFaces
            (
                const dictionary& motionDict,
                const labelList& surfaceToCellZone,
                labelList& namedSurfaceIndex,
                labelList& cellToZone
            ) const;

            //- Returns list with for every face the label of the nearest
            //  patch. Any unreached face (disconnected mesh?) becomes
            //  adaptPatchIDs[0]
            labelList nearestPatch(const labelList& adaptPatchIDs) const;

            //- Returns list with for every internal face -1 or the patch
            //  they should be baffled into. If removeEdgeConnectedCells is set
            //  removes cells based on perpendicularAngle.
            labelList markFacesOnProblemCells
            (
                const dictionary& motionDict,
                const labelList& surfaceToCellZone,
                labelList& namedSurfaceIndex,
                labelList& cellToZone
            ) const;

            //- Returns list with for every internal face -1 or the patch
            //  they should be baffled into.
            labelList markFacesOnProblemCellsGeometric
            (
                const snapParameters& snapParams,
                const dictionary& motionDict,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch
            ) const;

        // Split handling

           void walkEdgeDir
           (
               const label celli,
               const vector dir,
               const label start,
               const label end,
               DynamicList<label>& loop
            );

            void diagonalCut
            (
                const cellFeatures& cellFeat,
                const Map<label>& pointMap,
                const boolList&  featureEdgePoints,
                const label start,
                const label end,
                DynamicList<label>& loop,
                boolList& cutFaces
             );

            //Create cell loop for tetrahedral splitting
            labelList tetSplit
            (
                const cellFeatures& cellFeat,
                const boolList& cutNodes,
                const cell& c,
                boolList& cutFaces
             );

            //Create cell loop for prism splitting
            labelList prismSplit
            (
                const cellFeatures& cellFeat,
                const boolList& cutNodes,
                const label& celli,
                const cell& c,
                boolList& cutFaces
             );

            //Create cell loop for pyramid splitting
            //creates non-manifold tet converted to pyramid in splitPyramids
            labelList pyrSplit
            (
                const cellFeatures& cellFeat,
                const boolList& cutNodes,
                const label& celli,
                const cell& c,
                const label& pyrStart,
                const label& pyrEnd,
                boolList& cutFaces
             );

            labelList pyrSplitSideFace
            (
                const cellFeatures& cellFeat,
                const label& celli,
                const cell& c,
                const label& pyrStart,
                const label& pyrMid,
                const label& pyrEnd,
                boolList& cutFaces
             );

            //Split mesh pyramids
            void splitPyramids
            (
                const List<labelList>& pyrStartEnd,
                const List<scalarField>& cellWeightsPyr,
                List<labelList>& cellLoopsPyr,
                labelList& cutCellsPyr,
                pointField& movedPoints
            );

            //Split parallel faces not handled by meshCutter
            void splitProcessorFaces
            (
                List<label>& cutCells,
                List<label>& cutCellsPyr,
                List<labelList>& cellLoops,
                List<labelList>& cellLoopsPyr,
                List<labelList>& pyrStartEnd,
                pointField& movedPoints
            );

            //Select how to split the cells
            void selectCutType
            (
                const boolList& excludeTetSplit,
                const boolList& excludedCells,
                const boolList& cutNodes,
                const vectorField& cutNormals,
                DynamicList<label>& cutCells,
                DynamicList<label>& cutCellsPyr,
                DynamicList<labelList>& cellLoops,
                DynamicList<labelList>& cellLoopsPyr,
                DynamicList<labelList>& pyrStartEnd
             );

            //Filter out any problematic cell cuts
            void filterCutCells
            (
                List<label>& cutCells,
                List<label>& cutCellsPyr,
                List<labelList>& cellLoops,
                List<labelList>& cellLoopsPyr,
                List<labelList>& pyrStartEnd
             );

             // reset displacement of non cut points
             void resetCutDisplacements
             (
                 const List<labelList>& cellLoops,
                 pointField& movedPoints
             );

             //select which cells not to split
             boolList setRefinementCells() const;

             //select which nodes to cut
             boolList calculateCutNodes
             (
                 vectorField& cutNormals
             ) const;

             //trye an improve cut nodes form topological checks
             boolList updateCutNodes
             (
                 const vectorField& cutNormals,
                 boolList& cutNodes,
                 pointField& movedPoints
             ) const;

             //reallocate face owners for incorrect non manifold split
             void resetNonManifoldSplits
             (
                 const labelList& newFaceOwner,
                 const labelList& newFaceNeighbour,
                 List<labelList>& cellLoopsPyr,
                 List<label>& cutCellsPyr,
                 List<labelList>& pyrStartEnd,
                 pointField& movedPoints
              );

             //cut the selected cells
             void splitCells
             (
                 List<label>& cutCells,
                 List<label>& cutCellsPyr,
                 List<labelList>& cellLoops,
                 List<labelList>& cellLoopsPyr,
                 List<labelList>& pyrStartEnd,
                 pointField& movedPoints
             );

        // Wrapper handling

             autoPtr<regionSplit> splitWrappedMesh
             (
                 const boolList& cutoffCells
              ) const;

            //- Modify fvSchemes and fvSolution dictionaries to suit solver
            void setFvDicts() const;

            //- Create volumetric sources if required
            void sourceSetup
            (
                const scalar& volDist,
                const scalar& Vs,
                volScalarField& dVs
            ) const;

            //- Solve Euler equation to identify leaks and create wrap
            labelList solve
            (
                const dictionary& dict,
                const pointField& locationsInMesh,
                const labelList& keepCells
            ) const;


        // Baffle merging

            //- Extract those baffles (duplicate) faces that are on the edge
            //  of a baffle region. These are candidates for merging.
            List<labelPair> freeStandingBaffles
            (
                const List<labelPair>&,
                const scalar freeStandingAngle
            ) const;

        // Zone handling

            //- Finds zone per cell for cells inside closed named surfaces.
            void findCellZone
            (
                const labelList& namedSurfaces,
                const labelList& surfaceToCellZone,
                labelList& namedSurfaceIndex,
                const PackedList<1>& rebaffledCells,
                labelList& cellToZone
            ) const;



            //- Finds zone per cell for cells inside closed named surfaces.
            //  (uses geometric test for insideness)
            //  Adapts namedSurfaceIndex so all faces on boundary of cellZone
            //  have corresponding faceZone.
            void findCellZoneGeometric
            (
                const pointField& neiCc,
                const labelList& closedNamedSurfaces,
                labelList& namedSurfaceIndex,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            //- Finds zone per cell for cells inside region for which name
            //  is specified.
            void findCellZoneInsideWalk
            (
                const pointField& locationsInMesh,
                const labelList& zonesInMesh,
                const labelList& blockedFace, // per face -1 or some index >= 0
                labelList& namedSurfaceIndex,
                labelList& cellToZone
            ) const;

            //- Finds zone per cell for cells inside region for which name
            //  is specified.
            void findCellZoneInsideWalk
            (
                const pointField& locationsInMesh,
                const wordList& zoneNamesInMesh,
                const labelList& faceToZone, // per face -1 or some index >= 0
                labelList& namedSurfaceIndex,
                labelList& cellToZone
            ) const;

            //- Finds zone per cell for cells inside closed named surfaces.
            void findCellZoneClosed
            (
                const labelList& closedNamedSurfaces,   // indices of named surfaces
                const labelList& surfaceToCellZone,     // cell zone index per surface
                labelList& namedSurfaceIndex,     // per face index of named surface
                labelList& cellToZone,
                label minZoneRegionSize
            ) const;

            //- Finds zone per cell. Uses topological walk with all faces
            //  marked in unnamedSurfaceRegion (intersections with unnamed
            //  surfaces) and namedSurfaceIndex (intersections with named
            //  surfaces) regarded as blocked.
            void findCellZoneTopo
            (
                const label backgroundZoneID,
                const List<point>& locationsInMesh,
                const labelList& unnamedSurfaceRegion,
                const labelList& namedSurfaceIndex,
                const labelList& surfaceToCellZone,
                labelList& cellToZone
            ) const;

            //- Finds size of split regions and removes from zones if too small
            void excludeSmallRegions
            (
                const labelHashSet& namedLocationsSet,
                labelList& globalRegion1,
                labelList& globalRegion2,
                labelList& namedSurfaceIndex,
                labelList& cellToZone,
                const label minSize,
                const Switch namedLocationsRezone
             ) const;

            labelList selectPinchedEdges
            (
                 const labelList& cellToZone,
                 const labelList& keepCell
            ) const;

            label resetRegionFaces
            (
                const labelList& keepCell,
                labelList& cellToZone,
                labelList& globalRegion1,
                labelList& globalRegion2,
                labelList& testFaces,
                label defaultRegion
            ) const;

            //- Move single outside cells to inside in gaps
            void addOutsideContactCells
            (
                const refinementParameters& refineParams,
                const pointField& cellCentres,
                labelList& globalRegion1,
                labelList& globalRegion2,
                labelList& cellToZone,
                labelList& testFaces
            ) const;

            //- Try and find cells that can be kept in corner regions
            void reallocateCornerCells
            (
                const refinementParameters& refineParams,
                labelList& globalRegion1,
                labelList& globalRegion2,
                labelList& cellToZone,
                labelList& testFaces
            ) const;

            void reallocateSingleRemovalCells
            (
                labelList& globalRegion1,
                labelList& globalRegion2,
                labelList& cellToZone,
                labelList& testFaces
            ) const;

            void zoneProblemCellReallocate
            (
                const labelList& globalRegion1,
                const labelList& globalRegion2,
                labelList& namedSurfaceIndex,
                labelList& cellToZone
            ) const;

            //- Make namedSurfaceIndex consistent with cellToZone
            //  - clear out any blocked faces inbetween same cell zone.
            void makeConsistentFaceIndex
            (
                const labelList& cellToZone,
                labelList& namedSurfaceIndex
            ) const;

            //- Calculate cellZone allocation
            void zonify
            (
                const label backgroundZoneID,
                const refinementParameters& refineParams,
                const label minZoneRegionSize,

                labelList& testFaces,
                labelList& globalRegion1,
                labelList& globalRegion2,

                labelList& cellToZone,
                labelList& namedSurfaceIndex,
                labelList& namedIntersections,
                PackedBoolList& posOrientation
            ) const;

            //- Put cells into cellZone, faces into faceZone
            void zonify
            (
                const PackedBoolList& isMasterFace,
                const labelList& cellToZone,
                const labelList& neiCellZone,
                const labelList& faceToZone,
                const PackedBoolList& meshFlipMap,
                polyTopoChange& meshMod
            ) const;

            //- Check for leaks in cellZones
            void checkZoneLeakPaths
            (
                const refinementParameters& refineParams,
                const labelList& namedSurfaceIndex,
                const labelList& globalRegion1,
                const labelList& globalRegion2,
                const labelList& cellToZone,
                const writer<scalar>& leakPathFormatter
            ) const;

            //- Allocate faceZoneName
            void allocateInterRegionFaceZone
            (
                const label ownZone,
                const label neiZone,
                wordPairHashTable& zonesToFaceZone,
                LabelPairMap<word>& zoneIDsToFaceZone
            ) const;

            //- Remove any loose standing cells
            void handleSnapProblems
            (
                const snapParameters& snapParams,
                const bool useTopologicalSnapDetection,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const dictionary& motionDict,
                Time& runTime,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const refinementParameters& refineParams
            );

            // Some patch utilities

            //- Get all faces in faceToZone that have no cellZone on
            //  either side.
            labelList freeStandingBaffleFaces
            (
                const labelList& faceToZone,
                const labelList& cellToZone,
                const labelList& neiCellZone
            ) const;

            //- Determine per patch edge the number of master faces. Used
            //  to detect non-manifold situations.
            void calcPatchNumMasterFaces
            (
                const PackedBoolList& isMasterFace,
                const indirectPrimitivePatch& patch,
                labelList& nMasterFaces
            ) const;

            //- Determine per patch face the (singly-) connected zone it
            //  is in. Return overall number of zones.
            label markPatchZones
            (
                const indirectPrimitivePatch& patch,
                const labelList& nMasterFaces,
                labelList& faceToZone
            ) const;

            //- Make faces consistent.
            void consistentOrientation
            (
                const PackedBoolList& isMasterFace,
                const indirectPrimitivePatch& patch,
                const labelList& nMasterFaces,
                const labelList& faceToZone,
                const Map<label>& zoneToOrientation,
                PackedBoolList& meshFlipMap
            ) const;

        //- Disallow default bitwise copy construct
        meshRefinement(const meshRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const meshRefinement&);

public:

    //- Runtime type information
    ClassName("meshRefinement");


    // Constructors

        //- Construct from components
        meshRefinement
        (
            fvMesh& mesh,
            const scalar mergeDistance,
            const bool readLevel,
            const bool overwrite,
            const meshControl& controller,
            const Tuple2<bool, scalar> crackDetection,
            const bool adddRays,
            const refinementSurfaces&,
            const shellSurfaces&,
            const shellSurfaces&,
            autoPtr<refinementFeatures>&,
            const scalar manualLevel0,
            const bool useVDB = false
        );


    // Member Functions

        // Access

            void clearOldPoints()
            {
                deleteDemandDrivenData(oldPointsPtr_);
            };

            pointField& oldPoints() const
            {
                return *oldPointsPtr_;
            };

            //- Reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }
            fvMesh& mesh()
            {
                return mesh_;
            }

            const labelIOList& cellRegionID() const
            {
                return cellRegionID_;
            }

            void updateRegionID
            (
                const labelList& newCellRegion
            )
            {
                cellRegionID_ = newCellRegion;
            }

            scalar mergeDistance() const
            {
                return mergeDistance_;
            }

            //- Overwrite the mesh?
            bool overwrite() const
            {
                return overwrite_;
            }

            //- Whether to perform expensive crack detection
            bool checkForCracks() const
            {
                return crackDetection_.first();
            }

            //Activate/Deactivate cell level wrapping
            void setWrapActive(const bool active)
            {
                wrapActive_ = active;
            }

            //Activate/Deactivate gap field creation
            void setGapFieldCreation(const bool active)
            {
                createGapField_ = active;
            }

            //- Crack tolerance based on minimum edge length
            scalar crackTol() const
            {
                return crackTol_;
            }

            //- Crack tolerance based on minimum edge length
            bool addedRays() const
            {
                return addedRays_;
            }

            //- (points)instance of mesh upon construction
            const word& oldInstance() const
            {
                return oldInstance_;
            }

            //- Reference to surface search engines
            const refinementSurfaces& surfaces() const
            {
                return surfaces_;
            }

            //- Reference to feature edge mesh
            const autoPtr<refinementFeatures>& features() const
            {
                return featuresPtr_;
            }

            //- Reference to refinement shells (regions)
            const shellSurfaces& shells() const
            {
                return shells_;
            }

            autoPtr<refinementFeatures>& features()
            {
                return featuresPtr_;
            }

            //- Reference to meshcutting engine
            const hexRef8& meshCutter() const
            {
                return meshCutter_;
            }

            //- Reference to meshcutting engine
            hexRef8& meshCutter()
            {
                return meshCutter_;
            }

            //- Per mesh face whether face is wrapped
            const labelList& wrapIndex() const
            {
                return wrapIndex_;
            }

            labelList& wrapIndex()
            {
                return wrapIndex_;
            }

            //- Per start-end edge the index of the surface hit
            const labelList& surfaceIndex() const
            {
                return surfaceIndex_;
            }

            labelList& surfaceIndex()
            {
                return surfaceIndex_;
            }

            //- For faces originating from processor faces store the original
            //  patch
            const Map<label>& faceToCoupledPatch() const
            {
                return faceToCoupledPatch_;
            }

            //- Additional face data that is maintained across
            //  topo changes. Every entry is a list over all faces.
            //  Bit of a hack. Additional flag to say whether to maintain master
            //  only (false) or increase set to account for face-from-face.
            const List<Tuple2<mapType, labelList>>& userFaceData() const
            {
                return userFaceData_;
            }

            List<Tuple2<mapType, labelList>>& userFaceData()
            {
                return userFaceData_;
            }

            const meshControl& controller() const
            {
                return controller_;
            }

        // Other

            // Reset face zone flip map (excluding freestanding zones)
            void resetFaceZoneFlipMap();

            //- Remove cells. Put exposedFaces into exposedPatchIDs.
            autoPtr<mapPolyMesh> doRemoveCells
            (
                const labelList& cellsToRemove,
                const labelList& exposedFaces,
                const labelList& exposedPatchIDs,
                removeCells& cellRemover,
                const bool updateIntersection = true
             );

            //- Returns list with for every internal face -1 or the patch
            //  they should be baffled into. If removeEdgeConnectedCells is set
            //  removes cells based on perpendicularAngle.
            labelList markFacesOnProblemCells
            (
                const dictionary& motionDict,
                const refinementParameters& refineParams,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const labelList& globalToMasterPatch,
                const labelList& adaptPatchIDs,
                const labelList& checkSurfaces,
                const List<labelPair> removeFaces = List<labelPair>()
            ) const;

            //- reallocate zones to remove holes
            label removeZoneHoles(const bool external);

            //- Close gaps with a single cell across them
            void closeSingleCellGaps
            (
                const labelList& closureSurfaces,
                const pointField& locationsInMesh,
                const labelList& globalRegion1,
                const labelList& wrapLevel,
                labelList& wrappedFaces
            ) const;

            //- Count number of intersections (local)
            label countHits
            (
                const bool threaded = false
            ) const;

            //- Create non unit weight field for decomposer to estimate
            // work required during snapping
            void calcSnapWeights
            (
                const scalar& snapWeights,
                scalarField& weightField
            ) const;

            //- Redecompose according to cell count
            //  keepZoneFaces : find all faceZones from zoned surfaces and keep
            //                  owner and neighbour together
            //  keepBaffles   : find all baffles and keep them together
            autoPtr<mapDistributePolyMesh> balance
            (
                const bool applyDictionaryConstraints,
                const bool keepZoneFaces,
                const bool keepBaffles,
                const scalarField& cellWeights,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                bool updateSurf = true
            );

            //- Redecompose according to cell count with explicit blockFaces set
            autoPtr<mapDistributePolyMesh> balance
            (
                const boolList blockedFace,
                const scalarField& cellWeights,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                bool updateSurf = true
            );

            //- Redistribute to higher number of processors
            void redistributeToMany
            (
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor
            );

            //- Print statistics for list of meshes
            void printMeshData(PtrList<fvMesh>& meshes);

            //- Get faces with intersection.
            labelList intersectedFaces() const;

            //- Get faces with intersection and neigbouring faces.
            labelList intersectedAndNeighbouring() const;

            //- Get faces with intersection and nedge cell faces.
            labelList intersectedAndEdgeCellFaces() const;

            //- Get points on surfaces with intersection and boundary faces.
            labelList intersectedPoints() const;

            //- Get cells on surfaces with intersection and boundary faces.
            labelList intersectedCells() const;

            //- Populate face addressing array for faces in given patches
           static void calcPatchAddressing
           (
               const polyMesh&,
               const labelList&,
               labelList&
           );

            //- Create patch from set of patches
            static autoPtr<indirectPrimitivePatch> makePatch
            (
                const polyMesh& mesh,
                const labelList& includedPatchIDs
            );

            //- Helper function to make a pointVectorField with correct
            //  bcs for mesh movement:
            //  - adaptPatchIDs         : fixedValue
            //  - processor             : calculated (so free to move)
            //  - cyclic/wedge/symmetry : slip
            //  - other                 : slip
            static tmp<pointVectorField> makeDisplacementField
            (
                const pointMesh& pMesh,
                const labelList& adaptPatchIDs
            );

            //- Helper function: check that face zones are synced
            static void checkCoupledFaceZones(const polyMesh&);

            //- Helper: calculate edge weights (1/length)
            static void calculateEdgeWeights
            (
                const polyMesh& mesh,
                const PackedBoolList& isMasterEdge,
                const boolList& isBoundaryEdge,
                const boolList& isBoundaryPoint,
                const labelList& meshEdges,
                const labelList& meshPoints,
                const edgeList& edges,
                scalarField& edgeWeights,
                scalarField& invSumWeight
            );

            //- Helper: weighted sum (over all subset of mesh points) by
            //  summing contribution from (master) edges
            template<class Type>
            static void weightedSum
            (
                const polyMesh& mesh,
                const PackedBoolList& isMasterEdge,
                const boolList& isBoundaryEdge,
                const boolList& isBoundaryPoint,
                const labelList& meshEdges,
                const labelList& meshPoints,
                const edgeList& edges,
                const scalarField& edgeWeights,
                const Field<Type>& pointData,
                Field<Type>& sum
            );

            // Repatch region
            void repatchRegion(const dictionary& dict);

            // Repatch based on input surfaces
            void surfaceToPatch(const dictionary& dict);

            //- Return boundary edgeMesh from an input set of mesh patches
            static edgeMesh calcBoundaryEdgeMesh
            (
                const polyMesh& mesh,
                const labelHashSet& edgeSet
            );

            //- Source/Target boundary edge error calculation
            // from input patch lists
            static void calcSourceTargetMatch
            (
                const polyMesh& mesh,
                const wordReList& sourcePatches,
                const wordReList& targetPatches,
                const bool writeVTK
            );

            //- Source/Target boundary edge error calculation
            // from an input dictionary
            static void sourceTargetChecking
            (
                const polyMesh& mesh,
                const dictionary& dict
            );

            //- Triangulate problem boundary faces
            void triangulateProblemBoundaryFaces
            (
                const dictionary& triangulateDict
            );

        // Refinement

            //- Calculate maximum cell level
            label maxCellLevel() const;

            //- Is local topology a small gap?
            bool isGap
            (
                const scalar,
                const vector&,
                const vector&,
                const vector&,
                const vector&
            ) const;

            //- Is local topology a small gap normal to the test vector
            bool isNormalGap
            (
                const scalar,
                const vector&,
                const vector&,
                const vector&,
                const vector&
            ) const;

            //- Calculate list of gap cells to refine.
            labelList refineGaps();

            //- Calculate list of removal interface cells to refine
            labelList markRemovalInterfaceRefine
            (
                const refinementParameters& refineParams
            ) const;

            //- Calculate list of cells to refine.
            labelList refineCandidates
            (
                const label iterationNum,
                const refinementParameters& refineParams,

                const bool initialRefinement,
                const bool featureRefinement,
                const bool featureDistanceRefinement,
                const bool internalRefinement,
                const bool surfaceRefinement,
                const bool curvatureRefinement,
                const bool proximityRefinement,
                const bool smallFeatureRefinement,
                const bool gapRefinement,
                const bool bigGapRefinement,
                const bool spreadGapSize
            );

            //- Refine some cells
            autoPtr<mapPolyMesh> refine
            (
                const labelList& cellsToRefine,
                const pointField& locationsInMesh
            );

            //- Balance the mesh
            autoPtr<mapDistributePolyMesh> balance
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                labelList& cellsToRefine,
                const refinementParameters& refineParams
            );

            //- Refine some cells and rebalance
            autoPtr<mapDistributePolyMesh> refineAndBalance
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const labelList& cellsToRefine,
                const refinementParameters& refineParams
            );

            //- Balance before refining some cells
            autoPtr<mapDistributePolyMesh> balanceAndRefine
            (
                const string& msg,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const labelList& cellsToRefine,
                const refinementParameters& refineParams
            );

            // Refine cell other side of intersected face
            labelList refineVoids();


        // Baffle handling

            static void updateBaffles
            (
                const mapPolyMesh &map,
                List<labelPair>& baffles,
                const List<labelPair> added = List<labelPair>()
             );

            //- Split off unreachable areas of mesh.
            void baffleAndSplitMesh
            (
                const bool handleSnapProblems,

                // How to remove problem snaps
                const snapParameters& snapParams,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const dictionary& motionDict,
                Time& runTime,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const refinementParameters& refineParams,
                const writer<scalar>& leakPathFormatter,
                const bool threaded = false
            );

            //- Merge free-standing baffles
            void mergeFreeStandingBaffles
            (
                const snapParameters& snapParams,
                const bool useTopologicalSnapDetection,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const scalar planarAngle,
                const dictionary& motionDict,
                Time& runTime,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const refinementParameters& refineParams,
                const writer<scalar>& leakPathFormatter,
                const bool threaded = false
            );

            //- Look for free standing holes and try to
            // filter out based the surface orientation
            void filterFreeStandingHoles
            (
                const refinementParameters& refineParams,
                labelList& namedSurfaceIndex
            ) const;

            //- Merge single couples
            void mergeSingleCouples();

            //Remove cells is disconnected regions
            autoPtr<mapPolyMesh> removeDisconnectedRegions();

            void moveWrongSidedCells
            (
                const refinementParameters& refineParams,
                bool wrongSidedCorner
            );

            //remove cells whose cell centres fall other side of surface
            void removeOtherSide
            (
                const refinementParameters& refineParams,
                const labelList& globalToMasterPatch
            );

            //remove cells whose cell centres fall other side of surface
            // and also remove cells created from originating cell
            void removeOtherSide
            (
                const refinementParameters& refineParams,
                const labelList& globalToMasterPatch,
                const mapPolyMesh& map
            );

            //Refine cells if there removal might result in a void being created
            void checkRemovedAndRefine
            (
                const refinementParameters& refineParams,
                const labelList& globalToMasterPatch,
                const label maxIter
            );

            //Input cells to refine and if any resultant refined cell
            // falls wrong side of geometry remove cell selected for refinement
            bool checkRefinedAndRemove
            (
                const labelList& cellsToRefine,
                const bool checkPointCells,
                const bool unnamedOnly
            );

            //Remove cells from dual mesh that might produce snapping issues
            label removeProblemCells();

            //Remove cells based on cellSet selection
            label removeSelectedCellSets
            (
                const dictionary& dict
            );

            //remove holes by cell deletion
            label removeHoleCells();

            //remove cells that might cause problems for dualisation
            label removeProblemDualisationCells
            (
                const labelHashSet& noLayerPatchIDs,
                const bool additionalChecks
            );

            //Detect holes and close them since this might cause problems
            // adding layers for dualised mesh
            label baffleHoles
            (
                label closeHoleSize = labelMax
            );

            //- Split off unreachable areas
            //  of mesh for shell meshing.
            void splitMesh
            (
                const refinementParameters& refineParams
            );

            //- Split off (with optional buffer layers) unreachable areas
            //  of mesh. Does not introduce baffles.
            autoPtr<mapPolyMesh> splitMesh
            (
                const label nBufferLayers,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const refinementParameters& refineParams,
                const writer<scalar>& leakPathFormatter
            );

            //- Remove any loose standing cells after zoning type boundary
            void handleZoneSnapProblems
            (
                const snapParameters& snapParams,
                const dictionary& motionDict,
                Time& runTime,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const refinementParameters& refineParams,
                List<labelPair>& baffles
            );

            //- Remove cells that might cause issues for extrusion methods
            void removeExtrusionProblemCells
            (
                const dictionary& meshDict,
                const refinementParameters& refineParams
            );

            //- Find boundary points that connect to more than one cell
            //  region and split them.
            autoPtr<mapPolyMesh> dupNonManifoldPoints(const localPointRegion&);

            //- Find boundary points that connect to more than one cell
            //  region and split them.
            autoPtr<mapPolyMesh> dupNonManifoldPoints();

            //- Find boundary points that are on faceZones of type boundary
            //  and duplicate them
            autoPtr<mapPolyMesh> dupNonManifoldBoundaryPoints();

            //- Merge duplicate points
            autoPtr<mapPolyMesh> mergePoints
            (
                const labelList& pointToDuplicate
            );

            //- Create baffle for every internal face where ownPatch != -1.
            //  External faces get repatched according to ownPatch (nbrPatch
            //  should be -1 for these)
            autoPtr<mapPolyMesh> createBaffles
            (
                const labelList& ownPatch,
                const labelList& nbrPatch,
                const bool updateIntersection = true
            );

            //- Debug helper: check faceZones are not on processor patches
            void checkZoneFaces() const;

            //- Get zones of given type
            labelList getZones
            (
                const List<surfaceZonesInfo::faceZoneType>& fzTypes
            ) const;

            //- Subset baffles according to zones
            static List<labelPair> subsetBaffles
            (
                const polyMesh& mesh,
                const labelList& zoneIDs,
                const List<labelPair>& baffles
            );

            //- Create baffles for faces on faceZones. Return created baffles
            //  (= pairs of faces) and corresponding faceZone
            autoPtr<mapPolyMesh> createZoneBaffles
            (
                const labelList& zoneIDs,
                List<labelPair>& baffles,
                labelList& originatingFaceZone,
                const bool updateIntersection = true
            );

            //- Merge baffles. Gets pairs of faces and boundary faces to move
            //  onto (coupled) patches
            autoPtr<mapPolyMesh> mergeBaffles
            (
                const List<labelPair>&,
                const Map<label>& faceToPatch,
                const bool reallocateFaceZone,
                const bool updateIntersection = true
            );

            //- Merge all baffles on faceZones
            autoPtr<mapPolyMesh> mergeZoneBaffles
            (
                const bool doInternalZones,
                const bool doBaffleZones,
                const bool updateIntersection = true
            );

            //- Put faces/cells into zones according to surface specification.
            //  Returns null if no zone surfaces present. Regions containing
            //  locationsInMesh/regionsInMesh will be put in corresponding
            //  cellZone. keepPoints is for backwards compatibility and sets
            //  all yet unassigned cells to be non-zoned (zone = -1)
            void zonify
            (
                const refinementParameters& refineParams,
                const writer<scalar>& leakPathFormatter,
                wordPairHashTable& zonesToFaceZone,
                bool keepAlreadyZonedCells,
                bool resetUnsetZones,
                bool removeUnsetCells
            );

        // Geometric wrapper handling

            void markWrapFaces(const pointField& locationsInMesh);

        // Solver Wrapper handling

            //- Remove pressurised containers from mesh
            void splitWrappedMesh
            (
                const refinementParameters& refineParams,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const dictionary& dict,
                const dictionary& motionDict,
                const bool handleSnapProblems,
                const bool removeEdgeConnectedCells,
                const scalarField& perpendicularAngle,
                const writer<scalar>& leakPathFormatter
            );

        // Thin gap handling

            //- Select thin gap cells for keeping
            void markThinGapCells
            (
                const pointField& cellCentres,
                boolList& gapCells,
                boolList& protectedFaces,
                DynamicList<label>& baffleFaces
            ) const;

            //- Re-baffle selected thin gap cells
            void rebaffleThinGapCells
            (
                const pointField& locationsInMesh,
                const pointField& cellCentres,
                labelList& globalRegion1,
                labelList& globalRegion2,
                labelList& testFaces
            ) const;

            //- Ensure no internal points in gap
            void createSingleCellGap();

        // Cell cutting

            //- Perform edge snapping as precursor to splitting
            void snapAndCut();

        // Other topo changes

            //- Helper: keep the largest regions
            labelList keepLargestRegions
            (
                const regionSplit& cellRegion
            );

            //- Reset master region id
            void updateMasterRegions();

            //- Helper:append patch to end of mesh.
            static label appendPatch
            (
                fvMesh&,
                const word&,
                const dictionary&
            );

            //- Helper:add patch to mesh. Update all registered fields.
            //  Used by addMeshedPatch to add patches originating from surfaces.
            static label addPatch(fvMesh&, const word& name, const dictionary&);

            //- Add patch originating from meshing. Update meshedPatches_.
            label addMeshedPatch(const word& name, const dictionary&);

            //- Remove any zero-sized patches and accompanying fields
            static void removeZeroSizedPatches(fvMesh& mesh);

            //- Get patchIDs for patches added in addMeshedPatch.
            labelList meshedPatches() const;

            //- Get patchIDs for original patches.
            labelList unmeshedPatches() const;

            //- Add/lookup faceZone and update information. Return index of
            //  faceZone
            label addFaceZone
            (
                const word& fzName,
                const word& masterPatch,
                const word& slavePatch,
                const surfaceZonesInfo::faceZoneType& fzType
            );

            //- Lookup faceZone information. Return false if no information
            //  for faceZone
            bool getFaceZoneInfo
            (
                const word& fzName,
                label& masterPatchID,
                label& slavePatchID,
                surfaceZonesInfo::faceZoneType& fzType
            ) const;

            //- Select coupled faces that are not collocated
            void selectSeparatedCoupledFaces(boolList&) const;

            //- Find any intersection of surface. Store in surfaceIndex_.
            void updateIntersections
            (
                const labelList& changedFaces,
                const bool threaded = false
            );

            //- Find cell associated with keepPoint
            static label findCell
            (
                const point& keepPoint,
                const polyMesh& mesh,
                const hexRef8& meshCutter
            );

            //- Find region point is in. Uses optional perturbation to re-test.
            static label findRegion
            (
                const polyMesh&,
                const labelList& cellRegion,
                const vector& perturbVec,
                const point& p
            );

            static void findRegions
            (
                const polyMesh&,
                const hexRef8&,
                const vector& perturbVec,
                const pointField& locationsInMesh,
                const pointField& locationsOutsideMesh,
                const writer<scalar>& leakPathFormatter,
                const bool& calcLeakPath,
                const label nRegions,
                labelList& cellRegion,
                const boolList& blockedFace
            );

            static void calculateLeakPath
            (
                const word leakType,
                const word leakPathName,
                const polyMesh& mesh,
                const boolList& blockedFace,
                const pointField& locationsInMesh,
                const pointField& locationsOutsideMesh,
                const writer<scalar>& leakPathFormatter
            );

            //- Split mesh. Keep part containing point. Return empty map if
            //  no cells removed.
            autoPtr<mapPolyMesh> splitMeshRegions
            (
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const pointField& locationsInMesh,
                const pointField& locationsOutsideMesh,
                const writer<scalar>& leakPathFormatter,
                const bool& calcLeakPath
            );

            //- Split faces into two
            void doSplitFaces
            (
                const labelList& splitFaces,
                const labelPairList& splits,
                polyTopoChange& meshMod
            ) const;

            //- Split faces along diagonal. Maintain mesh quality. Return
            //  total number of faces split.
            label splitFacesUndo
            (
                const labelList& splitFaces,
                const labelPairList& splits,
                const dictionary& motionDict,

                labelList& duplicateFace,
                List<labelPair>& baffles
            );

            //- Update local numbering for mesh redistribution
            void distribute
            (
                const mapDistributePolyMesh&,
                const bool& updateSurf
            );

            //- Check mesh not too large for 32 bit label compilation
            void checkMeshSize();

            //- Update for external change to mesh. changedFaces are in new mesh
            //  face labels.
            void updateMesh
            (
                const mapPolyMesh&,
                const labelList& changedFaces,
                const bool updateIntersection = true
            );

            //- Helper: reorder list according to map.
            template<class T>
            static void updateList
            (
                const labelList& newToOld,
                const T& nullValue,
                List<T>& elems
            );


            // Restoring : is where other processes delete and reinsert data.

                //- Signal points/face/cells for which to store data
                void storeData
                (
                    const labelList& pointsToStore,
                    const labelList& facesToStore,
                    const labelList& cellsToStore
                );

                //- Update local numbering + undo
                //  Data to restore given as new pointlabel + stored pointlabel
                //  (i.e. what was in pointsToStore)
                void updateMesh
                (
                    const mapPolyMesh&,
                    const labelList& changedFaces,
                    const Map<label>& pointsToRestore,
                    const Map<label>& facesToRestore,
                    const Map<label>& cellsToRestore,
                    const bool updateIntersection = true
                );


                //- Select faces in proximity of surface
                // (in cellZone innerGrid) for mesh checking
                labelList selectInnerFaces() const;

            // Merging coplanar faces and edges

                //- Merge coplanar faces if sets are of size mergeSize
                //  (usually 4)
                label mergePatchFaces
                (
                    const scalar minCos,
                    const scalar concaveCos,
                    const label mergeSize,
                    const labelList& patchIDs
                );

                //- Merge coplanar faces. preserveFaces is != -1 for faces
                //  to be preserved
                label mergePatchFacesUndo
                (
                    const scalar minCos,
                    const scalar concaveCos,
                    const scalar maxAreaRatio,
                    const labelList& patchIDs,
                    const dictionary& motionDict,
                    const bool updateIntersection,
                    const bool maintainPatches,
                    const List<labelList> &featureShadows = List<labelList>(),
                    const vectorField* surfNormal = nullptr,
                    List<labelPair>* baffles = nullptr
                );

                autoPtr<mapPolyMesh> doRemovePoints
                (
                    removePoints& pointRemover,
                    const boolList& pointCanBeDeleted,
                    const bool updateIntersection = true
                );

                autoPtr<mapPolyMesh> doRestorePoints
                (
                    removePoints& pointRemover,
                    const labelList& facesToRestore,
                    const bool updateIntersection = true
                );

                labelList collectFaces
                (
                    const labelList& candidateFaces,
                    const labelHashSet& set
                ) const;

                // Pick up faces of cells of faces in set.
                labelList growFaceCellFace
                (
                    const UList<label>& set
                ) const;

                // Pick up faces of cells of faces in set.
                labelList growFaceCellFace
                (
                    const labelHashSet& set
                ) const;

                //- Merge edges, maintain mesh quality. Return global number
                //  of edges merged
                label mergeEdgesUndo
                (
                    const scalar minCos,
                    const dictionary& motionDict,
                    const bool updateIntersection,
                    List<labelPair>* baffles = nullptr
                );

                //- Remove any unused points in the mesh
                static void removeUnusedPoints(polyMesh& mesh);

        // Debug/IO

            //- Debugging: check that all faces still obey start()>end()
            void checkData();

            static void testSyncPointList
            (
                const string& msg,
                const polyMesh& mesh,
                const List<scalar>& fld
            );

            static void testSyncPointList
            (
                const string& msg,
                const polyMesh& mesh,
                const List<point>& fld
            );

            //- Compare two lists over all boundary faces
            template<class T>
            void testSyncBoundaryFaceList
            (
                const scalar mergeDistance,
                const string&,
                const UList<T>&,
                const UList<T>&
            ) const;

            //- Print list according to (collected and) sorted coordinate
            template<class T>
            static void collectAndPrint
            (
                const UList<point>& points,
                const UList<T>& data
            );

            //- Determine master point for subset of points. If coupled
            //  chooses only one
            static PackedBoolList getMasterPoints
            (
                const polyMesh& mesh,
                const labelList& meshPoints
            );

            //- Determine master edge for subset of edges. If coupled
            //  chooses only one
            static PackedBoolList getMasterEdges
            (
                const polyMesh& mesh,
                const labelList& meshEdges
            );

            //- Print some mesh stats.
            void printMeshInfo(const bool, const string&) const;

            //- Replacement for Time::timeName() : return oldInstance (if
            //  overwrite_)
            word timeName() const;

            //- Set instance of all local IOobjects
            void setInstance(const fileName&);

            //- Write mesh and all data
            bool write(bool writeSurf = true) const;

            //- Write refinement level as volScalarFields for postprocessing
            void dumpRefinementLevel() const;

            //- Debug: Write intersection information to OBJ format
            void dumpIntersections(const fileName& prefix) const;

            //- Do any one of above IO functions
            void write
            (
                const debugType debugFlags,
                const writeType writeFlags,
                const fileName&,
                bool writeSurf = true
            ) const;

            //- Helper: remove all relevant files from mesh instance
            static void removeFiles(const polyMesh&);

            //- Helper: calculate average
            template<class T>
            static T gAverage
            (
                const PackedBoolList& isMasterElem,
                const UList<T>& values
            );

            //- Get/set write level
            static writeType writeLevel();
            static void writeLevel(const writeType);

            //- Get/set output level
            static outputType outputLevel();
            static void outputLevel(const outputType);


            //- Helper: convert wordList into bit pattern using provided
            //  NamedEnum
            template<class Enum>
            static int readFlags(const Enum& namedEnum, const wordList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshRefinement/meshRefinementTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
