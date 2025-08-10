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
    (c) 2017 OpenCFD Ltd.
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::snappyLayerDriver

Description
    All to do with adding layers

SourceFiles
    snappyLayerDriver.C

\*---------------------------------------------------------------------------*/

#ifndef snappyLayerDriver_H
#define snappyLayerDriver_H

#include "meshRefinement/meshRefinement.H"
#include "hexReport/hexReport.H"
#include "meshes/primitiveShapes/plane/plane.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class removePoints;
class pointSet;
class motionSmoother;
class addPatchCellLayer;
class faceSet;
class layerParameters;

/*---------------------------------------------------------------------------*\
                           Class snappyLayerDriver Declaration
\*---------------------------------------------------------------------------*/

class snappyLayerDriver
{
public:

    // Public data types

        //- Extrusion controls
        enum extrudeMode
        {
            NOEXTRUDE,      //!< Do not extrude. No layers added.
            EXTRUDE,        //!< Extrude
            EXTRUDEREMOVE   //!< Extrude but afterwards remove added
                            //!< faces locally
        };

private:

    // Private classes

        //- Combine operator class for equalizing displacements.
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
        class nomalsCombine
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (y != vector(GREAT, GREAT, GREAT))
                {
                    if (x == vector(GREAT, GREAT, GREAT))
                    {
                        x = y;
                    }
                    else
                    {
                        x += y;
                        x /= mag(x) + SMALL;
                    }
                }
            }
        };

        //- Combine operator class to combine normal with other normal.
        class minusCombine
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (y != vector(GREAT, GREAT, GREAT))
                {
                    if (x == vector(GREAT, GREAT, GREAT))
                    {
                        x = y;
                    }
                    else
                    {
                        x -= y;
                    }
                }
            }
        };

        //- Find a unique edge direction.
        struct uniqueEdgeNormal
        {
            public:

            void operator()(vector& x, const vector& y) const
            {
                if (x != vector(GREAT, GREAT, GREAT))
                {
                    return;
                }
                else if (y != vector(GREAT, GREAT, GREAT))
                {
                    x = y;
                }
            }
        };

    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        autoPtr<searchableSurfaces>& allGeometryPtr_;

        //- Mesh report
        hexReport& stats_;

        //- Control of meshing process
        const meshControl& controller_;

        //- From surface region to patch
        const labelList globalToMasterPatch_;

        //- From surface region to patch
        const labelList globalToSlavePatch_;

        //- List of layer patches
        labelList layerPatchIDs_;

        //- List of zone layer patches
        labelList layerZonePatchIDs_;
    // Private Member Functions

        // Layers
            // return list of layer patches
            const labelList& layerPatchIDs() const
            {
                return layerPatchIDs_;
            }
            // return list of layer patches
            const labelList& layerZonePatchIDs() const
            {
                return layerZonePatchIDs_;
            }
            // Set list of layer patches
            void setLayerPatchIDs(const labelList& layerPatchIDs)
            {
                layerPatchIDs_ = layerPatchIDs;
            }
            // Set list of zone layer patches
            void setZoneLayerPatchIDs(const labelList& layerZonePatchIDs)
            {
                layerZonePatchIDs_ = layerZonePatchIDs;
            }

            //- For debugging: Dump displacement to .obj files
            static void dumpDisplacement
            (
                const fileName&,
                const indirectPrimitivePatch&,
                const vectorField&,
                const List<extrudeMode>&
            );

            //- Average point wise data to face wise
            static tmp<scalarField> avgPointData
            (
                const indirectPrimitivePatch&,
                const scalarField& pointFld
            );

            //- Check that primitivePatch is not multiply connected.
            //  Collect non-manifold points in pointSet.
            static void checkManifold
            (
                const indirectPrimitivePatch&,
                pointSet& nonManifoldPoints
            );

            //- Check that mesh outside is not multiply connected.
            void checkMeshManifold() const;

            // Static extrusion setup

                //- Unset extrusion on point. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const label patchPointi,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Unset extrusion on face. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const face& localFace,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Truncate index in face
                static label constrainFp(const label sz, const label fp);

                //- Count common points between face and its neighbours
                void countCommonPoints
                (
                    const indirectPrimitivePatch&,
                    const label facei,
                    Map<label>&
                ) const;

                //- Check if any common points form single string. Return
                //  false if not.
                bool checkCommonOrder
                (
                    const label nCommon,
                    const face& curFace,
                    const face& nbFace
                ) const;

                //- Check if any common points form single string; unmark
                //  points on face if not
                void checkCommonOrder
                (
                    const indirectPrimitivePatch& pp,
                    const label facei,
                    const Map<label>& nCommonPoints,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Check if any common points form single string; unmark
                //  points on face if not
                void handleNonStringConnected
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion at non-manifold points.
                void handleNonManifolds
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const labelListList& edgeGlobalFaces,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Prevent layers being generated on both sides of
                // single sided patches
                void handleSingleSidedPatches
                (
                    const indirectPrimitivePatch& pp,
                    const layerParameters& layerParams,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Prevent layers being generated in excluded regions
                void handleExcludedRegions
                (
                    const indirectPrimitivePatch& pp,
                    const layerParameters& layerParams,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on feature edges. Assumes non-manifold
                //  edges already handled.
                void handleFeatureAngle
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const scalar minCos,
                    const bool growConvexEdge,
                    const bool growConcaveEdge,
                    pointField& patchDisp,
                    vectorField& patchEdgeNormals,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus,
                    PackedList<1>& isConvexEdgePoint,
                    PackedList<1>& isConcaveEdgePoint
                ) const;

                //- No extrusion on warped faces
                void handleWarpedFaces
                (
                    const indirectPrimitivePatch& pp,
                    const scalar faceRatio,
                    const scalar edge0Len,
                    const labelList& cellLevel,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion at manifold edges
                void handleNonManifoldEdges
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on selected edges
                void handleMiscEdges
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion at grown up edges where angle is too large
                void handleGrownUpEdges
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const labelList& grownUpIDs,
                    const scalar grownUpAngleTerminate,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                 ) const;

                //- Determine the number of layers per point from the number of
                //  layers per surface.
                void setNumLayers
                (
                    const labelList& patchToNLayers,
                    const labelList& patchIDs,
                    const indirectPrimitivePatch& pp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus,
                    label& nIdealAddedCells
                ) const;

                //- Helper function to make a pointVectorField with correct
                //  bcs for layer addition:
                //  - numLayers > 0         : fixedValue
                //  - numLayers == 0        : fixedValue (always zero)
                //  - processor             : calculated (so free to move)
                //  - cyclic/wedge/symmetry : slip
                //  - other                 : slip
                static tmp<pointVectorField> makeLayerDisplacementField
                (
                    const pointMesh& pMesh,
                    const layerParameters& layerParams
                );

                //- Grow no-extrusion layer.
                void growNoExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- See what zones and patches edges should be extruded into
                void determineSidePatches
                (
                    fvMesh& mesh,

                    const globalIndex& globalFaces,
                    const labelListList& edgeGlobalFaces,
                    const indirectPrimitivePatch& pp,
                    const labelList grownUpPatchIDs,

                    labelList& edgePatchID,
                    labelList& edgeZoneID,
                    boolList& edgeFlip,
                    labelList& inflateFaceID
                );

                //- Calculate pointwise wanted and minimum thickness.
                //  thickness: wanted thickness
                //  minthickness: when to give up and not extrude
                //  Gets per patch parameters and determine pp pointwise
                //  parameters.
                void calculateLayerThickness
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& patchIDs,
                    const layerParameters& layerParams,
                    const labelList& cellLevel,
                    const scalar edge0Len,

                    labelList& patchNLayers,
                    scalarField& thickness,
                    scalarField& minThickness,
                    scalarField& targetExpansion,
                    scalarField& targetfch,
                    boolList& fixedfch
                ) const;


            // Extrusion execution

                //- Make sure that a face doesn't have two non-consecutive
                // areas not extruded
                label stopNonConsecutiveExtrusion
                (
                    const globalIndex& globalFaces,
                    const labelListList& edgeGlobalFaces,
                    const indirectPrimitivePatch& pp,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Synchronize displacement among coupled patches.
                void syncPatchDisplacement
                (
                    const indirectPrimitivePatch& pp,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Get nearest point on surface to snap to
                void getPatchDisplacement
                (
                    const indirectPrimitivePatch& pp,
                    const scalarField& thickness,
                    const scalarField& minThickness,
                    const labelList& grownUpIDs,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- For truncateDisplacement: find strings of edges
                bool sameEdgeNeighbour
                (
                    const labelListList& globalEdgeFaces,
                    const label myGlobalFacei,
                    const label nbrGlobFacei,
                    const label edgeI
                ) const;

                //- For truncateDisplacement: find strings of edges
                void getVertexString
                (
                    const indirectPrimitivePatch& pp,
                    const labelListList& globalEdgeFaces,
                    const label facei,
                    const label edgeI,
                    const label myGlobFacei,
                    const label nbrGlobFacei,
                    DynamicList<label>& vertices
                ) const;

                //- Truncates displacement
                // - for all patchFaces in the faceset displacement gets set
                //   to zero
                // - all displacement < minThickness gets set to zero
                // - all non-consecutive extrusions get set to 0
                label truncateDisplacement
                (
                    const globalIndex& globalFaces,
                    const labelListList& edgeGlobalFaces,
                    const indirectPrimitivePatch& pp,
                    const scalarField& minThickness,
                    const faceSet& illegalPatchFaces,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Setup layer information (at points and faces) to
                //  modify mesh topology in
                //  regions where layer mesh terminates. Guarantees an
                //  optional slow decreasing of the number of layers.
                //  Returns the number of layers per face and per point
                //  to go into the actual layer addition engine.
                void setupLayerInfoTruncation
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& patchNLayers,
                    const List<extrudeMode>& extrudeStatus,
                    const label& nBufferCellsNoExtrude,
                    const label& terminationStrategy,
                    const label& layerRecovery,
                    const PackedList<1>& isConvexEdgePoint,
                    const PackedList<1>& isConcaveEdgePoint,
                    labelList& nPatchPointLayers,
                    labelList& nPatchFaceLayers
                ) const;

                //- Does any of the cells use a face from faces?
                static bool cellsUseFace
                (
                    const polyMesh& mesh,
                    const labelList& cellLabels,
                    const labelList& faceLabels,
                    const labelHashSet& faces
                );

                //- Checks the newly added cells and locally unmarks points
                //  so they will not get extruded next time round. Returns
                //  global number of unmarked points (0 if all was fine)
                static label checkAndUnmark
                (
                    const addPatchCellLayer& addLayer,
                    const dictionary& motionDict,
                    const bool additionalReporting,
                    const List<labelPair>& baffles,
                    const indirectPrimitivePatch& pp,
                    const fvMesh&,
                    const labelList& faceMap,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus,

                    boolList& nonConvergedAfterNlayerIter,
                    const label nLayerIter
                );

                //- Count global number of extruded faces
                static label countExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    const List<extrudeMode>& extrudeStatus
                );

                //- After adding to mesh get the new baffles
                static List<labelPair> getBafflesOnAddedMesh
                (
                    const polyMesh& mesh,
                    const labelList& newToOldFaces,
                    const List<labelPair>& baffles
                );

                //- Collect layer faces and layer cells into bools
                //  for ease of handling
                static void getLayerCellsFaces
                (
                    const polyMesh&,
                    const addPatchCellLayer&,
                    const scalarField& oldRealThickness,

                    labelList& cellStatus,
                    scalarField& faceRealThickness
                );

            // Mesh shrinking (to create space for layers)

                //- Write volFields,cellSet,faceSet for layers depending
                //  on write level
        //- Disallow default bitwise copy construct
        snappyLayerDriver(const snappyLayerDriver&);

        //- Disallow default bitwise assignment
        void operator=(const snappyLayerDriver&);


public:

    //- Runtime type information
    ClassName("snappyLayerDriver");

    // Constructors

        //- Construct from components
        snappyLayerDriver
        (
            meshRefinement& meshRefiner,
            autoPtr<searchableSurfaces>& allGeometryPtr,
            hexReport& stats,
            const meshControl& controller,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch
        );


    // Member Functions

           //- Calculate and return split faces for each patch face
           List<face> calcSplitFace
           (
               const polyMesh& mesh,
               const indirectPrimitivePatch& pp,
               const label& patchFaceI,
               vectorField& patchDisp,
               labelList& patchNLayers,
               List<extrudeMode>& extrudeStatus,
               PackedList<1>& isTriSplit,
               label& nChanged
            );

            //- Split faces at layer terminations
            void splitLayerTerminationFaces
            (
                const layerParameters& layerParams,
                const globalIndex& globalFaces,
                const labelListList& edgeGlobalFaces,
                const indirectPrimitivePatch& pp,
                const scalarField& minThickness,
                vectorField& patchDisp,
                labelList& patchNLayers,
                List<extrudeMode>& extrudeStatus,
                polyTopoChange& meshMod,
                PackedList<1>& isTriSplit
            );


            //- Remove points.
            static autoPtr<mapPolyMesh> doRemovePoints
            (
                meshRefinement& meshRefiner,
                removePoints& pointRemover,
                const boolList& pointCanBeDeleted
            );

            //- Restore faces (which contain removed points)
            static autoPtr<mapPolyMesh> doRestorePoints
            (
                meshRefinement& meshRefiner,
                removePoints& pointRemover,
                const labelList& facesToRestore
            );

            //- Return candidateFaces that are also in set.
            static labelList collectFaces
            (
                meshRefinement& meshRefiner,
                const labelList& candidateFaces,
                const labelHashSet& set
            );

            //- Pick up faces of cells of faces in set.
            static labelList growFaceCellFace
            (
                meshRefinement& meshRefiner,
                const labelHashSet&
            );

            //- Add cell layers
            void addLayers
            (
                layerParameters& layerParams,
                const dictionary& motionDict,
                const labelList& patchIDs,
                const label multiLayerIter
            );

            //- Add layers according to the dictionary settings
            bool doLayers
            (
                const dictionary& shrinkDict,
                const dictionary& motionDict,
                decompositionMethod& decomposer,
                fvMeshDistribute& distributor,
                const label multiLayerIter
            );

            static PackedBoolList getMasterPPEdges
            (
                const polyMesh& mesh,
                const indirectPrimitivePatch& pp
            );

            PackedBoolList getMasterPPPoints
            (
                const indirectPrimitivePatch& pp
            ) const;

            List<edge> calculateFaceCuts
            (
                const List<face>& updatedFaces
            ) const;

            void updateCutFaces
            (
                const List<face>& updatedFaces,
                const label masterFaceI,
                polyTopoChange& meshMod
            ) const;

            //Calculate planes for grown up patches
            List<plane> calculateGrownUpPlanes
            (
                const labelList& grownUpIDs
            ) const;

            //Project grown up patch points onto plane
            void projectToGrownUpPatches
            (
                const layerParameters& layerParams,
                const labelList& grownUpIDs,
                const List<plane>& grownUpPlanes
            ) const;

            //Calculate layer stretching based on fch and layer height
            scalar calculateStretching
            (
                const scalar initStr,
                const scalar a1,
                const scalar a2,
                const label nLay
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "snappyHexMeshDriver/snappyLayerDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
