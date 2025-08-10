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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::snappySnapDriver

Description
    All to do with snapping to surface

SourceFiles
    snappySnapDriver.C
    snappySnapDriverFeature.C

\*---------------------------------------------------------------------------*/

#ifndef snappySnapDriver_H
#define snappySnapDriver_H

#include "meshRefinement/meshRefinement.H"
#include "fields/Fields/DynamicField/DynamicField.H"
#include "patchDist/wallPoint/wallPoint.H"
#include "snappyHexMeshDriver/refinementParameters/refinementParameters.H"
#include "snappyHexMeshDriver/pointData/pointData.H"
#include "meshControl/meshControl.H"
#include "edgeClassification/edgeClassification.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class motionSmoother;
class refinementParameters;
class snapParameters;
class pointConstraint;
class featureEdgeMesh;
class featureLineSnapper;

/*---------------------------------------------------------------------------*\
                           Class snappySnapDriver Declaration
\*---------------------------------------------------------------------------*/

class snappySnapDriver
{
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

    // Select iteration count to use optimisation during feature snap
    class optStages
    {
        public:

        enum optType
        {
            NONE = 0,
            SMOOTH = 1,
            ERROR = 2,
            RELAXEDERROR = 4,
            BOUNDARYRELAXEDERROR = 8,
            EXTRUDE = 16
        };

        optStages() : event_() {}

        void add(const label index, const direction opt)
        {
            event_.append(Tuple2<label, direction>(index,opt));
        }

        direction find(const label iter) const
        {
            forAll(event_, eventi)
            {
                if (event_[eventi].first() == iter)
                {
                    return event_[eventi].second();
                }
            }

            return optStages::NONE;
        }

        private:

        DynamicList<Tuple2<label,direction>> event_;

    };

    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        fvMeshDistribute& distributor_;

        //- Reference to mesh dictionary
        const dictionary& meshDict_;

        //- Reference to mesh process controller
        const meshControl& controller_;

        //- From global surface region to master side patch
        const labelList globalToMasterPatch_;

        //- From global surface region to slave side patch
        const labelList globalToSlavePatch_;


    // Private Member Functions


        // Snapping

            //- Calculates (geometric) shared points
            //  Requires PackedList to be sized and initalised
            static label getCollocatedPoints
            (
                const scalar tol,
                const pointField&,
                PackedBoolList&
            );

            //- Calculate displacement (over all mesh points) to move points
            //  to average of connected cell centres
            static tmp<pointField> smoothInternalDisplacement
            (
                const meshRefinement& meshRefiner,
                const motionSmoother&
            );

            //- Calculate displacement per patch point to smooth out patch.
            //  Quite complicated in determining which points to move where.
            static tmp<pointField> smoothPatchDisplacement
            (
                const motionSmoother&,
                const List<labelPair>&,
                const bool preSmoothBaffles
            );

            static tmp<pointField> avg
            (
                const indirectPrimitivePatch&,
                const pointField&
            );

            //- Calculate displacement per patch point. Wip.
            static pointField smoothLambdaMuPatchDisplacement
            (
                const motionSmoother& meshMover,
                const List<labelPair>& baffles
            );


            //- Per edge distance to patch
            static tmp<scalarField> edgePatchDist
            (
                const pointMesh&,
                const indirectPrimitivePatch&
            );

            //- Write displacement as .obj file.
            static void dumpMove
            (
                const fileName&,
                const pointField&,
                const pointField&
            );

            //- Check displacement is outwards pointing
            static bool outwardsDisplacement
            (
                const indirectPrimitivePatch&,
                const vectorField&
            );

            //- Detect warpage
            void detectWarpedFaces
            (
                const scalar featureCos,
                const indirectPrimitivePatch& pp,

                DynamicList<label>& splitFaces,
                DynamicList<labelPair>& splits
            ) const;

            //- Get per face -1 or label of opposite face if on internal/baffle
            //  faceZone
            labelList getInternalOrBaffleDuplicateFace() const;

            //- Get points both on patch and facezone.
            static labelList getZoneSurfacePoints
            (
                const fvMesh& mesh,
                const indirectPrimitivePatch&,
                const word& zoneName
            );

            //- Get points both on patch and facezone.
            template<class FaceList>
            static labelList getFacePoints
            (
                const indirectPrimitivePatch& pp,
                const FaceList& faces
            );

            //- Per patch point calculate point on nearest surface.
            //  Return displacement of patch points.
            static void calcNearestSurface
            (
                const refinementSurfaces& surfaces,

                const labelList& surfacesToTest,
                const labelListList& regionsToTest,

                const pointField& localPoints,
                const labelList& zonePointIndices,

                scalarField& minSnapDist,
                labelList& snapSurf,
                vectorField& patchDisp,

                // Optional: nearest point, normal
                pointField& nearestPoint,
                vectorField& nearestNormal
            );

            void splitNonplanarQuads
            (
                const indirectPrimitivePatch &pp,
                const pointField &newPoints,
                const vectorField &surfHitPatchPointNormals
            ) const;

            void updateAfterMeshMods
            (
                const mapPolyMesh &map,
                const fvMesh &mesh,
                const labelList &adaptPatchIDs,
                indirectPrimitivePatch &pp,
                featureLineSnapper &flSnap
             );

            //Clean up mesh as the result of any topology changes
            void cleanupMesh(const dictionary& motionDict);

            //Correct convex layers cells by moving internal points
            void correctConvexLayerCells();

            // Check how well boundary mesh faces conform to the surfaces
            void checkSurfaceConformity
            (
                const snapParameters& snapParams,
                const indirectPrimitivePatch& pp
            );


            // Check how well boundary mesh faces conform to the surfaces
            void refineLayerCells();

            void splitDegenerateCells
            (
                const labelList &adaptPatchIDs,
                const dictionary& motionDict,
                indirectPrimitivePatch &pp
             );


            // Feature line snapping

                //- Is point on two feature edges that make a largish angle?
                bool isFeaturePoint
                (
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const PackedBoolList& isFeatureEdge,
                    const label pointi
                ) const;


                void smoothAndConstrain
                (
                    const PackedBoolList& isMasterEdge,
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const List<pointConstraint>& constraints,
                    vectorField& disp
                ) const;
                //void smoothAndConstrain2
                //(
                //    const bool applyConstraints,
                //    const indirectPrimitivePatch& pp,
                //    const List<pointConstraint>& constraints,
                //    vectorField& disp
                //) const;
                void calcNearest
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    vectorField& pointDisp,
                    vectorField& pointSurfaceNormal,
                    vectorField& pointRotation
                ) const;
                void calcNearestFace
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,
                    const scalarField& faceSnapDist,
                    vectorField& faceDisp,
                    vectorField& faceSurfaceNormal,
                    labelList& faceSurfaceRegion
                    //vectorField& faceRotation
                ) const;
                void calcNearestFacePointProperties
                (
                    const label iter,
                    const indirectPrimitivePatch& pp,

                    const vectorField& faceDisp,
                    const vectorField& faceSurfaceNormal,
                    const labelList& faceSurfaceRegion,

                    List<List<point>>& pointFaceSurfNormals,
                    List<List<point>>& pointFaceDisp,
                    List<List<point>>& pointFaceCentres,
                    List<labelList>&    pointFacePatchID
                ) const;

                void correctAttraction
                (
                    const DynamicList<point>& surfacePoints,
                    const DynamicList<label>& surfaceCounts,
                    const point& edgePt,
                    const vector& edgeNormal,   // normalised normal
                    const point& pt,
                    vector& edgeOffset  // offset from pt to point on edge
                ) const;


                //- For any reverse (so from feature back to mesh) attraction:
                //  add attraction if diagonal points on face attracted
                void stringFeatureEdges
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Remove constraints of points next to multi-patch points
                //  to give a bit more freedom of the mesh to conform to the
                //  multi-patch points. Bit dodgy for simple cases.
                void releasePointsNextToMultiPatch
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Detect any diagonal attraction. Returns indices in face
                //  or (-1, -1) if none
                labelPair findDiagonalAttraction
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& patchAttraction,
                    const List<pointConstraint>& patchConstraints,
                    const label facei
                ) const;

                scalar pyrVol
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& featureAttraction,
                    const face& localF,
                    const point& cc
                ) const;
                void facePoints
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& featureAttraction,
                    const vectorField& nearestAttraction,
                    const face& f,
                    DynamicField<point>& points
                ) const;
                scalar pyrVol
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& featureAttraction,
                    const vectorField& nearestAttraction,
                    const face& localF,
                    const point& cc
                ) const;
                Tuple2<point, vector> centreAndNormal
                (
                    const indirectPrimitivePatch& pp,
                    const vectorField& featureAttraction,
                    const vectorField& nearestAttraction,
                    const face& localF
                ) const;
                bool isSplitAlignedWithFeature
                (
                    const scalar featureCos,
                    const point& newPt0,
                    const pointConstraint& pc0,
                    const point& newPt1,
                    const pointConstraint& pc1
                ) const;
                bool isConcave
                (
                    const point& c0,
                    const vector& area0,
                    const point& c1,
                    const vector& area1,
                    const scalar concaveCos
                ) const;
                labelPair findDiagonalAttraction
                (
                    const scalar featureCos,
                    const scalar concaveCos,
                    const scalar minAreaFraction,
                    const indirectPrimitivePatch& pp,
                    const vectorField& patchAttraction,
                    const List<pointConstraint>& patchConstraints,
                    const vectorField& nearestAttraction,
                    const vectorField& nearestNormal,
                    const label faceI,

                    DynamicField<point>& points0,
                    DynamicField<point>& points1
                ) const;

                //- Do all logic on whether to add face cut to diagonal
                //  attraction
                void splitDiagonals
                (
                    const scalar featureCos,
                    const scalar concaveCos,
                    const scalar minAreaFraction,

                    const indirectPrimitivePatch& pp,
                    const vectorField& nearestAttraction,
                    const vectorField& nearestNormal,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints,
                    DynamicList<label>& splitFaces,
                    DynamicList<labelPair>& splits
                ) const;

                //- Avoid attraction across face diagonal since would
                //  cause face squeeze
                void avoidDiagonalAttraction
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Write some stats about constraints
                void writeStats
                (
                    const indirectPrimitivePatch& pp,
                    const PackedBoolList& isMasterPoint,
                    const List<pointConstraint>& patchConstraints
                ) const;

                //- Return hit if on multiple points
                pointIndexHit findMultiPatchPoint
                (
                    const point& pt,
                    const labelList& patchIDs,
                    const List<point>& faceCentres
                ) const;

                //- Return hit if faces-on-the-same-normalplane are on multiple
                //  patches
                pointIndexHit findMultiPatchPoint
                (
                    const point& pt,
                    const labelList& pfPatchID,
                    const DynamicList<vector>& surfaceNormals,
                    const labelList& faceToNormalBin
                ) const;

                //- Return index of similar normal
                label findNormal
                (
                    const scalar featureCos,
                    const vector& faceSurfaceNormal,
                    const DynamicList<vector>& surfaceNormals
                ) const;

                //- Determine attraction and constraints for single point
                //  using sampled surrounding of the point
                void featureAttractionUsingReconstruction
                (
                    const label iter,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    const label pointi,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    DynamicList<point>& surfacePoints,
                    DynamicList<vector>& surfaceNormals,
                    labelList& faceToNormalBin,

                    vector& patchAttraction,
                    pointConstraint& patchConstraint
                ) const;

                //- Determine attraction and constraints for all points
                //  using sampled surrounding of the point
                void featureAttractionUsingReconstruction
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Determine geometric features and attraction to equivalent
                //  surface features
                void determineFeatures
                (
                    const label iter,
                    const scalar featureCos,
                    const bool multiRegionFeatureSnap,

                    const indirectPrimitivePatch&,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Determine features originating from bafles and
                //  and add attraction to equivalent surface features
                void determineBaffleFeatures
                (
                    const label iter,
                    const bool baffleFeaturePoints,
                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    // Feature-point to pp point
                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;
                void reverseAttractMeshPoints
                (
                    const label iter,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,

                    // Feature-point to pp point
                    const List<labelList>& pointAttractor,
                    const List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    const List<List<DynamicList<point>>>& edgeAttractors,
                    const List<List<DynamicList<pointConstraint>>>&,

                    const vectorField& rawPatchAttraction,
                    const List<pointConstraint>& rawPatchConstraints,

                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Find point on nearest feature edge (within searchDist).
                //  Return point and feature
                //  and store feature-edge to mesh-point and vice versa
                Tuple2<label, pointIndexHit> findNearFeatureEdge
                (
                    const bool isRegionEdge,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const label pointi,
                    const point& estimatedPt,

                    List<List<DynamicList<point>>>&,
                    List<List<DynamicList<pointConstraint>>>&,
                    vectorField&,
                    List<pointConstraint>&
                ) const;

                //- Find nearest feature point (within searchDist).
                //  Return feature point
                //  and store feature-point to mesh-point and vice versa.
                //  If another mesh point already referring to this feature
                //  point and further away, reset that one to a near feature
                //  edge (using findNearFeatureEdge above)
                Tuple2<label, pointIndexHit> findNearFeaturePoint
                (
                    const bool isRegionEdge,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const label pointi,
                    const point& estimatedPt,

                    // Feature-point to pp point
                    List<labelList>& pointAttractor,
                    List<List<pointConstraint>>& pointConstraints,
                    // Feature-edge to pp point
                    List<List<DynamicList<point>>>& edgeAttractors,
                    List<List<DynamicList<pointConstraint>>>& edgeConstraints,
                    // pp point to nearest feature
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                void featureAttractionUsingFeatureEdges
                (
                    const label iter,
                    const bool multiRegionFeatureSnap,

                    const bool detectBaffles,
                    const bool baffleFeaturePoints,
                    const bool releasePoints,
                    const bool stringFeatures,
                    const bool avoidDiagonal,

                    const scalar featureCos,

                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    const vectorField& nearestNormal,

                    const List<List<point>>& pointFaceSurfNormals,
                    const List<List<point>>& pointFaceDisp,
                    const List<List<point>>& pointFaceCentres,
                    const labelListList& pointFacePatchID,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                void preventFaceSqueeze
                (
                    const label iter,
                    const scalar featureCos,
                    const indirectPrimitivePatch& pp,
                    const scalarField& snapDist,
                    const vectorField& nearestAttraction,

                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints
                ) const;

                //- Top level feature attraction routine. Gets given
                //  displacement to nearest surface in nearestDisp
                //  and calculates new displacement taking into account
                //  features
                vectorField calcNearestSurfaceFeature
                (
                    const snapParameters& snapParams,
                    const bool alignMeshEdges,
                    const label iter,
                    const scalar featureCos,
                    const scalar featureAttract,
                    const scalarField& snapDist,
                    const vectorField& nearestDisp,
                    const vectorField& nearestNormal,
                    motionSmoother& meshMover,
                    vectorField& patchAttraction,
                    List<pointConstraint>& patchConstraints,

                    DynamicList<label>& splitFaces,
                    DynamicList<labelPair>& splits
                ) const;


        //- Disallow default bitwise copy construct
        snappySnapDriver(const snappySnapDriver&);

        //- Disallow default bitwise assignment
        void operator=(const snappySnapDriver&);


public:

    //- Runtime type information
    ClassName("snappySnapDriver");


    // Constructors

        //- Construct from components
        snappySnapDriver
        (
            meshRefinement& meshRefiner,
            decompositionMethod& decomposer,
            fvMeshDistribute& distributor,
            const dictionary& meshDict,
            const meshControl& controller,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch
        );


    // Member Functions

        // Snapping

            //Identify which faces are zoned (internal or baffle)
            static PackedList<1> setZonedFaces(const meshRefinement& meshRefiner);

            //Identify which faces are zoned excluding boundary zoned
            static PackedList<1> setZonedFacesExclBoundaryNamed
            (
                const meshRefinement& meshRefiner,
                const refinementParameters& refineParams
            );

            //Identify which faces are zoned and disbale baffle checks on
            static PackedList<1> setZonedBaffleCheckFaces
            (
                const meshRefinement& meshRefiner,
                const refinementParameters& refineParams
            );

            void resetPrimitivePatchAddressing
            (
                const labelList& adaptPatchIDs,
                indirectPrimitivePatch &pp
            );

            //- Merge baffles.
            autoPtr<mapPolyMesh> mergeZoneBaffles(const List<labelPair>&);

            //- Calculate edge length per patch point.
            static scalarField calcSnapDistance
            (
                const meshRefinement& meshRefiner,
                const snapParameters& snapParams,
                const indirectPrimitivePatch&
            );

            ////- Get patches generated for surfaces.
            //labelList getSurfacePatches() const;

            //- Smooth the mesh (patch and internal) to increase visibility
            //  of surface points (on castellated mesh) w.r.t. surface.
            static void preSmoothPatch
            (
                const meshRefinement& meshRefiner,
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );

            //- Smooth in vicinity of points not converged during feature snap
            void smoothNonConvergent
            (
                const indirectPrimitivePatch& pp,
                const boolList convergedPoints,
                pointField& newPatchPoints
            ) const;

            //- Helper: calculate average cell centre per point
            static tmp<pointField> avgCellCentres
            (
                const fvMesh& mesh,
                const indirectPrimitivePatch&
            );

            //- Per patch point override displacement if in gap situation
            void detectNearSurfaces
            (
                const scalar planarCos,
                const indirectPrimitivePatch&,
                const pointField& nearestPoint,
                const vectorField& nearestNormal,
                vectorField& disp
            ) const;

            //- Per patch point calculate point on nearest surface. Set as
            //  boundary conditions of motionSmoother displacement field. Return
            //  displacement of patch points.
            static vectorField calcNearestSurface
            (
                const bool strictRegionSnap,
                const meshRefinement& meshRefiner,
                const labelList& globalToMasterPatch,
                const labelList& globalToSlavePatch,
                const scalarField& snapDist,
                const indirectPrimitivePatch&,
                pointField& nearestPoint,
                vectorField& nearestNormal
            );

            //- Per patch point calculate point on nearest surface. Set as
            //  boundary conditions of motionSmoother displacement field. Return
            //  displacement of patch points.
            vectorField calcNearestSurface
            (
                const scalarField& snapDist,
                const List<labelPair>& baffles,
                motionSmoother& meshMover
            );

            //- Smooth the displacement field to the internal.
            void smoothDisplacement
            (
                const snapParameters& snapParams,
                motionSmoother&
            ) const;

            //- Interogate surface to find surface normals
            void findSurfaceNormals
            (
                const labelList& checkSurfaces,
                const indirectPrimitivePatch& pp,
                const scalarField& snapDist,
                const pointField& newPoints,
                const pointField& ppFaceCentres,
                const pointField& ppCellCentres,
                const vectorField& ppFaceAreas,
                const boolList& stationaryFaces,
                const bool& fastFeatureSnap,
                List<pointIndexHit>& finalHitInfo,
                labelList& pointIndex,
                vectorField& pNormals
            ) const;

            //- Perform feature snap initialisation
            void featureSnapInit
            (
                const bool initialSnap,
                pointField& newPoints,
                vectorField& pointDisp,
                scalarField& sumAreasMesh,
                vectorField& sumDispMesh,
                boolList& gapBoundaryPts
            ) const;

            //- Perform optimisation stages during feature snap
            void featureOptimise
            (
                const label iter,
                const optStages& ostages,
                const snapParameters& snapParams,
                const dictionary& motionDict,
                const indirectPrimitivePatch& pp,
                const labelList& snapSurfaces,
                const labelList& meshEdges,
                const boolList& boundaryPoints,
                const boolList& gapBoundaryPts,
                pointField& newPoints,
                motionSmoother& meshMover,
                boolList& stationaryPoints
            ) const;

            //- Recalculate surface point normal to snap to
            // based on an extended stencil
            void recalculatePointNormals
            (
                const snapParameters& snapParams,
                const indirectPrimitivePatch& pp,
                const labelList& snapSurfaces,
                const scalarField& snapDist,
                const pointField& newPoints,
                const List<pointIndexHit>& finalHitInfo,
                const labelList& pointIndex,
                const pointField& ppFaceCentres,
                const vectorField& ppFaceAreas,
                vectorField& pNormals
            ) const;

            //- Acute/reflex and weak feature snap
            void tightAngleCorrect
            (
                const bool weakSnap,
                const scalar weakAngleMin,
                const scalar weakAngleMax,
                const scalar tightAngle,
                const scalar minTightAngle,
                const pointField& newPoints,
                const indirectPrimitivePatch& pp,
                const labelList& meshEdges,
                const List<Tuple2<edgeClassification::edgeType,scalar>>& eType,
                const List<pointIndexHit>& finalHitInfo,
                const labelList& pointIndex,
                const vectorField& ppFaceAreas,
                const vectorField& pNormals,
                scalarField& sumAreasPP,
                vectorField& sumDispPP
            ) const;

            //- Perform surface smoothing during feature snap
            void featureSnapSurfaceSmooth
            (
                const label iter,
                const label nSurfSmooth,
                const indirectPrimitivePatch& pp,
                const List<Tuple2<edgeClassification::edgeType,scalar>>& eType,
                const PackedBoolList& isMasterEdge,
                const PackedBoolList& isMasterFace,
                const PackedBoolList& isMasterPPEdge,
                const scalarField& snapDist,
                const boolList& boundaryPoints,
                const scalarField& smoothingWeights,
                const boolList& excludeFaces,
                const boolList& zonedEdges,
                pointField& newPoints,
                vectorField& pointDisp,
                pointField& ppFaceCentres,
                vectorField& ppFaceAreas,
                boolList& stationaryPoints,
                scalarField& sumAreasMesh,
                vectorField& sumDispMesh,
                scalarField& sumAreasPP,
                vectorField& sumDispPP
            ) const;

            //- Snap faces to surface, trying to preserve surface features
            void featureSnap
            (
                const label maxIter,
                const label smoothStartIter,
                const labelList& snapSurfaces,
                const indirectPrimitivePatch& pp,
                const snapParameters& snapParams,
                const dictionary& motionDict,
                motionSmoother& meshMover,
                bool initialSnap,
                const optStages& ostages
            ) const;

            //Perform additional corrections during snap for extrude method
            void extrudeCorrect
            (
                const indirectPrimitivePatch& pp,
                const dictionary& motionDict,
                const boolList& boundaryPoints,
                motionSmoother& meshMover,
                pointField& newPoints
            ) const;

            //Perform additional corrections during snap for corner points
            // when using the extrude method
            void cornerCorrect
            (
                const indirectPrimitivePatch& pp,
                const labelList& meshEdges,
                const boolList& boundaryPoints,
                pointField& newPoints
            ) const;

            //- Do the hard work: move the mesh according to displacement,
            //  locally relax the displacement. Return true if ended up with
            //  correct mesh, false if not.
            bool scaleMesh
            (
                const snapParameters& snapParams,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );

            autoPtr<mapPolyMesh> surfaceToPatch
            (
                const snapParameters& snapParams,
                const refinementParameters& refineParams
            );

            //- Update baffle labels
            static void updateBaffles
            (
                const mapPolyMesh &map,
                List<labelPair>& baffles
            );

            //- Merge patch faces based on surface normals
            void mergePatchFacesUndo
            (
                const dictionary& motionDict,
                const scalarField& snapDist,
                indirectPrimitivePatch &pp,
                List<labelPair>& baffles,
                bool calcNormals = true
            );

            //- Final merge of patch faces
            void mergePatchFacesUndo
            (
                const dictionary& motionDict,
                const List<labelList> &featureShadows,
                const bool maintainPatches
            );

            //- Try splitting wraped/region cells
            // Used by extrude and dual mesher
            void splitCells
            (
                const List<labelList>& shadowChains,
                const refinementParameters& refineParams,
                const snapParameters& snapParams,
                const labelList& adaptPatchIDs,
                const PackedList<1>& isZonedFace,
                indirectPrimitivePatch &pp
            );

            //- Try smoothing sliver faces generated during snapping
            void smoothSliverFaces
            (
                const dictionary& motionDict,
                const snapParameters& snapParams,
                const label nInitErrors,
                const labelList& adaptPatchIDs,
                const List<labelPair>& baffles,
                indirectPrimitivePatch &pp
            );

            //- Try smoothing non-conformal faces generated during snapping
            void smoothNonConformalFaces
            (
                const dictionary& motionDict,
                const snapParameters& snapParams,
                const label nInitErrors,
                const labelList& adaptPatchIDs,
                const List<labelPair>& baffles,
                indirectPrimitivePatch &pp
            );

            //- Try and collapse small edges
            void collapseSmallEdges
            (
                const dictionary& motionDict,
                const labelList& adaptPatchIDs,
                const snapParameters& snapParams,
                indirectPrimitivePatch &pp
            );

            //- Try and correct concave faces
            void correctConcaveFaces
            (
                const dictionary& motionDict,
                const labelList& adaptPatchIDs,
                const snapParameters& snapParams,
                indirectPrimitivePatch &pp
            );

            //Add additional tet cells that may not have been added during splitting
            void addTetrahedralCellsToSplitMesh
            (
                const labelList& adaptPatchIDs
            );

            //Return list of zone patches (INTERNAL and BAFFLE)
            labelList zonePatchIDs() const;

            //Return list of zone slave patches (INTERNAL and BAFFLE)
            labelList zonePatchSlaveIDs() const;

            //Remove cells with no added layers and re-snap (extrude method)
            void removeLayersExtrude
            (
                const dictionary& motionDict,
                const snapParameters& snapParams,
                const label nInitErrors,
                const dictionary& layerDict,
                const labelList& adaptPatchIDs,
                indirectPrimitivePatch &pp
            );

            //Remove cells with no added layers and re-snap (dual method)
            void removeLayersDual
            (
                const dictionary& motionDict,
                const snapParameters& snapParams,
                const label nInitErrors,
                const dictionary& layerDict,
                const labelList& adaptPatchIDs,
                const labelList& featureIDs,
                indirectPrimitivePatch &pp
            );

            void snapAndMerge
            (
                snapParameters& snapParams
            );

            // Timing test per surface to find nearest point
            void testSnapPerformance
            (
                const indirectPrimitivePatch &pp,
                const snapParameters& snapParams
            );

            void doSnap
            (
                snapParameters& snapParams,
                const refinementParameters& refineParams
            );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "snappyHexMeshDriver/snappySnapDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
