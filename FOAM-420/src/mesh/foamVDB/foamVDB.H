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
    (c) 2020 ESI

Class
    Foam::foamVDB

Description

SourceFiles
    foamVDBI.H
    foamVDB.C
    foamVDBIO.C

\*---------------------------------------------------------------------------*/

#ifndef foamVDB_H
#define foamVDB_H

#include "include/OpenVDB.H"

//#include "global/clock/clock.H"
//#include <chrono>
//#include "cpuTime/cpuTime.H"

#include <vector>
#include "triSurface/triSurface.H"
#include "meshes/primitiveShapes/point/pointList.H"
#include "fvMesh/fvMesh.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchField.H"
#include "primitives/Vector/labelVector/labelVector.H"
#include "primitives/ints/lists/labelList.H"

#include "include/TBBTimer.H"

namespace Foam
{
    //struct Timer
    //{
    //    std::chrono::time_point<std::chrono::high_resolution_clock> start_, end_;
    //    std::chrono::duration<float> duration_;
    //    //clock timer_;
    //    cpuTime cpuTimer_; //https://pubs.opengroup.org/onlinepubs/9699919799/functions/times.html
    //    word text_;
    //    label threads_;

    //    Timer(word text="")
    //    {
    //        start_ = std::chrono::high_resolution_clock::now();
    //        //timer_ = clock();
    //        cpuTimer_ = cpuTime();
    //        text_ = text;
    //        threads_ = tbb::this_task_arena::max_concurrency();
    //    }

    //    ~Timer()
    //    {
    //        end_ = std::chrono::high_resolution_clock::now();
    //        duration_ = end_ - start_;

    //        scalar cpuTime = cpuTimer_.elapsedCpuTime();

    //        Info<< "=== Timer: "
    //            << setw(38) << text_;

    //        // print time of each processor
    //        if (Pstream::parRun())
    //        {
    //            List<scalar>   procCpuTime(Pstream::nProcs());
    //            List<scalar> procClockTime(Pstream::nProcs());

    //            procCpuTime[Pstream::myProcNo()]   = cpuTime;
    //            procClockTime[Pstream::myProcNo()] = duration_.count();

    //            IPstream::gatherList(procCpuTime);
    //            IPstream::gatherList(procClockTime);

    //            if (Pstream::master())
    //            {
    //                List<scalar> procEfficiency
    //                (
    //                    100.0 * (procCpuTime / threads_) /  procClockTime
    //                );

    //                Info<< " == clockTime (s) avg " << average(procClockTime) << setprecision(3)
    //                    << " - min " << min(procClockTime) << setprecision(3)
    //                    << " - max " << max(procClockTime) << setprecision(3) << nl
    //                    //<< "    cpuTime (s)    " << flatOutput(procCpuTime) << setprecision(3) << nl
    //                    //<< "    clockTime (s)  " << flatOutput(procClockTime) << setprecision(3) << nl
    //                    //<< "    efficiency (%) " << flatOutput(procEfficiency) << setprecision(3)
    //                    << "    efficiency (%) avg " << average(procEfficiency) << setprecision(3)
    //                    << " - min " << min(procEfficiency) << setprecision(3)
    //                    << " - max " << max(procEfficiency) << setprecision(3)
    //                    << endl;
    //            }
    //        }
    //        else
    //        {
    //            scalar eff = 100.0 * (cpuTime / threads_) /  duration_.count();

    //            Info<< " ====== cpuTime " << setw(5) << cpuTime << " s "
    //                //<< "===== clockTime " << setw(3) << timer_.elapsedClockTime() << " s "
    //                << "==== clockTime " << setw(9) << duration_.count() << setprecision(3) << " s "
    //                << "==== efficiency " << setw(4) << eff << setprecision(3) << "% "
    //                << "on " << threads_ << " threads" << endl;
    //        }
    //    }
    //}; //struct Timer

class meshRefinement;
class refinementSurfaces;
class refinementParameters;

struct vdbPoint
{
    point ijk_;
    label id_;
    bool active_;

    vdbPoint
    (
        const point& ijk,
        const label id
    )
    :
        ijk_(ijk),
        id_(id),
        active_(true)
    {}
};

struct vdbFace
{
    labelList vdbPoints_;
    label owner_;
    label neighbour_; //or patchID for boundary faces
    label cellLevel_;
    label id_;
    bool active_;

    vdbFace
    (
        const labelList& points,
        const label own,
        const label neigh,
        const label level,
        const label id
    )
    :
        vdbPoints_(points),
        owner_(own),
        neighbour_(neigh),
        cellLevel_(level),
        id_(id),
        active_(true)
    {}
};

// Coordinates of fine voxels inside a coarse one
extern const openvdb::Coord COARSE_TO_FINE[8];

// Voxel-face adjacent neighbours
extern const openvdb::Coord COORD_OFFSETS[6];

// edge-face connectivity in COORD_OFFSETS indices
extern const List<labelList> edgeFace;

// Voxel-face adjacent fine neighbours
extern const openvdb::Coord COORD_OFFSETS_SPLIT[24];

// edge-fine voxels connectivity in COORD_OFFSET_SPLITS indices
extern const List<labelList> fineEdgeFace;

// Voxel-edge adjacent neighbours
extern const openvdb::Coord EDGE_COORD_OFFSETS[12];

// Voxel-edge adjacent fine neighbours
extern const openvdb::Coord FINE_EDGE_COORD_OFFSETS[24];

extern const List<face> localFaces;
extern const List<word> faceNames;
extern const List<edge> localEdges;

// edges expressed as connection between start point
// of edges with lower index
// e.g. localEdge[8] = edge{localEdge[0].start(), localEdge[4].start()};
extern const List<labelList> edgeEdges;

// faces as edge loop
extern const List<labelList> edgeLoops;

/*---------------------------------------------------------------------------*\
                         Class foamVDB Declaration
\*---------------------------------------------------------------------------*/
class foamVDB
{

protected:

    // Protected data

        //const dictionary& openvdbDict_;

        const scalar voxelSize_;

        const label maxCellLevel_;

private:

    // Private data

    // Private Member Functions

        //- calculate voxel-voxel connectivity
        labelListList calcVoxelVoxels
        (
            std::vector<IndexGrid::Ptr>&  globalIDGrids,
            const label nVoxels
        );

        //- refine voxels which have more than 3 split faces
        void refineDanglingVoxels
        (
            std::vector<FloatGrid::Ptr>&  cellLevelGrids
        );

        //- sync grids across processors
        // topology difference of overlapping voxels
        // returns max cellID of sender processor
        label syncGrids
        (
            std::vector<FloatGrid::Ptr>&  cellLevelGrids,
            std::vector<IndexGrid::Ptr>&  globalIDGrids,
            const List<List<boundBox>>& procBBoxes,
            const label procN
        );

        //- Count voxelID and return ID offset for each cellLevel
        labelList countVoxels
        (
            const std::vector<FloatGrid::Ptr>&  cellLevelGrids,
            std::vector<IndexGrid::Ptr>&  globalIDGrids,
            pointField& voxelCentres
        );

        //- Count voxelID and return ID offset for each cellLevel
        // uses LeafManager
        labelList countVoxelsLeaf
        (
            const std::vector<FloatGrid::Ptr>&  cellLevelGrids,
            std::vector<IndexGrid::Ptr>&  globalIDGrids,
            pointField& voxelCentres
        );

        //- Returns vertices of a voxel
        const std::vector<openvdb::Coord> getVoxelVertices
        (
            const label pointsPerEdge,
            const openvdb::Coord& coord
        );

        //- Append to global list of points
        void addPoint
        (
            openvdb::Coord p,
            IndexGrid::Accessor& pointIDAcc,
            IndexGrid::Accessor& pointLevelAcc,
            std::vector<openvdb::Coord>& edgeCoords,
            labelList& edge,
            openvdb::Int32& nPoints,
            const label cellLevel,
            std::vector<point>& xyzPoints,
            std::vector<vdbPoint>& vdbPoints
        );

        //- Add mid point to edge and append to global list of points
        void addMidPoint
        (
            openvdb::Coord p,
            IndexGrid::Accessor& pointIDAcc,
            IndexGrid::Accessor& pointLevelAcc,
            std::vector<openvdb::Coord>& edgeCoords,
            labelList& edge,
            openvdb::Int32& nPoints,
            const label cellLevel,
            std::vector<point>& xyzPoints,
            std::vector<vdbPoint>& vdbPoints,
            openvdb::Int32& pointID
        );

        //- Count active points/faces in given list of vdbPoint/vdbFace
        template<class T>
        label activeCount
        (
            std::vector<T>& list
        );

        //- Add displacement to voxel face vertices
        // (snap to surface)
        // return true if face area is positive
        bool addDisplacement
        (
            label cellLevel,
            const face& f,
            std::vector<openvdb::Coord>& faceCoords,
            labelHashSet& displacedPoints,
            std::vector<point>& xyzPoints,
            std::vector<vdbPoint>& vdbPoints,
            Vec3dGrid::ConstAccessor& displAcc,
            Vec3dGrid::ConstAccessor& fineDisplAcc,
            const bool isSplitFace = false
        );

        //- Add empty patch to polyMesh
        void addEmptyPatch
        (
            fvMesh& mesh,
            const word& patchName,
            const dictionary& patchInfo
        );

        //- Group regions which have same minMax level and same
        // curvatureRefineAngle
        labelListList groupSameLevelRegions
        (
            const wordList& regNames,
            const refinementSurfaces& surfaces,
            const scalarField& featAngles
        );

        //- Counts active voxels inside a CoordBBox
        openvdb::Index64 onVoxelsInCoordBbox
        (
            FloatGrid::Ptr grid,
            const openvdb::CoordBBox& bBox
        );

        //- slice grids in one direction for hierarchical decomposition
        std::vector<FloatGrid::Ptr> sliceGrids
        (
            const std::vector<FloatGrid::Ptr>& grids,
            openvdb::math::BBox<openvdb::Vec3s>& bBoxSlice,
            const label& proci,
            const Vector<Foam::label>& n,
            const label& dir,
            const label& desiredVoxels,
            const scalar& maxLoadUnbalance
        );

        //- Recursive method to do a balanced hierarchical decomposition along 3 axis
        void doDecomposition
        (
            const std::vector<FloatGrid::Ptr>& grids,
            std::vector<std::vector<FloatGrid::Ptr>>& decomposedGrids,
            const Vector<Foam::label>& n,
            FixedList<direction, 3>& decompOrder,
            const label& nthDir,
            label procOffset,
            const label& nProcs,
            const scalar& maxLoadUnbalance
        );

        //- Calls the doDecomposition along slice in nthDir for first 2 directions
        //- if 3rd direction assigns grids to processors
        void fillProcGridsOrDecompose
        (
            const std::vector<FloatGrid::Ptr>& grids,
            std::vector<std::vector<FloatGrid::Ptr>>& decomposedGrids,
            const Vector<Foam::label>& n,
            FixedList<direction, 3>& decompOrder,
            const label nthDir,
            const label procI,
            const label offset,
            const label nProcs,
            const scalar maxLoadUnbalance
        );

        //- Disallow default bitwise copy construct and assignemnt
        foamVDB(const foamVDB&);
        void operator=(const foamVDB&);

public:

    // Constructors

        //- Construct from minimum voxelSize and max cell level
        foamVDB(const scalar& voxelSize, const label& maxCellLevel);


    //- Destructor
    ~foamVDB()
    {}

    // Member Functions

        //- autocreate decompositionDict for VDB grids
        dictionary createDecompDict
        (
            const Time& runTime,
            const word method
        );

        //- Check wanted number of subdomains and create proc folders if
        // they do not exist, return decomposeParDict and nDomainsPerNode
        dictionary checkDecomposeParDict
        (
            const Time& runTime,
            labelList& nDomainsPerNode,
            label& procOffset
        );

        //- Voxelize surfaces and refinement regions
        void voxelize
        (
            meshRefinement& meshRefiner,
            const refinementParameters& refineParams,
            const dictionary& vdbDict
        );

        //- Calculates face addressing and returns OpenFOAM mesh data
        void vdbGridsToPolyMesh
        (
            const std::vector<FloatGrid::Ptr>& cellLevelGrids,
            const std::vector<Vec3dGrid::Ptr>& displacementGrids,
            const std::vector<BoolGrid::Ptr>& innerGrids,
            const std::vector<BoolGrid::Ptr>& outerGrids,
            const label& nRegions,
            const point& levelSetOffset,
            const point& scale,
            const quaternion& R,
            const openvdb::CoordBBox& boundingBox,
            const dictionary& vdbDict,

            fvMesh&    mesh,
            labelList& hexRef8cellLevel,
            labelList& hexRef8pointLevel
        );

        //- Distribute triangles weighted with refinement level
        triSurface distributeSurface
        (
            const triSurface& globalSurf,
            const labelList& refinementLevels
        );

        //- Build a zero area triSurface out of feature edges
        //- to use as input for an Unsigned Distance Field
        triSurface featureTriSurface
        (
            const triSurface& surf,
            const scalar featAngle,
            const word name
        );

        //- Offset levelset generated by MultiResGrid starting from a fine grid
        //- by means of topDownRestrict
        void topDownAdvect
        (
            FloatGrid::Ptr grid,
            const label cellLevel,
            label maxLevel = -1
        );

        //- Offset levelset generated by MultiResGrid starting from a coarse grid
        void bottomUpAdvect
        (
            FloatGrid::Ptr grid,
            const label cellLevel
        );

        //- Advect a narrow band level set by delta
        void advectLevelSet
        (
            FloatGrid::Ptr grid,
            const scalar delta,
            const bool fromFineGrid = false
        );

        //- Returns coord from global index space to local index
        static const openvdb::Coord toGlobal
        (
            const openvdb::Coord& xyzLocal,
            openvdb::math::Transform::ConstPtr transform
        );

        //- Returns coord from local index space to global index
        static const openvdb::Coord toLocal
        (
            const openvdb::Coord& xyzGlobal,
            openvdb::math::Transform::ConstPtr transform
        );

        //- Splits one face of level N to 4 faces of level N+1
        void splitVoxelFace
        (
            const pointList& points,
            std::atomic<label>& pointCount,
            const openvdb::math::Vec4ui& f,
            const openvdb::Coord& coarseFaceCoord,
            openvdb::math::Transform::ConstPtr fineTransform,
            openvdb::math::Transform::ConstPtr coarseTransform,
            IndexGrid::Accessor& fineIndexGridAccessor,
            IndexGrid::Accessor& fineAddedIndexesGridAccessor,
            BoolGrid::ConstAccessor& cellLevelInterfaceFaceGridAccessor,
            IndexGrid::ConstAccessor& coarseOwnerGridAccessor,
            IndexGrid::Accessor& fineOwnerGridAccessor,
            IndexGrid::ConstAccessor& coarsePatchIDGridAccessor,
            IndexGrid::Accessor& finePatchIDGridAccessor,
            FaceGrid::Accessor& fineFaceGridAccessor,
            bool forceAdd
        );

        ////GGG no need to do this, as patchID comes from polygonIDGrids
        //// assign boundary voxels to closest region
        ////- Assign voxels to the closest patch
        //void checkPatchVoxelsProximity
        //(
        //    const std::vector<std::vector<FloatGrid::Ptr>>& boundaryCellLevels
        //);

        //- Slice boundBox along longest axis (returned) by the number of processors
        void sliceBoundBox
        (
            treeBoundBox& bb,
            const label coarserCellLevel,
            vector& axis
        );

        //- Converts cuttingPlane to triSurface
        const triSurface getVDBCuttingPlane
        (
            const edgeList& edges,
            const scalarField& dotProducts,
            const pointField& points,
            const plane& cutPlane
        );

        //- Closes a triSurface with bounding box coming
        //  from distributedTriSurface decomposition
        void clipClosedSurface
        (
            triSurface& surf,
            treeBoundBox& surfBb,
            const label coarserCellLevel
        );

        void decomposeVDBSurface
        (
            triSurface& surf,
            const label coarserCellLevel
        );

        //- Returns triangulation of 6 planes of a bounding box
        //const List<triSurface> boundBoxToPlanes
        const List<plane> boundBoxToPlanes
        (
            const treeBoundBox& box
        );

        //- Triangulates a bounding box
        const triSurface boundBoxToTriSurface
        (
            treeBoundBox& box,
            const label shellLevel = -1,
            const bool decompose = false
        );

        //- Splits a triSurface by patch or by same refinement level
        const List<triSurface> splitSurface
        (
            const triSurface& surf,
            List<labelList> sameMinMaxLevel = List<labelList>(0)
        );

        //- Add surfaces in List to one triSurface
        // all the surfs have the same patches with same ordering
        // similar to Foam::triSurfaceLoader::load with loadingOption::OFFSET_REGION
        const triSurface addSurfaces
        (
            List<triSurface>& surfs
        );

        //- Gather triSurfaces from all processors and combine them in a global triSurface
        const triSurface gatherTriSurface
        (
            triSurface& inputSurf
        );

        //- Counts active voxels in all grids
        template<typename GridType>
        label getActiveVoxels
        (
            const std::vector<typename GridType::Ptr>& grids
        );

        //- Finds bounding box of all grids expressed in level0 index space
        openvdb::math::BBox<openvdb::Vec3s> getBoundingBox
        (
            const std::vector<FloatGrid::Ptr>& grids
        );

        //- Decompose using OpenFOAM decomposer. Return per cellLevel nVoxels going to which proc
        labelListList decomposeGrids
        (
            const std::vector<FloatGrid::Ptr>& cellLevelGrids,
            std::vector<IndexGrid::Ptr>&  globalIDGrids,
            std::vector<std::vector<FloatGrid::Ptr>>& decomposedGrids,
            const dictionary& decompositionDict,
            const labelList& nVoxelsStart,
            const pointField& voxelCentres
        );

        //- Decompose all cellLevel grids with specified method
        // balancedCellLevels: every proc has the same number of voxels for each cellLevel
        // globalBalance:      every processor has the same number of voxels but
        //                         not necessarly same voxels for each cellLevel
        // hierarchical:
        void decomposeGrids
        (
            const std::vector<FloatGrid::Ptr>& grids,
            std::vector<std::vector<FloatGrid::Ptr>>& decomposedGrids,
            const dictionary& decompositionDict,
            const label& nProcs,
            const scalar& maxLoadUnbalance
        );

        //- Convert Foam::List<T> to OpenVDB vector (multithreaded)
        template<class Type, class OpenVDBType>
        void foamToVDB
        (
            const List<Type>& inputList,
            std::vector<OpenVDBType>& outputList
        );

        //- Create a narrow band level set from a closed triangulated surface
        //  or a unsigned distance field from an open triangulated surface
        FloatGrid::Ptr meshToLevelSet
        (
            const triSurface& inputSurf,
            const label cellLevel,
            const scalar halfWidth = 3.0,
            const int conversionFlags = 0, //i.e. Signed Distance Field
            typename FloatGrid::template ValueConverter<openvdb::Int32>::Type * polygonIndexGrid = nullptr,
            typename FloatGrid::template ValueConverter<bool>::Type * maskGrid = nullptr,
            typename FloatGrid::template ValueConverter<openvdb::Vec4i>::Type * polygonListGrid = nullptr,
            typename FloatGrid::template ValueConverter<openvdb::Vec3d>::Type * displacementGrid = nullptr
        );

        //- Offset the zero level set. This is needed when levelSetRebuild
        //  operates on Unsigned Distance Field
        //void offsetLevelSetValues
        //(
        //    const FloatGrid::Ptr& grid,
        //    const scalar& offset
        //);

        //- Checks if the given grid is a narrow-band Signed Distance Field
        bool isOpenSurface
        (
            const FloatGrid::Ptr& grid
        );

        //- Resample the input grid on the output grid index space
        void resample
        (
            const BoolGrid& inGrid,
            BoolGrid& outGrid
        );

        //- Checks if it is a valid pointer and the grid has active voxels
        template<typename GridType>
        bool isValid
        (
            typename GridType::ConstPtr grid
        );

        //- Send vector of grids to neighbour node
        template<typename GridType>
        static void sendGrids
        (
            const std::vector<typename GridType::Ptr>& inputGrids,
            PstreamBuffers& pBufSize,
            PstreamBuffers& pBufGrid,
            const label node
        );

        //- Send grid to neighbour node
        template<typename GridType>
        static void sendGrid
        (
            const typename GridType::Ptr grid,
            PstreamBuffers& pBufSize,
            PstreamBuffers& pBufGrid,
            const label node
        );

        //- Receive vector of grids from neighbour node
        template<typename GridType>
        static std::vector<typename GridType::Ptr> receiveGrids
        (
            PstreamBuffers& pBufSize,
            PstreamBuffers& pBufGrid,
            const label node
        );

        //- Receive grid from neighbour node
        template<typename GridType>
        static typename GridType::Ptr receiveGrid
        (
            PstreamBuffers& pBufSize,
            PstreamBuffers& pBufGrid,
            const label node
        );

        //- Combine (topology union) collection of grids in one and remove voxels inside geometry
        template<typename GridType>
        void combineGrids
        (
            const word name,
            std::vector<std::vector<typename GridType::Ptr>>&,
            std::vector<typename GridType::Ptr>&
        );

        //- Combine (topology union) collection of grids in one and remove voxels inside geometry
        void combineGrids
        (
            const word name,
            std::vector<std::vector<FloatGrid::Ptr>>&,
            std::vector<FloatGrid::Ptr>&,
            const std::vector<BoolGrid::Ptr>&
        );

        //- Combine (topology union) collection of grids in one and remove voxels inside geometry
        void combineGrids
        (
            const word name,
            std::vector<std::vector<FloatGrid::Ptr>>&,
            std::vector<FloatGrid::Ptr>&,
            const std::vector<FloatGrid::Ptr>&
        );

        //- Topoly union of valid boundaries level sets
        FloatGrid::Ptr topologyUnionLevelSets
        (
            const std::vector<std::vector<FloatGrid::Ptr>>& boundaryCellLevels,
            const label cellLevel
        );

        //- Computes csgUnion of boundary level sets
        // if convertToSDF == true tries to rebuild the level set to get inside/outside
        FloatGrid::Ptr csgUnionLevelSets
        (
            const std::vector<std::vector<FloatGrid::Ptr>>& boundaryCellLevels,
            const label cellLevel,
            const bool convetToSDF = false,
            const scalar halfWidth = 3
        );

        //- Populates with active voxels the inside of a narrowband levelSet
        std::vector<FloatGrid::Ptr> getInteriorGrids
        (
            const std::vector<FloatGrid::Ptr>& narrowBandLevelSets,
            const label maxShellLevel
        );

        FloatGrid::Ptr createCurvatureMask
        (
            const FloatGrid& inGrid,
            const scalar curvature
        );

        ////- Creates partial SDF from whole geometry and maskGrid
        //void createPartialNarrowBandLevelSetGrids
        //(
        //    const triSurface& inputSurf,
        //    const label minLevel,
        //    const label maxLevel,
        //    const word regionName,
        //    std::vector<FloatGrid::Ptr>& SDFCellLevels,
        //    std::vector<FloatGrid::Ptr>& boundaryN,
        //    std::vector<BoolGrid::Ptr>& interiorGrids,
        //    std::vector<IndexGrid::Ptr>& polygonIdCellLevels
        //);

        //- creates a vector of narrowband levelSets
        std::vector<FloatGrid::Ptr> multiResSDFandDispl
        (
            const triSurface& inputSurf,
            const label minLevel,
            const label maxLevel,
            std::vector<BoolGrid::Ptr>& interiorGrids,
            //std::vector<IndexGrid::Ptr>& polygonIdCellLevels
            std::vector<Vec3dGrid::Ptr>& displacementGrids
        );

        //- creates a vector of narrowband levelSets
        std::vector<FloatGrid::Ptr> multiResSDF
        (
            const triSurface& inputSurf,
            const label minLevel,
            const label maxLevel,
            std::vector<BoolGrid::Ptr>& interiorGrids
            //std::vector<IndexGrid::Ptr>& polygonIdCellLevels
            //std::vector<Vec3dGrid::Ptr>& displacementGrids
        );

        const triSurface volumeToMesh
        (
            const scalar iso,
            FloatGrid::Ptr refGrid,
            const word surfName
        );

        //- creates a vector of narrowband levelSets
        // in parallel
        std::vector<FloatGrid::Ptr> createNarrowBandLevelSetGridsPar
        (
            const triSurface& inputSurf
        );

        //- multithreaded gather scatter
        template<class T>
        void allGatherList
        (
            List<T>& list
        );

        //- sends and receive grids bounding boxes for each cellLevel
        template<typename GridType>
        List<List<boundBox>> sendReceiveBoundBoxes
        (
            const std::vector<typename GridType::Ptr>& grids
        );

        //- sends and receive grids across processors
        template<typename GridType>
        std::vector<openvdb::GridPtrVecPtr> sendReceiveGrids
        (
            std::vector<typename GridType::Ptr>& procGrids,
            const word& gridsName,
            bool GridPtrVecPtr
        );

        //- sends and receive grids across processors
        template<typename GridType>
        std::vector<std::vector<typename GridType::Ptr>> sendReceiveGrids
        (
            std::vector<typename GridType::Ptr>& procGrids,
            const word& gridsName
        );

        //- create mesh fields out of VDB fields (distance field, cellLevel, etc)
        template<class Type>
        void createVolField
        (
            fvMesh& mesh,
            const std::vector<Type>& field,
            const word& name
        );

        //- gather list of VDB grids and interior of refinement grading,
        //  do topology union of grids and remove interior voxels
        void intersectProcGrids
        (
            std::vector<FloatGrid::Ptr>& procGrids,
            std::vector<BoolGrid::Ptr>& interiorGrids
        );

        //- gather list of VDB grids
        template<typename GridType>
        void gatherProcGrids
        (
            std::vector<typename GridType::Ptr>& procGrids
        );

        //- remove from grid all voxels with negative value in SDFGrid
        void removeSolidVoxels
        (
            FloatGrid::Ptr SDFGrid,
            FloatGrid::Ptr grid
        );

        //- Generate Level Of Detail sequence of the input grid (at minLevel).
        // restrict (fine->coarse) from minLevel to level 0
        // prolongate (coarse->fine) from minLevel to maxLevel
        std::vector<FloatGrid::Ptr> generateLODSequence
        (
            const FloatGrid::Ptr grid,
            const label minLevel,
            const label maxLevel
        );

        //- creates a patch Unsigned Distance Fields
        // at coarser resolution to calculate grading
        // and at minLevel resolution
        // Fills bufferCellLevels and patchCellLevels for surfaceGrading
        void multiResUDF
        (
            const word& name,
            const triSurface& inputSurf,
            const triSurface& globalSurf,
            const label regioni,
            const refinementSurfaces& surfaces,
            const label nCellsBetweenLevels,
            std::vector<BoolGrid::Ptr> interiorGrids,
            std::vector<std::vector<FloatGrid::Ptr>>& bufferCellLevels,
            std::vector<std::vector<FloatGrid::Ptr>>& patchCellLevels,
            bool curvature = false
        );

        //- calculates a partial SDF from intersection between high resolution UDF
        // of a patch and low resolution SDF of whole geometry
        BoolGrid::Ptr partialInterior
        (
            FloatGrid::Ptr fineUDF,
            const label surfLevel,
            const label nCellsBetweenLevels,
            std::vector<FloatGrid::Ptr> SDFCellLevels,
            BoolGrid::Ptr mask
        );

        //- creates a dense levelSet in the interior of a searchableBox
        FloatGrid::Ptr searchableBoxToLevelSet
        (
            const point min,
            const point max,
            const label cellLevel
        );

        //- creates a vector of levelSets for a given surface
        //- (with refinement grading)
        //- with halfWidth from distances
        std::vector<FloatGrid::Ptr> distanceRefinementGrading
        (
            const word& name,
            const triSurface& inputSurf,
            const List<labelVector>& distanceLevels,
            const scalarField& distances,
            const label nCellsBetweenLevels,
            const openvdb::CoordBBox blockMeshBBox,
            std::vector<BoolGrid::Ptr>& refinementInterior
        );

        //- creates a vector of levelSets for a given surface
        //- with filled interior and refinement grading
        std::vector<FloatGrid::Ptr> refinementRegionGrading
        (
            const word& name,
            const triSurface& inputSurf,
            const label insideLevel,
            const label nCellsBetweenLevels,
            const openvdb::CoordBBox blockMeshBBox,
            std::vector<BoolGrid::Ptr>& refinementInterior
        );

        //- exclude voxels which intersect a coarser voxel
        void trimGrids
        (
            const std::vector<FloatGrid::Ptr>& grids,
            const label maxSurfLevel,
            const label nUpLevels
        );

        //- exclude voxels which intersect a coarser voxel
        template<typename GridType>
        void trimGrids
        (
            const std::vector<typename GridType::Ptr>& grids
        );

        //- topologyDifference of vector of grids
        void topologyDifference
        (
            const std::vector<FloatGrid::Ptr>& resultGrids,
            const std::vector<FloatGrid::Ptr>& grids
        );

        //- topologyUnion of vector of grids
        void topologyUnion
        (
            const std::vector<FloatGrid::Ptr>& resultGrids,
            const std::vector<FloatGrid::Ptr>& grids
        );

        //- Converts a boundingBox from the index space of cellLevel 0 to space index of any cellLevel
        openvdb::CoordBBox getLevelNBoundingBox
        (
            const openvdb::math::BBox<openvdb::Vec3s>& bBox,
            const label& cellLevel,
            const label& decomposeDir
        );

        //- deepCopy of vector of grids
        std::vector<FloatGrid::Ptr> deepCopy
        (
            const std::vector<FloatGrid::Ptr>& grids
        );

        //- exclude voxels outside given bounding box
        openvdb::Index64 activeVoxelCountInBBox
        (
            const std::vector<FloatGrid::Ptr>& grids,
            const label& maxSurfLevel,
            const openvdb::math::BBox<openvdb::Vec3s>& bBox,
            label decomposeDir = -1
        );

        //- exclude voxels outside given bounding box
        void clipGrids
        (
            const std::vector<FloatGrid::Ptr>& grids,
            const label& maxSurfLevel,
            const openvdb::math::BBox<openvdb::Vec3s>& bBox,
            label decomposeDir = -1
        );

        //- exclude voxels outside given bounding box
        void clipGrids
        (
            const std::vector<FloatGrid::Ptr>& grids,
            const label& maxSurfLevel,
            const openvdb::CoordBBox& bBox
        );

        //- creates a vector of levelSets with refinement grading at decreasing cellLevels
        // from surface to farfield
        std::vector<FloatGrid::Ptr> surfaceGrading
        (
            const word& name,
            const std::vector<std::vector<FloatGrid::Ptr>>& patchCellLevels,
            const std::vector<std::vector<FloatGrid::Ptr>>& bufferCellLevels,
            const label maxSurfLevel,
            const label nCellsBetweenLevels,
            const openvdb::CoordBBox blockMeshBBox,
            std::vector<BoolGrid::Ptr>& refinementInterior
        );

        //- calculate coarser refinement grading grids
        void calculateBufferGrids
        (
            const word& name,
            const label bufferLevel,
            const label nCellsBetweenLevels,
            FloatGrid::Ptr bufferGrid,
            label& voxelCount,
            std::vector<FloatGrid::Ptr>& gradingGrids,
            std::vector<BoolGrid::Ptr>& refinementInterior
        );

        //- calculate coarser refinement grading grids
        void calculateBufferGrids
        (
            const word& name,
            const label bufferLevel,
            const scalarField& isoValues,
            const scalarField& widths,
            FloatGrid::Ptr bufferGrid,
            label& voxelCount,
            std::vector<FloatGrid::Ptr>& gradingGrids,
            std::vector<BoolGrid::Ptr>& refinementInterior
        );

        //- Return the voxelSize
        const scalar&  voxelSize() const
        {
            return voxelSize_;
        }

        //- Return the voxelSize
        const label&  maxCellLevel() const
        {
            return maxCellLevel_;
        }

    // Member Operators

    // Friend Functions

    // Friend Operators

    // IOstream Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #include "foamVDBI.H"

#ifdef NoRepository
#   include "foamVDBTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
