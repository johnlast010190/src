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

\*---------------------------------------------------------------------------*/

#include "foamVDB.H"

//clash with macro expansion in <openvdb_dir>/include/openvdb/math/Vec3.h:663
#undef Log
#include <openvdb/tools/MeshToVolume.h> //"MeshToVolume.h"
#include "MultiResGrid.h"
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/LevelSetUtil.h> //for sdfInteriorMask
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/TopologyToLevelSet.h>
#include <openvdb/tools/Composite.h> // for csgUnion()
#include <openvdb/tools/Morphology.h> //for dilation and erosion
#include <openvdb/tools/FastSweeping.h> //for sdfToSdf
#include <openvdb/tools/Diagnostics.h> //for checkLevelSet
#include <openvdb/io/Stream.h>
//#include <openvdb/io/Archive.h> //for DEFAULT_COMPRESSION_FLAGS
//#include <openvdb/points/StreamCompression.h>
#include <openvdb/tools/GridOperators.h> //for meanCurvature
#include <openvdb/tools/LevelSetAdvect.h>
#include <openvdb/tree/LeafManager.h>
#include <tbb/global_control.h>
#include <atomic>
#define Log if (log) Info

#include "meshes/primitiveShapes/point/point.H"
#include "meshes/meshShapes/labelledTri/labelledTri.H"
#include "triSurface/patches/geometricSurfacePatchList.H"
//#include "triSurface/booleanOps/booleanSurface/booleanSurface.H"
#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "triSurface/surfaceFeatures/surfaceFeatures.H"
#include "include/parallelForAll.H"
#include "decompositionMethod/decompositionMethod.H"
#include "geomDecomp/geomDecomp.H"


using MultiResGrid    = openvdb::tools::MultiResGrid<FloatGrid::TreeType>;
using MultiResVecGrid = openvdb::tools::MultiResGrid<Vec3dGrid::TreeType>;

namespace Foam
{
    const openvdb::Coord COARSE_TO_FINE[8] =
    {
        openvdb::Coord(0, 0, 0),
        openvdb::Coord(1, 0, 0),
        openvdb::Coord(1, 1, 0),
        openvdb::Coord(0, 1, 0),
        openvdb::Coord(0, 0, 1),
        openvdb::Coord(1, 0, 1),
        openvdb::Coord(1, 1, 1),
        openvdb::Coord(0, 1, 1)
    };

    // Voxel-face adjacent neighbours
    const openvdb::Coord COORD_OFFSETS[6] =
    {
        openvdb::Coord(-1,  0,  0), // f0 mX
        openvdb::Coord( 1,  0,  0), // f1 pX
        openvdb::Coord( 0, -1,  0), // f2 mY
        openvdb::Coord( 0,  1,  0), // f3 pY
        openvdb::Coord( 0,  0, -1), // f4 mZ
        openvdb::Coord( 0,  0,  1)  // f5 pZ
    };

    // edge-face connectivity in COORD_OFFSETS indices
    const List<labelList> edgeFace
    ({
        labelList{0, 2}, // e0  mX mY
        labelList{0, 5}, // e1  mX pZ
        labelList{0, 3}, // e2  mX pY
        labelList{0, 4}, // e3  mX mZ
        labelList{1, 4}, // e4  pX mZ
        labelList{1, 3}, // e5  pX pY
        labelList{1, 5}, // e6  pX pZ
        labelList{1, 2}, // e7  pX mY
        labelList{2, 4}, // e8  mY mZ
        labelList{2, 5}, // e9  mY pZ
        labelList{3, 4}, // e10 pY mZ
        labelList{3, 5}  // e11 pY pZ
    });

    // Voxel-face adjacent fine neighbours
    const openvdb::Coord COORD_OFFSETS_SPLIT[24] =
    {
        openvdb::Coord(-1,  0,  0), //  0 mX 0
        openvdb::Coord(-1,  0,  1), //  1 mX 4
        openvdb::Coord(-1,  1,  1), //  2 mX 7
        openvdb::Coord(-1,  1,  0), //  3 mX 3
        openvdb::Coord( 2,  0,  0), //  4  pX 1
        openvdb::Coord( 2,  1,  0), //  5  pX 2
        openvdb::Coord( 2,  1,  1), //  6  pX 6
        openvdb::Coord( 2,  0,  1), //  7  pX 5
        openvdb::Coord( 0, -1,  0), //  8 mY 0
        openvdb::Coord( 1, -1,  0), //  9 mY 1
        openvdb::Coord( 1, -1,  1), // 10 mY 5
        openvdb::Coord( 0, -1,  1), // 11 mY 4
        openvdb::Coord( 1,  2,  0), // 12  pY 2
        openvdb::Coord( 0,  2,  0), // 13  pY 3
        openvdb::Coord( 0,  2,  1), // 14  pY 7
        openvdb::Coord( 1,  2,  1), // 15  pY 6
        openvdb::Coord( 0,  0, -1), // 16 mZ 0
        openvdb::Coord( 0,  1, -1), // 17 mZ 3
        openvdb::Coord( 1,  1, -1), // 18 mZ 2
        openvdb::Coord( 1,  0, -1), // 19 mZ 1
        openvdb::Coord( 0,  0,  2), // 20  pZ 4
        openvdb::Coord( 1,  0,  2), // 21  pZ 5
        openvdb::Coord( 1,  1,  2), // 22  pZ 6
        openvdb::Coord( 0,  1,  2)  // 23  pZ 7
    };

    // edge-fine voxels connectivity in COORD_OFFSETS_SPLIT indices
    const List<labelList> fineEdgeFace
    ({
        labelList{ 8, 11}, // e0 mX
        labelList{ 0,  1}, // e0 mY

        labelList{20, 23}, // e1 mX
        labelList{ 1,  2}, // e1 pZ

        labelList{13, 14}, // e2 mX
        labelList{ 2,  3}, // e2 pY

        labelList{16, 17}, // e3 mX
        labelList{ 0,  3}, // e3 mZ

        labelList{18, 19}, // e4 pX
        labelList{ 4,  5}, // e4 mZ

        labelList{12, 15}, // e5 pX
        labelList{ 5,  6}, // e5 pY

        labelList{21, 22}, // e6 pX
        labelList{ 6,  7}, // e6 pZ

        labelList{ 9, 10}, // e7 pX
        labelList{ 4,  7}, // e7 mY

        labelList{16, 19}, // e8 mY
        labelList{ 8,  9}, // e8 mZ

        labelList{20, 21}, // e9 mY
        labelList{10, 11}, // e9 pZ

        labelList{17, 18}, // e10 pY
        labelList{12, 13}, // e10 mZ

        labelList{22, 23}, // e11 pY
        labelList{14, 15}  // e11 pZ
    });

    // Voxel-edge adjacent neighbours
    const openvdb::Coord EDGE_COORD_OFFSETS[12] =
    {
        openvdb::Coord(-1, -1,  0), // mX e0  0-4
        openvdb::Coord(-1,  0,  1), // mX e1  4-7
        openvdb::Coord(-1,  1,  0), // mX e2  7-3
        openvdb::Coord(-1,  0, -1), // mX e3  3-0

        openvdb::Coord( 1,  0, -1), // pX e4  1-2
        openvdb::Coord( 1,  1,  0), // pX e5  2-6
        openvdb::Coord( 1,  0,  1), // pX e6  6-5
        openvdb::Coord( 1, -1,  0), // pX e7  5-1

        openvdb::Coord( 0, -1, -1), // mY e8  0-1
        openvdb::Coord( 0, -1,  1), // mY e9  5-4

        openvdb::Coord( 0,  1, -1), // pY e10 2-3
        openvdb::Coord( 0,  1,  1), // pY e11 7-6
    };

    // Voxel-edge adjacent fine neighbours
    const openvdb::Coord FINE_EDGE_COORD_OFFSETS[24] =
    {
        openvdb::Coord(-1, -1,  0), // mX e0  0-4
        openvdb::Coord(-1, -1,  1), // mX e0  0-4

        openvdb::Coord(-1,  0,  2), // mX e1  4-7
        openvdb::Coord(-1,  1,  2), // mX e1  4-7

        openvdb::Coord(-1,  2,  0), // mX e2  7-3
        openvdb::Coord(-1,  2,  1), // mX e2  7-3

        openvdb::Coord(-1,  0, -1), // mX e3  3-0
        openvdb::Coord(-1,  1, -1), // mX e3  3-0

        openvdb::Coord( 2,  0, -1), // pX e4  1-2
        openvdb::Coord( 2,  1, -1), // pX e4  1-2

        openvdb::Coord( 2,  2,  0), // pX e5  2-6
        openvdb::Coord( 2,  2,  1), // pX e5  2-6

        openvdb::Coord( 2,  0,  2), // pX e6  6-5
        openvdb::Coord( 2,  1,  2), // pX e6  6-5

        openvdb::Coord( 2, -1,  0), // pX e7  5-1
        openvdb::Coord( 2, -1,  1), // pX e7  5-1

        openvdb::Coord( 0, -1, -1), // mY e8  0-1
        openvdb::Coord( 1, -1, -1), // mY e8  0-1

        openvdb::Coord( 0, -1,  2), // mY e9  5-4
        openvdb::Coord( 1, -1,  2), // mY e9  5-4

        openvdb::Coord( 0,  2, -1), // pY e10 2-3
        openvdb::Coord( 1,  2, -1), // pY e10 2-3

        openvdb::Coord( 0,  2,  2), // pY e11 7-6
        openvdb::Coord( 1,  2,  2), // pY e11 7-6
    };

    const List<face> localFaces
    ({
        face{0, 4, 7, 3}, // x-min
        face{1, 2, 6, 5}, // x-max
        face{0, 1, 5, 4}, // y-min
        face{2, 3, 7, 6}, // y-max
        face{0, 3, 2, 1}, // z-min
        face{4, 5, 6, 7}  // z-max
    });

    const List<word> faceNames
    ({
        "minusX",
        " plusX",
        "minusY",
        " plusY",
        "minusZ",
        " plusZ"
    });

    const List<edge> localEdges
    ({
        edge{0, 4}, // e0  x-min
        edge{4, 7}, // e1
        edge{7, 3}, // e2
        edge{3, 0}, // e3
        edge{1, 2}, // e4  x-max
        edge{2, 6}, // e5
        edge{6, 5}, // e6
        edge{5, 1}, // e7
        edge{0, 1}, // e8  y-min
        edge{5, 4}, // e9
        edge{2, 3}, // e10 y-max
        edge{7, 6}  // e11
    });

    // edges expressed as connection between start point
    // of edges with lower index
    // e.g. localEdge[8] = edge{localEdge[0].start(), localEdge[4].start()};
    const List<labelList> edgeEdges
    ({
        labelList{0, 0}, //0 to 7 not used
        labelList{0, 0},
        labelList{0, 0},
        labelList{0, 0},
        labelList{0, 0},
        labelList{0, 0},
        labelList{0, 0},
        labelList{0, 0},
        labelList{0, 4}, // 8
        labelList{7, 1},
        labelList{5, 3},
        labelList{2, 6}
    });

    const List<labelList> edgeLoops
    ({
        labelList{     0,      1,      2,      3}, // x-min
        labelList{     4,      5,      6,      7}, // x-max
        labelList{     8, /*-*/7,      9, /*-*/0}, // y-min
        labelList{    10, /*-*/2,     11, /*-*/5}, // y-max
        labelList{/*-*/3,/*-*/10, /*-*/4, /*-*/8}, // z-min
        labelList{/*-*/9, /*-*/6,/*-*/11, /*-*/1}  // z-max
    });

    // TBB body object for threaded copy of points from VDB to OpenFOAM format
    struct PointListCopy
    {
        using PointList = std::unique_ptr<openvdb::Vec3s[]>;

        const PointList& pointsIn_;
        pointField& pointsOut_;

        PointListCopy
        (
            const PointList& pointsIn,
            pointField& pointsOut
        )
        :
            pointsIn_(pointsIn),
            pointsOut_(pointsOut)
        {
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t n = range.begin(); n < range.end(); ++n)
            {
                point p
                (
                    pointsIn_[n].x(),
                    pointsIn_[n].y(),
                    pointsIn_[n].z()
                );

                pointsOut_[n] = p;
            }
        }
    }; //PointListCopy


    // TBB body object for threaded copy of triangles from VDB to OpenFOAM format
    struct CopyTriangles
    {
        const openvdb::tools::PolygonPool& polygons_;
        List<labelledTri>& faces_;

        CopyTriangles
        (
            const openvdb::tools::PolygonPool& polygons,
            List<labelledTri>& faces
        )
        :
            polygons_(polygons),
            faces_(faces)
        {
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t n = range.begin(); n < range.end(); ++n)
            {
                labelledTri f
                {
                    label(polygons_.triangle(n)[0]),
                    label(polygons_.triangle(n)[1]),
                    label(polygons_.triangle(n)[2])
                };

                faces_[n] = f;
            }
        }
    }; //CopyTriangles


    // TBB body object for threaded split and copy of quads from VDB to OpenFOAM format
    struct CopyQuads
    {
        const openvdb::tools::PolygonPool& polygons_;
        List<labelledTri>& faces_;

        CopyQuads
        (
            const openvdb::tools::PolygonPool& polygons,
            List<labelledTri>& faces
        )
        :
            polygons_(polygons),
            faces_(faces)
        {
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t n = range.begin(); n < range.end(); ++n)
            {
                // split quad in 2 triangles
                labelledTri f
                {
                    label(polygons_.quad(n)[0]),
                    label(polygons_.quad(n)[1]),
                    label(polygons_.quad(n)[2])
                };

                faces_[2*n] = f;

                labelledTri f1
                {
                    label(polygons_.quad(n)[0]),
                    label(polygons_.quad(n)[2]),
                    label(polygons_.quad(n)[3])
                };

                faces_[2*n + 1] = f1;
            }
        }
    }; //CopyQuads
} //namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelListList Foam::foamVDB::calcVoxelVoxels
(
    std::vector<IndexGrid::Ptr>&  globalIDGrids,
    const label nVoxels
)
{
    Timer t("voxel-voxel connections");

    labelListList voxelVoxels(nVoxels);

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, maxCellLevel_+1),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t cellLevel = r.begin(); cellLevel < r.end(); ++cellLevel)
            {
                openvdb::tools::foreach
                (
                    globalIDGrids[cellLevel]->cbeginValueOn(),
                    [&](const IndexGrid::ValueOnCIter& iter)
                    {
                        DynamicList<label> nbrs;
                        nbrs.reserve(16);

                        openvdb::Coord ijk, nijk;
                        openvdb::Int32 nID;

                        IndexGrid::ConstAccessor acc
                        (
                            globalIDGrids[cellLevel]->getConstAccessor()
                        );

                        IndexGrid::ConstAccessor coarseAcc
                        (
                            //globalIDGrids[cellLevel-1]->getConstAccessor()
                            globalIDGrids[cellLevel - label(cellLevel > 0)]->getConstAccessor()
                        );

                        IndexGrid::ConstAccessor fineAcc
                        (
                            //globalIDGrids[cellLevel+1]->getConstAccessor()
                            globalIDGrids[cellLevel + label(cellLevel < maxCellLevel_)]->getConstAccessor()
                        );

                        ijk = iter.getCoord();

                        for (size_t i = 0; i < 6; ++i)
                        {
                            nijk = ijk + COORD_OFFSETS[i];

                            if (acc.probeValue(nijk, nID) && nID >= 0)
                            {
                                nbrs.append(nID);
                                continue;
                            }

                            //check coarser cellLevel
                            if (cellLevel > 0)
                            {
                                switch ( (nijk[0] & 1) | ((nijk[1] & 1) << 1) | ((nijk[2] & 1) << 2) )
                                {
                                    case 0:// all even
                                        nID = coarseAcc.getValue(nijk>>1);
                                        break;
                                    case 1:// x is odd
                                        nID = coarseAcc.getValue(nijk.offsetBy(-1,0,0)>>1);
                                        break;
                                    case 2:// y is odd
                                        nID = coarseAcc.getValue(nijk.offsetBy(0,-1,0)>>1);
                                        break;
                                    case 3:// x&y are odd
                                        nID = coarseAcc.getValue(nijk.offsetBy(-1,-1,0)>>1);
                                        break;
                                    case 4:// z is odd
                                        nID = coarseAcc.getValue(nijk.offsetBy(0,0,-1)>>1);
                                        break;
                                    case 5:// x&z are odd
                                        nID = coarseAcc.getValue(nijk.offsetBy(-1,0,-1)>>1);
                                        break;
                                    case 6:// y&z are odd
                                        nID = coarseAcc.getValue(nijk.offsetBy(0,-1,-1)>>1);
                                        break;
                                    case 7:// all are odd
                                        nID = coarseAcc.getValue(nijk.offsetBy(-1,-1,-1)>>1);
                                        break;
                                } //switch

                                if (nID >= 0)
                                {
                                    nbrs.append(nID);
                                    continue;
                                }
                            }

                            //check finer cellLevel
                            if (cellLevel < maxCellLevel_)
                            {
                                for (label j = 0; j < 4; j++)
                                {
                                    nijk = (ijk << 1) + COORD_OFFSETS_SPLIT[4*i + j]; //face adjacent

                                    if (fineAcc.probeValue(nijk, nID) && nID >= 0)
                                    {
                                        nbrs.append(nID);
                                    }
                                }
                            }
                        } // for 6 voxels neighbours

                        voxelVoxels[iter.getValue()] = std::move(nbrs);
                    },
                    /*threaded*/true,
                    /*shareOp*/false
                );
            } //for cellLevel
        } //TBB body voxel-voxel connections
    );

    return voxelVoxels;
} //calcVoxelVoxels


void Foam::foamVDB::refineDanglingVoxels
(
    std::vector<FloatGrid::Ptr>&  cellLevelGrids
)
{
    Timer t("Dangling cells refinement");

    // check only voxels at cellLevel interface
    std::vector<BoolGrid::Ptr> interfaceGrids(maxCellLevel_, nullptr);

    for (label cellLevel = 0; cellLevel < maxCellLevel_; ++cellLevel)
    {
        BoolGrid::Ptr interfaceGrid = BoolGrid::create();

        interfaceGrid->topologyUnion(*cellLevelGrids[cellLevel]);

        openvdb::tools::dilateActiveValues(interfaceGrid->tree(), 1);

        interfaceGrid->topologyDifference(*cellLevelGrids[cellLevel]);

        openvdb::tools::dilateActiveValues(interfaceGrid->tree(), 1);

        interfaceGrid->topologyIntersection(*cellLevelGrids[cellLevel]);

        interfaceGrids[cellLevel] = interfaceGrid;
    }

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, maxCellLevel_),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t cellLevel = r.begin(); cellLevel < r.end(); ++cellLevel)
            {
                openvdb::tools::foreach
                (
                    interfaceGrids[cellLevel]->cbeginValueOn(),
                    [&](const BoolGrid::ValueOnCIter& iter)
                    {
                        FloatGrid::Accessor coarseAcc
                        (
                            cellLevelGrids[cellLevel]->getAccessor()
                        );

                        FloatGrid::Accessor fineAcc
                        (
                            cellLevelGrids[cellLevel+1]->getAccessor()
                        );

                        label nSplitFaces = 0;

                        // get coord at finer level
                        const openvdb::Coord ijk = iter.getCoord() << 1;

                        for (size_t i = 0; i < 6; ++i)
                        {
                            const openvdb::Coord nijk = ijk + COORD_OFFSETS_SPLIT[4*i];

                            if (fineAcc.isValueOn(nijk))
                            {
                                nSplitFaces++;
                            }
                        }

                        if (nSplitFaces > 2)
                        {
                            coarseAcc.setValueOff(iter.getCoord());

                            for (label i = 0; i < 8; ++i)
                            {
                                fineAcc.setValueOn(ijk + COARSE_TO_FINE[i]);
                            }
                        }
                    },
                    /*threaded*/true,
                    /*shareOp*/false
                );
            } //for cellLevel
        } //TBB body dangling cells
    );
} // refineDanglingVoxels


Foam::label Foam::foamVDB::syncGrids
(
    std::vector<FloatGrid::Ptr>&  cellLevelGrids,
    std::vector<IndexGrid::Ptr>&  globalIDGrids,
    const List<List<boundBox>>& procBBoxes,
    const label procN
)
{
    Timer timer("sync globalIDGrids multithread");

    label nVoxels = 0;

    label minID, maxID;

    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
        PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

        const boundBox& myProcBB = procBBoxes[Pstream::myProcNo()][cellLevel];

        //send updated gloabalID grid
        if (procN == Pstream::myProcNo())
        {
            tbb::parallel_for
            (
                tbb::blocked_range<size_t>(0, Pstream::nProcs()),
                SendGridOp<IndexGrid>
                (
                    globalIDGrids[cellLevel],
                    pBufSize,
                    pBufGrid,
                    myProcBB,
                    procBBoxes[cellLevel],
                    /*isCellLevelGrid*/true
                )
            );
        }

        pBufSize.finishedSends();
        pBufGrid.finishedSends();

        //receive from procN
        if (procN < Pstream::myProcNo())
        {
            const boundBox& otherProcBB = procBBoxes[procN][cellLevel];

            if (myProcBB.overlaps(otherProcBB))
            {
                IndexGrid::Ptr globalIDProcI =
                    receiveGrid<IndexGrid>
                    (
                        pBufSize,
                        pBufGrid,
                        procN
                    );

                cellLevelGrids[cellLevel]->topologyDifference(*globalIDProcI);

                globalIDProcI->evalMinMax(minID, maxID);

                nVoxels = max(nVoxels, maxID);

                if (isValid<IndexGrid>(globalIDGrids[cellLevel]))
                {
                    //merge tree if already added by a previous proc
                    globalIDGrids[cellLevel]->tree().merge(globalIDProcI->tree()); //empties globalIDProcI tree
                }
                else
                {
                    globalIDGrids[cellLevel] = globalIDProcI;
                }
            } //if overlaps
        } //for proci (receive)
    } //for cellLevel

    return nVoxels;
} //syncGrids


Foam::labelList Foam::foamVDB::countVoxels
(
    const std::vector<FloatGrid::Ptr>&  cellLevelGrids,
    std::vector<IndexGrid::Ptr>&  globalIDGrids,
    pointField& voxelCentres
)
{
    Timer t("count voxels valueOn");

    labelList nVoxelsStart(maxCellLevel_ + 1, 0);

    forAll(nVoxelsStart, cellLevel)
    {
        IndexGrid::Ptr globalIDGrid = IndexGrid::create(-labelMax);
        globalIDGrid->setGridClass(openvdb::GRID_UNKNOWN);
        globalIDGrid->topologyUnion(*cellLevelGrids[cellLevel]);

        globalIDGrids[cellLevel] = globalIDGrid;

        if (cellLevel == 0) continue;
        nVoxelsStart[cellLevel] =
            nVoxelsStart[cellLevel-1] + cellLevelGrids[cellLevel-1]->activeVoxelCount();
    }

    label totalVoxels = nVoxelsStart[maxCellLevel_] + cellLevelGrids[maxCellLevel_]->activeVoxelCount();

    voxelCentres.setSize(totalVoxels);

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, maxCellLevel_+1),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t cellLevel = r.begin(); cellLevel < r.end(); ++cellLevel)
            {
                label nVoxels = nVoxelsStart[cellLevel];

                globalIDGrids[cellLevel]->tree().voxelizeActiveTiles();

                for (auto iter = globalIDGrids[cellLevel]->beginValueOn(); iter; ++iter)
                {
                    iter.setValue(nVoxels);

                    const openvdb::Coord ijk = iter.getCoord() << (maxCellLevel_ - cellLevel);

                    voxelCentres[nVoxels++] = point(ijk.x(), ijk.y(), ijk.z());
                }
            }
        }
    );

    return nVoxelsStart;
} //countVoxels


Foam::labelList Foam::foamVDB::countVoxelsLeaf
(
    const std::vector<FloatGrid::Ptr>&  cellLevelGrids,
    std::vector<IndexGrid::Ptr>&  globalIDGrids,
    pointField& voxelCentres
)
{
    Timer t("count voxels LeafManager");

    using LeafManager = openvdb::tree::LeafManager<IndexGrid::TreeType>;

    std::vector<LeafManager*> leafManagers;
    leafManagers.reserve(maxCellLevel_ + 1);

    labelList nVoxelsStart(maxCellLevel_ + 1, 0);

    forAll(nVoxelsStart, cellLevel)
    {
        IndexGrid::Ptr globalIDGrid = IndexGrid::create(-labelMax);
        globalIDGrid->setGridClass(openvdb::GRID_UNKNOWN);
        globalIDGrid->topologyUnion(*cellLevelGrids[cellLevel]);

        globalIDGrid->tree().voxelizeActiveTiles();

        globalIDGrids[cellLevel] = globalIDGrid;

        leafManagers.emplace_back(new LeafManager(globalIDGrid->tree()));

        if (cellLevel == 0) continue;
        nVoxelsStart[cellLevel] =
            nVoxelsStart[cellLevel-1] + leafManagers.end()[-2]->activeLeafVoxelCount();
    }

    label totalVoxels = nVoxelsStart[maxCellLevel_] + leafManagers.back()->activeLeafVoxelCount();

    voxelCentres.setSize(totalVoxels);

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, maxCellLevel_+1),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t cellLevel = r.begin(); cellLevel < r.end(); ++cellLevel)
            {
                label nVoxels = nVoxelsStart[cellLevel];

                LeafManager& mgr = *(leafManagers[cellLevel]);

                // Compute the prefix sum of offsets to active voxels
                size_t* mOffsets;
                size_t size=0, voxelCount = mgr.getPrefixSum(mOffsets, size, /*grainSize*/1);

                mgr.foreach
                (
                    [&](IndexGrid::TreeType::LeafNodeType& leaf, size_t leafIndex)
                    {
                        size_t nLeafVoxels = nVoxels + mOffsets[leafIndex];

                        for (auto iter = leaf.beginValueOn(); iter; ++iter)
                        {
                            iter.setValue(nLeafVoxels);

                            const openvdb::Coord ijk = iter.getCoord() << (maxCellLevel_ - cellLevel);

                            voxelCentres[nLeafVoxels++] = point(ijk.x(), ijk.y(), ijk.z());
                        }
                    }
                );
            }
        }
    );

    return nVoxelsStart;
} //countVoxelsLeaf


const std::vector<openvdb::Coord> Foam::foamVDB::getVoxelVertices
(
    const label pointsPerEdge,
    const openvdb::Coord& coord
)
{
    label xmin = pointsPerEdge * (coord.x());
    label xmax = pointsPerEdge * (coord.x() + 1);
    label ymin = pointsPerEdge * (coord.y());
    label ymax = pointsPerEdge * (coord.y() + 1);
    label zmin = pointsPerEdge * (coord.z());
    label zmax = pointsPerEdge * (coord.z() + 1);

    std::vector<openvdb::Coord> voxelVertices;

    voxelVertices.reserve(8);

    voxelVertices.emplace_back(openvdb::Coord(xmin, ymin, zmin));
    voxelVertices.emplace_back(openvdb::Coord(xmax, ymin, zmin));
    voxelVertices.emplace_back(openvdb::Coord(xmax, ymax, zmin));
    voxelVertices.emplace_back(openvdb::Coord(xmin, ymax, zmin));
    voxelVertices.emplace_back(openvdb::Coord(xmin, ymin, zmax));
    voxelVertices.emplace_back(openvdb::Coord(xmax, ymin, zmax));
    voxelVertices.emplace_back(openvdb::Coord(xmax, ymax, zmax));
    voxelVertices.emplace_back(openvdb::Coord(xmin, ymax, zmax));

    return voxelVertices;
} //getVoxelVertices


void Foam::foamVDB::addPoint
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
)
{
    openvdb::Int32 pointID;

    edgeCoords.emplace_back(p);

    // check if point already exists
    if (pointIDAcc.probeValue(p, pointID))
    {
        edge.append(pointID);
    }
    else // add point to pointList
    {
        pointLevelAcc.setValue(p, cellLevel);

        edge.append(nPoints);

        pointIDAcc.setValue(p, nPoints);

        //transform to global index space
        p <<= maxCellLevel_ - cellLevel;

        xyzPoints.emplace_back(p.x(), p.y(), p.z());
        vdbPoints.emplace_back(point(p.x(), p.y(), p.z()), nPoints);
        ++nPoints;
    }
} //addPoint


void Foam::foamVDB::addMidPoint
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
)
{
    labelList splitEdge;

    // duplicate last element
    edgeCoords.emplace_back(edgeCoords.back());
    // add midpoint
    edgeCoords[1] = p;

    splitEdge.append(edge.first());

    // check if point already exists
    if (pointIDAcc.probeValue(p, pointID))
    {
        splitEdge.append(pointID);
    }
    else // add point to pointList
    {
        pointLevelAcc.setValue(p, cellLevel);

        pointID = nPoints;

        splitEdge.append(nPoints);

        pointIDAcc.setValue(p, nPoints);

        //transform to global index space
        p <<= maxCellLevel_ - cellLevel;

        xyzPoints.emplace_back(p.x(), p.y(), p.z());
        vdbPoints.emplace_back(point(p.x(), p.y(), p.z()), nPoints);
        ++nPoints;
    }

    splitEdge.append(edge.last());

    edge = splitEdge;
} //addMidPoint


void Foam::foamVDB::addEmptyPatch
(
    fvMesh& mesh,
    const word& patchName,
    const dictionary& patchInfo
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches =
        const_cast<fvBoundaryMesh&>(mesh.boundary());

    label patchi = polyPatches.size();

    // Add polyPatch at the end
    polyPatches.setSize(patchi+1);
    polyPatches.set
    (
        patchi,
        polyPatch::New
        (
            patchName,
            patchInfo,
            patchi,
            polyPatches
        )
    );

    fvPatches.setSize(patchi+1);
    fvPatches.set
    (
        patchi,
        fvPatch::New
        (
            polyPatches[patchi],  // point to newly added polyPatch
            mesh.boundary()
        )
    );
} // addEmptyPatch

Foam::labelListList Foam::foamVDB::groupSameLevelRegions
(
    const wordList& regNames,
    const refinementSurfaces& surfaces,
    const scalarField& featAngles
)
{
    label nRegions = regNames.size();

    boolList alreadyAdded(nRegions, false);

    labelListList sameMinMaxLevel;

    forAll(regNames, regioni)
    {
        if (!alreadyAdded[regioni])
        {
            DynamicList<label> thisMinMaxLevel;

            thisMinMaxLevel.append(regioni);

            alreadyAdded[regioni] = true;

            for (label regionj = regioni+1 ; regionj < nRegions; regionj++)
            {
                if
                (
                    surfaces.minLevel()[regioni] == surfaces.minLevel()[regionj]
                 && surfaces.maxLevel()[regioni] == surfaces.maxLevel()[regionj]
                 && featAngles[regioni] == featAngles[regionj]
                )
                {
                    thisMinMaxLevel.append(regionj);

                    alreadyAdded[regionj] = true;
                }
            }
            sameMinMaxLevel.append(thisMinMaxLevel);
        }
    }

    return sameMinMaxLevel;
}

openvdb::Index64 Foam::foamVDB::onVoxelsInCoordBbox
(
    FloatGrid::Ptr grid,
    const openvdb::CoordBBox& bBox
)
{
    openvdb::Index64 activeLeafVoxelCount = 0;

    using TreeType = FloatGrid::TreeType;

    // Iterate over pointers to const LeafNodes.
    for (TreeType::LeafCIter iter = grid->tree().cbeginLeaf(); iter; ++iter)
    {
        const TreeType::LeafNodeType* leaf = iter.getLeaf();

        openvdb::CoordBBox leafBBox = leaf->getNodeBoundingBox();

        if (bBox.isInside(leafBBox))
        {
            activeLeafVoxelCount += leaf->onVoxelCount();
        }
        else if (bBox.hasOverlap(leafBBox))
        {
            for (auto vIter = leaf->cbeginValueOn(); vIter; ++vIter)
            {
                if (bBox.isInside(vIter.getCoord()))
                {
                    activeLeafVoxelCount++;
                }
            }
        }
    }

    return activeLeafVoxelCount;
}

std::vector<FloatGrid::Ptr> Foam::foamVDB::sliceGrids
(
    const std::vector<FloatGrid::Ptr>& grids,
    openvdb::math::BBox<openvdb::Vec3s>& bBoxSlice,
    const label& slicei,
    const Vector<label>& n,
    const label& dir,
    const label& desiredVoxels,
    const scalar& maxLoadUnbalance
)
{
    scalar minDelta = 1 / scalar(Foam::pow(2, maxCellLevel_));

    openvdb::math::BBox<openvdb::Vec3s> bBox = getBoundingBox(grids);

    point bBoxMin(bBox.min().x(), bBox.min().y(), bBox.min().z());
    point bBoxMax(bBox.max().x(), bBox.max().y(), bBox.max().z());

    reduce(
        std::tie(bBoxMin, bBoxMax),
        ParallelOp<minOp<point>, maxOp<point>>{}
    );

    scalar xMax = slicei == 0 ? bBoxMin.x() + ((bBoxMax.x() - bBoxMin.x()) / (n[dir]-slicei))
                              : bBoxSlice.max().x() + 2.*minDelta;

    scalar yMax = slicei == 0 ? bBoxMin.y() + ((bBoxMax.y() - bBoxMin.y()) / (n[dir]-slicei))
                              : bBoxSlice.max().y() + 2.*minDelta;

    scalar zMax = slicei == 0 ? bBoxMin.z() + ((bBoxMax.z() - bBoxMin.z()) / (n[dir]-slicei))
                              : bBoxSlice.max().z() + 2.*minDelta;

    bBoxSlice =
        openvdb::math::BBox<openvdb::Vec3s>
        (
            openvdb::Vec3s
            (
                bBoxMin.x(),
                bBoxMin.y(),
                bBoxMin.z()
            ),
            openvdb::Vec3s
            (
                dir == 0 ? xMax : bBoxMax.x(),
                dir == 1 ? yMax : bBoxMax.y(),
                dir == 2 ? zMax : bBoxMax.z()
            )
        );

    std::vector<FloatGrid::Ptr> sliceGrids = deepCopy(grids);

    for (size_t cellLevel=0; cellLevel < grids.size(); ++cellLevel)
    {
        sliceGrids[cellLevel]->tree().voxelizeActiveTiles();
    }

    label nVoxelsSlice =
        activeVoxelCountInBBox
        (
            sliceGrids,
            maxCellLevel_,
            bBoxSlice,
            dir
        );

    { //Timer timer("reduce nVoxelsSlice");
    reduce(nVoxelsSlice, sumOp<label>());
    }

    scalar unbalance = (desiredVoxels - nVoxelsSlice) / scalar(desiredVoxels);

    scalar unbalanceOld = unbalance;

    bool overshoot = false;

    bool triggeredOnce = false;

    if (Pstream::master())
    {
        std::cout<< "\nslice n" << slicei
            << " in direction " << (dir == 0 ? "x" : (dir == 1 ? "y" : "z"))
            << " desiredVoxels " << desiredVoxels
            << "\nnVoxelsSlice " <<nVoxelsSlice
            << " unbalance " << unbalance
            << " bBoxSlice " << bBoxSlice
            << " minDelta " << minDelta
            << std::endl;
    }

    scalar minUnbalance = 1;

    scalar delta = abs(unbalance) < 0.5 ? minDelta : 0.5;

    while (abs(unbalance) > maxLoadUnbalance)
    {
        if (sign(unbalance) == sign(unbalanceOld) && !overshoot)
        {
            delta *= 2.;
        }
        else
        {
            overshoot = true;
            delta /= 2.;
        }

        unbalanceOld = unbalance;

        delta =
            unbalance > 0
          ? Foam::max(Foam::min( abs(delta),  2),  minDelta)
          : Foam::min(Foam::max(-abs(delta), -2), -minDelta);

        bBoxSlice =
            openvdb::math::BBox<openvdb::Vec3s>
            (
                bBoxSlice.min(),
                openvdb::Vec3s
                (
                    dir == 0 ? bBoxSlice.max().x() + delta : bBoxSlice.max().x(),
                    dir == 1 ? bBoxSlice.max().y() + delta : bBoxSlice.max().y(),
                    dir == 2 ? bBoxSlice.max().z() + delta : bBoxSlice.max().z()
                )
            );

        nVoxelsSlice =
            activeVoxelCountInBBox
            (
                sliceGrids,
                maxCellLevel_,
                bBoxSlice,
                dir
            );

        { //Timer timer("reduce nVoxelsSlice");
        reduce(nVoxelsSlice, sumOp<label>());
        }

        unbalance = (desiredVoxels - nVoxelsSlice) / scalar(desiredVoxels);

        if (Pstream::master())
        {
            std::cout<< "nVoxelsSlice " << nVoxelsSlice
                << " unbalance " << unbalance
                << " bBoxSlice " << bBoxSlice
                << " delta to previous " << delta
                << std::endl;
        }

        // exit if trapped in infinite loop
        if (abs(unbalance) < minUnbalance)
        {
            minUnbalance = abs(unbalance);
        }
        else if (abs(unbalance) == minUnbalance && abs(delta) == minDelta)
        {
            if (triggeredOnce)
            {
                Info<< "breaking infinite loop..." << endl;
                break;
            }
            triggeredOnce = true;
        }
    } // while balance

    clipGrids
    (
        sliceGrids,
        maxCellLevel_,
        bBoxSlice,
        dir
    );

    return sliceGrids;
} // sliceGrids

void Foam::foamVDB::doDecomposition
(
    const std::vector<FloatGrid::Ptr>& grids,
    std::vector<std::vector<FloatGrid::Ptr>>& decomposedGrids,
    const Vector<label>& n,
    FixedList<direction, 3>& decompOrder,
    const label& nthDir,
    label procOffset,
    const label& nProcs,
    const scalar& maxLoadUnbalance
)
{
    const label decompDir = decompOrder[nthDir];

    if (n[decompDir] == 1)
    {
        fillProcGridsOrDecompose
        (
            grids,
            decomposedGrids,
            n,
            decompOrder,
            nthDir,
            procOffset, // procI
            procOffset, // offset
            nProcs,
            maxLoadUnbalance
        );

        return;
    }

    label nVoxels = getActiveVoxels<FloatGrid>(grids);

    reduce(nVoxels, sumOp<label>());

    const label desiredVoxels = nVoxels / n[decompDir];

    openvdb::math::BBox<openvdb::Vec3s> bBoxSlice;

    for (label i = 0; i < n[decompDir] - 1; ++i)
    {
        std::vector<FloatGrid::Ptr> aSliceGrid =
            sliceGrids
            (
                grids,
                bBoxSlice,
                i,              // slicei
                n,
                decompDir,      // dir
                desiredVoxels,
                maxLoadUnbalance
            );

        topologyDifference
        (
            grids,
            aSliceGrid
        );

        label procI = procOffset + i;

        label offset =
            (nthDir == 0)
          ? (i * (nProcs / n[decompDir]))
          : procOffset + (i * (nProcs / (n[decompDir] * n[decompOrder[nthDir-1]])));

        fillProcGridsOrDecompose
        (
            aSliceGrid,
            decomposedGrids,
            n,
            decompOrder,
            nthDir,
            procI,
            offset,
            nProcs,
            maxLoadUnbalance
        );

        if (i == n[decompDir] - 2)
        {
            procI++;

            offset =
                (nthDir == 0)
              ? ((i+1) * (nProcs / n[decompDir]))
              : procOffset + ((i+1) * (nProcs / (n[decompDir] * n[decompOrder[nthDir-1]])));

            fillProcGridsOrDecompose
            (
                grids,
                decomposedGrids,
                n,
                decompOrder,
                nthDir,
                procI,
                offset,
                nProcs,
                maxLoadUnbalance
            );
        } // if last slice
    } // for slice
} // doDecomposition

void Foam::foamVDB::fillProcGridsOrDecompose
(
    const std::vector<FloatGrid::Ptr>& grids,
    std::vector<std::vector<FloatGrid::Ptr>>& decomposedGrids,
    const Vector<label>& n,
    FixedList<direction, 3>& decompOrder,
    const label nthDir,
    const label procI,
    const label offset,
    const label nProcs,
    const scalar maxLoadUnbalance
)
{
    if (nthDir == 2) //decomposition along last axis
    {
        label nVoxels = getActiveVoxels<FloatGrid>(grids);

        reduce(nVoxels, sumOp<label>());

        openvdb::math::BBox<openvdb::Vec3s> bBox = getBoundingBox(grids);

        point bBoxMin(bBox.min().x(), bBox.min().y(), bBox.min().z());
        point bBoxMax(bBox.max().x(), bBox.max().y(), bBox.max().z());

        reduce(
            std::tie(bBoxMin, bBoxMax),
            ParallelOp<minOp<point>, maxOp<point>>{}
        );

        boundBox bb(bBoxMin, bBoxMax);

        if (Pstream::master())
        {
            Info<<"processor " << procI
                << " nVoxels " << nVoxels
                << " bBox " << bb
                << endl;
        }

        for (size_t cellLevel=0; cellLevel < grids.size(); ++cellLevel)
        {
            decomposedGrids[cellLevel][procI] = grids[cellLevel];
        }
    }
    else
    {
        doDecomposition
        (
            grids,
            decomposedGrids,
            n,
            decompOrder,
            nthDir + 1,
            offset,
            nProcs,
            maxLoadUnbalance
        );
    }
} // fillProcGridsOrDecompose

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVDB::foamVDB
(
    const scalar& voxelSize,
    const label& maxCellLevel
)
:
    voxelSize_(voxelSize),
    maxCellLevel_(maxCellLevel)
{

    // Initialize the OpenVDB library.  This must be called at least
    // once per program and may safely be called multiple times.
    openvdb::initialize();

    Info<< "OpenVDB using "
        << tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism)
        << " threads"
        << endl;
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


Foam::dictionary Foam::foamVDB::createDecompDict
(
    const Time& runTime,
    const word method
)
{
    dictionary decomposeDict;

    const label nProcs = Pstream::nProcs();

    decomposeDict.add("numberOfSubdomains", nProcs);

    if (nProcs == 1)
    {
        decomposeDict.add("method", "none");
        return decomposeDict;
    }

    if (method == "hierarchical")
    {
        decomposeDict.add("method", "hierarchical");
        decomposeDict.add(word("hierarchicalCoeffs"), dictionary());
        decomposeDict.subDict("hierarchicalCoeffs").add
        (
            word("order"),
            word("yxz")
        );
        decomposeDict.subDict("hierarchicalCoeffs").add
        (
            word("delta"),
            0.001
        );

        //Calculate hierarchical factors
        label np_3 = label(pow(nProcs, 1.0/3.0));

        label firstFactor = 0;
        for (label i = np_3; i <= nProcs; ++i)
        {
            if (nProcs % i == 0)
            {
                firstFactor = i;
                break;
            }
        }
        label remainder = nProcs/firstFactor;
        label np_2 = label(pow(remainder, 0.5));

        label secondFactor = 0;
        for (label i = np_2; i <= nProcs; ++i)
        {
            if (remainder % i == 0)
            {
                secondFactor = i;
                break;
            }
        }

        label thirdFactor = nProcs / firstFactor / secondFactor;

        labelVector n(thirdFactor, secondFactor, firstFactor);

        Info<< "Using hierarchical coefficients " << n << endl;

        decomposeDict.subDict("hierarchicalCoeffs").add("n", n);
    }
    else
    {
        decomposeDict.add("method", method);
    }

    return decomposeDict;
} //createDecompDict


Foam::dictionary Foam::foamVDB::checkDecomposeParDict
(
    const Time& runTime,
    labelList& nDomainsPerNode,
    label& procOffset
)
{
    // Always check for domains. I might run in serial on one node
    // but request decomposition on 8 procs
    dictionary decompositionDict =
        IOdictionary
        (
            IOobject
            (
                "decomposeParDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );

    const label nDomains = readLabel(decompositionDict.lookup("numberOfSubdomains"));

    const word method = decompositionDict.lookup("method");

    if (method == "hierarchical")
    {
        const dictionary& hierarchicalDict = decompositionDict.subDict("hierarchicalCoeffs");

        const Vector<label> n = hierarchicalDict.lookup("n");

        if (nDomains != n.x()*n.y()*n.z())
        {
            FatalErrorInFunction
                << "Wrong number of processor divisions:" << nl
                << "Number of domains    : " << nDomains << nl
                << "Wanted decomposition : " << n
                << exit(FatalError);
        }

        const word order(hierarchicalDict.lookup("order"));

        if (order.size() != 3)
        {
            FatalErrorInFunction
                << "decomposeParDict.hierarchicalCoeffs: "
                << "Number of characters in order (" << order << ") != 3"
                << exit(FatalError);
        }

        for (label i = 0; i < 3; ++i)
        {
            if
            (
                order[i] != 'x'
             && order[i] != 'y'
             && order[i] != 'z'
            )
            {
                FatalErrorInFunction
                    << "decomposeParDict.hierarchicalCoeffs: "
                    << "Illegal decomposition order " << order << endl
                    << "It should only contain x, y or z" << exit(FatalError);
            }
        }
    }
    else if (method != "balancedCellLevels" && method != "globalBalance")
    {
        //should not get here. Already checked in construction of refinementParameters
        FatalErrorInFunction
            << "Unknown decompositionMethod "
            << method << nl << nl
            << "Valid decompositionMethods for VDB are : " << endl
            << "\nbalancedCellLevels"
            << "\nglobalBalance"
            << "\nhierarchical"
            << exit(FatalError);
    }

    nDomainsPerNode =
        labelList
        (
            Pstream::nProcs(),
            nDomains / Pstream::nProcs()
        );

    label remainder = nDomains % Pstream::nProcs();

    if (remainder > 0)
    {
        // skip master
        for (label i = 1; i <= remainder; ++i)
        {
            nDomainsPerNode[i]++;
        }
    }

    // create processor folders
    procOffset = 0;

    forAll(nDomainsPerNode, i)
    {
        if (i == Pstream::myProcNo()) break;

        procOffset += nDomainsPerNode[i];
    }

    for (label i = 0; i < nDomainsPerNode[Pstream::myProcNo()]; ++i)
    {
        label procNo = procOffset + i;

        const word procDir("processor" + Foam::name(procNo));

        mkDir(procDir);
    }

    return decompositionDict;
}


Foam::triSurface Foam::foamVDB::distributeSurface
(
    const triSurface& globalSurf,
    const labelList& refinementLevels
)
{
    Timer timer("distributeSurface");

    labelList faceMap;

    surfacePatchList patches(globalSurf.calcPatches(faceMap));

    scalarField weightedArea(faceMap.size());

    forAll(faceMap, i)
    {
        const label facei = faceMap[i];
        const label regioni = globalSurf[facei].region();

        const label mul = Foam::pow(2, refinementLevels[regioni]);
        // face area * refinement^2
        weightedArea[i] = globalSurf[facei].mag(globalSurf.points()) * mul;// * mul;
    }

    const scalar overAllArea = sum(weightedArea);
    const scalar wantedArea = overAllArea / Pstream::nProcs();

    //if (debug)
    {
        Info<< "overAllArea " << overAllArea
            << " - wantedArea per proc " << wantedArea << endl;
    }

    surfacePatchList weightedPatches(Pstream::nProcs());

    label proc = 0;
    scalar area = 0;
    label startFacei = 0;
    label patchSize = 0;

    forAll(weightedArea, i)
    {
        area += weightedArea[i];

        weightedPatches[proc].size()++;

        if (area > wantedArea || i == weightedArea.size() - 1)
        {
            //Info<< "Area proc" << proc << ": " << area << endl;

            weightedPatches[proc].index() = proc;
            weightedPatches[proc].name() = word("proc") + Foam::name(proc);
            weightedPatches[proc].start() = startFacei;

            proc++;
            startFacei = i + 1;
            // reset scanned area
            area = 0;
        }
    }

    // every processor subsets its faces
    boolList includeMap(globalSurf.size(), false);

    const label& start = weightedPatches[Pstream::myProcNo()].start();

    for
    (
        label i = 0;
        i < weightedPatches[Pstream::myProcNo()].size();
        ++i
    )
    {
        const label facei = faceMap[i + start];
        includeMap[facei] = true;
    }

    labelList pMap, fMap;

    triSurface procSurf = globalSurf.subsetMesh(includeMap, pMap, fMap);

    //if (debug)
    {
        const fileName surfName =
            word("surfProc") + Foam::name(Pstream::myProcNo()) + word(".obj");
        Pout<< "Writing " << surfName << endl;
        procSurf.write(surfName);
        //Info<< "weightedPatches " << weightedPatches << endl;
    }

    return procSurf;
} // distributeSurface


Foam::triSurface Foam::foamVDB::featureTriSurface
(
    const triSurface& surf,
    const scalar featAngle,
    const word name
)
{
    //Timer timer("surfaceFeatures " + name);

    struct CreateFeatureTri
    {
        const surfaceFeatures& features_;
        const label externalStart_;
        const triSurface& surf_;
        List<labelledTri>& featTriangles_;

        CreateFeatureTri
        (
            const surfaceFeatures& features,
            const triSurface& surf,
            List<labelledTri>& featTriangles
        )
        :
            features_(features),
            externalStart_(features_.externalStart()),
            surf_(surf),
            featTriangles_(featTriangles)
        {}

        void operator()(label edgeI) const
        {
            const label& feIndex = features_.featureEdges()[edgeI + externalStart_];
            const edge& fe = surf_.edges()[feIndex];

            featTriangles_[edgeI] =
                labelledTri
                (
                    fe.start(),
                    fe.end(),
                    fe.start()
                );
        }
    };

    surfaceFeatures features
    (
        surf,
        180.0 - featAngle
    );

    // create zero area triangles from feature edges as input for openvdb
    List<labelledTri> featTriangles(features.featureEdges().size() - features.nRegionEdges());

    CreateFeatureTri op
    (
        features,
        surf,
        featTriangles
    );

    parallelForAll(featTriangles, op);

    triSurface featSurf(featTriangles, surf.localPoints());

    return featSurf;
} // featureTriSurface


void Foam::foamVDB::topDownAdvect
(
    FloatGrid::Ptr grid,
    const label cellLevel,
    label maxLevel /*= -1*/
)
{
    if (maxLevel == -1) maxLevel = maxCellLevel_;

    scalar delta = 0;

    // cumulative offset because all the coarse grids have been generated
    // at once with topDownRestrict
    for (label i = maxLevel; i > cellLevel; i--)
    {
        delta -= 0.5 * Foam::pow(2, maxLevel - i);
    }

    delta *= voxelSize_ * Foam::pow(2, maxCellLevel_ - maxLevel);

    advectLevelSet(grid, delta, true);
}

void Foam::foamVDB::bottomUpAdvect
(
    FloatGrid::Ptr grid,
    const label cellLevel
)
{
    scalar delta = 0.5 * Foam::pow(2, maxCellLevel_ - cellLevel);

    //scalar delta = 0;
    //// cumulative offset because all the fine grids have been generated
    //// at once with bottomUpProlongate
    //for (label i = 1; i <= cellLevel; ++i)
    //{
    //    delta += 0.5 * Foam::pow(2, maxCellLevel_ - i);
    //}

    delta *= voxelSize_;
    //Info<<"GGG  cellLevel " << cellLevel << " - delta " << delta<<endl;

    advectLevelSet(grid, delta);
}

void Foam::foamVDB::advectLevelSet
(
    FloatGrid::Ptr grid,
    const scalar delta,
    const bool fromFineGrid /*= false*/
)
{
    Timer timer("advectLevelSet");

    grid->setGridClass(openvdb::GRID_LEVEL_SET);

    //if (!fromFineGrid)
    //{
    //    openvdb::tools::erodeActiveValues(grid->tree(), 2);
    //}

    using VectT   = openvdb::Vec3fGrid;
    using FieldT  = openvdb::tools::DiscreteField<VectT>;
    using AdvectT = openvdb::tools::LevelSetAdvection<FloatGrid, FieldT>;

    VectT vect(openvdb::Vec3f(delta, delta, delta));

    FieldT field(vect);

    AdvectT advect(*grid, field);

    advect.setSpatialScheme(openvdb::math::HJWENO5_BIAS);

    advect.setTemporalScheme(openvdb::math::TVD_RK2);

    advect.advect(0, 1);

    if (fromFineGrid)
    {
        //resize narrowband halfwidth
        using TrackerT = openvdb::tools::LevelSetTracker<FloatGrid>;
        TrackerT tracker(*grid);
        tracker.setSpatialScheme(openvdb::math::FIRST_BIAS);
        tracker.setTemporalScheme(openvdb::math::TVD_RK1);

        tracker.resize(/*halfWidth*/3);
    }
} // advectLevelSet

const openvdb::Coord Foam::foamVDB::toGlobal
(
    const openvdb::Coord& xyzLocal,
    openvdb::math::Transform::ConstPtr transform
)
{
    const openvdb::Vec3d xyz =
        transform->indexToWorld(xyzLocal);

    const openvdb::Coord xyzCoord
    (
        xyz.x(),
        xyz.y(),
        xyz.z()
    );

    return xyzCoord;
}

const openvdb::Coord Foam::foamVDB::toLocal
(
    const openvdb::Coord& xyzGlobal,
    openvdb::math::Transform::ConstPtr transform
)
{
    const openvdb::Vec3d xyz =
        transform->worldToIndex(xyzGlobal.asVec3d());

    const openvdb::Coord xyzCoord
    (
        xyz.x(),
        xyz.y(),
        xyz.z()
    );

    return xyzCoord;
}

void Foam::foamVDB::splitVoxelFace
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
)
{
    const openvdb::Coord coarseFaceCentre = toGlobal(coarseFaceCoord, coarseTransform);

    for (unsigned pointI=0; pointI < 4; pointI++)
    {
        label next     = (pointI == 3) ? 0 : (pointI + 1);
        label previous = (pointI == 0) ? 3 : (pointI - 1);

        openvdb::Coord midNext
        (
            (points[f[pointI]].x() + points[f[next]].x()) / 2.,
            (points[f[pointI]].y() + points[f[next]].y()) / 2.,
            (points[f[pointI]].z() + points[f[next]].z()) / 2.
        );

        const openvdb::Coord midNextLocalCoord = toLocal(midNext, fineTransform);

        auto midNextIndex = fineIndexGridAccessor.getValue(midNextLocalCoord);

        //add midNext point if not already added
        if (midNextIndex < 0)
        {
            midNextIndex = pointCount.fetch_add(1);

            fineIndexGridAccessor.setValue(midNextLocalCoord, midNextIndex);
            fineAddedIndexesGridAccessor.setValue(midNextLocalCoord, midNextIndex);

            //if (debug)
            //{
            //    std::cout<<"Added new point: " << midNextLocalCoord
            //        << " (midNext) global coord: " << midNext
            //        << " at index " << midNextIndex << std::endl;
            //}
        }

        openvdb::Coord midPrevious
        (
            (points[f[pointI]].x() + points[f[previous]].x()) / 2.,
            (points[f[pointI]].y() + points[f[previous]].y()) / 2.,
            (points[f[pointI]].z() + points[f[previous]].z()) / 2.
        );

        const openvdb::Coord midPreviousLocalCoord = toLocal(midPrevious, fineTransform);

        auto midPreviousIndex = fineIndexGridAccessor.getValue(midPreviousLocalCoord);

        //add midPrevious point if not already added
        if (midPreviousIndex < 0)
        {
            midPreviousIndex = pointCount.fetch_add(1);

            fineIndexGridAccessor.setValue(midPreviousLocalCoord, midPreviousIndex);
            fineAddedIndexesGridAccessor.setValue(midPreviousLocalCoord, midPreviousIndex);

            //if (debug)
            //{
            //    std::cout<<"Added new point: " << midPreviousLocalCoord
            //        << " (midPrevious) global coord: " << midPrevious
            //        << " at index " << midPreviousIndex << std::endl;
            //}
        }

        const openvdb::Coord newFaceCentre
        (
            (coarseFaceCentre.x() + 2*points[f[pointI]].x()) / 2.,
            (coarseFaceCentre.y() + 2*points[f[pointI]].y()) / 2.,
            (coarseFaceCentre.z() + 2*points[f[pointI]].z()) / 2.
        );

        const openvdb::Coord newFaceCentreLocal = toLocal(newFaceCentre, fineTransform);

        openvdb::Coord xyzCoarseFaceCentre
        (
            coarseFaceCentre.x() / 2.,
            coarseFaceCentre.y() / 2.,
            coarseFaceCentre.z() / 2.
        );

        const openvdb::Coord xyzCoarseFaceCentreLocal = toLocal(xyzCoarseFaceCentre, fineTransform);

        openvdb::Int32 coarseFaceCentreIndex = fineIndexGridAccessor.getValue(xyzCoarseFaceCentreLocal);

        if (coarseFaceCentreIndex < 0)
        {
            coarseFaceCentreIndex = pointCount.fetch_add(1);

            fineIndexGridAccessor.setValue(xyzCoarseFaceCentreLocal, coarseFaceCentreIndex);
            fineAddedIndexesGridAccessor.setValue(xyzCoarseFaceCentreLocal, coarseFaceCentreIndex);
        }

        bool addNewFace =
            forceAdd
          ? true
          : !cellLevelInterfaceFaceGridAccessor.isValueOn(newFaceCentreLocal);

        // add newFace as boundary face if it does not exist yet
        if (addNewFace)
        {
            openvdb::math::Vec4ui newFace
            (
                f[pointI],
                midNextIndex,
                coarseFaceCentreIndex,
                midPreviousIndex
            );

            auto coarseOwner = coarseOwnerGridAccessor.getValue(coarseFaceCoord);

            auto coarsePatchID = coarsePatchIDGridAccessor.getValue(coarseFaceCoord);

            fineOwnerGridAccessor.setValueOn(newFaceCentreLocal, coarseOwner);

            finePatchIDGridAccessor.setValueOn(newFaceCentreLocal, coarsePatchID);

            fineFaceGridAccessor.setValueOn(newFaceCentreLocal, newFace);
        }
    } // for pointI < 4
} // splitVoxelFace

//void Foam::foamVDB::checkPatchVoxelsProximity
//(
//    const std::vector<std::vector<FloatGrid::Ptr>>& boundaryCellLevels
//)
//{
//    Timer timer("checkPatchVoxelProximity");
//
//    struct VoxelProximityOp
//    {
//        FloatGrid& surfGrid_;
//        //FloatGrid& otherSurfGrid_;
//        //FloatGrid::Accessor surfGridAcc_;
//        FloatGrid::ConstAccessor otherSurfGridAcc_;
//
//        VoxelProximityOp
//        (
//            FloatGrid& surfGrid,
//            FloatGrid& otherSurfGrid
//        )
//        :
//            surfGrid_(surfGrid),
//            //otherSurfGrid_(otherSurfGrid)
//            //surfGridAcc_(surfGrid.getAccessor()),
//            otherSurfGridAcc_(otherSurfGrid.getConstAccessor())
//        {}
//
//        inline void operator()(const FloatGrid::ValueOnCIter& iter) const
//        {
//            FloatGrid::Accessor surfGridAcc(surfGrid_.getAccessor());
//            //FloatGrid::ConstAccessor otherSurfGridAcc(otherSurfGrid_.getConstAccessor());
//
//            if
//            (
//                surfGridAcc.getValue(iter.getCoord())
//              > otherSurfGridAcc_.getValue(iter.getCoord())
//            )
//            {
//                surfGridAcc.setValueOff(iter.getCoord());
//                //surfGridAcc.setValue(iter.getCoord(), otherSurfGridAcc_.getValue(iter.getCoord())); //GGG
//            }
//        }
//    };
//
//    for (size_t patchI = 0; patchI < boundaryCellLevels.size(); patchI++)
//    {
//        for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
//        {
//            if (!boundaryCellLevels[patchI][cellLevel]) continue;
//
//            FloatGrid::Ptr subSurfGrid = boundaryCellLevels[patchI][cellLevel];
//
//            //if (debug)
//            //{
//            //    std::cout << "Neighbour patches of " << subSurfGrid->getName() <<": " <<std::endl;
//            //}
//            for (size_t patchJ = 0; patchJ < boundaryCellLevels.size(); patchJ++)
//            {
//                if (patchI == patchJ) continue;
//
//                if (!boundaryCellLevels[patchJ][cellLevel]) continue;
//
//                FloatGrid::Ptr otherSubSurfGrid = boundaryCellLevels[patchJ][cellLevel];
//
//                FloatGrid::Ptr intersection = otherSubSurfGrid->deepCopy();
//
//                intersection->topologyIntersection(*subSurfGrid);
//
//                if (intersection->activeVoxelCount() > 0)
//                {
//                    //if (debug)
//                    //{
//                    //    std::cout<< "    " << otherSubSurfGrid->getName() << std::endl;
//                    //}
//                    openvdb::tools::foreach
//                    (
//                        intersection->cbeginValueOn(),
//                        VoxelProximityOp
//                        (
//                            *subSurfGrid,
//                            *otherSubSurfGrid
//                        ),
//                        /*threaded*/true,
//                        /*shareOp*/false //true
//                    );
//                }
//            } // for patchJ
//            //if (debug)
//            //{
//            //    std::cout<< std::endl;
//            //}
//        } // for cellLevel
//    }//for patchI
//} // checkPatchVoxelsProximity


void Foam::foamVDB::sliceBoundBox
(
    treeBoundBox& bb,
    const label coarserCellLevel,
    vector& axis
)
{
    const vector& span = bb.span();

    //slice along longest axis
    label sliceDir = span[0] > span[1] ? 0 : 1;
    sliceDir = span[2] > span[sliceDir] ? 2 : sliceDir;

    scalar delta = span[sliceDir] / Pstream::nProcs();

    scalar surfVoxelSize = voxelSize_ * Foam::pow(2, maxCellLevel_ - coarserCellLevel);

    if (delta <= surfVoxelSize)
    {
        Pout<<"WARNING!! bounding box delta ("
            << delta
            << ") < surface voxelSize ("
            << surfVoxelSize
            << ")" << endl;
    }

    point& bbMin = bb.min();
    point& bbMax = bb.max();

    bbMin[sliceDir] =
        Foam::max
        (
            bbMin[sliceDir],
            bbMin[sliceDir] + delta*Pstream::myProcNo() - surfVoxelSize
        );

    bbMax[sliceDir] =
        Foam::min
        (
            bbMax[sliceDir],
            bbMin[sliceDir] + delta + surfVoxelSize
        );

    axis =
        vector
        (
            sliceDir == 0 ? 1 : 0,
            sliceDir == 1 ? 1 : 0,
            sliceDir == 2 ? 1 : 0
        );
} // sliceBoundBox


const Foam::triSurface Foam::foamVDB::getVDBCuttingPlane
(
    const edgeList& edges,
    const scalarField& dotProducts,
    const pointField& points,
    const plane& cutPlane
)
{
    DynamicList<point> cutPoints(points.size()/100);

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if
        (
            (dotProducts[e[0]] <= 0.0 && dotProducts[e[1]] > 0.0)
         || (dotProducts[e[1]] <= 0.0 && dotProducts[e[0]] > 0.0)
        )
        {
            const point& p0 = points[e[0]];
            const point& p1 = points[e[1]];

            scalar alpha = cutPlane.lineIntersect(linePointRef(p0, p1));

            if (alpha < SMALL)
            {
                cutPoints.append(p0);
            }
            else if (alpha >= (1.0 - SMALL))
            {
                cutPoints.append(p1);
            }
            else
            {
                cutPoints.append((1-alpha)*p0 + alpha*p1);
            }
        }
    } //forAll edges
    cutPoints.shrink();

    triSurface planeSurf;

    if (cutPoints.size())
    {
        treeBoundBox planeBB(cutPoints);

        List<labelledTri> triangles(12);

        //Vertex coordinates. In octant coding
        triangles[0]  = labelledTri(0, 1, 2);
        triangles[1]  = labelledTri(1, 3, 2);
        triangles[2]  = labelledTri(2, 3, 6);
        triangles[3]  = labelledTri(3, 7, 6);
        triangles[4]  = labelledTri(6, 7, 4);
        triangles[5]  = labelledTri(7, 5, 4);
        triangles[6]  = labelledTri(4, 5, 0);
        triangles[7]  = labelledTri(5, 1, 0);
        triangles[8]  = labelledTri(4, 0, 6);
        triangles[9]  = labelledTri(0, 2, 6);
        triangles[10] = labelledTri(1, 5, 3);
        triangles[11] = labelledTri(5, 7, 3);

        planeSurf = triSurface(triangles, planeBB.points());

        planeSurf.cleanup(/*verbose*/false);
    }
    else
    {
        pointField pf(1, cutPlane.refPoint());

        List<labelledTri> triangles(1, labelledTri(0,0,0));

        planeSurf = triSurface(triangles, pf);
    }

    return planeSurf;
} //getVDBCuttingPlane


void Foam::foamVDB::clipClosedSurface
(
    triSurface& surf,
    treeBoundBox& surfBb,
    const label coarserCellLevel
)
{
    Timer timer("clipClosedSurface");

    List<plane> boundPlanes = boundBoxToPlanes(surfBb);

    const edgeList& edges = surf.edges();
    const pointField& points = surf.localPoints();

    label nPlaneFaces  = 0;
    label nPlanePoints = 0;

    List<labelledTri> allFaces(surf);
    vectorField allPoints(surf.points());

    forAll(boundPlanes, I)
    {
        const plane& planeI = boundPlanes[I];

        DynamicList<point> cutPoints(surf.localPoints().size()/100);

        const scalarField dotProducts
            ((points - planeI.refPoint()) & planeI.normal());

        triSurface triPlane =
            getVDBCuttingPlane
            (
                edges,
                dotProducts,
                points,
                planeI
            );

        //add this plane faces to allFaces
        forAll(triPlane, faceI)
        {
            allFaces.append
            (
                labelledTri
                (
                    triPlane[faceI][0] + allPoints.size(),
                    triPlane[faceI][1] + allPoints.size(),
                    triPlane[faceI][2] + allPoints.size(),
                    surf.patches().size()
                )
            );
        }

        allPoints.append(triPlane.points());
    }

    geometricSurfacePatchList newPatches(surf.patches());

    newPatches.append
    (
        geometricSurfacePatch
        (
            "VDBPlane",
            surf.patches().size()
        )
    );

    surf = triSurface(allFaces, newPatches, allPoints);

    //triSurfaceSearch triSearch1(bbox);
    //triSurfaceSearch triSearch2(surf);

    //surfaceIntersection inter(triSearch1, triSearch2);

    //booleanSurface decomposedSurf(bbox, surf, inter, booleanSurface::INTERSECTION);

    //surf = decomposedSurf;

    //if (debug)
    {
        surf.write("closedInputSurf_proc" + std::to_string(Pstream::myProcNo()) + ".obj");
        surf.writeStats(Pout);
    }
} //clipClosedSurface

void Foam::foamVDB::decomposeVDBSurface
(
    triSurface& surf,
    const label coarserCellLevel
)
{
    Timer timer("DecomposeVDBSurface");

    treeBoundBox surfBBox(surf.localPoints());

    vector normal;

    sliceBoundBox
    (
        surfBBox,
        coarserCellLevel,
        normal
    );

    point& bbMin = surfBBox.min();
    point& bbMax = surfBBox.max();

    plane minPlane(bbMin, normal);
    plane maxPlane(bbMax, normal);

    const edgeList& edges = surf.edges();
    const pointField& points = surf.localPoints();

    DynamicList<point> cutPoints(surf.localPoints().size()/100);

    const scalarField dotProductsMinPlane
        ((points - minPlane.refPoint()) & minPlane.normal());

    const scalarField dotProductsMaxPlane
        ((points - maxPlane.refPoint()) & maxPlane.normal());

    triSurface minSurf =
        getVDBCuttingPlane
        (
            edges,
            dotProductsMinPlane,
            points,
            minPlane
        );

    triSurface maxSurf =
        getVDBCuttingPlane
        (
            edges,
            dotProductsMaxPlane,
            points,
            maxPlane
        );

    boolList includeMap(surf.size(), true);

    forAll(surf, faceI) //TODO tbb parallel loop
    {
        const labelledTri& f = surf[faceI];

        // exclude triangle if all of its points are outside bounding box
        label n = 0;
        forAll(f, pointI)
        {
            const point& p = surf.points()[f[pointI]];

            if
            (
                ((p - bbMin) & normal) < 0
             || ((p - bbMax) & normal) > 0
            )
            {
                n++;
            }
        } // forAll triangle points

        if (n == 3)
        {
            includeMap[faceI] = false;
        }
    } // forAll triangles

    // Subset triSurface
    labelList pointMap;
    labelList faceMap;

    triSurface subSurf =
        surf.subsetMesh
        (
            includeMap,
            pointMap,
            faceMap
        );

    // Make new storage
    List<labelledTri> facesAll(subSurf.size() + minSurf.size() + maxSurf.size());

    const pointField& points1 = subSurf.points();
    const pointField& points2 = minSurf.points();
    const pointField& points3 = maxSurf.points();

    vectorField pointsAll(points1.size() + points2.size() + points3.size());

    label pointi = 0;

    // Copy points1 into pointsAll
    forAll(points1, point1i)
    {
        pointsAll[pointi++] = points1[point1i];
    }
    // Add surface2 points
    forAll(points2, point2i)
    {
        pointsAll[pointi++] = points2[point2i];
    }
    // Add surface3 points
    forAll(points3, point3i)
    {
        pointsAll[pointi++] = points3[point3i];
    }

    label trianglei = 0;

    // Determine map for both regions
    labelList patch1Map(subSurf.patches().size());

    patch1Map = identity(subSurf.patches().size());

    label nNewPatches = patch1Map.size() + 1;

    // Copy triangles1 into trianglesAll
    forAll(subSurf, facei)
    {
        const labelledTri& tri = subSurf[facei];
        labelledTri& destTri = facesAll[trianglei++];

        destTri.triFace::operator=(tri);
        destTri.region() = patch1Map[tri.region()];
    }

    // Add (renumbered) surface2 triangles
    forAll(minSurf, facei)
    {
        const labelledTri& tri = minSurf[facei];

        labelledTri& destTri = facesAll[trianglei++];
        destTri[0] = tri[0] + points1.size();
        destTri[1] = tri[1] + points1.size();
        destTri[2] = tri[2] + points1.size();
        destTri.region() = patch1Map.size();
    }

    // Add (renumbered) surface3 triangles
    forAll(maxSurf, facei)
    {
        const labelledTri& tri = maxSurf[facei];

        labelledTri& destTri = facesAll[trianglei++];
        destTri[0] = tri[0] + points1.size() + points2.size();
        destTri[1] = tri[1] + points1.size() + points2.size();
        destTri[2] = tri[2] + points1.size() + points2.size();
        destTri.region() = patch1Map.size();
    }

    geometricSurfacePatchList newPatches(nNewPatches);

    forAll(subSurf.patches(), patchi)
    {
        newPatches[patch1Map[patchi]] = subSurf.patches()[patchi];
    }

    newPatches[patch1Map.size()] =
        geometricSurfacePatch
        (
            "VDBPlane",
            patch1Map.size()
        );

    // Construct new surface mesh
    surf = triSurface(facesAll, newPatches, pointsAll);

    //if (debug)
    {
        surf.write("VDBinputSurf_proc" + std::to_string(Pstream::myProcNo()) + ".obj");
        surf.writeStats(Pout);
    }
} //decomposeVDBSurface


//const Foam::List<Foam::triSurface> Foam::foamVDB::boundBoxToPlanes
const Foam::List<Foam::plane> Foam::foamVDB::boundBoxToPlanes
(
    const treeBoundBox& box
)
{
    pointField& boxPoints = box.points().ref();

    List<plane> planes(6);
    //List<triSurface> planes(6);

    //List<labelledTri> triangles(2);

    ////Vertex coordinates. In octant coding
    //triangles[0]  = labelledTri(0, 1, 2);
    //triangles[1]  = labelledTri(1, 3, 2);

    //pointField points(4);

    //bottom face
    //points[0] = boxPoints[0];
    //points[1] = boxPoints[1];
    //points[2] = boxPoints[2];
    //points[3] = boxPoints[3];
    //planes[0] = triSurface(triangles, points);
    planes[0] = plane(boxPoints[0], boxPoints[1], boxPoints[2]);

    //left face
    //points[0] = boxPoints[4];
    //points[1] = boxPoints[0];
    //points[2] = boxPoints[6];
    //points[3] = boxPoints[2];
    //planes[1] = triSurface(triangles, points);
    //
    //for some reason boxPoint[0][0] is always VERYSMALL ...I'll avoid that for now
    //boxPoints 8(( 1.36896e-316 -1.01496 -0.296862 ) ( 1.07265 -1.01496 -0.296862 ) ( 0.117375 1.01496 -0.296862 ) ( 1.07265 1.01496 -0.296862 ) ( 0.117375 -1.01496 0.965872 ) ( 1.07265 -1.01496 0.965872 ) ( 0.117375 1.01496 0.965872 ) ( 1.07265 1.01496 0.965872 ))
    //planes[1] = plane(boxPoints[4], boxPoints[0], boxPoints[6]);
    planes[1] = plane(boxPoints[4], boxPoints[2], boxPoints[6]);

    //right face
    //points[0] = boxPoints[1];
    //points[1] = boxPoints[5];
    //points[2] = boxPoints[3];
    //points[3] = boxPoints[7];
    //planes[2] = triSurface(triangles, points);
    planes[2] = plane(boxPoints[1], boxPoints[5], boxPoints[3]);

    //top face
    //points[0] = boxPoints[6];
    //points[1] = boxPoints[7];
    //points[2] = boxPoints[4];
    //points[3] = boxPoints[5];
    //planes[3] = triSurface(triangles, points);
    planes[3] = plane(boxPoints[6], boxPoints[7], boxPoints[4]);

    //front face
    //points[0] = boxPoints[4];
    //points[1] = boxPoints[5];
    //points[2] = boxPoints[0];
    //points[3] = boxPoints[1];
    //planes[4] = triSurface(triangles, points);
    planes[4] = plane(boxPoints[4], boxPoints[5], boxPoints[0]);

    //back face
    //points[0] = boxPoints[2];
    //points[1] = boxPoints[3];
    //points[2] = boxPoints[6];
    //points[3] = boxPoints[7];
    //planes[5] = triSurface(triangles, points);
    planes[5] = plane(boxPoints[2], boxPoints[3], boxPoints[6]);

    return planes;
} //boundBoxToPlanes


const Foam::triSurface Foam::foamVDB::boundBoxToTriSurface
(
    treeBoundBox& box,
    const label shellLevel /*= -1*/,
    const bool decompose /*= false*/
)
{
    List<labelledTri> triangles(12);

    //Vertex coordinates. In octant coding
    triangles[0]  = labelledTri(0, 1, 2);
    triangles[1]  = labelledTri(1, 3, 2);
    triangles[2]  = labelledTri(2, 3, 6);
    triangles[3]  = labelledTri(3, 7, 6);
    triangles[4]  = labelledTri(6, 7, 4);
    triangles[5]  = labelledTri(7, 5, 4);
    triangles[6]  = labelledTri(4, 5, 0);
    triangles[7]  = labelledTri(5, 1, 0);
    triangles[8]  = labelledTri(4, 0, 6);
    triangles[9]  = labelledTri(0, 2, 6);
    triangles[10] = labelledTri(1, 5, 3);
    triangles[11] = labelledTri(5, 7, 3);

    if (Pstream::parRun() && decompose)
    {
        vector normal;

        sliceBoundBox
        (
            box,
            shellLevel,
            normal
        );
    } // if parallel

    triSurface surf(triangles, box.points());

    //if (debug)
    //{
    //    surf.write("refinementBox_proc" + std::to_string(Pstream::myProcNo()) + ".obj");
    //}

    return surf;
}

const Foam::List<Foam::triSurface> Foam::foamVDB::splitSurface
(
    const triSurface& surf,
    List<labelList> sameMinMaxLevel /*= List<labelList>(0)*/
)
{
    //Timer timer("splitSurface");

    labelList faceMap;
    surfacePatchList patches(surf.calcPatches(faceMap));

    label nRegions = patches.size();

    if (sameMinMaxLevel.size())
    {
        nRegions = sameMinMaxLevel.size();
    }

    List<triSurface> subSurfs(nRegions);

    for (label patchI = 0; patchI < nRegions; patchI++)
    {
        boolList includeMap(surf.size(), false);

        if (sameMinMaxLevel.size())
        {
            forAll(surf, faceI)
            {
                const labelledTri& f = surf[faceI];

                forAll(sameMinMaxLevel[patchI], regioni)
                {
                    if (f.region() == sameMinMaxLevel[patchI][regioni])
                    {
                        includeMap[faceI] = true;
                    }
                }
            }
        }
        else if (patches[patchI].size())
        {
            const label& start = patches[patchI].start();

            for
            (
                label i = 0;
                i < patches[patchI].size();
                ++i
            )
            {
                const label facei = faceMap[i + start];
                includeMap[facei] = true;
            }
        }
        else
        {
            subSurfs[patchI] =
                triSurface
                (
                    List<labelledTri>(0),
                    surf.patches(),
                    pointField(0)
                );

            continue;
        }

        // Subset triSurface
        labelList pointMap;
        labelList faceMap0;

        subSurfs[patchI] =
            triSurface
            (
                surf.subsetMesh
                (
                    includeMap,
                    pointMap,
                    faceMap0
                )
            );
    }

    //if (debug)
    //{
    //    forAll(subSurfs, i)
    //    {
    //        subSurfs[i].write
    //        (
    //            word("subSurf") + Foam::name(i) + word("_proc") + Foam::name(Pstream::myProcNo()) + word(".obj")
    //        );
    //    }
    //}

    return subSurfs;
} //splitSurface


const Foam::triSurface Foam::foamVDB::addSurfaces
(
    List<triSurface>& surfs
)
{
    List<labelledTri> faces;
    pointField points;
    DynamicList<geometricSurfacePatch> patches(16*surfs.size());

    forAll(surfs, surfi)
    {
        triSurface addsurf = surfs[surfi];

        List<labelledTri> addfaces(addsurf.xferFaces());
        List<point> addpoints(addsurf.xferPoints());

        // region offset
        const label regoff = patches.size();
        patches.append(addsurf.patches());

        // Offset the points for all additional surfaces
        if (surfi)
        {
            const label ptoff = points.size();

            for (labelledTri& f : addfaces)
            {
                forAll(f, fi)
                {
                    f[fi] += ptoff;
                }

                f.region() += regoff;
            }

            faces.append(addfaces);
            points.append(addpoints);
        }
        else
        {
            faces.transfer(addfaces);
            points.transfer(addpoints);
        }
    }

    triSurface globalSurf
    (
        faces,
        patches,
        points,
        true
    );

    //globalSurf.cleanup(/*verbose*/ false);

    return globalSurf;
} //addSurfaces


const Foam::triSurface Foam::foamVDB::gatherTriSurface
(
    triSurface& inputSurf
)
{
    Timer timer("gatherTriSurface");

    List<triSurface> inputSurfs(Pstream::nProcs());

    inputSurfs[Pstream::myProcNo()] = inputSurf;

    PstreamBuffers pBufSurf(Pstream::commsTypes::nonBlocking);

    // myProcNo sends to all other procs
    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        if (proci == Pstream::myProcNo()) continue;

        //non-blocking communication
        UOPstream toProcSurf(proci, pBufSurf);
        toProcSurf << inputSurfs[Pstream::myProcNo()];
    }

    //Start sending and receiving and block
    pBufSurf.finishedSends();

    // myProcNo receives from all other procs
    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        if (proci == Pstream::myProcNo()) continue;

        UIPstream fromProcSurf(proci, pBufSurf);
        fromProcSurf >> inputSurfs[proci];
    }

    triSurface globalSurf = addSurfaces(inputSurfs);

    //if (debug)
    {
        inputSurf.write
        (
            "gatherTriSurface_proc" + Foam::name(Pstream::myProcNo()) + std::string(".obj")
        );
        globalSurf.write("gatherTriSurface_globalSurf.obj");
    }

    return globalSurf;
} //gatherTriSurface

openvdb::math::BBox<openvdb::Vec3s> Foam::foamVDB::getBoundingBox
(
    const std::vector<FloatGrid::Ptr>& grids
)
{
    openvdb::math::BBox<openvdb::Vec3s> bBox0
    (
        openvdb::Vec3s(0,0,0),
        openvdb::Vec3s(0,0,0)
    );

    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        scalar edgeLength = Foam::pow(2, cellLevel);

        scalar padding = edgeLength - 1;

        openvdb::CoordBBox bBox = grids[cellLevel]->evalActiveVoxelBoundingBox();

        openvdb::math::BBox<openvdb::Vec3s> bBoxToLevel0
        (
            openvdb::Vec3s
            (
                bBox.min().x() / edgeLength,
                bBox.min().y() / edgeLength,
                bBox.min().z() / edgeLength
            ),
            openvdb::Vec3s
            (
                (bBox.max().x() - padding) / edgeLength,
                (bBox.max().y() - padding) / edgeLength,
                (bBox.max().z() - padding) / edgeLength
            )
        );

        bBox0.expand(bBoxToLevel0);

        //std::cout<< "level "<< cellLevel
        //    <<" bBox " << bBox
        //    <<" bBoxToLevel0 " << bBoxToLevel0
        //    <<std::endl;
    }
    return bBox0;
} // getBoundingBox


Foam::labelListList Foam::foamVDB::decomposeGrids
(
    const std::vector<FloatGrid::Ptr>& cellLevelGrids,
    std::vector<IndexGrid::Ptr>&  globalIDGrids,
    std::vector<std::vector<FloatGrid::Ptr>>& decomposedGrids,
    const dictionary& decompositionDict,
    const labelList& nVoxelsStart,
    const pointField& voxelCentres
)
{
    const word method = decompositionDict.lookup("method");

    Timer timer("Decompose grids OF (" + method + ")");

    const label& nProcs = Pstream::nProcs();

    std::vector<std::vector<FloatGrid::Accessor>> decomposedGridAccessors;

    // initialize decomposedGrids
    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        std::vector<FloatGrid::Ptr> levelNprocGrids;
        std::vector<FloatGrid::Accessor> procAccessors;

        for (label i = 0; i < nProcs; ++i)
        {
            FloatGrid::Ptr procM = FloatGrid::create();

            procM->setTransform
            (
                cellLevelGrids[cellLevel]->transform().copy()
            );

            procM->setGridClass(openvdb::GRID_UNKNOWN);

            word gridName = "proc " + Foam::name(i) + " level " + Foam::name(cellLevel);
            procM->setName(gridName);

            levelNprocGrids.push_back(procM);

            procAccessors.push_back
            (
                procM->getAccessor()
            );
        }

        decomposedGrids.push_back(levelNprocGrids);
        decomposedGridAccessors.push_back(procAccessors);
    } // for cellLevel


    autoPtr<decompositionMethod> decomposer =
        decompositionMethod::New
        (
            decompositionDict
        );

    labelList distribution;

    if (isA<geomDecomp>(decomposer()))
    {
        Timer t("decomposer.decompose() " + method);
        distribution =
            decomposer().decompose
            (
                voxelCentres,
                scalarField(voxelCentres.size(), 1) //weights
            );
    }
    else
    {
        Timer t("decomposer.decompose() " + method);
        labelListList globalCellCells =
            calcVoxelVoxels(globalIDGrids, voxelCentres.size());

        distribution =
            decomposer().decompose
            (
                globalCellCells,
                voxelCentres,
                scalarField(globalCellCells.size(), 1) //weights
            );
    } //Timer decomposer

    //nVoxelsInProc[cellLevel][procI];
    labelListList nVoxelsInProc(maxCellLevel_ + 1, labelList(nProcs, 0));

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, maxCellLevel_+1),
        [&](const tbb::blocked_range<size_t>& r)
        {
            for (size_t cellLevel = r.begin(); cellLevel < r.end(); ++cellLevel)
            {
                label nVoxels = nVoxelsStart[cellLevel];

                cellLevelGrids[cellLevel]->tree().voxelizeActiveTiles();

                for (auto iter = cellLevelGrids[cellLevel]->cbeginValueOn(); iter; ++iter)
                {
                    label procI = distribution[nVoxels++];

                    nVoxelsInProc[cellLevel][procI]++;

                    decomposedGridAccessors[cellLevel][procI].setValue
                    (
                        iter.getCoord(),
                        iter.getValue()
                    );
                }
            }
        }
    );

    return nVoxelsInProc;
} //decomposeGrids voxelCentres


void Foam::foamVDB::decomposeGrids
(
    const std::vector<FloatGrid::Ptr>& grids,
    std::vector<std::vector<FloatGrid::Ptr>>& decomposedGrids,
    const dictionary& decompositionDict,
    const label& nProcs,
    const scalar& maxLoadUnbalance
)
{
    const word method = decompositionDict.lookup("method");

    Timer timer("Decompose grids (" + method + ")");

    std::vector<std::vector<FloatGrid::Accessor>> decomposedGridAccessors;

    // initialize decomposedGrids
    for (size_t cellLevel = 0; cellLevel < grids.size(); ++cellLevel)
    {
        std::vector<FloatGrid::Ptr> levelNprocGrids;
        std::vector<FloatGrid::Accessor> procAccessors;

        for (label i = 0; i < nProcs; ++i)
        {
            FloatGrid::Ptr procM = FloatGrid::create();

            procM->setTransform
            (
                grids[cellLevel]->transform().copy()
            );

            procM->setGridClass(openvdb::GRID_UNKNOWN);

            word gridName = "proc " + Foam::name(i) + " level " + Foam::name(cellLevel);
            procM->setName(gridName);

            levelNprocGrids.push_back(procM);

            procAccessors.push_back
            (
                procM->getAccessor()
            );
        }

        decomposedGrids.push_back(levelNprocGrids);
        decomposedGridAccessors.push_back(procAccessors);
    } // for cellLevel

    if (method == "balancedCellLevels")
    {
        for (size_t cellLevel=0; cellLevel < grids.size(); ++cellLevel)
        {
            int nVoxels = grids[cellLevel]->activeVoxelCount();

            const label voxelsPerProc = nVoxels / nProcs;

            const label leftOver = nProcs - (nVoxels % nProcs);

            //openvdb::tools::pruneLevelSet(cellLevelGrids[cellLevel]->tree());

            for (auto iter = grids[cellLevel]->cbeginValueOn(); iter; ++iter)
            {
                for (int i = nProcs - 1; i >= 0; i--)
                {
                    int threshold =
                        i <= leftOver
                      ? i * voxelsPerProc
                      : i * voxelsPerProc + i - leftOver;

                    if (nVoxels > threshold)
                    {
                        decomposedGridAccessors[cellLevel][i].setValue(iter.getCoord(), iter.getValue());
                        break;
                    }
                }
                nVoxels--;
            } // for valueOn
        } // for cellLevel
    }
    else if (method == "globalBalance")
    {
        label nVoxels = getActiveVoxels<FloatGrid>(grids);

        const label voxelsPerProc = nVoxels / nProcs;

        const label leftOver = nProcs - (nVoxels % nProcs);

        for (size_t cellLevel=0; cellLevel < grids.size(); ++cellLevel)
        {
            for (auto iter = grids[cellLevel]->cbeginValueOn(); iter; ++iter)
            {
                for (int i = nProcs - 1; i >= 0; i--)
                {
                    int threshold =
                        i <= leftOver
                      ? i * voxelsPerProc
                      : i * voxelsPerProc + i - leftOver;

                    if (nVoxels > threshold)
                    {
                        decomposedGridAccessors[cellLevel][i].setValue(iter.getCoord(), iter.getValue());
                        break;
                    }
                }
                nVoxels--;
            }
        }
    }
    else if (method == "hierarchical")
    {
        const dictionary& hierarchicalDict = decompositionDict.subDict("hierarchicalCoeffs");

        const Vector<label> n = hierarchicalDict.lookup("n");

        FixedList<direction, 3> decompOrder;

        const word order(hierarchicalDict.lookup("order"));

        for (label i = 0; i < 3; ++i)
        {
            if (order[i] == 'x')
            {
                decompOrder[i] = 0;
            }
            else if (order[i] == 'y')
            {
                decompOrder[i] = 1;
            }
            else if (order[i] == 'z')
            {
                decompOrder[i] = 2;
            }
        }

        doDecomposition
        (
            grids,
            decomposedGrids,
            n,
            decompOrder,
            0,  // nthDir
            0,  // procOffset
            nProcs,
            maxLoadUnbalance
        );
    }
    else //should not get here as checked already in createDecompDict
    {
        FatalErrorInFunction
            << "Unknown decompositionMethod "
            << method << nl << nl
            << "Supported decompositionMethods "
            << "for VDB grids are:"
            << "\n  balancedCellLevels"
            << "\n  globalBalance"
            << "\n  hierarchical"
            << exit(FatalError);
    }
} //decomposeGrids


std::vector<FloatGrid::Ptr> Foam::foamVDB::getInteriorGrids
(
    const std::vector<FloatGrid::Ptr>& narrowBandLevelSets,
    const label maxShellLevel
)
{
    struct SetOutsideValueOff
    {
        static inline void op(const FloatGrid::ValueOnCIter& iter, FloatGrid::Accessor& accessor)
        {
            if (iter.getValue() > -VSMALL)
            {
                accessor.setValueOff(iter.getCoord());
            }
        }
    };

    const label& cellLevels = narrowBandLevelSets.size();

    std::vector<FloatGrid::Ptr> interiorGrids(cellLevels, nullptr);

    for (label cellLevel=0; cellLevel < cellLevels; ++cellLevel)
    {
        if (!narrowBandLevelSets[cellLevel]) continue;

        Timer timer;

        FloatGrid::Ptr interiorGrid = narrowBandLevelSets[cellLevel]->deepCopy();

        if (cellLevel <= maxShellLevel)
        {
            auto interior = openvdb::tools::sdfInteriorMask(*interiorGrid);

            interior->tree().voxelizeActiveTiles();

            interiorGrid->topologyUnion(*interior);

            openvdb::tools::erodeActiveValues(interiorGrid->tree(), 1); //erode same number of narrowBand bandWidth
            ////openvdb::tools::erodeActiveValues(interiorGrid->tree(), 3);
        }
        else
        {
            //TODO try foreach as trasnformValues threaded not working as expected
            openvdb::tools::transformValues
            (
                narrowBandLevelSets[cellLevel]->cbeginValueOn(),
                *interiorGrid,
                SetOutsideValueOff::op,
                /*threaded*/false,//true,
                /*shareOp*/false//true
                /*merge=MERGE_ACTIVE_STATES*/
            );

            openvdb::tools::dilateActiveValues(interiorGrid->tree(), 1);
        }

        //if (debug)
        {
            Info<< "Geometry interior level " << cellLevel
                << ": active voxels ";
            if (Pstream::parRun())
            {
                List<label> procActiveVoxels(Pstream::nProcs());
                procActiveVoxels[Pstream::myProcNo()] = interiorGrid->activeVoxelCount();
                IPstream::gatherList(procActiveVoxels);
                Info<< procActiveVoxels << endl;
            }
            else
            {
                Info<< interiorGrid->activeVoxelCount() <<endl;
            }
        }

        //interiorGrid->setName(gridName);
        interiorGrids[cellLevel] = interiorGrid;
    }
    return interiorGrids;
} //getInteriorGrids


// Data structure to convert from OpenFOAM triSurface to OpenVDB format
struct MeshDataAdapter
{
    const Foam::triSurface surf_;

    MeshDataAdapter
    (
        const Foam::triSurface& surf
    )
    :
        surf_(surf)
    {}

    size_t polygonCount() const { return surf_.size(); }
    size_t pointCount() const { return surf_.points().size(); }
    size_t vertexCount(size_t n) const { return 3; }

    //- Returns position pos in local grid index space
    //  for polygon n and vertex v
    void getIndexSpacePoint(size_t n, size_t v, openvdb::Vec3d& pos) const
    {
        const Foam::label& idx = surf_[Foam::label(n)][Foam::label(v)];
        const Foam::point& p = surf_.points()[idx];
        pos[0] = double(p[0]);
        pos[1] = double(p[1]);
        pos[2] = double(p[2]);
    }
}; //struct MeshDataAdapter


std::vector<FloatGrid::Ptr> Foam::foamVDB::generateLODSequence
(
    const FloatGrid::Ptr grid,
    const label minLevel,
    const label maxLevel
)
{
    Timer timer("generateLODsequence");

    std::vector<FloatGrid::Ptr> grids(maxCellLevel_ + 1, nullptr);

    // Generate LOD sequence from minLevel to level 0 by restriction (fine->coarse)
    MultiResGrid mrgCoarse(minLevel + 1, *grid);

    // Generate LOD sequence from minLevel to maxLevel by prolongation (coarse->fine)
    const scalar maxLevelVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - maxLevel));
    MultiResGrid mrgFine
    (
        max((maxLevel - minLevel + 1), 2),
        grid->deepCopy(),
        maxLevelVoxelSize
    );

    for (label cellLevel = 0; cellLevel <= maxLevel; ++cellLevel)
    {
        if (cellLevel == minLevel)
        {
            grids[cellLevel] = grid;
        }
        else if (cellLevel < minLevel)
        {
            // in OpenVDB level 0 is the finest
            size_t vdbLevel = mrgCoarse.coarsestLevel() - cellLevel;
            grids[cellLevel] = mrgCoarse.grid(vdbLevel);
        }
        else
        {
            size_t vdbLevel = mrgFine.coarsestLevel() - cellLevel + minLevel;
            grids[cellLevel] = mrgFine.grid(vdbLevel);
        }

        Info<< "cellLevel " << cellLevel << " "
            << grids[cellLevel]->getName()
            << (cellLevel == minLevel ? "           " : "")
            << " voxels: " << grids[cellLevel]->activeVoxelCount()
            << endl;
    }

    return grids;
} //generateLODSequence


void Foam::foamVDB::multiResUDF
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
    bool curvature /*=false*/
)
{
    Timer timer("multiResUDF " + name);

    word prefix = "Boundary " + name;

    label surfLevel = surfaces.minLevel()[regioni];

    if (curvature)
    {
        prefix = name;
        surfLevel = surfaces.maxLevel()[regioni];
    }

    // coarse Unsigned Distance Field
    FloatGrid::Ptr udf =
        meshToLevelSet
        (
            inputSurf,
            surfLevel - 1,
            /*halfWidth*/nCellsBetweenLevels,
            /*conversioFlags*/openvdb::tools::UNSIGNED_DISTANCE_FIELD
        );
    word gridName = prefix + " level " + std::to_string(surfLevel - 1);
    udf->setName(gridName);

    std::vector<FloatGrid::Ptr> coarseUDF(maxCellLevel_ + 1, nullptr);
    coarseUDF[surfLevel - 1] = udf;
    bufferCellLevels[regioni] = coarseUDF;

    const scalar patchVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - surfLevel));
    // Prolongate (coarse->fine) bufferGrid to patchGrid (1 level finer)
    // Generate LOD sequence from coarsest level
    MultiResGrid mrg( /*levels*/2, udf->deepCopy(), patchVoxelSize);
    FloatGrid::Ptr highResGrid = mrg.grid(0);

    gridName = prefix + " level " + std::to_string(surfLevel);
    highResGrid->setName(gridName);

    std::vector<FloatGrid::Ptr> fineUDF(maxCellLevel_ + 1, nullptr);
    fineUDF[surfLevel] = highResGrid->deepCopy();
    patchCellLevels[regioni] = fineUDF;
} // multiResUDF


BoolGrid::Ptr Foam::foamVDB::partialInterior
(
    FloatGrid::Ptr fineUDF,
    const label surfLevel,
    const label nCellsBetweenLevels,
    std::vector<FloatGrid::Ptr> SDFCellLevels,
    BoolGrid::Ptr mask
)
{
    openvdb::tools::dilateActiveValues(mask->tree(), nCellsBetweenLevels);
    FloatGrid::Ptr topologyLs =
        openvdb::tools::topologyToLevelSet
        (
            *mask,
            /*halfWidth*/3,
            /*closingSteps*/0,
            /*dilation*/0,
            /*smoothingSteps*/0
        );

    //intersection between mask topologyLs and prolongated SDF
    topologyLs =
        openvdb::tools::csgIntersectionCopy<FloatGrid>
        (
            *SDFCellLevels[surfLevel],
            *topologyLs
        );

    //remove voxels with coarse distance
    topologyLs->topologyDifference(*fineUDF);
    //replace them with high-resolution dist
    topologyLs->tree().merge(fineUDF->tree()); //empties fineUDF tree
    //topologyLs->tree().voxelizeActiveTiles();

    //extract enclosed region
    BoolGrid::Ptr enclGrid = openvdb::tools::extractEnclosedRegion
    (
        *topologyLs,
        /*isovalue*/topologyLs->transform().voxelSize().x()
    );
    openvdb::tools::erodeActiveValues(enclGrid->tree(), 1);

    return enclGrid;
} //partialInterior


struct maskInsideVoxels
{
    static inline void op
    (
        const FloatGrid::ValueOnCIter& iter,
        BoolGrid::Accessor& accessor
    )
    {
        if (iter.getValue() < 0.0)
        {
            accessor.setValue(iter.getCoord(), true);
        }
        else
        {
            accessor.setValue(iter.getCoord(), false);
        }
    }
};

struct setInsideVoxelsOn
{
    static inline void op
    (
        const FloatGrid::ValueOnCIter& iter,
        BoolGrid::Accessor& accessor
    )
    {
        if (iter.getValue() < 0.0)
        {
            accessor.setValueOn(iter.getCoord(), true);
        }
    }
};

// see tools/Composite.h compReplace
struct SurfCurvLeaf
{
    FloatGrid::ConstAccessor distAcc_;
    Foam::scalar voxelSize_, curvature_;

    SurfCurvLeaf
    (
        const FloatGrid& distGrid,
        const Foam::scalar curvature
    )
    :
        voxelSize_(distGrid.transform().voxelSize().x()),
        curvature_(curvature),
        distAcc_(distGrid.getConstAccessor())
    {}

    inline void operator()
    (
        const openvdb::FloatTree::LeafIter& leafIter
    ) const
    {
        //FloatGrid::ConstAccessor distAcc(distGrid_.getConstAccessor);

        for (auto iter = leafIter->beginValueOn(); iter; ++iter)
        {
            if (abs(distAcc_.getValue(iter.getCoord())) > voxelSize_)
            {
                iter.setValueOff();
            }
            else if (abs(iter.getValue()) < curvature_)
            {
                iter.setValueOff();
            }
        }
    }
}; // SurfCurvLeaf

struct SurfCurv
{
    FloatGrid::ConstAccessor distAcc_;
    Foam::scalar voxelSize_, curvature_;

    SurfCurv
    (
        const FloatGrid& distGrid,
        const Foam::scalar curvature
    )
    :
        voxelSize_(distGrid.transform().voxelSize().x()),
        curvature_(curvature),
        distAcc_(distGrid.getConstAccessor())
    {}

    inline void operator()
    (
        const FloatGrid::ValueOnIter& iter
    ) const
    {
        if (abs(distAcc_.getValue(iter.getCoord())) > voxelSize_)
        {
            iter.setValueOff();
        }
        else if (abs(iter.getValue()) < curvature_)
        {
            iter.setValueOff();
        }
    }
}; // SurfCurv


FloatGrid::Ptr Foam::foamVDB::createCurvatureMask
(
    const FloatGrid& inGrid,
    const scalar curvature
)
{
    FloatGrid::Ptr curv = openvdb::tools::meanCurvature(inGrid);

    const scalar voxelSize = inGrid.transform().voxelSize().x();

    //SurfCurv op(inGrid, curvature);
    SurfCurvLeaf op(inGrid, curvature);

    openvdb::tools::foreach
    (
        //curv->beginValueOn(),
        curv->tree().beginLeaf(),
        op,
        /*threaded*/true,
        /*shareOp*/false
    );

    std::vector<FloatGrid::Ptr> segments;
    openvdb::tools::segmentActiveVoxels(*curv, segments);

    label segmentVoxels = 0;
    for (size_t i = 0; i < segments.size(); ++i)
    {
        segmentVoxels = segments[i]->activeVoxelCount();

        if (segmentVoxels < 2)
        {
            curv->topologyDifference(*segments[i]);
        }
        //else
        //{
        //Info<<"Segment " << i
        //    <<" active voxels " << segments[i]->activeVoxelCount()
        //    <<endl;
        //}
    }
    openvdb::tools::dilateActiveValues(curv->tree(), 2);

    Info<< "Active voxels in curvature mask (curv > " << curvature << "): "
        << curv->activeVoxelCount()
        << endl;

    return curv;
}

//void Foam::foamVDB::createPartialNarrowBandLevelSetGrids
//(
//    const triSurface& inputSurf,
//    const label minLevel,
//    const label maxLevel,
//    const word regionName,
//    std::vector<FloatGrid::Ptr>& SDFCellLevels,
//    std::vector<FloatGrid::Ptr>& boundaryN,
//    std::vector<BoolGrid::Ptr>& interiorGrids,
//    std::vector<IndexGrid::Ptr>& polygonIdCellLevels
//)
//{
//    BoolGrid::Ptr interior;
//
//    for (label cellLevel = minLevel; cellLevel <= maxLevel; ++cellLevel)
//    {
//        word gridName = "Partial narrow band SDF level "
//                        + std::to_string(cellLevel)
//                        + " region " + regionName;
//
//        Timer timer(gridName);
//
//        // create maskGrid
//        BoolGrid::Ptr boundaryMask = BoolGrid::create(false);
//        boundaryMask->setTransform(boundaryN[cellLevel]->transformPtr());
//
//        FloatGrid::Ptr coarseMask = SDFCellLevels[cellLevel-1]->deepCopy();
//        FloatGrid::Ptr boundaryGrid = boundaryN[cellLevel-1]->deepCopy();
//
//        openvdb::tools::dilateActiveValues
//        (
//            boundaryGrid->tree(),
//            3
//        );
//
//        coarseMask->topologyIntersection(*boundaryGrid);
//
//        FloatGrid::Ptr fineSDF = FloatGrid::create();
//        fineSDF->setTransform(boundaryN[cellLevel]->transformPtr());
//
//        coarseMask->setGridClass(openvdb::GRID_UNKNOWN);
//        //openvdb::tools::resampleToMatch<openvdb::tools::QuadraticSampler>(*coarseMask, *fineSDF);
//        //openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(*coarseMask, *fineSDF);
//        openvdb::tools::resampleToMatch<openvdb::tools::PointSampler>(*coarseMask, *fineSDF);
//
//        openvdb::tools::transformValues
//        (
//            fineSDF->cbeginValueOn(),
//            *boundaryMask,
//            maskInsideVoxels::op
//        );
//        // END create maskGrid
//
//        // offset surface to be cell-centered
//        // as voxel is located at min corner of a cell
//        const scalar worldEdgeLength = voxelSize_*Foam::pow(2, maxCellLevel_ - cellLevel);
//
//        const point surfOffset(-0.5*worldEdgeLength);
//
//        triSurface offsetSurf(inputSurf);
//
//        offsetSurf.translatePoints(surfOffset);
//
//        // transform points to local grid index space
//        offsetSurf.scalePoints(1.0 / (worldEdgeLength));
//
//        const MeshDataAdapter vdbSurfMesh(offsetSurf);
//
//        //IndexGrid::Ptr polygonIdGrid = IndexGrid::create();
//        using IntGridT = typename FloatGrid::template ValueConverter<openvdb::Int32>::Type;
//        typename IntGridT::Ptr polygonIdGrid;
//        polygonIdGrid.reset(new IntGridT(0));
//
//        FloatGrid::Ptr partialSDFGrid =
//            openvdb::tools::meshToVolume<FloatGrid>
//            (
//                vdbSurfMesh,
//                boundaryMask->transform(),
//                /*exteriorBandWidth =*/10.0f/*3.0f*/,
//                /*interiorBandWidth =*/10.0f/*3.0f*/,
//                openvdb::tools::PARTIAL_SIGNED_DISTANCE_FIELD,
//                polygonIdGrid.get(),
//                boundaryMask.get()
//            );
//
//        partialSDFGrid->setGridClass(openvdb::GRID_UNKNOWN);
//        partialSDFGrid->setName(gridName);
//
//        BoolGrid::Ptr interior = BoolGrid::create(false);
//
//        openvdb::tools::transformValues
//        (
//            partialSDFGrid->cbeginValueOn(),
//            *interior,
//            setInsideVoxelsOn::op
//        );
//
//        //interior->tree().voxelizeActiveTiles(); //save space for compressed streaming?
//
//        if (interiorGrids[cellLevel]->activeVoxelCount() > 0)
//        {
//            interiorGrids[cellLevel]->topologyUnion(*interior);
//        }
//        else
//        {
//            interior->setName("interiorGrids level " + std::to_string(cellLevel));
//
//            interiorGrids[cellLevel] = interior;
//        }
//
//        if (isValid<FloatGrid>(SDFCellLevels[cellLevel]))
//        {
//            openvdb::tools::csgUnion
//            (
//                *SDFCellLevels[cellLevel],
//                *partialSDFGrid // the operation always leaves this empty
//                ///*GGG debug*/*fineSDF
//            );
//        }
//        else
//        {
//            SDFCellLevels[cellLevel] = partialSDFGrid;//fineSDF;//GGG debug
//        }
//
//        //TODO cannot do csgUnion of polygonIdGrids, I have to merge them
//        polygonIdCellLevels[cellLevel] = polygonIdGrid;
//    } // for cellLevel
//} // createPartialNarrowBandLevelSetGrids


std::vector<FloatGrid::Ptr> Foam::foamVDB::multiResSDFandDispl
(
    const triSurface& inputSurf,
    const label minLevel,
    const label maxLevel,
    std::vector<BoolGrid::Ptr>& interiorGrids,
    //std::vector<IndexGrid::Ptr>& polygonIdCellLevels
    std::vector<Vec3dGrid::Ptr>& displacementGrids
)
{
    Timer timer("multiResSDFandDispl (maxLevel " + std::to_string(maxLevel) + ")");

    std::vector<FloatGrid::Ptr> SDFGrids(maxCellLevel_ + 1, nullptr);

    for (label cellLevel = minLevel-1; cellLevel <= maxLevel; ++cellLevel)
    {
        //Vec3dGrid::Ptr displacementGrid = Vec3dGrid::create();
        using DisplGrid = typename FloatGrid::template ValueConverter<openvdb::Vec3d>::Type;
        typename DisplGrid::Ptr displacementGrid;
        displacementGrid.reset(new DisplGrid(openvdb::Vec3d(0)));

        FloatGrid::Ptr sdf =
            meshToLevelSet
            (
                inputSurf,
                cellLevel,
                /*halfWidth*/3,
                /*conversionFlags*/0, //i.e. Signed Distance Field
                /*polygonIndexGrid*/nullptr,
                /*maskGrid*/nullptr,
                /*polygonListGrid*/nullptr,
                displacementGrid.get()
            );

        const word gridName = "Narrow band SDF level " + std::to_string(cellLevel);
        sdf->setName(gridName);

        SDFGrids[cellLevel] = sdf;
        displacementGrids[cellLevel] = displacementGrid;

        interiorGrids[cellLevel] = openvdb::tools::sdfInteriorMask(*sdf);

        Info<< "cellLevel " << cellLevel << " "
            << sdf->getName()
            << " voxels: " << sdf->activeVoxelCount()
            << endl;
    }
    return SDFGrids;
} //multiResSDFandDispl


std::vector<FloatGrid::Ptr> Foam::foamVDB::multiResSDF
(
    const triSurface& inputSurf,
    const label minLevel,
    const label maxLevel,
    std::vector<BoolGrid::Ptr>& interiorGrids
    //std::vector<IndexGrid::Ptr>& polygonIdCellLevels
    //std::vector<Vec3dGrid::Ptr>& displacementGrids
)
{
    Timer timer("multiResSDF (minLevel " + std::to_string(minLevel) + ")");

    //////IndexGrid::Ptr polygonIdGrid = IndexGrid::create();
    ////using IntGrid = typename FloatGrid::template ValueConverter<openvdb::Int32>::Type;
    ////typename IntGrid::Ptr polygonIdGrid;
    ////polygonIdGrid.reset(new IntGrid(0));

    ////Vec3dGrid::Ptr displacementGrid = Vec3dGrid::create();
    //using DisplGrid = typename FloatGrid::template ValueConverter<openvdb::Vec3d>::Type;
    //typename DisplGrid::Ptr displacementGrid;
    //displacementGrid.reset(new DisplGrid(openvdb::Vec3d(0)));

    FloatGrid::Ptr sdf =
        meshToLevelSet
        (
            inputSurf,
            minLevel,
            /*halfWidth*/3,
            /*conversionFlags*/0 //i.e. Signed Distance Field
            ///*polygonIndexGrid*/nullptr,
            ///*maskGrid*/nullptr,
            ///*polygonListGrid*/nullptr,
            //displacementGrid.get()
        );

    const word gridName = "Narrow band SDF level " + std::to_string(minLevel);
    sdf->setName(gridName);

    ////polygonIdCellLevels[minLevel] = polygonIdGrid;
    //displacementGrids[minLevel] = displacementGrid;

    std::vector<FloatGrid::Ptr> SDFCellLevels =
        generateLODSequence
        (
            sdf,
            minLevel,
            maxLevel
        );

    // fill interiorGrids of coarse grids
    for (label cellLevel = 0; cellLevel <= maxLevel; ++cellLevel)
    {
        interiorGrids[cellLevel] = openvdb::tools::sdfInteriorMask(*SDFCellLevels[cellLevel]);

        //const label edgeLength = Foam::pow(2, maxCellLevel_ - cellLevel);

        //scalar isoValue = 0.5 * edgeLength; //GGG

        //BoolGrid::Ptr interior = openvdb::tools::sdfInteriorMask(*SDFCellLevels[cellLevel], isoValue);

        //if (!isValid<BoolGrid>(interior)) continue;

        ////interior->tree().voxelizeActiveTiles(); //save space for compressed streaming?
        //interior->pruneGrid();

        //if (isValid<BoolGrid>(interiorGrids[cellLevel]))
        //{
        //    interiorGrids[cellLevel]->topologyUnion(*interior);
        //}
        //else
        //{
        //    interior->setName("interiorGrids level " + std::to_string(cellLevel));

        //    interiorGrids[cellLevel] = interior;
        //}
        //
    }

    return SDFCellLevels;
} // multiResSDF


const Foam::triSurface Foam::foamVDB::volumeToMesh
(
    const scalar iso,
    FloatGrid::Ptr refGrid,
    const word surfName
)
{
    Timer timer("volumeToMesh " + surfName);

    openvdb::tools::VolumeToMesh mesher(iso, /*adaptivity*/0);

    mesher(*refGrid);

    pointField points(mesher.pointListSize());

    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, points.size()),
        PointListCopy
        (
            mesher.pointList(),
            points
        )
    );

    openvdb::tools::PolygonPoolList& polygonPoolList = mesher.polygonPoolList();
    //
    //// Count triangles and quads
    //size_t numQuads = 0, numTriangles = 0;
    //for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n)
    //{
    //    openvdb::tools::PolygonPool& polygons = polygonPoolList[n];
    //    numTriangles += polygons.numTriangles();
    //    numQuads += polygons.numQuads();
    //}

    // Set triangles
    List<labelledTri> faces;
    //faces.resize(numTriangles + 2*numQuads);
    List<labelledTri> tmpFaces;

    for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n)
    {
        openvdb::tools::PolygonPool& polygons = polygonPoolList[n];

        tmpFaces.resize(polygons.numTriangles());

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, polygons.numTriangles()),
            CopyTriangles
            (
                polygons,
                tmpFaces
            )
        );

        faces.append(tmpFaces);

        tmpFaces.resize(2*polygons.numQuads());

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, polygons.numQuads()),
            CopyQuads
            (
                polygons,
                tmpFaces
            )
        );

        faces.append(tmpFaces);
    }

    triSurface surf
    (
        faces,
        points
    );

    return surf;
} //volumeToMesh


std::vector<FloatGrid::Ptr> Foam::foamVDB::createNarrowBandLevelSetGridsPar
(
    const triSurface& inputSurf
)
{
    std::vector<FloatGrid::Ptr> SDFCellLevels(maxCellLevel_ + 1, nullptr);
    std::vector<IndexGrid::Ptr> polygonIdCellLevels(maxCellLevel_ + 1, nullptr);

    openvdb::GridPtrVecPtr SDFProcGrids(new openvdb::GridPtrVec);

    for (label cellLevel = maxCellLevel_; cellLevel >= 0; --cellLevel)
    {
        // Example: proc - cellLevel
        //          0       9
        //          1       8
        //          2       7
        //          3       6 5 4 3 2 1
        if
        (
            Pstream::myProcNo() < Pstream::nProcs() - 1
         && Pstream::myProcNo() + cellLevel != maxCellLevel_
        )
        {
            continue;
        }

        if (Pstream::myProcNo() + cellLevel > maxCellLevel_)
        {
            continue;
        }

        word gridName = "Narrow band SDF level " + std::to_string(cellLevel);

        Timer timer(gridName);

        //Pout<< "Processing grid: " + gridName<<endl;

        // offset surface to be cell-centered
        // as voxel is located at min corner of a cell
        const label edgeLength = Foam::pow(2, maxCellLevel_ - cellLevel);

        const point surfOffset(-(voxelSize_ * edgeLength) / 2.);

        triSurface offsetSurf(inputSurf);

        offsetSurf.translatePoints(surfOffset);

        // transform points to local grid index space
        offsetSurf.scalePoints(1.0 / (voxelSize_ * edgeLength));

        const MeshDataAdapter vdbSurfMesh(offsetSurf);

        //IndexGrid::Ptr polygonIdGrid = IndexGrid::create();
        using IntGridT = typename FloatGrid::template ValueConverter<openvdb::Int32>::Type;
        typename IntGridT::Ptr polygonIdGrid;
        polygonIdGrid.reset(new IntGridT(0));

        openvdb::math::Transform::Ptr linearTransform =
            openvdb::math::Transform::createLinearTransform
            (
                voxelSize_ * edgeLength
            );

        FloatGrid::Ptr SDFGrid =
            openvdb::tools::meshToVolume<FloatGrid>
            (
                vdbSurfMesh,
                *linearTransform,
                /*exteriorBandWidth = */3.0f,
                /*interiorBandWidth = */3.0f,
                /*int flags = */0, //computes Signed Distance Field!
                polygonIdGrid.get() //typename GridType::template ValueConverter<Int32>::Type* polygonIndexGrid = nullptr
            );

        SDFGrid->setGridClass(openvdb::GRID_LEVEL_SET);
        SDFGrid->setName(gridName);

        SDFCellLevels[cellLevel] = SDFGrid;
        polygonIdCellLevels[cellLevel] = polygonIdGrid;
    }

    gatherProcGrids<FloatGrid>(SDFCellLevels);
    //gatherProcGrids<IndexGrid>(polygonIdCellLevels);

    return SDFCellLevels;
} // createNarrowBandLevelSetGridsPar


template<typename GridType>
std::vector<openvdb::GridPtrVecPtr> Foam::foamVDB::sendReceiveGrids
(
    std::vector<typename GridType::Ptr>& procGrids,
    const word& gridsName,
    bool GridPtrVecPtr
)
{
    Timer timer("sendReceiveGrids " + gridsName);

    PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

    // Stream the grids to a string.
    // The grids get compressed! see openvdb::io::Archive::DEFAULT_COMPRESSION_FLAGS
    std::ostringstream ostr(std::ios_base::binary);

    openvdb::io::Stream(ostr).write(procGrids);

    std::streamsize strSize = ostr.str().size();

    //if (debug)
    //{
    //    //find compression ratio
    //    std::ostringstream uncompressedOstr(std::ios_base::binary);

    //    openvdb::io::Stream gridStream(uncompressedOstr);

    //    gridStream.setCompression(openvdb::io::COMPRESS_NONE);

    //    gridStream.write(procGrids);

    //    std::streamsize uncompressedSize = uncompressedOstr.str().size();

    //    Pout<< gridsName << " initial size: "
    //        << uncompressedSize / 1e6 << " Mb; "
    //        << "compressed size : "
    //        << strSize / 1e6 << " Mb; "
    //        << "compression ratio: "
    //        << float(uncompressedSize) / float(strSize)
    //        << endl;
    //}

    const label nProcs = Pstream::nProcs();

    // myProcNo sends to all other procs
    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == Pstream::myProcNo()) continue;

        //non-blocking communication
        UOPstream toProcSize(proci, pBufSize);
        toProcSize << strSize;

        UOPstream toProc(proci, pBufGrid);
        toProc.write(ostr.str().data(), strSize);
    }

    //Start sending and receiving and block
    pBufSize.finishedSends();
    pBufGrid.finishedSends();

    std::vector<openvdb::GridPtrVecPtr> procGridList(nProcs, nullptr);

    // myProcNo receives from all other procs
    for (label proci = 0; proci < nProcs; ++proci)
    {
        if (proci == Pstream::myProcNo()) continue;

        // non-blocking communication
        UIPstream fromProcSize(proci, pBufSize);
        std::streamsize strSize;
        fromProcSize >> strSize;

        char * incoming = new char[strSize];

        UIPstream fromProc(proci, pBufGrid);
        fromProc.read(incoming, strSize);

        std::stringstream ss;
        ss.write(incoming, strSize);

        std::istringstream istr;
        istr.str(ss.str());

        openvdb::io::Stream strm(istr);
        openvdb::GridPtrVecPtr grids = strm.getGrids();
        procGridList[proci] = grids;
    }

    return procGridList;
} //sendReceiveGrids GridPtrVecPtr


template<typename GridType>
std::vector<std::vector<typename GridType::Ptr>> Foam::foamVDB::sendReceiveGrids
(
    std::vector<typename GridType::Ptr>& procGrids,
    const word& gridsName
)
{
    std::vector<openvdb::GridPtrVecPtr> procGridList =
        sendReceiveGrids<GridType>
        (
            procGrids,
            gridsName,
            /*GridPtrVecPtr*/true
        );

    using GridTypePtr = typename GridType::Ptr;

    std::vector<std::vector<GridTypePtr>> procGridCollection(Pstream::nProcs());

    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        if (proci == Pstream::myProcNo())
        {
            procGridCollection[proci] = procGrids;
        }
        else
        {
            openvdb::GridPtrVecPtr grids = procGridList[proci];

            std::vector<GridTypePtr> procICellLevelGrids(maxCellLevel_ + 1);

            for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
            {
                word gridName = gridsName + " level " + std::to_string(cellLevel)
                                + " proc" + std::to_string(proci);

                openvdb::GridBase::Ptr baseGrid = openvdb::findGridByName(*grids, gridName);

                GridTypePtr grid;

                if (baseGrid)
                {
                    grid = openvdb::gridPtrCast<GridType>(baseGrid);
                }
                else
                {
                    grid = GridType::create();
                }

                grid->setName(gridName);

                procICellLevelGrids[cellLevel] = grid;
            }

            procGridCollection[proci] = procICellLevelGrids;
        }
    }

    return procGridCollection;
} //sendReceiveGrids


void Foam::foamVDB::intersectProcGrids
(
    std::vector<FloatGrid::Ptr>& procGrids,
    std::vector<BoolGrid::Ptr>& interiorGrids
)
{
    std::vector<std::vector<FloatGrid::Ptr>> procGridCollection =
        sendReceiveGrids<FloatGrid>
        (
            procGrids,
            "cellLevelGrids"
        );

    std::vector<std::vector<BoolGrid::Ptr>> interiorGridCollection =
        sendReceiveGrids<BoolGrid>
        (
            interiorGrids,
            "interiorGrids"
        );

    // topology union of interior voxels of each processor
    combineGrids<BoolGrid>
    (
        "interiorGrids procs",
        interiorGridCollection,
        interiorGrids
    );

    // topology union of processor grids and remove interior voxels
    combineGrids
    (
        "procGrids",
        procGridCollection,
        procGrids,
        interiorGrids
    );
}// void intersectProcGrids interiorGrids


template<typename GridType>
void Foam::foamVDB::gatherProcGrids
(
    std::vector<typename GridType::Ptr>& procGrids
)
{
    std::vector<openvdb::GridPtrVecPtr> procGridList =
        sendReceiveGrids<FloatGrid>
        (
            procGrids,
            "SDF",
            /*GridPtrVecPtr*/true
        );

    for (label proci = 0; proci < Pstream::nProcs(); ++proci)
    {
        if (proci == Pstream::myProcNo()) continue;

        openvdb::GridPtrVecPtr grids = procGridList[proci];

        for (label cellLevel = maxCellLevel_; cellLevel >= 0; --cellLevel)
        {
            // Example:
            // proc0 = [_ _ _ _ _ _ _ _ _ 9]
            // proc1 = [_ _ _ _ _ _ _ _ 8 _]
            // proc2 = [_ _ _ _ _ _ _ 7 _ _]
            // proc3 = [0 1 2 3 4 5 6 _ _ _]

            if
            (
                proci < Pstream::nProcs() - 1
             && proci + cellLevel != maxCellLevel_
            )
            {
                continue;
            }

            if (proci + cellLevel > maxCellLevel_)
            {
                continue;
            }

            word gridName = "Narrow band SDF level " + std::to_string(cellLevel);

            //Pout<< "receiving " << gridName << " from proc " << proci <<endl;

            openvdb::GridBase::Ptr baseGrid = openvdb::findGridByName(*grids, gridName);
            typename GridType::Ptr grid = openvdb::gridPtrCast<GridType>(baseGrid);

            procGrids[cellLevel] = grid;
        } // for cellLevel
    } // for proci
} //gatherProcGrids


void Foam::foamVDB::removeSolidVoxels
(
    FloatGrid::Ptr SDFGrid,
    FloatGrid::Ptr grid
)
{
    struct removeSolidVoxelsOp
    {
        FloatGrid::ConstAccessor SDFAcc_;

        removeSolidVoxelsOp
        (
            FloatGrid& SDFGrid
        )
        :
            SDFAcc_(SDFGrid.getConstAccessor())
        {}

        inline void operator()(const FloatGrid::ValueOnIter& iter) const
        {
            if (SDFAcc_.getValue(iter.getCoord()) < 0.0f)
            {
                iter.setValueOff();
            }
        }
    };

    openvdb::tools::foreach
    (
        grid->beginValueOn(),
        removeSolidVoxelsOp
        (
            *SDFGrid
        ),
        /*threaded*/true,
        /*shareOp*/false
    );
} // removeSolidVoxels


//std::vector<FloatGrid::Ptr> Foam::foamVDB::createNarrowBandMultiResGrids
//(
//    const triSurface& inputSurf
//)
//{
//    std::vector<FloatGrid::Ptr> surfCellLevels(maxCellLevel_ + 1, nullptr);
//
//    point surfOffset(-(voxelSize_ / 2.));
//
//    triSurface offsetSurf(inputSurf);
//
//    offsetSurf.translatePoints(surfOffset);
//
//    {
//        Timer timer;
//        FloatGrid::Ptr levelNSurfGrid =
//            meshToLevelSet
//            (
//                offsetSurf,
//                maxCellLevel_,
//                /*halfWidth*/1
//            );
//
//        word gridName = "Boundary narrow band level " + std::to_string(maxCellLevel_) + " cells";
//
//        levelNSurfGrid->setName(gridName);
//        surfCellLevels[maxCellLevel_] = levelNSurfGrid;
//
//        //if (debug)
//        {
//            Info<< gridName << ": active voxels ";
//            if (Pstream::parRun())
//            {
//                List<label> procActiveVoxels(Pstream::nProcs());
//                procActiveVoxels[Pstream::myProcNo()] = levelNSurfGrid->activeVoxelCount();
//                IPstream::gatherList(procActiveVoxels);
//                Info<< procActiveVoxels << endl;
//            }
//            else
//            {
//                Info<< levelNSurfGrid->activeVoxelCount() <<endl;
//            }
//        }
//    }
//
//    const word name = "narrow band";
//
//    {
//        Timer timer("RestrictGrids");
//        restrictGrids
//        (
//            name,
//            surfCellLevels,
//            maxCellLevel_
//        );
//    }
//
//    return surfCellLevels;
//} //createNarrowBandMultiResGrids


void Foam::foamVDB::resample
(
    const BoolGrid& inGrid,
    BoolGrid& outGrid
)
{
    Timer timer("resample");

    const openvdb::math::Transform
        &sourceTransform = inGrid.transform(),
        &targetTransform = outGrid.transform();

    openvdb::Mat4R xform =
        sourceTransform.baseMap()->getAffineMap()->getMat4() *
        targetTransform.baseMap()->getAffineMap()->getMat4().inverse();

    openvdb::tools::GridTransformer transformer(xform);

    // Resample using nearest-neighbor interpolation.
    // Transform the input grid and store the result in the output grid.
    transformer.transformGrid<openvdb::tools::PointSampler, BoolGrid>
    (
        inGrid,
        outGrid
    );
} // resample


void Foam::foamVDB::combineGrids
(
    const word name,
    std::vector<std::vector<FloatGrid::Ptr>>& cellLevelGridsCollection,
    std::vector<FloatGrid::Ptr>& cellLevelGrids,
    const std::vector<BoolGrid::Ptr>& interiorGrids
)
{
    combineGrids<FloatGrid>
    (
        name,
        cellLevelGridsCollection,
        cellLevelGrids
    );

    // remove solid and intersecting voxels
    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        if (interiorGrids[cellLevel]->activeVoxelCount() > 0)
        {
            cellLevelGrids[cellLevel]->topologyDifference(*interiorGrids[cellLevel]);
            openvdb::tools::pruneInactive(cellLevelGrids[cellLevel]->tree());
        }
    } // for cellLevel
} // combineGrids interiorGrids

void Foam::foamVDB::combineGrids
(
    const word name,
    std::vector<std::vector<FloatGrid::Ptr>>& cellLevelGridsCollection,
    std::vector<FloatGrid::Ptr>& cellLevelGrids,
    const std::vector<FloatGrid::Ptr>& SDFCellLevels
)
{
    combineGrids<FloatGrid>
    (
        name,
        cellLevelGridsCollection,
        cellLevelGrids
    );

    Timer timer("remove solid voxels");

    // remove solid voxels
    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        //removeSolidVoxels
        //(
        //    SDFCellLevels[cellLevel],
        //    cellLevelGrids[cellLevel]
        //);

        auto interior = openvdb::tools::sdfInteriorMask(*SDFCellLevels[cellLevel]);

        interior->pruneGrid();

        cellLevelGrids[cellLevel]->topologyDifference(*interior);
        openvdb::tools::pruneInactive(cellLevelGrids[cellLevel]->tree());
    } // for cellLevel
} // combineGrids SDFCellLevels


std::vector<FloatGrid::Ptr> Foam::foamVDB::distanceRefinementGrading
(
    const word& name,
    const triSurface& inputSurf,
    const List<labelVector>& distanceLevels,
    const scalarField& distances,
    const label nCellsBetweenLevels,
    const openvdb::CoordBBox blockMeshBBox,
    std::vector<BoolGrid::Ptr>& refinementInterior
)
{
    Timer timer("distanceRefinementGrading " + name);
    Info<<nl;
    std::vector<FloatGrid::Ptr> surfCellLevels(maxCellLevel_ + 1, nullptr);

    label finestDistanceLevel = distanceLevels[0][0];

    label bufferLevel = finestDistanceLevel - 1;

    const scalar bufferVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - bufferLevel));

    const scalar halfWidth = distances[0]/bufferVoxelSize + 1;

    // coarse Unsigned Distance Field
    FloatGrid::Ptr bufferGrid =
        meshToLevelSet
        (
            inputSurf,
            bufferLevel,
            halfWidth,
            /*conversioFlags*/openvdb::tools::UNSIGNED_DISTANCE_FIELD
        );

    const scalar fineVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - distanceLevels[0][0]));

    // Prolongate (coarse->fine) bufferGrid to inside voxels (1 level finer)
    // Generate LOD sequence from coarsest level
    MultiResGrid mrg(/*levels*/2, bufferGrid->deepCopy(), fineVoxelSize);

    FloatGrid::Ptr insideGrid = mrg.grid(0);

    auto interior = openvdb::tools::sdfInteriorMask(*insideGrid, /*isoValue*/distances[0]);
    insideGrid->topologyUnion(*interior);
    insideGrid->topologyIntersection(*interior);

    // dilate to overlap (fill gaps) with coarser level
    openvdb::tools::dilateActiveValues(insideGrid->tree(), 2, openvdb::tools::NN_FACE_EDGE);

    const word gridName = name + " level " + Foam::name(distanceLevels[0][0]) + " cells";

    insideGrid->setName(gridName);

    label voxelCount = insideGrid->activeVoxelCount();
    //if (debug)
    {
        Info<< gridName << ": active voxels ";
        if (Pstream::parRun())
        {
            List<label> procActiveVoxels(Pstream::nProcs());
            procActiveVoxels[Pstream::myProcNo()] = insideGrid->activeVoxelCount();
            IPstream::gatherList(procActiveVoxels);
            Info<< procActiveVoxels << endl;
        }
        else
        {
            Info<< insideGrid->activeVoxelCount() <<endl;
        }
    }

    surfCellLevels[distanceLevels[0][0]] = insideGrid;

    scalar extWidth = nCellsBetweenLevels + 1;

    if (distances.size() > 1)
    {
        extWidth = (distances[1] - distances[0]) / bufferVoxelSize;
    }

    bufferGrid =
        openvdb::tools::levelSetRebuild<FloatGrid>
        (
            *bufferGrid,
            /*isoValue*/distances[0],
            extWidth,
            /*int*/0.1
        );

    voxelCount = 0;

    scalarField isoValues(bufferLevel+1);
    scalarField widths(bufferLevel+1);

    label i = 0;
    for (int cellLevel = bufferLevel; cellLevel >= 0; --cellLevel)
    {
        const scalar currentVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - cellLevel));

        if (i < distances.size()-1)
        {
            ++i;
            isoValues[cellLevel] = distances[i] - distances[i-1];
        }
        else
        {
            // + 0.5 so if for example nCellsBetweenLevels is 3, in particular areas will prefer 4 to 2
            // also, 0.5 more square refinement grading
            //       0.25 more rounded refinement grading
            isoValues[cellLevel] = currentVoxelSize * (nCellsBetweenLevels + 0.5);
        }

        if (i+1 < distances.size())
        {
            widths[cellLevel] = 1 + (distances[i+1] - distances[i]) / currentVoxelSize;
        }
        else
        {
            widths[cellLevel] = nCellsBetweenLevels + 1;
        }
    }

    calculateBufferGrids
    (
        name,
        bufferLevel,
        isoValues,
        widths,
        bufferGrid,
        voxelCount,
        surfCellLevels,
        refinementInterior
    );

    if (voxelCount > 0)
    {
        clipGrids
        (
            surfCellLevels,
            finestDistanceLevel,
            blockMeshBBox
        );

        trimGrids
        (
            surfCellLevels,
            finestDistanceLevel,
            /*nUpLevels*/1
        );
    }

    return surfCellLevels;
} // distanceRefinementGrading


std::vector<FloatGrid::Ptr> Foam::foamVDB::refinementRegionGrading
(
    const word& name,
    const triSurface& inputSurf,
    const label insideLevel,
    const label nCellsBetweenLevels,
    const openvdb::CoordBBox blockMeshBBox,
    std::vector<BoolGrid::Ptr>& refinementInterior
)
{
    Timer timer("refinementRegionGrading " + name);
    Info<<nl;
    std::vector<FloatGrid::Ptr> surfCellLevels(maxCellLevel_ + 1, nullptr);

    label bufferLevel = insideLevel - 1;

    FloatGrid::Ptr bufferGrid =
        meshToLevelSet
        (
            inputSurf,
            bufferLevel,
            /*halfWidth*/nCellsBetweenLevels + 1,
            openvdb::tools::DISABLE_RENORMALIZATION
        );

    const scalar insideVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - insideLevel));

    // Prolongate (coarse->fine) bufferGrid to inside voxels (1 level finer)
    // Generate LOD sequence from coarsest level
    MultiResGrid mrg( /*levels*/2, bufferGrid->deepCopy(), insideVoxelSize);
    FloatGrid::Ptr insideGrid = mrg.grid(0);

    auto interior = openvdb::tools::sdfInteriorMask(*insideGrid);
    insideGrid->topologyUnion(*interior);
    insideGrid->topologyIntersection(*interior);

    // dilate to overlap (fill gaps) with coarser level
    openvdb::tools::dilateActiveValues(insideGrid->tree(), 2, openvdb::tools::NN_FACE_EDGE);

    const word gridName = name + " level " + Foam::name(insideLevel) + " cells";

    insideGrid->setName(gridName);

    label voxelCount = insideGrid->activeVoxelCount();
    //if (debug)
    {
        Info<< gridName << ": active voxels ";
        if (Pstream::parRun())
        {
            List<label> procActiveVoxels(Pstream::nProcs());
            procActiveVoxels[Pstream::myProcNo()] = insideGrid->activeVoxelCount();
            IPstream::gatherList(procActiveVoxels);
            Info<< procActiveVoxels << endl;
        }
        else
        {
            Info<< insideGrid->activeVoxelCount() <<endl;
        }
    }

    surfCellLevels[insideLevel] = insideGrid;

    calculateBufferGrids
    (
        name,
        bufferLevel,
        nCellsBetweenLevels,
        bufferGrid,
        voxelCount,
        surfCellLevels,
        refinementInterior
    );

    if (voxelCount > 0)
    {
        clipGrids
        (
            surfCellLevels,
            insideLevel,
            blockMeshBBox
        );

        trimGrids
        (
            surfCellLevels,
            insideLevel,
            /*nUpLevels*/1
        );
    }

    return surfCellLevels;
} // refinementRegionGrading


void Foam::foamVDB::trimGrids
(
    const std::vector<FloatGrid::Ptr>& grids,
    const label maxSurfLevel,
    const label nUpLevels
)
{
    struct ExcludeFineVoxelsOp
    {
        label nUpLevel_;
        FloatGrid& fineGrid_;

        ExcludeFineVoxelsOp
        (
            const label& nUpLevel,
            FloatGrid& fineGrid
        )
        :
            nUpLevel_(nUpLevel),
            fineGrid_(fineGrid)
        {}

        //recursive function to touch all voxels at finer levels
        inline void setValuesOff
        (
            const openvdb::Coord& ijk,
            FloatGrid::Accessor& acc,
            label upLevel
        ) const
        {
            for (size_t i = 0; i < 8; ++i)
            {
                const openvdb::Coord ijkFine = (ijk << 1) + COARSE_TO_FINE[i];

                if (nUpLevel_ > upLevel)
                {
                    setValuesOff(ijkFine, acc, upLevel++);
                }
                else
                {
                    acc.setValueOff(ijkFine);
                }
            }
        }

        inline void operator()(const FloatGrid::ValueOnCIter& iter) const
        {
            FloatGrid::Accessor fineGridAcc(fineGrid_.getAccessor());

            const openvdb::Coord& ijk = iter.getCoord();

            setValuesOff(ijk, fineGridAcc, 1);
        }
    };

    for (int cellLevel = 0; cellLevel < maxSurfLevel; ++cellLevel)
    {
        if (!grids[cellLevel]) continue;

        FloatGrid::Ptr levelNSurfGrid = grids[cellLevel];

        for (label nUpLevel = 1; nUpLevel <= nUpLevels; nUpLevel++)
        {
            if (grids[cellLevel+nUpLevel] && maxSurfLevel - cellLevel >= nUpLevel)
            {
                //if (debug)
                //{
                //    std::cout <<"Excluding " << grids[cellLevel+nUpLevel]->getName()
                //          <<" which are inside level " << cellLevel << " cells"
                //          << std::endl;
                //}

                openvdb::tools::foreach
                (
                    levelNSurfGrid->cbeginValueOn(),
                    ExcludeFineVoxelsOp
                    (
                        nUpLevel,
                        *grids[cellLevel+nUpLevel]
                    ),
                    /*threaded*/true,
                    /*shareOp*/false //true
                );

                //if (debug)
                //{
                //    Info<<"Updated voxel count level "<<cellLevel+nUpLevel
                //        <<": "<<grids[cellLevel+nUpLevel]->activeVoxelCount()
                //        <<endl;
                //}
            } // if nUpLevel exists
        } // for nUpLevel
    } // for cellLevel
}


void Foam::foamVDB::topologyDifference
(
    const std::vector<FloatGrid::Ptr>& resultGrids,
    const std::vector<FloatGrid::Ptr>& grids
)
{
    for (size_t i = 0; i < grids.size(); ++i)
    {
        resultGrids[i]->topologyDifference(*grids[i]);
        openvdb::tools::pruneInactive(resultGrids[i]->tree());
    }
}

void Foam::foamVDB::topologyUnion
(
    const std::vector<FloatGrid::Ptr>& resultGrids,
    const std::vector<FloatGrid::Ptr>& grids
)
{
    for (size_t i = 0; i < grids.size(); ++i)
    {
        resultGrids[i]->topologyUnion(*grids[i]);
    }
}

std::vector<FloatGrid::Ptr> Foam::foamVDB::deepCopy
(
    const std::vector<FloatGrid::Ptr>& grids
)
{
    std::vector<FloatGrid::Ptr> gridsCopy;

    for (size_t i = 0; i < grids.size(); ++i)
    {
        FloatGrid::Ptr gridCopy = grids[i]->deepCopy();

        gridsCopy.push_back
        (
            gridCopy
        );
    }

    return gridsCopy;
}

openvdb::CoordBBox Foam::foamVDB::getLevelNBoundingBox
(
    const openvdb::math::BBox<openvdb::Vec3s>& bBox,
    const label& cellLevel,
    const label& decomposeDir
)
{
    scalar edgeLength = Foam::pow(2, cellLevel);

    scalar padding = edgeLength - 1;

    openvdb::CoordBBox levelNBBox
    (
        openvdb::Coord
        (
            bBox.min().x() * edgeLength,
            bBox.min().y() * edgeLength,
            bBox.min().z() * edgeLength
        ),
        openvdb::Coord
        (
            std::floor(bBox.max().x() * edgeLength + (decomposeDir == 0 ? 0 : padding)),
            std::floor(bBox.max().y() * edgeLength + (decomposeDir == 1 ? 0 : padding)),
            std::floor(bBox.max().z() * edgeLength + (decomposeDir == 2 ? 0 : padding))
        )
    );

    return levelNBBox;
}

openvdb::Index64 Foam::foamVDB::activeVoxelCountInBBox
(
    const std::vector<FloatGrid::Ptr>& grids,
    const label& maxSurfLevel,
    const openvdb::math::BBox<openvdb::Vec3s>& bBox,
    label decomposeDir /*= -1*/
)
{
    openvdb::Index64 activeVoxels = 0;

    for (label cellLevel = 0; cellLevel <= maxSurfLevel; ++cellLevel)
    {
        openvdb::CoordBBox levelNboundingBox =
            getLevelNBoundingBox
            (
                bBox,
                cellLevel,
                decomposeDir
            );

        activeVoxels += onVoxelsInCoordBbox(grids[cellLevel], levelNboundingBox);
    }

    return activeVoxels;
}

void Foam::foamVDB::clipGrids
(
    const std::vector<FloatGrid::Ptr>& grids,
    const label& maxSurfLevel,
    const openvdb::math::BBox<openvdb::Vec3s>& bBox,
    label decomposeDir /*= -1*/
)
{
    for (label cellLevel = 0; cellLevel <= maxSurfLevel; ++cellLevel)
    {
        openvdb::CoordBBox levelNboundingBox =
            getLevelNBoundingBox
            (
                bBox,
                cellLevel,
                decomposeDir
            );

        grids[cellLevel]->clip(levelNboundingBox);
        grids[cellLevel]->pruneGrid();
    }
}

void Foam::foamVDB::clipGrids
(
    const std::vector<FloatGrid::Ptr>& grids,
    const label& maxSurfLevel,
    const openvdb::CoordBBox& coordBBox
)
{
    openvdb::math::BBox<openvdb::Vec3s> bBox
    (
        openvdb::Vec3s
        (
            coordBBox.min().x(),
            coordBBox.min().y(),
            coordBBox.min().z()
        ),
        openvdb::Vec3s
        (
            coordBBox.max().x(),
            coordBBox.max().y(),
            coordBBox.max().z()
        )
    );

    clipGrids
    (
        grids,
        maxSurfLevel,
        bBox
    );
}

std::vector<FloatGrid::Ptr> Foam::foamVDB::surfaceGrading
(
    const word& name,
    const std::vector<std::vector<FloatGrid::Ptr>>& patchCellLevels,
    const std::vector<std::vector<FloatGrid::Ptr>>& bufferCellLevels,
    const label maxSurfLevel,
    const label nCellsBetweenLevels,
    const openvdb::CoordBBox blockMeshBBox,
    std::vector<BoolGrid::Ptr>& refinementInterior
)
{
    Info<<nl;
    std::vector<FloatGrid::Ptr> surfCellLevels(maxCellLevel_ + 1, nullptr);

    FloatGrid::Ptr patchUnion =
        topologyUnionLevelSets
        (
            patchCellLevels,
            maxSurfLevel
        );

    surfCellLevels[maxSurfLevel] = patchUnion;

    label bufferLevel = maxSurfLevel - 1;

    FloatGrid::Ptr coarsePatchUnion =
        csgUnionLevelSets
        (
            bufferCellLevels,
            bufferLevel
        );

    const scalar finerVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - maxSurfLevel));
    scalar isoValue = finerVoxelSize * (nCellsBetweenLevels + 0.5);

    FloatGrid::Ptr bufferGrid =
        openvdb::tools::levelSetRebuild<FloatGrid>
        (
            *coarsePatchUnion,
            isoValue,
            /*ext*/nCellsBetweenLevels + 1,
            /*int*/3*nCellsBetweenLevels //fill internal cavities
        );

    label voxelCount = 0;

    calculateBufferGrids
    (
        name,
        bufferLevel,
        nCellsBetweenLevels,
        bufferGrid,
        voxelCount,
        surfCellLevels,
        refinementInterior
    );

    clipGrids
    (
        surfCellLevels,
        maxSurfLevel,
        blockMeshBBox
    );

    trimGrids
    (
        surfCellLevels,
        maxSurfLevel,
        /*nUpLevels*/1
    );

    return surfCellLevels;
} // surfaceGrading boundaryCellLevel


void Foam::foamVDB::calculateBufferGrids
(
    const word& name,
    const label bufferLevel,
    const label nCellsBetweenLevels,
    FloatGrid::Ptr bufferGrid,
    label& voxelCount,
    std::vector<FloatGrid::Ptr>& gradingGrids,
    std::vector<BoolGrid::Ptr>& refinementInterior
)
{
    scalarField isoValues(bufferLevel+1);

    const scalarField widths(bufferLevel+1, nCellsBetweenLevels+1);

    for (int cellLevel = bufferLevel; cellLevel >= 0; --cellLevel)
    {
        const scalar currentVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - cellLevel));

        // + 0.5 so if for example nCellsBetweenLevels is 3, in particular areas will prefer 4 to 2
        // also, 0.5 more square refinement grading
        //       0.25 more rounded refinement grading
        isoValues[cellLevel] = currentVoxelSize * (nCellsBetweenLevels + 0.5);
    }

    calculateBufferGrids
    (
        name,
        bufferLevel,
        isoValues,
        widths,
        bufferGrid,
        voxelCount,
        gradingGrids,
        refinementInterior
    );

}


void Foam::foamVDB::calculateBufferGrids
(
    const word& name,
    const label bufferLevel,
    const scalarField& isoValues,
    const scalarField& widths,
    FloatGrid::Ptr bufferGrid,
    label& voxelCount,
    std::vector<FloatGrid::Ptr>& gradingGrids,
    std::vector<BoolGrid::Ptr>& refinementInterior
)
{
    scalar transformOffset = 0;

    for (int cellLevel = bufferLevel; cellLevel >= 0; --cellLevel)
    {
        //std::cout<< "bufferGrid Diagnostics:\n"
        //    << openvdb::tools::checkLevelSet(*bufferGrid, 9)
        //    << std::endl;

        FloatGrid::Ptr grid = bufferGrid->deepCopy();

        auto interior = openvdb::tools::sdfInteriorMask(*grid);

        //Do not remove coarse cells inside a cavity
        if (cellLevel == bufferLevel)
        {
            //extract enclosed region
            BoolGrid::Ptr cavity =
                openvdb::tools::extractEnclosedRegion
                (
                    *grid,
                    -isoValues[cellLevel]
                );
            interior->topologyDifference(*cavity);
            grid->topologyUnion(*cavity);
        }

        //interior->tree().voxelizeActiveTiles(); //save space for compressed streaming?
        interior->pruneGrid();

        // dilate to overlap (fill gaps) with coarser level
        openvdb::tools::dilateActiveValues(grid->tree(), 2);

        grid->topologyDifference(*interior);
        openvdb::tools::pruneInactive(grid->tree());

        refinementInterior[cellLevel] = interior;

        const word gridName = name + " level " + Foam::name(cellLevel) + " cells";

        grid->setName(gridName);

        grid->tree().prune();

        gradingGrids[cellLevel] = grid;

        //update zero level-set for next cellLevel
        //if (cellLevel > 0)
        //{
        //    Timer timer("fineToCoarse Mrg");
        //    const scalar currentVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - cellLevel));

        //    auto bufferGrid1 = openvdb::tools::dilateSdf
        //    (
        //        *bufferGrid,
        //        2*(nCellsBetweenLevels+1),
        //        openvdb::tools::NN_FACE_EDGE
        //    );

        //    bufferGrid1 = openvdb::tools::sdfToSdf
        //    (
        //        *bufferGrid1,
        //        currentVoxelSize*nCellsBetweenLevels
        //    );

        //    // Generate LOD sequence from finest level
        //    MultiResGrid fineToCoarse(/*levels*/2, *bufferGrid);

        //    bufferGrid1 = fineToCoarse.grid(1);
        //}
        if (cellLevel > 0)
        {
            //Timer timer("fineToCoarse levelSetRebuild");
            const scalar currentVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - cellLevel));
            const scalar coarseVoxelSize = voxelSize_ * (1 << (maxCellLevel_ - cellLevel + 1));

            openvdb::math::Transform::Ptr coarseTransform =
                openvdb::math::Transform::createLinearTransform
                (
                    coarseVoxelSize
                );

            transformOffset += 0.5 * (1.0/Foam::pow(2, bufferLevel - cellLevel + 1));

            openvdb::Vec3d offset(transformOffset, transformOffset, transformOffset);

            coarseTransform->preTranslate(offset);

            bufferGrid =
                openvdb::tools::levelSetRebuild<FloatGrid>
                (
                    *bufferGrid,
                    isoValues[cellLevel],
                    /*ext*/widths[cellLevel],
                    /*int*/0.1,
                    coarseTransform.get()
                );
        } //update zero level-set for grading

        //if (debug)
        {
            Info<< gridName << ": active voxels ";
            if (Pstream::parRun())
            {
                List<label> procActiveVoxels(Pstream::nProcs());
                procActiveVoxels[Pstream::myProcNo()] = grid->activeVoxelCount();
                IPstream::gatherList(procActiveVoxels);
                Info<< procActiveVoxels << endl;
            }
            else
            {
                Info<< grid->activeVoxelCount() <<endl;
            }
        }

        voxelCount += grid->activeVoxelCount();
    } //for cellLevel
} //calculateBufferGrids


FloatGrid::Ptr Foam::foamVDB::searchableBoxToLevelSet
(
    const point min,
    const point max,
    const label cellLevel
)
{
    //Timer timer;
    scalar edgeLength = Foam::pow(2, maxCellLevel_ - cellLevel);

    openvdb::math::Transform::Ptr linearTransform =
        openvdb::math::Transform::createLinearTransform
        (
            voxelSize_ * edgeLength
        );

    const openvdb::Vec3s pmin
    (
        linearTransform->worldToIndex
        (
            openvdb::Vec3s
            (
                min.x(),
                min.y(),
                min.z()
            )
        )
    );

    const openvdb::Vec3s pmax
    (
        linearTransform->worldToIndex
        (
            openvdb::Vec3s
            (
                max.x(),
                max.y(),
                max.z()
            )
        )
    );

    openvdb::CoordBBox boxBBox
    (
        openvdb::Coord
        (
            pmin.x(),
            pmin.y(),
            pmin.z()
        ),
        openvdb::Coord
        (
            pmax.x(),
            pmax.y(),
            pmax.z()
        )
    );

    FloatGrid::Ptr boxLevelSet = FloatGrid::create(0.0);
    boxLevelSet->setTransform(linearTransform);
    boxLevelSet->fill(boxBBox, /*value*/1.0);

    const word gridName = "background level " + std::to_string(cellLevel) + " cells";

    boxLevelSet->setName(gridName);

    Info<< nl << gridName<< ": active voxels ";
    if (Pstream::parRun())
    {
        List<label> procActiveVoxels(Pstream::nProcs());
        procActiveVoxels[Pstream::myProcNo()] = boxLevelSet->activeVoxelCount();
        IPstream::gatherList(procActiveVoxels);
        Info<< procActiveVoxels << endl;
    }
    else
    {
        Info<< boxLevelSet->activeVoxelCount() <<endl;
    }

    return boxLevelSet;
}

FloatGrid::Ptr Foam::foamVDB::meshToLevelSet
(
    const triSurface& inputSurf,
    const label cellLevel,
    const scalar halfWidth, /*= 3*/
    const int conversionFlags, /*= 0*/ //i.e. Signed Distance Field
    typename FloatGrid::template ValueConverter<openvdb::Int32>::Type * polygonIndexGrid,/*= nullptr*/
    typename FloatGrid::template ValueConverter<bool>::Type * maskGrid,                  /*= nullptr*/
    typename FloatGrid::template ValueConverter<openvdb::Vec4i>::Type * polygonListGrid, /*= nullptr*/
    typename FloatGrid::template ValueConverter<openvdb::Vec3d>::Type * displacementGrid /*= nullptr*/
)
{
    // The grid spacing
    const scalar delta = voxelSize_* Foam::pow(2, maxCellLevel_ - cellLevel);

    // A linear transform with the correct spacing
    openvdb::math::Transform::Ptr linearTransform =
        openvdb::math::Transform::createLinearTransform(delta);

    // The offset to cell-center points
    const scalar offset = -0.5*delta;

    triSurface idxSpaceSurf(inputSurf);
    // offset surface to be cell-centered
    idxSpaceSurf.translatePoints(point(offset));
    // transform points to local grid index space
    idxSpaceSurf.scalePoints(1.0 / delta);

    const MeshDataAdapter vdbSurfMesh(idxSpaceSurf);

    FloatGrid::Ptr surfGrid =
        //openvdb::tools::meshToVolume<FloatGrid>
        //(
        //    vdbSurfMesh,
        //    *linearTransform,
        //    /*exteriorBandWidth = */halfWidth,
        //    /*interiorBandWidth = */halfWidth,
        //    conversionFlags,
        //    polygonIndexGrid,
        //    maskGrid,
        //    polygonListGrid,
        //    displacementGrid
        //);
        openvdb::tools::meshToVolume<FloatGrid>
        (
            vdbSurfMesh,
            *linearTransform,
            /*exteriorBandWidth = */halfWidth,
            /*interiorBandWidth = */halfWidth,
            conversionFlags,
            polygonIndexGrid
        );

    return surfGrid;
}

//void Foam::foamVDB::offsetLevelSetValues
//(
//    const FloatGrid::Ptr& grid,
//    const scalar& offset
//)
//{
//    struct OffsetOp
//    {
//        scalar offset_;
//
//        OffsetOp
//        (
//            const scalar& offset
//        )
//        :
//            offset_(offset)
//        {}
//
//        inline void operator()(const FloatGrid::ValueOnIter& iter) const
//        {
//            iter.setValue(*iter + offset_);
//        }
//    };
//
//    openvdb::tools::foreach(grid->beginValueOn(), OffsetOp(offset));
//}

bool Foam::foamVDB::isOpenSurface
(
    const FloatGrid::Ptr& grid
)
{
    std::atomic<bool> openSurface;

    struct IsOpenOp
    {
        std::atomic<bool> * const openSurface_;

        IsOpenOp
        (
            std::atomic<bool>& openSurface
        )
        :
            openSurface_(&openSurface)
        {}

        inline void operator()(const FloatGrid::ValueOnCIter& iter) const
        {
            if (!*openSurface_) return;

            if (iter.getValue() < 0)
            {
                *openSurface_ = false;
                return;
            }
        }
    };

    openSurface = true;

    openvdb::tools::foreach(grid->cbeginValueOn(), IsOpenOp(openSurface));

    if (openSurface)
    {
        Info<<"Surface is open! "<<endl;
    }
    else
    {
        Info<<"Surface is closed! "<<endl;
    }

    return bool(openSurface);
} // is OpenSurface


FloatGrid::Ptr Foam::foamVDB::topologyUnionLevelSets
(
    const std::vector<std::vector<FloatGrid::Ptr>>& boundaryCellLevels,
    const label cellLevel
)
{
    label a = 0;

    //find first valid grid
    for (size_t i = 0; i < boundaryCellLevels.size(); ++i)
    {
        if (boundaryCellLevels[i][cellLevel])
        {
            a = i;
            break;
        }
    }

    FloatGrid::Ptr topoUnion = boundaryCellLevels[a][cellLevel];

    //add valid grids
    for (size_t i = a + 1; i < boundaryCellLevels.size(); ++i)
    {
        if (boundaryCellLevels[i][cellLevel])
        {
            topoUnion->topologyUnion(*boundaryCellLevels[i][cellLevel]);
        }
    }

    return topoUnion;
} // topologyUnionLevelSets


FloatGrid::Ptr Foam::foamVDB::csgUnionLevelSets
(
    const std::vector<std::vector<FloatGrid::Ptr>>& boundaryCellLevels,
    const label cellLevel,
    const bool convertToSDF, /*= false*/
    const scalar halfWidth /*= 3*/
)
{
    std::vector<FloatGrid::Ptr> validBoundaries;

    for (size_t i = 0; i < boundaryCellLevels.size(); ++i)
    {
        if (boundaryCellLevels[i][cellLevel])
        {
            validBoundaries.push_back(boundaryCellLevels[i][cellLevel]->deepCopy());
        }
    }

    label nUnion = ceil(validBoundaries.size() / 2.0);

    std::vector<FloatGrid::Ptr> unionGrids(nUnion);

    for (size_t i = 0; i < unionGrids.size(); ++i)
    {
        uLabel regioni = 2*i;

        if (regioni + 1 == validBoundaries.size())
        {
            //copy last if size is not even number
            unionGrids[i] =
                //boundaryCellLevels[regioni][cellLevel]->deepCopy();
                validBoundaries[regioni];
        }
        else
        {
            const label cellCountA = validBoundaries[regioni]->activeVoxelCount();
            const label cellCountB = validBoundaries[regioni+1]->activeVoxelCount();

            if (cellCountA > 0 && cellCountB > 0)
            {
                unionGrids[i] =
                    openvdb::tools::csgUnionCopy<FloatGrid>
                    (
                        *validBoundaries[regioni],
                        *validBoundaries[regioni + 1]
                    );
            }
            else if (cellCountA == 0) // or both have 0 voxels, will clean in the next loop
            {
                unionGrids[i] =
                    validBoundaries[regioni + 1];
            }
            else // cellCountB == 0
            {
                unionGrids[i] =
                    validBoundaries[regioni];
            }
        }
    }

    // find first unionGrid with nonzero active voxels
    label a = 0;
    for (size_t i = 0; i < unionGrids.size(); ++i)
    {
        if (unionGrids[i]->activeVoxelCount() > 0)
        {
            a = i;
            break;
        }
    }

    for (size_t i = a + 1; i < unionGrids.size(); ++i)
    {
        if (unionGrids[i]->activeVoxelCount() > 0)
        {
            openvdb::tools::csgUnion
            (
                *unionGrids[a],
                *unionGrids[i] // the operation always leaves this empty
            );
        }
    }

    if (convertToSDF)
    {
        label edgeLength = Foam::pow(2, maxCellLevel_ - cellLevel);

        const scalar newZero = voxelSize_ * edgeLength;

        FloatGrid::Ptr grid =
            openvdb::tools::sdfToSdf
            (
                *unionGrids[0],
                newZero
            );

        //std::cout << "SDF Diagnostics:\n" << openvdb::tools::checkLevelSet(*grid, 9) << std::endl;
        //std::cout << "unionGrid Diagnostics:\n" << openvdb::tools::checkLevelSet(*unionGrids[0], 9) << std::endl;

        BoolGrid::Ptr maskGrid = openvdb::tools::extractEnclosedRegion(*grid);

        maskGrid->tree().voxelizeActiveTiles();
        openvdb::tools::erodeActiveValues(maskGrid->tree(), 1);

        { Timer timer("csgUnionLevelSets: topologyToLevelSet");
            /*FloatGrid::Ptr*/ grid = openvdb::tools::topologyToLevelSet
                (
                    *maskGrid,
                    halfWidth,
                    /*closingSteps*/0,
                    /*dilation*/0,
                    /*smoothingSteps*/0
                );
        }
        //std::cout << "voxelSize " << grid->voxelSize().x()<< std::endl;
        //std::cout << "topologyToLevelSet Diagnostics:\n" << openvdb::tools::checkLevelSet(*grid, 9) << std::endl;

        isOpenSurface(grid);

        return grid;
    }
    else
    {
        return unionGrids[a];
    }
} // csgUnionLevelSets


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * Friend Operators * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// ************************************************************************* //
