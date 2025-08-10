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
    (c) 2022 ESI

\*---------------------------------------------------------------------------*/

#include "foamVDB.H"

#include <openvdb/tools/Morphology.h> //for dilation and erosion
#include <openvdb/tools/Composite.h>

#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"

namespace Foam
{
    static bool less(const point& a, const point& b)
    {
        if (a.x() == b.x())
        {
            if (a.y() == b.y())
            {
                return (a.z() < b.z());
            }
            return (a.y() < b.y());
        }
        return (a.x() < b.x());
    }

    //- To compare points
    class pointLess
    {
        const pointList& values_;

        public:

        pointLess(const pointList& values)
            :
                values_(values)
        {}

        bool operator()(const label a, const label b) const
        {
            return less(values_[a], values_[b]);
        }
    };

    //transform points from index space to world space
    struct TransformPoints
    {
        const std::vector<vdbPoint>& pointsIn_;
        pointField& pointsOut_;
        const point& levelSetOffset_;
        const point& scale_;
        const quaternion& R_;
        const scalar voxelSize_;

        TransformPoints
        (
            const std::vector<vdbPoint>& pointsIn,
            pointField& pointsOut,
            const point& levelSetOffset,
            const point& scale,
            const quaternion& R,
            const scalar voxelSize
        )
        :
            pointsIn_(pointsIn),
            pointsOut_(pointsOut),
            levelSetOffset_(levelSetOffset),
            scale_(scale),
            R_(R),
            voxelSize_(voxelSize)
        {
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t n = range.begin(); n < range.end(); ++n)
            {
                if (!pointsIn_[n].active_) continue;

                const label id = pointsIn_[n].id_;
                pointsOut_[id] =
                    R_.transform
                    (
                        cmptMultiply
                        (
                            (((pointsIn_[n].ijk_ / /*pointsPerEdge*/2) * voxelSize_) - levelSetOffset_),
                            scale_
                        )
                    );
            }
        }
    }; //TransformPoints


    struct SetFineVoxelsOff
    {
        IndexGrid::ConstAccessor acc_;

        SetFineVoxelsOff(IndexGrid& coarseIDGrid)
        :
            acc_(coarseIDGrid.getConstAccessor())
        {}

        inline void operator()(const BoolGrid::ValueOnIter& iter) const
        {
            const openvdb::Coord& ijk = iter.getCoord();

            openvdb::Int32 coarseCellID;

            switch ( (ijk[0] & 1) | ((ijk[1] & 1) << 1) | ((ijk[2] & 1) << 2) )
            {
                case 0:// all even
                    coarseCellID = acc_.getValue(ijk>>1);
                    break;
                case 1:// x is odd
                    coarseCellID = acc_.getValue(ijk.offsetBy(-1,0,0)>>1);
                    break;
                case 2:// y is odd
                    coarseCellID = acc_.getValue(ijk.offsetBy(0,-1,0)>>1);
                    break;
                case 3:// x&y are odd
                    coarseCellID = acc_.getValue(ijk.offsetBy(-1,-1,0)>>1);
                    break;
                case 4:// z is odd
                    coarseCellID = acc_.getValue(ijk.offsetBy(0,0,-1)>>1);
                    break;
                case 5:// x&z are odd
                    coarseCellID = acc_.getValue(ijk.offsetBy(-1,0,-1)>>1);
                    break;
                case 6:// y&z are odd
                    coarseCellID = acc_.getValue(ijk.offsetBy(0,-1,-1)>>1);
                    break;
                case 7:// all are odd
                    coarseCellID = acc_.getValue(ijk.offsetBy(-1,-1,-1)>>1);
                    break;
            } //switch

            if (coarseCellID >= 0)
            {
                iter.setValueOff();
            }
        }
    }; // SetFineVoxelsOff


    struct SetProcID
    {
        BoolGrid::ConstAccessor acc_;
        label neighProcID_;

        SetProcID
        (
            BoolGrid& coarseIDGrid,
            const label neighProcID
        )
        :
            acc_(coarseIDGrid.getConstAccessor()),
            neighProcID_(neighProcID)
        {}

        inline void operator()(const BoolGrid::ValueOnCIter& iter, IndexGrid::Accessor& myLevelNcellIDAccessor) const
        {
            const openvdb::Coord& ijk = iter.getCoord();

            bool neighProcIsOn = false;

            switch ( (ijk[0] & 1) | ((ijk[1] & 1) << 1) | ((ijk[2] & 1) << 2) )
            {
                case 0:// all even
                    neighProcIsOn = acc_.isValueOn(ijk>>1);
                    break;
                case 1:// x is odd
                    neighProcIsOn = acc_.isValueOn(ijk.offsetBy(-1,0,0)>>1);
                    break;
                case 2:// y is odd
                    neighProcIsOn = acc_.isValueOn(ijk.offsetBy(0,-1,0)>>1);
                    break;
                case 3:// x&y are odd
                    neighProcIsOn = acc_.isValueOn(ijk.offsetBy(-1,-1,0)>>1);
                    break;
                case 4:// z is odd
                    neighProcIsOn = acc_.isValueOn(ijk.offsetBy(0,0,-1)>>1);
                    break;
                case 5:// x&z are odd
                    neighProcIsOn = acc_.isValueOn(ijk.offsetBy(-1,0,-1)>>1);
                    break;
                case 6:// y&z are odd
                    neighProcIsOn = acc_.isValueOn(ijk.offsetBy(0,-1,-1)>>1);
                    break;
                case 7:// all are odd
                    neighProcIsOn = acc_.isValueOn(ijk.offsetBy(-1,-1,-1)>>1);
                    break;
            } //switch

            if (neighProcIsOn)
            {
                myLevelNcellIDAccessor.setValueOn(ijk, neighProcID_);
            }
        }
    }; // SetProcID


    struct SetFineProcID
    {
        BoolGrid::ConstAccessor acc_;
        label neighProcID_;

        SetFineProcID
        (
            BoolGrid& neighProcGrid,
            const label neighProcID
        )
        :
            acc_(neighProcGrid.getConstAccessor()),
            neighProcID_(neighProcID)
        {}

        inline void operator()(const BoolGrid::ValueOnCIter& iter, IndexGrid::Accessor& myLevelNcellIDAccessor) const
        {
            //scale coord to finer level
            openvdb::Coord ijk = iter.getCoord() << 1;

            openvdb::Coord ijkFine;

            for (label i = 0; i < 8; ++i)
            {
                ijkFine = ijk + COARSE_TO_FINE[i];

                if (acc_.isValueOn(ijkFine))
                {
                    myLevelNcellIDAccessor.setValueOn(ijkFine, neighProcID_);
                }
            }

        }
    }; // SetFineProcID
} //namespace Foam


bool Foam::foamVDB::addDisplacement
(
    label cellLevel,
    const face& f,
    std::vector<openvdb::Coord>& faceCoords,
    labelHashSet& displacedPoints,
    std::vector<point>& xyzPoints,
    std::vector<vdbPoint>& vdbPoints,
    Vec3dGrid::ConstAccessor& displAcc,
    Vec3dGrid::ConstAccessor& fineDisplAcc,
    const bool isSplitFace /*= false*/
)
{
    return true; //GGG

    if (isSplitFace || f.size() > 4)
    {
        cellLevel += 1;
    }

    openvdb::Vec3d displ;

    forAll(f, pointI)
    {
        const label& pointID = f[pointI];

        if (!displacedPoints.found(pointID))
        {
            const point& p = xyzPoints[pointID];

            openvdb::Coord ijk //divide by pointsPerEdge
            (
                (openvdb::Int32(p.x()/2) >> (maxCellLevel_ - cellLevel)),
                (openvdb::Int32(p.y()/2) >> (maxCellLevel_ - cellLevel)),
                (openvdb::Int32(p.z()/2) >> (maxCellLevel_ - cellLevel))
            );

            //displacement in indexSpace (times pointsPerEdge)
            if (isSplitFace || f.size() > 4)
            {
                displ = 2*(fineDisplAcc.getValue(ijk) - openvdb::Vec3d(0.5));
            }
            else
            {
                displ = 2*(displAcc.getValue(ijk) - openvdb::Vec3d(0.5));
            }

            //convert to world space
            displ *= Foam::pow(2, (maxCellLevel_ - cellLevel));

            xyzPoints[pointID] -= point(displ.x(), displ.y(), displ.z());
            vdbPoints[pointID].ijk_ -= point(displ.x(), displ.y(), displ.z());

            displacedPoints.insert(pointID);
        }
    }

    //GGG
    if (f.size() != 4) return true;

    // Calculate face area
    // Compute an estimate of the centre as the average of the points
    point pAvg = Zero;
    forAll(f, pI)
    {
        pAvg += xyzPoints[f[pI]];
    }
    pAvg /= f.size();

    // Compute the face area normal and unit normal by summing up the
    // normals of the triangles formed by connecting each edge to the
    // point average.
    vector sumA = Zero;
    forAll(f, pI)
    {
        const point& p =     xyzPoints[f[pI]];
        const point& pNext = xyzPoints[f[(pI == f.size()-1 ? 0 : pI+1)]];

        const vector a = (pNext - p)^(pAvg - p);

        sumA += a;
    }
    scalar faceArea = mag(0.5*sumA);

    scalar voxelEdge = 2 << (maxCellLevel_ - cellLevel);
    scalar voxelFaceArea = Foam::pow(voxelEdge, 2);

    if (faceArea < 0.5*voxelFaceArea)
    {
        Pout<<"GGG collapsing face!"<<endl;
        forAll(f, pI)
        {
            const point& p = xyzPoints[f[pI]];
            for (label next = 3; next > pI; --next)
            {
                const point& pNext = xyzPoints[f[next]];
                scalar delta = mag(p - pNext);

                if (delta < 0.5*scalar(voxelEdge))
                {
                    Pout<<"    GGG collapsing edge "<< f[pI]<<"-"<<f[next]<< endl;
                    point mid = (vdbPoints[f[pI]].ijk_ + vdbPoints[f[next]].ijk_)/2;

                    label id = vdbPoints[f[pI]].id_;
                    if (vdbPoints[f[next]].id_ < id)
                    {
                        id = vdbPoints[f[next]].id_;
                        vdbPoints[f[pI]].active_ = false;
                        vdbPoints[f[pI]].id_ = id;
                        vdbPoints[f[next]].ijk_ = mid;
                    }
                    else
                    {
                        vdbPoints[f[next]].active_ = false;
                        vdbPoints[f[next]].id_ = id;
                        vdbPoints[f[pI]].ijk_ = mid;
                    }

                    xyzPoints[f[pI]] = mid;
                    xyzPoints[f[next]] = mid;
                    break;
                }
            }
        }

        label nInactive = 0;
        forAll(f, pI)
        {
            if (!vdbPoints[f[pI]].active_)
            {
                nInactive++;
            }
        }
        if (nInactive > 1)
        {
            return false;
        }
    }
    return true;
} // addDisplacement

void Foam::foamVDB::vdbGridsToPolyMesh
(
    const std::vector<FloatGrid::Ptr>& cellLevelGrids,
    const std::vector<Vec3dGrid::Ptr>& displacementGrids,
    const std::vector<BoolGrid::Ptr>& innerGrids,
    const std::vector<BoolGrid::Ptr>& outerGrids,
    const label& nRegions,
    const point& levelSetOffset,
    const point& scale,
    const quaternion& R,
    const openvdb::CoordBBox& boundingBox, //level0 bounds
    const dictionary& vdbDict,

    fvMesh&    mesh,
    labelList& hexRef8cellLevel,
    labelList& hexRef8pointLevel
)
{
    std::vector<vdbPoint>             vdbPoints;
    std::vector<vdbFace>              vdbInternalFaces;
    std::vector<std::vector<vdbFace>> vdbBoundaryFaces(nRegions,          std::vector<vdbFace>());
    std::vector<std::vector<vdbFace>> vdbProcFaces    (Pstream::nProcs(), std::vector<vdbFace>());

    label nInternal = 0;
    label nBoundaryFace = 0;
    labelList nProcFace(Pstream::nProcs(), 0);


    pointField points;
    DynamicList<face> faces;
    DynamicList<label> owner;
    DynamicList<label> neighbour;
    List<label> boundaryMeshSizes(nRegions, 0);
    List<label> boundaryMeshStarts(nRegions, 0);

    const bool subsetInnerGrid(vdbDict.lookupOrDefault<bool>("subsetInnerGrid", false));
    const bool subsetOuterGrid(vdbDict.lookupOrDefault<bool>("subsetOuterGrid", false));
    const bool balanceSurfaceCells(vdbDict.lookupOrDefault<bool>("balanceSurfaceCells", true));
    const bool threaded(vdbDict.lookupOrDefault<bool>("threadedVDBtoPolyMesh", true));

    std::vector<label> innerSubset;
    std::vector<label> outerSubset;
    if (!subsetInnerGrid && !subsetOuterGrid)
    {
        const label nInner = getActiveVoxels<BoolGrid>(innerGrids);
        innerSubset.reserve(nInner);

        const label nOuter = getActiveVoxels<BoolGrid>(outerGrids);
        outerSubset.reserve(nOuter);
    }

    //fields
    std::vector<scalar> cellLevelField;
    std::vector<scalar> cellIDField;
    std::vector<scalar> procIDField;
    std::vector<scalar> distanceField;

    std::vector<IndexGrid::Ptr> cellIDGrids    (maxCellLevel_ + 1, nullptr);
    std::vector<IndexGrid::Ptr> pointIDGrids   (maxCellLevel_ + 1, nullptr);
    std::vector<IndexGrid::Ptr> pointLevelGrids(maxCellLevel_ + 1, nullptr);

    // topology of neighbour procs
    std::vector<BoolGrid::Ptr> coarseNeighProc(Pstream::nProcs(), nullptr);

    BoolGrid::Ptr coarseHaloGrid = BoolGrid::create(-labelMax);

    openvdb::Int32 nCells = 0;

    { Timer timer("populate cellIDGrids");
    // List[proci][cellLevel]
    List<List<boundBox>> procBBoxes =
        sendReceiveBoundBoxes<FloatGrid>(cellLevelGrids);

    // populate cellIDGrids
    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        const label edgeLength = 1 << (maxCellLevel_ - cellLevel); //in (smallest) voxel units

        openvdb::math::Transform::Ptr levelNTransform =
            openvdb::math::Transform::createLinearTransform
            (
                edgeLength
            );

        // pointIDs
        IndexGrid::Ptr levelNpointIDGrid = IndexGrid::create(-1);
        levelNpointIDGrid->setTransform(levelNTransform);
        levelNpointIDGrid->setGridClass(openvdb::GRID_UNKNOWN);
        pointIDGrids[cellLevel] = levelNpointIDGrid;

        // pointLevel
        IndexGrid::Ptr levelNpointLevelGrid = IndexGrid::create(-1);
        levelNpointLevelGrid->setTransform(levelNTransform);
        levelNpointLevelGrid->setGridClass(openvdb::GRID_UNKNOWN);
        pointLevelGrids[cellLevel] = levelNpointLevelGrid;

        // cellIDs
        IndexGrid::Ptr levelNcellIDGrid = IndexGrid::create(-labelMax);
        levelNcellIDGrid->setTransform(levelNTransform);
        levelNcellIDGrid->setGridClass(openvdb::GRID_UNKNOWN);

        IndexGrid::Accessor levelNcellIDAccessor = levelNcellIDGrid->getAccessor();

        FloatGrid::Ptr levelNGrid = cellLevelGrids[cellLevel];
        levelNGrid->tree().voxelizeActiveTiles();

        //TODO remove this loop
        // switch off values of finer grid in case they overlap
        if (cellLevel < maxCellLevel_)
        {
            FloatGrid::Ptr fineCellIDGrid = cellLevelGrids[cellLevel + 1];

            FloatGrid::Accessor acc = fineCellIDGrid->getAccessor();

            openvdb::Int32 neighID;
            for (auto iter = levelNGrid->cbeginValueOn(); iter; ++iter)
            {
                //scale coord to finer level
                openvdb::Coord ijk = iter.getCoord() << 1;

                for (label i = 0; i < 8; ++i)
                {
                    openvdb::Coord ijkFine = ijk + COARSE_TO_FINE[i];

                    if (acc.isValueOn(ijkFine))
                    {
                        acc.setValueOff(ijkFine);

                        ijkFine <<= (maxCellLevel_ - cellLevel);

                        point p(ijkFine.x(), ijkFine.y(), ijkFine.z());

                        Pout<< "WARNING! At cellevel " << cellLevel + 1
                            << " Found overlapping voxel at " << p;

                        p = (((p / /*pointsPerEdge*/2) * voxelSize_) - levelSetOffset);
                        p = cmptMultiply(scale, p);

                        Pout<< " - world coord " << p << endl;
                    }
                }
            }
        }

        // switch off voxels inside coarse neighbour proc
        if (cellLevel > 0)
        {
            for (label proci = 0; proci < Pstream::nProcs(); ++proci)
            {
                if (proci == Pstream::myProcNo()) continue;

                // coarseNeighProc saved from previous cellLevel iteration
                if (!coarseNeighProc[proci]) continue;

                if (threaded) //1
                {
                    using IterRange = typename openvdb::tree::IteratorRange<BoolGrid::ValueOnCIter>;
                    tbb::parallel_for
                    (
                        IterRange(coarseNeighProc[proci]->cbeginValueOn()),
                        [&](IterRange r)
                        {
                            openvdb::Coord ijk, ijkFine;

                            FloatGrid::Accessor acc = cellLevelGrids[cellLevel]->getAccessor();

                            for (; r; ++r)
                            {
                                auto iter = r.iterator();
                                const openvdb::Coord& ijk0 = iter.getCoord();

                                //scale coord to finer level
                                ijk = ijk0 << 1;

                                for (label i = 0; i < 8; ++i)
                                {
                                    ijkFine = ijk + COARSE_TO_FINE[i];

                                    acc.setValueOff(ijkFine);
                                }
                            }
                        }
                    );
                }
                else
                {
                    FloatGrid::Accessor acc = cellLevelGrids[cellLevel]->getAccessor();
                    for (auto iter = coarseNeighProc[proci]->cbeginValueOn(); iter; ++iter)
                    {
                        //scale coord to finer level
                        openvdb::Coord ijk = iter.getCoord() << 1;

                        for (label i = 0; i < 8; ++i)
                        {
                            openvdb::Coord ijkFine = ijk + COARSE_TO_FINE[i];

                            acc.setValueOff(ijkFine);
                        }
                    }
                }
            }
        }

        // set cellID (not thread-safe leave it serial)
        for (auto iter = levelNGrid->cbeginValueOn(); iter; ++iter)
        {
            levelNcellIDAccessor.setValueOn(iter.getCoord(), nCells++);
        }

        // Original copy of levelNcellIDGrid (without halo voxels)
        BoolGrid::Ptr interiorTopology = BoolGrid::create();

        interiorTopology->topologyUnion(*levelNcellIDGrid);

        //create halo grid
        BoolGrid::Ptr haloGrid = BoolGrid::create();

        haloGrid->topologyUnion(*interiorTopology);

        openvdb::tools::dilateActiveValues(haloGrid->tree(), 1, openvdb::tools::NN_FACE_EDGE);

        haloGrid->topologyDifference(*interiorTopology);
        openvdb::tools::pruneInactive(haloGrid->tree());

        //add neighbour procID of coarser cellLevel
        if (cellLevel > 0)
        {
            //remove haloGrid voxels overlapping with active coarse voxels
            //if (threaded) //2
            {
                SetFineVoxelsOff op(*cellIDGrids[cellLevel - 1]);
                openvdb::tools::foreach
                (
                    haloGrid->beginValueOn(),
                    op,
                    threaded,
                    false //DO NOT SHARE OPERATOR as grid accessor does caching...
                );
            }
            //else
            //{
            //    IndexGrid::Ptr coarseCellIDGrid = cellIDGrids[cellLevel - 1];

            //    IndexGrid::ConstAccessor acc = coarseCellIDGrid->getConstAccessor();

            //    for (auto iter = haloGrid->beginValueOn(); iter; ++iter)
            //    {
            //        const openvdb::Coord& ijk = iter.getCoord();

            //        openvdb::Int32 coarseCellID;

            //        switch ( (ijk[0] & 1) | ((ijk[1] & 1) << 1) | ((ijk[2] & 1) << 2) )
            //        {
            //            case 0:// all even
            //                coarseCellID = acc.getValue(ijk>>1);
            //                break;
            //            case 1:// x is odd
            //                coarseCellID = acc.getValue(ijk.offsetBy(-1,0,0)>>1);
            //                break;
            //            case 2:// y is odd
            //                coarseCellID = acc.getValue(ijk.offsetBy(0,-1,0)>>1);
            //                break;
            //            case 3:// x&y are odd
            //                coarseCellID = acc.getValue(ijk.offsetBy(-1,-1,0)>>1);
            //                break;
            //            case 4:// z is odd
            //                coarseCellID = acc.getValue(ijk.offsetBy(0,0,-1)>>1);
            //                break;
            //            case 5:// x&z are odd
            //                coarseCellID = acc.getValue(ijk.offsetBy(-1,0,-1)>>1);
            //                break;
            //            case 6:// y&z are odd
            //                coarseCellID = acc.getValue(ijk.offsetBy(0,-1,-1)>>1);
            //                break;
            //            case 7:// all are odd
            //                coarseCellID = acc.getValue(ijk.offsetBy(-1,-1,-1)>>1);
            //                break;
            //        } //switch

            //        if (coarseCellID >= 0)
            //        {
            //            iter.setValueOff();
            //        }
            //    } //for haloGrid voxel (turnoff fine voxels overlapping coarse ones)
            //}

            //store in halo neighbour procID of coarser cellLevel
            // (fine voxel checks for coarse processor neighbours)
            for (label proci = 0; proci < Pstream::nProcs(); ++proci)
            {
                if (proci == Pstream::myProcNo()) continue;

                // coarseNeighProc saved from previous cellLevel iteration
                if (!coarseNeighProc[proci]) continue;

                // neighbour processor ID is stored as negative procNo minus 1
                // e.g. proc0 -> -1
                //      proc1 -> -2
                //      proc15-> -16
                openvdb::Int32 neighProcID = -(proci + 1);

                //if (threaded) //3
                {
                    SetProcID op(*coarseNeighProc[proci], neighProcID);

                    IndexGrid::Ptr addedIDGrid = IndexGrid::create(-labelMax);

                    openvdb::tools::transformValues
                    (
                        haloGrid->cbeginValueOn(),
                        *addedIDGrid,
                        op,
                        threaded,
                        false //DO NOT SHARE OPERATOR as grid accessor does caching...
                    );

                    openvdb::tools::compReplace(levelNcellIDGrid->tree(), addedIDGrid->tree());
                }
                //else
                //{
                //    BoolGrid::ConstAccessor acc = coarseNeighProc[proci]->getConstAccessor();

                //    for (auto iter = haloGrid->cbeginValueOn(); iter; ++iter)
                //    {
                //        const openvdb::Coord& ijk = iter.getCoord();

                //        bool neighProcIsOn = false;

                //        switch ( (ijk[0] & 1) | ((ijk[1] & 1) << 1) | ((ijk[2] & 1) << 2) )
                //        {
                //            case 0:// all even
                //                neighProcIsOn = acc.isValueOn(ijk>>1);
                //                break;
                //            case 1:// x is odd
                //                neighProcIsOn = acc.isValueOn(ijk.offsetBy(-1,0,0)>>1);
                //                break;
                //            case 2:// y is odd
                //                neighProcIsOn = acc.isValueOn(ijk.offsetBy(0,-1,0)>>1);
                //                break;
                //            case 3:// x&y are odd
                //                neighProcIsOn = acc.isValueOn(ijk.offsetBy(-1,-1,0)>>1);
                //                break;
                //            case 4:// z is odd
                //                neighProcIsOn = acc.isValueOn(ijk.offsetBy(0,0,-1)>>1);
                //                break;
                //            case 5:// x&z are odd
                //                neighProcIsOn = acc.isValueOn(ijk.offsetBy(-1,0,-1)>>1);
                //                break;
                //            case 6:// y&z are odd
                //                neighProcIsOn = acc.isValueOn(ijk.offsetBy(0,-1,-1)>>1);
                //                break;
                //            case 7:// all are odd
                //                neighProcIsOn = acc.isValueOn(ijk.offsetBy(-1,-1,-1)>>1);
                //                break;
                //        } //switch

                //        if (neighProcIsOn)
                //        {
                //            levelNcellIDAccessor.setValueOn(ijk, neighProcID);
                //        }
                //    } // for voxel of dilated grid
                //}
            } //for neighbour proc
        } //add neigh procID of coarser cellLevel


        // send my cellIDGrid if myProcBB overlaps
        // with halo (+ coarse halo) of procN
        PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
        PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

        boundBox myProcBB = procBBoxes[Pstream::myProcNo()][cellLevel];

        // send to overlapping proc
        for (label proci = 0; proci < Pstream::nProcs(); ++proci)
        {
            if (proci == Pstream::myProcNo()) continue;

            boundBox otherProcBB = procBBoxes[proci][cellLevel];

            //expand otherProcBB to include coarser cellLevel
            if (cellLevel > 0)
            {
                boundBox coarseOtherProcBB = procBBoxes[proci][cellLevel - 1];
                //scale to current cellLevel
                coarseOtherProcBB.min() *= 2.0;
                coarseOtherProcBB.max() *= 2.0;

                otherProcBB.add(coarseOtherProcBB);
            }

            //expand otherProcBB to include finer cellLevel
            if (cellLevel < maxCellLevel_)
            {
                boundBox fineOtherProcBB = procBBoxes[proci][cellLevel + 1];
                //scale to current cellLevel
                fineOtherProcBB.min() /= 2.0;
                fineOtherProcBB.max() /= 2.0;

                otherProcBB.add(fineOtherProcBB);
            }

            //inflate to include halo and find overlap with myProc voxels
            otherProcBB.inflateAbs(1.0);

            if (myProcBB.overlaps(otherProcBB))
            {
                //send the original cellID grid (without halo values)
                sendGrid<BoolGrid>
                (
                    interiorTopology,
                    pBufSize,
                    pBufGrid,
                    proci
                );
            }
        } //for proci (send)

        pBufSize.finishedSends();
        pBufGrid.finishedSends();

        //expand myProcBB to include coarser cellLevel
        if (cellLevel > 0)
        {
            boundBox coarseMyProcBB = procBBoxes[Pstream::myProcNo()][cellLevel - 1];
            //scale to current cellLevel
            coarseMyProcBB.min() *= 2.0;
            coarseMyProcBB.max() *= 2.0;

            myProcBB.add(coarseMyProcBB);
        }
        //expand myProcBB to include finer cellLevel
        if (cellLevel < maxCellLevel_)
        {
            boundBox fineMyProcBB = procBBoxes[Pstream::myProcNo()][cellLevel + 1];
            //scale to current cellLevel
            fineMyProcBB.min() /= 2.0;
            fineMyProcBB.max() /= 2.0;

            myProcBB.add(fineMyProcBB);
        }

        //inflate to include halo and find overlap with procN voxels
        myProcBB.inflateAbs(1.0);

        //receive from overlapping procs
        for (label proci = 0; proci < Pstream::nProcs(); ++proci)
        {
            if (proci == Pstream::myProcNo()) continue;

            const boundBox& otherProcBB = procBBoxes[proci][cellLevel];

            if (coarseNeighProc[proci]) coarseNeighProc[proci]->clear();

            if (myProcBB.overlaps(otherProcBB))
            {
                BoolGrid::Ptr boundaryProcI =
                    receiveGrid<BoolGrid>
                    (
                        pBufSize,
                        pBufGrid,
                        proci
                    );

                openvdb::Int32 neighProcID = -(proci + 1);

                if (boundaryProcI->activeVoxelCount() > 0)
                {
                    if (threaded) //4
                    {
                        IndexGrid::Ptr addedIDGrid = IndexGrid::create(neighProcID);

                        addedIDGrid->topologyUnion(*haloGrid);

                        addedIDGrid->topologyIntersection(*boundaryProcI);

                        haloGrid->topologyDifference(*addedIDGrid);

                        openvdb::tools::pruneInactive(haloGrid->tree());

                        openvdb::tools::compReplace(levelNcellIDGrid->tree(), addedIDGrid->tree());
                    }
                    else
                    {
                        BoolGrid::Ptr myProcNeighbours = haloGrid->deepCopy();

                        myProcNeighbours->topologyIntersection(*boundaryProcI);

                        for (auto iter = myProcNeighbours->cbeginValueOn(); iter; ++iter)
                        {
                            levelNcellIDAccessor.setValueOn(iter.getCoord(), neighProcID);
                        }

                        haloGrid->topologyDifference(*myProcNeighbours);
                        openvdb::tools::pruneInactive(haloGrid->tree());
                    }

                }

                // check if received voxels from procN are neighbours with coarser voxels in myProc
                // if so, add procID to voxels at cellLevel interface
                // (corse voxel checks for fine processor neighbours)
                if (cellLevel > 0)
                {
                    //if (threaded) //5
                    {
                        SetFineProcID op(*boundaryProcI, neighProcID);

                        IndexGrid::Ptr addedIDGrid = IndexGrid::create(-labelMax);

                        // coarseHaloGrid saved from previous cellLevel iteration
                        openvdb::tools::transformValues
                        (
                            coarseHaloGrid->cbeginValueOn(),
                            *addedIDGrid,
                            op,
                            threaded,
                            false //DO NOT SHARE OPERATOR as grid accessor does caching...
                        );

                        openvdb::tools::compReplace(levelNcellIDGrid->tree(), addedIDGrid->tree());
                    }
                    //else
                    //{
                    //    BoolGrid::ConstAccessor neighProcAcc = boundaryProcI->getConstAccessor();

                    //    openvdb::Coord ijk, ijkFine;

                    //    // coarseHaloGrid saved from previous cellLevel iteration
                    //    for (auto iter = coarseHaloGrid->cbeginValueOn(); iter; ++iter)
                    //    {
                    //        //scale coord to fine level
                    //        ijk = iter.getCoord() << 1;

                    //        for (label i = 0; i < 8; ++i)
                    //        {
                    //            ijkFine = ijk + COARSE_TO_FINE[i];

                    //            if (neighProcAcc.isValueOn(ijkFine))
                    //            {
                    //                levelNcellIDAccessor.setValueOn(ijkFine, neighProcID);
                    //            }
                    //        }
                    //    }
                    //}
                }

                // save coarse procI to check at next iteration
                // (fine voxel checks for coarse processor neighbours)
                if (cellLevel < maxCellLevel_)
                {
                    coarseNeighProc[proci] = boundaryProcI;
                }
            } // if proci overlaps myProc
            else
            {
                coarseNeighProc[proci] = nullptr;
            }
        } //for proci (receive)

        cellIDGrids[cellLevel] = levelNcellIDGrid;

        coarseHaloGrid->clear();

        if (cellLevel < maxCellLevel_)
        {
            coarseHaloGrid = haloGrid;
        }
    } //populate cellIDGrids
    } // Timer

    std::vector<point> xyzPoints;
    std::vector<label> pointLevel;
    xyzPoints.reserve(6*nCells); //first estimate of point count
    pointLevel.reserve(6*nCells);
    labelHashSet displacedPoints;

    //processorBoundaries
    List<DynamicList<point>> procsFaceCentres(Pstream::nProcs());
    List<DynamicList<face>> procsFaces(Pstream::nProcs());
    List<DynamicList<label>> procOwners(Pstream::nProcs());
    List<word> procBoundaryNames(Pstream::nProcs());
    List<label> procBoundarySizes(Pstream::nProcs(), 0);

    forAll(procBoundaryNames, i)
    {
        procBoundaryNames[i] =
            "procBoundary" + Foam::name(Pstream::myProcNo()) + "to" + Foam::name(i);
    }

    //boundary patches
    List<DynamicList<face>> patchFaces(nRegions);
    List<DynamicList<label>> patchOwners(nRegions);

    label interfaceSize = 0;
    DynamicList<face>  interfaceFaces(0);
    DynamicList<label> interfaceOwner(0);

    label defaultSize = 0;
    DynamicList<face>  defaultFaces(0);
    DynamicList<label> defaultOwner(0);

    //fields
    cellLevelField.reserve(nCells);
    cellIDField.reserve(nCells);
    procIDField.reserve(nCells);
    distanceField.reserve(nCells);

    //std::vector<scalar> triangleIDField;
    //triangleIDField.reserve(nCells);

    openvdb::Int32 nPoints = 0;

    label voxelCount = 0;

    std::vector<BoolGrid::Ptr> otherSubsets(maxCellLevel_ + 1, nullptr);
    if (subsetInnerGrid)
    {
        otherSubsets = outerGrids;
    }
    else if (subsetOuterGrid)
    {
        otherSubsets = innerGrids;
    }
    else
    {
        for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
        {
            otherSubsets[cellLevel] = BoolGrid::create(false);
        }
    }

    { Timer timer("face addressing");

    // wind tunnel bounding box in global coord (2 pointsPerEdge)
    openvdb::Int32 xmin = (2* boundingBox.min().x())    << maxCellLevel_;
    openvdb::Int32 xmax = (2*(boundingBox.max().x()+1)) << maxCellLevel_;
    openvdb::Int32 ymin = (2* boundingBox.min().y())    << maxCellLevel_;
    openvdb::Int32 ymax = (2*(boundingBox.max().y()+1)) << maxCellLevel_;
    openvdb::Int32 zmin = (2* boundingBox.min().z())    << maxCellLevel_;
    openvdb::Int32 zmax = (2*(boundingBox.max().z()+1)) << maxCellLevel_;

    // loop through voxels to create face addressing
    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        Timer timer("    cellLevel " + std::to_string(cellLevel));

        IndexGrid::Accessor cellIDAcc     = cellIDGrids[cellLevel]->getAccessor();
        IndexGrid::Accessor pointIDAcc    = pointIDGrids[cellLevel]->getAccessor();
        IndexGrid::Accessor pointLevelAcc = pointLevelGrids[cellLevel]->getAccessor();

        Vec3dGrid::ConstAccessor displAcc = displacementGrids[cellLevel]->getConstAccessor();

        label fineLevel = Foam::min(cellLevel + 1, maxCellLevel_);
        IndexGrid::Accessor fineCellIDAcc  =  cellIDGrids[fineLevel]->getAccessor();
        IndexGrid::Accessor finePointIDAcc = pointIDGrids[fineLevel]->getAccessor();

        BoolGrid::ConstAccessor otherSubsetFineAcc = otherSubsets[fineLevel]->getConstAccessor();
        BoolGrid::ConstAccessor innerGridAcc = innerGrids[cellLevel]->getConstAccessor();

        Vec3dGrid::ConstAccessor fineDisplAcc = displacementGrids[fineLevel]->getConstAccessor();

        FloatGrid::Ptr levelNGrid = cellLevelGrids[cellLevel];
        levelNGrid->tree().voxelizeActiveTiles();

        for (auto iter = levelNGrid->cbeginValueOn(); iter; ++iter)
        {
            const openvdb::Coord& ijk = iter.getCoord();

            label patchID = 0;

            const std::vector<openvdb::Coord> voxelVertices =
                getVoxelVertices(/*pointsPerEdge*/2, ijk);

            const openvdb::Int32& ownID = cellIDAcc.getValue(ijk);

            cellIDField.emplace_back(ownID);

            distanceField.emplace_back(iter.getValue());

            if (!subsetInnerGrid && !subsetOuterGrid)
            {
                if (innerGridAcc.isValueOn(ijk))
                {
                    innerSubset.emplace_back(ownID);
                }
                else
                {
                    outerSubset.emplace_back(ownID);
                }
            }

            openvdb::Int32 neighID;
            openvdb::Int32 pointID;

            //save edges/splitEdges of this voxel
            List<labelList> voxelEdges(12);
            std::vector<std::vector<openvdb::Coord>> voxelEdgesCoord(12, std::vector<openvdb::Coord>());


            boolList splitEdges(12, false);

            // if face get split save pointID of face centre in order
            // to create split face in the next face loop
            labelList faceCentreIDs(6, -1);
            std::vector<openvdb::Coord> faceCentreCoords(6);

            // first loop through voxel faces to determine edges that need to be split
            // and add points
            forAll(localFaces, i)
            {
                openvdb::Coord neighCoord = ijk + COORD_OFFSETS[i];

                cellIDAcc.probeValue(neighCoord, neighID);

                bool doSplitFace  = false;

                //Info<<"    face " <<faceNames[i]<<endl;

                openvdb::Coord fineNeighCoord;
                openvdb::Int32 fineNeighID;

                // first check if face needs to be split in 4
                if
                (
                    neighID == -labelMax // not internal or processor face at this cellLevel
                 && cellLevel < maxCellLevel_
                )
                {
                    // loop across subfaces
                    for (label j = 0; j < 4; j++)
                    {
                        fineNeighCoord = (ijk << 1) + COORD_OFFSETS_SPLIT[4*i + j];

                        fineCellIDAcc.probeValue(fineNeighCoord, fineNeighID);

                        if
                        (
                            fineNeighID > -labelMax
                         || otherSubsetFineAcc.isValueOn(fineNeighCoord)
                        )
                        {
                            doSplitFace = true;
                            //Info<<"    face " <<faceNames[i]
                            //    <<" splitEdges " <<edgeLoops[i]<<endl;

                            for (label eI = 0; eI < 4; eI++)
                            {
                                splitEdges[edgeLoops[i][eI]] = true;
                            }

                            // calculate face centre
                            // sum of vertex coordinates divided by 4 in index space
                            openvdb::Coord fc =
                            (
                                voxelVertices[localFaces[i][0]] +
                                voxelVertices[localFaces[i][1]] +
                                voxelVertices[localFaces[i][2]] +
                                voxelVertices[localFaces[i][3]]
                            ) >> 2;

                            faceCentreIDs[i] = nPoints;
                            faceCentreCoords[i] = fc;

                            //Info<<"           faceCentreID "<<nPoints<<endl;

                            pointIDAcc.setValue(fc, nPoints);
                            pointLevelAcc.setValue(fc, cellLevel+1);
                            finePointIDAcc.setValue(fc<<1, nPoints);

                            //transform to global index space
                            fc <<= (maxCellLevel_ - cellLevel);

                            xyzPoints.emplace_back(fc.x(), fc.y(), fc.z());
                            vdbPoints.emplace_back(point(fc.x(), fc.y(), fc.z()), nPoints);
                            ++nPoints;
                            break;
                        }
                    } //for fine neighbour
                } //check if face needs to be split in 4

                // second check if any internal or boundary face
                // needs to split any edge
                // don't check edges of z-min (i==4) and z-max (i==5)
                if (!doSplitFace && cellLevel < maxCellLevel_ && i < 4)
                {
                    openvdb::Coord edgeNeighCoord;
                    openvdb::Int32 edgeNeighID;

                    for (label eI = 0; eI < 4; eI++)
                    {
                        // two edges only for y-min and y-max
                        if (i > 1 && eI%2 == 1) continue;

                        const label edgeID = edgeLoops[i][eI];

                        edgeNeighCoord = ijk + EDGE_COORD_OFFSETS[edgeID];

                        cellIDAcc.probeValue(edgeNeighCoord, edgeNeighID);

                        if (edgeNeighID == -labelMax) // no edge adjacent neighbour of same cellLevel
                        {
                            // check fine voxels adjacent to edge
                            //for (label k = 0; k < 8; k++)
                            for (label k = 0; k < 2; k++)
                            {
                                //// skip fine voxels out of plane
                                //if
                                //(
                                //    k == localEdges[edgeID].start()
                                // || k == localEdges[edgeID].end()
                                //) continue;

                                //fineNeighCoord = (edgeNeighCoord << 1) + COARSE_TO_FINE[k];

                                ////GGG
                                fineNeighCoord = (ijk << 1) + FINE_EDGE_COORD_OFFSETS[2*edgeID + k];

                                fineCellIDAcc.probeValue(fineNeighCoord, fineNeighID);

                                if
                                (
                                    fineNeighID > -labelMax
                                 || otherSubsetFineAcc.isValueOn(fineNeighCoord)
                                )
                                {
                                    //Info<<"    face " <<faceNames[i]
                                    //    <<" splitEdge " <<edgeID<<endl;;
                                    splitEdges[edgeID] = true;
                                    break;
                                }
                            } // check fine voxels adjacent to edge

                            //check fine voxels adjacent to active neighbour voxel
                            if (splitEdges[edgeID] == false)
                            //if
                            //(
                            //    ((i==1 && edgeID==4))// || (i==3 && edgeID==10))
                            // && splitEdges[edgeID] == false
                            //)
                            {
                                openvdb::Coord eijk0 = ijk + COORD_OFFSETS[edgeFace[edgeID][0]];
                                openvdb::Coord eijk1 = ijk + COORD_OFFSETS[edgeFace[edgeID][1]];

                                openvdb::Int32 eID0, eID1;
                                cellIDAcc.probeValue(eijk0, eID0);
                                cellIDAcc.probeValue(eijk1, eID1);

                                bool checkFineVoxels = false;

                                label id = 0;

                                if (eID0 == -labelMax && eID1 > -labelMax)
                                {
                                    checkFineVoxels = true;
                                    id = 1;
                                }
                                else if (eID1 == -labelMax && eID0 > -labelMax)
                                {
                                    checkFineVoxels = true;
                                    id = 0;
                                }


                                //for (label k = 0; k < 2; k++)
                                if (checkFineVoxels)
                                {
                                    openvdb::Coord edgeFaceCoord, fineEdgeFaceCoord;

                                    //edgeFaceCoord = ijk + COORD_OFFSETS[edgeFace[edgeID][k]];
                                    edgeFaceCoord = ijk + COORD_OFFSETS[edgeFace[edgeID][id]];

                                    cellIDAcc.probeValue(edgeFaceCoord, edgeNeighID);

                                    if (edgeNeighID > -labelMax)
                                    {
                                        for (label z = 0; z < 2; z++)
                                        {
                                            fineEdgeFaceCoord = (edgeFaceCoord << 1) + COORD_OFFSETS_SPLIT[fineEdgeFace[2*edgeID+id][z]];

                                            fineCellIDAcc.probeValue(fineEdgeFaceCoord, fineNeighID);

                                            if
                                            (
                                                fineNeighID > -labelMax
                                             || otherSubsetFineAcc.isValueOn(fineEdgeFaceCoord)
                                            )
                                            {
                                                //Info<<"    face " <<faceNames[i]
                                                //    <<" splitEdge " <<edgeID<<endl;;
                                                splitEdges[edgeID] = true;
                                                break;
                                            }
                                        }
                                    }
                                }
                            } //check fine voxels adjacent to active neighbour voxel
                        } // if there is no neighbour at same cellLevel
                    } //for face edge
                } //check edge-voxels adjacency

                //add points and split edges
                for (label eI = 0; eI < 4; eI++)
                {
                    const label edgeID = edgeLoops[i][eI];
                    const openvdb::Coord& edgeStart = voxelVertices[localEdges[edgeID].start()];
                    const openvdb::Coord& edgeEnd   = voxelVertices[localEdges[edgeID].end()];

                    //face x-min and x-max add all voxel points and the first 8 edges
                    if (i < 2)
                    {
                        // add edge START point
                        if (eI == 0)
                        {
                            addPoint
                            (
                                edgeStart,
                                pointIDAcc,
                                pointLevelAcc,
                                voxelEdgesCoord[edgeID],
                                voxelEdges[edgeID],
                                nPoints,
                                cellLevel,
                                xyzPoints,
                                vdbPoints
                            );
                        }
                        else
                        {
                            voxelEdges[edgeID].append
                            (
                                voxelEdges[edgeID - 1].last()
                            );

                            voxelEdgesCoord[edgeID].emplace_back
                            (
                                voxelEdgesCoord[edgeID - 1].back()
                            );
                        }

                        // add edge END point
                        if (eI == 3)
                        {
                            voxelEdges[edgeID].append
                            (
                                voxelEdges[edgeID - 3].first()
                            );

                            voxelEdgesCoord[edgeID].emplace_back
                            (
                                voxelEdgesCoord[edgeID - 3].front()
                            );
                        }
                        else
                        {
                            addPoint
                            (
                                edgeEnd,
                                pointIDAcc,
                                pointLevelAcc,
                                voxelEdgesCoord[edgeID],
                                voxelEdges[edgeID],
                                nPoints,
                                cellLevel,
                                xyzPoints,
                                vdbPoints
                            );
                        }
                    } // x-min and x-max
                    else if (i < 4 && eI%2 == 0)
                    {
                        // y-min(i==2) and y-max(i==3) add only 2 unique edges each
                        // e.g. y-min = edgeLoops[2] adds only edge 8 (i.e. edgeLoops[2][0])
                        //      and edge 9 (i.e. edgeLoops[2][2]

                        // add edge START point
                        voxelEdges[edgeID].append
                        (
                            voxelEdges[edgeEdges[edgeID].first()].first()
                        );
                        voxelEdgesCoord[edgeID].emplace_back
                        (
                            voxelEdgesCoord[edgeEdges[edgeID].first()].front()
                        );

                        // add edge END point
                        voxelEdges[edgeID].append
                        (
                            voxelEdges[edgeEdges[edgeID].last()].first()
                        );
                        voxelEdgesCoord[edgeID].emplace_back
                        (
                            voxelEdgesCoord[edgeEdges[edgeID].last()].front()
                        );
                    } // y-min and y-max


                    // add mid point if edge needs to be split
                    // this will also update edges added previously for example by x-min or x-max
                    if (splitEdges[edgeID] && voxelEdges[edgeID].size() == 2)
                    {
                        const openvdb::Coord edgeMid =
                        (
                            edgeStart
                          + edgeEnd
                        ) >> 1;

                        addMidPoint
                        (
                            edgeMid,
                            pointIDAcc,
                            pointLevelAcc,
                            voxelEdgesCoord[edgeID],
                            voxelEdges[edgeID],
                            nPoints,
                            cellLevel,
                            xyzPoints,
                            vdbPoints,
                            pointID
                        );

                        // update pointIDGrid of finer cellLevel
                        finePointIDAcc.setValue(edgeMid<<1, pointID);

                        pointID = pointIDAcc.getValue(edgeStart);
                        finePointIDAcc.setValue(edgeStart<<1, pointID);

                        pointID = pointIDAcc.getValue(edgeEnd);
                        finePointIDAcc.setValue(edgeEnd<<1, pointID);
                    }
                } // for face edge eI
            }// first loop through voxel faces to determine edges that need to be split and add points

            // add faces
            forAll(localFaces, i)
            {
                openvdb::Coord neighCoord = ijk + COORD_OFFSETS[i];

                // check neighbour of same cellLevel
                // if neighbour is active the face is either internal or is a processor face
                cellIDAcc.probeValue(neighCoord, neighID);

                // === 0 === internal face, already added
                if (neighID >= 0 && neighID < ownID)
                {
                    /*if (debug)*/
                    //Info<< "    face " << faceNames[i]
                    //    << " - own " << ownID
                    //    << " - neigh " << neighID
                    //    << " already added"
                    //    << endl;

                    continue;
                }

                face f;
                std::vector<openvdb::Coord> faceCoords;

                for (label eI = 0; eI < 4; eI++)
                {
                    const label& edgeID = edgeLoops[i][eI];

                    if ((i < 2) || (i < 4 && eI%2 == 0))
                    {
                        // x-min, x-max, y-min[0,2], y-max[0,2]
                        f.append(voxelEdges[edgeID].first());
                        faceCoords.emplace_back
                        (
                            voxelEdgesCoord[edgeID].front()
                        );
                    }
                    else
                    {
                        // z-min, z-max, y-min[1,3], y-max[1,3]
                        f.append(voxelEdges[edgeID].last());
                        faceCoords.emplace_back
                        (
                            voxelEdgesCoord[edgeID].back()
                        );
                    }


                    if (voxelEdges[edgeID].size() == 3)
                    {
                        f.append(voxelEdges[edgeID][1]);
                        faceCoords.emplace_back
                        (
                            voxelEdgesCoord[edgeID][1]
                        );
                    }
                }

                // === 1 === internal face
                if (neighID >= 0)
                {
                    /*if (debug)*/
                    //Info<< "    face " << faceNames[i]
                    //    << " - own " << ownID
                    //    << " - neigh " << neighID
                    //    << " - new face " << f
                    //    << " - internal"
                    //    << endl;

                    owner.append(ownID);
                    neighbour.append(neighID);

                    faces.append(f);

                    vdbInternalFaces.emplace_back
                    (
                        //faceCoords, ownID, neighID, cellLevel, nInternal++
                        f, ownID, neighID, cellLevel, nInternal++
                    );

                    continue;
                }

                // calculate face centre
                // (for sorting later if procBoundary
                // or assign to windTunnel patch if boundary)
                // sum of vertex coordinates divided by 4 in index space
                openvdb::Coord fc =
                (
                    voxelVertices[localFaces[i][0]] +
                    voxelVertices[localFaces[i][1]] +
                    voxelVertices[localFaces[i][2]] +
                    voxelVertices[localFaces[i][3]]
                ) >> 2;

                //transform to global index space
                fc <<= (maxCellLevel_ - cellLevel);

                // === 2 === processorBoundary
                if (neighID >= -(Pstream::nProcs() + 1))
                {
                    const label procID = -(neighID + 1);

                    procOwners[procID].append(ownID);

                    // match processor face anchor point f[0] of pY
                    // 2376 -> 3762
                    if (i == 3) //plusY
                    {
                        label pos = 1;
                        if (voxelEdges[10].size() == 3)
                        {
                            pos = 2;
                        }

                        inplaceRotateList<List, label>(f, -pos);

                        //rotate to the left
                        std::rotate
                        (
                            faceCoords.begin(),
                            faceCoords.begin() + pos,
                            faceCoords.end()
                        );
                    }

                    procsFaces[procID].append(f);

                    vdbProcFaces[procID].emplace_back
                    (
                        //faceCoords, ownID, /*neighID*/-1, cellLevel, nProcFace[procID]++
                        f, ownID, /*neighID*/-1, cellLevel, nProcFace[procID]++
                    );

                    /*if (debug)*/
                    //Pout<< "    face " << faceNames[i]
                    //    << " - own " << ownID
                    //    << " - neigh " << neighID
                    //    << " - toProc" << procID
                    //    << " - new face " << f
                    //    << endl;

                    procsFaceCentres[procID].append
                    (
                        point(fc.x(), fc.y(), fc.z())
                    );

                    continue;
                }

                // == 3 == cellLevel interface
                if (faceCentreIDs[i] > -1)
                {
                    const label& fcID = faceCentreIDs[i]; //face centre ID
                    const openvdb::Coord& fcCoord = faceCentreCoords[i]; //face centre ID

                    /*if (debug)*/
                    //Info<< "    face " << faceNames[i]
                    //    << " - own " << ownID
                    //    << " - neigh " << neighID
                    //    << " - new face " << f
                    //    << " - split face - (coarse faceID " << fcID
                    //    << ")" << endl;

                    for (label j = 0; j < 4; j++)
                    {
                        label prev = (j + 3)%4;

                        //Info<< " voxelEdges[edgeLoops[i]["<<j<<"]] " << voxelEdges[edgeLoops[i][j]] <<endl;

                        face subf(4);
                        std::vector<openvdb::Coord> subFaceCoords(4);

                        label pointID;
                        if ((i < 2) || (i < 4 && j%2 == 0))
                        {
                            subf[0] = voxelEdges[edgeLoops[i][j]].first();
                            subFaceCoords[0] = voxelEdgesCoord[edgeLoops[i][j]].front();
                        }
                        else
                        {
                            subf[0] = voxelEdges[edgeLoops[i][j]].last();
                            subFaceCoords[0] = voxelEdgesCoord[edgeLoops[i][j]].back();
                        }
                        subf[1] = voxelEdges[edgeLoops[i][j]][1];    //mid point of current edge
                        subf[2] = fcID;
                        subf[3] = voxelEdges[edgeLoops[i][prev]][1]; // mid point of previous edge

                        subFaceCoords[1] = voxelEdgesCoord[edgeLoops[i][j]][1];    //mid point of current edge
                        subFaceCoords[2] = fcCoord;
                        subFaceCoords[3] = voxelEdgesCoord[edgeLoops[i][prev]][1]; // mid point of previous edge

                        ///*if (debug)*/ Info<<"      subFace["<<j<<"]: "<<subf<<endl;

                        // check if there is a neighbour at a finer level
                        neighCoord = (ijk << 1) + COORD_OFFSETS_SPLIT[4*i + j];

                        fineCellIDAcc.probeValue(neighCoord, neighID);

                        // === 3a === internal face (at cellLevel interface)
                        if (neighID >= 0)
                        {
                            owner.append(ownID);
                            neighbour.append(neighID);

                            faces.append(subf);

                            vdbInternalFaces.emplace_back
                            (
                                //subFaceCoords, ownID, neighID, cellLevel, nInternal++
                                subf, ownID, neighID, cellLevel, nInternal++
                            );

                            // update fine grid with owner of this subface
                            fineCellIDAcc.setValue((neighCoord - COORD_OFFSETS[i]), ownID);

                            continue;
                        }

                        // find split face centre for sorting later
                        // or assign to wind tunnel boundaries
                        point subFcentre =
                        (
                            xyzPoints[subf[0]] +
                            xyzPoints[subf[1]] +
                            xyzPoints[subf[2]] +
                            xyzPoints[subf[3]]
                        ) / 4.0;

                        // === 3b === processorBoundary (at cellLevel interface)
                        if (neighID >= -(Pstream::nProcs() + 1))
                        {
                            const label procID = -(neighID + 1);

                            procOwners[procID].append(ownID);

                            // match processor face anchor point f[0]
                            if (i == 3) //plusY
                            {
                                // 2376 -> 3762
                                inplaceRotateList<List, label>(subf, -(1-j));

                                if (j < 2)
                                {
                                    //rotate to the left
                                    std::rotate
                                    (
                                        subFaceCoords.begin(),
                                        subFaceCoords.begin() + (1-j),
                                        subFaceCoords.end()
                                    );
                                }
                                else
                                {
                                    //rotate to the right
                                    std::rotate
                                    (
                                        subFaceCoords.rbegin(),
                                        subFaceCoords.rbegin() + (j-1),
                                        subFaceCoords.rend()
                                    );
                                }
                            }
                            else
                            {
                                inplaceRotateList<List, label>(subf, j);

                                //rotate to the right
                                std::rotate
                                (
                                    subFaceCoords.rbegin(),
                                    subFaceCoords.rbegin() + j,
                                    subFaceCoords.rend()
                                );
                            }

                            procsFaces[procID].append(subf);

                            vdbProcFaces[procID].emplace_back
                            (
                                //subFaceCoords, ownID, /*neighID*/-1, cellLevel, nProcFace[procID]++
                                subf, ownID, /*neighID*/-1, cellLevel, nProcFace[procID]++
                            );

                            procsFaceCentres[procID].append(subFcentre);

                            continue;
                        }

                        bool isValidFace = true;

                        // === 3c === boundary (split face)
                        label domainPatchID = 0;
                        if      (subFcentre.x() == xmin) domainPatchID = 0;
                        else if (subFcentre.x() == xmax) domainPatchID = 1;
                        else if (subFcentre.y() == ymin) domainPatchID = 2;
                        else if (subFcentre.y() == ymax) domainPatchID = 3;
                        else if (subFcentre.z() == zmin) domainPatchID = 4;
                        else if (subFcentre.z() == zmax) domainPatchID = 5;

                        if
                        (
                            subFcentre.x() == xmin || subFcentre.x() == xmax
                         || subFcentre.y() == ymin || subFcentre.y() == ymax
                         || subFcentre.z() == zmin || subFcentre.z() == zmax
                        )
                        {
                            patchFaces[domainPatchID].append(subf);
                            patchOwners[domainPatchID].append(ownID);
                            boundaryMeshSizes[domainPatchID]++;
                        }
                        else
                        {
                            defaultFaces.append(subf);
                            defaultOwner.append(ownID);
                            defaultSize++;
                        }
                        //else
                        //{
                        //    //add displacement to boundary points
                        //    isValidFace = addDisplacement
                        //        (
                        //            cellLevel,
                        //            subf,
                        //            subFaceCoords,
                        //            displacedPoints,
                        //            xyzPoints,
                        //            vdbPoints,
                        //            displAcc,
                        //            fineDisplAcc,
                        //            true //isSplitFace
                        //        );
                        ////GGG
                        //    interfaceFaces.append(subf);
                        //    interfaceOwner.append(ownID);
                        //    interfaceSize++;
                        ////GGG
                        //}

                        ////TODO get face region
                        //if (isValidFace)
                        //{
                        //    vdbBoundaryFaces[0].emplace_back
                        //    (
                        //        //subFaceCoords, ownID, /*neighID*/-1, cellLevel, nBoundaryFace++
                        //        subf, ownID, /*neighID*/-1, cellLevel, nBoundaryFace++
                        //    );
                        //}
                    } // for j<4 - add 4 split faces

                    continue;
                } // if splitFace

                /*if (debug)*/
                //Info<< "    face " << faceNames[i]
                //    << " - own " << ownID
                //    << " - neigh " << neighID
                //    << " - new face " << f;

                bool isValidFace = true;

                // === 4 === boundary face
                label domainPatchID = 0;
                if      (fc.x() == xmin) domainPatchID = 0;
                else if (fc.x() == xmax) domainPatchID = 1;
                else if (fc.y() == ymin) domainPatchID = 2;
                else if (fc.y() == ymax) domainPatchID = 3;
                else if (fc.z() == zmin) domainPatchID = 4;
                else if (fc.z() == zmax) domainPatchID = 5;

                if
                (
                    fc.x() == xmin || fc.x() == xmax
                 || fc.y() == ymin || fc.y() == ymax
                 || fc.z() == zmin || fc.z() == zmax
                )
                {
                    patchFaces[domainPatchID].append(f);
                    patchOwners[domainPatchID].append(ownID);
                    boundaryMeshSizes[domainPatchID]++;
                }
                else
                {
                    defaultFaces.append(f);
                    defaultOwner.append(ownID);
                    defaultSize++;
                }
                //else
                //{
                //    //add displacement to boundary points
                //    isValidFace = addDisplacement
                //        (
                //            cellLevel,
                //            f,
                //            faceCoords,
                //            displacedPoints,
                //            xyzPoints,
                //            vdbPoints,
                //            displAcc,
                //            fineDisplAcc
                //        );
                ////
                ////GGG
                //    interfaceFaces.append(f);
                //    interfaceOwner.append(ownID);
                //    interfaceSize++;
                ////GGG
                //}

                //if (isValidFace)
                //{
                //    vdbBoundaryFaces[0].emplace_back
                //    (
                //        //faceCoords, ownID, /*neighID*/-1, cellLevel, nBoundaryFace++
                //        f, ownID, /*neighID*/-1, cellLevel, nBoundaryFace++
                //    );
                //}
            } // forAll localFaces i

            voxelCount++;
        } // for cellID

        //write fields
        cellLevelField.resize(voxelCount, cellLevel);

        if (Pstream::parRun())
        {
            procIDField.resize(voxelCount, Pstream::myProcNo());
        }

        levelNGrid->clear();
        cellIDGrids[cellLevel]->clear();
        //pointIDGrids[cellLevel]->clear();
    } // for cellLevel
    } //Timer face addressing

    //GGG
    nPoints = activeCount<vdbPoint>(vdbPoints);

    //do not add collapsed faces
    for (size_t i = 0; i < vdbInternalFaces.size(); ++i)
    {
        label nActive = 0;
        for (size_t pointI = 0; pointI < vdbInternalFaces[i].vdbPoints_.size(); ++pointI)
        {
            label vdbPointID = vdbInternalFaces[i].vdbPoints_[pointI];

            if (vdbPoints[vdbPointID].active_)
            {
                ++nActive;
            }
        }
        if (nActive < 3)
        {
            vdbInternalFaces[i].active_= false;
        }
    }
    nInternal = activeCount<vdbFace>(vdbInternalFaces);

    for (size_t i = 0; i < vdbBoundaryFaces[0].size(); ++i)
    {
        label nActive = 0;
        for (size_t pointI = 0; pointI < vdbBoundaryFaces[0][i].vdbPoints_.size(); ++pointI)
        {
            label vdbPointID = vdbBoundaryFaces[0][i].vdbPoints_[pointI];

            if (vdbPoints[vdbPointID].active_)
            {
                ++nActive;
            }
        }
        if (nActive < 3)
        {
            vdbBoundaryFaces[0][i].active_= false;
        }
    }
    nBoundaryFace = activeCount<vdbFace>(vdbBoundaryFaces[0]);

    for (size_t i = 0; i < vdbProcFaces.size(); ++i)
    {
        nProcFace[i] = activeCount<vdbFace>(vdbProcFaces[i]);
    }

    label nFaces = nInternal + nBoundaryFace + sum(nProcFace);
    faceList       foamFaces(nFaces);
    labelList      foamOwner(nFaces);
    labelList      foamNeighbour(nInternal);

    //{
    //    Timer t("vdbFaces to OpenFOAM faces");

    //    //TODO multithreaded
    //    for (size_t i = 0; i < vdbInternalFaces.size(); ++i)
    //    {
    //        if (!vdbInternalFaces[i].active_) continue;

    //        const label& cellLevel = vdbInternalFaces[i].cellLevel_;
    //        const label& faceID    = vdbInternalFaces[i].id_;

    //        face f;
    //        for (size_t pointI = 0; pointI < vdbInternalFaces[i].vdbPoints_.size(); ++pointI)
    //        {
    //            label vdbPointID = vdbInternalFaces[i].vdbPoints_[pointI];

    //            if (vdbPoints[vdbPointID].active_)
    //            {
    //                label pointID = vdbPoints[vdbPointID].id_;

    //                f.append(pointID);
    //            }
    //        }

    //        foamFaces[faceID] = f;
    //        foamOwner[faceID] = vdbInternalFaces[i].owner_;
    //        foamNeighbour[faceID] = vdbInternalFaces[i].neighbour_;
    //    }

    //    //TODO multithreaded
    //    for (size_t i = 0; i < vdbBoundaryFaces[0].size(); ++i)    //TODO for each patch
    //    {
    //        if (!vdbBoundaryFaces[0][i].active_) continue;

    //        const label& cellLevel = vdbBoundaryFaces[0][i].cellLevel_;
    //        const label faceID    = nInternal + vdbBoundaryFaces[0][i].id_;

    //        face f;
    //        for (size_t pointI = 0; pointI < vdbBoundaryFaces[0][i].vdbPoints_.size(); ++pointI)
    //        {
    //            label vdbPointID = vdbBoundaryFaces[0][i].vdbPoints_[pointI];

    //            if (vdbPoints[vdbPointID].active_)
    //            {
    //                label pointID = vdbPoints[vdbPointID].id_;

    //                f.append(pointID);
    //            }
    //        }

    //        foamFaces[faceID]     = f;
    //        foamOwner[faceID]     = vdbBoundaryFaces[0][i].owner_;
    //        //foamNeighbour[faceID] = vdbBoundaryFaces[0][i].neighbour_;
    //    }

    //    for (size_t proci = 0; proci < Pstream::nProcs(); ++proci)
    //    {
    //        labelList sortedIndices;
    //        sortedOrder
    //        (
    //            procsFaceCentres[proci],
    //            sortedIndices,
    //            pointLess(procsFaceCentres[proci])
    //        );
    //    //TODO multithreaded
    //    for (size_t i = 0; i < vdbProcFaces[proci].size(); ++i)
    //    {
    //        if (!vdbProcFaces[proci][i].active_) continue;

    //        const label& cellLevel = vdbProcFaces[proci][i].cellLevel_;
    //        const label faceID = nInternal + nBoundaryFace + sortedIndices[vdbProcFaces[proci][i].id_];

    //        face f;
    //        for (size_t pointI = 0; pointI < vdbProcFaces[proci][i].vdbPoints_.size(); ++pointI)
    //        {
    //            label vdbPointID = vdbProcFaces[proci][i].vdbPoints_[pointI];

    //            //if (vdbPoints[vdbPointID].active_)
    //            {
    //                label pointID = vdbPoints[vdbPointID].id_;

    //                f.append(pointID);
    //            }
    //        }

    //        foamFaces[faceID]     = f;
    //        foamOwner[faceID]     = vdbProcFaces[proci][i].owner_;
    //        //foamNeighbour[faceID] = vdbProcFaces[proci][i].neighbour_;
    //    }
    //    }

    //}
    //GGG

    //append patch and processor faces to internal faces
    //label patchStart = faces.size();
    label patchStart = nInternal;

    forAll(patchFaces, patchI)
    {
        faces.append(patchFaces[patchI]);
        owner.append(patchOwners[patchI]);

        boundaryMeshStarts[patchI] = patchStart;

        patchStart += boundaryMeshSizes[patchI];
    }

    label allInterfaceSize = interfaceSize;
    label allDefaultSize = defaultSize;

    reduce(
        std::tie(allInterfaceSize, allDefaultSize),
        ParallelOp<maxOp<label>, maxOp<label>>{}
    );

    if (allDefaultSize)
    {
        faces.append(defaultFaces);
        owner.append(defaultOwner);

        dictionary patchInfo;
        patchInfo.add("type", wallPolyPatch::typeName);
        patchInfo.add("nFaces", 0);
        patchInfo.add("startFace", 0);

        addEmptyPatch(mesh, "defaultFaces", patchInfo);

        boundaryMeshStarts.append(patchStart);
        boundaryMeshSizes.append(defaultSize);

        patchStart += defaultSize;
    }

    if (allInterfaceSize)
    {
        faces.append(interfaceFaces);
        owner.append(interfaceOwner);

        dictionary patchInfo;
        patchInfo.add("type", wallPolyPatch::typeName);
        patchInfo.add("nFaces", 0);
        patchInfo.add("startFace", 0);

        word interfaceName = "allBoundaries";
        if (subsetInnerGrid)
        {
            interfaceName = "innerToOuter";
        }
        else if (subsetOuterGrid)
        {
            interfaceName = "outerToInner";
        }
        addEmptyPatch(mesh, interfaceName, patchInfo);

        boundaryMeshStarts.append(patchStart);
        //boundaryMeshSizes.append(interfaceSize);
        boundaryMeshSizes.append(nBoundaryFace);

        //patchStart += interfaceSize;
        patchStart += nBoundaryFace;
    }

    { Timer timer("proc face sorting");
    //sort processor faces
    forAll(procBoundarySizes, proci)
    {
        procBoundarySizes[proci] += procsFaces[proci].size();

        if (procBoundarySizes[proci])
        {
            labelList sortedIndices;
            sortedOrder
            (
                procsFaceCentres[proci],
                sortedIndices,
                pointLess(procsFaceCentres[proci])
            );

            forAll(procsFaceCentres[proci], i)
            {
                //if (debug)
                //{
                //    Pout<< i << " "
                //        << "face "
                //        <<procsFaces[proci][sortedIndices[i]]
                //        << " fc "
                //        << cmptMultiply((((procsFaceCentres[proci][sortedIndices[i]])/2)*voxelSize_)-levelSetOffset, scale)
                //        <<endl;
                //}
                faces.append(procsFaces[proci][sortedIndices[i]]);
                owner.append(procOwners[proci][sortedIndices[i]]);
            }
        }
    }
    } // Timer proc face sorting


    // ### Construct polyMesh ### //

    label totalCells = nCells;

    Info<< "\nTotal number of cells: "
        << returnReduce(totalCells, sumOp<label>())
        << endl;

    if (Pstream::parRun())
    {
        labelList cellsPerProc(Pstream::nProcs());
        cellsPerProc[Pstream::myProcNo()] = nCells;
        IPstream::gatherList(cellsPerProc);

        Info<< "\tcells per processor: " << endl;
        Info<< "\t" << cellsPerProc << nl <<endl;
    }

    //pointField points(xyzPoints.size());
    points.setSize(nPoints);

    //transform points from index space to world space
    const auto grainSize =
        std::max<size_t>
        (
            vdbPoints.size() / tbb::this_task_arena::max_concurrency(),
            1024
        );

    const tbb::blocked_range<size_t> range(0, vdbPoints.size(), grainSize);

    tbb::parallel_for
    (
        range,
        TransformPoints(vdbPoints, points, levelSetOffset, scale, R, voxelSize_),
        tbb::simple_partitioner()
    );

    //points = cmptMultiply(((xyzPoints / /*pointsPerEdge*/2) * voxelSize_) - levelSetOffset, scale);

    //add processor boundaries
    if (Pstream::parRun())
    {
        const bool procAsWall(vdbDict.lookupOrDefault<bool>("procAsWall", false));

        forAll(procBoundarySizes, neighI)
        {
            if (procBoundarySizes[neighI] == 0 && !procAsWall) continue;

            dictionary patchInfo;

            word patchName = procBoundaryNames[neighI];

            if (procAsWall)
            {
                Info<< "Debug: writing processor boundaries as wallPolyPatches"
                    << endl;

                patchInfo.add("type", wallPolyPatch::typeName);
                patchName = "procBoundaryMeto" + Foam::name(neighI);
            }
            else
            {
                patchInfo.add("type", processorPolyPatch::typeName);
                patchInfo.add("myProcNo", Pstream::myProcNo());
                patchInfo.add("neighbProcNo", neighI);
            }

            Pout<< "neigh " << neighI
                << " start " << patchStart
                << " size " << procBoundarySizes[neighI]
                << " patchName " << patchName
                << endl;

            patchInfo.add("nFaces", 0);
            patchInfo.add("startFace", 0);

            addEmptyPatch
            (
                mesh,
                patchName,
                patchInfo
            );

            boundaryMeshStarts.append(patchStart);
            boundaryMeshSizes.append(procBoundarySizes[neighI]);

            patchStart += procBoundarySizes[neighI];
        }
    } //if parRun add processor boundaries

    //if(debug)
    //{
    //    Info<< "points " << points.size() <<endl;
    //    Info<< "faces " << faces <<endl;
    //    Info<< "owner " << owner<<endl;
    //    Info<< "neighbour " << neighbour<<endl;
    //    Info<< "boundaryMeshStarts " << boundaryMeshStarts<<endl;
    //    Info<< "boundaryMeshSizes " << boundaryMeshSizes <<endl;
    //    Info<< "mesh.boundaryMesh.names() " << mesh.boundaryMesh().names() <<endl;
    //    Info<< "mesh.boundaryMesh.types() " << mesh.boundaryMesh().types() <<endl;
    //}

    //TODO update with correct pointIDs!
    //if (false)
    {
    // update hexRef8.cellLevel
    hexRef8cellLevel.setSize(nCells);

    forAll(hexRef8cellLevel, celli)
    {
        hexRef8cellLevel[celli] = cellLevelField[celli];
    }

    // update hexRef8.pointLevel
    hexRef8pointLevel.setSize(points.size());

    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        IndexGrid::ConstAccessor pointIDAcc    = pointIDGrids[cellLevel]->getConstAccessor();
        IndexGrid::ConstAccessor pointLevelAcc = pointLevelGrids[cellLevel]->getConstAccessor();

        for (auto iter = pointIDGrids[cellLevel]->cbeginValueOn(); iter; ++iter)
        {
            hexRef8pointLevel[iter.getValue()] = pointLevelAcc.getValue(iter.getCoord());
        }
    }
    }

    { Timer t("mesh.resetPrimitives");
    mesh.resetPrimitives
    (
        xferMove(points),
        //foamFaces.xfer(),       // xferMove(faces),
        //foamOwner.xfer(),       // xferMove(owner),
        //foamNeighbour.xfer(),   // xferMove(neighbour),
        faces.xfer(),       // xferMove(faces),
        owner.xfer(),       // xferMove(owner),
        neighbour.xfer(),   // xferMove(neighbour),
        boundaryMeshSizes,  // patchSizes
        boundaryMeshStarts, // patchStarts
        true                // validBoundary (do parallel sync)
    );
    }

    if (!subsetInnerGrid && !subsetOuterGrid && balanceSurfaceCells)
    {
        Timer t("inner/outer cellZones");
        cellZoneMesh& cellZones = mesh.cellZones();

        labelList innerCells(innerSubset.size());
        labelList outerCells(outerSubset.size());

        forAll(innerCells, i)
        {
            innerCells[i] = innerSubset[i];
        }
        forAll(outerCells, i)
        {
            outerCells[i] = outerSubset[i];
        }

        //create inner zone
        cellZones.setSize(2);
        cellZones.set
        (
            0,
            new cellZone
            (
                "innerGrid", //name
                innerCells,  //addressing
                0,           //index
                cellZones    //cellZoneMesh
            )
        );
        cellZones.set
        (
            1,
            new cellZone
            (
                "outerGrid", //name
                outerCells,  //addressing
                1,           //index
                cellZones    //cellZoneMesh
            )
        );
    }

    //add volFields
    createVolField<scalar>(mesh, cellLevelField, "cellLevelVDB");
    //if (debug)
    {
        //createVolField<scalar>(mesh, cellIDField,     "cellID");
        //createVolField<scalar>(mesh, distanceField,   "distance");
        //createVolField<scalar>(mesh, triangleIDField, "triangleID");

        if (Pstream::parRun())
        {
            createVolField<scalar>(mesh, procIDField, "procID");
        }
    }
} // vdbGridsToPolyMesh
