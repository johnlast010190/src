#ifndef GEOMETRY_UTIL
#define GEOMETRY_UTIL

#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h> // for openvdb::tools::MeshToVoxelEdgeData
#include <openvdb/tools/VolumeToMesh.h> // for findFeaturePoint
#include <openvdb/tree/LeafManager.h>
#include <openvdb/util/Util.h> // for openvdb::util::COORD_OFFSETS

label findRegion
(
    point faceCentre,
    vector faceNormal,
    const openvdb::math::Transform& transform,
    openvdb::tree::ValueAccessor<const IndexGrid::TreeType>& idxAcc,
    const triSurface& refSurf
)
{
    openvdb::Vec3d pos;
    openvdb::Coord ijk, nijk;
    label idx, tmpIdx, region = 0;

    pos[0] = faceCentre.x();
    pos[1] = faceCentre.y();
    pos[2] = faceCentre.z();

    pos = transform.worldToIndex(pos);

    ijk[0] = int(std::floor(pos[0]));
    ijk[1] = int(std::floor(pos[1]));
    ijk[2] = int(std::floor(pos[2]));

    idx = idxAcc.getValue(ijk);

    vector refNormal = refSurf.faceNormals()[idx];

    if (abs(faceNormal & refNormal) > 0.95)
    {
        region = refSurf[idx].region();
    }
    else
    {
        //for (size_t i = 0; i < 18; ++i)
        for (size_t i = 0; i < 26; ++i)
        {
            nijk = ijk + openvdb::util::COORD_OFFSETS[i];

            if (idxAcc.probeValue(nijk, tmpIdx) && tmpIdx != idx)
            {
                vector tmpRefNormal = refSurf.faceNormals()[tmpIdx];

                if (abs(faceNormal & tmpRefNormal) > 0.95)
                {
                    region = refSurf[tmpIdx].region();
                    break;
                }
            }
        }
    }

    return region;
}

class PrimCpy
{
public:
    PrimCpy(const openvdb::tools::PolygonPoolList& primsIn, const std::vector<size_t>& indexList,
        List<labelledTri>& primsOut)
        : mPrimsIn(primsIn)
        , mIndexList(indexList)
        , mPrimsOut(&primsOut)
        {
        }

    void runParallel()
    {
        tbb::parallel_for(tbb::blocked_range<size_t>(0, mIndexList.size()), *this);
    }

    void runSerial()
    {
        (*this)(tbb::blocked_range<size_t>(0, mIndexList.size()));
    }

    inline void operator()(const tbb::blocked_range<size_t>& range) const
    {
        labelledTri tri(0,0,0,0);
        List<labelledTri>& primsOut = *mPrimsOut;

        for (size_t n = range.begin(); n < range.end(); ++n) {
            size_t index = mIndexList[n];
            openvdb::tools::PolygonPool& polygons = mPrimsIn[n];

            // Copy quads
            for (size_t i = 0, I = polygons.numQuads(); i < I; ++i) {
                const openvdb::Vec4I& quad = polygons.quad(i);

                // split quad in 2 triangles
                tri[0] = label(quad[0]);
                tri[1] = label(quad[1]);
                tri[2] = label(quad[2]);
                primsOut[index++] = tri;

                tri[0] = label(quad[0]);
                tri[1] = label(quad[2]);
                tri[2] = label(quad[3]);
                primsOut[index++] = tri;
            }
            polygons.clearQuads();

            // Copy triangles (adaptive mesh)
            for (size_t i = 0, I = polygons.numTriangles(); i < I; ++i) {
                const openvdb::Vec3I& triangle = polygons.triangle(i);
                tri[0] = label(triangle[0]);
                tri[1] = label(triangle[1]);
                tri[2] = label(triangle[2]);
                primsOut[index++] = tri;
            }
            polygons.clearTriangles();
        }
    }

private:
    const openvdb::tools::PolygonPoolList& mPrimsIn;
    const std::vector<size_t>& mIndexList;
    List<labelledTri> * const mPrimsOut;
}; //PrimCpy


// TBB body object for threaded copy of points from VDB to OpenFOAM format
class PointListCopy
{
    using PointList = std::unique_ptr<openvdb::Vec3s[]>;

    public:
        PointListCopy
        (
            const PointList& pointsIn,
            pointField& pointsOut
        )
        :
            mPointsIn(pointsIn),
            mPointsOut(pointsOut)
        {
        }

        void operator()(const tbb::blocked_range<size_t>& range) const
        {
            for (size_t n = range.begin(); n < range.end(); ++n)
            {
                point p
                (
                    mPointsIn[n].x(),
                    mPointsIn[n].y(),
                    mPointsIn[n].z()
                );

                mPointsOut[n] = p;
            }
        }

    private:
        const PointList& mPointsIn;
        pointField& mPointsOut;
};

///// TBB body object for threaded world to voxel space transformation and copy of points
class TransformOp
{
public:
    TransformOp
    (
        triSurface const * const gdp,
        const openvdb::math::Transform& transform,
        std::vector<openvdb::Vec3s>& pointList
    );

    void operator()(const tbb::blocked_range<size_t>&) const;

private:
    triSurface const * const mGdp;
    const openvdb::math::Transform& mTransform;
    std::vector<openvdb::Vec3s>* const mPointList;
};

TransformOp::TransformOp
(
    triSurface const * const gdp,
    const openvdb::math::Transform& transform,
    std::vector<openvdb::Vec3s>& pointList
)
:
    mGdp(gdp),
    mTransform(transform),
    mPointList(&pointList)
{
}

void
TransformOp::operator()(const tbb::blocked_range<size_t>& r) const
{
    point pos;
    openvdb::Vec3d ipos;

    for (size_t i = r.begin(); i != r.end(); i++)
    {
        pos = mGdp->localPoints()[i];
        ipos = mTransform.worldToIndex(openvdb::Vec3d(pos.x(), pos.y(), pos.z()));

        (*mPointList)[i] = openvdb::Vec3s(ipos);
    }
}

///// @brief   TBB body object for threaded primitive copy
class PrimCpyOp
{
public:
    PrimCpyOp
    (
        triSurface const * const gdp,
        std::vector<openvdb::Vec4I>& primList
    );

    void operator()(const tbb::blocked_range<size_t>&) const;

private:
    triSurface const * const mGdp;
    std::vector<openvdb::Vec4I>* const mPrimList;
};

PrimCpyOp::PrimCpyOp
(
    triSurface const * const gdp,
    std::vector<openvdb::Vec4I>& primList
)
:
    mGdp(gdp),
    mPrimList(&primList)
{
}


void
PrimCpyOp::operator()(const tbb::blocked_range<size_t>& r) const
{
    openvdb::Vec4I prim;

    for (size_t i = r.begin(); i != r.end(); i++)
    {
        const labelledTri& tri = mGdp->localFaces()[i];

        forAll(tri, j)
        {
            prim[j] = tri[j];
        }

        prim[3] = openvdb::util::INVALID_IDX;

        (*mPrimList)[i] = prim;
    }
}

/// TBB body object for threaded sharp feature construction
class SharpenFeaturesOp
{
public:
    using EdgeData = openvdb::tools::MeshToVoxelEdgeData;

    SharpenFeaturesOp
    (
        pointField& meshGeo,
        const triSurface& refGeo,
        EdgeData& edgeData,
        const openvdb::math::Transform& xform,
        const openvdb::BoolTree* mask = nullptr
    );

    void operator()(const tbb::blocked_range<size_t>&) const;

private:
    pointField& mMeshGeo; // points of the resulting surface mesh
    const triSurface& mRefGeo; // is inputSurf
    EdgeData& mEdgeData;
    const openvdb::math::Transform& mXForm;
    const openvdb::BoolTree* mMaskTree;
};

SharpenFeaturesOp::SharpenFeaturesOp
(
    pointField& meshGeo,
    const triSurface& refGeo,
    EdgeData& edgeData,
    const openvdb::math::Transform& xform,
    const openvdb::BoolTree * mask
)
:
    mMeshGeo(meshGeo),
    mRefGeo(refGeo),
    mEdgeData(edgeData),
    mXForm(xform),
    mMaskTree(mask)
{
}

void
SharpenFeaturesOp::operator()
(
    const tbb::blocked_range<size_t>& range
) const
{
    openvdb::tools::MeshToVoxelEdgeData::Accessor acc = mEdgeData.getAccessor();

    using BoolAccessor = openvdb::tree::ValueAccessor<const openvdb::BoolTree>;
    std::unique_ptr<BoolAccessor> maskAcc;

    if (mMaskTree) {
        maskAcc.reset(new BoolAccessor(*mMaskTree));
    }

    label ptnOffset;

    point tmpN, tmpP, avgP;
    treeBoundBox cell;

    openvdb::Vec3d pos, normal;
    openvdb::Coord ijk;

    std::vector<openvdb::Vec3d> points(12), normals(12);
    std::vector<openvdb::Index32> primitives(12);

    for (size_t ptnOffset = range.begin(); ptnOffset != range.end(); ptnOffset++)
    {
        tmpP = mMeshGeo[ptnOffset];
        pos[0] = tmpP.x();
        pos[1] = tmpP.y();
        pos[2] = tmpP.z();

        pos = mXForm.worldToIndex(pos);

        ijk[0] = int(std::floor(pos[0]));
        ijk[1] = int(std::floor(pos[1]));
        ijk[2] = int(std::floor(pos[2]));

        if (maskAcc && !maskAcc->isValueOn(ijk)) continue;

        points.clear();
        normals.clear();
        primitives.clear();

        // get voxel-edge intersections
        mEdgeData.getEdgeData(acc, ijk, points, primitives);

        avgP = point(0.0, 0.0, 0.0);

        // get normal list
        for (size_t n = 0, N = points.size(); n < N; ++n) {

            avgP.x() = static_cast<float>(avgP.x() + points[n].x());
            avgP.y() = static_cast<float>(avgP.y() + points[n].y());
            avgP.z() = static_cast<float>(avgP.z() + points[n].z());

            tmpN = mRefGeo.faceNormals()[primitives[n]];

            normal[0] = tmpN.x();
            normal[1] = tmpN.y();
            normal[2] = tmpN.z();

            normals.push_back(normal);
        }

        // Calculate feature point position
        if (points.size() > 1)
        {

            pos = openvdb::tools::findFeaturePoint(points, normals);

            // Constrain points to stay inside their initial
            // coordinate cell.
            cell =
                treeBoundBox
                (
                    point(double(ijk[0]), double(ijk[1]), double(ijk[2])),
                    point(double(ijk[0]+1), double(ijk[1]+1), double(ijk[2]+1))
                );

            //cell.expandBounds(0.3, 0.3, 0.3);
            //cell.inflate(0.3);
            cell.inflate(0.6);

            if (!cell.contains(point(pos[0], pos[1], pos[2])))
            {
                point org
                (
                    static_cast<float>(pos[0]),
                    static_cast<float>(pos[1]),
                    static_cast<float>(pos[2])
                );

                avgP *= 1.f / float(points.size());
                //UT_Vector3 dir = avgP - org;
                //dir.normalize();

                //double distance;
                point intersection;

                //if (cell.intersectRay(org, dir, 1E17F, &distance) > 0)
                if (cell.intersects(org, avgP, intersection))
                {

                    //tmpP = org + dir * distance;
                    tmpP = intersection;

                    pos[0] = tmpP.x();
                    pos[1] = tmpP.y();
                    pos[2] = tmpP.z();
                }
            }

            pos = mXForm.indexToWorld(pos);

            tmpP.x() = static_cast<float>(pos[0]);
            tmpP.y() = static_cast<float>(pos[1]);
            tmpP.z() = static_cast<float>(pos[2]);

            mMeshGeo[ptnOffset] = tmpP;
        }
    } //for range.begin
} //SharpenFeaturesOp::operator()


/// TBB body object for threaded sharp feature construction
template<typename IndexTreeType, typename BoolTreeType>
class GenAdaptivityMaskOp
{
public:
    using BoolLeafManager = openvdb::tree::LeafManager<BoolTreeType>;

    GenAdaptivityMaskOp
    (
        const triSurface& refGeo,
        const IndexTreeType& indexTree,
        BoolLeafManager&,
        float edgetolerance = 0.0
    );

    void run(bool threaded = true);

    void operator()(const tbb::blocked_range<size_t>&) const;

private:
    const triSurface& mRefGeo;
    const IndexTreeType& mIndexTree;
    BoolLeafManager& mLeafs;
    float mEdgeTolerance;
};


template<typename IndexTreeType, typename BoolTreeType>
GenAdaptivityMaskOp<IndexTreeType, BoolTreeType>::GenAdaptivityMaskOp
(
    const triSurface& refGeo,
    const IndexTreeType& indexTree,
    BoolLeafManager& leafMgr,
    float edgetolerance
)
:
    mRefGeo(refGeo),
    mIndexTree(indexTree),
    mLeafs(leafMgr),
    mEdgeTolerance(edgetolerance)
{
    mEdgeTolerance = std::max(0.0f, mEdgeTolerance);
    mEdgeTolerance = std::min(1.0f, mEdgeTolerance);
}


template<typename IndexTreeType, typename BoolTreeType>
void
GenAdaptivityMaskOp<IndexTreeType, BoolTreeType>::run(bool threaded)
{
    if (threaded) {
        tbb::parallel_for(mLeafs.getRange(), *this);
    } else {
        (*this)(mLeafs.getRange());
    }
}


template<typename IndexTreeType, typename BoolTreeType>
void
GenAdaptivityMaskOp<IndexTreeType, BoolTreeType>::operator()
(
    const tbb::blocked_range<size_t>& range
) const
{
    using IndexAccessorType = typename openvdb::tree::ValueAccessor<const IndexTreeType>;
    IndexAccessorType idxAcc(mIndexTree);

    point tmpN, normal;
    label tmpIdx;

    openvdb::Coord ijk, nijk;
    typename BoolTreeType::LeafNodeType::ValueOnIter iter;

    for (size_t n = range.begin(); n < range.end(); ++n)
    {
        iter = mLeafs.leaf(n).beginValueOn();

        for (; iter; ++iter)
        {
            ijk = iter.getCoord();

            bool edgeVoxel = false;

            int idx = idxAcc.getValue(ijk);

            normal = mRefGeo.faceNormals()[idx];

            for (size_t i = 0; i < 18; ++i)
            {
                nijk = ijk + openvdb::util::COORD_OFFSETS[i];

                if (idxAcc.probeValue(nijk, tmpIdx) && tmpIdx != idx)
                {
                    tmpN = mRefGeo.faceNormals()[tmpIdx];

                    if ((normal & tmpN) < mEdgeTolerance)
                    {
                        edgeVoxel = true;
                        break;
                    }
                }
            }

            if (!edgeVoxel) iter.setValueOff();
        }
    }
}

#endif // GEOMETRY_UTIL
