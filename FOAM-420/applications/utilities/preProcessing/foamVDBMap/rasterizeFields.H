struct RasterizationData
{
    using Ptr = std::unique_ptr<RasterizationData>;
    using FloatTree    = FloatGrid::TreeType;
    using Vec3STree    = Vec3SGrid::TreeType;
    using FloatTreeAcc = openvdb::tree::ValueAccessor<FloatTree>;
    using Vec3STreeAcc = openvdb::tree::ValueAccessor<Vec3STree>;

    std::vector<FloatTree>    scalarTrees;
    std::vector<Vec3STree>    vectorTrees;

    //std::vector<FloatTreeAcc> scalarAccessors;
    //std::vector<Vec3STreeAcc> vectorAccessors;

    RasterizationData
    (
        const label nScalarFields,
        const label nVectorFields
    )
    :
        scalarTrees(nScalarFields, FloatTree(SMALL)),
        vectorTrees(nVectorFields, Vec3STree(Vec3SGrid::ValueType(SMALL)))//,
      //  scalarAccessors(),
      //  vectorAccessors()
    {
        //scalarAccessors.resize(nScalarFields);
        //for (size_t i = 0; i <= nScalarFields; i++)
        //{
        //    scalarAccessors.push_back
        //    (
        //        FloatTreeAcc(scalarTrees[i])
        //    );
        //}

        //vectorAccessors.resize(nVectorFields);
        //for (size_t i = 0; i <= nVectorFields; i++)
        //{
        //    vectorAccessors.push_back
        //    (
        //        Vec3STreeAcc(vectorTrees[i])
        //    );
        //}
    }
}; // RasterizationData

class Rasterizer
{
public:

    using DataTable = tbb::enumerable_thread_specific<RasterizationData::Ptr>;

    Rasterizer
    (
        DataTable& dataTable,
        const cellList& cells,
        const faceList& faces,
        const pointField& points,
        const List<scalarField>& scalarFields,
        const List<vectorField>& vectorFields,
        const openvdb::CoordBBox& targetBBox,
        openvdb::math::Transform::ConstPtr transformPtr
    )
    :
        dataTable_(&dataTable),
        cells_(cells),
        faces_(faces),
        points_(points),
        scalarFields_(scalarFields),
        vectorFields_(vectorFields),
        targetBBox_(targetBBox),
        transformPtr_(transformPtr)
    {
        ////if (debug)
        //if (Pstream::master())
        //{
        //    std::cout<< "    Rasterizer thread ID "
        //        << tbb::this_tbb_thread::get_id()
        //        << std::endl;
        //}
    }

    void operator()(const tbb::blocked_range<size_t>& range) const
    {
        typename RasterizationData::Ptr& dataPtr = dataTable_->local();
        if (!dataPtr) dataPtr.reset(new RasterizationData(scalarFields_.size(), vectorFields_.size()));

        RasterizationData& data = *dataPtr;

        for (size_t cellI = range.begin(); cellI < range.end(); cellI++)
        {
            //TODO optimization here is having boundbox of cells
            // already precalculated
            boundBox bb
            (
                cells_[cellI].points
                (
                    faces_,
                    points_
                ),
                /*doReduce*/false
            );

            openvdb::CoordBBox cellBBox =
                toCoordBBox
                (
                    bb,
                    transformPtr_
                );

            if (targetBBox_.hasOverlap(cellBBox))
            {
                // scalars
                forAll(scalarFields_, i)
                {
                    scalar value = scalarFields_[i][cellI];

                    //if (cellBBox.volume() == 1)
                    //{
                    //    data.scalarAccessors[i].setValue
                    //    (
                    //        cellBBox.min(),
                    //        value
                    //    );
                    //}
                    //else
                    //{
                        data.scalarTrees[i].fill
                        (
                            cellBBox,
                            value
                        );
                    //}
                }

                // vectors
                forAll(vectorFields_, i)
                {
                    const vector& v = vectorFields_[i][cellI];

                    Vec3SGrid::ValueType value (v.x(), v.y(), v.z());

                    //if (cellBBox.volume() == 1)
                    //{
                    //    data.vectorAccessors[i].setValue
                    //    (
                    //        cellBBox.min(),
                    //        value
                    //    );
                    //}
                    //else
                    //{
                        data.vectorTrees[i].fill
                        (
                            cellBBox,
                            value
                        );
                    //}
                }
            } // if targetBBox overlaps cell

            //GridType::Accessor acc(grid_.getAccessor());
            //openvdb::Coord xyz;
            //for (label x = cellBBox.min().x(); x <= cellBBox.max().x(); x++)
            //{
            //    xyz.setX(x);
            //    for (label y = cellBBox.min().y(); y <= cellBBox.max().y(); y++)
            //    {
            //        xyz.setY(y);
            //        for (label z = cellBBox.min().z(); z <= cellBBox.max().z(); z++)
            //        {
            //            xyz.setZ(z);
            //            acc.setValueOn(xyz, fieldSource_[cellI]);
            //        } //for z
            //    } // for y
            //} // for x
        } // for cellI
    } // operator ()

private:

    static openvdb::CoordBBox toCoordBBox
    (
        const boundBox& bb,
        openvdb::math::Transform::ConstPtr transform
    )
    {
        const openvdb::Coord bbMin
        (
            transform->worldToIndexNodeCentered
            (
                openvdb::Vec3s
                (
                    bb.min().x(),
                    bb.min().y(),
                    bb.min().z()
                )
            )
        );
        const openvdb::Coord bbMax
        (
            transform->worldToIndexNodeCentered
            (
                openvdb::Vec3s
                (
                    bb.max().x(),
                    bb.max().y(),
                    bb.max().z()
                )
            )
        );

        return openvdb::CoordBBox(bbMin, bbMax);
    }

    DataTable* const dataTable_;
    const cellList& cells_;
    const faceList& faces_;
    const pointField& points_;
    const List<scalarField>& scalarFields_;
    const List<vectorField>& vectorFields_;
    const openvdb::CoordBBox& targetBBox_;
    openvdb::math::Transform::ConstPtr transformPtr_;
}; // Rasterizer


template<typename GridType>
inline void mergeTrees
(
    std::vector<typename GridType::Ptr>& grids,
    std::vector<typename GridType::TreeType>& dataTrees
)
{
    using GridPtr = typename GridType::Ptr;
    using TreeType = typename GridType::TreeType;

    for (size_t i = 0; i < grids.size(); i++)
    {
        GridPtr gridPtr = grids[i];
        TreeType& dataTree = dataTrees[i];

        gridPtr->tree().merge(dataTree);
    }
} // mergeTrees


template<typename GridType>
inline void mergeTrees
(
    tbb::task_group& tasks,
    std::vector<typename GridType::Ptr>& grids,
    std::vector<typename GridType::TreeType>& dataTrees
)
{
    using GridPtr = typename GridType::Ptr;
    using TreeType = typename GridType::TreeType;

    for (size_t i = 0; i < grids.size(); i++)
    {
        GridPtr gridPtr = grids[i];
        TreeType& dataTree = dataTrees[i];

        tasks.run
        (
            [gridPtr, &dataTree]
            {
                gridPtr->tree().merge(dataTree);
            }
        );
    }
} // mergeTrees

template<typename polyMeshBasic>
inline void rasterizeFields
(
  const polyMeshBasic& mesh,
  const List<scalarField>& scalarFields,
  const List<vectorField>& vectorFields,
  const openvdb::CoordBBox& targetBBox,
  std::vector<FloatGrid::Ptr>& scalarGrids,
  std::vector<Vec3SGrid::Ptr>& vectorGrids
)
{
    using DataTable = tbb::enumerable_thread_specific<RasterizationData::Ptr>;
    using FloatTree = FloatGrid::TreeType;
    using Vec3STree = Vec3SGrid::TreeType;

    DataTable data;

    const tbb::blocked_range<size_t> cellRange(0, mesh.nCells());

    const openvdb::CoordBBox paddedBB = targetBBox.expandBy(4);

    tbb::tick_count t0 = tbb::tick_count::now();

    tbb::parallel_for
    (
        cellRange,
        Rasterizer
        (
            data,
            mesh.cells(),
            mesh.faces(),
            mesh.points(),
            scalarFields,
            vectorFields,
            paddedBB,
            (
                scalarGrids.size()
              ? scalarGrids[0]->transformPtr()
              : vectorGrids[0]->transformPtr()
            )
        )
    );

    //tbb::tick_count t1 = tbb::tick_count::now();
    //Info<<"Rasterized in " << (t1-t0).seconds() << endl;

    // merge trees
    for (typename DataTable::iterator dataI = data.begin(); dataI != data.end(); ++dataI)
    {
        RasterizationData& dataItem = **dataI;

        tbb::task_group tasks;

        mergeTrees<FloatGrid>(tasks, scalarGrids, dataItem.scalarTrees);
        mergeTrees<Vec3SGrid>(tasks, vectorGrids, dataItem.vectorTrees);

        tasks.wait();
    }

    //tbb::tick_count t2 = tbb::tick_count::now();
    //Info<< "Merged trees in " << (t2-t1).seconds()
    //    << " - total: " << (t2-t0).seconds()
    //    << endl;
}
