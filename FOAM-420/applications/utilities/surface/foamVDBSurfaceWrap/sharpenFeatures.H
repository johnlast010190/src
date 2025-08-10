// Setup level set mesher
openvdb::tools::VolumeToMesh mesher(iso, adaptivity);

//SOP_OpenVDB_To_Polygons.cc:619
//mesher.setSurfaceMask(maskGridPtr, invertmask);
//A boolean grid whose active topology defines the region to mesh.
//
//SOP_OpenVDB_To_Polygons.cc:641
//mesher.setSpatialAdaptivity(grid);
//A scalar grid used as a spatial multiplier for the adaptivity threshold.



//SOP_OpenVDB_To_Polygons.cc:693
//std::list<openvdb::GridBase::ConstPtr> grids;
//referenceMeshing<openvdb::FloatGrid>(grids, mesher, refGeo, boss, time);

//referenceMeshing
openvdb::tools::MeshToVoxelEdgeData edgeData;

openvdb::math::Transform::Ptr transform = surfLevelSet->transform().copy();
const FloatGrid::ValueType backgroundValue = surfLevelSet->background();
const openvdb::GridClass gridClass = surfLevelSet->getGridClass();

using IntGridT = typename FloatGrid::template ValueConverter<openvdb::Int32>::Type;
typename IntGridT::Ptr indexGrid; // replace
indexGrid.reset(new IntGridT(0));


std::vector<openvdb::Vec3s> pointList; //points of original surface in index space
std::vector<openvdb::Vec4I> primList;  //quads and tri of original shape

pointList.resize(globalSurf.localPoints().size());
primList.resize(globalSurf.localFaces().size());

tbb::parallel_for
(
    tbb::blocked_range<size_t>(0, pointList.size()),
    TransformOp(&globalSurf, *transform, pointList)
);

tbb::parallel_for
(
    tbb::blocked_range<size_t>(0, primList.size()),
    PrimCpyOp(&globalSurf, primList)
);

openvdb::tools::QuadAndTriangleDataAdapter<openvdb::Vec3s, openvdb::Vec4I>
    mesh(pointList, primList); // original mesh (globalSurf) in openVDB format

//float bandWidth = 3.0;
//
//if (gridClass != openvdb::GRID_LEVEL_SET)
//{
//    bandWidth = float(backgroundValue) / float(transform->voxelSize()[0]);
//}
////else
////{
////Info<<"GGG openvdb::GRID_LEVEL_SET " <<endl;
////}
//
//FloatGrid::ConstPtr refGrid = // levelSet of original mesh
//    openvdb::tools::meshToVolume<FloatGrid>
//    (
//        mesh, // pointList and primList (globalSurf) in index space
//        *transform,
//        bandWidth,
//        bandWidth,
//        0,              //  optional conversion flags defined in MeshToVolumeFlags
//        indexGrid.get() // optional grid output that will contain the closest-polygon index for each voxel in the active narrow band region
//    );


using TreeType = typename FloatGrid::TreeType;
using BoolTreeType = typename TreeType::template ValueConverter<bool>::Type;
typename BoolTreeType::Ptr maskTree;

//////if (sharpenFeatures) {
//    maskTree = typename BoolTreeType::Ptr(new BoolTreeType(false));
//    maskTree->topologyUnion(indexGrid->tree());
//    openvdb::tree::LeafManager<BoolTreeType> maskLeafs(*maskTree);
//
//    GenAdaptivityMaskOp<typename IntGridT::TreeType, BoolTreeType> op
//    (
//        globalSurf, //*refGeo,
//        indexGrid->tree(),
//        maskLeafs,
//        edgetolerance
//    );
//    op.run();
//
//    openvdb::tools::pruneInactive(*maskTree);
//
//    openvdb::tools::dilateActiveValues(*maskTree, 2,
//                openvdb::tools::NN_FACE, openvdb::tools::IGNORE_TILES);
//
//    // A boolean tree whose active topology defines the adaptivity mask.
//    mesher.setAdaptivityMask(maskTree);
//////}


const double iadaptivity = 0.1; //TODO //read from dictionary
// When surfacing fractured SDF fragments, the original unfractured SDF grid can be used to eliminate seam lines and tag polygons that are coincident with the reference surface with the POLYFLAG_EXTERIOR flag and polygons that are in proximity to the seam lines with the POLYFLAG_FRACTURE_SEAM flag. (The performance cost for using this reference based scheme compared to the regular meshing scheme is approximately 15% for the first fragment and neglect-able for subsequent fragments.)
//mesher.setRefGrid(refGrid, iadaptivity);
mesher.setRefGrid(surfLevelSet, iadaptivity);

//mesher(*refGrid); //VolumeToMesh
{ Timer timer("volumeToMesh");
mesher(*surfLevelSet); //VolumeToMesh
}

//copyMesh(*gdp, mesher, boss, usePolygonSoup,
//    keepVdbName ? grid->getName().c_str() : nullptr,
//    surfaceGroup, interiorGroup, seamGroup, seamPointGroup);

// Set patch name/index.
geometricSurfacePatchList patches(nRegions);
forAll(patches, patchi)
{
    patches[patchi] = geometricSurfacePatch(allRegNames[patchi], patchi);
}

// Set points
pointField of_points(mesher.pointListSize());

tbb::parallel_for
(
    tbb::blocked_range<size_t>(0, of_points.size()),
    PointListCopy
    (
        mesher.pointList(),
        of_points
    )
);

if (sharpenFeatures)
{
    Timer timer("level-set to mesh sharpenFeatures");

    edgeData.convert(pointList, primList);

    // move points to sharp features
    tbb::parallel_for
    (
        tbb::blocked_range<size_t>(0, of_points.size()),
        SharpenFeaturesOp
        (
            of_points,
            globalSurf,
            edgeData,
            *transform,
            maskTree.get()
        )
    );
}

List<labelledTri> faces;

{ Timer t("copy primitives");
// Copy faces
openvdb::tools::PolygonPoolList& polygonPoolList = mesher.polygonPoolList();

size_t numPrimitives = 0;
std::vector<size_t> indexlist(mesher.polygonPoolListSize());

for (size_t n = 0, N = mesher.polygonPoolListSize(); n < N; ++n)
{
    const openvdb::tools::PolygonPool& polygons = polygonPoolList[n];
    indexlist[n] = numPrimitives;
    numPrimitives += 2*polygons.numQuads();
    numPrimitives += polygons.numTriangles();
}

faces.resize(numPrimitives);
PrimCpy primCpy(polygonPoolList, indexlist, faces);
primCpy.runParallel();
}


/////////////////////////

//remeshed surface with sharpen features
triSurface wrapSurf
    (
        faces,
        patches,
        of_points
    );

{ Timer timer("mergePoints");
wrapSurf = triSurfaceTools::mergePoints(wrapSurf, 0.1*minEdge);
}

// offset to voxel centre
wrapSurf.translatePoints(point(0.5*minEdge));
