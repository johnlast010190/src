
    triSurface inputSurf;
    forAll(surfaceGeometry, surfi)
    {
        label geomi = surfaceGeometry[surfi];
        const searchableSurface& s = allGeometryPtr()[geomi];

        if (isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& triMesh = refCast<const triSurfaceMesh>(s);
            //triSurface& inputSurf = const_cast<triSurface&>
            inputSurf = const_cast<triSurface&>
            (
                refCast<const triSurface>(triMesh)
            );
        }
    }

    inputSurf.translatePoints(levelSetOffset);

    openvdb::GridPtrVec vdbGrids;

    FloatGrid::Ptr highResGrid =
        ovdb.meshToLevelSet
        (
            inputSurf,
            maxLevel,
            (nCellsBetweenLevels << maxLevel)
        );
    highResGrid->setName("highResLevelSet");

    openvdb::CoordBBox highResCoordBB =
        toCoordBBox
        (
            samplingBBox,
            highResGrid->transformPtr()
        );

    highResGrid->clip(highResCoordBB);

    highResGrid->insertMeta
    (
        "samplingBoxCoordMin",
        openvdb::Vec3IMetadata
        (
            highResCoordBB.min().asVec3i()
        )
    );
    highResGrid->insertMeta
    (
        "samplingBoxCoordMax",
        openvdb::Vec3IMetadata
        (
            highResCoordBB.max().asVec3i()
        )
    );

    if (saveHighResGrid || maxLevel == 0)
    {
        vdbGrids.push_back(highResGrid);
    }

    Info<< "\nActive voxels in high resolution narrow-band levelset "
        << setw(10) << highResGrid->activeVoxelCount()
        << endl;

    if (maxLevel > 0)
    {
        typedef openvdb::tools::MultiResGrid<openvdb::FloatTree> MultiResGrid;

        const size_t levels = maxLevel + 1;

        // Generate LOD sequence from finest level
        MultiResGrid mrg(levels, *highResGrid, /* reduction by injection */ false );

        FloatGrid::Ptr coarsestGrid = mrg.grid(mrg.coarsestLevel());

        Info<< "Active voxels in low resolution narrow-band levelset  "
            << setw(10) << coarsestGrid->activeVoxelCount()
            << endl;

        coarsestGrid->setName("lowResLevelSet");

        openvdb::CoordBBox coarsestGridCoordBB =
            toCoordBBox
            (
                samplingBBox,
                coarsestGrid->transformPtr()
            );

        coarsestGridCoordBB.max().x() = coarsestGridCoordBB.min().x() + outputDivisions[0]-1;
        coarsestGridCoordBB.max().y() = coarsestGridCoordBB.min().y() + outputDivisions[1]-1;
        coarsestGridCoordBB.max().z() = coarsestGridCoordBB.min().z() + outputDivisions[2]-1;

        coarsestGrid->clip(coarsestGridCoordBB);

        coarsestGrid->insertMeta
        (
            "samplingBoxCoordMin",
            openvdb::Vec3IMetadata
            (
                coarsestGridCoordBB.min().asVec3i()
            )
        );
        coarsestGrid->insertMeta
        (
            "samplingBoxCoordMax",
            openvdb::Vec3IMetadata
            (
                coarsestGridCoordBB.max().asVec3i()
            )
        );

        if (normalizeLevelset)
        {
            openvdb::tools::foreach
            (
                coarsestGrid->beginValueOn(),
                Norm(coarsestGrid->background())
            );
        }

        vdbGrids.push_back(coarsestGrid);

        cellLevelGrids[0] = coarsestGrid->deepCopy();
    }
    else
    {
        if (normalizeLevelset)
        {
            openvdb::tools::foreach
            (
                highResGrid->beginValueOn(),
                Norm(highResGrid->background())
            );
        }

        cellLevelGrids[0] = highResGrid->deepCopy();
    }

    if (visualizeLevelSet)
    {
        Info<< "\nWriting distanceField in OpenFOAM format for visualization" << endl;

        openvdb::CoordBBox bBox =
            toCoordBBox
            (
                samplingBBox,
                cellLevelGrids[0]->transformPtr()
            );

        //bBox.max().offset(0, 0, 1);

        openvdb::MaskGrid mask;
        mask.setTransform(cellLevelGrids[0]->transformPtr());
        mask.sparseFill(bBox, true);

        if (maskSdf)
        {
            cellLevelGrids[0] = openvdb::tools::maskSdf(*cellLevelGrids[0], mask);
        }
        else
        {
            cellLevelGrids[0]->topologyUnion(mask);
        }

        if (normalizeLevelset)
        {
            openvdb::tools::foreach
            (
                cellLevelGrids[0]->beginValueOn(),
                Norm(cellLevelGrids[0]->background())
            );
        }

        cellLevelGrids[0]->clip(bBox);

        cellLevelGrids[0]->setName("distanceField");

        Info<< "Active voxels in low resolution distanceField  . . . .  "
            << setw(10) << cellLevelGrids[0]->activeVoxelCount()
            << endl;

        distanceGrid = cellLevelGrids[0];

        // do not update processor boundaries
        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        //set output time to serial case
        Foam::Time runTime
        (
            Foam::Time::controlDictName,
            args.rootPath(),
            args.globalCaseName()
        );

        autoPtr<fvMesh> meshPtr =
            writeMesh
            (
                cellLevelGrids,
                ovdb,
                runTime,
                levelSetOffset
            );

        writeVolFields<scalar, FloatGrid>(meshPtr(), cellLevelGrids);

        UPstream::parRun() = oldParRun;
    }

    const word vdbOutput = sampleGeometryDict.lookupOrDefault<word>("levelsetOutput", "geometry.vdb");

    Info<< "\nSaving "
        << ((saveHighResGrid || maxLevel == 0) ? "highResLevelSet " : "")
        << ((saveHighResGrid && maxLevel > 0)  ? "and " : "")
        << (maxLevel > 0 ? "lowResLevelSet " : "")
        << "to " << vdbOutput << endl;

    openvdb::io::File file(vdbOutput);
    file.write(vdbGrids);
    file.close();

