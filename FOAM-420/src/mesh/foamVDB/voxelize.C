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
    (c) 2021 ESI

\*---------------------------------------------------------------------------*/

#include "foamVDB.H"

////clash with macro expansion in <openvdb_dir>/include/openvdb/math/Vec3.h:663
//#undef Log
#include "MeshToVolume.h"
#include "MultiResGrid.h"
#include <openvdb/tools/Morphology.h> //for dilation and erosion
#include <openvdb/tools/Composite.h>
//#define Log if (log) Info

#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "searchableSurfaces/searchableBox/searchableBox.H"
#include "distributedTriSurfaceMesh/distributedTriSurfaceMesh.H"
#include "meshRefinement/meshRefinement.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "shellSurfaces/shellSurfaces.H"
#include "snappyHexMeshDriver/refinementParameters/refinementParameters.H"

#ifdef FOAM_USE_OPENCASCADE
  #include "distributedCADSurface/distributedCADSurface.H"
#endif


void Foam::foamVDB::voxelize
(
    meshRefinement& meshRefiner,
    const refinementParameters& refineParams,
    const dictionary& vdbDict
)
{
    const refinementSurfaces& surfaces = meshRefiner.surfaces();
    const shellSurfaces&      shells   = meshRefiner.shells();
    //TODO read refinement features too?
    //const refinementFeatures& features = meshRefiner.features()();

    fvMesh& mesh = meshRefiner.mesh();

    wordList allRegNames;
    wordList allRegNamesAndSlaves = mesh.boundaryMesh().names();

    // remove "*_slave" and VDBdomain faces from region names
    forAll(allRegNamesAndSlaves, i)
    {
        if (i < 6) continue; //first 6 patches are VDBdomain

        const word& regName = allRegNamesAndSlaves[i];

        if (regName.endsWith("_slave")) continue;

        allRegNames.append(regName);
    }

    const label nRegions = allRegNames.size();

    const label nCellsBetweenLevels = max(3, refineParams.nBufferLayers());
    const word  vdbDecompMethod     = refineParams.vdbDecompMethod();

    //create hierarchicalDecomposition for VDB grids decomposition
    const dictionary decompositionDict =
        createDecompDict
        (
            mesh.time(),
            vdbDecompMethod //wanted VDB decomposition method
        );

    // ### Initialize storage ### //

    // for every boundary patch, vector of Unsigned Distance Field
    // of refinement grading for every cellLevel
    // (filled only with patch at minLevel-1)
    std::vector<std::vector<FloatGrid::Ptr>> bufferCellLevels
    (
        nRegions,
        std::vector<FloatGrid::Ptr>
        (
            maxCellLevel_ + 1,
            nullptr
        )
    );

    // filled only with patch at minLevel
    std::vector<std::vector<FloatGrid::Ptr>> patchCellLevels
    (
        nRegions,
        std::vector<FloatGrid::Ptr>
        (
            maxCellLevel_ + 1,
            nullptr
        )
    );

    // collection of all following grids for each cell level:
    // - background voxels at cellLevel 0
    // - refinementRegions (shells)
    // - levelSet of each patch with relative buffer (nCellsBetweenLevels)
    // - levelSet of curvature (edges) with buffer
    std::vector<std::vector<FloatGrid::Ptr>> cellLevelGridsCollection;

    // collection of interior of each refinement grading for each of the grids above
    std::vector<std::vector<BoolGrid::Ptr>> cellLevelInteriorCollection;

    // For each cell level, resulting grids from topology union/intersection
    // of the above collection of grids
    std::vector<FloatGrid::Ptr> cellLevelGrids(maxCellLevel_ + 1, nullptr);

    // For each cell level, resulting grids from topology union
    // of cellLevelInteriorCollection
    std::vector<BoolGrid::Ptr> interiorGrids(maxCellLevel_ + 1, nullptr);

    // For each cell level, displacement field of boundary points
    std::vector<Vec3dGrid::Ptr> surfaceDisplacementGrids(maxCellLevel_ + 1, nullptr);
    std::vector<Vec3dGrid::Ptr> featureDisplacementGrids(maxCellLevel_ + 1, nullptr);

    // Grids in proximity of surfaces to snap
    std::vector<BoolGrid::Ptr> innerGrids(maxCellLevel_ + 1, nullptr);
    std::vector<BoolGrid::Ptr> outerGrids(maxCellLevel_ + 1, nullptr);

    // ### Setup background voxels (level 0) ### //

    const vector scale(vdbDict.lookupOrDefault<vector>("scaleVoxel", vector::one));
    const vector invScale = cmptDivide(vector::one, scale);

    vector origin(Zero);
    quaternion R    = quaternion(vector::one, 0.0);
    quaternion invR = quaternion(vector::one, 0.0);

    if (vdbDict.found("coordinateSystem"))
    {
        const dictionary& coordDict = vdbDict.subDict("coordinateSystem");

        coordinateSystem coord("cartesian", coordDict);
        origin = coord.origin();
        R = quaternion(coord.R());

        if (coordDict.found("yawPitchRoll"))
        {
            vector v = coordDict.lookup("yawPitchRoll");
            v *= Foam::constant::mathematical::pi/180.0;

            scalar yaw   = v.x();
            scalar pitch = v.y();
            scalar roll  = v.z();

            R  = quaternion(vector(0, 0, 1), yaw);
            R *= quaternion(vector(0, 1, 0), pitch);
            R *= quaternion(vector(1, 0, 0), roll);

            invR  = quaternion(vector(-1, 0, 0), roll);
            invR *= quaternion(vector(0, -1, 0), pitch);
            invR *= quaternion(vector(0, 0, -1), yaw);
        }
        else if (coordDict.found("rollPitchYaw"))
        {
            vector v = coordDict.lookup("rollPitchYaw");
            v *= Foam::constant::mathematical::pi/180.0;

            R = quaternion(v.x(), v.y(), v.z());

            scalar yaw   = v.z();
            scalar pitch = v.y();
            scalar roll  = v.x();

            invR  = quaternion(vector(0, 0, -1), yaw);
            invR *= quaternion(vector(0, -1, 0), pitch);
            invR *= quaternion(vector(-1, 0, 0), roll);
        }

        Info<<"GGG coord " << coord<< endl;
        Info<<"GGG rotation " << R<< endl;
        Info<<"GGG invRotation " << invR<< endl;
    }

    point levelSetOffset;

    // boundBox used to clip all the other VDB grids
    openvdb::CoordBBox boundingBox;

    const scalar level0Edge = voxelSize_ * Foam::pow(2, maxCellLevel_);

    {
        //Timer t("boxToLevelSet");

        point boundMin = cmptMultiply(point(vdbDict.lookup("min")), invScale);
        point boundMax = cmptMultiply(point(vdbDict.lookup("max")), invScale);

        // shift xmin, ymin, zmin of wind tunnel to closest integer value
        // will shift back when converting points from index space to world space
        levelSetOffset =
            point
            (
                -std::fmod(boundMin.x(), level0Edge),
                -std::fmod(boundMin.y(), level0Edge),
                -std::fmod(boundMin.z(), level0Edge)
            );

        boundMin += levelSetOffset;
        boundMax += levelSetOffset;

        // cells are built using x+1 voxels, so reduce boundMax by one to
        // avoid having a final larger domain
        boundMax -= point(level0Edge);

        std::vector<FloatGrid::Ptr> boundsCellLevels(maxCellLevel_ + 1, nullptr);

        boundsCellLevels[0] = searchableBoxToLevelSet(boundMin, boundMax, /*cellLevel*/0);

        boundingBox = boundsCellLevels[0]->evalActiveVoxelBoundingBox();

        //if (debug)
        //{
        //    Info<< "\nlevelSetOffset " << levelSetOffset
        //        << nl << "boundMin " << boundMin
        //        << nl << "boundMax " << boundMax
        //        << endl;
        //    std::cout << "boundingBox " << boundingBox << std::endl;
        //}

        cellLevelGridsCollection.push_back(boundsCellLevels);
    }

    // ### Voxelize refinement regions ### //
    {
        Timer timer("Voxelize refinement regions");

        //add shell refinement to cellLevelGridsCollection
        forAll(shells.shells(), shelli)
        {
            label geomi = shells.shells()[shelli];

            const word& shellName = surfaces.geometry().names()[geomi];

            const searchableSurface& s = surfaces.geometry()[geomi];

            label shellLevel = shells.levels()[shelli][0][0]; //assuming iso refinement

            triSurface triShell;

            // distribute shells across processors
            if (shelli % Pstream::nProcs() == Pstream::myProcNo())
            {
                if (isA<triSurfaceMesh>(s))
                {
                    const triSurfaceMesh& triMesh = refCast<const triSurfaceMesh>(s);

                    triShell =
                        const_cast<triSurface&>
                        (
                            refCast<const triSurface>(triMesh)
                        );
                }
                else if (isA<searchableBox>(s))
                {
                    treeBoundBox sBox(refCast<const searchableBox>(s));

                    triShell = boundBoxToTriSurface(sBox, shellLevel, /*decompose*/false);
                }
            }

            triShell.rotate(invR); //TODO from dict if in global coordSystem
            triShell.scalePoints(invScale);
            triShell.translatePoints(levelSetOffset);

            if (shells.modes()[shelli] == shellSurfaces::INSIDE)
            {
                std::vector<BoolGrid::Ptr> shellRefinementInterior(maxCellLevel_ + 1, nullptr);

                std::vector<FloatGrid::Ptr> shellRefinementCellLevels =
                    refinementRegionGrading
                    (
                        shellName,
                        triShell,
                        shellLevel,
                        nCellsBetweenLevels,
                        boundingBox,
                        shellRefinementInterior
                    );

                cellLevelGridsCollection.push_back(shellRefinementCellLevels);
                cellLevelInteriorCollection.push_back(shellRefinementInterior);
            }
            else if (shells.modes()[shelli] == shellSurfaces::DISTANCE)
            {
                //TODO decompose triSurface
                std::vector<BoolGrid::Ptr> shellRefinementInterior(maxCellLevel_ + 1, nullptr);

                std::vector<FloatGrid::Ptr> shellRefinementCellLevels =
                    distanceRefinementGrading
                    (
                        shellName,
                        triShell,
                        shells.levels()[shelli],
                        shells.shellIsoDistances(shelli),
                        nCellsBetweenLevels,
                        boundingBox,
                        shellRefinementInterior
                    );

                cellLevelGridsCollection.push_back(shellRefinementCellLevels);
                cellLevelInteriorCollection.push_back(shellRefinementInterior);
            }
        } // shellCellLevels
    } // Timer voxelize refinement regions


    //// find geometry used for both refinementSurfaces
    //// and refinementRegions (distance mode)
    //const wordList& allGeomNames = surfaces.geometry().names();

    //labelHashSet distRefinementGeom(surfaces.surfaces());
    //distRefinementGeom.retain(labelHashSet(shells.shells()));
    //Info<< "Used for both distance and surface refinement:" << endl;
    //for (const label geomi : distRefinementGeom)
    //{
    //    Info<< "    " << allGeomNames[geomi] << nl;
    //}


    // ### Voxelize surfaces ### //

    const labelList surfaceGeometry = surfaces.surfaces();

    List<triSurface> surfList(surfaceGeometry.size());

    bool distributeSurf = true;

    forAll(surfaceGeometry, surfi)
    {
        label geomi = surfaceGeometry[surfi];

        const searchableSurface& s = surfaces.geometry()[geomi];

        triSurface inputSurf;

        if (isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& triMesh = refCast<const triSurfaceMesh>(s);
            //triSurface& inputSurf = const_cast<triSurface&>
            inputSurf = const_cast<triSurface&>
            (
                refCast<const triSurface>(triMesh)
            );
        }

        inputSurf.rotate(invR);
        inputSurf.scalePoints(invScale);
        inputSurf.translatePoints(levelSetOffset);

        surfList[surfi] = inputSurf;

        // do not distribute if at least one surface is distributed<Tri|CAD>Surface
        if
        (
            distributeSurf
         && (
                isA<distributedTriSurfaceMesh>(s)
#ifdef FOAM_USE_OPENCASCADE
             || isA<distributedCADSurface>(s)
#endif
            )
        )
        {
            distributeSurf = false;
        }
    }

    // every processor has globalSurf
    triSurface globalSurf = addSurfaces(surfList);

    /////////////////////////////
    ////TODO exclude from globalSurf surfaces needed only for cellZones (e.g. radiators)
    //std::vector<BoolGrid::Ptr> SDFinteriorGrids(maxCellLevel_ + 1, nullptr);
    //std::vector<FloatGrid::Ptr> SDFGrids =
    //    //multiResSDFandDispl
    //    //(
    //    //    globalSurf,
    //    //    min(surfaces.minLevel()),
    //    //    max(surfaces.maxLevel()),
    //    //    SDFinteriorGrids,
    //    //    surfaceDisplacementGrids
    //    //);
    //    multiResSDF
    //    (
    //        globalSurf,
    //        min(surfaces.minLevel()),
    //        max(surfaces.maxLevel()),
    //        SDFinteriorGrids
    //    );

    //cellLevelInteriorCollection.push_back(SDFinteriorGrids);

    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        if (!isValid<Vec3dGrid>(surfaceDisplacementGrids[cellLevel]))
        {
            surfaceDisplacementGrids[cellLevel] = Vec3dGrid::create();
        }
    }

    /////////////////////////////

    triSurface myProcSurf;
    if (Pstream::parRun() && distributeSurf)
    {
        // every processor distributes but retains only its portion
        myProcSurf = distributeSurface(globalSurf, surfaces.maxLevel());
    }
    else
    {
        myProcSurf = globalSurf;
    }

    // separate surface regions
    List<triSurface> subSurfs = splitSurface(myProcSurf);

    // find min and max level of existing regions in this processor
    label myProcMinLevel = max(surfaces.minLevel());
    label myProcMaxLevel = 0;

    forAll(subSurfs, i)
    {
        const triSurface& s = subSurfs[i];

        // skip surface which is not present on this processor
        if (s.size() == 0) continue;

        myProcMinLevel =
            Foam::min
            (
                myProcMinLevel,
                surfaces.minLevel()[i]
            );

        myProcMaxLevel =
            Foam::max
            (
                myProcMaxLevel,
                surfaces.maxLevel()[i]
            );
    }

    // calculate Unsigned Distance Field of every patch
    { Timer timer("boundary Unsigned Distance Field and refinement grading");

        { Timer timer("boundary Unsigned Distance Field");

            forAll(subSurfs, i)
            {
                multiResUDF
                (
                    allRegNames[i],
                    subSurfs[i],
                    globalSurf,
                    i,
                    surfaces,
                    nCellsBetweenLevels,
                    interiorGrids,
                    bufferCellLevels,
                    patchCellLevels
                );
            }
        } //Timer boundary UDF

        // Calculate Level Sets of surface regions and grading
        //  - minLevel on surface, coarser levels for grading (nCellsBetweenLevels)
        //
        // Calculate Level Sets of feature edges (curvature) and grading
        //  - maxLevel on edges, coarser levels for grading (nCellsBetweenLevels)
        //
        {
            Timer timer("Refinement grading of surfaces");

            scalar featAngle = 0;

            if (refineParams.curvature() > 0)
            {
                featAngle = Foam::radToDeg(Foam::acos(refineParams.curvature())); //Foam::cos(degToRad("resolveFeatureAngle"));
            }

            scalarField featAngles(nRegions, featAngle);

            forAll(featAngles, surfI)
            {
                if (surfaces.featureRefineAngle(surfI, /*region*/0) > 0)
                {
                    featAngles[surfI] = Foam::radToDeg(surfaces.featureRefineAngle(surfI, /*region*/0));
                }
            }

            // group regions which have same minMax level and same curvature refinement
            labelListList sameMinMaxLevel =
                groupSameLevelRegions
                (
                    allRegNames,
                    surfaces,
                    featAngles
                );

            const List<triSurface> sameLevelSurfs =
                splitSurface(myProcSurf, sameMinMaxLevel);

            forAll(sameMinMaxLevel, samei)
            {
                const label& minSurfLevel = surfaces.minLevel()[sameMinMaxLevel[samei][0]];
                const label& maxSurfLevel = surfaces.maxLevel()[sameMinMaxLevel[samei][0]];
                const scalar& sameFeatAngle = featAngles[sameMinMaxLevel[samei][0]];

                const label& nPatches = sameMinMaxLevel[samei].size();

                Info<< "\n~~~~~~~~~~~~~~~~~~\nProcessing "
                    << nPatches
                    << " patch" << (nPatches > 1 ? "es" : "")
                    << " with level ( " << minSurfLevel
                    << " " << maxSurfLevel
                    << " ) and featureAngle "
                    << sameFeatAngle << " deg :"
                    << endl;

                forAll(sameMinMaxLevel[samei], patchi)
                {
                    Info<< "    "
                        << allRegNames[sameMinMaxLevel[samei][patchi]]
                        << endl;
                }

                word patchName = "";

                if (sameMinMaxLevel[samei].size() == 1)
                {
                    patchName = allRegNames[sameMinMaxLevel[samei][0]];
                }
                else
                {
                    patchName = Foam::name(nPatches)
                                + "_regions_with_level_("
                                + Foam::name(minSurfLevel)
                                + "_"
                                + Foam::name(maxSurfLevel)
                                + ")";
                }

                { //Timer timer(patchName + word(" base"));
                std::vector<BoolGrid::Ptr> patchAndBufferInterior(maxCellLevel_ + 1, nullptr);

                //patches and buffer cells cellLevels
                std::vector<FloatGrid::Ptr> patchAndBufferCellLevels =
                    surfaceGrading
                    (
                        patchName,
                        patchCellLevels,
                        bufferCellLevels,
                        minSurfLevel,
                        nCellsBetweenLevels,
                        boundingBox,
                        patchAndBufferInterior
                    );

                cellLevelGridsCollection.push_back(patchAndBufferCellLevels);
                cellLevelInteriorCollection.push_back(patchAndBufferInterior);

                // save innerGrids
                for (label cellLevel = minSurfLevel; cellLevel >= minSurfLevel-1; --cellLevel)
                {
                    BoolGrid::Ptr innerGrid = BoolGrid::create(false);

                    innerGrid->topologyUnion(*patchAndBufferCellLevels[cellLevel]);
                    openvdb::tools::dilateActiveValues(innerGrid->tree(), 3);

                    if (isValid<BoolGrid>(innerGrids[cellLevel]))
                    {
                        innerGrids[cellLevel]->topologyUnion(*innerGrid);
                    }
                    else
                    {
                        const word gridName = "innerGrid_level " + Foam::name(cellLevel);
                        innerGrid->setName(gridName);
                        innerGrids[cellLevel] = innerGrid;
                    }
                }
                } // Timer timer(patchName + word(" base"));

                //curvature and buffer cells CellLevels
                if (maxSurfLevel != minSurfLevel)
                {
                    //Timer timer(patchName + word(" curvature"));

                    word curvName = patchName + "_curv_" + std::to_string(int(ceil(sameFeatAngle))) +"_deg";

                    //TODO consider reading refinementFeatures instead?
                    triSurface featSurf =
                        featureTriSurface
                        (
                            sameLevelSurfs[samei],
                            sameFeatAngle,
                            curvName
                        );

                    std::vector<std::vector<FloatGrid::Ptr>> curvBufferCellLevels
                    (
                        nRegions,
                        std::vector<FloatGrid::Ptr>
                        (
                            maxCellLevel_ + 1,
                            nullptr
                        )
                    );

                    std::vector<std::vector<FloatGrid::Ptr>> curvCellLevels
                    (
                        nRegions,
                        std::vector<FloatGrid::Ptr>
                        (
                            maxCellLevel_ + 1,
                            nullptr
                        )
                    );

                    multiResUDF
                    (
                        curvName,
                        featSurf,
                        globalSurf,
                        /*regioni*/sameMinMaxLevel[samei][0],
                        surfaces,
                        nCellsBetweenLevels,
                        interiorGrids,
                        curvBufferCellLevels,
                        curvCellLevels,
                        /*curvature*/true
                    );

                    std::vector<BoolGrid::Ptr> curvatureAndBufferInterior(maxCellLevel_ + 1, nullptr);

                    std::vector<FloatGrid::Ptr> curvatureAndBufferCellLevels =
                        surfaceGrading
                        (
                            curvName,
                            curvCellLevels,
                            curvBufferCellLevels,
                            maxSurfLevel,
                            nCellsBetweenLevels,
                            boundingBox,
                            curvatureAndBufferInterior
                        );

                    cellLevelGridsCollection.push_back(curvatureAndBufferCellLevels);
                    cellLevelInteriorCollection.push_back(curvatureAndBufferInterior);

                    // save innerGrids
                    for (label cellLevel = maxSurfLevel; cellLevel >= minSurfLevel+1; --cellLevel)
                    {
                        BoolGrid::Ptr innerGrid = BoolGrid::create(false);

                        innerGrid->topologyUnion(*curvatureAndBufferCellLevels[cellLevel]);
                        openvdb::tools::dilateActiveValues(innerGrid->tree(), 3);

                        if (isValid<BoolGrid>(innerGrids[cellLevel]))
                        {
                            innerGrids[cellLevel]->topologyUnion(*innerGrid);
                        }
                        else
                        {
                            const word gridName = "innerGrid_level " + Foam::name(cellLevel);
                            innerGrid->setName(gridName);
                            innerGrids[cellLevel] = innerGrid;
                        }
                    }
                } // if minLevel != maxLevel
            } //forAll sameMinMaxLevel
        } // Timer timer("Refinement grading of surfaces");
    } //Timer timer("boundary Unsigned Distance Field and refinement grading");

    Info<< nl;

    // local topology union of all interiorGrids
    combineGrids<BoolGrid>
    (
        "interiorGrids",
        cellLevelInteriorCollection,
        interiorGrids
    );

    Info<< nl;

    // local topology union of all cellLevelGrids and remove interiorGrids
    combineGrids
    (
        "cellLevelGrids",
        cellLevelGridsCollection,
        cellLevelGrids,
        interiorGrids
    );


    // ### Synchronize voxels across nodes ### //

    { Timer timer("wait for other procs");
        List<labelList> wait(Pstream::nProcs());
        wait[Pstream::myProcNo()] = labelList(1, Pstream::myProcNo());
        allGatherList<labelList>(wait);
    }

    // List[proci][cellLevel]
    List<List<boundBox>> procBBoxes =
        sendReceiveBoundBoxes<FloatGrid>(cellLevelGrids);

    List<List<boundBox>> procInteriorBBoxes =
        sendReceiveBoundBoxes<BoolGrid>(interiorGrids);

    // convert List[proci][cellLevel] to List[cellLevel][proci]
    List<List<boundBox>> procsBB
    (
        maxCellLevel_ + 1,
        List<boundBox>
        (
            Pstream::nProcs(),
            boundBox()
        )
    );

    List<List<boundBox>> procsInteriorBB
    (
        maxCellLevel_ + 1,
        List<boundBox>
        (
            Pstream::nProcs(),
            boundBox()
        )
    );

    for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
    {
        for (label proci = 0; proci < Pstream::nProcs(); ++proci)
        {
            procsBB[cellLevel][proci] = procBBoxes[proci][cellLevel];

            procsInteriorBB[cellLevel][proci] = procInteriorBBoxes[proci][cellLevel];
        }
    }

    //if (Pstream::parRun()) //count voxels also in serial
    {
        // send my interiorGrid to procs that overlap with myProcNo
        if (Pstream::haveThreads()) // && procs > threads per proc
        {
            Timer timer("sync interiorGrids multithread");

            for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
            {
                PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
                PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

                boundBox myProcBB = procInteriorBBoxes[Pstream::myProcNo()][cellLevel];

                //TODO better with tasks?
                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, Pstream::nProcs()),
                    SendGridOp<BoolGrid>
                    (
                        interiorGrids[cellLevel],
                        pBufSize,
                        pBufGrid,
                        myProcBB, //interior
                        procsBB[cellLevel],
                        /*isCellLevelGrid*/false
                    )
                );

                pBufSize.finishedSends();
                pBufGrid.finishedSends();

                myProcBB = procBBoxes[Pstream::myProcNo()][cellLevel];

                std::vector<BoolGrid::Ptr> gridsFromProcs(Pstream::nProcs(), nullptr);

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, Pstream::nProcs()),
                    ReceiveGridsOp<BoolGrid>
                    (
                        gridsFromProcs,
                        pBufSize,
                        pBufGrid,
                        myProcBB,
                        procsInteriorBB[cellLevel],
                        /*isCellLevelGrid*/false
                    )
                );

                for (label i = 0; i < Pstream::nProcs(); ++i)
                {
                    if (gridsFromProcs[i])
                    {
                        cellLevelGrids[cellLevel]->topologyDifference(*gridsFromProcs[i]);
                    }
                }
                openvdb::tools::pruneInactive(cellLevelGrids[cellLevel]->tree());
            } //for cellLevel
        } //Timer sync interiorGrids multithreaded
        else
        {
            Timer timer("sync interiorGrids");
            for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
            {
                boundBox myProcBB = procInteriorBBoxes[Pstream::myProcNo()][cellLevel];

                PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
                PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

                // send to overlapping proc
                for (label proci = 0; proci < Pstream::nProcs(); ++proci)
                {
                    if (proci == Pstream::myProcNo()) continue;

                    const boundBox& otherProcBB = procBBoxes[proci][cellLevel];

                    if (myProcBB.overlaps(otherProcBB))
                    {
                        sendGrid<BoolGrid>
                        (
                            interiorGrids[cellLevel],
                            pBufSize,
                            pBufGrid,
                            proci
                        );
                    }
                } //for proci (send)

                pBufSize.finishedSends();
                pBufGrid.finishedSends();

                myProcBB = procBBoxes[Pstream::myProcNo()][cellLevel];

                //receive from overlapping procs
                for (label proci = 0; proci < Pstream::nProcs(); ++proci)
                {
                    if (proci == Pstream::myProcNo()) continue;

                    const boundBox& otherProcBB = procInteriorBBoxes[proci][cellLevel];

                    if (myProcBB.overlaps(otherProcBB))
                    {
                        BoolGrid::Ptr interiorProcI =
                            receiveGrid<BoolGrid>
                            (
                                pBufSize,
                                pBufGrid,
                                proci
                            );

                        cellLevelGrids[cellLevel]->topologyDifference(*interiorProcI);
                    }
                } //for proci (receive)
                openvdb::tools::pruneInactive(cellLevelGrids[cellLevel]->tree());
            }// for cellLevel (share interiorGrid)
        } //syncInteriorGrids serial

        //// update procBBoxes? as they might be smaller after removing interiorGrids
        //Info<<"GGG procBBoxes pre " << procBBoxes <<endl;
        //procBBoxes = sendAndReceiveBoundBoxes<FloatGrid>(cellLevelGrids);
        //Info<<"GGG procBBoxes post " << procBBoxes <<endl;

        // send my cellLevelGrid to procs of lower rank that overlap with myProcNo
        if (Pstream::haveThreads()) // && procs > threads per proc
        {
            Timer timer("sync cellLevelGrids multithread");

            for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
            {
                PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
                PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

                const boundBox& myProcBB = procBBoxes[Pstream::myProcNo()][cellLevel];

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, Pstream::nProcs()),
                    SendGridOp<FloatGrid>
                    (
                        cellLevelGrids[cellLevel],
                        pBufSize,
                        pBufGrid,
                        myProcBB,
                        procsBB[cellLevel],
                        /*isCellLevelGrid*/true
                    )
                );

                pBufSize.finishedSends();
                pBufGrid.finishedSends();

                std::vector<FloatGrid::Ptr> gridsFromProcs(Pstream::nProcs(), nullptr);

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, Pstream::nProcs()),
                    ReceiveGridsOp<FloatGrid>
                    (
                        gridsFromProcs,
                        pBufSize,
                        pBufGrid,
                        myProcBB,
                        procsBB[cellLevel],
                        /*isCellLevelGrid*/true
                    )
                );

                for (label i = 0; i < Pstream::nProcs(); ++i)
                {
                    if (gridsFromProcs[i])
                    {
                        cellLevelGrids[cellLevel]->topologyDifference(*gridsFromProcs[i]);
                    }
                }
                openvdb::tools::pruneInactive(cellLevelGrids[cellLevel]->tree());
            } //for cellLevel
        } //Timer sync cellLevelGrids multithreaded
        else
        {
            Timer timer("sync cellLevelGrids");
            for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
            {
                const boundBox& myProcBB = procBBoxes[Pstream::myProcNo()][cellLevel];

                PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
                PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

                for (label proci = 0; proci < Pstream::nProcs(); ++proci)
                {
                    if (proci == Pstream::myProcNo()) continue;

                    const boundBox& otherProcBB = procBBoxes[proci][cellLevel];

                    if (myProcBB.overlaps(otherProcBB) && proci > Pstream::myProcNo())
                    {
                        sendGrid<FloatGrid>
                        (
                            cellLevelGrids[cellLevel],
                            pBufSize,
                            pBufGrid,
                            proci
                        );
                    } //if overlaps
                } //for proci (send)

                //Start sending and receiving and block
                pBufSize.finishedSends();
                pBufGrid.finishedSends();

                //receive from overlapping procs
                for (label proci = 0; proci < Pstream::nProcs(); ++proci)
                {
                    if (proci == Pstream::myProcNo()) continue;

                    const boundBox& otherProcBB = procBBoxes[proci][cellLevel];

                    if (myProcBB.overlaps(otherProcBB) && proci < Pstream::myProcNo())
                    {
                        FloatGrid::Ptr cellLevelProcI =
                            receiveGrid<FloatGrid>
                            (
                                pBufSize,
                                pBufGrid,
                                proci
                            );

                        cellLevelGrids[cellLevel]->topologyDifference(*cellLevelProcI);
                    } //if overlaps
                } //for proci (receive)
                openvdb::tools::pruneInactive(cellLevelGrids[cellLevel]->tree());
            }// for cellLevel (share interiorGrid)
        } //sync cellLevelGrids serial
    } // if parRun synchronize active voxels across processors

    //done also during grid decomposition for each procCellLevelGrids
    //if (false)
    {
        trimGrids<FloatGrid>(cellLevelGrids);
    }


    // ### Count active voxels and assign ID ### //
    std::vector<IndexGrid::Ptr>  globalIDGrids(maxCellLevel_ + 1, nullptr);

    pointField voxelCentres;

    labelList nVoxelsStart =
        countVoxels //GGG
        //countVoxelsLeaf //fix free(): invalid pointer on cirrus!
        (
            cellLevelGrids,
            globalIDGrids,
            voxelCentres
        );


    // ### Balance voxels across nodes ### //
    if (Pstream::parRun())
    {
        // decomposedCellLevelGrids[cellLevel][proci]
        std::vector<std::vector<FloatGrid::Ptr>> decomposedCellLevelGrids;

        // only finds which voxel goes to which processor
        // actual distribution of voxels is done later

        //nVoxelsForProc[cellLevel][procI];
        labelListList nVoxelsForProc =
            decomposeGrids
            (
                cellLevelGrids,
                globalIDGrids,
                decomposedCellLevelGrids,
                decompositionDict,
                nVoxelsStart,
                voxelCentres
            );

        //make it available to other procs
        //procVoxelsInProcNodes[myNode][cellLevel][toProc]
        List<labelListList> procVoxelsForProc(Pstream::nProcs());
        procVoxelsForProc[Pstream::myProcNo()] = nVoxelsForProc;

        { Timer timer("allGather procVoxelsInProcNodes");
            allGatherList<labelListList>(procVoxelsForProc);
        }


        // distribute decomposition
        // replace content of cellLevelGrids
        {
            Timer timer("distribute decomposition serial");

            for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
            {
                // add first voxels in myProc
                cellLevelGrids[cellLevel] = decomposedCellLevelGrids[cellLevel][Pstream::myProcNo()];

                PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
                PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

                // myProcNo sends to proci the processor voxels belonging to proci
                for (label proci = 0; proci < Pstream::nProcs(); ++proci)
                {
                    if (proci == Pstream::myProcNo()) continue;

                    if (nVoxelsForProc[cellLevel][proci] > 0)
                    {
                        sendGrid<FloatGrid>
                        (
                            decomposedCellLevelGrids[cellLevel][proci],
                            pBufSize,
                            pBufGrid,
                            proci
                        );
                    }
                } // for proci (send)

                pBufSize.finishedSends();
                pBufGrid.finishedSends();

                //receive from proci all voxels belonging to myProcNo
                for (label proci = 0; proci < Pstream::nProcs(); ++proci)
                {
                    if (proci == Pstream::myProcNo()) continue;

                    if (procVoxelsForProc[proci][cellLevel][Pstream::myProcNo()] > 0)
                    {
                        FloatGrid::Ptr grid =
                            receiveGrid<FloatGrid>
                            (
                                pBufSize,
                                pBufGrid,
                                proci
                            );

                        cellLevelGrids[cellLevel]->topologyUnion(*grid);
                    }
                } // for proci (receive)
            } // for cellLevel
        } // distribute decomposition serial

        //if (false)
        { Timer timer("trim decomposed grids");
            trimGrids<FloatGrid>(cellLevelGrids);
        } //Timer

        // List[proci][cellLevel]
        procBBoxes =
            sendReceiveBoundBoxes<FloatGrid>(cellLevelGrids);
    }

    // Refine cells to remove dangling cells
    //if (false)
    {
        refineDanglingVoxels(cellLevelGrids);

        trimGrids<FloatGrid>(cellLevelGrids);
    }

    // Add boundary cells which have 5 active neighbours (fill holes)
    {
        Timer timer("fill voxel holes");

        const label minSurfLevel = min(surfaces.minLevel());
        const label maxSurfLevel = max(surfaces.maxLevel());

        for (label cellLevel = maxSurfLevel; cellLevel >= minSurfLevel; --cellLevel)
        {
            openvdb::Coord nijk;
            label nInternal;

            FloatGrid::Accessor acc = cellLevelGrids[cellLevel]->getAccessor();

            //TODO multithreded
            for (auto iter = cellLevelGrids[cellLevel]->cbeginValueOff(); iter; ++iter)
            {
                const openvdb::Coord& ijk = iter.getCoord();

                nInternal = 0;

                for (label i = 0; i < 6; ++i)
                {
                    nijk = ijk + COORD_OFFSETS[i];

                    bool isOn = acc.isValueOn(nijk);

                    if (isOn)
                    {
                        nInternal++;
                    }
                }

                if (nInternal >= 5)
                {
                    acc.setValueOn(ijk, 1.0);
                }
            }
        }
    }

    // Remove cells with 2 or less internal faces
    if (false)
    {
        Timer timer("remove twoInternalFaces");

        openvdb::Coord ijk, nijk, nnijk;
        label nInternal;

        for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
        {
            FloatGrid::ConstAccessor acc = cellLevelGrids[cellLevel]->getConstAccessor();
            label one = (cellLevel < maxCellLevel_ ? 1 : 0 );
            FloatGrid::ConstAccessor fineAcc = cellLevelGrids[cellLevel+one]->getConstAccessor();

            label mOne = (cellLevel > 0 ? -1 : 0 );
            FloatGrid::ConstAccessor coarseAcc = cellLevelGrids[cellLevel+mOne]->getConstAccessor();

            label it = 0;
            //label nOn;
            //do
            //{
            //    nOn = 0;
            //    //TODO multithreaded
            //    for (auto iter = cellLevelGrids[cellLevel]->beginValueOff(); iter; ++iter)
            //    {
            //        const openvdb::Coord& ijk = iter.getCoord();

            //        nInternal = 0;

            //        for (label i = 0; i < 6; ++i)
            //        {
            //            nijk = ijk + COORD_OFFSETS[i];

            //            bool isOn = acc.isValueOn(nijk);

            //            if (isOn)
            //            {
            //                nInternal++;
            //            }

            //            //check if split face of coarser voxels are neighbour
            //            if (!isOn && cellLevel > 0)
            //            {
            //                bool coarseVoxelOn = false;

            //                switch ( (nijk[0] & 1) | ((nijk[1] & 1) << 1) | ((nijk[2] & 1) << 2) )
            //                {
            //                    case 0:// all even
            //                        coarseVoxelOn = coarseAcc.isValueOn(nijk>>1);
            //                        break;
            //                    case 1:// x is odd
            //                        coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(-1,0,0)>>1);
            //                        break;
            //                    case 2:// y is odd
            //                        coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(0,-1,0)>>1);
            //                        break;
            //                    case 3:// x&y are odd
            //                        coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(-1,-1,0)>>1);
            //                        break;
            //                    case 4:// z is odd
            //                        coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(0,0,-1)>>1);
            //                        break;
            //                    case 5:// x&z are odd
            //                        coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(-1,0,-1)>>1);
            //                        break;
            //                    case 6:// y&z are odd
            //                        coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(0,-1,-1)>>1);
            //                        break;
            //                    case 7:// all are odd
            //                        coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(-1,-1,-1)>>1);
            //                        break;
            //                } //switch

            //                if (coarseVoxelOn)
            //                {
            //                    nInternal++;
            //                    isOn = true;
            //                }
            //            }

            //            if (!isOn && cellLevel < maxCellLevel_)
            //            {
            //                for (label j = 0; j < 4; j++)
            //                {
            //                    nnijk = (ijk << 1) + COORD_OFFSETS_SPLIT[4*i + j]; //only face adjacent

            //                    if (fineAcc.isValueOn(nnijk))
            //                    {
            //                        nInternal++;
            //                        break;
            //                    }
            //                }
            //            }
            //        }

            //        if (nInternal >= 5)
            //        {
            //            iter.setValue(1.0);
            //            nOn++;
            //        }
            //    }
            //    ++it;
            //}
            //while ( nOn > 0 && it < 20);
            //Info<<"GGG nOn " << nOn
            //<<" it " << it
            //<<endl;

            label nOff;
            it = 0;
            do
            {
                nOff = 0;
                //TODO multithreaded
                for (auto iter = cellLevelGrids[cellLevel]->beginValueOn(); iter; ++iter)
                {
                    const openvdb::Coord& ijk = iter.getCoord();

                    nInternal = 0;

                    for (label i = 0; i < 6; ++i)
                    {
                        nijk = ijk + COORD_OFFSETS[i];

                        bool isOn = acc.isValueOn(nijk);

                        if (isOn)
                        {
                            nInternal++;
                        }

                        //check if split face of coarser voxels are neighbour
                        if (!isOn && cellLevel > 0)
                        {
                            bool coarseVoxelOn = false;

                            switch ( (nijk[0] & 1) | ((nijk[1] & 1) << 1) | ((nijk[2] & 1) << 2) )
                            {
                                case 0:// all even
                                    coarseVoxelOn = coarseAcc.isValueOn(nijk>>1);
                                    break;
                                case 1:// x is odd
                                    coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(-1,0,0)>>1);
                                    break;
                                case 2:// y is odd
                                    coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(0,-1,0)>>1);
                                    break;
                                case 3:// x&y are odd
                                    coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(-1,-1,0)>>1);
                                    break;
                                case 4:// z is odd
                                    coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(0,0,-1)>>1);
                                    break;
                                case 5:// x&z are odd
                                    coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(-1,0,-1)>>1);
                                    break;
                                case 6:// y&z are odd
                                    coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(0,-1,-1)>>1);
                                    break;
                                case 7:// all are odd
                                    coarseVoxelOn = coarseAcc.isValueOn(nijk.offsetBy(-1,-1,-1)>>1);
                                    break;
                            } //switch

                            if (coarseVoxelOn)
                            {
                                nInternal++;
                                isOn = true;
                            }
                        }

                        if (!isOn && cellLevel < maxCellLevel_)
                        {
                            for (label j = 0; j < 4; j++)
                            //for (label j = 0; j < 8; j++)
                            {
                                nnijk = (ijk << 1) + COORD_OFFSETS_SPLIT[4*i + j]; //only face adjacent
                                //nnijk = (nijk << 1) + COARSE_TO_FINE[j]; //check all 8 fine neighbour voxels

                                if (fineAcc.isValueOn(nnijk))
                                {
                                    nInternal++;
                                    break;
                                }
                            }
                        }
                    }

                    if (nInternal <= 2)
                    {
                        iter.setValueOff();
                        nOff++;
                    }
                }
                ++it;
            }
            while (nOff > 0 && it < 20);
        Info<<"GGG nOff " << nOff
            <<" it " << it
            <<endl;
        }
    } // remove twoInternalFaces


    // ### subset grids around car for snapping ### //

    const bool subsetInnerGrid(vdbDict.lookupOrDefault<bool>("subsetInnerGrid", false));
    const bool subsetOuterGrid(vdbDict.lookupOrDefault<bool>("subsetOuterGrid", false));

    {
        Timer t("subset <inner|outer>Grids");

        for (label cellLevel = maxCellLevel_; cellLevel >= 0; --cellLevel)
        {
            if (!isValid<BoolGrid>(innerGrids[cellLevel]))
            {
                BoolGrid::Ptr innerGrid = BoolGrid::create(false);
                const word gridName = "innerGrid_level " + Foam::name(cellLevel);
                innerGrid->setName(gridName);
                innerGrids[cellLevel] = innerGrid;
            }
        }

        // send my innerGrids to all procs for subsetting
        if (Pstream::haveThreads())
        {
            Timer timer("    sync innerGrids multithread");

            for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
            {
                PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
                PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, Pstream::nProcs()),
                    SendGridOp<BoolGrid>
                    (
                        innerGrids[cellLevel],
                        pBufSize,
                        pBufGrid,
                        boundBox(point(0), point(1)), //dummyBB
                        List<boundBox>(Pstream::nProcs(), boundBox(point(0), point(1))),
                        /*isCellLevelGrid*/false
                    )
                );

                pBufSize.finishedSends();
                pBufGrid.finishedSends();

                std::vector<BoolGrid::Ptr> gridsFromProcs(Pstream::nProcs(), nullptr);

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, Pstream::nProcs()),
                    ReceiveGridsOp<BoolGrid>
                    (
                        gridsFromProcs,
                        pBufSize,
                        pBufGrid,
                        boundBox(point(0), point(1)), //dummyBB
                        List<boundBox>(Pstream::nProcs(), boundBox(point(0), point(1))),
                        /*isCellLevelGrid*/false
                    )
                );

                for (label i = 0; i < Pstream::nProcs(); ++i)
                {
                    if (gridsFromProcs[i])
                    {
                        innerGrids[cellLevel]->topologyUnion(*gridsFromProcs[i]);
                    }
                }
                openvdb::tools::pruneInactive(innerGrids[cellLevel]->tree());
            } //for cellLevel
        } //Timer sync innerGrids multithreaded

        //subset cellLevelGrids
        for (label cellLevel = maxCellLevel_; cellLevel >= 0; --cellLevel)
        {
            BoolGrid::Ptr outerGrid = BoolGrid::create(false);
            const word gridName = "outerGrid_level " + Foam::name(cellLevel);
            outerGrid->setName(gridName);
            outerGrid->topologyUnion(*cellLevelGrids[cellLevel]);

            if (isValid<BoolGrid>(innerGrids[cellLevel]))
            {
                innerGrids[cellLevel]->topologyIntersection(*cellLevelGrids[cellLevel]);

                outerGrid->topologyDifference(*innerGrids[cellLevel]);
                openvdb::tools::pruneInactive(outerGrid->tree());
            }

            outerGrids[cellLevel] = outerGrid;

            if (subsetInnerGrid)
            {
                cellLevelGrids[cellLevel]->topologyIntersection(*innerGrids[cellLevel]);
            }
            else if (subsetOuterGrid)
            {
                cellLevelGrids[cellLevel]->topologyIntersection(*outerGrids[cellLevel]);
            }
        }
    } //subset innerGrids


    // write vdb grids to file
    const bool writeGrids(vdbDict.lookupOrDefault<bool>("writeGrids", false));

    if (writeGrids)
    {
        Timer timer("Write VDB grids");

        const word procDir
        (
            Pstream::parRun()
         ?  "processor" + Foam::name(Pstream::myProcNo())
         :  "./"
        );

        openvdb::io::File file(procDir + "/cellLevelGrids.vdb");

        // Add the grid pointer to a container.
        openvdb::GridPtrVec grids;

        for (label cellLevel = 0; cellLevel <= maxCellLevel_; ++cellLevel)
        {
            grids.push_back(cellLevelGrids[cellLevel]);
        }

        // Write out the contents of the container.
        file.write(grids);
        file.close();
    } //Timer timer("Write VDB grids");


    //hexRef8 data
    labelList hexRef8cellLevel;
    labelList hexRef8pointLevel;

    // ### Convert voxels to OpenFOAM mesh ### //
    vdbGridsToPolyMesh
    (
        cellLevelGrids,
        surfaceDisplacementGrids,
        innerGrids,
        outerGrids,
        /*nRegions*/ mesh.boundaryMesh().names().size(), // also include "*_slave"
        levelSetOffset,
        scale,
        R,
        boundingBox, //level0 bounds
        vdbDict,

        mesh,
        hexRef8cellLevel,
        hexRef8pointLevel
    );

    meshRefiner.meshCutter().updateLevels
    (
        hexRef8pointLevel,
        hexRef8cellLevel
    );
    meshRefiner.meshCutter().updateLevel0EdgeLength
    (
        level0Edge
    );

    //if (debug)
    {
        Timer t("mesh.write()");
        mesh.write();

        meshRefiner.meshCutter().syncLevel();
        meshRefiner.meshCutter().write();

    //    bool report = true;
    //    mesh.checkPoints(report);
    //    mesh.checkUpperTriangular(report);
    //    mesh.checkCellsZipUp(report);
    //    mesh.checkFaceVertices(report);
    //    mesh.checkClosedBoundary(report);
    //    mesh.checkClosedCells(report);
    //    mesh.checkFaceAreas(report);
    //    mesh.checkCellVolumes(report);
    //    mesh.checkFaceOrthogonality(report);
    //    mesh.checkFacePyramids(report);
    }
} //voxelize

