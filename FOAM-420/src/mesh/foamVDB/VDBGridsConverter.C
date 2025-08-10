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
    (c) 2020 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "VDBGridsConverter.H"
#include "fields/Fields/DynamicField/DynamicField.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"
#include "memInfo/memInfo.H"

//clash with macro expansion in <openvdb_dir>/include/openvdb/math/Vec3.h:657:16
#undef Log
#include <openvdb/tools/ValueTransformer.h> // for foreach()

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const std::vector<openvdb::Coord> Foam::VDBGridsConverter::getVoxelVertices
(
    const openvdb::Coord& coord,
    const label pointsPerEdge /* = 2 */
)
{
    label xmin = pointsPerEdge * (coord.x());
    label xmax = pointsPerEdge * (coord.x() + 1);
    label ymin = pointsPerEdge * (coord.y());
    label ymax = pointsPerEdge * (coord.y() + 1);
    label zmin = pointsPerEdge * (coord.z());
    label zmax = pointsPerEdge * (coord.z() + 1);

    std::vector<openvdb::Coord> voxelVertices(8);

    voxelVertices[0] = openvdb::Coord(xmin, ymin, zmin);
    voxelVertices[1] = openvdb::Coord(xmax, ymin, zmin);
    voxelVertices[2] = openvdb::Coord(xmax, ymax, zmin);
    voxelVertices[3] = openvdb::Coord(xmin, ymax, zmin);
    voxelVertices[4] = openvdb::Coord(xmin, ymin, zmax);
    voxelVertices[5] = openvdb::Coord(xmax, ymin, zmax);
    voxelVertices[6] = openvdb::Coord(xmax, ymax, zmax);
    voxelVertices[7] = openvdb::Coord(xmin, ymax, zmax);

    return voxelVertices;
}


const std::vector<openvdb::Coord> Foam::VDBGridsConverter::getVoxelFaceCentres
(
    const openvdb::Coord& coord,
    const label pointsPerEdge /* = 2 */
)
{
    // in index space the face centres are stored in a grid with double the size of vertex grid
    label xmin = pointsPerEdge_ * (coord.x());
    label ymin = pointsPerEdge_ * (coord.y());
    label zmin = pointsPerEdge_ * (coord.z());

    std::vector<openvdb::Coord> voxelFaceCentres(6);

    voxelFaceCentres[0] = openvdb::Coord(2*xmin - pointsPerEdge_, 2*ymin, 2*zmin).offset(pointsPerEdge_); // 2*x-min
    voxelFaceCentres[1] = openvdb::Coord(2*xmin + pointsPerEdge_, 2*ymin, 2*zmin).offset(pointsPerEdge_); // 2*x-max
    voxelFaceCentres[2] = openvdb::Coord(2*xmin, 2*ymin - pointsPerEdge_, 2*zmin).offset(pointsPerEdge_); // 2*y-min
    voxelFaceCentres[3] = openvdb::Coord(2*xmin, 2*ymin + pointsPerEdge_, 2*zmin).offset(pointsPerEdge_); // 2*y-max
    voxelFaceCentres[4] = openvdb::Coord(2*xmin, 2*ymin, 2*zmin - pointsPerEdge_).offset(pointsPerEdge_); // 2*z-min
    voxelFaceCentres[5] = openvdb::Coord(2*xmin, 2*ymin, 2*zmin + pointsPerEdge_).offset(pointsPerEdge_); // 2*z-max

    return voxelFaceCentres;
}


template<typename GridT>
Foam::VDBGridsConverter::InterfaceFace Foam::VDBGridsConverter::isCellLevelInterfaceFace
(
    const face& face,
    openvdb::math::Transform::ConstPtr coarseTransform,
    typename GridT::ConstAccessor& cCoarseFaceGridAccessor
)
{
   for (unsigned i = 0; i < 4; i++)
   {
       const openvdb::Coord faceVertexGlobal
       (
           xyzPointsList[face[i]].x(),
           xyzPointsList[face[i]].y(),
           xyzPointsList[face[i]].z()
       );

       const openvdb::Coord faceVertexCoarseLocal =
           hVDB_.toLocal
           (
               faceVertexGlobal,
               coarseTransform
           );

       openvdb::Coord xyzFaceCoarseLocal =
           openvdb::Coord
           (
               2. * faceVertexCoarseLocal.x(),
               2. * faceVertexCoarseLocal.y(),
               2. * faceVertexCoarseLocal.z()
           );

       if (cCoarseFaceGridAccessor.isValueOn(xyzFaceCoarseLocal))
       {
           //if (debug)
           //{
           //    std::cout<< "faceVertexGlobal " << faceVertexGlobal
           //        << " is at face centre of cellLevel " << cellLevel - 1
           //        <<std::endl;
           //}

           InterfaceFace iFace = {true, xyzFaceCoarseLocal};
           return iFace;
       }
   }

   openvdb::Coord xyz;
   InterfaceFace iFace = {false, xyz};
   return iFace;
}


void Foam::VDBGridsConverter::convert()
{
    memInfo mem;
    //Info<< "Initial memory " << mem.update().size() << " kB" << endl;

    struct PointCounterOp
    {
        //using AccessorT = openvdb::tree::ValueAccessorRW<IndexGrid::TreeType, true>;
        //not threaded but much faster than using thread-safe ValueAccessorRW
        using AccessorT = openvdb::tree::ValueAccessor<IndexGrid::TreeType, true>;

        AccessorT& accessor_;

        std::atomic<label> * const pointCount_;

        PointCounterOp
        (
            AccessorT& accessor,
            std::atomic<label>& pointCount

        )
        :
            accessor_(accessor),
            pointCount_(&pointCount)
        {}

        void operator()(const FloatGrid::ValueOnCIter& iter) const
        {
            // Note: this code must be thread-safe.
            const std::vector<openvdb::Coord> voxelVertices =
                getVoxelVertices(iter.getCoord());

            for (unsigned vertexI=0; vertexI < voxelVertices.size(); vertexI++)
            {
                const openvdb::Coord& xyz = voxelVertices[vertexI];

                if (accessor_.getValue(xyz) < 0)
                {
                    accessor_.setValue
                    (
                        xyz,
                        pointCount_->fetch_add(1)
                    );
                }
            } // for voxelVertex
        } //operator ()
    }; // PointCounterOp

    struct UpdatePointListOp
    {
        using TransformPtr = openvdb::math::Transform::Ptr;

        TransformPtr indexTransform_;
        pointList * const xyzPointsList_;

        UpdatePointListOp
        (
            TransformPtr indexTransform,
            pointList& xyzPointsList
        )
        :
            indexTransform_(indexTransform),
            xyzPointsList_(&xyzPointsList)
        {}

        inline void operator()(const IndexGrid::ValueOnCIter& iter) const
        {
            const openvdb::Coord& xyz = iter.getCoord();

            const openvdb::Coord xyzGlobal = foamVDB::toGlobal(xyz, indexTransform_);

            (*xyzPointsList_)[iter.getValue()] =
                (
                    point
                    (
                        xyzGlobal.x(),
                        xyzGlobal.y(),
                        xyzGlobal.z()
                    )
                );
        }
    }; // UpdatePointListOp

    struct UpdateFinerIndexGridOp
    {
        using AccessorI = openvdb::tree::ValueAccessorRW<IndexGrid::TreeType, true>;
        using cAccessorI = openvdb::tree::ValueAccessorRW<const IndexGrid::TreeType, true>;

        cAccessorI& coarseIndexGridAccessor_;
        AccessorI& fineIndexGridAccessor_;
        IndexGrid& coarseIndexGrid_;
        IndexGrid& fineIndexGrid_;
        pointList * const xyzPointsList_;

        UpdateFinerIndexGridOp
        (
            cAccessorI& coarseIndexGridAccessor,
            AccessorI& fineIndexGridAccessor,
            IndexGrid& coarseIndexGrid,
            IndexGrid& fineIndexGrid,
            pointList& xyzPointsList
        )
        :
            coarseIndexGridAccessor_(coarseIndexGridAccessor),
            fineIndexGridAccessor_(fineIndexGridAccessor),
            coarseIndexGrid_(coarseIndexGrid),
            fineIndexGrid_(fineIndexGrid),
            xyzPointsList_(&xyzPointsList)
        {}

        inline void operator()(const FaceGrid::ValueOnCIter& iter) const
        {
            //const openvdb::Coord coarseFaceCentreLocal = iter.getCoord();

            const openvdb::math::Vec4ui& f = iter.getValue();

            std::vector<openvdb::Coord> cellLevelInterfacePoints(4);

            for (unsigned n=0; n < 4; n++)
            {
                cellLevelInterfacePoints[n] =
                    openvdb::Coord
                    (
                        (*xyzPointsList_)[f[n]].x(),
                        (*xyzPointsList_)[f[n]].y(),
                        (*xyzPointsList_)[f[n]].z()
                    );
            }

            for (unsigned pointI=0; pointI < 4; pointI++)
            {
                const openvdb::Coord& xyzGlobal = cellLevelInterfacePoints[pointI];

                const openvdb::Coord xyzCoord = foamVDB::toLocal(xyzGlobal, coarseIndexGrid_.transformPtr());

                const openvdb::Coord xyzFineCoord = foamVDB::toLocal(xyzGlobal, fineIndexGrid_.transformPtr());

                const label pointIndex = coarseIndexGridAccessor_.getValue(xyzCoord);

                fineIndexGridAccessor_.setValue(xyzFineCoord, pointIndex);
            }
        } //operator ()
    }; // UpdateFinerIndexGridOp

    const label& maxCellLevel = hVDB_.maxCellLevel();

    std::vector<IndexGrid::Ptr> indexGrids;
    std::vector<IndexGrid::Ptr> ownerGrids;
    std::vector<IndexGrid::Ptr> patchIDGrids;
    std::vector<FaceGrid::Ptr>  faceGrids;
    std::vector<BoolGrid::Ptr>  procBoundaryGrids;
    std::vector<BoolGrid::Ptr>  cellLevelInterfaceFaceGrids;
    std::vector<BoolGrid::Ptr>  refinedBoundaryFaceGrids;

    // initialize indexGrids
    for (label cellLevel = 0; cellLevel <= maxCellLevel; cellLevel++)
    {
        label edgeLength = Foam::pow(2, maxCellLevel - cellLevel); //in (smallest) voxel units

        //if (debug)
        //{
        //    Info<< "edgeLength (in voxel units) "
        //        << edgeLength
        //        << " - cellLevel "
        //        << cellLevel << endl;
        //}

        openvdb::math::Transform::Ptr levelNTransform =
            openvdb::math::Transform::createLinearTransform
            (
                edgeLength
            );

        IndexGrid::Ptr levelNIndexGrid = IndexGrid::create(-1);
        levelNIndexGrid->setTransform(levelNTransform);
        indexGrids.push_back(levelNIndexGrid);

        // the index space of the following grids is doubled as values are stored at face centre
        IndexGrid::Ptr levelNOwnerGrid = IndexGrid::create(-1);
        levelNOwnerGrid->setTransform(levelNTransform);
        ownerGrids.push_back(levelNOwnerGrid);

        IndexGrid::Ptr levelNPatchIDGrid = IndexGrid::create(-1);
        levelNPatchIDGrid->setTransform(levelNTransform);
        patchIDGrids.push_back(levelNPatchIDGrid);

        FaceGrid::Ptr levelNFaceGrid = FaceGrid::create();
        levelNFaceGrid->setTransform(levelNTransform);
        faceGrids.push_back(levelNFaceGrid);

        BoolGrid::Ptr levelNCellLevelInterfaceFaceGrid = BoolGrid::create();
        levelNCellLevelInterfaceFaceGrid->setTransform(levelNTransform);
        cellLevelInterfaceFaceGrids.push_back(levelNCellLevelInterfaceFaceGrid);

        BoolGrid::Ptr levelNRefinedBoundaryFaceGrid = BoolGrid::create();
        levelNRefinedBoundaryFaceGrid->setTransform(levelNTransform);
        refinedBoundaryFaceGrids.push_back(levelNRefinedBoundaryFaceGrid);
    } // initialize indexGrids

    // initialize total number of points and cells
    std::atomic<label> pointCount{0};
    openvdb::Int32 cellCount = 0;

    DynamicList<face> faces;
    DynamicList<label> owner;
    DynamicList<label> neighbour;

    //fields
    DynamicField<scalar> cellLevelField;
    DynamicField<scalar> distanceField;
    DynamicField<scalar> patchIDField;
    DynamicField<scalar> procIDField;

    for (label cellLevel=0; cellLevel <= maxCellLevel; cellLevel++)
    {
        Timer timer("all cellLevel " + std::to_string(cellLevel));
        //FloatGrid::Ptr levelNGrid = *it;

        label coarseLevel = (cellLevel > 0) ? cellLevel - 1 : 0;

        // Get an accessor for coordinate-based access to voxels.
        IndexGrid::Accessor       levelNIndexGridAccessor = indexGrids[cellLevel]->getAccessor();
        IndexGrid::ConstAccessor cLevelNIndexGridAccessor = indexGrids[cellLevel]->getConstAccessor();

        IndexGrid::Accessor       levelNOwnerGridAccessor = ownerGrids[cellLevel]->getAccessor();
        IndexGrid::ConstAccessor cLevelNOwnerGridAccessor = ownerGrids[cellLevel]->getConstAccessor();
        IndexGrid::ConstAccessor cCoarseOwnerGridAccessor = ownerGrids[coarseLevel]->getConstAccessor();

        IndexGrid::Accessor       levelNPatchIDGridAccessor = patchIDGrids[cellLevel]->getAccessor();
        IndexGrid::ConstAccessor cLevelNPatchIDGridAccessor = patchIDGrids[cellLevel]->getConstAccessor();
        IndexGrid::ConstAccessor cCoarsePatchIDGridAccessor = patchIDGrids[coarseLevel]->getConstAccessor();

        FaceGrid::Accessor       levelNFaceGridAccessor = faceGrids[cellLevel]->getAccessor();
        FaceGrid::Accessor       coarseFaceGridAccessor = faceGrids[coarseLevel]->getAccessor();
        FaceGrid::ConstAccessor cLevelNFaceGridAccessor = faceGrids[cellLevel]->getConstAccessor();
        FaceGrid::ConstAccessor cCoarseFaceGridAccessor = faceGrids[coarseLevel]->getConstAccessor();

        BoolGrid::Accessor       levelNCellLevelInterfaceFaceGridAccessor = cellLevelInterfaceFaceGrids[cellLevel]->getAccessor();
        BoolGrid::ConstAccessor cLevelNCellLevelInterfaceFaceGridAccessor = cellLevelInterfaceFaceGrids[cellLevel]->getConstAccessor();

        BoolGrid::Accessor       coarseRefinedBoundaryFaceGridAccessor = refinedBoundaryFaceGrids[coarseLevel]->getAccessor();

        openvdb::math::Transform::Ptr indexTransform = indexGrids[cellLevel]->transformPtr();

        IndexGrid::Ptr levelNAddedIndexesGrid = IndexGrid::create(-1);
        levelNAddedIndexesGrid->setTransform(indexTransform);
        IndexGrid::Accessor levelNAddedIndexesGridAccessor = levelNAddedIndexesGrid->getAccessor();

        label nPointsCoarse = xyzPointsList.size();
        label cellLevelActiveVoxels = 0;

        for (size_t id = 0; id < allCellsGrids_.size(); id++)
        {
            FloatGrid::Ptr levelNGrid = allCellsGrids_[id][cellLevel];

            if (!levelNGrid) continue;
            if (levelNGrid->activeVoxelCount() == 0) continue;

            label levelNvoxelCount = 0;

            label patchID = id - 1;

            //Timer timer("patchID " + std::to_string(patchID) + " cellLevel " + std::to_string(cellLevel));
            //////////////

            levelNGrid->tree().voxelizeActiveTiles();

            std::string name = levelNGrid->getName();
            int n = name.length();
            char gridName[n+1];
            std::strcpy(gridName, name.c_str());

            label activeVoxels = levelNGrid->activeVoxelCount();
            label checkInterval = activeVoxels > 100000 ? 100000 : 1000;
            cellLevelActiveVoxels += activeVoxels;

            // populate xyzPointsList with coordinate of voxels
            // and add index of point to indexGrid
            {
                //Timer timer("create cells cellLevel " + std::to_string(cellLevel));

                //populate indexGrid and update pointCount
                //TODO
                { //Timer timer("PointCounterOp");
    ///////////////////////////                /////////////////////////
                //openvdb::tools::transformValues
                //(
                //    levelNGrid->cbeginValueOn(),
                //    *indexGrids[cellLevel],
                //    VoxelCounter::op,
                //    /*threaded*/false,//true,
                //    /*shareOp*/true
                //    /*merge=MERGE_ACTIVE_STATES*/
                //);
                //std::cout<<"GGG pointCount "<< pointCount <<std::endl;
    ///////////////////////////                /////////////////////////
                //openvdb::tree::ValueAccessorRW<IndexGrid::TreeType, true> acc(indexGrids[cellLevel]->tree());
                openvdb::tree::ValueAccessor<IndexGrid::TreeType, true> indexAcc(indexGrids[cellLevel]->tree());
                //FloatGrid::ConstAccessor gridAcc(levelNGrid->getConstAccessor());
                openvdb::tools::foreach
                (
                    levelNGrid->cbeginValueOn(),
                    PointCounterOp(indexAcc, pointCount),
                    /*threaded*/false,//true, // not threaded but much faster than using thread-safe ValueAccessorRW
                    /*shareOp*/false//true
                );
                } //Timer timer("PointCounterOp");

                ///////////////
                //openvdb::tree::ValueAccessorRW<IndexGrid::TreeType, true> acc(indexGrids[cellLevel]->tree());
                //PointCounter pointCounter(boundingBox, cellLevel, acc);
                //PointCounter::IterRange range(levelNGrid->cbeginValueOn());
                //tbb::parallel_for(range, pointCounter);

                //std::cout<<"GGG activeVoxels "<<levelNGrid->activeVoxelCount() <<" - " << activeVoxels<<std::endl;
                //std::cout<<"GGG pointCount "<< pointCount <<std::endl;

                xyzPointsList.resize(pointCount);
                //dump points to pointList
                openvdb::tools::foreach
                (
                    indexGrids[cellLevel]->cbeginValueOn(),
                    UpdatePointListOp(indexTransform, xyzPointsList),
                    /*threaded*/(Pstream::parRun() ? false : true),
                    /*shareOp*/true
                );

                for (auto iter = levelNGrid->cbeginValueOn(); iter; ++iter)
                {
                    const std::vector<openvdb::Coord> voxelVertices =
                        getVoxelVertices(iter.getCoord(), pointsPerEdge_);

                    List<label> cellPoints(8);

                    for (unsigned vertexI=0; vertexI < voxelVertices.size(); vertexI++)
                    {
                        const openvdb::Coord& xyz = voxelVertices[vertexI];

                        cellPoints[vertexI] = cLevelNIndexGridAccessor.getValue(xyz);
                    }

                    //populate boundary and internal face grids
                    List<face> voxelFaces(6);
                    voxelFaces[0] = face{cellPoints[0], cellPoints[4], cellPoints[7], cellPoints[3]}; // x-min
                    voxelFaces[1] = face{cellPoints[1], cellPoints[2], cellPoints[6], cellPoints[5]}; // x-max
                    voxelFaces[2] = face{cellPoints[0], cellPoints[1], cellPoints[5], cellPoints[4]}; // y-min
                    voxelFaces[3] = face{cellPoints[3], cellPoints[7], cellPoints[6], cellPoints[2]}; // y-max
                    voxelFaces[4] = face{cellPoints[0], cellPoints[3], cellPoints[2], cellPoints[1]}; // z-min
                    voxelFaces[5] = face{cellPoints[4], cellPoints[5], cellPoints[6], cellPoints[7]}; // z-max

                    const std::vector<openvdb::Coord> voxelFaceCentres =
                        getVoxelFaceCentres(iter.getCoord(), pointsPerEdge_);

                    forAll(voxelFaces, faceI)
                    {
                        label minIndex = Foam::min(voxelFaces[faceI]);

                        openvdb::Coord xyzFaceCoarseLocal;
                        InterfaceFace interfaceFace = {false, xyzFaceCoarseLocal};

                        // check if face is at cellLevelInterface,
                        // if so a face vertex at this level matches face centre of coarser level
                        if (cellLevel > 0 && minIndex < nPointsCoarse)
                        {
                            interfaceFace =
                                isCellLevelInterfaceFace<FaceGrid>
                                (
                                    voxelFaces[faceI],
                                    faceGrids[cellLevel-1]->transformPtr(),
                                    cCoarseFaceGridAccessor
                                );
                        }

                        const openvdb::Coord& xyz = voxelFaceCentres[faceI];

                        label ownerValue = cLevelNOwnerGridAccessor.getValue(xyz);

                        // if face has owner value already it means it was added previously from an other voxel,
                        // therefore is an internal face (shared by 2 voxel)
                        // and the current voxel is neighbour of this face
                        if (ownerValue > -1)
                        {
                            //if (debug)
                            //{
                            //    Info<< "\nFound internal face: "
                            //        << voxelFaces[faceI]
                            //        << endl;
                            //    std::cout<< "faceIndex " << xyz <<std::endl;
                            //}

                            levelNFaceGridAccessor.setValueOff(xyz);

                            owner.append(ownerValue);
                            neighbour.append(cellCount);
                            voxelFaces[faceI].flip();
                            faces.append(voxelFaces[faceI]);
                        }
                        else if (interfaceFace.isInterfaceFace)
                        {
                            //if (debug)
                            //{
                            //    Info<< "\nFound internal face at cellLevel interface: "
                            //        << voxelFaces[faceI]
                            //        << endl;
                            //    std::cout<< "faceIndex " << xyz
                            //        << " - xyzFaceCoarseLocal " << interfaceFace.xyz  <<std::endl;
                            //}

                            levelNCellLevelInterfaceFaceGridAccessor.setValueOn(xyz, true);

                            levelNFaceGridAccessor.setValueOff(xyz);

                            //find owner of coarse face
                            auto oldOwner = cCoarseOwnerGridAccessor.getValue(interfaceFace.xyz);

                            //mark coarse face as refined
                            coarseRefinedBoundaryFaceGridAccessor.setValueOn(interfaceFace.xyz);

                            owner.append(oldOwner);
                            neighbour.append(cellCount);
                            voxelFaces[faceI].flip();
                            faces.append(voxelFaces[faceI]);
                        }
                        else
                        {
                            //if (debug)
                            //{
                            //    Info<< "\nAdding boundary face: "
                            //        << voxelFaces[faceI]
                            //        << endl;
                            //    std::cout<< "faceIndex " << xyz <<std::endl;
                            //}

                            openvdb::math::Vec4ui f
                            (
                                voxelFaces[faceI][0],
                                voxelFaces[faceI][1],
                                voxelFaces[faceI][2],
                                voxelFaces[faceI][3]
                            );

                            levelNFaceGridAccessor.setValueOn(xyz, f);

                            levelNOwnerGridAccessor.setValueOn(xyz, cellCount);

                            levelNPatchIDGridAccessor.setValueOn(xyz, patchID);
                        }
                    }

                    if (levelNvoxelCount > 0 && levelNvoxelCount % checkInterval == 0)
                    {
                        label memoryLimit = 100000; //TODO
                        label memory = mem.update().size() / 1000;

                        //Pout<<"Memory : " << memory << "Mb"<<endl;

                        if (memory > memoryLimit)
                        {
                            FatalErrorInFunction
                                << "Exceeded memory limit ("
                                << memoryLimit
                                << " Mb), aborting"
                                << abort(FatalError);
                        }
                        //time += runTime.cpuTimeIncrement();
    	        	    //printf("Generated %i/%i level %i cells in %fs (mem %i Mb)\r", levelNvoxelCount, activeVoxels, cellLevel, time, memory);
    	        	    //fflush(stdout);
    	            }

                    cellCount++;
                    levelNvoxelCount++;
                }
            }//Timer create cells

            {//Timer
                //Timer timer("write fields cellLevel " + std::to_string(cellLevel));

                scalarField levelNCells(levelNvoxelCount, cellLevel);
                cellLevelField.append(levelNCells);

                scalarField patchIDCells(levelNvoxelCount, patchID);
                patchIDField.append(patchIDCells);

                if (Pstream::parRun())
                {
                    scalarField procIDCells(levelNvoxelCount, Pstream::myProcNo());
                    procIDField.append(procIDCells);
                }

                for (auto iter = levelNGrid->cbeginValueOn(); iter; ++iter)
                {
                    distanceField.append(iter.getValue());
                }
            }// Timer write fields
        } // for patchID

        // switch off boundary faces of coarser cellLevel
        // and make sure there are exactly 4 faces at finer level to replace coarse face
        if (cellLevel > 0)
        {
            //Timer timer("splitHex faces " + std::to_string(cellLevel));
            for (auto iter = refinedBoundaryFaceGrids[cellLevel-1]->cbeginValueOn(); iter; ++iter)
            {
                const openvdb::Coord& coarseFaceCoord = iter.getCoord();

                const openvdb::math::Vec4ui& f = cCoarseFaceGridAccessor.getValue(coarseFaceCoord);

                coarseFaceGridAccessor.setValueOff(coarseFaceCoord);

                hVDB_.splitVoxelFace
                (
                    xyzPointsList,
                    pointCount,
                    f,                                      // face
                    coarseFaceCoord,                        // coarseFaceCoord
                    indexTransform,                         // fineTransform
                    faceGrids[cellLevel-1]->transformPtr(), // coarseTransform
                    levelNIndexGridAccessor,                // fineIndexGridAccessor
                    levelNAddedIndexesGridAccessor,         // fineAddedIndexesGridAccessor
                    cLevelNCellLevelInterfaceFaceGridAccessor, // cellLevelInterfaceFaceGridAccessor
                    cCoarseOwnerGridAccessor,               // coarseOwnerGridAccessor
                    levelNOwnerGridAccessor,                // fineOwnerGridAccessor
                    cCoarsePatchIDGridAccessor,             // coarsePatchIDGridAccessor
                    levelNPatchIDGridAccessor,              // finePatchIDGridAccessor
                    levelNFaceGridAccessor,                 // fineFaceGridAccessor
                    false                                   // forceAdd
                );
            } // for refinedBoundaryFace valueOn

            xyzPointsList.resize(pointCount);
            openvdb::tools::foreach
            (
                levelNAddedIndexesGrid->cbeginValueOn(),
                UpdatePointListOp(indexTransform, xyzPointsList),
                /*threaded*/(Pstream::parRun() ? false : true),
                /*shareOp*/true
            );
        } // if cellLevel > 0

        //update indexGrid of next cellLevel
        if (cellLevel < maxCellLevel)
        {
            //Timer timer("update indexGrid next cellLevel " + std::to_string(cellLevel));

            //thread-safe accessors
            openvdb::tree::ValueAccessorRW<const IndexGrid::TreeType, true> coarseIndexAcc(indexGrids[cellLevel]->tree());
            openvdb::tree::ValueAccessorRW<IndexGrid::TreeType, true> fineIndexAcc(indexGrids[cellLevel+1]->tree());

            openvdb::tools::foreach
            (
                faceGrids[cellLevel]->cbeginValueOn(),
                UpdateFinerIndexGridOp
                (
                    coarseIndexAcc,
                    fineIndexAcc,
                    *indexGrids[cellLevel],
                    *indexGrids[cellLevel + 1],
                    xyzPointsList
                ),
                /*threaded*/(Pstream::parRun() ? false : true),
                /*shareOp*/true
            );
        } //if (cellLevel < maxCellLevel)

        //Info<< "Generated " << activeVoxels
        //    << " cells of level "<< cellLevel
        //    << " patchID "<< patchID
        //    << " in " << runTime.cpuTimeIncrement()
        //    << "s (clockTime "<< runTime.clockTimeIncrement()
        //    << "s)"
        //    << endl;

        //if (debug)
        if (false)
        {
            std::cout<<"level"<<cellLevel<<"IndexGrid : Memory usage " << indexGrids[cellLevel]->memUsage() <<std::endl;
                        indexGrids[cellLevel]->pruneGrid();
            std::cout<<"level"<<cellLevel<<"IndexGrid : Memory usage " << indexGrids[cellLevel]->memUsage() <<std::endl;

            std::cout<<"level"<<cellLevel<<"FaceGrid : Memory usage " << faceGrids[cellLevel]->memUsage() <<std::endl;
                        faceGrids[cellLevel]->pruneGrid();
            std::cout<<"level"<<cellLevel<<"FaceGrid : Memory usage " << faceGrids[cellLevel]->memUsage() <<std::endl;

            std::cout<<"level"<<cellLevel<<"CellLevelInterfaceFaceGrid : Memory usage " << cellLevelInterfaceFaceGrids[cellLevel]->memUsage() <<std::endl;
                        cellLevelInterfaceFaceGrids[cellLevel]->pruneGrid();
            std::cout<<"level"<<cellLevel<<"CellLevelInterfaceFaceGrid : Memory usage " << cellLevelInterfaceFaceGrids[cellLevel]->memUsage() <<std::endl;

            std::cout<<"level"<<cellLevel<<"refinedBoundaryFaceGrid : Memory usage " << refinedBoundaryFaceGrids[cellLevel]->memUsage() <<std::endl;
                        refinedBoundaryFaceGrids[cellLevel]->pruneGrid();
            std::cout<<"level"<<cellLevel<<"refinedBoundaryFaceGrid : Memory usage " << refinedBoundaryFaceGrids[cellLevel]->memUsage() <<std::endl;

            std::cout<<"level"<<cellLevel<<"OwnerGrid : Memory usage " << ownerGrids[cellLevel]->memUsage() <<std::endl;
                        ownerGrids[cellLevel]->pruneGrid();
            std::cout<<"level"<<cellLevel<<"OwnerGrid : Memory usage " << ownerGrids[cellLevel]->memUsage() <<std::endl;
        }

        if (Pstream::parRun())
        {
            BoolGrid::Ptr maskGrid = BoolGrid::create(false);

            maskGrid->setTransform(faceGrids[cellLevel]->transformPtr());

            maskGrid->topologyUnion(*faceGrids[cellLevel]);

            procBoundaryGrids.push_back(maskGrid);
        }

        Info<< "Generated " << cellLevelActiveVoxels
            << " cells of level "<< cellLevel
            //<< " in " << runTime.cpuTimeIncrement()
            //<< "s (clockTime "<< runTime.clockTimeIncrement()
            //<< "s)"
            << endl;
    } // for cellLevel

    //////////////////////////////////////
    // Assign faces to boundary patches
    //////////////////////////////////////

    // initialize defaultFaces
    label defaultSize = 0;
    DynamicList<face>  defaultFaces(0);
    DynamicList<label> defaultOwner(0);

    //initialize list of patchID faceLists
    List<DynamicList<face>> patchIDsFaces(boundaryPatchNames_.size());
    List<DynamicList<label>> ownersFaces(boundaryPatchNames_.size());

    for (label cellLevel=0; cellLevel <= maxCellLevel; cellLevel++)
    {
        FaceGrid::Ptr levelNFaceGrid = faceGrids[cellLevel];

        IndexGrid::ConstAccessor cLevelNOwnerGridAccessor = ownerGrids[cellLevel]->getConstAccessor();

        for (auto iter = levelNFaceGrid->cbeginValueOn(); iter; ++iter)
        {
            openvdb::math::Vec4ui f = iter.getValue();

            label patchID = patchIDGrids[cellLevel]->getConstAccessor().getValue(iter.getCoord());

            openvdb::Int32 ownerValue = cLevelNOwnerGridAccessor.getValue(iter.getCoord());

            //if (debug)
            //{
            //    std::cout<<"Boundary face "<< iter.getCoord()
            //            << " cellLevel " << cellLevel
            //            << " patchID " << patchID
            //            << std::endl;
            //}

            if (patchID >= 0)
            {
                patchIDsFaces[patchID].append
                (
                    face{label(f[0]), label(f[1]), label(f[2]), label(f[3])}
                );

                ownersFaces[patchID].append
                (
                    ownerValue
                );

                boundaryPatchSizes_[patchID]++;
            }
            else
            {
                defaultFaces.append
                (
                    face{label(f[0]), label(f[1]), label(f[2]), label(f[3])}
                );

                defaultOwner.append(ownerValue);

                defaultSize++;
            }
        }
    } // for cellLevel

    label patchStart = faces.size();

    forAll(patchIDsFaces, patchI)
    {
        faces.append(patchIDsFaces[patchI]);
        owner.append(ownersFaces[patchI]);

        boundaryPatchStarts_[patchI] = patchStart;

        patchStart += boundaryPatchSizes_[patchI];
    }

    if (defaultFaces.size())
    {
        faces.append(defaultFaces);
        owner.append(defaultOwner);

        word defaultFacesName = "defaultFaces";
        word defaultFacesType = emptyPolyPatch::typeName;

        boundaryPatchNames_.append(defaultFacesName);
        boundaryPatchTypes_.append(defaultFacesType);

        boundaryPatchStarts_.append(patchStart);
        boundaryPatchSizes_.append(defaultSize);

        patchStart += defaultSize;
    }

    //transform points from index space to world space
    points_ =
        pointField
        (
            ((xyzPointsList / pointsPerEdge_) * hVDB_.voxelSize())
        );

    faces_     = faces.xfer();
    owner_     = owner.xfer();
    neighbour_ = neighbour.xfer();

    cellLevelField_ = xferMoveTo<scalarField>(cellLevelField);
    distanceField_  = xferMoveTo<scalarField>(distanceField);
    patchIDField_   = xferMoveTo<scalarField>(patchIDField);
    procIDField_    = xferMoveTo<scalarField>(procIDField);
} // convert()

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VDBGridsConverter::VDBGridsConverter
(
    foamVDB& hVDB,
    gridLists& allCellsGrids,
    label pointsPerEdge /* = 2 */
)
:
    hVDB_(hVDB),
    allCellsGrids_(allCellsGrids),
    pointsPerEdge_(pointsPerEdge)
{
    convert();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::VDBGridsConverter::~VDBGridsConverter()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //


// ************************************************************************* //
