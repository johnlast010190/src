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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

//OpenVDB
#include <openvdb/openvdb.h>
#include <openvdb/io/Stream.h>
#include <openvdb/tools/GridTransformer.h> // for resampleToMatch()
#include <openvdb/tools/FastSweeping.h> //for maskSdf
#include "MultiResGrid.h"
#include <tbb/global_control.h>
#include <tbb/tick_count.h>

//OpenFOAM
#include "VDBGridsConverter.H"
#include "cfdTools/general/include/fvCFD.H"
#include "global/argList/argList.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "memInfo/memInfo.H"
#include "fields/Fields/DynamicField/DynamicField.H"
//#include "sampledSets.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "db/IOobjectList/IOobjectList.H"
//#include "include/parallelForAll.H"
#include "refinementSurfaces/refinementSurfaces.H"

using namespace Foam;

// * * *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "rasterizeFields.H"

template<typename GridType>
void sendGrids
(
    const std::vector<typename GridType::Ptr>& inputGrids,
    PstreamBuffers& pBufSize,
    PstreamBuffers& pBufGrid,
    const label node
)
{
    std::ostringstream ostr(std::ios_base::binary);

    openvdb::GridPtrVecPtr grids(new openvdb::GridPtrVec);

    for (size_t i = 0; i < inputGrids.size(); i++)
    {
        grids->push_back(inputGrids[i]);
    }

    openvdb::io::Stream(ostr).write(*grids);

    std::streamsize strSize = ostr.str().size();

    std::string gridName = inputGrids[0]->getName();

    UOPstream toProcSize(node, pBufSize);
    toProcSize << strSize;

    UOPstream toProc(node, pBufGrid);
    toProc.write(ostr.str().data(), strSize);
} // sendGrids


template<typename GridType>
std::vector<typename GridType::Ptr> receiveGrids
(
    PstreamBuffers& pBufSize,
    PstreamBuffers& pBufGrid,
    const label node
)
{
    // non-blocking communication
    UIPstream fromProcSize(node, pBufSize);
    std::streamsize strSize;
    fromProcSize >> strSize;

    char * incoming = new char[strSize];

    UIPstream fromProc(node, pBufGrid);
    fromProc.read(incoming, strSize);

    std::stringstream ss;
    ss.write(incoming, strSize);

    std::istringstream istr;
    istr.str(ss.str());

    openvdb::io::Stream strm(istr);
    openvdb::GridPtrVecPtr streamGrids = strm.getGrids();

    std::vector<typename GridType::Ptr> grids;

    for (openvdb::GridPtrVec::iterator it = streamGrids->begin(); it < streamGrids->end(); ++it)
    {
        typename GridType::Ptr grid = openvdb::gridPtrCast<GridType>(*it);

        grids.push_back(grid);
    }

    return grids;
} // receiveGrids


template<typename GridType>
void gatherGrids
(
    std::vector<typename GridType::Ptr>& grids
)
{
    // master collects grids from other processors
    PstreamBuffers pBufSize(Pstream::commsTypes::nonBlocking);
    PstreamBuffers pBufGrid(Pstream::commsTypes::nonBlocking);

    // send my grids to master
    if (Pstream::myProcNo() != Pstream::masterNo())
    {
        sendGrids<GridType>
        (
            grids,
            pBufSize,
            pBufGrid,
            Pstream::masterNo()
        );
    }

    pBufSize.finishedSends();
    pBufGrid.finishedSends();

    if (Pstream::myProcNo() == Pstream::masterNo())
    {
        for
        (
            int proci = Pstream::firstSlave();
            proci <= Pstream::lastSlave();
            proci++
        )
        {
            std::vector<typename GridType::Ptr> receivedGrids =
                receiveGrids<GridType>
                (
                    pBufSize,
                    pBufGrid,
                    proci
                );

            for (size_t i = 0; i < grids.size(); i++)
            {
                grids[i]->tree().merge(receivedGrids[i]->tree()); // empties receiving tree
            }
        }
    }
} // gatherGrids


template<typename GridType>
void restrictGrids
(
    const std::vector<typename GridType::Ptr>& fineGrids,
    std::vector<typename GridType::Ptr>& coarseGrids,
    const label maxLevel
)
{
    for (size_t i = 0; i < fineGrids.size(); i++)
    {
        const word name = fineGrids[i]->getName();

        typename GridType::Ptr coarseGrid = GridType::create(typename GridType::ValueType(SMALL));
        openvdb::math::Transform::Ptr xform = fineGrids[i]->transformPtr()->copy();
        xform->preScale( scalar(1 << maxLevel) );
        coarseGrid->setTransform(xform);
        coarseGrid->setName(name);

        openvdb::util::NullInterrupter interrupter;
        openvdb::tools::doResampleToMatch<openvdb::tools::PointSampler>(*fineGrids[i], *coarseGrid, interrupter);

        label fineVoxels = fineGrids[i]->activeVoxelCount();
        label coarseVoxels = coarseGrid->activeVoxelCount();
        //reduce(fineVoxels, sumOp<label>());
        //reduce(coarseVoxels, sumOp<label>());
        Info<< "\n[0] Active voxels in grid " << name << nl
            << setw(20) << "high resolution: "
            << setw(10) << fineVoxels << nl
            << setw(20) << "low resolution: "
            << setw(10) << coarseVoxels
            << endl;

        coarseGrids.push_back(coarseGrid);
    }
}

template<typename GridType>
inline void restrictGrids
(
    tbb::task_group& tasks,
    const std::vector<typename GridType::Ptr>& fineGrids,
    std::vector<typename GridType::Ptr>& coarseGrids,
    const label maxLevel
)
{
    tasks.run
    (
        [&fineGrids, &coarseGrids, maxLevel]
        {
            restrictGrids<GridType>(fineGrids, coarseGrids, maxLevel);
        }
    );
}

inline void setValue
(
    scalar& value,
    const FloatGrid::ValueType& v
)
{
    value = v;
}

inline void setValue
(
    vector& value,
    const Vec3SGrid::ValueType& v
)
{
    value.x() = v.x();
    value.y() = v.y();
    value.z() = v.z();
}

template<class Type, typename GridType>
void writeVolFields
(
    const fvMesh& mesh,
    std::vector<typename GridType::Ptr>& grids
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    for (size_t i = 0; i < grids.size(); i++)
    {
        const word name = grids[i]->getName();

        Info<< "Writing field: " << name << endl;

        fieldType volField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensioned<Type>("zero", dimless, Zero),
            zeroGradientFvPatchField<Type>::typeName
        );

        label cellI = 0;
        for (auto iter = grids[i]->cbeginValueOn(); iter; ++iter)
        {
            setValue(volField[cellI++], iter.getValue());
        }

        volField.write();
    }
}


autoPtr<fvMesh> writeMesh
(
    std::vector<FloatGrid::Ptr> cellLevelGrids,
    foamVDB& ovdb,
    Time& runTime,
    point levelSetOffset
)
{
    std::vector<std::vector<FloatGrid::Ptr>> allCellsGrids;

    allCellsGrids.push_back(cellLevelGrids);

    VDBGridsConverter vdbGrids
    (
        ovdb,
        allCellsGrids
    );

    vdbGrids.points() -= levelSetOffset;

    fileName path = runTime.constant();

    autoPtr<fvMesh> meshPtr
    (
        new fvMesh
        (
            IOobject
            (
                polyMesh::defaultRegion,
                path,
                runTime
            ),
            xferCopy(vdbGrids.points()),
            vdbGrids.faces().xfer(),
            vdbGrids.owner().xfer(),
            vdbGrids.neighbour().xfer()
        )
    );

    fvMesh& mesh = meshPtr();

    List<polyPatch*> polyPatches(vdbGrids.boundaryPatchStarts().size());

    forAll(polyPatches, patchI)
    {
        polyPatches[patchI] = polyPatch::New
        (
            vdbGrids.boundaryPatchTypes()[patchI],
            vdbGrids.boundaryPatchNames()[patchI],
            vdbGrids.boundaryPatchSizes()[patchI],
            vdbGrids.boundaryPatchStarts()[patchI],
            patchI,
            mesh.boundaryMesh()
        ).ptr();
    }

    mesh.addFvPatches(polyPatches);

    Info<< "\nWriting polyMesh" <<endl;
    mesh.write();

    return meshPtr;
} //writeMesh


struct Norm
{
    // assuming symmetric narrow-band levelset
    scalar max_, range_;
    Norm(const scalar& max): max_(max), range_(2*max) {}
    inline void operator()(const FloatGrid::ValueOnIter& iter) const
    {
        iter.setValue(-(*iter - max_) / range_);
    }
};


label checkProcessorFolders(const fileName& dir)
{
    label nProcs = 0;

    while
    (
        isDir
        (
            dir/"processor"
          + Foam::name(++nProcs)
        )
    )
    {}

    //wait for all procs to finish counting
    gMax(labelList(Pstream::nProcs(), 0));

    return nProcs;

    //// create dummy source folders
    //if (nProcs < Pstream::nProcs())
    //{
    //    instantList timeDirs;
    //    const word procFolder = "processor" + Foam::name(Pstream::myProcNo());
    //    if (Pstream::master())
    //    {
    //        const bool oldParRun = Pstream::parRun();
    //        Pstream::parRun() = false;
    //        timeDirs = Time::findTimes(dir/procFolder, "constant");
    //        Pstream::parRun() = oldParRun;
    //    }
    //    Pstream::scatter(timeDirs);
    //    forAll(timeDirs, i)
    //    {
    //        mkDir(dir/procFolder/timeDirs[i].name());
    //    }
    //}
    //else if (nProcs > Pstream::nProcs())
    //{
    //    FatalError
    //        << "Running with " << Pstream::nProcs()
    //        << " processors but case "
    //        << dir << " is decomposed in "
    //        << nProcs << " processors. "
    //        << nl << "Please run with maximum number of source and target processors"
    //        << exit(FatalError);
    //}
}

openvdb::CoordBBox toCoordBBox
(
    const boundBox& bb,
    openvdb::math::Transform::ConstPtr transform
)
{
    const openvdb::Coord bbMin
    (
        transform->worldToIndexCellCentered
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
        transform->worldToIndexCellCentered
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


template<class Type, typename GridType>
void initializeGrids
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields,
    const scalar voxelSize,
    wordList& names,
    std::vector<typename GridType::Ptr>& grids
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    IOobjectList fields = objects.lookupClass(fieldType::typeName);

    forAllIter(IOobjectList, fields, fieldIter)
    {
        const word& fieldName = fieldIter()->name();

        if (selectedFields.empty() || selectedFields.found(fieldName))
        {
            openvdb::math::Transform::Ptr linearTransform =
                openvdb::math::Transform::createLinearTransform(voxelSize);

            Info<<"Creating VDB grid for field " << fieldName << endl;

            names.append(fieldName);

            typename GridType::Ptr vdbGrid = GridType::create(typename GridType::ValueType(SMALL));

            vdbGrid->setTransform(linearTransform);

            vdbGrid->setName(fieldName);

            grids.push_back(vdbGrid);
        }
    }
}


// very lightweight struct for reading basic mesh data
// and multi-threaded calculation of face-cell addressing
struct polyMeshBasic
{
    pointIOField      points_;
    faceCompactIOList faces_;
    labelIOList       owner_;
    labelIOList       neighbour_;
    boundBox          bounds_;
    label             nCells_;
    cellList          cellFaceAddr_;

    polyMeshBasic
    (
        const Time& runTime,
        const fileName& pointsInstance,
        const fileName& facesInstance,
        const point levelSetOffset
    )
    :
        points_
        (
            IOobject
            (
                "points",
                pointsInstance,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            )
        ),
        faces_
        (
            IOobject
            (
                "faces",
                facesInstance,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            )
        ),
        owner_
        (
            IOobject
            (
                "owner",
                facesInstance,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            )
        ),
        neighbour_
        (
            IOobject
            (
                "neighbour",
                facesInstance,
                polyMesh::meshSubDir,
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            )
        ),
        bounds_(points_, /*doReduce*/false),
        nCells_(-1),
        cellFaceAddr_()
    {
        points_ += levelSetOffset;
        bounds_.min() += levelSetOffset;
        bounds_.max() += levelSetOffset;
    }

    inline label nCells() const
    {
        return nCells_;
    }

    const boundBox& bounds() const
    {
        return bounds_;
    }

    const pointField& points() const
    {
        return points_;
    }

    const faceList& faces() const
    {
        return faces_;
    }

    const cellList& cells() const
    {
        return cellFaceAddr_;
    }

    label reduce(const labelList& llist, label identity)
    {
        const auto nThreads =
            tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

        const auto grainSize =
            std::max<size_t>
            (
                llist.size() / nThreads,
                1024
            );

        //Info<<"grainSize : " << grainSize <<endl;
        const tbb::blocked_range<size_t> range(0, llist.size(), grainSize);

        return tbb::parallel_reduce
        (
            range,
            identity,
            /*func*/[&](const tbb::blocked_range<size_t>& r, const label init)
            {
                //if (Pstream::master())
                //{
                //    std::cout<< "    reduce thread ID "
                //        << std::this_thread::get_id()
                //        << " - range " << r.begin()
                //        << ":" << r.end()
                //        << std::endl;
                //}
                label res = init;
                for (size_t i = r.begin(); i != r.end(); ++i)
                {
                    res = Foam::max(res, llist[i]);
                }
                return res;
            },
            /*reduction*/[](const label x, const label y)
            {
                return Foam::max(x, y);
            },
            tbb::simple_partitioner()
        );
    }

    label countCells(const labelList& own, const labelList& nei)
    {
        label nCells = reduce(own, -1);

        nCells = reduce(nei, nCells);

        ++nCells;

        return nCells;
    }

    void calcCells()
    {
        nCells_ = countCells(owner_, neighbour_);

        calcCells(cellFaceAddr_, owner_, neighbour_);
    }


    inline void parallelCount
    (
        const labelList& list,
        std::vector<std::atomic<label>>& ncf
    )
    {
        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, list.size()),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t facei = r.begin(); facei < r.end(); facei++)
                {
                    const label& celli = list[facei];
                    ncf[celli]++;
                }
            }
        );
    }

    inline void parallelFill
    (
        const labelList& list,
        std::vector<std::atomic<label>>& ncf,
        cellList& cellFaceAddr
    )
    {
        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, list.size()),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t facei = r.begin(); facei < r.end(); facei++)
                {
                    const label& celli = list[facei];
                    cellFaceAddr[celli][ncf[celli]++] = facei;
                }
            }
        );
    }

    void calcCells
    (
        cellList& cellFaceAddr,
        const labelList& own,
        const labelList& nei
    )
    {
        const auto nThreads =
            tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

        cellFaceAddr.setSize(nCells_);

        std::vector<std::atomic<label>> ncf(nCells_);

        tbb::task_group tasks;
        tasks.run
        (
            [&]{ parallelCount(own, ncf); }
        );
        tasks.run_and_wait
        (
            [&]{ parallelCount(nei, ncf); }
        );

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, nCells_),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t celli = r.begin(); celli < r.end(); celli++)
                {
                    cellFaceAddr[celli].setSize(ncf[celli]);
                    ncf[celli] = 0;
                }
            }
        );

        tasks.run
        (
            [&]{ parallelFill(own, ncf, cellFaceAddr); }
        );
        tasks.run_and_wait
        (
            [&]{ parallelFill(nei, ncf, cellFaceAddr); }
        );
    } // calcCells
}; //polyMeshBasic


template<class Type>
void readFields
(
    const label fieldSize,
    const Time& runTime,
    const wordList& names,
    List<Field<Type>>& fields
)
{
    fields.setSize(names.size());

    // serial loop because read from stream
    // can be done by one thread at a time
    forAll(names, i)
    {
        const word& fieldName = names[i];

        localIOdictionary fieldDict
        (
            IOobject
            (
                fieldName,
                runTime.timeName(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                /*registerObject*/false
            ),
            GeometricField<Type, fvPatchField, volMesh>::typeName
        );

        Field<Type> fieldSource("internalField", fieldDict, fieldSize);

        fields[i] = fieldSource;
    }
} // readFields


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Rasterizes source fields and maps them to a target mesh or a sampleBox"
    );

    argList::noFunctionObjects();

    argList::noCheckProcessorDirectories();

    // Create argList. This will check for non-existing processor dirs.
    // Create processor directory if non-existing
    Foam::argList args
    (
        argc,
        argv,
        /*checkArgs*/ true,
        /*checkOpts*/ true,
        /*initialise*/ true,
        /*needsThread*/true
    );

    Info<<"MPI_THREAD_MULTIPLE "
        << (Pstream::haveThreads() ? "true" : "false")
        << endl;

    //if (Pstream::parRun() && !isDir(args.path()))
    //{
    //    mkDir(args.path());
    //}

    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

    Foam::Time runTimeTarget(Foam::Time::controlDictName, args);

    memInfo mem;
    Info<< "Initial memory " << mem.update().size() << " kB" << endl;

    const word dictName("foamVDBMapDict");

    Info<< "\nReading " << dictName << nl << endl;
    IOdictionary mapDict
    (
        IOobject
        (
            dictName,
            runTimeTarget.system(),
            "",
            runTimeTarget,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    const label threads = mapDict.lookupOrDefault<label>("threads", 0);
    std::unique_ptr<tbb::global_control> control;
    if (threads > 0)
    {
        // note, threads == 0 means use all threads (default), so don't
        // manually create a tbb::global_control in this case
        control.reset
        (
            new tbb::global_control
            (
                tbb::global_control::max_allowed_parallelism,
                threads
            )
        );
    }

    bool debug = mapDict.lookupOrDefault<bool>("debug", false);

    const bool samplingBox = mapDict.found("samplingBox");

    bool visualization = false;

    point samplingBoundMin;
    point samplingBoundMax;
    point levelSetOffset(0,0,0);
    // initialize samplingCoordBBox to LARGE
    openvdb::CoordBBox samplingCoordBBox
    (
        openvdb::Coord::min(),
        openvdb::Coord::max()
    );

    boundBox samplingBBox;

    scalar voxelSize = mapDict.lookupOrDefault<scalar>("voxelSize", -1.0);
    scalar outputVoxelSize = voxelSize;

    label samplingLevel = 0, outputLevel = 0;
    labelList outputDivisions(3);

    FloatGrid::Ptr distanceGrid = FloatGrid::create();

    if (samplingBox)
    {
        const dictionary& boxDict = mapDict.subDict("samplingBox");

        outputLevel = boxDict.lookupOrDefault<label>("outputLevel", 0);
        samplingLevel = boxDict.lookupOrDefault<label>("samplingLevel", 1);

        outputVoxelSize = voxelSize * Foam::pow(2, (samplingLevel - outputLevel));

        samplingBoundMin = boxDict.lookup("min");
        samplingBoundMax = boxDict.lookup("max");

        outputDivisions[0] = label((samplingBoundMax.x() - samplingBoundMin.x()) / outputVoxelSize);
        outputDivisions[1] = label((samplingBoundMax.y() - samplingBoundMin.y()) / outputVoxelSize);
        outputDivisions[2] = label((samplingBoundMax.z() - samplingBoundMin.y()) / outputVoxelSize);

        if (boxDict.found("outputDivisions"))
        {
            outputDivisions = labelList(boxDict.lookup("outputDivisions"));

            if (outputDivisions.size() != 3)
            {
                FatalErrorIn(args.executable())
                    << "outputDivisions not size 3!"
                    << exit(FatalError);
            }

            outputVoxelSize =
            (
                ((samplingBoundMax.x() - samplingBoundMin.x()) / outputDivisions[0])
              + ((samplingBoundMax.y() - samplingBoundMin.y()) / outputDivisions[1])
              + ((samplingBoundMax.z() - samplingBoundMin.z()) / outputDivisions[2])
            ) / 3.0;

            samplingBoundMax.x() = samplingBoundMin.x() + outputVoxelSize*(outputDivisions[0] - 1);
            samplingBoundMax.y() = samplingBoundMin.y() + outputVoxelSize*(outputDivisions[1] - 1);
            samplingBoundMax.z() = samplingBoundMin.z() + outputVoxelSize*(outputDivisions[2] - 1);

            voxelSize = outputVoxelSize / Foam::pow(2, (samplingLevel - outputLevel));
        }

        if (voxelSize <= 0)
        {
            FatalErrorIn(args.executable())
                << "Please specify a positive voxelSize!"
                << exit(FatalError);
        }

        Info<< "Sampling resolution: " << voxelSize << nl
            << "Output   resolution: " << outputVoxelSize << nl
            << "Output bounding box:" << nl
            << "    min " << samplingBoundMin << nl
            << "    max " << samplingBoundMax << nl
            << "Output divisions: " << outputDivisions << nl
            << endl;

        visualization =
            boxDict.lookupOrDefault<bool>("visualization", false);

        levelSetOffset =
            point
            (
                -std::fmod(samplingBoundMin.x(), 1.0),
                -std::fmod(samplingBoundMin.y(), 1.0),
                -std::fmod(samplingBoundMin.z(), 1.0)
            );

        // shift xmin, ymin, zmin of samplingBox to closest integer value
        // will shift back when converting points from index space to world space
        samplingBoundMin += levelSetOffset;
        samplingBoundMax += levelSetOffset;

        samplingBBox.add(samplingBoundMin);
        samplingBBox.add(samplingBoundMax);

        samplingCoordBBox.reset
        (
            openvdb::Coord
            (
                samplingBoundMin.x() / voxelSize,
                samplingBoundMin.y() / voxelSize,
                samplingBoundMin.z() / voxelSize
            ),
            openvdb::Coord
            (
                samplingBoundMax.x() / voxelSize,
                samplingBoundMax.y() / voxelSize,
                samplingBoundMax.z() / voxelSize
            )
        );

        // check for geometry input
        if (boxDict.found("sampleGeometry") && Pstream::master())
        {
            const dictionary& sampleGeometryDict = boxDict.subDict("sampleGeometry");

            const bool maskSdf           = sampleGeometryDict.lookupOrDefault<bool>("maskSdf", false);
            const bool normalizeLevelset = sampleGeometryDict.lookupOrDefault<bool>("normalizeLevelset", false);
            const bool saveHighResGrid   = sampleGeometryDict.lookupOrDefault<bool>("saveHighResGrid", false);
            const bool visualizeLevelSet = sampleGeometryDict.lookupOrDefault<bool>("visualizeLevelSet", false);

            const label nCellsBetweenLevels = sampleGeometryDict.lookupOrDefault<label>("nCellsBetweenLevels", 3);

            const dictionary& geometryDict = sampleGeometryDict.subDict("geometry");

            autoPtr<searchableSurfaces> allGeometryPtr
            (
                new searchableSurfaces
                (
                    IOobject
                    (
                        "abc",                      // dummy name
                        runTimeTarget.time().constant(),  // directory
                        "triSurface",               // instance
                        runTimeTarget.time(),             // registry
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    geometryDict,
                    /*singleRegionName*/true,
                    /*alternativeSurfacePath*/List<fileName>(0),
                    /*maxRegionSize*/labelMax
                )
            );

            autoPtr<refinementSurfaces> surfacesPtr
            (
                new refinementSurfaces
                (
                    allGeometryPtr(),
                    sampleGeometryDict.subDict("refinementSurfaces"),
                    sampleGeometryDict,
                    /*gapLevelIncrement*/0
                )
            );

            Info<< "Read refinement surfaces in = "
                << runTimeTarget.time().cpuTimeIncrement() << " s" << nl << endl;

            refinementSurfaces& surfaces = surfacesPtr();

            labelList surfaceGeometry = surfaces.surfaces();

            label maxLevel = Foam::max(surfaces.maxLevel()) - Foam::min(surfaces.minLevel());

            scalar minEdge = outputVoxelSize / Foam::pow(2, maxLevel);

            Info<< "geometry supersampling resolution: " << minEdge
                << nl << endl;

            foamVDB ovdb(minEdge, maxLevel);

            FloatGrid::Ptr dummyGrid = FloatGrid::create();
            std::vector<FloatGrid::Ptr> cellLevelGrids(ovdb.maxCellLevel() + 1, dummyGrid);

            runTimeTarget.time().clockTimeIncrement();

            #include "createNarrowBandLevelSet.H"

            Info<< "Generated levelset in = "
                << runTimeTarget.time().clockTimeIncrement() << " s" << nl << endl;

        } //if sampleGeometry
    } //if samplingBox

    HashSet<word> selectedFields;
    if (mapDict.found("fields"))
    {
        mapDict.lookup("fields")() >> selectedFields;
    }
    else
    {
        Info<< "Not mapping fields.\nEnd." << endl;
        return 0;
    }

    const fileName rootDirSource = fileName(mapDict.lookup("sourceCase")).toAbsolute();

    if (!isDir(rootDirSource.path()))
    {
        FatalErrorInFunction
            << "source case directory: "
            << rootDirSource
            << " does not exist!"
            << exit(FatalError);
    }

    //check sourceCase processors
    const label nProcsSource = checkProcessorFolders(rootDirSource);

    Info<< "\nsourceCase " << rootDirSource
        << " decomposed in " << nProcsSource
        << " processors." << endl;

    word sourceTimeName("latestTime");

    if (mapDict.found("sourceTime"))
    {
        ITstream sourceTimeEntry = mapDict.lookup("sourceTime");

        if (sourceTimeEntry[0] == word("latestTime"))
        {
            sourceTimeName = word("latestTime");
        }
        else
        {
            scalar sourceTimeScalar = readScalar(mapDict.lookup("sourceTime"));
            sourceTimeName = name(sourceTimeScalar);
        }
    }

    // assign a list of source processor folders to each node
    labelList nSourceMeshesInProc(Pstream::nProcs());
    labelList sourceMeshToProc(identity(nProcsSource));
    labelListList sourceMeshNoInProc =
        fvMeshDistribute::calcMeshToProcMap
        (
            nProcsSource,
            nSourceMeshesInProc,
            sourceMeshToProc
        );
    //e.g. 2 procs, 5 meshes
    // sourceMeshNoInProc
    // 2
    // (
    //   ( 0 2 3 )
    //   ( 1 4 )
    // )

    const labelList& mySourceMeshes =
        sourceMeshNoInProc[Pstream::myProcNo()];

    std::vector<FloatGrid::Ptr> scalarGrids;
    std::vector<Vec3SGrid::Ptr> vectorGrids;
    //TODO tensors, symmTensors etc

    runTimeTarget.time().clockTimeIncrement();

    // initialize scalar and vector grids
    Time runTimeSource0
    (
        Time::controlDictName,
        rootDirSource,
        (
            nProcsSource > 1
          ? fileName(word("processor0"))
          : fileName(word("./"))
        )
    );

    instantList sourceTimes = runTimeSource0.times();
    label sourceTimeIndex = runTimeSource0.timeIndex();

    if (sourceTimeName == "latestTime")
    {
        sourceTimeIndex = sourceTimes.size() - 1;
    }
    else
    {
        IStringStream is(sourceTimeName);
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            readScalar(is)
        );
    }

    runTimeSource0.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

    // Search for list of objects
    IOobjectList objects0(runTimeSource0, runTimeSource0.timeName());

    //Info<< "Objects found in sourceTime "
    //    << runTimeSource0.timeName()
    //    << objects0.toc() << endl;

    Info<< "\nMapping from source time "
        << runTimeSource0.timeName()
        << nl << endl;

    wordList scalarNames, vectorNames;

    initializeGrids<scalar, FloatGrid>
    (
         objects0,
         selectedFields,
         voxelSize,
         scalarNames,
         scalarGrids
    );

    initializeGrids<vector, Vec3SGrid>
    (
         objects0,
         selectedFields,
         voxelSize,
         vectorNames,
         vectorGrids
    );

    const fileName pointsInstance = runTimeSource0.findInstance(polyMesh::meshSubDir, "points");
    const fileName facesInstance  = runTimeSource0.findInstance(polyMesh::meshSubDir, "faces");

    Info<<"\nInitialized grids in "
        << runTimeTarget.time().clockTimeIncrement() << " (cpu "
        << runTimeTarget.time().cpuTimeIncrement() << ") s" << nl << endl;

    // do not update processor boundaries
    bool oldParRun = UPstream::parRun();
    UPstream::parRun() = false;

    scalar rasterTime = 0, allRasterTime = 0;

    forAll(mySourceMeshes, i)
    {
        label proci = sourceMeshNoInProc[Pstream::myProcNo()][i];

        //if (Pstream::master())
        //{
        //    std::cout<< "rasterizeProc thread ID "
        //        << tbb::this_tbb_thread::get_id()
        //        << std::endl;
        //}

        Time runTimeSourceI
        (
            runTimeSource0.controlDict(),
            rootDirSource,
            (
                nProcsSource > 1
              ? fileName(word("processor") + name(proci))
              : fileName(word("./"))
            )
        );

        // assuming all processor folders have same times
        runTimeSourceI.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

        polyMeshBasic meshI(runTimeSourceI, pointsInstance, facesInstance, levelSetOffset);

        Info<< "\nRead source mesh "
            << (
                   nProcsSource > 1
                 ? rootDirSource/fileName(word("processor") + name(proci))
                 : rootDirSource
               )
            << " in "
            << runTimeSourceI.time().clockTimeIncrement() << " (cpu "
            << runTimeSourceI.time().cpuTimeIncrement() << ") s"
            << endl;

        // continue if this processor boundBox is outside of samplingCoordBBox
        const boundBox& procBb = meshI.bounds();
        if (samplingBox && !procBb.overlaps(samplingBBox))
        {
            //Pout<< "Processor" << proci
            //    << " does not overlap with target, skipping..."
            //    << endl;
            continue;
        }

        // calculate face-cell addressing
        runTimeSourceI.time().clockTimeIncrement();
        runTimeSourceI.time().cpuTimeIncrement();

        meshI.calcCells();

        //Info<< "\nCalculated addressing in "
        //    << runTimeSourceI.time().clockTimeIncrement() << " (cpu "
        //    << runTimeSourceI.time().cpuTimeIncrement() << ") s"
        //    << nl << endl;

        // Rasterize volFields
        // ~~~~~~~~~~~~~
        List<scalarField> scalarFields;
        List<vectorField> vectorFields;

        runTimeSourceI.time().clockTimeIncrement();
        runTimeSourceI.time().cpuTimeIncrement();

        readFields<scalar>
        (
            meshI.nCells(),
            runTimeSourceI,
            scalarNames,
            scalarFields
        );

        readFields<vector>
        (
            meshI.nCells(),
            runTimeSourceI,
            vectorNames,
            vectorFields
        );
        Info<< "Read fields of "
            << (word("processor") + name(proci))
            << " in "
            << runTimeSourceI.time().clockTimeIncrement() << " s (cpu "
            << runTimeSourceI.time().cpuTimeIncrement() << " s)"
            << endl;

        rasterizeFields<polyMeshBasic>
        (
            meshI,
            scalarFields,
            vectorFields,
            samplingCoordBBox,
            scalarGrids,
            vectorGrids
        );

        rasterTime = runTimeSourceI.time().clockTimeIncrement();
        Info<< "Rasterized fields in "
            << rasterTime << " s (cpu "
            << runTimeSourceI.time().cpuTimeIncrement() << " s)"
            << nl << endl;

        allRasterTime += rasterTime;
    } //forAll mySourceMeshes

    UPstream::parRun() = oldParRun;

    Info<<"\nRead source mesh and rasterized fields in "
        << runTimeTarget.time().clockTimeIncrement() << " s (cpu "
        << runTimeTarget.time().cpuTimeIncrement() << " s)"
        << nl << "(total rasterization clockTime " << allRasterTime
        << ")" << endl;

    std::vector<FloatGrid::Ptr> coarseScalarGrids;
    std::vector<Vec3SGrid::Ptr> coarseVectorGrids;

    const label maxLevel = samplingLevel - outputLevel;

    if (maxLevel > 0)
    {
        tbb::tick_count t0 = tbb::tick_count::now();

        tbb::task_group tasks;

        restrictGrids<FloatGrid>(tasks, scalarGrids, coarseScalarGrids, maxLevel);
        restrictGrids<Vec3SGrid>(tasks, vectorGrids, coarseVectorGrids, maxLevel);

        tasks.wait();

        tbb::tick_count t1 = tbb::tick_count::now();
        Info<<"restrictGrids in: " << (t1-t0).seconds() << endl;
    }
    else
    {
        for (size_t i = 0; i < scalarGrids.size(); i++)
        {
            coarseScalarGrids.push_back(scalarGrids[i]);
        }
        for (size_t i = 0; i < vectorGrids.size(); i++)
        {
            coarseVectorGrids.push_back(vectorGrids[i]);
        }
    }

    Info<< "\nMerged and restricted grids in "
        << runTimeTarget.time().clockTimeIncrement() << " s (cpu "
        << runTimeTarget.time().cpuTimeIncrement() << " s)" << nl << endl;

    // master collects grids from other processors
    // TODO consider multi-threading
    gatherGrids<FloatGrid>(coarseScalarGrids);
    gatherGrids<Vec3SGrid>(coarseVectorGrids);

    Info<< "\nMaster collected grids in "
        << runTimeTarget.time().clockTimeIncrement() << " s (cpu "
        << runTimeTarget.time().cpuTimeIncrement() << " s)" << nl << endl;

    // write VDB grids to file
    if (Pstream::master())
    {
        openvdb::CoordBBox coordBB =
            toCoordBBox
            (
                samplingBBox,
                coarseScalarGrids[0]->transformPtr()
            );

        coordBB.max().x() = coordBB.min().x() + outputDivisions[0]-1;
        coordBB.max().y() = coordBB.min().y() + outputDivisions[1]-1;
        coordBB.max().z() = coordBB.min().z() + outputDivisions[2]-1;

        openvdb::GridPtrVec fieldGrids;

        const word vdbOutput = mapDict.lookupOrDefault<word>("fieldsOutput", "fields.vdb");

        Info<< "\nSaving";

        for (size_t i = 0; i < coarseScalarGrids.size(); i++)
        {
            Info<< " " << coarseScalarGrids[i]->getName();

            // clip to boundingBox
            coarseScalarGrids[i]->clip(coordBB);

            // explicitly set to SMALL inactive values to avoid conversion errors to python
            for (auto iter = coarseScalarGrids[i]->beginValueOff(); iter; ++iter)
            {
                iter.setValue(FloatGrid::ValueType(SMALL));
            }

            fieldGrids.push_back(coarseScalarGrids[i]);
        }
        for (size_t i = 0; i < coarseVectorGrids.size(); i++)
        {
            Info<< " " << coarseVectorGrids[i]->getName();

            // clip to boundingBox
            coarseVectorGrids[i]->clip(coordBB);

            // explicitly set to SMALL inactive values to avoid conversion errors to python
            for (auto iter = coarseVectorGrids[i]->beginValueOff(); iter; ++iter)
            {
                iter.setValue(Vec3SGrid::ValueType(SMALL));
            }

            fieldGrids.push_back(coarseVectorGrids[i]);
        }

        Info<< " to " << vdbOutput << endl;

        if (coarseScalarGrids[0])
        {
            Info<< "Active voxels: " << coarseScalarGrids[0]->activeVoxelCount() << endl;
        }
        else if (coarseVectorGrids[0])
        {
            Info<< "Active voxels: " << coarseVectorGrids[0]->activeVoxelCount() << endl;
        }

        openvdb::io::File file(vdbOutput);
        file.write(fieldGrids);
        file.close();
    }

    Info<< "\nWrote to disk in "
        << runTimeTarget.time().clockTimeIncrement() << " s (cpu "
        << runTimeTarget.time().cpuTimeIncrement() << " s)" << nl << endl;


    if (samplingBox && visualization && Pstream::master())
    {
        // do not update processor boundaries
        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        Info<< "\nWriting data in OpenFOAM format for visualization" << endl;

        foamVDB ovdb(outputVoxelSize, /*maxLevel*/0);

        std::vector<FloatGrid::Ptr> cellLevelGrids(1, coarseScalarGrids[0]);

        //push distanceField
        //coarseScalarGrids.push_back(distanceGrid);

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

        fvMesh& mesh = meshPtr();

        writeVolFields<scalar, FloatGrid>(mesh, coarseScalarGrids);
        writeVolFields<vector, Vec3SGrid>(mesh, coarseVectorGrids);

        UPstream::parRun() = oldParRun;
    }

    Info<< "\nelapsedCpuTime: " << runTimeTarget.elapsedCpuTime() << "s\n"
        << "elapsedClockTime: " << runTimeTarget.elapsedClockTime() << "s,\n"
        << mem.update().peak() << " kB (peak)"
        << nl << endl;

    return 0;
}

// ************************************************************************* //
