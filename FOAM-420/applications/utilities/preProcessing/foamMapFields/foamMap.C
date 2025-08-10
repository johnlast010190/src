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

Description
    Map two flow fields with arbitrary unstructured mesh.
    Must be run on maximum number of source and target processors.

\*---------------------------------------------------------------------------*/
#include "cfdTools/general/include/fvCFD.H"
#include "foamMap.H"
#include "meshes/polyMesh/polyPatches/derived/baseWall/wall.H"
#include "regionProperties/regionProperties.H"
#include "memInfo/memInfo.H"
#include "redistributePar/loadOrCreateMesh.H" //from redistributePar
#include "fvMeshSubset/fvMeshSubset.H"


void checkProcessorFolders(const fileName& dir)
{
    const fileName lck(dir/"foamMap.lock");

    while (Foam::isFile(lck))
    {
        Info<< "Found lock file " << lck
            << ". Waiting..."<< endl;
        Foam::sleep(10);
    }

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

    // create dummy source folders
    if (nProcs < Pstream::nProcs())
    {
        instantList timeDirs;
        const word procFolder = "processor" + Foam::name(Pstream::myProcNo());
        if (Pstream::master())
        {
            const bool oldParRun = Pstream::parRun();

            //Only create lock file if it doesn't already exist
            if (!Foam::isFile(lck))
            {
                Info<< "Creating lock file" << endl;

                OFstream os(lck);
                os << "foamMap: create dummy source folders\n";
                os.flush();
            }
            else
            {
                //otherwise check again
                Foam::sleep(1);
                checkProcessorFolders(dir);
            }
            Pstream::parRun() = false;
            timeDirs = Time::findTimes(dir/procFolder, "constant");
            Pstream::parRun() = oldParRun;
        }
        Pstream::scatter(timeDirs);
        forAll(timeDirs, i)
        {
            mkDir(dir/procFolder/timeDirs[i].name());
        }
    }
    else if (nProcs > Pstream::nProcs())
    {
        FatalError
            << "Running with " << Pstream::nProcs()
            << " processors but case "
            << dir << " is decomposed in "
            << nProcs << " processors. "
            << nl << "Please run with maximum number of source and target processors"
            << exit(FatalError);
    }
} //checkProcessorFolders


autoPtr<fvMesh> detectAndGetMesh
(
    const word& regionName,
    const Time& runTime,
    boolList& haveMesh
)
{
    fileName meshSubDir;

    if (regionName == polyMesh::defaultRegion)
    {
        meshSubDir = polyMesh::meshSubDir;
    }
    else
    {
        meshSubDir = regionName/polyMesh::meshSubDir;
    }

    fileName masterInstDir;

    if (Pstream::master())
    {
        masterInstDir = runTime.findInstance
        (
            meshSubDir,
            "faces",
            IOobject::READ_IF_PRESENT
        );
    }
    Pstream::scatter(masterInstDir);

    // Check who has a mesh
    const fileName meshPath =
        runTime.path()/masterInstDir/meshSubDir/"faces";

    Info<< "    Checking for mesh in " << meshPath << nl << endl;

    haveMesh = boolList(Pstream::nProcs(), false);
    haveMesh[Pstream::myProcNo()] = isFile(meshPath);
    Pstream::allGatherList(haveMesh);
    Info<< "    Per processor mesh availability : " << haveMesh << endl;

    autoPtr<fvMesh> meshPtr =
        loadOrCreateMesh
        (
            IOobject
            (
                regionName,
                masterInstDir,
                runTime,
                Foam::IOobject::MUST_READ
            ),
            /*removeDummyMesh*/false
        );

    return meshPtr;
} //detectAndGetMesh


void writeDummyFields
(
    const boolList& haveMesh,
    const fvMesh& mesh,
    const Time& runTime
)
{
    autoPtr<fvMeshSubset> subsetterPtr;

    // Find last non-processor patch.
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label nonProcI = -1;

    forAll(patches, patchI)
    {
        if (isA<processorPolyPatch>(patches[patchI]))
        {
            break;
        }
        nonProcI++;
    }

    if (nonProcI == -1)
    {
        FatalErrorInFunction
            << "Cannot find non-processor patch on processor "
            << Pstream::myProcNo() << endl
            << " Current patches:" << patches.names()
            << abort(FatalError);
    }

    // Subset 0 cells, no parallel comms. This is used to create
    // zero-sized fields.
    subsetterPtr.reset(new fvMeshSubset(mesh));
    subsetterPtr().setLargeCellSubset(labelHashSet(0), nonProcI, /*syncPar*/false);

    PtrList<volScalarField> volScalarFields;
    PtrList<volVectorField> volVectorFields;

    IOobjectList objects(mesh, runTime.timeName());

    readFields
    (
        haveMesh,
        mesh,
        subsetterPtr,
        objects,
        volScalarFields
    );

    readFields
    (
        haveMesh,
        mesh,
        subsetterPtr,
        objects,
        volVectorFields
    );

    //write empty field
    if (!haveMesh[Pstream::myProcNo()])
    {
        forAll(volScalarFields, i)
        {
            volScalarFields[i].write();
        }
        forAll(volVectorFields, i)
        {
            volVectorFields[i].write();
        }
    }

    subsetterPtr.reset();
} //writeDummyFields

void reportMemory
(
    memInfo& mem,
    const word& desc,
    scalar& maxUsage
)
{
    scalar memory = scalar(mem.update().size())*1.0e-3;
    if (memory > maxUsage)
    {
        scalar maxUsage = memory;

        Pout<< "Maximum memory usage now: " << maxUsage << " MB," << desc << endl;
    }
}

int sourceMapSeq
(
    int argc,
    char *argv[],
    foamMap& fieldmap,
    const Time& runTime,
    const fvMesh& mesh
)
{
    scalar maxUsage = 0;
    word desc = "before reading source mesh";
    memInfo mem;

    if (fieldmap.reportMemUsage())
    {
        reportMemory(mem, desc, maxUsage);
    }

    if (fieldmap.reportMemUsage())
    {
        word desc = "after reading source mesh";
        reportMemory(mem, desc, maxUsage);
    }

    fieldmap.getFieldTypes(mesh, runTime.timeName());

    fieldmap.createSource
    (
        &mesh,
        &runTime
    );

    fieldmap.source().createFields(runTime.timeName());

    if (fieldmap.reportMemUsage())
    {
        word desc = "after creating source fields";
        reportMemory(mem, desc, maxUsage);
    }

    fieldmap.source().setInputs();

    //if (fieldmap.wdistMap_)
    //{
    fieldmap.getWallDistPatchs
    (
        mesh
    );

    const volScalarField y =
        fieldmap.wallDistPatchs().size() == 0
      ? wallDist::New(mesh).y()
      : wallDist::New(mesh, fieldmap.wallDistPatchs(), "wall").y();

    if (fieldmap.reportMemUsage())
    {
        word desc = "after getting wall distance in source";
        reportMemory(mem, desc, maxUsage);
    }
    //}

    fieldmap.source().buildKdTrees(y);

    if (fieldmap.reportMemUsage())
    {
        word desc = "after building kdTree";
        reportMemory(mem, desc, maxUsage);
    }

    fieldmap.source().storeFields();
    if (fieldmap.reportMemUsage())
    {
        word desc = "after store fields";
        reportMemory(mem, desc, maxUsage);
    }

    return 0;
}

int sourceMap
(
    int argc,
    char *argv[],
    foamMap& fieldmap,
    const Time& runTime,
    const fvMesh& mesh,
    boolList& haveMesh
)
{
    Info<< "Reading source" << endl;

    scalar maxUsage = 0;
    word desc = "before reading source mesh";
    memInfo mem;

    if (fieldmap.reportMemUsage())
    {
        reportMemory(mem, desc, maxUsage);
    }

    const bool allHaveMesh = (findIndex(haveMesh, false) == -1);
    if (!allHaveMesh)
    {
        writeDummyFields
        (
            haveMesh,
            mesh,
            runTime
        );
    }

    fieldmap.getFieldTypes(mesh, runTime.timeName());

    fieldmap.createSource
    (
        &mesh,
        &runTime
    );

    fieldmap.source().createFields(runTime.timeName());

    fieldmap.source().setInputs();

    if (fieldmap.reportMemUsage())
    {
        desc = "after reading source fields";
        reportMemory(mem, desc, maxUsage);
    }

    if (fieldmap.wdistMap_)
    {
        fieldmap.getWallDistPatchs(mesh);

        fieldmap.wait();

        const volScalarField y =
            fieldmap.wallDistPatchs().size() == 0
          ? wallDist::New(mesh).y()
          : wallDist::New(mesh, fieldmap.wallDistPatchs(), "wall").y();

        scalar maxwdist = gMax(y);

        fieldmap.source().setWallDist(y, maxwdist);
    }

    if (fieldmap.reportMemUsage())
    {
        desc = "after calculating wall distance in source case";
        reportMemory(mem, desc, maxUsage);
    }

    fieldmap.source().storeFields();

    if (fieldmap.reportMemUsage())
    {
        desc = "before exiting source case";
        reportMemory(mem, desc, maxUsage);
    }

    if (!haveMesh[Pstream::myProcNo()])
    {
        // We created a dummy mesh file above. Delete it.
        const fileName meshFiles = runTime.path();
        //Pout<< "Removing dummy mesh " << meshFiles << endl;
        rmDir(meshFiles);
    }

    if (!allHaveMesh)
    {
        //wait for all proc to finish
        gMax(labelList(Pstream::nProcs(), 0));
        if (Pstream::master())
        {
            //Info<<"\nRemoving "<< fileName(runTime.rootPath()/"foamMap.lock")<<endl;
            rm(runTime.rootPath()/"foamMap.lock");
        }
    }

    return 0;
}

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "sourceCase",
        "dir",
        "Specify source case directory"
    );

    argList::addOption
    (
        "sourceTime",
        "scalar|latestTime",
        "Specify the time step/iteration of the source case to map from (default: latestTime)"
    );

    argList::addOption
    (
        "targetTime",
        "scalar",
        "Specify the time step/iteration of the target case to map to (default: 0)"
    );

    argList::addOption
    (
        "sourceRegion",
        "word",
        "Specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "Specify the target region"
    );
    argList::addBoolOption
    (
        "allRegions",
        "Operate on all regions in target regionProperties"
    );
    argList::addOption
    (
        "regionMaps",
        "List<Pair<word>>",
        "Specify source region to target region maps"
    );

    argList::addOption
    (
        "mapScalarFields",
        "list",
        "Specify the scalar field names to be mapped. E.g. '(p nut)'"
    );

    argList::addOption
    (
        "mapVectorFields",
        "list",
        "Specify the vector field names to be mapped. E.g. '(U)'"
    );

    argList::addOption
    (
        "nwdist",
        "label",
        "Specify the number of divisions in the wall distance range[0,1]"
    );

    argList::addOption
    (
        "rhoRefSource",
        "scalar",
        "Specify the reference density value of the source field"
    );

    argList::addOption
    (
        "rhoRefTarget",
        "scalar",
        "Specify the reference density value of the target field"
    );

    argList::addOption
    (
        "UrefSource",
        "vector",
        "Specify the reference velocity of the source field. E.g. '(1 0 0)'"
    );

    argList::addOption
    (
        "UrefTarget",
        "vector",
        "Specify the reference velocity of the target field. E.g. '(1 0 0)'"
    );

    argList::addOption
    (
        "UrotDegreeFromSource",
        "scalar",
        "Specify the in-flow velocity vector rotate degree from source to target"
    );

    argList::addOption
    (
        "fieldTypes",
        "HashTable<word>",
        "Specify the property types of each maping field. E.g. '(p pressure U velocity k turbEnergy)'"
    );

    argList::addOption
    (
        "alphaMax",
        "scalar",
        "Specify the upper wall-distance value to be included in the wall-distance-based search trees"
    );

    argList::addBoolOption
    (
        "mapBoundary",
        "Specify whether the map include boundaries"
    );

    argList::addBoolOption
    (
        "interpolation",
        "Specify whether interpolation will be used in the mapping"
    );

    argList::addOption
    (
        "function",
        "map|validate",
        "Specify the purpose of running foamMap: map is for mapping two flows or map a flow onto \
a grid, validate is to validate the method by mapping a field into itself."
    );

    argList::addNote
    (
"foamMap maps an arbitrary number of flow fields from an unstructed mesh domain\n\
onto another unstructured domain. The two domains must be similar in nature but\n\
do not have to be of the same size. The boundary conditions do not need to be \n\
the same. The in-flow velocity can be rotated and scaled from the source domain\n\
to the target domain."
    );

    argList::noCheckProcessorDirectories();

#if !defined( WIN32 ) && !defined( WIN64 )
#include "include/addProfilingOption.H"
#endif

    #include "include/setRootCase.H"

    //check target processors
    checkProcessorFolders(args.rootPath()/args.globalCaseName());

    Foam::Time runTimeTarget(Foam::Time::controlDictName, args);

    foamMap fieldmap;

    IOdictionary dict
    (
        IOobject
        (
            "foamMapDict",
            runTimeTarget.system(),
            "",
            runTimeTarget,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fieldmap.setInput(dict, args);

    fieldmap.myid_     = Pstream::myProcNo();
    fieldmap.masterid_ = Pstream::masterNo();
    fieldmap.nprocs_   = Pstream::nProcs();

    if (Pstream::parRun())
    {
        fieldmap.parInit(fieldmap.nprocs_);
    }

    fieldmap.setOptions(args);

    fieldmap.mapCase_ = cwd();

    const fileName rootDirSource = fileName(fieldmap.sourceCase_).toAbsolute();

    if (!isDir(rootDirSource.path()))
    {
        FatalErrorInFunction
            << "source case directory: "
            << rootDirSource
            << " does not exist!"
            << exit(FatalError);
    }

    fileName caseName = "./";
    if (Pstream::parRun())
    {
        caseName = caseName/fileName(word("processor")+name(fieldmap.myid_));
    }

    Info<< "sourceCase: " << rootDirSource << endl;

    //check sourceCase processors
    checkProcessorFolders(rootDirSource);

    Time runTimeSource
    (
        Time::controlDictName,
        rootDirSource,
        caseName
    );

    instantList sourceTimes = runTimeSource.times();
    label sourceTimeIndex = runTimeSource.timeIndex();

    if (fieldmap.mapTimeName_ == "latestTime")
    {
        sourceTimeIndex = sourceTimes.size() - 1;
    }
    else
    {
        IStringStream is(fieldmap.mapTimeName_);
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            readScalar(is)
        );
    }

    runTimeSource.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);
    fieldmap.mapTimeName_ = runTimeSource.timeName();

    Info<< "sourceTime: " << runTimeSource.timeName() << endl;

    const word oldTargetTime = runTimeTarget.timeName();

    instantList targetTimes = runTimeTarget.times();
    label targetTimeIndex = runTimeTarget.timeIndex();

    if (fieldmap.tgtTimeName_ == "latestTime")
    {
        targetTimeIndex = targetTimes.size() - 1;

        runTimeTarget.setTime
        (
            targetTimes[targetTimeIndex],
            targetTimes[targetTimeIndex].value()
        );
    }
    else
    {
        IStringStream is(fieldmap.tgtTimeName_);

        const instant targetTime(readScalar(is), fieldmap.tgtTimeName_);

        runTimeTarget.setTime(targetTime, 0);
    }

    fieldmap.tgtTimeName_ = runTimeTarget.timeName();

    if
    (
        oldTargetTime != fieldmap.tgtTimeName_
     && !isDir(fileName(args.path()/fieldmap.tgtTimeName_))
    )
    {
        if
        (
            !cp
            (
                fileName(args.path()/oldTargetTime),
                fileName(args.path()/fieldmap.tgtTimeName_),
                /*followLink*/false
            )
        )
        {
            FatalErrorInFunction
                << "Cannot copy time directory "
                << oldTargetTime
                << " to "
                << fieldmap.tgtTimeName_
                << exit(FatalError);
        }
    }

    Info<< "targetTime: " << runTimeTarget.timeName() << endl;

    //warn the user if targetTime is different from startTime in controlDict
    if (Foam::name(runTimeTarget.startTime().value()) != runTimeTarget.timeName())
    {
        Warning
            << "Mapping to time " << runTimeTarget.timeName()
            << " which is different from startTime " << runTimeTarget.startTime().value()
            << endl;
    }

    if (fieldmap.allRegions_)
    {
        Info<< "\nMapping for all regions in target regionProperties"
            << endl;

        fieldmap.sourceRegions_.resize(0);
        fieldmap.targetRegions_.resize(0);

        regionProperties rp(runTimeTarget);

        forAllConstIter(HashTable<wordList>, rp, iter)
        {
            const wordList& regions = iter();
            forAll(regions, i)
            {
                if (findIndex(fieldmap.targetRegions_, regions[i]) == -1)
                {
                    fieldmap.targetRegions_.append(regions[i]);
                    fieldmap.sourceRegions_.append(regions[i]);
                }
            }
        }
        fieldmap.targetRegionDirs_ = fieldmap.targetRegions_;
    }

    scalar maxUsage = 0;
    word desc = "before reading target mesh";
    memInfo mem;

    forAll(fieldmap.sourceRegions_, regioni)
    {
        const word& sourceRegion = fieldmap.sourceRegions_[regioni];
        const word& targetRegion = fieldmap.targetRegions_[regioni];

        Info<< "\nMapping source region: "<< sourceRegion
            << "\n     to target region: "<< targetRegion
            << nl<<endl;

        autoPtr<fvMesh> sourceMeshPtr;

        if (Pstream::parRun())
        {
            boolList haveMesh;

            sourceMeshPtr =
                detectAndGetMesh
                (
                    sourceRegion,
                    runTimeSource,
                    haveMesh
                );

            sourceMap
            (
                argc,
                argv,
                fieldmap,
                runTimeSource,
                sourceMeshPtr(),
                haveMesh
            );

            fieldmap.wait();
        }
        else
        {
            chDir(fieldmap.sourceCase_);

            sourceMeshPtr.reset
            (
                new fvMesh
                (
                    IOobject
                    (
                        sourceRegion,
                        runTimeSource.timeName(),
                        runTimeSource,
                        IOobject::MUST_READ
                    )
                )
            );

            sourceMapSeq
            (
                argc,
                argv,
                fieldmap,
                runTimeSource,
                sourceMeshPtr()
            );
        }

        Info<< "\nReading target" << endl;

        boolList targetHaveMesh;

        autoPtr<fvMesh> targetMeshPtr =
            detectAndGetMesh
            (
                targetRegion,
                runTimeTarget,
                targetHaveMesh
            );

        fvMesh& targetMesh = targetMeshPtr();

        if (fieldmap.reportMemUsage())
        {
            desc = "after reading target mesh";
            reportMemory(mem, desc, maxUsage);
        }

        const bool allHaveMesh = (findIndex(targetHaveMesh, false) == -1);
        if (!allHaveMesh)
        {
            writeDummyFields
            (
                targetHaveMesh,
                targetMesh,
                runTimeTarget
            );
        }

        fieldmap.createTarget
        (
            &targetMesh,
            &runTimeTarget
        );

        fieldmap.target().createFields(fieldmap.tgtTimeName_);
        Info<< "Target field created." << nl << endl;

        fieldmap.checkFieldDimensions();

        sourceMeshPtr.reset();


        if (fieldmap.reportMemUsage())
        {
            desc = "after reading target fields";
            reportMemory(mem, desc, maxUsage);
        }

        fieldmap.target().setInputs();

        if (fieldmap.scaleSource())
        {
            fieldmap.source().transformPoints
            (
                fieldmap.target().gbox(),
                fieldmap.source().xyz()
            );

            Info<< "point transformed..." << endl;
        }

        fieldmap.getDomainFields(); //build searchTree for each processor

        if (fieldmap.reportMemUsage())
        {
            desc = "after get domain fields";
            reportMemory(mem, desc, maxUsage);
        }

        Info<< "domain field obtained..." << endl;

        if (fieldmap.wdistMap_)
        {
            fieldmap.getWallDistPatchs
            (
                targetMesh
            );

            if (fieldmap.wallDistPatchs().size() == 0)
            {
                Info<< "set wall distance..." << endl;
                const volScalarField y = wallDist::New(targetMesh).y();
                scalar maxwdist = gMax(y);
                fieldmap.target().setWallDist(y, maxwdist);
                if (Pstream::parRun())
                {
                    fieldmap.wait();
                }
            }
            else
            {
                const volScalarField y = wallDist::New(targetMesh,fieldmap.wallDistPatchs(),"wall").y();
                scalar maxwdist = gMax(y);
                fieldmap.target().setWallDist(y, maxwdist);
                if (Pstream::parRun())
                {
                    fieldmap.wait();
                }
            }
        } // if wdistMap

        if (fieldmap.reportMemUsage())
        {
            desc = "after creating wall distance field in target";
            reportMemory(mem, desc, maxUsage);
        }

        if (Pstream::parRun())
        {
            Info<< "building parallel search trees...start." << endl;
            fieldmap.buildpKdTrees();

            if (fieldmap.reportMemUsage())
            {
                desc = "after building parallel search trees in target";
                reportMemory(mem, desc, maxUsage);
            }

            fieldmap.mapFields();    //parallel run map fields
            if (fieldmap.reportMemUsage())
            {
                desc = "after mapping fields in target";
                reportMemory(mem, desc, maxUsage);
            }
        }
        else
        {
            if (fieldmap.function_ == "validate")
            {
                Info<< "Validate foamMap..." << endl;
                fieldmap.search_error();
                fieldmap.map_error("p");
                fieldmap.map_error("nuTilda");
                fieldmap.Umap_error("U");
            }
            else
            {
                Info<< "mapping fields..." << endl;
                fieldmap.mapFields();
                if (fieldmap.reportMemUsage())
                {
                    desc = "after mapping fields in target";
                    reportMemory(mem, desc, maxUsage);
                }
            }
        }

        runTimeTarget.writeNow();

        fieldmap.clearSourceTarget();

        if (!targetHaveMesh[Pstream::myProcNo()])
        {
            // We created a dummy mesh file above. Delete it.
            const fileName meshFiles = runTimeTarget.path();
            //Pout<< "Removing dummy mesh " << meshFiles << endl;
            rmDir(meshFiles);
        }
        targetMeshPtr.reset();
    } //forAll regions

    Info<< "ExecutionTime = " << runTimeTarget.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTimeTarget.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
