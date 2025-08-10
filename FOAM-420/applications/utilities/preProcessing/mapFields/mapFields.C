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
    (c) 2011-2016 OpenFOAM Foundation

Application
    mapFields

Group
    grpPreProcessingUtilities

Description
    Maps volume fields from one mesh to another, reading and
    interpolating all fields present in the time directory of both cases.

    Parallel and non-parallel cases are handled without the need to reconstruct
    them first.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "meshToMesh0/meshToMesh0.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "MapMeshes.H"
#include "decompositionModel.H"
#include "regionProperties/regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int readNumProcs
(
    const argList& args,
    const word& optionName,
    const word& regionDir,
    const Time& runTime
)
{
    const word dictName = "decomposeParDict";
    fileName dictFile;
    if (args.optionReadIfPresent(optionName, dictFile) && isDir(dictFile))
    {
        dictFile = dictFile / dictName;
    }

    return readInt
    (
        IOdictionary
        (
            decompositionModel::selectIO
            (
                IOobject
                (
                    dictName,
                    runTime.system(),
                    regionDir,
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                dictFile
            )
        ).lookup("numberOfSubdomains")
    );
}


void mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const meshToMesh0::order& mapOrder,
    const bool subtract
)
{
    if (subtract)
    {
        MapConsistentMesh<minusEqOp>
        (
            meshSource,
            meshTarget,
            mapOrder
        );
    }
    else
    {
        MapConsistentMesh<eqOp>
        (
            meshSource,
            meshTarget,
            mapOrder
        );
    }
}


void mapSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches,
    const meshToMesh0::order& mapOrder,
    const bool subtract
)
{
    if (subtract)
    {
        MapSubMesh<minusEqOp>
        (
            meshSource,
            meshTarget,
            patchMap,
            cuttingPatches,
            mapOrder
        );
    }
    else
    {
        MapSubMesh<eqOp>
        (
            meshSource,
            meshTarget,
            patchMap,
            cuttingPatches,
            mapOrder
        );
    }
}


void mapConsistentSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const meshToMesh0::order& mapOrder,
    const bool subtract
)
{
    if (subtract)
    {
        MapConsistentSubMesh<minusEqOp>
        (
            meshSource,
            meshTarget,
            mapOrder
        );
    }
    else
    {
        MapConsistentSubMesh<eqOp>
        (
            meshSource,
            meshTarget,
            mapOrder
        );
    }
}


wordList addProcessorPatches
(
    const fvMesh& meshTarget,
    const wordList& cuttingPatches
)
{
    // Add the processor patches to the cutting list
    HashTable<label> cuttingPatchTable;
    forAll(cuttingPatches, i)
    {
        cuttingPatchTable.insert(cuttingPatches[i], i);
    }

    forAll(meshTarget.boundary(), patchi)
    {
        if (isA<processorFvPatch>(meshTarget.boundary()[patchi]))
        {
            if
            (
               !cuttingPatchTable.found
                (
                    meshTarget.boundaryMesh()[patchi].name()
                )
            )
            {
                cuttingPatchTable.insert
                (
                    meshTarget.boundaryMesh()[patchi].name(),
                    -1
                );
            }
        }
    }

    return cuttingPatchTable.toc();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "map volume fields from one mesh to another"
    );
    argList::noParallel();
    argList::validArgs.append("sourceCase");

    argList::addOption
    (
        "sourceTime",
        "scalar|'latestTime'",
        "specify the source time"
    );
    argList::addOption
    (
        "sourceRegion",
        "word",
        "specify the source region"
    );
    argList::addOption
    (
        "targetRegion",
        "word",
        "specify the target region"
    );
    argList::addBoolOption
    (
        "allRegions",
        "operate on all regions in target regionProperties"
    );
    argList::addOption
    (
        "regionMaps",
        "List<Pair<word>>",
        "Specify source region to target region maps"
    );
    argList::addBoolOption
    (
        "parallelSource",
        "the source is decomposed"
    );
    argList::addBoolOption
    (
        "parallelTarget",
        "the target is decomposed"
    );
    argList::addBoolOption
    (
        "consistent",
        "source and target geometry and boundary conditions identical"
    );
    argList::addOption
    (
        "mapMethod",
        "word",
        "specify the mapping method"
    );
    argList::addBoolOption
    (
        "subtract",
        "subtract mapped source from target"
    );
    argList::addOption
    (
        "sourceDecomposeParDict",
        "file",
        "read decomposePar dictionary from specified location"
    );
    argList::addOption
    (
        "targetDecomposeParDict",
        "file",
        "read decomposePar dictionary from specified location"
    );


    argList args(argc, argv);

    WarningInFunction << nl
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << nl
        << "mapFields is in a deprecated state and "
        << "is now superseded by foamMap." << nl
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

#if defined(WIN64) || defined(WIN32)
#include "include/forceLoadLibraries.H"
#endif
    if (!args.check())
    {
        FatalError.exit();
    }
     #include "include/createTime.H"
   Info<< "Time starts = "
        << runTime.cpuTimeIncrement()
        << " s\n" << endl << endl;
    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());

    fileName casePath = args[1];
    const fileName rootDirSource = casePath.path().toAbsolute();
    const fileName caseDirSource = casePath.name();

    const bool parallelSource = args.optionFound("parallelSource");
    const bool parallelTarget = args.optionFound("parallelTarget");
    const bool consistent = args.optionFound("consistent");

    meshToMesh0::order mapOrder = meshToMesh0::INTERPOLATE;
    if (args.optionFound("mapMethod"))
    {
        const word mapMethod(args["mapMethod"]);
        if (mapMethod == "mapNearest")
        {
            mapOrder = meshToMesh0::MAP;
        }
        else if (mapMethod == "interpolate")
        {
            mapOrder = meshToMesh0::INTERPOLATE;
        }
        else if (mapMethod == "cellPointInterpolate")
        {
            mapOrder = meshToMesh0::CELL_POINT_INTERPOLATE;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown mapMethod " << mapMethod << ". Valid options are: "
                << "mapNearest, interpolate and cellPointInterpolate"
                << exit(FatalError);
        }

        Info<< "Mapping method: " << mapMethod << endl;
    }

    const bool subtract = args.optionFound("subtract");
    if (subtract)
    {
        Info<< "Subtracting mapped source field from target" << endl;
    }


    #include "createTimes.H"

    const bool allRegions = args.optionFound("allRegions");
    const bool regionMaps = args.optionFound("regionMaps");

    wordList sourceRegions;
    wordList targetRegions;
    wordList sourceRegionDirs;
    wordList targetRegionDirs;
    if (regionMaps)
    {
        List<Pair<word>> regionMap;
        args.optionLookup("regionMaps")() >> regionMap;
        forAll(regionMap, regioni)
        {
            sourceRegions.append(regionMap[regioni].first());
            targetRegions.append(regionMap[regioni].second());
        }
        sourceRegionDirs = sourceRegions;
        targetRegionDirs = targetRegions;
    }
    else if (allRegions)
    {
        Info<< "Mapping for all regions in target regionProperties" << nl
            << endl;
        regionProperties rp(runTimeTarget);
        forAllConstIter(HashTable<wordList>, rp, iter)
        {
            const wordList& regions = iter();
            forAll(regions, i)
            {
                if (findIndex(targetRegions, regions[i]) == -1)
                {
                    targetRegions.append(regions[i]);
                    sourceRegions.append(regions[i]);
                }
            }
        }
        sourceRegionDirs = sourceRegions;
        targetRegionDirs = targetRegions;
    }
    else
    {
        Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;

        if (args.optionFound("sourceRegion"))
        {
            sourceRegions = wordList(1, args["sourceRegion"]);
            sourceRegionDirs = sourceRegions;
            Info<< "Source region: " << sourceRegions[0] << endl;
        }
        else
        {
            sourceRegions = wordList(1, fvMesh::defaultRegion);
            sourceRegionDirs = wordList(1, word::null);
        }

        Info<< "Target: " << rootDirTarget << " " << caseDirTarget << endl;

        if (args.optionFound("targetRegion"))
        {
            targetRegions = wordList(1, args["targetRegion"]);
            targetRegionDirs = targetRegions;
            Info<< "Target region: " << targetRegions[0] << endl;
        }
        else
        {
            targetRegions = wordList(1, fvMesh::defaultRegion);
            targetRegionDirs = wordList(1, word::null);
        }
    }

    forAll(sourceRegions, regioni)
    {
        const word& sourceRegion = sourceRegions[regioni];
        const word& targetRegion = targetRegions[regioni];
        const word& targetRegionDir = targetRegionDirs[regioni];
        const word& sourceRegionDir = sourceRegionDirs[regioni];

        Info<<"Maping source region : "<< sourceRegion
            <<" to target region: "<< targetRegion << endl;

        HashTable<word> patchMap;
        wordList cuttingPatches;

        if (!consistent)
        {
            IOdictionary mapFieldsDict
            (
                IOobject
                (
                    "mapFieldsDict",
                    runTimeTarget.system(),
                    targetRegionDir,
                    runTimeTarget,
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE,
                    false
                )
             );

            mapFieldsDict.lookup("patchMap") >> patchMap;
            mapFieldsDict.lookup("cuttingPatches") >>  cuttingPatches;
        }

        if (parallelSource && !parallelTarget)
        {
            const int nProcs = readNumProcs
            (
                args,
                "sourceDecomposeParDict",
                sourceRegionDir,
                runTimeSource
            );

            Info<< "Create target mesh\n" << endl;

            fvMesh meshTarget
            (
                IOobject
                (
                    targetRegion,
                    runTimeTarget.timeName(),
                    runTimeTarget
                )
            );

            Info<< "Target mesh size: " << meshTarget.nCells() << endl;

            for (int proci=0; proci<nProcs; proci++)
            {
                Info<< nl << "Source processor " << proci << endl;

                Time runTimeSource
                (
                    Time::controlDictName,
                    rootDirSource,
                    caseDirSource/fileName(word("processor") + name(proci))
                );

                #include "setTimeIndex.H"

                fvMesh meshSource
                (
                    IOobject
                    (
                        sourceRegion,
                        runTimeSource.timeName(),
                        runTimeSource
                    )
                );

                Info<< "mesh size: " << meshSource.nCells() << endl;

                if (consistent)
                {
                    mapConsistentSubMesh
                    (
                        meshSource,
                        meshTarget,
                        mapOrder,
                        subtract
                    );
                }
                else
                {
                    mapSubMesh
                    (
                        meshSource,
                        meshTarget,
                        patchMap,
                        cuttingPatches,
                        mapOrder,
                        subtract
                    );
                }
            }
        }
        else if (!parallelSource && parallelTarget)
        {
            const int nProcs = readNumProcs
            (
                args,
                "targetDecomposeParDict",
                targetRegionDir,
                runTimeTarget
            );


            Info<< "Create source mesh\n" << endl;

            #include "setTimeIndex.H"

            fvMesh meshSource
            (
                IOobject
                (
                    sourceRegion,
                    runTimeSource.timeName(),
                    runTimeSource
                )
            );

            Info<< "Source mesh size: " << meshSource.nCells() << endl;

            for (int proci=0; proci<nProcs; proci++)
            {
                Info<< nl << "Target processor " << proci << endl;

                Time runTimeTarget
                (
                    Time::controlDictName,
                    rootDirTarget,
                    caseDirTarget/fileName(word("processor") + name(proci))
                );

                fvMesh meshTarget
                (
                    IOobject
                    (
                        targetRegion,
                        runTimeTarget.timeName(),
                        runTimeTarget
                    )
                );

                Info<< "mesh size: " << meshTarget.nCells() << endl;

                if (consistent)
                {
                    mapConsistentSubMesh
                    (
                        meshSource,
                        meshTarget,
                        mapOrder,
                        subtract
                    );
                }
                else
                {
                    mapSubMesh
                    (
                        meshSource,
                        meshTarget,
                        patchMap,
                        addProcessorPatches(meshTarget, cuttingPatches),
                        mapOrder,
                        subtract
                     );
                }
            }
        }
        else if (parallelSource && parallelTarget)
        {
            const int nProcsSource = readNumProcs
            (
                args,
                "sourceDecomposeParDict",
                sourceRegionDir,
                runTimeSource
            );
            const int nProcsTarget = readNumProcs
            (
                args,
                "targetDecomposeParDict",
                targetRegionDir,
                runTimeTarget
            );

            List<boundBox> bbsTarget(nProcsTarget);
            List<bool> bbsTargetSet(nProcsTarget, false);

            for (int procISource=0; procISource<nProcsSource; procISource++)
            {
                Info<< nl << "Source processor " << procISource << endl;

                Time runTimeSource
                (
                    Time::controlDictName,
                    rootDirSource,
                    caseDirSource/fileName(word("processor")+name(procISource))
                );

                #include "setTimeIndex.H"

                fvMesh meshSource
                (
                    IOobject
                    (
                        sourceRegion,
                        runTimeSource.timeName(),
                        runTimeSource
                    )
                );

                Info<< "mesh size: " << meshSource.nCells() << endl;

                boundBox bbSource(meshSource.bounds());

                for (int procITarget=0; procITarget<nProcsTarget; procITarget++)
                {
                    if
                    (
                        !bbsTargetSet[procITarget]
                        || (
                            bbsTargetSet[procITarget]
                            && bbsTarget[procITarget].overlaps(bbSource)
                            )
                    )
                    {
                        Info<< nl << "Target processor " << procITarget << endl;

                        Time runTimeTarget
                        (
                            Time::controlDictName,
                            rootDirTarget,
                            caseDirTarget/fileName(word("processor")
                                                   + name(procITarget))
                        );

                        fvMesh meshTarget
                        (
                            IOobject
                            (
                                targetRegion,
                                runTimeTarget.timeName(),
                                runTimeTarget
                            )
                        );

                        Info<< "mesh size: " << meshTarget.nCells() << endl;

                        bbsTarget[procITarget] = meshTarget.bounds();
                        bbsTargetSet[procITarget] = true;

                        if (bbsTarget[procITarget].overlaps(bbSource))
                        {
                            if (consistent)
                            {
                                mapConsistentSubMesh
                                (
                                    meshSource,
                                    meshTarget,
                                    mapOrder,
                                    subtract
                                );
                            }
                            else
                            {
                                mapSubMesh
                                (
                                    meshSource,
                                    meshTarget,
                                    patchMap,
                                    addProcessorPatches
                                    (
                                        meshTarget,
                                        cuttingPatches
                                    ),
                                    mapOrder,
                                    subtract
                                );
                            }
                        }
                    }
                }
            }
        }
        else
        {
            #include "setTimeIndex.H"

            Info<< "Create meshes\n" << endl;

            fvMesh meshSource
            (
                IOobject
                (
                    sourceRegion,
                    runTimeSource.timeName(),
                    runTimeSource
                )
            );

            fvMesh meshTarget
            (
                IOobject
                (
                    targetRegion,
                    runTimeTarget.timeName(),
                    runTimeTarget
                )
            );

            Info<< "Source mesh size: " << meshSource.nCells() << tab
                << "Target mesh size: " << meshTarget.nCells() << nl << endl;

            if (consistent)
            {
                mapConsistentMesh(meshSource, meshTarget, mapOrder, subtract);
            }
            else
            {
                mapSubMesh
                (
                    meshSource,
                    meshTarget,
                    patchMap,
                    cuttingPatches,
                    mapOrder,
                    subtract
                );
            }
        }
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;
    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
