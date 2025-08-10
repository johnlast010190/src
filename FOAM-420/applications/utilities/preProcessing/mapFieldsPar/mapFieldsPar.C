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
    (c) 2017 Esi Ltd.
    (c) 2015 OpenCFD Ltd.

Application
    mapFieldsPar

Group
    grpPreProcessingUtilities

Description
    Maps volume fields from one mesh to another, reading and
    interpolating all fields present in the time directory of both cases.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "meshToMesh/meshToMesh.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "MapMeshes.H"
#include "regionProperties/regionProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void mapConsistentMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const dictionary& dict,
    const word& mapMethod,
    const word& AMIMapMethod,
    const bool subtract,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Consistently creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    meshToMesh interp
    (
        meshSource,
        meshTarget,
        mapMethod,
        AMIMapMethod,
        true,
        dict
    );

    if (subtract)
    {
        MapMesh<minusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
    else
    {
        MapMesh<plusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
}


void mapSubMesh
(
    const fvMesh& meshSource,
    const fvMesh& meshTarget,
    const dictionary& dict,
    const HashTable<word>& patchMap,
    const wordList& cuttingPatches,
    const word& mapMethod,
    const word& AMIMapMethod,
    const bool subtract,
    const HashSet<word>& selectedFields,
    const bool noLagrangian
)
{
    Info<< nl << "Creating and mapping fields for time "
        << meshSource.time().timeName() << nl << endl;

    meshToMesh interp
    (
        meshSource,
        meshTarget,
        mapMethod,
        AMIMapMethod,
        patchMap,
        cuttingPatches,
        true,
        dict
    );

    if (subtract)
    {
        MapMesh<minusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
    else
    {
        MapMesh<plusEqOp>
        (
            interp,
            selectedFields,
            noLagrangian
        );
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "map volume fields from one mesh to another"
    );

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
        "consistent",
        "source and target geometry and boundary conditions identical"
    );
    argList::addOption
    (
        "mapMethod",
        "word",
        "specify the mapping method "
        "(direct|mapNearest|cellVolumeWeight|correctedCellVolumeWeight)"
    );
    argList::addOption
    (
        "patchMapMethod",
        "word",
        "specify the patch mapping method (direct|mapNearest|faceAreaWeight)"
    );
    argList::addBoolOption
    (
        "subtract",
        "subtract mapped source from target"
    );
    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be mapped. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );
    argList::addBoolOption
    (
        "noLagrangian",
        "skip mapping lagrangian positions and fields"
    );

    argList args(argc, argv);

    WarningInFunction << nl
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << nl
        << "mapFieldsPar is in a deprecated state and "
        << "is now superseded by foamMap." << nl
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

    #if defined(WIN64) || defined(WIN32)
    fileName rootDirTarget(args.rootPath(), true);
    fileName caseDirTarget(args.globalCaseName(), true);
    #else
    fileName rootDirTarget(args.rootPath());
    fileName caseDirTarget(args.globalCaseName());
    #endif

    #if defined(WIN64) || defined(WIN32)
    const fileName casePath(args[1],true);
    #else
    const fileName casePath = args[1];
    #endif
    const fileName rootDirSource = casePath.path();
    const fileName caseDirSource = casePath.name();

    const bool consistent = args.optionFound("consistent");


    word mapMethod = meshToMesh::interpolationMethodNames_
    [
        meshToMesh::imCellVolumeWeight
    ];

    if  (args.optionReadIfPresent("mapMethod", mapMethod))
    {
        Info<< "Mapping method: " << mapMethod << endl;
    }


    word patchMapMethod;
    if (meshToMesh::interpolationMethodNames_.hasEnum(mapMethod))
    {
        // Lookup corresponding AMI method
        meshToMesh::interpolationMethod method =
            meshToMesh::interpolationMethodNames_[mapMethod];

        patchMapMethod = AMIPatchToPatchInterpolation::interpolationMethodToWord
        (
            meshToMesh::interpolationMethodAMI(method)
        );
    }

    // Optionally override
    if (args.optionFound("patchMapMethod"))
    {
        patchMapMethod = args["patchMapMethod"];

        Info<< "Patch mapping method: " << patchMapMethod << endl;
    }


    if (patchMapMethod.empty())
    {
        FatalErrorInFunction
            << "No valid patchMapMethod for method " << mapMethod
            << ". Please supply one through the 'patchMapMethod' option"
            << exit(FatalError);
    }

    const bool subtract = args.optionFound("subtract");
    if (subtract)
    {
        Info<< "Subtracting mapped source field from target" << endl;
    }

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }

    const bool noLagrangian = args.optionFound("noLagrangian");

    #if defined( WIN32 ) || defined( WIN64 )
    //- needed to force load ptscotch and othe libs at runtime -DC
    #include "include/forceLoadLibraries.H"
    #endif

    #include "createTimes.H"

    const bool allRegions = args.optionFound("allRegions");
    const bool regionMaps = args.optionFound("regionMaps");

    wordList sourceRegions;
    wordList targetRegions;
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
        targetRegionDirs = targetRegions;
    }
    else
    {
        Info<< "Source: " << rootDirSource << " " << caseDirSource << endl;

        if (args.optionFound("sourceRegion"))
        {
            sourceRegions = wordList(1, args["sourceRegion"]);
            Info<< "Source region: " << sourceRegions[0] << endl;
        }
        else
        {
            sourceRegions = wordList(1, fvMesh::defaultRegion);
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

    #include "setTimeIndex.H"

    forAll(sourceRegions, regioni)
    {
        const word& sourceRegion = sourceRegions[regioni];
        const word& targetRegion = targetRegions[regioni];
        const word& targetRegionDir = targetRegionDirs[regioni];

        Info<<"Mapping source region : "<< sourceRegion
            <<" to target region: "<< targetRegion << endl;

        HashTable<word> patchMap;
        wordList cuttingPatches;

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

        if (!consistent)
        {
            mapFieldsDict.lookup("patchMap") >> patchMap;
            mapFieldsDict.lookup("cuttingPatches") >>  cuttingPatches;
        }


        Info<< "\nCreate meshes\n" << endl;

        fvMesh meshSource
        (
            IOobject
            (
                sourceRegion,
                runTimeSource.timeName(),
                runTimeSource,
                IOobject::MUST_READ
            )
        );

        fvMesh meshTarget
        (
            IOobject
            (
                targetRegion,
                runTimeTarget.timeName(),
                runTimeTarget,
                IOobject::MUST_READ
             )
        );

        Info<< "Source mesh size: " << meshSource.globalData().nTotalCells()
            << tab << "Target mesh size: "
            << meshTarget.globalData().nTotalCells() << nl << endl;

        if (consistent)
        {
            mapConsistentMesh
            (
                meshSource,
                meshTarget,
                mapFieldsDict,
                mapMethod,
                patchMapMethod,
                subtract,
                selectedFields,
                noLagrangian
            );
        }
        else
        {
            mapSubMesh
            (
                meshSource,
                meshTarget,
                mapFieldsDict,
                patchMap,
                cuttingPatches,
                mapMethod,
                patchMapMethod,
                subtract,
                selectedFields,
                noLagrangian
            );
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
