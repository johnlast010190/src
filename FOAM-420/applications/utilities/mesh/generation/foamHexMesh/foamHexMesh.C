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
    (c) 2010-2022 Esi Ltd.
    (c) 2015-2017 OpenCFD Ltd.

Application
    foamHexMesh

Group
    grpMeshGenerationUtilities

Description
    Automatic split hex mesher.
    Refines and snaps to surface before adding layers.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"
#include "db/Time/Time.H"
#include "global/etcFiles/etcFiles.H"
#include "regionProperties/regionProperties.H"
#include "foamMeshGenerator/foamMeshGenerator.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void createRegionDictDefaults(const word& regionName)
{
    if (!Pstream::parRun() || Pstream::master())
    {
        fileName ctrlName("system/controlDict");

        IFstream cs(ctrlName);
        if (!cs)
        {
           Info<< "Creating 'system/controlDict'" << endl;

           cp
           (
              findEtcFile
              (
                  "caseDicts/preProcessing/caseSetup/settings/"
                  "defaultControlDict.cfg"
              ),
              "system/controlDict"
           );
        }

        fileName schName;
        fileName solName;
        if (regionName.size())
        {
            schName = "system/" + regionName + "/fvSchemes";
            solName = "system/" + regionName + "/fvSolution";
        }
        else
        {
           schName = "system/fvSchemes";
           solName = "system/fvSolution";
        }

        IFstream scs(schName);
        if (!scs)
        {
            Info<< "Creating " << schName <<endl;
            cp
            (
                findEtcFile
                (
                    "caseDicts/preProcessing/caseSetup/settings/"
                    "defaultFvSchemes.cfg"
                ),
                schName
            );
        }

        IFstream sos(solName);
        if (!sos)
        {
            Info<< "Creating " << solName <<endl;
            cp
            (
                findEtcFile
                (
                    "caseDicts/preProcessing/caseSetup/settings/"
                    "defaultFvSolution.cfg"
                ),
                solName
            );
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "checkGeometry",
        "check all surface geometry for quality"
    );
    argList::addOption
    (
        "time",
        "scalar",
        "Output mesh to specified time"
    );
    argList::addBoolOption
    (
        "constant",
        "Output mesh to constant"
    );
    argList::addBoolOption
    (
        "withZero",
        "Output mesh to 0 folder"
    );
    Foam::argList::addOption
    (
        "regions",
        "(name1 .. nameN)",
        "specify mesh regions"
    );
    argList::addOption
    (
        "region",
        "word",
        "Create a mesh region"
    );
    argList::addBoolOption
    (
        "writeDict",
        "Write out mesh dictionary file"
    );
    argList::addBoolOption
    (
        "writeIntermediateMeshes",
        "write meshes at intermediate stages"
    );
    argList::addBoolOption
    (
        "noBaseMeshCheck",
        "perform no checks on validity of base mesh"
    );
    argList::noCheckProcessorDirectories();

#if !defined( WIN32 ) && !defined( WIN64 )
#include "include/addProfilingOption.H"
#endif

#include "include/addDictOption.H"

#include "include/addOverwriteOption.H"
    // Create argList. This will check for non-existing processor dirs.
    // Create processor directory if non-existing
    Foam::argList args
    (
        argc,
        argv,
        /*checkArgs*/ true,
        /*checkOpts*/ true,
        /*initialise*/ true,
        /*needsThread*/false
    );

    //Info<<"MPI_THREAD_MULTIPLE "
    //    << (Pstream::haveThreads() ? "true" : "false")
    //    << endl;

#if defined( WIN32 ) || defined( WIN64 )
//- needed to force load ptscotch and othe libs at runtime -DC
#include "include/forceLoadLibraries.H"
// For some reason ptscotch isn't loaded properly by forceLoadLibraries
// We therefore do it manually here
    Foam::string ptScotcchString("ptscotchDecomp");
    Foam::dlOpen(ptScotcchString);
#endif

    bool rmCollated = false;
    if (Pstream::parRun() && Pstream::master())
    {
        fileName caseDir =
            args.path()/".."/"processors"+Foam::name(Pstream::nProcs());
        caseDir.clean();
        if (!isDir(caseDir))
        {
           rmCollated = true;
           mkDir(caseDir);
        }
    }

    if (Foam::fileHandler().type() == "uncollated")
    {
        if (rmCollated)
        {
            fileName caseDir =
                args.path()/".."/"processors"+Foam::name(Pstream::nProcs());
            caseDir.clean();
            rmDir(caseDir);
        }

        if (Pstream::parRun() && !isDir(args.path()))
        {
            mkDir(args.path());
        }
        if (!args.checkRootCase())
        {
            Foam::FatalError.exit();
        }
    }

    // Mesh writing options
    const bool writeIntermediate = args.optionFound("writeIntermediateMeshes");
    bool writeConstant = true;
    bool overwrite = false;
    scalar writeTime(0);
    if (args.optionFound("overwrite"))
    {
        writeConstant = false;
        overwrite = true;
    }
    else if (args.optionFound("constant"))
    {
        writeConstant = true;
    }
    else
    {
        if (args.optionFound("withZero"))
        {
            writeConstant = false;
        }
        else if (args.optionFound("time"))
        {
            writeConstant = false;
            writeTime = args.optionRead<scalar>("time");
        }
    }

    //Only check surface geometry
    const bool checkGeometry = args.optionFound("checkGeometry");

    fileName decompDictFile;
    args.optionReadIfPresent("decomposeParDict", decompDictFile);

    bool baseCheck = true;
    if (args.optionFound("noBaseMeshCheck"))
    {
        baseCheck = false;
    }

    bool writeDict = false;
    if (args.optionFound("writeDict"))
    {
        writeDict = true;
    }

    //Optional single region mesh dictionary location
    fileName dictPath = "";

    wordList regionNames;
    wordList regionDirs;
    if (args.optionFound("regions"))
    {
        regionNames = args.optionRead<wordList>("regions");
        forAll(regionNames, regioni)
        {
           createRegionDictDefaults(regionNames[regioni]);
        }
        regionDirs = regionNames;
    }
    else
    {
        word regionName;
        if (args.optionReadIfPresent("region", regionName))
        {
            regionNames = wordList(1, regionName);
            regionDirs = regionNames;
            createRegionDictDefaults(regionName);
        }
        else
        {
            regionNames = wordList(1, polyMesh::defaultRegion);
            regionDirs = wordList(1, word::null);
            createRegionDictDefaults(word());
        }

        if (args.optionFound("dict"))
        {
           dictPath = args["dict"];
           if (isDir(dictPath))
           {
              dictPath = dictPath / "foamHexMeshDict";
           }
        }
    }

    forAll(regionNames, regionI)
    {
        #include "include/createTime.H"
        runTime.functionObjects().off();

        const word& regionName = regionNames[regionI];
        const word& regionDir = regionDirs[regionI];

        Info<< "Generating mesh region : "<<regionName << nl;

        foamMeshGenerator foamMeshGen(runTime);
        foamMeshGen.generateMesh
        (
           regionName,
           regionDir,
           dictPath,
           decompDictFile,
           baseCheck,
           overwrite,
           writeIntermediate,
           writeConstant,
           writeTime,
           checkGeometry,
           writeDict,
           true //whether to write mesh
        );
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
