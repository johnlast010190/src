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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2009 OpenCFD Ltd.

Application
    foamMeshSetup

Description
    Geometric population of a foamHexMeshDict.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "primitives/strings/wordRes/wordRes.H"
#include "primitives/strings/lists/stringListOps.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "triSurface/triSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("overwrite", "");
    Foam::argList args(argc, argv);

#include "include/createTime.H"
    runTime.functionObjects().off();

    const bool overwrite = args.optionFound("overwrite");

    fileName baseDir
    (
        "${FOAM_PROJECT_DIR}/etc/caseDicts/preProcessing/createCase"
    );

    if (isFile(runTime.time().system()/"foamHexMeshDict") && overwrite)
    {
        //Make a backup of foamHexMeshDict
        cp
        (
            runTime.time().system()/"foamHexMeshDict",
            runTime.time().system()/"foamHexMeshDict.backup"
        );
        // copy a templated foamHexMeshDict

        cp
        (
            baseDir + "/createCase.foamHexMeshDict",
            runTime.time().system()/"foamHexMeshDict"
        );
    }
    else if (!isFile(runTime.time().system()/"foamHexMeshDict") || overwrite)
    {
        //foamHexMeshDict
        cp
        (
            baseDir + "/createCase.foamHexMeshDict",
            runTime.time().system()/"foamHexMeshDict"
        );
    }

    // Read meshing dictionary
    IOdictionary meshDict
    (
       IOobject
       (
            "foamHexMeshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       )
    );

    IOobject setupHeader
    (
        "foamMeshSetupDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    List<wordRe> patterns(0);
    Switch allRegions = false;

    if (setupHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary setupDict
        (
            IOobject
            (
                "foamMeshSetupDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
         );

        List<string> strList =
            setupDict.lookupOrDefault("volRefine", List<string>());
        patterns.setSize(strList.size());
        forAll(patterns, i)
        {
            patterns[i] = wordRe(strList[i], wordRe::REGEX);
        }

        allRegions = setupDict.lookupOrDefault("setAllPatches", false);
    }
    wordRes patternMatch(patterns);

    fileName triSurfaceDir = runTime.time().constant()/"triSurface";
    fileNameList triFiles(readDir(triSurfaceDir, fileName::FILE));

    // all surface geometry
    dictionary& geometryDict = meshDict.subDict("geometry");

    // refinement parameters
    dictionary& refineDict = meshDict.subDict("castellatedMeshControls");
    dictionary& refineSurfDict = refineDict.subDict("refinementSurfaces");
    dictionary& refineRegionsDict = refineDict.subDict("refinementRegions");

    // layer addition parameters
    dictionary& layerDict = meshDict.subDict("addLayersControls");
    dictionary& layerSpecificationDict = layerDict.subDict("layers");

    forAll(triFiles, i)
    {
        fileName name = triFiles[i];

        word ext = name.ext();

        if
        (
            ext == "stl" || ext == "nas" || ext == "obj"
            || ext == "ebs" || ext == "bts"
        )
        {
            fileName surf = name.lessExt();

            //If setting all region data read triSurface
            wordList regions(0);
            if (allRegions)
            {
                triSurface triSurf(triSurfaceDir/name);
                regions.setSize(triSurf.patches().size());
                forAll(regions, i)
                {
                    regions[i] = triSurf.patches()[i].name();
                }
            }

            if (!geometryDict.found(name))
            {
                geometryDict.add(word(name), dictionary(), false);

                char defaultType[] = "triSurfaceMesh";
                geometryDict.subDict(name).add
                (
                    word("type"),
                    defaultType,
                    false
                );

                geometryDict.subDict(name).add
                (
                    word("name"),
                    word(surf),
                    false
                );

                bool volSurface = false;

                if (patternMatch.match(string(name)))
                {
                    volSurface = true;
                    if (!refineRegionsDict.found(surf))
                    {
                        Info<<" Adding volume refinement for surface: "
                            <<name<<endl;

                        refineRegionsDict.add(word(surf), dictionary(), false);
                        char defaultMode[] = "inside";
                        refineRegionsDict.subDict(surf).add
                        (
                            word("mode"),
                            defaultMode,
                            false
                         );
                        char defaultLevel[] = "((1E15 0))";
                        refineRegionsDict.subDict(surf).add
                        (
                            word("levels"),
                            defaultLevel,
                            false
                         );
                    }
                }
                else
                {
                    if (!refineSurfDict.found(surf))
                    {
                        Info<<" Adding surface refinement for surface: "
                            << name <<endl;

                        refineSurfDict.add(word(surf), dictionary(), false);

                        dictionary& surfDict = refineSurfDict.subDict(surf);

                        char defaultLev[] = "(0 0)";

                        surfDict.add
                        (
                            word("level"),
                            defaultLev,
                            false
                         );
                        surfDict.add
                        (
                            word("regions"),
                            dictionary(),
                            false
                        );
                        if (allRegions && regions.size() > 1)
                        {
                            forAll(regions, i)
                            {
                                surfDict.subDict("regions").add
                                (
                                    regions[i],
                                    dictionary(),
                                    false
                                );
                                dictionary& regionDict =
                                    surfDict.subDict("regions");

                                regionDict.subDict(regions[i]).add
                                (
                                    word("level"),
                                    defaultLev,
                                    false
                                );
                            }
                        }

                    }
                }
                if (!volSurface)
                {
                    if (allRegions)
                    {
                        forAll(regions, i)
                        {
                            word regionPatch
                            (
                                regions.size() == 1
                                ? surf : surf + "_" + regions[i]
                            );
                            if (!layerSpecificationDict.found(regionPatch))
                            {
                                layerSpecificationDict.add
                                (
                                    regionPatch,
                                    dictionary(),
                                    false
                                );
                                char defaultLayers[] = "0";
                                layerSpecificationDict.subDict(regionPatch).add
                                (
                                    word("nSurfaceLayers"),
                                    defaultLayers,
                                    false
                                );
                            }
                        }
                    }
                    else
                    {
                        fileName wildName(surf + ".*");
                        if (!layerSpecificationDict.found(wildName))
                        {
                            layerSpecificationDict.add
                            (
                                wildName,
                                dictionary(),
                                false
                            );
                            char defaultLayers[] = "0";
                            layerSpecificationDict.subDict(wildName).add
                            (
                                word("nSurfaceLayers"),
                                defaultLayers,
                                false
                            );
                        }
                    }
                }
            }
        }
    }

    //Remove dictionary entries if triSurface removed
    wordList geomNames = geometryDict.toc();
    HashSet<fileName> allTriFiles;

    forAll(triFiles, i)
    {
        allTriFiles.insert(triFiles[i]);
    }

    forAll(geomNames, i)
    {
        const dictionary& geom = geometryDict.subDict(geomNames[i]);
        const word type = geom.lookup("type");
        if (type == "triSurfaceMesh")
        {
            if (!allTriFiles.found(geomNames[i]))
            {
                Info<<"Removing entries for surface: "
                    << geomNames[i]
                    <<" which no longer exists in triSurface folder"<<endl;

                //remove geometry entry
                word name = geom.lookup("name");

                geometryDict.remove(geomNames[i]);
                if (refineRegionsDict.found(name))
                {
                    refineRegionsDict.remove(name);
                }

                if (refineSurfDict.found(name))
                {
                    refineSurfDict.remove(name);
                }

                const List<keyType> wildCards =
                    layerSpecificationDict.keys(true);
                forAll(wildCards, i)
                {
                    regExp re(wildCards[i]);
                    if (re.match(name+"_"))
                    {
                        layerSpecificationDict.remove(wildCards[i]);
                    }
                }

                const List<keyType> nonWildCards =
                    layerSpecificationDict.keys(false);
                word wildName(name + ".*");
                forAll(nonWildCards, i)
                {
                    regExp re(wildName);
                    if (re.match(nonWildCards[i]))
                    {
                        layerSpecificationDict.remove(nonWildCards[i]);
                    }
                }
            }
        }
    }

    meshDict.regIOobject::writeObject
    (
        IOstream::ASCII,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED,
        true
    );

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
