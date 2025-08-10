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
    (c) 2016-2019 OpenCFD Ltd.

Application
    foamToVTK

Group
    grpPostProcessingUtilities

Description
    General OpenFOAM to VTK file writer.

    Other bits
    - Handles volFields, pointFields, surfaceScalarField, surfaceVectorField
      fields.
    - Mesh topo changes.
    - Output legacy or xml VTK format in ascii or binary.
    - Single time step writing.
    - Write subset only.
    - Optional decomposition of cells.

Usage
    \b foamToVTK [OPTION]

    Options:
      - \par -ascii
        Write VTK data in ASCII format instead of binary.

      - \par -legacy
        Write VTK data in legacy format instead of XML format

      - \par -fields \<fields\>
        Specify single or multiple fields to write (all by default)
        For example,
        \verbatim
          -fields T
          -fields '(p T U \"alpha.*\")'
        \endverbatim
        The quoting is required to avoid shell expansions and to pass the
        information as a single argument.

      - \par -surfaceFields
        Write surfaceScalarFields (e.g., phi)

      - \par -cellSet \<name\>
      - \par -cellZone \<name\>
        Restrict conversion to either the cellSet or the cellZone.

      - \par -faceSet \<name\>
      - \par -pointSet \<name\>
        Restrict conversion to the faceSet or pointSet.

      - \par -faceZones zone or zone list
        Specify single faceZone or or multiple faceZones (name or regex)
        to write

      - \par -nearCellValue
        Output cell value on patches instead of patch value itself

      - \par -no-boundary
        Suppress output for all boundary patches

      - \par -no-internal
        Suppress output for internal (volume) mesh

      - \par -no-lagrangian
        Suppress writing Lagrangian positions and fields.

      - \par -no-point-data
        Suppress conversion of pointFields. No interpolated PointData.

      - \par -with-ids
        Additional mesh id fields (cellID, procID, patchID)

      - \par -with-point-ids
        Additional pointID field for internal mesh

      - \par -poly-decomp
        Decompose polyhedral cells into tets/pyramids

      - \par -one-boundary
        Combine all patches into a single boundary file

      - \par -patches NAME | LIST
        Specify single patch or multiple patches (name or regex) to write
        For example,
        \verbatim
          -patches top
          -patches '( front \".*back\" )'
        \endverbatim

      - \par -excludePatches NAME | LIST
        Specify single or multiple patches (name or regex) not to convert.
        For example,
        \verbatim
          -excludePatches '( inlet_1 inlet_2 "proc.*")'
        \endverbatim

Note
    The mesh subset is handled by fvMeshSubsetProxy. Slight inconsistency in
    interpolation: on the internal field it interpolates the whole volField
    to the whole-mesh pointField and then selects only those values it
    needs for the subMesh (using the fvMeshSubset cellMap(), pointMap()
    functions). For the patches however it uses the
    fvMeshSubset.interpolate function to directly interpolate the
    whole-mesh values onto the subset patch.

\*---------------------------------------------------------------------------*/

// Use std::iota in place of 1912 implementation of identity()
#include <numeric>

#include "cfdTools/general/include/fvCFD.H"
#include "meshes/pointMesh/pointMesh.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"
#include "interpolation/volPointInterpolation/volPointInterpolation.H"
#include "meshes/polyMesh/zones/ZoneMesh/faceZoneMesh.H"
#include "fields/areaFields/areaFields.H"
#include "fvMeshSubsetProxy/fvMeshSubsetProxy.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"
// #include "HashOps.H"
#include "regionProperties/regionProperties.H"
#include "primitives/strings/lists/stringListOps.H"

#include "Cloud/Cloud.H"
#include "readFields.H"
#include "reportFields.H"

#include "vtk/file/foamVtmWriter.H"
#include "vtk/output/foamVtkInternalWriter.H"
#include "vtk/output/foamVtkPatchWriter.H"
#include "vtk/output/foamVtkSurfaceMeshWriter.H"
#include "vtk/output/foamVtkSurfaceFieldWriter.H"
#include "output/foamVtkWriteTopoSet.H"
#include "vtk/file/foamVtkSeriesWriter.H"

#include "writeAreaFields.H"
#include "writeDimFields.H"
#include "writeVolFields.H"
#include "writePointFields.H"

#include "memInfo/memInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

labelList getSelectedPatches
(
    const polyBoundaryMesh& patches,
    const wordRes& whitelist,
    const wordRes& blacklist
)
{
    DynamicList<label> patchIDs(patches.size());

    for (const polyPatch& pp : patches)
    {
        if (isType<emptyPolyPatch>(pp))
        {
            continue;
        }
        else if (Pstream::parRun() && bool(isA<processorPolyPatch>(pp)))
        {
            break; // No processor patches for parallel output
        }

        const word& patchName = pp.name();

        bool accept = false;

        if (whitelist.size())
        {
            const auto matched = whitelist.matched(patchName);

            accept =
            (
                matched == wordRe::LITERAL
              ? true
              : (matched == wordRe::REGEX && !blacklist.match(patchName))
            );
        }
        else
        {
            accept = !blacklist.match(patchName);
        }

        if (accept)
        {
            patchIDs.append(pp.index());
        }
    }

    return patchIDs.shrink();
}

bool foundType
(
    const IOobjectList& list,
    const word& type
)
{
    bool found(false);
    forAllConstIters(list, iter)
    {
        const IOobject* io = iter();
        if (type(io->headerClassName()))
        {
            found = true;
        }
    }
    return found;
}

//
// Process args for output options (-ascii, -legacy)
//
vtk::outputOptions getOutputOptions(const argList& args)
{
    // Default is inline ASCII xml
    vtk::outputOptions opts;

    if (args.optionFound("legacy"))
    {
        opts.legacy(true);

        if (!args.optionFound("ascii"))
        {
            if (sizeof(floatScalar) != 4 || sizeof(label) != 4)
            {
                opts.ascii(true);

                WarningInFunction
                    << "Using ASCII rather than legacy binary VTK format since "
                    << "floatScalar and/or label are not 4 bytes in size."
                    << nl << endl;
            }
            else
            {
                opts.ascii(false);
            }
        }
    }
    else
    {
        opts.ascii(args.optionFound("ascii"));
    }

    return opts;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "General OpenFOAM to VTK file writer"
    );
    timeSelector::addOptions();

    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of binary"
    );
    argList::addBoolOption
    (
        "legacy",
        "Write legacy format instead of xml"
    );
    argList::addBoolOption
    (
        "poly-decomp",
        "Decompose polyhedral cells into tets/pyramids"
    );
    // argList::ignoreOptionCompat
    // (
    //     {"xml", 1806},  // xml is now default, use -legacy to toggle
    //     false           // bool option, no argument
    // );
    // argList::ignoreOptionCompat
    // (
    //     {"poly", 1806}, // poly is now default, use -poly-decomp to toggle
    //     false           // bool option, no argument
    // );

    argList::addOption
    (
        "cellSet",
        "name",
        "Convert mesh subset corresponding to specified cellSet"
    );
    argList::addOption
    (
        "cellZone",
        "name",
        "Convert mesh subset corresponding to specified cellZone"
    );
    argList::addOption
    (
        "faceSet",
        "name",
        "Convert specified faceSet only"
    );
    argList::addOption
    (
        "pointSet",
        "name",
        "Convert specified pointSet only"
    );
    argList::addOption
    (
        "faceZones",
        "wordRes",
        "Specify single or multiple faceZones to write\n"
        "Eg, 'cells' or '( slice \"mfp-.*\" )'."
    );

    argList::addOption
    (
        "fields",
        "wordRes",
        "Specify single or multiple fields to write (all by default)\n"
        "Eg, 'T' or '(p T U \"alpha.*\")'"
    );

    argList::addBoolOption
    (
        "processor-fields",
        "Write field values on processor boundaries only"
    );
    argList::addBoolOption
    (
        "surfaceFields",
        "Write surfaceScalarFields (eg, phi)"
    );
    argList::addBoolOption
    (
        "finiteAreaFields",
        "Write finite area fields"
    );
    argList::addBoolOption
    (
        "nearCellValue",
        "Use cell value on patches instead of patch value itself"
    );
    argList::addBoolOption
    (
        "no-boundary",
        "Suppress output for boundary patches"
    );
    argList::addBoolOption
    (
        "no-internal",
        "Suppress output for internal volume mesh"
    );
    argList::addBoolOption
    (
        "no-lagrangian",  // noLagrangian
        "Suppress writing lagrangian positions and fields"
    );
    // argList::addOptionCompat("no-lagrangian", {"noLagrangian", 1806});

    argList::addBoolOption
    (
        "no-point-data",  // noPointValues
        "Suppress conversion of pointFields. No interpolated PointData."
    );
    // argList::addOptionCompat("no-point-data", {"noPointValues", 1806});

    argList::addBoolOption
    (
        "with-ids",
        "Additional mesh id fields (cellID, procID, patchID)"
    );

    argList::addBoolOption
    (
        "with-point-ids",
        "Additional pointID field for internal mesh"
    );

    argList::addBoolOption
    (
        "one-boundary",  // allPatches
        "Combine all patches into a single file"
    );
    // argList::addOptionCompat("one-boundary", {"allPatches", 1806});

    #include "include/addRegionOption.H"

    argList::addOption
    (
        "regions",
        "wordRes",
        "Operate on selected regions from regionProperties.\n"
        "Eg, '( gas \"solid.*\" )'"
    );
    argList::addBoolOption
    (
        "allRegions",
        "Operate on all regions in regionProperties"
    );

    argList::addOption
    (
        "patches",
        "wordRes",
        "Specify single patch or multiple patches to write\n"
        "Eg, 'top' or '( front \".*back\" )'"
    );
    argList::addOption
    (
        "excludePatches",
        "wordRes",
        "Specify single patch or multiple patches to exclude from writing."
        " Eg, 'outlet' or '( inlet \".*Wall\" )'"
    );
    // argList::ignoreOptionCompat
    // (
    //     {"noFaceZones", 1806},  // faceZones are only enabled on demand
    //     false                   // bool option, no argument
    // );
    // argList::ignoreOptionCompat
    // (
    //     {"noLinks", 1806},      // ignore never make any links
    //     false                   // bool option, no argument
    // );
    // argList::ignoreOptionCompat
    // (
    //     {"useTimeName", 1806},  // ignore - works poorly with VTM formats
    //     false                   // bool option, no argument
    // );
    argList::addBoolOption
    (
        "overwrite",
        "Remove any existing VTK output directory"
    );
    argList::addOption
    (
        "name",
        "subdir",
        "Directory name for VTK output (default: 'VTK')"
    );

    #include "include/setRootCase.H"

    const bool decomposePoly = args.optionFound("poly-decomp");
    const bool doBoundary    = !args.optionFound("no-boundary");
    const bool doInternal    = !args.optionFound("no-internal");
    const bool doLagrangian  = !args.optionFound("no-lagrangian");
    const bool doFiniteArea  = args.optionFound("finiteAreaFields");
    const bool doSurfaceFields = args.optionFound("surfaceFields");

    const bool oneBoundary   = args.optionFound("one-boundary") && doBoundary;
    const bool nearCellValue = args.optionFound("nearCellValue") && doBoundary;
    const bool allRegions    = args.optionFound("allRegions");

    const vtk::outputOptions writeOpts = getOutputOptions(args);

    bool processorFieldsOnly = false;

    if (args.optionFound("processor-fields"))
    {
        if (!Pstream::parRun())
        {
            Info<< "Ignoring processor patch writing in serial"
                << nl << endl;
        }
        else if (writeOpts.legacy())
        {
            Info<< "Ignoring processor patch writing in legacy format"
                << nl << endl;
        }
        else
        {
            processorFieldsOnly = true;

            Info<< "Writing processor patch fields only"
                << nl << endl;
        }
    }

    if (nearCellValue)
    {
        Info<< "Using neighbouring cell value instead of patch value"
            << nl << endl;
    }

    const bool noPointValues = args.optionFound("no-point-data");
    if (noPointValues)
    {
        Info<< "Point fields and interpolated point data"
            << " disabled with the '-no-point-data' option"
            << nl;
    }

    const bool withPointIds = args.optionFound("with-point-ids");
    if (withPointIds)
    {
        Info<< "Write point ids requested";

        if (noPointValues)
        {
            Info<< ", but ignored due to the '-no-point-data' option";
        }

        Info<< nl;
    }

    const bool withMeshIds = args.optionFound("with-ids");
    if (withMeshIds)
    {
        Info<< "Writing mesh ids (cell, patch, proc) requested" << nl;
    }

    wordRes includePatches, excludePatches;
    if (doBoundary)
    {
        if (args.optionFound("patches"))
        {
            includePatches = args.optionReadList<wordRe>("patches");
            Info<< "Including patches " << flatOutput(includePatches)
                << nl << endl;
        }

        if (args.optionFound("excludePatches"))
        {
            excludePatches = args.optionReadList<wordRe>("excludePatches");
            Info<< "Excluding patches " << flatOutput(excludePatches)
                << nl << endl;
        }
    }

    // Can be specified as empty (ie, no fields)
    wordRes selectedFields;
    const bool useFieldFilter = args.optionFound("fields");
    if (useFieldFilter)
    {
        selectedFields = args.optionReadList<wordRe>("fields");
    }

    // Non-mandatory
    List<wordRe> temp;
    if (args.optionFound("faceZones"))
    {
        temp = args.optionReadList<wordRe>("faceZones");
    }
    const wordRes selectedFaceZones(temp);

    #include "include/createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    // Information for file series
    HashTable<vtk::seriesWriter, fileName> vtkSeries;

    wordList regionNames;
    wordRes selectRegions;
    if (allRegions)
    {
        regionNames =
            regionProperties(runTime).groupNames();
            // regionProperties(runTime, IOobject::READ_IF_PRESENT).names();

        if (regionNames.empty())
        {
            Info<< "Warning: "
                << "No regionProperties - assuming default region"
                << nl << endl;

            regionNames.resize(1);
            regionNames.first() = fvMesh::defaultRegion;
        }
        else
        {
            Info<< "Using all regions in regionProperties" << nl
                << "    "<< flatOutput(regionNames) << nl;
        }
    }
    else if (args.optionFound("regions"))
    {
        selectRegions = args.optionReadList<wordRe>("regions");

        if (selectRegions.empty())
        {
            regionNames.resize(1);
            regionNames.first() = fvMesh::defaultRegion;
        }
        else if
        (
            selectRegions.size() == 1 && !selectRegions.first().isPattern()
        )
        {
            // Identical to -region NAME
            regionNames.resize(1);
            regionNames.first() = selectRegions.first();
        }
        else
        {
            regionNames =
                regionProperties(runTime).groupNames();
                // regionProperties(runTime, IOobject::READ_IF_PRESENT).names();

            if (regionNames.empty())
            {
                Info<< "Warning: "
                    << "No regionProperties - assuming default region"
                    << nl << endl;

                regionNames.resize(1);
                regionNames.first() = fvMesh::defaultRegion;
            }
            else
            {
                inplaceSubsetStrings(selectRegions, regionNames);

                if (regionNames.empty())
                {
                    Info<< "No matching regions ... stopping" << nl << endl;
                    return 1;
                }

                Info<< "Using matching regions: "
                    << flatOutput(regionNames) << nl;
            }
        }
    }
    else
    {
        regionNames.resize(1);
        regionNames.first() = args.optionLookupOrDefault<word>("region", fvMesh::defaultRegion);
    }


    // Names for sets and zones
    word cellSelectionName;
    word faceSetName;
    word pointSetName;

    fvMeshSubsetProxy::subsetType cellSubsetType = fvMeshSubsetProxy::NONE;

    string vtkName = args.globalCaseName();

    if (regionNames.size() == 1)
    {
        if (args.optionReadIfPresent("cellSet", cellSelectionName))
        {
            vtkName = cellSelectionName;
            cellSubsetType = fvMeshSubsetProxy::SET;

            Info<< "Converting cellSet " << cellSelectionName
                << " only. New outside faces as \"oldInternalFaces\"."
                << nl;
        }
        else if (args.optionReadIfPresent("cellZone", cellSelectionName))
        {
            vtkName = cellSelectionName;
            cellSubsetType = fvMeshSubsetProxy::ZONE;

            Info<< "Converting cellZone " << cellSelectionName
                << " only. New outside faces as \"oldInternalFaces\"."
                << nl;
        }

        args.optionReadIfPresent("faceSet", faceSetName);
        args.optionReadIfPresent("pointSet", pointSetName);
    }
    else
    {
        const wordList typeSetZones
        (
            { "cellSet", "cellZone", "faceSet", "pointSet"}
        );
        for
        (
            const word& opt : typeSetZones
        )
        {
            if (args.optionFound(opt))
            {
                Info<< "Ignoring -" << opt << " for multi-regions" << nl;
            }
        }
    }


    cpuTime timer;
    memInfo mem;
    Info<< "Initial memory " << mem.update().size() << " kB" << endl;

    #include "createMeshes.H"

    // Directory management

    // Sub-directory for output
    const word vtkDirName = args.optionLookupOrDefault<word>("name", "VTK");

    const fileName outputDir(args.rootPath()/args.globalCaseName()/vtkDirName);

    if (Pstream::master())
    {
        // Overwrite or create the VTK/regionName directories.
        // For the default region, this is simply "VTK/"

        fileName regionDir;
        for (const word& regionName : regionNames)
        {
            if (regionName != polyMesh::defaultRegion)
            {
                regionDir = outputDir / regionName;
            }
            else
            {
                regionDir = outputDir;
            }

            if (args.optionFound("overwrite") && isDir(regionDir))
            {
                Info<< "Removing old directory "
                    << regionDir
                    << nl << endl;
                rmDir(regionDir);
            }
            mkDir(regionDir);
        }
    }


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        const word timeDesc = "_" + Foam::name(runTime.timeIndex());
        const scalar timeValue = runTime.value();

        Info<< "Time: " << runTime.timeName() << endl;


        // Accumulate information for multi-region VTM
        vtk::vtmWriter vtmMultiRegion;

        // vtmMultiRegion.set(vtkDir/vtkName + timeDesc)

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];

            fileName regionPrefix;
            if (regionName != polyMesh::defaultRegion)
            {
                regionPrefix = regionName;
            }

            auto& meshProxy = meshProxies[regioni];
            auto& vtuMeshCells = vtuMappings[regioni];

            // polyMesh::readUpdateState meshState = mesh.readUpdate();

            // Check for new polyMesh/ and update mesh, fvMeshSubset
            // and cell decomposition.
            polyMesh::readUpdateState meshState =
                meshProxy.readUpdate();

            const fvMesh& mesh = meshProxy.mesh();

            if
            (
                meshState == polyMesh::TOPO_CHANGE
             || meshState == polyMesh::TOPO_PATCH_CHANGE
            )
            {
                // Trigger change for vtk cells too
                vtuMeshCells.clear();
            }

            // Write topoSets before attempting anything else
            {
                #include "convertTopoSet.H"
                if (wroteTopoSet)
                {
                    continue;
                }
            }

            // Search for list of objects for this time
            IOobjectList objects(meshProxy.baseMesh(), runTime.timeName());

            if (useFieldFilter)
            {
                objects.filterKeys(selectedFields);
            }

            // Prune restart fields
            objects.filterKeys
            (
                [](const word& k){ return k.endsWith("_0"); },
                true  // prune
            );

            const hashedWordList pointFieldTypes
            {
                pointScalarField::typeName,
                pointVectorField::typeName,
                pointSphericalTensorField::typeName,
                pointSymmTensorField::typeName,
                pointTensorField::typeName,
                pointScalarField::Internal::typeName,
                pointVectorField::Internal::typeName,
                pointSphericalTensorField::Internal::typeName,
                pointSymmTensorField::Internal::typeName,
                pointTensorField::Internal::typeName
            };

            if (noPointValues)
            {
                // Prune point fields unless specifically requested
                for (auto fieldType : pointFieldTypes)
                {
                    objects.filterKeys
                    (
                        fieldType,
                        true  // prune
                    );
                }
            }

            if (processorFieldsOnly)
            {
                // Processor-patches only and continue
                #include "convertProcessorPatches.H"
                continue;
            }

            // Volume, internal, point fields
            #include "convertVolumeFields.H"

            // Surface fields
            #include "convertSurfaceFields.H"

            // Finite-area mesh and fields - need not exist
            #include "convertAreaFields.H"

            // Write lagrangian data
            #include "convertLagrangian.H"
        }

        // Emit multi-region vtm
        if (Pstream::master() && regionNames.size() > 1)
        {
            fileName outputName
            (
                outputDir/vtkName + "-regions" + timeDesc + ".vtm"
            );

            vtmMultiRegion.setTime(timeValue);
            vtmMultiRegion.write(outputName);

            fileName seriesName(vtk::seriesWriter::base(outputName));

            vtk::seriesWriter& series = vtkSeries(seriesName);

            // First time?
            // Load from file, verify against filesystem,
            // prune time >= currentTime
            if (series.empty())
            {
                series.load(seriesName, true, timeValue);
            }

            series.append(timeValue, outputName);
            series.write(seriesName);
        }

        Info<< "Wrote in "
            << timer.cpuTimeIncrement() << " s, "
            << mem.update().size() << " kB" << endl;
    }


    Info<< "\nEnd: "
        << timer.elapsedCpuTime() << " s, "
        << mem.update().peak() << " kB (peak)\n" << endl;

    return 0;
}


// ************************************************************************* //
