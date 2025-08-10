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
    (c) 2017 OpenCFD Ltd.
    (c) 2010-2016 Esi Ltd.
    (c) 2020 Esi Ltd.

Application
    surfaceMeshTriangulate

Group
    grpSurfaceUtilities

Description
    Extracts surface from a polyMesh. Depending on output surface format
    triangulates faces.

    Region numbers on faces cannot be guaranteed to be the same as the patch
    indices.

    Optionally only triangulates named patches.

    If run in parallel the processor patches get filtered out by default and
    the mesh gets merged (based on topology).

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/timeSelector.H"
#include "meshToSurface/meshToSurface.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract surface from a polyMesh"
    );
    timeSelector::addOptions();

    argList::validArgs.append("output file");
    #include "include/addRegionOption.H"
    argList::addBoolOption
    (
        "excludeProcPatches",
        "exclude processor patches"
    );
    argList::addOption
    (
        "patches",
        "(patch0 .. patchN)",
        "only triangulate selected patches (wildcards supported)"
    );
    argList::addOption
    (
        "faceZones",
        "(fz0 .. fzN)",
        "triangulate selected faceZones (wildcards supported)"
    );
    argList::addBoolOption
    (
        "consistent",
        "keep surface patches in mesh order"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    const fileName userOutFileName(args[1]);

    if (!userOutFileName.hasExt())
    {
        FatalErrorInFunction
            << "Missing extension on output name " << userOutFileName
            << exit(FatalError);
    }

    Info<< "Extracting surface from boundaryMesh ..."
        << endl << endl;

    const bool includeProcPatches =
       !(
            args.optionFound("excludeProcPatches")
         || Pstream::parRun()
        );

    const bool consistent = args.optionFound("consistent");

    if (includeProcPatches)
    {
        Info<< "Including all processor patches." << nl << endl;
    }
    else if (Pstream::parRun())
    {
        Info<< "Excluding all processor patches." << nl << endl;
    }


    Info<< "Reading mesh from time " << runTime.value() << endl;

    #include "include/createNamedPolyMesh.H"

    // User specified times
    instantList timeDirs = timeSelector::select0(runTime, args);

    forAll(timeDirs, timeIndex)
    {
        runTime.setTime(timeDirs[timeIndex], timeIndex);
        Info<< "Time [" << timeIndex << "] = " << runTime.timeName();

        fileName outFileName;
        if (timeDirs.size() == 1)
        {
            outFileName = userOutFileName;
            Info<< nl;
        }
        else
        {
            polyMesh::readUpdateState meshState = mesh.readUpdate();
            if (timeIndex && meshState == polyMesh::UNCHANGED)
            {
                Info<<"  ... no mesh change." << nl;
                continue;
            }
            else
            {
                Info<< nl;
            }

            // The filename based on the original, but with additional
            // time information. The extension was previously checked that
            // it exists
            std::string::size_type dot = userOutFileName.rfind('.');

            outFileName =
                userOutFileName.substr(0, dot) + "_"
              + Foam::name(runTime.value()) + "."
              + userOutFileName.ext();
        }

        // Create local surface from:
        // - explicitly named patches only (-patches (at your option)
        // - all patches (default in sequential mode)
        // - all non-processor patches (default in parallel mode)
        // - all non-processor patches (sequential mode, -excludeProcPatches
        //   (at your option)

        // Construct table of patches to include.
        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        labelHashSet includePatches(bMesh.size());

        if (args.optionFound("patches"))
        {
            includePatches = bMesh.patchSet
            (
                wordReList(args.optionLookup("patches")())
            );
        }
        else
        {
            forAll(bMesh, patchi)
            {
                const polyPatch& patch = bMesh[patchi];

                if (includeProcPatches || !isA<processorPolyPatch>(patch))
                {
                    includePatches.insert(patchi);
                }
            }
        }


        const faceZoneMesh& fzm = mesh.faceZones();
        labelHashSet includeFaceZones(fzm.size());

        if (args.optionFound("faceZones"))
        {
            wordReList zoneNames(args.optionLookup("faceZones")());
            const wordList allZoneNames(fzm.names());
            forAll(zoneNames, i)
            {
                const wordRe& zoneName = zoneNames[i];

                labelList zoneIDs = findStrings(zoneName, allZoneNames);

                forAll(zoneIDs, j)
                {
                    includeFaceZones.insert(zoneIDs[j]);
                }

                if (zoneIDs.empty())
                {
                    WarningInFunction
                        << "Cannot find any faceZone name matching "
                        << zoneName << endl;
                }

            }
            Info<< "Additionally triangulating faceZones "
                << UIndirectList<word>
                  (
                      allZoneNames,
                      includeFaceZones.sortedToc()
                  )
                << endl;
        }

        meshToSurface surfaceGenerator
        (
            mesh,
            includePatches,
            includeFaceZones,
            consistent
        );

        surfaceGenerator.write(outFileName);

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
