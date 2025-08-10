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
    (c) 2010-2012 Esi Ltd.

Application
    surfaceSplitByPatch

Group
    grpSurfaceUtilities

Description
    Writes regions of triSurface to separate files.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "triSurface/triSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "write surface mesh regions to separate files"
    );

    argList::noParallel();
    argList::validArgs.append("input surfaceFile");
    argList::validOptions.insert("patches","regions to extract");
    argList::validOptions.insert("outputFile","output file name");
    argList::validOptions.insert("extension","output file extension");
    argList::validOptions.insert("outputDir","output file directory");
    argList args(argc, argv);

    const fileName surfName = args[1];

    Info<< "Reading surf from " << surfName << " ..." << nl << endl;

    fileName surfBase = surfName.lessExt();

    word extension = surfName.ext();
    if (args.optionFound("extension"))
    {
        extension = args.option("extension");
    }

    fileName outputDir("./");
    if (args.optionFound("outputDir"))
    {
        outputDir = args.option("outputDir");
    }

    triSurface surf(surfName);

    HashSet<word> includePatches;
    if (args.optionFound("patches"))
    {
        args.optionLookup("patches")() >> includePatches;

        Info<< "Finding regions " << includePatches << nl << endl;
    }

    word outFileName;
    if (args.optionFound("outputFile"))
    {
        outFileName = args.option("outputFile");
    }
    else
    {
        Info<< "Writing regions to separate files ..."
            << nl << endl;
    }

    const geometricSurfacePatchList& patches = surf.patches();

    boolList includeMap(surf.size(), false);

    forAll(patches, patchi)
    {
        const geometricSurfacePatch& pp = patches[patchi];

        word patchName = pp.name();

        if (patchName.empty())
        {
            patchName = "patch" + Foam::name(patchi);
        }

        if (!args.optionFound("patches") || includePatches.found(patchName))
        {
            // Collect faces of region
            if (!args.optionFound("outputFile"))
            {
                includeMap = false;
            }


            forAll(surf, facei)
            {
                const labelledTri& f = surf[facei];

                if (f.region() == patchi)
                {
                    includeMap[facei] = true;
                }
            }

            if (!args.optionFound("outputFile"))
            {
                fileName outFile
                (
                    outputDir + '/' + surfBase + '_' + patchName + '.' + extension
                );

                // Subset triSurface
                labelList pointMap;
                labelList faceMap;

                triSurface subSurf
                (
                    surf.subsetMesh
                    (
                        includeMap,
                        pointMap,
                        faceMap
                    )
                );

                Info<< "   Writing patch " << patchName
                    << " to file " << outFile
                    << endl;

                subSurf.write(outFile);
            }
        }
    }

    if (args.optionFound("outputFile"))
    {
        fileName outFile
        (
           outFileName
        );

        // Subset triSurface
        labelList pointMap;
        labelList faceMap;

        triSurface subSurf
        (
            surf.subsetMesh
            (
                includeMap,
                pointMap,
                faceMap
            )
        );

        Info<< "   Writing patches"
            << " to file " << outFile
            << endl;

        subSurf.write(outFile);
    }
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
