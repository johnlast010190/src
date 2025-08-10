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
    (c) 2020 Esi Ltd.

Class

Description
    read or create modifyGIBBoundaryDict on the fly

SourceFiles
    readOrCreateDictionary.H

\*---------------------------------------------------------------------------*/

#ifndef readOrCreateDictionary_H
#define readOrCreateDictionary_H

#include "db/Time/Time.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    IOdictionary readOrCreateDictionary
    (
        const Time& runTime,
        const argList& args
    )
    {
        word dictName = "modifyGIBBoundaryDict";
        IOobject dictIO
        (
            dictName,
            runTime.caseSystem(),
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        );

        if (args.optionFound("dict"))
        {
            fileName dictPath = args["dict"];
            if (isDir(dictPath))
            {
                dictPath = dictPath / dictName;
            }
            dictIO =
                IOobject
                (
                    dictPath,
                    runTime,
                    IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                );
        }

        IOdictionary cGIBDict
        (
            dictIO,
            dictionary()
        );

        if (args.optionFound("name"))
        {
            cGIBDict.clear();

            word gibName = args["name"];
            word faceZoneName = args["faceZone"];

            word masterName = gibName;
            dictionary mDict;
            mDict.add("type", "indirectWall");
            mDict.add("faceZone", faceZoneName);
            mDict.add("indirectPolyPatchType", "master");

            word slaveName = gibName+"_slave";
            dictionary sDict;
            sDict.add("type", "indirectWall");
            sDict.add("faceZone", faceZoneName);
            sDict.add("indirectPolyPatchType", "slave");

            PtrList<entry> boundaryList(2);
            boundaryList.set
            (
                0,
                new dictionaryEntry
                (
                    masterName,
                    dictionary::null,
                    mDict
                )
            );

            boundaryList.set
            (
                1,
                new dictionaryEntry
                (
                    slaveName,
                    dictionary::null,
                    sDict
                )
            );

            cGIBDict.add("boundary", boundaryList);

            cGIBDict.add("initialization", dictionary());
            dictionary& initDict = cGIBDict.subDict("initialization");

            if (args.optionFound("triSurfaceName"))
            {
                word triSurfaceName = args["triSurfaceName"];

                initDict.add("triSurfaceName", triSurfaceName);
                initDict.add("dynamicFvMesh", "solidBodyGIBFvMesh");
                initDict.add("initializationType", "dynamicGIBFvMesh");
                if (args.optionFound("orientByPatches"))
                {
                    const wordList region0Patch =
                        args.optionRead<wordList>("orientByPatches");
                    initDict.add("region0Patch", region0Patch);
                }
                else
                {
                    initDict.add("region0Patch", wordList());
                }
            }
            else
            {
                initDict.add("initializationType", "none");
                //- add orientation
                if (args.optionFound("orientByPoint"))
                {
                    const point outside = args.optionRead<point>("orientByPoint");
                    initDict.add("orientByPoint", outside);
                }
            }

            if (args.optionFound("exportGeometry"))
            {
                word surfaceName = args["exportGeometry"];
                initDict.add("exportGeometry", surfaceName);
            }

            Info<< "GIBs are set from options:" << endl;
            Info<< cGIBDict <<endl;
        }

//      Info<< exit(FatalError);

        return cGIBDict;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
