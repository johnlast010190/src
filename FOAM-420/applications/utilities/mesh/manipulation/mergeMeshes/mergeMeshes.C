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
    (c) 2016 OpenCFD Ltd.
    (c) 2017-2023 Esi Ltd.

Application
    mergeMeshes

Group
    grpMeshManipulationUtilities

Description
    Merges two meshes.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "mergePolyMesh.H"
#include "sets/topoSets/topoSet.H"
#include "regionProperties/regionProperties.H"
#include "processorMeshes.H"

using namespace Foam;

void getRootCase(fileName& casePath)
{
    casePath.clean();

    if (casePath.empty() || casePath == ".")
    {
        // handle degenerate form and '.'
        casePath = cwd();
    }
    else if (casePath[0] != '/' && casePath.name() == "..")
    {
        // avoid relative cases ending in '..' - makes for very ugly names
        casePath = cwd()/casePath;
        casePath.clean();
    }

    #if defined( WIN32 ) || defined( WIN64 )
    //- convert \ slash to / for windows only
    for (label i=0; i < casePath.size(); i++)
    {
        if (casePath[i] == '\\')
            casePath[i] ='/';
    }
    #endif
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "merge two meshes"
    );

    #include "include/addOverwriteOption.H"

    argList::addBoolOption
    (
        "allRegionsToMesh",
        "Merge all regions"
    );

    argList::validArgs.append("masterCase");

    argList::addOption
    (
        "masterRegion",
        "name"
    );

    argList::validArgs.append("addCase");
    argList::addOption
    (
        "addRegion",
        "name",
        "specify alternative mesh region for the additional mesh"
    );

    argList args(argc, argv);
    if (!args.check())
    {
         FatalError.exit();
    }

    const bool overwrite = args.optionFound("overwrite");

    fileName masterCase = args[1];
    word masterRegion = polyMesh::defaultRegion;
    args.optionReadIfPresent("masterRegion", masterRegion);

    fileName addCase = args[2];
    word addRegion = polyMesh::defaultRegion;
    args.optionReadIfPresent("addRegion", addRegion);

    // Since we don't use argList processor directory detection, add it to
    // the casename ourselves so it triggers the logic inside TimePath.
    const fileName& cName = args.caseName();
    std::string::size_type pos = cName.find("processor");
    if (pos != string::npos && pos != 0)
    {
        fileName processorName = cName.substr(pos, cName.size()-pos);
        masterCase += '/' + processorName;
        addCase += '/' + processorName;
    }


    getRootCase(masterCase);
    getRootCase(addCase);

    Info<< "Master:      " << masterCase << "  region " << masterRegion << nl
        << "mesh to add: " << addCase    << "  region " << addRegion << endl;

    #include "createTimes.H"

    regionProperties rp(runTimeMaster);

    bool allRegionsToMesh = false;
    if (args.optionFound("allRegionsToMesh"))
    {
        if (rp.found())
        {
            allRegionsToMesh = true;
        }
        else
        {
            WarningInFunction
                << "'allRegionsToMesh' is set but region properties file "
                << "does not exist. Deactivating 'allRegionsToMesh'" << endl;
        }
    }

    Info<< "Reading master mesh for time = " << runTimeMaster.timeName() << nl;

    Info<< "Create mesh\n" << endl;

    autoPtr<mergePolyMesh> masterMeshPtr;

    if (allRegionsToMesh)
    {
        masterMeshPtr.reset
        (
            new mergePolyMesh
            (
                IOobject
                (
                    masterRegion,
                    runTimeMaster.constant(),
                    runTimeMaster
                ),
                pointField(0),
                faceList(0),
                cellList(0)
            )
        );
    }
    else
    {
        masterMeshPtr.reset
        (
            new mergePolyMesh
            (
                IOobject
                (
                    masterRegion,
                    runTimeMaster.timeName(),
                    runTimeMaster
                 )
             )
        );
    }

    mergePolyMesh& masterMesh = masterMeshPtr();
    const word oldInstance = masterMesh.pointsInstance();

    wordList regionNames;
    wordList regionDirs;
    if (allRegionsToMesh)
    {
        Info<< "merging all regions in regionProperties" << nl
            << endl;
        forAllConstIter(HashTable<wordList>, rp, iter)
        {
            const wordList& regions = iter();
            forAll(regions, i)
            {
                if (findIndex(regionNames, regions[i]) == -1)
                {
                    regionNames.append(regions[i]);
                }
            }
        }
        regionDirs = regionNames;
    }
    else
    {
        regionNames = wordList(1, addRegion);
        regionDirs = regionNames;
    }

    forAll(regionNames, regionI)
    {
        const word& regionName = regionNames[regionI];

        Info<< "Reading mesh region : "<<regionName
            <<" to add for time = " << runTimeToAdd.timeName() << nl;

        Info<< "Create mesh\n" << endl;
        polyMesh meshToAdd
        (
            IOobject
            (
                regionName,
                runTimeToAdd.timeName(),
                runTimeToAdd
            )
        );

        masterMesh.addMesh(meshToAdd);
        masterMesh.merge();

        if (regionI < regionNames.size()-1)
        {
            masterMesh.updateMeshMod();
        }
    }

    if (!overwrite)
    {
        runTimeMaster++;
        Info<< "Writing combined mesh to " << runTimeMaster.timeName() << endl;
    }
    else
    {
        masterMesh.setInstance(oldInstance);
        Info<< "Writing combined mesh to " << oldInstance << endl;
    }

    masterMesh.write();
    topoSet::removeFiles(masterMesh);
    processorMeshes::removeFiles(masterMesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
