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
    (c) 2010-2023 Esi Ltd
    (c) 1991-2008 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "stateIndex.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::stateIndex::regionProperties()
{
    dictionary& constDict(global_.constant());
    constDict.add("regionProperties", dictionary());
    dictionary& regionProperties(constDict.subDict("regionProperties"));

    HashTable<wordList> regionPropertiesRegions(1);

    forAll(regions_, rI)
    {
        bool foundType = false;
        forAllIter(HashTable<wordList>, regionPropertiesRegions, iter)
        {
            if
            (
                stateFunction::regionTypeNames_[word(iter.key())]
             == regions_[rI].region()
            )
            {
                foundType = true;
                wordList& typeNames(iter());
                typeNames.setSize(typeNames.size()+1);
                typeNames[typeNames.size() - 1] = regions_[rI].regionName();
            }
        }

        if (!foundType)
        {
            word regTypeName
            (
                stateFunction::regionTypeNames_[regions_[rI].region()]
            );

            regionPropertiesRegions.insert
            (
                regTypeName,
                wordList(1, regions_[rI].regionName())
            );
        }
    }

    regionProperties.add("regions", regionPropertiesRegions);
    if (solutionSchedule_.active())
    {
        regionProperties.add
        (
            "solutionSchedule",
            solutionSchedule_.getSchedulerDictionary()
        );
        regionProperties.add
        (
            "solutionMesh",
            solutionSchedule_.getSolutionMeshDictionary()
        );
    }



}

Foam::Xfer<Foam::dictionary> Foam::stateIndex::assembleRegionInput
(
    const dictionary& input,
    const word& regionName,
    const word& rgID
) const
{
    //regionDefaults and regExp regions
    dictionary regionInput;

    //- type merge
    if (input.found("regions"))
    {
        const dictionary& regions(input.subDict("regions"));

        // List of pattern entries followed by explicit ones
        List<keyType> keys = regions.keys(true);
        keys.append(regions.keys(false));

        forAll(keys, i)
        {
            if (keys[i](rgID))
            {
                regionInput.merge
                (
                    regions.subDict(keys[i])
                );
            }
        }

        forAll(keys, i)
        {
            if (keys[i](regionName))
            {
                regionInput.merge
                (
                    regions.subDict(keys[i])
                );
            }
        }
    }

    if (regionInput.empty())
    {
        FatalErrorInFunction
            << "No input specified that applies to region " << regionName
            << " of group " << rgID << exit(FatalError);
    }

    return regionInput.xfer();
}

void Foam::stateIndex::assembleDictionaries
(
    const dictionary& input,
    dictionary& defaults
)
{
    Info<< nl << "Defining regions" << endl;
    Info       << "================" << endl;

    //insert region0 into regionGroups_ if it is empty
    if (regionGroups_.empty() && !solutionSchedule_.active())
    {
        regionGroups_.insert
        (
            "fluid",
            wordList(1, polyMesh::defaultRegion)
        );
    }

    // ***1***
    //Construct region states and setup state dependent entries
    if (solutionSchedule_.active())
    {
        forAll(solutionSchedule_.instances(), inI)
        {
            const solutionInstance& inst = solutionSchedule_.instances()[inI];
            forAll(inst.solutionRegions(), solRegI)
            {
                const wordList& regionNames = inst.solutionRegions()[solRegI];
                forAll(regionNames, rI)
                {
                    word regionGroupID = word::null;
                    forAllConstIter(HashTable<wordList>, regionGroups_, iter)
                    {
                        const wordList& regionGroupNames(iter());
                        forAll(regionGroupNames, rrI)
                        {
                            if (regionNames[rI]==regionGroupNames[rrI])
                            {
                                regionGroupID = iter.key();
                            }
                        }
                    }
                    dictionary regionInput
                    (
                        assembleRegionInput
                        (
                            input,
                            regionNames[rI],
                            regionGroupID
                        )
                    );
                    Info<< regionNames[rI] << " : ";
                    word meshName = inst.whichMesh(regionNames[rI]);
                    regions_.append
                    (
                        stateFunction::New
                        (
                            regionNames[rI],
                            regionInput,
                            defaults,
                            global_,
                            *this,
                            meshName
                        )
                    );
                }
            }
        }
    }
    else
    {
        forAllConstIter(HashTable<wordList>, regionGroups_, iter)
        {
            const word& regionGroupID = iter.key();
            const wordList& regionNames(iter());

            forAll(regionNames, rI)
            {
                dictionary regionInput
                (
                    assembleRegionInput
                    (
                        input,
                        regionNames[rI],
                        regionGroupID
                    )
                );
                Info<< regionNames[rI] << " : ";

                regions_.append
                (
                    stateFunction::New
                    (
                        regionNames[rI], regionInput, defaults, global_, *this
                    )
                );
            }
        }
    }

    // ***2***
    // construct regionProperties from region states
    regionProperties();


    // ***3***
    // procedural correction of state dictionaries
    global_.correct();
    forAll(regions_, ri)
    {
        regions_[ri].initialise(); //build field definitions
        regions_[ri].correct(); //perform other procedural generation
        regions_[ri].finalise(); //perform other procedural generation
        regions_[ri].updateMaster(global_); //update master from regions
    }

    // ***4***
    // Merge global direct user input for system/constant and inject functions
    //into controlDict
    global_.finalise();

}

void Foam::stateIndex::openLibs()
{
    if (global_.system().found("controlDict"))
    {
        const_cast<Time&>(runTime_).libs().open
        (
            global_.system().subDict("controlDict"),
            "libs"
        );
    }
    else
    {
        FatalError << "system/controlDict not found"
                   << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stateIndex::stateIndex
(
    const Time& runTime,
    const dictionary& input,
    dictionary& defaults,
    const bool distributed,
    const bool collated
)
:
    runTime_(runTime),
    regionGroups_
    (
        input.subDict("global")
        .lookupOrDefault<HashTable<wordList>>
        ("regionGroups", HashTable<wordList>(0))
    ),
    solutionSchedule_
    (
        runTime,
        IOdictionary
        (
            IOobject
            (
                "caseSetupDict",
                runTime.time().system(),
                runTime.db()
            ),
            input.subDict("global")
        ),
        true
    ),
    global_
    (
        runTime,
        input.subDict("global"),
        defaults,
        distributed,
        collated,
        *this
    ),
    regions_()
{
    // Assemble dictionaries
    assembleDictionaries(input, defaults);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::stateIndex::writeDictionaries()
{
    Info<< nl << "Writing dictionaries" << endl;
    Info       << "====================" << endl;

    bool writeMasterSystemDicts = true;
    bool writeMasterConstantDicts = true;

    // Write constant/system dictionaries
    forAll(regions_, ri)
    {
        //special handling of single region cases to prevent multiple overwrites
        if (regions_[ri].meshName() == polyMesh::defaultRegion)
        {
            //merge master onto region
            regions_[ri].system().merge(global_.system());
            writeMasterSystemDicts = false;
            if (regions_[ri].regionName() == polyMesh::defaultRegion)
            {
                regions_[ri].constant().merge(global_.constant());

                writeMasterConstantDicts = false;
            }
        }

        // system/f*Schemes and system/f*Solution are only written once per
        // mesh, currently. As a temporary solution we handle this by merging
        // and only writing once.
        // Check ahead for multiple contexts using the same mesh, and
        // back for those we have already merged
        wordList singleDicts
        (
            {
                "controlDict",
                "fvSolution",
                "fvSchemes",
                "faSolution",
                "faSchemes"
            }
        );
        bool exclude = false;
        forAll(regions_, rj)
        {
            if (regions_[rj].meshName() == regions_[ri].meshName())
            {
                if (rj < ri)
                {
                    exclude = true;
                    break;
                }
                else if (rj > ri)
                {
                    for (const word& subDict : regions_[rj].system().toc())
                    {
                        for (const word& singleDict : singleDicts)
                        {
                            if (subDict == singleDict)
                            {
                                regions_[ri].system().add
                                (
                                    subDict, dictionary(), false
                                );
                                regions_[ri].system().subDict(subDict).merge
                                (
                                    regions_[rj].system().subDict(subDict)
                                );
                            }
                            break;
                        }
                    }
                }
            }
        }

        regions_[ri].writeAllDictionaries(runTime_, singleDicts, exclude);
    }

    if (writeMasterSystemDicts)
    {
        global_.writeSystemDictionaries(runTime_);
    }
    if (writeMasterConstantDicts)
    {
        global_.writeConstantDictionaries(runTime_);
    }
}

void Foam::stateIndex::createRegionObjects()
{
    // load libraries so run-time selectable field classes can load
    openLibs();

    //generate meshes, store fields
    forAll(regions_, ri)
    {
        regions_[ri].createObjects(runTime_);
    }
}


void Foam::stateIndex::modifyMesh(scalar writePause)
{
    //modify mesh and update on disk
    Info<< nl << nl << "Modifying boundary mesh" << endl;
    Info<< "=======================" << endl;

    forAll(regions_, ri)
    {
        regions_[ri].modifyMesh(writePause);
    }

}


void Foam::stateIndex::createFields()
{
    //initialise fields
    Info<< nl << nl << "Creating fields" << endl;
    Info             << "===============" << endl;

    forAll(regions_, ri)
    {
        regions_[ri].createFields();
    }
}

void Foam::stateIndex::resetBoundaries(const bool &reset)
{
    global_.resetBoundarySwitch(reset);
    forAll(regions_, ri)
    {
        regions_[ri].resetBoundarySwitch(reset);
    }
}
// ************************************************************************* //
