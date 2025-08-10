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
    (c) 2012 OpenFOAM Foundation
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "regionProperties/regionProperties.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "primitives/UnbracketedTuple2/UnbracketedTuple2.H"
#include "fvMesh/fvMesh.H"
#include "global/argList/argList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionProperties::regionProperties(const Time& runTime)
:
    dict_
    (
        IOobject
        (
            "regionProperties",
            runTime.time().constant(),
            runTime.db(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    )
{
    // Initialise with defaults. One group ('default') containing only the
    // global ('region0'). Overwritten below if entries found.
    groupNames_.append("default");
    this->insert("default", wordList(1, fvMesh::defaultRegion));

    if (found())
    {
        dict_.readOpt() = IOobject::MUST_READ_IF_MODIFIED;

        if (dict_.found("regions"))
        {
            // Read like this (rather than straight into the HashTable) to
            // preserve ordering while keeping file format backward-compatible
            List<UnbracketedTuple2<word, wordList>> regionGroups
            (
                dict_.lookup("regions")
            );

            wordList groupNames(regionGroups.size());
            forAll(groupNames, i)
            {
                groupNames[i] = regionGroups[i].first();
            }
            groupNames_ = hashedWordList(groupNames);

            HashTable<wordList> htTmp(dict_.lookup("regions"));
            this->HashTable<wordList>::transfer(htTmp);
        }

        consolidatedGroupNames_ =
            hashedWordList
            (
                dict_.lookupOrDefault("consolidatedGroups", hashedWordList())
            );

        inactiveRegionsNames_ =
            hashedWordList
            (
                dict_.lookupOrDefault("inactiveRegions", hashedWordList())
            );
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionProperties::~regionProperties()
{}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::regionDir(const word& regionName)
{
    return
        regionName == polyMesh::defaultRegion
      ? word::null
      : regionName;
}


Foam::wordList Foam::selectRegionNames
(
    const argList& args,
    const Time& runTime
)
{
    const bool allRegions = args.optionFound("allRegions");

    wordList regionNames;

    if (allRegions)
    {
        const regionProperties rp(runTime);
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
    }
    else
    {
        word regionName;
        if (args.optionReadIfPresent("region", regionName))
        {
            regionNames = wordList(1, regionName);
        }
        else
        {
            regionNames = wordList(1, polyMesh::defaultRegion);
        }
    }

    return regionNames;
}


// ************************************************************************* //
