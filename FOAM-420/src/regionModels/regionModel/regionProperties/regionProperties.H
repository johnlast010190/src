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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2019-2022 Esi Ltd.

Class
    Foam::regionProperties

Description
    Simple class to hold region information for coupled region simulations.

    The 'regions' list gives per group (e.g. 'fluid', 'solid') the names of
    the regions. The group names have no particular meaning at this level.

    Groups that should be treated by domain decomposition as one effective
    region are listed in the 'consolidatedGroups' list.

    Class derives from HashTable with key as the group name. Group names in
    original order are available from groupNames(), and consolidated
    groups from consolidatedGroupNames().

SourceFiles
    regionProperties.C

\*---------------------------------------------------------------------------*/

#ifndef regionProperties_H
#define regionProperties_H

#include "db/Time/Time.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "primitives/strings/lists/hashedWordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class regionProperties Declaration
\*---------------------------------------------------------------------------*/

class regionProperties
:
    public HashTable<wordList>
{

        //- Disallow default bitwise copy construct
        regionProperties(const regionProperties&);

        //- Disallow default bitwise assignment
        void operator=(const regionProperties&);

        IOdictionary dict_;

        hashedWordList groupNames_;
        hashedWordList consolidatedGroupNames_;
        hashedWordList inactiveRegionsNames_;

public:

    // Constructors

        //- Construct from components
        regionProperties(const Time& runTime);


    //- Destructor
    ~regionProperties();


    // Member Functions

        bool found()
        {
            return dict_.typeHeaderOk<IOdictionary>(true);
        }

        const dictionary& dict()
        {
            return dict_;
        }

        const hashedWordList& groupNames()
        {
            return groupNames_;
        }

        const hashedWordList& consolidatedGroupNames()
        {
            return consolidatedGroupNames_;
        }

        const hashedWordList& inactiveRegions() const
        {
            return inactiveRegionsNames_;
        }

        wordList allRegionNames() const
        {
            wordList regionNames;
            forAllConstIter(HashTable<wordList>, *this, iter)
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
            return regionNames;
        }
};


const word& regionDir(const word& regionName);

wordList selectRegionNames(const argList& args, const Time& runTime);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
