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
    (c) 2012-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Thermo, class Table>
typename Table::iterator Foam::basicThermo::lookupThermo
(
    const dictionary& thermoTypeDict,
    Table& tablePtr,
    const int nCmpt,
    const char* cmptNames[],
    const word& thermoTypeName,
    const word& phaseName
)
{
    // Lookup the thermo package
    typename Table::iterator cstrIter = tablePtr.find(thermoTypeName);

    // Print error message if package not found in the table
    if (cstrIter == tablePtr.end())
    {
        FatalErrorInFunction
            << "Unknown " << Thermo::typeName << " type " << nl
            << "thermoType" << thermoTypeDict << nl << nl
            << "Valid " << Thermo::typeName << " types are:"
            << nl << nl;

        // Get the list of all the suitable thermo packages available
        wordList validThermoTypeNames;
        for (const auto& p : tablePtr) {
            validThermoTypeNames.append(p.first);
        }

        // Build a table of the thermo packages constituent parts
        // Note: row-0 contains the names of constituent parts
        List<wordList> validThermoTypeNameCmpts
        (
            validThermoTypeNames.size() + 1
        );

        validThermoTypeNameCmpts[0].setSize(nCmpt);
        forAll(validThermoTypeNameCmpts[0], j)
        {
            validThermoTypeNameCmpts[0][j] = cmptNames[j];
        }

        // Split the thermo package names into their constituent parts
        // Removing incompatible entries from the list
        label j = 0;
        forAll(validThermoTypeNames, i)
        {
            wordList names
            (
                Thermo::splitThermoName(validThermoTypeNames[i], nCmpt)
            );

            if (names.size())
            {
                validThermoTypeNameCmpts[j++] = names;
            }
        }
        validThermoTypeNameCmpts.setSize(j);

        // Print the table of available packages
        // in terms of their constituent parts
        printTable(validThermoTypeNameCmpts, FatalError);

        FatalError<< exit(FatalError);
    }

    return cstrIter;
}


template<class Thermo, class Table>
typename Table::iterator Foam::basicThermo::lookupThermo
(
    const dictionary& thermoDict,
    Table& tablePtr,
    const word& phaseName
)
{
    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));

        Info<< "Selecting thermodynamics package " << thermoTypeDict << endl;

        if (thermoTypeDict.found("properties"))
        {
            const int nCmpt = 4;
            const char* cmptNames[nCmpt] =
            {
                "type",
                "mixture",
                "properties",
                "energy"
            };

            // Construct the name of the thermo package from the components
            const word thermoTypeName
            (
                word(thermoTypeDict.lookup("type")) + '<'
              + word(thermoTypeDict.lookup("mixture")) + '<'
              + word(thermoTypeDict.lookup("properties")) + ','
              + word(thermoTypeDict.lookup("energy")) + ">>"
            );

            return lookupThermo<Thermo, Table>
            (
                thermoTypeDict,
                tablePtr,
                nCmpt,
                cmptNames,
                thermoTypeName
            );
        }
        else
        {
            const int nCmpt = 7;
            const char* cmptNames[nCmpt] =
            {
                "type",
                "mixture",
                "transport",
                "thermo",
                "equationOfState",
                "specie",
                "energy"
            };

            // Construct the name of the thermo package from the components
            const word thermoTypeName
            (
                word(thermoTypeDict.lookup("type")) + '<'
              + word(thermoTypeDict.lookup("mixture")) + '<'
              + word(thermoTypeDict.lookup("transport")) + '<'
              + word(thermoTypeDict.lookup("thermo")) + '<'
              + word(thermoTypeDict.lookup("equationOfState")) + '<'
              + word(thermoTypeDict.lookup("specie")) + ">>,"
              + word(thermoTypeDict.lookup("energy")) + ">>>"
            );

            return lookupThermo<Thermo, Table>
            (
                thermoTypeDict,
                tablePtr,
                nCmpt,
                cmptNames,
                thermoTypeName
            );
        }
    }
    else
    {
        if (thermoDict.found("thermoType"))
        {
            FatalErrorInFunction
                << "thermoType specification isn't supported option." << nl
                << exit(FatalError);
        }

        word materialType_
        (
            thermoDict.optionalSubDict(phaseName).lookup("materialType")
        );
        word extraPhaseMessage;
        word phaseDashes;
        if (phaseName != word::null)
        {
            extraPhaseMessage = " \"" + phaseName + "\"";
            phaseDashes = std::string(phaseName.size() + 3, '-');
        }
        Info<< nl
            << "Material settings" << extraPhaseMessage << nl
            << "-----------------" << phaseDashes << nl << nl
            << "Selecting materials package \"" << materialType_
            << "\"." << endl;

        //TODO: In future even without multiple species this could be reacting
        // mixture. That way no "hack" is required to load up mat. props
        if
        (
            materialType_ == "fluid"
         && thermoDict.optionalSubDict(phaseName).found("species")
        )
        {
            materialType_ = "reactingFluid";
        }
        else if
        (
            materialType_ == "psiFluid"
         && thermoDict.optionalSubDict(phaseName).found("species")
        )
        {
            materialType_ = "psiReactingFluid";
        }

        typename Table::iterator cstrIter = tablePtr.find(materialType_);

        if (cstrIter == tablePtr.end())
        {
            FatalErrorInFunction
                << "Unknown " << Thermo::typeName << " materialType "
                << materialType_ << nl << nl
                << exit(FatalError);
        }

        return cstrIter;
    }
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    IOdictionary thermoDict(loadDictionary(obr, phaseName, false));

    typename Thermo::objectRegistryConstructorTable::iterator cstrIter =
        lookupThermo<Thermo, typename Thermo::objectRegistryConstructorTable>
        (
            thermoDict,
            Thermo::objectRegistryConstructorTable_(),
            phaseName
        );

    return autoPtr<Thermo>(cstrIter->second(obr, phaseName));
}


template<class Thermo>
Foam::autoPtr<Thermo> Foam::basicThermo::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName
)
{
    typename Thermo::dictionaryConstructorTable::iterator cstrIter =
        lookupThermo<Thermo, typename Thermo::dictionaryConstructorTable>
        (
            dict,
            Thermo::dictionaryConstructorTable_(),
            phaseName
        );

    return autoPtr<Thermo>(cstrIter->second(obr, dict, phaseName));
}


// ************************************************************************* //
