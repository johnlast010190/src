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
    (c) 2012-2016 OpenFOAM Foundation
    (c) 2021-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "chemistryModel/basicChemistryModel/basicChemistryModel.H"
#include "include/dummyThermo.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ChemistryModel>
Foam::autoPtr<ChemistryModel> Foam::basicChemistryModel::New
(
    const objectRegistry& obr,
    const word& phaseName
)
{
    IOdictionary chemistryDict
    (
        IOobject
        (
            IOobject::groupName("chemistryProperties", phaseName),
            obr.time().constant(),
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word chemistryTypeName;

    if (chemistryDict.isDict("chemistryType"))
    {
        const dictionary& chemistryTypeDict
        (
            chemistryDict.subDict("chemistryType")
        );

        Info<< "Selecting chemistry type " << chemistryTypeDict << endl;

        const int nCmpt = 8;
        const char* cmptNames[nCmpt] =
        {
            "chemistrySolver",
            "chemistryThermo",
            "reactionThermo",
            "transport",
            "thermo",
            "equationOfState",
            "specie",
            "energy"
        };
        IOobject thermoHeader
        (
            IOobject::groupName(basicThermo::dictName, phaseName),
            obr.time().caseConstant(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        );
        // Old style dict will be replaced by new style
        IOobject materialsHeader
        (
            basicThermo::matDictName,
            obr.time().caseSystem(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        );
        const bool isMaterials(materialsHeader.typeHeaderOk<IOdictionary>(true));

        IOdictionary thermoDict(isMaterials ? materialsHeader : thermoHeader);

        word thermoTypeName;

        dictionary thermoTypeDict;
        if (thermoDict.isDict("thermoType"))
        {
            thermoTypeDict = thermoDict.subDict("thermoType");
            thermoTypeName =
                word(thermoTypeDict.lookup("transport")) + '<'
              + word(thermoTypeDict.lookup("thermo")) + '<'
              + word(thermoTypeDict.lookup("equationOfState")) + '<'
              + word(thermoTypeDict.lookup("specie")) + ">>,"
              + word(thermoTypeDict.lookup("energy")) + ">";
        }
        else if (thermoDict.optionalSubDict(phaseName).found("materialType"))
        {
            const word thermoType
            (
                thermoDict.optionalSubDict(phaseName).lookup<word>("materialType")
            );
            const word newThermo = dummyThermo::typeName();
            if (thermoType.find(newThermo))
            {
                thermoTypeName = newThermo;
            }
            else
            {
                FatalIOErrorInFunction(thermoDict)
                    << "unsupported materialType " << newThermo
                    << exit(FatalIOError);
            }
        }
        else
        {
            FatalIOErrorInFunction(thermoDict)
                << "thermoType is in the old format and must be upgraded"
                << exit(FatalIOError);
        }

        Switch isTDAC(chemistryTypeDict.lookupOrDefault("TDAC", false));

        // Construct the name of the chemistry type from the components
        if (isTDAC)
        {
            chemistryTypeName =
                word(chemistryTypeDict.lookup("chemistrySolver")) + '<'
              + "TDACChemistryModel<"
              + word(chemistryTypeDict.lookup("chemistryThermo")) + ','
              + thermoTypeName + ">>";
        }
        else
        {
            chemistryTypeName =
                word(chemistryTypeDict.lookup("chemistrySolver")) + '<'
              + "chemistryModel<"
              + word(chemistryTypeDict.lookup("chemistryThermo")) + ','
              + thermoTypeName + ">>";
        }

        typename ChemistryModel::objectRegistryConstructorTable::iterator cstrIter =
            ChemistryModel::objectRegistryConstructorTable_().find(chemistryTypeName);

        if (cstrIter == ChemistryModel::objectRegistryConstructorTable_().end())
        {
            FatalErrorInFunction
                << "Unknown " << ChemistryModel::typeName << " type " << nl
                << "chemistryType" << chemistryTypeDict << nl
                << "thermoType" << thermoTypeDict << nl << nl
                << "Valid " << ChemistryModel ::typeName << " types are:"
                << nl << nl;

            // Get the list of all the suitable chemistry packages available
            wordList validChemistryTypeNames;
            for (const auto& p : ChemistryModel::objectRegistryConstructorTable_()) {
                validChemistryTypeNames.append(p.first);
            }

            // Build a table of the thermo packages constituent parts
            // Note: row-0 contains the names of constituent parts
            List<wordList> validChemistryTypeNameCmpts
            (
                validChemistryTypeNames.size() + 1
            );

            validChemistryTypeNameCmpts[0].setSize(nCmpt);
            forAll(validChemistryTypeNameCmpts[0], j)
            {
                validChemistryTypeNameCmpts[0][j] = cmptNames[j];
            }

            // Split the thermo package names into their constituent parts
            forAll(validChemistryTypeNames, i)
            {
                validChemistryTypeNameCmpts[i+1] = basicThermo::splitThermoName
                (
                    validChemistryTypeNames[i],
                    nCmpt
                );
            }

            // Print the table of available packages
            // in terms of their constituent parts
            printTable(validChemistryTypeNameCmpts, FatalError);

            FatalError<< exit(FatalError);
        }

        return autoPtr<ChemistryModel>(cstrIter->second(obr, phaseName));
    }
    else
    {
        chemistryTypeName =
            word(chemistryDict.lookup("chemistryType"));

        Info<< "Selecting chemistry type " << chemistryTypeName << endl;

        const auto ctor = ctorTableLookup(std::string{ChemistryModel::typeName} + " type", ChemistryModel::objectRegistryConstructorTable_(), chemistryTypeName);
        return autoPtr<ChemistryModel>(ctor(obr, phaseName));
    }
}

// ************************************************************************* //
