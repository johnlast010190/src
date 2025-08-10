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
    (c) 2013-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "basicSolidChemistryModel/basicSolidChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicSolidChemistryModel> Foam::basicSolidChemistryModel::
New
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

    const dictionary& chemistryTypeDict
    (
        chemistryDict.subDict("chemistryType")
    );

    Info<< "Selecting chemistry type " << chemistryTypeDict << endl;

    const int nCmpt = 13;
    const char* cmptNames[nCmpt] =
    {
        "chemistrySolver",
        "chemistryThermo",
        "baseChemistry",
        "transport",
        "thermo",
        "equationOfState",
        "specie",
        "energy",
        "transport",
        "thermo",
        "equationOfState",
        "specie",
        "energy"
    };

    IOdictionary thermoDict
    (
        IOobject
        (
            basicThermo::dictName,
            obr.time().constant(),
            obr,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    const dictionary& solidThermoTypeDict(thermoDict.subDict("thermoType"));
    word solidThermoTypeName
    (
        word(solidThermoTypeDict.lookup("transport")) + '<'
      + word(solidThermoTypeDict.lookup("thermo")) + '<'
      + word(solidThermoTypeDict.lookup("equationOfState")) + '<'
      + word(solidThermoTypeDict.lookup("specie")) + ">>,"
      + word(solidThermoTypeDict.lookup("energy")) + ">"
    );

    const dictionary& gasThermoTypeDict(thermoDict.subDict("gasThermoType"));
    word gasThermoTypeName
    (
        word(gasThermoTypeDict.lookup("transport")) + '<'
      + word(gasThermoTypeDict.lookup("thermo")) + '<'
      + word(gasThermoTypeDict.lookup("equationOfState")) + '<'
      + word(gasThermoTypeDict.lookup("specie")) + ">>,"
      + word(gasThermoTypeDict.lookup("energy")) + ">"
    );

    // Construct the name of the chemistry type from the components
    word chemistryTypeName
    (
        word(chemistryTypeDict.lookup("chemistrySolver")) + '<'
      + word(chemistryTypeDict.lookup("chemistryThermo")) + '<'
      + typeName + ','
      + solidThermoTypeName + ',' + gasThermoTypeName + ">>"
    );

    Info<< "chemistryTypeName " << chemistryTypeName << endl;

    objectRegistryConstructorTable::iterator cstrIter =
        objectRegistryConstructorTable_().find(chemistryTypeName);

    if (cstrIter == objectRegistryConstructorTable_().end())
    {
        FatalErrorInFunction
            << "Unknown " << typeName << " type " << nl
            << "chemistryType" << chemistryTypeDict << nl << nl
            << "Valid " << typeName << " types are:"
            << nl << nl;

        // Get the list of all the suitable chemistry packages available
        wordList validChemistryTypeNames;
        for (const auto& p : objectRegistryConstructorTable_()) {
            validChemistryTypeNames.append(p.first);
        }
        Info<< validChemistryTypeNames << endl;

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

    return autoPtr<basicSolidChemistryModel>(cstrIter->second(obr, phaseName));
}


// ************************************************************************* //
