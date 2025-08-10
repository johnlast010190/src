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
    (c) 2016 OpenFOAM Foundation
    (c) 2021-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "chemistryModel/TDACChemistryModel/reduction/chemistryReductionMethod/chemistryReductionMethod.H"
#include "include/dummyThermo.H"
#include "db/Time/Time.H"
#include "basicThermo/basicThermo.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
Foam::autoPtr<Foam::chemistryReductionMethod<CompType, ThermoType>>
Foam::chemistryReductionMethod<CompType, ThermoType>::New
(
    const IOdictionary& dict,
    TDACChemistryModel<CompType, ThermoType>& chemistry
)
{
    IOobject thermoHeader
    (
        basicThermo::dictName,
        dict.db().time().constant(),
        dict.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );
    IOobject materialsHeader
    (
        basicThermo::matDictName,
        dict.db().time().system(),
        dict.db(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );
    const bool isMaterials(materialsHeader.typeHeaderOk<IOdictionary>(true));

    IOdictionary thermoDict(isMaterials ? materialsHeader : thermoHeader);

    word thermoTypeName;

    if (thermoDict.isDict("thermoType"))
    {
        const dictionary& thermoTypeDict(thermoDict.subDict("thermoType"));
        thermoTypeName =
            word(thermoTypeDict.lookup("transport")) + '<'
          + word(thermoTypeDict.lookup("thermo")) + '<'
          + word(thermoTypeDict.lookup("equationOfState")) + '<'
          + word(thermoTypeDict.lookup("specie")) + ">>,"
          + word(thermoTypeDict.lookup("energy")) + ">";
    }
    else if (thermoDict.found("materialType"))
    {
        const word thermoType(thermoDict.lookup<word>("materialType"));
        const word newThermo = dummyThermo::typeName();
        if (thermoType.find(newThermo))
        {
            thermoTypeName = newThermo;
        }
        else
        {
            FatalIOErrorInFunction(thermoDict)
                << "unsupported materialType"
                << exit(FatalIOError);
        }
    }
    else
    {
        FatalIOErrorInFunction(thermoDict)
            << "thermoType is in the old format and must be upgraded"
            << exit(FatalIOError);
    }

    dictionary MRdict(dict.subDict("reduction"));

    word chemistryReductionMethodTypeName =
        word(MRdict.lookup("method")) + '<'
      + word(dict.subDict("chemistryType").lookup("chemistryThermo")) + ','
      + thermoTypeName + '>';

    const auto ctor = ctorTableLookup("chemistryReductionMethodType type", dictionaryConstructorTable_(), chemistryReductionMethodTypeName);
    return autoPtr<chemistryReductionMethod<CompType, ThermoType>>
    (
        ctor(dict, chemistry)
    );
}


// ************************************************************************* //
