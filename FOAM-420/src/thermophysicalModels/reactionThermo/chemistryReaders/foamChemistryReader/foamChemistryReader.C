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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "chemistryReaders/foamChemistryReader/foamChemistryReader.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
Foam::speciesTable& Foam::foamChemistryReader<ThermoType>::setSpecies
(
    const dictionary& dict,
    speciesTable& species
)
{
    wordList s(dict.lookup("species"));
    species.transfer(s);
    return species;
}


template<class ThermoType>
Foam::dictionary Foam::foamChemistryReader<ThermoType>::readThermoDict
(
    const dictionary& thermoDict
)
{
    if (thermoDict.found("foamChemistryThermoFile"))
    {
        return
            dictionary
            (
                IFstream
                (
                    fileName(thermoDict.lookup("foamChemistryThermoFile")).expand()
                )()
            );
    }
    else
    {
        wordList species = thermoDict.lookup<wordList>("species");
        dictionary newThermoDict;
        forAll(species, speciei)
        {
            newThermoDict.add(species[speciei], thermoDict.subDict(species[speciei]));
        }
        return newThermoDict;
    }
    return dictionary();
}


template<class ThermoType>
void Foam::foamChemistryReader<ThermoType>::readSpeciesComposition()
{
    if (!chemDict_.found("elements"))
    {
        return;
    }

    wordList e(chemDict_.lookup("elements"));
    label currentElementIndex(0);

    DynamicList<word> elementNames_;
    HashTable<label> elementIndices_;

    forAll(e, ei)
    {
        if (!elementIndices_.found(e[ei]))
        {
            elementIndices_.insert(e[ei], currentElementIndex++);
            elementNames_.append(e[ei]);
        }
        else
        {
            IOWarningInFunction(chemDict_)
                << "element " << e[ei] << " already in table." << endl;
        }
    }

    // Loop through all species in thermoDict to retrieve
    // the species composition
    forAll(speciesTable_, si)
    {
        if (thermoDict_.subDict(speciesTable_[si]).isDict("elements"))
        {
            dictionary currentElements
            (
                thermoDict_.subDict(speciesTable_[si]).subDict("elements")
            );

            wordList currentElementsName(currentElements.toc());
            List<specieElement> currentComposition(currentElementsName.size());

            forAll(currentElementsName, eni)
            {
                currentComposition[eni].name() = currentElementsName[eni];

                currentComposition[eni].nAtoms() =
                    currentElements.lookupOrDefault
                    (
                        currentElementsName[eni],
                        0
                    );
            }

            // Add current specie composition to the hash table
            speciesCompositionTable::iterator specieCompositionIter
            (
                speciesComposition_.find(speciesTable_[si])
            );

            if (specieCompositionIter != speciesComposition_.end())
            {
                speciesComposition_.erase(specieCompositionIter);
            }

            speciesComposition_.insert(speciesTable_[si], currentComposition);
        }
        else
        {
            FatalIOErrorInFunction(thermoDict_)
                << "Specie " << speciesTable_[si]
                << " does not contain element description."
                << exit(FatalIOError);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::foamChemistryReader<ThermoType>::foamChemistryReader
(
    const fileName& reactionsFileName,
    speciesTable& species,
    const fileName& thermoFileName,
    const objectRegistry& obr
)
:
    chemistryReader<ThermoType>(),
    chemDict_
    (
        IFstream
        (
            fileName(reactionsFileName).expand()
        )()
    ),
    thermoDict_
    (
        IFstream
        (
            fileName(thermoFileName).expand()
        )()
    ),
    speciesTable_(setSpecies(chemDict_, species)),
    speciesThermo_(thermoDict_, typename ThermoType::iNew(obr)),
    reactions_(speciesTable_, speciesThermo_, chemDict_)
{
    readSpeciesComposition();
}


template<class ThermoType>
Foam::foamChemistryReader<ThermoType>::foamChemistryReader
(
    const dictionary& thermoDict,
    speciesTable& species,
    const objectRegistry& obr
)
:
    chemistryReader<ThermoType>(),
    chemDict_
    (
        thermoDict.found("foamChemistryFile")
      ? IFstream
         (
             fileName(thermoDict.lookup("foamChemistryFile")).expand()
         )()
      : thermoDict
    ),
    thermoDict_(readThermoDict(thermoDict)),
    speciesTable_(setSpecies(chemDict_, species)),
    speciesThermo_(thermoDict_, typename ThermoType::iNew(obr)),
    reactions_(speciesTable_, speciesThermo_, chemDict_)
{
    readSpeciesComposition();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
