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
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "mixtures/multiComponentMixtureBase/multiComponentMixtureBase.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixtureBase<ThermoType>::constructSpeciesData
(
    const objectRegistry& obr,
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(obr, thermoDict.subDict(species_[i]))
        );
    }

    return speciesData_[0];
}


template<class ThermoType>
void Foam::multiComponentMixtureBase<ThermoType>::correctMassFractions()
{
    // Multiplication by 1.0 changes Yt patches to "calculated"
    volScalarField Yt("Yt", 1.0*Y_[0]);

    for (label n=1; n<Y_.size(); n++)
    {
        Yt += Y_[n];
    }

    if (min(Yt).value() < ROOTVSMALL)
    {
        FatalErrorInFunction
            << "Sum of species mass fractions is zero or less";
        if (min(Yt.primitiveField()) < ROOTVSMALL)
        {
            FatalErrorInFunction
                << " in internal field";
        }
        else
        {
            forAll(Yt.boundaryField(), patchi)
            {
                if (min(Yt.boundaryField()[patchi]) < ROOTVSMALL)
                {
                    FatalErrorInFunction
                        << " in patch field "
                        << Yt.boundaryField()[patchi].patch().name();
                }
            }
        }
        FatalErrorInFunction
            << " (fields:";
        forAll(Y_, Yi)
        {
            FatalErrorInFunction << " " << Y_[Yi].name();
        }
        FatalErrorInFunction << ")." << nl << exit(FatalError);

    }

    forAll(Y_, n)
    {
        Y_[n] /= Yt;
    }
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::autoPtr<Foam::multiComponentMixtureBase<ThermoType>>
Foam::multiComponentMixtureBase<ThermoType>::New
(
    const dictionary& dict,
    const objectRegistry& obr,
    const speciesTable& species,
    PtrList<volScalarField>& Y,
    const word& phaseName
)
{
    const word modelType(dict.lookupOrDefault<word>("multiComponentMixtureBase", "multiComponent"));

    // Info<< "Selecting incompressible transport model " << modelType << endl;

    const auto ctor = ctorTableLookup("multiComponentMixtureBase type", dictionaryConstructorTable_(), modelType);
    return autoPtr<multiComponentMixtureBase<ThermoType>>
        (ctor(dict, obr, species, Y, phaseName));
}


template<class ThermoType>
Foam::autoPtr<Foam::multiComponentMixtureBase<ThermoType>>
Foam::multiComponentMixtureBase<ThermoType>::New
(
    const dictionary& dict,
    const wordList& specieNames,
    const HashPtrTable<ThermoType>& thermoData,
    const objectRegistry& obr,
    const speciesTable& species,
    PtrList<volScalarField>& Y,
    const word& phaseName
)
{
    const word modelType(dict.lookupOrDefault<word>("multiComponentMixtureBase", "multiComponent"));

    // Info<< "Selecting incompressible transport model " << modelType << endl;

    const auto ctor = ctorTableLookup("multiComponentMixtureBase type", objectRegistryConstructorTable_(), modelType);
    return autoPtr<multiComponentMixtureBase<ThermoType>>
        (ctor(dict, specieNames, thermoData, obr, species, Y, phaseName));
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentMixtureBase<ThermoType>::multiComponentMixtureBase
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const HashPtrTable<ThermoType>& thermoData,
    const objectRegistry& obr,
    const speciesTable& species,
    PtrList<volScalarField>& Y,
    const word& phaseName
)
:
    species_(species),
    Y_(Y),
    speciesData_(species_.size()),
    mixture_("mixture", *thermoData[specieNames[0]]),
    mixtureVol_("volMixture", *thermoData[specieNames[0]])
{
    forAll(species_, i)
    {
        speciesData_.set
        (
            i,
            new ThermoType(*thermoData[species_[i]])
        );
    }

    correctMassFractions();
}


template<class ThermoType>
Foam::multiComponentMixtureBase<ThermoType>::multiComponentMixtureBase
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const speciesTable& species,
    PtrList<volScalarField>& Y,
    const word& phaseName
)
:
    species_(species),
    Y_(Y),
    speciesData_(species_.size()),
    mixture_("mixture", constructSpeciesData(obr, thermoDict)),
    mixtureVol_("volMixture", speciesData_[0])
{
    correctMassFractions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
void Foam::multiComponentMixtureBase<ThermoType>::read
(
    const objectRegistry& obr,
    const dictionary& thermoDict
)
{
    forAll(species_, i)
    {
        speciesData_[i] = ThermoType(obr, thermoDict.subDict(species_[i]));
    }
}


// ************************************************************************* //
