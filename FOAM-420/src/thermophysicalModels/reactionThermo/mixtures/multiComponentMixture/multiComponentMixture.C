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

#include "mixtures/multiComponentMixture/multiComponentMixture.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::constructSpeciesData
(
    const objectRegistry& obr,
    const dictionary& thermoDict
)
{
    return componentMixture_->constructSpeciesData(obr, thermoDict);
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::correctMassFractions()
{
    componentMixture_->correctMassFractions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const HashPtrTable<ThermoType>& thermoData,
    const objectRegistry& obr,
    const word& phaseName
)
:
    basicSpecieMixture(thermoDict, specieNames, obr, phaseName),
    componentMixture_
    (
        multiComponentMixtureBase<ThermoType>::New
        (
            thermoDict, specieNames, thermoData, obr,
            this->species_, this->Y_, phaseName
        )
    )
{}


template<class ThermoType>
Foam::multiComponentMixture<ThermoType>::multiComponentMixture
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const word& phaseName
)
:
    basicSpecieMixture
    (
        thermoDict,
        thermoDict.lookup("species"),
        obr,
        phaseName
    ),
    componentMixture_
    (
        multiComponentMixtureBase<ThermoType>::New
        (
            thermoDict, obr, this->species_, this->Y_, phaseName
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellMixture
(
    const label celli
) const
{
    return componentMixture_->cellMixture(celli);
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    return componentMixture_->patchFaceMixture(patchi, facei);
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    return componentMixture_->cellVolMixture(p, T, celli);
}


template<class ThermoType>
const ThermoType& Foam::multiComponentMixture<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    return componentMixture_->patchFaceVolMixture(p, T, patchi, facei);
}


template<class ThermoType>
void Foam::multiComponentMixture<ThermoType>::read
(
    const objectRegistry& obr,
    const dictionary& thermoDict
)
{
    componentMixture_->read(obr, thermoDict);
}


// ************************************************************************* //
