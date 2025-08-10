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

#include "mixtures/multiComponentMixtureBase/multiComponent/multiComponent.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoType>
Foam::multiComponent<ThermoType>::multiComponent
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
    multiComponentMixtureBase<ThermoType>(thermoDict, specieNames, thermoData,
    obr, species, Y, phaseName)
{}


template<class ThermoType>
Foam::multiComponent<ThermoType>::multiComponent
(
    const dictionary& thermoDict,
    const objectRegistry& obr,
    const speciesTable& species,
    PtrList<volScalarField>& Y,
    const word& phaseName
)
:
    multiComponentMixtureBase<ThermoType>(thermoDict, obr, species, Y, phaseName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
const ThermoType& Foam::multiComponent<ThermoType>::cellMixture
(
    const label celli
) const
{
    this->mixture_ = this->Y_[0][celli]*this->speciesData_[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        this->mixture_ += this->Y_[n][celli]*this->speciesData_[n];
    }

    return this->mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponent<ThermoType>::patchFaceMixture
(
    const label patchi,
    const label facei
) const
{
    this->mixture_ = this->Y_[0].boundaryField()[patchi][facei]*this->speciesData_[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        this->mixture_ += this->Y_[n].boundaryField()[patchi][facei]*this->speciesData_[n];
    }

    return this->mixture_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponent<ThermoType>::cellVolMixture
(
    const scalar p,
    const scalar T,
    const label celli
) const
{
    scalar rhoInv = 0.0;
    forAll(this->speciesData_, i)
    {
        rhoInv += this->Y_[i][celli]/this->speciesData_[i].rho(p, T);
    }

    this->mixtureVol_ =
        this->Y_[0][celli]/this->speciesData_[0].rho(p, T)/rhoInv*this->speciesData_[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        this->mixtureVol_ +=
            this->Y_[n][celli]/this->speciesData_[n].rho(p, T)/rhoInv*this->speciesData_[n];
    }

    return this->mixtureVol_;
}


template<class ThermoType>
const ThermoType& Foam::multiComponent<ThermoType>::
patchFaceVolMixture
(
    const scalar p,
    const scalar T,
    const label patchi,
    const label facei
) const
{
    scalar rhoInv = 0.0;
    forAll(this->speciesData_, i)
    {
        rhoInv +=
            this->Y_[i].boundaryField()[patchi][facei]/this->speciesData_[i].rho(p, T);
    }

    this->mixtureVol_ =
        this->Y_[0].boundaryField()[patchi][facei]/this->speciesData_[0].rho(p, T)/rhoInv
      * this->speciesData_[0];

    for (label n=1; n<this->Y_.size(); n++)
    {
        this->mixtureVol_ +=
            this->Y_[n].boundaryField()[patchi][facei]/this->speciesData_[n].rho(p,T)
          / rhoInv*this->speciesData_[n];
    }

    return this->mixtureVol_;
}



// ************************************************************************* //
