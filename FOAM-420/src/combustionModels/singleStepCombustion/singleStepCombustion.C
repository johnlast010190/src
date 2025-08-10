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

\*---------------------------------------------------------------------------*/

#include "singleStepCombustion/singleStepCombustion.H"
#include "finiteVolume/fvm/fvmSup.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
singleStepCombustion<CombThermoType, ThermoType>::singleStepCombustion
(
    const word& modelType,
    const objectRegistry& obr,
    const word& combustionProperties,
    const word& phaseName
)
:
    CombThermoType(modelType, obr, phaseName),
    singleMixturePtr_(nullptr),
    wFuel_
    (
        IOobject
        (
            IOobject::groupName("wFuel", phaseName),
            this->mesh().time().timeName(),
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh(),
        dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
    ),
    semiImplicit_(readBool(this->coeffs_.lookup("semiImplicit")))
{
    if (isA<singleStepReactingMixture<ThermoType>>(this->thermo()))
    {
        singleMixturePtr_ =
            &dynamic_cast<singleStepReactingMixture<ThermoType>&>
            (
                this->thermo()
            );
    }
    else
    {
        FatalErrorInFunction
            << "Inconsistent thermo package for " << this->type() << " model:\n"
            << "    " << this->thermo().type() << nl << nl
            << "Please select a thermo package based on "
            << "singleStepReactingMixture" << exit(FatalError);
    }

    if (semiImplicit_)
    {
        Info<< "Combustion mode: semi-implicit" << endl;
    }
    else
    {
        Info<< "Combustion mode: explicit" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
singleStepCombustion<CombThermoType, ThermoType>::~singleStepCombustion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
tmp<fvScalarMatrix> singleStepCombustion<CombThermoType, ThermoType>::R
(
    volScalarField& Y
) const
{
    const label specieI =
        this->thermoPtr_->composition().species()[Y.member()];

    volScalarField wSpecie
    (
        wFuel_*singleMixturePtr_->specieStoichCoeffs()[specieI]
    );

    if (semiImplicit_)
    {
        const label fNorm = singleMixturePtr_->specieProd()[specieI];
        const volScalarField fres(singleMixturePtr_->fres(specieI));
        wSpecie /= max(fNorm*(Y - fres), scalar(1e-2));

        return -fNorm*wSpecie*fres + fNorm*fvm::Sp(wSpecie, Y);
    }
    else
    {
        return wSpecie + fvm::Sp(0.0*wSpecie, Y);
    }
}


template<class CombThermoType, class ThermoType>
tmp<volScalarField>
singleStepCombustion<CombThermoType, ThermoType>::Qdot() const
{
    const label fuelI = singleMixturePtr_->fuelIndex();
    volScalarField& YFuel =
        const_cast<volScalarField&>(this->thermoPtr_->composition().Y(fuelI));

    return -singleMixturePtr_->qFuel()*(R(YFuel) & YFuel);
}


template<class CombThermoType, class ThermoType>
bool singleStepCombustion<CombThermoType, ThermoType>::read()
{
    if (CombThermoType::read())
    {
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
