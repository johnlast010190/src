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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "diffusion/diffusion.H"
#include "finiteVolume/fvc/fvcGrad.H"

namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
diffusion<CombThermoType, ThermoType>::diffusion
(
    const word& modelType,
    const objectRegistry& obr,
    const word& combustionProperties,
    const word& phaseName
)
:
    singleStepCombustion<CombThermoType, ThermoType>
    (
        modelType,
        obr,
        combustionProperties,
        phaseName
    ),
    C_(readScalar(this->coeffs().lookup("C"))),
    oxidantName_(this->coeffs().template lookupOrDefault<word>("oxidant", "O2"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
diffusion<CombThermoType, ThermoType>::~diffusion()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
void diffusion<CombThermoType, ThermoType>::correct()
{
    this->wFuel_.forceAssign(
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0)
    );

    if (this->active())
    {
        this->singleMixturePtr_->fresCorrect();

        const label fuelI = this->singleMixturePtr_->fuelIndex();

        const volScalarField& YFuel =
            this->thermoPtr_->composition().Y()[fuelI];

        if (this->thermoPtr_->composition().contains(oxidantName_))
        {
            const volScalarField& YO2 =
                this->thermoPtr_->composition().Y(oxidantName_);

            this->wFuel_.forceAssign(
                C_*this->turbulence().muEff()
               *mag(fvc::grad(YFuel) & fvc::grad(YO2))
               *pos0(YFuel)*pos0(YO2)
            );
        }
    }
}


template<class CombThermoType, class ThermoType>
bool diffusion<CombThermoType, ThermoType>::read()
{
    if (singleStepCombustion<CombThermoType, ThermoType>::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        this->coeffs().readIfPresent("oxidant", oxidantName_);
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
