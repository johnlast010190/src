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
    (c) 2017-2022 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "mixtures/multiComponentMixtureBase/multiComponentMixtureBase.H"
#include "mixtures/multiComponentMixtureBase/multiComponent/multiComponent.H"

#include "include/thermoPhysicsTypes.H"
#include "include/dummyThermo.H"
#include "include/solidThermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
makeMultiComponentMixture(constGasHThermoPhysics);
makeMultiComponentMixture(gasHThermoPhysics);
makeMultiComponentMixture(constIncompressibleGasHThermoPhysics);
makeMultiComponentMixture(incompressibleGasHThermoPhysics);
makeMultiComponentMixture(constIcoFluidHThermoPhysics);
makeMultiComponentMixture(icoPoly8HThermoPhysics);
makeMultiComponentMixture(tabulatedGasHThermoPhysics);
makeMultiComponentMixture(constGasEThermoPhysics);
makeMultiComponentMixture(gasEThermoPhysics);
makeMultiComponentMixture(constIncompressibleGasEThermoPhysics);
makeMultiComponentMixture(incompressibleGasEThermoPhysics);
makeMultiComponentMixture(constIcoFluidEThermoPhysics);
makeMultiComponentMixture(icoPoly8EThermoPhysics);
makeMultiComponentMixture(tabulatedGasEThermoPhysics);
makeMultiComponentMixture(hConstSolidThermoPhysics);
makeMultiComponentMixture(hPowerSolidThermoPhysics);
makeMultiComponentMixture(hTransportThermoPoly8SolidThermoPhysics);
makeMultiComponentMixture(hExpKappaConstSolidThermoPhysics);

makeMultiComponentMixtureType(multiComponent, constGasHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, gasHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constIncompressibleGasHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, incompressibleGasHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constIcoFluidHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, icoPoly8HThermoPhysics);
makeMultiComponentMixtureType(multiComponent, tabulatedGasHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constGasEThermoPhysics);
makeMultiComponentMixtureType(multiComponent, gasEThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constIncompressibleGasEThermoPhysics);
makeMultiComponentMixtureType(multiComponent, incompressibleGasEThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constIcoFluidEThermoPhysics);
makeMultiComponentMixtureType(multiComponent, icoPoly8EThermoPhysics);
makeMultiComponentMixtureType(multiComponent, tabulatedGasEThermoPhysics);
makeMultiComponentMixtureType(multiComponent, hConstSolidThermoPhysics);
makeMultiComponentMixtureType(multiComponent, hPowerSolidThermoPhysics);
makeMultiComponentMixtureType(multiComponent, hTransportThermoPoly8SolidThermoPhysics);
makeMultiComponentMixtureType(multiComponent, hExpKappaConstSolidThermoPhysics);

// Materials
makeMultiComponentMixture(dummyThermo);
makeMultiComponentMixtureType(multiComponent, dummyThermo);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
