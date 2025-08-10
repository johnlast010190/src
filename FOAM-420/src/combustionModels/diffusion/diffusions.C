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

\*---------------------------------------------------------------------------*/

#include "combustionModel/makeCombustionTypes.H"

#include "include/thermoPhysicsTypes.H"
#include "psiCombustionModel/psiThermoCombustion/psiThermoCombustion.H"
#include "rhoCombustionModel/rhoThermoCombustion/rhoThermoCombustion.H"
#include "diffusion/diffusion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Combustion models based on sensibleEnthalpy
makeCombustionTypesThermo
(
    diffusion,
    psiThermoCombustion,
    gasHThermoPhysics,
    psiCombustionModel
);

makeCombustionTypesThermo
(
    diffusion,
    psiThermoCombustion,
    constGasHThermoPhysics,
    psiCombustionModel
);

makeCombustionTypesThermo
(
    diffusion,
    rhoThermoCombustion,
    gasHThermoPhysics,
    rhoCombustionModel
);

makeCombustionTypesThermo
(
    diffusion,
    rhoThermoCombustion,
    constGasHThermoPhysics,
    rhoCombustionModel
);


// Combustion models based on sensibleInternalEnergy

makeCombustionTypesThermo
(
    diffusion,
    psiThermoCombustion,
    gasEThermoPhysics,
    psiCombustionModel
);

makeCombustionTypesThermo
(
    diffusion,
    psiThermoCombustion,
    constGasEThermoPhysics,
    psiCombustionModel
);

makeCombustionTypesThermo
(
    diffusion,
    rhoThermoCombustion,
    gasEThermoPhysics,
    rhoCombustionModel
);

makeCombustionTypesThermo
(
    diffusion,
    rhoThermoCombustion,
    constGasEThermoPhysics,
    rhoCombustionModel
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
