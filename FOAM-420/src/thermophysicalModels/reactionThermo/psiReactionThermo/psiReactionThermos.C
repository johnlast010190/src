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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "makeReactionThermo.H"

#include "psiReactionThermo/psiReactionThermo.H"
#include "psiThermo/hePsiThermo.H"

#include "specie/specie.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "thermo/hConst/hConstThermo.H"
#include "thermo/janaf/janafThermo.H"
#include "thermo/hTabulated/hTabulatedThermo.H"
#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"
#include "thermo/thermo/thermo.H"
#include "transport/const/constTransport.H"
#include "transport/sutherland/sutherlandTransport.H"
#include "transport/tabulated/tabulatedTransport.H"

#include "mixtures/homogeneousMixture/homogeneousMixture.H"
#include "mixtures/inhomogeneousMixture/inhomogeneousMixture.H"
#include "mixtures/veryInhomogeneousMixture/veryInhomogeneousMixture.H"
#include "mixtures/multiComponentMixture/multiComponentMixture.H"
#include "mixtures/reactingMixture/reactingMixture.H"
#include "mixtures/singleStepReactingMixture/singleStepReactingMixture.H"

#include "include/thermoPhysicsTypes.H"
#include "include/dummyThermo.H"


#include "materials/matHePsiThermo/matHePsiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// constTransport, hConstThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);


// sutherlandTransport, hConstThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

// tabulatedTransport, hConstThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hTabulatedThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hTabulatedThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hTabulatedThermo,
    perfectGas,
    specie
);

// sutherlandTransport, janafThermo

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    homogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeReactionThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);


// Multi-component thermo for sensible enthalpy

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    constGasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    gasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    tabulatedGasHThermoPhysics
);


// Multi-component thermo for internal energy

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    constGasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    multiComponentMixture,
    tabulatedGasEThermoPhysics
);

// Multi-component reaction thermo for sensible enthalpy

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    constGasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    tabulatedGasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    gasHThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    singleStepReactingMixture,
    gasHThermoPhysics
);

// Multi-component reaction thermo for internal energy

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    constGasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    singleStepReactingMixture,
    gasEThermoPhysics
);

makeReactionMixtureThermo
(
    psiThermo,
    psiReactionThermo,
    hePsiThermo,
    reactingMixture,
    tabulatedGasEThermoPhysics
);


// New material library thermo
makeMatReactionThermo
(
    psiThermo,
    psiReactionThermo,
    matHePsiThermo,
    reactingMixture,
    dummyThermo,
    psiReactingFluid
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
