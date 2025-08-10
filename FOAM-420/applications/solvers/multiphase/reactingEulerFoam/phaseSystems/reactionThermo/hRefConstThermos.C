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
    (c) 2015-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "makeReactionThermo.H"
#include "fluidThermo/makeThermo.H"

#include "rhoReactionThermo/rhoReactionThermo.H"
#include "rhoThermo/heRhoThermo.H"

#include "specie/specie.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "equationOfState/perfectFluid/perfectFluid.H"
#include "equationOfState/rhoConst/rhoConst.H"

#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"

#include "thermo/hRefConst/hRefConstThermo.H"

#include "transport/const/constTransport.H"

#include "mixtures/pureMixture/pureMixture.H"
#include "mixtures/multiComponentMixture/multiComponentMixture.H"

#include "include/thermoPhysicsTypes.H"

#include "mixtures/multiComponentMixtureBase/multiComponentMixtureBase.H"
#include "mixtures/multiComponentMixtureBase/multiComponent/multiComponent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Thermo type typedefs:

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectGas<specie>
        >,
        sensibleEnthalpy
    >
> constRefGasHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectFluid<specie>
        >,
        sensibleEnthalpy
    >
> constRefFluidHThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectGas<specie>
        >,
        sensibleInternalEnergy
    >
> constRefGasEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            perfectFluid<specie>
        >,
        sensibleInternalEnergy
    >
> constRefFluidEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            rhoConst<specie>
        >,
        sensibleInternalEnergy
    >
> constRefRhoConstEThermoPhysics;

typedef
constTransport
<
    species::thermo
    <
        hRefConstThermo
        <
            rhoConst<specie>
        >,
        sensibleEnthalpy
    >
> constRefRhoConstHThermoPhysics;


// pureMixture, sensibleEnthalpy:

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hRefConstThermo,
    rhoConst,
    specie
);


// pureMixture, sensibleInternalEnergy:

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    perfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hRefConstThermo,
    rhoConst,
    specie
);


// multiComponentMixture, sensibleInternalEnergy:

makeReactionMixtureThermo
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefGasEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefFluidEThermoPhysics
);

makeReactionMixtureThermo
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefRhoConstEThermoPhysics
);


// multiComponentMixture, sensibleEnthalpy:

makeReactionMixtureThermo
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefRhoConstHThermoPhysics
);

makeReactionMixtureThermo
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefFluidHThermoPhysics
);

makeReactionMixtureThermo
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    constRefGasHThermoPhysics
);

makeMultiComponentMixture(constRefGasHThermoPhysics);
makeMultiComponentMixture(constRefGasEThermoPhysics);
makeMultiComponentMixture(constRefFluidHThermoPhysics);
makeMultiComponentMixture(constRefFluidEThermoPhysics);
makeMultiComponentMixture(constRefRhoConstHThermoPhysics);
makeMultiComponentMixture(constRefRhoConstEThermoPhysics);

makeMultiComponentMixtureType(multiComponent, constRefGasHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constRefGasEThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constRefFluidHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constRefFluidEThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constRefRhoConstHThermoPhysics);
makeMultiComponentMixtureType(multiComponent, constRefRhoConstEThermoPhysics);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
