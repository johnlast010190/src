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
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interfacialCompositionModels/interfaceCompositionModels/interfaceCompositionModel/interfaceCompositionModel.H"
#include "InterfaceCompositionModel.H"
#include "interfacialCompositionModels/interfaceCompositionModels/Henry/Henry.H"
#include "interfacialCompositionModels/interfaceCompositionModels/NonRandomTwoLiquid/NonRandomTwoLiquid.H"
#include "interfacialCompositionModels/interfaceCompositionModels/Raoult/Raoult.H"
#include "interfacialCompositionModels/interfaceCompositionModels/Saturated/Saturated.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "makeReactionThermo.H"

#include "include/thermoPhysicsTypes.H"

#include "equationOfState/rhoConst/rhoConst.H"
#include "equationOfState/perfectFluid/perfectFluid.H"

#include "mixtures/pureMixture/pureMixture.H"
#include "mixtures/multiComponentMixture/multiComponentMixture.H"
#include "mixtures/reactingMixture/reactingMixture.H"
#include "mixtures/SpecieMixture/SpecieMixture.H"

#include "rhoThermo/rhoThermo.H"
#include "rhoReactionThermo/rhoReactionThermo.H"
#include "rhoThermo/heRhoThermo.H"

#include "mixtures/multiComponentMixtureBase/multiComponentMixtureBase.H"
#include "mixtures/multiComponentMixtureBase/multiComponent/multiComponent.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        constTransport
        <
            species::thermo
            <
                hConstThermo
                <
                    perfectFluid<specie>
                >,
                sensibleInternalEnergy
            >
        > constFluidEThermoPhysics;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // multi-component liquid
    makeReactionThermo
    (
        rhoThermo,
        rhoReactionThermo,
        heRhoThermo,
        multiComponentMixture,
        constTransport,
        sensibleInternalEnergy,
        hConstThermo,
        perfectFluid,
        specie
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    using namespace interfaceCompositionModels;

    // multi-component gas in the presence of a pure liquid
    makeInterfaceCompositionType
    (
        Saturated,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        gasEThermoPhysics,
        heRhoThermo,
        rhoThermo,
        pureMixture,
        constFluidEThermoPhysics
    );

    // reacting gas in the presence of a pure liquid
    makeInterfaceCompositionType
    (
        Saturated,
        heRhoThermo,
        rhoReactionThermo,
        reactingMixture,
        gasEThermoPhysics,
        heRhoThermo,
        rhoThermo,
        pureMixture,
        constFluidEThermoPhysics
    );

    // multi-component gas in the presence of a multi-component liquid
    makeSpecieInterfaceCompositionType
    (
        Saturated,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constGasEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics
    );

    // multi-component liquid in the presence of a multi-component gas
    makeSpecieInterfaceCompositionType
    (
        Henry,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constFluidEThermoPhysics,
        heRhoThermo,
        rhoReactionThermo,
        multiComponentMixture,
        constGasEThermoPhysics
    );

    makeMultiComponentMixture(constFluidEThermoPhysics);
    makeMultiComponentMixtureType(multiComponent, constFluidEThermoPhysics);
}

// ************************************************************************* //
