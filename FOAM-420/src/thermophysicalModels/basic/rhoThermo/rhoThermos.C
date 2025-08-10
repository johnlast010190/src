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

#include "rhoThermo/rhoThermo.H"
#include "fluidThermo/makeThermo.H"

#include "specie/specie.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "equationOfState/incompressiblePerfectGas/incompressiblePerfectGas.H"
#include "equationOfState/Boussinesq/Boussinesq.H"
#include "equationOfState/rhoConst/rhoConst.H"
#include "equationOfState/perfectFluid/perfectFluid.H"
#include "equationOfState/PengRobinsonGas/PengRobinsonGas.H"
#include "equationOfState/adiabaticPerfectFluid/adiabaticPerfectFluid.H"

#include "thermo/hConst/hConstThermo.H"
#include "thermo/hFunctionT/hFunctionTThermo.H"
#include "thermo/janaf/janafThermo.H"
#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"
#include "thermo/sensibleInternalEnergy/sensibleInternalEnergy.H"
#include "thermo/thermo/thermo.H"

#include "transport/const/constTransport.H"
#include "transport/functionT/functionTTransport.H"
#include "transport/sutherland/sutherlandTransport.H"
#include "transport/tabulated/tabulatedTransport.H"

#include "equationOfState/icoPolynomial/icoPolynomial.H"
#include "thermo/hPolynomial/hPolynomialThermo.H"
#include "transport/polynomial/polynomialTransport.H"

#include "rhoThermo/heRhoThermo.H"
#include "mixtures/pureMixture/pureMixture.H"

#include "materialModels/materialMacros.H"
#include "materials/matHeRhoThermo/matHeRhoThermo.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
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
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
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
    hConstThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hConstThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hConstThermo,
    icoPolynomial,
    specie
);


makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
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
    hConstThermo,
    rhoConst,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
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
    hConstThermo,
    adiabaticPerfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hConstThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    perfectGas,
    specie
);


/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
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
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
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
    hConstThermo,
    rhoConst,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
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
    hConstThermo,
    adiabaticPerfectFluid,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    polynomialTransport,
    sensibleInternalEnergy,
    hPolynomialThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    icoPolynomial,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    incompressiblePerfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hConstThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    Boussinesq,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    tabulatedTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    tabulatedTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeMatThermo
(
    rhoThermo,
    speciesMassFractions,
    matHeRhoThermo,
    fluid
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
