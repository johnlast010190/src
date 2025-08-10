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
    (c) 2015-2022 Esi Ltd.
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "psiThermo/psiThermo.H"
#include "fluidThermo/makeThermo.H"

#include "specie/specie.H"
#include "equationOfState/perfectGas/perfectGas.H"
#include "equationOfState/PengRobinsonGas/PengRobinsonGas.H"
#include "equationOfState/adiabaticPerfectFluid/adiabaticPerfectFluid.H"
#include "thermo/hConst/hConstThermo.H"
#include "thermo/hFunctionT/hFunctionTThermo.H"
#include "thermo/eConst/eConstThermo.H"
#include "thermo/janaf/janafThermo.H"
#include "thermo/hTabulated/hTabulatedThermo.H"
#include "thermo/sensibleEnthalpy/sensibleEnthalpy.H"
#include "thermo/sensibleInternalEnergy/sensibleInternalEnergy.H"
#include "thermo/thermo/thermo.H"

#include "transport/const/constTransport.H"
#include "transport/functionT/functionTTransport.H"
#include "transport/sutherland/sutherlandTransport.H"
#include "transport/tabulated/tabulatedTransport.H"

#include "thermo/hPolynomial/hPolynomialThermo.H"
#include "transport/polynomial/polynomialTransport.H"

#include "psiThermo/hePsiThermo.H"
#include "mixtures/pureMixture/pureMixture.H"

#include "materialModels/materialMacros.H"
#include "materials/matHePsiThermo/matHePsiThermo.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    adiabaticPerfectFluid,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hConstThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    hPolynomialThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    polynomialTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    functionTTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleEnthalpy,
    janafThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hTabulatedThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hTabulatedThermo,
    adiabaticPerfectFluid,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleEnthalpy,
    hFunctionTThermo,
    adiabaticPerfectFluid,
    specie
);


/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    constTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    functionTTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    hFunctionTThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    sutherlandTransport,
    sensibleInternalEnergy,
    janafThermo,
    PengRobinsonGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleInternalEnergy,
    eConstThermo,
    perfectGas,
    specie
);

makeThermo
(
    psiThermo,
    hePsiThermo,
    pureMixture,
    tabulatedTransport,
    sensibleInternalEnergy,
    hTabulatedThermo,
    perfectGas,
    specie
);

makeMatThermo
(
    psiThermo,
    speciesMassFractions,
    matHePsiThermo,
    psiFluid
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
