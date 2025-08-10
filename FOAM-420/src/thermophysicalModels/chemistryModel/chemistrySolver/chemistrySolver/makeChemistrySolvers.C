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

#include "chemistrySolver/chemistrySolver/makeChemistrySolverTypes.H"

#include "include/thermoPhysicsTypes.H"
#include "include/dummyThermo.H"
#include "chemistryModel/psiChemistryModel/psiChemistryModel.H"
#include "chemistryModel/rhoChemistryModel/rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistrySolverTypes(psiChemistryModel, constGasHThermoPhysics);
    makeChemistrySolverTypes(psiChemistryModel, gasHThermoPhysics);
    makeChemistrySolverTypes
    (
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypes
    (
        psiChemistryModel,
        incompressibleGasHThermoPhysics)
    ;
    makeChemistrySolverTypes(psiChemistryModel, constIcoFluidHThermoPhysics);
    makeChemistrySolverTypes(psiChemistryModel, icoPoly8HThermoPhysics);
    makeChemistrySolverTypes(psiChemistryModel, tabulatedGasHThermoPhysics);

    makeChemistrySolverTypes(rhoChemistryModel, constGasHThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, gasHThermoPhysics);
    makeChemistrySolverTypes
    (
        rhoChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypes
    (
        rhoChemistryModel,
        incompressibleGasHThermoPhysics
    );
    makeChemistrySolverTypes(rhoChemistryModel, constIcoFluidHThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, icoPoly8HThermoPhysics);

    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistrySolverTypes(psiChemistryModel, constGasEThermoPhysics);
    makeChemistrySolverTypes(psiChemistryModel, gasEThermoPhysics);
    makeChemistrySolverTypes
    (
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypes
    (
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypes(psiChemistryModel, constIcoFluidEThermoPhysics);
    makeChemistrySolverTypes(psiChemistryModel, icoPoly8EThermoPhysics);
    makeChemistrySolverTypes(psiChemistryModel, tabulatedGasEThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, constGasEThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, gasEThermoPhysics);
    makeChemistrySolverTypes
    (
        rhoChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypes
    (
        rhoChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistrySolverTypes(rhoChemistryModel, constIcoFluidEThermoPhysics);
    makeChemistrySolverTypes(rhoChemistryModel, icoPoly8EThermoPhysics);

    // New thermo models
    makeChemistrySolverTypes(psiChemistryModel, dummyThermo);
    makeChemistrySolverTypes(rhoChemistryModel, dummyThermo);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
