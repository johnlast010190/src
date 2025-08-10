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
    (c) 2016 OpenFOAM Foundation
    (c) 2021-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "chemistryModel/TDACChemistryModel/reduction/makeChemistryReductionMethods.H"

#include "include/thermoPhysicsTypes.H"

#include "chemistryModel/psiChemistryModel/psiChemistryModel.H"
#include "chemistryModel/rhoChemistryModel/rhoChemistryModel.H"

#include "include/dummyThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry solvers based on sensibleEnthalpy
    makeChemistryReductionMethods(psiChemistryModel, constGasHThermoPhysics);
    makeChemistryReductionMethods(psiChemistryModel, gasHThermoPhysics);
    makeChemistryReductionMethods
    (
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods
    (
        psiChemistryModel,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods(psiChemistryModel, constIcoFluidHThermoPhysics);
    makeChemistryReductionMethods(psiChemistryModel, icoPoly8HThermoPhysics);
    makeChemistryReductionMethods(psiChemistryModel, tabulatedGasHThermoPhysics);

    makeChemistryReductionMethods(rhoChemistryModel, constGasHThermoPhysics);
    makeChemistryReductionMethods(rhoChemistryModel, gasHThermoPhysics);
    makeChemistryReductionMethods
    (
        rhoChemistryModel,
        constIncompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods
    (
        rhoChemistryModel,
        incompressibleGasHThermoPhysics
    );
    makeChemistryReductionMethods(rhoChemistryModel, constIcoFluidHThermoPhysics);
    makeChemistryReductionMethods(rhoChemistryModel, icoPoly8HThermoPhysics);


    // Chemistry solvers based on sensibleInternalEnergy
    makeChemistryReductionMethods(psiChemistryModel, constGasEThermoPhysics);
    makeChemistryReductionMethods(psiChemistryModel, gasEThermoPhysics);
    makeChemistryReductionMethods
    (
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods
    (
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods(psiChemistryModel, constIcoFluidEThermoPhysics);
    makeChemistryReductionMethods(psiChemistryModel, icoPoly8EThermoPhysics);
    makeChemistryReductionMethods(psiChemistryModel, tabulatedGasEThermoPhysics);

    makeChemistryReductionMethods(rhoChemistryModel, constGasEThermoPhysics);
    makeChemistryReductionMethods(rhoChemistryModel, gasEThermoPhysics);
    makeChemistryReductionMethods
    (
        rhoChemistryModel,
        constIncompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods
    (
        rhoChemistryModel,
        incompressibleGasEThermoPhysics
    );
    makeChemistryReductionMethods(rhoChemistryModel, constIcoFluidEThermoPhysics);
    makeChemistryReductionMethods(rhoChemistryModel, icoPoly8EThermoPhysics);

    // New thermo models
    makeChemistryReductionMethods(psiChemistryModel, dummyThermo);
    makeChemistryReductionMethods(rhoChemistryModel, dummyThermo);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
