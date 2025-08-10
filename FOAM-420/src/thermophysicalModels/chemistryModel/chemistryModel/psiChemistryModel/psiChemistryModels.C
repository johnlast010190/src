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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

InClass
    Foam::psiChemistryModel

Description
    Creates chemistry model instances templated on the type of thermodynamics

\*---------------------------------------------------------------------------*/

#include "chemistryModel/makeChemistryModel.H"

#include "chemistryModel/psiChemistryModel/psiChemistryModel.H"
#include "chemistryModel/chemistryModel/chemistryModel.H"
#include "chemistryModel/TDACChemistryModel/TDACChemistryModel.H"
#include "include/thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Chemistry moldels based on sensibleEnthalpy
    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        gasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constIcoFluidHThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        icoPoly8HThermoPhysics
    );


    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        gasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constIncompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        incompressibleGasHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constIcoFluidHThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        icoPoly8HThermoPhysics
    );


    // Chemistry moldels based on sensibleInternalEnergy
    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        gasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        constIcoFluidEThermoPhysics
    );

    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        icoPoly8EThermoPhysics
    );


    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        gasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constIncompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        incompressibleGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        constIcoFluidEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        icoPoly8EThermoPhysics
    );

    // tabulated
    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        tabulatedGasEThermoPhysics
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        tabulatedGasEThermoPhysics
    );
    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        tabulatedGasHThermoPhysics
    );
    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        tabulatedGasHThermoPhysics
    );

    // New thermo models
    makeChemistryModel
    (
        chemistryModel,
        psiChemistryModel,
        dummyThermo
    );

    makeChemistryModel
    (
        TDACChemistryModel,
        psiChemistryModel,
        dummyThermo
    );

}

// ************************************************************************* //
