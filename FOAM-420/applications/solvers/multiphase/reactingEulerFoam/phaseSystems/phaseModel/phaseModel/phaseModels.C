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

#include "rhoThermo/rhoThermo.H"
#include "rhoReactionThermo/rhoReactionThermo.H"

#include "rhoCombustionModel/rhoCombustionModel/rhoCombustionModel.H"

#include "phaseModel.H"
#include "phaseModel/ThermoPhaseModel/ThermoPhaseModel.H"
#include "phaseModel/IsothermalPhaseModel/IsothermalPhaseModel.H"
#include "phaseModel/AnisothermalPhaseModel/AnisothermalPhaseModel.H"
#include "phaseModel/PurePhaseModel/PurePhaseModel.H"
#include "phaseModel/MultiComponentPhaseModel/MultiComponentPhaseModel.H"
#include "phaseModel/InertPhaseModel/InertPhaseModel.H"
#include "phaseModel/ReactingPhaseModel/ReactingPhaseModel.H"
#include "phaseModel/MovingPhaseModel/MovingPhaseModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    typedef
        MovingPhaseModel
        <
            AnisothermalPhaseModel
            <
                PurePhaseModel
                <
                    InertPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        purePhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        purePhaseModel,
        phaseSystem,
        purePhaseModel
    );

    typedef
        MovingPhaseModel
        <
            IsothermalPhaseModel
            <
                PurePhaseModel
                <
                    InertPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoThermo>
                    >
                >
            >
        >
        pureIsothermalPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        pureIsothermalPhaseModel,
        phaseSystem,
        pureIsothermalPhaseModel
    );

    typedef
        MovingPhaseModel
        <
            AnisothermalPhaseModel
            <
                MultiComponentPhaseModel
                <
                    InertPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>
                    >
                >
            >
        >
        multiComponentPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        multiComponentPhaseModel,
        phaseSystem,
        multiComponentPhaseModel
    );

    typedef
        MovingPhaseModel
        <
            AnisothermalPhaseModel
            <
                MultiComponentPhaseModel
                <
                    ReactingPhaseModel
                    <
                        ThermoPhaseModel<phaseModel, rhoReactionThermo>,
                        combustionModels::rhoCombustionModel
                    >
                >
            >
        >
        reactingPhaseModel;

    addNamedToRunTimeSelectionTable
    (
        phaseModel,
        reactingPhaseModel,
        phaseSystem,
        reactingPhaseModel
    );
}

// ************************************************************************* //
