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
    (c) 2013-2015 OpenFOAM Foundation

InClass
    Foam::basicSolidChemistryModel

Description
    Creates solid chemistry model instances templated on the type of
    solid thermodynamics

\*---------------------------------------------------------------------------*/

#include "makeSolidChemistryModel.H"

#include "pyrolysisChemistryModel/pyrolysisChemistryModel.H"
#include "basicSolidChemistryModel/basicSolidChemistryModel.H"
#include "solidChemistryModel/solidChemistryModel.H"
#include "include/solidThermoPhysicsTypes.H"
#include "include/thermoPhysicsTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// hConstSolidThermoPhysics

makeSolidChemistryModel
(
    solidChemistryModel,
    pyrolysisChemistryModel,
    basicSolidChemistryModel,
    hConstSolidThermoPhysics
);

makeSolidGasChemistryModel
(
    solidChemistryModel,
    pyrolysisChemistryModel,
    basicSolidChemistryModel,
    hConstSolidThermoPhysics,
    gasHThermoPhysics
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// hPowerSolidThermoPhysics

makeSolidChemistryModel
(
    solidChemistryModel,
    pyrolysisChemistryModel,
    basicSolidChemistryModel,
    hPowerSolidThermoPhysics
);

makeSolidGasChemistryModel
(
    solidChemistryModel,
    pyrolysisChemistryModel,
    basicSolidChemistryModel,
    hPowerSolidThermoPhysics,
    gasHThermoPhysics
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// hExpKappaConstSolidThermoPhysics

makeSolidChemistryModel
(
    solidChemistryModel,
    pyrolysisChemistryModel,
    basicSolidChemistryModel,
    hExpKappaConstSolidThermoPhysics
);

makeSolidGasChemistryModel
(
    solidChemistryModel,
    pyrolysisChemistryModel,
    basicSolidChemistryModel,
    hExpKappaConstSolidThermoPhysics,
    gasHThermoPhysics
);


// ************************************************************************* //
