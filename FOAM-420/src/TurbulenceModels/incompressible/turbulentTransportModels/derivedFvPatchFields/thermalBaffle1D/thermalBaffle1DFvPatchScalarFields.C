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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "thermalBaffle1DFvPatchScalarFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTemplateTypeNameAndDebugWithName
(
    Foam::incompressible::constSolid_thermalBaffle1DFvPatchScalarField,
    "incompressible::thermalBaffle1D<hConstSolidThermoPhysics>",
    0
);

defineTemplateTypeNameAndDebugWithName
(
    Foam::incompressible::expoSolid_thermalBaffle1DFvPatchScalarField,
    "incompressible::thermalBaffle1D<hPowerSolidThermoPhysics>",
    0
);

namespace Foam
{
namespace incompressible
{
    addToPatchFieldRunTimeSelection
    (
        fvPatchScalarField,
        constSolid_thermalBaffle1DFvPatchScalarField
    );

    addToPatchFieldRunTimeSelection
    (
        fvPatchScalarField,
        expoSolid_thermalBaffle1DFvPatchScalarField
    );
}
}


// ************************************************************************* //
