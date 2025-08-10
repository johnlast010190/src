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

\*---------------------------------------------------------------------------*/

#include "Cole.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "compressibleTurbulenceModel.H"
#include "ThermalDiffusivity/ThermalDiffusivity.H"
#include "PhaseCompressibleTurbulenceModel/PhaseCompressibleTurbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace wallBoilingModels
{
namespace departureFrequencyModels
{
    defineTypeNameAndDebug(Cole, 0);
    addToRunTimeSelectionTable
    (
        departureFrequencyModel,
        Cole,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureFrequencyModels::
Cole::Cole(const dictionary& dict)
:
    departureFrequencyModel()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallBoilingModels::departureFrequencyModels::
Cole::~Cole()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::wallBoilingModels::departureFrequencyModels::
Cole::fDeparture
(
    const phaseModel& liquid,
    const phaseModel& vapor,
    const label patchi,
    const scalarField& dDep
) const
{
    // Gravitational acceleration
    const uniformDimensionedVectorField& g =
        liquid.mesh().lookupObject<uniformDimensionedVectorField>("g");

    const fvPatchScalarField& rhoLiquid =
        liquid.turbulence().rho().boundaryField()[patchi];

    const fvPatchScalarField& rhoVapor =
        vapor.turbulence().rho().boundaryField()[patchi];

    return sqrt
    (
        4*mag(g).value()
       *max(rhoLiquid - rhoVapor, scalar(0.1))
       /(3*dDep*rhoLiquid)
    );
}


// ************************************************************************* //
