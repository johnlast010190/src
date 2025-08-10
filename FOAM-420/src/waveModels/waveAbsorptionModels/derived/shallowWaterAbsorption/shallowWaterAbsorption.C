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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2015 IH-Cantabria

\*---------------------------------------------------------------------------*/

#include "waveAbsorptionModels/derived/shallowWaterAbsorption/shallowWaterAbsorption.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(shallowWaterAbsorption, 0);
    addToRunTimeSelectionTable(waveModel, shallowWaterAbsorption, patch);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::waveModels::shallowWaterAbsorption::setLevel
(
    const scalar t,
    const scalar tCoeff,
    scalarField& level
) const
{
    level = waterDepthRef_;
}


void Foam::waveModels::shallowWaterAbsorption::setVelocity
(
    const scalar t,
    const scalar tCoeff,
    const scalarField& level
)
{
    // Apply zero-gradient condition to z-component of velocity only
    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);
    U_ = U.boundaryField()[patch_.index()].patchInternalField();
    U_.replace(0, 0);
    U_.replace(1, 0);
}


void Foam::waveModels::shallowWaterAbsorption::setAlpha
(
    const scalarField& level
)
{
    // Set alpha as zero-gradient
    const volScalarField& alpha =
        mesh_.lookupObject<volScalarField>(alphaName_);

    alpha_ = alpha.boundaryField()[patch_.index()].patchInternalField();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::shallowWaterAbsorption::shallowWaterAbsorption
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    waveAbsorptionModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::shallowWaterAbsorption::~shallowWaterAbsorption()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::shallowWaterAbsorption::readDict
(
    const dictionary& overrideDict
)
{
    return waveAbsorptionModel::readDict(overrideDict);
}


// ************************************************************************* //
