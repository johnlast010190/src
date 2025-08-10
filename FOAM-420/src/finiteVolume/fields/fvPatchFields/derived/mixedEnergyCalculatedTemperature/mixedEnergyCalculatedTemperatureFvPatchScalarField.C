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
    (c) 2023 Esi Ltd.
    (c) 2019-2023 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mixedEnergyCalculatedTemperature/mixedEnergyCalculatedTemperatureFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug
    (
        mixedEnergyCalculatedTemperatureFvPatchScalarField,
        0
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixedEnergyCalculatedTemperatureFvPatchScalarField::
mixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(p, iF),
    heRefValue_(p.size()),
    heRefGrad_(p.size()),
    heValueFraction_(p.size())
{}


Foam::mixedEnergyCalculatedTemperatureFvPatchScalarField::
mixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    calculatedFvPatchScalarField(p, iF),
    heRefValue_(p.size()),
    heRefGrad_(p.size()),
    heValueFraction_(p.size())
{
    calculatedFvPatchScalarField::evaluate();
}


Foam::mixedEnergyCalculatedTemperatureFvPatchScalarField::
mixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const mixedEnergyCalculatedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    calculatedFvPatchScalarField(ptf, p, iF, mapper),
    heRefValue_(mapper(ptf.heRefValue_)),
    heRefGrad_(mapper(ptf.heRefGrad_)),
    heValueFraction_(mapper(ptf.heValueFraction_))
{}


Foam::mixedEnergyCalculatedTemperatureFvPatchScalarField::
mixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const mixedEnergyCalculatedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    calculatedFvPatchScalarField(ptf, iF),
    heRefValue_(ptf.heRefValue_),
    heRefGrad_(ptf.heRefGrad_),
    heValueFraction_(ptf.heValueFraction_)
{}



Foam::mixedEnergyCalculatedTemperatureFvPatchScalarField::
mixedEnergyCalculatedTemperatureFvPatchScalarField
(
    const mixedEnergyCalculatedTemperatureFvPatchScalarField& ptf
)
:
    calculatedFvPatchScalarField(ptf),
    heRefValue_(ptf.heRefValue_),
    heRefGrad_(ptf.heRefGrad_),
    heValueFraction_(ptf.heValueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixedEnergyCalculatedTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    calculatedFvPatchScalarField::autoMap(m);
    m(heRefValue_, heRefValue_);
    m(heRefGrad_, heRefGrad_);
    m(heValueFraction_, heValueFraction_);
}


void Foam::mixedEnergyCalculatedTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    calculatedFvPatchScalarField::rmap(ptf, addr);

    const mixedEnergyCalculatedTemperatureFvPatchScalarField& mptf =
        refCast<const mixedEnergyCalculatedTemperatureFvPatchScalarField>(ptf);

    heRefValue_.rmap(mptf.heRefValue_, addr);
    heRefGrad_.rmap(mptf.heRefGrad_, addr);
    heValueFraction_.rmap(mptf.heValueFraction_, addr);
}

// ************************************************************************* //
