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
    (c) 2019-2023 Esi Ltd.
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fixedElectricCurrentFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedElectricCurrentFvPatchScalarField::
fixedElectricCurrentFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    electricalBoundaryBase(p)
{
    gradient() = 0;
}


Foam::fixedElectricCurrentFvPatchScalarField::
fixedElectricCurrentFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    electricalBoundaryBase(p, dict),
    I_(Function1<scalar>::New("I", dict))
{
    fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

    if (dict.found("gradient"))
    {
        // Full restart
        gradient() = scalarField("gradient", dict, p.size());
    }
    else
    {
        // Start from user entered data. Assume zeroGrad initially since sigma
        // not available yet.
        gradient() = 0;
    }
}


Foam::fixedElectricCurrentFvPatchScalarField::
fixedElectricCurrentFvPatchScalarField
(
    const fixedElectricCurrentFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    electricalBoundaryBase(p),
    I_(ptf.I_, false)
{}


Foam::fixedElectricCurrentFvPatchScalarField::
fixedElectricCurrentFvPatchScalarField
(
    const fixedElectricCurrentFvPatchScalarField& ecpsf
)
:
    fixedGradientFvPatchScalarField(ecpsf),
    electricalBoundaryBase(ecpsf.patch()),
    I_(ecpsf.I_, false)
{}


Foam::fixedElectricCurrentFvPatchScalarField::
fixedElectricCurrentFvPatchScalarField
(
    const fixedElectricCurrentFvPatchScalarField& ecpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ecpsf, iF),
    electricalBoundaryBase(ecpsf.patch()),
    I_(ecpsf.I_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedElectricCurrentFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
}


void Foam::fixedElectricCurrentFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);
}


void Foam::fixedElectricCurrentFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    gradient() =
        I_().value(this->internalField().mesh().time().timeOutputValue())
       /(gSum(patch().magSf())*(nfSigma() & patch().nf()));

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedElectricCurrentFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeEntry("I", I_());
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedElectricCurrentFvPatchScalarField
    );
}

// ************************************************************************* //
