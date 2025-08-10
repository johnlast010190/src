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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/wallVelocity/wallVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF)
{}


Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const wallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper)
{}


Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false)
{
    // Make sure inputValue also has correct value without frame
    if (dict.found("value") && !dict.found("inputValue"))
    {
        setInputValue(vectorField("value", dict, p.size()));
    }
    wallVelocityFvPatchVectorField::updateCoeffs();
}


Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const wallVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf)
{}


Foam::wallVelocityFvPatchVectorField::wallVelocityFvPatchVectorField
(
    const wallVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    vectorField Up(inputValue());
    frameFieldUpdate(Up, false);

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::CoeffField<Foam::vector>>
Foam::wallVelocityFvPatchVectorField::gradientInternalBCoeffs() const
{
    typedef CoeffField<vector> TypeCoeffField;
    typedef typename TypeCoeffField::squareTypeField squareTypeField;

    tmp<TypeCoeffField> bct(new TypeCoeffField(this->size()));
    TypeCoeffField& bc = bct.ref();

    squareTypeField& bcSq = bc.asSquare();

    const fvPatch& p = this->patch();

    bcSq = -p.deltaCoeffs()*(I - p.nf()*p.nf());

    return bct;
}


Foam::tmp<Foam::vectorField>
Foam::wallVelocityFvPatchVectorField::gradientBoundaryCoeffs() const
{
    const fvPatch& p = this->patch();

    tmp<vectorField> nHat = p.nf();
    tmp<vectorField> gbv
    (
        *this - (*this & nHat())*nHat()
      + (this->patchInternalField() &  nHat())*nHat()
    );

    //not fully implicit, but lets see
    return p.deltaCoeffs()*gbv();
}


Foam::tmp<Foam::vectorField>
Foam::wallVelocityFvPatchVectorField::gradientBoundaryBCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


Foam::tmp<Foam::vectorField>
Foam::wallVelocityFvPatchVectorField::snGrad() const
{
    tmp<vectorField> nHat = this->patch().nf();
    vectorField dfield(*this - this->patchInternalField());
    dfield -= (dfield & nHat())*nHat();
    return (dfield*this->patch().deltaCoeffs());
}


void Foam::wallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField(fvPatchVectorField, wallVelocityFvPatchVectorField);
}

// ************************************************************************* //
