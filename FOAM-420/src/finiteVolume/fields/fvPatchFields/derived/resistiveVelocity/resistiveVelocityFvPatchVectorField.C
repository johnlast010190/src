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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/resistiveVelocity/resistiveVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/fvPatchFields/derived/resistivePressure/resistivePressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

resistiveVelocityFvPatchVectorField::
resistiveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(p, iF),
    pName_("p"),
    phiName_("phi")
{
    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = symmTensor::zero;
}


resistiveVelocityFvPatchVectorField::
resistiveVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    directionMixedFvPatchVectorField(p, iF),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));

    refValue() = vector::zero;
    refGrad() = vector::zero;
    valueFraction() = symmTensor::zero;
}


resistiveVelocityFvPatchVectorField::
resistiveVelocityFvPatchVectorField
(
    const resistiveVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    directionMixedFvPatchVectorField(ptf, p, iF, mapper),
    pName_(ptf.pName_),
    phiName_(ptf.phiName_)
{}


resistiveVelocityFvPatchVectorField::
resistiveVelocityFvPatchVectorField
(
    const resistiveVelocityFvPatchVectorField& ptf
)
:
    directionMixedFvPatchVectorField(ptf),
    pName_(ptf.pName_),
    phiName_(ptf.phiName_)
{}


resistiveVelocityFvPatchVectorField::
resistiveVelocityFvPatchVectorField
(
    const resistiveVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    directionMixedFvPatchVectorField(ptf, iF),
    pName_(ptf.pName_),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void resistiveVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>(db(), phiName_);


    const fvPatchField<scalar>& pp =
        patch().lookupPatchFieldInDb<volScalarField, scalar>(this->db(), pName_);

    if (isA<resistivePressureFvPatchScalarField>(pp))
    {
        //patch normal gradient
        const resistivePressureFvPatchScalarField& ppres
            = dynamic_cast<const resistivePressureFvPatchScalarField&>(pp);

        valueFraction()
            = neg(phip)*(I - sqr(patch().nf()))
            + ppres.Ct()*pos0(phip)*(I - sqr(patch().nf()));
    }
    else
    {
        FatalErrorInFunction
            << "Can only be used in conjunction with a "
            << "resistivePressureFvPatchScalarField for pressure."
            << exit(FatalError);
    }

    directionMixedFvPatchVectorField::updateCoeffs();
    directionMixedFvPatchVectorField::evaluate();
}

tmp<Field<vector>> resistiveVelocityFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>(db(), phiName_);

    return
        (
            neg(phip) *
            (*this
          - cmptMultiply
            (
                valueInternalCoeffs(this->patch().weights()),
                this->patchInternalField()
            ))
        );
}


void resistiveVelocityFvPatchVectorField::write(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    directionMixedFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    resistiveVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
