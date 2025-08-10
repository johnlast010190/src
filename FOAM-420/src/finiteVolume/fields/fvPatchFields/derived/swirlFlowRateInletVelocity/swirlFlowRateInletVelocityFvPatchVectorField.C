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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/swirlFlowRateInletVelocity/swirlFlowRateInletVelocityFvPatchVectorField.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    flowRateBase(patch(), db(), false),
    rpm_(),
    outlet_(true)
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    flowRateBase(patch(), db(), ptf, false),
    rpm_(ptf.rpm_, false),
    outlet_(ptf.outlet_)
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    flowRateBase(patch(), db(), dict, false),
    rpm_(Function1<scalar>::New("rpm", dict)),
    outlet_(dict.lookupOrDefault<Switch>("implicitOutlet", true))
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    flowRateBase(patch(), db(), ptf, false),
    rpm_(ptf.rpm_, false),
    outlet_(ptf.outlet_)
{}


Foam::swirlFlowRateInletVelocityFvPatchVectorField::
swirlFlowRateInletVelocityFvPatchVectorField
(
    const swirlFlowRateInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    flowRateBase(patch(), db(), ptf, false),
    rpm_(ptf.rpm_, false),
    outlet_(ptf.outlet_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::swirlFlowRateInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    const scalar flowRate = currentFlowRate();
    const scalar rpm = rpm_->value(t);

    const scalar totArea   = gSum(patch().magSf());
    const scalar avgU = -flowRate/totArea;

    const vector avgCenter = gSum(patch().Cf()*patch().magSf())/totArea;
    const vector avgNormal = gSum(patch().Sf())/totArea;

    // Update angular velocity - convert [rpm] to [rad/s]
    tmp<vectorField> tangentialVelocity
        (
            (rpm*constant::mathematical::pi/30.0)
          * (patch().Cf() - avgCenter) ^ avgNormal
        );

    tmp<vectorField> n = patch().nf();

    forceAssign
    (
        tangentialVelocity
            + n*avgU/multiplyByRhoOrOne(scalarField(this->size(),1), false)
    );

    /*const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        // volumetric flow-rate
        forceAssign(tangentialVelocity + n*avgU);
    }
    else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                this->db(),
                rhoName_
            );

        // mass flow-rate
        forceAssign(tangentialVelocity + n*avgU/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl << exit(FatalError);
    }*/

    fixedValueFvPatchField<vector>::updateCoeffs();
}


Foam::tmp<Foam::vectorField>
Foam::swirlFlowRateInletVelocityFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{

    if (outlet_)
    {
        tmp<vectorField> vic
        (
            new vectorField(this->size(), pTraits<vector>::one)
        );

        vic.ref() *= pos0((*this) & patch().nf()());

        return vic;
    }
    else
    {
        return fixedValueFvPatchField<vector>::valueInternalCoeffs(w);
    }

}


Foam::tmp<Foam::vectorField>
Foam::swirlFlowRateInletVelocityFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{

    if (outlet_)
    {
        tmp<vectorField> vbc
        (
            new vectorField(*this)
        );

        vbc.ref() *= neg((*this) & patch().nf()());

        return vbc;
    }
    else
    {
        return fixedValueFvPatchField<vector>::valueBoundaryCoeffs(w);
    }
}


void Foam::swirlFlowRateInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    flowRateBase::write(os);
    rpm_->writeData(os);
    writeEntryIfDifferent<Switch>(os, "implicitOutlet", true, outlet_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       swirlFlowRateInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
