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
    (c) 2006-2009 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/nutTurbulentIntensityLengthScaleInlet/nutTurbulentIntensityLengthScaleInletFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField::
nutTurbulentIntensityLengthScaleInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    intensity_(0.05),
    length_(0.01),
    Cmu_(0.09),
    phiName_("phi"),
    UName_("U")
{}

Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField::
nutTurbulentIntensityLengthScaleInletFvPatchScalarField
(
    const nutTurbulentIntensityLengthScaleInletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    intensity_(ptf.intensity_),
    length_(ptf.length_),
    Cmu_(ptf.Cmu_),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_)
{}

Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField::
nutTurbulentIntensityLengthScaleInletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict, false),
    intensity_(readScalar(dict.lookup("intensity"))),
    length_(readScalar(dict.lookup("length"))),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    UName_(dict.lookupOrDefault<word>("U", "U"))
{
    if (intensity_ < 0 || intensity_ > 1)
    {
        FatalErrorInFunction
            << "Turbulence intensity should be specified as a fraction 0-1 "
               "of the mean velocity\n"
               "    value given is " << intensity_
            << "\n    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }

    fvPatchField<scalar>::operator=
    (
        Field<scalar>("value", dict, p.size())
    );

    this->refValue() = *this;
    this->refGrad() = pTraits<scalar>::zero;
    this->valueFraction() = SMALL;
}

Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField::
nutTurbulentIntensityLengthScaleInletFvPatchScalarField
(
    const nutTurbulentIntensityLengthScaleInletFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    intensity_(ptf.intensity_),
    length_(ptf.length_),
    Cmu_(ptf.Cmu_),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_)
{}


Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField::
nutTurbulentIntensityLengthScaleInletFvPatchScalarField
(
    const nutTurbulentIntensityLengthScaleInletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    intensity_(ptf.intensity_),
    length_(ptf.length_),
    Cmu_(ptf.Cmu_),
    phiName_(ptf.phiName_),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        patch().lookupPatchFieldInDb<volVectorField, vector>(db(), UName_);

    scalarField k( max(SMALL,(1.5*sqr(intensity_)*magSqr(Up))) );

    scalar CmuPr4 = ::pow(Cmu_, 0.25);

    this->refValue() = CmuPr4 * sqrt(k) * length_;

    const Field<scalar>& phip = this->patch().lookupPatchFieldInDb<surfaceScalarField>
    (
        db(),
        phiName_,
        reinterpret_cast<const surfaceScalarField*>(0),
        reinterpret_cast<const scalar*>(0)
    );

    this->valueFraction() = 1.0 - pos0(phip);

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("intensity", intensity_);
    os.writeEntry("length", length_);
    writeEntryIfDifferent<scalar>(os, "Cmu", 0.09, Cmu_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::nutTurbulentIntensityLengthScaleInletFvPatchScalarField::operator=
(
    const fvPatchField<scalar>& ptf
)
{
    fvPatchField<scalar>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutTurbulentIntensityLengthScaleInletFvPatchScalarField
    );
}

// ************************************************************************* //
