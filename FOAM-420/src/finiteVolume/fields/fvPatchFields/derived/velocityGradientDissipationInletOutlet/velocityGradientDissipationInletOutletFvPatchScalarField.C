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

#include "fields/fvPatchFields/derived/velocityGradientDissipationInletOutlet/velocityGradientDissipationInletOutletFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/volFields.H"
#include "fields/fvPatchFields/derived/windProfileDirectionVelocity/windProfileDirectionVelocityFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityGradientDissipationInletOutletFvPatchScalarField::
velocityGradientDissipationInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    phiName_("phi"),
    Uname_("U"),
    Cmu_(0.09),
    Lmax_(GREAT)
{}

Foam::velocityGradientDissipationInletOutletFvPatchScalarField::
velocityGradientDissipationInletOutletFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF, dict, false),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    Uname_(dict.lookupOrDefault<word>("U", "U")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    Lmax_(dict.lookupOrDefault<scalar>("Lmax", GREAT))
{

    fvPatchField<scalar>::operator=
    (
        Field<scalar>("value", dict, p.size())
    );

    this->refValue() = *this;
    this->refGrad() = 0.0;
    this->valueFraction() = 1;
}

Foam::velocityGradientDissipationInletOutletFvPatchScalarField::
velocityGradientDissipationInletOutletFvPatchScalarField
(
    const velocityGradientDissipationInletOutletFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    Uname_(ptf.Uname_),
    Cmu_(ptf.Cmu_),
    Lmax_(ptf.Lmax_)
{}


Foam::velocityGradientDissipationInletOutletFvPatchScalarField::
velocityGradientDissipationInletOutletFvPatchScalarField
(
    const velocityGradientDissipationInletOutletFvPatchScalarField& ptf
)
:
    mixedFvPatchField<scalar>(ptf),
    phiName_(ptf.phiName_),
    Uname_(ptf.Uname_),
    Cmu_(ptf.Cmu_),
    Lmax_(ptf.Lmax_)
{}


Foam::velocityGradientDissipationInletOutletFvPatchScalarField::
velocityGradientDissipationInletOutletFvPatchScalarField
(
    const velocityGradientDissipationInletOutletFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(ptf, iF),
    phiName_(ptf.phiName_),
    Uname_(ptf.Uname_),
    Cmu_(ptf.Cmu_),
    Lmax_(ptf.Lmax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::velocityGradientDissipationInletOutletFvPatchScalarField::
updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        patch().lookupPatchFieldInDb<volVectorField, vector>
        (
            this->db(),
            Uname_
        );

    if (isA<windProfileDirectionVelocityFvPatchVectorField>(Up))
    {
        tmp<scalarField> dUdz
        (
            refCast<const windProfileDirectionVelocityFvPatchVectorField>(Up)
            .gradient()
        );

        const fvPatchField<scalar>& kp =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                this->db(),
                "k"
            );

        scalarField sqrtk(Foam::sqrt(kp));
        scalarField Lturb( ::sqrt(Cmu_)*sqrtk / max(dUdz(), SMALL) );
        Lturb = min(Lturb, Lmax_);


        if (internalField().name() == "omega")
        {
            this->refValue() = sqrtk / Lturb;
        }
        else
        {
            this->refValue() = Cmu_ * kp * sqrtk / Lturb;
        }

    }
    else
    {
        FatalErrorInFunction
            << Up.type() << " U boundary is not compatible with "
            << type() << exit(FatalError);
    }

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            this->db(),
            phiName_
         );

    this->valueFraction() = 1.0 - pos0(phip);

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::velocityGradientDissipationInletOutletFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "U", "U", Uname_);
    writeEntryIfDifferent<scalar>(os, "Cmu", 0.09, Cmu_);
    writeEntryIfDifferent<scalar>(os, "Lmax", GREAT, Lmax_);
    writeEntry("value", os);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::velocityGradientDissipationInletOutletFvPatchScalarField::operator=
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
        velocityGradientDissipationInletOutletFvPatchScalarField
    );
}

// ************************************************************************* //
