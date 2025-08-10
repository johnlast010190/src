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

#include "fields/fvPatchFields/derived/windProfileDirectionVelocity/windProfileDirectionVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "global/unitConversion/unitConversion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

vector windProfileDirectionVelocityFvPatchVectorField::direction() const
{
    if (normalDirection_)
    {
        return normalised(-gAverage(patch().nf()));
    }
    else
    {
        // also calculate wind direction unit vector
        const scalar t = this->db().time().timeOutputValue();
        const scalar windDirectionDeg = windDirectionDeg_->value(t);

        vector windDirection(Zero);

        windDirection.x() = -cos(degToRad(90 - windDirectionDeg));
        windDirection.y() = -sin(degToRad(90 - windDirectionDeg));

        const tensor T
        (
            rotationTensor
            (
                vector(0, 0, 1),
                normalised(distribution_->axis())
            )
        );

        return normalised(transform(T, windDirection));
    }

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

windProfileDirectionVelocityFvPatchVectorField::
windProfileDirectionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(p, iF),
    distribution_(),
    phiName_("phi"),
    normalDirection_(false),
    windDirectionDeg_()
{
    refValue() = *this;
    refGrad() = vector::zero;
    valueFraction() = 0.0;
}


windProfileDirectionVelocityFvPatchVectorField::
windProfileDirectionVelocityFvPatchVectorField
(
    const windProfileDirectionVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchVectorField(ptf, p, iF, mapper),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_),
    normalDirection_(ptf.normalDirection_),
    windDirectionDeg_(ptf.windDirectionDeg_, false)
{}


windProfileDirectionVelocityFvPatchVectorField::
windProfileDirectionVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchVectorField(p, iF, dict, false),
    distribution_(new fieldProfile<scalar>(p, dict)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    normalDirection_(dict.lookupOrDefault<Switch>("normalDirection", false)),
    windDirectionDeg_(Function1<scalar>::New("windDirection", dict))
{
    if (dict.found("value"))
    {
        forceAssign
        (
            Field<vector>("value", dict, p.size())
        );
    }
    else
    {
        forceAssign
        (
            fvPatchField<vector>::patchInternalField()()
        );
    }

    this->refGrad() = vector::zero;
    this->refValue() = vector::zero;
    this->valueFraction() = 1.0;
}


windProfileDirectionVelocityFvPatchVectorField::
windProfileDirectionVelocityFvPatchVectorField
(
    const windProfileDirectionVelocityFvPatchVectorField& ptf
)
:
    mixedFvPatchVectorField(ptf),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_),
    normalDirection_(ptf.normalDirection_),
    windDirectionDeg_(ptf.windDirectionDeg_, false)
{}


windProfileDirectionVelocityFvPatchVectorField::
windProfileDirectionVelocityFvPatchVectorField
(
    const windProfileDirectionVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchVectorField(ptf, iF),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_),
    normalDirection_(ptf.normalDirection_),
    windDirectionDeg_(ptf.windDirectionDeg_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void windProfileDirectionVelocityFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    vector dir(direction());
    if (!normalDirection_)
    {
        makeVectorGlobal(dir);
    }
    this->refValue() = dir*distribution_->value();

    if (phiName_ != "none")
    {
        const fvsPatchField<scalar>& phip =
            patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
            (
                this->db(),
                phiName_
            );

        scalarField negDirDotSf(1 - pos0(dir & patch().Sf()));

        tmp<vectorField> nf(patch().nf());

        forAll(negDirDotSf, facei)
        {
            if (!negDirDotSf[facei])
            {
                this->refValue()[facei] -=
                    nf()[facei]*(nf()[facei] & this->refValue()[facei]);
            }
        }

        scalarField negPhi(1 - pos0(phip));
        this->valueFraction() = pos0(negDirDotSf + negPhi - 0.5);
    }
    else
    {
        this->valueFraction() = 1.0;
    }

    mixedFvPatchVectorField::updateCoeffs();
}


tmp<Field<scalar>> windProfileDirectionVelocityFvPatchVectorField::gradient
(
) const
{
    return distribution_->gradient();
}


void windProfileDirectionVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);

    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }

    writeEntryIfDifferent<Switch>
    (
        os,
        "normalDirection",
        false,
        normalDirection_
    );
    windDirectionDeg_->writeData(os);
    distribution_->write(os);

    writeEntry("value", os);
    referenceFrameFvPatch<vector>::write(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void windProfileDirectionVelocityFvPatchVectorField::operator=
(
    const fvPatchField<vector>& ptf
)
{
    forceAssign
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    windProfileDirectionVelocityFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
