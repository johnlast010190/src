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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/timeVaryingMappedTotalPressure/timeVaryingMappedTotalPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "addStaticHead/staticHead.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
    timeVaryingMappedTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    psiName_("none"),
    gamma_(0.0),
    staticHead_(*this),
    freeStreamVelocity_(nullptr)
{}


Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
    timeVaryingMappedTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    psiName_(dict.lookupOrDefault<word>("psi", "none")),
    gamma_(psiName_ != "none" ? readScalar(dict.lookup("gamma")) : 1),
    staticHead_(*this, p, dict),
    freeStreamVelocity_
    (
        dict.found("freeStreamVelocity")
        ? new point(dict.lookup("freeStreamVelocity"))
        : nullptr
    )
{}


Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
    timeVaryingMappedTotalPressureFvPatchScalarField
(
    const timeVaryingMappedTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    psiName_(ptf.psiName_),
    gamma_(ptf.gamma_),
    staticHead_(*this, ptf.staticHead_, mapper),
    freeStreamVelocity_
    (
         ptf.freeStreamVelocity_.valid()
       ? new vector(ptf.freeStreamVelocity_())
       : nullptr
    )
{}


Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
    timeVaryingMappedTotalPressureFvPatchScalarField
(
    const timeVaryingMappedTotalPressureFvPatchScalarField& tppsf
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(tppsf),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    staticHead_(*this, tppsf.staticHead_),
    freeStreamVelocity_
    (
         tppsf.freeStreamVelocity_.valid()
       ? new vector(tppsf.freeStreamVelocity_())
       : nullptr
    )
{}


Foam::timeVaryingMappedTotalPressureFvPatchScalarField::
    timeVaryingMappedTotalPressureFvPatchScalarField
(
    const timeVaryingMappedTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    timeVaryingMappedFixedValueFvPatchScalarField(tppsf, iF),
    UName_(tppsf.UName_),
    phiName_(tppsf.phiName_),
    psiName_(tppsf.psiName_),
    gamma_(tppsf.gamma_),
    staticHead_(*this, tppsf.staticHead_),
    freeStreamVelocity_
    (
         tppsf.freeStreamVelocity_.valid()
       ? new vector(tppsf.freeStreamVelocity_())
       : nullptr
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingMappedTotalPressureFvPatchScalarField::updateCoeffs
(
    const scalarField& p0p,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }

    tmp<scalarField> pp;

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            db(),
            phiName_
        );

    if (internalField().dimensions() == dimPressure)
    {
        if (psiName_ == "none")
        {
            // Variable density and low-speed compressible flow

            const fvPatchField<scalar>& rho =
                patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), staticHead_.getRhoName());

            pp = p0p - 0.5*rho*(1.0 - pos0(phip))*magSqr(Up);
        }
        else
        {
            // High-speed compressible flow

            const fvPatchField<scalar>& psip =
                patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), psiName_);

            if (gamma_ > 1)
            {
                scalar gM1ByG = (gamma_ - 1)/gamma_;

                pp =
                    p0p
                   /pow
                    (
                        (1.0 + 0.5*psip*gM1ByG*(1.0 - pos0(phip))*magSqr(Up)),
                        1.0/gM1ByG
                    );
            }
            else
            {
                pp = p0p/(1.0 + 0.5*psip*(1.0 - pos0(phip))*magSqr(Up));
            }
        }

    }
    else if (internalField().dimensions() == dimPressure/dimDensity)
    {
        // Incompressible flow
        pp = p0p - 0.5*neg(phip)*magSqr(Up);
    }
    else
    {
        FatalErrorInFunction
            << " Incorrect pressure dimensions " << internalField().dimensions()
            << nl
            << "    Should be " << dimPressure
            << " for compressible/variable density flow" << nl
            << "    or " << dimPressure/dimDensity
            << " for incompressible flow," << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << exit(FatalError);
    }


    staticHead_.addStaticHead(pp.ref());


    forceAssign(pp);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::timeVaryingMappedTotalPressureFvPatchScalarField::updateCoeffs()
{
    scalarField p0(this->size(), 0.0);

    checkTable();

    // Interpolate between the sampled data

    scalar wantedAverage;

    if (endSampleTime_ == -1)
    {
        // Only start value
        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, non-interpolated values"
                << " from start time:"
                << sampleTimes_[startSampleTime_].name() << nl;
        }

        p0 = startSampledValues_;
        wantedAverage = startAverage_;
    }
    else
    {
        scalar start = sampleTimes_[startSampleTime_].value();
        scalar end = sampleTimes_[endSampleTime_].value();

        scalar s = (this->db().time().value() - start)/(end - start);

        if (debug)
        {
            Pout<< "updateCoeffs : Sampled, interpolated values"
                << " between start time:"
                << sampleTimes_[startSampleTime_].name()
                << " and end time:" << sampleTimes_[endSampleTime_].name()
                << " with weight:" << s << endl;
        }

        p0=((1 - s)*startSampledValues_ + s*endSampledValues_);
        wantedAverage = (1 - s)*startAverage_ + s*endAverage_;
    }

    // Enforce average. Either by scaling (if scaling factor > 0.5) or by
    // offsetting.
    if (setAverage_)
    {
        const scalarField& fld = p0;

        scalar averagePsi =
            gSum(this->patch().magSf()*fld)
           /gSum(this->patch().magSf());

        if (debug)
        {
            Pout<< "updateCoeffs :"
                << " actual average:" << averagePsi
                << " wanted average:" << wantedAverage
                << endl;
        }

        if (mag(averagePsi) < VSMALL)
        {
            // Field too small to scale. Offset instead.
            const scalar offset = wantedAverage - averagePsi;
            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " offsetting with:" << offset << endl;
            }
            p0=(fld + offset);
        }
        else
        {
            const scalar scale = mag(wantedAverage)/mag(averagePsi);

            if (debug)
            {
                Pout<< "updateCoeffs :"
                    << " scaling with:" << scale << endl;
            }
            p0=(scale*fld);
        }
    }

    // Apply offset to mapped values
    if (offset_.valid())
    {
        const scalar t = this->db().time().timeOutputValue();
        p0=(p0 + offset_->value(t));
    }

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(p0)
            << " max:" << gMax(p0)
            << " avg:" << gAverage(p0) << endl;
    }

    updateCoeffs
    (
        p0,
        patch().lookupPatchFieldInDb<volVectorField, vector>(db(), UName())
    );
}


void Foam::timeVaryingMappedTotalPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    timeVaryingMappedFixedValueFvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    os.writeEntry("psi", psiName_);
    os.writeEntry("gamma", gamma_);
    staticHead_.write(os);

    if (freeStreamVelocity_.valid())
    {
        os.writeEntry("freeStreamVelocity", freeStreamVelocity_());
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        timeVaryingMappedTotalPressureFvPatchScalarField
    );
}

// ************************************************************************* //
