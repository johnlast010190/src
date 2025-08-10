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
    (c) 2017 OpenCFD Ltd
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "derivedFvPatchFields/fanPressure/fanPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "addStaticHead/staticHead.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* NamedEnum
    <
        fanPressureFvPatchScalarField::fanFlowDirection,
        2
    >::names[] =
    {
        "in",
        "out"
    };

    const NamedEnum
    <
        fanPressureFvPatchScalarField::fanFlowDirection,
        2
    > fanPressureFvPatchScalarField::fanFlowDirectionNames_;

    template<>
    const char* NamedEnum
    <
        fanPressureFvPatchScalarField::pressureMode,
        2
    >::names[] =
    {
        "static",
        "total"
    };

    const NamedEnum
    <
        fanPressureFvPatchScalarField::pressureMode,
        2
    > fanPressureFvPatchScalarField::pressureModeNames_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::fanPressureFvPatchScalarField::addDynamicAndStaticHead
(
    const scalarField& p0p,
    const vectorField& Up
)
{
    if (updated())
    {
        return;
    }
    scalarField pp(p0p);

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            this->db(),
            phiName_
        );

    if (pressureMode_ == pTot)
    {
        if (internalField().dimensions() == dimPressure)
        {
            const fvPatchField<scalar>& rho =
                patch().lookupPatchFieldInDb<volScalarField, scalar>
                (
                    this->db(),
                    staticHead_.getRhoName()
                );

            pp -= 0.5*rho*magSqr(Up)*neg(phip);
        }
        else if (internalField().dimensions() == dimPressure/dimDensity)
        {
            pp -= 0.5*magSqr(Up)*neg(phip);
        }
    }

    staticHead_.addStaticHead(pp);

    // Relaxation
    const label currentTimeIndex = this->db().time().timeIndex();
    if (!updateOnesPerTimestep_ || currentTimeIndex != timeIndex_)
    {
        // Don't relax first time (old pressure not yet set)
        if (relax_ < 1.0)
        {
            if (!oldPressure_.valid())
            {
                oldPressure_.reset(new scalarField(pp));
            }
            else
            {
                oldPressure_() = internalField();
                pp = oldPressure_() + relax_*(pp - oldPressure_());
            }
        }

        timeIndex_ = currentTimeIndex;

        forceAssign(pp);
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    timeIndex_(0),
    updateOnesPerTimestep_(false),
    relax_(1.0),
    isTimeCurve_(false),
    fanCurve_(new Function1Types::Constant<scalar>("fanCurve", 0.0)),
    timeFanCurve_(nullptr),
    direction_(ffdOut),
    nonDimensional_(false),
    rpm_(0.0),
    dm_(0.0),
    pressureMode_(pSt),
    UName_(word::null),
    phiName_(word::null),
    p0_(p.size(), 0.0),
    staticHead_(*this)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    timeIndex_(ptf.timeIndex_),
    updateOnesPerTimestep_(ptf.updateOnesPerTimestep_),
    relax_(ptf.relax_),
    isTimeCurve_(ptf.isTimeCurve_),
    fanCurve_(!isTimeCurve_ ? ptf.fanCurve_().clone().ptr() : nullptr),
    timeFanCurve_(isTimeCurve_ ? ptf.timeFanCurve_().clone().ptr() : nullptr),
    direction_(ptf.direction_),
    nonDimensional_(ptf.nonDimensional_),
    rpm_(ptf.rpm_),
    dm_(ptf.dm_),
    pressureMode_(ptf.pressureMode_),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    p0_(mapper(ptf.p0_)),
    staticHead_(*this, ptf.staticHead_, mapper)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict, false),
    timeIndex_(0),
    updateOnesPerTimestep_(dict.lookupOrDefault("updateOnesPerTimestep", false)),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    direction_(fanFlowDirectionNames_.read(dict.lookup("direction"))),
    nonDimensional_(dict.lookupOrDefault<Switch>("nonDimensional", false)),
    rpm_(dict.lookupOrDefault<scalar>("rpm", 0.0)),
    dm_(dict.lookupOrDefault<scalar>("dm", 0.0)),
    pressureMode_(pSt),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    p0_("p0", dict, p.size()),
    staticHead_(*this, p, dict)
{
    if (dict.found("pressureMode"))
    {
        pressureMode_ = pressureModeNames_.read(dict.lookup("pressureMode"));
    }
    if (nonDimensional_)
    {
        dict.lookup("rpm") >> rpm_;
        dict.lookup("dm") >> dm_;
    }
    isTimeCurve_ = dict.lookupOrDefault("timeFanCurve", false);
    if (isTimeCurve_)
    {
        timeFanCurve_ = Function2<scalar>::New("fanCurve", dict);
    }
    else
    {
        if (dict.found("fanCurve"))
        {
            fanCurve_ = Function1<scalar>::New("fanCurve", dict);
        }
        else
        {
            FatalIOErrorIn
            (
                "fanPressureFvPatchScalarField::"
                "fanPressureFvPatchScalarField"
                "(const fvPatch&, const DimensionedField<vector, volMesh>&,"
                " const dictionary&)",
                dict
            )   << "Please supply 'fanCurve'"
                << exit(FatalIOError);
        }
    }


    if (dict.found("value"))
    {
        fvPatchField<scalar>::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<scalar>::operator=(p0_);
    }

    if (direction_ == ffdOut && pressureMode_ == pTot)
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " totalPressure curve is specified"
            << " and fan flowing at the outlet." << nl
            << " Switching the behaviour to static pressure."
            << endl;
    }
}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& pfopsf
)
:
    fixedValueFvPatchScalarField(pfopsf),
    timeIndex_(pfopsf.timeIndex_),
    updateOnesPerTimestep_(pfopsf.updateOnesPerTimestep_),
    relax_(pfopsf.relax_),
    isTimeCurve_(pfopsf.isTimeCurve_),
    fanCurve_(!isTimeCurve_ ? pfopsf.fanCurve_().clone().ptr() : nullptr),
    timeFanCurve_(isTimeCurve_ ? pfopsf.timeFanCurve_().clone().ptr() : nullptr),
    direction_(pfopsf.direction_),
    nonDimensional_(pfopsf.nonDimensional_),
    rpm_(pfopsf.rpm_),
    dm_(pfopsf.dm_),
    pressureMode_(pfopsf.pressureMode_),
    UName_(pfopsf.UName_),
    phiName_(pfopsf.phiName_),
    p0_(pfopsf.p0_),
    staticHead_(*this, pfopsf.staticHead_)
{}


Foam::fanPressureFvPatchScalarField::fanPressureFvPatchScalarField
(
    const fanPressureFvPatchScalarField& pfopsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(pfopsf, iF),
    timeIndex_(pfopsf.timeIndex_),
    updateOnesPerTimestep_(pfopsf.updateOnesPerTimestep_),
    relax_(pfopsf.relax_),
    isTimeCurve_(pfopsf.isTimeCurve_),
    fanCurve_(!isTimeCurve_ ? pfopsf.fanCurve_().clone().ptr() : nullptr),
    timeFanCurve_(isTimeCurve_ ? pfopsf.timeFanCurve_().clone().ptr() : nullptr),
    direction_(pfopsf.direction_),
    nonDimensional_(pfopsf.nonDimensional_),
    rpm_(pfopsf.rpm_),
    dm_(pfopsf.dm_),
    pressureMode_(pfopsf.pressureMode_),
    UName_(pfopsf.UName_),
    phiName_(pfopsf.phiName_),
    p0_(pfopsf.p0_),
    staticHead_(*this, pfopsf.staticHead_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve flux field
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    int dir = 2*direction_ - 1;

    // Average volumetric flow rate
    scalar volFlowRate = 0;

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        volFlowRate = dir*gSum(phip);
    }
    else if (phi.dimensions() == dimVelocity*dimArea*dimDensity)
    {
        const scalarField& rhop =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                this->db(),
                staticHead_.getRhoName()
            );
        volFlowRate = dir*gSum(phip/rhop);
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
                << "\n    on patch " << patch().name()
                << " of field " << internalField().name()
                << " in file " << internalField().objectPath() << nl
                << exit(FatalError);
    }

    if (nonDimensional_)
    {
        // Create an adimensional flow rate
        volFlowRate =
            120.0*volFlowRate/pow3(constant::mathematical::pi)/pow3(dm_);
    }

    // Pressure drop for this flow rate
    scalar pdFan;
    if (timeFanCurve_)
    {
        pdFan =
            timeFanCurve_->value
            (
                this->db().time().timeOutputValue(),
                max(volFlowRate, 0.0)
            );
    }
    else
    {
        pdFan = fanCurve_->value(max(volFlowRate, 0.0));
    }

    if (nonDimensional_)
    {
        // Convert the adimensional deltap from curve into deltaP
        pdFan = pdFan*pow4(constant::mathematical::pi)*rpm_*sqr(dm_)/1800;
    }
    if (debug)
    {
        Info<< "volFlowrate:" << volFlowRate << endl;
        Info<< "pdfan:" << pdFan << endl;
    }

    this->addDynamicAndStaticHead
    (
        p0_ - dir*pdFan,
        patch().lookupPatchFieldInDb<volVectorField, vector>
        (
            this->db(),
            UName_
        )
    );
}


void Foam::fanPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("timeFanCurve", isTimeCurve_);
    if (isTimeCurve_)
    {
        timeFanCurve_->writeData(os);
    }
    else
    {
        fanCurve_->writeData(os);
    }
    if (updateOnesPerTimestep_)
    {
        os.writeEntry("updateOnesPerTimestep", updateOnesPerTimestep_);
    }
    if (relax_ != 1.0)
    {
        os.writeEntry("relax", relax_);
    }
    os.writeEntry("direction", fanFlowDirectionNames_[direction_]);
    os.writeEntry("nonDimensional", nonDimensional_);
    os.writeEntry("rpm", rpm_);
    os.writeEntry("dm", dm_);
    os.writeEntry("pressureMode", pressureModeNames_[pressureMode_]);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    if (!staticHead_.active())
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", staticHead_.getRhoName());
    }
    p0_.writeEntry("p0", os);

    staticHead_.write(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fanPressureFvPatchScalarField
    );
};


// ************************************************************************* //
