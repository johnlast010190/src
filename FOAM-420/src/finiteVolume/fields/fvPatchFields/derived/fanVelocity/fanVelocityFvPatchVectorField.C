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
    (c) 2010-2017 Esi Ltd.
    (c) 2006 Mark Olesen

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/fanVelocity/fanVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "primitives/functions/Function1/Constant/Constant.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
const char* NamedEnum
<
    fanVelocityFvPatchVectorField::fanFlowDirection,
    2
>::names[] =
{
    "in",
    "out"
};

const NamedEnum
<
    fanVelocityFvPatchVectorField::fanFlowDirection,
    2
> fanVelocityFvPatchVectorField::fanFlowDirectionNames_;


template<>
const char* NamedEnum
<
    fanVelocityFvPatchVectorField::fanPressureInputType,
    2
>::names[] =
{
    "fixed",
    "patch"
};

const NamedEnum
<
    fanVelocityFvPatchVectorField::fanPressureInputType,
    2
> fanVelocityFvPatchVectorField::fanPressureInputTypeNames_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scalar fanVelocityFvPatchVectorField::patchMeanPressure(label pi) const
{
    const fvPatch& cPatch(patch().boundaryMesh()[pi]);

    // Total patch area
    scalar Apatch  = gSum(cPatch.magSf());

    // patch density
    scalarField rhoField(cPatch.size(), 1.0);

    const volScalarField& pressure(db().lookupObject<volScalarField>(pName_));

    // check if pressure needs to be converted to Pascal

    if
    (
        convertPressureToPascal_
     && pressure.dimensions() == dimPressure/dimDensity
    )
    {
        //first check for density field
        if (db().foundObject<volScalarField>(rhoName_))
        {
            const volScalarField& density
            (
                db().lookupObject<volScalarField>(rhoName_)
            );

            rhoField = density.boundaryField()[pi];
        }
        else if (rhoInlet_ > 0)
        {
            rhoField = rhoInlet_;
        }
    }

    scalar Pavg
        = gSum(pressure.boundaryField()[pi] * rhoField *cPatch.magSf())/Apatch;

    return Pavg;
}


scalar fanVelocityFvPatchVectorField::filterInput(scalar val, scalar dt)
{
    //pop old element if size is equal to sampleSize_

    while
    (
        storedDp_.size()
     >= sampleSize_->value(this->db().time().timeOutputValue())
    )
    {
        storedDp_.pop();
        storedDt_.pop();
    }

    //push element onto FIFO stack
    storedDp_.push(val);
    storedDt_.push(dt);


    //loop over stacks to collect mean values
    scalar dtTotal = 0;
    scalar dpTotal = 0;

    FIFOStack<scalar>::iterator dtIter = storedDt_.begin();
    scalar storeSize = storedDp_.size();
    scalar counter = 1;

    forAllConstIter(FIFOStack<scalar>, storedDp_, iter)
    {
        scalar weight = counter/(storeSize + 1);
        dpTotal += iter()*dtIter()*weight;
        dtTotal += dtIter()*weight;
        ++dtIter;
        ++counter;
    }

    return dpTotal/dtTotal;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fanVelocityFvPatchVectorField::
fanVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    flowRateOutletVelocityFvPatchVectorField(p, iF),
    fanCurve_(),
    direction_(ffdIn),
    pName_("p"),
    pressureInput_(fpiPatch),
    fanNbPatchPressure_(0.0),
    pressurePatch_(""),
    convertPressureToPascal_(true),
    sampleSize_(),
    storedDp_(),
    storedDt_(),
    curTimeIndex_(-1)
{}

fanVelocityFvPatchVectorField::
fanVelocityFvPatchVectorField
(
    const fanVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    flowRateOutletVelocityFvPatchVectorField(ptf, p, iF, mapper),
    fanCurve_(ptf.fanCurve_().clone().ptr()),
    direction_(ptf.direction_),
    pName_(ptf.pName_),
    pressureInput_(ptf.pressureInput_),
    fanNbPatchPressure_(ptf.fanNbPatchPressure_),
    pressurePatch_(ptf.pressurePatch_),
    convertPressureToPascal_(ptf.convertPressureToPascal_),
    sampleSize_(ptf.sampleSize_->clone().ptr()),
    storedDp_(ptf.storedDp_),
    storedDt_(ptf.storedDt_),
    curTimeIndex_(ptf.curTimeIndex_)
{}

fanVelocityFvPatchVectorField::
fanVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    flowRateOutletVelocityFvPatchVectorField(p, iF),
    fanCurve_(),
    direction_(ffdOut),
    pName_(dict.lookupOrDefault<word>("p", "p")),
    pressureInput_(fpiFixed),
    fanNbPatchPressure_(dict.lookupOrDefault<scalar>("p0", 0.0)),
    pressurePatch_(""),
    convertPressureToPascal_
    (
        dict.lookupOrDefault<Switch>("convertPressureToPascal", true)
    ),
    sampleSize_(),
    storedDp_(),
    storedDt_(),
    curTimeIndex_(-1)

{
    // read filter size function
    if (dict.found("filterSize"))
    {
        sampleSize_ = Function1<scalar>::New("filterSize", dict);
    }
    else
    {
        sampleSize_.reset
        (
            new Function1Types::Constant<scalar>("filterSize", 200)
        );
    }

    // set inlet/outlet status
    if (dict.found("direction"))
    {
        direction_ = fanFlowDirectionNames_.read
        (
            dict.lookup("direction")
        );
    }

    // specify whether to use fixed or mapped "otherside" pressure
    if (dict.found("pressureInputType"))
    {
        pressureInput_ = fanPressureInputTypeNames_.read
        (
            dict.lookup("pressureInputType")
        );
    }

    // volumetric fan curve
    if (dict.found("volumetricFanCurve"))
    {
        fanCurve_ = Function1<scalar>::New("volumetricFanCurve", dict);
        volumetric_ = true;
    }
    else if (dict.found("fanCurve")) // volumetric fan curve
    {
        fanCurve_ = Function1<scalar>::New("fanCurve", dict);
        volumetric_ = true;
    }
    else //"massFanCurve"
    {
        fanCurve_ = Function1<scalar>::New("massFanCurve", dict);
        volumetric_ = false;
    }

    //get patch name for other side of fan boundary pair
    if (pressureInput_ == fpiPatch)
    {
        if (direction_ == ffdIn)
        {
            pressurePatch_ = word(dict.lookup("fanOutletPatch"));
        }
        else if (direction_ == ffdOut)
        {
            pressurePatch_ = word(dict.lookup("fanInletPatch"));
        }

        if (patch().boundaryMesh().findPatchID(pressurePatch_) == -1)
        {
            //patch not found, throw exception
            FatalErrorInFunction
              << "Neighbour patch " << pressurePatch_ << " not found." << nl
              << exit(FatalError);
        }
    }

    //data used for filtering
    // historic pressure data
    if (dict.found("dpData"))
    {
        storedDp_ = FIFOStack<scalar>(dict.lookup("dpData"));
    }
    // historic time step size
    if (dict.found("dtData"))
    {
        storedDp_ = FIFOStack<scalar>(dict.lookup("dtData"));
    }



    // set base class members of flowRateOutletVelocity BC
    currentFlowRate_ = dict.lookupOrDefault<scalar>("currentFlowRate", 0.0);

    rhoName_ = dict.lookupOrDefault<word>("rho", "none");

    rhoInlet_ = dict.lookupOrDefault<scalar>("rhoInlet", -VGREAT);
    rhoInlet_ = dict.lookupOrDefault<scalar>("rhoRef", rhoInlet_);

    rhoAverage_ = dict.lookupOrDefault<Switch>("rhoAverage", true);

    maxDelta_ = dict.lookupOrDefault<scalar>("maxDelta", 0.2);

    phiName_ = dict.lookupOrDefault<word>("phi", "phi");

    if (direction_ == ffdOut)
    {
        uniform_ = dict.lookupOrDefault<Switch>("uniform", false);
    }
    else
    {
        uniform_ = dict.lookupOrDefault<Switch>("uniform", true);
    }

    fluxValueSwitch_ = dict.lookupOrDefault<Switch>("fluxValueSwitch_", true);

    fanVelocityFvPatchVectorField::initialise(dict);

}

fanVelocityFvPatchVectorField::
fanVelocityFvPatchVectorField
(
    const fanVelocityFvPatchVectorField& ptf
)
:
    flowRateOutletVelocityFvPatchVectorField(ptf),
    fanCurve_(ptf.fanCurve_().clone().ptr()),
    direction_(ptf.direction_),
    pName_(ptf.pName_),
    pressureInput_(ptf.pressureInput_),
    fanNbPatchPressure_(ptf.fanNbPatchPressure_),
    pressurePatch_(ptf.pressurePatch_),
    convertPressureToPascal_(ptf.convertPressureToPascal_),
    sampleSize_(ptf.sampleSize_->clone().ptr()),
    storedDp_(ptf.storedDp_),
    storedDt_(ptf.storedDt_),
    curTimeIndex_(ptf.curTimeIndex_)
{}

fanVelocityFvPatchVectorField::
fanVelocityFvPatchVectorField
(
    const fanVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    flowRateOutletVelocityFvPatchVectorField(ptf, iF),
    fanCurve_(ptf.fanCurve_().clone().ptr()),
    direction_(ptf.direction_),
    pName_(ptf.pName_),
    pressureInput_(ptf.pressureInput_),
    fanNbPatchPressure_(ptf.fanNbPatchPressure_),
    pressurePatch_(ptf.pressurePatch_),
    convertPressureToPascal_(ptf.convertPressureToPascal_),
    sampleSize_(ptf.sampleSize_->clone().ptr()),
    storedDp_(ptf.storedDp_),
    storedDt_(ptf.storedDt_),
    curTimeIndex_(ptf.curTimeIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void fanVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated() || curTimeIndex_ == this->db().time().timeIndex())
    {
        return;
    }
    curTimeIndex_ = this->db().time().timeIndex();

    // set flow direction multiplier
    int dir = 2*direction_ - 1;

    // neighbour patch pressure
    if (pressureInput_ == fpiPatch)
    {
        label neibI = patch().boundaryMesh().findPatchID(pressurePatch_);
        fanNbPatchPressure_ = patchMeanPressure(neibI);
    }

    // pressure change
    scalar deltaP
        = dir * (fanNbPatchPressure_ - patchMeanPressure(patch().index()));

    // filter input pressure change using historical data
    scalar filteredDeltaP
        = filterInput(deltaP, this->db().time().deltaT().value());

    // get flow rate from fan curve
    if (debug)
    {
        InfoInFunction << "Fan curve input pressure difference "
            << "(instantaneous/filtered): "
            << deltaP << " " << filteredDeltaP << endl;
    }

    scalar flowRate = fanCurve_->value(filteredDeltaP);

    // set direction based on direction_
    flowRate *= dir;

    flowRateOutletVelocityFvPatchVectorField::updateCoeffs(flowRate);
}


void fanVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);

    fanCurve_->writeData(os);
    os.writeEntry("direction", fanFlowDirectionNames_[direction_]);
    os.writeEntry("currentFlowRate", currentFlowRate_);
    writeEntryIfDifferent<word>(os, "p", "p", pName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>
    (
        os,
        "pressureInputType",
        fanPressureInputTypeNames_[fpiFixed],
        fanPressureInputTypeNames_[pressureInput_]
    );

    if (pressureInput_ == fpiPatch)
    {
        if (direction_ == ffdIn)
        {
            os.writeEntry("fanOutletPatch", pressurePatch_);
        }
        else if (direction_ == ffdOut)
        {
            os.writeEntry("fanInletPatch", pressurePatch_);
        }
    }

    writeEntryIfDifferent<scalar>(os, "p0", 0.0, fanNbPatchPressure_);

    writeEntryIfDifferent<Switch>
    (
        os, "convertPressureToPascal", true, convertPressureToPascal_
    );

    sampleSize_->writeData(os);
    os.writeEntry("dpData", storedDp_);
    os.writeEntry("dtData", storedDt_);

    if (!volumetric_ || convertPressureToPascal_)
    {
        writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
        writeEntryIfDifferent<scalar>(os, "rhoRef", -VGREAT, rhoInlet_);
        writeEntryIfDifferent<Switch>(os, "rhoAverage", true, rhoAverage_);
    }

    writeEntryIfDifferent<scalar>(os, "maxDelta", 0.2, maxDelta_);
    writeEntryIfDifferent<Switch>(os, "uniform", false, uniform_);
    writeEntryIfDifferent<Switch>(os, "fluxValueSwitch", true, fluxValueSwitch_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    fanVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
