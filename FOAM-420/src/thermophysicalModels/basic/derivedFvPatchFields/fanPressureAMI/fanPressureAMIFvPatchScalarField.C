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

#include "derivedFvPatchFields/fanPressureAMI/fanPressureAMIFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "addStaticHead/staticHead.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char*
    NamedEnum
    <
        fanPressureAMIFvPatchScalarField::directionMode,
        2
    >::names[] =
    {
        "plane",
        "cylindrical"
    };
}

const Foam::NamedEnum
<
    Foam::fanPressureAMIFvPatchScalarField::directionMode,
    2
> Foam::fanPressureAMIFvPatchScalarField::directionModeNames_;


// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

void Foam::fanPressureAMIFvPatchScalarField::initializeFlowDirection()
{
    reverseFlow_ = 1;
    switch (directionMode_)
    {
        case plane:
        {
            vector patchSf = gSum(this->patch().Sf());
            scalar flowDirDotPatchSf = flowDir_&patchSf;

            if (flowDirDotPatchSf<0)
            {
                reverseFlow_=-1;
            }

            break;
        }
        case cylindrical:
        {
            if (radialFlowDir_=="in")
            {
                reverseFlow_ = -1;
            }

            vectorField cfFromOrigin(this->patch().Cf() - origin_);

            scalar sfDotlocalCoor = gSum(this->patch().Sf()&cfFromOrigin);

            if (sfDotlocalCoor<0)
            {
                reverseFlow_*=-1;
            }

            break;
        }
    }
    //Info<< this->patch().name() << tab  << volFlowRate << tab << flowDir_<<  endl;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fanPressureAMIFvPatchScalarField::fanPressureAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF),
    referenceFrameFvPatch<vector>(p),
    directionMode_(plane),
    phiName_("phi"),
    fanCurve_(new Function1Types::Constant<scalar>("fanCurve", 0.0)),
    flowDir_(Zero),
    origin_(Zero),
    radialFlowDir_(word::null),
    reverseFlow_(Zero),
    nonDimensional_(false),
    rpm_(0.0),
    dm_(0.0),
    staticHead_(*this),
    relax_(1.0),
    minJump_(-GREAT),
    maxJump_(GREAT)
{}


Foam::fanPressureAMIFvPatchScalarField::fanPressureAMIFvPatchScalarField
(
    const fanPressureAMIFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, p, iF, mapper),
    referenceFrameFvPatch<vector>
    (
        p,
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        mapper
    ),
    directionMode_(ptf.directionMode_),
    phiName_(ptf.phiName_),
    fanCurve_(ptf.fanCurve_().clone().ptr()),
    flowDir_(ptf.flowDir_),
    origin_(ptf.origin_),
    radialFlowDir_(ptf.radialFlowDir_),
    reverseFlow_(ptf.reverseFlow_),
    nonDimensional_(ptf.nonDimensional_),
    rpm_(ptf.rpm_),
    dm_(ptf.dm_),
    staticHead_(*this, ptf, mapper),
    relax_(ptf.relax_),
    minJump_(ptf.minJump_),
    maxJump_(ptf.maxJump_)
{}


Foam::fanPressureAMIFvPatchScalarField::fanPressureAMIFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF),
    referenceFrameFvPatch<vector>(dict, p),
    directionMode_
    (
        directionModeNames_.readOrDefault
        (
            dict.lookupOrDefault<word>
            (
                "directionMode",
                "plane"
            ),
            plane
        )
    ),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    flowDir_(Zero),
    origin_(Zero),
    radialFlowDir_(word::null),
    reverseFlow_(Zero),
    nonDimensional_(dict.lookupOrDefault<Switch>("nonDimensional", false)),
    rpm_(dict.lookupOrDefault<scalar>("rpm", 0.0)),
    dm_(dict.lookupOrDefault<scalar>("dm", 0.0)),
    staticHead_(*this, p, dict),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    minJump_(dict.lookupOrDefault<scalar>("minJump", -GREAT)),
    maxJump_(dict.lookupOrDefault<scalar>("maxJump", GREAT))
{
    if (nonDimensional_)
    {
        dict.lookup("rpm") >> rpm_;
        dict.lookup("dm") >> dm_;
    }

    if (dict.found("fanCurve"))
    {
        fanCurve_ = Function1<scalar>::New("fanCurve", dict);
    }
    else
    {
        FatalIOErrorIn
        (
            "fanPressureAMIFvPatchScalarField::"
            "fanPressureAMIFvPatchScalarField"
            "(const fvPatch&, const DimensionedField<vector, volMesh>&,"
            " const dictionary&)",
            dict
        )   << "Please supply 'fanCurve'"
            << exit(FatalIOError);
    }

    //- read information for flow direction
    switch (directionMode_)
    {
        case plane:
        {
            flowDir_ = dict.lookup("flowDirection");

            break;
        }
        case cylindrical:
        {
            if (dict.found("origin"))
            {
                origin_ = dict.lookup("origin");
            }
            else
            {
                origin_ = gAverage(this->patch().patch().points());
            }

            dict.lookup("radialFlowDirection") >> radialFlowDir_;
            if (radialFlowDir_!="in" && radialFlowDir_!="out")
            {
                FatalErrorInFunction
                    << "radialFlowDirection is defined as: "
                    << radialFlowDir_
                    << ". Only in or out modes are allowed."
                        << exit(FatalError);
            }

            break;
        }
    }
}


Foam::fanPressureAMIFvPatchScalarField::fanPressureAMIFvPatchScalarField
(
    const fanPressureAMIFvPatchScalarField& pfopsf
)
:
    fixedJumpAMIFvPatchField<scalar>(pfopsf),
    referenceFrameFvPatch<vector>
    (
        pfopsf.patch(),
        pfopsf.coorFramePtr_,
        pfopsf.inputValue(),
        pfopsf.inletFlux_,
        pfopsf.frameName_
    ),
    directionMode_(pfopsf.directionMode_),
    phiName_(pfopsf.phiName_),
    fanCurve_(pfopsf.fanCurve_().clone().ptr()),
    flowDir_(pfopsf.flowDir_),
    origin_(pfopsf.origin_),
    radialFlowDir_(pfopsf.radialFlowDir_),
    reverseFlow_(pfopsf.reverseFlow_),
    nonDimensional_(pfopsf.nonDimensional_),
    rpm_(pfopsf.rpm_),
    dm_(pfopsf.dm_),
    staticHead_(*this, pfopsf),
    relax_(pfopsf.relax_),
    minJump_(pfopsf.minJump_),
    maxJump_(pfopsf.maxJump_)
{}


Foam::fanPressureAMIFvPatchScalarField::fanPressureAMIFvPatchScalarField
(
    const fanPressureAMIFvPatchScalarField& pfopsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(pfopsf, iF),
    referenceFrameFvPatch<vector>
    (
        pfopsf.patch(),
        pfopsf.coorFramePtr_,
        pfopsf.inputValue(),
        pfopsf.inletFlux_,
        pfopsf.frameName_
    ),
    directionMode_(pfopsf.directionMode_),
    phiName_(pfopsf.phiName_),
    fanCurve_(pfopsf.fanCurve_().clone().ptr()),
    flowDir_(pfopsf.flowDir_),
    origin_(pfopsf.origin_),
    radialFlowDir_(pfopsf.radialFlowDir_),
    reverseFlow_(pfopsf.reverseFlow_),
    nonDimensional_(pfopsf.nonDimensional_),
    rpm_(pfopsf.rpm_),
    dm_(pfopsf.dm_),
    staticHead_(*this, pfopsf),
    relax_(pfopsf.relax_),
    minJump_(pfopsf.minJump_),
    maxJump_(pfopsf.maxJump_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fanPressureAMIFvPatchScalarField::updateCoeffs()
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

    if (reverseFlow_ == 0)
    {
        initializeFlowDirection();
    }


    // Average volumetric flow rate
    scalar volFlowRate = 0;

    if (phi.dimensions() == dimVelocity*dimArea)
    {
        volFlowRate = gSum(phip);
    }
    else if (phi.dimensions() == dimVelocity*dimArea*dimDensity)
    {
        const scalarField& rhop =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                this->db(),
                staticHead_.getRhoName()
            );
        volFlowRate = gSum(phip/rhop);
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
            120.0*volFlowRate/(pow3(constant::mathematical::pi)*dm_*rpm_);
    }
    volFlowRate *= reverseFlow_;

    // Pressure drop for this flow rate
    scalar pdFan = fanCurve_->value(max(volFlowRate, 0.0));

    if (nonDimensional_)
    {
        // Convert the adimensional deltap from curve into deltaP
        pdFan = pdFan*pow4(constant::mathematical::pi)*sqr(rpm_*dm_)/1800;
    }

    tmp<scalarField> pp(new scalarField(this->size(), reverseFlow_*pdFan));


    staticHead_.addStaticHead(pp.ref());

    min(pp.ref(), maxJump_);
    max(pp.ref(), minJump_);
    jump_ = relax_*pp.ref() + (1-relax_)*jump_;

    fvPatchScalarField::updateCoeffs();
}


void Foam::fanPressureAMIFvPatchScalarField::write(Ostream& os) const
{
    fixedJumpAMIFvPatchField<scalar>::write(os);
    if (!this->cyclicAMIPatch().owner())
    {
        jump_.writeEntry("jump", os);
    }
    referenceFrameFvPatch::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    if (!staticHead_.active())
    {
        os.writeEntry("rho", staticHead_.getRhoName());
    }
    os.writeEntry("directionMode", directionModeNames_[directionMode_]);
    switch (directionMode_)
    {
        case plane:
        {
            os.writeEntry("flowDirection", flowDir_);

            break;
        }
        case cylindrical:
        {
            os.writeEntry("origin", origin_);
            os.writeEntry("radialFlowDirection", radialFlowDir_);

            break;
        }
    }
    fanCurve_->writeData(os);
    os.writeEntry("nonDimensional", nonDimensional_);
    if (nonDimensional_)
    {
        os.writeEntry("rpm", rpm_);
        os.writeEntry("dm", dm_);
    }

    staticHead_.write(os);

    writeEntryIfDifferent<scalar>(os, "relax", 1.0, relax_);
    writeEntryIfDifferent<scalar>(os, "minJump", -GREAT, minJump_);
    writeEntryIfDifferent<scalar>(os, "maxJump", GREAT, maxJump_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fanPressureAMIFvPatchScalarField
    );
};


// ************************************************************************* //
