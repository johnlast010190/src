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
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/flowRateOutletVelocity/flowRateOutletVelocityFvPatchVectorField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "primitives/functions/Function1/Constant/Constant.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


namespace Foam
{

void flowRateOutletVelocityFvPatchVectorField::scaleFlux
(
    const scalar flowRate,
    scalarField& phiPN
) const
{
    //relax velocity boundary update

    //we will block normal reversed flow, scale accordingly
    phiPN *= pos0(sign(flowRate)*phiPN);

    const scalar cFlowRate = gSum(phiPN);

    const scalar flowFrac
        = mag(cFlowRate)/max(mag(flowRate), SMALL);

    scalar Atotal = gSum(patch().magSf());

    scalar addFlux(0);
    scalar multiplyFlux(1.0);

    //ramp from current to target velocity
    if (flowFrac < 0.5 || sign(cFlowRate) != sign(flowRate))
    {
        //use additive when far below target
        addFlux = maxDelta_ * (flowRate - cFlowRate);
    }
    else
    {
        //use scaling when close to or above target
        multiplyFlux = min
        (
            1 + maxDelta_,
            max
            (
                1 - maxDelta_,
                flowRate/(sign(cFlowRate)*max(mag(cFlowRate), SMALL))
            )
        );
    }

    // modify fluxes closer to target value
    phiPN = multiplyFlux*phiPN + addFlux*patch().magSf()/Atotal;

}

tmp<scalarField> flowRateOutletVelocityFvPatchVectorField::nonUniformVelocity
(
    const scalar flowRate
) const
{
    // zero-gradient derived flux field
    tmp<scalarField> phiPN;

    // patch area
    scalar Atotal = gSum(patch().magSf());

    phiPN = multiplyByRhoOrOne(patch().Sf() & this->patchInternalField(), true);

    scaleFlux(flowRate, phiPN.ref());

    if (debug && phaseName() != "none")
    {
        const fvPatchField<scalar>& phasep =
            patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                db(), phaseName()
            );

        Info<< "current flow rate: " << gSum(phiPN()) << nl
            << "current phase flow rate: " << gSum(phiPN()*phasep) << nl
            << "phase wetting: " << gSum(phasep*patch().magSf())/Atotal << endl;
    }

    //TODO: Preserving prior behaviour by not dividing here by average rho if
    // selected - is this really intended?
    phiPN.ref() /= multiplyByRhoOrOne(patch().magSf(), false);

    return phiPN;
}

tmp<scalarField> flowRateOutletVelocityFvPatchVectorField::uniformVelocity
(
    const scalar flowRate
) const
{
    // zero-gradient derived flux field
    tmp<scalarField> phiPN;

    // patch area
    const scalar Atotal = gSum(patch().magSf());
    const scalar avgRhoU = flowRate/Atotal;

    phiPN = avgRhoU/multiplyByRhoOrOne(scalarField(this->size(), 1));

    return phiPN;
}

void flowRateOutletVelocityFvPatchVectorField::initialise
(
    const dictionary& dict
)
{
    maxDelta_ = min(1, max(0, maxDelta_));

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, patch().size())
        );

        tmp<vectorField> n = patch().nf();
        fixedNormalSlipFvPatchVectorField::fixedValue_
            = (*this & n())*n();

        fixedNormalSlipFvPatchVectorField::updateCoeffs();
    }
    else
    {
        evaluate(Pstream::commsTypes::blocking);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

flowRateOutletVelocityFvPatchVectorField::
flowRateOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedNormalSlipFvPatchVectorField(p, iF),
    flowRateBase(patch(), db(), false),
    currentFlowRate_(0),
    maxDelta_(0.2),
    phiName_("phi"),
    uniform_(false),
    fluxValueSwitch_(true),
    relax_(0.2)
{}


flowRateOutletVelocityFvPatchVectorField::
flowRateOutletVelocityFvPatchVectorField
(
    const flowRateOutletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedNormalSlipFvPatchVectorField(ptf, p, iF, mapper),
    flowRateBase(patch(), db(), ptf, false),
    currentFlowRate_(ptf.currentFlowRate_),
    maxDelta_(ptf.maxDelta_),
    phiName_(ptf.phiName_),
    uniform_(ptf.uniform_),
    fluxValueSwitch_(ptf.fluxValueSwitch_),
    relax_(ptf.relax_)
{}


flowRateOutletVelocityFvPatchVectorField::
flowRateOutletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedNormalSlipFvPatchVectorField(p, iF),
    flowRateBase(patch(), db(), dict, false),
    currentFlowRate_(dict.lookupOrDefault<scalar>("currentFlowRate", 0.0)),
    maxDelta_(dict.lookupOrDefault<scalar>("maxDelta", 0.2)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    uniform_(dict.lookupOrDefault<Switch>("uniform", false)),
    fluxValueSwitch_(dict.lookupOrDefault<Switch>("fluxValueSwitch", true)),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.2))
{
    initialise(dict);
}




flowRateOutletVelocityFvPatchVectorField::
flowRateOutletVelocityFvPatchVectorField
(
    const flowRateOutletVelocityFvPatchVectorField& ptf
)
:
    fixedNormalSlipFvPatchVectorField(ptf),
    flowRateBase(patch(), db(), ptf, false),
    flowRate_(ptf.flowRate_().clone().ptr()),
    currentFlowRate_(ptf.currentFlowRate_),
    maxDelta_(ptf.maxDelta_),
    phiName_(ptf.phiName_),
    uniform_(ptf.uniform_),
    fluxValueSwitch_(ptf.fluxValueSwitch_),
    relax_(ptf.relax_)
{}


flowRateOutletVelocityFvPatchVectorField::
flowRateOutletVelocityFvPatchVectorField
(
    const flowRateOutletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedNormalSlipFvPatchVectorField(ptf, iF),
    flowRateBase(patch(), db(), ptf, false),
    currentFlowRate_(ptf.currentFlowRate_),
    maxDelta_(ptf.maxDelta_),
    phiName_(ptf.phiName_),
    uniform_(ptf.uniform_),
    fluxValueSwitch_(ptf.fluxValueSwitch_),
    relax_(ptf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<Field<vector>> flowRateOutletVelocityFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{

    tmp<Field<vector>> vic
    (
        new Field<vector>(this->size(), pTraits<vector>::one)
    );

    if (fluxValueSwitch_)
    {
        const Field<scalar>& phip =
            this->patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
            (
                this->db(),
                phiName_
            );

        vic.ref() *= pos0(phip);
    }
    else
    {
        vic.ref() *= pos0(currentFlowRate_);
    }

    return vic;
}

tmp<Field<vector>> flowRateOutletVelocityFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{

    tmp<Field<vector>> vbc
    (
        new Field<vector>(*this)
    );

    if (fluxValueSwitch_)
    {
        const Field<scalar>& phip =
            this->patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
            (
                this->db(),
                phiName_
            );

        vbc.ref() *= neg(phip);
    }
    else
    {
        vbc.ref() *= neg(currentFlowRate_);
    }

    return vbc;
}

Foam::tmp<Foam::Field<Foam::vector>>
Foam::flowRateOutletVelocityFvPatchVectorField
::valueDivInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<Field<vector>> vict
    (
        new Field<vector>(this->size(), Zero)
    );

    return vict;
}

Foam::tmp<Foam::Field<Foam::vector>>
Foam::flowRateOutletVelocityFvPatchVectorField
::valueDivBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    tmp<Field<vector>> vbct
    (
        new Field<vector>(*this)
    );

    return vbct;
}




void flowRateOutletVelocityFvPatchVectorField::updateCoeffs
(
    const scalar flowRate
)
{
    //store current flow rate for switching schemes
    currentFlowRate_ = flowRate;

    //normal velocity magnitude
    tmp<scalarField> Upn;

    //populate the zero-gradient derived flux field
    //non-uniform (scaled zero gradient), very slow for inlet
    if (!uniform_)
    {
        Upn = nonUniformVelocity(flowRate);
    }
    else //slug flow
    {
        Upn = uniformVelocity(flowRate);
    }

    fixedNormalSlipFvPatchVectorField::fixedValue_ = Upn()*this->patch().nf();
    fixedNormalSlipFvPatchVectorField::updateCoeffs();
}


void flowRateOutletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    updateCoeffs(currentFlowRate());
}


void Foam::flowRateOutletVelocityFvPatchVectorField::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>& bEq
) const
{
    if (relax_==1)
    {
        return;
    }
    bEq.boundaryRelax(relax_, this->patch().index());
}



void flowRateOutletVelocityFvPatchVectorField::write(Ostream& os) const
{
    transformFvPatchField<vector>::write(os);
    flowRateBase::write(os);

    os.writeEntry("currentFlowRate", currentFlowRate_);

    writeEntryIfDifferent<scalar>(os, "maxDelta", 0.2, maxDelta_);

    writeEntryIfDifferent<Switch>(os, "uniform", false, uniform_);
    writeEntryIfDifferent<Switch>(os, "fluxValueSwitch", true, fluxValueSwitch_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);

    writeEntryIfDifferent<scalar>(os, "relax", 0.2, relax_);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    flowRateOutletVelocityFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
