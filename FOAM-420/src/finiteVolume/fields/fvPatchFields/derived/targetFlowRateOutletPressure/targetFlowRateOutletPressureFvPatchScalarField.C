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
    (c) 2020-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/targetFlowRateOutletPressure/targetFlowRateOutletPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::targetFlowRateOutletPressureFvPatchScalarField::checkPhiDimensions
(
    const fvsPatchField<scalar>& phiP
) const
{
    const dimensionSet& phiDims = phiP.internalField().dimensions();
    const bool isPhiMassFlow = (phiDims == dimDensity*dimVelocity*dimArea);
    if (!isPhiMassFlow && (phiDims != dimVelocity*dimArea))
    {
        FatalErrorInFunction
            << phiName_ << " field dimensions are not consistent"
            << abort(FatalError);
    }
    return isPhiMassFlow;
}


Foam::scalar
Foam::targetFlowRateOutletPressureFvPatchScalarField::calcTargetFlowRate
(
    const bool isPhiMassFlow
) const
{
    scalar tarFlux = currentFlowRate();
    if (!isSplitBC())
    {
        if (isPhiMassFlow)
        {
            if (volumetric())
            {
                tarFlux *= rho_;
            }
        }
        else
        {
            if (!volumetric())
            {
                tarFlux /= rho_;
            }
        }
    }
    return tarFlux;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::targetFlowRateOutletPressureFvPatchScalarField::
targetFlowRateOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    flowRateBase(patch(), db(), true),
    phiName_("phi"),
    relax_(0.05),
    minP_(-GREAT),
    maxP_(GREAT),
    rho_(0.0)
{
    this->refValue() = 0.0;
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::targetFlowRateOutletPressureFvPatchScalarField::
targetFlowRateOutletPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF, dict, false),
    flowRateBase(patch(), db(), dict, true),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.05)),
    minP_(dict.lookupOrDefault<scalar>("minP", -GREAT)),
    maxP_(dict.lookupOrDefault<scalar>("maxP",  GREAT)),
    rho_(0.0)
{
    frameName_ = dict.lookupOrDefault<word>("referenceFrame", "");
    inletFlux_ = true;
    if (dict.found("value"))
    {
        this->refValue() = scalarField("value", dict, p.size());
        fvPatchField<scalar>::operator= (this->refValue());
    }
    else
    {
        this->refValue() = 0.0;
    }
    this->refGrad() = 0.0;
    this->valueFraction() = 0.0;
}


Foam::targetFlowRateOutletPressureFvPatchScalarField::
targetFlowRateOutletPressureFvPatchScalarField
(
    const targetFlowRateOutletPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    flowRateBase(patch(), db(), ptf, true),
    phiName_(ptf.phiName_),
    relax_(ptf.relax_),
    minP_(ptf.minP_),
    maxP_(ptf.maxP_),
    rho_(0.0)
{}


Foam::targetFlowRateOutletPressureFvPatchScalarField::
targetFlowRateOutletPressureFvPatchScalarField
(
    const targetFlowRateOutletPressureFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    flowRateBase(patch(), db(), ptf, true),
    phiName_(ptf.phiName_),
    relax_(ptf.relax_),
    minP_(ptf.minP_),
    maxP_(ptf.maxP_),
    rho_(0.0)
{}


Foam::targetFlowRateOutletPressureFvPatchScalarField::
targetFlowRateOutletPressureFvPatchScalarField
(
    const targetFlowRateOutletPressureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    flowRateBase(patch(), db(), ptf, true),
    phiName_(ptf.phiName_),
    relax_(ptf.relax_),
    minP_(ptf.minP_),
    maxP_(ptf.maxP_),
    rho_(0.0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::targetFlowRateOutletPressureFvPatchScalarField::updateCoeffs()
{

    if (updated())
    {
        return;
    }

    if (this->internalField().name() == "pPot")
    {
        //- hack to by-pass the caseSetup potential solver.
        //  transform the boundary to zero-Neumann (artificial wall pressure BC)
        //  Otherwise the pot solver will never converge and air will be
        //  sucked/blowed from the boundary
        //  Caution: if potential solving field changes, this must be modified.
        //  Update 2022: the coeffs reversed back. So now in caseSetup the BC
        //  behaves as a fixedPressure with a value of the specified one.
        //  Careful for the value compressible and incompressible solvers
//        this->valueFraction() = 0.0;
//        this->refGrad() = 0.0;
        this->valueFraction() = 1.0;
        this->refGrad() = 0.0;
    }
    else
    {
        fvsPatchField<scalar> phiP =
            patch().lookupPatchFieldInDb<surfaceScalarField, scalar>
            (
                this->db(),
                phiName_
            );
        const bool isPhiMassFlow = checkPhiDimensions(phiP);

        // Load rho for compressible
        const fvPatchField<scalar>* rhoP =
            isPhiMassFlow
          ? &patch().lookupPatchFieldInDb<volScalarField, scalar>
            (
                this->db(),
                this->rhoName()
            )
          : nullptr;

        if (!this->db().foundObject<foamCoupledControl>(solutionControl::typeName))
        {
            if (isPhiMassFlow)
            {
                scalarField relPhiP(phiP.size(), 0.0);
                makeRelative(relPhiP);
                phiP += (*rhoP)*relPhiP;
            }
            else
            {
                makeRelative(phiP);
            }
        }

        const scalar curFlux = gSum(phiP);
        const scalar surface = gSum(patch().magSf());
        rho_ = calcRho(this->rhoName(), rhoP, surface);

        scalar tarFlux = calcTargetFlowRate(isPhiMassFlow);

        scalar dp =
            0.5*(sign(curFlux)*sqr(curFlux) - sqr(tarFlux))/sqr(surface);

        if (isPhiMassFlow)
        {
            dp /= rho_;
        }

        scalarField& p = *this;
        p += relax_*dp;

        //- limiters to bound values (if required)
        if (minP_ != -GREAT)
        {
            forAll(p, fI)
            {
                p[fI] = max(minP_, p[fI]);
            }
        }
        if (maxP_ != GREAT)
        {
            forAll(p, fI)
            {
                p[fI] = min(maxP_, p[fI]);
            }
        }

        this->refValue() = p;
        this->valueFraction() = 1.0;
    }

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::targetFlowRateOutletPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchScalarField::write(os);
    flowRateBase::write(os);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<scalar>(os, "relax", 0.05, relax_);
    writeEntryIfDifferent<scalar>(os, "minP", -GREAT, minP_);
    writeEntryIfDifferent<scalar>(os, "maxP", GREAT, maxP_);
}


Foam::scalar Foam::targetFlowRateOutletPressureFvPatchScalarField::calcRho
(
    const word& rhoName,
    const fvPatchField<scalar>* rhop,
    const scalar surface
)
{
    if (rhop)
    {
        return gSum((*rhop)*(*rhop).patch().magSf())/surface;
    }
    else if (rho_ < SMALL)
    {
        IOdictionary tp
        (
            IOobject
            (
                "transportProperties",
                db().time().constant(),
                db().time(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        if (!tp.found(rhoName))
        {
            WarningInFunction
                << "Cannot find rho in dictionary " << tp.name()
                << ".\n Disable rho multiplication. "
                << "This message should only be shown when "
                << "using buoyantBoussinesq(S/P)impleFoam!" << endl;
            return 1.0;
        }
        else
        {
            const scalar rho =
                dimensionedScalar("rho", tp.lookup(rhoName)).value();
            if (rho < SMALL)
            {
                return 1.0;
            }
            return rho;
        }
    }

    return rho_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        targetFlowRateOutletPressureFvPatchScalarField
    );
}

// ************************************************************************* //
