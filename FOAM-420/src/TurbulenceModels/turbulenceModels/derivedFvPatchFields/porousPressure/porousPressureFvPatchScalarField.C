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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "porousPressureFvPatchScalarField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char*
    NamedEnum
    <
        porousPressureFvPatchScalarField::porousMode,
        3
    >::names[] =
    {
        "darcyForchheimer",
        "alphaBeta",
        "powerLaw"
    };
}

const Foam::NamedEnum
<
    Foam::porousPressureFvPatchScalarField::porousMode,
    3
> Foam::porousPressureFvPatchScalarField::porousModeNames_;


// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

void Foam::porousPressureFvPatchScalarField::updateSegBasedCoeffs()
{
    if (updated())
    {
        return;
    }

    word rAUlookup = rAUName_;
    if (!db().foundObject<volScalarField>(rAUlookup))
    {
        rAUlookup = "(1|A(U))";
        if (!db().foundObject<volScalarField>(rAUlookup))
        {
            word addMessage = word::null;
            if (rAUlookup != rAUlookup)
            {
                addMessage = ", neither ";
                addMessage += rAUlookup;
            }
            WarningInFunction
                << "Could not find specified rUA field: " << rAUName_
                << addMessage << ". "
                << "Boundary update skipped..."
                << endl;
            fixedValueFvPatchScalarField::updateCoeffs();
            return;
        }
        else
        {
            WarningInFunction
                << "Could not find specified rUA field: " << rAUName_
                << ", defaulting to " << rAUlookup << " for this update."
                << endl;
        }
    }

    const fvPatchField<scalar>& rAp =
        patch().lookupPatchFieldInDb<volScalarField, scalar>
        (
            db(),
            rAUlookup
        );

    vectorField Up
    (
        db().lookupObject<volVectorField>
        (
            UName()
        ).boundaryField()[patch().index()]
    );

    makeRelative(Up);

    scalarField Un(Up & patch().nf()); //(phip / patch().magSf());

    scalarField sqrUt(magSqr(Up - Un*patch().nf()));

    scalarField& p = *this;

    scalarField pi(this->patchInternalField());

    // calculate base gradient before modifying pressure
    scalarField dp0(this->snGrad());

    const scalarField& rd = patch().deltaCoeffs();

    scalarField dUn(this->size(), 0.0);

    //iterative pressure solver based on
    // Un_new = Un_old - 1/A * (grad(p) - grad(p_old))
    // and
    // p = p0 + (C1 + 0.5*C2*|Un|)Un + 0.5*Ct*sqr(Ut)
    // rho is embedded in the C1 C2 for compressible
    tmp<scalarField> tC1 = C1Field();
    const scalarField& C1 = tC1();

    tmp<scalarField> tC2 = C2Field();
    scalarField& C2 = tC2.ref();

    tmp<scalarField> tC3 = C3Field();
    const scalarField& C3 = tC3();

    const scalarField& p0f = p0();
    forAll(p, facei)
    {
        label n = 0;
        scalar err = GREAT;

        do
        {
            scalar pold = p[facei];

            dUn[facei] =
                - (p[facei] - pi[facei])*rAp[facei]*rd[facei]
                + dp0[facei]*rAp[facei];

            scalar UN = Un[facei] + dUn[facei];
            scalar f =
                p[facei] - p0f[facei]
              - (C1[facei] + C2[facei]*mag(UN))*UN
              - C3[facei]*sqrUt[facei];

            scalar df =
                1
              + C2[facei]*UN*(UN*rd[facei]*rAp[facei])/max(mag(UN), VSMALL)
              + C1[facei]*rAp[facei]*rd[facei]
              + C2[facei]*rAp[facei]*rd[facei]*mag(UN);

            p[facei] = pold - f/(sign(df)*max(mag(df), VSMALL));

            err =
                mag(pold - p[facei])/max(SMALL, max(mag(p[facei]), mag(pold)));

            n++;

            if (n == 100)
            {
                WarningInFunction
                    << "Boundary pressure solver did not converge. error ="
                    << err << endl;
            }

        } while (err > 0.000001 && n < 100);

    }

    //convergence acceleration
    //correct outlet pressure and patch adjacent pressure using means
    //to accelerate convergence of boundary
    const scalar patchArea = gSum(patch().magSf());
    const scalar pMean = gSum(p*patch().magSf())/patchArea;
    const scalar tpMean =
        gSum
        (
            (p0f + (C1 + C2*mag(Un))*Un + C3*sqrUt)*patch().magSf()
        )/patchArea;

    p += - 0.1 * (pMean - tpMean);

    fixedValueFvPatchScalarField::updateCoeffs();
}


Foam::scalar Foam::porousPressureFvPatchScalarField::calcLinearValue
(
    const scalar& T
) const
{
    return (aT1_ + aT2_*TRef_ + aT3_*sqr(TRef_))*(TRef_/T);
}


Foam::scalar Foam::porousPressureFvPatchScalarField::calcNonLinearValue
(
    const scalar& T
) const
{
    return (bT1_ + bT2_*TRef_ + bT3_*sqr(TRef_))*sqr(TRef_/T);
}


void Foam::porousPressureFvPatchScalarField::
correctTemperatureDepLinearTerm
(
    scalarField& linearCoeff
) const
{
    if (temperatureDependence_)
    {
        if (patch().boundaryMesh().mesh().foundObject<volScalarField>(TName_))
        {
            const fvPatchField<scalar>& Tp =
                patch().lookupPatchFieldInDb<volScalarField, scalar>
                (
                    db(),
                    TName_
                );

            forAll(Tp, facei)
            {
                linearCoeff[facei] *= calcLinearValue(Tp[facei]);
            }
        }
    }
}


void Foam::porousPressureFvPatchScalarField::correctTemperatureDepNonLinearTerm
(
    scalarField& nonLinearCoeff
) const
{
    if (temperatureDependence_)
    {
        if (patch().boundaryMesh().mesh().foundObject<volScalarField>(TName_))
        {
            const fvPatchField<scalar>& Tp =
                patch().lookupPatchFieldInDb<volScalarField, scalar>
                (
                    db(),
                    TName_
                );

            forAll(Tp, facei)
            {
                nonLinearCoeff[facei] *= calcNonLinearValue(Tp[facei]);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousPressureFvPatchScalarField::porousPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    basePressureFvPatchScalarField(p, iF),
    porousMode_(darcyForchheimer),
    d_(),
    f_(),
    tangentf_(),
    alpha_(),
    beta_(),
    tangentBeta_(),
    C0_(),
    C1_(),
    length_(0),
    relax_(0.7),
    useAveragedRho_(true),
    oldRho_(patch().size(), 1.0),
    isOldRhoInit_(false),
    rho_(patch().size(), 1.0),
    rhoRelax_(0.1),
    rAUName_("rAU"),
    temperatureDependence_(false),
    TName_(word::null),
    TRef_(0.0),
    aT1_(0.0),
    aT2_(0.0),
    aT3_(0.0),
    bT1_(0.0),
    bT2_(0.0),
    bT3_(0.0)
{}


Foam::porousPressureFvPatchScalarField::porousPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    basePressureFvPatchScalarField(p, iF, dict),
    porousMode_
    (
        porousModeNames_.readOrDefault
        (
            dict.lookupOrDefault<word>
            (
                "porousMode",
                "darcyForchheimer"
            ),
            darcyForchheimer
        )
    ),
    d_(),
    f_(),
    tangentf_(),
    alpha_(),
    beta_(),
    tangentBeta_(),
    C0_(),
    C1_(),
    length_(readScalar(dict.lookup("length"))),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.7)),
    useAveragedRho_(dict.lookupOrDefault<Switch>("useAveragedRho", true)),
    oldRho_(patch().size(), 1.0),
    isOldRhoInit_(false),
    rho_
    (
        dict.found("rhoInit")
      ? scalarField("rhoInit", dict, patch().size())
      : scalarField(patch().size(), 1.0)
    ),
    rhoRelax_(dict.lookupOrDefault<scalar>("rhoRelax", 0.1)),
    rAUName_(dict.lookupOrDefault<word>("rAU", "(1|A(U))")),
    temperatureDependence_
    (
        dict.lookupOrDefault<Switch>("temperatureDependence", "false")
    ),
    TName_(word::null),
    TRef_(0.0),
    aT1_(0.0),
    aT2_(0.0),
    aT3_(0.0),
    bT1_(0.0),
    bT2_(0.0),
    bT3_(0.0)
{
    inletFlux_ = true;
    switch (porousMode_)
    {
        case darcyForchheimer:
        {
            d_ = Function1<scalar>::New("d", dict);
            f_ = Function1<scalar>::New("f", dict);
            if (dict.found("tangentf"))
            {
                FatalErrorInFunction
                    << "Tangent \"tangentf\" componnent isn't supported yet."
                    << exit(FatalError);
                tangentf_ = Function1<scalar>::New("tangentf", dict);
            }

            break;
        }
        case alphaBeta:
        {
            alpha_ = Function1<scalar>::New("alpha", dict);
            beta_ = Function1<scalar>::New("beta", dict);
            if (dict.found("tangentBeta"))
            {
                FatalErrorInFunction
                    << "Tangent \"tangentBeta\" "
                    << "componnent isn't supported yet."
                    << exit(FatalError);
                tangentBeta_ = Function1<scalar>::New("tangentBeta", dict);
            }

            break;
        }
        case powerLaw:
        {
            C0_ = Function1<scalar>::New("C0", dict);
            C1_ = Function1<scalar>::New("C1", dict);

            break;
        }
    }
    if (temperatureDependence_)
    {
        TName_ = dict.lookupOrDefault<word>("T", "T");
        TRef_ = dict.lookupOrDefault<scalar>("TRef", 0.0);
        aT1_ = dict.lookupOrDefault<scalar>("a1", 0.0);
        aT2_ = dict.lookupOrDefault<scalar>("a2", 0.0);
        aT3_ = dict.lookupOrDefault<scalar>("a3", 0.0);
        bT1_ = dict.lookupOrDefault<scalar>("b1", 0.0);
        bT2_ = dict.lookupOrDefault<scalar>("b2", 0.0);
        bT3_ = dict.lookupOrDefault<scalar>("b3", 0.0);
    }
}


Foam::porousPressureFvPatchScalarField::porousPressureFvPatchScalarField
(
    const porousPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basePressureFvPatchScalarField(ptf, p, iF, mapper),
    porousMode_(ptf.porousMode_),
    d_(ptf.d_, false),
    f_(ptf.f_, false),
    tangentf_(ptf.tangentf_, false),
    alpha_(ptf.alpha_, false),
    beta_(ptf.beta_, false),
    tangentBeta_(ptf.tangentBeta_, false),
    C0_(ptf.C0_, false),
    C1_(ptf.C1_, false),
    length_(ptf.length_),
    relax_(ptf.relax_),
    useAveragedRho_(ptf.useAveragedRho_),
    oldRho_(ptf.oldRho_),
    isOldRhoInit_(false),
    rho_(ptf.rho_),
    rhoRelax_(ptf.rhoRelax_),
    rAUName_(ptf.rAUName_),
    temperatureDependence_(ptf.temperatureDependence_),
    TName_(ptf.TName_),
    TRef_(ptf.TRef_),
    aT1_(ptf.aT1_),
    aT2_(ptf.aT2_),
    aT3_(ptf.aT3_),
    bT1_(ptf.bT1_),
    bT2_(ptf.bT2_),
    bT3_(ptf.bT3_)
{}


Foam::porousPressureFvPatchScalarField::porousPressureFvPatchScalarField
(
    const porousPressureFvPatchScalarField& tppsf
)
:
    basePressureFvPatchScalarField(tppsf),
    porousMode_(tppsf.porousMode_),
    d_(tppsf.d_, false),
    f_(tppsf.f_, false),
    tangentf_(tppsf.tangentf_, false),
    alpha_(tppsf.alpha_, false),
    beta_(tppsf.beta_, false),
    tangentBeta_(tppsf.tangentBeta_, false),
    C0_(tppsf.C0_, false),
    C1_(tppsf.C1_, false),
    length_(tppsf.length_),
    relax_(tppsf.relax_),
    useAveragedRho_(tppsf.useAveragedRho_),
    oldRho_(tppsf.oldRho_),
    isOldRhoInit_(false),
    rho_(tppsf.rho_),
    rhoRelax_(tppsf.rhoRelax_),
    rAUName_(tppsf.rAUName_),
    temperatureDependence_(tppsf.temperatureDependence_),
    TName_(tppsf.TName_),
    TRef_(tppsf.TRef_),
    aT1_(tppsf.aT1_),
    aT2_(tppsf.aT2_),
    aT3_(tppsf.aT3_),
    bT1_(tppsf.bT1_),
    bT2_(tppsf.bT2_),
    bT3_(tppsf.bT3_)
{}


Foam::porousPressureFvPatchScalarField::porousPressureFvPatchScalarField
(
    const porousPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    basePressureFvPatchScalarField(tppsf, iF),
    porousMode_(tppsf.porousMode_),
    d_(tppsf.d_, false),
    f_(tppsf.f_, false),
    tangentf_(tppsf.tangentf_, false),
    alpha_(tppsf.alpha_, false),
    beta_(tppsf.beta_, false),
    tangentBeta_(tppsf.tangentBeta_, false),
    C0_(tppsf.C0_, false),
    C1_(tppsf.C1_, false),
    length_(tppsf.length_),
    relax_(tppsf.relax_),
    useAveragedRho_(tppsf.useAveragedRho_),
    oldRho_(tppsf.oldRho_),
    isOldRhoInit_(false),
    rho_(tppsf.rho_),
    rhoRelax_(tppsf.rhoRelax_),
    rAUName_(tppsf.rAUName_),
    temperatureDependence_(tppsf.temperatureDependence_),
    TName_(tppsf.TName_),
    TRef_(tppsf.TRef_),
    aT1_(tppsf.aT1_),
    aT2_(tppsf.aT2_),
    aT3_(tppsf.aT3_),
    bT1_(tppsf.bT1_),
    bT2_(tppsf.bT2_),
    bT3_(tppsf.bT3_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousPressureFvPatchScalarField::multiplyByRho
(
    scalarField& coeff
) const
{
    if (internalField().dimensions() == dimPressure)
    {
        const scalarField& rhop =
            db().lookupObject<volScalarField>
            (
                rhoName()
            ).boundaryField()[patch().index()];

        // First time the rho is used initialise old rho to current rho
        if (!isOldRhoInit_)
        {
            if (useAveragedRho_)
            {
                oldRho_ = gSum(rhop*patch().magSf())/gSum(patch().magSf());
            }
            else
            {
                oldRho_ = rhop;
            }

            isOldRhoInit_ = true;
        }

        if (useAveragedRho_)
        {
            rho_ = gSum(rhop*patch().magSf())/gSum(patch().magSf());
        }
        else
        {
            rho_ = rhop;
        }

        if (rhoRelax_ < 1.0)
        {
            rho_ = (oldRho_ + (rho_ - oldRho_)*rhoRelax_);
        }

        coeff *= rho_;
    }
}


Foam::tmp<Foam::scalarField>
Foam::porousPressureFvPatchScalarField::C0Field() const
{
    return tmp<scalarField>(new scalarField(this->size(), Zero));
}


Foam::tmp<Foam::scalarField>
Foam::porousPressureFvPatchScalarField::C1Field() const
{
    const scalar t = db().time().timeOutputValue();
    tmp<scalarField> tC1(new scalarField(this->size(), length_));
    scalarField& C1 = tC1.ref();
    switch (porousMode_)
    {
        case darcyForchheimer:
        {
            const word turbName
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    internalField().group()
                )
            );
            if (db().found(turbName))
            {
                const turbulenceModel& turbModel =
                    db().lookupObject<turbulenceModel>(turbName);
                C1 *= d_->value(t)*turbModel.nu(patch().index());
                correctTemperatureDepLinearTerm(C1);
            }
            else
            {
                WarningInFunction
                    << "Not using viscosity for porosity on the boundary "
                    << this->patch().name() << " should be reported only when "
                    << "turbulence model isn't available." << endl;
                C1 *= 0.0;
            }
            break;
        }
        case alphaBeta:
        {
            C1 *= alpha_->value(t);
            correctTemperatureDepLinearTerm(C1);
            break;
        }
        case powerLaw:
        {
            const vectorField Up =
                patch().lookupPatchField<volVectorField, scalar>(UName());
            C1 *= C0_->value(t)*pow(magSqr(Up), C1_->value(t) - 1.0);
            break;
        }
    }

    multiplyByRho(C1);

    return tC1;
}


Foam::tmp<Foam::scalarField>
Foam::porousPressureFvPatchScalarField::C2Field() const
{
    const scalar t = db().time().timeOutputValue();
    tmp<scalarField> tC2
    (
        new scalarField(this->size(), length_)
    );
    scalarField& C2 = tC2.ref();
    switch (porousMode_)
    {
        case darcyForchheimer:
        {
            C2 *= 0.5*f_->value(t);
            correctTemperatureDepNonLinearTerm(C2);
            break;
        }
        case alphaBeta:
        {
            C2 *= beta_->value(t);
            correctTemperatureDepNonLinearTerm(C2);
            break;
        }
        case powerLaw:
        {
            C2 = Zero;
            break;
        }
    }

    multiplyByRho(C2);

    return tC2;
}


Foam::tmp<Foam::scalarField>
Foam::porousPressureFvPatchScalarField::C3Field() const
{
    // TODO Proper validation of tangent component is needed that includes
    // check if the temperature correction is correct
    // Most likely will be removed in future
    tmp<scalarField> tC3(new scalarField(this->size(), 0.0));
    // scalarField& C3 = tC3.ref();
    // const scalar t = db().time().timeOutputValue();
    // switch (porousMode_)
    // {
    //     case darcyForchheimer:
    //     {
    //         if (tangentf_.valid())
    //         {
    //             C3 = 0.5*length_*tangentf_->value(t);
    //             correctTemperatureDepNonLinearTerm(C3);
    //         }
    //         break;
    //     }
    //     case alphaBeta:
    //     {
    //         if (tangentBeta_.valid())
    //         {
    //             C3 = length_*tangentBeta_->value(t);
    //             correctTemperatureDepNonLinearTerm(C3);
    //         }
    //         break;
    //     }
    //     case powerLaw:
    //     {
    //         break;
    //     }
    // }

    // multiplyByRho(C3);

    return tC3;
}


void Foam::porousPressureFvPatchScalarField::updateCoeffs()
{
    // Update old rho used for relaxation of rho
    if (!updated())
    {
        if (internalField().dimensions() == dimPressure)
        {
            oldRho_ = rho_;
        }
    }

    if (isCoupledSolver())
    {
        basePressureFvPatchScalarField::updateCoeffs();
    }
    else
    {
        updateSegBasedCoeffs();
    }
}


void Foam::porousPressureFvPatchScalarField::boundaryRelaxMatrix
(
    fvBlockMatrix<vector>& bEq
) const
{
    if (relax_ < 1)
    {
        bEq.boundaryRelax(relax_, this->patch().index());
    }
}


void Foam::porousPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    referenceFrameFvPatch::write(os);
    p0().writeEntry("p0", os);
    os.writeEntry("porousMode", porousModeNames_[porousMode_]);

    writeEntryIfDifferent<word>(os, "U", "U", UName());
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName());
    switch (porousMode_)
    {
        case darcyForchheimer:
        {
            d_->writeData(os);
            f_->writeData(os);
            if (tangentf_.valid())
            {
                tangentf_->writeData(os);
            }
            break;
        }
        case alphaBeta:
        {
            alpha_->writeData(os);
            beta_->writeData(os);
            if (tangentBeta_.valid())
            {
                tangentBeta_->writeData(os);
            }
            break;
        }
        case powerLaw:
        {
            C0_->writeData(os);
            C1_->writeData(os);
            break;
        }
    }
    os.writeEntry("length", length_);
    writeEntryIfDifferent<scalar>(os, "relax", 0.7, relax_);

    // Write rho related variables
    writeEntryIfDifferent<Switch>(os, "useAverageRho", true, useAveragedRho_);
    writeEntryIfDifferent<scalar>(os, "rhoRelax", 0.1, rhoRelax_);
    rho_.writeEntry("rhoInit", os);

    writeEntryIfDifferent<word>(os, "rAU", "(1|A(U))", rAUName_);
    writeEntryIfDifferent<Switch>
    (
        os,
        "temperatureDependence",
        false,
        temperatureDependence_
    );
    if (temperatureDependence_)
    {
        if (porousMode_ == powerLaw)
        {
            FatalErrorInFunction
                << "Temperature dependency is not supported in  "
                << powerLaw
                << " model."
                << exit(FatalError);
        }
        os.writeEntry("T", TName_);
        os.writeEntry("TRef", TRef_);
        os.writeEntry("a1", aT1_);
        os.writeEntry("a2", aT2_);
        os.writeEntry("a3", aT3_);
        os.writeEntry("b1", bT1_);
        os.writeEntry("b2", bT2_);
        os.writeEntry("b3", bT3_);
    }
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        porousPressureFvPatchScalarField
    );
}

// ************************************************************************* //
