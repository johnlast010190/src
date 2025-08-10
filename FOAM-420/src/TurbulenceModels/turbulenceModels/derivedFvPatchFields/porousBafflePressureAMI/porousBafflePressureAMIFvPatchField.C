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
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "porousBafflePressureAMIFvPatchField.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "turbulenceModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "finiteVolume/fvc/fvcMeshPhi.H"
#include "finiteVolume/fvc/fvcFlux.H"
#include "finiteVolume/fvc/fvcReconstruct.H"
#include "cfdTools/general/solutionControl/foamCoupledControl/foamCoupledControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char*
    NamedEnum
    <
        porousBafflePressureAMIFvPatchField::porousMode,
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
    Foam::porousBafflePressureAMIFvPatchField::porousMode,
    3
> Foam::porousBafflePressureAMIFvPatchField::porousModeNames_;


// * * * * * * * * * * * * * Private Functions   * * * * * * * * * * * * * * //

Foam::scalar Foam::porousBafflePressureAMIFvPatchField::calcDcorrValue
(
    const scalar& T
) const
{
    return (aT1_ + aT2_*TRef_ + aT3_*sqr(TRef_))*(TRef_/T);
}


Foam::scalar Foam::porousBafflePressureAMIFvPatchField::calcFcorrValue
(
    const scalar& T
) const
{
    return (bT1_ + bT2_*TRef_ + bT3_*sqr(TRef_))*sqr(TRef_/T);
}

void Foam::porousBafflePressureAMIFvPatchField::stabiliseJump
(
    const scalarField& jump0
)
{
    if (!uniformJump_ && (maxJump_!=GREAT))
    {
        min(jump_, maxJump_);
    }
    if (!uniformJump_ && (minJump_!=-GREAT))
    {
        max(jump_, minJump_);
    }
    if (relax_<1.0)
    {
       jump_ = relax_*jump_+ (1-relax_)*jump0;
    }
}

Foam::tmp<Foam::scalarField>
Foam::porousBafflePressureAMIFvPatchField::computeUn() const
{
    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    const fvsPatchField<scalar>& phip =
        patch().patchField<surfaceScalarField, scalar>(phi);

    tmp<scalarField> Unt
    (
        new scalarField(phip/patch().magSf())
    );
    scalarField& Un = Unt.ref();

    // for moving mesh cases subtract velocity due to meshPhi
    if (patch().boundaryMesh().mesh().moving())
    {
        Un -= fvc::meshPhi
        (
            db().lookupObject<volVectorField>("U")
        )->boundaryField()[patch().index()]/patch().magSf();
    }
    // for patch inside MRF, get relative velocity by subtracting
    // the frame velocity
    if
    (
        coorFramePtr() &&
        !this->db().foundObject<foamCoupledControl>(solutionControl::typeName)
    )
    {
        Un -= coorFramePtr()->frameVelocity(patch().Cf()) & patch().nf();
    }

    if (phi.dimensions() == dimMass/dimTime)
    {
        Un /= patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), rhoName_);
    }

    return Unt;
}


void Foam::porousBafflePressureAMIFvPatchField::computeUniformVelocity
(
    scalarField& Un
) const
{
    if (uniformJump_)
    {
        const scalarField& magSf = patch().magSf();
        scalar patchArea = gSum(magSf);

        Un = gSum(Un*magSf);

        if (patchArea>SMALL)
        {
            Un /= patchArea;
        }
    }
}


void Foam::porousBafflePressureAMIFvPatchField::correctTemperatureDependency
(
    scalarField& linerCoeff,
    scalarField& nonLinearCoeff
) const
{
    if (temperatureDependence_)
    {
        if (patch().boundaryMesh().mesh().foundObject<volScalarField>(TName_))
        {
            const volScalarField& T = patch().boundaryMesh().mesh().
                lookupObject<volScalarField>(TName_);

            const fvPatchField<scalar>& Tp =
                patch().patchField<volScalarField, scalar>(T);

            forAll(Tp, fI)
            {
                linerCoeff[fI] *= calcDcorrValue(Tp[fI]);
                nonLinearCoeff[fI] *= calcFcorrValue(Tp[fI]);
            }
        }
    }
}


void Foam::porousBafflePressureAMIFvPatchField::computeDFJump()
{
    const scalar t = db().time().timeOutputValue();
    scalarField d(this->size(), d_->value(t));
    scalarField f(this->size(), f_->value(t));
    correctTemperatureDependency(d, f);

    tmp<scalarField> Unt(computeUn());
    scalarField Un = Unt.ref();

    //- store previous jump. It might be used later in the stabilization
    scalarField jump0 = jump_;

// for alternative mode,full U needed
    if (alternativeModel_)
    {
        vectorField Up =
            patch().lookupPatchFieldInDb<volVectorField, vector>(db(), "U");

        const fvMesh& mesh = patch().boundaryMesh().mesh();

        // for moving mesh cases subtract velocity due to meshPhi
        if (mesh.moving())
        {
            Up -= fvc::reconstruct
            (
                (fvc::meshPhi(db().lookupObject<volVectorField>("U")))
            )->boundaryField()[patch().index()];
        }
        // for patch inside MRF, get relative velocity by subtracting
        // the frame velocity
        if (coorFramePtr())
        {
            Up -= coorFramePtr()->frameVelocity(patch().Cf());
        }

        const scalar Ct1 = Ct1_->value(t);
        const scalar Ct2 = Ct2_->value(t);

        tmp<vectorField> Ut(Up - (Up & patch().nf())*patch().nf());

        if (fluxNormalVelocity_)
        {
            Up -= ((Up & patch().nf()) - Un)*patch().nf();
        }

        if (uniformJump_)
        {
            const scalarField& magSf = patch().magSf();
            scalar patchArea = gSum(magSf);

            Up = gSum(Up*magSf);
            Ut.ref() = gSum(Ut()*magSf);
            Un = gSum(Un*magSf);

            if (patchArea>SMALL)
            {
                Up /= patchArea;
                Ut.ref() /= patchArea;
                Un /= patchArea;
            }
        }

        scalarField magU(mag(Up));
        scalarField magUt(mag(Ut));

        if
        (
            !db().foundObject<turbulenceModel>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    internalField().group()
                )
            )
        )
        {
            jump_ = -sign(Un)*(f*0.5*magU)*magU*length_;
        }
        else
        {
            const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    internalField().group()
                )
            );

            jump_ =
                -sign(Un)
                *max(
                    d*turbModel.nu(patch().index())
                  + f*0.5*magU,
                    scalar(0)
                 )*magU*length_;
        }

        jump_ *= (1 + Ct1*magUt/magU + Ct2*sqr(magUt/magU));
    }
    else
    {
        computeUniformVelocity(Un);

        scalarField magUn(mag(Un));

        if
        (
            !db().foundObject<turbulenceModel>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    internalField().group()
                )
            )
        )
        {
            jump_ = -sign(Un)*(f*0.5*magUn)*magUn*length_;
        }
        else
        {
            const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
            (
                IOobject::groupName
                (
                    turbulenceModel::propertiesName,
                    internalField().group()
                )
            );

            jump_ =
                -sign(Un)
                *max(
                    d*turbModel.nu(patch().index())
                  + f*0.5*magUn,
                    scalar(0)
                 )*magUn*length_;
        }
    }

    if (internalField().dimensions() == dimPressure)
    {
        jump_ *= patch().lookupPatchFieldInDb<volScalarField, scalar>(db(), rhoName_);
    }

    stabiliseJump(jump0);

    printDebugInfo(Un);
}


void Foam::porousBafflePressureAMIFvPatchField::computeABJump()
{
    const scalar t = db().time().timeOutputValue();
    scalarField alpha(this->size(), alpha_->value(t));
    scalarField beta(this->size(), beta_->value(t));
    correctTemperatureDependency(alpha, beta);

    tmp<scalarField> Unt(computeUn());
    scalarField Un = Unt.ref();

    //- store previous jump. It might be used later in the stabilization
    scalarField jump0 = jump_;

    computeUniformVelocity(Un);
    scalarField magUn(mag(Un));

    jump_ = -sign(Un)*
    (
        alpha + beta*magUn
    )*magUn*length_;

    if (internalField().dimensions() == dimPressure)
    {
        jump_ *= patch().lookupPatchFieldInDb<volScalarField, scalar>
        (
            db(),
            rhoName_
        );
    }

    stabiliseJump(jump0);

    printDebugInfo(Un);
}


void Foam::porousBafflePressureAMIFvPatchField::computePLJump()
{
    const scalar t = db().time().timeOutputValue();
    const scalar c0 = C0_->value(t);
    const scalar c1 = C1_->value(t);
    const scalar C1m1 = c1 - 1.0;

    tmp<scalarField> Unt(computeUn());
    scalarField Un = Unt.ref();

    //- store previous jump. It might be used later in the stabilization
    scalarField jump0 = jump_;

    computeUniformVelocity(Un);

    scalarField magUn(mag(Un));

    jump_ = -sign(Un)*c0*pow(magSqr(Un), C1m1)*magUn*length_;

    if (internalField().dimensions() == dimPressure)
    {
        jump_ *= patch().lookupPatchFieldInDb<volScalarField, scalar>
        (
            db(),
            rhoName_
        );
    }

    stabiliseJump(jump0);

    printDebugInfo(Un);
}


void Foam::porousBafflePressureAMIFvPatchField::printDebugInfo
(
    const scalarField& Un
) const
{
    if (debug)
    {
        scalar avePressureJump = gAverage(jump_);
        scalar aveVelocity = gAverage(Un);

        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << " Average pressure drop :" << avePressureJump
            << " Average velocity :" << aveVelocity
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousBafflePressureAMIFvPatchField::porousBafflePressureAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF),
    referenceFrameFvPatch<vector>(p),
    porousMode_(darcyForchheimer),
    phiName_("phi"),
    rhoName_("rho"),
    d_(),
    f_(),
    Ct1_(),
    Ct2_(),
    alpha_(),
    beta_(),
    C0_(),
    C1_(),
    length_(0),
    relax_(1.0),
    minJump_(-GREAT),
    maxJump_(GREAT),
    uniformJump_(false),
    alternativeModel_(false),
    fluxNormalVelocity_(true),
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


Foam::porousBafflePressureAMIFvPatchField::porousBafflePressureAMIFvPatchField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpAMIFvPatchField<scalar>(p, iF),
    referenceFrameFvPatch<vector>(dict, p),
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
    )
    ,
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    d_(),
    f_(),
    Ct1_(),
    Ct2_(),
    alpha_(),
    beta_(),
    C0_(),
    C1_(),
    length_(readScalar(dict.lookup("length"))),
    relax_(dict.lookupOrDefault<scalar>("relax", 1.0)),
    minJump_(dict.lookupOrDefault<scalar>("minJump", -GREAT)),
    maxJump_(dict.lookupOrDefault<scalar>("maxJump", GREAT)),
    uniformJump_(dict.lookupOrDefault<Switch>("uniformJump", false)),
    alternativeModel_(dict.lookupOrDefault<Switch>("alternativeModel", false)),
    fluxNormalVelocity_
    (
        dict.lookupOrDefault<Switch>("fluxNormalVelocity", true)
    ),
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
    switch (porousMode_)
    {
        case darcyForchheimer:
        {
            d_ = Function1<scalar>::New("d", dict);
            f_ = Function1<scalar>::New("f", dict);
            if (alternativeModel_)
            {
                Ct1_ = Function1<scalar>::New("Ct1", dict);
                Ct2_ = Function1<scalar>::New("Ct2", dict);
            }

            break;
        }
        case alphaBeta:
        {
            alpha_ = Function1<scalar>::New("alpha", dict);
            beta_ = Function1<scalar>::New("beta", dict);

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
    if (dict.found("jump"))
    {
        jump_ = Field<scalar>("jump", dict, p.size());
    }
    fvPatchField<scalar>::operator=
    (
        Field<scalar>("value", dict, p.size())
    );
}


Foam::porousBafflePressureAMIFvPatchField::porousBafflePressureAMIFvPatchField
(
    const porousBafflePressureAMIFvPatchField& ptf,
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
    porousMode_(ptf.porousMode_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    d_(ptf.d_, false),
    f_(ptf.f_, false),
    Ct1_(ptf.Ct1_, false),
    Ct2_(ptf.Ct2_, false),
    alpha_(ptf.alpha_, false),
    beta_(ptf.beta_, false),
    C0_(ptf.C0_, false),
    C1_(ptf.C1_, false),
    length_(ptf.length_),
    relax_(ptf.relax_),
    minJump_(ptf.minJump_),
    maxJump_(ptf.maxJump_),
    uniformJump_(ptf.uniformJump_),
    alternativeModel_(ptf.alternativeModel_),
    fluxNormalVelocity_(ptf.fluxNormalVelocity_),
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


Foam::porousBafflePressureAMIFvPatchField::porousBafflePressureAMIFvPatchField
(
    const porousBafflePressureAMIFvPatchField& ptf
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf),
    referenceFrameFvPatch<vector>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_
    ),
    porousMode_(ptf.porousMode_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    d_(ptf.d_, false),
    f_(ptf.f_, false),
    Ct1_(ptf.Ct1_, false),
    Ct2_(ptf.Ct2_, false),
    alpha_(ptf.alpha_, false),
    beta_(ptf.beta_, false),
    C0_(ptf.C0_, false),
    C1_(ptf.C1_, false),
    length_(ptf.length_),
    relax_(ptf.relax_),
    minJump_(ptf.minJump_),
    maxJump_(ptf.maxJump_),
    uniformJump_(ptf.uniformJump_),
    alternativeModel_(ptf.alternativeModel_),
    fluxNormalVelocity_(ptf.fluxNormalVelocity_),
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


Foam::porousBafflePressureAMIFvPatchField::porousBafflePressureAMIFvPatchField
(
    const porousBafflePressureAMIFvPatchField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedJumpAMIFvPatchField<scalar>(ptf, iF),
    referenceFrameFvPatch<vector>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_
    ),
    porousMode_(ptf.porousMode_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    d_(ptf.d_, false),
    f_(ptf.f_, false),
    Ct1_(ptf.Ct1_, false),
    Ct2_(ptf.Ct2_, false),
    alpha_(ptf.alpha_, false),
    beta_(ptf.beta_, false),
    C0_(ptf.C0_, false),
    C1_(ptf.C1_, false),
    length_(ptf.length_),
    relax_(ptf.relax_),
    minJump_(ptf.minJump_),
    maxJump_(ptf.maxJump_),
    uniformJump_(ptf.uniformJump_),
    alternativeModel_(ptf.alternativeModel_),
    fluxNormalVelocity_(ptf.fluxNormalVelocity_),
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::porousBafflePressureAMIFvPatchField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedJumpAMIFvPatchField<scalar>::autoMap(m);
    referenceFrameFvPatch<vector>::autoMap(m);
}


void Foam::porousBafflePressureAMIFvPatchField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedJumpAMIFvPatchField<scalar>::rmap(ptf, addr);
    const porousBafflePressureAMIFvPatchField& pbptf =
        refCast<const porousBafflePressureAMIFvPatchField>(ptf);
    referenceFrameFvPatch<vector>::rmap(pbptf, addr);
}


void Foam::porousBafflePressureAMIFvPatchField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fixedJumpAMIFvPatchField<scalar>::autoMapGIB(mapper);
    referenceFrameFvPatch<vector>::autoMapGIB(mapper);
}


void Foam::porousBafflePressureAMIFvPatchField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!db().foundObject<surfaceScalarField>(phiName_))
    {
        return;
    }

    switch (porousMode_)
    {
        case darcyForchheimer:
        {
            computeDFJump();

            break;
        }
        case alphaBeta:
        {
            computeABJump();

            break;
        }
        case powerLaw:
        {
            computePLJump();

            break;
        }
    }

    fixedJumpAMIFvPatchField<scalar>::updateCoeffs();
}


void Foam::porousBafflePressureAMIFvPatchField::write(Ostream& os) const
{
    fixedJumpAMIFvPatchField<scalar>::write(os);
    if (!this->cyclicAMIPatch().owner())
    {
        jump_.writeEntry("jump", os);
    }
    referenceFrameFvPatch::write(os);

    os.writeEntry("porousMode", porousModeNames_[porousMode_]);

    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    switch (porousMode_)
    {
        case darcyForchheimer:
        {
            d_->writeData(os);
            f_->writeData(os);
            if (alternativeModel_)
            {
                Ct1_->writeData(os);
                Ct2_->writeData(os);
            }

            os.writeEntry("alternativeModel", alternativeModel_);

            if (alternativeModel_)
            {
                os.writeEntry("fluxNormalVelocity", fluxNormalVelocity_);
            }
            break;
        }
        case alphaBeta:
        {
            alpha_->writeData(os);
            beta_->writeData(os);

            break;
        }
        case powerLaw:
        {
            C0_->writeData(os);
            C1_->writeData(os);

            break;
        }
    }
    writeEntryIfDifferent<Switch>
    (
        os, "temperatureDependence", false, temperatureDependence_
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
    os.writeEntry("length", length_);
    writeEntryIfDifferent<scalar>(os, "relax", 1.0, relax_);
    writeEntryIfDifferent<scalar>(os, "minJump", -GREAT, minJump_);
    writeEntryIfDifferent<scalar>(os, "maxJump", GREAT, maxJump_);
    os.writeEntry("uniformJump", uniformJump_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        porousBafflePressureAMIFvPatchField
    );
}

// ************************************************************************* //
