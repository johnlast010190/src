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
    (c) 2012-2013 OpenFOAM Foundation
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvPatchFields/coupledWallEnergy/coupledWallEnergyFvPatchScalarField.H"
#include "db/Time/Time.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "solidThermo/solidThermo.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "fields/volFields/volFields.H"
#include "fvPatchFields/regionCoupledFlux/regionCoupledEnergyFluxFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Private members  * * * * * * * * * * * * *//

void coupledWallEnergyFvPatchScalarField::setThermos() const
{
    if (!nbrThermoPtr_)
    {
        nbrThermoPtr_ =
        (
            &regionCoupledPatch_.nbrMesh().lookupObject<basicThermo>
            (
                IOobject::groupName
                (
                    basicThermo::dictName,
                    neighbourPhaseName_
                )
            )
        );
    }

    if (!thermoPtr_)
    {
        thermoPtr_ =
        (
            &this->db().lookupObject<basicThermo>
            (
                IOobject::groupName
                (
                    basicThermo::dictName,
                    internalField().group()
                )
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledWallEnergyFvPatchScalarField::
coupledWallEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchField<scalar>(p, iF),
    boundaryKappa
    (
        db(), patch(), "lookup", "undefined", "undefined", iF.group()
    ),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr),
    lWeights_(p.weights()),
    fCorrEval_(p.size(), 0.0),
    contactInsulance_(0),
    QrName_("Qr"),
    QrNbrName_("Qr"),
    Qsum_(p.size(), 0.0),
    qURF_(0.3)
{
    // The he field is constructed with this constructor, not read from
    // dictionary. Look up the dictionary settings from the corresponding
    // temperature field.

    DeprecationWarningInFunction
    (
        this->typeName, "boundary condition", 40000
    )
        << "Please use the "
        << regionCoupledEnergyFluxFvPatchScalarField::typeName
        << " boundary instead." << nl << endl;

    if (iF.member() != "T")
    {
        const basicThermo& thermo =
            iF.db().lookupObject<basicThermo>
            (
                iF.groupName(basicThermo::dictName, iF.group())
            );
        const fvPatchField& Tpf = thermo.T().boundaryField()[p.index()];
        const coupledWallEnergyFvPatchScalarField& Tepf =
            refCast<const coupledWallEnergyFvPatchScalarField>(Tpf);
        method_ = Tepf.method_;
        kappaName_ = Tepf.kappaName_;
        alphaAniName_ = Tepf.alphaAniName_;
        thicknessLayers_ = Tepf.thicknessLayers_;
        kappaLayers_ = Tepf.kappaLayers_;
        fCorrEval_ = Tepf.fCorrEval_;
        contactInsulance_ = Tepf.contactInsulance_;
        QrName_ = Tepf.QrName_;
        QrNbrName_ = Tepf.QrNbrName_;
        Qsum_ = Tepf.Qsum_;
        qURF_ = Tepf.qURF_;

        if (Tepf.neighbourFieldName_ != Tepf.internalField().name())
        {
            neighbourFieldName_ = Tepf.neighbourFieldName_;
        }
        else
        {
            neighbourFieldName_ = iF.name();
        }
        neighbourPhaseName_ = Tepf.neighbourPhaseName_;
    }

}


coupledWallEnergyFvPatchScalarField::
coupledWallEnergyFvPatchScalarField
(
    const coupledWallEnergyFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCoupledFvPatchField<scalar>(ptf, p, iF, mapper),
    boundaryKappa(db(), patch(), ptf),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr),
    lWeights_(mapper(ptf.lWeights_)),
    fCorr_(mapper(ptf.fCorr_)),
    fCorrEval_(mapper(ptf.fCorrEval_)),
    pCpEffByNbrCpEff_(mapper(ptf.pCpEffByNbrCpEff_)),
    kappaLayers_(ptf.kappaLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_),
    QrName_(ptf.QrName_),
    QrNbrName_(ptf.QrNbrName_),
    Qsum_(mapper(ptf.Qsum_)),
    qURF_(ptf.qURF_)
{}


coupledWallEnergyFvPatchScalarField::
coupledWallEnergyFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCoupledFvPatchField<scalar>(p, iF, dict),
    boundaryKappa(db(), patch(), dict, iF.group()),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr),
    lWeights_(p.weights()),
    fCorrEval_
    (
        dict.found("fCorr")
      ? scalarField(word("fCorr"), dict, p.size()).xfer()
      : scalarField(p.size(), 0.0).xfer()
    ),
    contactInsulance_(0.0),
    QrName_(dict.lookupOrDefault("Qr", word("Qr"))),
    QrNbrName_(dict.lookupOrDefault("QrNbr", word("Qr"))),
    Qsum_
    (
        dict.found("Qsum")
      ? scalarField(word("Qsum"), dict, p.size()).xfer()
      : scalarField(p.size(), 0.0).xfer()
    ),
    qURF_(dict.lookupOrDefault<scalar>("Qurf", 0.3))
{
    if (dict.found("kappaLayers") || dict.found("thicknessLayers"))
    {
        dict.lookup("kappaLayers") >> kappaLayers_;
        dict.lookup("thicknessLayers") >> thicknessLayers_;

        if (thicknessLayers_.size() != kappaLayers_.size())
        {
            FatalIOErrorInFunction(dict)
                << "kappaLayers and thicknessLayers must be lists of the same "
                << "length" << nl << endl;
        }

        if (kappaLayers_.size() > 0)
        {
            // Calculate effective thermal resistance by harmonic averaging
            forAll(kappaLayers_, iLayer)
            {
                contactInsulance_ += thicknessLayers_[iLayer] / kappaLayers_[iLayer];
            }
        }
    }
}


coupledWallEnergyFvPatchScalarField::
coupledWallEnergyFvPatchScalarField
(
    const coupledWallEnergyFvPatchScalarField& ptf
)
:
    regionCoupledFvPatchField<scalar>(ptf),
    boundaryKappa(db(), patch(), ptf),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr),
    lWeights_(ptf.lWeights_),
    fCorrEval_(ptf.fCorrEval_),
    kappaLayers_(ptf.kappaLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_),
    QrName_(ptf.QrName_),
    QrNbrName_(ptf.QrNbrName_),
    Qsum_(ptf.Qsum_),
    qURF_(ptf.qURF_)
{}


coupledWallEnergyFvPatchScalarField::
coupledWallEnergyFvPatchScalarField
(
    const coupledWallEnergyFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchField<scalar>(ptf, iF),
    boundaryKappa(db(), patch(), ptf),
    nbrThermoPtr_(nullptr),
    thermoPtr_(nullptr),
    lWeights_(ptf.lWeights_),
    fCorrEval_(ptf.fCorrEval_),
    kappaLayers_(ptf.kappaLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_),
    QrName_(ptf.QrName_),
    QrNbrName_(ptf.QrNbrName_),
    Qsum_(ptf.Qsum_),
    qURF_(ptf.qURF_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledWallEnergyFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    regionCoupledFvPatchField<scalar>::autoMap(m);
    m(fCorrEval_, fCorrEval_);
    m(Qsum_, Qsum_);
}


void coupledWallEnergyFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    regionCoupledFvPatchField<scalar>::rmap(ptf, addr);

    const coupledWallEnergyFvPatchScalarField& tiptf =
        refCast<const coupledWallEnergyFvPatchScalarField>(ptf);

    fCorr_.rmap(tiptf.fCorr_, addr);
    fCorrEval_.rmap(tiptf.fCorrEval_, addr);
    pCpEffByNbrCpEff_.rmap(tiptf.pCpEffByNbrCpEff_, addr);
    kappaLayers_ = tiptf.kappaLayers_;
    thicknessLayers_ = tiptf.thicknessLayers_;
    Qsum_.rmap(tiptf.Qsum_, addr);
}


void coupledWallEnergyFvPatchScalarField::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    regionCoupledFvPatchField<scalar>::autoMapGIB(mapper);
    mapper.map(fCorrEval_, scalar(0));
    mapper.map(Qsum_, scalar(0));
}


tmp<scalarField> coupledWallEnergyFvPatchScalarField::snGrad() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::snGrad();
    }

    return patch().deltaCoeffs() * (*this - patchInternalField());
}


tmp<scalarField>
coupledWallEnergyFvPatchScalarField::snGrad(const scalarField&) const
{
    return snGrad();
}


void coupledWallEnergyFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
    if (!regionCoupled())
    {
        regionCoupledFvPatchField<scalar>::evaluate();
        return;
    }

    updateCoeffs();
    if (ishe())
    {
        scalarField::operator=
        (
            lWeights_ * patchInternalField()
          + (1.0-lWeights_) * pCpEffByNbrCpEff_ * patchNeighbourField()
          + fCorrEval_
        );
    }
    else
    {
        scalarField::operator=
        (
            lWeights_ * patchInternalField()
          + (1.0-lWeights_) * patchNeighbourField()
          + fCorrEval_
        );
    }

    fvPatchScalarField::evaluate();
}


void coupledWallEnergyFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    else if (!regionCoupled())
    {
        fvPatchScalarField::updateCoeffs();
        return;
    }

    // Change the tag in case we are inside initEvaluate/evaluate and there are
    // processor comms underway.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag+1;

    setThermos();

    const fvPatch& patch = this->patch();
    const label patchi = patch.index();
    const scalarField& deltaCoeffs = patch.deltaCoeffs();
    const scalarField& pp = thermoPtr_->p().boundaryField()[patchi];
    const scalarField& pT0 = thermoPtr_->T().boundaryField()[patchi];
    scalarField pCpEff;
    if (ishe())
    {
        pCpEff = thermoPtr_->he(pp, pT0, patchi)/pT0;
    }
    else
    {
        pCpEff = thermoPtr_->Cp(pp, pT0, patchi);
    }
    const scalarField deltaByKappa
    (
        scalar(1)/stabilise((nfAlphaEff()*pCpEff & patch.nf())*deltaCoeffs, SMALL)
    );

    const fvPatch& nbrPatch = regionCoupledPatch_.nbrPatch();
    label nbrPatchi = nbrPatch.index();
    const coupledWallEnergyFvPatchScalarField& nbrPatchField =
        refCast<const coupledWallEnergyFvPatchScalarField>
        (
            neighbourFvPatchField()
        );
    const fvMesh& nbrMesh = nbrPatch.boundaryMesh().mesh();
    const scalarField& nbrDeltaCoeffs =
        nbrMesh.deltaCoeffs().boundaryField()[nbrPatch.index()];
    nbrPatchField.setThermos();
    const scalarField& pNbrp = nbrThermoPtr_->p().boundaryField()[nbrPatchi];
    const scalarField& pNbrT0 = nbrThermoPtr_->T().boundaryField()[nbrPatchi];
    scalarField pNbrCpEff;
    if (ishe())
    {
        pNbrCpEff = nbrThermoPtr_->he(pNbrp, pNbrT0, nbrPatchi)/pNbrT0;
    }
    else
    {
        pNbrCpEff = nbrThermoPtr_->Cp(pNbrp, pNbrT0, nbrPatchi);
    }
    scalarField nbrDeltaByKappa
    (
        scalar(1)
      / stabilise
        (
            (nbrPatchField.nfAlphaEff()*pNbrCpEff & nbrPatch.nf())*nbrDeltaCoeffs,
            SMALL
        )
    );

    // Interpolate across region boundary. Use this-side value as default in
    // case of missing mapping (low weight) in AMI
    scalarField nbrDeltaByKappaInterp( interpolateFromNeighbour(nbrDeltaByKappa, deltaByKappa) );

    scalarField pNbrCpEffInterp( interpolateFromNeighbour(pNbrCpEff, pCpEff) );
    pCpEffByNbrCpEff_ = pCpEff/pNbrCpEffInterp;

    // Set weights
    lWeights_ =
        scalar(1)
      - (
            deltaByKappa
           /(
               deltaByKappa+nbrDeltaByKappaInterp
             + contactInsulance_+nbrPatchField.contactInsulance_
            )
        );

    // Explicit correction to the face interpolated value
    fCorr_ = scalarField(this->size(), 0);

    // Use the saved gradient to (1) ensure consistency with that used in
    // Laplacian operator, (2) avoid unnecessary recalculation of the whole
    // field and (3) prevent synchronisation issues with consolidated regions
    // if we were to calculate our neighbour region's gradient.

    const fvMesh& mesh(patch.boundaryMesh().mesh());
    const volScalarField& field =
        mesh.lookupObject<volScalarField>(internalField().name());
    word gradName = "grad(" + field.name() + ")";

    const volScalarField& nbrField =
        nbrMesh.lookupObject<volScalarField>(coupledField());
    word nbrGradName = "grad(" + nbrField.name() + ")";

    if
    (
        isA<solidThermo>(*thermoPtr_)
     && !refCast<const solidThermo>(*thermoPtr_).isotropic()
    )
    {
        autoPtr<volVectorField> pgrad;
        if (!mesh.foundObject<volVectorField>(gradName))
        {
            serialThreads::pauseSwitching();
            pgrad.set(fvc::grad(field).ptr());
            pgrad().rename(gradName);
            serialThreads::resumeSwitching();
        }

        const volVectorField& grad =
            mesh.lookupObject<volVectorField>(gradName);

        scalarField corr
        (
            deltaByKappa
           *(nfAlphaEff()*pCpEff & grad.boundaryField()[patchi])
          - (*this - this->patchInternalField())
        );
        if (ishe())
        {
            corr /= pCpEff;
        }
        fCorr_ -= lWeights_*corr;
    }
    if
    (
        isA<solidThermo>(*nbrThermoPtr_)
     && !refCast<const solidThermo>(*nbrThermoPtr_).isotropic()
    )
    {
        autoPtr<volVectorField> pnbrGrad;
        if (!nbrMesh.foundObject<volVectorField>(nbrGradName))
        {
            serialThreads::pauseSwitching();
            pnbrGrad.set
            (
                fvc::grad(nbrField).ptr()
            );
            pnbrGrad().rename(nbrGradName);
            serialThreads::resumeSwitching();
        }

        const volVectorField& nbrGrad =
            nbrMesh.lookupObject<volVectorField>(nbrGradName);

        scalarField nbrCorr
        (
            nbrDeltaByKappa
           *(
                nbrPatchField.nfAlphaEff()*pNbrCpEff
              & nbrGrad.boundaryField()[nbrPatchi]
            )
          - (nbrPatchField - nbrPatchField.patchInternalField())
        );
        if (ishe())
        {
            nbrCorr /= pNbrCpEff;
        }
        fCorr_ -=
            (scalar(1)-lWeights_)
           *interpolateFromNeighbour(nbrCorr, scalarField(patch.size(), 0.0));
    }
    if (ishe())
    {
        fCorr_ *= pCpEff;
    }

    // Update heat sources, using Newton's method to account for possible
    // temperature dependence
    // We are solving
    // he_f = lWeight*he + (1-lWeight)*heNbr*pCpEffByNbrCpEff
    //        + fCorr + sourceCorr(he_f)

    scalarField C
    (
        lWeights_*patchInternalField()
      + (1-lWeights_)*patchNeighbourField()*pCpEffByNbrCpEff_
      + fCorr_
    );

    scalarField corrCoeff
    (
        (scalar(1)-lWeights_)
       *(nbrDeltaByKappaInterp+nbrPatchField.contactInsulance_)
    );
    if (ishe())
    {
        corrCoeff *= pCpEff;
    }

    label i = 0;
    scalar phe_err = GREAT;
    scalar maxErr = 1e-5;
    label maxLoops = 100;

    scalarField phe_new = *this;
    scalarField phe_old;

    scalarField Qopt;
    do
    {
        phe_old = phe_new;
        scalarField bSourceDeriv, bSourceDerivNbr;
        Qopt = boundarySources(phe_old, bSourceDeriv);
        Qopt +=
            interpolateFromNeighbour
            (
                neighbourFvPatchField().boundarySources
                (
                    interpolateToNeighbour
                    (
                        phe_old/pCpEffByNbrCpEff_,
                        neighbourFvPatchField()
                    ),
                    bSourceDerivNbr
                ),
                scalarField(this->size(), scalar(0))
            );
        bSourceDeriv +=
            interpolateFromNeighbour
            (
                bSourceDerivNbr,
                scalarField(this->size(), scalar(0))
            )
           *pCpEffByNbrCpEff_;

        scalarField sourceCorr( Qopt*corrCoeff );
        scalarField sourceCorrDeriv( bSourceDeriv*corrCoeff );

        scalarField f( C + sourceCorr - phe_old );
        scalarField df( sourceCorrDeriv - 1 );

        phe_new = phe_old - f/df;

        phe_err = gMax(mag(phe_new-phe_old)/phe_old);
        i++;

    }
    while (i < maxLoops && phe_err > maxErr);

    if (i == maxLoops)
    {
        WarningInFunction
            << "Non-convergence in Newton's method for "
            << "temperature-dependent boundary source, patch: "
            << this->patch().name() << nl << endl;
    }

    // Apply under-relaxation to Qsum
    Qsum_ = (1-qURF_)*Qsum_ + qURF_*Qopt;

    // Compute final fCorr_
    fCorr_ += Qsum_*corrCoeff;

    // Restore tag
    UPstream::msgType() = oldTag;

    fvPatchScalarField::updateCoeffs();
}


tmp<scalarField>
coupledWallEnergyFvPatchScalarField::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::valueInternalCoeffs(w);
    }

    return tmp<scalarField>(new scalarField(patch().size(), 0));
}


tmp<scalarField>
coupledWallEnergyFvPatchScalarField::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::valueBoundaryCoeffs(w);
    }

    // Note - for interfaces, this is a coefficient of the boundary value, not
    // the value itself
    return tmp<scalarField>(new scalarField(patch().size(), 1));
}


tmp<scalarField>
coupledWallEnergyFvPatchScalarField::gradientInternalCoeffs() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::gradientInternalCoeffs();
    }

    return patch().deltaCoeffs()
      * (valueInternalCoeffs(patch().weights()) - 1.0);
}


tmp<scalarField>
coupledWallEnergyFvPatchScalarField::gradientBoundaryCoeffs() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::gradientBoundaryCoeffs();
    }

    // Note - for interfaces, this is a coefficient of the boundary value,
    // not a source value
    return patch().deltaCoeffs() * valueBoundaryCoeffs(patch().weights());
}


void coupledWallEnergyFvPatchScalarField::regionCoupledBoundaryCoeffs
(
    const fvScalarMatrix& matrix,
    const scalarField& bouCoeffs,
    const scalarField& intCoeffs,
    scalarField& coupledBouCoeffs,
    scalarField& coupledIntCoeffs
)
{
    if (!regionCoupled())
    {
        return;
    }

    // Existing boundary coeffs are the (negative of) the coefficient of the
    // boundary value. Multiply these by the coefficients of the boundary
    // values as a function of internal and neighbour values
    coupledIntCoeffs = intCoeffs - lWeights_*bouCoeffs;
    if (ishe())
    {
        coupledBouCoeffs = (1.0-lWeights_) * pCpEffByNbrCpEff_ * bouCoeffs;
    }
    else
    {
        coupledBouCoeffs = (1.0-lWeights_) * bouCoeffs;
    }

    fCorrEval_ = fCorr_;
}


tmp<scalarField> coupledWallEnergyFvPatchScalarField::boundarySources
(
    const scalarField& pf, scalarField& df
) const
{
    // Forward to the companion temperature boundary
    if (ishe())
    {
        const label patchi = this->patch().index();
        const fvPatchScalarField& pT0 =
            thermoPtr_->T().boundaryField()[patchi];
        const fvPatchScalarField& pp =
            thermoPtr_->p().boundaryField()[patchi];
        scalarField pT( thermoPtr_->THE(pf, pp, pT0, patchi) );

        tmp<scalarField> f = pT0.boundarySources(pT, df);
        // Return derivative with respect to patch field value,
        // i.e. he not T
        scalarField dTdhe( 1.0/thermoPtr_->Cpv(pp, pT, patchi) );
        df *= dTdhe;
        return f;
    }
    else
    {
        return regionCoupledFvPatchField<scalar>::boundarySources(pf, df);
    }
}


void coupledWallEnergyFvPatchScalarField::write(Ostream& os) const
{
    regionCoupledFvPatchField<scalar>::write(os);
    boundaryKappa::write(os);
    if (QrName_ != "Qr")
    {
        os.writeEntry("Qr", QrName_);
    }
    if (QrNbrName_ != "Qr")
    {
        os.writeEntry("QrNbr", QrNbrName_);
    }
    if (kappaLayers_.size())
    {
        kappaLayers_.writeEntry("kappaLayers", os);
        thicknessLayers_.writeEntry("thicknessLayers", os);
    }

    //- Recover Qsum and fCorr
    scalarField QsumOut(Qsum_);
    scalarField fCorrOut(fCorrEval_);
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    if
    (
        internalField().member() == "T"
     && mesh.foundObject<basicThermo>(basicThermo::dictName)
    )
    {
        const basicThermo& thermo =
            mesh.lookupObject<basicThermo>
            (
                basicThermo::dictName
            );

        const fvPatchField& hep = thermo.he().boundaryField()[patch().index()];
        const coupledWallEnergyFvPatchScalarField& cwep =
            refCast<const coupledWallEnergyFvPatchScalarField>(hep);

        QsumOut = cwep.Qsum_;
        fCorrOut = cwep.fCorrEval_;
    }


    QsumOut.writeEntry("Qsum", os);
    fCorrOut.writeEntry("fCorr", os);

    writeEntryIfDifferent<scalar>
    (
        os,
        "Qurf",
        0.3,
        qURF_
    );

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


makePatchTypeField
(
    fvPatchScalarField,
    coupledWallEnergyFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
