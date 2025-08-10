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
    (c) 2010-2019 Esi Ltd.
    (c) 2012-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvPatchFields/coupledElectricPotential/coupledElectricPotentialFvPatchScalarField.H"
#include "db/Time/Time.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "solidThermo/solidThermo.H"
#include "derivedFvPatchFields/phaseChangeHumidity/phaseChangeHumidityFvPatchScalarField.H"
#include "fields/volFields/volFields.H"
#include "fvPatchFields/regionCoupledFlux/regionCoupledFluxFvPatchScalarField.H"


// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchField<scalar>(p, iF),
    electricalBoundaryBase(p),
    lWeights_(p.weights()),
    contactResistance_(0),
    contactInsulance_(0)
{}


coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCoupledFvPatchField<scalar>(ptf, p, iF, mapper),
    electricalBoundaryBase(p),
    lWeights_(ptf.lWeights_),
    contactResistance_(ptf.contactResistance_),
    sigmaLayers_(ptf.sigmaLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_)
{}


coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCoupledFvPatchField<scalar>(p, iF, dict),
    electricalBoundaryBase(p, dict),
    lWeights_(p.weights()),
    contactResistance_(scalar(0)),
    contactInsulance_(scalar(0))
{
    DeprecationWarningInFunction
    (
        this->typeName, "boundary condition", 40000
    )
        << "Please use the "
        << regionCoupledFluxFvPatchScalarField::typeName
        << " boundary instead." << nl << endl;

    if (dict.found("contactResistance") && dict.found("sigmaLayers"))
    {
        FatalErrorInFunction
            << "Either contact resistance or contact layers should be "
            << "specified in boundary "
            << p.name() << " of field " << iF.name() << endl;
    }

    if (dict.found("contactResistance"))
    {
        dict.lookup("contactResistance") >> contactResistance_;
        contactInsulance_ = contactResistance_*gSum(p.magSf());
    }
    else if (dict.found("sigmaLayers") || dict.found("thicknessLayers"))
    {
        dict.lookup("sigmaLayers") >> sigmaLayers_;
        dict.lookup("thicknessLayers") >> thicknessLayers_;

        if (thicknessLayers_.size() != sigmaLayers_.size())
        {
            FatalIOErrorInFunction(dict)
                << "sigmaLayers and thicknessLayers must be lists of the same "
                << "length" << nl << endl;
        }

        if (sigmaLayers_.size() > 0)
        {
            // Calculate effective electrical contact insulance by harmonic averaging
            forAll(sigmaLayers_, iLayer)
            {
                contactInsulance_ += thicknessLayers_[iLayer] / sigmaLayers_[iLayer];
            }
        }
    }
}


coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& ptf
)
:
    regionCoupledFvPatchField<scalar>(ptf),
    electricalBoundaryBase(ptf.patch()),
    lWeights_(ptf.lWeights_),
    contactResistance_(ptf.contactResistance_),
    sigmaLayers_(ptf.sigmaLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_)
{}


coupledElectricPotentialFvPatchScalarField::
coupledElectricPotentialFvPatchScalarField
(
    const coupledElectricPotentialFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchField<scalar>(ptf, iF),
    electricalBoundaryBase(ptf.patch()),
    lWeights_(ptf.lWeights_),
    contactResistance_(ptf.contactResistance_),
    sigmaLayers_(ptf.sigmaLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    contactInsulance_(ptf.contactInsulance_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledElectricPotentialFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    regionCoupledFvPatchField<scalar>::autoMap(m);
}


void coupledElectricPotentialFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    regionCoupledFvPatchField<scalar>::rmap(ptf, addr);
}


tmp<scalarField> coupledElectricPotentialFvPatchScalarField::snGrad() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::snGrad();
    }

    return patch().deltaCoeffs() * (*this - patchInternalField());
}


tmp<scalarField>
coupledElectricPotentialFvPatchScalarField::snGrad(const scalarField&) const
{
    return snGrad();
}


void coupledElectricPotentialFvPatchScalarField::evaluate
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

    scalarField::operator=
    (
        lWeights_ * patchInternalField()
        + (1.0-lWeights_) * patchNeighbourField()
        + fCorr_
    );

    fvPatchScalarField::evaluate();
}


void coupledElectricPotentialFvPatchScalarField::updateCoeffs()
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

    const fvPatch& patch = this->patch();
    const label patchi = patch.index();
    const scalarField& deltaCoeffs = patch.deltaCoeffs();
    const scalarField deltaBySigma
    (
        scalar(1)/stabilise((nfSigma() & patch.nf())*deltaCoeffs, SMALL)
    );

    const fvPatch& nbrPatch = regionCoupledPatch_.nbrPatch();
    label nbrPatchi = nbrPatch.index();
    const coupledElectricPotentialFvPatchScalarField& nbrPatchField =
        refCast<const coupledElectricPotentialFvPatchScalarField>
        (
            neighbourFvPatchField()
        );
    const fvMesh& nbrMesh = nbrPatch.boundaryMesh().mesh();
    const scalarField& nbrDeltaCoeffs =
        nbrMesh.deltaCoeffs().boundaryField()[nbrPatch.index()];
    scalarField nbrDeltaBySigma
    (
        scalar(1)
      / stabilise
        (
            (nbrPatchField.nfSigma() & nbrPatch.nf())*nbrDeltaCoeffs,
            SMALL
        )
    );
    // Interpolate across region boundary. Use this-side value as default in
    // case of missing mapping (low weight) in AMI
    scalarField nbrDeltaBySigmaInterp( interpolateFromNeighbour(nbrDeltaBySigma, deltaBySigma) );

    // Set weights
    lWeights_ =
        scalar(1)
      - (
            deltaBySigma
           /(
               deltaBySigma+nbrDeltaBySigmaInterp
             + contactInsulance_+nbrPatchField.contactInsulance_
            )
        );

    // Explicit correction to the face interpolated value

    fCorr_ = scalarField(patch.size(), 0);

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

    if (anisotropic() && mesh.foundObject<volVectorField>(gradName))
    {
        volVectorField& grad =
            mesh.lookupObjectRef<volVectorField>(gradName);

        scalarField corr
        (
            deltaBySigma
           *(nfSigma() & grad.boundaryField()[patchi])
          - (*this - this->patchInternalField())
        );
        fCorr_ -= lWeights_*corr;
    }
    if
    (
        nbrPatchField.anisotropic()
     && nbrMesh.foundObject<volVectorField>(nbrGradName)
    )
    {
        volVectorField& nbrGrad =
            nbrMesh.lookupObjectRef<volVectorField>(nbrGradName);

        scalarField nbrCorr
        (
            nbrDeltaBySigma
           *(
                nbrPatchField.nfSigma()
              & nbrGrad.boundaryField()[nbrPatchi]
            )
          - (nbrPatchField - nbrPatchField.patchInternalField())
        );
        fCorr_ -=
            (scalar(1)-lWeights_)
           *interpolateFromNeighbour(nbrCorr, scalarField(patch.size(), 0));
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    fvPatchScalarField::updateCoeffs();
}


tmp<scalarField>
coupledElectricPotentialFvPatchScalarField::valueInternalCoeffs
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
coupledElectricPotentialFvPatchScalarField::valueBoundaryCoeffs
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
coupledElectricPotentialFvPatchScalarField::valueCorrection() const
{
    return tmp<scalarField>(new scalarField(patch().size(), 0));
}


tmp<scalarField>
coupledElectricPotentialFvPatchScalarField::gradientInternalCoeffs() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::gradientInternalCoeffs();
    }

    return patch().deltaCoeffs()
      * (valueInternalCoeffs(patch().weights()) - 1.0);
}


tmp<scalarField>
coupledElectricPotentialFvPatchScalarField::gradientBoundaryCoeffs() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::gradientBoundaryCoeffs();
    }

    // Note - for interfaces, this is a coefficient of the boundary value,
    // not a source value
    return patch().deltaCoeffs() * valueBoundaryCoeffs(patch().weights());
}


tmp<scalarField>
coupledElectricPotentialFvPatchScalarField::gradientCorrection() const
{
    return tmp<scalarField>(new scalarField(patch().size(), 0));
}


void coupledElectricPotentialFvPatchScalarField::regionCoupledBoundaryCoeffs
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
    coupledBouCoeffs = (1.0-lWeights_) * bouCoeffs;
}


void coupledElectricPotentialFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    if (sigmaLayers_.size())
    {
        sigmaLayers_.writeEntry("sigmaLayers", os);
        thicknessLayers_.writeEntry("thicknessLayers", os);
    }
    else
    {
        os.writeEntry("contactResistance", contactResistance_);
    }

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


makePatchTypeField
(
    fvPatchScalarField,
    coupledElectricPotentialFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
