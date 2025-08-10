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
    (c) 2010-2020 Esi Ltd.
    (c) 2012-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvPatchFields/coupledDiffusion/coupledDiffusionFvPatchScalarField.H"
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

coupledDiffusionFvPatchScalarField::
coupledDiffusionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchField<scalar>(p, iF),
    lWeights_(p.weights()),
    fCorr_(p.size(), scalar(0)),
    diffusivityName_("rhoD_"+iF.name()),
    rContactDiffusivity_(0)
{}


coupledDiffusionFvPatchScalarField::
coupledDiffusionFvPatchScalarField
(
    const coupledDiffusionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    regionCoupledFvPatchField<scalar>(ptf, p, iF, mapper),
    lWeights_(ptf.lWeights_),
    fCorr_(ptf.fCorr_),
    diffusivityName_(ptf.diffusivityName_),
    diffusivityLayers_(ptf.diffusivityLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    rContactDiffusivity_(ptf.rContactDiffusivity_)
{}


coupledDiffusionFvPatchScalarField::
coupledDiffusionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    regionCoupledFvPatchField<scalar>(p, iF, dict),
    lWeights_(p.weights()),
    fCorr_(p.size(), scalar(0)),
    diffusivityName_
    (
        dict.lookupOrDefault("diffusivityName", word("rhoD_"+iF.name()))
    ),
    rContactDiffusivity_(scalar(0))
{
    DeprecationWarningInFunction
    (
        this->typeName, "boundary condition", 40000
    )
        << "Please use the "
        << regionCoupledFluxFvPatchScalarField::typeName
        << " boundary instead." << nl << endl;

    if (dict.found("diffusivityLayers") || dict.found("thicknessLayers"))
    {
        dict.lookup("diffusivityLayers") >> diffusivityLayers_;
        dict.lookup("thicknessLayers") >> thicknessLayers_;

        if (thicknessLayers_.size() != diffusivityLayers_.size())
        {
            FatalIOErrorInFunction(dict)
                << "diffusivityLayers and thicknessLayers must be lists of the "
                << "same length" << nl << endl;
        }

        if (diffusivityLayers_.size() > 0)
        {
            // Calculate effective contact insulance by harmonic averaging
            forAll(diffusivityLayers_, iLayer)
            {
                rContactDiffusivity_ += thicknessLayers_[iLayer] / diffusivityLayers_[iLayer];
            }
        }
    }
    forceAssign(patchInternalField());
}


coupledDiffusionFvPatchScalarField::
coupledDiffusionFvPatchScalarField
(
    const coupledDiffusionFvPatchScalarField& ptf
)
:
    regionCoupledFvPatchField<scalar>(ptf),
    lWeights_(ptf.lWeights_),
    fCorr_(ptf.fCorr_),
    diffusivityName_(ptf.diffusivityName_),
    diffusivityLayers_(ptf.diffusivityLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    rContactDiffusivity_(ptf.rContactDiffusivity_)
{}


coupledDiffusionFvPatchScalarField::
coupledDiffusionFvPatchScalarField
(
    const coupledDiffusionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    regionCoupledFvPatchField<scalar>(ptf, iF),
    lWeights_(ptf.lWeights_),
    fCorr_(ptf.fCorr_),
    diffusivityName_(ptf.diffusivityName_),
    diffusivityLayers_(ptf.diffusivityLayers_),
    thicknessLayers_(ptf.thicknessLayers_),
    rContactDiffusivity_(ptf.rContactDiffusivity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void coupledDiffusionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    regionCoupledFvPatchField<scalar>::autoMap(m);
}


void coupledDiffusionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    regionCoupledFvPatchField<scalar>::rmap(ptf, addr);
}


bool Foam::coupledDiffusionFvPatchScalarField::anisotropic() const
{
    return
        patch().boundaryMesh().mesh().foundObject<volSymmTensorField>
        (
            diffusivityName_
        );
}


Foam::tmp<Foam::vectorField>
Foam::coupledDiffusionFvPatchScalarField::nfDiffusivity() const
{
    if (anisotropic())
    {
        const symmTensorField& diffusivityWall =
            patch().lookupPatchField<volSymmTensorField, scalar>
            (
                diffusivityName_
            );

        return patch().nf() & diffusivityWall;
    }
    else
    {
        return
            patch().nf()*patch().lookupPatchField<volScalarField, scalar>
            (
                diffusivityName_
            );
    }
}


tmp<scalarField> coupledDiffusionFvPatchScalarField::snGrad() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::snGrad();
    }

    return patch().deltaCoeffs() * (*this - patchInternalField());
}


tmp<scalarField>
coupledDiffusionFvPatchScalarField::snGrad(const scalarField&) const
{
    return snGrad();
}


void coupledDiffusionFvPatchScalarField::evaluate
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
        valueInternalCoeffs(patch().weights())*patchInternalField()
      + valueBoundaryCoeffs(patch().weights())*patchNeighbourField()
      + fCorr_
    );

    fvPatchScalarField::evaluate();
}


void coupledDiffusionFvPatchScalarField::updateCoeffs()
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
    const scalarField deltaByDiffusivity
    (
        scalar(1)/stabilise((nfDiffusivity() & patch.nf())*deltaCoeffs, SMALL)
    );

    // Defer the update if neighbour's diffusivity field is not yet available
    const fvPatch& nbrPatch = regionCoupledPatch_.nbrPatch();
    label nbrPatchi = nbrPatch.index();
    const fvMesh& nbrMesh = nbrPatch.boundaryMesh().mesh();
    const coupledDiffusionFvPatchScalarField& nbrPatchField =
        refCast<const coupledDiffusionFvPatchScalarField>
        (
            neighbourFvPatchField()
        );
    if
    (
        (
            nbrPatchField.anisotropic()
         && !nbrMesh.foundObject<volSymmTensorField>
            (
                nbrPatchField.diffusivityName_
            )
        )
     || (
            !nbrPatchField.anisotropic()
         && !nbrMesh.foundObject<volScalarField>(nbrPatchField.diffusivityName_)
        )
    )
    {
        // Restore tag
        UPstream::msgType() = oldTag;
        return;
    }

    const scalarField& nbrDeltaCoeffs =
        nbrMesh.deltaCoeffs().boundaryField()[nbrPatch.index()];
    scalarField nbrDeltaByDiffusivity
    (
        scalar(1)
      / stabilise
        (
            (nbrPatchField.nfDiffusivity() & nbrPatch.nf())*nbrDeltaCoeffs,
            SMALL
        )
    );
    // Interpolate across region boundary. Use this-side value as default in
    // case of missing mapping (low weight) in AMI
    scalarField nbrDeltaByDiffusivityInterp( interpolateFromNeighbour(nbrDeltaByDiffusivity, deltaByDiffusivity) );

    // Set weights
    lWeights_ =
        scalar(1)
      - (
            deltaByDiffusivity
           /(
               deltaByDiffusivity+nbrDeltaByDiffusivityInterp
             + rContactDiffusivity_+nbrPatchField.rContactDiffusivity_
            )
        );

    // Explicit correction to the face interpolated value

    fCorr_ = scalarField(patch.size(), 0);

    // Use the saved gradient to (1) ensure consistency with that used in
    // Laplacian operator, (2) avoid unnecessary recalculation of the whole
    // field and (3) prevent synchronisation issues with consolidated regions
    // if we were to calculate our neighbour region's gradient.

    const fvMesh& mesh(patch.boundaryMesh().mesh());
    if (anisotropic())
    {
        const volScalarField& field =
            mesh.lookupObject<volScalarField>(internalField().name());
        word gradName = "grad(" + field.name() + ")";

        if (mesh.foundObject<volVectorField>(gradName))
        {
            volVectorField& grad =
                mesh.lookupObjectRef<volVectorField>(gradName);

            scalarField corr
            (
                deltaByDiffusivity
                *(nfDiffusivity() & grad.boundaryField()[patchi])
                - (*this - this->patchInternalField())
            );
            fCorr_ -= lWeights_*corr;
        }
    }
    if (nbrPatchField.anisotropic())
    {
        const volScalarField& nbrField =
            nbrMesh.lookupObject<volScalarField>(coupledField());
        word nbrGradName = "grad(" + nbrField.name() + ")";

        if (nbrMesh.foundObject<volVectorField>(nbrGradName))
        {
            volVectorField& nbrGrad =
                nbrMesh.lookupObjectRef<volVectorField>(nbrGradName);

            scalarField nbrCorr
            (
                nbrDeltaByDiffusivity
               *(
                    nbrPatchField.nfDiffusivity()
                  & nbrGrad.boundaryField()[nbrPatchi]
                )
              - (nbrPatchField - nbrPatchField.patchInternalField())
            );
            fCorr_ -=
                (scalar(1)-lWeights_)
               *interpolateFromNeighbour(nbrCorr, scalarField(patch.size(), 0));
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    fvPatchScalarField::updateCoeffs();

    // Update the other side coefficients in case it was deferred
    const_cast<coupledDiffusionFvPatchScalarField&>(nbrPatchField).updateCoeffs();
}


tmp<scalarField>
coupledDiffusionFvPatchScalarField::valueInternalCoeffs
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
coupledDiffusionFvPatchScalarField::valueBoundaryCoeffs
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
coupledDiffusionFvPatchScalarField::valueCorrection() const
{
    return tmp<scalarField>(new scalarField(patch().size(), 0));
}


tmp<scalarField>
coupledDiffusionFvPatchScalarField::gradientInternalCoeffs() const
{
    if (!regionCoupled())
    {
        return regionCoupledFvPatchField<scalar>::gradientInternalCoeffs();
    }

    return patch().deltaCoeffs()
      * (valueInternalCoeffs(patch().weights()) - 1.0);
}


tmp<scalarField>
coupledDiffusionFvPatchScalarField::gradientBoundaryCoeffs() const
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
coupledDiffusionFvPatchScalarField::gradientCorrection() const
{
    return tmp<scalarField>(new scalarField(patch().size(), 0));
}


void coupledDiffusionFvPatchScalarField::regionCoupledBoundaryCoeffs
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


void coupledDiffusionFvPatchScalarField::write(Ostream& os) const
{
    regionCoupledFvPatchField<scalar>::write(os);

    if (diffusivityName_ != "D_"+this->internalField().name())
    {
        os.writeEntry("diffusivityName", diffusivityName_);
    }

    if (diffusivityLayers_.size())
    {
        diffusivityLayers_.writeEntry("diffusivityLayers", os);
        thicknessLayers_.writeEntry("thicknessLayers", os);
    }

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


makePatchTypeField
(
    fvPatchScalarField,
    coupledDiffusionFvPatchScalarField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
