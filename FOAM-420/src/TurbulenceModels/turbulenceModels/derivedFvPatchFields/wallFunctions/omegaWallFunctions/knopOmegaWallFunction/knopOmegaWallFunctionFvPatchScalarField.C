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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2010-2024 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "knopOmegaWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutkWallFunction/nutkWallFunctionFvPatchScalarField.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar knopOmegaWallFunctionFvPatchScalarField::tolerance_ = 0.1;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void knopOmegaWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void knopOmegaWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("E", E_);
    os.writeEntry("beta1", beta1_);
}

void knopOmegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const auto& omega =
        dynamic_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<knopOmegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            knopOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}

void knopOmegaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    bool registerField = true;
    const auto& omega =
        dynamic_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const fvMesh& mesh = omega.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            !registerField
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<knopOmegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCellsList = bf[patchi].patch().faceCells();
            for (int celli : faceCellsList)
            {
                weights[celli]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(omegaPatches, i)
    {
        label patchi = omegaPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), 0.0);
    omega_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}

knopOmegaWallFunctionFvPatchScalarField&
knopOmegaWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const auto& omega =
        dynamic_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const auto& opf =
        refCast<const knopOmegaWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<knopOmegaWallFunctionFvPatchScalarField&>(opf);
}

void knopOmegaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            knopOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            knopOmegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf.forceAssign(scalarField(omega0, opf.patch().faceCells()));
        }
    }
}

void knopOmegaWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbulenceModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchI = patch.index();
    const scalarField& y = turbulenceModel.y()[patchI];

    const scalar Cmu25 = pow025(Cmu_);

    typedef DimensionedField<scalar, volMesh> FieldType;
    const auto& G =
      db().lookupObject<FieldType>(turbulenceModel.GName());

    const tmp<volScalarField> tk = turbulenceModel.k();
    const volScalarField& k = tk();

    const tmp<volScalarField> tnu = turbulenceModel.nu();
    const scalarField& nuw = tnu().boundaryField()[patchI];

    const tmp<volScalarField> tnut = turbulenceModel.nut();
    const volScalarField& nut = tnut();
    const scalarField& nutw = nut.boundaryField()[patchI];
    const fvPatchScalarField& nutwpf = nut.boundaryField()[patchI];

    const fvPatchVectorField& Uw = turbulenceModel.U().boundaryField()[patchI];
    const scalarField magGradUw(mag(Uw.snGrad()));

    const auto& nutWF =
        dynamic_cast<const nutWallFunctionFvPatchScalarField&>(nutwpf);

    tmp<scalarField> yp  = nutWF.yPlus();
    vectorField n ( patch.nf() );

    // Set omega and G
    scalar omegaWallValue;
    scalar roughnessHeight = nutWF.roughWallCoeffs().roughnessHeight();
    forAll(nutw, facei)
    {
        label faceCellI = patch.faceCells()[facei];
        const scalar w = cornerWeights[facei];

        scalar omegaVis = 6.0*nuw[facei]/(beta1_*sqr(y[facei]));
        scalar uStar = yp()[facei]*nuw[facei]/y[facei];
        scalar KsPlus = uStar*roughnessHeight/nuw[facei];

        scalar Fr2 =
          min(1, pow((KsPlus/30), (2./3.)))
         *min(1, pow((KsPlus/45), (1./4.)))
         *min(1, pow((KsPlus/60), (1./4.)));
        scalar d0 = Fr2 * 0.03 * roughnessHeight;

        scalar omegaLog = sqrt(k[faceCellI])/(Cmu25*kappa_*y[facei]);

        scalar omegab1 = omegaVis + omegaLog;

        scalar omegab2 = pow(pow(omegaVis,1.2) + pow(omegaLog,1.2),1/1.2);

        scalar yPlus = yp()[facei];
        scalar arg = 0.1*yPlus;

        scalar phi = tanh(pow(arg,4));

        omegaWallValue =  phi*omegab1 + (1-phi)*omegab2;

        omegaWallValue = min
        (
            uStar/(sqrt(Cmu_)*kappa_*max(d0, ROOTVSMALL)),
            omegaWallValue
        );
        omega0[faceCellI] += w*omegaWallValue;

        const scalar Rey = y[facei]*sqrt(k[faceCellI])/nuw[facei];
        const scalar lamFrac = exp(-Rey/11);
        const scalar turbFrac = 1 - lamFrac;

        scalar lamPart = G[faceCellI];
        scalar turbPart = 0;

        //- Sanity checks, they can be triggered in SP
        if (yPlus!=0)
        {
            const scalar uPlus =  (1/kappa_)*log(E_*yPlus);
            if ((uPlus!=0))
            {
                lamPart *= lamFrac;
                turbPart = turbFrac * sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                    /(nuw[facei]*kappa_*yPlus);
            }
        }
        G0[faceCellI] += w*(lamPart + turbPart);
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

knopOmegaWallFunctionFvPatchScalarField::knopOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    beta1_(0.075),
    yPlusLam_(nutkWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


knopOmegaWallFunctionFvPatchScalarField::knopOmegaWallFunctionFvPatchScalarField
(
    const knopOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    beta1_(ptf.beta1_),
    yPlusLam_(ptf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


knopOmegaWallFunctionFvPatchScalarField::knopOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    beta1_(dict.lookupOrDefault<scalar>("beta1", 0.075)),
    yPlusLam_(nutkWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


knopOmegaWallFunctionFvPatchScalarField::knopOmegaWallFunctionFvPatchScalarField
(
    const knopOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    yPlusLam_(owfpsf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


knopOmegaWallFunctionFvPatchScalarField::knopOmegaWallFunctionFvPatchScalarField
(
    const knopOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    yPlusLam_(owfpsf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
scalarField& knopOmegaWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return omegaPatch(master_).G();
}


scalarField& knopOmegaWallFunctionFvPatchScalarField::omega(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            omega_ = 0.0;
        }

        return omega_;
    }

    return omegaPatch(master_).omega(init);
}

void knopOmegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), omega(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& omega0 = this->omega();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>(turbModel.GName())
        );

    auto& omega = const_cast<FieldType&>(internalField());

    const scalarField::subField& magFA = patch().patch().magFaceAreas();
    const scalarField& magSf = patch().magSf();

    scalarField weights(patch().size());
    forAll(weights, facei)
    {
        scalar& w = weights[facei];

        w = magFA[facei] > SMALL ? magSf[facei]/magFA[facei] : 1;

        w = w <= tolerance_ ? 0 : (w - tolerance_)/(1 - tolerance_);
    }

    forAll(weights, facei)
    {
        const scalar w = weights[facei];
        const label celli = patch().faceCells()[facei];

        G[celli] = (1 - w)*G[celli] + w*G0[celli];
        omega[celli] = (1 - w)*omega[celli] + w*omega0[celli];
    }

    this->forceAssign(scalarField(omega, patch().faceCells()));

    fvPatchField<scalar>::updateCoeffs();
}


void knopOmegaWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    const scalarField::subField& magFA = patch().patch().magFaceAreas();
    const scalarField& magSf = patch().magSf();

    scalarField weights(patch().size());
    forAll(weights, facei)
    {
        scalar& w = weights[facei];

        w = magFA[facei] > SMALL ? magSf[facei]/magFA[facei] : 1;

        w = w <= tolerance_ ? 0 : (w - tolerance_)/(1 - tolerance_);
    }

    matrix.setValues(patch().faceCells(), patchInternalField()(), weights);

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void knopOmegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    knopOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
