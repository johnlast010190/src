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
    (c) 2017 OpenCFD Ltd.
    (c) 2011-2019 OpenFOAM Foundation
    (c) 2010-2024 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "omegaWallFunctionFvPatchScalarField.H"
#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar omegaWallFunctionFvPatchScalarField::tolerance_ = 0.1;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void omegaWallFunctionFvPatchScalarField::checkType()
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


void omegaWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("E", E_);
    os.writeEntry("beta1", beta1_);
    os.writeEntry("blended", blended_);
    os.writeEntry("legacy", legacy_);
}


void omegaWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            opf.master() = master;
        }
    }
}


void omegaWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

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
            false // do not register
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    );

    DynamicList<label> omegaPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<omegaWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            omegaPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                label celli = faceCells[i];
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


omegaWallFunctionFvPatchScalarField&
omegaWallFunctionFvPatchScalarField::omegaPatch(const label patchi)
{
    const volScalarField& omega =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = omega.boundaryField();

    const omegaWallFunctionFvPatchScalarField& opf =
        refCast<const omegaWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<omegaWallFunctionFvPatchScalarField&>(opf);
}


void omegaWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbModel,
    scalarField& G0,
    scalarField& omega0
)
{
    // Ask for the wall distance in advance. Some processes may not have any
    // corner weights, so we cannot allow functions inside the if statements
    // below to trigger the wall distance calculation.
    turbModel.y();

    // Accumulate all of the G and omega contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            opf.calculate(turbModel, w, opf.patch(), G0, omega0);
        }
    }

    // Apply zero-gradient condition for omega
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            omegaWallFunctionFvPatchScalarField& opf = omegaPatch(patchi);

            opf.forceAssign(scalarField(omega0, opf.patch().faceCells()));
        }
    }
}


void omegaWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& omega0
)
{
    const label patchi = patch.index();

    const scalarField& y = turbModel.y()[patchi];

    const scalar Cmu25 = pow025(Cmu_);
    const scalar Cmu5 = sqrt(Cmu_);

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    typedef DimensionedField<scalar, volMesh> FieldType;
    const FieldType& G =
      db().lookupObject<FieldType>(turbModel.GName());

    // Set omega and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];
        const scalar w = cornerWeights[facei];

        const scalar Rey = y[facei]*sqrt(k[celli])/nuw[facei];
        const scalar yPlus = Cmu25*Rey;
        const scalar uPlus = (1/kappa_)*log(E_*yPlus);

        if (blended_)
        {
            if (!legacy_)
            {
                const scalar lamFrac = exp(-Rey/11);
                const scalar turbFrac = 1 - lamFrac;

                const scalar uStar = sqrt
                (
                    lamFrac*nuw[facei]*magGradUw[facei] + turbFrac*Cmu5*k[celli]
                );

                const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));
                const scalar omegaLog = uStar/(Cmu5*kappa_*y[facei]);

                omega0[celli] += w*(lamFrac*omegaVis + turbFrac*omegaLog);

                scalar lamPart = G[celli];
                scalar turbPart = 0;

                //- Sanity checks, they can be triggered in SP
                if ((uPlus!=0) && (yPlus!=0))
                {
                    lamPart *= lamFrac;
                    turbPart =
                        turbFrac
                        *sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                        /(nuw[facei]*kappa_*yPlus);
                }

                G0[celli] += w*(lamPart + turbPart);
            }
            else
            {
                const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));
                const scalar omegaLog = sqrt(k[celli])/(Cmu25*kappa_*y[facei]);
                omega0[celli] += w*sqrt(sqr(omegaVis) + sqr(omegaLog));

                G0[celli] +=
                    w
                    *(nutw[facei] + nuw[facei])
                    *magGradUw[facei]
                    *Cmu25*sqrt(k[celli])
                    /(kappa_*y[facei]);
            }
        }
        else // not blended
        {
            if (yPlus < yPlusLam_) // viscous layer
            {
                const scalar omegaVis = 6*nuw[facei]/(beta1_*sqr(y[facei]));

                omega0[celli] += w*omegaVis;

                G0[celli] += w*G[celli];
            }
            else
            {
                const scalar uStar = sqrt(Cmu5*k[celli]);
                const scalar omegaLog = uStar/(Cmu5*kappa_*y[facei]);

                omega0[celli] += w*omegaLog;

                G0[celli] +=
                    w*
                    sqr(uStar*magGradUw[facei]*y[facei]/uPlus)
                    /(nuw[facei]*kappa_*yPlus);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
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
    blended_(true),
    legacy_(false),
    yPlusLam_(nutWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& ptf,
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
    blended_(ptf.blended_),
    legacy_(ptf.legacy_),
    yPlusLam_(ptf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
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
    blended_(dict.lookupOrDefault<Switch>("blended", true)),
    legacy_(dict.lookupOrDefault<Switch>("legacy", false)),
    yPlusLam_(nutWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // apply zero-gradient condition on start-up
    this->forceAssign(patchInternalField());
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedValueFvPatchField<scalar>(owfpsf),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    legacy_(owfpsf.legacy_),
    yPlusLam_(owfpsf.yPlusLam_),
    G_(),
    omega_(),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


omegaWallFunctionFvPatchScalarField::omegaWallFunctionFvPatchScalarField
(
    const omegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(owfpsf, iF),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    beta1_(owfpsf.beta1_),
    blended_(owfpsf.blended_),
    legacy_(owfpsf.legacy_),
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

scalarField& omegaWallFunctionFvPatchScalarField::G(bool init)
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


scalarField& omegaWallFunctionFvPatchScalarField::omega(bool init)
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


void omegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
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

    FieldType& omega = const_cast<FieldType&>(internalField());

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


void omegaWallFunctionFvPatchScalarField::manipulateMatrix
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


void omegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    omegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
