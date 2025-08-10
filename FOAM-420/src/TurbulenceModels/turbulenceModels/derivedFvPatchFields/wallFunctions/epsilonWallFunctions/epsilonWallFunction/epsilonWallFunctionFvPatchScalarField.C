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
    (c) 2008-2010, 2017 OpenCFD Ltd.
    (c) 2023-2024 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "epsilonWallFunctionFvPatchScalarField.H"
#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::epsilonWallFunctionFvPatchScalarField::tolerance_ = 0.1;


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::epsilonWallFunctionFvPatchScalarField::checkType()
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


void Foam::epsilonWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("E", E_);
    os.writeEntry("lowReCorrection", lowReCorrection_);
}


void Foam::epsilonWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void Foam::epsilonWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

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

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                weights[faceCells[i]]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(epsilonPatches, i)
    {
        label patchi = epsilonPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), Zero);
    epsilon_.setSize(internalField().size(), Zero);

    initialised_ = true;
}


Foam::epsilonWallFunctionFvPatchScalarField&
Foam::epsilonWallFunctionFvPatchScalarField::epsilonPatch(const label patchi)
{
    const volScalarField& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const epsilonWallFunctionFvPatchScalarField& epf =
        refCast<const epsilonWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<epsilonWallFunctionFvPatchScalarField&>(epf);
}


void Foam::epsilonWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& G0,
    scalarField& epsilon0
)
{
    // Ask for the wall distance in advance. Some processes may not have any
    // corner weights, so we cannot allow functions inside the if statements
    // below to trigger the wall distance calculation.
    turbulence.y();

    // Accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // Apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            epf.forceAssign(scalarField(epsilon0, epf.patch().faceCells()));
        }
    }
}


void Foam::epsilonWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();

    const scalarField& y = turbModel.y()[patchi];

    const scalar Cmu25 = pow025(Cmu_);
    const scalar Cmu75 = pow(Cmu_, 0.75);

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    // Set epsilon and G
    forAll(nutw, facei)
    {
        const label celli = patch.faceCells()[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];

        const scalar w = cornerWeights[facei];

        // Default high-Re form
        scalar epsilonc = w*Cmu75*pow(k[celli], 1.5)/(kappa_*y[facei]);
        scalar Gc =
            w
           *(nutw[facei] + nuw[facei])
           *magGradUw[facei]
           *Cmu25*sqrt(k[celli])
           /(kappa_*y[facei]);

        if (lowReCorrection_ && yPlus < yPlusLam_)
        {
            epsilonc = w*2.0*k[celli]*nuw[facei]/sqr(y[facei]);
            Gc = 0;
        }

        epsilon0[celli] += epsilonc;

        G0[celli] += Gc;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    yPlusLam_(nutWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    epsilon_(),
    lowReCorrection_(false),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    yPlusLam_(ptf.yPlusLam_),
    G_(),
    epsilon_(),
    lowReCorrection_(ptf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
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
    yPlusLam_(nutWallFunctionFvPatchScalarField::yPlusLam(kappa_, E_)),
    G_(),
    epsilon_(),
    lowReCorrection_(dict.lookupOrDefault<Switch>("lowReCorrection", false)),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();

    // Apply zero-gradient condition on start-up
    this->forceAssign(patchInternalField());
}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_),
    yPlusLam_(ewfpsf.yPlusLam_),
    G_(),
    epsilon_(),
    lowReCorrection_(ewfpsf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    Cmu_(ewfpsf.Cmu_),
    kappa_(ewfpsf.kappa_),
    E_(ewfpsf.E_),
    yPlusLam_(ewfpsf.yPlusLam_),
    G_(),
    epsilon_(),
    lowReCorrection_(ewfpsf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    cornerWeights_()
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::epsilonWallFunctionFvPatchScalarField::G(bool init)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return epsilonPatch(master_).G();
}


Foam::scalarField& Foam::epsilonWallFunctionFvPatchScalarField::epsilon
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void Foam::epsilonWallFunctionFvPatchScalarField::updateCoeffs()
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
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& epsilon = const_cast<FieldType&>(internalField());

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
        epsilon[celli] = (1 - w)*epsilon[celli] + w*epsilon0[celli];
    }

    this->forceAssign(scalarField(epsilon, patch().faceCells()));

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonWallFunctionFvPatchScalarField::manipulateMatrix
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


void Foam::epsilonWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        epsilonWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
