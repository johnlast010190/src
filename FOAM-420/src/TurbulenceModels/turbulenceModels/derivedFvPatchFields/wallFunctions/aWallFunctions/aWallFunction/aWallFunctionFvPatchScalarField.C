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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "aWallFunctionFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "turbulenceModel.H"
#include "matrices/RectangularMatrix/RectangularMatrix.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "derivedFvPatchFields/wallFunctions/analyticalProfiles/FkEBLagAnalyticalProfile/FkEBLagAnalyticalProfile.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

scalar aWallFunctionFvPatchScalarField::tolerance_ = 1e-5;
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void aWallFunctionFvPatchScalarField::checkType()
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

void aWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& a =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = a.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<aWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            aWallFunctionFvPatchScalarField& epf = aPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}

void aWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& a =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = a.boundaryField();

    const fvMesh& mesh = a.mesh();

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

    DynamicList<label> aPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<aWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            aPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                weights[faceCells[i]]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(aPatches, i)
    {
        label patchi = aPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    a_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}

aWallFunctionFvPatchScalarField&
aWallFunctionFvPatchScalarField::aPatch(const label patchi)
{
    const volScalarField& a =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = a.boundaryField();

    const aWallFunctionFvPatchScalarField& epf =
        refCast<const aWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<aWallFunctionFvPatchScalarField&>(epf);
}

void aWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& a0
)
{
    // Accumulate all of a contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            aWallFunctionFvPatchScalarField& epf = aPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), a0);
        }
    }

    // Apply zero-gradient condition for a
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            aWallFunctionFvPatchScalarField& epf = aPatch(patchi);

            epf.forceAssign(scalarField(a0, epf.patch().faceCells()));
        }
    }
}

void aWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& a0
)
{
    const label patchi = patch.index();

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    FkEBLagAnalyticalProfile Fk(lag_);

    forAll(*this, facei)
    {
        label celli = patch.faceCells()[facei];
        scalar x = min
        (
           k[celli]*y[facei]*y[facei]/nuw[facei]/nuw[facei],
           Fk.maxFkValue()
        );
        scalar yPlus = Fk.getYPlus(x);
        if (yPlus >= 17.0)
        {
           a0[celli] = 1./(1. + pow(17./yPlus, 4./3.));
        }
        else
        {
            a0[celli] = 1.0 - exp(-yPlus/24.52);
        }
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

aWallFunctionFvPatchScalarField::aWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    a_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(true)
{
    checkType();
}


aWallFunctionFvPatchScalarField::aWallFunctionFvPatchScalarField
(
    const aWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    a_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(ptf.lag_)
{
    checkType();
}


aWallFunctionFvPatchScalarField::aWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    a_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(dict.lookupOrDefault<Switch>("lag", true))
{
    checkType();
}


aWallFunctionFvPatchScalarField::aWallFunctionFvPatchScalarField
(
    const aWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchField<scalar>(awfpsf),
    a_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(awfpsf.lag_)
{
    checkType();
}


aWallFunctionFvPatchScalarField::aWallFunctionFvPatchScalarField
(
    const aWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(awfpsf, iF),
    a_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(awfpsf.lag_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalarField& aWallFunctionFvPatchScalarField::a
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            a_ = 0.0;
        }

        return a_;
    }

    return aPatch(master_).a(init);
}

void aWallFunctionFvPatchScalarField::updateCoeffs()
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
        calculateTurbulenceFields(turbModel, a(true));
    }

    const scalarField& a0 = this->a();
    typedef DimensionedField<scalar, volMesh> FieldType;
    FieldType& a = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];
        a[celli] = a0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void aWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
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
        calculateTurbulenceFields(turbModel, a(true));
    }

    const scalarField& a0 = this->a();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& a = const_cast<FieldType&>(internalField());

    scalarField& af = *this;

    // Only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            a[celli] = (1.0 - w)*a[celli] + w*a0[celli];
            af[facei] = a[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}

void aWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}

void aWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField()());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}

void aWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }
    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constrainta(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& a
        = internalField();

    label nConstrainedCells = 0;


    forAll(weights, facei)
    {
        // Only set the values if the weights are > tolerance
        if (weights[facei] > 1e-05)
        {
            nConstrainedCells++;

            label celli = faceCells[facei];

            constraintCells.append(celli);
            constrainta.append(a[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << nConstrainedCells
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues
    (
        constraintCells,
        scalarField(constrainta.xfer())
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void aWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    aWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
