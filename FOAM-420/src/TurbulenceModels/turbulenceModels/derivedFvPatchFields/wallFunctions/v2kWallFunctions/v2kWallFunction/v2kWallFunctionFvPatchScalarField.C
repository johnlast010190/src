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

#include "v2kWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "matrices/RectangularMatrix/RectangularMatrix.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "derivedFvPatchFields/wallFunctions/analyticalProfiles/v2kEBLagAnalyticalProfile/v2kEBLagAnalyticalProfile.H"
#include "derivedFvPatchFields/wallFunctions/analyticalProfiles/FkEBLagAnalyticalProfile/FkEBLagAnalyticalProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

scalar v2kWallFunctionFvPatchScalarField::tolerance_ = 1e-5;
// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void v2kWallFunctionFvPatchScalarField::checkType()
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


void v2kWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
}

void v2kWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const volScalarField& v2k =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = v2k.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<v2kWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            v2kWallFunctionFvPatchScalarField& epf = v2kPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void v2kWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const volScalarField& v2k =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = v2k.boundaryField();

    const fvMesh& mesh = v2k.mesh();

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

    DynamicList<label> v2kPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<v2kWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            v2kPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            forAll(faceCells, i)
            {
                weights[faceCells[i]]++;
            }
        }
    }

    cornerWeights_.setSize(bf.size());
    forAll(v2kPatches, i)
    {
        label patchi = v2kPatches[i];
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    v2k_.setSize(internalField().size(), 0.0);

    initialised_ = true;
}


v2kWallFunctionFvPatchScalarField&
v2kWallFunctionFvPatchScalarField::v2kPatch(const label patchi)
{
    const volScalarField& v2k =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = v2k.boundaryField();

    const v2kWallFunctionFvPatchScalarField& epf =
        refCast<const v2kWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<v2kWallFunctionFvPatchScalarField&>(epf);
}

void v2kWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& v2k0
)
{
    // Accumulate all of a contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            v2kWallFunctionFvPatchScalarField& epf = v2kPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), v2k0);
        }
    }

    // Apply zero-gradient condition for a
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            v2kWallFunctionFvPatchScalarField& epf = v2kPatch(patchi);

            epf.forceAssign(scalarField(v2k0, epf.patch().faceCells()));
        }
    }
}

void v2kWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& v2k0
)
{
    const label patchi = patch.index();

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();


    v2kEBLagAnalyticalProfile v2kProf;
    FkEBLagAnalyticalProfile Fk(lag_);

    // Set v2/k wall values
    forAll(*this, facei)
    {
        label celli = patch.faceCells()[facei];
        scalar x = min
        (
           k[celli]*y[facei]*y[facei]/nuw[facei]/nuw[facei],
           Fk.maxFkValue()
        );

        scalar yPlus = Fk.getYPlus(x);
        v2k0[celli] = v2kProf.getv2k(yPlus);
    }
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

v2kWallFunctionFvPatchScalarField::v2kWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    v2k_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(true)
{
    checkType();
}


v2kWallFunctionFvPatchScalarField::v2kWallFunctionFvPatchScalarField
(
    const v2kWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    v2k_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(ptf.lag_)
{
    checkType();
}


v2kWallFunctionFvPatchScalarField::v2kWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    v2k_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(dict.lookupOrDefault<Switch>("lag", true))
{
    checkType();
}


v2kWallFunctionFvPatchScalarField::v2kWallFunctionFvPatchScalarField
(
    const v2kWallFunctionFvPatchScalarField& v2kwfpsf
)
:
    fixedValueFvPatchField<scalar>(v2kwfpsf),
    v2k_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(v2kwfpsf.lag_)
{
    checkType();
}


v2kWallFunctionFvPatchScalarField::v2kWallFunctionFvPatchScalarField
(
    const v2kWallFunctionFvPatchScalarField& v2wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(v2wfpsf, iF),
    v2k_(),
    initialised_(false),
    master_(-1),
    cornerWeights_(),
    lag_(v2wfpsf.lag_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::scalarField& v2kWallFunctionFvPatchScalarField::v2k
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            v2k_ = 0.0;
        }

        return v2k_;
    }

    return v2kPatch(master_).v2k(init);
}

void v2kWallFunctionFvPatchScalarField::updateCoeffs()
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
        calculateTurbulenceFields(turbModel, v2k(true));
    }

    const scalarField& v2k0 = this->v2k();
    typedef DimensionedField<scalar, volMesh> FieldType;
    FieldType& v2k = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];
        v2k[celli] = v2k0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}

void v2kWallFunctionFvPatchScalarField::updateWeightedCoeffs
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
        calculateTurbulenceFields(turbModel, v2k(true));
    }

    const scalarField& v2k0 = this->v2k();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& v2k = const_cast<FieldType&>(internalField());

    scalarField& v2kf = *this;

    // Only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        scalar w = weights[facei];

        if (w > tolerance_)
        {
            label celli = patch().faceCells()[facei];

            v2k[celli] = (1.0 - w)*v2k[celli] + w*v2k0[celli];
            v2kf[facei] = v2k[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}

void v2kWallFunctionFvPatchScalarField::manipulateMatrix
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

void v2kWallFunctionFvPatchScalarField::manipulateMatrix
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
    DynamicList<scalar> constraintV2k(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& v2k
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
            constraintV2k.append(v2k[celli]);
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
        scalarField(constraintV2k.xfer())
    );

    fvPatchField<scalar>::manipulateMatrix(matrix);
}

void v2kWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void v2kWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    v2kWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
