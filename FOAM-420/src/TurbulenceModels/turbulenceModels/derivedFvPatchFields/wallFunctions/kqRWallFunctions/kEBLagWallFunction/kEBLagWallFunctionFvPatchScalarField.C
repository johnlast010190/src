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

#include "derivedFvPatchFields/wallFunctions/nutWallFunctions/nutWallFunction/nutWallFunctionFvPatchScalarField.H"
#include "kEBLagWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "derivedFvPatchFields/wallFunctions/analyticalProfiles/FkEBLagAnalyticalProfile/FkEBLagAnalyticalProfile.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void kEBLagWallFunctionFvPatchScalarField::checkType()
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


void kEBLagWallFunctionFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kEBLagWallFunctionFvPatchScalarField::kEBLagWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchField<scalar>(p, iF)
{
    checkType();
}


kEBLagWallFunctionFvPatchScalarField::kEBLagWallFunctionFvPatchScalarField
(
    const kEBLagWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    zeroGradientFvPatchField<scalar>(ptf, p, iF, mapper)
{
    checkType();
}


kEBLagWallFunctionFvPatchScalarField::kEBLagWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    zeroGradientFvPatchField<scalar>(p, iF, dict)
{
    checkType();
}


kEBLagWallFunctionFvPatchScalarField::kEBLagWallFunctionFvPatchScalarField
(
    const kEBLagWallFunctionFvPatchScalarField& v2wfpsf
)
:
    zeroGradientFvPatchField<scalar>(v2wfpsf)
{
    checkType();
}


kEBLagWallFunctionFvPatchScalarField::kEBLagWallFunctionFvPatchScalarField
(
    const kEBLagWallFunctionFvPatchScalarField& v2wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    zeroGradientFvPatchField<scalar>(v2wfpsf, iF)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void kEBLagWallFunctionFvPatchScalarField::updateCoeffs()
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

    FkEBLagAnalyticalProfile Fk(true);

    const label patchi = patch().index();
    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    typedef DimensionedField<scalar, volMesh> FieldType;
    FieldType& Sk =
        const_cast<FieldType&>
        (
            db().lookupObject<FieldType>("Sk")
        );

    const volScalarField& nut =
        db().lookupObject<volScalarField>("nut");

    const volScalarField::Boundary& bf = nut.boundaryField();

    const nutWallFunctionFvPatchScalarField& nutpf =
        refCast<const nutWallFunctionFvPatchScalarField>(bf[patchi]);

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));

    forAll(*this, facei)
    {
        label celli = patch().faceCells()[facei];

        scalar x = min
        (
           k[celli]*y[facei]*y[facei]/nuw[facei]/nuw[facei],
           Fk.maxFkValue()
        );

        scalar yPlus = Fk.getYPlus(x);
        scalar rH = nutpf.roughWallCoeffs().roughnessHeight();
        scalar uk = yPlus*nuw[facei]/y[facei];

        scalar KsPlus = uk*rH/nuw[facei];

        scalar Uplus = 1./0.41*log(1+0.41*yPlus)
          + 7.8*(1.-exp(-yPlus/11.)-yPlus/11.*exp(-yPlus/3.));
        scalar Du = 1./0.41*log(1 + KsPlus/exp(3.25*0.41));
        Uplus -= Du;
        scalar tw = nuw[facei]*magGradU[facei]*yPlus/Uplus;

        scalar nutPlus = 3*0.41/2.*yPlus*(1-exp(-yPlus/26.))*(1-exp(-yPlus/26.));
        scalar PPlus = 0.95*nutPlus/(1+nutPlus)/(1+nutPlus);
        Sk[celli] = PPlus*uk*uk*tw/nuw[facei];
    }

    zeroGradientFvPatchField<scalar>::updateCoeffs();
}

void kEBLagWallFunctionFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    zeroGradientFvPatchField<scalar>::evaluate(commsType);
}


void kEBLagWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    zeroGradientFvPatchField<scalar>::write(os);
    this->writeEntry("value", os);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    makePatchTypeField
    (
        fvPatchScalarField,
        kEBLagWallFunctionFvPatchScalarField
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
