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
    (c) 2010-2016, Esi Ltd

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatKaderWallFunction/alphatKaderWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void alphatKaderWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatKaderWallFunctionFvPatchScalarField::
alphatKaderWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Cmu_(0.09)
{
    checkType();
}


alphatKaderWallFunctionFvPatchScalarField::
alphatKaderWallFunctionFvPatchScalarField
(
    const alphatKaderWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_)
{
    checkType();
}


alphatKaderWallFunctionFvPatchScalarField::
alphatKaderWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09))
{
    checkType();
}


alphatKaderWallFunctionFvPatchScalarField::
alphatKaderWallFunctionFvPatchScalarField
(
    const alphatKaderWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Cmu_(awfpsf.Cmu_)
{
    checkType();
}


alphatKaderWallFunctionFvPatchScalarField::
alphatKaderWallFunctionFvPatchScalarField
(
    const alphatKaderWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Cmu_(awfpsf.Cmu_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatKaderWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const compressible::turbulenceModel& turbModel =
        db().lookupObject<compressible::turbulenceModel>
    (
        IOobject::groupName
        (
            compressible::turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    // Field data
    const label patchi = patch().index();

    const tmp<scalarField> talphaw = turbModel.alpha(patchi);
    const scalarField& alphaw = talphaw();

    //now calc Tplus
    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    scalarField Pr( muw/alphaw );

    scalarField yPlus( max(turbModel.yPlus(patchi, Cmu_), SMALL) ); //NOTE: This was 0.1 instead of SMALL

    scalarField Beta
    (
        sqr(3.85*pow(Pr, 0.333) - 1.3)
        + 2.12*log(Pr)
    );

    scalarField Gamma
    (
        -0.01 * pow4(Pr*yPlus)
        /(1.0 + 5.0*yPlus*pow3(Pr))
    );

    scalarField Tplus
    (
        Pr*yPlus*exp(Gamma)
        +(2.12*log(1.0 + yPlus) + Beta)* exp(1.0/Gamma)
    );

    scalarField alphaEffw( yPlus*muw/Tplus );

    forceAssign(max(scalar(0.0), alphaEffw - alphaw));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void alphatKaderWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<scalar>(os, "Cmu", 0.09, Cmu_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatKaderWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
