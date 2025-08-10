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

#include "turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtKaderWallFunction/DtKaderWallFunctionFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DtKaderWallFunctionFvPatchScalarField::
DtKaderWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Cmu_(0.09)
{}


DtKaderWallFunctionFvPatchScalarField::
DtKaderWallFunctionFvPatchScalarField
(
    const DtKaderWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_)
{}


DtKaderWallFunctionFvPatchScalarField::
DtKaderWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09))
{}


DtKaderWallFunctionFvPatchScalarField::
DtKaderWallFunctionFvPatchScalarField
(
    const DtKaderWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Cmu_(awfpsf.Cmu_)
{}


DtKaderWallFunctionFvPatchScalarField::
DtKaderWallFunctionFvPatchScalarField
(
    const DtKaderWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Cmu_(awfpsf.Cmu_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DtKaderWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Retrieve turbulence properties from model
    const compressible::turbulenceModel& turbModel =
        db().lookupObject<compressible::turbulenceModel>
        (
            IOobject::groupName
            (
                compressible::turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    string pofix(string(this->internalField().name())(2, 100));
    const word subDictName(pofix + "concentrationTransport");
    bool isMat = (basicThermo::dictName == basicThermo::matDictName);

    //read laminar mass diffusivity
    dimensionedScalar D
    (
        isMat ? "D" : "D" + pofix,
        isMat
      ? turbModel.transport().optionalSubDict(subDictName).lookup("D")
      : turbModel.transport().lookup("D" + pofix)
    );

    //read turbulent Schmidt number
    dimensionedScalar Sct
    (
        isMat ? "Sct" : "Sct" + pofix,
        isMat
      ? turbModel.transport().optionalSubDict(subDictName).lookup("Sct")
      : turbModel.transport().lookup("Sct" + pofix)
    );


    label patchI = patch().index();

    //now calc Tplus
    tmp<scalarField> tmuw = turbModel.mu(patchI);
    const scalarField& muw = tmuw();
    const scalarField& rhow = turbModel.rho().boundaryField()[patchI];

    scalarField Sc( muw/(rhow*D.value()) );

    scalarField yPlus( max(turbModel.yPlus(patchI, Cmu_), scalar(0.1)) );

    scalarField Beta
    (
        sqr(3.85*pow(Sc, 0.333) - 1.3)
        + 2.12*log(Sc)
    );

    scalarField Gamma
    (
        -0.01 * pow4(Sc*yPlus)
        /(1.0 + 5.0*yPlus*pow3(Sc))
    );

    scalarField Mplus
    (
        Sc*yPlus*exp(Gamma)
        +(2.12*log(1.0 + yPlus) + Beta)* exp(1.0/Gamma)
    );

    scalarField Deffw( yPlus*muw/Mplus );

    forceAssign(max(scalar(0.0), Deffw - rhow*D.value()));

    fixedValueFvPatchScalarField::updateCoeffs();
}


void DtKaderWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<scalar>(os, "Cmu", 0.09, Cmu_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, DtKaderWallFunctionFvPatchScalarField);

// New name
addNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    DtKaderWallFunctionFvPatchScalarField,
    DtKaderWallFunction
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
