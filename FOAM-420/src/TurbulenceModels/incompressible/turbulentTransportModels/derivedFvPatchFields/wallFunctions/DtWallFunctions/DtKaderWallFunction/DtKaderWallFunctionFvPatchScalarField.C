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
    (c) 2010-2016 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "DtKaderWallFunctionFvPatchScalarField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "incompressibleTurbulenceModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
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
    const incompressibleTurbulenceModel& turbModel =
        db().lookupObject<incompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    label patchI = patch().index();

    const IOdictionary& transportProperties
        = db().lookupObject<IOdictionary>("transportProperties");


    string pofix(string(this->internalField().name())(2, 100));

    //read laminar vapour/mixture diffusion coefficient
    dimensionedScalar D
    (
        "D" + pofix,
        transportProperties.lookup("D" + pofix)
    );

    //now calc Tplus
    fvPatchField<scalar> nuw
        = turbModel.nu()().boundaryField()[patchI];

    scalarField Sc( nuw/D.value() );

    scalarField yPlus( max(turbModel.yPlus(patchI, Cmu_), SMALL) );

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

    scalarField Tplus( Sc*yPlus*exp(Gamma)+(2.12*log(yPlus) + Beta)* exp(1.0/Gamma) );

    scalarField Deffw( yPlus*nuw/Tplus );

    forceAssign(max(scalar(0.0), Deffw - D.value()));
}


void DtKaderWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntryIfDifferent<scalar>(os, "Cmu", 0.09, Cmu_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchScalarField, DtKaderWallFunctionFvPatchScalarField);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
