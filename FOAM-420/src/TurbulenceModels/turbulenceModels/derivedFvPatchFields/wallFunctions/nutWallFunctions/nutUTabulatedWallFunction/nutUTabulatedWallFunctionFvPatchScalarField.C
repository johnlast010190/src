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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nutUTabulatedWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "global/etcFiles/etcFiles.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
Foam::fileName Foam::nutUTabulatedWallFunctionFvPatchScalarField::getTableFileName
(
    const dictionary& dict
)
{
    if (fileName_ != "default")
    {
        return ("constant"/fileName_);
    }
    return findEtcFile("turbulenceData/kOmegaSST.data");
}

tmp<scalarField> nutUTabulatedWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    vectorField nf(patch().nf()());
    scalarField magGradU(this->size(), 0.0);
    {
        vectorField gradU(Uw.snGrad());
        gradU -= nf*(gradU & nf);
        magGradU = mag(gradU);
        magGradU = max(magGradU, ROOTVSMALL);
    }
    const scalarField magUp ( magGradU/this->patch().deltaCoeffs() );
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return
        max
        (
            scalar(0),
            (
                    sqr(calcuTau(yPlus(), magUp))/(magGradU+ROOTVSMALL)
                    - nuw
            )
        );
}


tmp<scalarField> nutUTabulatedWallFunctionFvPatchScalarField::calcuTau
(
    const scalarField& yPlus,
    const scalarField& Up
) const
{
    tmp<scalarField> tuTau(new scalarField(patch().size(), 0.0));
    scalarField& uTau = tuTau.ref();

    forAll(uTau, fi)
    {
        uTau[fi] = Up[fi]/uPlusTable_(yPlus[fi]);
    }
    return tuTau;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUTabulatedWallFunctionFvPatchScalarField::
nutUTabulatedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    fileName_("undefined-Name"),
    uPlusTableFileData_("undefined-uPlusTableName"),
    uPlusTable_(uPlusTableFileData_)
{}


nutUTabulatedWallFunctionFvPatchScalarField::
nutUTabulatedWallFunctionFvPatchScalarField
(
    const nutUTabulatedWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    fileName_(ptf.fileName_),
    uPlusTableFileData_(ptf.uPlusTableFileData_),
    uPlusTable_(ptf.uPlusTable_)
{}


nutUTabulatedWallFunctionFvPatchScalarField::
nutUTabulatedWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    fileName_(dict.lookupOrDefault<word>("dataFile", "default")),
    uPlusTableFileData_(getTableFileName(dict)),
    uPlusTable_(uPlusTableFileData_)
{}


nutUTabulatedWallFunctionFvPatchScalarField::
nutUTabulatedWallFunctionFvPatchScalarField
(
    const nutUTabulatedWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    fileName_(wfpsf.fileName_),
    uPlusTableFileData_(wfpsf.uPlusTableFileData_),
    uPlusTable_(wfpsf.uPlusTable_)
{}


nutUTabulatedWallFunctionFvPatchScalarField::
nutUTabulatedWallFunctionFvPatchScalarField
(
    const nutUTabulatedWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    fileName_(wfpsf.fileName_),
    uPlusTableFileData_(wfpsf.uPlusTableFileData_),
    uPlusTable_(wfpsf.uPlusTable_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutUTabulatedWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    vectorField nf(patch().nf()());

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchi];

    scalarField magGradU(this->size(), 0.0);
    {
        vectorField gradU(Uw.snGrad());
        gradU -= nf*(gradU & nf);
        magGradU = mag(gradU);
        magGradU = max(magGradU, ROOTVSMALL);
    }

    return y*sqrt((nuw + *this) * magGradU)/nuw;
}


void nutUTabulatedWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeEntry("value", os);
    writeEntryIfDifferent
    (
        os,
        "dataFile",
        word("default"),
        fileName_
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutUTabulatedWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
