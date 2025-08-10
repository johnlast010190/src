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

#include "DtJayatillekeWallFunctionFvPatchScalarField.H"
#include "incompressibleTurbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar DtJayatillekeWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar DtJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label DtJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void DtJayatillekeWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


scalar DtJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Scat
) const
{
    return 9.24*(pow(Scat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Scat));
}


scalar DtJayatillekeWallFunctionFvPatchScalarField::yPlusSpec
(
    const scalar P,
    const scalar Scat
) const
{
    scalar ypt = 11.0;

    for (int i=0; i<maxIters_; i++)
    {
        scalar f = ypt - (log(E_*ypt)/kappa_ + P)/Scat;
        scalar df = 1.0 - 1.0/(ypt*kappa_*Scat);
        scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance_)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
    }

    return ypt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

DtJayatillekeWallFunctionFvPatchScalarField::
DtJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


DtJayatillekeWallFunctionFvPatchScalarField::
DtJayatillekeWallFunctionFvPatchScalarField
(
    const DtJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


DtJayatillekeWallFunctionFvPatchScalarField::
DtJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{
    checkType();
}


DtJayatillekeWallFunctionFvPatchScalarField::
DtJayatillekeWallFunctionFvPatchScalarField
(
    const DtJayatillekeWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


DtJayatillekeWallFunctionFvPatchScalarField::
DtJayatillekeWallFunctionFvPatchScalarField
(
    const DtJayatillekeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DtJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const IOdictionary& transportProperties
        = db().lookupObject<IOdictionary>("transportProperties");

    string pofix(string(this->internalField().name())(2, 100));

    //read laminar mass diffusivity
    dimensionedScalar D
    (
        "D" + pofix,
        transportProperties.lookup("D" + pofix)
    );

    //read turbulent Schmidt number
    dimensionedScalar Sct
    (
        "Sct" + pofix,
        transportProperties.lookup("Sct" + pofix)
    );

    const incompressibleTurbulenceModel& turbModel =
        db().lookupObject<incompressibleTurbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );

    // Field data
    const label patchI = patch().index();

    const scalarField rhow
    (
        turbModel.rho()().boundaryField()[patchI]
    );

    tmp<scalarField> yPlus(turbModel.yPlus(patchI));

    const scalarField& nuw = turbModel.nu()().boundaryField()[patchI];
    scalarField Sc( nuw/D.value() );

    const scalarField nutw = turbModel.nut()().boundaryField()[patchI];

    scalarField& Dtw = *this;

    const scalarField& ry = patch().deltaCoeffs();

    // Populate boundary values
    forAll(Dtw, faceI)
    {
        scalar uTau = yPlus()[faceI] * ry[faceI] * nuw[faceI];

        if (uTau > ROOTVSMALL)
        {

            scalar yPlusf = yPlus()[faceI];

            // Molecular-to-turbulenbt Prandtl number ratio
            scalar Scat = Sc[faceI]/Sct.value();

            // Species sublayer thickness
            scalar P = Psmooth(Scat);
            scalar yPlusSpec = this->yPlusSpec(P, Scat);

            // Evaluate new effective species diffusivity
            scalar Deff = 0.0;
            if (yPlusf < yPlusSpec)
            {
                Deff = D.value();
            }
            else
            {
                Deff = (rhow[faceI]*uTau/ry[faceI])
                    / (Sct.value() * (1.0/kappa_*log(E_*yPlusf) + P));
            }

            // Update turbulent species diffusivity
            Dtw[faceI] = max(0.0, Deff - D.value());

        }
        else
        {
            Dtw[faceI] = 0.0;
        }
    }

    fixedValueFvPatchScalarField::updateCoeffs();
}


void DtJayatillekeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("E", E_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    DtJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
