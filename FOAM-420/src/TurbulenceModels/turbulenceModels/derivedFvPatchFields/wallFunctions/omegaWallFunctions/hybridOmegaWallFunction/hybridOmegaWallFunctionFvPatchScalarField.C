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
    (c) 2010-2012 Esi Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "hybridOmegaWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void hybridOmegaWallFunctionFvPatchScalarField::checkType()
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hybridOmegaWallFunctionFvPatchScalarField::hybridOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF),
    UName_("U"),
    kName_("k"),
    nuName_("nu"),
    nutName_("nut"),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8),
    tw_(0.057)
{
    checkType();
}


hybridOmegaWallFunctionFvPatchScalarField::hybridOmegaWallFunctionFvPatchScalarField
(
    const hybridOmegaWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedInternalValueFvPatchField<scalar>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    kName_(ptf.kName_),
    nuName_(ptf.nuName_),
    nutName_(ptf.nutName_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_),
    tw_(ptf.tw_)
{
    checkType();
}


hybridOmegaWallFunctionFvPatchScalarField::hybridOmegaWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedInternalValueFvPatchField<scalar>(p, iF, dict),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    kName_(dict.lookupOrDefault<word>("k", "k")),
    nuName_(dict.lookupOrDefault<word>("nu", "nu")),
    nutName_(dict.lookupOrDefault<word>("nut", "nut")),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8)),
    tw_(dict.lookupOrDefault<scalar>("tw", 0.057))
{
    checkType();
}


hybridOmegaWallFunctionFvPatchScalarField::hybridOmegaWallFunctionFvPatchScalarField
(
    const hybridOmegaWallFunctionFvPatchScalarField& owfpsf
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    tw_(owfpsf.tw_)
{
    checkType();
}


hybridOmegaWallFunctionFvPatchScalarField::hybridOmegaWallFunctionFvPatchScalarField
(
    const hybridOmegaWallFunctionFvPatchScalarField& owfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedInternalValueFvPatchField<scalar>(owfpsf, iF),
    UName_(owfpsf.UName_),
    kName_(owfpsf.kName_),
    nuName_(owfpsf.nuName_),
    nutName_(owfpsf.nutName_),
    Cmu_(owfpsf.Cmu_),
    kappa_(owfpsf.kappa_),
    E_(owfpsf.E_),
    tw_(owfpsf.tw_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hybridOmegaWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (!this->updated())
    {
        const fvMesh& mesh(patch().boundaryMesh().mesh());

        const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
             )
        );

        const scalar yPlusLam = turbModel.yPlusLam(kappa_, E_);
        const scalarField& y = turbModel.y()[patch().index()];

        const scalar Cmu25 = pow(Cmu_, 0.25);
        const scalar Cmu75 = pow(Cmu_, 0.75);

        typedef DimensionedField<scalar, volMesh> FieldType;
        FieldType& G =
            const_cast<FieldType&>
            (
                db().lookupObject<FieldType>
                (
                    turbModel.GName()
                 )
            );

        FieldType& omega = const_cast<FieldType&>(internalField());

        const scalarField& k = db().lookupObject<volScalarField>(kName_);

        const scalarField& nuw =
            patch().lookupPatchField<volScalarField, scalar>(nuName_);

        const tmp<scalarField> tnutw = turbModel.nut(this->patch().index());
        const scalarField& nutw = tnutw();

        const fvPatchVectorField& Uw =
            patch().lookupPatchField<volVectorField, vector>(UName_);

        const scalarField magGradUw
        (
            mag(Uw.snGrad() - patch().nf() * (Uw.snGrad() & patch().nf()))
        );

        // Set omega and G
        forAll(nutw, faceI)
        {
            label faceCellI = patch().faceCells()[faceI];

            //yPlus from wall shear to allow for laminar region
            scalar tauw = (nutw[faceI]+nuw[faceI])*magGradUw[faceI];
            scalar utau = max(SMALL,::sqrt(tauw));

            scalar yPlus = utau*y[faceI]/nuw[faceI];

            //star hybrid formulation for k production and omega

            //phi blending
            scalar phiHL = sqr(1.0-exp(-yPlus/yPlusLam));

             //inverse apparent cell thickness
            scalar sv = patch().magSf()[faceI]/mesh.V()[faceCellI];

            if (yPlus == 0)
            {
                G[faceCellI] = 0;
            }
            else
            {
                G[faceCellI] =
                    phiHL *
                    (
                        tauw*magGradUw[faceI]*y[faceI]*sv
                        -Cmu75 * sv/kappa_*log(E_*yPlus)
                        *::sqrt(k[faceCellI]) * k[faceCellI]
                    )
                    - (1.0-phiHL) *
                    (
                        4.0*nuw[faceI]*sqr(sv)*k[faceCellI]
                    );
            }

            //STAR continuous WF formulation for Omega
            scalar Rek = y[faceI]*::sqrt(k[faceCellI])/nuw[faceI];
            scalar Lomega = y[faceI]*(1.0-exp(-tw_*Rek));
            omega[faceCellI] = ::sqrt(k[faceCellI])
                /max((kappa_*Cmu25*Lomega), SMALL);
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void hybridOmegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "k", "k", kName_);
    writeEntryIfDifferent<word>(os, "nu", "nu", nuName_);
    writeEntryIfDifferent<word>(os, "nut", "nut", nutName_);
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("E", E_);
    os.writeEntry("tw", tw_);

    fixedInternalValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    hybridOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
