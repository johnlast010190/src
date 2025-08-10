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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/DtWallFunctions/DtJayatillekeWallFunction/DtJayatillekeWallFunctionFvPatchScalarField.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar DtJayatillekeWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar DtJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label DtJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

scalar DtJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar DtJayatillekeWallFunctionFvPatchScalarField::yPlusMass
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
{}


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
{}


DtJayatillekeWallFunctionFvPatchScalarField::
DtJayatillekeWallFunctionFvPatchScalarField
(
    const DtJayatillekeWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{}


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
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void DtJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
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

    // Field data
    const label patchI = patch().index();

    tmp<scalarField> tmuw = turbModel.mu(patchI);
    const scalarField& muw = tmuw();
    const tmp<scalarField> tmutw = turbModel.mut(patchI);
    const scalarField& mutw = tmutw();

    scalarField& Dtw = *this;

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    const scalarField magGradUw( mag(Uw.snGrad()) );
    const scalarField magUp( magGradUw/this->patch().deltaCoeffs() );

    const scalarField& rhow = turbModel.rho().boundaryField()[patchI];
    const fvPatchScalarField& speciesw =
        patch().lookupPatchField<volScalarField, scalar>("pofix");

    const scalarField& ry = patch().deltaCoeffs();

    // Mass flux [W/m2] - lagging Dtw
    const scalarField mDot( (Dtw + rhow*D.value())*speciesw.snGrad() );

    // Populate boundary values
    forAll(Dtw, faceI)
    {
        // Calculate uTau using Newton-Raphson iteration
        scalar uTau =
            sqrt((mutw[faceI] + muw[faceI])/rhow[faceI]*magGradUw[faceI]);

        if (uTau > ROOTVSMALL)
        {
            label iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUp[faceI]/uTau, maxExp_);
                scalar fkUu = exp(kUu) - 1.0 - kUu*(1.0 + 0.5*kUu);

                scalar f =
                    - uTau/(ry[faceI]*muw[faceI]/rhow[faceI])
                    + magUp[faceI]/uTau
                    + 1.0/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    - 1.0/(ry[faceI]*muw[faceI]/rhow[faceI])
                    - magUp[faceI]/sqr(uTau)
                    - 1.0/E_*kUu*fkUu/uTau;

                scalar uTauNew = uTau - f/df;
                err = mag((uTau - uTauNew)/uTau);
                uTau = uTauNew;

            } while (uTau>VSMALL && err>tolerance_ && ++iter<maxIters_);

            scalar yPlus = uTau/ry[faceI]/(muw[faceI]/rhow[faceI]);

            // Molecular Schmidt number
            scalar Sc = muw[faceI]/(rhow[faceI]*D.value());

            // Molecular-to-turbulenbt Schmidt number ratio
            scalar Scat = Sc/Sct.value();

            // Mass diffusivity sublayer thickness
            scalar P = Psmooth(Scat);
            scalar yPlusMass = this->yPlusMass(P, Scat);

            // Evaluate new effective mass diffusivity
            scalar Deff = 0.0;
            if (yPlus < yPlusMass)
            {
                scalar A = mDot[faceI]*rhow[faceI]*uTau/ry[faceI];
                scalar B = mDot[faceI]*Sc*yPlus;
                scalar C = Sc*0.5*rhow[faceI]*uTau*sqr(magUp[faceI]);
                Deff = A/(B + C + VSMALL);
            }
            else
            {
                scalar A = mDot[faceI]*rhow[faceI]*uTau/ry[faceI];
                scalar B = mDot[faceI]*Sct.value()*(1.0/kappa_*log(E_*yPlus) + P);
                scalar magUc = uTau/kappa_*log(E_*yPlusMass) - mag(Uw[faceI]);
                scalar C =
                    0.5*rhow[faceI]*uTau
                   *(Sct.value()*sqr(magUp[faceI])
                   + (Sc - Sct.value())*sqr(magUc));
                Deff = A/(B + C + VSMALL);
            }

            // Update turbulent mass diffusivity
            Dtw[faceI] = max(0.0, Deff - D.value());

            if (debug)
            {
                Info<< "    uTau           = " << uTau << nl
                    << "    Sc             = " << Sc << nl
                    << "    Sct            = " << Sct.value() << nl
                    << "    mDot           = " << mDot[faceI] << nl
                    << "    yPlus          = " << yPlus << nl
                    << "    yPlusMass      = " << yPlusMass << nl
                    << "    Deff           = " << Deff << nl
                    << "    D              = " << D.value() << nl
                    << "    Dtw            = " << Dtw[faceI] << nl
                    << endl;
            }
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

// New name
addNamedToPatchFieldRunTimeSelection
(
    fvPatchScalarField,
    DtJayatillekeWallFunctionFvPatchScalarField,
    DtJayatillekeWallFunction
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
