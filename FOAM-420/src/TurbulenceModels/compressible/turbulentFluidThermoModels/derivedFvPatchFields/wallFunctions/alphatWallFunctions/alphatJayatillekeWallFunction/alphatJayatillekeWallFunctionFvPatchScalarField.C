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
    (c) 2010-2012 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "turbulentFluidThermoModels/derivedFvPatchFields/wallFunctions/alphatWallFunctions/alphatJayatillekeWallFunction/alphatJayatillekeWallFunctionFvPatchScalarField.H"
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

scalar alphatJayatillekeWallFunctionFvPatchScalarField::maxExp_ = 50.0;
scalar alphatJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Patch type for patch " << patch().name() << " must be wall\n"
            << "Current patch type is " << patch().type() << nl
            << exit(FatalError);
    }
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11.0;

    for (int i=0; i<maxIters_; i++)
    {
        scalar f = ypt - (log(E_*ypt)/kappa_ + P)/Prat;
        scalar df = 1.0 - 1.0/(ypt*kappa_*Prat);
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

alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85),
    Cmu_(0.09),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.lookupOrDefault<scalar>("Prt", 0.85)),
    Cmu_(dict.lookupOrDefault<scalar>("Cmu", 0.09)),
    kappa_(dict.lookupOrDefault<scalar>("kappa", 0.41)),
    E_(dict.lookupOrDefault<scalar>("E", 9.8))
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    Prt_(awfpsf.Prt_),
    Cmu_(awfpsf.Cmu_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    Prt_(awfpsf.Prt_),
    Cmu_(awfpsf.Cmu_),
    kappa_(awfpsf.kappa_),
    E_(awfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

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

    const tmp<scalarField> tmuw = turbModel.mu(patchi);
    const scalarField& muw = tmuw();

    const tmp<scalarField> tmutw = turbModel.mut(patchi);
    const scalarField& mutw = tmutw();

    const tmp<scalarField> talphaw = turbModel.alpha(patchi);
    const scalarField& alphaw = talphaw();

    scalarField& alphatw = *this;

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw( mag(Uw.snGrad()) );
    const scalarField magUp( magGradUw/this->patch().deltaCoeffs() );

    const scalarField& rhow = turbModel.rho().boundaryField()[patchi];
    const fvPatchScalarField& hew =
        turbModel.transport().he().boundaryField()[patchi];

    const scalarField& ry = patch().deltaCoeffs();

    // Heat flux [W/m2] - lagging alphatw
    const scalarField qDot
    (
        turbModel.transport().alphaEff(alphatw, patchi)*hew.snGrad()
    );

    // Populate boundary values
    forAll(alphatw, facei)
    {
        // Calculate uTau using Newton-Raphson iteration
        scalar uTau =
            sqrt((mutw[facei] + muw[facei])/rhow[facei]*magGradUw[facei]);

        if (uTau > ROOTVSMALL)
        {
            label iter = 0;
            scalar err = GREAT;

            do
            {
                scalar kUu = min(kappa_*magUp[facei]/uTau, maxExp_);
                scalar fkUu = exp(kUu) - 1.0 - kUu*(1.0 + 0.5*kUu);

                scalar f =
                    - uTau/(ry[facei]*muw[facei]/rhow[facei])
                    + magUp[facei]/uTau
                    + 1.0/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    - 1.0/(ry[facei]*muw[facei]/rhow[facei])
                    - magUp[facei]/sqr(uTau)
                    - 1.0/E_*kUu*fkUu/uTau;

                scalar uTauNew = uTau - f/df;
                err = mag((uTau - uTauNew)/uTau);
                uTau = uTauNew;

            } while (uTau>VSMALL && err>tolerance_ && ++iter<maxIters_);

            scalar yPlus = uTau/ry[facei]/(muw[facei]/rhow[facei]);

            // Molecular Prandtl number
            scalar Pr = muw[facei]/alphaw[facei];

            // Molecular-to-turbulenbt Prandtl number ratio
            scalar Prat = Pr/Prt_;

            // Thermal sublayer thickness
            scalar P = Psmooth(Prat);
            scalar yPlusTherm = this->yPlusTherm(P, Prat);

            // Evaluate new effective thermal diffusivity
            scalar alphaEff = 0.0;
            if (yPlus < yPlusTherm)
            {
                scalar A = qDot[facei]*rhow[facei]*uTau/ry[facei];
                scalar B = qDot[facei]*Pr*yPlus;
                scalar C = Pr*0.5*rhow[facei]*uTau*sqr(magUp[facei]);
                alphaEff = A/(B + C + VSMALL);
            }
            else
            {
                scalar A = qDot[facei]*rhow[facei]*uTau/ry[facei];
                scalar B = qDot[facei]*Prt_*(1.0/kappa_*log(E_*yPlus) + P);
                scalar magUc = uTau/kappa_*log(E_*yPlusTherm) - mag(Uw[facei]);
                scalar C =
                    0.5*rhow[facei]*uTau
                   *(Prt_*sqr(magUp[facei]) + (Pr - Prt_)*sqr(magUc));
                alphaEff = A/(B + C + VSMALL);
            }

            // Update turbulent thermal diffusivity
            alphatw[facei] = max(0.0, alphaEff - alphaw[facei]);

            if (debug)
            {
                Info<< "    uTau           = " << uTau << nl
                    << "    Pr             = " << Pr << nl
                    << "    Prt            = " << Prt_ << nl
                    << "    qDot           = " << qDot[facei] << nl
                    << "    yPlus          = " << yPlus << nl
                    << "    yPlusTherm     = " << yPlusTherm << nl
                    << "    alphaEff       = " << alphaEff << nl
                    << "    alphaw         = " << alphaw[facei] << nl
                    << "    alphatw      = " << alphatw[facei] << nl
                    << endl;
            }
        }
        else
        {
            alphatw[facei] = 0.0;
        }
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void alphatJayatillekeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("Prt", Prt_);
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("E", E_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
