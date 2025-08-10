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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nutUKnopSAWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar nutUKnopSAWallFunctionFvPatchScalarField::calcYPlus
(
    const scalar y,
    const scalar utau,
    const scalar nu
)  const
{
    return (y*utau/nu);
}


tmp<scalarField> nutUKnopSAWallFunctionFvPatchScalarField::uTauReichardt
(
    const scalarField& magGradU,
    const scalarField& utauInit
) const
{
    // stuff needed everywhere
    const label patchI = patch().index();
    vectorField nf(patch().nf()());

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchI];

    const tmp<scalarField> tnuw = turbModel.nu(patchI);
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    scalarField magUp(mag(Uw.snGrad()/this->patch().deltaCoeffs()));
    magUp = max(magUp, ROOTVSMALL);

    tmp<scalarField> tuTau(new scalarField(utauInit));
    scalarField& uTau = tuTau.ref();

    label maxIter(20);
    scalar maxErr(1e-3);

    forAll(uTau, fi)
    {
        scalar& ut(uTau[fi]);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar uplus(magUp[fi]/ut);
                scalar duplus(-magUp[fi]/sqr(ut));

                scalar yp(calcYPlus(y[fi], ut, nuw[fi]));
                scalar dyp(y[fi]/nuw[fi]);

                // log law blending factor
                scalar argBl(yp/CyR_);
                scalar phiBl = argBl*argBl*argBl*argBl; //pow4
                phiBl = tanh(phiBl);

                // Reichardt
                scalar Frei = log(1+Cr1_*yp)/kappa_
                    + Cr2_*(1-exp(-yp/Cr3_) - yp/Cr3_*exp(-yp/Cr4_));

                // Std log law
                scalar Flog = log(yp)/kappa_ + Br1_;

                // Combined functional
                scalar Freim = (1-phiBl)*Frei + phiBl*Flog - uplus;

                // Derivative of Freim

                // dFreim/dut = d/dut(1-phiBl)*Frei + (1-phiBl)*d/dut(Frei)
                //              d/dut(phiBl)*Flog + phiBl*d/dut(Flog)

                scalar dphiBl = (1 - sqr(phiBl))*4*(argBl*argBl*argBl)*dyp/CyR_;

                scalar dFrei = 1/(1+Cr1_*yp)/kappa_*(Cr1_*dyp)
                             + Cr2_*(-exp(-yp/Cr3_))*(-dyp/Cr3_)
                             + Cr2_*(dyp/Cr3_)*exp(-yp/Cr4_)
                             + Cr2_*(yp/Cr3_)*exp(-yp/Cr4_)*(-dyp/Cr4_);

                scalar dFlog = 1/yp/kappa_*dyp;

                scalar dFreim = -dphiBl*Frei + (1-phiBl)*dFrei
                                + dphiBl*Flog + phiBl*dFlog - duplus;

                scalar uto = max(ut, ROOTVSMALL);
                ut -= Freim/dFreim;
                if (ut < ROOTVSMALL) ut = 0.5*uto;

                err = mag((ut - uto)/uto);

            } while (ut > ROOTVSMALL && err > maxErr && ++iter < maxIter);

            if (iter >= maxIter)
            {
                WarningInFunction
                    << "Knop SA wall function (Reichardt section) did not "
                    << "converge within max iteration limit:" << maxIter << nl
                    << tab << "patch: " << this->patch().name()
                    << ", face: " << fi << ", error[%]: " << err
                    << ", utau: " << ut << ", y+: "<< calcYPlus(y[fi], ut, nuw[fi])
                    << endl;
            }

        }
    }

    return tuTau;
}

tmp<scalarField> nutUKnopSAWallFunctionFvPatchScalarField::uTauSpaldingN
(
    const scalarField& magGradU,
    const label N,
    const scalarField& utauInit
) const
{

    // stuff needed everywhere
    const label patchI = patch().index();
    vectorField nf(patch().nf()());

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchI];
    const tmp<scalarField> tnuw = turbModel.nu(patchI);
    const scalarField& nuw = tnuw();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    scalarField magUp(mag(Uw.snGrad()/this->patch().deltaCoeffs()));
    magUp = max(magUp, ROOTVSMALL);

    tmp<scalarField> tuTau(new scalarField(utauInit));
    scalarField& uTau = tuTau.ref();

    label maxIter(10);
    scalar maxErr(1e-3);
    scalar FsnC(exp(-kappa_*Bsa1_));
    label N1 = N+1;

    forAll(uTau, fi)
    {
        scalar& ut(uTau[fi]);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err = GREAT;

            do
            {
                scalar uplus(magUp[fi]/ut);
                scalar duplus(-magUp[fi]/sqr(ut));

                scalar kUu = kappa_*uplus;
                scalar dkUu = kappa_*duplus;

                scalar yp(calcYPlus(y[fi], ut, nuw[fi]));
                scalar dyp(y[fi]/nuw[fi]);

                scalar Fsp = -1;
                scalar dFsp = 0;

                for (label n=1; n<N1; n++)
                {
                    Fsp -= pow(kUu, n)/scalar(factorial(n));
                    dFsp -= n*pow(kUu, n-1)/scalar(factorial(n)) * dkUu;
                }

                Fsp += exp(kUu);
                Fsp *= FsnC;
                Fsp += uplus - yp;

                dFsp += exp(kUu)*dkUu;
                dFsp *= FsnC;
                dFsp += duplus - dyp;

                scalar uto = max(ut, ROOTVSMALL);
                ut -= Fsp/dFsp;

                err = mag((ut - uto)/uto);

            } while (ut > ROOTVSMALL && err > maxErr && ++iter < maxIter);

            if (iter >= maxIter)
            {
                WarningInFunction
                    << "Knop SA wall function (Spalding section) did not "
                    << "converge within max iteration limit:" << maxIter << nl
                    << tab << "patch: " << this->patch().name()
                    << ", face: " << fi << ", error[%]: " << err
                    << ", utau: " << ut << ", y+: "<< calcYPlus(y[fi], ut, nuw[fi])
                    << endl;
            }

        }
    }

    return tuTau;
}

tmp<scalarField> nutUKnopSAWallFunctionFvPatchScalarField::calcNut() const
{
    // stuff needed everywhere
    const label patchI = patch().index();
    vectorField nf(patch().nf()());

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchI];
    scalarField magGradU(this->size(), 0.0);
    {
        vectorField gradU(Uw.snGrad());
        gradU -= nf*(gradU & nf);
        magGradU = mag(gradU);
        magGradU = max(magGradU, ROOTVSMALL);
    }

    const tmp<scalarField> tnuw = turbModel.nu(patchI);
    const scalarField& nuw = tnuw();

    return max
    (
        scalar(0),
        sqr(calcUTau(magGradU))/magGradU - nuw
    );
}


tmp<scalarField> nutUKnopSAWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    // stuff needed everywhere
    const label patchI = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchI];
    const tmp<scalarField> tnuw = turbModel.nu(patchI);
    const scalarField& nuw = tnuw();
    const scalarField& nutw = *this;

    //scalarField utauInit(magGradU);
    scalarField utauInit(sqrt((nutw + nuw)*magGradU));

    // calc Reichardt utau
    tmp<scalarField> utauR(uTauReichardt(magGradU, utauInit));


    // calc Spalding N utau
    utauInit = utauR();
    tmp<scalarField> utauSN(uTauSpaldingN(magGradU, 5, utauInit));


    // iterate to get y+ dependent blending function
    tmp<scalarField> tuTau(new scalarField(utauInit));
    scalarField& uTau(tuTau.ref());

    label maxIter(10);
    scalar maxErr(1e-3);

    forAll(uTau, fi)
    {
        int iter = 0;
        scalar ypErr(GREAT);
        scalar yp(calcYPlus(y[fi], uTau[fi], nuw[fi]));


        do
        {
            scalar phiSA(yp/CySA_);
            phiSA *= phiSA*phiSA;
            phiSA = tanh(phiSA);

            uTau[fi] = (1-phiSA)*utauSN()[fi] + phiSA*utauR()[fi];

            scalar ypold(max(yp, ROOTVSMALL));
            yp = calcYPlus(y[fi], uTau[fi], nuw[fi]);
            ypErr = mag(yp-ypold)/ypold;

        } while (ypErr > maxErr && ++iter < maxIter);

        if (iter >= maxIter)
        {
            WarningInFunction
                << "Knop SA wall function (Blending section) did not "
                << "converge within max iteration limit:" << maxIter << nl
                << tab << "patch: " << this->patch().name()
                << ", face: " << fi << ", y+ error[%]: " << ypErr
                << ", utau: " << uTau[fi] << ", y+: " << yp
                << endl;
        }
    }

    return tuTau;
}

void nutUKnopSAWallFunctionFvPatchScalarField::
writeLocalEntries(Ostream& os) const
{
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("E", E_);

    os.writeEntry("Cr1", Cr1_);
    os.writeEntry("Cr2", Cr2_);
    os.writeEntry("Cr3", Cr3_);
    os.writeEntry("Cr4", Cr4_);

    os.writeEntry("Br1", Br1_);
    os.writeEntry("CyR", CyR_);
    os.writeEntry("Bsa1", Bsa1_);
    os.writeEntry("CySA", CySA_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUKnopSAWallFunctionFvPatchScalarField::
nutUKnopSAWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    Cr1_(0.4),
    Cr2_(7.8),
    Cr3_(11.0),
    Cr4_(3.0),
    Br1_(5.1),
    CyR_(27.0),
    Bsa1_(5.2),
    CySA_(24.0)
{}


nutUKnopSAWallFunctionFvPatchScalarField::
nutUKnopSAWallFunctionFvPatchScalarField
(
    const nutUKnopSAWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    Cr1_(ptf.Cr1_),
    Cr2_(ptf.Cr2_),
    Cr3_(ptf.Cr3_),
    Cr4_(ptf.Cr4_),
    Br1_(ptf.Br1_),
    CyR_(ptf.CyR_),
    Bsa1_(ptf.Bsa1_),
    CySA_(ptf.CySA_)
{}


nutUKnopSAWallFunctionFvPatchScalarField::
nutUKnopSAWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    Cr1_(dict.lookupOrDefault<scalar>("Cr1", 0.4)),
    Cr2_(dict.lookupOrDefault<scalar>("Cr2", 7.8)),
    Cr3_(dict.lookupOrDefault<scalar>("Cr3", 11.0)),
    Cr4_(dict.lookupOrDefault<scalar>("Cr4", 3.0)),
    Br1_(dict.lookupOrDefault<scalar>("Br1", 5.1)),
    CyR_(dict.lookupOrDefault<scalar>("CyR", 27.0)),
    Bsa1_(dict.lookupOrDefault<scalar>("Bsa1", 5.2)),
    CySA_(dict.lookupOrDefault<scalar>("CySA", 24.0))
{}


nutUKnopSAWallFunctionFvPatchScalarField::
nutUKnopSAWallFunctionFvPatchScalarField
(
    const nutUKnopSAWallFunctionFvPatchScalarField& ptf
)
:
    nutWallFunctionFvPatchScalarField(ptf),
    Cr1_(ptf.Cr1_),
    Cr2_(ptf.Cr2_),
    Cr3_(ptf.Cr3_),
    Cr4_(ptf.Cr4_),
    Br1_(ptf.Br1_),
    CyR_(ptf.CyR_),
    Bsa1_(ptf.Bsa1_),
    CySA_(ptf.CySA_)
{}


nutUKnopSAWallFunctionFvPatchScalarField::
nutUKnopSAWallFunctionFvPatchScalarField
(
    const nutUKnopSAWallFunctionFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(ptf, iF),
    Cr1_(ptf.Cr1_),
    Cr2_(ptf.Cr2_),
    Cr3_(ptf.Cr3_),
    Cr4_(ptf.Cr4_),
    Br1_(ptf.Br1_),
    CyR_(ptf.CyR_),
    Bsa1_(ptf.Bsa1_),
    CySA_(ptf.CySA_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutUKnopSAWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();
    vectorField nf(patch().nf()());

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    scalarField magGradU(this->size(), 0.0);
    {
        vectorField gradU(Uw.snGrad());
        gradU -= nf*(gradU & nf);
        magGradU = mag(gradU);
        magGradU = max(magGradU, ROOTVSMALL);
    }

    return y * sqrt((nuw + *this) * magGradU) / nuw;
}


void nutUKnopSAWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutUKnopSAWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
