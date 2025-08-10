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

#include "nutUKnopKOmegaWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

scalar nutUKnopKOmegaWallFunctionFvPatchScalarField::calcYPlus
(
    scalar y,
    scalar utau,
    scalar nu
) const
{
    return max(min((y*utau/nu), GREAT), ROOTVSMALL);
}


tmp<scalarField> nutUKnopKOmegaWallFunctionFvPatchScalarField::uTauReichardt
(
    const scalarField& magGradU,
    const scalarField& utauInit
) const
{
    // stuff needed everywhere
    const label patchI = patch().index();
    vectorField nf(patch().nf()());

    const auto& turbModel = db().lookupObject<turbulenceModel>
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

    bool printWarningMessage = false;

    forAll(uTau, fi)
    {
        scalar roughnessHeight = this->roughWallCoeffs().roughnessHeight();
        scalar& ut(uTau[fi]);

        // Limit initial solution for ut to cope with tricky flow initilizations
        ut = max(1e-06, ut);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err;

            do
            {
                scalar uplus(magUp[fi]/ut);
                scalar duplus(-magUp[fi]/sqr(ut));

                scalar yp(calcYPlus(y[fi], ut, nuw[fi]));
                scalar dyp(y[fi]/nuw[fi]);

                scalar coeffLim = 0.5;
                roughnessHeight = min(roughnessHeight, y[fi]/0.3*exp(kappa_*(Br1_-coeffLim)));
                scalar KsPlus = ut*roughnessHeight/nuw[fi];
                scalar DU = 1./kappa_*log(1+0.3*KsPlus);
                scalar dKsPlus = 0.3*roughnessHeight/(kappa_*(1.+0.3*KsPlus))/nuw[fi];

                // log law blending factor
                scalar argBl(yp/CyR_);
                scalar phiBl = pow4(argBl);
                phiBl = tanh(phiBl);

                // Reichardt
                scalar Frei = log(1+Cr1_*yp)/kappa_
                    + Cr2_*(1-exp(-yp/Cr3_) - yp/Cr3_*exp(-yp/Cr4_))-DU;
                // Std log law
                scalar Flog = log(yp)/kappa_ + Br1_ - DU;
                // Combined functional
                scalar Freim = (1-phiBl)*Frei + phiBl*Flog - uplus;

                // Derivative of Freim

                // dFreim/dut = d/dut(1-phiBl)*Frei + (1-phiBl)*d/dut(Frei)
                //              d/dut(phiBl)*Flog + phiBl*d/dut(Flog)

                scalar dphiBl = (1 - sqr(phiBl))*4*(argBl*argBl*argBl)*dyp/CyR_;

                scalar dFrei = 1/(1+Cr1_*yp)/kappa_*(Cr1_*dyp)
                             + Cr2_*(-exp(-yp/Cr3_))*(-dyp/Cr3_)
                             + Cr2_*(dyp/Cr3_)*exp(-yp/Cr4_)
                             + Cr2_*(yp/Cr3_)*exp(-yp/Cr4_)*(-dyp/Cr4_)-dKsPlus;

                scalar dFlog = 1/yp/kappa_*dyp - dKsPlus;
                scalar dFreim = -dphiBl*Frei + (1-phiBl)*dFrei
                                + dphiBl*Flog + phiBl*dFlog - duplus;

                scalar uto = max(ut, ROOTVSMALL);
                ut -= Freim/dFreim;
                if (ut < ROOTVSMALL) ut = 0.5*uto;

                err = mag((ut - uto)/uto);

            } while (ut > ROOTVSMALL && err > maxErr && ++iter < maxIter);

            if (iter >= maxIter)
            {
                printWarningMessage = true;
                if (debug)
                {
                    WarningInFunction
                    << "Knop kOmega wall function (Reichardt section) did not "
                    << "converge within max iteration limit:" << maxIter << nl
                    << tab << "patch: " << this->patch().name()
                    << ", face: " << fi << ", error[%]: " << err
                    << ", utau: " << ut << ", y+: "<< calcYPlus(y[fi], ut, nuw[fi])
                    << endl;
                }
            }
        }
    }

    if (printWarningMessage)
    {
        WarningInFunction
            << "Knop kOmega wall function (Reichardt section) did not "
            << "converge within max iteration limit:" << maxIter << endl;
    }

    return tuTau;
}

tmp<scalarField> nutUKnopKOmegaWallFunctionFvPatchScalarField::uTauSpaldingN
(
    const scalarField& magGradU,
    label N,
    const scalarField& utauInit
) const
{

    // stuff needed everywhere
    const label patchI = patch().index();
    vectorField nf(patch().nf()());

    const auto& turbModel = db().lookupObject<turbulenceModel>
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
    scalar roughnessHeight = roughWallCoeffs().roughnessHeight();
    if (roughnessHeight>0)
    {
        //With roughness correction it may take longer to converge
        maxIter = 25;
    }

    scalar maxErr(1e-3);
    scalar FsnC(exp(-kappa_*Bsa1_));
    label N1 = N+1;
    bool printWarningMessage = false;
    forAll(uTau, fi)
    {
        scalar& ut(uTau[fi]);

        if (ut > ROOTVSMALL)
        {
            int iter = 0;
            scalar err;

            do
            {
                scalar uplus = min((magUp[fi]/ut), 50.);
                scalar duplus(-magUp[fi]/sqr(ut));
                scalar KsPlus = ut*roughnessHeight/nuw[fi];

                scalar DU = 1./kappa_*log(1+0.3*KsPlus);
                scalar dKsPlus = 0.3*roughnessHeight/(kappa_*(1.+0.3*KsPlus))/nuw[fi];

                scalar kUu = kappa_*( uplus-DU);
                scalar dkUu = kappa_*(duplus-dKsPlus);

                scalar yp(calcYPlus(y[fi], ut, nuw[fi]));
                scalar dyp(y[fi]/nuw[fi]);

                scalar Fsp = -1;
                scalar dFsp = 0;

                for (label n=1; n<N1; n++)
                {
                    Fsp -= pow(kUu, n)/scalar(factorial(n));
                    dFsp -= n*pow(kUu, n-1)/scalar(factorial(n)) * dkUu;
                }

                Fsp += exp(kappa_*uplus);
                Fsp *= FsnC;
                Fsp += uplus - yp;

                dFsp += exp(kappa_*uplus)*kappa_*duplus;
                dFsp *= FsnC;
                dFsp += duplus - dyp;

                scalar uto = max(ut, ROOTVSMALL);
                ut = (1-newtonsRelax_)*uto + newtonsRelax_*(uto - Fsp/dFsp);

                err = mag((ut - uto)/uto);

            } while (ut > ROOTVSMALL && err > maxErr && ++iter < maxIter);

            if (iter >= maxIter)
            {
                printWarningMessage = true;
                if (debug)
                {
                    WarningInFunction
                    << "Knop kOmega wall function (Spalding section) did not "
                    << "converge within max iteration limit:" << maxIter << nl
                    << tab << "patch: " << this->patch().name()
                    << ", face: " << fi << ", error[%]: " << err
                    << ", utau: " << ut << ", y+: "<< calcYPlus(y[fi], ut, nuw[fi])
                    << endl;
                }
            }

        }
        if (uTau[fi]<ROOTVSMALL && roughnessHeight>0)
        {
            uTau[fi] = utauInit[fi];
        }
    }

    if (printWarningMessage)
    {
        WarningInFunction
            << "Knop kOmega wall function (Spalding section) did not "
            << "converge within max iteration limit:" << maxIter << endl;
    }

    return tuTau;
}

tmp<scalarField> nutUKnopKOmegaWallFunctionFvPatchScalarField::calcNut() const
{
    // stuff needed everywhere
    const label patchI = patch().index();
    vectorField nf(patch().nf()());

    const auto& turbModel = db().lookupObject<turbulenceModel>
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

    tmp<scalarField> tnutCurrent(sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw);
    scalarField& nutCurrent = tnutCurrent.ref();
    tmp<scalarField> tnutw(new scalarField(*this));
    scalarField& nutw = tnutw.ref();

    if (roughWallCoeffs().roughnessIsActive())
    {
        forAll(nutCurrent, fi)
        {
            nutCurrent[fi] = max
            (
                min
                (
                    nutCurrent[fi],
                    max(2*nutw[fi], nuw[fi])
                ),
                0.5*nutw[fi]
            );
        }
    }

    return max
    (
        scalar(0),
        tnutCurrent
    );
}


tmp<scalarField> nutUKnopKOmegaWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    // stuff needed everywhere
    const label patchI = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
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
    label N1 = 3;
    tmp<scalarField> utauSN(uTauSpaldingN(magGradU, N1, utauInit));


    // iterate to get y+ dependent blending function
    tmp<scalarField> tuTau(new scalarField(utauInit));
    scalarField& uTau(tuTau.ref());

    label maxIter(10);
    scalar roughnessHeight = this->roughWallCoeffs().roughnessHeight();
    if (roughnessHeight>0)
    {
        maxIter=20;
    }
    scalar maxErr(3e-3);

    bool printWarningMessage = false;

    forAll(uTau, fi)
    {
        int iter = 0;
        scalar ypErr;
        scalar yp(calcYPlus(y[fi], uTau[fi], nuw[fi]));

        do
        {
            scalar phiKO(sqr(yp/CyKOmega_));
            phiKO = tanh(phiKO);

            uTau[fi] = (1-phiKO)*utauSN()[fi] + phiKO*utauR()[fi];

            scalar ypold(max(yp, ROOTVSMALL));
            yp = calcYPlus(y[fi], uTau[fi], nuw[fi]);
            ypErr = mag(yp-ypold)/ypold;

        } while (ypErr > maxErr && ++iter < maxIter);

        if (iter >= maxIter)
        {
            if (debug)
            {
                WarningInFunction
                << "Knop kOmega wall function (Blending section) did not "
                << "converge within max iteration limit:" << maxIter << nl
                << tab << "patch: " << this->patch().name()
                << ", face: " << fi << ", y+ error[%]: " << ypErr
                << ", utau: " << uTau[fi] << ", y+: " << yp
                << endl;
            }
            printWarningMessage = true;
        }
    }

    if (printWarningMessage)
    {
        WarningInFunction
            << "Knop kOmega wall function (Blending section) did not "
            << "converge within max iteration limit:" << maxIter << nl
            << endl;
    }

    return tuTau;
}

void nutUKnopKOmegaWallFunctionFvPatchScalarField::
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
    os.writeEntry("CyKOmega", CyKOmega_);
    os.writeEntry("NewtonsRelax", newtonsRelax_);
    roughWallCoeffs_.writeLocalEntries(os);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

nutUKnopKOmegaWallFunctionFvPatchScalarField::
nutUKnopKOmegaWallFunctionFvPatchScalarField
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
    CyKOmega_(50.0),
    newtonsRelax_(0.9)
{}


nutUKnopKOmegaWallFunctionFvPatchScalarField::
nutUKnopKOmegaWallFunctionFvPatchScalarField
(
    const nutUKnopKOmegaWallFunctionFvPatchScalarField& ptf,
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
    CyKOmega_(ptf.CyKOmega_),
    newtonsRelax_(ptf.newtonsRelax_)
{}


nutUKnopKOmegaWallFunctionFvPatchScalarField::
nutUKnopKOmegaWallFunctionFvPatchScalarField
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
    CyKOmega_(dict.lookupOrDefault<scalar>("CyKOmega", 50.0)),
    newtonsRelax_(dict.lookupOrDefault<scalar>("NewtonsRelax", 0.9))
{}


nutUKnopKOmegaWallFunctionFvPatchScalarField::
nutUKnopKOmegaWallFunctionFvPatchScalarField
(
    const nutUKnopKOmegaWallFunctionFvPatchScalarField& ptf
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
    CyKOmega_(ptf.CyKOmega_),
    newtonsRelax_(ptf.newtonsRelax_)
{}


nutUKnopKOmegaWallFunctionFvPatchScalarField::
nutUKnopKOmegaWallFunctionFvPatchScalarField
(
    const nutUKnopKOmegaWallFunctionFvPatchScalarField& ptf,
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
    CyKOmega_(ptf.CyKOmega_),
    newtonsRelax_(ptf.newtonsRelax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> nutUKnopKOmegaWallFunctionFvPatchScalarField::yPlus() const
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
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    scalarField magGradU(this->size(), 0.0);
    {
        vectorField gradU(Uw.snGrad());
        gradU -= nf*(gradU & nf);
        magGradU = mag(gradU);
        magGradU = max(magGradU, SMALL);
    }

    return min((y*sqrt((nuw + *this) * magGradU)/nuw), GREAT);
}


void nutUKnopKOmegaWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    nutUKnopKOmegaWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
