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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "IAPWSIF97EOS.H"
#include "primitives/strings/fileName/fileName.H"
#include "db/IOobject/IOobject.H"

#include "global/unitConversion/unitConversion.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary Foam::IAPWSIF97EOS::loadDataDict
(
    const objectRegistry& obr
) const
{
    fileName etcFileName = "thermoData/IAPWS-IF97EOS";
    fileName waterSteamTablesFileName = findEtcFile(etcFileName);
    waterSteamTablesFileName.expand();
    waterSteamTablesFileName.toAbsolute();

    Info<< "Reading IAPWS-IF97EOS from: "
        << waterSteamTablesFileName << nl << endl;

    IOdictionary readDict
    (
        IOobject
        (
            waterSteamTablesFileName,
            obr,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );
    return readDict;
}


Foam::Field<Foam::scalarField> Foam::IAPWSIF97EOS::loadLambda1Coeffs
(
    const dictionary& dict
) const
{
    Field<scalarField> lambda1n(6);
    wordList nCoeffsNames({"n1", "n2", "n3", "n4", "n5", "n6"});
    forAll(lambda1n, ni)
    {
        lambda1n[ni] = kappaDict_.subDict("LAMBDA1").lookup<scalarField>(nCoeffsNames[ni]);
    }
    return lambda1n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IAPWSIF97EOS::IAPWSIF97EOS(const objectRegistry& obr)
:
    dict_(loadDataDict(obr)),
    n_(dict_.subDict("saturationCurve").lookup<scalarField>("n")),
    n23_(dict_.subDict("boundary23").lookup<scalarField>("n")),
    region1_(dict_.subDict("GibbsFreeEnergyRegion1")),
    region2_(dict_.subDict("GibbsFreeEnergyRegion2")),
    region3_(dict_.subDict("GibbsFreeEnergyRegion3")),
    region5_(dict_.subDict("GibbsFreeEnergyRegion5")),
    muDict_(dict_.subDict("muModelCoeffs")),
    muTmax1_(muDict_.lookup<scalar>("Tmax1")),
    muTmax2_(muDict_.lookup<scalar>("Tmax2")),
    muN0_(muDict_.lookup<scalarField>("n0")),
    muI_(muDict_.lookup<scalarField>("I")),
    muJ_(muDict_.lookup<scalarField>("J")),
    muN_(muDict_.lookup<scalarField>("n")),
    kappaDict_(dict_.subDict("kappaModelCoeffs")),
    lambda0n_(kappaDict_.subDict("LAMBDA0").lookup<scalarField>("n")),
    lambda1n_(loadLambda1Coeffs(kappaDict_)),
    lambda2n_(kappaDict_.subDict("LAMBDA2").lookup<scalarField>("n")),
    lambda2n1_(kappaDict_.subDict("LAMBDA2").lookup<scalarField>("n1")),
    delta1_(kappaDict_.subDict("LAMBDA2").lookup<scalar>("delta1")),
    lambda2n2_(kappaDict_.subDict("LAMBDA2").lookup<scalarField>("n2")),
    delta2_(kappaDict_.subDict("LAMBDA2").lookup<scalar>("delta2")),
    lambda2n3_(kappaDict_.subDict("LAMBDA2").lookup<scalarField>("n3")),
    delta3_(kappaDict_.subDict("LAMBDA2").lookup<scalar>("delta3")),
    lambda2n4_(kappaDict_.subDict("LAMBDA2").lookup<scalarField>("n4")),
    delta4_(kappaDict_.subDict("LAMBDA2").lookup<scalar>("delta4")),
    lambda2n5_(kappaDict_.subDict("LAMBDA2").lookup<scalarField>("n5")),
    pc_(dict_.lookup<scalar>("pc")),
    kappaR_(kappaDict_.lookup<scalar>("R"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::IAPWSIF97EOS::saturation
(
    scalar ps,
    scalar Ts,
    bool calcTemp
) const
{
    // This can be done either from pressure or from temperature
    const scalar pStar = 1e6;
    const scalar TStar = 1;

    if (!calcTemp)
    {
        // The Saturation-Pressure Equation (Basic Equation)
        // Range: 273.15 K <= T <= 647.096 K
        // Validation data:
        // T[K] ps [MPa]
        // 300  0.353658941e-2
        // 500  0.263889776e1
        // 600  0.123443146e2
        const scalar theta = (Ts/TStar) + n_[8]/((Ts/TStar) - n_[9]);
        const scalar A = sqr(theta) + n_[0]*theta + n_[1];
        const scalar B = n_[2]*sqr(theta) + n_[3]*theta + n_[4];
        const scalar C = n_[5]*sqr(theta) + n_[6]*theta + n_[7];
        return pStar*pow(2*C/(-B + sqrt(sqr(B) - 4*A*C)), 4);
    }

    // The Saturation-Temperature Equation (Backward Equation)
    // Range: 611.212677 Pa <= p <= 22.064 MPa
    // Validation data:
    // p    T
    // 0.1  372.755919
    // 1    453.035632
    // 10   584.149488
    const scalar beta = pow(ps/pStar, 0.25);
    const scalar E = sqr(beta) + n_[2]*beta + n_[5];
    const scalar F = n_[0]*sqr(beta) + n_[3]*beta + n_[6];
    const scalar G = n_[1]*sqr(beta) + n_[4]*beta + n_[7];
    const scalar D = (2*G)/(-F - sqrt(sqr(F) - 4*E*G));

    // The saturation temperature for given pressure is represented by the equation
    return TStar*(n_[9] + D - sqrt(sqr(n_[9] + D) - 4*(n_[8] + n_[9]*D)))*0.5;
}


Foam::scalar Foam::IAPWSIF97EOS::boundary23
(
    scalar p,
    scalar T,
    bool calcTemp
) const
{
    // This can be done either from pressure or from temperature
    const scalar pStar = 1e6;
    const scalar TStar = 1;

    if (!calcTemp)
    {
        // Calculate pressure [Pa]
        const scalar THETA(T/TStar);
        return pStar*(n23_[0] + n23_[1]*THETA + n23_[2]*sqr(THETA));
    }

    // Calculate temperature [K]
    const scalar phi(p/pStar);
    return TStar*(n23_[3] + sqrt((phi - n23_[4])/n23_[2]));
}


Foam::label Foam::IAPWSIF97EOS::whichRegion
(
    const scalar p,
    const scalar T
) const
{
    const scalar Tmin = 273.15;
    const scalar T13 = 623.15;
    const scalar T25 = 1073.15;
    const scalar T23max = 863.15;
    const scalar Tmax = 2273.15;
    const scalar pmax = 100e6;
    const scalar pmin = 0;
    const scalar pmax5 = 50e6;
    const scalar p123 = 16.5292e6;

    // Check if it is outside range
    bool isOutOfRange = false;
    if (p > pmax || T > Tmax || T < Tmin || p < pmin)
    {
        isOutOfRange = true;
    }
    else if ((T > T25 && T <= Tmax) && p > pmax5)
    {
        isOutOfRange = true;
    }
    if (isOutOfRange)
    {
        FatalErrorInFunction
            << "Combination of p: " << p << " [Pa], T: " << T
            << " [K] is outside of the range for model \"IAPWSIF97EOS\"." << nl
            << abort(FatalError);
    }

    // Now we have to handle only boundary regions

    // Region 1 or 2
    if (T <= T13)
    {
        if (T >= saturation(p))
        {
            return region2;
        }
        else
        {
            return region1;
        }
    }

    // Region 2 or 3
    else if (T > T13 && T <= T25)
    {
        if
        (
            p <= p123
         || (T > T23max || T >= boundary23(p))
        )
        {
            return region2;
        }
        else
        {
            return region3;
        }
    }
    else
    {
        return region5;
    }

    return -1;
}


Foam::scalar Foam::IAPWSIF97EOS::u(const scalar p, const scalar T) const
{
    switch (whichRegion(p, T))
    {
        case region1:
            return region1_.u(p, T);

        case region2:
            return region2_.u(p, T);

        case region3:
            return region3_.u(p, T);

        case region5:
            return region5_.u(p, T);
    }

    return Zero;
}


Foam::scalar Foam::IAPWSIF97EOS::s(const scalar p, const scalar T) const
{
    switch (whichRegion(p, T))
    {
        case region1:
            return region1_.s(p, T);

        case region2:
            return region2_.s(p, T);

        case region3:
            return region3_.s(p, T);

        case region5:
            return region5_.s(p, T);
    }

    return Zero;
}


Foam::scalar Foam::IAPWSIF97EOS::h(const scalar p, const scalar T) const
{
    switch (whichRegion(p, T))
    {
        case region1:
            return region1_.h(p, T);

        case region2:
            return region2_.h(p, T);

        case region3:
            return region3_.h(p, T);

        case region5:
            return region5_.h(p, T);
    }

    return Zero;
}


Foam::scalar Foam::IAPWSIF97EOS::Cv(const scalar p, const scalar T) const
{
    switch (whichRegion(p, T))
    {
        case region1:
            return region1_.Cv(p, T);

        case region2:
            return region2_.Cv(p, T);

        case region3:
            return region3_.Cv(p, T);

        case region5:
            return region5_.Cv(p, T);
    }

    return Zero;
}


Foam::scalar Foam::IAPWSIF97EOS::Cp(const scalar p, const scalar T) const
{
    switch (whichRegion(p, T))
    {
        case region1:
            return region1_.Cp(p, T);

        case region2:
            return region2_.Cp(p, T);

        case region3:
            return region3_.Cp(p, T);

        case region5:
            return region5_.Cp(p, T);
    }

    return Zero;
}


Foam::scalar Foam::IAPWSIF97EOS::w(const scalar p, const scalar T) const
{
    switch (whichRegion(p, T))
    {
        case region1:
            return region1_.w(p, T);

        case region2:
            return region2_.w(p, T);

        case region3:
            return region3_.w(p, T);

        case region5:
            return region5_.w(p, T);
    }

    return Zero;
}


Foam::scalar Foam::IAPWSIF97EOS::rho(const scalar p, const scalar T) const
{
    switch (whichRegion(p, T))
    {
        case region1:
            return region1_.rho(p, T);

        case region2:
            return region2_.rho(p, T);

        case region3:
            return region3_.rho(p, T);

        case region5:
            return region5_.rho(p, T);
    }

    return Zero;
}


Foam::scalar Foam::IAPWSIF97EOS::mu(const scalar p, const scalar T) const
{
    // Check range
    bool outOfRange = false;
    if (T < 273.15 || p < 0 || T > 1273.15 || p > 100e6)
    {
        outOfRange = true;
    }
    else if (T > 1073 && T <= 1173 && p > 50e6)
    {
        outOfRange = true;
    }
    if (outOfRange)
    {
        FatalErrorInFunction
            << "Combination of p: " << p << " [Pa], T: " << T
            << " [K] is outside of the range for "
            << "viscosity \"IAPWSIF97EOS\" model." << nl
            << abort(FatalError);
    }

    const scalar theta = T/647.096;
    scalar PSITheta = 0;
    forAll(muN0_, i)
    {
        PSITheta += muN0_[i]*pow(theta, -i);
    }
    PSITheta = sqrt(theta)/PSITheta;

    const scalar rhoR(rho(p, T));
    const scalar delta = rhoR/322.0;

    scalar PSIDeltaTheta = 0;
    forAll(muI_, i)
    {
        PSIDeltaTheta +=
            muN_[i]*pow(delta - 1, muI_[i])*pow(1.0/theta - 1, muJ_[i]);
    }
    PSIDeltaTheta = Foam::exp(delta*PSIDeltaTheta);

    const scalar muStar = 1e-6;

    return muStar*PSITheta*PSIDeltaTheta;
}


Foam::scalar Foam::IAPWSIF97EOS::KT(const scalar p, const scalar T) const
{
    const scalar rhoR(rho(p, T));
    scalar CpR(Cp(p, T));

    const label regionNumber = whichRegion(p, T);
    if (regionNumber == region3 && (CpR < 0 || CpR > 1e16))
    {
        CpR = 1e16;
    }

    scalar CvR(Cv(p, T));
    const scalar gamma(CpR/CvR);
    const scalar speedOfSound(w(p, T));

    // Isothermal compressibility
    scalar kt(1/(rhoR*sqr(speedOfSound))*gamma);
    if (regionNumber == region3 && (kt < 0 || kt > 1e16))
    {
        kt = 1e16;
    }

    return kt;
}


Foam::scalar Foam::IAPWSIF97EOS::kappa(const scalar p, const scalar T) const
{
    // Check range
    bool outOfRange = false;
    if (T < 273.15 || p < 0 || T > 1273.15 || p > 100e6)
    {
        outOfRange = true;
    }
    else if (T > 1073 && T <= 1173 && p > 50e6)
    {
        outOfRange = true;
    }
    if (outOfRange)
    {
        FatalErrorInFunction
            << "Combination of p: " << p << " [Pa], T: " << T
            << " [K] is outside of the range for "
            << "conductivity \"IAPWSIF97EOS\" model." << nl
            << abort(FatalError);
    }

    const scalar theta = T/647.096;
    scalar LAMBDA0 = 0;
    forAll(lambda0n_, i)
    {
        LAMBDA0 += lambda0n_[i]*pow(theta, -i);
    }
    LAMBDA0 = sqrt(theta)/LAMBDA0;

    const scalar rhoR(rho(p, T));
    const scalar delta = rhoR/322.0;

    scalar LAMBDA1 = 0;
    forAll(lambda1n_.first(), i)
    {
        scalar sum = 0;
        forAll(lambda1n_, j)
        {
            sum += lambda1n_[j][i]*pow(delta - 1, j);
        }
        LAMBDA1 += pow(1/theta - 1, i)*sum;
    }
    LAMBDA1 = Foam::exp(delta*LAMBDA1);

    scalar CpR(Cp(p, T));

    const label regionNumber = whichRegion(p, T);
    if (regionNumber == region3 && (CpR < 0 || CpR > 1e16))
    {
        CpR = 1e16;
    }

    scalar CvR(Cv(p, T));
    const scalar gamma(CpR/CvR);
    const scalar b = gamma;

    // Isothermal compressibility
    scalar Kt(KT(p, T));

    const scalarField* CnPtr = nullptr;
    if (delta <= delta1_)
    {
        CnPtr = &lambda2n1_;
    }
    else if (delta > delta1_ && delta <= delta2_)
    {
        CnPtr = &lambda2n2_;
    }
    else if (delta > delta2_ && delta <= delta3_)
    {
        CnPtr = &lambda2n3_;
    }
    else if (delta > delta3_ && delta <= delta4_)
    {
        CnPtr = &lambda2n4_;
    }
    else
    {
        CnPtr = &lambda2n5_;
    }
    const scalarField& Cn = *CnPtr;

    scalar C = 0;
    forAll(Cn, i)
    {
        C += Cn[i]*pow(delta, i);
    }
    C = 1.0/C;

    scalar B = pc_*delta*Kt - lambda2n_[4]*(1/theta)*C;
    B = max(B, 0);

    const scalar a = lambda2n_[2]*pow(delta*B, lambda2n_[3]);
    scalar A;
    if (a < 1e-7)
    {
        A = 0;
    }
    else
    {
        A =
            (lambda2n_[1]/a)
           *(
                (1 - (1/b))*Foam::atan(a)
              + (a/b)
              - 1 + Foam::exp(-1/(1/a + sqr(a)/(3*sqr(delta))))
            );
    }

    const scalar muStar = 1e-6;
    const scalar PSI(mu(p, T)/muStar);

    // Starnard said that we should use slightly different R = 416.51805
    // however in validation data they provide they still use R = 461.526
    // It isn't clear which gas constant should be used in this case
    // LAMBDA2 still causes slight deviation
    const scalar LAMBDA2 = (lambda2n_[0]*delta*theta*CpR*A)/(PSI*461.526);

    const scalar kappaStar = 1e-3;
    const scalar kappa = kappaStar*(LAMBDA0*LAMBDA1 + LAMBDA2);

    return kappa;
}


// ************************************************************************* //
