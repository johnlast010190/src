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

#include "region3IAPWSIF97EOS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::region3IAPWSIF97EOS::region3IAPWSIF97EOS(const dictionary& dict)
:
    I_(dict.lookup<scalarField>("I")),
    J_(dict.lookup<scalarField>("J")),
    n_(dict.lookup<scalarField>("n")),
    Tc_(dict.lookup<scalar>("Tc")),
    rhoc_(dict.lookup<scalar>("rhoc")),
    R_(dict.lookup<scalar>("R"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::region3IAPWSIF97EOS::deltaRatio(const scalar rho) const
{
    return rho/rhoc_;
}


Foam::scalar Foam::region3IAPWSIF97EOS::tauRatio(const scalar T) const
{
    return Tc_/T;
}


Foam::scalar Foam::region3IAPWSIF97EOS::PHI
(
    const scalar delta,
    const scalar tau
) const
{
    scalar PHI = n_[0]*log(delta);
    for (label i = 1; i < I_.size(); i++)
    {
        PHI += n_[i]*pow(delta, I_[i])*pow(tau, J_[i]);
    }
    return PHI;
}


Foam::scalar Foam::region3IAPWSIF97EOS::PHIDelta
(
    const scalar delta,
    const scalar tau
) const
{
    scalar PHIDelta = (n_[0]/delta);
    for (label i = 1; i < I_.size(); i++)
    {
        PHIDelta += n_[i]*I_[i]*pow(delta, I_[i] - 1)*pow(tau, J_[i]);
    }
    return PHIDelta;
}


Foam::scalar Foam::region3IAPWSIF97EOS::PHIDeltaDelta
(
    const scalar delta,
    const scalar tau
) const
{
    scalar PHIDeltaDelta(-n_[0]/sqr(delta));
    for (label i = 1; i < I_.size(); i++)
    {
        PHIDeltaDelta +=
            n_[i]*I_[i]*(I_[i] - 1)*pow(delta, I_[i] - 2)*pow(tau, J_[i]);
    }
    return PHIDeltaDelta;
}


Foam::scalar Foam::region3IAPWSIF97EOS::PHITau
(
    const scalar delta,
    const scalar tau
) const
{
    scalar PHITau = 0;
    for (label i = 1; i < I_.size(); i++)
    {
        PHITau += n_[i]*J_[i]*pow(delta, I_[i])*pow(tau, J_[i] - 1);
    }
    return PHITau;
}


Foam::scalar Foam::region3IAPWSIF97EOS::PHITauTau
(
    const scalar delta,
    const scalar tau
) const
{
    scalar PHITauTau = 0;
    for (label i = 1; i < I_.size(); i++)
    {
        PHITauTau +=
            n_[i]*J_[i]*(J_[i] - 1)*pow(delta, I_[i])*pow(tau, J_[i] - 2);
    }
    return PHITauTau;
}


Foam::scalar Foam::region3IAPWSIF97EOS::PHIDeltaTau
(
    const scalar delta,
    const scalar tau
) const
{
    scalar PHIDeltaTau = 0;
    for (label i = 1; i < I_.size(); i++)
    {
        PHIDeltaTau +=
            n_[i]*I_[i]*J_[i]*pow(delta, I_[i] - 1)*pow(tau, J_[i] - 1);
    }
    return PHIDeltaTau;
}


Foam::scalar Foam::region3IAPWSIF97EOS::u
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));
    return tau*PHITau(delta, tau)*R_*T;
}


Foam::scalar Foam::region3IAPWSIF97EOS::s
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));
    return (tau*PHITau(delta, tau) - PHI(delta, tau))*R_;
}


Foam::scalar Foam::region3IAPWSIF97EOS::h
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));
    return (tau*PHITau(delta, tau) + delta*PHIDelta(delta, tau))*R_*T;
}


Foam::scalar Foam::region3IAPWSIF97EOS::Cv
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));
    return -sqr(tau)*PHITauTau(delta, tau)*R_;
}


Foam::scalar Foam::region3IAPWSIF97EOS::Cp
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));
    return
        (
          - sqr(tau)*PHITauTau(delta, tau)
          + sqr(delta*PHIDelta(delta, tau)
          - delta*tau*PHIDeltaTau(delta, tau))
           /(
                2*delta*PHIDelta(delta, tau)
              + sqr(delta)*PHIDeltaDelta(delta, tau)
            )
        )*R_;
}


Foam::scalar Foam::region3IAPWSIF97EOS::w
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));

    return
        sqrt
        (
            (
                2*delta*PHIDelta(delta, tau)
              + sqr(delta)*PHIDeltaDelta(delta, tau)
              - sqr(delta*PHIDelta(delta, tau)
              - delta*tau*PHIDeltaTau(delta, tau))
               /(sqr(tau)*PHITauTau(delta, tau))
            )*R_*T
        );
}


Foam::scalar Foam::region3IAPWSIF97EOS::rho
(
    const scalar p,
    const scalar T
) const
{
    const scalar tau(tauRatio(T));

    scalar rho = 10;

    // Reset initial rho
    // Make initial rho closer to actual rho by marching
    // through with deltaRho = 50 kg/m^3
    // TODO This is quite expensive is there better option?
    scalar deltaRho = 50;
    scalar PHIDelta = 0;
    for (label i = 1; i < 100; i++)
    {
        if (rho >= 1000)
        {
            break;
        }

        // Derivitive with respect to delta
        PHIDelta = (n_[0]/(rho/rhoc_));
        for (label i = 1; i < I_.size(); i++)
        {
            PHIDelta +=
                n_[i]*I_[i]*pow(rho/rhoc_, I_[i] - 1)*pow(tau, J_[i]);
        }

        const scalar pEst = PHIDelta*(rho/rhoc_)*rho*R_*T;
        if (pEst > p)
        {
            break;
        }
        rho += deltaRho;
    }

    // Newton-Raphson method to match p
    label iter = 0;
    const label maxIter = 100;
    scalar tol = 1e-8; // [kg/m^3]
    scalar rhoEst;
    do
    {
        // Derivitive with respect to delta
        PHIDelta = (n_[0]/(rho/rhoc_));
        for (label i = 1; i < I_.size(); i++)
        {
            PHIDelta +=
                n_[i]*I_[i]*pow(rho/rhoc_, I_[i] - 1)*pow(tau, J_[i]);
        }

        // Pressure
        scalar frho = PHIDelta*(rho/rhoc_)*rho*R_*T - p;

        // Derivitive of pressure function in rho
        scalar dfrho = n_[0]*R_*T;
        for (label i = 1; i < I_.size(); i++)
        {
            dfrho +=
                R_*T*(I_[i] + 1)*n_[i]*I_[i]
               *pow(tau, J_[i])*pow(rho/rhoc_, I_[i]);
        }
        rhoEst = rho;
        rho = rhoEst - frho/dfrho;

        if (iter++ > maxIter)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter
                << abort(FatalError);
        }
    } while (mag(rhoEst - rho) > tol);
    return rho;
}


Foam::scalar Foam::region3IAPWSIF97EOS::mu
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));

    return
        (
          - (
                delta*PHIDelta(delta, tau)
              + sqr(delta)*PHIDeltaDelta(delta, tau)
              + delta*tau*PHIDeltaTau(delta, tau)
            )
           /(
                sqr
                (
                    delta*PHIDelta(delta, tau)
                  - delta*tau*PHIDeltaTau(delta, tau)
                )
              - sqr(tau)
               *PHITauTau(delta, tau)
               *(
                    2*delta*PHIDelta(delta, tau)
                  + sqr(delta)*PHIDeltaDelta(delta, tau)
                )
            )
        )
       *(1/(R_*rhoR));
}


Foam::scalar Foam::region3IAPWSIF97EOS::deltaT
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));

    return
        (
            1
          - (
                delta*PHIDelta(delta, tau)
              - delta*tau*PHIDeltaTau(delta, tau)
            )
           /(
                2*delta*PHIDelta(delta, tau)
              + sqr(delta)*PHIDeltaDelta(delta, tau)
            )
        )
       *(1/rhoR);
}


Foam::scalar Foam::region3IAPWSIF97EOS::betaS
(
    const scalar p,
    const scalar T
) const
{
    scalar rhoR = rho(p, T);
    scalar delta(deltaRatio(rhoR));
    scalar tau(tauRatio(T));

    return
        (
            (
                delta*PHIDelta(delta, tau)
              - delta*tau*PHIDeltaTau(delta, tau)
            )
           /(
                sqr
                (
                    delta*PHIDelta(delta, tau)
                  - delta*tau*PHIDeltaTau(delta, tau)
                )
              - sqr(tau)
               *PHITauTau(delta, tau)
               *(
                    2*delta*PHIDelta(delta, tau)
                  + sqr(delta)*PHIDeltaDelta(delta, tau)
                )
            )
        )
       *(1/(R_*rhoR));
}


// ************************************************************************* //
