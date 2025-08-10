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

#include "region5IAPWSIF97EOS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::region5IAPWSIF97EOS::region5IAPWSIF97EOS(const dictionary& dict)
:
    I0_(dict.lookup<scalarField>("I0")),
    J0_(dict.lookup<scalarField>("J0")),
    n0_(dict.lookup<scalarField>("n0")),
    I_(dict.lookup<scalarField>("I")),
    J_(dict.lookup<scalarField>("J")),
    n_(dict.lookup<scalarField>("n")),
    pStar_(dict.lookup<scalar>("pStar")),
    TStar_(dict.lookup<scalar>("TStar")),
    R_(dict.lookup<scalar>("R"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::region5IAPWSIF97EOS::phiRatio(const scalar p) const
{
    return p/pStar_;
}


Foam::scalar Foam::region5IAPWSIF97EOS::tauRatio(const scalar T) const
{
    return TStar_/T;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gamma0
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gamma0 = log(phi);
    forAll(I0_, i)
    {
        gamma0 += n0_[i]*pow(tau, J0_[i]);
    }
    return gamma0;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gamma0Phi
(
    const scalar phi,
    const scalar tau
) const
{
    return 1/phi;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gamma0PhiPhi
(
    const scalar phi,
    const scalar tau
) const
{
    return -1/sqr(phi);
}


Foam::scalar Foam::region5IAPWSIF97EOS::gamma0Tau
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gamma0Tau = 0;
    forAll(J0_, i)
    {
        gamma0Tau += n0_[i]*J0_[i]*pow(tau, J0_[i] - 1);
    }
    return gamma0Tau;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gamma0TauTau
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gamma0TauTau = 0;
    forAll(J0_, i)
    {
        gamma0TauTau += n0_[i]*J0_[i]*(J0_[i] - 1)*pow(tau, J0_[i] - 2);
    }
    return gamma0TauTau;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gamma0PhiTau
(
    const scalar phi,
    const scalar tau
) const
{
    return 0;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gammar
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammar = 0;
    forAll(I_, i)
    {
        gammar += n_[i]*pow(phi, I_[i])*pow(tau, J_[i]);
    }
    return gammar;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gammarPhi
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammarPhi = 0;
    forAll(I_, i)
    {
        gammarPhi += n_[i]*I_[i]*pow(phi, I_[i] - 1)*pow(tau, J_[i]);
    }
    return gammarPhi;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gammarPhiPhi
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammarPhiPhi = 0;
    forAll(I_, i)
    {
        gammarPhiPhi +=
            n_[i]*I_[i]*(I_[i] - 1)*pow(phi, I_[i] - 2)*pow(tau, J_[i]);
    }
    return gammarPhiPhi;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gammarTau
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammarTau = 0;
    forAll(I_, i)
    {
        gammarTau += n_[i]*J_[i]*pow(phi, I_[i])*pow(tau, J_[i] - 1);
    }
    return gammarTau;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gammarTauTau
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammarTauTau = 0;
    forAll(I_, i)
    {
        gammarTauTau +=
            n_[i]*J_[i]*(J_[i] - 1)*pow(phi, I_[i])*pow(tau, J_[i] - 2);
    }
    return gammarTauTau;
}


Foam::scalar Foam::region5IAPWSIF97EOS::gammarPhiTau
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammarPhiTau = 0;
    forAll(I_, i)
    {
        gammarPhiTau +=
            n_[i]*J_[i]*I_[i]*pow(phi, I_[i] - 1)*pow(tau, J_[i] - 1);
    }
    return gammarPhiTau;
}


Foam::scalar Foam::region5IAPWSIF97EOS::v
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return phi*(gamma0Phi(phi, tau) + gammarPhi(phi, tau))*R_*T/p;
}


Foam::scalar Foam::region5IAPWSIF97EOS::u
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return
        R_*T
       *(
            tau*(gamma0Tau(phi, tau) + gammarTau(phi, tau))
          - phi*(gamma0Phi(phi, tau) + gammarPhi(phi, tau))
        );
}


Foam::scalar Foam::region5IAPWSIF97EOS::s
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return
        R_
       *(
            tau*(gamma0Tau(phi, tau) + gammarTau(phi, tau))
          - (gamma0(phi, tau) + gammar(phi, tau))
        );
}


Foam::scalar Foam::region5IAPWSIF97EOS::h
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return R_*T*tau*(gamma0Tau(phi, tau) + gammarTau(phi, tau));
}


Foam::scalar Foam::region5IAPWSIF97EOS::Cv
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return
       R_
      *(
          - sqr(tau)*(gamma0TauTau(phi, tau) + gammarTauTau(phi, tau))
          - sqr(1.0 + phi*gammarPhi(phi, tau) - tau*phi*gammarPhiTau(phi, tau))
           /(1.0 - sqr(phi)*gammarPhiPhi(phi, tau))
       );
}


Foam::scalar Foam::region5IAPWSIF97EOS::Cp
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return -R_*sqr(tau)*(gamma0TauTau(phi, tau) + gammarTauTau(phi, tau));
}


Foam::scalar Foam::region5IAPWSIF97EOS::w
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return
        sqrt
        (
            R_*T*
            (
                (
                    1
                  + 2*phi*gammarPhi(phi, tau)
                  + sqr(phi)*sqr(gammarPhi(phi, tau))
                )
               /(
                    1
                  - sqr(phi)*gammarPhiPhi(phi, tau)
                  + sqr
                    (
                        1
                      + phi*gammarPhi(phi, tau)
                      - tau*phi*gammarPhiTau(phi, tau)
                    )
                   /(
                        sqr(tau)
                      *(gamma0TauTau(phi, tau) + gammarTauTau(phi, tau))
                    )
                )
            )
        );
}


Foam::scalar Foam::region5IAPWSIF97EOS::rho
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return p/(R_*T*phi*(gamma0Phi(phi, tau) + gammarPhi(phi, tau)));
}


// ************************************************************************* //
