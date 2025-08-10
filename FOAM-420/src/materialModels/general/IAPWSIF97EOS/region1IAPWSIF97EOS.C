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

#include "region1IAPWSIF97EOS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::region1IAPWSIF97EOS::region1IAPWSIF97EOS(const dictionary& dict)
:
    I_(dict.lookup<scalarField>("I")),
    J_(dict.lookup<scalarField>("J")),
    n_(dict.lookup<scalarField>("n")),
    pStar_(dict.lookup<scalar>("pStar")),
    TStar_(dict.lookup<scalar>("TStar")),
    R_(dict.lookup<scalar>("R"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::region1IAPWSIF97EOS::phiRatio(const scalar p) const
{
    return p/pStar_;
}


Foam::scalar Foam::region1IAPWSIF97EOS::tauRatio(const scalar T) const
{
    return TStar_/T;
}


Foam::scalar Foam::region1IAPWSIF97EOS::gamma
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gamma = 0;
    forAll(I_, i)
    {
        gamma += n_[i]*pow(7.1 - phi, I_[i])*pow(tau - 1.222, J_[i]);
    }
    return gamma;
}


Foam::scalar Foam::region1IAPWSIF97EOS::gammaPhi
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammaPhi = 0;
    forAll(I_, i)
    {
        gammaPhi -=
            n_[i]*I_[i]*pow(7.1 - phi, I_[i] - 1)*pow(tau - 1.222, J_[i]);
    }
    return gammaPhi;
}


Foam::scalar Foam::region1IAPWSIF97EOS::gammaPhiPhi
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammaPhiPhi = 0;
    forAll(I_, i)
    {
        gammaPhiPhi +=
            n_[i]*I_[i]*(I_[i] - 1)
           *pow(7.1 - phi, I_[i] - 2)
           *pow(tau - 1.222, J_[i]);
    }
    return gammaPhiPhi;
}


Foam::scalar Foam::region1IAPWSIF97EOS::gammaTau
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammaTau = 0;
    forAll(I_, i)
    {
        gammaTau +=
            n_[i]*J_[i]*pow(7.1 - phi, I_[i])*pow(tau - 1.222, J_[i] - 1);
    }
    return gammaTau;
}


Foam::scalar Foam::region1IAPWSIF97EOS::gammaTauTau
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammaTauTau = 0;
    forAll(I_, i)
    {
        gammaTauTau +=
            n_[i]*J_[i]*(J_[i] - 1)
           *pow(7.1 - phi, I_[i])
           *pow(tau - 1.222, J_[i] - 2);
    }
    return gammaTauTau;
}


Foam::scalar Foam::region1IAPWSIF97EOS::gammaPhiTau
(
    const scalar phi,
    const scalar tau
) const
{
    scalar gammaPhiTau = 0;
    forAll(I_, i)
    {
        gammaPhiTau -=
            n_[i]*J_[i]*I_[i]
           *pow(7.1 - phi, I_[i] - 1)
           *pow(tau - 1.222, J_[i] - 1);
    }
    return gammaPhiTau;
}


Foam::scalar Foam::region1IAPWSIF97EOS::v
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return phi*gammaPhi(phi, tau)*R_*T/p;
}


Foam::scalar Foam::region1IAPWSIF97EOS::u
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return R_*T*(tau*gammaTau(phi, tau) - phi*gammaPhi(phi, tau));
}


Foam::scalar Foam::region1IAPWSIF97EOS::s
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return R_*(tau*gammaTau(phi, tau) - gamma(phi, tau));
}


Foam::scalar Foam::region1IAPWSIF97EOS::h
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return R_*T*tau*gammaTau(phi, tau);
}


Foam::scalar Foam::region1IAPWSIF97EOS::Cv
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
          - sqr(tau)*gammaTauTau(phi, tau)
          + sqr
            (
                gammaPhi(phi, tau) - tau*gammaPhiTau(phi, tau)
            )/gammaPhiPhi(phi, tau)
        );
}


Foam::scalar Foam::region1IAPWSIF97EOS::Cp
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return -R_*sqr(tau)*gammaTauTau(phi, tau);
}


Foam::scalar Foam::region1IAPWSIF97EOS::w
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    const scalar tgammaPhi(gammaPhi(phi, tau));

    return
        sqrt
        (
            R_*T*
            (
                sqr(tgammaPhi)
               /(
                    sqr(tgammaPhi - tau*gammaPhiTau(phi, tau))
                   /(sqr(tau)*gammaTauTau(phi, tau))
                  - gammaPhiPhi(phi, tau)
                )
            )
        );
}


Foam::scalar Foam::region1IAPWSIF97EOS::rho
(
    const scalar p,
    const scalar T
) const
{
    const scalar phi(phiRatio(p));
    const scalar tau(tauRatio(T));
    return p/(R_*T*phi*gammaPhi(phi, tau));
}


// ************************************************************************* //
