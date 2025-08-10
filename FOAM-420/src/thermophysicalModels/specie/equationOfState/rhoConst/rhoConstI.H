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
    (c) 2012-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "equationOfState/rhoConst/rhoConst.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoConst<Specie>::rhoConst
(
    const Specie& sp,
    const fvMesh& mesh,
    const scalar rho
)
:
    Specie(sp),
    rho_(rho)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoConst<Specie>::rhoConst
(
    const word& name,
    const rhoConst<Specie>& ico
)
:
    Specie(name, ico),
    rho_(ico.rho_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::rhoConst<Specie>>
Foam::rhoConst<Specie>::clone() const
{
    return autoPtr<rhoConst<Specie>>(new rhoConst<Specie>(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::rho(scalar p, scalar T) const
{
    return rho_;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::H(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::Cp(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::S(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::psi(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::Z(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::rhoConst<Specie>::CpMCv(scalar p, scalar T) const
{
    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::rhoConst<Specie>::operator+=(const rhoConst<Specie>& ico)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ico);

    if (mag(this->Y()) > SMALL)
    {
        rho_ = this->Y()/(Y1/rho_ + ico.Y()/ico.rho_);
    }
}


template<class Specie>
inline void Foam::rhoConst<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator+
(
    const rhoConst<Specie>& ico1,
    const rhoConst<Specie>& ico2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ico1)
      + static_cast<const Specie&>(ico2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return rhoConst<Specie>
        (
            sp,
            ico1.mesh(),
            ico1.rho_
        );
    }
    else
    {
        return rhoConst<Specie>
        (
            sp,
            ico1.mesh(),
            sp.Y()/(ico1.Y()/ico1.rho_ + ico2.Y()/ico2.rho_)
        );
    }
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator*
(
    const scalar s,
    const rhoConst<Specie>& ico
)
{
    return rhoConst<Specie>(s*static_cast<const Specie&>(ico), ico.mesh(), ico.rho_);
}


template<class Specie>
inline Foam::rhoConst<Specie> Foam::operator==
(
    const rhoConst<Specie>& ico1,
    const rhoConst<Specie>& ico2
)
{
    return rhoConst<Specie>
    (
        static_cast<const Specie&>(ico1) == static_cast<const Specie&>(ico2),
        ico1.mesh(),
        NAN
    );
}


// ************************************************************************* //
