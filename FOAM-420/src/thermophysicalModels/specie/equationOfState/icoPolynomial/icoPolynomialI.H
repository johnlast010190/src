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
    (c) 2016 Esi Ltd.
    (c) 2011-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "equationOfState/icoPolynomial/icoPolynomial.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoPolynomial<Specie, PolySize>::icoPolynomial
(
    const Specie& sp,
    const fvMesh& mesh,
    const Polynomial<PolySize>& rhoCoeffs,
    scalar rhoMin,
    scalar rhoMax
)
:
    Specie(sp),
    rhoCoeffs_(rhoCoeffs),
    rhoMin_(rhoMin),
    rhoMax_(rhoMax)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::icoPolynomial<Specie, PolySize>::icoPolynomial
(
    const word& name,
    const icoPolynomial<Specie, PolySize>& ip
)
:
    Specie(name, ip),
    rhoCoeffs_(ip.rhoCoeffs_),
    rhoMin_(ip.rhoMin_),
    rhoMax_(ip.rhoMax_)
{}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<Specie, PolySize>>
Foam::icoPolynomial<Specie, PolySize>::clone() const
{
    return autoPtr<icoPolynomial<Specie, PolySize>>
    (
        new icoPolynomial<Specie, PolySize>(*this)
    );
}


template<class Specie, int PolySize>
inline Foam::autoPtr<Foam::icoPolynomial<Specie, PolySize>>
Foam::icoPolynomial<Specie, PolySize>::New(const objectRegistry& obr, const dictionary& dict)
{
    return autoPtr<icoPolynomial<Specie, PolySize>>
    (
        new icoPolynomial<Specie, PolySize>(obr, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::rho
(
    scalar p,
    scalar T
) const
{
    return max(rhoMin_, min(rhoMax_, rhoCoeffs_.value(T)));
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::H
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::Cp
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::S
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::psi
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::Z
(
    scalar p,
    scalar T
) const
{
    return 0;
}


template<class Specie, int PolySize>
inline Foam::scalar Foam::icoPolynomial<Specie, PolySize>::CpMCv
(
    scalar p,
    scalar T
) const
{
    scalar drdt = rhoCoeffs_.derivative(T);
    scalar rw = rhoCoeffs_.value(T);

    if (rw > rhoMax_)
    {
        rw = rhoMax_;
        drdt = 0;
    }
    else if (rw < rhoMin_)
    {
        rw = rhoMin_;
        drdt = 0;
    }

    return -(p/sqr(rw))*drdt;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
inline void Foam::icoPolynomial<Specie, PolySize>::operator=
(
    const icoPolynomial<Specie, PolySize>& ip
)
{
    Specie::operator=(ip);

    rhoCoeffs_ = ip.rhoCoeffs_;

    rhoMin_ = ip.rhoMin_;
    rhoMax_ = ip.rhoMax_;
}


template<class Specie, int PolySize>
inline void Foam::icoPolynomial<Specie, PolySize>::operator+=
(
    const icoPolynomial<Specie, PolySize>& ip
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ip);

    if (mag(this->Y()) > SMALL)

    {
        Y1 /= this->Y();
        const scalar Y2 = ip.Y()/this->Y();

        rhoCoeffs_ = Y1*rhoCoeffs_ + Y2*ip.rhoCoeffs_;
        rhoMin_ = Y1*rhoMin_ + Y2*ip.rhoMin_;
        rhoMax_ = Y1*rhoMax_ + Y2*ip.rhoMax_;
    }
}


template<class Specie, int PolySize>
inline void Foam::icoPolynomial<Specie, PolySize>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie, int PolySize>
Foam::icoPolynomial<Specie, PolySize> Foam::operator+
(
    const icoPolynomial<Specie, PolySize>& ip1,
    const icoPolynomial<Specie, PolySize>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
      + static_cast<const Specie&>(ip2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return icoPolynomial<Specie, PolySize>
        (
            sp,
            ip1.mesh(),
            ip1.rhoCoeffs_
        );
    }
    else
    {
        const scalar Y1 = ip1.Y()/sp.Y();
        const scalar Y2 = ip2.Y()/sp.Y();

        return icoPolynomial<Specie, PolySize>
        (
            sp,
            ip1.mesh(),
            Y1*ip1.rhoCoeffs_ + Y2*ip2.rhoCoeffs_,
            Y1*ip1.rhoMin_ + Y2*ip2.rhoMin_,
            Y1*ip1.rhoMax_ + Y2*ip2.rhoMax_

        );
    }
}


template<class Specie, int PolySize>
Foam::icoPolynomial<Specie, PolySize> Foam::operator*
(
    const scalar s,
    const icoPolynomial<Specie, PolySize>& ip
)
{
    return icoPolynomial<Specie, PolySize>
    (
        s*static_cast<const Specie&>(ip),
        ip.mesh(),
        ip.rhoCoeffs_,
        ip.rhoMin_,
        ip.rhoMax_
    );
}


template<class Specie, int PolySize>
Foam::icoPolynomial<Specie, PolySize> Foam::operator==
(
    const icoPolynomial<Specie, PolySize>& ip1,
    const icoPolynomial<Specie, PolySize>& ip2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ip1)
     == static_cast<const Specie&>(ip2)
    );

    const scalar Y1 = ip1.Y()/sp.Y();
    const scalar Y2 = ip2.Y()/sp.Y();

    return icoPolynomial<Specie, PolySize>
    (
        sp,
        ip1.mesh(),
        Y2*ip2.rhoCoeffs_ - Y1*ip1.rhoCoeffs_,
        Y2*ip2.rhoMin_ - Y1*ip1.rhoMin_,
        Y2*ip2.rhoMax_ - Y1*ip1.rhoMax_
    );
}


// ************************************************************************* //
