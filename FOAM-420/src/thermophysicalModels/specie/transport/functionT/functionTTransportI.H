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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Zero/ZeroConstant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::functionTTransport<Thermo>::functionTTransport
(
    const Thermo& t,
    const Function1<scalar>& mu,
    const Function1<scalar>& kappa
)
:
    Thermo(t),
    mu_(mu.clone().ptr()),
    kappa_(kappa.clone().ptr())
{}


template<class Thermo>
inline Foam::functionTTransport<Thermo>::functionTTransport
(
    const word& name,
    const functionTTransport& ct
)
:
    Thermo(name, ct),
    mu_(ct.mu_.clone()),
    kappa_(ct.kappa_.clone())
{}


template<class Thermo>
inline Foam::autoPtr<Foam::functionTTransport<Thermo>>
Foam::functionTTransport<Thermo>::clone() const
{
    return autoPtr<functionTTransport<Thermo>>
    (
        new functionTTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::functionTTransport<Thermo>>
Foam::functionTTransport<Thermo>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<functionTTransport<Thermo>>
    (
        new functionTTransport<Thermo>(obr, dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::functionTTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    return mu_().value(T);
}


template<class Thermo>
inline Foam::scalar Foam::functionTTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    return kappa_().value(T);
}


template<class Thermo>
inline Foam::scalar Foam::functionTTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return kappa(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::functionTTransport<Thermo>::operator=
(
    const functionTTransport<Thermo>& ct
)
{
    Thermo::operator=(ct);

    mu_.reset(ct.mu_.clone());
    kappa_.reset(ct.kappa_.clone());
}


template<class Thermo>
inline void Foam::functionTTransport<Thermo>::operator+=
(
    const functionTTransport<Thermo>& st
)
{
    //scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        //Y1 /= this->Y();
        //scalar Y2 = st.Y()/this->Y();

        // TODO: Implement with a new Function1 which is a linear combination
        // of other Function1s
        NotImplemented;
        //mu_ = Y1*mu_ + Y2*st.mu_;
        //kappa_ = Y1*kappa_ + Y2*st.kappa_;
    }
}


template<class Thermo>
inline void Foam::functionTTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::functionTTransport<Thermo> Foam::operator+
(
    const functionTTransport<Thermo>& ct1,
    const functionTTransport<Thermo>& ct2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ct1) + static_cast<const Thermo&>(ct2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return functionTTransport<Thermo>
        (
            t,
            Function1Types::ZeroConstant<scalar>("mu", dictionary()),
            Function1Types::ZeroConstant<scalar>("kappa", dictionary())
        );
    }
    else
    {
        //scalar Y1 = ct1.Y()/t.Y();
        //scalar Y2 = ct2.Y()/t.Y();

        NotImplemented;
    }
}


template<class Thermo>
inline Foam::functionTTransport<Thermo> Foam::operator*
(
    const scalar s,
    const functionTTransport<Thermo>& ct
)
{
    return functionTTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.mu_,
        ct.kappa_
    );
}


// ************************************************************************* //
