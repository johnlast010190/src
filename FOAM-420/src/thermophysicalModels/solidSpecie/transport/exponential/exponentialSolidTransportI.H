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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::exponentialSolidTransport<Thermo>::exponentialSolidTransport
(
    const Thermo& t,
    const scalar kappa0,
    const scalar n0,
    const scalar Tref
)
:
    Thermo(t),
    kappa0_(kappa0),
    n0_(n0),
    Tref_(Tref)
{}


template<class Thermo>
inline Foam::exponentialSolidTransport<Thermo>::exponentialSolidTransport
(
    const word& name,
    const exponentialSolidTransport& ct
)
:
    Thermo(name, ct),
    kappa0_(ct.kappa0_),
    n0_(ct.n0_),
    Tref_(ct.Tref_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::exponentialSolidTransport<Thermo>>
Foam::exponentialSolidTransport<Thermo>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<exponentialSolidTransport<Thermo>>
    (
        new exponentialSolidTransport<Thermo>(obr, dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::exponentialSolidTransport<Thermo>::kappa
(
    const scalar p, const scalar T
) const
{
    return (kappa0_*pow(T/Tref_, n0_));
}


template<class Thermo>
inline Foam::vector Foam::exponentialSolidTransport<Thermo>::Kappa
(
    const scalar p, const scalar T
) const
{
    const scalar kappa(kappa0_*pow(T/Tref_, n0_));
    return vector(kappa, kappa, kappa);
}


template<class Thermo>
inline Foam::scalar Foam::exponentialSolidTransport<Thermo>::
mu(const scalar p, const scalar T) const
{
    NotImplemented;
}


template<class Thermo>
inline Foam::scalar Foam::exponentialSolidTransport<Thermo>::
alphah(const scalar p, const scalar T) const
{
    return kappa(p, T)/this->Cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


template<class Thermo>
inline void Foam::exponentialSolidTransport<Thermo>::operator=
(
    const exponentialSolidTransport<Thermo>& ct
)
{
    kappa0_ = ct.kappa0_;
    n0_ = ct.n0_;
    Tref_ = ct.Tref_;
}


template<class Thermo>
inline void Foam::exponentialSolidTransport<Thermo>::operator+=
(
    const exponentialSolidTransport<Thermo>& ct
)
{
    scalar Y1 = this->Y();

    Y1 /= this->Y();
    scalar Y2 = ct.Y()/this->Y();

    kappa0_ = Y1*kappa0_ + Y2*ct.kappa0_;
    n0_ = (Y1*n0_ + Y2*ct.n0_);
    Tref_ = (Y1*Tref_ + Y2*ct.Tref_);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::exponentialSolidTransport<Thermo> Foam::operator*
(
    const scalar s,
    const exponentialSolidTransport<Thermo>& ct
)
{
    return exponentialSolidTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.kappa0_,
        ct.n0_,
        ct.Tref_
    );
}

// ************************************************************************* //
