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
inline Foam::constAnIsoSolidTransport<Thermo>::constAnIsoSolidTransport
(
    const Thermo& t,
    const vector kappa
)
:
    Thermo(t),
    kappa_(kappa)
{}


template<class Thermo>
inline Foam::constAnIsoSolidTransport<Thermo>::constAnIsoSolidTransport
(
    const word& name,
    const constAnIsoSolidTransport& ct
)
:
    Thermo(name, ct),
    kappa_(ct.kappa_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::constAnIsoSolidTransport<Thermo>>
Foam::constAnIsoSolidTransport<Thermo>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<constAnIsoSolidTransport<Thermo>>
    (
        new constAnIsoSolidTransport<Thermo>(obr, dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::constAnIsoSolidTransport<Thermo>::
kappa(const scalar p, const scalar T) const
{
    return mag(kappa_);
}

template<class Thermo>
inline Foam::vector Foam::constAnIsoSolidTransport<Thermo>::
Kappa(const scalar p, const scalar T) const
{
    return kappa_;
}


template<class Thermo>
inline Foam::scalar Foam::constAnIsoSolidTransport<Thermo>::
mu(const scalar p, const scalar T) const
{
    NotImplemented;
}


template<class Thermo>
inline Foam::vector Foam::constAnIsoSolidTransport<Thermo>::
alphah(const scalar p, const scalar T) const
{
    return kappa_/this->Cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::constAnIsoSolidTransport<Thermo>::operator=
(
    const constAnIsoSolidTransport<Thermo>& ct
)
{
    kappa_ = ct.kappa_;
}


template<class Thermo>
inline void Foam::constAnIsoSolidTransport<Thermo>::operator+=
(
    const constAnIsoSolidTransport<Thermo>& ct
)
{
    scalar Y1 = this->Y();

    Y1 /= this->Y();
    scalar Y2 = ct.Y()/this->Y();

    kappa_ = Y1*kappa_ + Y2*ct.kappa_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Thermo>
inline Foam::constAnIsoSolidTransport<Thermo> Foam::operator*
(
    const scalar s,
    const constAnIsoSolidTransport<Thermo>& ct
)
{
    return constAnIsoSolidTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.kappa_
    );
}

// ************************************************************************* //
