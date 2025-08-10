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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::functionTAnIsoSolidTransport<Thermo>::functionTAnIsoSolidTransport
(
    const Thermo& t,
    const Function1<vector>& kappa
)
:
    Thermo(t),
    kappa_(kappa.clone().ptr())
{}


template<class Thermo>
inline Foam::functionTAnIsoSolidTransport<Thermo>::functionTAnIsoSolidTransport
(
    const word& name,
    const functionTAnIsoSolidTransport& ct
)
:
    Thermo(name, ct),
    kappa_(ct.kappa_.clone())
{}


template<class Thermo>
inline Foam::autoPtr<Foam::functionTAnIsoSolidTransport<Thermo>>
Foam::functionTAnIsoSolidTransport<Thermo>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<functionTAnIsoSolidTransport<Thermo>>
    (
        new functionTAnIsoSolidTransport<Thermo>(obr, dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::functionTAnIsoSolidTransport<Thermo>::
kappa(const scalar p, const scalar T) const
{
    return mag(kappa_().value(T));
}

template<class Thermo>
inline Foam::vector Foam::functionTAnIsoSolidTransport<Thermo>::
Kappa(const scalar p, const scalar T) const
{
    return kappa_().value(T);
}


template<class Thermo>
inline Foam::scalar Foam::functionTAnIsoSolidTransport<Thermo>::
mu(const scalar p, const scalar T) const
{
    NotImplemented;
}


template<class Thermo>
inline Foam::vector Foam::functionTAnIsoSolidTransport<Thermo>::
alphah(const scalar p, const scalar T) const
{
    return Kappa()/this->Cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::functionTAnIsoSolidTransport<Thermo>::operator=
(
    const functionTAnIsoSolidTransport<Thermo>& ct
)
{
    kappa_.reset(ct.kappa_.clone());
}


template<class Thermo>
inline void Foam::functionTAnIsoSolidTransport<Thermo>::operator+=
(
    const functionTAnIsoSolidTransport<Thermo>& ct
)
{
    scalar Y1 = this->Y();

    Y1 /= this->Y();
    //scalar Y2 = ct.Y()/this->Y();

    //kappa_() = Y1*kappa_ + Y2*ct.kappa_();
    NotImplemented;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class Thermo>
inline Foam::functionTAnIsoSolidTransport<Thermo> Foam::operator*
(
    const scalar s,
    const functionTAnIsoSolidTransport<Thermo>& ct
)
{
    return functionTAnIsoSolidTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ct),
        ct.kappa_()
    );
}

// ************************************************************************* //
