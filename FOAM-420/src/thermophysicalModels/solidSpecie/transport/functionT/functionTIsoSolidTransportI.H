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

template<class thermo>
inline Foam::functionTIsoSolidTransport<thermo>::functionTIsoSolidTransport
(
    const thermo& t,
    const Function1<scalar>& kappa
)
:
    thermo(t),
    kappa_(kappa.clone().ptr())
{}


template<class thermo>
inline Foam::functionTIsoSolidTransport<thermo>::functionTIsoSolidTransport
(
    const word& name,
    const functionTIsoSolidTransport& ct
)
:
    thermo(name, ct),
    kappa_(ct.kappa_.clone())
{}


template<class Thermo>
inline Foam::autoPtr<Foam::functionTIsoSolidTransport<Thermo>>
Foam::functionTIsoSolidTransport<Thermo>::New
(
    const objectRegistry& obr,
    const dictionary& dict
)
{
    return autoPtr<functionTIsoSolidTransport<Thermo>>
    (
        new functionTIsoSolidTransport<Thermo>(obr, dict)
    );
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class thermo>
inline Foam::scalar Foam::functionTIsoSolidTransport<thermo>::
kappa(const scalar p, const scalar T) const
{
    return kappa_().value(T);
}

template<class thermo>
inline Foam::vector Foam::functionTIsoSolidTransport<thermo>::
Kappa(const scalar p, const scalar T) const
{
    scalar kappa = kappa_().value(T);
    return vector(kappa, kappa, kappa);
}


template<class thermo>
inline Foam::scalar Foam::functionTIsoSolidTransport<thermo>::
mu(const scalar p, const scalar T) const
{
    NotImplemented;
}


template<class thermo>
inline Foam::scalar Foam::functionTIsoSolidTransport<thermo>::
alphah(const scalar p, const scalar T) const
{
    return kappa()/this->Cp(p, T);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class thermo>
inline void Foam::functionTIsoSolidTransport<thermo>::operator=
(
    const functionTIsoSolidTransport<thermo>& ct
)
{
    thermo::operator=(ct);
    kappa_ = ct.kappa_;
}


template<class thermo>
inline void Foam::functionTIsoSolidTransport<thermo>::operator+=
(
    const functionTIsoSolidTransport<thermo>& ct
)
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


template<class thermo>
inline Foam::functionTIsoSolidTransport<thermo> Foam::operator*
(
    const scalar s,
    const functionTIsoSolidTransport<thermo>& ct
)
{
    return functionTIsoSolidTransport<thermo>
    (
        s*static_cast<const thermo&>(ct),
        ct.kappa_
    );
}


// ************************************************************************* //
