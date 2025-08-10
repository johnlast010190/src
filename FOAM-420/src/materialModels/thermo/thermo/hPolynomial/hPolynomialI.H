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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hPolynomial> Foam::hPolynomial::clone() const
{
    return autoPtr<hPolynomial>
    (
        new hPolynomial(*this)
    );
}


inline Foam::autoPtr<Foam::hPolynomial> Foam::hPolynomial::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hPolynomial>
    (
        new hPolynomial(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hPolynomial::CpThermo(const scalar T) const
{
    return CpCoeffs_.value(T);
}


inline Foam::scalar Foam::hPolynomial::HaThermo(const scalar T) const
{
    return hCoeffs_.value(T) + Hf_;
}


inline Foam::scalar Foam::hPolynomial::HsThermo(const scalar T) const
{
    return hCoeffs_.value(T);
}


inline Foam::scalar Foam::hPolynomial::Hc() const
{
    return Hf_;
}


inline Foam::scalar Foam::hPolynomial::SThermo(const scalar T) const
{
    return sCoeffs_.value(T);
}


// ************************************************************************* //
