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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hTabulated> Foam::hTabulated::clone() const
{
    return autoPtr<hTabulated>
    (
        new hTabulated(*this)
    );
}


inline Foam::autoPtr<Foam::hTabulated> Foam::hTabulated::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hTabulated>
    (
        new hTabulated(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hTabulated::CpThermo
(
    const scalar p,
    const scalar T
) const
{
    return CpTable_()(p, T);
}


inline Foam::scalar Foam::hTabulated::HaThermo
(
    const scalar p,
    const scalar T
) const
{
    return CpTable_().integrateOverY(p, 0, T) + Hf_;
}


inline Foam::scalar Foam::hTabulated::HsThermo
(
    const scalar p,
    const scalar T
) const
{
    return CpTable_().integrateOverY(p, 0, T);
}


inline Foam::scalar Foam::hTabulated::Hc() const
{
    return Hf_;
}


inline Foam::scalar Foam::hTabulated::SThermo
(
    const scalar p,
    const scalar T
) const
{
    NotImplemented;
}


// ************************************************************************* //
