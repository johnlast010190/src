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
    (c) 2015-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hRefConst> Foam::hRefConst::clone() const
{
    return autoPtr<hRefConst>
    (
        new hRefConst(*this)
    );
}


inline Foam::autoPtr<Foam::hRefConst> Foam::hRefConst::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hRefConst>
    (
        new hRefConst(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hRefConst::CpThermo() const
{
    return Cp_;
}


inline Foam::scalar Foam::hRefConst::HaThermo(const scalar T) const
{
    return Cp_*(T - Tref_) + Href_ + Hf_;
}


inline Foam::scalar Foam::hRefConst::HsThermo(const scalar T) const
{
    return Cp_*(T - Tref_) + Href_;
}


inline Foam::scalar Foam::hRefConst::Hc() const
{
    return Hf_;
}


inline Foam::scalar Foam::hRefConst::SThermo(const scalar T) const
{
    return Cp_*log(T/Tstd);
}


// ************************************************************************* //
