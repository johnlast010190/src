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

inline Foam::autoPtr<Foam::hConst> Foam::hConst::clone() const
{
    return autoPtr<hConst>
    (
        new hConst(*this)
    );
}


inline Foam::autoPtr<Foam::hConst> Foam::hConst::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hConst>
    (
        new hConst(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hConst::limit(const scalar T) const
{
    return T;
}


inline Foam::scalar Foam::hConst::CpThermo() const
{
    return Cp_;
}


inline Foam::scalar Foam::hConst::HaThermo(const scalar T) const
{
    return Cp_*T + Hf_;
}


inline Foam::scalar Foam::hConst::HsThermo(const scalar T) const
{
    return Cp_*T;
}


inline Foam::scalar Foam::hConst::Hc() const
{
    return Hf_;
}


inline Foam::scalar Foam::hConst::SThermo(const scalar T) const
{
    return Sf_ + Cp_*log(T/Tstd);
}


// ************************************************************************* //
