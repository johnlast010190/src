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
    (c) 2019-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::hFunctionT> Foam::hFunctionT::clone() const
{
    return autoPtr<hFunctionT>
    (
        new hFunctionT(*this)
    );
}


inline Foam::autoPtr<Foam::hFunctionT> Foam::hFunctionT::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<hFunctionT>
    (
        new hFunctionT(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::hFunctionT::CpThermo(const scalar T) const
{
    return Cp_().value(T);
}


inline Foam::scalar Foam::hFunctionT::HaThermo(const scalar T) const
{
    return Cp_().integrate(0, T) + Hf_;
}


inline Foam::scalar Foam::hFunctionT::Hc() const
{
    return Hf_;
}


inline Foam::scalar Foam::hFunctionT::SThermo(const scalar T) const
{
    return Cp_().integrateYoverX(Tstd, T) + Sf_;
}


// ************************************************************************* //
