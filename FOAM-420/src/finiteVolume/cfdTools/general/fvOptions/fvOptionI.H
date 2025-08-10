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
    (c) 2021 Esi Ltd

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::fv::option::name() const
{
    return name_;
}


inline const Foam::fvMesh& Foam::fv::option::mesh() const
{
    return mesh_;
}


inline const Foam::dictionary& Foam::fv::option::coeffs() const
{
    return coeffs_;
}


inline bool Foam::fv::option::active() const
{
    return active_;
}


inline bool Foam::fv::option::MRF() const
{
    return MRF_;
}


inline void Foam::fv::option::setApplied(const label fieldi)
{
    applied_[fieldi] = true;
}


inline void Foam::fv::option::setBoundaryApplied
(
    const word& fieldName, const label fieldPatchi
)
{
    boundaryApplied_[fieldName][fieldPatchi] = true;
}


inline const Foam::word Foam::fv::option::execHookName() const
{
    return execHookTypeNames_[execHook_];
}


inline Foam::Switch& Foam::fv::option::active()
{
    return active_;
}


inline Foam::Switch& Foam::fv::option::MRF()
{
    return MRF_;
}


// ************************************************************************* //
