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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "IAPWSIF97EOSModels.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::IAPWSIF97EOSModels>
Foam::IAPWSIF97EOSModels::clone() const
{
    return autoPtr<IAPWSIF97EOSModels>(new IAPWSIF97EOSModels(*this));
}


inline Foam::autoPtr<Foam::IAPWSIF97EOSModels> Foam::IAPWSIF97EOSModels::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<IAPWSIF97EOSModels>
    (
        new IAPWSIF97EOSModels(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::IAPWSIF97EOSModels::checkPhaseChange
(
    const scalar p,
    const scalar T
) const
{
    if (!supportPhaseChange_)
    {
        const label region = waterSteamModel_.whichRegion(p, T);
        const label region1 = waterSteamModel_.region1;
        const label region2 = waterSteamModel_.region2;
        const label region3 = waterSteamModel_.region3;

        if
        (
            regionNumber_ != region
         && (
                (regionNumber_ == region1  && region == region2)
             || (regionNumber_ == region2  && region == region1)
             || (regionNumber_ == region2  && region == region3)
             || (regionNumber_ == region3  && region == region2)
            )
        )
        {
            FatalErrorInFunction
                << "The model " << IAPWSIF97EOSModels::typeName
                << " switched phase region from " << (regionNumber_ + 1)
                << " to " << (region + 1) << " if the phase change "
                << " is accounted for in the solver, please use a flag "
                << "\"supportPhaseChange true;\" in the model subdict."
                << exit(FatalError);
        }
    }
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::rho
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.rho(p, T);
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::CpMCv
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return (waterSteamModel_.Cp(p, T) - waterSteamModel_.Cv(p, T));
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::Hs
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.h(p, T);
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::Ha
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.h(p, T) + Hf_;
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::S
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.s(p, T) + Sf_;
}



inline Foam::scalar Foam::IAPWSIF97EOSModels::Cp
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.Cp(p, T);
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::Cv
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.Cv(p, T);
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::Hc() const
{
    return Hf_;
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::mu
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.mu(p, T);
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::kappa
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.kappa(p, T);
}


inline Foam::scalar Foam::IAPWSIF97EOSModels::psi
(
    const scalar p,
    const scalar T
) const
{
    checkPhaseChange(p, T);
    return waterSteamModel_.KT(p, T);
}


// ************************************************************************* //
