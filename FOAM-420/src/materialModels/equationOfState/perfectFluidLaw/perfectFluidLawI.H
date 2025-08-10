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
    (c) 2012-2017 OpenFOAM Foundation
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "perfectFluidLaw.H"
#include "global/constants/thermodynamic/thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::autoPtr<Foam::perfectFluidLaw>
Foam::perfectFluidLaw::clone() const
{
    return autoPtr<perfectFluidLaw>(new perfectFluidLaw(*this));
}


inline Foam::autoPtr<Foam::perfectFluidLaw> Foam::perfectFluidLaw::New
(
    const objectRegistry& obr,
    const dictionary& dict,
    const word& phaseName,
    const word& specieName,
    const word& name
)
{
    return autoPtr<perfectFluidLaw>
    (
        new perfectFluidLaw(obr, dict, phaseName, specieName, name)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::perfectFluidLaw::rho
(
    const scalar p,
    const scalar T
) const
{
    return rho0_ + p/(R_*T);
}


inline Foam::scalar Foam::perfectFluidLaw::SDeparture(const scalar p) const
{
    return -R_*log(p/Pstd);
}


inline Foam::scalar Foam::perfectFluidLaw::psi(const scalar T) const
{
    return 1.0/(R_*T);
}


inline Foam::scalar Foam::perfectFluidLaw::Z() const
{
    return 1;
}


// ************************************************************************* //
