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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/RHWindProfile/RHWindDissipationRateProfile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    makeScalarFunction1(RHWindDissipationRateProfile);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::RHWindDissipationRateProfile::RHWindDissipationRateProfile
(
    const word& entryName,
    const dictionary& dict
)
:
    RHWindKineticEnergyProfile(entryName,dict)
{
}


Foam::Function1Types::RHWindDissipationRateProfile::RHWindDissipationRateProfile
(
     const RHWindDissipationRateProfile& wp
)
:
    RHWindKineticEnergyProfile(wp)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::RHWindDissipationRateProfile::~RHWindDissipationRateProfile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::RHWindDissipationRateProfile::value(const scalar z_) const
{

    scalar z0min(max(z0_, 0.00001));
    scalar Ustar = Kappa_ * Uref_ / (log((Href_  + z0min)/z0min));

    return pow3(mag(Ustar))/(Kappa_*(z_ + z0min));
}
// ************************************************************************* //
