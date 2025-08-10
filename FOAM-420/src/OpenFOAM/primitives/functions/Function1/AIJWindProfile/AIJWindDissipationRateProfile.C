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

#include "primitives/functions/Function1/AIJWindProfile/AIJWindDissipationRateProfile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    makeScalarFunction1(AIJWindDissipationRateProfile);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::AIJWindDissipationRateProfile::
AIJWindDissipationRateProfile
(
    const word& entryName,
    const dictionary& dict
)
:
    AIJWindKineticEnergyProfile(entryName, dict)
{
}


Foam::Function1Types::AIJWindDissipationRateProfile::
AIJWindDissipationRateProfile
(
    const AIJWindDissipationRateProfile& aijp
)
:
    AIJWindKineticEnergyProfile(aijp)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::AIJWindDissipationRateProfile::
~AIJWindDissipationRateProfile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::AIJWindDissipationRateProfile::value
(
    const scalar z
) const
{
    return
        sqrt(Cmu_)*AIJWindKineticEnergyProfile::value(z)*Uref_
       /Href_*alpha_*pow(z/Href_, alpha_ - 1.0);

}


// ************************************************************************* //
