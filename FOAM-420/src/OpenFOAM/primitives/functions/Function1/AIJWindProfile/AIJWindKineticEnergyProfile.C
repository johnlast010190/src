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

#include "primitives/functions/Function1/AIJWindProfile/AIJWindKineticEnergyProfile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    makeScalarFunction1(AIJWindKineticEnergyProfile);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::AIJWindKineticEnergyProfile::AIJWindKineticEnergyProfile
(
    const word& entryName,
    const dictionary& dict
)
:
    AIJWindProfile(entryName,dict),
    AIJPaper_(dict.lookupOrDefault<bool>("consistentTKE", false))
{
}


Foam::Function1Types::AIJWindKineticEnergyProfile::AIJWindKineticEnergyProfile
(
    const AIJWindKineticEnergyProfile& aijp
)
:
    AIJWindProfile(aijp),
    AIJPaper_(aijp.AIJPaper_)
{
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::AIJWindKineticEnergyProfile::~AIJWindKineticEnergyProfile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::AIJWindKineticEnergyProfile::value
(
    const scalar z
) const
{
    if (AIJPaper_)
    {
        return sqr(0.1*pow(z/zG_, -alpha_-0.05)*AIJWindProfile::value(z));
    }
    else
    {
        return
            3.0/2.0
           *sqr(0.1*pow(z/zG_, -alpha_ - 0.05)*AIJWindProfile::value(z));
    }
}

void Foam::Function1Types::AIJWindKineticEnergyProfile::writeData(Ostream& os)
const
{
    Function1<scalar>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os.beginBlock(word(this->name() + "Coeffs"));
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("Uref", Uref_);
    os.writeEntry("Href", Href_);
    os.writeEntry("alpha", alpha_);
    os.writeEntry("zG", zG_);
    os.writeEntry("consistentTKE", AIJPaper_);
    os.endBlock();
}
// ************************************************************************* //
