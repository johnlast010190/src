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

#include "primitives/functions/Function1/AIJWindProfile/AIJWindProfile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    makeScalarFunction1(AIJWindProfile);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::AIJWindProfile::AIJWindProfile
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<scalar>(entryName),
    Cmu_(0.0),
    Uref_(0.0),
    Href_(0.0),
    alpha_(0.0),
    zG_(0.0)
{
    read(dict);
}

Foam::Function1Types::AIJWindProfile::AIJWindProfile
(
    const AIJWindProfile& wp
)
:
    Function1<scalar>(wp),
    Cmu_(wp.Cmu_),
    Uref_(wp.Uref_),
    Href_(wp.Href_),
    alpha_(wp.alpha_),
    zG_(wp.zG_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::AIJWindProfile::~AIJWindProfile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::AIJWindProfile::value(const scalar z) const
{
    return mag(Uref_)*pow(z/Href_, alpha_);
}


void Foam::Function1Types::AIJWindProfile::read(const dictionary& coeffs)
{
    Cmu_ = coeffs.lookupOrDefault<scalar>("Cmu", 0.09);
    Uref_ = coeffs.lookup<scalar>("Uref");
    Href_ = coeffs.lookup<scalar>("Href");
    alpha_ = coeffs.lookup<scalar>("alpha");
    zG_ = coeffs.lookup<scalar>("zG");
}


void Foam::Function1Types::AIJWindProfile::writeData(Ostream& os) const
{
    Function1<scalar>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os.beginBlock(word(this->name() + "Coeffs"));
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("Uref", Uref_);
    os.writeEntry("Href", Href_);
    os.writeEntry("alpha", alpha_);
    os.writeEntry("zG", zG_);
    os.endBlock();
}

// ************************************************************************* //
