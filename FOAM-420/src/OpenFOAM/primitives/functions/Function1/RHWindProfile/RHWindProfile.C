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

#include "primitives/functions/Function1/RHWindProfile/RHWindProfile.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    makeScalarFunction1(RHWindProfile);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Function1Types::RHWindProfile::RHWindProfile
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<scalar>(entryName),
    Cmu_(0.0),
    Kappa_(0.0),
    Uref_(0.0),
    Href_(0.0),
    z0_(0.0)
{
    read(dict);
}

Foam::Function1Types::RHWindProfile::RHWindProfile
(
    const RHWindProfile& wp
)
:
    Function1<scalar>(wp),
    Cmu_(wp.Cmu_),
    Kappa_(wp.Kappa_),
    Uref_(wp.Uref_),
    Href_(wp.Href_),
    z0_(wp.z0_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::RHWindProfile::~RHWindProfile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::RHWindProfile::Ufrictional() const
{
    scalar z0min(max(z0_, 0.00001));
    return Kappa_ * Uref_ / (log((Href_  + z0min)/z0min));
}

Foam::scalar Foam::Function1Types::RHWindProfile::value(const scalar z_) const
{

  //limit surfacxe roughness to 0.01mm to prevent singularity
    scalar z0min(max(z0_, 0.00001));
    scalar Ustar = Ufrictional();

    return ((Ustar/Kappa_)*log((z_ + z0min)/z0min));
}

void Foam::Function1Types::RHWindProfile::read(const dictionary& coeffs)
{
    Cmu_ = coeffs.lookupOrDefault<scalar>("Cmu", 0.09);
    Kappa_ = coeffs.lookupOrDefault<scalar>("Kappa", 0.41);
    Uref_ = coeffs.lookup<scalar>("Uref");
    Href_ = coeffs.lookup<scalar>("Href");
    z0_ = coeffs.lookup<scalar>("z0");
}

void Foam::Function1Types::RHWindProfile::writeData(Ostream& os) const
{
    Function1<scalar>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os.beginBlock(word(this->name() + "Coeffs"));
    os.writeEntry("Cmu", Cmu_);
    os.writeEntry("Kappa", Kappa_);
    os.writeEntry("Uref", Uref_);
    os.writeEntry("Href", Href_);
    os.writeEntry("z0", z0_);
    os.endBlock();
}

// ************************************************************************* //
