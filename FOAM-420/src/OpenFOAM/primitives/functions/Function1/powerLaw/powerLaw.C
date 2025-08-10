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

#include "primitives/functions/Function1/powerLaw/powerLaw.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace Function1Types
{
    makeScalarFunction1(powerLaw);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::Function1Types::powerLaw::read(const dictionary& coeffs)
{
    a_ = coeffs.lookup<scalar>("a");
    b_ = coeffs.lookup<scalar>("b");
    c_ = coeffs.lookup<scalar>("c");
    d_ = coeffs.lookup<scalar>("d");
    maxValue_ = coeffs.lookupOrDefault("maxValue", floatScalarVGREAT);
    minValue_ = coeffs.lookupOrDefault("minValue", -floatScalarVGREAT);
    clip_ = coeffs.lookupOrDefault("clip", true);
}


Foam::Function1Types::powerLaw::powerLaw
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<scalar>(entryName)
{
    read(dict);
}


Foam::Function1Types::powerLaw::powerLaw(const powerLaw& se)
:
    Function1<scalar>(se),
    a_(se.a_),
    b_(se.b_),
    c_(se.c_),
    d_(se.d_),
    maxValue_(se.maxValue_),
    minValue_(se.minValue_),
    clip_(se.clip_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::powerLaw::~powerLaw()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::powerLaw::value(const scalar x) const
{
    if ((x - b_) < 0.0 && mag(c_) < 1.0)
    {
        if (clip_)
        {
            return (c_ < 0.0) ? maxValue_ : max(d_, minValue_);
        }
        else
        {
            FatalErrorInFunction
                << "The value " << x
                << " is outside of validity for: " << type()
                << exit(FatalError);
        }
    }

    return max(min(a_*pow(x - b_, c_) + d_, maxValue_), minValue_);
}


Foam::scalar Foam::Function1Types::powerLaw::derivative(const scalar x) const
{
    if ((x - b_) < 0.0 && mag(c_ - 1.0) < 1.0)
    {
        FatalErrorInFunction
            << "The derivative at point " << x
            << " is outside of validity for: " << type()
            << exit(FatalError);
    }

    // Check for the clip => constant function
    const scalar val = value(x);
    if (val == maxValue_ || val == minValue_)
    {
        return 0.0;
    }

    return a_*c_*pow(x - b_, c_ - 1);
}


void Foam::Function1Types::powerLaw::writeData(Ostream& os) const
{
    Function1<scalar>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));

    os.writeEntry("a", a_);
    os.writeEntry("b", b_);
    os.writeEntry("c", c_);
    os.writeEntry("d", d_);

    if (maxValue_ != floatScalarVGREAT)
    {
        os.writeEntry("maxValue", maxValue_);
    }

    if (maxValue_ != -floatScalarVGREAT)
    {
        os.writeEntry("minValue", maxValue_);
    }
    os.endBlock() << flush;
}


// ************************************************************************* //
