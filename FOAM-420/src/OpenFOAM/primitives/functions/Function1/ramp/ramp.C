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
    (c) 2017 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/ramp/ramp.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

void Foam::Function1Types::ramp::read(const dictionary& coeffs)
{
    start_ = coeffs.lookupOrDefault<scalar>("start", 0);
    duration_ = coeffs.lookup<scalar>("duration");
}


Foam::Function1Types::ramp::ramp
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<scalar>(entryName)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Function1Types::ramp::~ramp()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::Function1Types::ramp::linearRampDerivative(scalar t) const
{
    const scalar rampValue = ramp::linearRamp(t);

    // Make sure we can check values
    const scalar valMax = rampValue + doubleScalarVSMALL;
    const scalar valMin = rampValue - doubleScalarVSMALL;

    // Check if the limits were reached (derivative of const function is zero)
    if (valMax >= 1.0 || valMin <= 0)
    {
        return 0.0;
    }

    return 1.0/duration_;
}


void Foam::Function1Types::ramp::writeData(Ostream& os) const
{
    Function1<scalar>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os.beginBlock(word(this->name() + "Coeffs"));
    os.writeEntry("start", start_);
    os.writeEntry("duration", duration_);
    os.endBlock();
}


// ************************************************************************* //
