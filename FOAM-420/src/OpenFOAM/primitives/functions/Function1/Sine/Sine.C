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
    (c) 2016 OpenCFD Ltd.
    (c) 2016-2017 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Sine/Sine.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Sine<Type>::read(const dictionary& coeffs)
{
    t0_ = coeffs.lookupOrDefault<scalar>("t0", 0);
    amplitude_ = Function1<scalar>::New("amplitude", coeffs);
    frequency_ = Function1<scalar>::New("frequency", coeffs);
    scale_ = Function1<Type>::New("scale", coeffs);
    level_ = Function1<Type>::New("level", coeffs);
}


template<class Type>
Foam::Function1Types::Sine<Type>::Sine
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName)
{
    read(dict);
}


template<class Type>
Foam::Function1Types::Sine<Type>::Sine(const Sine<Type>& se)
:
    Function1<Type>(se),
    t0_(se.t0_),
    amplitude_(se.amplitude_, false),
    frequency_(se.frequency_, false),
    scale_(se.scale_, false),
    level_(se.level_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Sine<Type>::~Sine()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::Sine<Type>::value(const scalar t) const
{
    return
        amplitude_->value(t)
       *sin(constant::mathematical::twoPi*frequency_->value(t)*(t - t0_))
       *scale_->value(t)
      + level_->value(t);
}


template<class Type>
Type Foam::Function1Types::Sine<Type>::derivative(const scalar t) const
{
    const scalar twoPi = constant::mathematical::twoPi;

    // Values and derivitives
    const scalar f = frequency_->value(t);
    const scalar fDer = frequency_->derivative(t);
    const scalar a = amplitude_->value(t);
    const scalar aDer = amplitude_->derivative(t);
    const Type scale = scale_->value(t);
    const Type scaleDer = scale_->derivative(t);
    const Type levelDer = level_->derivative(t);
    const scalar tMinT0 = (t - t0_);

    // Using chain rule twice
    return
        (aDer*scale + a*scaleDer)*sin(twoPi*f*tMinT0)
      + a*scale*twoPi*(fDer*tMinT0 + f)*cos(twoPi*f*tMinT0)
      + levelDer;
}


template<class Type>
void Foam::Function1Types::Sine<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));

    os.writeEntry("t0", t0_);
    amplitude_->writeData(os);
    frequency_->writeData(os);
    scale_->writeData(os);
    level_->writeData(os);

    os.endBlock() << flush;
}


// ************************************************************************* //
