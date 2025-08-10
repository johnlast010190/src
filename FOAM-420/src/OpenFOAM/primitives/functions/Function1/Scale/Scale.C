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

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Scale/Scale.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
void Foam::Function1Types::Scale<Type>::read(const dictionary& coeffs)
{
    scale_ = Function1<scalar>::New("scale", coeffs);
    value_ = Function1<Type>::New("value", coeffs);
}


template<class Type>
Foam::Function1Types::Scale<Type>::Scale
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
Foam::Function1Types::Scale<Type>::Scale(const Scale<Type>& se)
:
    Function1<Type>(se),
    scale_(se.scale_, false),
    value_(se.value_, false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Scale<Type>::~Scale()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::Scale<Type>::value(const scalar x) const
{
    return scale_->value(x)*value_->value(x);
}


template<class Type>
Type Foam::Function1Types::Scale<Type>::derivative
(
    const scalar x
) const
{
    // Derivative chain rule
    return
        scale_->derivative(x)*value_->value(x)
      + scale_->value(x)*value_->derivative(x);
}



template<class Type>
void Foam::Function1Types::Scale<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os  << token::END_STATEMENT << nl;
    os.beginBlock(word(this->name() + "Coeffs"));
    scale_->writeData(os);
    value_->writeData(os);
    os.endBlock();
}


// ************************************************************************* //
