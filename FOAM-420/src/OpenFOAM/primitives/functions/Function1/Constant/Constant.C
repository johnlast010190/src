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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015 OpenCFD Ltd.
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "primitives/functions/Function1/Constant/Constant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Constant<Type>::Constant
(
    const word& entryName,
    const Type& value
)
:
    Function1<Type>(entryName),
    value_(value)
{}


template<class Type>
Foam::Function1Types::Constant<Type>::Constant
(
    const word& entryName,
    const dictionary& dict
)
:
    Function1<Type>(entryName),
    value_(Zero)
{
    Istream& is(dict.lookup(entryName));
    word entryType(is);
    is  >> value_;
}


template<class Type>
Foam::Function1Types::Constant<Type>::Constant
(
    const word& entryName,
    Istream& is
)
:
    Function1<Type>(entryName),
    value_(pTraits<Type>(is))
{}


template<class Type>
Foam::Function1Types::Constant<Type>::Constant(const Constant<Type>& cnst)
:
    Function1<Type>(cnst),
    value_(cnst.value_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Constant<Type>::~Constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::Constant<Type>::value(const scalar x) const
{
    return value_;
}


template<class Type>
Type Foam::Function1Types::Constant<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    return (x2 - x1)*value_;
}


template<class Type>
Type Foam::Function1Types::Constant<Type>::integrateYoverX
(
    const scalar x1,
    const scalar x2
) const
{
    return (log(x2) - log(x1))*value_;
}


template<class Type>
Type Foam::Function1Types::Constant<Type>::derivative
(
    const scalar x
) const
{
    return pTraits<Type>::zero;
}


template<class Type>
void Foam::Function1Types::Constant<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);

    os  << token::SPACE << value_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
