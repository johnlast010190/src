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
    (c) 2020 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "Constant2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Constant<Type>::Constant
(
    const word& name,
    const Type& val
)
:
    FieldFunction2<Type, Constant<Type>>(name),
    value_(val)
{}


template<class Type>
Foam::Function2s::Constant<Type>::Constant
(
    const word& name,
    const dictionary& dict
)
:
    FieldFunction2<Type, Constant<Type>>(name),
    value_(Zero)
{
    if (!dict.found(name))
    {
        dict.lookup("value") >> value_;
    }
    else
    {
        Istream& is(dict.lookup(name));
        word entryType(is);
        if (is.eof())
        {
            dict.lookup("value") >> value_;
        }
        else
        {
            is  >> value_;
        }
    }
}


template<class Type>
Foam::Function2s::Constant<Type>::Constant
(
    const word& name,
    Istream& is
)
:
    FieldFunction2<Type, Constant<Type>>(name),
    value_(pTraits<Type>(is))
{}


template<class Type>
Foam::Function2s::Constant<Type>::Constant(const Constant<Type>& cnst)
:
    FieldFunction2<Type, Constant<Type>>(cnst),
    value_(cnst.value_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function2s::Constant<Type>::~Constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::Function2s::Constant<Type>::writeData(Ostream& os) const
{
    Function2<Type>::writeData(os);
    os  << token::SPACE << value_ << token::END_STATEMENT << nl;
}


// ************************************************************************* i/
