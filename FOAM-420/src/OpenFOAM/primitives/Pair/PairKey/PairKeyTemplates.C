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
    (c) 2019 Esi Ltd.
    (c) 2014-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "primitives/Pair/PairKey/PairKey.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PairKey<Type>::hash::hash()
{}


template<class Type>
Foam::PairKey<Type>::PairKey()
{}


template<class Type>
Foam::PairKey<Type>::PairKey
(
    const word& name1,
    const word& name2,
    const bool ordered
)
:
    Pair<word>(name1, name2),
    ordered_(ordered)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::PairKey<Type>::~PairKey()
{}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class Type>
Foam::label Foam::PairKey<Type>::hash::operator()
(
    const PairKey<Type>& key
) const
{
    if (key.ordered_)
    {
        return
            typename Type::hash()
            (
                key.first(),
                typename Type::hash()(key.second())
            );
    }
    else
    {
        return
            typename Type::hash()(key.first())
          + typename Type::hash()(key.second());
    }
}


// * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * * //

template<class Type>
bool Foam::operator==
(
    const PairKey<Type>& a,
    const PairKey<Type>& b
)
{
    const label c = Pair<Type>::compare(a,b);

    return
        (a.ordered_ == b.ordered_)
     && (
            (a.ordered_ && (c == 1))
         || (!a.ordered_ && (c != 0))
        );
}


template<class Type>
bool Foam::operator!=
(
    const PairKey<Type>& a,
    const PairKey<Type>& b
)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * Istream Operator  * * * * * * * * * * * * * * //

template<class Type>
Foam::Istream& Foam::operator>>(Istream& is, PairKey<Type>& key)
{
    const FixedList<word, 3> temp(is);

    key.first() = temp[0];

    if (temp[1] == "and")
    {
        key.ordered_ = false;
    }
    else if (temp[1] == "in")
    {
        key.ordered_ = true;
    }
    else
    {
        FatalErrorInFunction
            << "Pair type is not recognised. "
            << temp
            << "Use (name1 in name2) for an ordered"
            << "pair, or (name1 and name2) for an unordered pair."
            << exit(FatalError);
    }

    key.second() = temp[2];

    return is;
}


// * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const PairKey<Type>& key)
{
    os  << token::BEGIN_LIST
        << key.first()
        << token::SPACE
        << (key.ordered_ ? "in" : "and")
        << token::SPACE
        << key.second()
        << token::END_LIST;

    return os;
}


// ************************************************************************* //
