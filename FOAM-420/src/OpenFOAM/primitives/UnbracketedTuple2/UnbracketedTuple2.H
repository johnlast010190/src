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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::UnbracketedTuple2

Description
    Same as Tuple2, but read and stored without list delimiters (brackets)

See also
    Foam::Tuple2

\*---------------------------------------------------------------------------*/

#ifndef UnbracketedTuple2_H
#define UnbracketedTuple2_H

#include "primitives/Tuple2/Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type1, class Type2>
class UnbracketedTuple2;

template<class Type1, class Type2>
inline Istream& operator>>(Istream& is, UnbracketedTuple2<Type1, Type2>& t2);

template<class Type1, class Type2>
inline Ostream& operator<<(Ostream& os, const UnbracketedTuple2<Type1, Type2>& t2);


/*---------------------------------------------------------------------------*\
                    class UnbracketedTuple2 Declaration
\*---------------------------------------------------------------------------*/

template<class Type1, class Type2>
class UnbracketedTuple2
:
    public Tuple2<Type1, Type2>
{

public:

    // Constructors

        //- Null constructor for lists and hashes
        inline UnbracketedTuple2()
        :
            Tuple2<Type1, Type2>()
        {}

        //- Construct from components
        inline UnbracketedTuple2(const Type1& f, const Type2& s)
        :
            Tuple2<Type1, Type2>(f, s)
        {}

        //- Construct from Istream
        inline UnbracketedTuple2(Istream& is)
        :
            Tuple2<Type1, Type2>(is)
        {}


    // IOstream operators

        //- Read UnbracketedTuple2 from Istream, discarding contents of existing
        friend Istream& operator>> <Type1, Type2>
        (
            Istream& is,
            UnbracketedTuple2<Type1, Type2>& t2
        );

        // Write UnbracketedTuple2 to Ostream.
        friend Ostream& operator<< <Type1, Type2>
        (
            Ostream& os,
            const UnbracketedTuple2<Type1, Type2>& t2
        );
};


template<class Type1, class Type2>
inline Istream& operator>>(Istream& is, UnbracketedTuple2<Type1, Type2>& t2)
{
    is >> t2.first() >> t2.second();
    is.check(FUNCTION_NAME);
    return is;
}


template<class Type1, class Type2>
inline Ostream& operator<<(Ostream& os, const UnbracketedTuple2<Type1, Type2>& t2)
{
    os  << t2.first() << token::SPACE << t2.second();
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
