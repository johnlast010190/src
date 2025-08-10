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
    (c) 2014-2016 OpenFOAM Foundation

Primitive
    int64

Description
    64bit integer

SourceFiles
    int64.C
    int64IO.C

\*---------------------------------------------------------------------------*/

#ifndef int64_H
#define int64_H

#include <cstdint>
#include <climits>
#include <cstdlib>

#include "primitives/strings/word/word.H"
#include "primitives/pTraits/pTraits.H"
#include "primitives/direction/direction.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a word representation of an int64
inline word name(const int64_t val)
{
    // no stripping required
    return word(std::to_string(val), false);
}


//- Return a word representation of an int64, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const char* fmt, const int64_t);


//- Return a word representation of an int64, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const std::string& fmt, const int64_t);


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

int64_t readInt64(Istream&);
bool read(const char*, int64_t&);
Istream& operator>>(Istream&, int64_t&);
Ostream& operator<<(Ostream&, const int64_t);

//- Template specialization for pTraits<int64_t>
template<>
class pTraits<int64_t>
{
    int64_t p_;

public:

    //- Component type
    typedef int64_t cmptType;


    // Member constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Rank of int64_t is 0
        static const direction rank = 0;

        //- Number of components in int64_t is 1
        static const direction nComponents = 1;


    // Static data members

        static const char* const typeName;
        static const char* const componentNames[];
        static constexpr int64_t zero = 0;
        static constexpr int64_t one = 1;
        static constexpr int64_t min = INT64_MIN;
        static constexpr int64_t max = INT64_MAX;
        static constexpr int64_t rootMax = pTraits<int64_t>::max;
        static constexpr int64_t rootMin = pTraits<int64_t>::min;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const int64_t& p):
            p_(p)
        {}

        //- Construct from Istream
        pTraits(Istream& is) {
            is >> p_;
        }


    // Member Functions

        //- Access to the int64_t value
        operator int64_t() const
        {
            return p_;
        }

        //- Access to the int value
        operator int64_t&()
        {
            return p_;
        }
};


inline int64_t mag(const int64_t l)
{
    return ::labs(l);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
