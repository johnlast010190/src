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
    int32

Description
    32bit integer

SourceFiles
    int32.C
    int32IO.C

\*---------------------------------------------------------------------------*/

#ifndef int32_H
#define int32_H

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

//- Return a word representation of an int32
inline word name(const int32_t val)
{
    // no stripping required
    return word(std::to_string(val), false);
}


//- Return a word representation of an int32, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const char* fmt, const int32_t);


//- Return a word representation of an int32, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const std::string&, const int32_t);


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

int32_t readInt32(Istream&);
bool read(const char*, int32_t&);
Istream& operator>>(Istream&, int32_t&);
Ostream& operator<<(Ostream&, const int32_t);

//  32-bit operating systems are no longer supported by FOAM
// On 32bit OSs long is not unambiguously int32_t (or int64_t) causing problems
// for IO operator resolution.
// This problem is avoided by explicitly defining the following operators:
// #if FOAM_ARCH_OPTION == 32
//     Istream& operator>>(Istream&, long&);
//     Ostream& operator<<(Ostream&, const long);
// #endif


//- Template specialization for pTraits<int32_t>
template<>
class pTraits<int32_t>
{
    int32_t p_;

public:

    //- Component type
    typedef int32_t cmptType;


    // Member constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of int32_t is 0
        static constexpr direction rank = 0;

        //- Number of components in int32_t is 1
        static constexpr direction nComponents = 1;


    // Static data members

        static const char* const typeName;
        static const char* const componentNames[];
        static constexpr int32_t zero = 0;
        static constexpr int32_t one = 1;
        static constexpr int32_t min = INT32_MIN;
        static constexpr int32_t max = INT32_MAX;
        static constexpr int32_t rootMax = pTraits<int32_t>::max;
        static constexpr int32_t rootMin = pTraits<int32_t>::min;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const int32_t& p):
            p_(p)
        {}

        //- Construct from Istream
        pTraits(Istream& is) {
            is >> p_;
        }


    // Member Functions

        //- Access to the int32_t value
        operator int32_t() const
        {
            return p_;
        }

        //- Access to the int value
        operator int32_t&()
        {
            return p_;
        }
};


inline int32_t mag(const int32_t l)
{
    return ::abs(l);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
