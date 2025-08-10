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
    uint64

Description
    64bit uinteger

SourceFiles
    uint64.C
    uint64IO.C

\*---------------------------------------------------------------------------*/

#ifndef uint64_H
#define uint64_H

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

//- Return a word representation of a uint64
inline word name(const uint64_t val)
{
    // no stripping required
    return word(std::to_string(val), false);
}


//- Return a word representation of a uint64_t, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const char* fmt, const uint64_t);


//- Return a word representation of a uint64_t, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const std::string& fmt, const uint64_t);


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

uint64_t readUint64(Istream&);
bool read(const char*, uint64_t&);
Istream& operator>>(Istream&, uint64_t&);
Ostream& operator<<(Ostream&, const uint64_t);

//- Template specialization for pTraits<uint64_t>
template<>
class pTraits<uint64_t>
{
    uint64_t p_;

public:

    //- Component type
    typedef uint64_t cmptType;


    // Member constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Rank of uint64_t is 0
        static const direction rank = 0;

        //- Number of components in uint64_t is 1
        static const direction nComponents = 1;


    // Static data members

        static const char* const typeName;
        static const char* const componentNames[];
        static const uint64_t zero;
        static const uint64_t one;
        static const uint64_t min;
        static const uint64_t max;
        static const uint64_t rootMax;
        static const uint64_t rootMin;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const uint64_t&);

        //- Construct from Istream
        pTraits(Istream&);


    // Member Functions

        //- Access to the uint64_t value
        operator uint64_t() const
        {
            return p_;
        }

        //- Access to the uint64_t value
        operator uint64_t&()
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
