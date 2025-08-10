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
    uint32

Description
    32bit uinteger

SourceFiles
    uint32.C
    uint32IO.C

\*---------------------------------------------------------------------------*/

#ifndef uint32_H
#define uint32_H

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

//- Return a word representation of a uint32
inline word name(const uint32_t val)
{
    // no stripping required
    return word(std::to_string(val), false);
}


//- Return a word representation of a uint32, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const char* fmt, const uint32_t);


//- Return a word representation of a uint32, using printf-style formatter.
//  The representation is not checked for valid word characters.
word name(const std::string& fmt, const uint32_t);


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

uint32_t readUint32(Istream&);
bool read(const char*, uint32_t&);
Istream& operator>>(Istream&, uint32_t&);
Ostream& operator<<(Ostream&, const uint32_t);

//- Template specialization for pTraits<uint32_t>
template<>
class pTraits<uint32_t>
{
    uint32_t p_;

public:

    //- Component type
    typedef uint32_t cmptType;


    // Member constants

        //- Dimensionality of space
        static const direction dim = 3;

        //- Rank of uint32_t is 0
        static const direction rank = 0;

        //- Number of components in uint32_t is 1
        static const direction nComponents = 1;


    // Static data members

        static const char* const typeName;
        static const char* const componentNames[];
        static const uint32_t zero;
        static const uint32_t one;
        static const uint32_t min;
        static const uint32_t max;
        static const uint32_t rootMax;
        static const uint32_t rootMin;


    // Constructors

        //- Construct from primitive
        explicit pTraits(const uint32_t&);

        //- Construct from Istream
        pTraits(Istream&);


    // Member Functions

        //- Access to the uint32_t value
        operator uint32_t() const
        {
            return p_;
        }

        //- Access to the uint32_t value
        operator uint32_t&()
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
