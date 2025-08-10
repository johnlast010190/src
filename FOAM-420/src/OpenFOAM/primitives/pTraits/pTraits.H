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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::pTraits

Description
    Traits class for primitives.

    All primitives need a specialised version of this class. The
    specialised version will normally also require a conversion
    method.

\*---------------------------------------------------------------------------*/

#ifndef pTraits_H
#define pTraits_H

#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;

/*---------------------------------------------------------------------------*\
                           Class pTraits Declaration
\*---------------------------------------------------------------------------*/

template<typename T>
constexpr bool IsIStream() {
    constexpr bool out = std::is_same_v<Istream, T> || std::is_base_of_v<Istream, T>;
    static_assert(out, "IStream expected");
    return out;
}

template<class PrimitiveType>
class pTraits
:
    public PrimitiveType
{

public:

    // Constructors

        //- Construct from primitive
        pTraits(const PrimitiveType& p)
        :
            PrimitiveType(p)
        {}

        //- Construct from Istream
        template<typename T, std::enable_if_t<IsIStream<T>() && std::is_constructible_v<PrimitiveType, Istream&>, bool> = true>
        pTraits(T& is)
            : PrimitiveType(is)
        {}

        template<typename T, std::enable_if_t<IsIStream<T>() && !std::is_constructible_v<PrimitiveType, Istream&>, bool> = true>
        pTraits(T& is)
            : PrimitiveType()
        {
            is >> *(static_cast<std::remove_const_t<PrimitiveType>*>(this));
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
