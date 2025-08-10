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

InClass
    Foam::contiguous

Description
    Template function to specify if the data of a type are contiguous.

    The default function specifies that data are not contiguous.
    This is specialised for the types (eg, primitives) with contiguous data.

\*---------------------------------------------------------------------------*/

#ifndef contiguous_H
#define contiguous_H

#include "primitives/ints/int/int.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Forward declaration of friend functions and operators
template<typename T, unsigned Size> class FixedList;
template<typename T> class Pair;
template<typename Type1, class Type2> class Tuple2;

template<typename T>
constexpr bool contiguous();

// TODO: Figure out why the wheels fall off in the serialiser when we simply do
//       `std::is_trivially_copyable_V<T>` here. That ought to work, but I thiink
//       something somewhere is checking `contiguous` when it should be checking
//       ostream::format()==BINARY to determine whether or not to use binary mode.
template<typename T>
constexpr bool is_contiguous = std::is_fundamental_v<T> || std::is_enum_v<T>;


// Pairs of contiguous things are contiguous.
template<typename T>
constexpr bool is_contiguous<Pair<T>> = contiguous<T>();

// Landmine:
// TODO: Why does uncommenting this break everything? Probably because the serialiser code is nuts.
//template<typename T, typename Q>
//constexpr bool is_contiguous<Tuple2<T, Q>> = contiguous<T>() && contiguous<Q>();

// Tuples of contiguous things are contiguous.
template<typename... T>
constexpr bool is_contiguous<std::tuple<T...>> = (contiguous<T>() && ...);

// FixedList of contiguous things are contiguous.
template<typename T, unsigned N>
constexpr bool is_contiguous<FixedList<T, N>> = contiguous<T>();

// Existing code specialises _this_, instead of the above (hence why the above keeps
// calling this function). Existing code also reads this.
template<class T>
constexpr bool contiguous() {
    return is_contiguous<T>;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
