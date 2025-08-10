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
    (c) 2017 OpenCFD Ltd.

Description
    Test miscellaneous C++ templates/functionality.

\*---------------------------------------------------------------------------*/

#include "primitives/strings/string/string.H"
#include "db/IOstreams/IOstreams.H"
#include "containers/Lists/UList/UList.H"
#include "containers/HashTables/HashSet/HashSet.H"

#include <typeinfo>
#include <type_traits>
#include <utility>

using namespace Foam;

// Macros to stringify macro contents.
#define STRINGIFY(content)      #content
#define STRING_QUOTE(input)     STRINGIFY(input)

#define PRINT_TYPEID(arg)       \
    Info<< typeid(arg).name() << " <= typeid of " << STRING_QUOTE(arg) << nl



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "various declaration types" << nl << nl;

    PRINT_TYPEID(label);
    PRINT_TYPEID(decltype(UList<label>::value_type()));
    PRINT_TYPEID(decltype(std::declval<UList<label>>().cbegin()));
    PRINT_TYPEID(decltype(*(std::declval<UList<label>>().cbegin())));
    Info<< nl;

    PRINT_TYPEID(decltype(HashTable<label>::key_type()));
    PRINT_TYPEID(decltype(HashTable<label>::value_type()));
    // Not yet: PRINT_TYPEID(decltype(HashTable<label>::mapped_type()));
    PRINT_TYPEID(decltype(std::declval<HashTable<label>>().begin()));
    PRINT_TYPEID(decltype(std::declval<const HashTable<label>>().begin()));
    PRINT_TYPEID(decltype(*(std::declval<HashTable<label>>().begin())));
    PRINT_TYPEID(decltype(*(std::declval<const HashTable<label>>().begin())));

    PRINT_TYPEID(decltype(std::declval<const HashTable<label>>().keys()));
    Info<< nl;

    PRINT_TYPEID(decltype(HashSet<label>::key_type()));
    PRINT_TYPEID(decltype(HashSet<label>::value_type()));
    // Not yet: PRINT_TYPEID(decltype(HashSet<label>::mapped_type()));
    PRINT_TYPEID(decltype(std::declval<HashSet<label>>().begin()));
    PRINT_TYPEID(decltype(std::declval<const HashSet<label>>().begin()));
    PRINT_TYPEID(decltype(*(std::declval<HashSet<label>>().begin())));
    PRINT_TYPEID(decltype(*(std::declval<const HashSet<label>>().begin())));
    Info<< nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
