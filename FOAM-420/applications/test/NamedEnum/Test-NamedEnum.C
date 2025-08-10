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

Description

\*---------------------------------------------------------------------------*/

#include "primitives/enums/NamedEnum.H"
#include "primitives/enums/Enum.H"
#include "db/IOstreams/IOstreams.H"

using namespace Foam;

class namedEnumTest
{
public:

    enum class option
    {
        A,
        B,
        C,
        D
    };

    enum class otherOption
    {
        A,
        B,
        C,
        D
    };

    static const Foam::NamedEnum<option, 4> optionNamed;

    static const Foam::Enum<otherOption> optionEnum;

    static const Foam::Enum<option> optionEnum2;
};


template<>
const char* Foam::NamedEnum<namedEnumTest::option, 4>::names[] =
{
    "a",
    "b",
    "c",
    "d",
};

const Foam::NamedEnum<namedEnumTest::option, 4> namedEnumTest::optionNamed;

const Foam::Enum<namedEnumTest::otherOption> namedEnumTest::optionEnum
{
    { namedEnumTest::otherOption::A, "a" },
    { namedEnumTest::otherOption::B, "b" },
    { namedEnumTest::otherOption::C, "c" },
    { namedEnumTest::otherOption::D, "d" },
};


const Foam::Enum<namedEnumTest::option> namedEnumTest::optionEnum2
(
    namedEnumTest::option::C,
    { "c", "d" }
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<<"NamedEnum: " << namedEnumTest::optionNamed << nl;
    Info<<"Enum: " << namedEnumTest::optionEnum << nl;
    Info<<"Enum: " << namedEnumTest::optionEnum2 << nl;

    dictionary testDict;
    testDict.add("lookup1", "c");

    Info<< nl
        << int(namedEnumTest::optionNamed["a"]) << nl
        << namedEnumTest::optionNamed[namedEnumTest::option::A] << nl;

    Info<< nl
        << int(namedEnumTest::optionEnum["a"]) << nl
        << namedEnumTest::optionEnum[namedEnumTest::otherOption::A] << nl;

    Info<< "--- test dictionary lookup ---" << endl;
    {
        Info<< "dict: " << testDict << endl;

        Info<< "got: "
            <<  int
                (
                    namedEnumTest::optionNamed.lookupOrDefault
                    (
                        "notFound",
                        testDict,
                        namedEnumTest::option::A
                    )
                )
            << nl;

        Info<< "got: "
            <<  int
                (
                    namedEnumTest::optionNamed.lookupOrDefault
                    (
                        "lookup1",
                        testDict,
                        namedEnumTest::option::A
                    )
                )
            << nl;

        Info<< "got: "
            <<  int
                (
                    namedEnumTest::optionEnum2.lookupOrDefault
                    (
                        "lookup1",
                        testDict,
                        namedEnumTest::option::A
                    )
                )
            << nl;
    }

    Info<< "--- test read ---" << endl;

    namedEnumTest::option dummy(namedEnumTest::optionNamed.read(Sin));
    Info<< namedEnumTest::optionNamed[dummy] << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
