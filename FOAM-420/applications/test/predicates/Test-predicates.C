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

Application
    Test-predicates

Description
    Simple tests using predicates

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/IOstreams.H"
#include "primitives/ints/lists/labelList.H"
#include "primitives/strings/lists/wordList.H"
#include "primitives/predicates/predicates.H"
#include "containers/Lists/ListOps/FlatOutput.H"
#include "regExp.H"

using namespace Foam;


template<class ListType, class UnaryPredicate>
label printMatching(const ListType& list, const UnaryPredicate& pred)
{
    label count = 0;

    Info<< "(";

    for (const auto& val : list)
    {
        if (pred(val))
        {
            if (count) Info<< ' ';
            Info<< val;
            ++count;
        }
    }

    Info<< ") => " << count << nl;

    return count;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    wordList words
    {
        "abc",
        "def",
        "hij",
        "abc_",
        "def_",
        "hij_",
    };

    labelRange range(-10, 40);
    labelList values(range.begin(), range.end());

    Info<<"words:  " << flatOutput(words) << endl;
    Info<<"values: " << flatOutput(values)  << endl;

    regExp matcher(".*_.*");

    Info<<"With '_': ";
    printMatching(words, matcher);

    Info<<"All: ";
    printMatching(words, predicates::always());

    Info<<"None: ";
    printMatching(words, predicates::never());

    Info<<"Neg values: ";
    printMatching(values, [](const label v) { return v < 0; });

    Info<<"Even values: ";
    printMatching(values, [](const label v) { return !(v % 2); });

    Info<<"All: ";
    printMatching(values, predicates::always());

    Info<<"None: ";
    printMatching(values, predicates::never());

    return 0;
}

// ************************************************************************* //
