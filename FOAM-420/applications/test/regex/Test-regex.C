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

Description
    Tests for regular expressions

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/IOstreams.H"
#include "db/IOobject/IOobject.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "regExp.H"
#include "containers/Lists/List/List.H"
#include "primitives/Tuple2/Tuple2.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    List<Tuple2<string, string>> rawList(IFstream("testRegexps")());
    Info<< "Test expressions:" << rawList << endl;
    IOobject::writeDivider(Info) << endl;

    List<string> groups;

    // Report matches:
    forAll(rawList, elemI)
    {
        const string& pat = rawList[elemI].first();
        const string& str = rawList[elemI].second();
        regExp re(pat);

        Info<< str << " =~ m/" << pat.c_str() << "/ == ";

        if (re.match(str, groups))
        {
            Info<< "true";
            if (re.ngroups())
            {
                Info<< nl << "groups: " << groups;
            }
        }
        else
        {
            if (re.search(str))
            {
                Info<< " partial match";
            }
            else
            {
                Info<< "false";
            }
        }
        Info<< endl;
    }

    Info<< nl << "test regExp(const char*) ..." << endl;
    string me("Mark");

    // Handling of null strings
    if (regExp(nullptr).match(me))
    {
        Info<< "fail - matched: " << me << endl;
    }
    else
    {
        Info<< "pass - null pointer is no expression" << endl;
    }

    // Normal match
    if (regExp("[Mm]ar[ck]").match(me))
    {
        Info<< "pass - matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    // Match ignore case
    if (regExp("mar[ck]", true).match(me))
    {
        Info<< "pass - matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    // Embedded prefix for match ignore case
    if (regExp("(?i)mar[ck]").match(me))
    {
        Info<< "pass - matched: " << me << endl;
    }
    else
    {
        Info<< "no match" << endl;
    }

    // Handling of empty expression
    if (regExp("").match(me))
    {
        Info<< "fail - matched: " << me << endl;
    }
    else
    {
        Info<< "pass - no match on empty expression" << endl;
    }

    // Embedded prefix - but expression is empty
    if (regExp("(?i)").match(me))
    {
        Info<< "fail - matched: " << me << endl;
    }
    else
    {
        Info<< "pass - no match on empty expression" << endl;
    }


    Info<< endl;

    return 0;
}


// ************************************************************************* //
