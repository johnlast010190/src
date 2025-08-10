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

Application
    Test-codeStream

Description

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/IOstreams/IOstreams.H"
#include "db/IOobject/IOobject.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/dictionary/dictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("dict .. dictN");
    argList args(argc, argv, false, true);

    Info<< nl
        << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
        << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
        << endl;

    if (args.size() <= 1)
    {
        Info<<"specify dictionaries to test\n";
    }
    else
    {
        IOobject::writeDivider(Info);
        for (label argI=1; argI < args.size(); ++argI)
        {
            const string& dictFile = args[argI];
            IFstream is(dictFile);

            dictionary dict(is);

            Info<< dict << endl;
        }
    }

    return 0;
}


// ************************************************************************* //
