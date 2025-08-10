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
    (c) 2011-2017 OpenFOAM Foundation

Application
    fileNameCleanTest

Description


\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "primitives/strings/fileName/fileName.H"
#include "containers/Lists/SubList/SubList.H"
#include "db/IOobject/IOobject.H"
#include "db/IOstreams/IOstreams.H"
#include "include/OSspecific.H"


using namespace Foam;

void printCleaning(fileName& pathName)
{
    Info<< "fileName = " << pathName << nl
        << "  path() = " << pathName.path() << nl
        << "  name() = " << pathName.name() << nl
        << "  joined = " << pathName.path()/pathName.name() << nl << nl;

    pathName.clean();

    Info<< "cleaned  = " << pathName << nl
        << "  path() = " << pathName.path() << nl
        << "  name() = " << pathName.name() << nl
        << "  joined = " << pathName.path()/pathName.name() << nl << nl;

    IOobject::writeDivider(Info);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::validArgs.insert("fileName .. fileNameN");
    argList::addOption("istream", "file", "test Istream values");

    argList args(argc, argv, false, true);

    if (args.size() <= 1 && args.options().empty())
    {
        args.printUsage();
    }

    fileName pathName;
    if (args.optionReadIfPresent("case", pathName))
    {
        Info<< nl
            << "-case" << nl
            << "path = " << args.path() << nl
            << "root = " << args.rootPath() << nl
            << "case = " << args.caseName() << nl
            << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
            << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
            << endl;

        printCleaning(pathName);
    }

    for (label argI=1; argI < args.size(); ++argI)
    {
        pathName = args[argI];
        printCleaning(pathName);
    }

    if (args.optionFound("istream"))
    {
        args.optionLookup("istream")() >> pathName;

        Info<< nl
            << "-case" << nl
            << "path = " << args.path() << nl
            << "root = " << args.rootPath() << nl
            << "case = " << args.caseName() << nl
            << "FOAM_CASE=" << getEnv("FOAM_CASE") << nl
            << "FOAM_CASENAME=" << getEnv("FOAM_CASENAME") << nl
            << endl;

        printCleaning(pathName);
    }

    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
