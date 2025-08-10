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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2007 OpenCFD Ltd.

Description
    Reads specified dictionary from constant or system and writes it out again

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "db/dictionary/dictionary.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/dictionary/dictionaryEntry/dictionaryEntry.H"
#include "db/dictionary/entry/entry.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "writeIODictToFile.H"

int main(int argc, char *argv[])
{
    argList::validArgs.append("dictionary");

#include "include/setRootCase.H"
#include "include/createTime.H"
    runTime.functionObjects().off();

    //will simply check constant and system for specified dicts
    //read and then write them

    word dictName(args.args()[1]);


    IOobject systemDictHeader
    (
        dictName,
        runTime.caseSystem(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
    IOobject constantDictHeader
    (
        dictName,
        runTime.caseConstant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (systemDictHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(systemDictHeader);
        Info<< "Writing " << dictName << "..." << endl;
        writeDictFile(dict);

        Info<< "End\n" << endl;
        return 0;
    }
    else if (constantDictHeader.typeHeaderOk<IOdictionary>(true))
    {
        IOdictionary dict(constantDictHeader);
        Info<< "Writing " << dictName << "..." << endl;
        writeDictFile(dict);

        Info<< "End\n" << endl;
        return 0;
    }
    else
    {
        Info<< "Could not find dictionary "
             << dictName
             << " in system or constant directories."
             << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
