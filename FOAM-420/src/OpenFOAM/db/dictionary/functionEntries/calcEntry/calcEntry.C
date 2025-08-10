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
    (c) 2011-2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/dictionary/functionEntries/calcEntry/calcEntry.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"
#include "db/dictionary/dictionary.H"
#include "db/dynamicLibrary/dynamicCode/dynamicCode.H"
#include "db/dictionary/functionEntries/codeStream/codeStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{
    defineTypeNameAndDebug(calcEntry, 0);

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        calcEntry,
        execute,
        dictionaryIstream
    );

    addToMemberFunctionSelectionTable
    (
        functionEntry,
        calcEntry,
        execute,
        primitiveEntryIstream
    );

}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::calcEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& thisEntry,
    Istream& is
)
{
    Info<< "Using #calcEntry at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::calcEntry::execute(..)",
        parentDict
    );

    // Read string
    string s(is);
    // Make sure we stop this entry
    //is.putBack(token(token::END_STATEMENT, is.lineNumber()));

    // Construct codeDict for codeStream
    // must reference parent for stringOps::expand to work nicely.
    dictionary codeSubDict;
    codeSubDict.add("code", "os << (" + s + ");");
    dictionary codeDict(parentDict, codeSubDict);

    codeStream::streamingFunctionType function = codeStream::getFunction
    (
        parentDict,
        codeDict
    );

    // use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // get the entry from this stream
    IStringStream resultStream(os.str());
    thisEntry.read(parentDict, resultStream);

    return true;
}


bool Foam::functionEntries::calcEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    Info<< "Using #calcEntry at line " << is.lineNumber()
        << " in file " <<  parentDict.name() << endl;

    dynamicCode::checkSecurity
    (
        "functionEntries::calcEntry::execute(..)",
        parentDict
    );

    // Read string
    string s(is);
    // Make sure we stop this entry
    //is.putBack(token(token::END_STATEMENT, is.lineNumber()));

    // Construct codeDict for codeStream
    // must reference parent for stringOps::expand to work nicely.
    dictionary codeSubDict;
    codeSubDict.add("code", "os << (" + s + ");");
    dictionary codeDict(parentDict, codeSubDict);

    codeStream::streamingFunctionType function = codeStream::getFunction
    (
        parentDict,
        codeDict
    );

    // use function to write stream
    OStringStream os(is.format());
    (*function)(os, parentDict);

    // get the entry from this stream
    IStringStream resultStream(os.str());
    parentDict.read(resultStream);

    return true;
}


// ************************************************************************* //
