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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/dictionary/functionEntries/inputModeEntry/inputModeEntry.H"
#include "db/dictionary/dictionary.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::inputModeEntry::typeName
(
    Foam::functionEntries::inputModeEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include inputModeEntries
int Foam::functionEntries::inputModeEntry::debug(0);

Foam::functionEntries::inputModeEntry::inputMode
    Foam::functionEntries::inputModeEntry::mode_(MERGE);

namespace Foam
{
namespace functionEntries
{
    addToMemberFunctionSelectionTable
    (
        functionEntry,
        inputModeEntry,
        execute,
        dictionaryIstream
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// we could combine this into execute() directly, but leave it here for now
void Foam::functionEntries::inputModeEntry::setMode(Istream& is)
{
    clear();

    const word mode(is);
    if (mode == "merge" || mode == "default")
    {
        mode_ = MERGE;
    }
    else if (mode == "overwrite")
    {
        mode_ = OVERWRITE;
    }
    else if (mode == "protect")
    {
        mode_ = PROTECT;
    }
    else if (mode == "warn")
    {
        mode_ = WARN;
    }
    else if (mode == "error")
    {
        mode_ = ERROR;
    }
    else
    {
        WarningInFunction
            << "unsupported input mode '" << mode
            << "' ... defaulting to 'merge'"
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::inputModeEntry::execute
(
    dictionary& unused,
    Istream& is
)
{
    setMode(is);
    return true;
}


void Foam::functionEntries::inputModeEntry::clear()
{
    mode_ = MERGE;
}


bool Foam::functionEntries::inputModeEntry::merge()
{
    return mode_ == MERGE;
}


bool Foam::functionEntries::inputModeEntry::overwrite()
{
    return mode_ == OVERWRITE;
}


bool Foam::functionEntries::inputModeEntry::protect()
{
    return mode_ == PROTECT;
}

bool Foam::functionEntries::inputModeEntry::error()
{
    return mode_ == ERROR;
}


// ************************************************************************* //
