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

Class
    Foam::functionEntries::includeIfPresentEntry

Description
    Specify a file to include if it exists. Expects a single string to follow.

    The \c \#includeIfPresent directive is similar to the \c \#include
    directive, but does not generate an error if the file does not exist.

See also
    Foam::functionEntries::includeEntry

SourceFiles
    includeIfPresentEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeIfPresentEntry_H
#define includeIfPresentEntry_H

#include "db/dictionary/functionEntries/includeEntry/includeEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                    Class includeIfPresentEntry Declaration
\*---------------------------------------------------------------------------*/

class includeIfPresentEntry
:
    public includeEntry
{
public:

    //- Include (if present) in a sub-dict context
    static bool execute(dictionary& parentDict, Istream& is);

    //- Include (if present) in a primitiveEntry context
    static bool execute
    (
        const dictionary& parentDict,
        primitiveEntry& entry,
        Istream& is
    );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
