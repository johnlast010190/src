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
    (c) 2017 Esi Ltd.
    (c) 2015-2016 OpenFOAM Foundation

Class
    Foam::functionEntries::includeEtcEntry

Description
    Specify an etc file to include when reading dictionaries, expects a
    single string to follow.

    Searches for files from user/group/shipped directories.
    The search scheme allows for version-specific and
    version-independent files using the following hierarchy:
    - \b user settings:
      - ~/.FOAMcore/\<VERSION\>
      - ~/.FOAMcore/
    - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is set):
      - $FOAM_CUSTOM_DATA_DIR/\<VERSION\>
      - $FOAM_CUSTOM_DATA_DIR
    - \b group (site) settings (when $FOAM_CUSTOM_DATA_DIR is not set):
      - $FOAM_PROJECT_DIR/../site/\<VERSION\>
      - $FOAM_PROJECT_DIR/../site/
    - \b other (shipped) settings:
      - $FOAM_PROJECT_DIR/etc/

    An example of the \c \#includeEtc directive:
    \verbatim
        #includeEtc "etcFile"
    \endverbatim

    The usual expansion of environment variables and other constructs is
    retained.

See also
    findEtcFile, fileName, string::expand()

SourceFiles
    includeEtcEntry.C

\*---------------------------------------------------------------------------*/

#ifndef includeEtcEntry_H
#define includeEtcEntry_H

#include "db/dictionary/functionEntries/functionEntry/functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class includeEtcEntry Declaration
\*---------------------------------------------------------------------------*/

class includeEtcEntry
:
    public functionEntry
{
    //- Expand include fileName and return
    static fileName resolveFile(const fileName& f, const dictionary& dict);


public:

    // Static data members

        //- Report to stdout which file is included
        static bool log;


    // Member Functions

        //- Merge file into dictionary
        static bool mergeFile(dictionary& parentDict, fileName file);

        //- Include file in a sub-dict context
        static bool execute(dictionary& parentDict, Istream& is);

        //- Include file in a primitiveEntry context
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
