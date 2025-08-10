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
    (c) 2016 OpenCFD Ltd.

Class
    Foam::ensight::FileName

Description
    Specification of a valid Ensight file-name.

    Spaces must be quoted,
    no '*' wildcards, no '%' (structured block continuation).

    Overall line length within case file is limited to 1024, but this is not
    yet addresssed.

\*---------------------------------------------------------------------------*/

#ifndef ensightFileName_H
#define ensightFileName_H

#include "primitives/strings/fileName/fileName.H"
#include "primitives/strings/word/word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace ensight
{

/*---------------------------------------------------------------------------*\
                      Class ensight::FileName Declaration
\*---------------------------------------------------------------------------*/

class FileName
:
    public fileName
{
    // Private Member Functions

        //- Strip invalid characters
        inline void stripInvalid();

public:

    // Constructors

        //- Construct as copy
        inline explicit FileName(const FileName&);

        //- Construct as copy of character array
        inline explicit FileName(const char*);

        //- Construct as copy of std::string
        inline explicit FileName(const std::string&);


    // Member functions

        //- Is this character valid for an ensight file-name
        inline static bool valid(char);


    // Member operators

    // Assignment (disabled)

        void operator=(const fileName&) = delete;
        void operator=(const word&) = delete;
        void operator=(const string&) = delete;
        void operator=(const std::string&) = delete;
        void operator=(const char*) = delete;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace ensight
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ensight/name/ensightFileNameI.H"

#endif

// ************************************************************************* //
