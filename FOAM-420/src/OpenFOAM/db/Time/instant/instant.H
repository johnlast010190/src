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
    (c) 2018-2022 OpenCFD Ltd
Class
    Foam::instant

Description
    An instant of time. Contains the time value and name.
    Uses Foam::Time when formatting the name.

SourceFiles
    instants.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_instant_H
#define Foam_instant_H

#include "db/Time/instant/Instant.H"
#include "primitives/strings/word/word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class instant Declaration
\*---------------------------------------------------------------------------*/

class instant
:
    public Instant<word>
{
public:

    // Static Data Members

        //- The type name (eg, for pTraits)
        static const char* const typeName;


    // Generated Methods

        //- Default construct, with value = 0 and empty name
        instant() = default;

        //- Copy construct
        instant(const instant&) = default;

        //- Move construct
        instant(instant&&) = default;

        //- Copy assignment
        instant& operator=(const instant&) = default;

        //- Move assignment
        instant& operator=(instant&&) = default;


    // Constructors

        //- Inherit all constructors (eg, copy/move construct from components)
        using Instant<word>::Instant;

        //- Construct from timeValue, auto generating the name
        explicit instant(scalar timeValue);

        //- Copy construct from timeName, parsing timeName for a value
        explicit instant(const word& timeName);

        //- Move construct from timeName, parsing timeName for a value
        explicit instant(word&& timeName);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
