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
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/IOobjects/IOdictionary/baseIOdictionary.H"
#include "db/objectRegistry/objectRegistry.H"
#include "db/IOstreams/Pstreams/Pstream.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(baseIOdictionary, 0);

bool baseIOdictionary::writeDictionaries
(
    debug::infoSwitch("writeDictionaries", 0)
);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::baseIOdictionary::baseIOdictionary(const IOobject& io)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();
}


Foam::baseIOdictionary::baseIOdictionary
(
    const IOobject& io,
    const dictionary& dict
)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();
}


Foam::baseIOdictionary::baseIOdictionary
(
    const IOobject& io,
    Istream& is
)
:
    regIOobject(io)
{
    dictionary::name() = IOobject::objectPath();
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::baseIOdictionary::~baseIOdictionary()
{}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

const Foam::word& Foam::baseIOdictionary::name() const
{
    return regIOobject::name();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::baseIOdictionary::operator=(const baseIOdictionary& rhs)
{
    dictionary::operator=(rhs);
}


// ************************************************************************* //
