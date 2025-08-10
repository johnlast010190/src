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
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "surfZone/surfZone/surfZone.H"
#include "db/dictionary/dictionary.H"
#include "primitives/strings/word/word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(surfZone, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZone::surfZone()
:
    surfZoneIdentifier(),
    size_(0),
    start_(0)
{}


Foam::surfZone::surfZone
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& geometricType
)
:
    surfZoneIdentifier(name, index, geometricType),
    size_(size),
    start_(start)
{}


Foam::surfZone::surfZone(Istream& is, const label index)
:
    surfZoneIdentifier(),
    size_(0),
    start_(0)
{
    word name(is);
    dictionary dict(is);

    operator=(surfZone(name, dict, index));
}


Foam::surfZone::surfZone
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    surfZoneIdentifier(name, dict, index),
    size_(readLabel(dict.lookup("nFaces"))),
    start_(readLabel(dict.lookup("startFace")))
{}


Foam::surfZone::surfZone(const surfZone& zone, const label index)
:
    surfZoneIdentifier(zone, index),
    size_(zone.size()),
    start_(zone.start())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfZone::write(Ostream& os) const
{
    writeDict(os);
}


void Foam::surfZone::writeDict(Ostream& os) const
{
    os.beginBlock(name());
    surfZoneIdentifier::write(os);
    os.writeEntry("nFaces", size());
    os.writeEntry("startFace", start());
    os.endBlock();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::surfZone::operator!=(const surfZone& rhs) const
{
    return !(*this == rhs);
}


bool Foam::surfZone::operator==(const surfZone& rhs) const
{
    return
    (
        size()  == rhs.size()
     && start() == rhs.start()
     && geometricType() == rhs.geometricType()
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, surfZone& zone)
{
    zone = surfZone(is, 0);

    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const surfZone& zone)
{
    zone.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
