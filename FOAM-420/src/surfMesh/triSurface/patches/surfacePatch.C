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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "triSurface/patches/surfacePatch.H"
#include "surfZone/surfZone/surfZone.H"
#include "db/dictionary/dictionary.H"
#include "primitives/strings/word/word.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfacePatch, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfacePatch::surfacePatch()
:
    geometricSurfacePatch(word::null, word::null, -1),
    size_(0),
    start_(0)
{}


Foam::surfacePatch::surfacePatch(const label index)
:
    geometricSurfacePatch(word::null, word::null, index),
    size_(0),
    start_(0)
{}


Foam::surfacePatch::surfacePatch
(
    const word& geometricType,
    const word& name,
    const label size,
    const label start,
    const label index
)
:
    geometricSurfacePatch(geometricType, name, index),
    size_(size),
    start_(start)
{}


Foam::surfacePatch::surfacePatch(Istream& is, const label index)
:
    geometricSurfacePatch(is, index),
    size_(0),
    start_(0)
{
    size_ = readLabel(is);
    start_ = readLabel(is);
}


Foam::surfacePatch::surfacePatch
(
    const word& name,
    const dictionary& dict,
    const label index
)
:
    geometricSurfacePatch(name, dict, index),
    size_(readLabel(dict.lookup("nFaces"))),
    start_(readLabel(dict.lookup("startFace")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfacePatch::write(Ostream& os) const
{
    os  << nl
        << static_cast<const geometricSurfacePatch&>(*this) << endl
        << size() << tab << start();
}

void Foam::surfacePatch::writeDict(Ostream& os) const
{
    os  << nl << name() << nl << token::BEGIN_BLOCK << nl;

    geometricSurfacePatch::writeDict(os);

    os  << "    nFaces " << size() << ';' << nl
        << "    startFace " << start() << ';' << nl
        << token::END_BLOCK << endl;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

Foam::surfacePatch::operator Foam::surfZone() const
{
    return surfZone
    (
        this->name(),
        this->size(),
        this->start(),
        this->index(),
        this->geometricType()
    );
}


bool Foam::surfacePatch::operator!=(const surfacePatch& p) const
{
    return !(*this == p);
}


bool Foam::surfacePatch::operator==(const surfacePatch& p) const
{
    return
    (
        (geometricType() == p.geometricType())
     && (size() == p.size())
     && (start() == p.start())
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfacePatch& p)
{
    p.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
