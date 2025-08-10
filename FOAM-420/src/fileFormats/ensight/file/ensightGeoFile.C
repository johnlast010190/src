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
    (c) 2011 OpenFOAM Foundation
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "ensight/file/ensightGeoFile.H"
#include "global/foamVersion.H"

// Macros to stringify macro contents.
#define STRINGIFY(content)      #content
#define STRING_QUOTE(input)     STRINGIFY(input)


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightGeoFile::initialize()
{
    writeBinaryHeader();

    // Description line 1
    write("Ensight Geometry File");
    newline();

    // Description line 2
    write(string("Written by FOAM-" + string(Foam::FOAMversion)));

    newline();
    write("node id assign");
    newline();

    write("element id assign");
    newline();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightGeoFile::ensightGeoFile
(
    const fileName& pathname,
    IOstream::streamFormat format
)
:
    ensightFile(pathname, format)
{
    initialize();
}


Foam::ensightGeoFile::ensightGeoFile
(
    const fileName& path,
    const fileName& name,
    IOstream::streamFormat format
)
:
    ensightFile(path, name, format)
{
    initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightGeoFile::~ensightGeoFile()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::ensightGeoFile::writeKeyword(const keyType& key)
{
    // ensure we get ensightFile::write(const string&)
    write(static_cast<const string&>(key));
    newline();

    return *this;
}


//
// Convenience Output Methods
//

void Foam::ensightGeoFile::beginPart
(
    const label index,
    const string& description
)
{
    beginPart(index);
    write(description);
    newline();
}


void Foam::ensightGeoFile::beginCoordinates(const label npoints)
{
    write("coordinates");
    newline();
    write(npoints);
    newline();
}


// ************************************************************************* //
