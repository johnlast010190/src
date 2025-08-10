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
    (c) 2021-2021 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "primitives/Vector/floatVector/floatVector.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/StringStreams/IStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The number of bytes in the STL binary header
static const unsigned BTSHeaderSize = 80;

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//Read of Esi Binary Triangulated Surface (BTS) format
bool Foam::triSurface::readBTS(const fileName& BTSfileName)
{
    autoPtr<istream> streamPtr
    (
        new ifstream(BTSfileName.c_str(), std::ios::binary)
    );

    istream& is = streamPtr();

    if (!is.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << BTSfileName
             << exit(FatalError);
    }

    // Read the BTS header
    char header[BTSHeaderSize];
    is.read(header, BTSHeaderSize);
    unsigned pos = 0;
    while (std::isspace(header[pos]) && pos < BTSHeaderSize)
    {
        ++pos;
    }

    bool validHeader = false;
    if
    (
        pos < (BTSHeaderSize-12)  // At least 12 chars remaining
        && std::toupper(header[pos+0]) == 'E'
        && std::toupper(header[pos+1]) == 'N'
        && std::toupper(header[pos+2]) == 'G'
        && std::toupper(header[pos+3]) == 'Y'
        && std::toupper(header[pos+4]) == 'S'
        && std::toupper(header[pos+6]) == 'B'
        && std::toupper(header[pos+7]) == 'I'
        && std::toupper(header[pos+8]) == 'N'
        && std::toupper(header[pos+9]) == 'A'
        && std::toupper(header[pos+10]) == 'R'
        && std::toupper(header[pos+11]) == 'Y'
    )
    {
        validHeader = true;
    }

    if (!validHeader)
    {
        FatalErrorInFunction
            << "File " << BTSfileName
            << " has invalid header : " << header
            << exit(FatalError);
    }

    int32_t nRegions;
    is.read(reinterpret_cast<char*>(&nRegions), sizeof(int32_t));
    geometricSurfacePatchList patches(nRegions);

    const unsigned regionHeaderSize = 160;
    char regionHeader[regionHeaderSize];
    for (label regioni = 0; regioni < nRegions; regioni++)
    {
        is.read(regionHeader, regionHeaderSize);
        word regionName(regionHeader);
        patches[regioni] = geometricSurfacePatch
        (
            regionName,
            regioni
        );
    }

    pointField points;
    List<labelledTri> faces;

    label regionFaceStart = 0;
    label regionPointStart = 0;
    for (label regioni = 0; regioni < nRegions; regioni++)
    {
        uint32_t sz;
        is.read(reinterpret_cast<char*>(&sz), sizeof(uint32_t));
        label nPts(sz);
        points.setSize(nPts+regionPointStart);
        for (label lp = 0; lp < nPts; lp++)
        {
            label pointi = lp + regionPointStart;
            floatVector pt;
            is.read(reinterpret_cast<char*>(&pt), sizeof(floatVector));
            points[pointi] = pt;
        }

        is.read(reinterpret_cast<char*>(&sz), sizeof(uint32_t));
        label nFaces(sz);
        faces.setSize(nFaces+regionFaceStart);
        for (label lf = 0; lf < nFaces; lf++)
        {
            label facei = lf + regionFaceStart;
            uint32_t a;
            is.read(reinterpret_cast<char*>(&a), sizeof(uint32_t));
            uint32_t b;
            is.read(reinterpret_cast<char*>(&b), sizeof(uint32_t));
            uint32_t c;
            is.read(reinterpret_cast<char*>(&c), sizeof(uint32_t));

            a += regionPointStart;
            b += regionPointStart;
            c += regionPointStart;
            faces[facei] = labelledTri(a, b, c, regioni);
        }
        regionPointStart += nPts;
        regionFaceStart += nFaces;
    }

    // Create triSurface
    *this = triSurface(faces, patches, points, true);

    stitchTriangles(SMALL, false);
    clearOut();
    checkDegenerate(false);

    return true;
}


// ************************************************************************* //
