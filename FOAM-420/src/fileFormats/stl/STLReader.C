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
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "stl/STLReader.H"
#include "containers/HashTables/Map/Map.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "meshes/meshTools/mergePoints.H"
#include "db/IOstreams/gzstream/gzstream.h"
#include "include/OSspecific.H"

#undef DEBUG_STLBINARY

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::fileFormats::STLReader::readBlockBINARY
(
    const fileName& filename
)
{
    sorted_ = true;
    format_ = UNKNOWN;

    autoPtr<istream> streamPtr
    (
        new ifstream(filename.c_str(), std::ios::binary)
    );

    // If the file is compressed, decompress it before reading.
    if (!streamPtr->good() && isFile(filename + ".gz", false))
    {
        streamPtr.reset(new igzstream((filename + ".gz").c_str()));
    }
    istream& is = streamPtr();

    DynamicList<label> dynSizes;

    label nTris = 0;
    label facei = 0;
    label ptI = 0;

    label nBlocks = 0;

    const unsigned STLHeaderSize = 80;
    // Read the STL header
    char header[STLHeaderSize];
    is.read(header, STLHeaderSize);

    const unsigned patchHeaderSize = 160;

    while (true)
    {
        // Read the STL header
        char patchHeader[patchHeaderSize];
        is.read(patchHeader, patchHeaderSize);

        if (is.eof())
        {
            break;
        }
        word blockName(patchHeader);

        // Read the number of triangles in the STl file
        // (note: read as int so we can check whether >2^31)
        int32_t nBlockTris;
        is.read(reinterpret_cast<char*>(&nBlockTris), sizeof(int32_t));

        nTris += nBlockTris;

        points_.setSize(3*nTris);
        zoneIds_.setSize(nTris);
        names_.setSize(nBlocks+1);
        names_[nBlocks] = blockName;

        dynSizes.append(0);

        for (label i = 0; i < nBlockTris; i++)
        {
            // Read STL triangle
            STLtriangle stlTri(is);

            // transcribe the vertices of the STL triangle -> points
            points_[ptI++] = stlTri.a();
            points_[ptI++] = stlTri.b();
            points_[ptI++] = stlTri.c();

            zoneIds_[facei++] = nBlocks;
            dynSizes[nBlocks]++;
        }
        nBlocks++;
    }

    sizes_.transfer(dynSizes);
    format_ = BLOCKBINARY;

    return true;
}


bool Foam::fileFormats::STLReader::readBINARY
(
    const fileName& filename
)
{
    sorted_ = true;
    format_ = UNKNOWN;

    label nTris = 0;
    autoPtr<istream> streamPtr = readBinaryHeader(filename, nTris);

    if (!streamPtr.valid())
    {
        FatalErrorInFunction
            << "Error reading file " << filename
            << " or file " << filename + ".gz"
            << exit(FatalError);
    }

    istream& is = streamPtr();

#ifdef DEBUG_STLBINARY
    Info<< "# " << nTris << " facets" << endl;
    label prevZone = -1;
#endif

    points_.setSize(3*nTris);
    zoneIds_.setSize(nTris);

    Map<label> lookup;
    DynamicList<label> dynSizes;

    label ptI = 0;
    label zoneI = -1;
    forAll(zoneIds_, facei)
    {
        // Read STL triangle
        STLtriangle stlTri(is);

        // transcribe the vertices of the STL triangle -> points
        points_[ptI++] = stlTri.a();
        points_[ptI++] = stlTri.b();
        points_[ptI++] = stlTri.c();

        // interpret STL attribute as a zone
        const label origId = stlTri.attrib();

        Map<label>::const_iterator fnd = lookup.find(origId);
        if (fnd != lookup.end())
        {
            if (zoneI != fnd())
            {
                // group appeared out of order
                sorted_ = false;
            }
            zoneI = fnd();
        }
        else
        {
            zoneI = dynSizes.size();
            lookup.insert(origId, zoneI);
            dynSizes.append(0);
        }

        zoneIds_[facei] = zoneI;
        dynSizes[zoneI]++;

#ifdef DEBUG_STLBINARY
        if (prevZone != zoneI)
        {
            if (prevZone != -1)
            {
                Info<< "endsolid zone" << prevZone << nl;
            }
            prevZone = zoneI;

            Info<< "solid zone" << prevZone << nl;
        }

        stlTri.print(Info);
#endif
    }

#ifdef DEBUG_STLBINARY
    if (prevZone != -1)
    {
        Info<< "endsolid zone" << prevZone << nl;
    }
#endif

    names_.clear();
    sizes_.transfer(dynSizes);

    format_ = BINARY;
    return true;
}


bool Foam::fileFormats::STLReader::readFile
(
    const fileName& filename,
    const STLFormat& format
)
{
    bool blockBinary(false);
    if
    (
        format == UNKNOWN ?
        detectBinaryHeader(filename,blockBinary) : format == BINARY
    )
    {
        if (blockBinary)
        {
            return readBlockBINARY(filename);
        }
        else
        {
            return readBINARY(filename);
        }
    }
    else
    {
        return readASCII(filename);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::STLReader::STLReader
(
    const fileName& filename
)
:
    sorted_(true),
    points_(),
    zoneIds_(),
    names_(),
    sizes_(),
    format_(STLCore::UNKNOWN)
{
    // Auto-detect ASCII/BINARY format
    readFile(filename, STLCore::UNKNOWN);
}


Foam::fileFormats::STLReader::STLReader
(
    const fileName& filename,
    const STLFormat& format
)
:
    sorted_(true),
    points_(),
    zoneIds_(),
    names_(),
    sizes_(),
    format_(STLCore::UNKNOWN)
{
    // Manually specified ASCII/BINARY format
    readFile(filename, format);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileFormats::STLReader::~STLReader()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fileFormats::STLReader::clear()
{
    sorted_ = true;
    points_.clear();
    zoneIds_.clear();
    names_.clear();
    sizes_.clear();
    format_ = UNKNOWN;
}


Foam::label Foam::fileFormats::STLReader::mergePointsMap
(
    labelList& pointMap
) const
{
    // With the merge distance depending on the input format (ASCII | BINARY),
    // but must be independent of FOAM_SP or FOAM_DP flag.
    // - floatScalarSMALL  = 1e-6
    // - doubleScalarSMALL = 1e-15

    return mergePointsMap
    (
        (format_ == BINARY ? 10 : 100) * doubleScalarSMALL,
        pointMap
    );
}


Foam::label Foam::fileFormats::STLReader::mergePointsMap
(
    const scalar mergeTol,
    labelList& pointMap
) const
{
    return Foam::mergePoints
    (
        points_,
        mergeTol,
        false, // verbose
        pointMap
    );
}


// ************************************************************************* //
