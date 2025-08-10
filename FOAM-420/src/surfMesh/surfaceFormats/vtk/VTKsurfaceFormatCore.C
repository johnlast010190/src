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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "surfaceFormats/vtk/VTKsurfaceFormatCore.H"
#include "global/clock/clock.H"
#include "vtk/output/foamVtkOutput.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::VTKsurfaceFormatCore::writeHeader
(
    vtk::formatter& format,
    const pointField& pts
)
{
    vtk::legacy::fileHeader
    (
        format,
        ("surface written " + clock::dateTime()),
        vtk::fileTag::POLY_DATA
    );

    vtk::legacy::beginPoints(format.os(), pts.size());

    vtk::writeList(format, pts);
    format.flush();
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeCellData
(
    vtk::formatter& format,
    const UList<surfZone>& zones
)
{
    // Zone ids as CellData

    // Number of faces covered by the zones
    label nFaces = 0;
    for (const auto& z : zones)
    {
        nFaces += z.size();
    }

    vtk::legacy::dataHeader
    (
        format.os(),
        vtk::fileTag::CELL_DATA,
        nFaces,
        1  // Only one field
    );

    vtk::legacy::intField
    (
        format.os(),
        "region",
        1, // nComponent
        nFaces
    );

    label zoneId = 0;
    for (const surfZone& zone : zones)
    {
        forAll(zone, i)
        {
            format.write(zoneId);
        }
        ++zoneId;
    }
    format.flush();
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeCellData
(
    vtk::formatter& format,
    const labelUList& zoneIds
)
{
    // Zone ids as CellData

    // Number of faces
    const label nFaces = zoneIds.size();

    vtk::legacy::dataHeader
    (
        format.os(),
        vtk::fileTag::CELL_DATA,
        nFaces,
        1  // Only one field
    );

    vtk::legacy::intField
    (
        format.os(),
        "region",
        1, // nComponent
        nFaces
    );

    vtk::writeList(format, zoneIds);
    format.flush();
}


// ************************************************************************* //
