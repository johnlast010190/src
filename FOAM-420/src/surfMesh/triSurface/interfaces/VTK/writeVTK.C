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
    (c) 2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "vtk/output/foamVtkOutput.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// File-scope constant.
//
// TODO: make this run-time selectable (ASCII | BINARY)
// - Legacy mode only

static const Foam::vtk::formatType fmtType =
    Foam::vtk::formatType::LEGACY_ASCII;
    // Foam::vtk::formatType::LEGACY_BINARY;


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triSurface::writeVTK
(
    const bool writeSorted,
    std::ostream& os
) const
{
    autoPtr<vtk::formatter> format =
        vtk::newFormatter(os, fmtType);

    // Header
    vtk::legacy::fileHeader
    (
        format(),
        "triSurface",
        vtk::fileTag::POLY_DATA
    );

    const pointField& pts = this->points();
    const auto& faceLst = this->surfFaces();

    const label nFaces = faceLst.size();

    // Points
    vtk::legacy::beginPoints(os, pts.size());

    vtk::writeList(format(), pts);
    format().flush();

    // connectivity count (without additional storage)
    // is simply 3 * nFaces

    vtk::legacy::beginPolys(os, nFaces, 3*nFaces);

    labelList faceMap;
    surfacePatchList patches(calcPatches(faceMap));

    const bool useFaceMap = (writeSorted && patches.size() > 1);

    if (useFaceMap)
    {
        label faceIndex = 0;
        for (const surfacePatch& patch : patches)
        {
            forAll(patch, i)
            {
                const Face& f = faceLst[faceMap[faceIndex++]];

                format().write(label(3));   // The size prefix
                vtk::writeList(format(), f);
            }
        }
        format().flush();


        // Write regions (zones) as CellData
        if (patches.size() > 1)
        {
            vtk::legacy::dataHeader
            (
                os,
                vtk::fileTag::CELL_DATA,
                nFaces,
                1  // Only one field
            );

            vtk::legacy::intField
            (
                os,
                "region",
                1, // nComponent
                nFaces
            );

            faceIndex = 0;
            for (const surfacePatch& patch : patches)
            {
                forAll(patch, i)
                {
                    const Face& f = faceLst[faceMap[faceIndex++]];
                    format().write(f.region());
                }
            }
            format().flush();
        }
    }
    else
    {
        // No faceMap (unsorted)

        for (const Face& f : faceLst)
        {
            format().write(label(3));   // The size prefix
            vtk::writeList(format(), f);
        }
        format().flush();


        // Write regions (zones) as CellData

        vtk::legacy::dataHeader
        (
            os,
            vtk::fileTag::CELL_DATA,
            faceLst.size(),
            1  // Only one field
        );

        vtk::legacy::intField
        (
            os,
            "region",
            1, // nComponent
            faceLst.size()
        );

        for (const Face& f : faceLst)
        {
            format().write(f.region());
        }
        format().flush();
    }
}


// ************************************************************************* //
