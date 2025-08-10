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

#include "surfaceFormats/vtp/VTPsurfaceFormat.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "vtk/output/foamVtkOutput.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// File-scope constant.
//
// TODO: make this run-time selectable
// - No append mode supported
// - Legacy mode is dispatched via 'VTKsurfaceFormat' instead

static const Foam::vtk::formatType fmtType =
    Foam::vtk::formatType::INLINE_ASCII;
    // Foam::vtk::formatType::INLINE_BASE64;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::writePolys
(
    vtk::formatter& format,
    const UList<Face>& faces
)
{
    format.tag(vtk::fileTag::POLYS);

    //
    // 'connectivity'
    //
    {
        uint64_t payLoad = 0;
        for (const auto& f : faces)
        {
            payLoad += f.size();
        }

        format.openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
            .closeTag();

        format.writeSize(payLoad * sizeof(label));

        for (const Face& f : faces)
        {
            vtk::writeList(format, f);
        }

        format.flush();
        format.endDataArray();
    }


    //
    // 'offsets'  (connectivity offsets)
    //
    {
        const uint64_t payLoad(faces.size() * sizeof(label));

        format
            .openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
            .closeTag();

        format.writeSize(payLoad);

        label off = 0;
        for (const auto& f : faces)
        {
            off += f.size();

            format.write(off);
        }

        format.flush();
        format.endDataArray();
    }

    format.endTag(vtk::fileTag::POLYS);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::VTPsurfaceFormat<Face>::VTPsurfaceFormat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurfaceProxy<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>&  faceLst = surf.surfFaces();
    const List<label>& faceMap = surf.faceMap();

    const List<surfZone>& zones =
    (
        surf.surfZones().empty()
      ? surfaceFormatsCore::oneZone(faceLst)
      : surf.surfZones()
    );

    const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

    std::ofstream os(filename.c_str(), std::ios::binary);

    autoPtr<vtk::formatter> format =
        vtk::newFormatter(os, fmtType);

    writeHeader(format(), pointLst, faceLst.size());

    if (useFaceMap)
    {
        format().tag(vtk::fileTag::POLYS);

        //
        // 'connectivity'
        //
        {
            uint64_t payLoad = 0;
            for (const auto& f : faceLst)
            {
                payLoad += f.size();
            }

            format().openDataArray<label>(vtk::dataArrayAttr::CONNECTIVITY)
                .closeTag();

            format().writeSize(payLoad * sizeof(label));

            label faceIndex = 0;
            for (const surfZone& zone : zones)
            {
                forAll(zone, i)
                {
                    const Face& f = faceLst[faceMap[faceIndex++]];

                    vtk::writeList(format(), f);
                }
            }

            format().flush();
            format().endDataArray();
        }


        //
        // 'offsets'  (connectivity offsets)
        //
        {
            const uint64_t payLoad(faceLst.size() * sizeof(label));

            format()
                .openDataArray<label>(vtk::dataArrayAttr::OFFSETS)
                    .closeTag();

            format().writeSize(payLoad);

            label off = 0, faceIndex = 0;
            for (const surfZone& zone : zones)
            {
                forAll(zone, i)
                {
                    const Face& f = faceLst[faceMap[faceIndex++]];

                    off += f.size();

                    format().write(off);
                }
            }

            format().flush();
            format().endDataArray();
        }

        format().endTag(vtk::fileTag::POLYS);
    }
    else
    {
        // Easy to write polys without a faceMap
        writePolys(format(), faceLst);
    }

    // Write regions (zones) as CellData
    if (zones.size() > 1)
    {
        writeCellData(format(), zones);
    }

    writeFooter(format());
}


template<class Face>
void Foam::fileFormats::VTPsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    std::ofstream os(filename.c_str(), std::ios::binary);

    autoPtr<vtk::formatter> format =
        vtk::newFormatter(os, fmtType);

    const List<Face>& faceLst = surf.surfFaces();

    writeHeader(format(), surf.points(), faceLst.size());

    // Easy to write polys without a faceMap
    writePolys(format(), faceLst);

    // Write regions (zones) as CellData
    writeCellData(format(), surf.zoneIds());

    writeFooter(format());
}


// ************************************************************************* //
