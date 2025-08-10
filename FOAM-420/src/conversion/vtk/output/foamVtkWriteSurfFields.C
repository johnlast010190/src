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

#include "vtk/output/foamVtkWriteSurfFields.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "fields/fvsPatchFields/constraint/empty/emptyFvsPatchFields.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "vtk/output/foamVtkOutput.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void Foam::vtk::writeSurfFields
(
    const fvMesh& mesh,
    const fileName& baseName,
    const vtk::outputOptions outOpts,
    const UPtrList<const surfaceVectorField>& surfVectorFields
)
{
    outputOptions opts(outOpts);
    opts.append(false);  // No append supported

    const bool legacy_(opts.legacy());

    std::ofstream os(baseName + (legacy_ ? ".vtk" : ".vtp"));
    autoPtr<vtk::formatter> format = opts.newFormatter(os);

    // Same payload size for points and vector fields!
    const int nCmpt(3);  // vector
    const uint64_t payLoad(mesh.nFaces() * 3 * sizeof(float));

    if (legacy_)
    {
        legacy::fileHeader(format(), "surfaceFields", vtk::fileTag::POLY_DATA);
        legacy::beginPoints(os, mesh.nFaces());
    }
    else
    {
        // XML (inline)

        format()
            .xmlHeader()
            .xmlComment("surfaceFields")
            .beginVTKFile(vtk::fileTag::POLY_DATA, "0.1");

        // Tricky - hide in beginPiece()
        format()
            .openTag(vtk::fileTag::PIECE)
            .xmlAttr(vtk::fileAttr::NUMBER_OF_POINTS, mesh.nFaces())
            .closeTag();

        format().tag(vtk::fileTag::POINTS)
            .openDataArray<float,3>(vtk::dataArrayAttr::POINTS)
            .closeTag();
    }

    const pointField& fc = mesh.faceCentres();

    format().writeSize(payLoad);
    vtk::writeList(format(), fc);
    format().flush();

    if (!legacy_)
    {
        format()
            .endDataArray()
            .endTag(vtk::fileTag::POINTS);
    }


    // Fields
    if (legacy_)
    {
        legacy::dataHeader
        (
            os,
            vtk::fileTag::POINT_DATA,
            mesh.nFaces(),
            surfVectorFields.size()
        );
    }
    else
    {
        format().tag(vtk::fileTag::POINT_DATA);
    }

    // surfVectorFields
    forAll(surfVectorFields, fieldi)
    {
        const auto& fld = surfVectorFields[fieldi];

        if (legacy_)
        {
            legacy::floatField(os, fld.name(), nCmpt, mesh.nFaces());
        }
        else
        {
            format().openDataArray<float, nCmpt>(fld.name())
                .closeTag();
        }

        format().writeSize(payLoad);

        for (label facei=0; facei < mesh.nInternalFaces(); ++facei)
        {
            vtk::write(format(), fld[facei]);
        }

        forAll(fld.boundaryField(), patchi)
        {
            const fvPatch& pp = mesh.boundary()[patchi];
            const auto& pf = fld.boundaryField()[patchi];

            if (isA<emptyFvsPatchVectorField>(pf))
            {
                // Note: loop over polypatch size, not fvpatch size.
                forAll(pp.patch(), i)
                {
                    vtk::write(format(), vector::zero);
                }
            }
            else
            {
                vtk::writeList(format(), pf);
            }
        }

        format().flush();

        if (!legacy_)
        {
            format().endDataArray();
        }
    }

    if (!legacy_)
    {
        format().endTag(vtk::fileTag::POINT_DATA);

        // slight cheat. </Piece> too
        format().endTag(vtk::fileTag::PIECE);

        format().endTag(vtk::fileTag::POLY_DATA)
            .endVTKFile();
    }
}


// ************************************************************************* //
