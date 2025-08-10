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
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "sampledSurface/writers/vtkb/vtkbSurfaceWriter.H"
#include "sampledSurface/writers/makeSurfaceWriterMethods.H"
#include "vtk/output/foamVtkOutput.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(vtkbSurfaceWriter);
    addToRunTimeSelectionTable(surfaceWriter, vtkbSurfaceWriter, wordDict);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtkbSurfaceWriter::writeGeometry
(
    std::ofstream& os,
    autoPtr<vtk::formatter>& format,
    const meshedSurf& surf
)
{
    const pointField& points = surf.points();
    const faceList& faces = surf.faces();

    vtk::legacy::fileHeader(format(), "sampleSurface", vtk::fileTag::POLY_DATA);

    // Write points and faces as polygons
    vtk::legacy::beginPoints(os, points.size());

    vtk::writeList(format(), points);
    format().flush();

    label nNodes = 0;
    forAll(faces, facei)
    {
        nNodes += faces[facei].size();
    }
    vtk::legacy::beginPolys(os, faces.size(), nNodes);
    forAll(faces, facei)
    {
        face f = faces[facei];
        format().write(f.size());  // The size prefix
        vtk::writeList(format(), f);
    }
    format().flush();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtkbSurfaceWriter::vtkbSurfaceWriter()
:
    surfaceWriter()
{}


Foam::vtkbSurfaceWriter::vtkbSurfaceWriter(const dictionary& dict)
:
    surfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtkbSurfaceWriter::~vtkbSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::vtkbSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    fileName fName(surfaceName + ".vtk");
    std::ofstream os(outputDir/fName);

    //Set for legacy/binary writing
    vtk::outputOptions opts;
    if (sizeof(floatScalar) != 4)
    {
        opts.ascii(true);
    }
    else
    {
        opts.ascii(false);
    }
    opts.legacy(true);

    autoPtr<vtk::formatter> format = opts.newFormatter(os);

    if (verbose)
    {
        Info<< "Writing geometry to " << fName << endl;
    }

    writeGeometry(os, format, surf);

    return fName;
}


// create write methods
defineSurfaceWriterWriteFields(Foam::vtkbSurfaceWriter);


// ************************************************************************* //
