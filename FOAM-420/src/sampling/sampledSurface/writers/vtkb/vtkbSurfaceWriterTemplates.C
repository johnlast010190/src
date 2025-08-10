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
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "db/IOstreams/Fstreams/OFstream.H"
#include "include/OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::fileName Foam::vtkbSurfaceWriter::writeTemplate
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const word& fieldName,
    const Field<Type>& values,
    const bool isNodeValues,
    const bool verbose
) const
{
    if (!isDir(outputDir))
    {
        mkDir(outputDir);
    }

    fileName fName(fieldName + '_' + surfaceName + ".vtk");

    std::ofstream os(outputDir/fName);
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

    writeGeometry(os, format, surf);

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << fName << endl;
    }

    // start writing data
    if (isNodeValues)
    {
        vtk::legacy::dataHeader(os, vtk::fileTag::POINT_DATA, values.size(), 1);
    }
    else
    {
        vtk::legacy::dataHeader(os, vtk::fileTag::CELL_DATA, values.size(), 1);
    }

    const int nCmpt(pTraits<Type>::nComponents);
    vtk::legacy::floatField(os, fieldName, nCmpt, values.size());
    const uint64_t payLoad(values.size() * nCmpt * sizeof(float));
    format().writeSize(payLoad);
    vtk::writeList(format(), values);
    format().flush();

    return fName;
}


// ************************************************************************* //
