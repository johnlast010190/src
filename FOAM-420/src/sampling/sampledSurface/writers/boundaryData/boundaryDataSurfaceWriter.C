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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "sampledSurface/writers/boundaryData/boundaryDataSurfaceWriter.H"
#include "sampledSurface/writers/makeSurfaceWriterMethods.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    makeSurfaceWriterType(boundaryDataSurfaceWriter);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryDataSurfaceWriter::boundaryDataSurfaceWriter()
:
    surfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boundaryDataSurfaceWriter::~boundaryDataSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::boundaryDataSurfaceWriter::write
(
    const fileName& outputDir,
    const fileName& surfaceName,
    const meshedSurf& surf,
    const bool verbose
) const
{
    const fileName baseDir(outputDir.path()/surfaceName);
    const fileName timeName(outputDir.name());

    const pointField& points = surf.points();

    // Construct dummy time to use as an objectRegistry
    const fileName caseDir(getEnv("FOAM_CASE"));
    Time dummyTime
    (
        caseDir.path(), //rootPath,
        caseDir.name(), //caseName,
        "system",       //systemName,
        "constant",     //constantName,
        false           //enableFunctionObjects
    );


    // Write points
    if (verbose)
    {
        Info<< "Writing points to " << baseDir/"points" << endl;
    }

    pointIOField pts
    (
        IOobject
        (
            baseDir/"points",
            dummyTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        points
    );

    {
        // Do like regIOobject::writeObject but don't do instance() adaptation
        // since this would write to e.g. 0/ instead of postProcessing/

        // Try opening an OFstream for object
        mkDir(pts.path());
        OFstream os(pts.objectPath());

        pts.writeData(os);
    }

    return baseDir;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// create write methods
defineSurfaceWriterWriteFields(Foam::boundaryDataSurfaceWriter);


// ************************************************************************* //
