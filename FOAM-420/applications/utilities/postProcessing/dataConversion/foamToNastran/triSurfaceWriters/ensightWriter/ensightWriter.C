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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2007 OpenCFD Ltd.

Description

\*---------------------------------------------------------------------------*/

#include "ensightWriter.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "include/OSspecific.H"
#include "meshes/GeoMesh/GeoMesh.H"
#include "fields/DimensionedFields/DimensionedField/DimensionedField.H"
#include "db/IOstreams/Fstreams/OFstream.H"


namespace Foam
{
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(ensightWriter, 0);
addToRunTimeSelectionTable(triSurfaceWriter, ensightWriter, word);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ensightWriter::swapWord(label& word32)
{
    char* mem = reinterpret_cast<char*>(&word32);

    char a = mem[0];
    mem[0] = mem[3];
    mem[3] = a;

    a = mem[1];
    mem[1] = mem[2];
    mem[2] = a;
}


void Foam::ensightWriter::swapWords(const label nWords, label* words32)
{
    for (label i = 0; i < nWords; i++)
    {
        swapWord(words32[i]);
    }
}


void Foam::ensightWriter::write
(
    std::ostream& os,
    const bool binary,
    List<floatScalar>& fField
)
{
    if (binary)
    {
#       ifdef LITTLEENDIAN
        swapWords(fField.size(), reinterpret_cast<label*>(fField.begin()));
#       endif
        os.write
        (
            reinterpret_cast<char*>(fField.begin()),
            fField.size()*sizeof(float)
        );

        os << std::endl;
    }
    else
    {
        forAll(fField, i)
        {
            os << fField[i] << ' ';

            if (i > 0 && (i % 10) == 0)
            {
                os << std::endl;
            }
        }
        os << std::endl;
    }
}


void Foam::ensightWriter::write
(
    std::ostream& os,
    const bool binary,
    labelList& elems
)
{
    if (binary)
    {
#       ifdef LITTLEENDIAN
        swapWords(elems.size(), reinterpret_cast<label*>(elems.begin()));
#       endif
        os.write
        (
            reinterpret_cast<char*>(elems.begin()),
            elems.size()*sizeof(label)
        );

        os << std::endl;
    }
    else
    {
        forAll(elems, i)
        {
            os << elems[i] << ' ';

            if (i > 0 && (i % 10) == 0)
            {
                os << std::endl;
            }
        }
        os << std::endl;
    }
}


// Store vector in dest.
void Foam::ensightWriter::insert(const point& pt, DynamicList<floatScalar>& dest)
{
    dest.append(float(pt.x()));
    dest.append(float(pt.y()));
    dest.append(float(pt.z()));
}


// Store labelList in dest.
void Foam::ensightWriter::insert(const labelList& source, DynamicList<label>& dest)
{
    forAll(source, i)
    {
        dest.append(source[i]);
    }
}


// Store scalarField in dest
void Foam::ensightWriter::insert
(
    const scalarField& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(source, i)
    {
       dest.append(float(source[i]));
    }
}


// Store pointField in dest
void Foam::ensightWriter::insert
(
    const List<point>& source,
    DynamicList<floatScalar>& dest
)
{
    forAll(source, i)
    {
       insert(source[i], dest);
    }
}


void Foam::ensightWriter::writeGeometry
(
    const bool binary,
    const triSurface& s,
    std::ostream& os
)
{
    // Write vertex coordinates
    const pointField& points = s.localPoints();
    const List<labelledTri>& faces = s.localFaces();

    os  << "POINTS " << points.size() << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*points.size());

    insert(points, ptField);

    write(os, binary, ptField.shrink());


    os << "POLYGONS " << faces.size() << ' ' << 4*faces.size()
        << std::endl;

    DynamicList<label> vertLabels(4*faces.size());

    forAll(faces, faceI)
    {
        const labelledTri& f = faces[faceI];

        vertLabels.append(f.size());

        insert(f, vertLabels);
    }
    write(os, binary, vertLabels.shrink());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::ensightWriter::ensightWriter()
:
    triSurfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightWriter::~ensightWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightWriter::write
(
    const fileName& path,           // <root>/<case>/sampleSurfaces
    const fileName& timeDir,        // time name
    const fileName& surfaceName,    // name of surface
    const triSurface& surface,
    PtrList<triSurfaceScalarField>& scalarFields,
    PtrList<triSurfaceVectorField>& vectorFields,
    PtrList<triSurfaceTensorField>& tensorFields
) const
{
    // Most of this is copied from foamToEnsight so if it doesn't make any
    // sense ...

    const word postProcDir = "EnSight";
    const word prepend = args.globalCaseName() + '.';

    fileName postProcPath = args.rootPath()/args.globalCaseName()/postProcDir;

    if (dir(postProcPath))
    {
        rmDir(postProcPath);
    }

    mkDir(postProcPath);

    // Open the Case file
    fileName ensightWriterCaseFileName = prepend + "case";

    OFstream ensightWriterCaseFile
    (
        postProcPath/ensightWriterCaseFileName,
        runTime.writeFormat(),
        runTime.writeVersion(),
        runTime.writeCompression()
    );

    Info<< nl << "Case file is " << ensightWriterCaseFileName << endl;

//XXX
    bool lagrangianExist = false;
    word sprayName;

#   include "ensightWriterCaseHeader.H"

    label nTimeSteps = 0;
    Times.setSize(1);
    label timeIndex = 0;

    word timeName = itoa(timeIndex);
    word timeFile = prepend + timeName;

    if (timeIndex == 0)
    {
        ensightWriterCaseFile << nl << "VARIABLE" << nl;
    }

}


// ************************************************************************* //
