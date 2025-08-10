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

#include "stl/STLReader.H"
#include "triSurface/triSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readSTL(const fileName& STLfileName, bool forceBinary)
{
    // Read in the values
    fileFormats::STLReader reader
    (
        STLfileName,
        (
            forceBinary
          ? fileFormats::STLCore::BINARY
          : fileFormats::STLCore::UNKNOWN
        )
    );

    // Get the map for stitched surface points, with merge tolerance depending
    // on the input format
    labelList pointMap;
    const label nUniquePoints = reader.mergePointsMap(pointMap);

    const auto& readpts = reader.points();
    const labelList& zoneIds = reader.zoneIds();

    pointField& pointLst = storedPoints();
    List<Face>& faceLst  = storedFaces();

    // Sizing
    pointLst.setSize(nUniquePoints);
    faceLst.setSize(zoneIds.size());

    // Assign points
    forAll(readpts, pointi)
    {
        pointLst[pointMap[pointi]] = readpts[pointi];
    }

    // Assign triangles
    label pointi = 0;
    forAll(faceLst, i)
    {
        Face& f = faceLst[i];

        f[0] = pointMap[pointi++];
        f[1] = pointMap[pointi++];
        f[2] = pointMap[pointi++];
        f.region() = zoneIds[i];
    }

    // Set patch name/index.
    if
    (
        reader.stlFormat() == fileFormats::STLCore::ASCII
        || reader.stlFormat() == fileFormats::STLCore::BLOCKBINARY
    )
    {
        const List<word>& names = reader.names();

        patches_.setSize(names.size());
        forAll(patches_, patchi)
        {
            patches_[patchi] = geometricSurfacePatch(names[patchi], patchi);
        }
    }

    checkDegenerate(false);

    return true;
}

// ************************************************************************* //
