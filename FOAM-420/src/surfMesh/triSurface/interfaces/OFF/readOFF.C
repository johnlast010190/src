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
    (c) 2022 Esi Ltd.

Description
    Geomview OFF polyList format. Does triangulation.

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "db/IOstreams/StringStreams/IStringStream.H"
#include "algorithms/polygonTriangulate/polygonTriangulate.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::triSurface::readOFF(const fileName& OFFfileName)
{
    IFstream OFFfile(OFFfileName);

    if (!OFFfile.good())
    {
        FatalErrorInFunction
            << "Cannot read file " << OFFfileName
            << exit(FatalError);
    }

    // Read header
    string hdr = getLineNoComment(OFFfile);
    if (hdr != "OFF")
    {
        FatalErrorInFunction
            << "OFF file " << OFFfileName
            << " does not start with 'OFF'"
            << exit(FatalError);
    }


    label nPoints, nEdges, nElems;

    string line = getLineNoComment(OFFfile);
    {
        IStringStream lineStream(line);
        lineStream >> nPoints >> nElems >> nEdges;
    }

    // Read points
    pointField points(nPoints);

    forAll(points, pointi)
    {
        scalar x, y, z;
        line = getLineNoComment(OFFfile);
        {
            IStringStream lineStream(line);
            lineStream >> x >> y >> z;
        }
        points[pointi] = point(x, y, z);
    }

    // Create a triangulation engine
    polygonTriangulate triEngine;

    // Read faces & triangulate them
    DynamicList<labelledTri> tris(nElems);

    for (label facei = 0; facei < nElems; facei++)
    {
        line = getLineNoComment(OFFfile);
        {
            IStringStream lineStream(line);

            label nVerts;
            lineStream >> nVerts;

            face f(nVerts);

            forAll(f, fp)
            {
                lineStream >> f[fp];
            }

            // Triangulate.
            if (nVerts == 3)
            {
                tris.append(labelledTri(f[0], f[1], f[2], 0));
            }
            else if (nVerts == 4)
            {
                tris.append(labelledTri(f[0], f[1], f[2], 0));
                tris.append(labelledTri(f[2], f[3], f[0], 0));
            }
            else
            {
                triEngine.triangulate(UIndirectList<point>(points, f));

                forAll(triEngine.triPoints(), triI)
                {
                    tris.append(labelledTri(triEngine.triPoints(triI, f), 0));
                }
            }
        }
    }

    tris.shrink();

    *this = triSurface(tris, points);

    return true;
}


// ************************************************************************* //
