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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

Application
    surfaceInertia

Description
    Calculates the eMesh file of intersction of surface with a plane

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "triSurface/triSurface.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshTools/meshTools.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "meshes/primitiveShapes/line/linePointRef.H"
#include "edgeMesh/featureEdgeMesh/featureEdgeMesh.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "searchableSurfaces/searchableSurface/searchableSurface.H"
#include "meshes/primitiveShapes/plane/plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Calculates the intersection of a given surface with a plane"
    );

    argList::noParallel();

#include "include/addDictOption.H"

#include "include/setRootCase.H"
#include "include/createTime.H"

    const word dictName("surfaceIntersectionDict");
#include "include/setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;

    const IOdictionary dict(dictIO);

    forAllConstIter(dictionary, dict, iter)
    {
        const dictionary& surfaceDict = iter().dict();
        const fileName surfFileName(word(surfaceDict.lookup("fileName")));

        const fileName sFeatFileName = surfFileName.lessExt().name();

        Info<< "Surface            : " << surfFileName << nl << endl;

        const word eName =
            surfaceDict.lookupOrDefault<word>("eMeshName", sFeatFileName);

        const Switch writeVTK =
            surfaceDict.lookupOrDefault<Switch>("writeVTK", "off");

        const word intersectionMethod = surfaceDict.lookup("intersectionMethod");

        triSurface surf(runTime.constantPath()/"triSurface"/surfFileName);
        const edgeList& edges = surf.edges();
        const pointField& points = surf.localPoints();
        const labelListList& faceEdges = surf.faceEdges();

        DynamicList<point> cutPoints(surf.localPoints().size()/100);
        DynamicList<edge> cutEdges(surf.edges().size()/100);

        if (intersectionMethod == "plane")
        {
            const dictionary planeDict = surfaceDict.subDict("planeCoeffs");
            plane cutPlane(planeDict);

            const scalarField dotProducts
                ((points - cutPlane.refPoint()) & cutPlane.normal());

            labelList edgePoint(edges.size(), -1);

            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];

                if
                (
                    (dotProducts[e[0]] <= 0.0 && dotProducts[e[1]] > 0.0)
                    || (dotProducts[e[1]] <= 0.0 && dotProducts[e[0]] > 0.0)
                 )
                {
                    edgePoint[edgeI] = cutPoints.size();

                    const point& p0 = points[e[0]];
                    const point& p1 = points[e[1]];

                    scalar alpha = cutPlane.lineIntersect(linePointRef(p0, p1));

                    if (alpha < SMALL)
                    {
                        cutPoints.append(p0);
                    }
                    else if (alpha >= (1.0 - SMALL))
                    {
                        cutPoints.append(p1);
                    }
                    else
                    {
                        cutPoints.append((1-alpha)*p0 + alpha*p1);
                    }
                }
            }
            cutPoints.shrink();

            forAll(surf, faceI)
            {
                const labelList& fEdges = faceEdges[faceI];

                edge cutEdge(-1,-1);
                label nFound = 0;
                forAll(fEdges, fEI)
                {
                    label edgeI = fEdges[fEI];

                    if (edgePoint[edgeI] != -1)
                    {
                        label cutPt = edgePoint[edgeI];
                        if (nFound == 0)
                        {
                            cutEdge[0] = cutPt;
                            nFound++;
                        }
                        else
                        {
                            if
                            (
                                mag(cutPoints[cutEdge[0]] - cutPoints[cutPt])
                                > SMALL
                            )
                            {
                                cutEdge[1] = cutPt;
                                nFound++;
                                break;
                            }
                        }
                    }
                }
                if (nFound == 2)
                {
                    cutEdges.append(cutEdge);
                }
            }
            cutEdges.shrink();
        }
        else if (intersectionMethod == "surface")
        {
            const dictionary sDict = surfaceDict.subDict("surfaceCoeffs");

            const word name  = sDict.lookup
            (
                "name"
             );

            autoPtr<IOobject> namedIO
            (
                new IOobject
                (
                    "dummy",                    // dummy name
                    runTime.time().constant(),  // directory
                    "triSurface",               // instance
                    runTime.time(),             // registry
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                 )
            );
            namedIO().rename(name);

            autoPtr<searchableSurface> searchablePtr
            (
                searchableSurface::New
                (
                    sDict.lookup("type"),
                    namedIO(),
                    sDict
                 )
             );


            scalar maxLength = -GREAT;
            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];
                scalar edgeLength = e.mag(points);

                if (edgeLength > maxLength)
                {
                    maxLength = edgeLength;
                }
            }
            boundBox bb = searchablePtr().bounds();
            //expand bounding box
            bb.min() -= point(maxLength,maxLength,maxLength);
            bb.max() += point(maxLength,maxLength,maxLength);

            DynamicList<label> markedEdges(edges.size()/10);
            forAll(edges, edgeI)
            {
                const edge& e = edges[edgeI];
                if (bb.contains(points[e[0]]) || bb.contains(points[e[1]]))
                {
                    markedEdges.append(edgeI);
                }
            }
            markedEdges.shrink();

            pointField start(markedEdges.size());
            pointField end(markedEdges.size());

            forAll(markedEdges, i)
            {
                label edgeI = markedEdges[i];
                const edge& e = edges[edgeI];
                start[i] = points[e[0]];
                end[i] = points[e[1]];
            }

            List<pointIndexHit> pHits;

            searchablePtr().findLineAny
            (
                start,
                end,
                pHits
             );

            labelList edgePoint(edges.size(), -1);
            forAll(markedEdges, i)
            {
                if (pHits[i].hit())
                {
                    label edgeI = markedEdges[i];
                    edgePoint[edgeI] = cutPoints.size();

                    const edge& e = edges[edgeI];

                    const point& p0 = points[e[0]];
                    const point& p1 = points[e[1]];

                    scalar num = mag(pHits[i].hitPoint() - p0);
                    scalar den = mag(p1 - p0);
                    scalar alpha = num/den;

                    if (alpha < SMALL)
                    {
                        cutPoints.append(p0);
                    }
                    else if (alpha >= (1.0 - SMALL))
                    {
                        cutPoints.append(p1);
                    }
                    else
                    {
                        cutPoints.append(pHits[i].hitPoint());
                    }
                }
            }
            cutPoints.shrink();

            forAll(surf, faceI)
            {
                const labelList& fEdges = faceEdges[faceI];

                edge cutEdge(-1,-1);
                label nFound = 0;
                forAll(fEdges, fEI)
                {
                    label edgeI = fEdges[fEI];

                    if (edgePoint[edgeI] != -1)
                    {
                        label cutPt = edgePoint[edgeI];
                        if (nFound == 0)
                        {
                            cutEdge[0] = cutPt;
                            nFound++;
                        }
                        else
                        {
                            if
                            (
                                mag(cutPoints[cutEdge[0]] - cutPoints[cutPt])
                                > SMALL
                            )
                            {
                                cutEdge[1] = cutPt;
                                nFound++;
                                break;
                            }
                        }
                    }
                }
                if (nFound == 2)
                {
                    cutEdges.append(cutEdge);
                }
            }
            cutEdges.shrink();
        }

        Info<<"Calculated "<<cutEdges.size()<<" cut edges"<<endl;

        featureEdgeMesh feMesh
        (
            IOobject
            (
                eName + ".eMesh",   // name
                runTime.constant(), // instance
                "triSurface",
                runTime,            // registry
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
             ),
            pointField(cutPoints),
            cutEdges
         );

        Info<< nl << "Writing intersection featureEdgeMesh to "
            << feMesh.objectPath() << endl;

        feMesh.regIOobject::write();

        if (writeVTK)
        {
            simpleVTKWriter
            (
                cutEdges,
                cutPoints
             ).write(eName + ".vtk");
        }
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
