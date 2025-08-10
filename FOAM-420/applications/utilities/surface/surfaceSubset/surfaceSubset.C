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
    (c) 2015 OpenCFD Ltd.
    (c) 2010-2016 Esi Ltd.

Application
    surfaceSubset

Group
    grpSurfaceUtilities

Description
    A surface analysis tool which sub-sets the triSurface
    to choose only a part of interest. Based on subsetMesh.

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"
#include "global/argList/argList.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "primitives/bools/Switch/Switch.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"
#include "meshes/boundBox/boundBox.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataTriSurface.H"
#include "primitives/random/Random/Random.H"
#include "primitives/strings/wordRes/wordReListMatcher.H"
#include "algorithms/indexedOctree/volumeType.H"
#include "meshes/primitiveShapes/plane/plane.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("surfaceSubsetDict");
    argList::validArgs.append("surfaceFile");
    argList::validArgs.append("output surfaceFile");
    argList args(argc, argv);

    Info<< "Reading dictionary " << args[1] << " ..." << endl;
    IFstream dictFile(args[1]);
    dictionary meshSubsetDict(dictFile);

    Info<< "Reading surface " << args[2] << " ..." << endl;
    triSurface surf1(args[2]);

    const fileName outFileName(args[3]);


    Info<< "Original:" << endl;
    surf1.writeStats(Info);
    Info<< endl;

    List<wordRe> markedPatches
    (
        meshSubsetDict.lookup("patches")
    );

    labelList markedPoints
    (
        meshSubsetDict.lookup("localPoints")
    );

    labelList markedEdges
    (
        meshSubsetDict.lookup("edges")
    );

    labelList markedFaces
    (
        meshSubsetDict.lookup("faces")
    );

    pointField markedZone
    (
        meshSubsetDict.lookup("zone")
    );

    if (markedZone.size() && markedZone.size() != 2)
    {
        FatalErrorInFunction
            << "zone specification should be two points, min and max of "
            << "the boundingbox" << endl
            << "zone:" << markedZone
            << exit(FatalError);
    }

    Switch addFaceNeighbours
    (
        meshSubsetDict.lookup("addFaceNeighbours")
    );

    const bool invertSelection =
        meshSubsetDict.lookupOrDefault("invertSelection", false);

    // Mark the cells for the subset

    // Faces to subset
    boolList facesToSubset(surf1.size(), false);


    // pick up faces in named patches
    if (markedPatches.size())
    {

        wordReListMatcher regionMatch(markedPatches);
        const geometricSurfacePatchList& tsps(surf1.patches());

        DynamicList<label> markRegionIndices(0);

        // collect marked region indices
        Info<< "Selecting patches:";
        forAll(tsps, patchI)
        {
            if (regionMatch.match(tsps[patchI].name()))
            {
                markRegionIndices.append(tsps[patchI].index());
                Info<< " " << tsps[patchI].name();
            }
        }

        Info<< endl;
        markRegionIndices.shrink();

        // mark faces with matching indices
        // by adding faces to facesToSubset
        forAll(surf1, faceI)
        {
            const labelledTri& ct(surf1[faceI]);

            forAll(markRegionIndices, ii)
            {
                if (ct.region() == markRegionIndices[ii])
                {
                    facesToSubset[faceI] = true;
                }
            }
        }
    }

    //
    // pick up faces connected to "localPoints"
    //

    if (markedPoints.size())
    {
        Info<< "Found " << markedPoints.size() << " marked point(s)." << endl;

        // pick up cells sharing the point

        forAll(markedPoints, pointi)
        {
            if
            (
                markedPoints[pointi] < 0
             || markedPoints[pointi] >= surf1.nPoints()
            )
            {
                FatalErrorInFunction
                    << "localPoint label " << markedPoints[pointi]
                    << "out of range."
                    << " The mesh has got "
                    << surf1.nPoints() << " localPoints."
                    << exit(FatalError);
            }

            const labelList& curFaces =
                surf1.pointFaces()[markedPoints[pointi]];

            forAll(curFaces, i)
            {
                facesToSubset[curFaces[i]] =  true;
            }
        }
    }



    //
    // pick up faces connected to "edges"
    //

    if (markedEdges.size())
    {
        Info<< "Found " << markedEdges.size() << " marked edge(s)." << endl;

        // pick up cells sharing the edge

        forAll(markedEdges, edgeI)
        {
            if
            (
                markedEdges[edgeI] < 0
             || markedEdges[edgeI] >= surf1.nEdges()
            )
            {
                FatalErrorInFunction
                    << "edge label " << markedEdges[edgeI]
                    << "out of range."
                    << " The mesh has got "
                    << surf1.nEdges() << " edges."
                    << exit(FatalError);
            }

            const labelList& curFaces = surf1.edgeFaces()[markedEdges[edgeI]];

            forAll(curFaces, i)
            {
                facesToSubset[curFaces[i]] =  true;
            }
        }
    }


    //
    // pick up faces with centre inside "zone"
    //

    if (markedZone.size() == 2)
    {
        const point& min = markedZone[0];
        const point& max = markedZone[1];

        Info<< "Using zone min:" << min << " max:" << max << endl;

        forAll(surf1, facei)
        {
            const point centre = surf1[facei].centre(surf1.points());

            if
            (
                (centre.x() >= min.x())
             && (centre.y() >= min.y())
             && (centre.z() >= min.z())
             && (centre.x() <= max.x())
             && (centre.y() <= max.y())
             && (centre.z() <= max.z())
            )
            {
                facesToSubset[facei] = true;
            }
        }
    }


    //
    // pick up faces on certain side of surface
    //

    if (meshSubsetDict.found("surface"))
    {
        const dictionary& surfDict = meshSubsetDict.subDict("surface");

        fileName surfName(surfDict.lookup("name"));

        Switch outside(surfDict.lookup("outside"));

        if (outside)
        {
            Info<< "Selecting all triangles with centre outside surface "
                << surfName << endl;
        }
        else
        {
            Info<< "Selecting all triangles with centre inside surface "
                << surfName << endl;
        }

        // Read surface to select on
        triSurface selectSurf(surfName);

        triSurfaceSearch searchSelectSurf
        (
            selectSurf,
            indexedOctree<treeDataTriSurface>::perturbTol(),
            8
        );

        const indexedOctree<treeDataTriSurface>& selectTree =
            searchSelectSurf.tree();

        // Check if face (centre) is in outside or inside.
        forAll(facesToSubset, facei)
        {
            if (!facesToSubset[facei])
            {
                const point fc(surf1[facei].centre(surf1.points()));

                volumeType t = selectTree.getVolumeType(fc);

                if
                (
                    outside
                  ? (t == volumeType::OUTSIDE)
                  : (t == volumeType::INSIDE)
                )
                {
                    facesToSubset[facei] = true;
                }
            }
        }
    }


    if (meshSubsetDict.found("plane"))
    {
        const dictionary& planeDict = meshSubsetDict.subDict("plane");

        const plane pl(planeDict);
        const scalar distance(readScalar(planeDict.lookup("distance")));
        const scalar cosAngle(readScalar(planeDict.lookup("cosAngle")));

        // Select all triangles that are close to the plane and
        // whose normal aligns with the plane as well.

        forAll(surf1.faceCentres(), facei)
        {
            const point& fc = surf1.faceCentres()[facei];
            const point& nf = surf1.faceNormals()[facei];

            if (pl.distance(fc) < distance && mag(pl.normal() & nf) > cosAngle)
            {
                facesToSubset[facei] = true;
            }
        }
    }



    //
    // pick up specified "faces"
    //

    // Number of additional faces picked up because of addFaceNeighbours
    label nFaceNeighbours = 0;

    if (markedFaces.size())
    {
        Info<< "Found " << markedFaces.size() << " marked face(s)." << endl;

        // Check and mark faces to pick up
        forAll(markedFaces, facei)
        {
            if
            (
                markedFaces[facei] < 0
             || markedFaces[facei] >= surf1.size()
            )
            {
                FatalErrorInFunction
                    << "Face label " << markedFaces[facei] << "out of range."
                    << " The mesh has got "
                    << surf1.size() << " faces."
                    << exit(FatalError);
            }

            // Mark the face
            facesToSubset[markedFaces[facei]] = true;

            // mark its neighbours if requested
            if (addFaceNeighbours)
            {
                const labelList& curFaces =
                    surf1.faceFaces()[markedFaces[facei]];

                forAll(curFaces, i)
                {
                    label facei = curFaces[i];

                    if (!facesToSubset[facei])
                    {
                        facesToSubset[facei] =  true;
                        nFaceNeighbours++;
                    }
                }
            }
        }
    }

    if (addFaceNeighbours)
    {
        Info<< "Added " << nFaceNeighbours
            << " faces because of addFaceNeighbours" << endl;
    }


    if (invertSelection)
    {
        Info<< "Inverting selection." << endl;

        forAll(facesToSubset, i)
        {
            facesToSubset[i] = facesToSubset[i] ? false : true;
        }
    }


    // Create subsetted surface
    labelList pointMap;
    labelList faceMap;
    bool compactPatches
    (
        meshSubsetDict.lookupOrDefault<bool>("compactEmptyPatches", false)
    );

    triSurface surf2
    (
        surf1.subsetMesh(facesToSubset, pointMap, faceMap, compactPatches)
    );

    Info<< "Subset:" << endl;
    surf2.writeStats(Info);
    Info<< endl;

    Info<< "Writing surface to " << outFileName << endl;

    surf2.write(outFileName);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
