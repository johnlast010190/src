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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2010-2016 Esi Ltd.

Description
    Finds a point inside closed manifold surface

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "triSurface/triSurface.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "triSurface/triSurfaceSearch/triSurfaceSearch.H"

using namespace Foam;

// Tolerance used for defining insideness
static const scalar defaultMinTol = 1E-3;
static const label defaultMaxIter = 7;

void getNextIntersections
(
    const indexedOctree<treeDataTriSurface>& octree,
    const point& start,
    const point& end,
    const vector& smallVec,
    DynamicList<pointIndexHit, 1, 1>& hits
)
{
    const vector dirVec(end-start);
    const scalar magSqrDirVec(magSqr(dirVec));

    // Initial perturbation amount
    vector perturbVec(smallVec);

    while (true)
    {
        // Start tracking from last hit.
        point pt = hits[hits.size()-1].hitPoint() + perturbVec;

        if (((pt-start)&dirVec) > magSqrDirVec)
        {
            return;
        }

        // See if any intersection between pt and end
        pointIndexHit inter = octree.findLine(pt, end);

        if (!inter.hit())
        {
            return;
        }

        // Check if already found this intersection
        bool duplicateHit = false;
        forAllReverse(hits, i)
        {
            if (hits[i].index() == inter.index())
            {
                duplicateHit = true;
                break;
            }
        }


        if (duplicateHit)
        {
            // Hit same triangle again. Increase perturbVec and try again.
            perturbVec *= 2;
        }
        else
        {
            // Proper hit
            hits.append(inter);
            // Restore perturbVec
            perturbVec = smallVec;
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("surface file");
    argList::validOptions.insert("minTol", "inside/outside tolerance");
    argList::validOptions.insert("outside", "check for outside point");
    argList::validOptions.insert("maxIter", "Maximum number search iterations");

    argList args(argc, argv);

    scalar minTol = defaultMinTol;
    args.optionReadIfPresent("minTol", minTol);
    label maxIter = defaultMaxIter;
    args.optionReadIfPresent("maxIter", maxIter);

    bool inside = true;
    if (args.optionFound("outside"))
    {
        inside = false;
    }

    Info<<"Reading surface: " << args[1] << endl;
    triSurface surf(args[1]);
    const labelListList& eFaces = surf.edgeFaces();

    label nSingleEdges = 0;
    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() == 1)
        {
            nSingleEdges++;
        }
    }

    label nMultEdges = 0;
    forAll(eFaces, edgeI)
    {
        const labelList& myFaces = eFaces[edgeI];

        if (myFaces.size() > 2)
        {
            nMultEdges++;
        }
    }

    Info<<"Checking of surface" << endl;
    if ((nSingleEdges != 0) || (nMultEdges != 0))
    {
        Info<< "Surface is not closed since not all edges connected to "
            << "two faces:" << endl
            << "    connected to one face : " << nSingleEdges << endl
            << "    connected to >2 faces : " << nMultEdges << endl;
    }
    else
    {
        Info<< "Surface is closed. All edges connected to two faces." << endl;
    }
    labelList faceZone;
    label numZones = surf.markZones(boolList(surf.nEdges(), false), faceZone);
    Pout<< "Number of unconnected parts : " << numZones << endl;

    if ((nSingleEdges != 0) || (nMultEdges != 0) || (numZones > 1))
    {
        WarningInFunction
            <<"Surface non-manifold which might produce "
            <<" unexpected behavour when finding a point inside the surface "
            <<endl;
    }

    treeBoundBox bb
    (
        vector(GREAT, GREAT, GREAT),
        vector(-GREAT, -GREAT, -GREAT)
     );

    forAll(surf, triI)
    {
        const labelledTri& f = surf[triI];

        forAll(f, fp)
        {
            label pointI = f[fp];
            bb.min() = ::Foam::min(bb.min(), surf.points()[pointI]);
            bb.max() = ::Foam::max(bb.max(), surf.points()[pointI]);
        }
    }

    scalar minDim = bb.minDim();
    if (minDim < SMALL)
    {
        FatalErrorIn("surfaceInsidePoint")
            <<" Surface is 2D "
            <<" with bounding box: "<<bb
            << abort(FatalError);
    }
    else if (minTol >= 0.5*minDim)
    {
        WarningInFunction
            <<"Default tolerance: " << minTol
            <<" is greater than  bounding box minimum dimension: " << minDim
            <<endl;
    }

    // Random number generator.
    Random rndGen(65431);

    DynamicList<treeBoundBox> subBoxes(1);
    subBoxes.append(bb);

    // Slightly extended bb. Slightly off-centred just so on symmetric
    // geometry there are less face/edge aligned items.
    bb = bb.extend(rndGen, 1E-4);
    bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    triSurfaceSearch searchSelectSurf
    (
        surf,
        indexedOctree<treeDataTriSurface>::perturbTol(),
        10
     );

    const indexedOctree<treeDataTriSurface>& octree =
        searchSelectSurf.tree();

    // Work array
    DynamicList<pointIndexHit, 1, 1> hits;
    point midBbPt = bb.midpoint();
    scalar magBb = bb.mag();
    pointField end(3);

    end[0] = midBbPt + point(magBb, 0., 0.);
    end[1] = midBbPt + point(0., magBb, 0.);
    end[2] = midBbPt + point(0., 0., magBb);

    pointField start(3);

    bool foundPt = false;
    bool terminate = false;
    label nIter = 0;
    point insidePt = point(GREAT, GREAT, GREAT);

    while (!foundPt)
    {
        Info<<"Search iteration: "<<nIter
            <<" Searching using sub bounding box size: "
            << subBoxes[0].minDim()<<endl;

        scalar furthestDist = -GREAT;

        forAll(subBoxes, boxI)
        {
            const treeBoundBox& sB = subBoxes[boxI];

            const vectorField dirVec(end-start);
            const scalarField magSqrDirVec(magSqr(dirVec));
            const vectorField smallVec
            (
                indexedOctree<treeDataTriSurface>::perturbTol()*dirVec
                + vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL)
             );
            start = sB.midpoint();

            List<List<pointIndexHit>> info;
            info.setSize(start.size());
            label nIn =0;
//            label nOut = 0;

            forAll(start, pointI)
            {
                // See if any intersection between pt and end
                pointIndexHit inter =
                    octree.findLine(start[pointI], end[pointI]);

                if (inter.hit())
                {
                    hits.clear();
                    hits.append(inter);

                    getNextIntersections
                    (
                        octree,
                        start[pointI],
                        end[pointI],
                        smallVec[pointI],
                        hits
                     );

                    info[pointI].transfer(hits);
                }
                else
                {
                    info[pointI].clear();
                }
                if (info[pointI].size() % 2 == 0)
                {
//                    nOut++;
                }
                else
                {
                    nIn++;
                }
            }
            if
            (
                (nIn == 3 && inside)
               || (nIn == 0 && !inside)
            )
            {
                pointIndexHit hit = octree.findNearest(start[0], magBb);
                if (hit.hit())
                {
                    scalar hitDist = mag(hit.hitPoint() - start[0]);
                    if (hitDist > minTol)
                    {
                        foundPt = true;
                        if (hitDist > furthestDist)
                        {
                            insidePt = start[0];
                            furthestDist = hitDist;
                        }
                    }
                }
            }
        }

        if (foundPt)
        {
            if (inside)
            {
                Info<<"Found Inside point: "<<insidePt<<endl;
            }
            else
            {
                Info<<"Found outside point: "<<insidePt<<endl;
            }
            Info<<"Distance to surface: "<<furthestDist<<endl;
        }


        DynamicList<treeBoundBox> newSubBoxes(subBoxes.size()*8);
        forAll(subBoxes, boxI)
        {
            for (direction oct = 0; oct < 8; oct++)
            {
                treeBoundBox ssBox = subBoxes[boxI].subBbox(oct);
                newSubBoxes.append(ssBox);
                if (ssBox.minDim() < minTol)
                {
                    terminate = true;
                }
            }
        }
        if (terminate || nIter == maxIter)
        {
            Info<<"Could not find inside/outside point. "<<endl;
            if (!inside)
            {
                Info<<"Projecting from bounding box dimensions. "
                    <<"Outside point: "<<end[0]<<endl;
            }
            break;
        }

        subBoxes.clear();
        newSubBoxes.shrink();
        subBoxes.transfer(newSubBoxes);
        nIter++;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
