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
    (c) 2019-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "surface/cuttingTriSurface/cuttingTriSurface.H"
#include "meshes/primitiveMesh/primitiveMesh.H"
#include "meshes/primitiveShapes/line/linePointRef.H"
#include "meshTools/meshTools.H"
#include "algorithms/polygonTriangulate/polygonTriangulate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Set values for what is close to zero and what is considered to
// be positive (and not just rounding noise)
//! \cond localScope
const Foam::scalar zeroish  = Foam::SMALL;
//! \endcond

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::indexedOctree<Foam::treeDataTriSurface>&
Foam::cuttingTriSurface::tree() const
{
    treePtr_.reset();

    if (treePtr_.empty())
    {
        // Calculate bb without constructing local point numbering.
        treeBoundBox bb
        (
            surface_.triSurface::points(),
            surface_.triSurface::meshPoints()
        );

        if (surface_.size())
        {
            label nPoints;
            PatchTools::calcBounds(surface_, bb, nPoints);

            if (nPoints != surface_.points().size())
            {
                WarningInFunction
                    << "Surface does not have compact point numbering."
                    << " Of " << surface().points().size()
                    << " only " << nPoints
                    << " are used. This might give problems in some routines."
                    << endl;
            }

            // Random number generator. Bit dodgy since not exactly random ;-)
            Random rndGen(65431);

            // Slightly extended bb. Slightly off-centred just so on symmetric
            // geometry there are less face/edge aligned items.
            bb = bb.extend(rndGen, 1e-4);
            bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        }

        const scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
        indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

        treePtr_.reset
        (
            new indexedOctree<treeDataTriSurface>
            (
                treeDataTriSurface(false, surface_, tolerance_),
                bb,
                10,             // maxLevel
                10,             // leafsize
                3.0             // duplicity
            )
        );

        indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
    }

    return treePtr_();
}


void Foam::cuttingTriSurface::calcEdgeCuts
(
    const primitiveMesh& mesh,
    const bool planarCut,
    scalarField& edgeCuts
)
{
    if (surface_.size() == 0)
    {
        return;
    }

    const indexedOctree<treeDataTriSurface>& octree = tree();

    if (planarCut)
    {
        //calculate face 0 direction
        vector fN0 = surface_.faceAreas()[0];

        vector aveCtr = Zero;
        vector aveNorm = Zero;
        scalar totArea = Zero;

        forAll(surface_, facei)
        {
            point fN = surface_.faceAreas()[facei];
            point fC = surface_.faceCentres()[facei];

            scalar area = mag(fN);
            aveCtr +=  area*fC;
            totArea += area;

            if ((fN & fN0) > 0)
            {
                aveNorm +=  fN;
            }
            else
            {
                aveNorm -=  fN;
            }
        }

        if (totArea > SMALL)
        {
            aveCtr /= totArea;
            aveNorm /= totArea;
        }

        plane cutPlane(aveCtr, aveNorm);

        forAll(mesh.edges(), edgei)
        {
            const edge e = mesh.edges()[edgei];

            point start = mesh.points()[e[0]];
            point end = mesh.points()[e[1]];

            scalar dP0 = (start - cutPlane.refPoint()) & cutPlane.normal();
            scalar dP1 = (end - cutPlane.refPoint()) & cutPlane.normal();
            if ((dP0 <= 0.0 && dP1 > 0.0) || (dP1 <= 0.0 && dP0 > 0.0))
            {
                scalar alpha = cutPlane.lineIntersect(linePointRef(start, end));
                edgeCuts[edgei] = alpha;
            }
        }

        scalar minEdgeLength = GREAT;
        forAll(mesh.cells(), celli)
        {
            const labelList& cEdges = mesh.cellEdges()[celli];

            bool possibleCut = false;
            forAll(cEdges, cEI)
            {
                label edgei = cEdges[cEI];
                if (edgeCuts[edgei] > -1)
                {
                    possibleCut = true;
                }
            }
            if (possibleCut)
            {
                forAll(cEdges, cEI)
                {
                    label edgei = cEdges[cEI];
                    const edge e = mesh.edges()[edgei];
                    scalar eLen = e.mag(mesh.points());
                    minEdgeLength = min(minEdgeLength, eLen);
                }
            }
        }
        reduce(minEdgeLength, minOp<scalar>());

        //set search distance
        scalar searchDist = 0.2*pow(minEdgeLength,2);

        forAll(mesh.edges(), edgei)
        {
            scalar alpha = edgeCuts[edgei];
            if (alpha > -1)
            {
                const edge e = mesh.edges()[edgei];

                point start = mesh.points()[e[0]];
                point end = mesh.points()[e[1]];
                point hPoint = (1-alpha)*start + alpha*end;

                pointIndexHit info = octree.findNearest(hPoint, searchDist);
                if (!info.hit())
                {
                    edgeCuts[edgei] = -1;
                }
            }
        }
    }
    else
    {
        forAll(mesh.edges(), edgei)
        {
            const edge e = mesh.edges()[edgei];

            point start = mesh.points()[e[0]];
            point end = mesh.points()[e[1]];
            pointIndexHit info = octree.findLineAny(start, end);
            if (info.hit())
            {
                scalar eLen = e.mag(mesh.points());
                edgeCuts[edgei] = mag(info.hitPoint()-start)/eLen;
            }
        }
    }
}


void Foam::cuttingTriSurface::calcCutCells
(
    const primitiveMesh& mesh,
    const scalarField& edgeCuts,
    const labelUList& cellIdLabels
)
{
    const labelListList& cellEdges = mesh.cellEdges();

    label listSize = cellEdges.size();
    if (notNull(cellIdLabels))
    {
        listSize = cellIdLabels.size();
    }

    meshCells_.setSize(listSize);
    label cutcelli(0);

    for (label listI = 0; listI < listSize; ++listI)
    {
        label celli = listI;

        if (notNull(cellIdLabels))
        {
            celli = cellIdLabels[listI];
        }

        const labelList& cEdges = cellEdges[celli];
        label nCutEdges = 0;
        forAll(cEdges, i)
        {
            label edgei = cEdges[i];

            if (edgeCuts[edgei] > -1)
            {
                nCutEdges++;

                if (nCutEdges > 2)
                {
                    meshCells_[cutcelli++] = celli;
                    break;
                }
            }
        }
    }

    // Set correct list size
    meshCells_.setSize(cutcelli);
}


void Foam::cuttingTriSurface::intersectEdges
(
    const primitiveMesh& mesh,
    const scalarField& edgeCuts,
    List<label>& edgePoint
)
{
    const edgeList& edges = mesh.edges();
    const pointField& points = mesh.points();

    // Per edge -1 or the label of the intersection point
    edgePoint.setSize(edges.size());

    DynamicList<point> dynCuttingPoints(4*meshCells_.size());

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];

        if (edgeCuts[edgeI] > -1)
        {
            // Edge is cut
            edgePoint[edgeI] = dynCuttingPoints.size();

            const point& p0 = points[e[0]];
            const point& p1 = points[e[1]];

            scalar alpha = edgeCuts[edgeI];

            if (alpha < zeroish)
            {
                dynCuttingPoints.append(p0);
            }
            else if (alpha >= (1.0 - zeroish))
            {
                dynCuttingPoints.append(p1);
            }
            else
            {
                dynCuttingPoints.append((1-alpha)*p0 + alpha*p1);
            }
        }
        else
        {
            edgePoint[edgeI] = -1;
        }
    }

    this->storedPoints().transfer(dynCuttingPoints);
}


bool Foam::cuttingTriSurface::walkCell
(
    const primitiveMesh& mesh,
    const labelUList& edgePoint,
    const label celli,
    label startEdgeI,
    DynamicList<label>& faceVerts
)
{
    label facei = -1;
    label edgeI = startEdgeI;
    label startFaceI = -1;

    label nIter = 0;
    labelHashSet visitedEdges;
    labelHashSet triedStartEdges;

    bool firstPass = true;

    const cell& c  = mesh.cells()[celli];
    edgeList cEdges = c.edges(mesh.faces());
    label nCutEdges = 0;

    forAll(cEdges, eI)
    {
        edge e = cEdges[eI];
        label meshEdgeI = meshTools::findEdge
              (
                  mesh.edges(),
                  mesh.pointEdges()[e[0]],
                  e[0],
                  e[1]
               );
        if (edgePoint[meshEdgeI] != -1)
        {
            nCutEdges++;
        }
    }

    faceVerts.clear();
    do
    {
        faceVerts.append(edgePoint[edgeI]);
        visitedEdges.insert(edgeI);

        // Cross edge to other face
        facei = meshTools::otherFace(mesh, celli, facei, edgeI);

        if (startFaceI == -1)
        {
            startFaceI = facei;
        }

        // Find next cut edge on face.
        const labelList& fEdges = mesh.faceEdges()[facei];

        label nextEdgeI = -1;

        //Note: here is where we should check for whether there are more
        // than 2 intersections with the face (warped/non-convex face).
        // If so should e.g. decompose the cells on both faces and redo
        // the calculation.

        forAll(fEdges, i)
        {
            label edge2I = fEdges[i];

            if
            (
                edge2I != edgeI && edgePoint[edge2I] != -1
                && !visitedEdges.found(edge2I)
            )
            {
                nextEdgeI = edge2I;
                break;
            }
        }

        if (nextEdgeI == -1)
        {
            //could not find edge (assuming finished loop or try new loop)
            nextEdgeI = startEdgeI;
        }

        edgeI = nextEdgeI;

        nIter++;

        if (nIter > 1000)
        {
            WarningInFunction
                << "Did not find closed walk along surface of cell " << celli
                << " starting from edge " << startEdgeI
                << " in " << nIter << " iterations." << nl
                << "Collected cutPoints so far:" << faceVerts
                << endl;

            return false;
        }

        if (edgeI == startEdgeI && faceVerts.size() != nCutEdges && firstPass)
        {
            firstPass = false;
            faceVerts.clear();
            visitedEdges.clear();
            facei = startFaceI;
            edgeI = startEdgeI;
        }
        else if (edgeI == startEdgeI)
        {
            if (faceVerts.size() != nCutEdges)
            {
                triedStartEdges.insert(edgeI);
                //try another starting edge
                const labelList& cEdges = mesh.cellEdges()[celli];

                label newEdgeI = -1;
                forAll(cEdges, cEdgeI)
                {
                    label eI = cEdges[cEdgeI];

                    if (edgePoint[eI] != -1 && !triedStartEdges.found(eI))
                    {
                        newEdgeI = eI;
                        break;
                    }
                }

                if (newEdgeI == -1)
                {
                    // Did not find another cut edge on faceI. Do what?
                    WarningInFunction
                        << "Did not find closed walk along surface of cell "
                        << celli
                        << " starting from edge " << startEdgeI
                        << " in " << nIter << " iterations." << nl
                        << "Collected cutPoints so far:" << faceVerts
                        << endl;

                    return false;
                }
                else
                {
                    facei = -1;
                    edgeI = startEdgeI;
                    startEdgeI = newEdgeI;
                    startFaceI = -1;
                    firstPass = true;
                    faceVerts.clear();
                    visitedEdges.clear();
                }
            }
            else
            {
                break;
            }
        }
    } while (true);

    if (faceVerts.size() >= 3)
    {
        return true;
    }
    else
    {
        WarningInFunction
            << "Did not find closed walk along surface of cell " << celli
            << " starting from edge " << startEdgeI << nl
            << "Collected cutPoints so far:" << faceVerts
            << endl;

        return false;
    }
}


void Foam::cuttingTriSurface::walkCellCuts
(
    const primitiveMesh& mesh,
    const bool triangulate,
    const labelUList& edgePoint
)
{
    const pointField& cutPoints = this->points();

    // use dynamic lists to handle triangulation and/or missed cuts
    DynamicList<face>  dynCutFaces(meshCells_.size());
    DynamicList<label> dynCutCells(meshCells_.size());

    // scratch space for calculating the face vertices
    DynamicList<label> faceVerts(10);

    // Create a triangulation engine
    polygonTriangulate triEngine;

    forAll(meshCells_, i)
    {
        label celli = meshCells_[i];

        // Find the starting edge to walk from.
        const labelList& cEdges = mesh.cellEdges()[celli];

        label startEdgeI = -1;

        forAll(cEdges, cEdgeI)
        {
            label edgeI = cEdges[cEdgeI];

            if (edgePoint[edgeI] != -1)
            {
                startEdgeI = edgeI;
                break;
            }
        }

        // Check for the unexpected ...
        if (startEdgeI == -1)
        {
            FatalErrorInFunction
                << "Cannot find cut edge for cut cell " << celli
                << abort(FatalError);
        }

        // Walk from starting edge around the circumference of the cell.
        bool okCut = walkCell
        (
            mesh,
            edgePoint,
            celli,
            startEdgeI,
            faceVerts
        );

        if (okCut)
        {
            face f(faceVerts);

            // the cut faces are usually quite ugly, so optionally triangulate
            if (triangulate)
            {
                triEngine.triangulate
                (
                    UIndirectList<point>(cutPoints, f)
                );

                forAll(triEngine.triPoints(), triI)
                {
                    dynCutFaces.append(triEngine.triPoints(triI, f));
                    dynCutCells.append(celli);
                }
            }
            else
            {
                dynCutFaces.append(f);
                dynCutCells.append(celli);
            }
        }
    }

    // Orient face to point in the same direction as the triSurface normal
    const indexedOctree<treeDataTriSurface>& octree = tree();
    boolList pointsDone(cutPoints.size(), false);
    forAll(dynCutFaces, i)
    {
        face& f = dynCutFaces[i];
        forAll(f, fp)
        {
            pointsDone[f[fp]] = true;
        }
        point fC = f.centre(cutPoints);
        pointIndexHit hit = octree.findNearest(fC, GREAT);
        if (hit.hit())
        {
            point fN = surface_.faceAreas()[hit.index()];
            if ((f.areaNormal(cutPoints) & fN) > 0)
            {
                f.flip();
            }
        }
    }

    //Check for unused points and remove
    labelList pointMap(cutPoints.size(),-1);
    DynamicList<point> dynCutPoints(cutPoints.size());

    label nValidPts = 0;
    forAll(cutPoints, pointi)
    {
        if (pointsDone[pointi])
        {
            dynCutPoints.append(cutPoints[pointi]);
            pointMap[pointi] = nValidPts;
            nValidPts++;
        }
    }

    if (nValidPts != cutPoints.size())
    {
        forAll(dynCutFaces, facei)
        {
            face& f = dynCutFaces[facei];
            forAll(f,fp)
            {
                f[fp] = pointMap[f[fp]];
            }
        }
        this->storedPoints().transfer(dynCutPoints);
    }

    this->storedFaces().transfer(dynCutFaces);
    meshCells_.transfer(dynCutCells);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cuttingTriSurface::cuttingTriSurface(const triSurface& tri)
:
    surface_(tri),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol())
{}


Foam::cuttingTriSurface::cuttingTriSurface
(
    const triSurface& tri,
    const primitiveMesh& mesh,
    const bool triangulate,
    const bool planarCut,
    const labelUList& cellIdLabels
)
:
    surface_(tri),
    tolerance_(indexedOctree<treeDataTriSurface>::perturbTol())
{
    reCut(mesh, triangulate, planarCut, cellIdLabels);
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cuttingTriSurface::reCut
(
    const primitiveMesh& mesh,
    const bool triangulate,
    const bool planarCut,
    const labelUList& cellIdLabels
)
{
    MeshStorage::clear();
    meshCells_.clear();

    scalarField edgeCuts(mesh.nEdges(), -1);
    calcEdgeCuts(mesh, planarCut, edgeCuts);

    // Determine cells that are (probably) cut.
    calcCutCells(mesh, edgeCuts, cellIdLabels);

    // Determine cutPoints and return list of edge cuts.
    // per edge -1 or the label of the intersection point
    labelList edgePoint;
    intersectEdges(mesh, edgeCuts, edgePoint);

    // Do topological walk around cell to find closed loop.
    walkCellCuts(mesh, triangulate, edgePoint);
}


void Foam::cuttingTriSurface::remapFaces
(
    const labelUList& faceMap
)
{
    // Recalculate the cells cut
    if (notNull(faceMap) && faceMap.size())
    {
        MeshStorage::remapFaces(faceMap);

        List<label> newCutCells(faceMap.size());
        forAll(faceMap, facei)
        {
            newCutCells[facei] = meshCells_[faceMap[facei]];
        }
        meshCells_.transfer(newCutCells);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cuttingTriSurface::operator=(const cuttingTriSurface& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    static_cast<MeshStorage&>(*this) = rhs;
    meshCells_ = rhs.meshCells();
}


// ************************************************************************* //
