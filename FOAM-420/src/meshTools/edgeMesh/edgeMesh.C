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
    (c) 2010-2012 Esi Ltd.
    (c) 2015-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "edgeMesh/edgeMesh.H"
#include "meshes/meshTools/mergePoints.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/runTimeSelection/memberFunctions/addToMemberFunctionSelectionTable.H"
#include "containers/HashTables/Map/Map.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "meshes/meshTools/matchPoints.H"
#include "meshes/meshShapes/edge/EdgeMap.H"
#include "containers/Lists/PackedList/PackedBoolList.H"
#include "fields/Fields/transformField/transformField.H"
#include "global/constants/mathematical/mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(edgeMesh, 0);
    defineRunTimeSelectionTable(edgeMesh, fileExtension);
    defineMemberFunctionSelectionTable(edgeMesh,write,fileExtension);
}


Foam::wordHashSet Foam::edgeMesh::readTypes()
{
    return wordHashSet(fileExtensionConstructorTable_());
}


Foam::wordHashSet Foam::edgeMesh::writeTypes()
{
    return wordHashSet(*writefileExtensionMemberFunctionTablePtr_);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::edgeMesh::canReadType
(
    const word& ext,
    const bool verbose
)
{
    return checkSupport
    (
        readTypes(),
        ext,
        verbose,
        "reading"
   );
}


bool Foam::edgeMesh::canWriteType
(
    const word& ext,
    const bool verbose
)
{
    return checkSupport
    (
        writeTypes(),
        ext,
        verbose,
        "writing"
    );
}


bool Foam::edgeMesh::canRead
(
    const fileName& name,
    const bool verbose
)
{
    word ext = name.ext();
    if (ext == "gz")
    {
        ext = name.lessExt().ext();
    }
    return canReadType(ext, verbose);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::edgeMesh::calcPointEdges() const
{
    if (pointEdgesPtr_.valid())
    {
        FatalErrorInFunction
            << "pointEdges already calculated." << abort(FatalError);
    }

    pointEdgesPtr_.reset(new labelListList(points_.size()));
    labelListList& pointEdges = pointEdgesPtr_();

    invertManyToMany(pointEdges.size(), edges_, pointEdges);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::edgeMesh::edgeMesh()
:
    fileFormats::edgeMeshFormatsCore(),
    points_(0),
    edges_(0),
    pointEdgesPtr_(nullptr)
{}


Foam::edgeMesh::edgeMesh
(
    const pointField& points,
    const edgeList& edges
)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(points),
    edges_(edges),
    pointEdgesPtr_(nullptr)
{}


Foam::edgeMesh::edgeMesh
(
    const Xfer<pointField>& pointLst,
    const Xfer<edgeList>& edgeLst
)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(0),
    edges_(0),
    pointEdgesPtr_(nullptr)
{
    points_.transfer(pointLst());
    edges_.transfer(edgeLst());
}

Foam::edgeMesh::edgeMesh(Istream& is)
:
    fileFormats::edgeMeshFormatsCore(),
    points_(0),
    edges_(0),
    pointEdgesPtr_(nullptr)
{
    read(is);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::edgeMesh::~edgeMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::edgeMesh::clear()
{
    points_.clear();
    edges_.clear();
    pointEdgesPtr_.clear();
}


void Foam::edgeMesh::reset
(
    const Xfer<pointField>& pointLst,
    const Xfer<edgeList>& edgeLst
)
{
    // Take over new primitive data.
    // Optimized to avoid overwriting data at all
    if (notNull(pointLst))
    {
        points_.transfer(pointLst());
    }

    if (notNull(edgeLst))
    {
        edges_.transfer(edgeLst());

        // connectivity likely changed
        pointEdgesPtr_.clear();
    }
}


void Foam::edgeMesh::transfer(edgeMesh& mesh)
{
    points_.transfer(mesh.points_);
    edges_.transfer(mesh.edges_);
    pointEdgesPtr_ = mesh.pointEdgesPtr_;
}


Foam::Xfer<Foam::edgeMesh> Foam::edgeMesh::xfer()
{
    return xferMove(*this);
}


Foam::label Foam::edgeMesh::regions(labelList& edgeRegion) const
{
    edgeRegion.setSize(edges_.size());
    edgeRegion = -1;

    label startEdgeI = 0;
    label currentRegion = 0;

    while (true)
    {
        while (startEdgeI < edges_.size() && edgeRegion[startEdgeI] != -1)
        {
            startEdgeI++;
        }

        if (startEdgeI == edges_.size())
        {
            break;
        }

        // Found edge that has not yet been assigned a region.
        // Mark connected region with currentRegion starting at startEdgeI.

        edgeRegion[startEdgeI] = currentRegion;
        labelList edgesToVisit(1, startEdgeI);

        while (edgesToVisit.size())
        {
            // neighbours of current edgesToVisit
            DynamicList<label> newEdgesToVisit(edgesToVisit.size());

            // Mark all point connected edges with current region.
            forAll(edgesToVisit, i)
            {
                label edgeI = edgesToVisit[i];

                // Mark connected edges
                const edge& e = edges_[edgeI];

                forAll(e, fp)
                {
                    const labelList& pEdges = pointEdges()[e[fp]];

                    forAll(pEdges, pEdgeI)
                    {
                        label nbrEdgeI = pEdges[pEdgeI];

                        if (edgeRegion[nbrEdgeI] == -1)
                        {
                            edgeRegion[nbrEdgeI] = currentRegion;
                            newEdgesToVisit.append(nbrEdgeI);
                        }
                    }
                }
            }

            edgesToVisit.transfer(newEdgesToVisit);
        }

        currentRegion++;
    }
    return currentRegion;
}


void Foam::edgeMesh::scalePoints(const scalar scaleFactor)
{
    // avoid bad scaling
    if (scaleFactor > 0 && scaleFactor != 1.0)
    {
        storedPoints() *= scaleFactor;
    }
}


void Foam::edgeMesh::scalePoints(const vector scaleFactor)
{
    // avoid bad scaling
    if (mag(scaleFactor) > 0)
    {
        storedPoints().replace
        (
            vector::X,
            scaleFactor.x()*storedPoints().component(vector::X)
        );
        storedPoints().replace
        (
            vector::Y,
            scaleFactor.y()*storedPoints().component(vector::Y)
        );
        storedPoints().replace
        (
            vector::Z,
            scaleFactor.z()*storedPoints().component(vector::Z)
        );
    }
}


void Foam::edgeMesh::translatePoints(const vector translateVector)
{
    // avoid bad scaling
    if (mag(translateVector) > 0)
    {
        storedPoints() += translateVector;
    }
}


void Foam::edgeMesh::rotate(quaternion R)
{
    storedPoints() = transform(R, storedPoints());
}


void Foam::edgeMesh::rotate(Pair<vector> n1n2)
{
    n1n2[0] /= mag(n1n2[0]);
    n1n2[1] /= mag(n1n2[1]);

    if (mag(n1n2[0] ^ n1n2[1]) > SMALL)
    {
        tensor T = rotationTensor(n1n2[0], n1n2[1]);
        storedPoints() = transform(T, storedPoints());
    }
    else
    {
        WarningInFunction
            << " n1n2 vectors parallel, disabling rotation"
            << endl;
    }
}


void Foam::edgeMesh::doTransforms(const dictionary& dict)
{
    //Perform geometric transforms
    if (dict.found("transforms"))
    {
        PtrList<dictionary> transforms(dict.lookup("transforms"));

        forAll(transforms, dictI)
        {
            const dictionary& transformDict = transforms[dictI];

            const word type(transformDict.lookup("type"));
            if (type == "scale")
            {
                const vector scale(transformDict.lookup("scaleVec"));
                const point pt
                (
                    transformDict.lookupOrDefault("aboutPoint", vector(0, 0, 0))
                );
                translatePoints(-pt);
                scalePoints(scale);
                translatePoints(pt);
            }
            else if (type == "translate")
            {
                const vector v(transformDict.lookup("translateVec"));
                translatePoints(v);
            }
            else if (type == "rotate")
            {
                const point pt
                (
                    transformDict.lookupOrDefault("aboutPoint", vector(0, 0, 0))
                );
                translatePoints(-pt);
                if (transformDict.found("n1n2"))
                {
                    Info<< "Applying rotation method n1n2" << endl;
                    Pair<vector> n1n2
                    (
                        transformDict.lookup("n1n2")
                    );
                    rotate(n1n2);
                }
                else if (transformDict.found("rollPitchYaw"))
                {
                    Info<< "Applying rotation method rollPitchYaw" << endl;
                    vector v = transformDict.lookup("rollPitchYaw");
                    v *= Foam::constant::mathematical::pi/180.0;

                    quaternion R(quaternion::rotationSequence::XYZ, v);
                    rotate(R);
                }
                else if (transformDict.found("yawPitchRoll"))
                {
                    Info<< "Applying rotation method yawPitchRoll" << endl;
                    vector v = transformDict.lookup("yawPitchRoll");
                    v *= Foam::constant::mathematical::pi/180.0;

                    scalar yaw = v.x();
                    scalar pitch = v.y();
                    scalar roll = v.z();

                    quaternion R = quaternion(vector(0, 0, 1), yaw);
                    R *= quaternion(vector(0, 1, 0), pitch);
                    R *= quaternion(vector(1, 0, 0), roll);
                    rotate(R);
                }
                else
                {
                    WarningInFunction
                        << "Surface transform: " << type
                        << " not performed. Requires n1n2, rollPitchYaw or"
                        << " yawPitchRoll keyword to be present" << endl;
                }
                translatePoints(pt);
            }
        }
    }

    return;
}


void Foam::edgeMesh::mergePoints(const scalar mergeDist)
{
    pointField newPoints;
    labelList pointMap;

    const bool hasMerged = Foam::mergePoints
    (
        points_,
        mergeDist,
        false,
        pointMap,
        newPoints,
        vector::zero
    );

    if (hasMerged)
    {
        pointEdgesPtr_.clear();   // connectivity change

        points_.transfer(newPoints);

        forAll(edges_, edgeI)
        {
            edge& e = edges_[edgeI];

            e[0] = pointMap[e[0]];
            e[1] = pointMap[e[1]];
        }
    }

    this->mergeEdges();
}


void Foam::edgeMesh::mergeEdges()
{
    HashSet<edge, Hash<edge>> uniqEdges(2*edges_.size());
    PackedBoolList pointIsUsed(points_.size());

    label nUniqEdges = 0;
    label nUniqPoints = 0;
    forAll(edges_, edgeI)
    {
        const edge& e = edges_[edgeI];

        // Remove degenerate and repeated edges
        // - reordering (e[0] < e[1]) is not really necessary
        if (e[0] != e[1] && uniqEdges.insert(e))
        {
            if (nUniqEdges != edgeI)
            {
                edges_[nUniqEdges] = e;
            }
            edges_[nUniqEdges].sort();
            ++nUniqEdges;

            if (pointIsUsed.set(e[0], 1))
            {
                ++nUniqPoints;
            }
            if (pointIsUsed.set(e[1], 1))
            {
                ++nUniqPoints;
            }
        }
    }

    if (debug)
    {
        Info<< "Merging duplicate edges: "
            << (edges_.size() - nUniqEdges)
            << " edges will be deleted, "
            << (points_.size() - nUniqPoints)
            << " unused points will be removed." << endl;
    }

    if (nUniqEdges < edges_.size())
    {
        pointEdgesPtr_.clear(); // connectivity change
        edges_.setSize(nUniqEdges);  // truncate
    }

    if (nUniqPoints < points_.size())
    {
        pointEdgesPtr_.clear(); // connectivity change

        // build a oldToNew point-map and rewrite the points.
        // We can do this simultaneously since the point order is unchanged
        // and we are only effectively eliminating some entries.
        labelList pointMap(points_.size(), -1);

        label newId = 0;
        forAll(pointMap, pointi)
        {
            if (pointIsUsed[pointi])
            {
                pointMap[pointi] = newId;

                if (newId < pointi)
                {
                    // copy down
                    points_[newId] = points_[pointi];
                }
                ++newId;
            }
        }
        points_.setSize(newId);

        // Renumber edges - already sorted (above)
        forAll(edges_, edgeI)
        {
            edge& e = edges_[edgeI];

            e[0] = pointMap[e[0]];
            e[1] = pointMap[e[1]];
        }
    }

}


Foam::label Foam::edgeMesh::findEdge
(
    const List<edge>& allEdges,
    const labelListList& allPointEdges,
    const edge& otherE
)
{
    // all edges connected to otherE[0]
    const labelList& pEdges = allPointEdges[otherE[0]];

    forAll(pEdges, i)
    {
        const edge& e = allEdges[pEdges[i]];

        if (otherE[0] == e[0])
        {
            if (otherE[1] == e[1])
            {
                return pEdges[i];
            }

        }
        else if (otherE[0] == e[1])
        {
            if (otherE[1] == e[0])
            {
                return pEdges[i];
            }
        }
    }
    return -1;
}


// Merge into allEdges.
void Foam::edgeMesh::merge
(
    const scalar mergeDist,
    const List<edge>& subEdges,
    const pointField& subPoints,

    List<edge>& allEdges,
    pointField& allPoints,

    labelList& edgeConstructMap,
    labelList& pointConstructMap
)
{
    labelList subToAll;
    matchPoints
    (
        subPoints,
        allPoints,
        scalarField(subPoints.size(), mergeDist),   // match distance
        false,                                      // verbose
        pointConstructMap
    );

    label nOldAllPoints = allPoints.size();


    // Add all unmatched points
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    label allPointI = nOldAllPoints;
    forAll(pointConstructMap, pointI)
    {
        if (pointConstructMap[pointI] == -1)
        {
            pointConstructMap[pointI] = allPointI++;
        }
    }

    if (allPointI > nOldAllPoints)
    {
        allPoints.setSize(allPointI);

        forAll(pointConstructMap, pointI)
        {
            if (pointConstructMap[pointI] >= nOldAllPoints)
            {
                allPoints[pointConstructMap[pointI]] = subPoints[pointI];
            }
        }
    }


    // To check whether triangles are same we use pointEdges.
    labelListList allPointEdges;
    invertManyToMany(nOldAllPoints, allEdges, allPointEdges);


    // Add all unmatched triangles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label allEdgeI = allEdges.size();
    allEdges.setSize(allEdgeI + subEdges.size());

    edgeConstructMap.setSize(subEdges.size());

    forAll(subEdges, edgeI)
    {
        const edge& subEdge = subEdges[edgeI];

        // Get triangle in new numbering
        edge mappedEdge
        (
            pointConstructMap[subEdge[0]],
            pointConstructMap[subEdge[1]]
        );


        // Check if all points were already in edgeMesh
        bool fullMatch = true;

        if (mappedEdge[0] >= nOldAllPoints)
        {
            fullMatch = false;
        }
        else if (mappedEdge[1] >= nOldAllPoints)
        {
            fullMatch = false;
        }

        if (fullMatch)
        {
            // All three points are mapped to old points. See if same
            // triangle.
            label i = findEdge
            (
                allEdges,
                allPointEdges,
                mappedEdge
            );

            if (i == -1)
            {
                // Add
                edgeConstructMap[edgeI] = allEdgeI;
                allEdges[allEdgeI] = mappedEdge;
                allEdgeI++;
            }
            else
            {
                edgeConstructMap[edgeI] = i;
            }
        }
        else
        {
            // Add
            edgeConstructMap[edgeI] = allEdgeI;
            allEdges[allEdgeI] = mappedEdge;
            allEdgeI++;
        }
    }
    allEdges.setSize(allEdgeI);
}


// Does any part of edge overlap bb.
bool Foam::edgeMesh::overlaps
(
    const List<treeBoundBox>& bbs,
    const point& p0,
    const point& p1
)
{
    forAll(bbs, bbI)
    {
        const treeBoundBox& bb = bbs[bbI];

        treeBoundBox edgeBb(p0, p0);
        edgeBb.min() = min(edgeBb.min(), p1);
        edgeBb.max() = max(edgeBb.max(), p1);

        //- Exact test of triangle intersecting bb

        // Quick rejection. If whole bounding box of edge is outside cubeBb then
        // there will be no intersection.
        if (bb.overlaps(edgeBb))
        {
            // Check if one or both edge point inside
            if (bb.contains(p0) || bb.contains(p1))
            {
                // One or more points inside
                return true;
            }

            // Now we have the difficult case: all points are outside but
            // connecting edges might go through cube. Use fast intersection
            // of bounding box.
            point pInter;
            bool intersect = bb.intersects(p0, p1, pInter);

            if (intersect)
            {
                return true;
            }
        }
    }
    return false;
}


// reconstruct edges meshes
void Foam::edgeMesh::reconstruct(const scalar& mergeDist)
{
    boundBox bb
    (
        point(-GREAT,-GREAT,-GREAT),
        point(GREAT,GREAT,GREAT)
    );

    const List<treeBoundBox> tbb(1, treeBoundBox(bb));

    distribute(tbb,mergeDist);

    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            for (int slave=1; slave<Pstream::nProcs(); slave++)
            {
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                toSlave << *this;
            }
        }
        else
        {
            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
             );
             edgeMesh::operator=(edgeMesh(fromMaster));
        }
    }
    pointEdgesPtr_.clear();
}


// Trim edge mesh ouside bb
void Foam::edgeMesh::trim
(
    const List<treeBoundBox>& bbs
)
{
    Map<label> pointMap(points().size());
    pointField pts(points().size());
    edgeList edg(edges().size());
    label nPoints = 0;
    label nEdges = 0;

    forAll(edges(), edgeI)
    {
        const edge& e = edges()[edgeI];
        label v0 = e[0];
        label v1 = e[1];

        if (overlaps(bbs, points()[v0], points()[v1]))
        {
            edg[nEdges] = e;
            nEdges++;

            if (pointMap.insert(v0,nPoints))
            {
                pts[nPoints] = points()[v0];
                nPoints++;
            }
            if (pointMap.insert(v1,nPoints))
            {
                pts[nPoints] = points()[v1];
                nPoints++;
            }
        }
    }
    pts.setSize(nPoints);
    edg.setSize(nEdges);

    forAll(edg, edgeI)
    {
        edg[edgeI][0] = pointMap[edg[edgeI][0]];
        edg[edgeI][1] = pointMap[edg[edgeI][1]];
    }

    edges_.transfer(edg);
    points_.transfer(pts);

    pointEdgesPtr_.clear();
}

void Foam::edgeMesh::subsetEdgeMeshMap
(
    const edgeMesh& eMesh,
    const boolList& include,
    const label nIncluded,
    labelList& newToOldPoints,
    labelList& oldToNewPoints,
    labelList& newToOldEdges
)
{
    newToOldEdges.setSize(nIncluded);
    newToOldPoints.setSize(eMesh.points().size());
    oldToNewPoints.setSize(eMesh.points().size());
    oldToNewPoints = -1;
    {
        label edgeI = 0;
        label pointI = 0;

        forAll(include, oldEdgei)
        {
            if (include[oldEdgei])
            {
                // Store new edge compact
                newToOldEdges[edgeI++] = oldEdgei;

                // Renumber labels for edge
                const edge& e = eMesh.edges()[oldEdgei];

                label oldPointI = e[0];
                if (oldToNewPoints[oldPointI] == -1)
                {
                    oldToNewPoints[oldPointI] = pointI;
                    newToOldPoints[pointI++] = oldPointI;
                }
                oldPointI = e[1];
                if (oldToNewPoints[oldPointI] == -1)
                {
                    oldToNewPoints[oldPointI] = pointI;
                    newToOldPoints[pointI++] = oldPointI;
                }
            }
        }
        newToOldPoints.setSize(pointI);
    }
}


Foam::edgeMesh Foam::edgeMesh::subsetEdgeMesh
(
    const edgeMesh& eMesh,
    const labelList& newToOldPoints,
    const labelList& oldToNewPoints,
    const labelList& newToOldEdges
)
{
    // Extract points
    pointField newPoints(newToOldPoints.size());
    forAll(newToOldPoints, i)
    {
        newPoints[i] = eMesh.points()[newToOldPoints[i]];
    }
    // Extract edge
    List<edge> newEdges(newToOldEdges.size());

    forAll(newToOldEdges, i)
    {
        // Get old vertex labels
        const edge& e = eMesh.edges()[newToOldEdges[i]];

        newEdges[i][0] = oldToNewPoints[e[0]];
        newEdges[i][1] = oldToNewPoints[e[1]];
    }

    // Reuse storage.
    return edgeMesh(newPoints, newEdges);
}


Foam::edgeMesh Foam::edgeMesh::subsetEdgeMesh
(
    const edgeMesh& eMesh,
    const boolList& include,
    labelList& newToOldPoints,
    labelList& newToOldEdges
)
{
    label n = 0;

    forAll(include, i)
    {
        if (include[i])
        {
            n++;
        }
    }

    labelList oldToNewPoints;
    subsetEdgeMeshMap
    (
        eMesh,
        include,
        n,
        newToOldPoints,
        oldToNewPoints,
        newToOldEdges
    );

    return subsetEdgeMesh
    (
        eMesh,
        newToOldPoints,
        oldToNewPoints,
        newToOldEdges
    );
}

Foam::edgeMesh Foam::edgeMesh::subsetEdgeMesh
(
    const edgeMesh& eMesh,
    const labelList& newToOldEdges,
    labelList& newToOldPoints
)
{
    const boolList include
    (
        createWithValues<boolList>
        (
            eMesh.edges().size(),
            false,
            newToOldEdges,
            true
        )
    );

    newToOldPoints.setSize(eMesh.points().size());
    labelList oldToNewPoints(eMesh.points().size(), -1);
    {
        label pointI = 0;

        forAll(include, oldEdgei)
        {
            if (include[oldEdgei])
            {
                // Renumber labels for edge
                const edge& e = eMesh.edges()[oldEdgei];
                label oldPointI = e[0];

                if (oldToNewPoints[oldPointI] == -1)
                {
                    oldToNewPoints[oldPointI] = pointI;
                    newToOldPoints[pointI++] = oldPointI;
                }

                oldPointI = e[1];

                if (oldToNewPoints[oldPointI] == -1)
                {
                    oldToNewPoints[oldPointI] = pointI;
                    newToOldPoints[pointI++] = oldPointI;
                }
            }
        }
        newToOldPoints.setSize(pointI);
    }

    return subsetEdgeMesh
    (
        eMesh,
        newToOldPoints,
        oldToNewPoints,
        newToOldEdges
    );
}


// Subset the part of edgeMesh that is overlapping bb.
Foam::edgeMesh Foam::edgeMesh::overlappingEdgeMesh
(
    const edgeMesh& eMesh,
    const List<treeBoundBox>& bbs,
    labelList& subPointMap,
    labelList& subEdgeMap
)
{
    const edgeList edges = eMesh.edges();

    // Determine what triangles to keep.
    boolList includedEdges(edges.size(), false);

    // Create a slightly larger bounding box.
    List<treeBoundBox> bbsX(bbs.size());
    const scalar eps = 1.0e-4;
    forAll(bbs, i)
    {
        const point mid = 0.5*(bbs[i].min() + bbs[i].max());
        const vector halfSpan = (1.0+eps)*(bbs[i].max() - mid);

        bbsX[i].min() = mid - halfSpan;
        bbsX[i].max() = mid + halfSpan;
    }

    forAll(edges, edgeI)
    {
        const edge& e = edges[edgeI];
        const point& p0 = eMesh.points()[e[0]];
        const point& p1 = eMesh.points()[e[1]];

        if (overlaps(bbsX, p0, p1))
        {
            includedEdges[edgeI] = true;
        }
    }

    return subsetEdgeMesh(eMesh, includedEdges, subPointMap, subEdgeMap);
}


void Foam::edgeMesh::distribute
(
    const List<treeBoundBox>& bbs,
    const scalar& mergeDist
)
{
    // Get bbs of all domains
    // ~~~~~~~~~~~~~~~~~~~~~~

    List<List<treeBoundBox>> procBb(Pstream::nProcs());

    procBb[Pstream::myProcNo()].setSize(bbs.size());
    forAll(bbs, i)
    {
        procBb[Pstream::myProcNo()][i] = bbs[i];
    }
    Pstream::allGatherList(procBb);

    labelListList edgeSendMap(Pstream::nProcs());
    labelListList pointSendMap(Pstream::nProcs());

    forAll(procBb, procI)
    {
        overlappingEdgeMesh
        (
            *this,
            procBb[procI],
            pointSendMap[procI],
            edgeSendMap[procI]
        );
    }

    // Send over how many edges/points I need to receive
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListList edgeSendSizes(Pstream::nProcs());
    edgeSendSizes[Pstream::myProcNo()].setSize(Pstream::nProcs());
    forAll(edgeSendMap, procI)
    {
        edgeSendSizes[Pstream::myProcNo()][procI] = edgeSendMap[procI].size();
    }
    Pstream::allGatherList(edgeSendSizes);

    // Exchange edgeMesh
    // ~~~~~~~~~~~~~~~~~

    // Storage for resulting edgeMesh
    List<edge> allEdges;
    pointField allPoints;

    labelListList edgeConstructMap(Pstream::nProcs());
    labelListList pointConstructMap(Pstream::nProcs());

    // My own edgeMesh first
    // ~~~~~~~~~~~~~~~~~~~~~

    {
        labelList pointMap;
        edgeMesh subEdgeMesh
        (
            subsetEdgeMesh
            (
                *this,
                edgeSendMap[Pstream::myProcNo()],
                pointMap
            )
        );

        allEdges = subEdgeMesh.edges();
        allPoints = subEdgeMesh.points();

        edgeConstructMap[Pstream::myProcNo()] = identity
        (
            edgeSendMap[Pstream::myProcNo()].size()
        );
        pointConstructMap[Pstream::myProcNo()] = identity
        (
            pointSendMap[Pstream::myProcNo()].size()
        );
    }

    // Send all
    // ~~~~~~~~

    forAll(edgeSendSizes, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            if (edgeSendSizes[Pstream::myProcNo()][procI] > 0)
            {
                OPstream str(Pstream::commsTypes::blocking, procI);

                labelList pointMap;
                edgeMesh subEdgeMesh
                (
                    subsetEdgeMesh
                    (
                        *this,
                        edgeSendMap[procI],
                        pointMap
                    )
                );
                str << subEdgeMesh;
           }
        }
    }

    // Receive and merge all
    // ~~~~~~~~~~~~~~~~~~~~~

    forAll(edgeSendSizes, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            if (edgeSendSizes[procI][Pstream::myProcNo()] > 0)
            {
                IPstream str(Pstream::commsTypes::blocking, procI);
                // Receive
                edgeMesh subEdgeMesh(str);

                // Merge into allEdges
                merge
                (
                    mergeDist,
                    subEdgeMesh.edges(),
                    subEdgeMesh.points(),

                    allEdges,
                    allPoints,
                    edgeConstructMap[procI],
                    pointConstructMap[procI]
                );
           }
        }
    }

    // Construct edgeMesh. Reuse storage.
    edgeMesh::operator=(edgeMesh(allPoints, allEdges));

    pointEdgesPtr_.clear();

}


// ************************************************************************* //
