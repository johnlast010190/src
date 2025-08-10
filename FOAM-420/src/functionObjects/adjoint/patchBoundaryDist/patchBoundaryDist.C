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

\*---------------------------------------------------------------------------*/

#include "patchBoundaryDist/patchBoundaryDist.H"
#include "fvMesh/fvMesh.H"
#include "memory/autoPtr/autoPtr.H"
#include "indexedOctree/treeDataEdge.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "algorithms/indexedOctree/indexedOctree.H"

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



void patchBoundaryDist::calculate()
{
//    autoPtr<indexedOctree<treeDataEdge>> tree(buildEdgeTree());
    //create index of edges

    labelList edgeOnPatch(mesh_.nEdges(), 0);

    const fvPatchList& patches = mesh_.boundary();
    const fvPatch& patch = patches[index_];
    const labelList& pMeshEdges = patch.patch().meshEdges();
    const edgeList& edges = mesh_.edges();
    const pointField& points = mesh_.points();

    forAll(pMeshEdges, leI)
    {
        edgeOnPatch[pMeshEdges[leI]] = 1;
    }

    // parallel sync of all edges belonging to target patch
    syncTools::syncEdgeList
    (
        mesh_,
        edgeOnPatch,
        maxEqOp<label>(),
        label(0)      // initial value
    );

    //now loop over edges from all non-constrained boundaries
    //(i.e. walls, inlets, outlets, patches)
    //contained bc's = empty, symmetryplane, processor, wedge, cyclic

    // list of collected edges
    edgeList patchBoundaryEdges(pMeshEdges.size());

    // new list of points - need a seperate list since they will be needed
    // on all processors
    // We will not attempt to rationalise point list
    // Each edge will have its own associated points, thus no sharing
    pointField patchBoundaryPoints(2*pMeshEdges.size());

    label bEdgeI = 0;

    forAll(patches, pI)
    {
        const fvPatch& cpatch = patches[pI];
        if (!cpatch.constrained() && pI != index_)
        {
            const labelList& cMeshEdges = cpatch.patch().meshEdges();
            forAll(cMeshEdges, leI)
            {
                if (edgeOnPatch[cMeshEdges[leI]] == 1)
                {
                    patchBoundaryEdges[bEdgeI] = edge(2*bEdgeI, 2*bEdgeI+1);
                    patchBoundaryPoints[2*bEdgeI] =
                    (
                        points[edges[cMeshEdges[leI]].start()]
                    );
                    patchBoundaryPoints[2*bEdgeI + 1] =
                    (
                        points[edges[cMeshEdges[leI]].end()]
                    );

                    bEdgeI++;
                }
            }
        }
    }

    patchBoundaryEdges.setSize(bEdgeI);
    patchBoundaryPoints.setSize(2*bEdgeI);
    reduce(bEdgeI, sumOp<label>());

    // do a list gather/scatter so all procs have a complete set of data
    if (Pstream::parRun())
    {
        List<List<edge>> gPatchEdges(Pstream::nProcs());
        List<pointField > gPatchPoints(Pstream::nProcs());

        gPatchEdges[Pstream::myProcNo()] = patchBoundaryEdges;
        gPatchPoints[Pstream::myProcNo()] = patchBoundaryPoints;

        Pstream::allGatherList(gPatchEdges);
        Pstream::allGatherList(gPatchPoints);

        //rebuild unified list from processor pieces
        patchBoundaryEdges.setSize(bEdgeI);
        patchBoundaryPoints.setSize(2*bEdgeI);
        bEdgeI = 0;

        forAll(gPatchEdges, procI)
        {
            forAll(gPatchEdges[procI], peI)
            {
                patchBoundaryEdges[bEdgeI] = edge(2*bEdgeI, 2*bEdgeI+1);
                patchBoundaryPoints[2*bEdgeI] = gPatchPoints[procI][2*peI];
                patchBoundaryPoints[2*bEdgeI+1] = gPatchPoints[procI][2*peI+1];
                bEdgeI++;
            }
        }
    }

    //identity addressing required for octree constructor
    labelList edgeIndices(identity(patchBoundaryEdges.size()));

    // need to compute bounding box boundary patch edges
    treeBoundBox bb(patchBoundaryPoints);

    //ensure there are no zero dimensions because the tree construction
    //can't handle this
    scalar smallNumber = 0.001;
    scalar mySmall = bb.avgDim() < 1.0 ? smallNumber : bb.avgDim()*smallNumber;
    //mySmall *= 1000.0;
    for (label i = 0; i < 3; ++i)
    {
        if (bb.min()[i] + mySmall >= bb.max()[i])
        {
            bb.min()[i] -= mySmall;
            bb.max()[i] += mySmall;
        }
    }

    //build octree search engine
    indexedOctree<treeDataEdge> tree
    (
        treeDataEdge
        (
            false,                     // cachebb
            patchBoundaryEdges,                     //list of edges
            patchBoundaryPoints,                    // points
            edgeIndices                    // selected edges
        ),
        bb,
        10,    // maxLevel
        1,     // leafsize
        3.0
    );

    //calculate distances from patch face centres to
    //un-constrained patch bondary edges

    forAll(patch, fI)
    {
        pointIndexHit hit
        = tree.findNearest(patch.Cf()[fI],GREAT);

        this->operator[](fI) = mag(hit.hitPoint()-patch.Cf()[fI]);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

patchBoundaryDist::patchBoundaryDist
(
    const fvMesh& mesh,
    label pI
)
:
    scalarField(mesh.boundary()[pI].size(), 1.0),
    mesh_(mesh),
    index_(pI)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

patchBoundaryDist::~patchBoundaryDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void patchBoundaryDist::correct()
{
    calculate();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
