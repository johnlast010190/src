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

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::labelList Foam::PrimitivePatch<FaceList, PointField>::meshEdges
(
    const edgeList& allEdges,
    const labelListList& cellEdges,
    const labelList& faceCells
) const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating labels of patch edges in mesh edge list" << endl;
    }

    // get reference to the list of edges on the patch
    const edgeList& PatchEdges = edges();

    const labelListList& EdgeFaces = edgeFaces();

    // create the storage
    labelList meshEdges(PatchEdges.size());

    bool found = false;

    // get reference to the points on the patch
    const labelList& pp = meshPoints();

    // WARNING: Remember that local edges address into local point list;
    // local-to-global point label translation is necessary
    forAll(PatchEdges, edgeI)
    {
        const edge curEdge
            (pp[PatchEdges[edgeI].start()], pp[PatchEdges[edgeI].end()]);

        found = false;

        // get the patch faces sharing the edge
        const labelList& curFaces = EdgeFaces[edgeI];

        forAll(curFaces, facei)
        {
            // get the cell next to the face
            label curCell = faceCells[curFaces[facei]];

            // get reference to edges on the cell
            const labelList& ce = cellEdges[curCell];

            forAll(ce, cellEdgeI)
            {
                if (allEdges[ce[cellEdgeI]] == curEdge)
                {
                    found = true;

                    meshEdges[edgeI] = ce[cellEdgeI];

                    break;
                }
            }

            if (found) break;
        }
    }

    return meshEdges;
}


template<class FaceList, class PointField>
Foam::labelList Foam::PrimitivePatch<FaceList, PointField>::meshEdges
(
    const edgeList& allEdges,
    const labelListList& pointEdges
) const
{
    if (debug)
    {
        InfoInFunction
            << "Calculating labels of patch edges in mesh edge list" << endl;
    }

    // get reference to the list of edges on the patch
    const edgeList& PatchEdges = edges();

    // create the storage
    labelList meshEdges(PatchEdges.size());

    // get reference to the points on the patch
    const labelList& pp = meshPoints();

    // WARNING: Remember that local edges address into local point list;
    // local-to-global point label translation is necessary
    forAll(PatchEdges, edgeI)
    {
        const label globalPointi = pp[PatchEdges[edgeI].start()];
        const edge curEdge(globalPointi, pp[PatchEdges[edgeI].end()]);

        const labelList& pe = pointEdges[globalPointi];

        forAll(pe, i)
        {
            if (allEdges[pe[i]] == curEdge)
            {
                meshEdges[edgeI] = pe[i];
                break;
            }
        }
    }

    return meshEdges;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::label Foam::PrimitivePatch<FaceList, PointField>::whichEdge
(
    const edge& e
) const
{
    // Get pointEdges from the starting point and search all the candidates
    const edgeList& Edges = edges();

    if (e.start() > -1 && e.start() < nPoints())
    {
        const labelList& pe = pointEdges()[e.start()];

        forAll(pe, peI)
        {
            if (e == Edges[pe[peI]])
            {
                return pe[peI];
            }
        }
    }

    // Edge not found.  Return -1
    return -1;
}


// ************************************************************************* //
