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
    This function calculates the list of patch edges, defined on the list of
    points supporting the patch. The edges are ordered:
    - 0..nInternalEdges-1 : upper triangular order
    - nInternalEdges..    : boundary edges (no particular order)

    Other patch addressing information is also calculated:
    - faceFaces with neighbour faces in ascending order
    - edgeFaces with ascending face order
    - faceEdges sorted according to edges of a face

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/PrimitivePatch/PrimitivePatch.H"
#include "containers/Lists/DynamicList/DynamicList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void Foam::PrimitivePatch<FaceList, PointField>::calcAddressing() const
{
    if (debug)
    {
        InfoInFunction << "Calculating patch addressing" << endl;
    }

    if (edgesPtr_ || faceFacesPtr_ || edgeFacesPtr_ || faceEdgesPtr_)
    {
        // it is considered an error to attempt to recalculate
        // if already allocated
        FatalErrorInFunction
            << "addressing already calculated"
            << abort(FatalError);
    }

    // get reference to localFaces
    const List<FaceType>& locFcs = localFaces();

    // get reference to pointFaces
    const labelListList& pf = pointFaces();

    // Guess the max number of edges and neighbours for a face
    label maxEdges = 0;
    forAll(locFcs, facei)
    {
        maxEdges += locFcs[facei].size();
    }

    // create the lists for the various results. (resized on completion)
    edgesPtr_ = new edgeList(maxEdges);
    edgeList& edges = *edgesPtr_;

    edgeFacesPtr_ = new labelListList(maxEdges);
    labelListList& edgeFaces = *edgeFacesPtr_;

    // faceFaces created using a dynamic list.  Cannot guess size because
    // of multiple connections
    List<DynamicList<label>> ff(locFcs.size());

    faceEdgesPtr_ = new labelListList(locFcs.size());
    labelListList& faceEdges = *faceEdgesPtr_;

    forAll(locFcs, facei)
    {
        labelList& curFaceEdges = faceEdges[facei];
        curFaceEdges.setSize(locFcs[facei].nEdges());

        forAll(curFaceEdges, faceEdgeI)
        {
            curFaceEdges[faceEdgeI] = -1;
        }
    }

    // This algorithm will produce a separated list of edges, internal edges
    // starting from 0 and boundary edges starting from the top and
    // growing down.

    label nEdges = 0;

    bool found = false;

    // For all local faces ...
    forAll(locFcs, facei)
    {
        // Get reference to vertices of current face and corresponding edges.
        const FaceType& curF = locFcs[facei];

        // Record the neighbour face.  Multiple connectivity allowed
        List<DynamicList<label>> neiFaces(curF.size());
        List<DynamicList<label>> edgeOfNeiFace(curF.size());

        label nNeighbours = 0;

        // For all edges ...
        int edgeI = 0;
        for (const edge& e : locFcs[facei].walkEdges())
        {
            // If the edge is already detected, skip
            if (faceEdges[facei][edgeI] >= 0) {
                edgeI++;
                continue;
            }

            found = false;

            // Collect neighbours for the current face vertex.

            const labelList& nbrFaces = pf[e.start()];

            forAll(nbrFaces, nbrFacei)
            {
                // set reference to the current neighbour
                label curNei = nbrFaces[nbrFacei];

                // Reject neighbours with the lower label
                if (curNei > facei)
                {
                    int neiEdgeI = 0;
                    for (const edge& searchEdge : locFcs[curNei].walkEdges())
                    {
                        if (searchEdge == e)
                        {
                            // Match
                            found = true;

                            neiFaces[edgeI].append(curNei);
                            edgeOfNeiFace[edgeI].append(neiEdgeI);

                            // Record faceFaces both ways
                            ff[facei].append(curNei);
                            ff[curNei].append(facei);

                            // Keep searching due to multiple connectivity
                        }

                        neiEdgeI++;
                    }
                }
            } // End of neighbouring faces

            if (found)
            {
                // Register another detected internal edge
                nNeighbours++;
            }

            edgeI++;
        } // End of current edges

        // Add the edges in increasing number of neighbours.
        // Note: for multiply connected surfaces, the lower index neighbour for
        // an edge will come first.

        // Add the faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = locFcs.size();

            forAll(neiFaces, nfI)
            {
                if (neiFaces[nfI].size() && neiFaces[nfI][0] < minNei)
                {
                    nextNei = nfI;
                    minNei = neiFaces[nfI][0];
                }
            }

            if (nextNei > -1)
            {
                // Add the face to the list of faces
                edges[nEdges] = locFcs[facei].faceEdge(nextNei);

                // Set face-edge and face-neighbour-edge to current face label
                faceEdges[facei][nextNei] = nEdges;

                DynamicList<label>& cnf = neiFaces[nextNei];
                DynamicList<label>& eonf = edgeOfNeiFace[nextNei];

                // Set edge-face addressing
                labelList& curEf = edgeFaces[nEdges];
                curEf.setSize(cnf.size() + 1);
                curEf[0] = facei;

                forAll(cnf, cnfI)
                {
                    faceEdges[cnf[cnfI]][eonf[cnfI]] = nEdges;

                    curEf[cnfI + 1] = cnf[cnfI];
                }

                // Stop the neighbour from being used again
                cnf.clear();
                eonf.clear();

                // Increment number of faces counter
                nEdges++;
            }
            else
            {
                FatalErrorInFunction
                    << "Error in internal edge insertion"
                    << abort(FatalError);
            }
        }
    }

    nInternalEdges_ = nEdges;

    // Do boundary faces

    forAll(faceEdges, facei)
    {
        labelList& curEdges = faceEdges[facei];

        forAll(curEdges, edgeI)
        {
            if (curEdges[edgeI] < 0)
            {
                // Grab edge and faceEdge
                edges[nEdges] = locFcs[facei].faceEdge(edgeI);
                curEdges[edgeI] = nEdges;

                // Add edgeFace
                labelList& curEf = edgeFaces[nEdges];
                curEf.setSize(1);
                curEf[0] = facei;

                nEdges++;
            }
        }
    }

    // edges
    edges.setSize(nEdges);

    // edgeFaces list
    edgeFaces.setSize(nEdges);

    // faceFaces list
    faceFacesPtr_ = new labelListList(locFcs.size());
    labelListList& faceFaces = *faceFacesPtr_;

    forAll(faceFaces, facei)
    {
        faceFaces[facei].transfer(ff[facei]);
    }

    if (debug)
    {
        InfoInFunction << "Finished calculating patch addressing" << endl;
    }
}


// ************************************************************************* //
