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
    (c) 2011-2015 OpenFOAM Foundation

Description
    Contains fix for PrimitivePatch addressing (which doesn't work if surface
    is non-manifold). Should be moved into PrimitivePatch.

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurface.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "primitives/transform/transform.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::triSurface::calcSortedEdgeFaces() const
{
    if (sortedEdgeFacesPtr_)
    {
        FatalErrorInFunction
            << "sortedEdgeFacesPtr_ already set"
            << abort(FatalError);
    }

    const labelListList& eFaces = edgeFaces();

    sortedEdgeFacesPtr_ = new labelListList(eFaces.size());
    labelListList& sortedEdgeFaces = *sortedEdgeFacesPtr_;

    sortedEdgeFaces = PatchTools::sortedEdgeFaces(*this);
}


void Foam::triSurface::calcEdgeOwner() const
{
    if (edgeOwnerPtr_)
    {
        FatalErrorInFunction
            << "edgeOwnerPtr_ already set"
            << abort(FatalError);
    }

    // create the owner list
    edgeOwnerPtr_ = new labelList(nEdges());
    labelList& edgeOwner = *edgeOwnerPtr_;

    forAll(edges(), edgeI)
    {
        const edge& e = edges()[edgeI];

        const labelList& myFaces = edgeFaces()[edgeI];

        if (myFaces.size() == 1)
        {
            edgeOwner[edgeI] = myFaces[0];
        }
        else
        {
            // Find the first face whose vertices are aligned with the edge.
            // (in case of multiply connected edge the best we can do)
            edgeOwner[edgeI] = -1;

            forAll(myFaces, i)
            {
                const labelledTri& f = localFaces()[myFaces[i]];

                if
                (
                    ((f[0] == e.start()) && (f[1] == e.end()))
                 || ((f[1] == e.start()) && (f[2] == e.end()))
                 || ((f[2] == e.start()) && (f[0] == e.end()))
                )
                {
                    edgeOwner[edgeI] = myFaces[i];

                    break;
                }
            }

            if (edgeOwner[edgeI] == -1)
            {
                FatalErrorInFunction
                    << "Edge " << edgeI << " vertices:" << e
                    << " is used by faces " << myFaces
                    << " vertices:"
                    << UIndirectList<labelledTri>(localFaces(), myFaces)()
                    << " none of which use the edge vertices in the same order"
                    << nl << "I give up" << abort(FatalError);
            }
        }
    }
}


// ************************************************************************* //
