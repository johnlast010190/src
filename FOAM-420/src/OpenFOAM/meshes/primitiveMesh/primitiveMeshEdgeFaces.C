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

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/primitiveMesh.H"
#include "containers/Lists/ListOps/ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::primitiveMesh::calcEdgeFaces() const
{
    FOAM_ASSERT(!efPtr_)
    {
        FatalErrorInFunction << "Should not re-calculate edge faces." << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "primitiveMesh::edgeFaces() : calculating edgeFaces" << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }

    // Invert faceEdges
    efPtr_ = new labelListList(nEdges());
    invertManyToMany(nEdges(), faceEdges(), *efPtr_);
}


const Foam::labelList& Foam::primitiveMesh::edgeFaces
(
    const label edgeI,
    DynamicList<label>& storage
) const
{
    if (hasEdgeFaces())
    {
        return edgeFaces()[edgeI];
    }
    else
    {
        // Use the fact that pointEdges are sorted in incrementing edge order
        const edge& e = edges()[edgeI];
        const labelList& pFaces0 = pointFaces()[e[0]];
        const labelList& pFaces1 = pointFaces()[e[1]];

        label i0 = 0;
        label i1 = 0;

        storage.clear();

        while (i0 < pFaces0.size() && i1 < pFaces1.size())
        {
            if (pFaces0[i0] < pFaces1[i1])
            {
                ++i0;
            }
            else if (pFaces0[i0] > pFaces1[i1])
            {
                ++i1;
            }
            else
            {
                //check face does contain edge (non-manifold cell)
                face f = faces()[pFaces0[i0]];
                label i3 = findIndex(f, e[0]);
                if (f[f.fcIndex(i3)] == e[1] || f[f.rcIndex(i3)] == e[1])
                {
                    // Equal. Append.
                    storage.append(pFaces0[i0]);
                }
                ++i0;
                ++i1;
            }
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::edgeFaces(const label edgeI) const
{
    return edgeFaces(edgeI, labels_);
}


// ************************************************************************* //
