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

\*---------------------------------------------------------------------------*/

#include "meshes/primitiveMesh/primitiveMesh.H"
#include "containers/Lists/ListOps/ListOps.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::primitiveMesh::calcEdgeCells() const
{
    FOAM_ASSERT(!ecPtr_)
    {
        FatalErrorInFunction << "Should not re-calculate edge faces." << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "primitiveMesh::edgeCells() : calculating edgeCells" << endl;

        if (debug == -1)
        {
            // For checking calls:abort so we can quickly hunt down
            // origin of call
            FatalErrorInFunction
                << abort(FatalError);
        }
    }
    // Invert cellEdges
    ecPtr_ = new labelListList(nEdges());
    invertManyToMany(nEdges(), cellEdges(), *ecPtr_);
}


const Foam::labelList& Foam::primitiveMesh::edgeCells
(
    const label edgeI,
    DynamicList<label>& storage
) const
{
    if (hasEdgeCells())
    {
        return edgeCells()[edgeI];
    }
    else
    {
        const labelList& own = faceOwner();
        const labelList& nei = faceNeighbour();

        // Construct edgeFaces
        DynamicList<label> eFacesStorage;
        const labelList& eFaces = edgeFaces(edgeI, eFacesStorage);

        storage.clear();

        // Do quadratic insertion.
        forAll(eFaces, i)
        {
            label facei = eFaces[i];

            {
                label ownCelli = own[facei];

                // Check if not already in storage
                forAll(storage, j)
                {
                    if (storage[j] == ownCelli)
                    {
                        ownCelli = -1;
                        break;
                    }
                }

                if (ownCelli != -1)
                {
                    storage.append(ownCelli);
                }
            }

            if (isInternalFace(facei))
            {
                label neiCelli = nei[facei];

                forAll(storage, j)
                {
                    if (storage[j] == neiCelli)
                    {
                        neiCelli = -1;
                        break;
                    }
                }

                if (neiCelli != -1)
                {
                    storage.append(neiCelli);
                }
            }
        }

        return storage;
    }
}


const Foam::labelList& Foam::primitiveMesh::edgeCells(const label edgeI) const
{
    return edgeCells(edgeI, labels_);
}


// ************************************************************************* //
