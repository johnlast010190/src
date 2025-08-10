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
    (c) 2016 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

Description
    region renumbering functions from renumberMesh

\*---------------------------------------------------------------------------*/

#include "renumberMethod/renumberMethod.H"


using namespace Foam;


// Determine face order such that inside region faces are sorted
// upper-triangular but inbetween region faces are handled like boundary faces.
labelList getRegionFaceOrder
(
    const primitiveMesh& mesh,
    const labelList& cellOrder,     // New to old cell
    const labelList& cellToRegion   // Old cell to region
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));

    labelList oldToNewFace(mesh.nFaces(), -1);

    label newFacei = 0;

    label prevRegion = -1;

    forAll(cellOrder, newCelli)
    {
        label oldCelli = cellOrder[newCelli];

        if (cellToRegion[oldCelli] != prevRegion)
        {
            prevRegion = cellToRegion[oldCelli];
        }

        const cell& cFaces = mesh.cells()[oldCelli];

        SortableList<label> nbr(cFaces.size());

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh.isInternalFace(facei))
            {
                // Internal face. Get cell on other side.
                label nbrCelli = reverseCellOrder[mesh.faceNeighbour()[facei]];
                if (nbrCelli == newCelli)
                {
                    nbrCelli = reverseCellOrder[mesh.faceOwner()[facei]];
                }

                if (cellToRegion[oldCelli] != cellToRegion[cellOrder[nbrCelli]])
                {
                    // Treat like external face. Do later.
                    nbr[i] = -1;
                }
                else if (newCelli < nbrCelli)
                {
                    // Celli is master
                    nbr[i] = nbrCelli;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        nbr.sort();

        forAll(nbr, i)
        {
            if (nbr[i] != -1)
            {
                oldToNewFace[cFaces[nbr.indices()[i]]] = newFacei++;
            }
        }
    }

    // Do region interfaces
    label nRegions = max(cellToRegion)+1;
    {
        // Sort in increasing region
        SortableList<label> sortKey(mesh.nFaces(), labelMax);

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            label ownRegion = cellToRegion[mesh.faceOwner()[facei]];
            label neiRegion = cellToRegion[mesh.faceNeighbour()[facei]];

            if (ownRegion != neiRegion)
            {
                sortKey[facei] =
                    min(ownRegion, neiRegion)*nRegions
                   +max(ownRegion, neiRegion);
            }
        }
        sortKey.sort();

        // Extract.
        label prevKey = -1;
        forAll(sortKey, i)
        {
            label key = sortKey[i];

            if (key == labelMax)
            {
                break;
            }

            if (prevKey != key)
            {
                prevKey = key;
            }

            oldToNewFace[sortKey.indices()[i]] = newFacei++;
        }
    }

    // Leave patch faces intact.
    for (label facei = newFacei; facei < mesh.nFaces(); facei++)
    {
        oldToNewFace[facei] = facei;
    }


    // Check done all faces.
    forAll(oldToNewFace, facei)
    {
        if (oldToNewFace[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position"
                << " for face " << facei
                << abort(FatalError);
        }
    }

    return invert(mesh.nFaces(), oldToNewFace);
}


// cellOrder: old cell for every new cell
// faceOrder: old face for every new face. Ordering of boundary faces not
//     changed.
autoPtr<mapPolyMesh> reorderMesh
(
    polyMesh& mesh,
    const labelList& cellOrder,
    const labelList& faceOrder
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));
    labelList reverseFaceOrder(invert(faceOrder.size(), faceOrder));

    faceList newFaces(reorder(reverseFaceOrder, mesh.faces()));
    labelList newOwner
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceOwner())
        )
    );
    labelList newNeighbour
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceNeighbour())
        )
    );

    // Check if any faces need swapping.
    labelHashSet flipFaceFlux(newOwner.size());
    forAll(newNeighbour, facei)
    {
        label own = newOwner[facei];
        label nei = newNeighbour[facei];

        if (nei < own)
        {
            newFaces[facei].flip();
            Swap(newOwner[facei], newNeighbour[facei]);
            flipFaceFlux.insert(facei);
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelList patchSizes(patches.size());
    labelList patchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    labelListList patchPointMap(patches.size());

    forAll(patches, patchi)
    {
        patchSizes[patchi] = patches[patchi].size();
        patchStarts[patchi] = patches[patchi].start();
        oldPatchNMeshPoints[patchi] = patches[patchi].nPoints();
        patchPointMap[patchi] = identity(patches[patchi].nPoints());
    }

    mesh.resetPrimitives
    (
        Xfer<pointField>::null(),
        xferMove(newFaces),
        xferMove(newOwner),
        xferMove(newNeighbour),
        patchSizes,
        patchStarts,
        true
    );


    // Re-do the faceZones
    {
        faceZoneMesh& faceZones = mesh.faceZones();
        faceZones.clearAddressing();
        forAll(faceZones, zoneI)
        {
            faceZone& fZone = faceZones[zoneI];
            labelList newAddressing(fZone.size());
            boolList newFlipMap(fZone.size());
            forAll(fZone, i)
            {
                label oldFacei = fZone[i];
                newAddressing[i] = reverseFaceOrder[oldFacei];
                if (flipFaceFlux.found(newAddressing[i]))
                {
                    newFlipMap[i] = !fZone.flipMap()[i];
                }
                else
                {
                    newFlipMap[i] = fZone.flipMap()[i];
                }
            }
            labelList newToOld;
            sortedOrder(newAddressing, newToOld);
            fZone.resetAddressing
            (
                UIndirectList<label>(newAddressing, newToOld)(),
                UIndirectList<bool>(newFlipMap, newToOld)()
            );
        }
    }
    // Re-do the cellZones
    {
        cellZoneMesh& cellZones = mesh.cellZones();
        cellZones.clearAddressing();
        forAll(cellZones, zoneI)
        {
            cellZones[zoneI] = UIndirectList<label>
            (
                reverseCellOrder,
                cellZones[zoneI]
            )();
            Foam::sort(cellZones[zoneI]);
        }
    }


    return autoPtr<mapPolyMesh>
    (
        new mapPolyMesh
        (
            mesh,                       // const polyMesh& mesh,
            mesh.nPoints(),             // nOldPoints,
            mesh.nFaces(),              // nOldFaces,
            mesh.nCells(),              // nOldCells,
            identity(mesh.nPoints()),   // pointMap,
            List<objectMap>(0),         // pointsFromPoints,
            faceOrder,                  // faceMap,
            List<objectMap>(0),         // facesFromPoints,
            List<objectMap>(0),         // facesFromEdges,
            List<objectMap>(0),         // facesFromFaces,
            cellOrder,                  // cellMap,
            List<objectMap>(0),         // cellsFromPoints,
            List<objectMap>(0),         // cellsFromEdges,
            List<objectMap>(0),         // cellsFromFaces,
            List<objectMap>(0),         // cellsFromCells,
            identity(mesh.nPoints()),   // reversePointMap,
            reverseFaceOrder,           // reverseFaceMap,
            reverseCellOrder,           // reverseCellMap,
            flipFaceFlux,               // flipFaceFlux,
            patchPointMap,              // patchPointMap,
            labelListList(0),           // pointZoneMap,
            labelListList(0),           // faceZonePointMap,
            labelListList(0),           // faceZoneFaceMap,
            labelListList(0),           // cellZoneMap,
            pointField(0),              // preMotionPoints,
            patchStarts,                // oldPatchStarts,
            oldPatchNMeshPoints,        // oldPatchNMeshPoints
            autoPtr<scalarField>()      // oldCellVolumes
        )
    );
}


// Return new to old cell numbering
labelList regionRenumber
(
    const renumberMethod& method,
    const fvMesh& mesh,
    const labelList& cellToRegion
)
{
    Info<< "Determining cell order:" << endl;

    labelList cellOrder(cellToRegion.size());

    label nRegions = max(cellToRegion)+1;

    labelListList regionToCells(invertOneToMany(nRegions, cellToRegion));

    label celli = 0;

    forAll(regionToCells, regionI)
    {
        //Info<< "    region " << regionI << " starts at " << celli << endl;

        // Make sure no parallel comms
        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        // Per region do a reordering.
        fvMeshSubset subsetter(mesh);
        subsetter.setLargeCellSubset(cellToRegion, regionI);

        const fvMesh& subMesh = subsetter.subMesh();

        labelList subCellOrder = method.renumber
        (
            subMesh,
            subMesh.cellCentres()
        );

        // Restore state
        UPstream::parRun() = oldParRun;

        const labelList& cellMap = subsetter.cellMap();

        forAll(subCellOrder, i)
        {
            cellOrder[celli++] = cellMap[subCellOrder[i]];
        }
    }
    Info<< endl;

    return cellOrder;
}
