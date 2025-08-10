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
    (c) 2017-2020 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::polyTopoChange::reorder
(
    const labelList& oldToNew,
    DynamicList<T>& lst
)
{
    // Create a new allocation, full of garbage (the only heap allocation this function performs).
    DynamicList<T> newList;
    newList.setSize(lst.size());

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        // Move every element into the new list at the proper position.
        if (newElemI != -1) {
            newList[newElemI] = std::move(lst[elemI]);
        }
    }

    // Move the new list on top of the output.
    lst = std::move(newList);
}


template<class T>
void Foam::polyTopoChange::reorder
(
    const labelList& oldToNew,
    List<DynamicList<T>>& lst
)
{
    // Create copy
    List<DynamicList<T>> oldLst(lst);

    forAll(oldToNew, elemI)
    {
        label newElemI = oldToNew[elemI];

        if (newElemI != -1)
        {
            lst[newElemI].transfer(oldLst[elemI]);
        }
    }
}


template<class T>
void Foam::polyTopoChange::renumberKey
(
    const labelList& oldToNew,
    Map<T>& elems
)
{
    Map<T> newElems(elems.size());

    forAllConstIter(typename Map<T>, elems, iter)
    {
        label newElem = oldToNew[iter.key()];

        if (newElem >= 0)
        {
            newElems.insert(newElem, iter());
        }
    }

    elems.transfer(newElems);
}


template<class Type>
Foam::autoPtr<Foam::mapPolyMesh> Foam::polyTopoChange::makeMesh
(
    autoPtr<Type>& newMeshPtr,
    const IOobject& io,
    const polyMesh& mesh,
    const bool syncParallel,
    const bool orderCells,
    const bool orderPoints
)
{
    if (debug)
    {
        Pout<< "polyTopoChange::changeMesh"
            << "(autoPtr<fvMesh>&, const IOobject&, const fvMesh&"
            << ", const bool, const bool, const bool)"
            << endl;
    }

    if (debug)
    {
        Pout<< "Old mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    // new mesh points
    pointField newPoints;
    // number of internal points
    label nInternalPoints;
    // patch slicing
    labelList patchSizes;
    labelList patchStarts;
    // inflate maps
    List<objectMap> pointsFromPoints;
    List<objectMap> facesFromPoints;
    List<objectMap> facesFromEdges;
    List<objectMap> facesFromFaces;
    List<objectMap> cellsFromPoints;
    List<objectMap> cellsFromEdges;
    List<objectMap> cellsFromFaces;
    List<objectMap> cellsFromCells;

    // old mesh info
    List<Map<label>> oldPatchMeshPointMaps;
    labelList oldPatchNMeshPoints;
    labelList oldPatchStarts;
    List<Map<label>> oldFaceZoneMeshPointMaps;

    // Compact, reorder patch faces and calculate mesh/patch maps.
    compactAndReorder
    (
        mesh,
        syncParallel,
        orderCells,
        orderPoints,

        nInternalPoints,
        newPoints,
        patchSizes,
        patchStarts,
        pointsFromPoints,
        facesFromPoints,
        facesFromEdges,
        facesFromFaces,
        cellsFromPoints,
        cellsFromEdges,
        cellsFromFaces,
        cellsFromCells,
        oldPatchMeshPointMaps,
        oldPatchNMeshPoints,
        oldPatchStarts,
        oldFaceZoneMeshPointMaps
    );

    const label nOldPoints(mesh.nPoints());
    const label nOldFaces(mesh.nFaces());
    const label nOldCells(mesh.nCells());
    autoPtr<scalarField> oldCellVolumes(new scalarField(mesh.cellVolumes()));


    // Create the mesh
    // ~~~~~~~~~~~~~~~

    IOobject noReadIO(io);
    noReadIO.readOpt() = IOobject::NO_READ;
    newMeshPtr.reset
    (
        new Type
        (
            noReadIO,
            xferMove(newPoints),
            faces_.xfer(),
            faceOwner_.xfer(),
            faceNeighbour_.xfer()
        )
    );
    Type& newMesh = newMeshPtr();

    // Clear out primitives
    {
        retiredPoints_.clearStorage();
        region_.clearStorage();
    }


    if (debug)
    {
        // Some stats on changes
        label nAdd, nInflate, nMerge, nRemove;
        countMap(pointMap_, reversePointMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Points:"
            << "  added(from point):" << nAdd
            << "  added(from nothing):" << nInflate
            << "  merged(into other point):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(faceMap_, reverseFaceMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Faces:"
            << "  added(from face):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other face):" << nMerge
            << "  removed:" << nRemove
            << nl;

        countMap(cellMap_, reverseCellMap_, nAdd, nInflate, nMerge, nRemove);
        Pout<< "Cells:"
            << "  added(from cell):" << nAdd
            << "  added(inflated):" << nInflate
            << "  merged(into other cell):" << nMerge
            << "  removed:" << nRemove
            << nl
            << endl;
    }


    {
        const polyBoundaryMesh& oldPatches = mesh.boundaryMesh();

        List<polyPatch*> newBoundary(oldPatches.size());

        forAll(oldPatches, patchi)
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(oldPatches[patchi]);
            newBoundary[patchi] = dpp.clone
            (
                newMesh.boundaryMesh(),
                patchi,
                patchSizes[patchi],
                patchStarts[patchi]
            ).ptr();
        }
        newMesh.addFvPatches(newBoundary);
    }


    // Zones
    // ~~~~~

    // Start off from empty zones.
    const pointZoneMesh& oldPointZones = mesh.pointZones();
    List<pointZone*> pZonePtrs(oldPointZones.size());
    {
        forAll(oldPointZones, i)
        {
            pZonePtrs[i] = new pointZone
            (
                oldPointZones[i].name(),
                labelList(0),
                i,
                newMesh.pointZones()
            );
        }
    }

    const faceZoneMesh& oldFaceZones = mesh.faceZones();
    List<faceZone*> fZonePtrs(oldFaceZones.size());
    {
        forAll(oldFaceZones, i)
        {
            fZonePtrs[i] = new faceZone
            (
                oldFaceZones[i].name(),
                labelList(0),
                boolList(0),
                i,
                newMesh.faceZones()
            );
        }
    }

    const cellZoneMesh& oldCellZones = mesh.cellZones();
    List<cellZone*> cZonePtrs(oldCellZones.size());
    {
        forAll(oldCellZones, i)
        {
            cZonePtrs[i] = new cellZone
            (
                oldCellZones[i].name(),
                labelList(0),
                i,
                newMesh.cellZones()
            );
        }
    }

    newMesh.addZones(pZonePtrs, fZonePtrs, cZonePtrs);

    // Inverse of point/face/cell zone addressing.
    // For every preserved point/face/cells in zone give the old position.
    // For added points, the index is set to -1
    labelListList pointZoneMap(mesh.pointZones().size());
    labelListList faceZoneFaceMap(mesh.faceZones().size());
    labelListList cellZoneMap(mesh.cellZones().size());

    resetZones(mesh, newMesh, pointZoneMap, faceZoneFaceMap, cellZoneMap);

    // Clear zone info
    {
        pointZone_.clearStorage();
        faceZone_.clearStorage();
        faceZoneFlip_.clearStorage();
        cellZone_.clearStorage();
    }

    // Patch point renumbering
    // For every preserved point on a patch give the old position.
    // For added points, the index is set to -1
    labelListList patchPointMap(newMesh.boundaryMesh().size());
    calcPatchPointMap
    (
        oldPatchMeshPointMaps,
        newMesh.boundaryMesh(),
        patchPointMap
    );

    // Create the face zone mesh point renumbering
    labelListList faceZonePointMap(newMesh.faceZones().size());
    calcFaceZonePointMap(newMesh, oldFaceZoneMeshPointMaps, faceZonePointMap);

    if (debug)
    {
        Pout<< "New mesh:" << nl;
        writeMeshStats(mesh, Pout);
    }

    labelHashSet flipFaceFluxSet(getSetIndices(flipFaceFlux_));

    return autoPtr<mapPolyMesh>
    (
        new mapPolyMesh
        (
            newMesh,
            nOldPoints,
            nOldFaces,
            nOldCells,

            pointMap_,
            pointsFromPoints,

            faceMap_,
            facesFromPoints,
            facesFromEdges,
            facesFromFaces,

            cellMap_,
            cellsFromPoints,
            cellsFromEdges,
            cellsFromFaces,
            cellsFromCells,

            reversePointMap_,
            reverseFaceMap_,
            reverseCellMap_,

            flipFaceFluxSet,

            patchPointMap,

            pointZoneMap,

            faceZonePointMap,
            faceZoneFaceMap,
            cellZoneMap,

            newPoints,          // if empty signals no inflation.
            oldPatchStarts,
            oldPatchNMeshPoints,
            oldCellVolumes,
            true                // steal storage.
        )
    );

    // At this point all member DynamicList (pointMap_, cellMap_ etc.) will
    // be invalid.
}


// ************************************************************************* //
