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

\*---------------------------------------------------------------------------*/

#include "mergePolyMesh.H"
#include "db/Time/Time.H"
#include "polyTopoChange/polyTopoChanger/polyTopoChanger.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(mergePolyMesh, 1);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mergePolyMesh::setPatchIndices
(
    const polyBoundaryMesh& bm,
    labelList& patchIndices,
    DynamicList<labelPair>& movedPatches
)
{
    const polyBoundaryMesh& oldPatches = boundaryMesh();

    DynamicList<word> newPatchNames(2*patchNames_.size());
    DynamicList<dictionary> newPatchDicts(2*patchDicts_.size());

    forAll(patchNames_, patchi)
    {
        if (word(patchDicts_[patchi]["type"]) != "processor")
        {
            newPatchNames.append(patchNames_[patchi]);
            newPatchDicts.append(patchDicts_[patchi]);
        }
    }

    boolList renamedBC(bm.size(), false);
    //Add new non-processor faces
    forAll(bm, i)
    {
        const polyPatch& p = bm[i];
        const word& pType = p.type();
        const word& pName = p.name();

        if (pType != "processor")
        {
            label oldPatchID = oldPatches.findPatchID(pName);
            if (oldPatchID != -1)
            {
                if (word(patchDicts_[oldPatchID]["type"]) == pType)
                {
                    patchIndices[i] = oldPatchID;
                }
                else
                {
                    {
                        OStringStream os;
                        p.write(os);
                        newPatchDicts.append(dictionary(IStringStream(os.str())()));
                    }

                    // Duplicate name is not allowed.  Create a composite name from the
                    // patch name and case name
                    const word& caseName = p.boundaryMesh().mesh().time().caseName();
                    word newPatchName(pName + "_" + caseName);
                    if
                    (
                        oldPatches.findPatchID(newPatchName) != -1
                        || bm.findPatchID(newPatchName) != -1
                    )
                    {
                        newPatchName = word(pName + "__" + caseName);
                    }
                    newPatchNames.append(newPatchName);

                    Info<< "label patchIndex(const polyPatch& p) : "
                        << "Patch " << p.index() << " named "
                        << pName << " in mesh " << caseName
                        << " already exists, but patch types "
                        << " do not match.\nCreating a composite name as "
                        << newPatchNames.last() << endl;

                    patchIndices[i] = newPatchNames.size() - 1;
                    renamedBC[i] = true;
                }
            }
            else
            {
                // Patch not found.  Append to the list
                {
                    OStringStream os;
                    p.write(os);
                    newPatchDicts.append(dictionary(IStringStream(os.str())()));
                }
                newPatchNames.append(pName);
                patchIndices[i] = newPatchNames.size() - 1;
            }
        }
    }

    //Rename cyclicAMI neigbour references if re-named
    forAll(bm, i)
    {
        const polyPatch& p = bm[i];
        const word& pType = p.type();
        if
        (
            pType == "cyclic"
            || pType == "cyclicAMI"
            || pType == "cyclicSlip"
            || pType == "cyclicPeriodicAMI"
        )
        {
            label newPatchID = patchIndices[i];
            dictionary& pDict = newPatchDicts[newPatchID];
            if (pDict.found("neighbourPatch"))
            {
                word neighbourName = pDict.lookup("neighbourPatch");
                label patchID = bm.findPatchID(neighbourName);
                if (patchID != -1 && renamedBC[patchID])
                {
                    word newPatchName = newPatchNames[patchIndices[patchID]];
                    pDict.add("neighbourPatch",newPatchName,true);
                }
            }
        }
    }

    //Add old processor faces
    forAll(patchNames_, patchI)
    {
        if (word(patchDicts_[patchI]["type"]) == "processor")
        {
            newPatchNames.append(patchNames_[patchI]);
            newPatchDicts.append(patchDicts_[patchI]);
            label newPatchI = newPatchNames.size() -1;
            if (patchI != newPatchI)
            {
                movedPatches.append(labelPair(patchI, newPatchI));
            }
        }
    }

    //Add new processor faces
    forAll(bm, i)
    {
        const polyPatch& p = bm[i];
        const word& pType = p.type();
        const word& pName = p.name();

        if (pType == "processor")
        {
            label oldPatchID = oldPatches.findPatchID(pName);
            if (oldPatchID != -1)
            {
                if (word(patchDicts_[oldPatchID]["type"]) == pType)
                {
                    forAll(newPatchNames, patchI)
                    {
                        if (pName == newPatchNames[patchI])
                        {
                            patchIndices[i] = patchI;
                            break;
                        }
                    }
                }
                else
                {
                    {
                        OStringStream os;
                        p.write(os);
                        newPatchDicts.append(dictionary(IStringStream(os.str())()));
                    }

                    // Duplicate name is not allowed.  Create a composite name from the
                    // patch name and case name
                    const word& caseName = p.boundaryMesh().mesh().time().caseName();

                    newPatchNames.append(pName + "_" + caseName);

                    Info<< "label patchIndex(const polyPatch& p) : "
                        << "Patch " << p.index() << " named "
                        << pName << " in mesh " << caseName
                        << " already exists, but patch types "
                        << " do not match.\nCreating a composite name as "
                        << newPatchNames.last() << endl;

                    patchIndices[i] = newPatchNames.size() - 1;
                }
            }
            else
            {
                // Patch not found.  Append to the list
                {
                    OStringStream os;
                    p.write(os);
                    newPatchDicts.append(dictionary(IStringStream(os.str())()));
                }
                newPatchNames.append(pName);

                patchIndices[i] = newPatchNames.size() - 1;
            }
        }
    }

    patchNames_.clear();
    patchDicts_.clear();
    patchNames_ = newPatchNames;
    patchDicts_ = newPatchDicts;
}


Foam::label Foam::mergePolyMesh::zoneIndex
(
    DynamicList<word>& names,
    const word& curName
)
{
    forAll(names, zonei)
    {
        if (names[zonei] == curName)
        {
            return zonei;
        }
    }

    // Not found.  Add new name to the list
    names.append(curName);

    return names.size() - 1;
}


void Foam::mergePolyMesh::sortProcessorPatches()
{
    Info<< "Reordering processor patches last" << endl;

    // Updates boundaryMesh() and meshMod_ to guarantee processor patches
    // are last. This could be done inside the merge() but it is far easier
    // to do separately.


    // 1. Shuffle the patches in the boundaryMesh

    const polyBoundaryMesh& oldPatches = boundaryMesh();

    DynamicList<polyPatch*> newPatches(oldPatches.size());

    labelList oldToSorted(oldPatches.size());

    forAll(oldPatches, patchi)
    {
        const polyPatch& pp = oldPatches[patchi];
        if (isA<directPolyPatch>(pp))
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(pp);

            if (!isA<processorPolyPatch>(dpp))
            {
                oldToSorted[patchi] = newPatches.size();
                newPatches.append
                (
                    dpp.clone
                    (
                        oldPatches,
                        oldToSorted[patchi],
                        0,
                        nInternalFaces()
                    ).ptr()
                );
            }
        }
    }
    forAll(oldPatches, patchi)
    {
        const polyPatch& pp = oldPatches[patchi];
        if (isA<directPolyPatch>(pp))
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(pp);

            if (isA<processorPolyPatch>(dpp))
            {
                oldToSorted[patchi] = newPatches.size();
                newPatches.append
                (
                    dpp.clone
                    (
                        oldPatches,
                        oldToSorted[patchi],
                        0,
                        nInternalFaces()
                    ).ptr()
                );
            }
        }
    }


    removeBoundary();
    addPatches(newPatches);


    // Update the polyTopoChange
    DynamicList<label>& patchID = const_cast<DynamicList<label>&>
    (
        meshMod_.region()
    );

    forAll(patchID, facei)
    {
        label patchi = patchID[facei];
        if (patchi != -1)
        {
            patchID[facei] = oldToSorted[patchID[facei]];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mergePolyMesh::mergePolyMesh(const IOobject& io)
:
    polyMesh(io),
    meshMod_(*this),
    patchNames_(2*boundaryMesh().size()),
    patchDicts_(2*boundaryMesh().size()),
    pointZoneNames_(),
    faceZoneNames_(),
    cellZoneNames_()
{
    // Insert the original patches into the list
    wordList curPatchNames = boundaryMesh().names();

    forAll(boundaryMesh(), patchi)
    {
        patchNames_.append(boundaryMesh()[patchi].name());

        OStringStream os;
        boundaryMesh()[patchi].write(os);
        patchDicts_.append(dictionary(IStringStream(os.str())()));
    }

    // Insert point, face and cell zones into the list

    // Point zones
    wordList curPointZoneNames = pointZones().names();
    if (curPointZoneNames.size())
    {
        pointZoneNames_.setCapacity(2*curPointZoneNames.size());
    }

    forAll(curPointZoneNames, zoneI)
    {
        pointZoneNames_.append(curPointZoneNames[zoneI]);
    }

    // Face zones
    wordList curFaceZoneNames = faceZones().names();

    if (curFaceZoneNames.size())
    {
        faceZoneNames_.setCapacity(2*curFaceZoneNames.size());
    }
    forAll(curFaceZoneNames, zoneI)
    {
        faceZoneNames_.append(curFaceZoneNames[zoneI]);
    }

    // Cell zones
    wordList curCellZoneNames = cellZones().names();

    if (curCellZoneNames.size())
    {
        cellZoneNames_.setCapacity(2*curCellZoneNames.size());
    }
    forAll(curCellZoneNames, zoneI)
    {
        cellZoneNames_.append(curCellZoneNames[zoneI]);
    }
}


Foam::mergePolyMesh::mergePolyMesh
(
    const IOobject& io,
    const pointField& points,
    const faceList& faces,
    const cellList& cells
)
:
    polyMesh
    (
        io,
        xferCopy(points),
        xferCopy(faces),
        xferCopy(cells)
    ),
    meshMod_(*this),
    patchNames_(0),
    patchDicts_(0),
    pointZoneNames_(0),
    faceZoneNames_(0),
    cellZoneNames_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mergePolyMesh::updateMeshMod(bool addMesh)
{
    // Clear topo change for the next operation
    meshMod_.clear();

    if (addMesh)
    {
        // Add current mesh as starting one for the next operation
        meshMod_.addMesh
        (
            *this,
            identity(boundaryMesh().size()),
            identity(pointZones().size()),
            identity(faceZones().size()),
            identity(cellZones().size())
         );
    }
}


void Foam::mergePolyMesh::addMesh(const polyMesh& m)
{
    // Add all the points, faces and cells of the new mesh

    // Add points

    label zoneID = -1;

    const pointField& p = m.points();
    labelList renumberPoints(p.size());

    const pointZoneMesh& pz = m.pointZones();
    labelList pointZoneIndices(pz.size());

    forAll(pz, zoneI)
    {
        pointZoneIndices[zoneI] = zoneIndex(pointZoneNames_, pz[zoneI].name());
    }

    forAll(p, pointi)
    {
        // Grab zone ID.  If a point is not in a zone, it will return -1
        zoneID = pz.whichZone(pointi);

        if (zoneID >= 0)
        {
            // Translate zone ID into the new index
            zoneID = pointZoneIndices[zoneID];
        }

        renumberPoints[pointi] =
            meshMod_.setAction
            (
                polyAddPoint
                (
                    p[pointi],            // Point to add
                    -1,                   // Master point (straight addition)
                    zoneID,               // Zone for point
                    pointi < m.nPoints()  // Is in cell?
                )
            );
    }

    // Add cells

    const cellList& c = m.cells();
    labelList renumberCells(c.size());

    const cellZoneMesh& cz = m.cellZones();
    labelList cellZoneIndices(cz.size());

    forAll(cz, zoneI)
    {
        cellZoneIndices[zoneI] = zoneIndex(cellZoneNames_, cz[zoneI].name());
    }

    forAll(c, celli)
    {
        // Grab zone ID.  If a cell is not in a zone, it will return -1
        zoneID = cz.whichZone(celli);

        if (zoneID >= 0)
        {
            // Translate zone ID into the new index
            zoneID = cellZoneIndices[zoneID];
        }

        renumberCells[celli] =
            meshMod_.setAction
            (
                polyAddCell
                (
                    -1,                   // Master point
                    -1,                   // Master edge
                    -1,                   // Master face
                    -1,                   // Master cell
                    zoneID                // Zone for cell
                )
            );
    }

    // Add faces
    const polyBoundaryMesh& bm = m.boundaryMesh();

    // Gather the patch indices
    labelList patchIndices(bm.size());
    DynamicList<labelPair> movedPatches(bm.size());
    setPatchIndices(bm, patchIndices, movedPatches);

    // Temporary: update number of allowable patches. This should be
    // determined at the top - before adding anything.
    meshMod_.setNumPatches(patchNames_.size());

    const faceZoneMesh& fz = m.faceZones();
    labelList faceZoneIndices(fz.size());

    forAll(fz, zoneI)
    {
        faceZoneIndices[zoneI] = zoneIndex(faceZoneNames_, fz[zoneI].name());
    }

    const faceList& f = m.faces();
    labelList renumberFaces(f.size());

    const labelList& own = m.faceOwner();
    const labelList& nei = m.faceNeighbour();

    label newOwn, newNei, newPatch, newZone;
    bool newZoneFlip;

    forAll(f, facei)
    {
        const face& curFace = f[facei];

        face newFace(curFace.size());

        forAll(curFace, pointi)
        {
            newFace[pointi] = renumberPoints[curFace[pointi]];
        }

        if (debug)
        {
            // Check that the face is valid
            if (min(newFace) < 0)
            {
                FatalErrorInFunction
                    << "Error in point mapping for face " << facei
                    << ".  Old face: " << curFace << " New face: " << newFace
                    << abort(FatalError);
            }
        }

        if (facei < m.nInternalFaces() || facei >= m.nFaces())
        {
            newPatch = -1;
        }
        else
        {
            newPatch = patchIndices[bm.whichPatch(facei)];
        }

        newOwn = own[facei];
        if (newOwn > -1) newOwn = renumberCells[newOwn];

        if (newPatch > -1)
        {
            newNei = -1;
        }
        else
        {
            newNei = nei[facei];
            newNei = renumberCells[newNei];
        }


        newZone = fz.whichZone(facei);
        newZoneFlip = false;

        if (newZone >= 0)
        {
            newZoneFlip = fz[newZone].flipMap()[fz[newZone].whichFace(facei)];

            // Grab the new zone
            newZone = faceZoneIndices[newZone];
        }

        renumberFaces[facei] =
            meshMod_.setAction
            (
                polyAddFace
                (
                    newFace,
                    newOwn,
                    newNei,
                    -1,
                    -1,
                    -1,
                    false,
                    newPatch,
                    newZone,
                    newZoneFlip
                )
            );
    }

    //If exiting processor patch has been moved update face patch ID's
    if (movedPatches.size())
    {
        const polyBoundaryMesh& oldPatches = boundaryMesh();
        const faceList& f = faces();
        const labelList& own = faceOwner();

        forAll(movedPatches, mPI)
        {
            label oldPatchI = movedPatches[mPI].first();
            label newPatchI = movedPatches[mPI].second();

            const polyPatch& pp = oldPatches[oldPatchI];
            label start = pp.start();
            forAll(pp, i)
            {
                label faceI = start+i;
                label zoneID = faceZones().whichZone(faceI);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones()[zoneID];
                    zoneFlip = fZone.flipMap()[zoneID];
                }

                meshMod_.setAction
                (
                    polyModifyFace
                    (
                        f[faceI],
                        faceI,
                        own[faceI],
                        -1,
                        false,
                        newPatchI,
                        false,
                        zoneID,
                        zoneFlip
                     )
                 );
            }
        }
    }
}


void Foam::mergePolyMesh::merge()
{
    Info<< "patch names: " << patchNames_ << nl
        << "patch dicts: " << patchDicts_ << nl
        << "point zone names: " << pointZoneNames_ << nl
        << "face zone names: " << faceZoneNames_ << nl
        << "cell zone names: " << cellZoneNames_ << endl;

    // Add the patches if necessary

    bool newPatchesFound = false;
    const polyBoundaryMesh& oldPatches = boundaryMesh();
    forAll(patchNames_, patchi)
    {
        word name = patchNames_[patchi];
        label oldPatchID = oldPatches.findPatchID(name);
        if (oldPatchID == -1)
        {
            newPatchesFound = true;
            break;
        }
    }
    reduce(newPatchesFound, orOp<bool>());

    if (newPatchesFound)
    {
        Info<< "Copying old patches" << endl;

        List<polyPatch*> newPatches(patchNames_.size());
        label endOfLastPatch = 0;
        const polyBoundaryMesh& oldPatches = boundaryMesh();

        forAll(patchNames_, patchi)
        {
            word name = patchNames_[patchi];

            label oldPatchID = oldPatches.findPatchID(name);

            if (oldPatchID != -1)
            {

                if (isA<directPolyPatch>(oldPatches[oldPatchID]))
                {
                    const directPolyPatch& dpp =
                        refCast<const directPolyPatch>(oldPatches[oldPatchID]);

                    newPatches[patchi] = dpp.clone(oldPatches).ptr();
                    endOfLastPatch = dpp.start() + dpp.size();
                }
            }
            else
            {
                // Add a new patch
                dictionary dict(patchDicts_[patchi]);
                dict.set("nFaces", 0);
                dict.set("startFace", endOfLastPatch);

                newPatches[patchi] =
                (
                    polyPatch::New
                    (
                        patchNames_[patchi],
                        dict,
                        patchi,
                        oldPatches
                    ).ptr()
                );
            }
        }

        removeBoundary();
        addPatches(newPatches);
    }

    // Add the zones if necessary
    if (pointZoneNames_.size() > pointZones().size())
    {
        Info<< "Adding new pointZones. " << endl;
        label nZones = pointZones().size();

        pointZones().setSize(pointZoneNames_.size());

        for (label zoneI = nZones; zoneI < pointZoneNames_.size(); zoneI++)
        {
            pointZones().set
            (
                zoneI,
                new pointZone
                (
                    pointZoneNames_[zoneI],
                    labelList(),
                    zoneI,
                    pointZones()
                )
            );
        }
    }
    if (cellZoneNames_.size() > cellZones().size())
    {
        Info<< "Adding new cellZones. " << endl;

        label nZones = cellZones().size();

        cellZones().setSize(cellZoneNames_.size());

        for (label zoneI = nZones; zoneI < cellZoneNames_.size(); zoneI++)
        {
            cellZones().set
            (
                zoneI,
                new cellZone
                (
                    cellZoneNames_[zoneI],
                    labelList(),
                    zoneI,
                    cellZones()
                )
            );
        }
    }
    if (faceZoneNames_.size() > faceZones().size())
    {
        Info<< "Adding new faceZones. " << endl;

        label nZones = faceZones().size();

        faceZones().setSize(faceZoneNames_.size());

        for (label zoneI = nZones; zoneI < faceZoneNames_.size(); zoneI++)
        {
            faceZones().set
            (
                zoneI,
                new faceZone
                (
                    faceZoneNames_[zoneI],
                    labelList(),
                    boolList(),
                    zoneI,
                    faceZones()
                )
            );
        }
    }


    // Shuffle the processor patches to be last
    sortProcessorPatches();

    // Change mesh. No inflation
    meshMod_.changeMesh(*this, false);
}


// ************************************************************************* //
