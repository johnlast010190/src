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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2023 Esi Ltd.

Application
    mergeOrSplitBaffles

Group
    grpMeshManipulationUtilities

Description
    Detects boundary faces that share points (baffles). Either merges them or
    duplicate the points.

Usage
    \b mergeOrSplitBaffles [OPTION]

    Options:
      - \par -detect
        Detect baffles and write to faceSet duplicateFaces.

      - \par -merge
        Detect baffles and convert to internal faces.

      - \par -split
        Detect baffles and duplicate the points (used so the two sides
        can move independently)

      - \par -noFields
        Do not update fields.

      - \par -dict \<dictionary\>
        Specify a dictionary to read actions from.


Note
    - can only handle pairwise boundary faces. So three faces using
      the same points is not handled (is illegal mesh anyway)

    - surfaces consisting of duplicate faces can be topologically split
    if the points on the interior of the surface cannot walk to all the
    cells that use them in one go.

    - Parallel operation (where duplicate face is perpendicular to a coupled
    boundary) is supported but not really tested.
    (Note that coupled faces themselves are not seen as duplicate faces)

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"
#include "meshTools/meshTools.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "regionSplit/localPointRegion.H"
#include "polyTopoChange/duplicatePoints/duplicatePoints.H"
#include "polyTopoChange/hexRef8/hexRef8Data.H"
#include "fields/ReadFields/ReadFields.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "processorMeshes.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void insertDuplicateMerge
(
    const polyMesh& mesh,
    const labelList& boundaryFaces,
    const labelList& duplicates,
    polyTopoChange& meshMod
)
{
    const faceList& faces = mesh.faces();
    const labelList& faceOwner = mesh.faceOwner();
    const faceZoneMesh& faceZones = mesh.faceZones();

    forAll(duplicates, bFacei)
    {
        label otherFacei = duplicates[bFacei];

        if (otherFacei != -1 && otherFacei > bFacei)
        {
            // Two duplicate faces. Merge.

            label face0 = boundaryFaces[bFacei];
            label face1 = boundaryFaces[otherFacei];

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (own0 < own1)
            {
                // Use face0 as the new internal face.
                label zoneID = faceZones.whichZone(face0);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                }

                meshMod.setAction(polyRemoveFace(face1));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face0],           // modified face
                        face0,                  // label of face being modified
                        own0,                   // owner
                        own1,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
            else
            {
                // Use face1 as the new internal face.
                label zoneID = faceZones.whichZone(face1);
                bool zoneFlip = false;

                if (zoneID >= 0)
                {
                    const faceZone& fZone = faceZones[zoneID];
                    zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                }

                meshMod.setAction(polyRemoveFace(face0));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face1],           // modified face
                        face1,                  // label of face being modified
                        own1,                   // owner
                        own0,                   // neighbour
                        false,                  // face flip
                        -1,                     // patch for face
                        false,                  // remove from zone
                        zoneID,                 // zone for face
                        zoneFlip                // face flip in zone
                    )
                );
            }
        }
    }
}


label patchSize(const polyMesh& mesh, const labelList& patchIDs)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    label sz = 0;
    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];
        sz += pp.size();
    }
    return sz;
}


labelList patchFaces(const polyMesh& mesh, const labelList& patchIDs)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    labelList faceIDs(patchSize(mesh, patchIDs));
    label sz = 0;
    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];

        forAll(pp, ppi)
        {
            faceIDs[sz++] = pp.start()+ppi;
        }
    }

if (faceIDs.size() != sz)
{
    FatalErrorInFunction << exit(FatalError);
}

    return faceIDs;
}


labelList findBaffles(const polyMesh& mesh, const labelList& boundaryFaces)
{
    // Get all duplicate face labels (in boundaryFaces indices!).
    labelList duplicates = localPointRegion::findDuplicateFaces
    (
        mesh,
        boundaryFaces
    );


    // Check that none are on processor patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(duplicates, bFacei)
    {
        if (duplicates[bFacei] != -1)
        {
            label facei = boundaryFaces[bFacei];
            label patchi = patches.whichPatch(facei);

            if (isA<processorPolyPatch>(patches[patchi]))
            {
                FatalErrorInFunction
                    << "Duplicate face " << facei
                    << " is on a processorPolyPatch."
                    << "This is not allowed." << nl
                    << "Face:" << facei
                    << " is on patch:" << patches[patchi].name()
                    << abort(FatalError);
            }
        }
    }


    // Write to faceSet for ease of postprocessing.
    {
        faceSet duplicateSet
        (
            mesh,
            "duplicateFaces",
            (mesh.nFaces() - mesh.nInternalFaces())/256
        );

        forAll(duplicates, bFacei)
        {
            label otherFacei = duplicates[bFacei];

            if (otherFacei != -1 && otherFacei > bFacei)
            {
                duplicateSet.insert(boundaryFaces[bFacei]);
                duplicateSet.insert(boundaryFaces[otherFacei]);
            }
        }

        Info<< "Writing " << returnReduce(duplicateSet.size(), sumOp<label>())
            << " duplicate faces to faceSet " << duplicateSet.objectPath()
            << nl << endl;
        duplicateSet.write();
    }

    return duplicates;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Detect faces that share points (baffles).\n"
        "Merge them or duplicate the points."
    );

    #include "include/addOverwriteOption.H"
    #include "include/addRegionOption.H"
    #include "include/addDictOption.H"
    argList::addBoolOption
    (
        "detectOnly",
        "find baffles only, but do not merge or split them"
    );
    argList::addBoolOption
    (
        "split",
        "topologically split duplicate surfaces"
    );
    argList::addBoolOption
    (
        "noFields",
        "do not update fields"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    runTime.functionObjects().off();
    #include "include/createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const bool readDict   = args.optionFound("dict");
    const bool split      = args.optionFound("split");
    const bool overwrite  = args.optionFound("overwrite");
    const bool detectOnly = args.optionFound("detectOnly");
    const bool fields     = !args.optionFound("noFields");

    if (readDict && (split || detectOnly))
    {
        FatalErrorInFunction
            << "Use of dictionary for settings not compatible with"
            << " using command line arguments for \"split\""
            << " or \"detectOnly\"" << exit(FatalError);
    }


    labelList detectPatchIDs;
    labelList splitPatchIDs;
    labelList mergePatchIDs;

    if (readDict)
    {
        const word dictName;
        #include "include/setSystemMeshDictionaryIO.H"

        Info<< "Reading " << dictName << "\n" << endl;
        IOdictionary dict(dictIO);

        if (dict.found("detect"))
        {
            wordReList patchNames(dict.subDict("detect").lookup("patches"));
            detectPatchIDs = patches.patchSet(patchNames).sortedToc();
            Info<< "Detecting baffles on " << detectPatchIDs.size()
                << " patches with "
                << returnReduce(patchSize(mesh, detectPatchIDs), sumOp<label>())
                << " faces" << endl;
        }
        if (dict.found("merge"))
        {
            wordReList patchNames(dict.subDict("merge").lookup("patches"));
            mergePatchIDs = patches.patchSet(patchNames).sortedToc();
            Info<< "Detecting baffles on " << mergePatchIDs.size()
                << " patches with "
                << returnReduce(patchSize(mesh, mergePatchIDs), sumOp<label>())
                << " faces" << endl;
        }
        if (dict.found("split"))
        {
            wordReList patchNames(dict.subDict("split").lookup("patches"));
            splitPatchIDs = patches.patchSet(patchNames).sortedToc();
            Info<< "Detecting baffles on " << splitPatchIDs.size()
                << " patches with "
                << returnReduce(patchSize(mesh, splitPatchIDs), sumOp<label>())
                << " faces" << endl;
        }
    }
    else
    {
        if (detectOnly)
        {
            detectPatchIDs = identity(patches.size());
        }
        else if (split)
        {
            splitPatchIDs = identity(patches.size());
        }
        else
        {
            mergePatchIDs = identity(patches.size());
        }
    }


    if (detectPatchIDs.size())
    {
        findBaffles(mesh, patchFaces(mesh, detectPatchIDs));

        if (detectOnly)
        {
            return 0;
        }
    }

    // Read objects in time directory
    IOobjectList objects(mesh, runTime.timeName());

    if (fields) Info<< "Reading geometric fields" << nl << endl;

    PtrList<volScalarField> vsFlds;
    if (fields) ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    if (fields) ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    if (fields) ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    if (fields) ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    if (fields) ReadFields(mesh, objects, vtFlds);


    PtrList<surfaceScalarField> ssFlds;
    if (fields) ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    if (fields) ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    if (fields) ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    if (fields) ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    if (fields) ReadFields(mesh, objects, stFlds);


    if (mergePatchIDs.size())
    {
        Info<< "Merging duplicate faces" << nl << endl;

        // Mesh change engine
        polyTopoChange meshMod(mesh);

        const labelList boundaryFaces(patchFaces(mesh, mergePatchIDs));

        // Get all duplicate face pairs (in boundaryFaces indices!).
        labelList duplicates(findBaffles(mesh, boundaryFaces));

        // Merge into internal faces.
        insertDuplicateMerge(mesh, boundaryFaces, duplicates, meshMod);

        if (!overwrite)
        {
            runTime++;
        }

        // Change the mesh. No inflation.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }

        if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }
        Info<< "Writing mesh to time " << runTime.timeName() << endl;
        mesh.write();
    }


    if (splitPatchIDs.size())
    {
        Info<< "Topologically splitting duplicate surfaces"
            << ", i.e. duplicating points internal to duplicate surfaces"
            << nl << endl;

        // Determine points on split patches
        DynamicList<label> candidates;
        {
            label sz = 0;
            forAll(splitPatchIDs, i)
            {
                sz += patches[splitPatchIDs[i]].nPoints();
            }
            candidates.setCapacity(sz);

            PackedBoolList isCandidate(mesh.nPoints());
            forAll(splitPatchIDs, i)
            {
                const labelList& mp = patches[splitPatchIDs[i]].meshPoints();
                forAll(mp, mpi)
                {
                    label pointi = mp[mpi];
                    if (isCandidate.set(pointi))
                    {
                        candidates.append(pointi);
                    }
                }
            }
        }

        // Analyse which points need to be duplicated
        localPointRegion regionSide(mesh, candidates);

        // Point duplication engine
        duplicatePoints pointDuplicator(mesh);

        // Mesh change engine
        polyTopoChange meshMod(mesh);

        // Insert topo changes
        pointDuplicator.setRefinement(regionSide, meshMod);

        hexRef8Data refData
        (
            IOobject
            (
                "dummy",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        if (!overwrite)
        {
            runTime++;
        }

        // Change the mesh. No inflation.
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false);

        // Update fields
        mesh.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }

        if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }
        Info<< "Writing mesh to time " << runTime.timeName() << endl;
        mesh.write();

        topoSet::removeFiles(mesh);
        processorMeshes::removeFiles(mesh);

        // Dump duplicated points (if any)
        const labelList& pointMap = map().pointMap();

        labelList nDupPerPoint(map().nOldPoints(), 0);

        pointSet dupPoints(mesh, "duplicatedPoints", 100);

        forAll(pointMap, pointi)
        {
            label oldPointi = pointMap[pointi];

            nDupPerPoint[oldPointi]++;

            if (nDupPerPoint[oldPointi] > 1)
            {
                dupPoints.insert(map().reversePointMap()[oldPointi]);
                dupPoints.insert(pointi);
            }
        }

        Info<< "Writing " << returnReduce(dupPoints.size(), sumOp<label>())
            << " duplicated points to pointSet "
            << dupPoints.objectPath() << nl << endl;

        dupPoints.write();

        // Update refinement data
        refData.updateMesh(map);
        refData.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
