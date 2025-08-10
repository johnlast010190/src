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
    (c) 2010-2012 Esi Ltd.
    (c) 2008 Icon CG Ltd.

Description
    Makes internal faces into boundary faces. Does not duplicate points. Use
    mergeOrSplitBaffles if you want this.

    Note: if any coupled patch face is selected for baffling automatically
    the opposite member is selected for baffling as well. Note that this
    is the same as repatching. This was added only for convenience so
    you don't have to filter coupled boundary out of your set.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/syncTools/syncTools.H"
#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "sets/topoSets/faceSet.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "fields/ReadFields/ReadFields.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicPolyPatch.H"

#include "decompositionMethod/decompositionMethod.H"
#include "regionSplit/regionSplit.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void setFaces
(
    fvMesh& mesh,
    boolList& blockedFace,
    labelList& newPatch,
    labelList& swapMap,
    const dictionary& dict
)
{
    word jumpName = word(dict.lookup("name")) + "Master";

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    label wantedPatchI = patches.findPatchID(jumpName);

    Info<< "Using patch " << jumpName << " at index " << wantedPatchI << endl;

    if (wantedPatchI == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << jumpName << exit(FatalError);
    }

    bool geomCut = true;
    vector basePoint(vector::zero);
    vector jumpDir(vector::zero);
    scalar maxRadSqr(0.);
    word zoneName("empty");
    label zoneI = -1;

    //dimensionedVector dimBase("base", dimLength, basePoint);
    if (dict.found("basePoint") && dict.found("jumpDir") && dict.found("radius"))
    {
        basePoint = dict.lookup("basePoint");
        jumpDir = dict.lookup("jumpDir");
        maxRadSqr = sqr(readScalar(dict.lookup("radius")));
    }
    else if (dict.found("faceZone") && dict.found("jumpDir"))
    {
        zoneName = word(dict.lookup("faceZone"));
        jumpDir = dict.lookup("jumpDir");
        zoneI = mesh.faceZones().findZoneID(zoneName);

        if (zoneI == -1)
        {
            FatalErrorInFunction
                << "Could not find faceZone: " << zoneName << exit(FatalError);
        }
        else
        {
            Info<<"Using zone name: "<<mesh.faceZones()[zoneI].name()
                <<" to construct cyclic jump"<<endl;
        }

        geomCut = false;
    }
    else
    {
        FatalErrorInFunction
            << "Could not find method A requring keywords faceZone and jumpDir "
            << "Or method B requiring keywords basePoint, jumpDir and radius"
            << exit(FatalError);
    }

    jumpDir /= mag(jumpDir);

    vectorField disp( mesh.cellCentres()- basePoint );
    scalarField radSqr( magSqr(mesh.faceCentres() - basePoint) );

    // Calculate coupled cell centre

    vectorField neiCc(mesh.nFaces()-mesh.nInternalFaces());

    for (label faceI = mesh.nInternalFaces(); faceI < mesh.nFaces(); faceI++)
    {
        neiCc[faceI-mesh.nInternalFaces()] =
            mesh.cellCentres()[mesh.faceOwner()[faceI]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiCc);

    for (label faceI = 0; faceI < mesh.nFaces(); faceI++)
    {
        bool foundFace = false;
        if (geomCut && (radSqr[faceI] < maxRadSqr))
        {
            foundFace  = true;
        }
        else if (mesh.faceZones().whichZone(faceI) == zoneI)
        {
            foundFace = true;
        }

        if (foundFace)
        {
            if (mesh.isInternalFace(faceI))
            {
                label own = mesh.faceOwner()[faceI];
                label nei = mesh.faceNeighbour()[faceI];

                scalar oside = disp[own] & jumpDir;
                scalar nside = disp[nei] & jumpDir;
                if (!geomCut)
                {
                    oside = (mesh.cellCentres()[own] - mesh.faceCentres()[faceI])
                         & jumpDir;
                    nside = (mesh.cellCentres()[nei] - mesh.faceCentres()[faceI])
                        & jumpDir;
                }

                if (sign(oside) != sign(nside))
                {
                    newPatch[faceI] = wantedPatchI;
                    blockedFace[faceI] = false;
                    if (sign(oside) == sign(label(-1)))
                    {
                        swapMap[faceI] = 1;
                    }
                }
            }
            else if (patches[patches.whichPatch(faceI)].coupled())
            {
                label own = mesh.faceOwner()[faceI];

                scalar oside = disp[own] & jumpDir;
                vector neiDisp = neiCc[faceI-mesh.nInternalFaces()]
                    - basePoint;
                scalar nside =  neiDisp & jumpDir;
                if (!geomCut)
                {
                    oside = (mesh.cellCentres()[own] - mesh.faceCentres()[faceI])
                        & jumpDir;
                    nside = (neiCc[faceI-mesh.nInternalFaces()] - mesh.faceCentres()[faceI])
                        & jumpDir;
                }


                if (sign(oside) != sign(nside))
                {
                    newPatch[faceI] = wantedPatchI;
                    blockedFace[faceI] = false;
                    if (sign(oside) == sign(label(-1)))
                    {
                        swapMap[faceI] = 1;
                    }
                }
            }
        }
    }
}


// Main program:

// takes internal faces detected via cell-centre intersections of a circular
// disc and turns them into internal cyclic boundaries for use with jumpCyclic
// boundary conditions

int main(int argc, char *argv[])
{
#include "include/addOverwriteOption.H"

#include "include/setRootCase.H"
#include "include/createTime.H"
#include "include/createMesh.H"

    const bool overwrite = args.optionFound("overwrite");
    const word oldInstance = mesh.pointsInstance();

    Info<< "Reading createPatchDict\n" << endl;

    PtrList<dictionary> jumpDicts
    (
        IOdictionary
        (
            IOobject
            (
                "createJumpDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("patches")
    );

    // first add new patches
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Old and new patches.
    DynamicList<polyPatch*> allPatches(patches.size() + 2*jumpDicts.size());

    label pStart = 0;
    // Copy old patches.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if (isA<directPolyPatch>(pp))
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(pp);

            if (!isA<processorPolyPatch>(dpp))
            {
                allPatches.append
                (
                    dpp.clone
                    (
                        patches,
                        patchI,
                        dpp.size(),
                        dpp.start()
                    ).ptr()
                );
                pStart = dpp.start() + dpp.size();
            }
        }
    }

    forAll(jumpDicts, jumpI)
    {
        const dictionary& dict = jumpDicts[jumpI];
        word jumpName = dict.lookup("name");
        word masterJumpName = jumpName + word("Master");
        word slaveJumpName = jumpName + word("Slave");

        {
            label destPatchI = patches.findPatchID(masterJumpName);

            if (destPatchI == -1)
            {
                destPatchI = allPatches.size();

                Info<< "Adding new cyclic patch " << masterJumpName
                    << " as patch " << destPatchI << endl;

                dictionary masterPatchDict;
                masterPatchDict.add("type", "cyclic");
                masterPatchDict.add("startFace", pStart);
                masterPatchDict.add("nFaces", 0);
                masterPatchDict.add("matchTolerance", 0.0001);
                masterPatchDict.add("neighbourPatch", slaveJumpName);

                allPatches.append
                (
                    polyPatch::New
                    (
                        masterJumpName,
                        masterPatchDict,
                        destPatchI,
                        patches
                     ).ptr()
                );
            }
            else
            {
                Warning << "Patch name " << masterJumpName << " already exists." << nl
                        << masterJumpName << " will not be added as a new jump cyclic."
                        << endl;
            }
        }

        {
            label destPatchI = patches.findPatchID(slaveJumpName);

            if (destPatchI == -1)
            {
                destPatchI = allPatches.size();

                Info<< "Adding new cyclic patch " << slaveJumpName
                    << " as patch " << destPatchI << endl;

                dictionary slavePatchDict;
                slavePatchDict.add("type", "cyclic");
                slavePatchDict.add("startFace", pStart);
                slavePatchDict.add("nFaces", 0);
                slavePatchDict.add("matchTolerance", 0.0001);
                slavePatchDict.add("neighbourPatch", masterJumpName);

                allPatches.append
                (
                    polyPatch::New
                    (
                        slaveJumpName,
                        slavePatchDict,
                        destPatchI,
                        patches
                     ).ptr()
                );
            }
            else
            {
                Warning << "Patch name " << slaveJumpName << " already exists." << nl
              << slaveJumpName << " will not be added as a new jump cyclic."
              << endl;
            }
        }
    }

    // Copy processor patches.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<directPolyPatch>(pp))
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(pp);

            if (isA<processorPolyPatch>(dpp))
            {
                allPatches.append
                (
                    dpp.clone
                    (
                        patches,
                        patchI,
                        dpp.size(),
                        dpp.start()
                    ).ptr()
                );
            }
        }
    }

    allPatches.shrink();

    mesh.removeFvBoundary();
    mesh.addFvPatches(allPatches);

    //split the mesh and add the newly created faces to correct cyclic
    //boundaries in correct order

    labelList newPatch(mesh.nFaces(), -1);
    //wap owner-neighbour so jump direction is consistent
    labelList swapMap(mesh.nFaces(), 0);
    boolList blockedFace(mesh.nFaces(), true);

    // main loop
    forAll(jumpDicts, jumpI)
    {
        const dictionary& dict = jumpDicts[jumpI];

        // Creating cyclics:
        // - internal faces         : converted into boundary faces.

        setFaces(mesh, blockedFace, newPatch, swapMap, dict);
    }

    syncTools::syncFaceList(mesh, blockedFace, orEqOp<bool>());
    syncTools::syncFaceList(mesh, swapMap, maxEqOp<label>());
    syncTools::syncFaceList(mesh, newPatch, maxEqOp<label>());

    if (Pstream::parRun())
    {
        // Read decomposePar dictionary
        IOdictionary decomposeDict
        (
            IOobject
            (
                "decomposeParDict",
                runTime.system(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        // Decomposition
        autoPtr<decompositionMethod> decomposerPtr = decompositionMethod::New
        (
            decomposeDict
        );
        decompositionMethod& decomposer = decomposerPtr();

        //- Mesh distribution engine
        autoPtr<fvMeshDistribute> distributorPtr;
        // Mesh distribution engine (uses tolerance to reconstruct meshes)
        const scalar mergeDist = 1e-8;
        distributorPtr.reset(new fvMeshDistribute(mesh, mergeDist));

        // Wanted distribution
        labelList distribution;
        List<labelPair> explicitConnections;

        // Faces that move as block onto single processor
        PtrList<labelList> specifiedProcessorFaces;
        labelList specifiedProcessor;

        distribution = decomposer.decompose
        (
            mesh,
            scalarField(mesh.nCells(), 1),
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            explicitConnections   // explicit connections
        );

        // Do actual sending/receiving of mesh
        autoPtr<mapDistributePolyMesh> map
            = distributorPtr().distribute(distribution);
        map().distributeFaceData(newPatch);
        map().distributeFaceData(swapMap);
    }


    // Mesh change container
    polyTopoChange meshMod(mesh);

    //label nSplit = 0;

    //first add jump-start faces
    forAll(newPatch, faceI)
    {
        if (newPatch[faceI] != -1)
        {
            const face& f = mesh.faces()[faceI];
            label zoneID = mesh.faceZones().whichZone(faceI);
            bool zoneFlip = false;
            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            if (swapMap[faceI])
            {
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        f,                          // modified face
                        faceI,                      // label of face
                        mesh.faceOwner()[faceI],    // owner
                        -1,                         // neighbour
                        false,                      // face flip
                        newPatch[faceI],            // patch for face
                        false,                      // remove from zone
                        zoneID,                     // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }
            else
            {
                meshMod.setAction
                (
                    polyAddFace
                    (
                        f.reverseFace(),            // modified face
                        mesh.faceNeighbour()[faceI],// owner
                        -1,                         // neighbour
                        -1,                         // masterPointID
                        -1,                         // masterEdgeID
                        faceI,                      // masterFaceID,
                        false,                      // face flip
                        newPatch[faceI],            // patch for face
                        zoneID,                     // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }
        }

        //nSplit++;
    }

    //now add jump-destination faces
    forAll(newPatch, faceI)
    {
        if (newPatch[faceI] != -1)
        {
            const face& f = mesh.faces()[faceI];
            label zoneID = mesh.faceZones().whichZone(faceI);
            bool zoneFlip = false;
            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(faceI)];
            }

            if (swapMap[faceI])
            {
                meshMod.setAction
                (
                    polyAddFace
                    (
                        f.reverseFace(),            // modified face
                        mesh.faceNeighbour()[faceI],// owner
                        -1,                         // neighbour
                        -1,                         // masterPointID
                        -1,                         // masterEdgeID
                        faceI,                      // masterFaceID,
                        false,                      // face flip
                        newPatch[faceI]+1,          // patch for face
                        zoneID,                     // zone for face
                        !zoneFlip                    // face flip in zone
                    )
                );
            }
            else
            {
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        f,                          // modified face
                        faceI,                      // label of face
                        mesh.faceOwner()[faceI],    // owner
                        -1,                         // neighbour
                        false,                      // face flip
                        newPatch[faceI]+1,          // patch for face
                        false,                      // remove from zone
                        zoneID,                     // zone for face
                        zoneFlip                    // face flip in zone
                    )
                );
            }
        }

        //nSplit++;
    }

    // Change the mesh. Change points directly (no inflation).
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

    if (!overwrite)
    {
        runTime++;
    }
    else
    {
       mesh.setInstance(oldInstance);
    }

    Info<< "Writing mesh." << endl;

    mesh.write();

    Info<< "End\n" << endl;

    return 0;

}


// ************************************************************************* //
