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
    (c) 2011-2018 OpenFOAM Foundation
    (c) 2015-2018 OpenCFD Ltd.
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMeshDistribute/fvMeshDistribute.H"
#include "db/IOstreams/Pstreams/PstreamCombineReduceOps.H"
#include "fvMeshAdder/fvMeshAdder.H"
#include "polyMeshAdder/faceCoupleInfo.H"
#include "fields/fvPatchFields/constraint/processor/processorFvPatchField.H"
#include "fields/fvsPatchFields/constraint/processor/processorFvsPatchField.H"
#include "meshes/polyMesh/polyPatches/constraint/processorCyclic/processorCyclicPolyPatch.H"
#include "fields/fvPatchFields/constraint/processorCyclic/processorCyclicFvPatchField.H"
#include "fields/fvPatchFields/constraint/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvPatchField.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "containers/Lists/CompactListList/CompactListList.H"
#include "fvMeshTools/fvMeshTools.H"
#include "primitives/Pair/labelPairHashes.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"
#include "meshes/meshTools/matchPoints.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMeshEntries.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshDistribute, 0);

    //- Less function class that can be used for sorting processor patches
    class lessProcPatches
    {
        const labelList& nbrProc_;
        const labelList& referPatchID_;
        const labelList& referNbrPatchID_;

    public:

        lessProcPatches
        (
            const labelList& nbrProc,
            const labelList& referPatchID,
            const labelList& referNbrPatchID
        )
        :
            nbrProc_(nbrProc),
            referPatchID_(referPatchID),
            referNbrPatchID_(referNbrPatchID)
        {}

        bool operator()(const label a, const label b)
        {
            // Lower processor ID's go first
            if (nbrProc_[a] < nbrProc_[b])
            {
                return true;
            }
            else if (nbrProc_[a] > nbrProc_[b])
            {
                return false;
            }

            // Non-cyclics go next
            else if (referPatchID_[a] == -1)
            {
                return true;
            }
            else if (referPatchID_[b] == -1)
            {
                return false;
            }

            // Cyclics should be ordered by refer patch ID if this is the owner
            // (lower processor ID), and by the neighbour refer patch ID if
            // this is the neighbour
            else if (Pstream::myProcNo() < nbrProc_[a])
            {
                return referPatchID_[a] < referPatchID_[b];
            }
            else
            {
                return referNbrPatchID_[a] < referNbrPatchID_[b];
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistribute::inplaceRenumberWithFlip
(
    const labelUList& oldToNew,
    const bool oldToNewHasFlip,
    const bool lstHasFlip,
    labelUList& lst
)
{
    if (!lstHasFlip && !oldToNewHasFlip)
    {
        Foam::inplaceRenumber(oldToNew, lst);
    }
    else
    {
        // Either input data or map encodes sign so result encodes sign

        forAll(lst, elemI)
        {
            // Extract old value and sign
            label val = lst[elemI];
            label sign = 1;
            if (lstHasFlip)
            {
                if (val > 0)
                {
                    val = val-1;
                }
                else if (val < 0)
                {
                    val = -val-1;
                    sign = -1;
                }
                else
                {
                    FatalErrorInFunction
                        << "Problem : zero value " << val
                        << " at index " << elemI << " out of " << lst.size()
                        << " list with flip bit" << exit(FatalError);
                }
            }


            // Lookup new value and possibly change sign
            label newVal = oldToNew[val];

            if (oldToNewHasFlip)
            {
                if (newVal > 0)
                {
                    newVal = newVal-1;
                }
                else if (newVal < 0)
                {
                    newVal = -newVal-1;
                    sign = -sign;
                }
                else
                {
                    FatalErrorInFunction
                        << "Problem : zero value " << newVal
                        << " at index " << elemI << " out of "
                        << oldToNew.size()
                        << " list with flip bit" << exit(FatalError);
                }
            }


            // Encode new value and sign
            lst[elemI] = sign*(newVal+1);
        }
    }
}


Foam::labelList Foam::fvMeshDistribute::select
(
    const bool selectEqual,
    const labelList& values,
    const label value
)
{
    label n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            n++;
        }
    }

    labelList indices(n);
    n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            indices[n++] = i;
        }
    }
    return indices;
}


void Foam::fvMeshDistribute::checkEqualWordList
(
    const string& msg,
    const wordList& lst
)
{
#ifndef RUNTIME_DEBUG_FLAGS
    return;
#endif

    // Check all procs have same names and in exactly same order.

    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = lst;
    Pstream::allGatherList(allNames);

    for (label proci = 1; proci < Pstream::nProcs(); proci++)
    {
        if (allNames[proci] != allNames[0])
        {
            FatalErrorInFunction
                << "When checking for equal " << msg.c_str() << " :" << endl
                << "processor0 has:" << allNames[0] << endl
                << "processor" << proci << " has:" << allNames[proci] << endl
                << msg.c_str() << " need to be synchronised on all processors."
                << exit(FatalError);
        }
    }
}


Foam::wordList Foam::fvMeshDistribute::mergeWordList(const wordList& procNames)
{
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = procNames;
    Pstream::allGatherList(allNames);

    wordHashSet mergedNames;
    forAll(allNames, proci)
    {
        mergedNames.insert(allNames[proci]);
    }
    return mergedNames.toc();
}


void Foam::fvMeshDistribute::printMeshInfo(const fvMesh& mesh)
{
    Pout<< nl << "Mesh: " << mesh.name() << nl
        << "Primitives:" << nl
        << "    points       :" << mesh.nPoints() << nl
        << "    bb           :" << boundBox(mesh.points(), false) << nl
        << "    internalFaces:" << mesh.nInternalFaces() << nl
        << "    faces        :" << mesh.nFaces() << nl
        << "    cells        :" << mesh.nCells() << nl;

    const fvBoundaryMesh& patches = mesh.boundary();

    Pout<< "Patches:" << endl;
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi].patch();

        Pout<< "    " << patchi << " name:" << pp.name()
            << " size:" << pp.size()
            << " start:" << pp.start()
            << " type:" << pp.type()
            << endl;
    }

    if (mesh.pointZones().size())
    {
        Pout<< "PointZones:" << endl;
        forAll(mesh.pointZones(), zoneI)
        {
            const pointZone& pz = mesh.pointZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << pz.name()
                << " size:" << pz.size()
                << endl;
        }
    }
    if (mesh.faceZones().size())
    {
        Pout<< "FaceZones:" << endl;
        forAll(mesh.faceZones(), zoneI)
        {
            const faceZone& fz = mesh.faceZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << fz.name()
                << " size:" << fz.size()
                << endl;
        }
    }
    if (mesh.cellZones().size())
    {
        Pout<< "CellZones:" << endl;
        forAll(mesh.cellZones(), zoneI)
        {
            const cellZone& cz = mesh.cellZones()[zoneI];
            Pout<< "    " << zoneI << " name:" << cz.name()
                << " size:" << cz.size()
                << endl;
        }
    }
}


void Foam::fvMeshDistribute::printCoupleInfo
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNewNbrProc
)
{
    Pout<< nl
        << "Current coupling info:"
        << endl;

    forAll(sourceFace, bFacei)
    {
        label meshFacei = mesh.nInternalFaces() + bFacei;

        Pout<< "    meshFace:" << meshFacei
            << " fc:" << mesh.faceCentres()[meshFacei]
            << " connects to proc:" << sourceProc[bFacei]
            << "/face:" << sourceFace[bFacei]
            << " which will move to proc:" << sourceNewNbrProc[bFacei]
            << endl;
    }
}


Foam::label Foam::fvMeshDistribute::findNonEmptyPatch() const
{
    // Finds (non-empty) patch that exposed internal and proc faces can be
    // put into.
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    label nonEmptyPatchi = -1;

    forAllReverse(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
            !isA<emptyPolyPatch>(pp)
         && !isA<cyclicAMIPolyPatch>(pp)
         && !isA<nonConformalPolyPatch>(pp)
         && !pp.coupled()
        )
        {
            nonEmptyPatchi = patchi;
            break;
        }
    }

    if (nonEmptyPatchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find a patch which is neither of type empty, cyclicAMI"
            << " nor coupled in patches " << patches.names() << endl
            << "There has to be at least one such patch for"
            << " distribution to work" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "findNonEmptyPatch : using patch " << nonEmptyPatchi
            << " name:" << patches[nonEmptyPatchi].name()
            << " type:" << patches[nonEmptyPatchi].type()
            << " to put exposed faces into." << endl;
    }


    // Do additional test for processor patches intermingled with non-proc
    // patches.
    label procPatchi = -1;

    forAll(patches, patchi)
    {
        if (isA<processorPolyPatch>(patches[patchi]))
        {
            procPatchi = patchi;
        }
        else if (procPatchi != -1)
        {
            FatalErrorInFunction
                << "Processor patches should be at end of patch list."
                << endl
                << "Have processor patch " << procPatchi
                << " followed by non-processor patch " << patchi
                << " in patches " << patches.names()
                << abort(FatalError);
        }
    }

    return nonEmptyPatchi;
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvMeshDistribute::generateTestField
(
    const fvMesh& mesh
)
{
    vector testNormal(1, 1, 1);
    testNormal /= mag(testNormal);

    tmp<surfaceScalarField> tfld
    (
        new surfaceScalarField
        (
            IOobject
            (
                "myFlux",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimless, Zero)
        )
    );
    surfaceScalarField& fld = tfld.ref();

    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    forAll(fld, facei)
    {
        fld[facei] = (n[facei] & testNormal);
    }

    surfaceScalarField::Boundary& fluxBf = fld.boundaryFieldRef();
    const surfaceVectorField::Boundary& nBf = n.boundaryField();

    forAll(fluxBf, patchi)
    {
        fvsPatchScalarField& fvp = fluxBf[patchi];

        scalarField newPfld(fvp.size());
        forAll(newPfld, i)
        {
            newPfld[i] = (nBf[patchi][i] & testNormal);
        }
        fvp.forceAssign(newPfld);
    }

    return tfld;
}


void Foam::fvMeshDistribute::testField(const surfaceScalarField& fld)
{
    const fvMesh& mesh = fld.mesh();

    vector testNormal(1, 1, 1);
    testNormal /= mag(testNormal);

    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    forAll(fld, facei)
    {
        scalar cos = (n[facei] & testNormal);

        if (mag(cos - fld[facei]) > 1e-6)
        {
            //FatalErrorInFunction
            WarningInFunction
                << "On internal face " << facei << " at "
                << mesh.faceCentres()[facei]
                << " the field value is " << fld[facei]
                << " whereas cos angle of " << testNormal
                << " with mesh normal " << n[facei]
                << " is " << cos
                //<< exit(FatalError);
                << endl;
        }
    }
    forAll(fld.boundaryField(), patchi)
    {
        const fvsPatchScalarField& fvp = fld.boundaryField()[patchi];
        const fvsPatchVectorField& np = n.boundaryField()[patchi];

        forAll(fvp, i)
        {
            scalar cos = (np[i] & testNormal);

            if (mag(cos - fvp[i]) > 1e-6)
            {
                label facei = fvp.patch().start()+i;
                //FatalErrorInFunction
                WarningInFunction
                    << "On face " << facei
                    << " on patch " << fvp.patch().name()
                    << " at " << mesh.faceCentres()[facei]
                    << " the field value is " << fvp[i]
                    << " whereas cos angle of " << testNormal
                    << " with mesh normal " << np[i]
                    << " is " << cos
                    //<< exit(FatalError);
                    << endl;
            }
        }
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::deleteProcPatches
(
    const label destinationPatch
)
{
    // Delete all processor patches. Move any processor faces into the last
    // non-processor patch.

    // New patchID per boundary faces to be repatched. Is -1 (no change)
    // or new patchID
    labelList newPatchID(mesh_.nBoundaryFaces(), -1);

//    label nProcPatches = 0;

    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            if (debug)
            {
                Pout<< "Moving all faces of patch " << pp.name()
                    << " into patch " << destinationPatch
                    << endl;
            }

            label offset = pp.start() - mesh_.nInternalFaces();

            forAll(pp, i)
            {
                newPatchID[offset+i] = destinationPatch;
            }

//            nProcPatches++;
        }
    }

    // Note: order of boundary faces been kept the same since the
    // destinationPatch is at the end and we have visited the patches in
    // incremental order.
    labelListList dummyFaceMaps;
    autoPtr<mapPolyMesh> map = repatch(newPatchID, mesh_, dummyFaceMaps);


    // Delete (now empty) processor patches.
    {
        labelList oldToNew(identity(mesh_.boundaryMesh().size()));
        label newi = 0;
        // Non processor patches first
        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (!isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
            {
                oldToNew[patchi] = newi++;
            }
        }
        label nNonProcPatches = newi;

        // Processor patches as last
        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
            {
                oldToNew[patchi] = newi++;
            }
        }
        fvMeshTools::reorderPatches(mesh_, oldToNew, nNonProcPatches, false);
    }

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::repatch
(
    const labelList& newPatchID,         // per boundary face -1 or new patchID
    fvMesh& mesh,
    labelListList& constructFaceMap
)
{
    polyTopoChange meshMod(mesh);

    forAll(newPatchID, bFacei)
    {
        if (newPatchID[bFacei] != -1)
        {
            label facei = mesh.nInternalFaces() + bFacei;

            label zoneID = mesh.faceZones().whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh.faces()[facei],       // modified face
                    facei,                      // label of face
                    mesh.faceOwner()[facei],   // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    newPatchID[bFacei],         // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                )
            );
        }
    }


    // Do mapping of fields from one patchField to the other ourselves since
    // is currently not supported by updateMesh.

    // Store boundary fields (we only do this for surfaceFields)
    PtrList<FieldField<fvsPatchField, scalar>> sFlds;
    saveBoundaryFields<scalar, surfaceMesh>(sFlds, mesh);
    PtrList<FieldField<fvsPatchField, vector>> vFlds;
    saveBoundaryFields<vector, surfaceMesh>(vFlds, mesh);
    PtrList<FieldField<fvsPatchField, sphericalTensor>> sptFlds;
    saveBoundaryFields<sphericalTensor, surfaceMesh>(sptFlds, mesh);
    PtrList<FieldField<fvsPatchField, symmTensor>> sytFlds;
    saveBoundaryFields<symmTensor, surfaceMesh>(sytFlds, mesh);
    PtrList<FieldField<fvsPatchField, tensor>> tFlds;
    saveBoundaryFields<tensor, surfaceMesh>(tFlds, mesh);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    //
    // NOTE: there is one very particular problem with this ordering.
    // We first create the processor patches and use these to merge out
    // shared points (see mergeSharedPoints below). So temporarily points
    // and edges do not match!

    bool oldParRun = UPstream::parRun();

    bool syncPar = true;

    if (nFinalProcs_ > Pstream::nProcs())
    {
        UPstream::parRun() = false;
        syncPar = false;
    }

    // Change the mesh (no inflation). Note: parallel comms allowed (if one mesh pr proc).
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, syncPar);

    // Update fields. No inflation, parallel sync (if one mesh per proc).
    mesh.updateMesh(map);
    UPstream::parRun() = oldParRun;

    // Map patch fields using stored boundary fields. Note: assumes order
    // of fields has not changed in object registry!
    mapBoundaryFields<scalar, surfaceMesh>(map, sFlds, mesh);
    mapBoundaryFields<vector, surfaceMesh>(map, vFlds, mesh);
    mapBoundaryFields<sphericalTensor, surfaceMesh>(map, sptFlds, mesh);
    mapBoundaryFields<symmTensor, surfaceMesh>(map, sytFlds, mesh);
    mapBoundaryFields<tensor, surfaceMesh>(map, tFlds, mesh);


    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }

    // Adapt constructMaps.

    if (debug)
    {
        label index = findIndex(map().reverseFaceMap(), -1);

        if (index != -1)
        {
            FatalErrorInFunction
                << "reverseFaceMap contains -1 at index:"
                << index << endl
                << "This means that the repatch operation was not just"
                << " a shuffle?" << abort(FatalError);
        }
    }

    forAll(constructFaceMap, proci)
    {
        inplaceRenumberWithFlip
        (
            map().reverseFaceMap(),
            false,
            true,
            constructFaceMap[proci]
        );
    }

    return map;
}


// Detect shared points. Need processor patches to be present.
// Background: when adding bits of mesh one can get points which
// share the same position but are only detectable to be topologically
// the same point when doing parallel analysis. This routine will
// merge those points.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::mergeSharedPoints
(
    const labelList& pointToGlobalMaster,
    fvMesh& mesh,
    labelListList& constructPointMap
)
{
    // Find out which sets of points get merged and create a map from
    // mesh point to unique point.

    label nShared = 0;
    forAll(pointToGlobalMaster, pointi)
    {
        if (pointToGlobalMaster[pointi] != -1)
        {
            nShared++;
        }
    }

    Map<label> globalMasterToLocalMaster(2*nShared);
    Map<label> pointToMaster(2*nShared);

    forAll(pointToGlobalMaster, pointi)
    {
        label globali = pointToGlobalMaster[pointi];
        if (globali != -1)
        {
            Map<label>::const_iterator iter = globalMasterToLocalMaster.find
            (
                globali
            );

            if (iter == globalMasterToLocalMaster.end())
            {
                // Found first point. Designate as master
                globalMasterToLocalMaster.insert(globali, pointi);
                pointToMaster.insert(pointi, pointi);
            }
            else
            {
                pointToMaster.insert(pointi, iter());
            }
        }
    }

    if
    (
        nFinalProcs_ > Pstream::nProcs()
     && pointToMaster.size() == 0
    )
    {
        return autoPtr<mapPolyMesh>(nullptr);
    }
    else if
    (
        nFinalProcs_ <= Pstream::nProcs()
     && (returnReduce(pointToMaster.size(), sumOp<label>()) == 0)
    )
    {
        return autoPtr<mapPolyMesh>(nullptr);
    }

    polyTopoChange meshMod(mesh);

    fvMeshAdder::mergePoints(mesh, pointToMaster, meshMod);

    bool oldParRun = UPstream::parRun();

    bool syncPar = true;

    if (nFinalProcs_ > Pstream::nProcs())
    {
        UPstream::parRun() = false;
        syncPar = false;
    }

    // Change the mesh (no inflation). Note: parallel comms allowed (if one mesh pr proc).
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, syncPar);

    // Update fields. No inflation, parallel sync (if one mesh per proc).
    mesh.updateMesh(map);
    UPstream::parRun() = oldParRun;

    // Adapt constructMaps for merged points.
    forAll(constructPointMap, meshi)
    {
        labelList& constructMap = constructPointMap[meshi];

        forAll(constructMap, i)
        {
            label oldPointi = constructMap[i];

            label newPointi = map().reversePointMap()[oldPointi];

            if (newPointi < -1)
            {
                constructMap[i] = -newPointi-2;
            }
            else if (newPointi >= 0)
            {
                constructMap[i] = newPointi;
            }
            else
            {
                FatalErrorInFunction
                    << "Problem. oldPointi:" << oldPointi
                    << " newPointi:" << newPointi << abort(FatalError);
            }
        }
    }
    return map;
}


void Foam::fvMeshDistribute::getCouplingData
(
    const labelList& distribution,
    labelList& sourceFace,
    labelList& sourceProc,
    labelList& sourcePatch,
    labelList& sourceNbrPatch,
    labelList& sourceNewNbrProc,
    labelList& sourcePointMaster
) const
{
    // Construct the coupling information for all (boundary) faces and
    // points

    const label nBnd = mesh_.nBoundaryFaces();
    sourceFace.setSize(nBnd);
    sourceProc.setSize(nBnd);
    sourcePatch.setSize(nBnd);
    sourceNbrPatch.setSize(nBnd);
    sourceNewNbrProc.setSize(nBnd);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get neighbouring meshFace labels and new processor of coupled boundaries.
    labelList nbrFaces(nBnd, -1);
    labelList nbrNewNbrProc(nBnd, -1);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            label offset = pp.start() - mesh_.nInternalFaces();

            // Mesh labels of faces on this side
            forAll(pp, i)
            {
                label bndI = offset + i;
                nbrFaces[bndI] = pp.start()+i;
            }

            // Which processor they will end up on
            SubList<label>(nbrNewNbrProc, pp.size(), offset) =
                UIndirectList<label>(distribution, pp.faceCells())();
        }
    }


    // Exchange the boundary data
    syncTools::swapBoundaryFaceList(mesh_, nbrFaces);
    syncTools::swapBoundaryFaceList(mesh_, nbrNewNbrProc);


    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        label offset = pp.start() - mesh_.nInternalFaces();

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Check which of the two faces we store.

            if (procPatch.owner())
            {
                // Use my local face labels
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = pp.start()+i;
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
            else
            {
                // Use my neighbours face labels
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = nbrFaces[bndI];
                    sourceProc[bndI] = procPatch.neighbProcNo();
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }


            label patchi = -1, nbrPatchi = -1;
            if (isA<processorCyclicPolyPatch>(pp))
            {
                patchi =
                    refCast<const processorCyclicPolyPatch>(pp)
                   .referPatchID();
                nbrPatchi =
                    refCast<const cyclicPolyPatch>(patches[patchi])
                   .nbrPatchID();
            }

            forAll(pp, i)
            {
                label bndI = offset + i;
                sourcePatch[bndI] = patchi;
                sourceNbrPatch[bndI] = nbrPatchi;
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cpp = refCast<const cyclicPolyPatch>(pp);

            if (cpp.owner())
            {
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = pp.start()+i;
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourcePatch[bndI] = patchi;
                    sourceNbrPatch[bndI] = cpp.nbrPatchID();
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
            else
            {
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = nbrFaces[bndI];
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourcePatch[bndI] = patchi;
                    sourceNbrPatch[bndI] = cpp.nbrPatchID();
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
        }
        else
        {
            // Normal (physical) boundary
            forAll(pp, i)
            {
                label bndI = offset + i;
                sourceFace[bndI] = -1;
                sourceProc[bndI] = -1;
                sourcePatch[bndI] = patchi;
                sourceNbrPatch[bndI] = -1;
                sourceNewNbrProc[bndI] = -1;
            }
        }
    }


    // Collect coupled (collocated) points
    sourcePointMaster.setSize(mesh_.nPoints());
    sourcePointMaster = -1;
    {
        // Assign global master point
        const globalIndex globalPoints(mesh_.nPoints());

        const globalMeshData& gmd = mesh_.globalData();
        const indirectPrimitivePatch& cpp = gmd.coupledPatch();
        const labelList& meshPoints = cpp.meshPoints();
        const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
        const labelListList& slaves = gmd.globalCoPointSlaves();

        labelList elems(slavesMap.constructSize(), -1);
        forAll(meshPoints, pointi)
        {
            const labelList& slots = slaves[pointi];

            if (slots.size())
            {
                // pointi is a master. Assign a unique label.

                label globalPointi = globalPoints.toGlobal(meshPoints[pointi]);
                elems[pointi] = globalPointi;
                forAll(slots, i)
                {
                    label sloti = slots[i];
                    if (sloti >= meshPoints.size())
                    {
                        // Filter out local collocated points. We don't want
                        // to merge these
                        elems[slots[i]] = globalPointi;
                    }
                }
            }
        }

        // Push slave-slot data back to slaves
        slavesMap.reverseDistribute(elems.size(), elems, false);

        // Extract back onto mesh
        forAll(meshPoints, pointi)
        {
            sourcePointMaster[meshPoints[pointi]] = elems[pointi];
        }
    }
}


// Subset the neighbourCell/neighbourProc fields
void Foam::fvMeshDistribute::subsetCouplingData
(
    const fvMesh& mesh,
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap,

    const labelList& oldDistribution,
    const labelList& oldFaceOwner,
    const labelList& oldFaceNeighbour,
    const label oldInternalFaces,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNbrPatch,
    const labelList& sourceNewNbrProc,
    const labelList& sourcePointMaster,

    labelList& subFace,
    labelList& subProc,
    labelList& subPatch,
    labelList& subNbrPatch,
    labelList& subNewNbrProc,
    labelList& subPointMaster
)
{
    subFace.setSize(mesh.nBoundaryFaces());
    subProc.setSize(mesh.nBoundaryFaces());
    subPatch.setSize(mesh.nBoundaryFaces());
    subNbrPatch.setSize(mesh.nBoundaryFaces());
    subNewNbrProc.setSize(mesh.nBoundaryFaces());

    forAll(subFace, newBFacei)
    {
        label newFacei = newBFacei + mesh.nInternalFaces();

        label oldFacei = faceMap[newFacei];

        // Was oldFacei internal face? If so which side did we get.
        if (oldFacei < oldInternalFaces)
        {
            subFace[newBFacei] = oldFacei;
            subProc[newBFacei] = Pstream::myProcNo();
            subPatch[newBFacei] = -1;

            label oldOwn = oldFaceOwner[oldFacei];
            label oldNei = oldFaceNeighbour[oldFacei];

            if (oldOwn == cellMap[mesh.faceOwner()[newFacei]])
            {
                // We kept the owner side. Where does the neighbour move to?
                subNewNbrProc[newBFacei] = oldDistribution[oldNei];
            }
            else
            {
                // We kept the neighbour side.
                subNewNbrProc[newBFacei] = oldDistribution[oldOwn];
            }
        }
        else
        {
            // Was boundary face. Take over boundary information
            label oldBFacei = oldFacei - oldInternalFaces;

            subFace[newBFacei] = sourceFace[oldBFacei];
            subProc[newBFacei] = sourceProc[oldBFacei];
            subPatch[newBFacei] = sourcePatch[oldBFacei];
            subNbrPatch[newBFacei] = sourceNbrPatch[oldBFacei];
            subNewNbrProc[newBFacei] = sourceNewNbrProc[oldBFacei];
        }
    }


    subPointMaster = UIndirectList<label>(sourcePointMaster, pointMap);
}


// Find cells on mesh whose faceID/procID match the neighbour cell/proc of
// domainMesh. Store the matching face.
void Foam::fvMeshDistribute::findCouples
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,

    const label domain,
    const primitiveMesh& domainMesh,
    const labelList& domainFace,
    const labelList& domainProc,
    const labelList& domainPatch,

    labelList& masterCoupledFaces,
    labelList& slaveCoupledFaces
)
{
    // Store domain neighbour as map so we can easily look for pair
    // with same face+proc.
    labelPairLookup map(domainFace.size());

    forAll(domainProc, bFacei)
    {
        if (domainProc[bFacei] != -1 && domainPatch[bFacei] == -1)
        {
            map.insert
            (
                labelPair(domainFace[bFacei], domainProc[bFacei]),
                bFacei
            );
        }
    }


    // Try to match mesh data.

    masterCoupledFaces.setSize(domainFace.size());
    slaveCoupledFaces.setSize(domainFace.size());
    label coupledI = 0;

    forAll(sourceFace, bFacei)
    {
        if (sourceProc[bFacei] != -1 && sourcePatch[bFacei] == -1)
        {
            labelPair myData(sourceFace[bFacei], sourceProc[bFacei]);

            labelPairLookup::const_iterator iter = map.find(myData);

            if (iter != map.end())
            {
                label nbrBFacei = iter();

                masterCoupledFaces[coupledI] = mesh.nInternalFaces() + bFacei;
                slaveCoupledFaces[coupledI] =
                    domainMesh.nInternalFaces()
                  + nbrBFacei;

                coupledI++;
            }
        }
    }

    masterCoupledFaces.setSize(coupledI);
    slaveCoupledFaces.setSize(coupledI);

    if (debug)
    {
        Pout<< "findCouples : found " << coupledI
            << " faces that will be stitched" << nl << endl;
    }
}


// Map data on boundary faces to new mesh (resulting from adding two meshes)
Foam::labelList Foam::fvMeshDistribute::mapBoundaryData
(
    const primitiveMesh& mesh,      // mesh after adding
    const mapAddedPolyMesh& map,
    const labelList& boundaryData0, // on mesh before adding
    const label nInternalFaces1,
    const labelList& boundaryData1  // on added mesh
)
{
    labelList newBoundaryData(mesh.nBoundaryFaces());

    forAll(boundaryData0, oldBFacei)
    {
        label newFacei = map.oldFaceMap()[oldBFacei + map.nOldInternalFaces()];

        // Face still exists (is necessary?) and still boundary face
        if (newFacei >= 0 && newFacei >= mesh.nInternalFaces())
        {
            newBoundaryData[newFacei - mesh.nInternalFaces()] =
                boundaryData0[oldBFacei];
        }
    }

    forAll(boundaryData1, addedBFacei)
    {
        label newFacei = map.addedFaceMap()[addedBFacei + nInternalFaces1];

        if (newFacei >= 0 && newFacei >= mesh.nInternalFaces())
        {
            newBoundaryData[newFacei - mesh.nInternalFaces()] =
                boundaryData1[addedBFacei];
        }
    }

    return newBoundaryData;
}


Foam::labelList Foam::fvMeshDistribute::mapPointData
(
    const primitiveMesh& mesh,      // mesh after adding
    const mapAddedPolyMesh& map,
    const labelList& boundaryData0, // on mesh before adding
    const labelList& boundaryData1  // on added mesh
)
{
    labelList newBoundaryData(mesh.nPoints());

    forAll(boundaryData0, oldPointi)
    {
        label newPointi = map.oldPointMap()[oldPointi];

        // Point still exists (is necessary?)
        if (newPointi >= 0)
        {
            newBoundaryData[newPointi] = boundaryData0[oldPointi];
        }
    }

    forAll(boundaryData1, addedPointi)
    {
        label newPointi = map.addedPointMap()[addedPointi];

        if (newPointi >= 0)
        {
            newBoundaryData[newPointi] = boundaryData1[addedPointi];
        }
    }

    return newBoundaryData;
}


// Remove cells. Add all exposed faces to patch oldInternalPatchi
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::doRemoveCells
(
    const labelList& cellsToRemove,
    const label oldInternalPatchi
)
{
    // Mesh change engine
    polyTopoChange meshMod(mesh_);

    // Cell removal topo engine. Do NOT synchronize parallel since
    // we are doing a local cell removal.
    removeCells cellRemover(mesh_, false);

    // Get all exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));

    // Insert the topo changes
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        labelList(exposedFaces.size(), oldInternalPatchi),  // patch for exposed
                                                            // faces.
        meshMod
    );


    //// Generate test field
    //tmp<surfaceScalarField> sfld(generateTestField(mesh_));

    // Save internal fields (note: not as DimensionedFields since would
    // get mapped)
    PtrList<Field<scalar>> sFlds;
    saveInternalFields(sFlds);
    PtrList<Field<vector>> vFlds;
    saveInternalFields(vFlds);
    PtrList<Field<sphericalTensor>> sptFlds;
    saveInternalFields(sptFlds);
    PtrList<Field<symmTensor>> sytFlds;
    saveInternalFields(sytFlds);
    PtrList<Field<tensor>> tFlds;
    saveInternalFields(tFlds);

    // Change the mesh. No inflation. Note: no parallel comms allowed.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, false);

    // Update fields
    mesh_.updateMesh(map);


    // Any exposed faces in a surfaceField will not be mapped. Map the value
    // of these separately (until there is support in all PatchFields for
    // mapping from internal faces ...)

    mapExposedFaces(map(), sFlds);
    mapExposedFaces(map(), vFlds);
    mapExposedFaces(map(), sptFlds);
    mapExposedFaces(map(), sytFlds);
    mapExposedFaces(map(), tFlds);


    //// Test test field
    //testField(sfld);


    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }


    return map;
}


// Delete and add processor patches. Changes mesh and returns per neighbour proc
// the processor patchID.
void Foam::fvMeshDistribute::addProcPatches
(
    const labelList& nbrProc,         // Processor that neighbour is now on
    const labelList& referPatchID,    // Original patch ID (or -1)
    const labelList& referNbrPatchID, // Original neighbour patch ID (or -1)
    const label myMeshNo,
    fvMesh& mesh,
    List<Map<label>>& procPatchID
)
{
    // Now use the neighbourFace/Proc to repatch the mesh. These lists
    // contain for all current boundary faces the global patchID (for non-proc
    // patch) or the processor.

    // Determine a visit order such that the processor patches get added
    // in order of increasing neighbour processor (and for same neighbour
    // processor (in case of processor cyclics) in order of increasing
    // 'refer' patch)
    labelList indices;
    sortedOrder
    (
        nbrProc,
        indices,
        lessProcPatches(nbrProc, referPatchID, referNbrPatchID)
    );

    procPatchID.setSize(nFinalProcs_);

    forAll(indices, i)
    {
        label bFacei = indices[i];
        label proci = nbrProc[bFacei];

        if (proci != -1 && proci != myMeshNo)
        {
            if (!procPatchID[proci].found(referPatchID[bFacei]))
            {
                // No patch for neighbour yet. Is either a normal processor
                // patch or a processorCyclic patch.

                if (referPatchID[bFacei] == -1)
                {
                    // Ordinary processor boundary

                    processorPolyPatch pp
                    (
                        0,              // size
                        mesh.nFaces(),
                        mesh.boundaryMesh().size(),
                        mesh.boundaryMesh(),
                        myMeshNo,
                        proci
                    );

                    procPatchID[proci].insert
                    (
                        referPatchID[bFacei],
                        fvMeshTools::addPatch
                        (
                            mesh,
                            pp,
                            dictionary(),   // optional per field patchField
                            processorFvPatchField<scalar>::typeName,
                            false           // not parallel sync
                        )
                    );
                }
                else
                {
                    const coupledPolyPatch& pcPatch =
                        refCast<const coupledPolyPatch>
                        (
                            mesh.boundaryMesh()[referPatchID[bFacei]]
                        );
                    processorCyclicPolyPatch pp
                    (
                        0,              // size
                        mesh.nFaces(),
                        mesh.boundaryMesh().size(),
                        mesh.boundaryMesh(),
                        myMeshNo,
                        proci,
                        pcPatch.name()
                    );

                    procPatchID[proci].insert
                    (
                        referPatchID[bFacei],
                        fvMeshTools::addPatch
                        (
                            mesh,
                            pp,
                            dictionary(),   // optional per field patchField
                            processorCyclicFvPatchField<scalar>::typeName,
                            false           // not parallel sync
                        )
                    );
                }
            }
        }
    }
} // addProcPatches


void Foam::fvMeshDistribute::addNccProcPatches
(
    PtrList<fvMesh>& myProcMeshes,
    const List<labelList>& meshNoInProc
)
{
    List<List<boolList>> finalProcHasOrig(Pstream::nProcs());
    List<List<boolList>> finalProcHasNbrOrig(Pstream::nProcs());

    finalProcHasOrig[Pstream::myProcNo()].setSize(myProcMeshes.size());
    finalProcHasNbrOrig[Pstream::myProcNo()].setSize(myProcMeshes.size());

    DynamicList<label> ownNccPatches;

    forAll(myProcMeshes, meshi)
    {
        fvMesh& procMesh = myProcMeshes[meshi];

        finalProcHasOrig[Pstream::myProcNo()][meshi] =
            boolList(procMesh.boundaryMesh().size(), false);

        finalProcHasNbrOrig[Pstream::myProcNo()][meshi] =
            boolList(procMesh.boundaryMesh().size(), false);

        forAll(procMesh.boundaryMesh(), nccPatchi)
        {
            const polyPatch& pp = procMesh.boundaryMesh()[nccPatchi];

            if (!isA<nonConformalCyclicPolyPatch>(pp)) continue;

            const nonConformalCyclicPolyPatch& nccPp =
                refCast<const nonConformalCyclicPolyPatch>(pp);
            const polyPatch& origPp = nccPp.origPatch();

            const nonConformalCyclicPolyPatch& nbrNccPp = nccPp.nbrPatch();
            const polyPatch& nbrOrigPp = nbrNccPp.origPatch();

            if (!nccPp.owner()) continue;

            // Append patch ID to the list of owner NCC patches
            if (Pstream::master() && meshi == 0)
            {
                ownNccPatches.append(nccPatchi);
            }

            finalProcHasOrig[Pstream::myProcNo()][meshi][nccPatchi] =
                !origPp.empty();

            finalProcHasNbrOrig[Pstream::myProcNo()][meshi][nccPatchi] =
                !nbrOrigPp.empty();
        }
    }

    Pstream::allGatherList(finalProcHasOrig);
    Pstream::allGatherList(finalProcHasNbrOrig);
    Pstream::scatter(ownNccPatches);

    // Ensure proper processor communication by ordering the patches so that
    // for the lower indexed processor, the owner interface comes first;
    // and for the higher indexed processor, the neighbour comes first.
    auto add = [&](const bool owner, const bool first)
    {
        forAll(ownNccPatches, nccPatchi)
        {
            const label nccPatchID = ownNccPatches[nccPatchi];

            List<boolList> finalProcHasPatchA(Pstream::nProcs());
            List<boolList> finalProcHasPatchB(Pstream::nProcs());

            finalProcHasPatchA[Pstream::myProcNo()].setSize
            (
                myProcMeshes.size()
            );

            finalProcHasPatchB[Pstream::myProcNo()].setSize
            (
                myProcMeshes.size()
            );

            forAll(myProcMeshes, meshi)
            {
                finalProcHasPatchA[Pstream::myProcNo()][meshi] =
                    owner
                  ? finalProcHasOrig[Pstream::myProcNo()][meshi][nccPatchID]
                  : finalProcHasNbrOrig[Pstream::myProcNo()][meshi][nccPatchID];

                finalProcHasPatchB[Pstream::myProcNo()][meshi] =
                    owner
                  ? finalProcHasNbrOrig[Pstream::myProcNo()][meshi][nccPatchID]
                  : finalProcHasOrig[Pstream::myProcNo()][meshi][nccPatchID];
            }

            Pstream::allGatherList(finalProcHasPatchA);
            Pstream::allGatherList(finalProcHasPatchB);

            forAll(myProcMeshes, meshi)
            {
                fvMesh& procMesh = myProcMeshes[meshi];

                const polyPatch& pp = procMesh.boundaryMesh()[nccPatchID];

                const nonConformalCyclicPolyPatch& nccPp =
                    refCast<const nonConformalCyclicPolyPatch>(pp);
                const polyPatch& origPp = nccPp.origPatch();

                const nonConformalCyclicPolyPatch& nbrNccPp = nccPp.nbrPatch();
                const polyPatch& nbrOrigPp = nbrNccPp.origPatch();

                const label nProcMeshA =
                    meshNoInProc[Pstream::myProcNo()][meshi];

                if (finalProcHasPatchA[Pstream::myProcNo()][meshi])
                {
                    forAll(finalProcHasPatchB, proci)
                    {
                        forAll(finalProcHasPatchB[proci], procMeshi)
                        {
                            const label nProcMeshB =
                                meshNoInProc[proci][procMeshi];

                            if
                            (
                                (
                                    (first && nProcMeshB > nProcMeshA)
                                 || (!first && nProcMeshB < nProcMeshA)
                                )
                             && finalProcHasPatchB[proci][procMeshi]
                            )
                            {
                                fvMeshTools::addPatch
                                (
                                    procMesh,
                                    nonConformalProcessorCyclicPolyPatch
                                    (
                                        0,
                                        procMesh.nFaces(),
                                        procMesh.boundaryMesh().size(),
                                        procMesh.boundaryMesh(),
                                        nProcMeshA,
                                        nProcMeshB,
                                        owner ? nccPp.name() : nbrNccPp.name(),
                                        owner ? origPp.name() : nbrOrigPp.name()
                                    ),
                                    dictionary(),
                                    nonConformalProcessorCyclicFvPatchField
                                    <scalar>::typeName,
                                    false    // not parallel synced
                                );
                            }
                        }
                    }
                }
            }
        }
    };

    add(true, true);
    add(false, true);
    add(false, false);
    add(true, false);
}


Foam::List<Foam::List<Foam::DynamicList<Foam::label>>>
Foam::fvMeshDistribute::procPolyFaceAddr
(
    const List<labelList>& subMap,
    const surfaceLabelField::Boundary origPolyFacesBf,
    const labelList& origFaceOwner,
    const labelList& cellToProc,
    const label nccPatchID,
    const label nccNbrPatchID
) const
{
    List<List<DynamicList<label>>> finalProcOrigPolyFaces(Pstream::nProcs());
    finalProcOrigPolyFaces[Pstream::myProcNo()].setSize(nFinalProcs_);

    // Loop through all finite-volume faces from original meshes for a given
    // non-conformal patch and add them to the addressing, if both owner and
    // neighbour face pair are on the same redistributed processor mesh.
    forAll(origPolyFacesBf[nccPatchID], polyFacei)
    {
        const label facei = origPolyFacesBf[nccPatchID][polyFacei];
        const label celli = origFaceOwner[facei];
        const label proci = cellToProc[celli];

        const label nbrFacei = origPolyFacesBf[nccNbrPatchID][polyFacei];
        const label nbrCelli = origFaceOwner[nbrFacei];
        const label nbrProci = cellToProc[nbrCelli];

        if (proci != nbrProci) continue;

        // Find face label on the distribution map from subsetted data back
        // to the original mesh. This index will be mapped later to the
        // constructed (redistributed) mesh.
        const label fIndex =
            subMap[proci].find(origPolyFacesBf[nccPatchID][polyFacei] + 1);

        if (fIndex == -1)
        {
            FatalErrorInFunction
                << "Poly-face corresponding to the non-conformal face "
                << polyFacei << " (from patchID " << nccPatchID
                << ") not found on the respective distribution map from the"
                << " original mesh." << exit(FatalError);
        }

        finalProcOrigPolyFaces[Pstream::myProcNo()][proci].append(fIndex);
    }

    return finalProcOrigPolyFaces;
}


void Foam::fvMeshDistribute::unconform
(
    PtrList<fvMesh>& myProcMeshes,
    const List<labelList>& meshNoInProc,
    const PtrList<mapDistributePolyMesh>& distMaps,
    const surfaceLabelField::Boundary origPolyFacesBf,
    const labelList& origFaceOwner,
    const labelList& cellToProc
)
{
    Info<< "\nUnconforming processor meshes" << endl;

    // Construct labels of finite-volume boundary faces on redistributed
    // processor meshes. Only needed for non-conformal cyclic patches.
    List<List<List<DynamicList<label>>>> polyFacesBAddr(Pstream::nProcs());
    polyFacesBAddr[Pstream::myProcNo()].setSize(myProcMeshes.size());
    Pstream::gatherList(polyFacesBAddr);

    forAll(origPolyFacesBf, patchi)
    {
        // Stop when a processor-like patch field is found
        if
        (
            origPolyFacesBf[patchi].patchType() == "processor"
         || origPolyFacesBf[patchi].patchType() == "processorCyclic"
         || origPolyFacesBf[patchi].patchType() == "nonConformalProcessorCyclic"
        ) break;

        const fvPatch& fvp = myProcMeshes[0].boundary()[patchi];

        if (!isA<nonConformalCyclicFvPatch>(fvp)) continue;

        const nonConformalCyclicFvPatch& nccFvp =
            refCast<const nonConformalCyclicFvPatch>(fvp);

        const label nbrPatchi = nccFvp.nbrPatchID();

        // Construct list of poly-face labels for this non-conformal cyclic
        // patch across all redistributed processor meshes, mapped from the
        // distribution maps from the original meshes.
        List<List<DynamicList<label>>> finalProcOrigPolyFaces =
            procPolyFaceAddr
            (
                distMaps[0].faceMap().subMap(),
                origPolyFacesBf,
                origFaceOwner,
                cellToProc,
                patchi,
                nbrPatchi
            );

        Pstream::gatherList(finalProcOrigPolyFaces);

        // Distribution maps from subsetted data to the redistributed meshes
        List<List<List<labelList>>> constructMaps(Pstream::nProcs());
        constructMaps[Pstream::myProcNo()].setSize(myProcMeshes.size());

        forAll(myProcMeshes, procMeshi)
        {
            constructMaps[Pstream::myProcNo()][procMeshi] =
                distMaps[procMeshi].faceMap().constructMap();
        }

        Pstream::gatherList(constructMaps);

        // Map the poly-face labels from their location on the final meshes
        // to the process number mesh distribution
        if (Pstream::master())
        {
            forAll(finalProcOrigPolyFaces, proci)
            {
                forAll(finalProcOrigPolyFaces[proci], destProci)
                {
                    const labelList& patchPolyFaces =
                        finalProcOrigPolyFaces[proci][destProci];

                    if (patchPolyFaces.size())
                    {
                        label p = -1, m = -1;
                        forAll(meshNoInProc, processori)
                        {
                            m = meshNoInProc[processori].find(destProci);
                            if (m != -1)
                            {
                                p = processori;
                                break;
                            }
                        }

                        if (!polyFacesBAddr[p][m].size())
                        {
                            polyFacesBAddr[p][m].setSize
                            (
                                origPolyFacesBf.size()
                            );
                        }

                        polyFacesBAddr[p][m][patchi].append
                        (
                            labelField
                            (
                                constructMaps[p][m][proci],
                                patchPolyFaces
                            ) - 1
                        );
                    }
                }
            }
        }
    }

    Pstream::scatter(polyFacesBAddr);

    // Distribute the poly-face boundary addressing from the original processors
    // and unconform all processor meshes. Set dummy data for the face geometry,
    // since it should not be used during redistribution.
    forAll(myProcMeshes, procMeshi)
    {
        fvMesh& procMesh = myProcMeshes[procMeshi];

        surfaceLabelField::Boundary polyFacesBf
        (
            surfaceLabelField::null(),
            procMesh.polyFacesBf()
        );
        surfaceVectorField Sf(procMesh.Sf().cloneUnSliced());
        surfaceVectorField Cf(procMesh.Cf().cloneUnSliced());

        // Insert finite-volume faces into non-conformal cyclics, for the
        // redistributed processor meshes with non-conformal faces.
        if (polyFacesBAddr[Pstream::myProcNo()][procMeshi].size())
        {
            forAll(polyFacesBAddr[Pstream::myProcNo()][procMeshi], patchi)
            {
                const labelList& pPolyFaces =
                    polyFacesBAddr[Pstream::myProcNo()][procMeshi][patchi];

                if (!pPolyFaces.size()) continue;

                const fvPatch& fvp = procMesh.boundary()[patchi];
                if (!isA<nonConformalCyclicFvPatch>(fvp))
                {
                    FatalErrorInFunction
                        << "Wrong patch ordering for the redistributed meshes."
                        << exit(FatalError);
                }

                polyFacesBf[patchi] = pPolyFaces;

                const label size = polyFacesBf[patchi].size();

                Sf.boundaryFieldRef()[patchi].resize(size, Zero);
                Cf.boundaryFieldRef()[patchi].resize(size, Zero);
            }
        }

        // Unset parallel flag to avoid problems with processor boundaries
        // communication when constructing the (dummy) mesh geometry, during
        // the unconforming operation.
        const bool oldParRun = Pstream::parRun();
        Pstream::parRun() = false;

        procMesh.unconform
        (
            polyFacesBf,
            Sf,
            Cf,
            NullObjectRef<surfaceScalarField>(),
            false
        );

        Pstream::parRun() = oldParRun;

        // Sanity check for matching faces between owner and neighbour
        // non-conformal cyclic patches
        forAll(procMesh.boundary(), nccPatchi)
        {
            const fvPatch& fvp = procMesh.boundary()[nccPatchi];
            if (!isA<nonConformalCyclicFvPatch>(fvp)) continue;

            const nonConformalCyclicFvPatch& nccFvp =
                refCast<const nonConformalCyclicFvPatch>(fvp);

            if (!nccFvp.owner()) continue;

            if (nccFvp.size() != nccFvp.nbrPatch().size())
            {
                FatalErrorInFunction
                    << "The number of non-conformal faces for patch "
                    << nccFvp.name() << " does not match its neighbour patch "
                    << nccFvp.nbrPatch().name() << exit(FatalError);
            }
        }
    }
}


// Get boundary faces to be repatched. Is -1 or new patchID
Foam::labelList Foam::fvMeshDistribute::getBoundaryPatch
(
    const labelList& nbrProc,               // new processor per boundary face
    const labelList& referPatchID,          // patchID (or -1) I originated from
    const List<Map<label>>& procPatchID,    // per proc the new procPatches
    const label myMeshNo
)
{
    labelList patchIDs(nbrProc);

    forAll(nbrProc, bFacei)
    {
        if (nbrProc[bFacei] == myMeshNo)
        {
            label origPatchi = referPatchID[bFacei];
            patchIDs[bFacei] = origPatchi;
        }
        else if (nbrProc[bFacei] != -1)
        {
            label origPatchi = referPatchID[bFacei];
            patchIDs[bFacei] = procPatchID[nbrProc[bFacei]][origPatchi];
        }
        else
        {
            patchIDs[bFacei] = -1;
        }
    }
    return patchIDs;
}


// Send mesh and coupling data.
void Foam::fvMeshDistribute::sendMesh
(
    const label domain,
    const fvMesh& mesh,

    const wordList& pointZoneNames,
    const wordList& faceZoneNames,
    const wordList& cellZoneNames,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNbrPatch,
    const labelList& sourceNewNbrProc,
    const labelList& sourcePointMaster,
    Ostream& toDomain
)
{
    if (debug)
    {
        Pout<< "Sending to domain " << domain << nl
            << "    nPoints:" << mesh.nPoints() << nl
            << "    nFaces:" << mesh.nFaces() << nl
            << "    nCells:" << mesh.nCells() << nl
            << "    nPatches:" << mesh.boundaryMesh().size() << nl
            << endl;
    }

    // Assume sparse, possibly overlapping point zones. Get contents
    // in merged-zone indices.
    CompactListList<label> zonePoints;
    {
        const pointZoneMesh& pointZones = mesh.pointZones();

        labelList rowSizes(pointZoneNames.size(), 0);

        forAll(pointZoneNames, nameI)
        {
            label myZoneID = pointZones.findZoneID(pointZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = pointZones[myZoneID].size();
            }
        }
        zonePoints.setSize(rowSizes);

        forAll(pointZoneNames, nameI)
        {
            label myZoneID = pointZones.findZoneID(pointZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zonePoints[nameI].deepCopy(pointZones[myZoneID]);
            }
        }
    }

    // Assume sparse, possibly overlapping face zones
    CompactListList<label> zoneFaces;
    CompactListList<bool> zoneFaceFlip;
    {
        const faceZoneMesh& faceZones = mesh.faceZones();

        labelList rowSizes(faceZoneNames.size(), 0);

        forAll(faceZoneNames, nameI)
        {
            label myZoneID = faceZones.findZoneID(faceZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = faceZones[myZoneID].size();
            }
        }

        zoneFaces.setSize(rowSizes);
        zoneFaceFlip.setSize(rowSizes);

        forAll(faceZoneNames, nameI)
        {
            label myZoneID = faceZones.findZoneID(faceZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zoneFaces[nameI].deepCopy(faceZones[myZoneID]);
                zoneFaceFlip[nameI].deepCopy(faceZones[myZoneID].flipMap());
            }
        }
    }

    // Assume sparse, possibly overlapping cell zones
    CompactListList<label> zoneCells;
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        labelList rowSizes(cellZoneNames.size(), 0);

        forAll(cellZoneNames, nameI)
        {
            label myZoneID = cellZones.findZoneID(cellZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = cellZones[myZoneID].size();
            }
        }

        zoneCells.setSize(rowSizes);

        forAll(cellZoneNames, nameI)
        {
            label myZoneID = cellZones.findZoneID(cellZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zoneCells[nameI].deepCopy(cellZones[myZoneID]);
            }
        }
    }
    ////- Assume full cell zones
    //labelList cellZoneID;
    //if (hasCellZones)
    //{
    //    cellZoneID.setSize(mesh.nCells());
    //    cellZoneID = -1;
    //
    //    const cellZoneMesh& cellZones = mesh.cellZones();
    //
    //    forAll(cellZones, zoneI)
    //    {
    //        labelUIndList(cellZoneID, cellZones[zoneI]) = zoneI;
    //    }
    //}

    // Send
    toDomain
        << mesh.points()
        << CompactListList<label, face>(mesh.faces())
        << mesh.faceOwner()
        << mesh.faceNeighbour()
        << mesh.boundaryMesh()

        << zonePoints
        << zoneFaces
        << zoneFaceFlip
        << zoneCells

        << sourceFace
        << sourceProc
        << sourcePatch
        << sourceNbrPatch
        << sourceNewNbrProc
        << sourcePointMaster;


    if (debug)
    {
        Pout<< "Started sending mesh to domain " << domain
            << endl;
    }
}


// Receive mesh. Opposite of sendMesh
Foam::autoPtr<Foam::fvMesh> Foam::fvMeshDistribute::receiveMesh
(
    const label domain,
    const wordList& pointZoneNames,
    const wordList& faceZoneNames,
    const wordList& cellZoneNames,
    const Time& runTime,
    labelList& domainSourceFace,
    labelList& domainSourceProc,
    labelList& domainSourcePatch,
    labelList& domainSourceNbrPatch,
    labelList& domainSourceNewNbrProc,
    labelList& domainSourcePointMaster,
    Istream& fromNbr
)
{
    pointField domainPoints(fromNbr);
    faceList domainFaces = CompactListList<label, face>(fromNbr)();
    labelList domainAllOwner(fromNbr);
    labelList domainAllNeighbour(fromNbr);
    PtrList<entry> patchEntries(fromNbr);

    CompactListList<label> zonePoints(fromNbr);
    CompactListList<label> zoneFaces(fromNbr);
    CompactListList<bool> zoneFaceFlip(fromNbr);
    CompactListList<label> zoneCells(fromNbr);

    fromNbr
        >> domainSourceFace
        >> domainSourceProc
        >> domainSourcePatch
        >> domainSourceNbrPatch
        >> domainSourceNewNbrProc
        >> domainSourcePointMaster;

    // Construct fvMesh
    autoPtr<fvMesh> domainMeshPtr
    (
        new fvMesh
        (
            IOobject
            (
                fvMesh::defaultRegion,
                runTime.timeName(),
                runTime,
                IOobject::NO_READ
            ),
            xferMove(domainPoints),
            xferMove(domainFaces),
            xferMove(domainAllOwner),
            xferMove(domainAllNeighbour),
            false                   // no parallel comms
        )
    );
    fvMesh& domainMesh = domainMeshPtr();

    List<polyPatch*> patches(patchEntries.size());

    forAll(patchEntries, patchi)
    {
        patches[patchi] = polyPatch::New
        (
            patchEntries[patchi].keyword(),
            patchEntries[patchi].dict(),
            patchi,
            domainMesh.boundaryMesh()
        ).ptr();
    }
    // Add patches; no parallel comms
    domainMesh.addFvPatches(patches, false);

    // Construct zones
    List<pointZone*> pZonePtrs(pointZoneNames.size());
    forAll(pZonePtrs, i)
    {
        pZonePtrs[i] = new pointZone
        (
            pointZoneNames[i],
            zonePoints[i],
            i,
            domainMesh.pointZones()
        );
    }

    List<faceZone*> fZonePtrs(faceZoneNames.size());
    forAll(fZonePtrs, i)
    {
        const boolList fzfI(zoneFaceFlip[i]);
        fZonePtrs[i] = new faceZone
        (
            faceZoneNames[i],
            zoneFaces[i],
            fzfI,
            i,
            domainMesh.faceZones()
        );
    }

    List<cellZone*> cZonePtrs(cellZoneNames.size());
    forAll(cZonePtrs, i)
    {
        cZonePtrs[i] = new cellZone
        (
            cellZoneNames[i],
            zoneCells[i],
            i,
            domainMesh.cellZones()
        );
    }
    domainMesh.addZones(pZonePtrs, fZonePtrs, cZonePtrs);

    return domainMeshPtr;
}


Foam::label Foam::fvMeshDistribute::getRotation
(
    const pointField& points,
    const face& f,
    const point& anchor,
    const scalar tol
)
{
    label anchorFp = -1;
    scalar minDistSqr = GREAT;

    forAll(f, fp)
    {
        scalar distSqr = magSqr(anchor - points[f[fp]]);

        if (distSqr < minDistSqr)
        {
            minDistSqr = distSqr;
            anchorFp = fp;
        }
    }

    if (anchorFp == -1 || Foam::sqrt(minDistSqr) > tol)
    {
        return -1;
    }
    else
    {
        // Check that anchor is unique
        forAll(f, fp)
        {
            scalar distSqr = magSqr(anchor - points[f[fp]]);

            if (distSqr == minDistSqr && fp != anchorFp)
            {
                WarningInFunction
                    << "Cannot determine unique anchor point on face "
                    << UIndirectList<point>(points, f) << endl
                    << "Both at index " << anchorFp << " and " << fp
                    << " the vertices have the same distance "
                    << Foam::sqrt(minDistSqr) << " to the anchor " << anchor
                    << ". Continuing but results might be wrong."
                    << nl << endl;
            }
        }

        // Positive rotation
        return (f.size() - anchorFp) % f.size();
    }
}


void Foam::fvMeshDistribute::reorderCoupledFaces
(
    PtrList<fvMesh>& myProcMeshes,
    const labelList& meshToProc,
    const labelListList& meshNoInProc
)
{
    // Allocate buffers (in worst case scenario all meshes in myProc are all neighbours)
    label bufSize = 0;
    forAll(meshNoInProc[Pstream::myProcNo()], i)
    {
        bufSize += i;
    }

    std::vector<PstreamBuffers> pBufs
    (
        nFinalProcs_,
        PstreamBuffers(Pstream::commsTypes::nonBlocking)
    );

    // Allocate internal streams
    std::vector<OStringStream> toLocalNeighbour
    (
        bufSize,
        OStringStream(IOstream::BINARY)
    );
    label nStr = 0;

    // Send ordering
    forAll(meshToProc, recvMesh)
    {
        forAll(myProcMeshes, i)
        {
            if (recvMesh <= meshNoInProc[Pstream::myProcNo()][i]) continue;

            fvMesh& mesh = myProcMeshes[i];

            const polyBoundaryMesh& boundary = mesh.boundaryMesh();

            forAll(boundary, patchi)
            {
                if (!isA<processorPolyPatch>(boundary[patchi])) continue;

                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(boundary[patchi]);

                if (procPatch.neighbProcNo() != recvMesh) continue;

                if (procPatch.owner())
                {
                    //initOrder
                    primitivePatch pp
                    (
                        SubList<face>
                        (
                            mesh.faces(),
                            procPatch.size(),
                            procPatch.start()
                        ),
                        mesh.points()
                    );

                    const pointField& ppPoints = pp.points();

                    pointField anchors(pp.size());

                    // Anchor point is the first point
                    forAll(pp, facei)
                    {
                        anchors[facei] = ppPoints[pp[facei][0]];
                    }

                    // Get the average of the points of each face. This is needed in
                    // case the face centroid calculation is incorrect due to the face
                    // having a very high aspect ratio.
                    pointField facePointAverages(pp.size(), Zero);
                    forAll(pp, fI)
                    {
                        const labelList& facePoints = pp[fI];

                        forAll(facePoints, pI)
                        {
                            facePointAverages[fI] += ppPoints[facePoints[pI]];
                        }

                        facePointAverages[fI] /= facePoints.size();
                    }

                    // findout if neighbour mesh is on this processor
                    const label recvProc = meshToProc[procPatch.neighbProcNo()];

                    if (debug)
                        Pout<< "Sending " << procPatch.name();

                    // Now send all info over to the neighbour
                    if (recvProc != Pstream::myProcNo())
                    {
                        if (debug)
                        {
                            Pout<< " to processor (node) " << recvProc
                                << endl;
                        }

                        UOPstream toNeighbour(recvProc, pBufs[procPatch.neighbProcNo()]);
                        toNeighbour << pp.faceCentres() << pp.faceNormals()
                                    << anchors << facePointAverages;
                    }
                    else
                    {
                        if (debug)
                        {
                            Pout<< " to internal stream"
                                << " (nStr " << nStr << ")" << endl;
                        }

                        toLocalNeighbour[nStr] << pp.faceCentres() << pp.faceNormals()
                                    << anchors << facePointAverages;
                        nStr++;
                    }
                } // if proc owner
            } // forAll procBoundaries
        } // forAll myProcMeshes
    } // forAll receiving meshes

    // Start sending&receiving from buffers
    for (auto i = 0; i < nFinalProcs_; ++i)
    {
        pBufs[i].finishedSends();
    }

    nStr = 0;

    // Receiving ordering
    forAll(myProcMeshes, i)
    {
        fvMesh& mesh = myProcMeshes[i];

        // Mapping for faces (old to new). Extends over all mesh faces for
        // convenience (could be just the external faces)
        labelList oldToNew(identity(mesh.faces().size()));

        // Rotation on new faces.
        labelList rotation(mesh.faces().size(), 0);

        const polyBoundaryMesh& boundary = mesh.boundaryMesh();

        labelList patchSizes(boundary.size());
        labelList patchStarts(boundary.size());

        bool anyChanged = false;

        forAll(boundary, patchi)
        {
            patchSizes[patchi]  = boundary[patchi].size();
            patchStarts[patchi] = boundary[patchi].start();

            if (!isA<processorPolyPatch>(boundary[patchi])) continue;

            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(boundary[patchi]);

            bool changed = false;

            //order
            primitivePatch pp
            (
                SubList<face>
                (
                    mesh.faces(),
                    procPatch.size(),
                    procPatch.start()
                ),
                mesh.points()
            );

            labelList faceMap(pp.size(), -1);
            labelList faceRotation(pp.size(), 0);

            if (procPatch.owner())
            {
                forAll(faceMap, patchFacei)
                {
                    faceMap[patchFacei] = patchFacei;
                }

                const pointField& ppPoints = pp.points();

                pointField anchors(pp.size());

                // Anchor point is the first point
                forAll(pp, facei)
                {
                    anchors[facei] = ppPoints[pp[facei][0]];
                }

                // Calculate typical distance from face centre
                scalarField tols
                (
                    procPatch.matchTolerance()
                   *procPatch.calcFaceTol(pp, pp.points(), pp.faceCentres())
                );

                forAll(faceMap, patchFacei)
                {
                    const point& wantedAnchor = anchors[patchFacei];

                    faceRotation[patchFacei] =
                        getRotation
                        (
                            ppPoints,
                            pp[patchFacei],
                            wantedAnchor,
                            tols[patchFacei]
                        );

                    if (faceRotation[patchFacei] > 0)
                    {
                        changed = true;
                    }
                }
            } // if procPatch.owner
            else
            {
                vectorField masterCtrs;
                vectorField masterNormals;
                vectorField masterAnchors;
                vectorField masterFacePointAverages;

                // findout if neighbour mesh is on this processor
                const label sendProc = meshToProc[procPatch.neighbProcNo()];

                if (debug) Pout<< "Receiving " << procPatch.name();

                // Receive data from neighbour
                if (sendProc != Pstream::myProcNo())
                {
                    if (debug)
                    {
                        Pout<< " from processor (node) " << sendProc
                            << endl;
                    }

                    UIPstream fromNeighbour(sendProc, pBufs[procPatch.myProcNo()]);
                    fromNeighbour >> masterCtrs >> masterNormals
                                  >> masterAnchors >> masterFacePointAverages;
                }
                else
                {
                    if (debug)
                    {
                        Pout<< " from internal stream"
                            << " (nStr " << nStr << ")" << endl;
                    }

                    IStringStream fromNeighbour
                    (
                        toLocalNeighbour[nStr].str(),
                        IOstream::BINARY
                    );

                    fromNeighbour >> masterCtrs >> masterNormals
                                  >> masterAnchors >> masterFacePointAverages;

                    nStr++;
                }

                // Calculate typical distance from face centre
                scalarField tols
                (
                   procPatch.matchTolerance()
                  *procPatch.calcFaceTol(pp, pp.points(), pp.faceCentres())
                );

                if (masterCtrs.size() != pp.size())
                {
                    FatalErrorInFunction
                        << "in patch:" << procPatch.name() << " : "
                        << "Local size of patch is " << pp.size() << " (faces)."
                        << endl
                        << "Received from neighbour " << masterCtrs.size()
                        << " faceCentres!"
                        << abort(FatalError);
                }

                // Geometric match of face centre vectors
                // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                // 1. Try existing ordering and transformation
                bool matchedAll =
                    matchPoints
                    (
                       pp.faceCentres(),
                       masterCtrs,
                       pp.faceNormals(),
                       masterNormals,
                       tols,
                       false,
                       faceMap
                    );

                //// Fallback: try using face point average for matching
                //if (!matchedAll)
                //{
                //    const pointField& ppPoints = pp.points();

                //    pointField facePointAverages(pp.size(), Zero);
                //    forAll(pp, fI)
                //    {
                //        const labelList& facePoints = pp[fI];

                //        forAll(facePoints, pI)
                //        {
                //            facePointAverages[fI] += ppPoints[facePoints[pI]];
                //        }

                //        facePointAverages[fI] /= facePoints.size();
                //    }

                //    scalarField tols2
                //    (
                //        procPatch.matchTolerance()
                //       *procPatch.calcFaceTol(pp, pp.points(), facePointAverages)
                //    );

                //    // Note that we do not use the faceNormals anymore for
                //    // comparison. Since we're
                //    // having problems with the face centres (e.g. due to extreme
                //    // aspect ratios) we will probably also have problems with
                //    // reliable normals calculation
                //    labelList faceMap2(faceMap.size(), -1);
                //    matchedAll = matchPoints
                //    (
                //        facePointAverages,
                //        masterFacePointAverages,
                //        tols2,
                //        true,
                //        faceMap2
                //    );

                //    forAll(faceMap, oldFacei)
                //    {
                //        if (faceMap[oldFacei] == -1)
                //        {
                //            faceMap[oldFacei] = faceMap2[oldFacei];
                //        }
                //    }

                //    matchedAll = true;

                //    forAll(faceMap, oldFacei)
                //    {
                //        if (faceMap[oldFacei] == -1)
                //        {
                //            matchedAll = false;
                //        }
                //    }
                //} // if !matchedAll

                if (!matchedAll)
                {
                    SeriousErrorInFunction
                        << "in patch:" << procPatch.name() << " : "
                        << "Cannot match vectors to faces on both sides of patch"
                        << endl
                        << "    masterCtrs[0]:" << masterCtrs[0] << endl
                        << "    ctrs[0]:" << pp.faceCentres()[0] << endl
                        << "    Check your topology changes or maybe you have"
                        << " multiple separated (from cyclics) processor patches"
                        << endl
                        << "    Continuing with incorrect face ordering from now on"
                        << endl;

                    changed = false;
                }

                // Set rotation.
                forAll(faceMap, oldFacei)
                {
                    // The face f will be at newFacei (after morphing) and we want
                    // its anchorPoint (= f[0]) to align with the anchorpoint for
                    // the corresponding face on the other side.

                    label newFacei = faceMap[oldFacei];

                    const point& wantedAnchor = masterAnchors[newFacei];

                    faceRotation[newFacei] =
                        getRotation
                        (
                            pp.points(),
                            pp[oldFacei],
                            wantedAnchor,
                            tols[oldFacei]
                        );

                    if (faceRotation[newFacei] == -1)
                    {
                        SeriousErrorInFunction
                            << "in patch " << procPatch.name()
                            << " : "
                            << "Cannot find point on face " << pp[oldFacei]
                            << " with vertices "
                            << UIndirectList<point>(pp.points(), pp[oldFacei])()
                            << " that matches point " << wantedAnchor
                            << " when matching the halves of processor patch "
                            << procPatch.name()
                            << "Continuing with incorrect face ordering from now on"
                            << endl;

                        changed = false;
                    }
                }

                forAll(faceMap, facei)
                {
                    if (faceMap[facei] != facei || faceRotation[facei] != 0)
                    {
                        changed = true;
                    }
                }
            } // if not procPatch.owner

            if (changed)
            {
                // Merge patch face reordering into mesh face reordering table
                label start = procPatch.start();

                forAll(faceMap, patchFacei)
                {
                    oldToNew[patchFacei + start] =
                        start + faceMap[patchFacei];
                }

                forAll(faceRotation, patchFacei)
                {
                    rotation[patchFacei + start] =
                        faceRotation[patchFacei];
                }

                anyChanged = true;
            }
        } //forAll procBoundaries

        if (anyChanged)
        {
            //- Current point set
            pointField oldPoints(mesh.points());

            //- Current faceList
            faceList faces(mesh.faces());

            //- Owner for all faces
            labelList faceOwner(mesh.faceOwner());

            //- Neighbour for internal faces (-1 for external faces)
            labelList faceNeighbour(mesh.faceNeighbour());

            // Reorder faces according to oldToNew.
            //reorderCompactFaces
            inplaceReorder(oldToNew, faces);
            inplaceReorder(oldToNew, faceOwner);

            // Reorder faceZone flipMap
            labelList newZoneID(faces.size(), -1);
            boolList zoneFlip(faces.size(), false);

            forAll(mesh.faceZones(), zoneI)
            {
                const labelList& faceLabels = mesh.faceZones()[zoneI];
                const boolList& flipMap = mesh.faceZones()[zoneI].flipMap();

                forAll(faceLabels, i)
                {
                    newZoneID[faceLabels[i]] = zoneI;
                    zoneFlip[faceLabels[i]] = flipMap[i];
                }
            }
            inplaceReorder(oldToNew, zoneFlip);
            inplaceReorder(oldToNew, newZoneID);

            forAll(mesh.faceZones(), zoneI)
            {
                labelList faceLabels(mesh.faceZones()[zoneI].size());
                boolList flipMap(mesh.faceZones()[zoneI].size());

                label index = 0;
                forAll(newZoneID, i)
                {
                    if (newZoneID[i] == zoneI)
                    {
                        faceLabels[index] = i;
                        flipMap[index++] = zoneFlip[i];
                    }
                }

                mesh.faceZones()[zoneI].clearAddressing();
                mesh.faceZones()[zoneI].resetAddressing
                (
                    faceLabels,
                    flipMap
                );
            }

            // Rotate faces (rotation is already in new face indices).
            forAll(rotation, facei)
            {
                if (rotation[facei] != 0)
                {
                    inplaceRotateList<List, label>(faces[facei], rotation[facei]);
                }
            }

            mesh.resetPrimitives
            (
                xferMove(oldPoints),
                faces.xfer(),
                faceOwner.xfer(),
                faceNeighbour.xfer(),
                patchSizes,
                patchStarts,
                /*syncParallel*/false
            );
        } //if anyChanged
    } // forAll myProcMeshes receive ordering
} //reorderCoupledFaces


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::fvMeshDistribute::fvMeshDistribute(fvMesh& mesh, const scalar mergeTol)
:
    mesh_(mesh),
    mergeTol_(mergeTol),
    nFinalProcs_(Pstream::nProcs())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvMeshDistribute::countCells
(
    const labelList& distribution
)
{
    labelList nCells(Pstream::nProcs(), 0);
    forAll(distribution, celli)
    {
        label newProc = distribution[celli];

        if (newProc < 0 || newProc >= Pstream::nProcs())
        {
            FatalErrorInFunction
                << "Distribution should be in range 0.." << Pstream::nProcs()-1
                << endl
                << "At index " << celli << " distribution:" << newProc
                << abort(FatalError);
        }
        nCells[newProc]++;
    }
    return nCells;
}

Foam::labelList Foam::fvMeshDistribute::countCells
(
    const labelList& distribution,
    const label nFinalProcs
)
{
    labelList nCells(nFinalProcs, 0);
    forAll(distribution, celli)
    {
        label newProc = distribution[celli];

        if (newProc < 0 || newProc >= nFinalProcs)
        {
            FatalErrorInFunction
                << "Distribution should be in range 0.." << nFinalProcs-1
                << endl
                << "At index " << celli << " distribution:" << newProc
                << abort(FatalError);
        }
        nCells[newProc]++;
    }
    return nCells;
}

template<typename... Args, size_t... Is>
void Foam::fvMeshDistribute::receiveAllFieldsImpl
(
    const label sendProc,
    fvMesh& domainMesh,
    const dictionary& fieldDicts,
    std::index_sequence<Is...>,

    Args&... args
)
{
    std::tuple<Args&...> argsT{args...};

    ((receiveFields
    (
        sendProc,
        std::get<Is * 2 + 1>(argsT),
        domainMesh,
        std::get<Is * 2>(argsT),
        fieldDicts.subDict(std::remove_reference_t<decltype(std::get<Is * 2>(argsT))>::value_type::typeName)
    )), ...);
}

template<typename... Args>
void Foam::fvMeshDistribute::receiveAllFields
(
    const label sendProc,
    fvMesh& domainMesh,
    const dictionary& fieldDicts,

    Args&... args
)
{
    receiveAllFieldsImpl(sendProc, domainMesh, fieldDicts, std::make_integer_sequence<size_t, sizeof...(Args) / 2>{}, args...);
} //receiveAllFields


Foam::labelListList Foam::fvMeshDistribute::calcMeshToProcMap
(
    const label nFinalProcs,
    labelList& nMeshesInProc,
    labelList& meshToProc
)
{
    label meshesPerProc = nFinalProcs/Pstream::nProcs();

    nMeshesInProc = meshesPerProc;

    for (label i = 0; i < nFinalProcs%Pstream::nProcs(); i++)
    {
        nMeshesInProc[i]++;
    }

    // map which final mesh goes to which processor
    label proc = 0, count = 1;
    for (label i = Pstream::nProcs(); i < nFinalProcs; i++)
    {
        if (count == nMeshesInProc[proc])
        {
            proc++;
            count = 1;
        }
        meshToProc[i] = proc;
        count++;
    }
    //e.g. 2 procs, 5 meshes
    // meshToProc 5(0 1 0 0 1)

    labelListList meshNoInProc
    (
        invertOneToMany
        (
            Pstream::nProcs(),
            meshToProc
        )
    );
    //e.g. 2 procs, 5 meshes
    // meshNoInProc
    // 2
    // (
    //   ( 0 2 3 )
    //   ( 1 4 )
    // )

    return meshNoInProc;
} //calcMeshToProcMap


void Foam::fvMeshDistribute::setNFinalProcs
(
    const label nDestProcs
)
{
    nFinalProcs_ = nDestProcs;
}


void Foam::fvMeshDistribute::initMeshList
(
    const fileName& meshSubDir,
    const label nDestProcs,
    fvMesh& baseMesh,

    PtrList<Time>&                 myProcTimes,
    PtrList<fvMesh>&               myProcMeshes,
    List<PtrList<volScalarField>>& volScalarFields,
    List<PtrList<volVectorField>>& volVectorFields
)
{
    bool oldParRun = UPstream::parRun();
    UPstream::parRun() = false;

    labelList nMeshesInProc(Pstream::nProcs());
    labelList meshToProc(identity(nDestProcs));
    labelListList meshNoInProc =
        calcMeshToProcMap
        (
            nDestProcs,
            nMeshesInProc,
            meshToProc
        );

    const label nMyMeshes = nMeshesInProc[Pstream::myProcNo()];

    myProcTimes.setSize(nMyMeshes);
    myProcMeshes.setSize(nMyMeshes);

    List<PtrList<volScalarField>>             volScalarFields0(nMyMeshes);
    List<PtrList<volVectorField>>             volVectorFields0(nMyMeshes);
    //List<PtrList<volSphericalTensorField>>    volSphericalTensorFields0(nMyMeshes);
    //List<PtrList<volSymmTensorField>>         volSymmTensorFields0(nMyMeshes);
    //List<PtrList<volTensorField>>             volTensorFields0(nMyMeshes);
    //List<PtrList<volVector1Field>>            volVector1Fields0(nMyMeshes);
    //List<PtrList<volVector4Field>>            volVector4Fields0(nMyMeshes);
    //List<PtrList<volTensor4Field>>            volTensor4Fields0(nMyMeshes);
    //List<PtrList<surfaceScalarField>>         surfaceScalarFields0(nMyMeshes);
    //List<PtrList<surfaceVectorField>>         surfaceVectorFields0(nMyMeshes);
    //List<PtrList<surfaceSphericalTensorField>>surfaceSphericalTensorFields0(nMyMeshes);
    //List<PtrList<surfaceSymmTensorField>>     surfaceSymmTensorFields0(nMyMeshes);
    //List<PtrList<surfaceTensorField>>         surfaceTensorFields0(nMyMeshes);

    const wordList pointZoneNames(baseMesh.pointZones().names());
    const wordList faceZoneNames (baseMesh.faceZones().names());
    const wordList cellZoneNames (baseMesh.cellZones().names());

    const polyBoundaryMesh& patchEntries = baseMesh.boundaryMesh();

    instantList timeDirs = Time::findTimes(baseMesh.time().path(), "constant");

    forAll(myProcMeshes, i)
    {
        fileName procFolder
        (
            word("processor")
          + Foam::name(meshNoInProc[Pstream::myProcNo()][i])
        );

        if (i > 0)
        {
            forAll(timeDirs, j)
            {
                mkDir(procFolder/timeDirs[j].name());
            }
        }

        fileName processorCasePath
        (
            baseMesh.time().globalCaseName()/procFolder
        );

        autoPtr<Time> procTimePtr
        (
            new Time
            (
                Time::controlDictName,
                baseMesh.time().rootPath(),
                processorCasePath,
                word("system"),
                word("constant"),
                /*enableFunctionObjects*/ false
            )
        );

        procTimePtr->TimeState::operator=(baseMesh.time());
        procTimePtr->setTime(baseMesh.time());

        myProcTimes.set(i, procTimePtr);

        Time& runTime = myProcTimes[i];

        word regionName =
            meshSubDir.component(0) == "polyMesh"
          ? fvMesh::defaultRegion
          : meshSubDir.component(0);

        autoPtr<fvMesh> meshPtr
        (
            new fvMesh
            (
                IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    IOobject::NO_READ
                ),
                xferCopy(pointField(0)),
                xferCopy(faceList(0)),
                xferCopy(cellList(0)),
                false                   // no parallel comms
            )
        );

        // Add patches
        DynamicList<polyPatch*> patches(patchEntries.size());

        forAll(patchEntries, patchi)
        {
            const polyPatch& pp = patchEntries[patchi];

            if (isA<directPolyPatch>(pp))
            {
                const directPolyPatch& dpp =
                    refCast<const directPolyPatch>(pp);

                if (!isA<processorPolyPatch>(dpp))
                {
                    patches.append
                    (
                        dpp.clone
                        (
                            meshPtr().boundaryMesh(),
                            patchi,
                            0, //newSize
                            0  //newStart
                        ).ptr()
                    );
                }
            }
        }

        // Add patches; no parallel comms
        meshPtr().addFvPatches(patches, false);

        // Add zones
        List<pointZone*> pz(pointZoneNames.size());
        forAll(pointZoneNames, i)
        {
            pz[i] =
                new pointZone
                (
                    pointZoneNames[i],
                    labelList(0),
                    i,
                    meshPtr().pointZones()
                );
        }

        List<faceZone*> fz(faceZoneNames.size());
        forAll(faceZoneNames, i)
        {
            fz[i] =
                new faceZone
                (
                    faceZoneNames[i],
                    labelList(0),
                    boolList(0),
                    i,
                    meshPtr().faceZones()
                );
        }

        List<cellZone*> cz(cellZoneNames.size());
        forAll(cellZoneNames, i)
        {
            cz[i] =
                new cellZone
                (
                    cellZoneNames[i],
                    labelList(0),
                    i,
                    meshPtr().cellZones()
                );
        }

        meshPtr().addZones(pz, fz, cz);

        // Subset 0 cells, no parallel comms. This is used to create
        // zero-sized fields.
        autoPtr<fvMeshSubset> subsetterPtr(new fvMeshSubset(baseMesh));
        subsetterPtr().setLargeCellSubset(labelHashSet(0), 0, false);

        // initialize empty fields
        initFields<volScalarField>             (baseMesh, meshPtr(), subsetterPtr, volScalarFields0[i]);
        initFields<volVectorField>             (baseMesh, meshPtr(), subsetterPtr, volVectorFields0[i]);
        //initFields<volSphericalTensorField>    (baseMesh, meshPtr(), subsetterPtr, volSphericalTensorFields0[i]);
        //initFields<volSymmTensorField>         (baseMesh, meshPtr(), subsetterPtr, volSymmTensorFields0[i]);
        //initFields<volTensorField>             (baseMesh, meshPtr(), subsetterPtr, volTensorFields0[i]);
        //initFields<volVector1Field>            (baseMesh, meshPtr(), subsetterPtr, volVector1Fields0[i]);
        //initFields<volVector4Field>            (baseMesh, meshPtr(), subsetterPtr, volVector4Fields0[i]);
        //initFields<volTensor4Field>            (baseMesh, meshPtr(), subsetterPtr, volTensor4Fields0[i]);
        //initFields<surfaceScalarField>         (baseMesh, meshPtr(), subsetterPtr, surfaceScalarFields0[i]);
        //initFields<surfaceVectorField>         (baseMesh, meshPtr(), subsetterPtr, surfaceVectorFields0[i]);
        //initFields<surfaceSphericalTensorField>(baseMesh, meshPtr(), subsetterPtr, surfaceSphericalTensorFields0[i]);
        //initFields<surfaceSymmTensorField>     (baseMesh, meshPtr(), subsetterPtr, surfaceSymmTensorFields0[i]);
        //initFields<surfaceTensorField>         (baseMesh, meshPtr(), subsetterPtr, surfaceTensorFields0[i]);

        myProcMeshes.set(i, meshPtr);
    }

    volScalarFields.transfer             (volScalarFields0);
    volVectorFields.transfer             (volVectorFields0);
    //volSphericalTensorFields.transfer    (volSphericalTensorFields0);
    //volSymmTensorFields.transfer         (volSymmTensorFields0);
    //volTensorFields.transfer             (volTensorFields0);
    //volVector1Fields.transfer            (volVector1Fields0);
    //volVector4Fields.transfer            (volVector4Fields0);
    //volTensor4Fields.transfer            (volTensor4Fields0);
    //surfaceScalarFields.transfer         (surfaceScalarFields0);
    //surfaceVectorFields.transfer         (surfaceVectorFields0);
    //surfaceSphericalTensorFields.transfer(surfaceSphericalTensorFields0);
    //surfaceSymmTensorFields.transfer     (surfaceSymmTensorFields0);
    //surfaceTensorFields.transfer         (surfaceTensorFields0);


    UPstream::parRun() = oldParRun;
} // initMeshList

Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::fvMeshDistribute::distribute
(
    const labelList& distribution,
    PtrList<fvMesh>& myProcMeshes, /*= dummyMeshList*/
    PtrList<mapDistributePolyMesh>& mapsDist /*= dummyMapList*/
)
{
    // Some checks on distribution
    if (distribution.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Size of distribution:"
            << distribution.size() << " mesh nCells:" << mesh_.nCells()
            << abort(FatalError);
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Check all processors have same non-proc patches in same order.
    if (patches.checkParallelSync(true))
    {
        FatalErrorInFunction
            << "This application requires all non-processor patches"
            << " to be present in the same order on all patches" << nl
            << "followed by the processor patches (which of course are unique)."
            << nl
            << "Local patches:" << mesh_.boundaryMesh().names()
            << abort(FatalError);
    }

    // Save some data for mapping later on
    const label nOldPoints(mesh_.nPoints());
    const label nOldFaces(mesh_.nFaces());
    const label nOldCells(mesh_.nCells());
    labelList oldPatchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    forAll(patches, patchi)
    {
        oldPatchStarts[patchi] = patches[patchi].start();
        oldPatchNMeshPoints[patchi] = patches[patchi].nPoints();
    }

    // Save a copy of the poly-face and face-cell-owner addressing to map
    // non-conformal meshes
    const bool origMeshConformal = mesh_.conformal();
    const surfaceLabelField::Boundary origPolyFacesBf
    (
        surfaceLabelField::null(),
        mesh_.polyFacesBf()
    );

    const labelList origFaceOwner(mesh_.faceOwner());

    // Short circuit trivial case.
    if (!Pstream::parRun() && nFinalProcs_ == 1)
    {
        // Collect all maps and return
        return autoPtr<mapDistributePolyMesh>
        (
            new mapDistributePolyMesh
            (
                mesh_,

                nOldPoints,
                nOldFaces,
                nOldCells,
                oldPatchStarts.xfer(),
                oldPatchNMeshPoints.xfer(),

                labelListList(1, identity(mesh_.nPoints())).xfer(),//subPointMap
                labelListList(1, identity(mesh_.nFaces())).xfer(), //subFaceMap
                labelListList(1, identity(mesh_.nCells())).xfer(), //subCellMap
                labelListList(1, identity(patches.size())).xfer(), //subPatchMap

                labelListList(1, identity(mesh_.nPoints())).xfer(),//pointMap
                labelListList(1, identity(mesh_.nFaces())).xfer(), //faceMap
                labelListList(1, identity(mesh_.nCells())).xfer(), //cellMap
                labelListList(1, identity(patches.size())).xfer()  //patchMap
            )
        );
    }


    // Collect any zone names
    const wordList pointZoneNames(mergeWordList(mesh_.pointZones().names()));
    const wordList faceZoneNames(mergeWordList(mesh_.faceZones().names()));
    const wordList cellZoneNames(mergeWordList(mesh_.cellZones().names()));



    // Local environment of all boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // A face is uniquely defined by
    //  - proc
    //  - local face no
    //
    // To glue the parts of meshes which can get sent from anywhere we
    // need to know on boundary faces what the above tuple on both sides is.
    // So we need to maintain:
    //  - original face
    //  - original processor id (= trivial)
    // For coupled boundaries (where the faces are 'duplicate') we take the
    // lowest numbered processor as the data to store.
    //
    // Additionally to create the procboundaries we need to know where the owner
    // cell on the other side now is: newNeighbourProc.
    //

    // physical boundary:
    //     sourceProc = -1
    //     sourceNewNbrProc = -1
    //     sourceFace = -1
    //     sourcePatch = patchID
    // processor boundary:
    //     sourceProc = proc (on owner side)
    //     sourceNewNbrProc = distribution of coupled cell
    //     sourceFace = face (on owner side)
    //     sourcePatch = -1
    // ?cyclic:
    // ?    sourceProc = proc
    // ?    sourceNewNbrProc = distribution of coupled cell
    // ?    sourceFace = face (on owner side)
    // ?    sourcePatch = patchID
    // processor-cyclic boundary:
    //     sourceProc = proc (on owner side)
    //     sourceNewNbrProc = distribution of coupled cell
    //     sourceFace = face (on owner side)
    //     sourcePatch = patchID

    labelList sourceFace;
    labelList sourceProc;
    labelList sourcePatch;
    labelList sourceNbrPatch;
    labelList sourceNewNbrProc;
    labelList sourcePointMaster;
    getCouplingData
    (
        distribution,
        sourceFace,
        sourceProc,
        sourcePatch,
        sourceNbrPatch,
        sourceNewNbrProc,
        sourcePointMaster
    );


    // Remove meshPhi. Since this would otherwise disappear anyway
    // during topo changes and we have to guarantee that all the fields
    // can be sent.
    mesh_.clearOut();
    mesh_.resetMotion();

    // Get data to send. Make sure is synchronised
    const wordList volScalars(mesh_.names(volScalarField::typeName));
    checkEqualWordList("volScalarFields", volScalars);
    const wordList volVectors(mesh_.names(volVectorField::typeName));
    checkEqualWordList("volVectorFields", volVectors);
    const wordList volSphereTensors
    (
        mesh_.names(volSphericalTensorField::typeName)
    );
    checkEqualWordList("volSphericalTensorFields", volSphereTensors);
    const wordList volSymmTensors(mesh_.names(volSymmTensorField::typeName));
    checkEqualWordList("volSymmTensorFields", volSymmTensors);
    const wordList volTensors(mesh_.names(volTensorField::typeName));
    checkEqualWordList("volTensorField", volTensors);

    const wordList volVectors1(mesh_.names(volVector1Field::typeName));
    checkEqualWordList("volVector1Fields", volVectors1);
    const wordList volVectors4(mesh_.names(volVector4Field::typeName));
    checkEqualWordList("volVector4Fields", volVectors4);
    const wordList volTensors4(mesh_.names(volTensor4Field::typeName));
    checkEqualWordList("volTensor4Field", volTensors4);

    const wordList surfScalars(mesh_.names(surfaceScalarField::typeName));
    checkEqualWordList("surfaceScalarFields", surfScalars);
    const wordList surfVectors(mesh_.names(surfaceVectorField::typeName));
    checkEqualWordList("surfaceVectorFields", surfVectors);
    const wordList surfSphereTensors
    (
        mesh_.names(surfaceSphericalTensorField::typeName)
    );
    checkEqualWordList("surfaceSphericalTensorFields", surfSphereTensors);
    const wordList surfSymmTensors
    (
        mesh_.names(surfaceSymmTensorField::typeName)
    );
    checkEqualWordList("surfaceSymmTensorFields", surfSymmTensors);
    const wordList surfTensors(mesh_.names(surfaceTensorField::typeName));
    checkEqualWordList("surfaceTensorFields", surfTensors);

    typedef volScalarField::Internal dimScalType;
    const wordList dimScalars(mesh_.names(dimScalType::typeName));
    checkEqualWordList("volScalarField::Internal", dimScalars);

    typedef volVectorField::Internal dimVecType;
    const wordList dimVectors(mesh_.names(dimVecType::typeName));
    checkEqualWordList("volVectorField::Internal", dimVectors);

    typedef volSphericalTensorField::Internal dimSphereType;
    const wordList dimSphereTensors(mesh_.names(dimSphereType::typeName));
    checkEqualWordList
    (
        "volSphericalTensorField::Internal",
        dimSphereTensors
    );

    typedef volSymmTensorField::Internal dimSymmTensorType;
    const wordList dimSymmTensors(mesh_.names(dimSymmTensorType::typeName));
    checkEqualWordList
    (
        "volSymmTensorField::Internal",
        dimSymmTensors
    );

    typedef volTensorField::Internal dimTensorType;
    const wordList dimTensors(mesh_.names(dimTensorType::typeName));
    checkEqualWordList("volTensorField::Internal", dimTensors);

    typedef volVector1Field::Internal dimVec1Type;
    const wordList dimVectors1(mesh_.names(dimVec1Type::typeName));
    checkEqualWordList("volVector1Field::Internal", dimVectors1);

    typedef volVector4Field::Internal dimVec4Type;
    const wordList dimVectors4(mesh_.names(dimVec4Type::typeName));
    checkEqualWordList("volVector4Field::Internal", dimVectors4);

    typedef volTensor4Field::Internal dimTensor4Type;
    const wordList dimTensors4(mesh_.names(dimTensor4Type::typeName));
    checkEqualWordList("volTensor4Field::Internal", dimTensors4);




    // Find patch to temporarily put exposed and processor faces into.
    label oldInternalPatchi = findNonEmptyPatch();



    // Delete processor patches, starting from the back. Move all faces into
    // oldInternalPatchi.
    labelList repatchFaceMap;
    {
        autoPtr<mapPolyMesh> repatchMap = deleteProcPatches(oldInternalPatchi);

        // Store face map (only face ordering that changed)
        repatchFaceMap = repatchMap().faceMap();

        // Reorder all boundary face data (sourceProc, sourceFace etc.)
        labelList bFaceMap
        (
            SubList<label>
            (
                repatchMap().reverseFaceMap(),
                mesh_.nBoundaryFaces(),
                mesh_.nInternalFaces()
            )
          - mesh_.nInternalFaces()
        );

        inplaceReorder(bFaceMap, sourceFace);
        inplaceReorder(bFaceMap, sourceProc);
        inplaceReorder(bFaceMap, sourcePatch);
        inplaceReorder(bFaceMap, sourceNbrPatch);
        inplaceReorder(bFaceMap, sourceNewNbrProc);
    }



    // Print a bit.
    if (debug)
    {
        Pout<< nl << "MESH WITH PROC PATCHES DELETED:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        printFieldInfo<volVector1Field>(mesh_);
        printFieldInfo<volVector4Field>(mesh_);
        printFieldInfo<volTensor4Field>(mesh_);
        Pout<< nl << endl;
    }

    // case is distributed in more processors than
    // the application is running on (i.e. every processor has to
    // handle one or more meshes)
    if (nFinalProcs_ > Pstream::nProcs())
    {
        // map which final mesh goes to which processor
        labelList nMeshesInProc(Pstream::nProcs());
        labelList meshToProc(identity(nFinalProcs_));
        labelListList meshNoInProc =
            calcMeshToProcMap
            (
                nFinalProcs_,
                nMeshesInProc,
                meshToProc
            );

        Info<< "Distributing from " << Pstream::nProcs()
            << " to " << nFinalProcs_ << " processors." << endl;

        if (debug)
        {
            Info<<"\nFor final redistribution, "
                <<" number of mesh domains for each processor: "
                << nMeshesInProc << endl;
            //e.g. 2 procs, 5 meshes
            // 2(3 2)

            Info<<"\nmeshToProc " << meshToProc <<endl;
            //e.g. 2 procs, 5 meshes
            // 5(0 1 0 0 1)

            Info<<"\nmeshNoInProc: " << meshNoInProc <<endl;
            //e.g. 2 procs, 5 meshes
            // 2
            // (
            //   ( 0 2 3 )
            //   ( 1 4 )
            // )
        }

        const label myProcNo = Pstream::myProcNo();
        const label nMyMeshes = nMeshesInProc[myProcNo];

        mapsDist.setSize(nMyMeshes);

        labelListList sourcesFace       (nMyMeshes, labelList(0));
        labelListList sourcesProc       (nMyMeshes, labelList(0));
        labelListList sourcesPatch      (nMyMeshes, labelList(0));
        labelListList sourcesNbrPatch   (nMyMeshes, labelList(0));
        labelListList sourcesNewNbrProc (nMyMeshes, labelList(0));
        labelListList sourcesPointMaster(nMyMeshes, labelList(0));

        // Maps from subsetted mesh (that is sent) back to original maps
        List<labelListList> subCellMaps (nMyMeshes, labelListList(nFinalProcs_));
        List<labelListList> subFaceMaps (nMyMeshes, labelListList(nFinalProcs_));
        List<labelListList> subPointMaps(nMyMeshes, labelListList(nFinalProcs_));
        List<labelListList> subPatchMaps(nMyMeshes, labelListList(nFinalProcs_));

        // Maps from subsetted mesh to reconstructed mesh
        List<labelListList> constructCellMaps (nMyMeshes, labelListList(nFinalProcs_));
        List<labelListList> constructFaceMaps (nMyMeshes, labelListList(nFinalProcs_));
        List<labelListList> constructPointMaps(nMyMeshes, labelListList(nFinalProcs_));
        List<labelListList> constructPatchMaps(nMyMeshes, labelListList(nFinalProcs_));


        // Find out schedule
        // ~~~~~~~~~~~~~~~~~

        labelListList nSendCells(Pstream::nProcs());
        nSendCells[myProcNo] = countCells(distribution, nFinalProcs_);
        Pstream::allGatherList(nSendCells);
        if (debug)
        {
            Info<< "\nnSendCells " << nSendCells << endl;
        }

        // Allocate buffers (every receiving mesh has its own buffer)
        std::vector<PstreamBuffers> pBufs
        (
            nFinalProcs_,
            PstreamBuffers(Pstream::commsTypes::nonBlocking)
        );

        // Allocate internal streams
        std::vector<OStringStream> meshStr
        (
            nMyMeshes,
            OStringStream(IOstream::BINARY)
        );
        std::vector<OStringStream> fieldStr
        (
            nMyMeshes,
            OStringStream(IOstream::BINARY)
        );
        label nStr = 0;

        // What to send to neighbouring domains
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        forAll(nSendCells[myProcNo], recvMesh)
        {
            if (nSendCells[myProcNo][recvMesh] > 0)
            {
                if (debug)
                {
                    Pout<< nl
                        << "SUBSETTING FOR DOMAIN " << recvMesh
                        << ". Cells to send:"
                        << nSendCells[myProcNo][recvMesh]
                        << nl << endl;
                }

                // Mesh subsetting engine
                fvMeshSubset subsetter(mesh_);

                // Subset the cells of the current domain.
                subsetter.setLargeCellSubset
                (
                    distribution,
                    recvMesh,
                    oldInternalPatchi,  // oldInternalFaces patch
                    false               // no parallel sync
                );

                subCellMaps[0][recvMesh] = subsetter.cellMap();
                subFaceMaps[0][recvMesh] = subsetter.faceFlipMap();
                inplaceRenumberWithFlip
                (
                    repatchFaceMap,
                    false,      // oldToNew has flip
                    true,       // subFaceMap has flip
                    subFaceMaps[0][recvMesh]
                );
                subPointMaps[0][recvMesh] = subsetter.pointMap();
                subPatchMaps[0][recvMesh] = subsetter.patchMap();


                // Subset the boundary fields (owner/neighbour/processor)
                labelList procSourceFace;
                labelList procSourceProc;
                labelList procSourcePatch;
                labelList procSourceNbrPatch;
                labelList procSourceNewNbrProc;
                labelList procSourcePointMaster;

                subsetCouplingData
                (
                    subsetter.subMesh(),
                    subsetter.pointMap(),       // from subMesh to mesh
                    subsetter.faceMap(),        //      ,,      ,,
                    subsetter.cellMap(),        //      ,,      ,,

                    distribution,               // old mesh distribution
                    mesh_.faceOwner(),          // old owner
                    mesh_.faceNeighbour(),
                    mesh_.nInternalFaces(),

                    sourceFace,
                    sourceProc,
                    sourcePatch,
                    sourceNbrPatch,
                    sourceNewNbrProc,
                    sourcePointMaster,

                    procSourceFace,
                    procSourceProc,
                    procSourcePatch,
                    procSourceNbrPatch,
                    procSourceNewNbrProc,
                    procSourcePointMaster
                );

                const label recvProc = meshToProc[recvMesh];

                if (debug)
                    Pout<< "Sending domain " << recvMesh;

                if (recvProc != myProcNo)
                {
                    if (debug)
                    {
                        Pout<< " to processor " << recvProc << nl
                            << "    nCells:" << subsetter.subMesh().nCells() << nl
                            << endl;
                    }

                    // Pstream for sending mesh and fields
                    UOPstream str(recvProc, pBufs[recvMesh]);

                    // Send to neighbour processor
                    sendMesh
                    (
                        recvProc,
                        subsetter.subMesh(),

                        pointZoneNames,
                        faceZoneNames,
                        cellZoneNames,

                        procSourceFace,
                        procSourceProc,
                        procSourcePatch,
                        procSourceNbrPatch,
                        procSourceNewNbrProc,
                        procSourcePointMaster,

                        str
                    );

                    // volFields
                    sendFields<volScalarField>(recvProc, volScalars, subsetter, str);
                    sendFields<volVectorField>(recvProc, volVectors, subsetter, str);
                    sendFields<volSphericalTensorField>(recvProc, volSphereTensors, subsetter, str);
                    sendFields<volSymmTensorField>(recvProc, volSymmTensors, subsetter, str);
                    sendFields<volTensorField>(recvProc, volTensors, subsetter, str);

                    sendFields<volVector1Field>(recvProc, volVectors1, subsetter, str);
                    sendFields<volVector4Field>(recvProc, volVectors4, subsetter, str);
                    sendFields<volTensor4Field>(recvProc, volTensors4, subsetter, str);

                    // surfaceFields
                    sendFields<surfaceScalarField>(recvProc, surfScalars, subsetter, str);
                    sendFields<surfaceVectorField>(recvProc, surfVectors, subsetter, str);
                    sendFields<surfaceSphericalTensorField>(recvProc, surfSphereTensors, subsetter, str);
                    sendFields<surfaceSymmTensorField>(recvProc, surfSymmTensors, subsetter, str);
                    sendFields<surfaceTensorField>(recvProc, surfTensors, subsetter, str);

                    // Dimensioned fields
                    sendFields<volScalarField::Internal>(recvProc, dimScalars, subsetter, str);
                    sendFields<volVectorField::Internal>(recvProc, dimVectors, subsetter, str);
                    sendFields<volSphericalTensorField::Internal>(recvProc, dimSphereTensors, subsetter, str);
                    sendFields<volSymmTensorField::Internal>(recvProc, dimSymmTensors, subsetter, str);
                    sendFields<volTensorField::Internal>(recvProc, dimTensors, subsetter, str);

                    sendFields<volVector1Field::Internal>(recvProc, dimVectors1, subsetter, str);
                    sendFields<volVector4Field::Internal>(recvProc, dimVectors4, subsetter, str);
                    sendFields<volTensor4Field::Internal>(recvProc, dimTensors4, subsetter, str);
                }
                else
                {
                    if (debug)
                    {
                        Pout<< "    nCells:" << subsetter.subMesh().nCells() << nl
                            << " (internal stream nStr " << nStr << ")" << nl
                            << endl;
                    }

                    // send to local stream
                    sendMesh
                    (
                        recvProc,
                        subsetter.subMesh(),

                        pointZoneNames,
                        faceZoneNames,
                        cellZoneNames,

                        procSourceFace,
                        procSourceProc,
                        procSourcePatch,
                        procSourceNbrPatch,
                        procSourceNewNbrProc,
                        procSourcePointMaster,

                        meshStr[nStr]
                    );

                    // volFields
                    sendFields<volScalarField>(recvProc, volScalars, subsetter, fieldStr[nStr]);
                    sendFields<volVectorField>(recvProc, volVectors, subsetter, fieldStr[nStr]);
                    sendFields<volSphericalTensorField>(recvProc, volSphereTensors, subsetter, fieldStr[nStr]);
                    sendFields<volSymmTensorField>(recvProc, volSymmTensors, subsetter, fieldStr[nStr]);
                    sendFields<volTensorField>(recvProc, volTensors, subsetter, fieldStr[nStr]);

                    sendFields<volVector1Field>(recvProc, volVectors1, subsetter, fieldStr[nStr]);
                    sendFields<volVector4Field>(recvProc, volVectors4, subsetter, fieldStr[nStr]);
                    sendFields<volTensor4Field>(recvProc, volTensors4, subsetter, fieldStr[nStr]);

                    // surfaceFields
                    sendFields<surfaceScalarField>(recvProc, surfScalars, subsetter, fieldStr[nStr]);
                    sendFields<surfaceVectorField>(recvProc, surfVectors, subsetter, fieldStr[nStr]);
                    sendFields<surfaceSphericalTensorField>(recvProc, surfSphereTensors, subsetter, fieldStr[nStr]);
                    sendFields<surfaceSymmTensorField>(recvProc, surfSymmTensors, subsetter, fieldStr[nStr]);
                    sendFields<surfaceTensorField>(recvProc, surfTensors, subsetter, fieldStr[nStr]);

                    // Dimensioned fields
                    sendFields<volScalarField::Internal>(recvProc, dimScalars, subsetter, fieldStr[nStr]);
                    sendFields<volVectorField::Internal>(recvProc, dimVectors, subsetter, fieldStr[nStr]);
                    sendFields<volSphericalTensorField::Internal>(recvProc, dimSphereTensors, subsetter, fieldStr[nStr]);
                    sendFields<volSymmTensorField::Internal>(recvProc, dimSymmTensors, subsetter, fieldStr[nStr]);
                    sendFields<volTensorField::Internal>(recvProc, dimTensors, subsetter, fieldStr[nStr]);

                    sendFields<volVector1Field::Internal>(recvProc, dimVectors1, subsetter, fieldStr[nStr]);
                    sendFields<volVector4Field::Internal>(recvProc, dimVectors4, subsetter, fieldStr[nStr]);
                    sendFields<volTensor4Field::Internal>(recvProc, dimTensors4, subsetter, fieldStr[nStr]);

                    nStr++;
                } // send internally
            } //if there are cells for other domains
        } //forAll nSendCells


        UPstream::parRun() = oldParRun;


        // Start sending&receiving from buffers
        for (auto i = 0; i < nFinalProcs_; ++i)
        {
            pBufs[i].finishedSends();
        }

        oldParRun = UPstream::parRun();
        UPstream::parRun() = false;


        // Receive and add what was sent
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(nSendCells, sendProc)
        {
            nStr = 0;

            for (label i = 0; i < nMyMeshes; i++)
            {
                label sendMesh = meshNoInProc[myProcNo][i];

                // Did processor sendProc send anything to me?
                if (nSendCells[sendProc][sendMesh] == 0) continue;

                // Receive from sendProc
                labelList domainSourceFace;
                labelList domainSourceProc;
                labelList domainSourcePatch;
                labelList domainSourceNbrPatch;
                labelList domainSourceNewNbrProc;
                labelList domainSourcePointMaster;

                PtrList<volScalarField> vsf;
                PtrList<volVectorField> vvf;
                PtrList<volSphericalTensorField> vsptf;
                PtrList<volSymmTensorField> vsytf;
                PtrList<volTensorField> vtf;

                PtrList<volVector1Field> vvf1;
                PtrList<volVector4Field> vvf4;
                PtrList<volTensor4Field> vtf4;

                PtrList<surfaceScalarField> ssf;
                PtrList<surfaceVectorField> svf;
                PtrList<surfaceSphericalTensorField> ssptf;
                PtrList<surfaceSymmTensorField> ssytf;
                PtrList<surfaceTensorField> stf;

                PtrList<volScalarField::Internal> dsf;
                PtrList<volVectorField::Internal> dvf;
                PtrList<volSphericalTensorField::Internal> dstf;
                PtrList<volSymmTensorField::Internal> dsytf;
                PtrList<volTensorField::Internal> dtf;

                PtrList<volVector1Field::Internal> dvf1;
                PtrList<volVector4Field::Internal> dvf4;
                PtrList<volTensor4Field::Internal> dtf4;

                autoPtr<fvMesh> domainMeshPtr;

                {
                    dictionary fieldDicts;

                    if (debug)
                    {
                        Pout<< "Receiving domain " << sendMesh;
                    }

                    // Opposite of sendMesh
                    if (sendProc != myProcNo)
                    {
                        if (debug)
                        {
                            Pout<< " from processor " << sendProc << nl
                                << "    nCells:" << nSendCells[sendProc][sendMesh];
                        }

                        // Pstream for receiving mesh and fields
                        UIPstream str(sendProc, pBufs[sendMesh]);

                        domainMeshPtr =
                            receiveMesh
                            (
                                sendProc,
                                pointZoneNames,
                                faceZoneNames,
                                cellZoneNames,

                                const_cast<Time&>(mesh_.time()),
                                domainSourceFace,
                                domainSourceProc,
                                domainSourcePatch,
                                domainSourceNbrPatch,
                                domainSourceNewNbrProc,
                                domainSourcePointMaster,
                                str
                            );
                        if (debug)
                        {
                            Pout<< " - received: " << domainMeshPtr().nCells() << nl
                                << endl;
                        }

                        // Receive fields. Read as single dictionary because
                        // of problems reading consecutive fields from single stream.
                        fieldDicts = dictionary(str);
                    }
                    else
                    {
                        if (debug)
                        {
                            Pout<< "    nCells:" << nSendCells[sendProc][sendMesh] << nl
                                << " (internal stream nStr " << nStr << ")" << nl
                                << endl;
                        }
                        // receive mesh from local OStringStream meshStr
                        // and fields from fieldStr
                        domainMeshPtr =
                            receiveMesh
                            (
                                sendProc,
                                pointZoneNames,
                                faceZoneNames,
                                cellZoneNames,

                                const_cast<Time&>(mesh_.time()),
                                domainSourceFace,
                                domainSourceProc,
                                domainSourcePatch,
                                domainSourceNbrPatch,
                                domainSourceNewNbrProc,
                                domainSourcePointMaster,
                                IStringStream
                                (
                                    meshStr[nStr].str(),
                                    IOstream::BINARY
                                )()
                            );

                        fieldDicts =
                            dictionary
                            (
                                IStringStream
                                (
                                    fieldStr[nStr].str(),
                                    IOstream::BINARY
                                )()
                            );

                        nStr++;
                    }

                    fvMesh& domainMesh = domainMeshPtr();
                    // Force construction of various on mesh.
                    //(void)domainMesh.globalData();

                    receiveAllFields
                    (
                        sendProc,
                        domainMesh,
                        fieldDicts,

                        vsf,   volScalars,
                        vvf,   volVectors,
                        vsptf, volSphereTensors,
                        vsytf, volSymmTensors,
                        vtf,   volTensors,
                        vvf1,  volVectors1,
                        vvf4,  volVectors4,
                        vtf4,  volTensors4,
                        ssf,   surfScalars,
                        svf,   surfVectors,
                        ssptf, surfSphereTensors,
                        ssytf, surfSymmTensors,
                        stf,   surfTensors,
                        dsf,   dimScalars,
                        dvf,   dimVectors,
                        dstf,  dimSphereTensors,
                        dsytf, dimSymmTensors,
                        dtf,   dimTensors,
                        dvf1,  dimVectors1,
                        dvf4,  dimVectors4,
                        dtf4,  dimTensors4
                    );
                }
                const fvMesh& domainMesh = domainMeshPtr();

                // Print a bit.
                if (debug)
                {
                    Pout<< nl << "RECEIVED DOMAIN " << sendMesh
                        << " FROM PROC " << sendProc << endl;
                    printMeshInfo(domainMesh);
                    printFieldInfo<volScalarField>(domainMesh);
                    printFieldInfo<volVectorField>(domainMesh);
                    printFieldInfo<volSphericalTensorField>(domainMesh);
                    printFieldInfo<volSymmTensorField>(domainMesh);
                    printFieldInfo<volTensorField>(domainMesh);
                    printFieldInfo<surfaceScalarField>(domainMesh);
                    printFieldInfo<surfaceVectorField>(domainMesh);
                    printFieldInfo<surfaceSphericalTensorField>(domainMesh);
                    printFieldInfo<surfaceSymmTensorField>(domainMesh);
                    printFieldInfo<surfaceTensorField>(domainMesh);

                    printFieldInfo<volVector1Field>(domainMesh);
                    printFieldInfo<volVector4Field>(domainMesh);
                    printFieldInfo<volTensor4Field>(domainMesh);
                }

                constructCellMaps [i][sendProc] = identity(domainMesh.nCells());
                constructFaceMaps [i][sendProc] = identity(domainMesh.nFaces()) + 1;
                constructPointMaps[i][sendProc] = identity(domainMesh.nPoints());
                constructPatchMaps[i][sendProc] =
                    identity(domainMesh.boundaryMesh().size());

                fvMesh& currentMesh = myProcMeshes[i];

                // Now this mesh we received (from sendProc) needs to be merged
                // with the current mesh. On the current mesh we have for all
                // boundaryfaces the original face and processor. See if we can
                // match these up to the received domainSourceFace and
                // domainSourceProc.
                labelList masterCoupledFaces;
                labelList slaveCoupledFaces;

                findCouples
                (
                    currentMesh,

                    sourcesFace[i],
                    sourcesProc[i],
                    sourcesPatch[i],

                    sendProc,
                    domainMesh,
                    domainSourceFace,
                    domainSourceProc,
                    domainSourcePatch,

                    masterCoupledFaces,
                    slaveCoupledFaces
                );

                // Generate additional coupling info (points, edges) from
                // faces-that-match
                faceCoupleInfo couples
                (
                    currentMesh,
                    masterCoupledFaces,
                    domainMesh,
                    slaveCoupledFaces,
                    mergeTol_,              // merge tolerance
                    true,                   // faces align
                    true,                   // couples are ordered already
                    false
                );

                // Add domainMesh to mesh
                // ~~~~~~~~~~~~~~~~~~~~~~

                autoPtr<mapAddedPolyMesh> map =
                    fvMeshAdder::add
                    (
                        currentMesh,
                        domainMesh,
                        couples,
                        false,      // no parallel comms
                        true        // fullyMapped (suppress mapper warning)
                    );

                // Update mesh data: sourceFace,sourceProc for added
                // mesh.

                sourcesFace[i] = mapBoundaryData
                (
                    currentMesh,
                    map(),
                    sourcesFace[i],
                    domainMesh.nInternalFaces(),
                    domainSourceFace
                );
                sourcesProc[i] = mapBoundaryData
                (
                    currentMesh,
                    map(),
                    sourcesProc[i],
                    domainMesh.nInternalFaces(),
                    domainSourceProc
                );
                sourcesPatch[i] = mapBoundaryData
                (
                    currentMesh,
                    map(),
                    sourcesPatch[i],
                    domainMesh.nInternalFaces(),
                    domainSourcePatch
                );
                sourcesNbrPatch[i] = mapBoundaryData
                (
                    currentMesh,
                    map(),
                    sourcesNbrPatch[i],
                    domainMesh.nInternalFaces(),
                    domainSourceNbrPatch
                );
                sourcesNewNbrProc[i] = mapBoundaryData
                (
                    currentMesh,
                    map(),
                    sourcesNewNbrProc[i],
                    domainMesh.nInternalFaces(),
                    domainSourceNewNbrProc
                );
                // Update pointMaster data
                sourcesPointMaster[i] = mapPointData
                (
                    currentMesh,
                    map(),
                    sourcesPointMaster[i],
                    domainSourcePointMaster
                );

                // Update all addressing so xxProcAddressing points to correct
                // item in masterMesh.
                const labelList& oldCellMap = map().oldCellMap();
                const labelList& oldFaceMap = map().oldFaceMap();
                const labelList& oldPointMap = map().oldPointMap();
                const labelList& oldPatchMap = map().oldPatchMap();

                //Note: old mesh faces never flipped!
                forAll(constructPatchMaps[i], meshi)
                {
                    if (meshi != sendProc && constructPatchMaps[i][meshi].size())
                    {
                        // Processor already in mesh (either myProcNo or received)
                        inplaceRenumber(oldCellMap, constructCellMaps[i][meshi]);
                        inplaceRenumberWithFlip
                        (
                            oldFaceMap,
                            false,
                            true,
                            constructFaceMaps[i][meshi]
                        );
                        inplaceRenumber(oldPointMap, constructPointMaps[i][meshi]);
                        inplaceRenumber(oldPatchMap, constructPatchMaps[i][meshi]);
                    }
                }

                labelHashSet flippedAddedFaces;
                {
                    // Find out if any faces of domain mesh were flipped (boundary
                    // faces becoming internal)
                    const label nBnd = domainMesh.nBoundaryFaces();
                    flippedAddedFaces.resize(nBnd/4);

                    for
                    (
                        label domainFaceI = domainMesh.nInternalFaces();
                        domainFaceI < domainMesh.nFaces();
                        domainFaceI++
                    )
                    {
                        label newFaceI = map().addedFaceMap()[domainFaceI];
                        label newCellI = currentMesh.faceOwner()[newFaceI];

                        label domainCellI = domainMesh.faceOwner()[domainFaceI];

                        if (newCellI != map().addedCellMap()[domainCellI])
                        {
                            flippedAddedFaces.insert(domainFaceI);
                        }
                    }
                }

                // Added processor
                inplaceRenumber(map().addedCellMap(), constructCellMaps[i][sendProc]);
                // Add flip
                for (const label domainFaceI : flippedAddedFaces)
                {
                    label& val = constructFaceMaps[i][sendProc][domainFaceI];
                    val = -val;
                }
                inplaceRenumberWithFlip
                (
                    map().addedFaceMap(),
                    false,
                    true,           // constructFaceMap has flip sign
                    constructFaceMaps[i][sendProc]
                );
                inplaceRenumber(map().addedPointMap(), constructPointMaps[i][sendProc]);
                inplaceRenumber(map().addedPatchMap(), constructPatchMaps[i][sendProc]);

                if (debug)
                {
                    Pout<< nl << "MERGED DOMAIN " << sendMesh
                        << " FROM PROC " << sendProc << endl;
                    printMeshInfo(currentMesh);
                    printFieldInfo<volScalarField>(currentMesh);
                    printFieldInfo<volVectorField>(currentMesh);
                    printFieldInfo<volSphericalTensorField>(currentMesh);
                    printFieldInfo<volSymmTensorField>(currentMesh);
                    printFieldInfo<volTensorField>(currentMesh);
                    printFieldInfo<surfaceScalarField>(currentMesh);
                    printFieldInfo<surfaceVectorField>(currentMesh);
                    printFieldInfo<surfaceSphericalTensorField>(currentMesh);
                    printFieldInfo<surfaceSymmTensorField>(currentMesh);
                    printFieldInfo<surfaceTensorField>(currentMesh);

                    printFieldInfo<volVector1Field>(currentMesh);
                    printFieldInfo<volVector4Field>(currentMesh);
                    printFieldInfo<volTensor4Field>(currentMesh);
                    Pout<< nl << endl;
                }
            } // forAll domains in proc (nSendCells[sendProc])
        } // forAll procs (nSendCells)

        UPstream::parRun() = oldParRun;

        // Print a bit.
        if (debug)
        {
            forAll(myProcMeshes, i)
            {
                fvMesh& mesh = myProcMeshes[i];

                Pout<< nl << "REDISTRIBUTED MESH "
                    << meshNoInProc[myProcNo][i] << endl;
                printMeshInfo(mesh);
                printFieldInfo<volScalarField>(mesh);
                printFieldInfo<volVectorField>(mesh);
                printFieldInfo<volSphericalTensorField>(mesh);
                printFieldInfo<volSymmTensorField>(mesh);
                printFieldInfo<volTensorField>(mesh);
                printFieldInfo<surfaceScalarField>(mesh);
                printFieldInfo<surfaceVectorField>(mesh);
                printFieldInfo<surfaceSphericalTensorField>(mesh);
                printFieldInfo<surfaceSymmTensorField>(mesh);
                printFieldInfo<surfaceTensorField>(mesh);

                printFieldInfo<volVector1Field>(mesh);
                printFieldInfo<volVector4Field>(mesh);
                printFieldInfo<volTensor4Field>(mesh);
                Pout<< nl << endl;
            }
        }


        // See if any originally shared points need to be merged. Note: does
        // parallel comms. After this points and edges should again be consistent.
        //mergeSharedPoints(sourcePointMaster, constructPointMap);
        forAll(myProcMeshes, i)
        {
            if (debug)
            {
                Pout<<"Merging points of mesh "
                    << meshNoInProc[myProcNo][i] << endl;
            }


            fvMesh& mesh = myProcMeshes[i];
            mergeSharedPoints(sourcesPointMaster[i], mesh, constructPointMaps[i]);


            // Add processorPatches
            // ~~~~~~~~~~~~~~~~~~~~

            // Per neighbour processor, per originating patch, the patchID
            // For faces resulting from internal faces or normal processor patches
            // the originating patch is -1. For cyclics this is the cyclic patchID.
            List<Map<label>> procPatchID;

            // Add processor and processorCyclic patches.
            addProcPatches
            (
                sourcesNewNbrProc[i],
                sourcesPatch[i],
                sourcesNbrPatch[i],
                meshNoInProc[myProcNo][i],
                mesh,
                procPatchID
            );


            // Put faces into correct patch. Note that we now have proper
            // processorPolyPatches again so repatching will take care of coupled face
            // ordering.

            // Get boundary faces to be repatched. Is -1 or new patchID
            labelList newPatchID
            (
                getBoundaryPatch
                (
                    sourcesNewNbrProc[i],
                    sourcesPatch[i],
                    procPatchID,
                    meshNoInProc[myProcNo][i]
                )
            );

            // Change patches. Since this might change ordering of coupled faces
            // we also need to adapt our constructMaps.
            repatch(newPatchID, mesh, constructFaceMaps[i]);

            // Bit of hack: processorFvPatchField does not get reset since created
            // from nothing so explicitly reset.
            initPatchFields<volScalarField, processorFvPatchField<scalar>>
            (
                Zero, mesh
            );
            initPatchFields<volVectorField, processorFvPatchField<vector>>
            (
                Zero, mesh
            );
            initPatchFields
            <
                volSphericalTensorField,
                processorFvPatchField<sphericalTensor>
            >
            (
                Zero, mesh
            );
            initPatchFields<volSymmTensorField, processorFvPatchField<symmTensor>>
            (
                Zero, mesh
            );
            initPatchFields<volTensorField, processorFvPatchField<tensor>>
            (
                Zero, mesh
            );

            initPatchFields<volVector1Field, processorFvPatchField<vector1>>
            (
                Zero, mesh
            );
            initPatchFields<volVector4Field, processorFvPatchField<vector4>>
            (
                Zero, mesh
            );
            initPatchFields<volTensor4Field, processorFvPatchField<tensor4>>
            (
                Zero, mesh
            );


            mesh.setInstance(mesh_.time().timeName());

            // Print a bit
            if (debug)
            {
                Pout<< nl << "FINAL MESH "
                    << meshNoInProc[myProcNo][i] << endl;
                printMeshInfo(mesh);
                printFieldInfo<volScalarField>(mesh);
                printFieldInfo<volVectorField>(mesh);
                printFieldInfo<volSphericalTensorField>(mesh);
                printFieldInfo<volSymmTensorField>(mesh);
                printFieldInfo<volTensorField>(mesh);
                printFieldInfo<surfaceScalarField>(mesh);
                printFieldInfo<surfaceVectorField>(mesh);
                printFieldInfo<surfaceSphericalTensorField>(mesh);
                printFieldInfo<surfaceSymmTensorField>(mesh);
                printFieldInfo<surfaceTensorField>(mesh);

                printFieldInfo<volVector1Field>(mesh);
                printFieldInfo<volVector4Field>(mesh);
                printFieldInfo<volTensor4Field>(mesh);
                Pout<< nl << endl;
            }
        } // forAll myProcMeshes add processor patches

        // Reorder processor faces
        reorderCoupledFaces
        (
            myProcMeshes,
            meshToProc,
            meshNoInProc
        );

        forAll(myProcMeshes, i)
        {
            fvMesh& mesh = myProcMeshes[i];

            if (debug)
            {
                OFstream str
                (
                    "mapDistributePolyMesh_proc"
                   +Foam::name(meshNoInProc[myProcNo][i])
                );

                str << "mesh.nCells() " << mesh.nCells() << nl
                    //<< "nOldPoints " << nOldPoints << nl
                    //<< "nOldFaces " << nOldFaces << nl
                    << "nOldCells " << (i > 0 ? 0 : nOldCells) << nl
                    //<< "oldPatchStarts " << oldPatchStarts << nl
                    //<< "oldPatchNMeshPoints " << oldPatchNMeshPoints << nl

                    //<< "subPointMap " << subPointMap << nl
                    //<< "subFaceMap " << subFaceMap << nl
                    << "subCellMap " << subCellMaps[i] << nl
                    //<< "subPatchMap " << subPatchMap << nl

                    //<< "constructPointMap " << constructPointMap << nl
                    //<< "constructFaceMap " << constructFaceMap << nl
                    << "constructCellMap " << constructCellMaps[i] << nl
                    //<< "constructPatchMap " << constructPatchMap << nl
                    << endl;
            }

            labelList emptyList(patches.size(), 0);

            // Collect this mesh map
            autoPtr<mapDistributePolyMesh> mapPtr
            (
                new mapDistributePolyMesh
                (
                    mesh,

                    (i > 0 ? 0 : nOldPoints),
                    (i > 0 ? 0 : nOldFaces),
                    (i > 0 ? 0 : nOldCells),
                    (i > 0 ? emptyList.xfer() : oldPatchStarts.xfer()),
                    (i > 0 ? emptyList.xfer() : oldPatchNMeshPoints.xfer()),

                    subPointMaps[i].xfer(),
                    subFaceMaps[i].xfer(),
                    subCellMaps[i].xfer(),
                    subPatchMaps[i].xfer(),

                    constructPointMaps[i].xfer(),
                    constructFaceMaps[i].xfer(),
                    constructCellMaps[i].xfer(),
                    constructPatchMaps[i].xfer(),

                    true,           // subFaceMap has flip
                    true            // constructFaceMap has flip
                )
            );

            mapsDist.set(i, mapPtr);
        }

        if (!origMeshConformal)
        {
            // Add any nonConformalProcessorCyclic patches. These are empty
            // at the moment, as the mesh should not be stitched when
            // redistributing. So, we don't need to do any re-patching here
            // like with the standard processor patches.
            addNccProcPatches(myProcMeshes, meshNoInProc);

            // Unconform the processor meshes, so the finite-volume addressing
            // gets written down.
            unconform
            (
                myProcMeshes,
                meshNoInProc,
                mapsDist,
                origPolyFacesBf,
                origFaceOwner,
                distribution
            );
        }

        // reset mesh_
        {
            pointField newPoints(0);
            faceList newFaces(0);
            labelList newOwner(0);
            labelList newNeighbour(0);

            const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

            labelList patchSizes(boundary.size());
            labelList patchStarts(boundary.size());

            mesh_.resetPrimitives
            (
                xferMove(newPoints),
                newFaces.xfer(),
                newOwner.xfer(),
                newNeighbour.xfer(),
                patchSizes,
                patchStarts,
                /*syncParallel*/false
            );
        }

        // reset fvPatchFields updated_ flag
        resetUpdate<volScalarField>();
        resetUpdate<volVectorField>();
        resetUpdate<volSphericalTensorField>();
        resetUpdate<volSymmTensorField>();
        resetUpdate<volTensorField>();
        resetUpdate<volVector1Field>();
        resetUpdate<volVector4Field>();
        resetUpdate<volTensor4Field>();

        // Return empty map for backward compatibility
        return autoPtr<mapDistributePolyMesh>
        (
            //mapsDist.set(0, nullptr)
            new mapDistributePolyMesh()
        );
    }
    else // standard behaviour where every proc has one mesh
    {
        // Maps from subsetted mesh (that is sent) back to original maps
        labelListList subCellMap(Pstream::nProcs());
        labelListList subFaceMap(Pstream::nProcs());
        labelListList subPointMap(Pstream::nProcs());
        labelListList subPatchMap(Pstream::nProcs());
        // Maps from subsetted mesh to reconstructed mesh
        labelListList constructCellMap(Pstream::nProcs());
        labelListList constructFaceMap(Pstream::nProcs());
        labelListList constructPointMap(Pstream::nProcs());
        labelListList constructPatchMap(Pstream::nProcs());


        // Find out schedule
        // ~~~~~~~~~~~~~~~~~

        labelListList nSendCells(Pstream::nProcs());
        nSendCells[Pstream::myProcNo()] = countCells(distribution);
        Pstream::allGatherList(nSendCells);


        // Allocate buffers
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);


        // What to send to neighbouring domains
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        forAll(nSendCells[Pstream::myProcNo()], recvProc)
        {
            if
            (
                recvProc != Pstream::myProcNo()
             && nSendCells[Pstream::myProcNo()][recvProc] > 0
            )
            {
                // Send to recvProc

                if (debug)
                {
                    Pout<< nl
                        << "SUBSETTING FOR DOMAIN " << recvProc
                        << " cells to send:"
                        << nSendCells[Pstream::myProcNo()][recvProc]
                        << nl << endl;
                }

                // Pstream for sending mesh and fields
                //OPstream str(Pstream::commsTypes::blocking, recvProc);
                UOPstream str(recvProc, pBufs);

                // Mesh subsetting engine
                fvMeshSubset subsetter(mesh_);

                // Subset the cells of the current domain.
                subsetter.setLargeCellSubset
                (
                    distribution,
                    recvProc,
                    oldInternalPatchi,  // oldInternalFaces patch
                    false               // no parallel sync
                );

                subCellMap[recvProc] = subsetter.cellMap();
                subFaceMap[recvProc] = subsetter.faceFlipMap();
                inplaceRenumberWithFlip
                (
                    repatchFaceMap,
                    false,      // oldToNew has flip
                    true,       // subFaceMap has flip
                    subFaceMap[recvProc]
                );
                subPointMap[recvProc] = subsetter.pointMap();
                subPatchMap[recvProc] = subsetter.patchMap();


                // Subset the boundary fields (owner/neighbour/processor)
                labelList procSourceFace;
                labelList procSourceProc;
                labelList procSourcePatch;
                labelList procSourceNbrPatch;
                labelList procSourceNewNbrProc;
                labelList procSourcePointMaster;

                subsetCouplingData
                (
                    subsetter.subMesh(),
                    subsetter.pointMap(),       // from subMesh to mesh
                    subsetter.faceMap(),        //      ,,      ,,
                    subsetter.cellMap(),        //      ,,      ,,

                    distribution,               // old mesh distribution
                    mesh_.faceOwner(),          // old owner
                    mesh_.faceNeighbour(),
                    mesh_.nInternalFaces(),

                    sourceFace,
                    sourceProc,
                    sourcePatch,
                    sourceNbrPatch,
                    sourceNewNbrProc,
                    sourcePointMaster,

                    procSourceFace,
                    procSourceProc,
                    procSourcePatch,
                    procSourceNbrPatch,
                    procSourceNewNbrProc,
                    procSourcePointMaster
                );


                // Send to neighbour
                sendMesh
                (
                    recvProc,
                    subsetter.subMesh(),

                    pointZoneNames,
                    faceZoneNames,
                    cellZoneNames,

                    procSourceFace,
                    procSourceProc,
                    procSourcePatch,
                    procSourceNbrPatch,
                    procSourceNewNbrProc,
                    procSourcePointMaster,

                    str
                );

                // volFields
                sendFields<volScalarField>(recvProc, volScalars, subsetter, str);
                sendFields<volVectorField>(recvProc, volVectors, subsetter, str);
                sendFields<volSphericalTensorField>
                (
                    recvProc,
                    volSphereTensors,
                    subsetter,
                    str
                );
                sendFields<volSymmTensorField>
                (
                    recvProc,
                    volSymmTensors,
                    subsetter,
                    str
                );
                sendFields<volTensorField>(recvProc, volTensors, subsetter, str);


                sendFields<volVector1Field>(recvProc, volVectors1, subsetter, str);
                sendFields<volVector4Field>(recvProc, volVectors4, subsetter, str);
                sendFields<volTensor4Field>(recvProc, volTensors4, subsetter, str);


                // surfaceFields
                sendFields<surfaceScalarField>
                (
                    recvProc,
                    surfScalars,
                    subsetter,
                    str
                );
                sendFields<surfaceVectorField>
                (
                    recvProc,
                    surfVectors,
                    subsetter,
                    str
                );
                sendFields<surfaceSphericalTensorField>
                (
                    recvProc,
                    surfSphereTensors,
                    subsetter,
                    str
                );
                sendFields<surfaceSymmTensorField>
                (
                    recvProc,
                    surfSymmTensors,
                    subsetter,
                    str
                );
                sendFields<surfaceTensorField>
                (
                    recvProc,
                    surfTensors,
                    subsetter,
                    str
                );

                // Dimensioned fields
                sendFields<volScalarField::Internal>
                (
                    recvProc,
                    dimScalars,
                    subsetter,
                    str
                );
                sendFields<volVectorField::Internal>
                (
                    recvProc,
                    dimVectors,
                    subsetter,
                    str
                );
                sendFields<volSphericalTensorField::Internal>
                (
                    recvProc,
                    dimSphereTensors,
                    subsetter,
                    str
                );
                sendFields<volSymmTensorField::Internal>
                (
                    recvProc,
                    dimSymmTensors,
                    subsetter,
                    str
                );
                sendFields<volTensorField::Internal>
                (
                    recvProc,
                    dimTensors,
                    subsetter,
                    str
                );


                sendFields<volVector1Field::Internal>
                (
                    recvProc,
                    dimVectors1,
                    subsetter,
                    str
                );

                sendFields<volVector4Field::Internal>
                (
                    recvProc,
                    dimVectors4,
                    subsetter,
                    str
                );

                sendFields<volTensor4Field::Internal>
                (
                    recvProc,
                    dimTensors4,
                    subsetter,
                    str
                );
            }
        }


        UPstream::parRun() = oldParRun;


        // Start sending&receiving from buffers
        pBufs.finishedSends();


        // Subset the part that stays
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        {
            // Save old mesh maps before changing mesh
            const labelList oldFaceOwner(mesh_.faceOwner());
            const labelList oldFaceNeighbour(mesh_.faceNeighbour());
            const label oldInternalFaces = mesh_.nInternalFaces();

            // Remove cells.
            autoPtr<mapPolyMesh> subMap
            (
                doRemoveCells
                (
                    select(false, distribution, Pstream::myProcNo()),
                    oldInternalPatchi
                )
            );

            // Addressing from subsetted mesh
            subCellMap[Pstream::myProcNo()] = subMap().cellMap();
            subFaceMap[Pstream::myProcNo()] = renumber
            (
                repatchFaceMap,
                subMap().faceMap()
            );
            // Insert the sign bit from face flipping
            labelList& faceMap = subFaceMap[Pstream::myProcNo()];
            forAll(faceMap, faceI)
            {
                faceMap[faceI] += 1;
            }
            const labelHashSet& flip = subMap().flipFaceFlux();
            for (const label facei : flip)
            {
                faceMap[facei] = -faceMap[facei];
            }
            subPointMap[Pstream::myProcNo()] = subMap().pointMap();
            subPatchMap[Pstream::myProcNo()] = identity(patches.size());

            // Initialize all addressing into current mesh
            constructCellMap[Pstream::myProcNo()] = identity(mesh_.nCells());
            constructFaceMap[Pstream::myProcNo()] = identity(mesh_.nFaces()) + 1;
            constructPointMap[Pstream::myProcNo()] = identity(mesh_.nPoints());
            constructPatchMap[Pstream::myProcNo()] = identity(patches.size());

            // Subset the mesh data: neighbourCell/neighbourProc
            // fields
            labelList domainSourceFace;
            labelList domainSourceProc;
            labelList domainSourcePatch;
            labelList domainSourceNbrPatch;
            labelList domainSourceNewNbrProc;
            labelList domainSourcePointMaster;

            subsetCouplingData
            (
                mesh_,                          // new mesh
                subMap().pointMap(),            // from new to original mesh
                subMap().faceMap(),             // from new to original mesh
                subMap().cellMap(),

                distribution,                   // distribution before subsetting
                oldFaceOwner,                   // owner before subsetting
                oldFaceNeighbour,               // neighbour        ,,
                oldInternalFaces,               // nInternalFaces   ,,

                sourceFace,
                sourceProc,
                sourcePatch,
                sourceNbrPatch,
                sourceNewNbrProc,
                sourcePointMaster,

                domainSourceFace,
                domainSourceProc,
                domainSourcePatch,
                domainSourceNbrPatch,
                domainSourceNewNbrProc,
                domainSourcePointMaster
            );

            sourceFace.transfer(domainSourceFace);
            sourceProc.transfer(domainSourceProc);
            sourcePatch.transfer(domainSourcePatch);
            sourceNbrPatch.transfer(domainSourceNbrPatch);
            sourceNewNbrProc.transfer(domainSourceNewNbrProc);
            sourcePointMaster.transfer(domainSourcePointMaster);
        }


        // Print a bit.
        if (debug)
        {
            Pout<< nl << "STARTING MESH:" << endl;
            printMeshInfo(mesh_);
            printFieldInfo<volScalarField>(mesh_);
            printFieldInfo<volVectorField>(mesh_);
            printFieldInfo<volSphericalTensorField>(mesh_);
            printFieldInfo<volSymmTensorField>(mesh_);
            printFieldInfo<volTensorField>(mesh_);
            printFieldInfo<surfaceScalarField>(mesh_);
            printFieldInfo<surfaceVectorField>(mesh_);
            printFieldInfo<surfaceSphericalTensorField>(mesh_);
            printFieldInfo<surfaceSymmTensorField>(mesh_);
            printFieldInfo<surfaceTensorField>(mesh_);

            printFieldInfo<volVector1Field>(mesh_);
            printFieldInfo<volVector4Field>(mesh_);
            printFieldInfo<volTensor4Field>(mesh_);
            Pout<< nl << endl;
        }



        // Receive and add what was sent
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        forAll(nSendCells, sendProc)
        {
            // Did processor sendProc send anything to me?
            if
            (
                sendProc != Pstream::myProcNo()
             && nSendCells[sendProc][Pstream::myProcNo()] > 0
            )
            {
                if (debug)
                {
                    Pout<< nl
                        << "RECEIVING FROM DOMAIN " << sendProc
                        << " cells to receive:"
                        << nSendCells[sendProc][Pstream::myProcNo()]
                        << nl << endl;
                }


                // Pstream for receiving mesh and fields
                UIPstream str(sendProc, pBufs);


                // Receive from sendProc
                labelList domainSourceFace;
                labelList domainSourceProc;
                labelList domainSourcePatch;
                labelList domainSourceNbrPatch;
                labelList domainSourceNewNbrProc;
                labelList domainSourcePointMaster;

                autoPtr<fvMesh> domainMeshPtr;

                PtrList<volScalarField> vsf;
                PtrList<volVectorField> vvf;
                PtrList<volSphericalTensorField> vsptf;
                PtrList<volSymmTensorField> vsytf;
                PtrList<volTensorField> vtf;

                PtrList<volVector1Field> vvf1;
                PtrList<volVector4Field> vvf4;
                PtrList<volTensor4Field> vtf4;

                PtrList<surfaceScalarField> ssf;
                PtrList<surfaceVectorField> svf;
                PtrList<surfaceSphericalTensorField> ssptf;
                PtrList<surfaceSymmTensorField> ssytf;
                PtrList<surfaceTensorField> stf;

                PtrList<volScalarField::Internal> dsf;
                PtrList<volVectorField::Internal> dvf;
                PtrList<volSphericalTensorField::Internal> dstf;
                PtrList<volSymmTensorField::Internal> dsytf;
                PtrList<volTensorField::Internal> dtf;

                PtrList<volVector1Field::Internal> dvf1;
                PtrList<volVector4Field::Internal> dvf4;
                PtrList<volTensor4Field::Internal> dtf4;


                // Opposite of sendMesh
                {
                    domainMeshPtr = receiveMesh
                    (
                        sendProc,
                        pointZoneNames,
                        faceZoneNames,
                        cellZoneNames,

                        const_cast<Time&>(mesh_.time()),
                        domainSourceFace,
                        domainSourceProc,
                        domainSourcePatch,
                        domainSourceNbrPatch,
                        domainSourceNewNbrProc,
                        domainSourcePointMaster,
                        str
                    );
                    fvMesh& domainMesh = domainMeshPtr();
                    // Force construction of various on mesh.
                    //(void)domainMesh.globalData();


                    // Receive fields. Read as single dictionary because
                    // of problems reading consecutive fields from single stream.
                    dictionary fieldDicts(str);

                    receiveAllFields
                    (
                        sendProc,
                        domainMesh,
                        fieldDicts,

                        vsf,   volScalars,
                        vvf,   volVectors,
                        vsptf, volSphereTensors,
                        vsytf, volSymmTensors,
                        vtf,   volTensors,
                        vvf1,  volVectors1,
                        vvf4,  volVectors4,
                        vtf4,  volTensors4,
                        ssf,   surfScalars,
                        svf,   surfVectors,
                        ssptf, surfSphereTensors,
                        ssytf, surfSymmTensors,
                        stf,   surfTensors,
                        dsf,   dimScalars,
                        dvf,   dimVectors,
                        dstf,  dimSphereTensors,
                        dsytf, dimSymmTensors,
                        dtf,   dimTensors,
                        dvf1,  dimVectors1,
                        dvf4,  dimVectors4,
                        dtf4,  dimTensors4
                    );
                }
                const fvMesh& domainMesh = domainMeshPtr();


                constructCellMap[sendProc] = identity(domainMesh.nCells());
                constructFaceMap[sendProc] = identity(domainMesh.nFaces()) + 1;
                constructPointMap[sendProc] = identity(domainMesh.nPoints());
                constructPatchMap[sendProc] =
                    identity(domainMesh.boundaryMesh().size());


                // Print a bit.
                if (debug)
                {
                    Pout<< nl << "RECEIVED MESH FROM:" << sendProc << endl;
                    printMeshInfo(domainMesh);
                    printFieldInfo<volScalarField>(domainMesh);
                    printFieldInfo<volVectorField>(domainMesh);
                    printFieldInfo<volSphericalTensorField>(domainMesh);
                    printFieldInfo<volSymmTensorField>(domainMesh);
                    printFieldInfo<volTensorField>(domainMesh);
                    printFieldInfo<surfaceScalarField>(domainMesh);
                    printFieldInfo<surfaceVectorField>(domainMesh);
                    printFieldInfo<surfaceSphericalTensorField>(domainMesh);
                    printFieldInfo<surfaceSymmTensorField>(domainMesh);
                    printFieldInfo<surfaceTensorField>(domainMesh);

                    printFieldInfo<volVector1Field>(domainMesh);
                    printFieldInfo<volVector4Field>(domainMesh);
                    printFieldInfo<volTensor4Field>(domainMesh);
                }


                // Now this mesh we received (from sendProc) needs to be merged
                // with the current mesh. On the current mesh we have for all
                // boundaryfaces the original face and processor. See if we can
                // match these up to the received domainSourceFace and
                // domainSourceProc.
                labelList masterCoupledFaces;
                labelList slaveCoupledFaces;
                findCouples
                (
                    mesh_,

                    sourceFace,
                    sourceProc,
                    sourcePatch,

                    sendProc,
                    domainMesh,
                    domainSourceFace,
                    domainSourceProc,
                    domainSourcePatch,

                    masterCoupledFaces,
                    slaveCoupledFaces
                );

                // Generate additional coupling info (points, edges) from
                // faces-that-match
                faceCoupleInfo couples
                (
                    mesh_,
                    masterCoupledFaces,
                    domainMesh,
                    slaveCoupledFaces,
                    mergeTol_,              // merge tolerance
                    true,                   // faces align
                    true,                   // couples are ordered already
                    false
                );


                // Add domainMesh to mesh
                // ~~~~~~~~~~~~~~~~~~~~~~

                autoPtr<mapAddedPolyMesh> map = fvMeshAdder::add
                (
                    mesh_,
                    domainMesh,
                    couples,
                    false,      // no parallel comms
                    true        // fullyMapped (suppress mapper warning)
                );

                // Update mesh data: sourceFace,sourceProc for added
                // mesh.

                sourceFace = mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceFace,
                    domainMesh.nInternalFaces(),
                    domainSourceFace
                );
                sourceProc = mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceProc,
                    domainMesh.nInternalFaces(),
                    domainSourceProc
                );
                sourcePatch = mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourcePatch,
                    domainMesh.nInternalFaces(),
                    domainSourcePatch
                );
                sourceNbrPatch = mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceNbrPatch,
                    domainMesh.nInternalFaces(),
                    domainSourceNbrPatch
                );
                sourceNewNbrProc = mapBoundaryData
                (
                    mesh_,
                    map(),
                    sourceNewNbrProc,
                    domainMesh.nInternalFaces(),
                    domainSourceNewNbrProc
                );
                // Update pointMaster data
                sourcePointMaster = mapPointData
                (
                    mesh_,
                    map(),
                    sourcePointMaster,
                    domainSourcePointMaster
                );


                // Update all addressing so xxProcAddressing points to correct
                // item in masterMesh.
                const labelList& oldCellMap = map().oldCellMap();
                const labelList& oldFaceMap = map().oldFaceMap();
                const labelList& oldPointMap = map().oldPointMap();
                const labelList& oldPatchMap = map().oldPatchMap();

                //Note: old mesh faces never flipped!
                forAll(constructPatchMap, proci)
                {
                    if (proci != sendProc && constructPatchMap[proci].size())
                    {
                        // Processor already in mesh (either myProcNo or received)
                        inplaceRenumber(oldCellMap, constructCellMap[proci]);
                        inplaceRenumberWithFlip
                        (
                            oldFaceMap,
                            false,
                            true,
                            constructFaceMap[proci]
                        );
                        inplaceRenumber(oldPointMap, constructPointMap[proci]);
                        inplaceRenumber(oldPatchMap, constructPatchMap[proci]);
                    }
                }


                labelHashSet flippedAddedFaces;
                {
                    // Find out if any faces of domain mesh were flipped (boundary
                    // faces becoming internal)
                    const label nBnd = domainMesh.nBoundaryFaces();
                    flippedAddedFaces.resize(nBnd/4);

                    for
                    (
                        label domainFaceI = domainMesh.nInternalFaces();
                        domainFaceI < domainMesh.nFaces();
                        domainFaceI++
                    )
                    {
                        label newFaceI = map().addedFaceMap()[domainFaceI];
                        label newCellI = mesh_.faceOwner()[newFaceI];

                        label domainCellI = domainMesh.faceOwner()[domainFaceI];

                        if (newCellI != map().addedCellMap()[domainCellI])
                        {
                            flippedAddedFaces.insert(domainFaceI);
                        }
                    }
                }


                // Added processor
                inplaceRenumber(map().addedCellMap(), constructCellMap[sendProc]);
                // Add flip
                for (const label domainFaceI : flippedAddedFaces)
                {
                    label& val = constructFaceMap[sendProc][domainFaceI];
                    val = -val;
                }
                inplaceRenumberWithFlip
                (
                    map().addedFaceMap(),
                    false,
                    true,           // constructFaceMap has flip sign
                    constructFaceMap[sendProc]
                );
                inplaceRenumber(map().addedPointMap(), constructPointMap[sendProc]);
                inplaceRenumber(map().addedPatchMap(), constructPatchMap[sendProc]);

                if (debug)
                {
                    Pout<< nl << "MERGED MESH FROM:" << sendProc << endl;
                    printMeshInfo(mesh_);
                    printFieldInfo<volScalarField>(mesh_);
                    printFieldInfo<volVectorField>(mesh_);
                    printFieldInfo<volSphericalTensorField>(mesh_);
                    printFieldInfo<volSymmTensorField>(mesh_);
                    printFieldInfo<volTensorField>(mesh_);
                    printFieldInfo<surfaceScalarField>(mesh_);
                    printFieldInfo<surfaceVectorField>(mesh_);
                    printFieldInfo<surfaceSphericalTensorField>(mesh_);
                    printFieldInfo<surfaceSymmTensorField>(mesh_);
                    printFieldInfo<surfaceTensorField>(mesh_);

                    printFieldInfo<volVector1Field>(mesh_);
                    printFieldInfo<volVector4Field>(mesh_);
                    printFieldInfo<volTensor4Field>(mesh_);
                    Pout<< nl << endl;
                }
            }
        }

        UPstream::parRun() = oldParRun;

        // Print a bit.
        if (debug)
        {
            Pout<< nl << "REDISTRIBUTED MESH:" << endl;
            printMeshInfo(mesh_);
            printFieldInfo<volScalarField>(mesh_);
            printFieldInfo<volVectorField>(mesh_);
            printFieldInfo<volSphericalTensorField>(mesh_);
            printFieldInfo<volSymmTensorField>(mesh_);
            printFieldInfo<volTensorField>(mesh_);
            printFieldInfo<surfaceScalarField>(mesh_);
            printFieldInfo<surfaceVectorField>(mesh_);
            printFieldInfo<surfaceSphericalTensorField>(mesh_);
            printFieldInfo<surfaceSymmTensorField>(mesh_);
            printFieldInfo<surfaceTensorField>(mesh_);

            printFieldInfo<volVector1Field>(mesh_);
            printFieldInfo<volVector4Field>(mesh_);
            printFieldInfo<volTensor4Field>(mesh_);
            Pout<< nl << endl;
        }


        // See if any originally shared points need to be merged. Note: does
        // parallel comms. After this points and edges should again be consistent.
        mergeSharedPoints(sourcePointMaster, mesh_, constructPointMap);


        // Add processorPatches
        // ~~~~~~~~~~~~~~~~~~~~

        // Per neighbour processor, per originating patch, the patchID
        // For faces resulting from internal faces or normal processor patches
        // the originating patch is -1. For cyclics this is the cyclic patchID.
        List<Map<label>> procPatchID;

        // Add processor and processorCyclic patches.
        addProcPatches
        (
            sourceNewNbrProc,
            sourcePatch,
            sourceNbrPatch,
            Pstream::myProcNo(),
            mesh_,
            procPatchID
        );

        // Put faces into correct patch. Note that we now have proper
        // processorPolyPatches again so repatching will take care of coupled face
        // ordering.

        // Get boundary faces to be repatched. Is -1 or new patchID
        labelList newPatchID
        (
            getBoundaryPatch
            (
                sourceNewNbrProc,
                sourcePatch,
                procPatchID,
                Pstream::myProcNo()
            )
        );

        // Change patches. Since this might change ordering of coupled faces
        // we also need to adapt our constructMaps.
        repatch(newPatchID, mesh_, constructFaceMap);

        // Bit of hack: processorFvPatchField does not get reset since created
        // from nothing so explicitly reset.
        initPatchFields<volScalarField, processorFvPatchField<scalar>>
        (
            Zero, mesh_
        );
        initPatchFields<volVectorField, processorFvPatchField<vector>>
        (
            Zero, mesh_
        );
        initPatchFields
        <
            volSphericalTensorField,
            processorFvPatchField<sphericalTensor>
        >
        (
            Zero, mesh_
        );
        initPatchFields<volSymmTensorField, processorFvPatchField<symmTensor>>
        (
            Zero, mesh_
        );
        initPatchFields<volTensorField, processorFvPatchField<tensor>>
        (
            Zero, mesh_
        );

        initPatchFields<volVector1Field, processorFvPatchField<vector1>>
        (
            Zero, mesh_
        );
        initPatchFields<volVector4Field, processorFvPatchField<vector4>>
        (
            Zero, mesh_
        );
        initPatchFields<volTensor4Field, processorFvPatchField<tensor4>>
        (
            Zero, mesh_
        );


        mesh_.setInstance(mesh_.time().timeName());

        // reset fvPatchFields updated_ flag
        resetUpdate<volScalarField>();
        resetUpdate<volVectorField>();
        resetUpdate<volSphericalTensorField>();
        resetUpdate<volSymmTensorField>();
        resetUpdate<volTensorField>();
        resetUpdate<volVector1Field>();
        resetUpdate<volVector4Field>();
        resetUpdate<volTensor4Field>();

        // Print a bit
        if (debug)
        {
            Pout<< nl << "FINAL MESH:" << endl;
            printMeshInfo(mesh_);
            printFieldInfo<volScalarField>(mesh_);
            printFieldInfo<volVectorField>(mesh_);
            printFieldInfo<volSphericalTensorField>(mesh_);
            printFieldInfo<volSymmTensorField>(mesh_);
            printFieldInfo<volTensorField>(mesh_);
            printFieldInfo<surfaceScalarField>(mesh_);
            printFieldInfo<surfaceVectorField>(mesh_);
            printFieldInfo<surfaceSphericalTensorField>(mesh_);
            printFieldInfo<surfaceSymmTensorField>(mesh_);
            printFieldInfo<surfaceTensorField>(mesh_);

            printFieldInfo<volVector1Field>(mesh_);
            printFieldInfo<volVector4Field>(mesh_);
            printFieldInfo<volTensor4Field>(mesh_);
            Pout<< nl << endl;

            OFstream str
                (
                    "mapDistributePolyMesh_proc"
                   +Foam::name(Pstream::myProcNo())
                   +"old"
                );

            str << "mesh_.nCells() " << mesh_.nCells() << nl
                 //<< "nOldPoints " << nOldPoints << nl
                 //<< "nOldFaces " << nOldFaces << nl
                 << "nOldCells " << nOldCells << nl
                 //<< "oldPatchStarts " << oldPatchStarts << nl
                 //<< "oldPatchNMeshPoints " << oldPatchNMeshPoints << nl

                 //<< "subPointMap " << subPointMap << nl
                 //<< "subFaceMap " << subFaceMap << nl
                 << "subCellMap " << subCellMap << nl
                 //<< "subPatchMap " << subPatchMap << nl

                 //<< "constructPointMap " << constructPointMap << nl
                 //<< "constructFaceMap " << constructFaceMap << nl
                 << "constructCellMap " << constructCellMap << nl
                 //<< "constructPatchMap " << constructPatchMap << nl
                 << endl;
        }

        // Collect the distribution maps
        autoPtr<mapDistributePolyMesh> mapPtr
        (
            new mapDistributePolyMesh
            (
                mesh_,

                nOldPoints,
                nOldFaces,
                nOldCells,
                oldPatchStarts.xfer(),
                oldPatchNMeshPoints.xfer(),

                subPointMap.xfer(),
                subFaceMap.xfer(),
                subCellMap.xfer(),
                subPatchMap.xfer(),

                constructPointMap.xfer(),
                constructFaceMap.xfer(),
                constructCellMap.xfer(),
                constructPatchMap.xfer(),

                true,           // subFaceMap has flip
                true            // constructFaceMap has flip
            )
        );

        if (!origMeshConformal)
        {
            myProcMeshes.setSize(1);
            myProcMeshes.set(0, &mesh_);

            mapsDist.setSize(1);
            mapsDist.set(0, mapPtr);

            labelListList meshNoInProc(Pstream::nProcs());
            for (label i = 0; i < Pstream::nProcs(); i++)
            {
                meshNoInProc[i] = labelList(1, i);
            }

            // Add any nonConformalProcessorCyclic patches
            addNccProcPatches(myProcMeshes, meshNoInProc);

            // Unconform the processor meshes
            unconform
            (
                myProcMeshes,
                meshNoInProc,
                mapsDist,
                origPolyFacesBf,
                origFaceOwner,
                distribution
            );

            myProcMeshes.set(0, nullptr).release();

            return autoPtr<mapDistributePolyMesh>(mapsDist.set(0, nullptr));
        }

        return mapPtr;
    } // else Pstream::nProcs() == nFinalProcs_
}


// ************************************************************************* //
