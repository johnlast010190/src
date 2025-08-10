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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshRefinement/meshRefinement.H"
#include "fvMesh/fvMesh.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/Time/Time.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "sets/topoSets/pointSet.H"
#include "sets/topoSets/faceSet.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshTools/meshTools.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyCell.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemoveFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "regionSplit/localPointRegion.H"
#include "polyTopoChange/duplicatePoints/duplicatePoints.H"
#include "regionSplit/regionSplit.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "global/unitConversion/unitConversion.H"
#include "surfaceFormats/obj/OBJstream.H"
#include "meshRefinement/patchFaceOrientation.H"
#include "algorithms/PatchEdgeFaceWave/PatchEdgeFaceWave.H"
#include "algorithms/PatchEdgeFaceWave/patchEdgeFaceRegion.H"
#include "polyMeshAdder/polyMeshAdder.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "snappyHexMeshDriver/refinementParameters/refinementParameters.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fields/volFields/volFields.H"
#include "sets/topoSets/cellSet.H"
#include "primitives/Vector/labelVector/labelVector.H"
#include "meshes/polyMesh/syncTools/dummyTransform.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::meshRefinement::createBaffle
(
    const label facei,
    const label ownPatch,
    const label nbrPatch,
    polyTopoChange& meshMod
) const
{
    const face& f = mesh_.faces()[facei];
    label zoneID = mesh_.faceZones().whichZone(facei);
    bool zoneFlip = false;

    if (zoneID >= 0)
    {
        const faceZone& fZone = mesh_.faceZones()[zoneID];
        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
    }

    meshMod.setAction
    (
        polyModifyFace
        (
            f,                          // modified face
            facei,                      // label of face
            mesh_.faceOwner()[facei],   // owner
            -1,                         // neighbour
            false,                      // face flip
            ownPatch,                   // patch for face
            false,                      // remove from zone
            zoneID,                     // zone for face
            zoneFlip                    // face flip in zone
        )
    );


    label dupFaceI = -1;

    if (mesh_.isInternalFace(facei))
    {
        if (nbrPatch == -1)
        {
            FatalErrorInFunction
                << "No neighbour patch for internal face " << facei
                << " fc:" << mesh_.faceCentres()[facei]
                << " ownPatch:" << ownPatch << abort(FatalError);
        }

        bool reverseFlip = false;
        if (zoneID >= 0)
        {
            reverseFlip = !zoneFlip;
        }

        dupFaceI = meshMod.setAction
        (
            polyAddFace
            (
                f.reverseFace(),            // modified face
                mesh_.faceNeighbour()[facei],// owner
                -1,                         // neighbour
                -1,                         // masterPointID
                -1,                         // masterEdgeID
                facei,                      // masterFaceID,
                true,                       // face flip
                nbrPatch,                   // patch for face
                zoneID,                     // zone for face
                reverseFlip                 // face flip in zone
            )
        );
    }
    return dupFaceI;
}


// Get an estimate for the patch the internal face should get. Bit heuristic.
Foam::label Foam::meshRefinement::getBafflePatch
(
    const labelList& mPatches,
    const labelList& facePatch,
    const label facei
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Loop over face points
    // for each point check all faces patch IDs
    // as soon as an ID >= 0 is found, break and assign that ID
    // to the current face.
    // Check first for real patch (so proper surface intersection and then
    // in facePatch array for patches to block off faces

    label possibleID = -2;
    labelHashSet addedPatchIDSet(mPatches);

    forAll(mesh_.faces()[facei], fp)
    {
        label pointI = mesh_.faces()[facei][fp];

        forAll(mesh_.pointFaces()[pointI], pf)
        {
            label pFaceI = mesh_.pointFaces()[pointI][pf];

            label patchI = patches.whichPatch(pFaceI);

            if
            (
                patchI != -1 && !patches[patchI].coupled()
                && addedPatchIDSet.found(patchI)
            )
            {
                return patchI;
            }
            else if (patchI != -1 && !patches[patchI].coupled())
            {
                possibleID = patchI;
            }

            if
            (
                facePatch[pFaceI] > -1
                && addedPatchIDSet.found(facePatch[pFaceI])
            )
            {
                return facePatch[pFaceI];
            }
            else if (facePatch[pFaceI] > -1)
            {
                possibleID = facePatch[pFaceI];
            }
        }
    }

    // Loop over owner and neighbour cells, looking for the first face with a
    // valid patch number
    const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[facei]];

    forAll(ownFaces, i)
    {
        label cFaceI = ownFaces[i];

        label patchI = patches.whichPatch(cFaceI);

        if
        (
            patchI != -1 && !patches[patchI].coupled()
            && addedPatchIDSet.found(patchI)
        )
        {
            return patchI;
        }
        else if (patchI != -1 && !patches[patchI].coupled() && possibleID == -2)
        {
            possibleID = patchI;
        }

        if (facePatch[cFaceI] > -1 && addedPatchIDSet.found(facePatch[cFaceI]))
        {
            return facePatch[cFaceI];
        }
        else if (facePatch[cFaceI] > -1 && possibleID == -2)
        {
            possibleID = facePatch[cFaceI];
        }
    }

    if (mesh_.isInternalFace(facei))
    {
        const cell& neiFaces = mesh_.cells()[mesh_.faceNeighbour()[facei]];

        forAll(neiFaces, i)
        {
            label cFaceI = neiFaces[i];

            label patchI = patches.whichPatch(cFaceI);

            if
            (
                patchI != -1 && !patches[patchI].coupled()
                && addedPatchIDSet.found(patchI)
            )
            {
                return patchI;
            }
            else if
            (
                patchI != -1 && !patches[patchI].coupled() && possibleID == -2
            )
            {
                possibleID = patchI;
            }


            if
            (
                facePatch[cFaceI] > -1
                && addedPatchIDSet.found(facePatch[cFaceI])
            )
            {
                return facePatch[cFaceI];
            }
            else if (facePatch[cFaceI] > -1  &&  possibleID == -2)
            {
                possibleID = facePatch[cFaceI];
            }
        }
    }

    //Set unmeshed patch otherwise unset (-2) and reset outside of function
    return possibleID;
}


//// Check if we are a boundary face and normal of surface does
//// not align with test vector. In this case there'd probably be
//// a freestanding 'baffle' so we might as well not create it.
//// Note that since it is not a proper baffle we cannot detect it
//// afterwards so this code cannot be merged with the
//// filterDuplicateFaces code.
//bool Foam::meshRefinement::validBaffleTopology
//(
//    const label facei,
//    const vector& n1,
//    const vector& n2,
//    const vector& testDir
//) const
//{
//
//    label patchI = mesh_.boundaryMesh().whichPatch(facei);
//    if (patchI == -1 || mesh_.boundaryMesh()[patchI].coupled())
//    {
//        return true;
//    }
//    else if (mag(n1&n2) > cos(degToRad(30)))
//    {
//        // Both normals aligned. Check that test vector perpendicularish to
//        // surface normal
//        scalar magTestDir = mag(testDir);
//        if (magTestDir > VSMALL)
//        {
//            if (mag(n1&(testDir/magTestDir)) < cos(degToRad(45)))
//            {
//                //Pout<< "** disabling baffling face "
//                //    << mesh_.faceCentres()[facei] << endl;
//                return false;
//            }
//        }
//    }
//    return true;
//}


void Foam::meshRefinement::getIntersections
(
    const pointField& cellCentres,
    const labelList& surfacesToTest,
    const labelList& testFaces,

    labelList& globalRegion1,
    labelList& globalRegion2,
    const bool checkCracks
) const
{
    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(cellCentres, neiLevel, neiCc);

    autoPtr<OBJstream> str;
    if (debug&OBJINTERSECTIONS)
    {
        mkDir(mesh_.time().path()/timeName());
        str.reset
        (
            new OBJstream
            (
                mesh_.time().path()/timeName()/"intersections.obj"
            )
        );

        Pout<< "getIntersections : Writing surface intersections to file "
            << str().name() << nl << endl;
    }

    globalRegion1.setSize(mesh_.nFaces());
    globalRegion1 = -1;
    globalRegion2.setSize(mesh_.nFaces());
    globalRegion2 = -1;


    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    {
        labelList minLevel;
        calcCellCellRays
        (
            cellCentres,
            neiCc,
            labelList(neiCc.size(), -1),
            testFaces,
            start,
            end,
            minLevel
        );
    }

    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;

    if (checkCracks)
    {
        surfaces_.findNearestPerturbedIntersection
        (
            surfacesToTest,
            start,
            end,
            crackTol(),
            addedRays(),

            surface1,
            hit1,
            region1,

            surface2,
            hit2,
            region2
         );
    }
    else
    {
        surfaces_.findNearestIntersection
        (
            surfacesToTest,
            start,
            end,

            surface1,
            hit1,
            region1,

            surface2,
            hit2,
            region2
         );
    }

    labelHashSet unMeshedIDs(unmeshedPatches());

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    boolList checkPatches(patches.size(), false);

    forAll(patches, patchI)
    {
        if (unMeshedIDs.found(patchI))
        {
            checkPatches[patchI] = true;
        }
    }
    DynamicList<label> checkFaces(mesh_.nFaces()/100);

    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        if (hit1[i].hit() && hit2[i].hit())
        {
            label patchI = patches.whichPatch(facei);
            if (patchI != -1 && checkPatches[patchI])
            {
                checkFaces.append(i);
            }

            if (str.valid())
            {
                str().write(linePointRef(start[i], hit1[i].rawPoint()));
                str().write
                (
                    linePointRef(hit1[i].rawPoint(), hit2[i].rawPoint())
                );
                str().write(linePointRef(hit2[i].rawPoint(), end[i]));
            }

            // Pick up the patches
            globalRegion1[facei] =
                surfaces_.globalRegion(surface1[i], region1[i]);
            globalRegion2[facei] =
                surfaces_.globalRegion(surface2[i], region2[i]);

            if (globalRegion1[facei] == -1 || globalRegion2[facei] == -1)
            {
                FatalErrorInFunction
                    << "problem." << abort(FatalError);
            }
        }
    }

    if
    (
        !controller_.dual()
        && returnReduce(checkFaces.size(), sumOp<label>()) > 0
    )
    {
        pointField sNormals(checkFaces.size(), vector(GREAT, GREAT, GREAT));

        forAll(surfacesToTest, sI)
        {
            DynamicList<pointIndexHit> localHits;
            label surfI = surfacesToTest[sI];

            forAll(checkFaces, cFI)
            {
                label i = checkFaces[cFI];
                if (surface1[i] == surfI)
                {
                    localHits.append(hit1[i]);
                }
            }
            localHits.shrink();

            pointField localNormals;
            label geomI = surfaces_.surfaces()[surfI];
            surfaces_.geometry()[geomI].getNormal
            (
                localHits,
                localNormals
            );

            label localI = 0;
            forAll(checkFaces, cFI)
            {
                label i = checkFaces[cFI];
                if (surface1[i] == surfI)
                {
                    sNormals[cFI] = localNormals[localI];
                    localI++;
                }
            }
        }

        label nReset = 0;
        forAll(checkFaces, cFI)
        {
            label facei = testFaces[checkFaces[cFI]];
            vector fN = mesh_.faceAreas()[facei];
            fN /= (mag(fN) + SMALL);

            if (mag(fN & sNormals[cFI]) < 0.707)
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
                nReset++;
            }
        }

        Info<< "Prevented resetting  of "
             << returnReduce(nReset, sumOp<label>())
             << " faces on blockMesh patches because not aligned with surface."
             << endl;
    }

    syncTools::syncFaceList
    (
        mesh_,
        globalRegion1,
        maxEqOp<label>()
    );

    syncTools::syncFaceList
    (
        mesh_,
        globalRegion2,
        maxEqOp<label>()
    );
}


void Foam::meshRefinement::keepBoundaryPointCells
(
    const labelList& testFaces,
    labelList& ownPatch,
    labelList& nbrPatch
)
{
    Info<<"Checking whether to keep boundary point cells"<<endl;

    boolList blockedFace(mesh_.nFaces(),false);
    forAll(testFaces, i)
    {
        label facei = testFaces[i];
        if (ownPatch[facei] != -1 || nbrPatch[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);

    labelHashSet keepRegionSet(keepLargestRegions(cellRegion));
    labelList cellToZone(mesh_.nCells(), -1);

    boolList markedCells(mesh_.nCells(), false);
    forAll(mesh_.cells(), celli)
    {
        if (!keepRegionSet.found(cellRegion[celli]))
        {
            markedCells[celli] = true;
            cellToZone[celli] = -2;
        }
    }

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    boolList neiMarkedCells;
    syncTools::swapBoundaryCellList(mesh_, markedCells, neiMarkedCells);

    blockedFace = false;
//    label nBlocked = 0;
    forAll(blockedFace,facei)
    {
        bool markedOwn = markedCells[mesh_.faceOwner()[facei]];
        bool markedNei = false;
        if (mesh_.isInternalFace(facei))
        {
            markedNei = markedCells[mesh_.faceNeighbour()[facei]];
        }
        else
        {
            markedNei = neiMarkedCells[facei-mesh_.nInternalFaces()];
        }

        if
        (
            (markedOwn && !markedNei)
            || (!markedOwn && markedNei)
        )
        {
            blockedFace[facei] = true;
//            nBlocked++;
        }
    }

    boolList blockedPoints(mesh_.nPoints(),false);
    forAll(mesh_.faces(), facei)
    {
        if (blockedFace[facei])
        {
            const face& f = mesh_.faces()[facei];
            forAll(f,fp)
            {
                label pointI = f[fp];
                if (!blockedPoints[pointI])
                {
                    blockedPoints[pointI] = true;
                }
            }
        }
    }

    // Sync
    syncTools::syncPointList
    (
        mesh_,
        blockedPoints,
        orEqOp<bool>(),
        false           // null value
    );

    label nChecks = 0;
    forAll(mesh_.points(), pointI)
    {
        if (blockedPoints[pointI])
        {
            const labelList& pCells = mesh_.pointCells()[pointI];
            forAll(pCells, pCI)
            {
                if (!markedCells[pCells[pCI]])
                {
                    nChecks++;
                }
            }
        }
    }

    //check point to cc intersection
    pointField start(nChecks);
    pointField end(nChecks);
    nChecks = 0;
    forAll(mesh_.points(), pointI)
    {
        if (blockedPoints[pointI])
        {
            const labelList& pCells = mesh_.pointCells()[pointI];
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                if (!markedCells[pCells[pCI]])
                {
                    start[nChecks] = mesh_.points()[pointI];
                    end[nChecks] = mesh_.cellCentres()[celli];
                    nChecks++;
                }
            }
        }
    }

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;

    const labelList& unamedAndBoundarySurfaces =
        surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
        (
            surfaces_.surfZones()
        );

    surfaces_.findNearestIntersection
    (
        unamedAndBoundarySurfaces,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
    );

    nChecks = 0;

    boolList removePoints(mesh_.nPoints(),false);
    forAll(mesh_.points(), pointI)
    {
        if (blockedPoints[pointI])
        {
            const labelList& pCells = mesh_.pointCells()[pointI];
            bool foundHit = false;
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                if (!markedCells[celli])
                {
                    if (hit1[nChecks].hit() && hit2[nChecks].hit())
                    {
                        foundHit = true;
                    }
                    nChecks++;
                }
            }
            if (foundHit)
            {
                removePoints[pointI] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        removePoints,
        orEqOp<bool>(),
        false           // null value
    );

    PackedBoolList keepCells(mesh_.nCells(),false);
//    label nKeep = 0;
    forAll(mesh_.points(), pointI)
    {
        if (blockedPoints[pointI] && !removePoints[pointI])
        {
            const labelList& pCells = mesh_.pointCells()[pointI];
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                if (markedCells[celli] && !keepCells[celli])
                {
                    keepCells[celli] = true;
//                    nKeep++;
                }
            }
        }
    }

    boolList unsetFaces(mesh_.nFaces(), false);
    markedCells = false;

    bool firstPass = true;
    labelList newCellToZone(mesh_.nCells(), -2);
    while (true)
    {
        label nSet = 0;

        forAll(mesh_.cells(), celli)
        {
            if (cellToZone[celli] < -1 && keepCells[celli])
            {
                const labelList& cFaces = mesh_.cells()[celli];
                label gReg1 = -1;
                label gReg2 = -1;

                forAll(cFaces, cFI)
                {
                    label facei = cFaces[cFI];
                    label own = mesh_.faceOwner()[facei];
                    label neiZone = -1;

                    if (own == celli)
                    {
                        if (mesh_.isInternalFace(facei))
                        {
                            neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
                        }
                        else
                        {
                            neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                        }
                    }
                    else
                    {
                        neiZone = cellToZone[own];
                    }

                    if (neiZone >= -1)
                    {
                        if (!firstPass)
                        {
                            cellToZone[celli] = neiZone;
                        }
                        else
                        {
                            newCellToZone[celli] = neiZone;
                        }
                        markedCells[celli] = true;
                        unsetFaces[facei] = true;
                        gReg1 = ownPatch[facei];
                        gReg2 = nbrPatch[facei];
                        nSet++;
                        break;
                    }
                }
                if (gReg1 != -1)
                {
                    forAll(cFaces, cFI)
                    {
                        label facei = cFaces[cFI];
                        if (blockedFace[facei])
                        {
                            unsetFaces[facei] = true;
                        }
                        else if (ownPatch[facei] == -1)
                        {
                            ownPatch[facei] = gReg1;
                            nbrPatch[facei] = gReg2;
                        }
                    }
                }
            }
        }

        if (firstPass || returnReduce(nSet, orOp<label>()) == 0)
        {
            if (firstPass)
            {
                firstPass = false;
                forAll(newCellToZone,celli)
                {
                    if (newCellToZone[celli] != -2)
                    {
                        cellToZone[celli] = newCellToZone[celli];
                    }
                }
            }
            else
            {
                break;
            }
        }

        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
        syncTools::syncFaceList(mesh_, nbrPatch, maxEqOp<label>());
    }

    syncTools::syncFaceList(mesh_, unsetFaces, orEqOp<bool>());

    syncTools::swapBoundaryCellList(mesh_, markedCells, neiMarkedCells);

    forAll(mesh_.faces(), facei)
    {
        if (unsetFaces[facei])
        {
            ownPatch[facei] = -1;
            nbrPatch[facei] = -1;
        }
        else
        {
            bool markedOwn = markedCells[mesh_.faceOwner()[facei]];
            bool markedNei = false;
            if (mesh_.isInternalFace(facei))
            {
                markedNei = markedCells[mesh_.faceNeighbour()[facei]];
            }
            else
            {
                markedNei = neiMarkedCells[facei-mesh_.nInternalFaces()];
            }

            if
            (
                markedOwn && markedNei
            )
            {
                ownPatch[facei] = -1;
                nbrPatch[facei] = -1;
            }
        }
    }

    syncTools::syncFaceList(mesh_, ownPatch, minEqOp<label>());
    syncTools::syncFaceList(mesh_, nbrPatch, minEqOp<label>());
}

void Foam::meshRefinement::getBafflePatches
(
    const refinementParameters& refineParams,
    const labelList& globalToMasterPatch,
    labelList& ownPatch,
    labelList& nbrPatch
) const
{
    // This determines the patches for the intersected faces to
    // - remove the outside of the mesh
    // - introduce baffles for (non-faceZone) intersections
    // Any baffles for faceZones (faceType 'baffle'/'boundary') get introduced
    // later

    // 1. Determine cell zones
    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Note that this does not determine the surface+region that was intersected
    // so that is done in step 2 below.

    //re-use intersection of unnamed surfaces calculated in zonify
    labelList testFaces;
    labelList globalRegion1;
    labelList globalRegion2;

    const label minZoneRegionSize = refineParams.minZoneRegionSize();
    const Switch namedLocationsRezone = refineParams.namedLocationsRezone();

    labelList cellToZone(mesh_.nCells(),-2);
    {
        labelList namedSurfaceIndex;  //filtered named intersections
        labelList namedIntersections; //all named intersections

        PackedBoolList posOrientation;
        zonify
        (
            -2,                 // zone to put unreached cells into
            refineParams,
            (
                namedLocationsRezone
                ? minZoneRegionSize
                : label(0)
            ),
            testFaces,
            globalRegion1,
            globalRegion2,

            cellToZone,
            namedSurfaceIndex,
            namedIntersections,
            posOrientation
        );
    }

    // The logic is quite complicated and depends on the cell zone allocation
    // (see zonify). Cells can have zone:
    //  -2  : unreachable
    //  -1  : in background zone
    //  >=0 : in named cellZone
    // Faces can be intersected by a
    //  - unnamed surface (no faceZone defined for it)
    //  - named surface (a faceZone defined for it)
    // Per intersected faces, depending on the cellToZone on either side of
    // the face we need to:
    //
    // surface type    |   cellToZone      |   action
    // ----------------+-------------------+---------
    // unnamed         |  -2  | same       |   -
    // unnamed         |  -2  | different  |   baffle
    //                 |      |            |
    // unnamed         |  -1  | same       |   baffle
    // unnamed         |  -1  | different  |   -
    //                 |      |            |
    // unnamed         | >=0  | same       |   baffle
    // unnamed         | >=0  | different  |   -
    //
    // named           |  -2  | same       |   -
    // named           |  -2  | different  |   see note
    //
    // named           |  -1  | same       |   -
    // named           |  -1  | different  |   -
    //
    // named           | >=0  | same       |   -
    // named           | >=0  | different  |   -
    //
    // So the big difference between surface with a faceZone and those
    // without is that 'standing-up-baffles' are not supported. Note that
    // they are still in a faceZone so can be split etc. later on.
    // Note: this all depends on whether we allow named surfaces
    //       to be outside the unnamed geometry. 2.3 does not allow this
    //       so we do the same. We could implement it but it would require
    //       re-testing of the intersections with the named surfaces to
    //       obtain the surface and region number.

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);


    // 3. Baffle all boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Baffle all boundary faces except those on outside of unvisited cells
    // Use patch according to globalRegion1,2
    // Note: the behaviour is slightly different from version3.0 and earlier
    //       in that it will not create baffles which are outside the meshed
    //       domain. On a medium testcase (motorBike tutorial geometry) this
    //       selects about 20 cells less (out of 120000). These cells are where
    //       there might e.g. be a single cell which is fully unreachable.

    ownPatch.setSize(mesh_.nFaces());
    ownPatch = -1;
    nbrPatch.setSize(mesh_.nFaces());
    nbrPatch = -1;
    label nFreeStanding = 0;

    labelList wrapSurfaces(surfaces().wrapLevelSurfaces());
    label wrapPatchID = -1;
    if (wrapSurfaces.size() > 0)
    {
        dictionary patchInfo;
        patchInfo.set("type", wallPolyPatch::typeName);
        wrapPatchID = addPatch(mesh_, "wrapped", patchInfo);
    }

    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const wordList& zonesInMesh = refineParams.zonesInMesh();

    //Check for location named zones
    DynamicList<label> namedLocations(locationsInMesh.size());
    forAll(locationsInMesh, i)
    {
        const word& name = zonesInMesh[i];
        if (name != "none")
        {
            label zoneID = mesh_.cellZones().findZoneID(name);
            if (zoneID != -1)
            {
                namedLocations.append(zoneID);
            }
        }
    }
    bool interZoneBaffles = false;
    if
    (
        namedLocations.size() < 2
        || refineParams.interZoneBaffles()
    )
    {
        interZoneBaffles = true;
    }

    labelHashSet namedLocationsSet(namedLocations);

    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
        {
            label ownMasterPatch = -1;
            bool wrapFace = false;
            if (globalRegion1[facei] != -1)
            {
                if (wrapPatchID != -1 && wrapIndex_[facei] != -1)
                {
                    wrapFace = true;
                    ownMasterPatch = wrapPatchID;
                }
                else
                {
                    ownMasterPatch = globalToMasterPatch[globalRegion1[facei]];
                }
            }
            label neiMasterPatch = -1;
            if (globalRegion2[facei] != -1)
            {
                if (wrapPatchID != -1 && wrapIndex_[facei] != -1)
                {
                    wrapFace = true;
                    neiMasterPatch = wrapPatchID;
                }
                else
                {
                    neiMasterPatch = globalToMasterPatch[globalRegion2[facei]];
                }
            }

            label ownZone = cellToZone[mesh_.faceOwner()[facei]];
            label neiZone = -1;

            if (mesh_.isInternalFace(facei))
            {
                neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
            }
            else
            {
                neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
            }

            if (ownZone == -2)
            {
                if (neiZone != -2)
                {
                    ownPatch[facei] = ownMasterPatch;
                    nbrPatch[facei] = neiMasterPatch;
                }
            }
            else if (neiZone == -2)
            {
                ownPatch[facei] = ownMasterPatch;
                nbrPatch[facei] = neiMasterPatch;
            }
            else if (ownZone == neiZone && !wrapFace)
            {
                // Free-standing baffle
                ownPatch[facei] = ownMasterPatch;
                nbrPatch[facei] = neiMasterPatch;
                nFreeStanding++;
            }
            else if (ownZone >= 0)
            {
                if
                (
                    interZoneBaffles
                    || (neiZone == -1 && !namedLocationsSet.found(ownZone))
                )
                {
                    ownPatch[facei] = ownMasterPatch;
                    nbrPatch[facei] = neiMasterPatch;
                }
            }
            else if (neiZone >= 0)
            {
                if
                (
                    interZoneBaffles
                    || (ownZone == -1 && !namedLocationsSet.found(neiZone))
                )
                {
                    ownPatch[facei] = ownMasterPatch;
                    nbrPatch[facei] = neiMasterPatch;
                }
            }
        }
    }

    // No need to parallel sync since intersection data (surfaceIndex_ etc.)
    // already guaranteed to be synced...
    // However:
    // - owncc and neicc are reversed on different procs so might pick
    //   up different regions reversed? No problem. Neighbour on one processor
    //   might not be owner on the other processor but the neighbour is
    //   not used when creating baffles from proc faces.
    // - tolerances issues occasionally crop up.
    syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    syncTools::syncFaceList(mesh_, nbrPatch, maxEqOp<label>());

    // Recheck free standing faces if present with no crack tests
    if
    (
        checkForCracks()
        && returnReduce(nFreeStanding, sumOp<label>()) != 0
    )
    {
        //Mark all removal cells to -2
        boolList blockedFace(mesh_.nFaces(),false);

        forAll(blockedFace, facei)
        {
            if
            (
                ownPatch[facei] != -1 || nbrPatch[facei] != -1
            )
            {
                 blockedFace[facei] = true;
            }
        }

        // Set region per cell based on walking
        regionSplit cellRegion(mesh_, blockedFace);
        blockedFace.clear();

        // Force calculation of face decomposition (used in findCell)
        (void)mesh_.tetBasePtIs();

        DynamicList<label> keepRegions(cellRegion.nRegions());
        // For all locationsInMesh find the cell
        forAll(locationsInMesh, i)
        {
            const point& insidePoint = locationsInMesh[i];

            label keepRegionI = -1;

            label celli = findCell
            (
                insidePoint,
                mesh_,
                meshCutter_
            );

            if (celli != -1)
            {
                keepRegionI = cellRegion[celli];
            }
            reduce(keepRegionI, maxOp<label>());
            keepRegions.append(keepRegionI);
        }

        labelHashSet keepRegionsSet(keepRegions);

        forAll(cellToZone, celli)
        {
            if (cellToZone[celli] != -2)
            {
                label regionI = cellRegion[celli];
                if (!keepRegionsSet.found(regionI))
                {
                    cellToZone[celli] = -2;
                }
            }
        }
        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

        labelList freeStandingFaces(nFreeStanding);
        nFreeStanding = 0;
        forAll(testFaces, i)
        {
            label facei = testFaces[i];

            label ownZone = cellToZone[mesh_.faceOwner()[facei]];
            label neiZone = -1;

            if (mesh_.isInternalFace(facei))
            {
                neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
            }
            else
            {
                neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
            }

            if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
            {
                bool wrapFace = false;
                if (wrapPatchID != -1 && wrapIndex_[facei] != -1)
                {
                    wrapFace = true;
                }

                if
                (
                    ownZone != -2 && neiZone != -2
                    && ownZone == neiZone && !wrapFace
                )
                {
                    freeStandingFaces[nFreeStanding++] = facei;
                }
            }
        }
        freeStandingFaces.setSize(nFreeStanding);

        labelList globalRegionBaffle1;
        labelList globalRegionBaffle2;
        getIntersections
        (
            mesh_.cellCentres(),
            surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
            (
                surfaces_.surfZones()
            ),
            freeStandingFaces,
            globalRegionBaffle1,
            globalRegionBaffle2,
            false
        );

        forAll(freeStandingFaces, i)
        {
            label facei = freeStandingFaces[i];
            if (globalRegionBaffle1[facei] == -1 || globalRegionBaffle2[facei] == -1)
            {
                ownPatch[facei] = -1;
                nbrPatch[facei] = -1;
            }
        }
        syncTools::syncFaceList(mesh_, ownPatch, minEqOp<label>());
        syncTools::syncFaceList(mesh_, nbrPatch, minEqOp<label>());
    }
}


void Foam::meshRefinement::updateBaffles
(
    const mapPolyMesh &map,
    List<labelPair>& baffles,
    const List<labelPair> added
)
{
    forAll(baffles, i)
    {
        label& face0 = baffles[i].first();
        label& face1 = baffles[i].second();

        if (face0 != -1)
        {
            face0 = map.reverseFaceMap()[face0];
        }
        if (face1 != -1)
        {
            face1 = map.reverseFaceMap()[face1];
        }
    }

    if (added.size())
    {
        baffles.append(added);
    }
}


Foam::Map<Foam::labelPair> Foam::meshRefinement::getZoneBafflePatches
(
    const bool allowBoundary,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
) const
{
    Map<labelPair> bafflePatch(mesh_.nFaces()/1000);

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    const faceZoneMesh& fZones = mesh_.faceZones();

    forAll(surfZones, surfI)
    {
        const word& faceZoneName = surfZones[surfI].faceZoneName();

        if (faceZoneName.size())
        {
            // Get zone
            label zoneI = fZones.findZoneID(faceZoneName);

            const faceZone& fZone = fZones[zoneI];

            // Get patch allocated for zone
            label globalRegionI = surfaces_.globalRegion(surfI, 0);
            labelPair zPatches
            (
                globalToMasterPatch[globalRegionI],
                globalToSlavePatch[globalRegionI]
            );

            Info<< "For zone " << fZone.name() << " found patches "
                << mesh_.boundaryMesh()[zPatches[0]].name() << " and "
                << mesh_.boundaryMesh()[zPatches[1]].name()
                << endl;

            forAll(fZone, i)
            {
                label facei = fZone[i];

                if (allowBoundary || mesh_.isInternalFace(facei))
                {
                    labelPair patches = zPatches;
                    if (fZone.flipMap()[i])
                    {
                       patches = reverse(patches);
                    }

                    if (!bafflePatch.insert(facei, patches))
                    {
                        FatalErrorInFunction
                            << "Face " << facei
                            << " fc:" << mesh_.faceCentres()[facei]
                            << " in zone " << fZone.name()
                            << " is in multiple zones!"
                            << abort(FatalError);
                    }
                }
            }
        }
    }
    return bafflePatch;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createBaffles
(
    const labelList& ownPatch,
    const labelList& nbrPatch,
    const bool updateIntersection /*= true*/
)
{
    if
    (
        ownPatch.size() != mesh_.nFaces()
     || nbrPatch.size() != mesh_.nFaces()
    )
    {
        FatalErrorInFunction
            << "Illegal size :"
            << " ownPatch:" << ownPatch.size()
            << " nbrPatch:" << nbrPatch.size()
            << ". Should be number of faces:" << mesh_.nFaces()
            << abort(FatalError);
    }

    if (debug)
    {
        labelList syncedOwnPatch(ownPatch);
        syncTools::syncFaceList(mesh_, syncedOwnPatch, maxEqOp<label>());
        labelList syncedNeiPatch(nbrPatch);
        syncTools::syncFaceList(mesh_, syncedNeiPatch, maxEqOp<label>());

        forAll(syncedOwnPatch, facei)
        {
            if
            (
                (ownPatch[facei] == -1 && syncedOwnPatch[facei] != -1)
             || (nbrPatch[facei] == -1 && syncedNeiPatch[facei] != -1)
            )
            {
                FatalErrorInFunction
                    << "Non synchronised at face:" << facei
                    << " on patch:" << mesh_.boundaryMesh().whichPatch(facei)
                    << " fc:" << mesh_.faceCentres()[facei] << endl
                    << "ownPatch:" << ownPatch[facei]
                    << " syncedOwnPatch:" << syncedOwnPatch[facei]
                    << " nbrPatch:" << nbrPatch[facei]
                    << " syncedNeiPatch:" << syncedNeiPatch[facei]
                    << abort(FatalError);
            }
        }
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nBaffles = 0;

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (ownPatch[facei] != -1)
        {
            // Create baffle or repatch face. Return label of inserted baffle
            // face.
            createBaffle
            (
                facei,
                ownPatch[facei],   // owner side patch
                nbrPatch[facei],   // neighbour side patch
                meshMod
            );
            nBaffles++;
        }
    }
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        label coupledPatchI = -1;
        if
        (
            pp.coupled()
        && !refCast<const coupledPolyPatch>(pp).transform().translates()
        )
        {
            coupledPatchI = patchI;
        }

        forAll(pp, i)
        {
            label facei = pp.start()+i;

            if (ownPatch[facei] != -1)
            {
                createBaffle
                (
                    facei,
                    ownPatch[facei],    // owner side patch
                    nbrPatch[facei],    // neighbour side patch
                    meshMod
                );

                if (coupledPatchI != -1)
                {
                    faceToCoupledPatch_.insert(facei, coupledPatchI);
                }

                nBaffles++;
            }
        }
    }


    autoPtr<mapPolyMesh> map;
    if (returnReduce(nBaffles, sumOp<label>()))
    {
        // Change the mesh (no inflation, parallel sync)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }


        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        //- Redo the intersections on the newly create baffle faces. Note that
        //  this changes also the cell centre positions.
        faceSet baffledFacesSet(mesh_, "baffledFacesSet", 2*nBaffles);

        const labelList& reverseFaceMap = map().reverseFaceMap();
        const labelList& faceMap = map().faceMap();

        // Pick up owner side of baffle
        forAll(ownPatch, oldFaceI)
        {
            label facei = reverseFaceMap[oldFaceI];

            if (ownPatch[oldFaceI] != -1 && facei >= 0)
            {
                const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[facei]];

                forAll(ownFaces, i)
                {
                    baffledFacesSet.insert(ownFaces[i]);
                }
            }
        }
        // Pick up neighbour side of baffle (added faces)
        forAll(faceMap, facei)
        {
            label oldFaceI = faceMap[facei];

            if (oldFaceI >= 0 && reverseFaceMap[oldFaceI] != facei)
            {
                const cell& ownFaces = mesh_.cells()[mesh_.faceOwner()[facei]];

                forAll(ownFaces, i)
                {
                    baffledFacesSet.insert(ownFaces[i]);
                }
            }
        }
        baffledFacesSet.sync(mesh_);

        updateMesh(map, baffledFacesSet.toc(), updateIntersection);
    }

    return map;
}


Foam::labelList Foam::meshRefinement::getZones
(
    const List<surfaceZonesInfo::faceZoneType>& fzTypes
) const
{
    const faceZoneMesh& faceZones = mesh_.faceZones();

    DynamicList<label> zoneIDs(faceZones.size());

    forAll(faceZones, zoneI)
    {
        const faceZone& fZone = faceZones[zoneI];

        label mpI, spI;
        surfaceZonesInfo::faceZoneType fzType;
        bool hasInfo = getFaceZoneInfo(fZone.name(), mpI, spI, fzType);

        if (hasInfo && findIndex(fzTypes, fzType) != -1)
        {
            zoneIDs.append(zoneI);
        }
    }
    return labelList(zoneIDs, true);
}


// Subset those baffles where both faces are on the same zone
Foam::List<Foam::labelPair> Foam::meshRefinement::subsetBaffles
(
    const polyMesh& mesh,
    const labelList& zoneIDs,
    const List<labelPair>& baffles
)
{
    const faceZoneMesh& faceZones = mesh.faceZones();

    // Mark zone per face
    labelList faceToZone(mesh.nFaces(), -1);

    forAll(zoneIDs, i)
    {
        label zoneID = zoneIDs[i];
        UIndirectList<label>(faceToZone, faceZones[zoneID]) = zoneID;
    }

    // Subset baffles
    DynamicList<labelPair> newBaffles(baffles.size());
    forAll(baffles, i)
    {
        const labelPair& p = baffles[i];
        label ftz0 = faceToZone[p[0]];
        label ftz1 = faceToZone[p[1]];
        if (ftz0 != -1 && (ftz0 == ftz1))
        {
            newBaffles.append(p);
        }
    }

    return List<labelPair>(newBaffles, true);
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::createZoneBaffles
(
    const labelList& zoneIDs,
    List<labelPair>& baffles,
    labelList& originatingFaceZone,
    const bool updateIntersection
)
{
    autoPtr<mapPolyMesh> map;

    if (zoneIDs.size() > 0)
    {
        const faceZoneMesh& faceZones = mesh_.faceZones();

        // Split internal faces on interface surfaces
        Info<< "Converting zoned faces into baffles ..." << endl;

        // Per (internal) face the patch it should go into
        labelList ownPatch(mesh_.nFaces(), -1);
        labelList nbrPatch(mesh_.nFaces(), -1);
        labelList faceZoneID(mesh_.nFaces(), -1);

        labelList nBaffles(zoneIDs.size(), 0);

        forAll(zoneIDs, j)
        {
            label zoneI = zoneIDs[j];
            const faceZone& fz = faceZones[zoneI];

            const word& masterName = faceZoneToMasterPatch_[fz.name()];
            label masterPatchI = mesh_.boundaryMesh().findPatchID(masterName);
            const word& slaveName = faceZoneToSlavePatch_[fz.name()];
            label slavePatchI = mesh_.boundaryMesh().findPatchID(slaveName);

            if (masterPatchI == -1 || slavePatchI == -1)
            {
                FatalErrorInFunction
                    << "Problem: masterPatchI:" << masterPatchI
                    << " slavePatchI:" << slavePatchI << exit(FatalError);
            }

            forAll(fz, i)
            {
                label facei = fz[i];

                if (mesh_.isInternalFace(facei))
                {
                    if (fz.flipMap()[i])
                    {
                        ownPatch[facei] = slavePatchI;
                        nbrPatch[facei] = masterPatchI;
                    }
                    else
                    {
                        ownPatch[facei] = masterPatchI;
                        nbrPatch[facei] = slavePatchI;
                    }
                    faceZoneID[facei] = zoneI;

                    nBaffles[j]++;
                }
            }
        }

        label nLocalBaffles = sum(nBaffles);


        label nTotalBaffles = returnReduce(nLocalBaffles, sumOp<label>());

        if (nTotalBaffles > 0)
        {
            Pstream::listCombineReduce(nBaffles, plusOp<label>());

            Info<< nl
                << setf(ios_base::left)
                << setw(30) << "FaceZone"
                << setw(10) << "FaceType"
                << setw(10) << "nBaffles"
                << nl
                << setw(30) << "--------"
                << setw(10) << "--------"
                << setw(10) << "--------"
                << endl;

            forAll(zoneIDs, j)
            {
                label zoneI = zoneIDs[j];
                const faceZone& fz = faceZones[zoneI];

                label mpI, spI;
                surfaceZonesInfo::faceZoneType fzType;
                bool hasInfo = getFaceZoneInfo(fz.name(), mpI, spI, fzType);

                if (hasInfo)
                {
                    Info<< setf(ios_base::left)
                        << setw(30) << fz.name()
                        << setw(10)
                        << surfaceZonesInfo::faceZoneTypeNames[fzType]
                        << setw(10) << nBaffles[j]
                        << nl;
                }
            }
            Info<< endl;

            // Create baffles.
            map = createBaffles(ownPatch, nbrPatch, updateIntersection);

            // Get pairs of faces created.
            // Just loop over faceMap and store baffle if we encounter a slave
            // face.

            baffles.setSize(nLocalBaffles);
            originatingFaceZone.setSize(nLocalBaffles);
            label baffleI = 0;

            const labelList& faceMap = map().faceMap();
            const labelList& reverseFaceMap = map().reverseFaceMap();

            for
            (
                label facei = mesh_.nInternalFaces();
                facei < mesh_.nFaces();
                facei++
            )
            {
                label oldFaceI = faceMap[facei];
                label masterFaceI = reverseFaceMap[oldFaceI];
                if (masterFaceI != facei && ownPatch[oldFaceI] != -1)
                {
                    baffles[baffleI] = labelPair(masterFaceI, facei);
                    originatingFaceZone[baffleI] = faceZoneID[oldFaceI];
                    baffleI++;
                }
            }

            if (baffleI != baffles.size())
            {
                FatalErrorInFunction
                    << "Had " << baffles.size() << " baffles to create "
                    << " but encountered " << baffleI
                    << " slave faces originating from patcheable faces."
                    << abort(FatalError);
            }

            if (debug&MESH)
            {
                const_cast<Time&>(mesh_.time())++;
                Pout<< "Writing zone-baffled mesh to time " << timeName()
                    << endl;
                write
                (
                    debugType(debug),
                    writeType(writeLevel() | WRITEMESH),
                    mesh_.time().path()/"baffles"
                );
            }
        }
        Info<< "Created " << nTotalBaffles << " baffles in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
    else
    {
        baffles.clear();
        originatingFaceZone.clear();
    }

    return map;
}


Foam::List<Foam::labelPair> Foam::meshRefinement::freeStandingBaffles
(
    const List<labelPair>& couples,
    const scalar planarAngle
) const
{
    // Done by counting the number of baffles faces per mesh edge. If edge
    // has 2 boundary faces and both are baffle faces it is the edge of a baffle
    // region.

    // All duplicate faces on edge of the patch are to be merged.
    // So we count for all edges of duplicate faces how many duplicate
    // faces use them.
    labelList nBafflesPerEdge(mesh_.nEdges(), 0);


    // This algorithm is quite tricky. We don't want to use edgeFaces and
    // also want it to run in parallel so it is now an algorithm over
    // all (boundary) faces instead.
    // We want to pick up any edges that are only used by the baffle
    // or internal faces but not by any other boundary faces. So
    // - increment count on an edge by 1 if it is used by any (uncoupled)
    //   boundary face.
    // - increment count on an edge by 1000000 if it is used by a baffle face
    // - sum in parallel
    //
    // So now any edge that is used by baffle faces only will have the
    // value 2*1000000+2*1.


    const label baffleValue = 1000000;


    // Count number of boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        // Count number of boundary faces. Discard coupled boundary faces.
        if (!pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                const labelList& fEdges = mesh_.faceEdges(facei);

                forAll(fEdges, fEdgeI)
                {
                    nBafflesPerEdge[fEdges[fEdgeI]]++;
                }
                facei++;
            }
        }
    }


    DynamicList<label> fe0;
    DynamicList<label> fe1;


    // Count number of duplicate boundary faces per edge
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(couples, i)
    {
        {
            label f0 = couples[i].first();
            const labelList& fEdges0 = mesh_.faceEdges(f0, fe0);
            forAll(fEdges0, fEdgeI)
            {
                nBafflesPerEdge[fEdges0[fEdgeI]] += baffleValue;
            }
        }

        {
            label f1 = couples[i].second();
            const labelList& fEdges1 = mesh_.faceEdges(f1, fe1);
            forAll(fEdges1, fEdgeI)
            {
                nBafflesPerEdge[fEdges1[fEdgeI]] += baffleValue;
            }
        }
    }

    // Add nBaffles on shared edges
    syncTools::syncEdgeList
    (
        mesh_,
        nBafflesPerEdge,
        plusEqOp<label>(),  // in-place add
        label(0)            // initial value
    );


    // Baffles which are not next to other boundaries and baffles will have
    // nBafflesPerEdge value 2*baffleValue+2*1 (from 2 boundary faces which
    // are both baffle faces)

    List<labelPair> filteredCouples(couples.size());
    label filterI = 0;

    forAll(couples, i)
    {
        const labelPair& couple = couples[i];

        if
        (
            patches.whichPatch(couple.first())
         == patches.whichPatch(couple.second())
        )
        {
            const labelList& fEdges = mesh_.faceEdges(couple.first());

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (nBafflesPerEdge[edgeI] == 2*baffleValue+2*1)
                {
                    filteredCouples[filterI++] = couple;
                    break;
                }
            }
        }
    }
    filteredCouples.setSize(filterI);


    label nFiltered = returnReduce(filteredCouples.size(), sumOp<label>());

    Info<< "freeStandingBaffles : detected "
        << nFiltered
        << " free-standing baffles out of "
        << returnReduce(couples.size(), sumOp<label>())
        << nl << endl;


    if (nFiltered > 0)
    {
        // Collect segments
        // ~~~~~~~~~~~~~~~~

        pointField start(filteredCouples.size());
        pointField end(filteredCouples.size());

        const pointField& cellCentres = mesh_.cellCentres();

        forAll(filteredCouples, i)
        {
            const labelPair& couple = filteredCouples[i];
            start[i] = cellCentres[mesh_.faceOwner()[couple.first()]];
            end[i] = cellCentres[mesh_.faceOwner()[couple.second()]];
        }

        // Extend segments a bit
        {
            const vectorField smallVec(ROOTSMALL*(end-start));
            start -= smallVec;
            end += smallVec;
        }


        // Do test for intersections
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;
        vectorField normal2;

        surfaces_.findNearestIntersection
        (
            identity(surfaces_.surfaces().size()),
            start,
            end,

            surface1,
            hit1,
            region1,
            normal1,

            surface2,
            hit2,
            region2,
            normal2
        );

        //mkDir(mesh_.time().path()/timeName());
        //OBJstream str
        //(
        //    mesh_.time().path()/timeName()/"flatBaffles.obj"
        //);

        const scalar planarAngleCos = Foam::cos(degToRad(planarAngle));

        label filterI = 0;
        forAll(filteredCouples, i)
        {
            const labelPair& couple = filteredCouples[i];

            // Note: for a baffle-surface we do not want to merge the baffle.
            // We could either check for hitting the same triangle (but you
            // might hit same point on neighbouring triangles due to tolerance
            // issues) or better just to compare the hit point.
            // This might still go wrong for a ray in the plane of the triangle
            // which might hit two different points on the same triangle due
            // to tolerances...

            if
            (
                hit1[i].hit()
             && hit2[i].hit()
             && mag(hit1[i].hitPoint()-hit2[i].hitPoint()) > mergeDistance_
            )
            {
                // Two different hits. Check angle.
                //str.write
                //(
                //    linePointRef(hit1[i].hitPoint(), hit2[i].hitPoint()),
                //    normal1[i],
                //    normal2[i]
                //);

                if ((normal1[i]&normal2[i]) > planarAngleCos)
                {
                    // Both normals aligned
                    vector n = end[i]-start[i];
                    scalar magN = mag(n);
                    if (magN > VSMALL)
                    {
                        filteredCouples[filterI++] = couple;
                    }
                }
            }
            else if (hit1[i].hit() || hit2[i].hit())
            {
                // Single hit. Do not include in freestanding baffles.
            }
        }

        filteredCouples.setSize(filterI);

        Info<< "freeStandingBaffles : detected "
            << returnReduce(filterI, sumOp<label>())
            << " planar (within " << planarAngle
            << " degrees) free-standing baffles out of "
            << nFiltered
            << nl << endl;
    }

    return filteredCouples;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeBaffles
(
    const List<labelPair>& couples,
    const Map<label>& faceToPatch,
    const bool reallocateFaceZone,
    const bool updateIntersection
)
{
    autoPtr<mapPolyMesh> map;

    if (returnReduce(couples.size()+faceToPatch.size(), sumOp<label>()))
    {
        HashTable<label, labelPair, labelPair::Hash<>> zoneIDsToFaceZone;
        labelList cellToZone(mesh_.nCells(), -1);
        if (reallocateFaceZone)
        {
            const cellZoneMesh& cellZones = mesh_.cellZones();

            forAll(cellZones, cellZoneI)
            {
                const cellZone& cZone = mesh_.cellZones()[cellZoneI];
                forAll(cZone, czi)
                {
                    label celli = cZone[czi];
                    cellToZone[celli] = cellZoneI;
                }
            }

            labelList neiCellZone;
            syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

            const faceZoneMesh& faceZones = mesh_.faceZones();
            forAll(faceZones, faceZoneI)
            {
                const faceZone& fZone = mesh_.faceZones()[faceZoneI];
                forAll(fZone, fzi)
                {
                    label facei = fZone[fzi];
                    label own = mesh_.faceOwner()[facei];
                    label ownZone = cellToZone[own];

                    label neiZone = -1;
                    if (facei < mesh_.nInternalFaces())
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        neiZone = cellToZone[nei];
                    }
                    else
                    {
                        neiZone = neiCellZone[facei - mesh_.nInternalFaces()];
                    }

                    if (ownZone != neiZone)
                    {
                        // Make sure lowest number cellZone is master.
                        // Non-cellZone areas are slave
                        bool swap = (ownZone > neiZone);

                        // Quick check whether we already have pair of zones
                        labelPair key(ownZone, neiZone);
                        if (swap)
                        {
                            Swap(key.first(), key.second());
                        }

                        HashTable<label, labelPair, labelPair::Hash<>>::
                            const_iterator zoneFnd = zoneIDsToFaceZone.find
                            (
                                key
                            );

                        if (zoneFnd == zoneIDsToFaceZone.end())
                        {
                            zoneIDsToFaceZone.insert(key, faceZoneI);
                            break;
                        }
                    }
                }
            }
            Pstream::mapCombineGather(zoneIDsToFaceZone, eqOp<label>());
            Pstream::mapCombineScatter(zoneIDsToFaceZone);
        }

        // Mesh change engine
        polyTopoChange meshMod(mesh_);

        const faceList& faces = mesh_.faces();
        const labelList& faceOwner = mesh_.faceOwner();
        const faceZoneMesh& faceZones = mesh_.faceZones();

        labelList edgeFaceZones(mesh_.nEdges(), -1);
        if (reallocateFaceZone)
        {
            forAll(mesh_.edges(), edgei)
            {
                const labelList& eFaces = mesh_.edgeFaces()[edgei];
                forAll(eFaces, eFI)
                {
                    label facei = eFaces[eFI];
                    label zoneID = faceZones.whichZone(facei);
                    edgeFaceZones[edgei] = max(edgeFaceZones[edgei],zoneID);
                }
            }
            syncTools::syncEdgeList
            (
                mesh_,
                edgeFaceZones,
                maxEqOp<label>(),  // in-place add
                label(-1)            // initial value
             );
        }

        forAll(couples, i)
        {
            label face0 = couples[i].first();
            label face1 = couples[i].second();

            // face1 < 0 signals a coupled face that has been converted to
            // baffle

            label own0 = faceOwner[face0];
            label own1 = faceOwner[face1];

            if (face1 < 0 || own0 < own1)
            {
                // Use face0 as the new internal face.
                label zoneID = -1;
                bool zoneFlip = false;

                label nei = (face1 < 0 ? -1 : own1);

                if (reallocateFaceZone && nei > -1)
                {
                    label ownZone = cellToZone[own0];
                    label neiZone = cellToZone[nei];
                    if (ownZone != neiZone)
                    {
                        zoneFlip = (ownZone > neiZone);
                        // Quick check whether we already have pair of zones
                        labelPair key(ownZone, neiZone);
                        if (zoneFlip)
                        {
                            Swap(key.first(), key.second());
                        }
                        HashTable<label, labelPair, labelPair::Hash<>>::
                            const_iterator zoneFnd = zoneIDsToFaceZone.find
                            (
                                key
                            );

                        if (zoneFnd != zoneIDsToFaceZone.end())
                        {
                            label foundEdge = zoneFnd();
                            const labelList& fEdges = mesh_.faceEdges()[face0];
                            label maxEdgeZone = -1;
                            forAll(fEdges,fEI)
                            {
                                label edgei = fEdges[fEI];
                                label edgeZone = edgeFaceZones[edgei];

                                if (foundEdge == edgeZone)
                                {
                                    zoneID = foundEdge;
                                    break;
                                }
                                maxEdgeZone = max(maxEdgeZone,edgeZone);
                            }
                            if (zoneID != foundEdge)
                            {
                                zoneID =  maxEdgeZone;
                            }
                        }
                        else
                        {
                            zoneFlip = false;
                        }
                    }
                    else
                    {
                        zoneID = faceZones.whichZone(face0);
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = faceZones[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                        }
                    }
                }
                else
                {
                    zoneID = faceZones.whichZone(face0);
                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = faceZones[zoneID];
                        zoneFlip = fZone.flipMap()[fZone.whichFace(face0)];
                    }
                }

                meshMod.setAction(polyRemoveFace(face1));
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        faces[face0],           // modified face
                        face0,                  // label of face being modified
                        own0,                   // owner
                        nei,                    // neighbour
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
                label zoneID = -1;
                bool zoneFlip = false;
                if (reallocateFaceZone)
                {
                    label ownZone = cellToZone[own1];
                    label neiZone = cellToZone[own0];
                    if (ownZone != neiZone)
                    {
                        zoneFlip = (ownZone > neiZone);
                        // Quick check whether we already have pair of zones
                        labelPair key(ownZone, neiZone);
                        if (zoneFlip)
                        {
                            Swap(key.first(), key.second());
                        }

                        HashTable<label, labelPair, labelPair::Hash<>>::
                            const_iterator zoneFnd = zoneIDsToFaceZone.find
                            (
                                key
                            );

                        if (zoneFnd != zoneIDsToFaceZone.end())
                        {
                            label foundEdge = zoneFnd();
                            const labelList& fEdges = mesh_.faceEdges()[face1];
                            label maxEdgeZone = -1;
                            forAll(fEdges,fEI)
                            {
                                label edgei = fEdges[fEI];
                                label edgeZone = edgeFaceZones[edgei];

                                if (foundEdge == edgeZone)
                                {
                                    zoneID = foundEdge;
                                    break;
                                }
                                maxEdgeZone = max(maxEdgeZone,edgeZone);
                            }
                            if (zoneID != foundEdge)
                            {
                                zoneID =  maxEdgeZone;
                            }
                        }
                        else
                        {
                            zoneFlip = false;
                        }
                    }
                    else
                    {
                        zoneID = faceZones.whichZone(face1);
                        if (zoneID >= 0)
                        {
                            const faceZone& fZone = faceZones[zoneID];
                            zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                        }
                    }
                }
                else
                {
                    zoneID = faceZones.whichZone(face1);

                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = faceZones[zoneID];
                        zoneFlip = fZone.flipMap()[fZone.whichFace(face1)];
                    }
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

        forAllConstIter(Map<label>, faceToPatch, iter)
        {
            label facei = iter.key();
            label patchI = iter();

            if (!mesh_.isInternalFace(facei))
            {
                FatalErrorInFunction
                    << "problem: face:" << facei
                    << " at:" << mesh_.faceCentres()[facei]
                    << "(wanted patch:" << patchI
                    << ") is an internal face" << exit(FatalError);
            }

            label zoneID = faceZones.whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = faceZones[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    faces[facei],           // modified face
                    facei,                  // label of face being modified
                    faceOwner[facei],       // owner
                    -1,                     // neighbour
                    false,                  // face flip
                    patchI,                 // patch for face
                    false,                  // remove from zone
                    zoneID,                 // zone for face
                    zoneFlip                // face flip in zone
                )
            );
        }


        // Change the mesh (no inflation)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing does not do this)
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Recalculate intersections on merged faces since
        // this seems to give problems? Note: should not be necessary since
        // baffles preserve intersections from when they were created.
        labelList newExposedFaces(2*couples.size());
        label newI = 0;

        forAll(couples, i)
        {
            label newFace0 = map().reverseFaceMap()[couples[i].first()];
            if (newFace0 != -1)
            {
                newExposedFaces[newI++] = newFace0;
            }

            label newFace1 = map().reverseFaceMap()[couples[i].second()];
            if (newFace1 != -1)
            {
                newExposedFaces[newI++] = newFace1;
            }
        }
        newExposedFaces.setSize(newI);
        updateMesh(map, newExposedFaces, updateIntersection);
    }

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergeZoneBaffles
(
    const bool doInternalZones,
    const bool doBaffleZones,
    const bool updateIntersection
)
{
    labelList zoneIDs;
    {
        DynamicList<surfaceZonesInfo::faceZoneType> fzTypes;
        if (doInternalZones)
        {
            fzTypes.append(surfaceZonesInfo::INTERNAL);
        }
        if (doBaffleZones)
        {
            fzTypes.append(surfaceZonesInfo::BAFFLE);
        }
        zoneIDs = getZones(fzTypes);
    }

    List<labelPair> zoneBaffles
    (
        subsetBaffles
        (
            mesh_,
            zoneIDs,
            localPointRegion::findDuplicateFacePairs(mesh_)
        )
    );

    autoPtr<mapPolyMesh> mapPtr;
    if (returnReduce(zoneBaffles.size(), sumOp<label>()))
    {
        mapPtr = mergeBaffles
        (
            zoneBaffles,
            Map<label>(0),
            false, //perform faceZone checks for new faces
            updateIntersection
        );
    }
    return mapPtr;
}


// Finds region per cell for cells inside closed named surfaces
void Foam::meshRefinement::findCellZoneGeometric
(
    const pointField& neiCc,
    const labelList& closedNamedSurfaces,   // indices of closed surfaces
    labelList& namedSurfaceIndex,           // per face index of named surface
    const labelList& surfaceToCellZone,     // cell zone index per surface

    labelList& cellToZone
) const
{
    const pointField& cellCentres = mesh_.cellCentres();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Check if cell centre is inside
    labelList insideSurfaces;
    surfaces_.findInside
    (
        closedNamedSurfaces,
        cellCentres,
        insideSurfaces
    );

    forAll(insideSurfaces, celli)
    {
        label surfI = insideSurfaces[celli];

        if (surfI != -1)
        {
            if (cellToZone[celli] == -2)
            {
                cellToZone[celli] = surfaceToCellZone[surfI];
            }
            else if (cellToZone[celli] == -1)
            {
                // ? Allow named surface to override background zone (-1)
                // This is used in the multiRegionHeater tutorial where the
                // locationInMesh is inside a named surface.
                cellToZone[celli] = surfaceToCellZone[surfI];
            }
        }
    }


    // Some cells with cell centres close to surface might have
    // had been put into wrong surface. Recheck with perturbed cell centre.


    // 1. Collect points

    // Count points to test.
    label nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            if (mesh_.isInternalFace(facei))
            {
                nCandidates += 2;
            }
            else
            {
                nCandidates += 1;
            }
        }
    }

    // Collect points.
    pointField candidatePoints(nCandidates);
    nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            label own = faceOwner[facei];
            const point& ownCc = cellCentres[own];

            if (mesh_.isInternalFace(facei))
            {
                label nei = faceNeighbour[facei];
                const point& neiCc = cellCentres[nei];
                // Perturbed cc
                const vector d = 1e-4*(neiCc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
                candidatePoints[nCandidates++] = neiCc+d;
            }
            else
            {
                //const point& neiFc = mesh_.faceCentres()[faceI];
                const point& neiFc = neiCc[facei-mesh_.nInternalFaces()];

                // Perturbed cc
                const vector d = 1e-4*(neiFc - ownCc);
                candidatePoints[nCandidates++] = ownCc-d;
            }
        }
    }


    // 2. Test points for inside

    surfaces_.findInside
    (
        closedNamedSurfaces,
        candidatePoints,
        insideSurfaces
    );


    // 3. Update zone information

    nCandidates = 0;
    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];

        if (surfI != -1)
        {
            label own = faceOwner[facei];

            if (mesh_.isInternalFace(facei))
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }

                label neiSurfI = insideSurfaces[nCandidates++];
                if (neiSurfI != -1)
                {
                    label nei = faceNeighbour[facei];

                    cellToZone[nei] = surfaceToCellZone[neiSurfI];
                }
            }
            else
            {
                label ownSurfI = insideSurfaces[nCandidates++];
                if (ownSurfI != -1)
                {
                    cellToZone[own] = surfaceToCellZone[ownSurfI];
                }
            }
        }
    }


    // Adapt the namedSurfaceIndex
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // for if any cells were not completely covered.

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label ownZone = cellToZone[mesh_.faceOwner()[facei]];
        label neiZone = cellToZone[mesh_.faceNeighbour()[facei]];

        if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
        {
            // Give face the zone of min cell zone (but only if the
            // cellZone originated from a closed, named surface)

            label minZone;
            if (ownZone == -1)
            {
                minZone = neiZone;
            }
            else if (neiZone == -1)
            {
                minZone = ownZone;
            }
            else
            {
                minZone = min(ownZone, neiZone);
            }

            // Make sure the cellZone originated from a closed surface
            label geomSurfI = findIndex(surfaceToCellZone, minZone);

            if (geomSurfI != -1)
            {
                namedSurfaceIndex[facei] = geomSurfI;
            }
        }
    }

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
                {
                    // Give face the min cell zone
                    label minZone;
                    if (ownZone == -1)
                    {
                        minZone = neiZone;
                    }
                    else if (neiZone == -1)
                    {
                        minZone = ownZone;
                    }
                    else
                    {
                        minZone = min(ownZone, neiZone);
                    }

                    // Make sure the cellZone originated from a closed surface
                    label geomSurfI = findIndex(surfaceToCellZone, minZone);

                    if (geomSurfI != -1)
                    {
                        namedSurfaceIndex[facei] = geomSurfI;
                    }
                }
            }
        }
    }

    // Sync
    syncTools::syncFaceList(mesh_, namedSurfaceIndex, maxEqOp<label>());
}


void Foam::meshRefinement::findCellZoneInsideWalk
(
    const pointField& locationsInMesh,
    const labelList& zonesInMesh,
    const labelList& faceToZone, // per face -1 or some index >= 0
    labelList& namedSurfaceIndex,

    labelList& cellToZone
) const
{
    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());
    //selectSeparatedCoupledFaces(blockedFace);

    forAll(blockedFace, facei)
    {
        if
        (
            namedSurfaceIndex.size() && namedSurfaceIndex[facei] == -1
           && faceToZone[facei] == -1
        )
        {
            blockedFace[facei] = false;
        }
        else if (namedSurfaceIndex.size() == 0 && faceToZone[facei] == -1)
        {
            blockedFace[facei] = false;
        }
        else
        {
            blockedFace[facei] = true;
        }
    }

    // No need to sync since faceToZone already is synced

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);
    blockedFace.clear();

    // Mark off which regions correspond to a zone
    // (note zone is -1 for the non-zoned bit so initialise to -2)
    labelList regionToZone(cellRegion.nRegions(), -2);


    // Force calculation of face decomposition (used in findCell)
    (void)mesh_.tetBasePtIs();

    // For all locationsInMesh find the cell
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];
        label zoneID = zonesInMesh[i];

        label keepRegionI = -1;

        label celli = findCell
        (
            insidePoint,
            mesh_,
            meshCutter_
        );

        if (celli != -1)
        {
            keepRegionI = cellRegion[celli];
        }
        reduce(keepRegionI, maxOp<label>());
/*
        // Find the region containing the insidePoint
        label keepRegionI = findRegion
        (
            mesh_,
            cellRegion,
            mergeDistance_*vector(1,1,1),
            insidePoint
        );
*/
        Info<< "For cellZone "
            <<  (
                    zoneID == -1
                  ? "none"
                  : mesh_.cellZones()[zoneID].name()
                )
            << " found point " << insidePoint
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorInFunction
                << "Point " << insidePoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }

        // Mark correspondence to zone
        regionToZone[keepRegionI] = zoneID;


        // Set all cells with this region to the zoneID
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        label nWarnings = 0;

        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == keepRegionI)
            {
                label& cZone = cellToZone[celli];

                if (cZone == -2)
                {
                    // First visit of cell
                    cZone = zoneID;
                }
                else if (cZone != zoneID)
                {
                    if (nWarnings < 10 && zoneID != -1)
                    {
                        WarningInFunction
                            << "Cell " << celli
                            << " at " << mesh_.cellCentres()[celli]
                            << " is inside cellZone "
                            <<  (
                                    zoneID == -1
                                  ? "none"
                                  : mesh_.cellZones()[zoneID].name()
                                )
                            << " from locationInMesh " << insidePoint
                            << " but already marked as being in zone : "
                            << (cZone < 0 ? "none" : mesh_.cellZones()[cZone].name())
                            << endl
                            << "This can happen if your surfaces are not "
                            << "(sufficiently) closed or locationsInMesh "
                            << "defined more than once for same region."
                            << endl;
                        nWarnings++;
                    }
                }
            }
        }
    }
}


void Foam::meshRefinement::findCellZoneInsideWalk
(
    const pointField& locationsInMesh,
    const wordList& zoneNamesInMesh,
    const labelList& faceToZone,    // per face -1 or some index >= 0
    labelList& namedSurfaceIndex,

    labelList& cellToZone
) const
{
    const cellZoneMesh& czs = mesh_.cellZones();

    labelList zoneIDs(zoneNamesInMesh.size());
    forAll(zoneNamesInMesh, i)
    {
        zoneIDs[i] = czs.findZoneID(zoneNamesInMesh[i]);
    }
    findCellZoneInsideWalk
    (
        locationsInMesh,
        zoneIDs,
        faceToZone,
        namedSurfaceIndex,

        cellToZone
    );
}


// Finds region per cell for cells inside named surfaces
void Foam::meshRefinement::findCellZoneClosed
(
    const labelList& closedNamedSurfaces,   // indices of named surfaces
    const labelList& surfaceToCellZone,     // cell zone index per surface
    labelList& namedSurfaceIndex,     // per face index of named surface
    labelList& cellToZone,
    label minZoneRegionSize
) const
{
    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    boolList freeStanding(surfZones.size(), false);
    forAll(surfZones, surfI)
    {
        freeStanding[surfI] = surfZones[surfI].freeStanding();
    }

    boolList blockedFaces(mesh_.nFaces(), false);

    forAll(mesh_.faces(), facei)
    {
        if (namedSurfaceIndex[facei] != -1)
        {
            blockedFaces[facei] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFaces, orEqOp<bool>());

    regionSplit globalRegion(mesh_, blockedFaces);

    const pointField& pts = mesh_.cellCentres();
    surfaces_.findCellZone
    (
        pts,
        globalRegion,
        surfaceToCellZone,
        cellToZone,
        minZoneRegionSize
    );

    // Adapt the namedSurfaceIndex
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // for if any cells were not completely covered.

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label ownZone = cellToZone[mesh_.faceOwner()[facei]];
        label neiZone = cellToZone[mesh_.faceNeighbour()[facei]];

        if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
        {
            label namedSurfI = findIndex
            (
                surfaceToCellZone,
                max(ownZone, neiZone)
            );

            if
            (
                namedSurfI != -1 && surfaceToCellZone[namedSurfI] != -1
                && !freeStanding[namedSurfI]
            )
            {
                // Give face the zone of max cell zone
                namedSurfaceIndex[facei] = namedSurfI;
            }
        }
    }

    labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces());
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                neiCellZone[facei-mesh_.nInternalFaces()] = ownZone;
            }
        }
    }
    syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
                {
                    label namedSurfI = findIndex
                    (
                        surfaceToCellZone,
                        max(ownZone, neiZone)
                     );

                    if
                    (
                        namedSurfI != -1 && surfaceToCellZone[namedSurfI] != -1
                        && !freeStanding[namedSurfI]
                    )
                    {
                        // Give face the zone of max cell zone
                        namedSurfaceIndex[facei] = namedSurfI;
                    }
                }
            }
        }
    }

    // Sync
    syncTools::syncFaceList
    (
        mesh_,
        namedSurfaceIndex,
        maxEqOp<label>()
    );
}


void Foam::meshRefinement::findCellZoneTopo
(
    const label backgroundZoneID,
    const List<point>& locationsInMesh,
    const labelList& unnamedSurfaceRegion,
    const labelList& namedSurfaceIndex,
    const labelList& surfaceToCellZone,
    labelList& cellToZone
) const
{
    // This routine splits the mesh into regions, based on the intersection
    // with a surface. The problem is that we know the surface which
    // (intersected) face belongs to (in namedSurfaceIndex) but we don't
    // know which side of the face it relates to. So all we are doing here
    // is get the correspondence between surface/cellZone and regionSplit
    // region.

    // Assumes:
    // - region containing keepPoint does not go into a cellZone
    // - all other regions can be found by crossing faces marked in
    //   namedSurfaceIndex.

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces());

    forAll(unnamedSurfaceRegion, facei)
    {
        if (unnamedSurfaceRegion[facei] == -1 && namedSurfaceIndex[facei] == -1)
        {
            blockedFace[facei] = false;
        }
        else
        {
            blockedFace[facei] = true;
        }
    }
    // No need to sync since namedSurfaceIndex already is synced

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);

    // Per mesh region the zone the cell should be put in.
    // -2   : not analysed yet
    // -1   : keepPoint region. Not put into any cellzone.
    // >= 0 : index of cellZone
    labelList regionToCellZone(cellRegion.nRegions(), -2);

    boolList alreadySet(surfaceToCellZone.size(), false);

    // See which cells already are set in the cellToZone (from geometric
    // searching) and use these to take over their zones.
    // Note: could be improved to count number of cells per region.
    forAll(cellToZone, celli)
    {
        label zoneI = cellToZone[celli];
        if (zoneI != -2)
        {
            label regionI = cellRegion[celli];
            if (regionToCellZone[regionI] == -2)
            {
                regionToCellZone[regionI] = zoneI;
            }
            if (zoneI >= 0)
            {
                alreadySet[zoneI] = true;
            }
        }
    }
    Pstream::listCombineReduce(alreadySet, orOp<bool>());

    // Synchronise regionToCellZone.
    // Note:
    // - region numbers are identical on all processors
    // - keepRegion is identical ,,
    // - cellZones are identical ,,
    Pstream::listCombineReduce(regionToCellZone, maxOp<label>());


    labelList keepCell(locationsInMesh.size(), -1);

    // Find the region containing the keepPoint
    forAll(locationsInMesh, i)
    {
        const point& keepPoint = locationsInMesh[i];

        // Find the region containing the keepPoint
        label keepRegionI = -1;

        label celli = findCell
        (
            keepPoint,
            mesh_,
            meshCutter_
        );

        if (celli != -1)
        {
            keepCell[i] = celli;
            keepRegionI = cellRegion[celli];
        }
        reduce(keepRegionI, maxOp<label>());

        Info<< "Found point " << keepPoint
            << " in global region " << keepRegionI
            << " out of " << cellRegion.nRegions() << " regions." << endl;

        if (keepRegionI == -1)
        {
            FatalErrorInFunction
                << "Point " << keepPoint
                << " is not inside the mesh." << nl
                << "Bounding box of the mesh:" << mesh_.bounds()
                << exit(FatalError);
        }

        // Mark default region with zone -1. Note that all regions should
        // already be matched to a cellZone through the loop over cellToZone.
        // This is just to mop up any remaining regions. Not sure whether
        // this is needed?
        if (regionToCellZone[keepRegionI] == -2)
        {
            regionToCellZone[keepRegionI] = -1;
        }
    }

    label nChangedEnter = 0;
    label nChangedLeave = 0;
    bool enter = true;

    boolList walkReset(surfaceToCellZone.size(), false);

    // Find correspondence between cell zone and surface
    // by changing cell zone every time we cross a surface.
    while (true)
    {
        // Synchronise regionToCellZone.
        // Note:
        // - region numbers are identical on all processors
        // - keepRegion is identical ,,
        // - cellZones are identical ,,
        // This done at top of loop to account for geometric matching
        // not being synchronised.
        Pstream::listCombineReduce(regionToCellZone, maxOp<label>());

        bool changed = false;

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        // Get coupled neighbour cellRegion
        labelList neiCellRegion;
        syncTools::swapBoundaryCellList(mesh_, cellRegion, neiCellRegion);

        // Internal faces
        forAll(mesh_.faces(), facei)
        {
            label surfI = namedSurfaceIndex[facei];
            if (unnamedSurfaceRegion[facei] != -1 || surfI == -1)
            {
                continue;
            }

            label patchi = patches.whichPatch(facei);

            bool coupled = false;
            if (patchi != -1)
            {
                if (patches[patchi].coupled())
                {
                    coupled = true;
                }
                else
                {
                    continue;
                }
            }

            label surfZoneI = surfaceToCellZone[surfI];
            if (!enter && surfZoneI > -1 && alreadySet[surfZoneI])
            {
                label ownRegion = cellRegion[mesh_.faceOwner()[facei]];
                label neiRegion;
                if (coupled)
                {
                    neiRegion = neiCellRegion[facei-mesh_.nInternalFaces()];
                }
                else
                {
                    neiRegion = cellRegion[mesh_.faceNeighbour()[facei]];
                }

                label ownRegionCellZone = regionToCellZone[ownRegion];
                label neiRegionCellZone = regionToCellZone[neiRegion];
                if (ownRegion != neiRegion)
                {
                    if
                    (
                        ownRegionCellZone == -2 && neiRegionCellZone > -1
                        && alreadySet[neiRegionCellZone]
                    )
                    {
                        regionToCellZone[ownRegion] = -1;
                        nChangedLeave++;
                        changed = true;
                    }
                    else if
                    (
                        neiRegionCellZone == -2 && ownRegionCellZone > -1
                        && alreadySet[ownRegionCellZone]
                    )
                    {
                        regionToCellZone[neiRegion] = -1;
                        nChangedLeave++;
                        changed = true;
                    }
                }
            }
            else if
            (
                surfZoneI == -1 || !alreadySet[surfZoneI]
            )
            {
                label ownRegion = cellRegion[mesh_.faceOwner()[facei]];
                label neiRegion;
                if (coupled)
                {
                    neiRegion = neiCellRegion[facei-mesh_.nInternalFaces()];
                }
                else
                {
                    neiRegion = cellRegion[mesh_.faceNeighbour()[facei]];
                }

                label ownRegionCellZone = regionToCellZone[ownRegion];
                label neiRegionCellZone = regionToCellZone[neiRegion];

                if (ownRegion != neiRegion)
                {
                    if (enter && surfZoneI != -1)
                    {
                        if (ownRegionCellZone == -2)
                        {
                            if
                            (
                                neiRegionCellZone == -1
                                ||
                                (
                                    neiRegionCellZone > -1
                                    &&
                                    (
                                        alreadySet[neiRegionCellZone]
                                        ||
                                        (
                                            walkReset[neiRegionCellZone]
                                            && surfZoneI != neiRegionCellZone
                                        )
                                    )
                                )
                            )
                            {
                                regionToCellZone[ownRegion] = surfZoneI;
                                nChangedEnter++;
                                changed = true;
                            }
                        }
                        else if (neiRegionCellZone == -2)
                        {
                            if
                            (
                                ownRegionCellZone == -1
                                ||
                                (
                                    ownRegionCellZone > -1
                                    &&
                                    (
                                        alreadySet[ownRegionCellZone]
                                        ||
                                        (
                                            walkReset[ownRegionCellZone]
                                            && surfZoneI != ownRegionCellZone
                                        )
                                    )
                                )
                            )
                            {
                                regionToCellZone[neiRegion] = surfZoneI;
                                nChangedEnter++;
                                changed = true;
                            }
                        }
                    }
                    else if (!enter && surfZoneI == -1)
                    {
                        if
                        (
                            ownRegionCellZone >= 0
                            && !alreadySet[ownRegionCellZone]
                            && neiRegionCellZone == -2
                        )
                        {
                            regionToCellZone[neiRegion] = -1;
                            nChangedLeave++;
                            changed = true;
                        }
                        else if
                        (
                            neiRegionCellZone >= 0
                            && !alreadySet[neiRegionCellZone]
                            && ownRegionCellZone == -2
                        )
                        {
                            regionToCellZone[ownRegion] = -1;
                            nChangedLeave++;
                            changed = true;
                        }
                    }
                    else if
                    (
                        !enter && surfZoneI > -1 && walkReset[surfZoneI]
                    )
                    {
                        if
                        (
                            ownRegionCellZone == surfZoneI
                            && neiRegionCellZone == -2
                        )
                        {
                            regionToCellZone[neiRegion] = -1;
                            nChangedLeave++;
                            changed = true;
                        }
                        else if
                        (
                            neiRegionCellZone == surfZoneI
                            && ownRegionCellZone == -2
                        )
                        {
                            regionToCellZone[ownRegion] = -1;
                            nChangedLeave++;
                            changed = true;
                        }
                    }
                }
            }
        }

        if (!returnReduce(changed, orOp<bool>()))
        {
            if (!enter)
            {
                label nReset = 0;
                forAll(regionToCellZone, regioni)
                {
                    label zi = regionToCellZone[regioni];
                    if (zi > -1 && !alreadySet[zi] && !walkReset[zi])
                    {
                        walkReset[zi] = true;
                        nReset++;
                    }
                }
                reduce(nReset,sumOp<label>());
                if (nReset > 0)
                {
                    Pstream::listCombineReduce(walkReset, orOp<bool>());
                }

                label totalChanged = nChangedLeave + nChangedEnter;
                reduce(totalChanged,sumOp<label>());
                if (totalChanged == 0 && nReset == 0)
                {
                    break;
                }
                else
                {
                    nChangedLeave = 0;
                    nChangedEnter = 0;
                }
            }
            enter = !enter;
        }
    }

    if (debug)
    {
        forAll(regionToCellZone, regionI)
        {
            Pout<< "Region " << regionI
                << " becomes cellZone:" << regionToCellZone[regionI]
                << endl;
        }
    }

    // Rework into cellToZone
    forAll(cellToZone, celli)
    {
        label regionI = cellRegion[celli];
        if (cellToZone[celli] == -2 && regionToCellZone[regionI] != -2)
        {
            cellToZone[celli] = regionToCellZone[regionI];
        }
    }

    // Reset unvisited cells (-2) connected to keepPoints
    forAll(mesh_.faces(), facei)
    {
        if (unnamedSurfaceRegion[facei] == -1)
        {
            blockedFace[facei] = false;
        }
        else
        {
            blockedFace[facei] = true;
        }
    }
    // No need to sync since namedSurfaceIndex already is synced

    // Set region per cell based on walking
    regionSplit cellUnnamedRegion(mesh_, blockedFace);

    boolList keepRegions(cellUnnamedRegion.nRegions(), false);

    // Find the region containing the keepPoint
    forAll(locationsInMesh, i)
    {
        label celli = keepCell[i];

        label keepRegionI = -1;
        if (celli != -1)
        {
            keepRegionI = cellUnnamedRegion[celli];
        }
        reduce(keepRegionI, maxOp<label>());

        if (keepRegionI >=0)
        {
            keepRegions[keepRegionI] = true;
        }
    }

    forAll(mesh_.cells(), celli)
    {
        if (cellToZone[celli] == -2 && keepRegions[cellUnnamedRegion[celli]])
        {
            cellToZone[celli] = -1;
        }
    }
}


// Finds size of split regions and removes from zones if too small
void Foam::meshRefinement::excludeSmallRegions
(
    const labelHashSet& namedLocationsSet,
    labelList& globalRegion1,
    labelList& globalRegion2,
    labelList& namedSurfaceIndex,
    labelList& cellToZone,
    const label minSize,
    const Switch namedLocationsRezone
) const
{

    Info<< nl
        << "Excluding small regions from zone" << nl
        << "---------------------------------" << nl
        << endl;

    const labelList origCellToZone = cellToZone;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    bool foundNamed
    (
        returnReduce(namedSurfaceIndex.size(), sumOp<label>()) > 0
        ? true : false
    );

    while (true)
    {
        label nChanged = 0;
        boolList blockedFaces(mesh_.nFaces(), false);

        // Get coupled neighbour cellZone. Set to -2 on non-coupled patches.
        labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces(), -2);
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    neiCellZone[facei-mesh_.nInternalFaces()] =
                        cellToZone[mesh_.faceOwner()[facei]];
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

        forAll(mesh_.faces(), facei)
        {
            label own = mesh_.faceOwner()[facei];
            label ownZone = cellToZone[own];

            label neiZone = -1;
            if (facei < mesh_.nInternalFaces())
            {
                label nei = mesh_.faceNeighbour()[facei];
                neiZone = cellToZone[nei];
            }
            else
            {
                neiZone = neiCellZone[facei - mesh_.nInternalFaces()];
            }

            if (ownZone != neiZone)
            {
                blockedFaces[facei] = true;
            }
        }
        syncTools::syncFaceList(mesh_, blockedFaces, orEqOp<bool>());

        regionSplit globalRegion(mesh_, blockedFaces);

        boolList interZoneRegions(globalRegion.nRegions(), false);
        if (namedLocationsRezone)
        {
            forAll(globalRegion, celli)
            {
                if (cellToZone[celli] == -2)
                {
                    label regioni = globalRegion[celli];
                    interZoneRegions[regioni] = true;
                }
            }
            Pstream::listCombineReduce(interZoneRegions, orOp<bool>());

            labelList minNbrRegionZone(globalRegion.nRegions(),labelMax);
            labelList maxNbrRegionZone(globalRegion.nRegions(),labelMin);
            forAll(globalRegion, celli)
            {
                if (cellToZone[celli] == -2)
                {
                    label regioni = globalRegion[celli];
                    const cell& c = mesh_.cells()[celli];
                    forAll(c, cfi)
                    {
                        label facei = c[cfi];
                        label own = mesh_.faceOwner()[facei];
                        label neiZone = -1;
                        label patchi = patches.whichPatch(facei);
                        if (patchi == -1)
                        {
                            label nei = mesh_.faceNeighbour()[facei];
                            neiZone = (own == celli)
                                ? cellToZone[nei] : cellToZone[own];
                        }
                        else if (patches[patchi].coupled())
                        {
                            neiZone =
                                neiCellZone[facei - mesh_.nInternalFaces()];
                        }

                        if
                        (
                            neiZone == -1
                            ||
                            (
                                neiZone > -1
                                && !namedLocationsSet.found(neiZone)
                            )
                        )
                        {
                            interZoneRegions[regioni] = false;
                        }

                        if (neiZone > -1)
                        {
                            minNbrRegionZone[regioni] = min
                            (
                                minNbrRegionZone[regioni],
                                neiZone
                            );
                            maxNbrRegionZone[regioni] = max
                            (
                                maxNbrRegionZone[regioni],
                                neiZone
                            );
                        }
                    }
                }
            }

            Pstream::listCombineReduce(interZoneRegions, andOp<bool>());
            Pstream::listCombineReduce(maxNbrRegionZone, maxOp<label>());
            Pstream::listCombineReduce(minNbrRegionZone, minOp<label>());

            forAll(interZoneRegions, regioni)
            {
                label minRegionZone = minNbrRegionZone[regioni];
                label maxRegionZone = maxNbrRegionZone[regioni];
                if
                (
                    interZoneRegions[regioni]
                    && minRegionZone > -1 && maxRegionZone > -1
                    && minRegionZone == maxRegionZone
                )
                {
                    interZoneRegions[regioni] = false;
                }
            }
            Pstream::listCombineReduce(interZoneRegions, andOp<bool>());
        }

        labelList regionSize(globalRegion.nRegions(),0);
        boolList validRegion(globalRegion.nRegions(),true);
        forAll(globalRegion, celli)
        {
            label regioni = globalRegion[celli];
            if (cellToZone[celli] == -2)
            {
                if (namedLocationsRezone)
                {
                    if (!interZoneRegions[regioni])
                    {
                        validRegion[regioni] = false;
                    }
                }
                else
                {
                    validRegion[regioni] = false;
                }
            }
            regionSize[regioni]++;
        }
        Pstream::listCombineReduce(regionSize, plusOp<label>());
        Pstream::listCombineReduce(validRegion, andOp<bool>());

        DynamicList<label> smallRegionID(globalRegion.nRegions());
        forAll(regionSize, regioni)
        {
            if (validRegion[regioni] && regionSize[regioni] <= minSize)
            {
                smallRegionID.append(regioni);
            }
        }

        labelHashSet excludeRegions(smallRegionID);
        List<label> newZoneID(globalRegion.nRegions(), -1);
        forAll(cellToZone, celli)
        {
            label globalRegionI = globalRegion[celli];
            if
            (
                excludeRegions.found(globalRegionI)
                && newZoneID[globalRegionI] == -1
            )
            {
                const cell& cellFaces = mesh_.cells()[celli];
                forAll(cellFaces, j)
                {
                    label ownZone = cellToZone[celli];
                    label meshFaceI = cellFaces[j];
                    label own = mesh_.faceOwner()[meshFaceI];
                    label neiZone = -1;

                    if (meshFaceI < mesh_.nInternalFaces())
                    {
                        label nei = mesh_.faceNeighbour()[meshFaceI];
                        neiZone =
                            (own == celli) ? cellToZone[nei] : cellToZone[own];
                    }
                    else
                    {
                        neiZone = neiCellZone[meshFaceI
                                              - mesh_.nInternalFaces()];
                    }
                    if (neiZone >= 0 && ownZone != neiZone)
                    {
                        newZoneID[globalRegionI] = neiZone;
                    }
                }
            }
        }
        Pstream::listCombineReduce(newZoneID, maxOp<label>());
        boolList markedCells(mesh_.nCells(), false);
        forAll(cellToZone, celli)
        {
            label globalRegionI = globalRegion[celli];
            bool foundExcluded = excludeRegions.found(globalRegionI);
            if (foundExcluded)
            {
                if (newZoneID[globalRegionI] != -1)
                {
                    cellToZone[celli] = newZoneID[globalRegionI];
                    markedCells[celli] = true;
                    nChanged++;
                }
                else if (cellToZone[celli] != -1)
                {
                    cellToZone[celli] = -1;
                    markedCells[celli] = true;
                    nChanged++;
                }
            }
        }

        if (returnReduce(nChanged, sumOp<label>()) == 0)
        {
            break;
        }

        if (namedLocationsRezone)
        {
            syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
            forAll(mesh_.cells(), celli)
            {
                label regioni = globalRegion[celli];
                if (markedCells[celli] && interZoneRegions[regioni])
                {
                    const labelList& c = mesh_.cells()[celli];
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];

                        label ownZone = cellToZone[celli];
                        label neiZone = -2;

                        label patchi = patches.whichPatch(facei);
                        if (patchi == -1)
                        {
                            label own = mesh_.faceOwner()[facei];
                            label nei
                            (
                                own == celli
                                ? mesh_.faceNeighbour()[facei] : own
                             );
                            neiZone = cellToZone[nei];
                        }
                        else if
                        (
                            patchi != -1
                            && isA<processorPolyPatch>(patches[patchi])
                        )
                        {
                            neiZone =
                                neiCellZone[facei-mesh_.nInternalFaces()];
                        }

                        if (ownZone == neiZone)
                        {
                            globalRegion1[facei] = -1;
                            globalRegion2[facei] = -1;
                        }
                    }
                }
            }

            syncTools::syncFaceList
            (
                mesh_,
                globalRegion1,
                minEqOp<label>()
            );

            syncTools::syncFaceList
            (
                mesh_,
                globalRegion2,
                minEqOp<label>()
             );
        }

        if (foundNamed)
        {
            syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
            forAll(mesh_.cells(), celli)
            {
                if (markedCells[celli])
                {
                    const labelList& cFaces = mesh_.cells()[celli];
                    label namedRegion1 = -1;
                    forAll(cFaces, cFI)
                    {
                        label facei = cFaces[cFI];
                        if (namedSurfaceIndex[facei] != -1)
                        {
                            namedRegion1 = namedSurfaceIndex[facei];
                            break;
                        }
                    }

                    if (namedRegion1 != -1)
                    {
                        forAll(cFaces, cFI)
                        {
                            label facei = cFaces[cFI];

                            if (globalRegion1[facei] != -1)
                            {
                                continue;
                            }

                            label ownZone = cellToZone[celli];
                            label neiZone = -2;

                            label patchI = patches.whichPatch(facei);
                            if (patchI == -1)
                            {
                                label own = mesh_.faceOwner()[facei];
                                label nei
                                (
                                    own == celli
                                    ? mesh_.faceNeighbour()[facei] : own
                                );
                                neiZone = cellToZone[nei];
                            }
                            else if
                            (
                                patchI != -1
                                && isA<processorPolyPatch>(patches[patchI])
                            )
                            {
                                neiZone =
                                    neiCellZone[facei-mesh_.nInternalFaces()];
                            }

                            if (ownZone == neiZone)
                            {
                                namedSurfaceIndex[facei] = -1;
                            }
                            else if (namedSurfaceIndex[facei] == -1)
                            {
                                namedSurfaceIndex[facei] = namedRegion1;
                            }
                        }
                    }
                }
            }

            syncTools::syncFaceList
            (
                mesh_,
                namedSurfaceIndex,
                maxEqOp<label>()
            );

            makeConsistentFaceIndex
            (
                cellToZone,
                namedSurfaceIndex
            );
        }
    }

    if (namedLocationsRezone)
    {
        const labelList& pointLevel = meshCutter_.pointLevel();
        const labelList& cellLevel = meshCutter_.cellLevel();
        boolList blockedFaces(mesh_.nFaces(), false);
        boolList blockedPoints(mesh_.nPoints(), false);
        labelList neiCellZone(mesh_.nFaces()-mesh_.nInternalFaces(), -2);
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label facei = pp.start()+i;
                    neiCellZone[facei-mesh_.nInternalFaces()] =
                        cellToZone[mesh_.faceOwner()[facei]];
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiCellZone);

        forAll(mesh_.faces(), facei)
        {
            label own = mesh_.faceOwner()[facei];
            label ownZone = cellToZone[own];
            label neiZone = -1;
            label patchi = patches.whichPatch(facei);
            if (patchi == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                neiZone = cellToZone[nei];
            }
            else if (patches[patchi].coupled())
            {
                neiZone =
                    neiCellZone[facei - mesh_.nInternalFaces()];
            }

            if (ownZone > -1 && neiZone > -1 && ownZone != neiZone)
            {
                if
                (
                    namedLocationsSet.found(neiZone)
                    && namedLocationsSet.found(ownZone)
                )
                {
                    blockedFaces[facei] = true;
                    const face& f = mesh_.faces()[facei];
                    forAll(f,fp)
                    {
                        blockedPoints[f[fp]] = true;
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, blockedFaces, orEqOp<bool>());
        syncTools::syncPointList
        (
            mesh_,
            blockedPoints,
            orEqOp<bool>(),
            false           // null value
        );

        forAll(blockedPoints, pointi)
        {
            if (blockedPoints[pointi])
            {
                const labelList& pFaces = mesh_.pointFaces()[pointi];
                forAll(pFaces, pfi)
                {
                    label facei = pFaces[pfi];
                    label own = mesh_.faceOwner()[facei];
                    label ownZone = cellToZone[own];
                    label neiZone = -1;
                    label patchi = patches.whichPatch(facei);
                    if (patchi == -1)
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        neiZone = cellToZone[nei];
                    }
                    else if (patches[patchi].coupled())
                    {
                        neiZone =
                            neiCellZone[facei - mesh_.nInternalFaces()];
                    }

                    if (ownZone == neiZone && namedLocationsSet.found(ownZone))
                    {
                        blockedFaces[facei] = true;
                        const face& f = mesh_.faces()[facei];
                        forAll(f,fp)
                        {
                            blockedPoints[f[fp]] = true;
                        }
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, blockedFaces, orEqOp<bool>());
        syncTools::syncPointList
        (
            mesh_,
            blockedPoints,
            orEqOp<bool>(),
            false           // null value
        );

        // Does cell have exactly 7 of its 8 anchor points on the boundary?
        PackedBoolList hasSevenBoundaryAnchorPoints(mesh_.nCells());
        // If so what is the remaining non-boundary anchor point?
        labelHashSet nonBoundaryAnchors(mesh_.nCells()/10000);
        DynamicList<label> dynCPoints;

        forAll(cellLevel, celli)
        {
            const labelList& cPoints = mesh_.cellPoints(celli, dynCPoints);

            // Get number of anchor points (pointLevel <= cellLevel)

            label nBoundaryAnchors = 0;
            label nonBoundaryAnchor = -1;

            forAll(cPoints, i)
            {
                label pointi = cPoints[i];

                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    // Anchor point
                    if (blockedPoints[pointi])
                    {
                        nBoundaryAnchors++;
                    }
                    else
                    {
                        // Anchor point which is not on the surface
                        nonBoundaryAnchor = pointi;
                    }
                }
            }

            if (nBoundaryAnchors == 8)
            {
                const cell& cFaces = mesh_.cells()[celli];
                forAll(cFaces, cfi)
                {
                    label facei = cFaces[cfi];
                    label own = mesh_.faceOwner()[facei];
                    label ownZone = cellToZone[own];
                    label neiZone = -1;
                    label patchi = patches.whichPatch(facei);
                    if (patchi == -1)
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        neiZone = cellToZone[nei];
                    }
                    else if (patches[patchi].coupled())
                    {
                        neiZone =
                            neiCellZone[facei - mesh_.nInternalFaces()];
                    }

                    if (ownZone == neiZone && namedLocationsSet.found(ownZone))
                    {
                        globalRegion1[facei] = -1;
                        globalRegion2[facei] = -1;
                    }
                }
            }
            else if (nBoundaryAnchors == 7)
            {
                // Mark the cell. Store the (single!) non-boundary anchor point.
                hasSevenBoundaryAnchorPoints.set(celli, 1u);
                nonBoundaryAnchors.insert(nonBoundaryAnchor);
            }
        }

        // Loop over all points. If a point is connected to 4 or more cells
        // with 7 anchor points on the boundary set those cell's non-boundary
        // faces to baffles
        labelList nCellsSevenAnchors(mesh_.nPoints(), 0);
        DynamicList<label> dynPCells;

        forAllConstIter(labelHashSet, nonBoundaryAnchors, iter)
        {
            label pointi = iter.key();

            //in VDB pointLevel might not be accurate at procBoundary interface
            if (pointi == -1) continue;

            const labelList& pCells = mesh_.pointCells(pointi, dynPCells);

            // Count number of 'hasSevenBoundaryAnchorPoints' cells.

            forAll(pCells, i)
            {
                if (hasSevenBoundaryAnchorPoints.get(pCells[i]) == 1u)
                {
                    nCellsSevenAnchors[pointi]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            nCellsSevenAnchors,
            plusEqOp<label>(),
            label(0)              // null value
        );

        forAllConstIter(labelHashSet, nonBoundaryAnchors, iter)
        {
            label pointI = iter.key();
            const labelList& pCells = mesh_.pointCells()[pointI];
            if (nCellsSevenAnchors[pointI] > 3)
            {
               // Point in danger of being what? Remove all 7-cells.
                forAll(pCells, i)
                {
                    label celli = pCells[i];

                    if (hasSevenBoundaryAnchorPoints.get(celli) == 1u)
                    {
                        const cell& cFaces = mesh_.cells()[celli];
                        forAll(cFaces, cfi)
                        {
                            label facei = cFaces[cfi];
                            label own = mesh_.faceOwner()[facei];
                            label ownZone = cellToZone[own];
                            label neiZone = -1;
                            label patchi = patches.whichPatch(facei);
                            if (patchi == -1)
                            {
                                label nei = mesh_.faceNeighbour()[facei];
                                neiZone = cellToZone[nei];
                            }
                            else if (patches[patchi].coupled())
                            {
                                neiZone =
                                    neiCellZone[facei - mesh_.nInternalFaces()];
                            }

                            if
                            (
                                ownZone == neiZone
                                && namedLocationsSet.found(ownZone)
                            )
                            {
                                globalRegion1[facei] = -1;
                                globalRegion2[facei] = -1;
                            }
                        }
                    }
                }
            }
        }
        syncTools::syncFaceList
        (
            mesh_,
            globalRegion1,
            minEqOp<label>()
        );

        syncTools::syncFaceList
        (
            mesh_,
            globalRegion2,
            minEqOp<label>()
        );

    }
}


Foam::labelList Foam::meshRefinement::selectPinchedEdges
(
    const labelList& cellToZone,
    const labelList& keepCell
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    labelList neiKeepCell;
    syncTools::swapBoundaryCellList
    (
        mesh_,
        keepCell,
        neiKeepCell
    );

    labelList boundFaces(mesh_.nFaces(), -2);

    forAll(mesh_.faces(), facei)
    {
        label patchi = patches.whichPatch(facei);
        label own = mesh_.faceOwner()[facei];
        label ownZone = cellToZone[own];
        label ownKept = keepCell[own];

        label neiZone = -2;
        label neiKept = -2;
        bool boundaryFace = false;
        if (patchi == -1)
        {
            label nei = mesh_.faceNeighbour()[facei];
            neiZone = cellToZone[nei];
            neiKept = keepCell[nei];
        }
        else if (patches[patchi].coupled())
        {
            neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
            neiKept = neiKeepCell[facei-mesh_.nInternalFaces()];
        }
        else
        {
            boundaryFace = true;
        }

        label keepOwn = max(ownZone, ownKept);
        label keepNei = max(neiZone, neiKept);

        if (boundaryFace && (ownZone != -2 || ownKept != -2))
        {
            boundFaces[facei] = max(ownZone,ownKept);
        }
        else if
        (
            (keepOwn == -2 && keepNei != -2)
            || (keepOwn != -2 && keepNei == -2)
        )
        {
            boundFaces[facei] = max(keepOwn,keepNei);
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        boundFaces,
        maxEqOp<label>()
    );

    labelList nEdgeBoundFaces(mesh_.nEdges(), 0);
    labelList nKeptCells(mesh_.nEdges(), 0);

    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));
    forAll(mesh_.edges(), edgei)
    {
        const labelList& eCells = mesh_.edgeCells()[edgei];
        const labelList& eFaces = mesh_.edgeFaces()[edgei];
        forAll(eCells, eCI)
        {
            label celli = eCells[eCI];
            if (cellToZone[celli] != -2 || keepCell[celli] != -2)
            {
                nKeptCells[edgei]++;
            }
        }
        forAll(eFaces, eFI)
        {
            label facei = eFaces[eFI];
            if (isMasterFace[facei] && boundFaces[facei] != -2)
            {
                nEdgeBoundFaces[edgei]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        nEdgeBoundFaces,
        plusEqOp<label>(),
        label(0)                  // initial value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        nKeptCells,
        plusEqOp<label>(),
        label(0)                  // initial value
    );

    DynamicList<label> pinchedEdges(mesh_.nEdges()/10);
    forAll(mesh_.edges(), edgei)
    {
        if (nKeptCells[edgei] == 2 && nEdgeBoundFaces[edgei] == 4)
        {
            pinchedEdges.append(edgei);
        }
    }

    return labelList(pinchedEdges, true);
}


Foam::label Foam::meshRefinement::resetRegionFaces
(
    const labelList& keepCell,
    labelList& cellToZone,
    labelList& globalRegion1,
    labelList& globalRegion2,
    labelList& testFaces,
    label defaultRegion
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    labelList neiKeepCell;
    syncTools::swapBoundaryCellList
    (
        mesh_,
        keepCell,
        neiKeepCell
    );

    DynamicList<label> newBaffleFaces(mesh_.nFaces()/10);
    forAll(mesh_.faces(), facei)
    {
        label patchi = patches.whichPatch(facei);
        label own = mesh_.faceOwner()[facei];
        label ownZone = cellToZone[own];

        if (patchi == -1)
        {
            label nei = mesh_.faceNeighbour()[facei];
            if
            (
                (keepCell[own] != -2 && cellToZone[nei] != -2)
                || (keepCell[nei] != -2 && cellToZone[own] != -2)
            )
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
            }
            else if (keepCell[own] != -2 && keepCell[nei] != -2)
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
            }
            else if
            (
                globalRegion1[facei] == -1
                && (keepCell[own] != -2 || keepCell[nei] != -2)
            )
            {
                globalRegion1[facei] = defaultRegion;
                globalRegion2[facei] = defaultRegion;
                newBaffleFaces.append(facei);
            }
        }
        else if (patches[patchi].coupled())
        {
            label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
            label keepNei = neiKeepCell[facei-mesh_.nInternalFaces()];
            label keepOwn = keepCell[own];
            if
            (
                (keepOwn != -2 && neiZone != -2)
                || (keepNei != -2 && ownZone != -2)
            )
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
            }
            else if (keepOwn != -2 && keepNei != -2)
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
            }
            else if
            (
                globalRegion1[facei] == -1
                && (keepOwn != -2 || keepNei != -2)
            )
            {
                globalRegion1[facei] = defaultRegion;
                globalRegion2[facei] = defaultRegion;
                newBaffleFaces.append(facei);
            }
        }
    }

    label nUpdated = 0;
    forAll(keepCell, celli)
    {
        if (keepCell[celli] != -2)
        {
            nUpdated++;
            cellToZone[celli] = keepCell[celli];
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        globalRegion1,
        maxEqOp<label>()
     );

    syncTools::syncFaceList
    (
        mesh_,
        globalRegion2,
        maxEqOp<label>()
    );

    //Ensure testFaces includes updates
    if (newBaffleFaces.size() > 0)
    {
        label sz = testFaces.size();
        testFaces.setSize(sz+newBaffleFaces.size());
        forAll(newBaffleFaces, i)
        {
            testFaces[sz+i] = newBaffleFaces[i];
        }
    }

    return nUpdated;
}


void Foam::meshRefinement::addOutsideContactCells
(
    const refinementParameters& refineParams,
    const pointField& cellCentres,
    labelList& globalRegion1,
    labelList& globalRegion2,
    labelList& cellToZone,
    labelList& testFaces
) const
{
    const labelList surfacesToTest(identity(surfaces_.surfaces().size()));

    if (surfacesToTest.size() == 0  || !refineParams.addContactCells())
    {
        return;
    }
    else
    {
        Info<<"Checking for outside contact cells to keep"<<endl;
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& pointLevel = meshCutter_.pointLevel();
    const labelList& cellLevel = meshCutter_.cellLevel();

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    labelList neiCellLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCellCentres(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(cellCentres, neiCellLevel, neiCellCentres);

    labelList nKeptEdgeCells(mesh_.nEdges(), 0);
    forAll(mesh_.edges(), edgei)
    {
        const labelList& eCells = mesh_.edgeCells()[edgei];
        forAll(eCells, eCI)
        {
            if (cellToZone[eCells[eCI]] != -2)
            {
                nKeptEdgeCells[edgei]++;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh_,
        nKeptEdgeCells,
        plusEqOp<label>(),
        label(0) // initial value
    );

    boolList gapPts(mesh_.nPoints(), false);
    forAll(nKeptEdgeCells, edgei)
    {
        if (nKeptEdgeCells[edgei] == 1)
        {
            const edge& e = mesh_.edges()[edgei];
            gapPts[e[0]] = true;
            gapPts[e[1]] = true;
        }
    }
    syncTools::syncPointList
    (
        mesh_,
        gapPts,
        orEqOp<bool>(),
        false                 // initial value
    );

    DynamicList<label> checkCells(mesh_.nFaces()/10);
    DynamicList<label> checkFaces(mesh_.nFaces()/10);
    DynamicField<point> start(mesh_.nFaces()/10);
    DynamicField<point> end(mesh_.nFaces()/10);

    labelList pointZone(mesh_.nPoints(), -2);
    forAll(mesh_.faces(), facei)
    {
        label region1 = globalRegion1[facei];
        label region2 = globalRegion2[facei];
        const face& f = mesh_.faces()[facei];

        label nGapPts = 0;
        forAll(f,fp)
        {
            if (gapPts[f[fp]])
            {
                nGapPts++;
            }
        }

        if
        (
            f.size() != 4 || nGapPts != f.size()
            || region1 == -1 || region2 == -1
        )
        {
            continue;
        }

        label patchi = patches.whichPatch(facei);
        label own = mesh_.faceOwner()[facei];
        label ownLevel = cellLevel[own];
        label ownZone = cellToZone[own];

        label neiZone = labelMin;
        label neiLevel = labelMin;

        label splitCellI = -1;
        point nbrCC = vector::zero;
        if (patchi == -1)
        {
            label nei = mesh_.faceNeighbour()[facei];
            neiZone = cellToZone[nei];
            neiLevel = cellLevel[nei];
            if (ownZone == -2 && neiZone != -2)
            {
                nbrCC = cellCentres[nei];
                splitCellI = own;
            }
            else if (ownZone != -2 && neiZone == -2)
            {
                nbrCC = cellCentres[own];
                splitCellI = nei;
            }
        }
        else if (patches[patchi].coupled())
        {
            neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
            neiLevel = neiCellLevel[facei-mesh_.nInternalFaces()];
            if (ownZone == -2 && neiZone != -2)
            {
                nbrCC = neiCellCentres[facei-mesh_.nInternalFaces()];
                splitCellI = own;
            }
        }
        else
        {
            continue;
        }

        if (splitCellI != -1 && ownLevel == neiLevel)
        {
            label maxZone = max(ownZone,neiZone);
            forAll(f,fp)
            {
                label pointi = f[fp];
                pointZone[pointi] = max(pointZone[pointi],maxZone);
            }

            const labelList& cPts =  mesh_.cellPoints()[splitCellI];
            label cLevel = cellLevel[splitCellI];
            const point& cc = cellCentres[splitCellI];
            forAll(cPts, cPtI)
            {
                label pointi = cPts[cPtI];
                if (pointLevel[pointi] <= cLevel)
                {
                    checkCells.append(splitCellI);
                    checkFaces.append(facei);
                    point midPt = 0.5*(cc + mesh_.points()[pointi]);

                    start.append(midPt);
                    end.append(nbrCC);
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pointZone,
        maxEqOp<label>(),
        label(-2)
    );

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    vectorField normal1;

    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;
    vectorField normal2;

    surfaces_.findNearestIntersection
    (
        surfacesToTest,
        start,
        end,

        surface1,
        hit1,
        region1,
        normal1,

        surface2,
        hit2,
        region2,
        normal2
    );

    labelList nNonHitsPerFace(mesh_.nFaces(), 0);
    labelList nHitsPerFace(mesh_.nFaces(), 0);
    pointField hitNormals(mesh_.nFaces(), vector::zero);
    forAll(start, i)
    {
        label facei = checkFaces[i];
        if (!hit1[i].hit())
        {
            nNonHitsPerFace[facei]++;
        }
        else
        {
            nHitsPerFace[facei]++;
            vector norm = normal1[i];
            norm /= (mag(norm) + SMALL);
            if (mag(hitNormals[facei]) > SMALL)
            {
                if ((norm & hitNormals[facei]) < 0)
                {
                    norm = -norm;
                }
            }
            hitNormals[facei] += norm;
        }
    }

    labelList keepCell(mesh_.nCells(), -2);
    forAll(start, i)
    {
        label facei = checkFaces[i];
        if (nNonHitsPerFace[facei] > 3 && nHitsPerFace[facei] > 0)
        {
            vector fN = mesh_.faceAreas()[facei];
            fN /= (mag(fN) + SMALL);
            hitNormals[facei] /= nHitsPerFace[facei];
            hitNormals[facei] /= (mag(hitNormals[facei]) + SMALL);
            if (mag(fN&hitNormals[facei]) < 0.34)
            {
                label celli = checkCells[i];
                label maxPZone = labelMin;
                const face& f = mesh_.faces()[facei];
                forAll(f,fp)
                {
                    maxPZone = max(maxPZone, pointZone[f[fp]]);
                }
                keepCell[celli] = maxPZone;
            }
        }
    }

    label defaultRegion = surfaces_.globalRegion(surfacesToTest[0],label(0));

    resetRegionFaces
    (
        keepCell,
        cellToZone,
        globalRegion1,
        globalRegion2,
        testFaces,
        defaultRegion
    );
}


void Foam::meshRefinement::reallocateCornerCells
(
    const refinementParameters& refineParams,
    labelList& globalRegion1,
    labelList& globalRegion2,
    labelList& cellToZone,
    labelList& testFaces
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh_));
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    const labelList& cellLevel = meshCutter_.cellLevel();
    const scalar edge0Len = meshCutter_.level0EdgeLength();

    labelList surfacesToTest(surfaces_.cornerCellSurfaces());

    if (surfacesToTest.size() == 0)
    {
        return;
    }
    else
    {
        Info<<"Checking for corner cells to keep"<<endl;
    }

    label nGrowthIter = 0;

    labelList keepCell(mesh_.nCells(), -2);

    if (refineParams.cornerRemoveBaffles())
    {
        const boolList& cornerCheckRegions = surfaces_.cornerCellRegions();
        labelList neiCellZone;
        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
        forAll(mesh_.faces(), facei)
        {
            label region1 = globalRegion1[facei];
            label region2 = globalRegion2[facei];
            if (region1 == -1 || region2 == -1)
            {
                continue;
            }

            label patchi = patches.whichPatch(facei);
            label own = mesh_.faceOwner()[facei];

            label ownZone = cellToZone[own];
            label neiZone = labelMin;

            if (patchi == -1)
            {
                neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
            }
            else if (patches[patchi].coupled())
            {
                neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
            }

            if
            (
                (ownZone != -2 && ownZone == neiZone)
                && (cornerCheckRegions[region1] || cornerCheckRegions[region2])
            )
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
            }
        }
    }

    vector axisI = vector(1, 0, 0);
    vector axisJ = vector(0, 1, 0);

    List<direction> faceDir(mesh_.nFaces(),direction(0));
    forAll(mesh_.faces(), facei)
    {
        vector fA = mesh_.faceAreas()[facei];
        fA /= (mag(fA) + SMALL);
        if (mag(axisI & fA) > 0.707)
        {
            faceDir[facei] = direction(0);
        }
        else if (mag(axisJ & fA) > 0.707)
        {
            faceDir[facei] = direction(1);
        }
        else
        {
            faceDir[facei] = direction(2);
        }
    }
    label defaultRegion = surfaces_.globalRegion(surfacesToTest[0],label(0));
    label maxCornerIter = refineParams.maxCornerIter();
    while (nGrowthIter < maxCornerIter)
    {
        nGrowthIter++;

        // Get coupled neighbour cellZone
        labelList neiCellZone;
        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

        labelList faceType(mesh_.nFaces(), -1);

        faceSet baffledFacesSet(mesh_, "baffledFacesSet", mesh_.nFaces()/10);
        labelList pointZone(mesh_.nPoints(), -2);

        pointField shellOuterDir(mesh_.nPoints(), vector::zero);
        labelList nShellPtFaces(mesh_.nPoints(), 0);

        forAll(mesh_.faces(), facei)
        {
            label region1 = globalRegion1[facei];
            label region2 = globalRegion2[facei];
            if (region1 == -1 || region2 == -1)
            {
                continue;
            }

            label patchi = patches.whichPatch(facei);
            label own = mesh_.faceOwner()[facei];

            label ownZone = cellToZone[own];
            label neiZone = labelMin;

            if (patchi == -1)
            {
                neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
            }
            else if (patches[patchi].coupled())
            {
                neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
            }

            if
            (
                (ownZone == -2 && neiZone > -2)
                || (neiZone == -2 && ownZone > -2)
            )
            {
                baffledFacesSet.insert(facei);
                faceType[facei] = 0;
                const face& f = mesh_.faces()[facei];
                label maxZone = max(ownZone,neiZone);

                point fN = mesh_.faceAreas()[facei];
                fN /= mag(fN);
                if (ownZone == -2)
                {
                    fN = -fN;
                }

                forAll(f,fp)
                {
                    label pointi = f[fp];
                    pointZone[pointi] = max(pointZone[pointi],maxZone);
                    if (isMasterFace[facei])
                    {
                        shellOuterDir[pointi] += fN;
                        nShellPtFaces[pointi]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            shellOuterDir,
            plusEqOp<point>(),
            vector::zero
        );

        syncTools::syncPointList
        (
            mesh_,
            nShellPtFaces,
            plusEqOp<label>(),
            label(0)
        );

        syncTools::syncPointList
        (
            mesh_,
            pointZone,
            maxEqOp<label>(),
            label(-2)
        );
        labelHashSet usedPoints(mesh_.nPoints()/100);
        labelHashSet faceLabels(baffledFacesSet.toc());

        forAllConstIter(labelHashSet, faceLabels, iter)
        {
            usedPoints.insert(mesh_.faces()[iter.key()]);
        }

        pointSet baffledPtsSet
        (
            mesh_,
            "baffledPtsSet",
            usedPoints
        );
        baffledPtsSet.sync(mesh_);
        label nMarked = 0;

        boolList bufferCells(mesh_.nCells(), false);
        forAllConstIter(pointSet, baffledPtsSet, iter)
        {
            const labelList& pCells = mesh_.pointCells(iter.key());

            forAll(pCells, pCelli)
            {
                label celli = pCells[pCelli];

                if (cellToZone[celli] != -2)
                {
                    continue;
                }
                bufferCells[celli] = true;
                const cell& cFaces = mesh_.cells()[celli];

                forAll(cFaces, cFacei)
                {
                    label facei = cFaces[cFacei];
                    if (faceType[facei] == -1)
                    {
                        faceType[facei] = 1;
                        nMarked++;
                     }
                }
            }
        }

        forAllConstIter(pointSet, baffledPtsSet, iter)
        {
            const labelList& pFaces = mesh_.pointFaces()[iter.key()];

            forAll(pFaces, pFacei)
            {
                label facei = pFaces[pFacei];
                if (faceType[facei] == 1)
                {
                    faceType[facei] = 2;
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh_,
            faceType,
            maxEqOp<label>()
        );

        DynamicList<label> markedFaces(nMarked);
        forAll(faceType, facei)
        {
            if (faceType[facei] > 0)
            {
                markedFaces.append(facei);
            }
        }

        pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());

        pointField perturbCellCentres = mesh_.cellCentres();

        forAll(mesh_.cells(), celli)
        {
            if (bufferCells[celli])
            {
                const labelList& cPts = mesh_.cellPoints()[celli];

                label nOuterPts = 0;
                point aveDir = vector::zero;
                forAll(cPts, cPI)
                {
                    label pointi = cPts[cPI];
                    if (nShellPtFaces[pointi] > 0)
                    {
                        nOuterPts++;
                        aveDir += (shellOuterDir[pointi]/nShellPtFaces[pointi]);
                    }
                }

                if (nOuterPts > 0)
                {
                    aveDir /= nOuterPts;
                    scalar magDisp = mag(aveDir);
                    if (magDisp > SMALL)
                    {
                        point cc = mesh_.cellCentres()[celli];
                        aveDir /= magDisp;
                        scalar clength = edge0Len / pow(2., cellLevel[celli]);
                        perturbCellCentres[celli] = cc + 0.01*clength*aveDir;
                    }
                }
            }
        }

        calcPerturbedNeighbourData
        (
            mesh_.cellCentres(),
            perturbCellCentres,
            neiLevel,
            neiCc
        );

        pointField start(markedFaces.size());
        pointField end(markedFaces.size());
        {
            labelList minLevel;
            calcCellCellRays
            (
                perturbCellCentres,
                neiCc,
                labelList(neiCc.size(), -1),
                markedFaces,
                start,
                end,
                minLevel
             );
        }

        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        vectorField normal1;

        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;
        vectorField normal2;

        surfaces_.findNearestIntersection
        (
            surfacesToTest,
            start,
            end,

            surface1,
            hit1,
            region1,
            normal1,

            surface2,
            hit2,
            region2,
            normal2
        );

        DynamicList<labelPair> gapFaces(markedFaces.size());
        List<DynamicList<Tuple2<point,vector>>> singleHit(mesh_.nCells());

        scalarField minCCFaceDirDistance(mesh_.nFaces(), GREAT);
        scalarField minHitDistance(mesh_.nFaces(), GREAT);

        forAll(markedFaces, i)
        {
            label facei = markedFaces[i];

            if (surface1[i] != -1 && surface2[i] != -1)
            {
                label sRegionToGlobal1 =
                    surfaces_.globalRegion(surface1[i],region1[i]);
                label sRegionToGlobal2 =
                    surfaces_.globalRegion(surface2[i],region2[i]);

                if
                (
                    !surfaces_.cornerCellRegions()[sRegionToGlobal1]
                    && !surfaces_.cornerCellRegions()[sRegionToGlobal2]
                )
                {
                    continue;
                }

                vector hitVec = hit2[i].hitPoint() - hit1[i].hitPoint();

                if (mag(hitVec) > SMALL)
                {
                    if
                    (
                        faceType[facei] != 2
                        ||
                        (
                            nGrowthIter > 1
                            && mag(normal1[i] & normal2[i]) < 0.707
                        )
                    )
                    {
                        continue;
                    }

                    point hitMid = 0.5*(hit2[i].hitPoint()+hit1[i].hitPoint());
                    scalar nDist = -1;

                    if (nGrowthIter > 1)
                    {
                        nDist = 0.5*
                            (mag(normal1[i] & hitVec)+mag(normal2[i] & hitVec));
                    }
                    else
                    {
                        nDist = max
                        (
                            mag(normal1[i] & hitVec),
                            mag(normal2[i] & hitVec)
                        );
                    }

                    label patchi = patches.whichPatch(facei);
                    label own = mesh_.faceOwner()[facei];

                    label minLevel = labelMax;
                    minLevel = min(minLevel, cellLevel[own]);

                    if (patchi == -1)
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        minLevel = min(minLevel, cellLevel[nei]);
                    }
                    else if (patches[patchi].coupled())
                    {
                        minLevel = min
                        (
                            minLevel,
                            neiLevel[facei-mesh_.nInternalFaces()]
                        );
                    }

                    scalar clength = edge0Len / pow(2., minLevel);

                    if (nDist > 0.2*clength)
                    {
                        gapFaces.append(labelPair(facei,i));

                        DynamicList<label> fCells(2);
                        fCells.append(own);
                        if (patchi == -1)
                        {
                            label nei = mesh_.faceNeighbour()[facei];
                            fCells.append(nei);
                        }

                        forAll(fCells, fCI)
                        {
                            label celli = fCells[fCI];
                            const cell& c = mesh_.cells()[celli];
                            const point& cc = mesh_.cellCentres()[celli];

                            scalar ccDir = -1;
                            if (faceDir[facei] == 0)
                            {
                                ccDir = cc.x();
                            }
                            else if (faceDir[facei] == 1)
                            {
                                ccDir = cc.y();
                            }
                            else
                            {
                                ccDir = cc.z();
                            }
                            scalar hDist = mag(cc-hitMid);
                            forAll(c, cFI)
                            {
                                label cfacei = c[cFI];

                                if (faceType[cfacei] == 0)
                                {
                                    minCCFaceDirDistance[facei] = min
                                    (
                                        minCCFaceDirDistance[facei],
                                        ccDir
                                    );
                                    minHitDistance[facei] = min
                                    (
                                        minHitDistance[facei],
                                        hDist
                                    );

                                    break;
                                }
                            }
                        }
                    }
                }
                else
                {
                    label patchi = patches.whichPatch(facei);
                    label own = mesh_.faceOwner()[facei];
                    Tuple2<point,vector> hInfo(hit1[i].hitPoint(),normal1[i]);

                    singleHit[own].append(hInfo);
                    if (patchi == -1)
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        singleHit[nei].append(hInfo);
                    }
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh_,
            minCCFaceDirDistance,
            minEqOp<scalar>()
        );

        syncTools::syncFaceList
        (
            mesh_,
            minHitDistance,
            minEqOp<scalar>()
        );

        labelList nClosestCells(mesh_.nFaces(), 0);

        forAll(gapFaces, i)
        {
            label facei = gapFaces[i].first();

            label markedIndex = gapFaces[i].second();
            label patchi = patches.whichPatch(facei);
            DynamicList<label> fCells(2);
            label own = mesh_.faceOwner()[facei];
            fCells.append(own);
            if (patchi == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                fCells.append(nei);
            }

            point hitMid =
                0.5*(hit2[markedIndex].hitPoint()+hit1[markedIndex].hitPoint());
            forAll(fCells, fCI)
            {
                label celli = fCells[fCI];
                const cell& c = mesh_.cells()[celli];

                forAll(c, cFI)
                {
                    label cfacei = c[cFI];

                    if (faceType[cfacei] == 0)
                    {
                        const point& cc = mesh_.cellCentres()[celli];
                        scalar hDist = mag(cc-hitMid);
                        if (hDist < (minHitDistance[facei]+SMALL))
                        {
                            nClosestCells[facei]++;
                        }
                        break;
                    }

                }
            }
        }
        syncTools::syncFaceList
        (
            mesh_,
            nClosestCells,
            plusEqOp<label>()
        );

        keepCell = -2;
        forAll(gapFaces, i)
        {
            label facei = gapFaces[i].first();
            label markedIndex = gapFaces[i].second();
            label patchi = patches.whichPatch(facei);
            const face& f = mesh_.faces()[facei];
            label maxPZone = labelMin;

            forAll(f,fp)
            {
                maxPZone = max(maxPZone, pointZone[f[fp]]);
            }

            DynamicList<label> fCells(2);
            label own = mesh_.faceOwner()[facei];
            fCells.append(own);
            if (patchi == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                fCells.append(nei);
            }

            scalar minCCDir =  minCCFaceDirDistance[facei];

            point hitMid =
                0.5*(hit2[markedIndex].hitPoint()+hit1[markedIndex].hitPoint());

            forAll(fCells, fCI)
            {
                label celli = fCells[fCI];
                const cell& c = mesh_.cells()[celli];
                const point& cc = mesh_.cellCentres()[celli];
                scalar ccDir = -1;
                if (faceDir[facei] == 0)
                {
                    ccDir = cc.x();
                }
                else if (faceDir[facei] == 1)
                {
                    ccDir = cc.y();
                }
                else
                {
                    ccDir = cc.z();
                }

                forAll(c, cFI)
                {
                    label cfacei = c[cFI];
                    if (nClosestCells[facei] == 2)
                    {
                        if (faceType[cfacei] == 0)
                        {
                            if (ccDir < (minCCDir + SMALL))
                            {
                                keepCell[celli] = maxPZone;
                            }
                            break;
                        }
                    }
                    else
                    {
                        if (faceType[cfacei] == 0)
                        {
                            scalar hDist = mag(cc-hitMid);
                            if (hDist < (minHitDistance[facei]+SMALL))
                            {
                                keepCell[celli] = maxPZone;
                            }
                            break;
                        }
                    }
                }
            }
        }

        forAll(mesh_.cells(), celli)
        {
            if (singleHit[celli].size() > 1 && cellToZone[celli] == -2)
            {
                scalar clength = edge0Len / pow(2., cellLevel[celli]);
                const DynamicList<Tuple2<point,vector>>& hits
                    = singleHit[celli];

                bool foundGap = false;
                for (label i = 0; i < hits.size(); i++)
                {
                    for (label j = i+1; j < hits.size(); j++)
                    {
                        vector hp1 = hits[i].first();
                        vector hp2 = hits[j].first();

                        vector np1 = hits[i].second();
                        vector np2 = hits[j].second();

                        vector separationVec = hp1-hp2;
                        scalar separationDist = -1;

                        if (nGrowthIter > 1)
                        {
                            separationDist = 0.5*
                            (
                                mag(separationVec & np1) +
                                mag(separationVec & np2)
                            );
                        }
                        else
                        {
                            separationDist = max
                            (
                                mag(separationVec & np1),
                                mag(separationVec & np2)
                            );
                        }


                        if
                        (
                            separationDist > 0.2*clength
                            && separationDist < clength
                        )
                        {
                            if (mag(np1&np2) > 0.707)
                            {
                                foundGap = true;
                            }
                            else if (nGrowthIter == 1)
                            {
                                foundGap = true;
                            }
                        }
                    }
                }
                if (foundGap)
                {
                    const labelList& cPts = mesh_.cellPoints()[celli];

                    label maxCellZone = -2;
                    forAll(cPts, cPI)
                    {
                        label pointi = cPts[cPI];
                        maxCellZone = max(maxCellZone, pointZone[pointi]);
                    }
                    keepCell[celli] = maxCellZone;
                }
            }
        }

        labelList checkEdges = selectPinchedEdges(cellToZone,keepCell);
        List<pointField> eCC(mesh_.nEdges(), pointField(0));
        forAll(checkEdges, cEI)
        {
            label edgei = checkEdges[cEI];
            const labelList& eCells = mesh_.edgeCells()[edgei];
            pointField checkPts(2);
            label nPts = 0;
            forAll(eCells, eCI)
            {
                label celli = eCells[eCI];
                if (keepCell[celli] == -2 && cellToZone[celli] == -2)
                {
                    checkPts[nPts++] = mesh_.cellCentres()[celli];
                }
            }
            checkPts.setSize(nPts);
            eCC[edgei] = checkPts;
        }

        syncTools::syncEdgeList
        (
            mesh_,
            eCC,
            pointFieldCombine(),
            pointField() // initial value
        );

        forAll(mesh_.cells(), celli)
        {
            if (keepCell[celli] != -2)
            {
                const labelList& cPts = mesh_.cellPoints()[celli];
                forAll(cPts, cPtI)
                {
                    label pointi = cPts[cPtI];
                    pointZone[pointi] = max(pointZone[pointi],keepCell[celli]);
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            pointZone,
            maxEqOp<label>(),
            label(-2)
        );

        forAll(checkEdges, cEI)
        {
            label edgei = checkEdges[cEI];
            const labelList& eCells = mesh_.edgeCells()[edgei];
            const pointField& nCC = eCC[edgei];

            point minPt(GREAT, GREAT, GREAT);

            if (nCC.size() > 1)
            {
                if ((nCC[0].x()+SMALL) < nCC[1].x())
                {
                    minPt = nCC[0];
                }
                else if ((nCC[1].x()+SMALL) < nCC[0].x())
                {
                    minPt = nCC[1];
                }
                else if ((nCC[0].y()+SMALL) < nCC[1].y())
                {
                    minPt = nCC[0];
                }
                else if ((nCC[1].y()+SMALL) < nCC[0].y())
                {
                    minPt = nCC[1];
                }
                else if ((nCC[0].z()+SMALL) < nCC[1].z())
                {
                    minPt = nCC[0];
                }
                else if ((nCC[1].z()+SMALL) < nCC[0].z())
                {
                    minPt = nCC[1];
                }

                forAll(eCells, eCI)
                {
                    label celli = eCells[eCI];
                    {
                        point cc = mesh_.cellCentres()[celli];
                        if (mag(minPt-cc) < SMALL)
                        {
                            const edge& e = mesh_.edges()[edgei];
                            keepCell[celli] =
                                max(pointZone[e[0]],pointZone[e[1]]);
                        }
                    }
                }
            }
            else
            {
                forAll(eCells, eCI)
                {
                    label celli = eCells[eCI];
                    if (keepCell[celli] == -2 && cellToZone[celli] == -2)
                    {
                        const edge& e = mesh_.edges()[edgei];
                        keepCell[celli] =
                            max(pointZone[e[0]],pointZone[e[1]]);
                    }
                }
            }
        }

        label nUpdated = resetRegionFaces
        (
            keepCell,
            cellToZone,
            globalRegion1,
            globalRegion2,
            testFaces,
            defaultRegion
        );

        if (returnReduce(nUpdated, sumOp<label>()) == 0)
        {
            break;
        }
    }

    keepCell = -2;

    labelList nEdgeCells(mesh_.nEdges(), 0);
    labelList nPointCells(mesh_.nPoints(), 0);
    forAll(mesh_.edges(), edgei)
    {
        const labelList& eCells = mesh_.edgeCells()[edgei];
        nEdgeCells[edgei] += eCells.size();
    }
    forAll(mesh_.points(), pointi)
    {
        const labelList& pCells = mesh_.pointCells()[pointi];
        nPointCells[pointi] += pCells.size();
    }
    syncTools::syncEdgeList
    (
        mesh_,
        nEdgeCells,
        plusEqOp<label>(),
        label(0) // initial value
    );
    syncTools::syncPointList
    (
        mesh_,
        nPointCells,
        plusEqOp<label>(),
        label(0) // initial value
    );

    labelList nEdgeCellsOutside(mesh_.nEdges(), 0);
    forAll(mesh_.edges(), edgei)
    {
        const labelList& eCells = mesh_.edgeCells()[edgei];
        forAll(eCells, eCI)
        {
            label celli = eCells[eCI];
            if (cellToZone[celli] == -2)
            {
                nEdgeCellsOutside[edgei]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        nEdgeCellsOutside,
        plusEqOp<label>(),
        label(0) // initial value
    );

    forAll(mesh_.edges(), edgei)
    {
        nEdgeCellsOutside[edgei] += (label(4)-nEdgeCells[edgei]);
    }

    labelList nPointCellsOutside(mesh_.nPoints(), 0);
    forAll(mesh_.points(), pointi)
    {
        const labelList& pCells = mesh_.pointCells()[pointi];
        forAll(pCells, pCI)
        {
            label celli = pCells[pCI];
            if (cellToZone[celli] == -2)
            {
                nPointCellsOutside[pointi]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        nPointCellsOutside,
        plusEqOp<label>(),
        label(0) // initial value
    );

    forAll(mesh_.points(), pointi)
    {
        nPointCellsOutside[pointi] += (label(8)-nPointCells[pointi]);
    }

    boolList edgeConnectedOutside(mesh_.nPoints(), false);

    forAll(mesh_.points(), pointi)
    {
        const labelList& pEdges = mesh_.pointEdges()[pointi];
        label nFound = 0;
        forAll(pEdges, pEI)
        {
            label edgei = pEdges[pEI];
            if (nEdgeCellsOutside[edgei] == 1)
            {
                nFound++;
            }
        }

        if (nFound == pEdges.size())
        {
            edgeConnectedOutside[pointi] = true;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        edgeConnectedOutside,
        andEqOp<bool>(),
        false                 // initial value
    );

    List<pointField> pCC(mesh_.nPoints(), pointField(0));
    labelList pointZoneID(mesh_.nPoints(), -2);
    forAll(mesh_.points(), pointi)
    {
        if
        (
            nPointCellsOutside[pointi] == 2
            && edgeConnectedOutside[pointi]
        )
        {
            const labelList& pCells = mesh_.pointCells()[pointi];
            DynamicList<point> checkPts(2);
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                pointZoneID[pointi] = max(pointZoneID[pointi],cellToZone[celli]);
                if (keepCell[celli] == -2 && cellToZone[celli] == -2)
                {
                    checkPts.append(mesh_.cellCentres()[celli]);
                }
            }
            pCC[pointi] = checkPts;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pointZoneID,
        maxEqOp<label>(),
        label(-2) // initial value
    );

    syncTools::syncPointList
    (
        mesh_,
        pCC,
        pointFieldCombine(),
        pointField() // initial value
    );

    forAll(mesh_.points(), pointi)
    {
        if
        (
            nPointCellsOutside[pointi] == 2
            && edgeConnectedOutside[pointi]
        )
        {
            const pointField& nCC = pCC[pointi];

            if (nCC.size() > 1)
            {
                point minPt(GREAT, GREAT, GREAT);

                if ((nCC[0].x()+SMALL) < nCC[1].x())
                {
                    minPt = nCC[0];
                }
                else if ((nCC[1].x()+SMALL) < nCC[0].x())
                {
                    minPt = nCC[1];
                }
                else if ((nCC[0].y()+SMALL) < nCC[1].y())
                {
                    minPt = nCC[0];
                }
                else if ((nCC[1].y()+SMALL) < nCC[0].y())
                {
                    minPt = nCC[1];
                }
                else if ((nCC[0].z()+SMALL) < nCC[1].z())
                {
                    minPt = nCC[0];
                }
                else if ((nCC[1].z()+SMALL) < nCC[0].z())
                {
                    minPt = nCC[1];
                }

                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];

                    if (keepCell[celli] == -2 && cellToZone[celli] == -2)
                    {
                        point cc = mesh_.cellCentres()[celli];
                        if (mag(minPt-cc) < SMALL)
                        {
                            keepCell[celli] = pointZoneID[pointi];
                        }
                    }
                }
            }
            else
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];

                    if (keepCell[celli] == -2 && cellToZone[celli] == -2)
                    {
                        keepCell[celli] = pointZoneID[pointi];
                    }
                }
            }
        }
    }

    resetRegionFaces
    (
        keepCell,
        cellToZone,
        globalRegion1,
        globalRegion2,
        testFaces,
        defaultRegion
    );

    if (debug)
    {
        const_cast<Time&>(mesh_.time())++;
        mesh_.write();

        volScalarField volCellToRegion
        (
            IOobject
            (
                "cellToRegion",
                timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
             ),
            mesh_,
            dimensionedScalar("zero", dimless, -1),
            zeroGradientFvPatchScalarField::typeName
         );

        forAll(mesh_.cells(), celli)
        {
            if (cellToZone[celli] != -2)
            {
                volCellToRegion[celli] = scalar(0);
            }
        }

        volCellToRegion.write();
    }

    return;
}


// Switch cells from removed to kept if all anchors boundary.
void Foam::meshRefinement::reallocateSingleRemovalCells
(
    labelList& globalRegion1,
    labelList& globalRegion2,
    labelList& cellToZone,
    labelList& testFaces
) const
{
    const labelList& pointLevel = meshCutter_.pointLevel();
    const labelList& cellLevel = meshCutter_.cellLevel();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get coupled neighbour cellZone
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    boolList isBoundaryPoint(mesh_.nPoints(), false);
    forAll(mesh_.faces(), facei)
    {
        if (globalRegion1[facei] != -1)
        {
            face f = mesh_.faces()[facei];
            forAll(f, fp)
            {
                isBoundaryPoint[f[fp]] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false               // null value
    );

    // Get per face whether is it master (of a coupled set of faces)
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    labelList nEdgeBoundaryFaces(mesh_.nEdges(), 0);
    forAll(mesh_.edges(), edgei)
    {
        const labelList& eFaces = mesh_.edgeFaces()[edgei];
        forAll(eFaces, eFI)
        {
            label facei = eFaces[eFI];
            if (isMasterFace[facei] && globalRegion1[facei] != -1)
            {
                label patchi = patches.whichPatch(facei);
                label own = mesh_.faceOwner()[facei];
                label ownZone = cellToZone[own];
                label neiZone = -2;
                if (patchi == -1)
                {
                    neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
                }
                else if (patches[patchi].coupled())
                {
                    neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                }

                if
                (
                    (ownZone == -2 && neiZone >= -1)
                    || (ownZone >= -1 && neiZone == -2)
                )
                {
                    nEdgeBoundaryFaces[edgei]++;
                }
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        nEdgeBoundaryFaces,
        plusEqOp<label>(),
        label(0)                  // initial value
    );

    boolList markedKeeping(mesh_.nCells(), false);
    forAll(cellLevel, celli)
    {
        const cell& c = mesh_.cells()[celli];

        if (cellToZone[celli] != -2 || c.size() != 6)
        {
            continue;
        }

        const labelList& cPoints = mesh_.cellPoints()[celli];

        // Get number of anchor points (pointLevel <= cellLevel)

        label nBoundaryAnchors = 0;
        forAll(cPoints, i)
        {
            label pointI = cPoints[i];
            if (pointLevel[pointI] <= cellLevel[celli])
            {
                // Anchor point
                if (isBoundaryPoint[pointI])
                {
                    nBoundaryAnchors++;
                }
            }
        }

        if (nBoundaryAnchors == 8)
        {
            label nkeptNbrs = 0;
            label openFace = -1;

            label minZone = labelMax;
            label maxZone = labelMin;
            forAll(c,cfi)
            {
                label facei = c[cfi];
                if (globalRegion1[facei] == -1)
                {
                    openFace = facei;
                    continue;
                }
                label patchi = patches.whichPatch(facei);
                label own = mesh_.faceOwner()[facei];
                if (patchi == -1)
                {
                    label neiZone
                    (
                        own == celli ?
                        cellToZone[mesh_.faceNeighbour()[facei]]
                        : cellToZone[own]
                    );

                    if (neiZone >= -1)
                    {
                        minZone = min(minZone,neiZone);
                        maxZone = max(maxZone,neiZone);
                        nkeptNbrs++;
                    }
                    else
                    {
                        openFace = facei;
                    }
                }
                else if (patches[patchi].coupled())
                {
                    label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                    if (neiZone >= -1)
                    {
                        minZone = min(minZone,neiZone);
                        maxZone = max(maxZone,neiZone);
                        nkeptNbrs++;
                    }
                    else
                    {
                        openFace = facei;
                    }
                }
                else
                {
                    openFace = facei;
                }
            }

            if (nkeptNbrs == 5 && openFace != -1 && minZone == maxZone)
            {
                vector openNorm = mesh_.faceAreas()[openFace];
                openNorm /= (mag(openNorm) + SMALL);
                label topFace = -1;

                forAll(c,cfi)
                {
                    label facei = c[cfi];
                    if (facei != openFace)
                    {
                        vector norm = mesh_.faceAreas()[facei];
                        norm /= (mag(norm) + SMALL);
                        if (mag(norm&openNorm) > 0.707)
                        {
                            topFace = facei;
                            break;
                        }
                    }
                }

                if (topFace != -1)
                {
                    const labelList& topEdges = mesh_.faceEdges()[topFace];
                    bool keepCell = true;

                    forAll(topEdges, fEI)
                    {
                        if (nEdgeBoundaryFaces[topEdges[fEI]] != 2)
                        {
                            keepCell = false;
                            break;
                        }
                    }

                    if (keepCell)
                    {
                        const labelList& openEdges =
                            mesh_.faceEdges()[openFace];
                        const labelList& cEdges = mesh_.cellEdges()[celli];

                        labelHashSet openEdgesSet(openEdges);
                        labelHashSet topEdgesSet(topEdges);
                        DynamicList<label> sideEdges(cEdges.size());

                        forAll(cEdges, cEI)
                        {
                            label edgei = cEdges[cEI];
                            if (!openEdgesSet.found(edgei)&&!topEdgesSet(edgei))
                            {
                                sideEdges.append(edgei);
                            }
                        }
                        label nSideMan = 0;
                        forAll(sideEdges, sei)
                        {
                            if (nEdgeBoundaryFaces[sideEdges[sei]] == 2)
                            {
                                nSideMan++;
                            }
                        }
                        if (nSideMan == sideEdges.size())
                        {
                            markedKeeping[celli] = true;
                        }
                    }
                }
            }
        }
    }

    labelList resetFaces(mesh_.nFaces(), -2);
    boolList neiMarkedKeeping;
    syncTools::swapBoundaryCellList(mesh_, markedKeeping, neiMarkedKeeping);

    label nMarked = 0;

    forAll(cellLevel, celli)
    {
        label nNeiMarked = 0;
        if (markedKeeping[celli])
        {
            const cell& c = mesh_.cells()[celli];
            label maxZone = -2;
            forAll(c,cfi)
            {
                label facei = c[cfi];

                label patchi = patches.whichPatch(facei);
                label own = mesh_.faceOwner()[facei];
                if (patchi == -1)
                {
                    label neiZone
                    (
                        own == celli ?
                        cellToZone[mesh_.faceNeighbour()[facei]]
                        : cellToZone[own]
                    );
                    maxZone = max(maxZone, neiZone);

                    if (own == celli)
                    {
                        if (markedKeeping[mesh_.faceNeighbour()[facei]])
                        {
                            nNeiMarked++;
                        }
                    }
                    else
                    {
                        if (markedKeeping[own])
                        {
                            nNeiMarked++;
                        }
                    }
                }
                else if (patches[patchi].coupled())
                {
                    label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                    maxZone = max(maxZone, neiZone);

                    if (neiMarkedKeeping[facei-mesh_.nInternalFaces()])
                    {
                        nNeiMarked++;
                    }
                }
            }

            if (nNeiMarked == 0 && maxZone > -2)
            {
                nMarked++;
                cellToZone[celli] = maxZone;
                forAll(c,cfi)
                {
                    label facei = c[cfi];

                    label patchi = patches.whichPatch(facei);
                    label own = mesh_.faceOwner()[facei];
                    if (patchi == -1)
                    {
                        label neiZone
                        (
                            own == celli ?
                            cellToZone[mesh_.faceNeighbour()[facei]]
                            : cellToZone[own]
                        );

                        if (neiZone >= -1)
                        {
                            resetFaces[facei] = -1;
                        }
                        else
                        {
                            resetFaces[facei] = 0;
                        }
                    }
                    else if (patches[patchi].coupled())
                    {
                        label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                        if (neiZone >= -1)
                        {
                            resetFaces[facei] = -1;
                        }
                        else
                        {
                            resetFaces[facei] = 0;
                        }
                    }
                }
            }
        }
    }

    nMarked = returnReduce(nMarked, sumOp<label>());

    if (nMarked > 0)
    {
        Info<<"Re-allocated " << nMarked << " isolated cells "<<endl;
        syncTools::syncFaceList
        (
            mesh_,
            resetFaces,
            maxEqOp<label>()
        );

        DynamicList<label> newBaffleFaces(nMarked);
        forAll(resetFaces, facei)
        {
            if (resetFaces[facei] == -1)
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
            }
            else if (resetFaces[facei] == 0)
            {
                newBaffleFaces.append(facei);
                globalRegion1[facei] = 0;
                globalRegion2[facei] = 0;
            }
        }

        syncTools::syncFaceList
        (
            mesh_,
            globalRegion1,
            maxEqOp<label>()
        );

        syncTools::syncFaceList
        (
            mesh_,
            globalRegion2,
            maxEqOp<label>()
         );

         //Ensure testFaces includes updates
        if (newBaffleFaces.size() > 0)
        {
            label sz = testFaces.size();
            testFaces.setSize(sz+newBaffleFaces.size());
            forAll(newBaffleFaces, i)
            {
                testFaces[sz+i] = newBaffleFaces[i];
            }
        }
    }
}


// Zoned cells with no or single boundary anchor points can result in
// poor snapping. This function reallocates such cells to the neighboring zone
void Foam::meshRefinement::zoneProblemCellReallocate
(
    const labelList& globalRegion1,
    const labelList& globalRegion2,
    labelList& namedSurfaceIndex,
    labelList& cellToZone
) const
{
    if (returnReduce(namedSurfaceIndex.size(), sumOp<label>()) == 0)
    {
        return;
    }

    const labelList& pointLevel = meshCutter_.pointLevel();
    const labelList& cellLevel = meshCutter_.cellLevel();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get coupled neighbour cellZone
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    boolList isBoundaryPoint(mesh_.nPoints(), false);
    boolList isBoundaryEdge(mesh_.nEdges(), false);
    forAll(mesh_.faces(), facei)
    {
        label patchi = patches.whichPatch(facei);
        if
        (
            (patchi != -1 && !patches[patchi].coupled())
            || globalRegion1[facei] != -1
            || namedSurfaceIndex[facei] != -1
        )
        {
            face f = mesh_.faces()[facei];
            forAll(f, fp)
            {
                isBoundaryPoint[f[fp]] = true;
            }
            const labelList& fEdges = mesh_.faceEdges()[facei];
            forAll(fEdges, fEI)
            {
                isBoundaryEdge[fEdges[fEI]] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false               // null value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false               // null value
    );

    //Mark 7 and 8 anchor points
    boolList markedCells(mesh_.nCells(), false);

    // Does cell have exactly 7 of its 8 anchor points on the boundary?
    PackedBoolList hasSevenBoundaryAnchorPoints(mesh_.nCells());
    // Does cell have a single internal edge
    PackedBoolList hasFewInternalEdges(mesh_.nCells());
    // If so what is the remaining non-boundary anchor point?
    labelHashSet nonBoundaryAnchors(mesh_.nCells()/10000);
    labelHashSet nonBoundaryEdges(mesh_.nEdges()/10000);
    DynamicList<label> dynCPoints;

    forAll(cellLevel, celli)
    {
        if (cellToZone[celli] == -2)
        {
            continue;
        }

        const labelList& cPoints = mesh_.cellPoints(celli, dynCPoints);

        // Get number of anchor points (pointLevel <= cellLevel)

        label nBoundaryAnchors = 0;
//        label nNonAnchorBoundary = 0;
        label nonBoundaryAnchor = -1;

        forAll(cPoints, i)
        {
            label pointI = cPoints[i];
            if (pointLevel[pointI] <= cellLevel[celli])
            {
                // Anchor point
                if (isBoundaryPoint[pointI])
                {
                    nBoundaryAnchors++;
                }
                else
                {
                    // Anchor point which is not on the surface
                    nonBoundaryAnchor = pointI;
                }
            }
            else if (isBoundaryPoint[pointI])
            {
//                nNonAnchorBoundary++;
            }
        }

        const labelList& cEdges = mesh_.cellEdges()[celli];
        label numBoundaryEdges = 0;
        label internalEdge = -1;
        forAll(cEdges, cEI)
        {
            label edgei = cEdges[cEI];
            if (isBoundaryEdge[edgei])
            {
                numBoundaryEdges++;
            }
            else
            {
                internalEdge = edgei;
            }
        }

        if (numBoundaryEdges == cEdges.size())
        {
            markedCells[celli] = true;
        }
        else if
        (
            nBoundaryAnchors == 8
            && numBoundaryEdges > cEdges.size()-3
        )
        {
            hasFewInternalEdges.set(celli, 1u);
            nonBoundaryEdges.insert(internalEdge);
         }
        else if (nBoundaryAnchors == 7)
        {
            hasSevenBoundaryAnchorPoints.set(celli, 1u);
            nonBoundaryAnchors.insert(nonBoundaryAnchor);
        }
    }

    // Loop over all points. If a point is connected to 4 or more cells
    // with 7 anchor points on the boundary set those cell's non-boundary faces
    // to baffles
    labelList nCellsSevenAnchors(mesh_.nPoints(), 0);
    DynamicList<label> dynPCells;

    forAllConstIter(labelHashSet, nonBoundaryAnchors, iter)
    {
        label pointI = iter.key();

        const labelList& pCells = mesh_.pointCells(pointI, dynPCells);
        // Count number of 'hasSevenBoundaryAnchorPoints' cells.
        forAll(pCells, i)
        {
            if (hasSevenBoundaryAnchorPoints.get(pCells[i]) == 1u)
            {
                nCellsSevenAnchors[pointI]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        nCellsSevenAnchors,
        plusEqOp<label>(),
        label(0)              // null value
    );

    forAllConstIter(labelHashSet, nonBoundaryAnchors, iter)
    {
        label pointI = iter.key();
        if (nCellsSevenAnchors[pointI] > 3)
        {
            const labelList& pCells = mesh_.pointCells()[pointI];
            forAll(pCells, i)
            {
                label celli = pCells[i];
                if (hasSevenBoundaryAnchorPoints.get(celli) == 1u)
                {
                    markedCells[celli] = true;
                }
            }
        }
    }

    labelList nCellsFewInternalEdge(mesh_.nEdges(), 0);
    forAllConstIter(labelHashSet, nonBoundaryEdges, iter)
    {
        label edgeI = iter.key();
        const labelList& eCells = mesh_.edgeCells()[edgeI];
        forAll(eCells, i)
        {
            if (hasFewInternalEdges.get(eCells[i]) == 1u)
            {
                nCellsFewInternalEdge[edgeI]++;
            }
        }
    }
    syncTools::syncEdgeList
    (
        mesh_,
        nCellsFewInternalEdge,
        plusEqOp<label>(),
        label(0)              // null value
    );

    forAllConstIter(labelHashSet, nonBoundaryEdges, iter)
    {
        label edgeI = iter.key();
        if
        (
            nCellsFewInternalEdge[edgeI] > 0
            || nCellsFewInternalEdge[edgeI] < 3
        )
        {
            const labelList& eCells = mesh_.edgeCells()[edgeI];
            forAll(eCells, i)
            {
                label celli = eCells[i];
                if (hasFewInternalEdges.get(celli) == 1u)
                {
                    markedCells[celli] = true;
                }
            }
        }
    }

    label nSwitched = 0;

    labelList oldCellToZone = cellToZone;
    boolList allMarked = markedCells;
    boolList allNeiMarked;
    syncTools::swapBoundaryCellList(mesh_, allMarked, allNeiMarked);

    forAll(mesh_.cells(), celli)
    {
        if (markedCells[celli])
        {
            markedCells[celli] = false;
            const labelList& cFaces = mesh_.cells()[celli];

            bool foundNamed = false;
            forAll(cFaces, cFI)
            {
                label facei = cFaces[cFI];
                if (namedSurfaceIndex[facei] != -1)
                {
                    foundNamed = true;
                    break;
                }
            }

            if (!foundNamed)
            {
                continue;
            }

            label switchID = -2;
            label newZonedFaces = 0;
            label oldZonedFaces = 0;
            bool foundHigherMarkedNei = false;
            forAll(cFaces, cFI)
            {
                label facei = cFaces[cFI];

                label ownZone = oldCellToZone[celli];
                label neiZone = -2;
                label patchI = patches.whichPatch(facei);
                bool neiMarked = false;
                bool intFace = false;
                if (patchI == -1)
                {
                    label own = mesh_.faceOwner()[facei];
                    label nei
                    (
                        own == celli
                        ? mesh_.faceNeighbour()[facei] : own
                    );
                    neiZone = oldCellToZone[nei];
                    neiMarked = allMarked[nei];
                    intFace = true;
                }
                else if
                (
                    patchI != -1
                    && isA<processorPolyPatch>(patches[patchI])
                )
                {
                    neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                    neiMarked = allNeiMarked[facei-mesh_.nInternalFaces()];
                    intFace = true;
                }

                if
                (
                    intFace
                    && ownZone > -2 && neiZone > -2
                )
                {
                    if (ownZone != neiZone)
                    {
                        if (!neiMarked)
                        {
                            switchID = neiZone;
                        }
                        else if (neiZone > ownZone)
                        {
                            foundHigherMarkedNei = true;
                        }

                        if
                        (
                            globalRegion1[facei] == -1
                            && namedSurfaceIndex[facei] != -1
                         )
                        {
                            oldZonedFaces++;
                        }
                    }
                    else if (ownZone == neiZone)
                    {
                        if
                        (
                            globalRegion1[facei] == -1
                            && namedSurfaceIndex[facei] == -1
                            && !neiMarked
                        )
                        {
                            newZonedFaces++;
                        }
                    }
                }
            }
            if
            (
                switchID != -2 && newZonedFaces < oldZonedFaces
                && !foundHigherMarkedNei
            )
            {
                markedCells[celli] = true;
                cellToZone[celli] = switchID;
                nSwitched++;
            }
        }
    }

    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
    forAll(mesh_.cells(), celli)
    {
        if (markedCells[celli])
        {
            const labelList& cFaces = mesh_.cells()[celli];
            label namedRegion1 = -1;
            forAll(cFaces, cFI)
            {
                label facei = cFaces[cFI];
                if (namedSurfaceIndex[facei] != -1)
                {
                    namedRegion1 = namedSurfaceIndex[facei];
                    break;
                }
            }

            if (namedRegion1 != -1)
            {
                forAll(cFaces, cFI)
                {
                    label facei = cFaces[cFI];

                    if (globalRegion1[facei] != -1)
                    {
                        continue;
                    }

                    label ownZone = cellToZone[celli];
                    label neiZone = -2;

                    label patchI = patches.whichPatch(facei);
                    if (patchI == -1)
                    {
                        label own = mesh_.faceOwner()[facei];
                        label nei
                        (
                            own == celli
                            ? mesh_.faceNeighbour()[facei] : own
                        );
                        neiZone = cellToZone[nei];
                    }
                    else if
                    (
                        patchI != -1
                        && isA<processorPolyPatch>(patches[patchI])
                    )
                    {
                        neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                    }

                    if (ownZone == neiZone)
                    {
                        namedSurfaceIndex[facei] = -1;
                    }
                    else if (namedSurfaceIndex[facei] == -1)
                    {
                        namedSurfaceIndex[facei] = namedRegion1;
                    }
                }
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        namedSurfaceIndex,
        maxEqOp<label>()
    );


    makeConsistentFaceIndex
    (
        cellToZone,
        namedSurfaceIndex
    );

    reduce(nSwitched, sumOp<label>());

    Info<<"Reallocated "<< nSwitched
        <<" zoned cells that might produce snapping issues "<<endl;
}


void Foam::meshRefinement::makeConsistentFaceIndex
(
    const labelList& cellToZone,
    labelList& namedSurfaceIndex
) const
{
    // Make namedSurfaceIndex consistent with cellToZone - clear out any
    // blocked faces inbetween same cell zone (or background (=-1))
    // Do not do this for surfaces relating to 'pure' faceZones i.e.
    // faceZones without a cellZone. Note that we cannot check here
    // for different cellZones on either side but no namedSurfaceIndex
    // since cellZones can now originate from locationsInMesh as well
    // (instead of only through named surfaces)

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    boolList freeStanding(surfZones.size(), false);
    forAll(surfZones, surfI)
    {
        freeStanding[surfI] = surfZones[surfI].freeStanding();
    }

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label ownZone = cellToZone[faceOwner[facei]];
        label neiZone = cellToZone[faceNeighbour[facei]];

        if
        (
            ownZone == neiZone
         && namedSurfaceIndex[facei] != -1
         && !freeStanding[namedSurfaceIndex[facei]]
        )
        {
            namedSurfaceIndex[facei] = -1;
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get coupled neighbour cellZone
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    // Use coupled cellZone to do check
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                if
                (
                    ownZone == neiZone
                 && namedSurfaceIndex[facei] != -1
                 && !freeStanding[namedSurfaceIndex[facei]]
                )
                {
                    namedSurfaceIndex[facei] = -1;
                }
            }
        }
        else
        {
            // Unzonify boundary faces
            forAll(pp, i)
            {
                label facei = pp.start()+i;

                if
                (
                    namedSurfaceIndex[facei] != -1
                 && !freeStanding[namedSurfaceIndex[facei]]
                )
                {
                    namedSurfaceIndex[facei] = -1;
                }
            }
        }
    }
}


void Foam::meshRefinement::getIntersections
(
    const pointField& cellCentres,
    const labelList& surfacesToTest,
    const labelList& testFaces,

    labelList& namedSurfaceIndex,
    PackedBoolList& posOrientation
) const
{
    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(cellCentres, neiLevel, neiCc);

    namedSurfaceIndex.setSize(mesh_.nFaces());
    namedSurfaceIndex = -1;

    posOrientation.setSize(mesh_.nFaces());
    posOrientation = false;

    // Statistics: number of faces per faceZone
    labelList nSurfFaces(surfaces_.surfZones().size(), 0);

    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());
    {
        labelList minLevel;
        calcCellCellRays
        (
            cellCentres,
            neiCc,
            labelList(neiCc.size(), -1),
            testFaces,
            start,
            end,
            minLevel
        );
    }

    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note that we intersect all intersected faces again. Could reuse
    // the information already in surfaceIndex_.

    labelList surface1;
    List<pointIndexHit> hit1;
    vectorField normal1;
    labelList surface2;
    List<pointIndexHit> hit2;
    vectorField normal2;
    {
        labelList region1;
        labelList region2;

        if (checkForCracks())
        {
            surfaces_.findNearestPerturbedIntersection
            (
                surfacesToTest,
                start,
                end,
                crackTol(),
                addedRays(),

                surface1,
                hit1,
                region1,
                normal1,

                surface2,
                hit2,
                region2,
                normal2
             );
        }
        else
        {
            surfaces_.findNearestIntersection
            (
                surfacesToTest,
                start,
                end,

                surface1,
                hit1,
                region1,
                normal1,

                surface2,
                hit2,
                region2,
                normal2
             );
        }
    }

    forAll(testFaces, i)
    {
        label facei = testFaces[i];
        const vector& area = mesh_.faceAreas()[facei];

        if (surface1[i] != -1)
        {
            // If both hit should probably choose 'nearest'
            if
            (
                surface2[i] != -1
             && (
                    magSqr(hit2[i].hitPoint())
                  < magSqr(hit1[i].hitPoint())
                )
            )
            {
                namedSurfaceIndex[facei] = surface2[i];
                posOrientation[facei] = ((area&normal2[i]) > 0);
                nSurfFaces[surface2[i]]++;
            }
            else
            {
                namedSurfaceIndex[facei] = surface1[i];
                posOrientation[facei] = ((area&normal1[i]) > 0);
                nSurfFaces[surface1[i]]++;
            }
        }
        else if (surface2[i] != -1)
        {
            namedSurfaceIndex[facei] = surface2[i];
            posOrientation[facei] = ((area&normal2[i]) > 0);
            nSurfFaces[surface2[i]]++;
        }
    }

    // surfaceIndex might have different surfaces on both sides if
    // there happen to be a (obviously thin) surface with different
    // regions between the cell centres. If one is on a named surface
    // and the other is not this might give problems so sync.
    syncTools::syncFaceList
    (
        mesh_,
        namedSurfaceIndex,
        maxEqOp<label>()
    );

    // Print a bit
    if (debug)
    {
        forAll(nSurfFaces, surfI)
        {
            Pout<< "Surface:"
                << surfaces_.names()[surfI]
                << "  nZoneFaces:" << nSurfFaces[surfI] << nl;
        }
        Pout<< endl;
    }
}


void Foam::meshRefinement::zonify
(
    const label backgroundZoneID,
    const refinementParameters& refineParams,
    const label minZoneRegionSize,

    labelList& testFaces,
    labelList& globalRegion1,
    labelList& globalRegion2,

    labelList& cellToZone,
    labelList& namedSurfaceIndex,
    labelList& namedIntersections,
    PackedBoolList& posOrientation
) const
{
    // Determine zones for cells and faces
    // cellToZone:
    // -2  : unset
    // -1  : not in any zone (zone 'none' or background zone)
    // >=0 : zoneID
    // namedSurfaceIndex, posOrientation:
    // -1  : face not intersected by named surface
    // >=0 : index of named surface
    //       (and posOrientation: surface normal v.s. face normal)

    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const wordList& zonesInMesh = refineParams.zonesInMesh();
    const Switch rezoneSnapProblems = refineParams.rezoneProblemCells();
    const Switch namedLocationsRezone = refineParams.namedLocationsRezone();

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();

    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));
    labelList unnamedSurfaces(surfaceZonesInfo::getUnnamedSurfaces(surfZones));
    labelList moveCentroidSurfaces(surfaces().moveCentroidSurfaces());
    labelList closureSurfaces(surfaces().singleCellClosureSurfaces());
    labelList wrapSurfaces(surfaces().wrapLevelSurfaces());

    // Get map from surface to cellZone (or -1)
    labelList surfaceToCellZone;
    if (namedSurfaces.size())
    {
        // Get/add cellZones corresponding to surface names
        surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
    }

    namedSurfaceIndex.clear();
    posOrientation.clear();

    if (moveCentroidSurfaces.size())
    {
        testFaces = intersectedAndEdgeCellFaces();
    }
    else if (checkForCracks() || closureSurfaces.size())
    {
        testFaces = intersectedAndNeighbouring();
    }
    else
    {
        testFaces = intersectedFaces();
    }

    pointField cellCentres = mesh_.cellCentres();

    // 1. Test all (unnamed & named) surfaces
    {
        // Whether to move cell centroids on surface geometry for
        // intersection checking
        if (moveCentroidSurfaces.size() > 0)
        {
            if (wrapSurfaces.size() > 0)
            {
                WarningInFunction
                    << "Keyword moveCentroidsTol is set to non-zero value "
                    << " and wrapLevel is active. This might cause issues "
                    << "to cell level wrapping."
                    << endl;
            }

            const scalarField& moveCentroidTol = surfaces().moveCentroids();
            boolList cellsToCheck(mesh_.nCells(), false);
            scalar tol = max(moveCentroidTol);
            scalar tolSqr = sqr(tol);

            label nSet = 0;
            forAll(testFaces, tFI)
            {
                label facei = testFaces[tFI];
                label patchI = patches.whichPatch(facei);
                label own = mesh_.faceOwner()[facei];
                if (!cellsToCheck[own])
                {
                    cellsToCheck[own] = true;
                    nSet++;
                }

                if (patchI == -1)
                {
                    label nei = mesh_.faceNeighbour()[facei];
                    if (!cellsToCheck[nei])
                    {
                        cellsToCheck[nei] = true;
                        nSet++;
                    }
                }
            }

            labelList checkCells(nSet, -1);
            pointField checkCCs(nSet, vector::zero);
            nSet = 0;
            forAll(mesh_.cells(), celli)
            {
                if (cellsToCheck[celli])
                {
                    checkCells[nSet] = celli;
                    checkCCs[nSet] = mesh_.cellCentres()[celli];
                    nSet++;
                }
            }

            boolList markedCCs(checkCCs.size(), false);
            forAll(moveCentroidSurfaces, i)
            {
                label surfI = moveCentroidSurfaces[i];
                labelList movedSurface(1, surfI);

                label nMarked = 0;
                forAll(markedCCs, mCCI)
                {
                    if (!markedCCs[mCCI])
                    {
                        nMarked++;
                    }
                }

                pointField notFoundCCs(nMarked);
                labelList localToMaster(nMarked);
                nMarked = 0;
                forAll(markedCCs, mCCI)
                {
                    if (!markedCCs[mCCI])
                    {
                        notFoundCCs[nMarked] = checkCCs[mCCI];
                        localToMaster[nMarked] = mCCI;
                        nMarked++;
                    }
                }

                List<pointIndexHit> nearestInfo;
                labelList nearestSurface;
                surfaces_.findNearest
                (
                    movedSurface,
                    notFoundCCs,
                    scalarField(notFoundCCs.size(), tolSqr),    // sqr of attraction
                    nearestSurface,
                    nearestInfo
                );

                DynamicList<pointIndexHit> localHits;
                DynamicList<label> cellIDs;
                forAll(nearestSurface, i)
                {
                    if (nearestInfo[i].hit())
                    {
                        point hPoint = nearestInfo[i].hitPoint();
                        label celli = checkCells[localToMaster[i]];
                        point cc = mesh_.cellCentres()[celli];
                        point ccToHPt = hPoint - cc;
                        scalar hitDist(magSqr(ccToHPt));

                        if (hitDist < tolSqr)
                        {
                            localHits.append(nearestInfo[i]);
                            cellIDs.append(celli);
                            markedCCs[localToMaster[i]] = true;
                        }
                    }
                }

                label geomI = surfaces_.surfaces()[surfI];

                pointField localNormals;
                surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

                forAll(localHits, localI)
                {
                    point lNormal = localNormals[localI];
                    point hPoint = localHits[localI].hitPoint();

                    label celli = cellIDs[localI];
                    cellCentres[celli] = hPoint + lNormal*tol;
                }
            }
        }

        getIntersections
        (
            cellCentres,
            unnamedSurfaces,  // surfacesToTest,
            testFaces,
            globalRegion1,
            globalRegion2,
            checkForCracks()
        );

        if (wrapSurfaces.size() > 0)
        {
            boolList markedFaces(mesh_.nFaces(), false);
            forAll(testFaces, i)
            {
                label facei = testFaces[i];
                markedFaces[facei] = true;
            }

            DynamicList<label> addedFaces(mesh_.nFaces()/10);
            forAll(wrapIndex_,facei)
            {
                if (wrapIndex_[facei] != -1 && globalRegion1[facei] == -1)
                {
                    globalRegion1[facei] = wrapIndex_[facei];
                    globalRegion2[facei] = wrapIndex_[facei];
                    if (!markedFaces[facei])
                    {
                        addedFaces.append(facei);
                    }
                }
            }

            if (addedFaces.size() > 0)
            {
                label sz = testFaces.size();
                testFaces.setSize(sz+addedFaces.size());
                forAll(addedFaces, i)
                {
                    testFaces[sz+i] = addedFaces[i];
                }
            }

            syncTools::syncFaceList
            (
                mesh_,
                globalRegion1,
                maxEqOp<label>()
            );

            syncTools::syncFaceList
            (
                mesh_,
                globalRegion2,
                maxEqOp<label>()
            );
        }
        else if (closureSurfaces.size() > 0)
        {
            labelList wrapFaces(mesh_.nFaces(), -1);
            const boolList singleCellClosure = surfaces_.singleCellClosure();
            labelList wrapLevel(singleCellClosure.size(), -1);
            forAll(singleCellClosure, regioni)
            {
                if (singleCellClosure[regioni])
                {
                    wrapLevel[regioni] = 0;
                }
            }

            closeSingleCellGaps
            (
                closureSurfaces,
                locationsInMesh,
                globalRegion1,
                wrapLevel,
                wrapFaces
            );

            label nWrapped = 0;
            forAll(wrapFaces, facei)
            {
                if (wrapFaces[facei] != -1 && globalRegion1[facei] == -1)
                {
                    nWrapped++;
                }
            }

            label sz = testFaces.size();
            testFaces.setSize(sz+nWrapped);
            nWrapped = 0;
            forAll(wrapFaces, facei)
            {
                if (wrapFaces[facei] != -1 && globalRegion1[facei] == -1)
                {
                    globalRegion1[facei] = wrapFaces[facei];
                    globalRegion2[facei] = wrapFaces[facei];
                    testFaces[sz+nWrapped] = facei;
                    nWrapped++;
                }
            }

            syncTools::syncFaceList
            (
                mesh_,
                globalRegion1,
                maxEqOp<label>()
            );

            syncTools::syncFaceList
            (
                mesh_,
                globalRegion2,
                maxEqOp<label>()
            );
        }
    }

    //Re-baffle thin gaps cells
    rebaffleThinGapCells
    (
        locationsInMesh,
        cellCentres,
        globalRegion1,
        globalRegion2,
        testFaces
    );

    if (namedSurfaces.size())
    {
        getIntersections
        (
            cellCentres,
            namedSurfaces,
            testFaces,
            namedSurfaceIndex,
            posOrientation
        );
    }
    namedIntersections = namedSurfaceIndex;

    //If added locationsInMesh regionise and filter out duplicates
    refineParams.filterLocations
    (
        mesh_,
        meshCutter_,
        globalRegion1,
        namedSurfaceIndex
    );

    // 2. Walk from locationsInMesh. Hard set cellZones.

    if (locationsInMesh.size())
    {
        Info<< "Setting cellZones according to locationsInMesh:" << endl;

        labelList locationsZoneIDs(zonesInMesh.size(), -1);
        forAll(locationsInMesh, i)
        {
            const word& name = zonesInMesh[i];

            Info<< "Location : " << locationsInMesh[i] << nl
                << "    cellZone : " << name << endl;

            if (name != "none")
            {
                label zoneID = mesh_.cellZones().findZoneID(name);
                if (zoneID == -1)
                {
                    FatalErrorInFunction << "problem" << abort(FatalError);
                }
                locationsZoneIDs[i] = zoneID;
            }
        }
        Info<< endl;

        // Assign cellZone according to seed points
        findCellZoneInsideWalk
        (
            locationsInMesh,    // locations
            locationsZoneIDs,   // index of cellZone (or -1)
            globalRegion1,      // per face -1 (unblocked) or >= 0 (blocked)
            namedSurfaceIndex,

            cellToZone
        );
    }

    // 3. Mark named-surfaces-with-geometric faces. Do geometric test. Soft set
    // cellZones. Correct through making consistent.

    // Closed surfaces with cellZone specified.
    labelList closedNamedSurfaces
    (
        surfaceZonesInfo::getClosedNamedSurfaces
        (
            surfZones,
            surfaces_.geometry(),
            surfaces_.surfaces()
        )
    );

    if (closedNamedSurfaces.size())
    {
        Info<< "Found " << closedNamedSurfaces.size()
            << " closed, named surfaces. Assigning cells in/outside"
            << " these surfaces to the corresponding cellZone."
            << nl << endl;

        findCellZoneClosed
        (
            closedNamedSurfaces,    // indices of closed surfaces
            surfaceToCellZone,      // cell zone index per surface
            namedSurfaceIndex,      // per face index of named surface

            cellToZone,
            minZoneRegionSize
        );
    }

    // 4. Mark named-surfaces-with-insidePoint. Hard set cellZones.

    labelList locationSurfaces
    (
        surfaceZonesInfo::getInsidePointNamedSurfaces(surfZones)
    );

    if (locationSurfaces.size())
    {
        Info<< "Found " << locationSurfaces.size()
            << " named surfaces with a provided inside point."
            << " Assigning cells inside these surfaces"
            << " to the corresponding cellZone."
            << nl << endl;

        // Collect per surface the -insidePoint -cellZone
        pointField insidePoints(locationSurfaces.size());
        labelList insidePointCellZoneIDs(locationSurfaces.size(), -1);
        forAll(locationSurfaces, i)
        {
            label surfI = locationSurfaces[i];
            insidePoints[i] = surfZones[surfI].zoneInsidePoint();

            const word& name = surfZones[surfI].cellZoneName();
            if (name != "none")
            {
                label zoneID = mesh_.cellZones().findZoneID(name);
                if (zoneID == -1)
                {
                    FatalErrorInFunction
                        << "problem"
                        << abort(FatalError);
                }
                insidePointCellZoneIDs[i] = zoneID;
            }
        }

        findCellZoneInsideWalk
        (
            insidePoints,           // locations
            insidePointCellZoneIDs, // index of cellZone
            globalRegion1,          // per face -1 (unblocked) or >= 0 (blocked)
            namedSurfaceIndex,

            cellToZone
        );
    }

    if (debug&MESH)
    {
        Pout<< "Writing cell zone allocation on mesh to time "
            << timeName() << endl;
        mesh_.write();

        volScalarField volCellToZone
        (
            IOobject
            (
                "cellToZoneBeforeTopo",
                timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(cellToZone, celli)
        {
            volCellToZone[celli] = cellToZone[celli];
        }
        volCellToZone.write();
    }

    // 5. Find any unassigned regions (from regionSplit)

    if (namedSurfaces.size())
    {
        Info<< "Walking from known cellZones; crossing a faceZone "
            << "face changes cellZone" << nl << endl;

        findCellZoneTopo
        (
            backgroundZoneID,
            locationsInMesh,
            globalRegion1,       // Intersections with unnamed surfaces
            namedSurfaceIndex,   // Intersections with named surfaces
            surfaceToCellZone,
            cellToZone
        );
    }

    if (namedSurfaces.size())
    {
        makeConsistentFaceIndex
        (
            cellToZone,
            namedSurfaceIndex
        );
    }

    // Merge small regions and make sure namedSurfaceIndex is unset
    // inbetween same cell zones
    if (minZoneRegionSize > 0)
    {
        DynamicList<label> namedLocations(locationsInMesh.size());
        forAll(locationsInMesh, i)
        {
            const word& name = zonesInMesh[i];
            if (name != "none")
            {
                label zoneID = mesh_.cellZones().findZoneID(name);
                if (zoneID != -1)
                {
                    namedLocations.append(zoneID);
                }
            }
        }
        labelHashSet namedLocationsSet(namedLocations);

        excludeSmallRegions
        (
            namedLocationsSet,
            globalRegion1,
            globalRegion2,
            namedSurfaceIndex,
            cellToZone,
            minZoneRegionSize,
            namedLocationsRezone
        );
    }

    if (rezoneSnapProblems && namedSurfaces.size())
    {
        zoneProblemCellReallocate
        (
            globalRegion1,
            globalRegion2,
            namedSurfaceIndex,
            cellToZone
        );
    }

    if (controller_.algorithm() == meshControl::EXTRUDE)
    {
        addOutsideContactCells
        (
            refineParams,
            cellCentres,
            globalRegion1,
            globalRegion2,
            cellToZone,
            testFaces
        );

        reallocateCornerCells
        (
            refineParams,
            globalRegion1,
            globalRegion2,
            cellToZone,
            testFaces
        );

        reallocateSingleRemovalCells
        (
            globalRegion1,
            globalRegion2,
            cellToZone,
            testFaces
        );
    }

    // Some stats
    if (debug)
    {
        label nZones = gMax(cellToZone)+1;

        label nUnvisited = 0;
        label nBackgroundCells = 0;
        labelList nZoneCells(nZones, 0);
        forAll(cellToZone, celli)
        {
            label zoneI = cellToZone[celli];
            if (zoneI >= 0)
            {
                nZoneCells[zoneI]++;
            }
            else if (zoneI == -1)
            {
                nBackgroundCells++;
            }
            else if (zoneI == -2)
            {
                nUnvisited++;
            }
            else
            {
                FatalErrorInFunction
                    << "problem" << exit(FatalError);
            }
        }
        reduce(
            std::tie(nUnvisited, nBackgroundCells),
            UniformParallelOp<sumOp<label>, 2>{}
        );
        Pstream::listCombineReduce(nZoneCells, sumOp<label>{});

        Info<< "nUnvisited      :" << nUnvisited << endl;
        Info<< "nBackgroundCells:" << nBackgroundCells << endl;
        Info<< "nZoneCells      :" << nZoneCells << endl;
    }

    if (debug&MESH)
    {
        Time& runTime = const_cast<Time&>(mesh_.time());
        runTime++;

        Pout<< "Writing cell zone allocation on mesh to time "
            << timeName() << endl;

        faceSet intersectedFaces(mesh_, "intersectedFaces", mesh_.nFaces()/100);

        forAll(globalRegion1, facei)
        {
            if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
            {
                intersectedFaces.insert(facei);
            }
        }
        intersectedFaces.instance() = timeName();
        Pout<< "Dumping " << intersectedFaces.size()
            << " intersected faces to "
            << intersectedFaces.objectPath() << endl;
        intersectedFaces.write();

        mesh_.write();

        volScalarField volCellToZone
        (
            IOobject
            (
                "cellToZone",
                timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(cellToZone, celli)
        {
            volCellToZone[celli] = cellToZone[celli];
        }
        volCellToZone.write();
    }
}


void Foam::meshRefinement::handleSnapProblems
(
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const refinementParameters& refineParams
)
{
    Info<< nl
        << "Introducing baffles to block off problem cells" << nl
        << "----------------------------------------------" << nl
        << endl;

    labelList facePatch;
    if (useTopologicalSnapDetection)
    {
        // Surfaces that need to be baffled
        const labelList surfacesToBaffle
        (
            surfaceZonesInfo::getUnnamedSurfaces(surfaces_.surfZones())
        );

        facePatch = markFacesOnProblemCells
        (
            motionDict,
            refineParams,
            removeEdgeConnectedCells,
            perpendicularAngle,
            globalToMasterPatch,
            meshedPatches(),
            surfacesToBaffle
        );
    }
    else
    {
        facePatch = markFacesOnProblemCellsGeometric
        (
            snapParams,
            motionDict,
            globalToMasterPatch,
            globalToSlavePatch
        );
    }
    Info<< "Analyzed problem cells in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    if (debug&MESH)
    {
        faceSet problemFaces(mesh_, "problemFaces", mesh_.nFaces()/100);

        forAll(facePatch, facei)
        {
            if (facePatch[facei] != -1)
            {
                problemFaces.insert(facei);
            }
        }
        problemFaces.instance() = timeName();
        Pout<< "Dumping " << problemFaces.size()
            << " problem faces to " << problemFaces.objectPath() << endl;
        problemFaces.write();
    }

    Info<< "Introducing baffles to delete problem cells." << nl << endl;

    if (debug)
    {
        runTime++;
    }

    // Create baffles with same owner and neighbour for now.
    createBaffles(facePatch, facePatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles");

    if (debug&MESH)
    {
        Pout<< "Writing extra baffled mesh to time "
            << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/"extraBaffles"
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


void Foam::meshRefinement::handleZoneSnapProblems
(
    const snapParameters& snapParams,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const refinementParameters& refineParams,
    List<labelPair>& baffles
)
{
    Info<< nl
        << "Introducing baffles to block off problem zone cells" << nl
        << "---------------------------------------------------" << nl
        << endl;

    label defaultPatch = 0;
    if (meshedPatches().size())
    {
        defaultPatch = meshedPatches()[0];
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const label minZoneRegionSize = refineParams.minZoneRegionSize();
    bool firstPass = true;
    while (true)
    {
        removeCells cellRemover(mesh_);

        labelList facePatch
        (
            markFacesOnProblemCells
            (
                motionDict,
                refineParams,
                false,  // perpendicular edge connected cells
                scalarField(0), // per region perpendicular angle
                globalToMasterPatch,
                meshedPatches(),
                identity(surfaces().surfaces().size())
             )
         );

        DynamicList<label> cellsToRemove(mesh_.nCells());
        forAll(mesh_.cells(), celli)
        {
            const cell c = mesh_.cells()[celli];
            bool leak = false;
            label patchI =  -1;

            forAll(c, cfI)
            {
                label facei = c[cfI];

                label patchFaceI = patches.whichPatch(facei);
                if (patchFaceI == -1)
                {
                    if (facePatch[facei] == -1)
                    {
                        leak = true;
                    }
                }
                else
                {
                    if (!patches[patchFaceI].coupled())
                    {
                        patchI = patchFaceI;
                    }
                    else
                    {
                        if (facePatch[facei] == -1)
                        {
                            leak = true;
                        }
                    }
                }
            }
            if (!leak)
            {
                cellsToRemove.append(celli);
                forAll(c, cfI)
                {
                    label facei = c[cfI];

                    label patchFaceI = patches.whichPatch(facei);
                    if (patchFaceI == -1)
                    {
                        facePatch[facei] = patchI;
                    }
                    else if
                    (
                        patchFaceI != -1
                        && patches[patchFaceI].coupled()
                    )
                    {
                        facePatch[facei] = patchI;
                    }
                }
            }
        }
        cellsToRemove.shrink();

        if (returnReduce(cellsToRemove.size(), sumOp<label>()) != 0)
        {
            // Pick up patches for exposed faces
            labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
            labelList exposedPatches(exposedFaces.size());

            forAll(exposedFaces, i)
            {
                label facei = exposedFaces[i];

                if (facePatch[facei] != -1)
                {
                    exposedPatches[i] = facePatch[facei];
                }
                else
                {
                    exposedPatches[i] = defaultPatch;
                }
            }

            autoPtr<mapPolyMesh> map = doRemoveCells
            (
                cellsToRemove,
                exposedFaces,
                exposedPatches,
                cellRemover
             );
            updateBaffles(map, baffles);
        }
        else if (!firstPass)
        {
            break;
        }

        if (firstPass || returnReduce(cellsToRemove.size(), sumOp<label>()) != 0)
        {
            firstPass = false;
            DynamicList<label> regionCellsToRemove(mesh_.nCells());

            regionSplit cellRegion(mesh_);
            boolList keepRegions(cellRegion.nRegions(), false);
            labelList regionSize(cellRegion.nRegions(),0);
            forAll(cellRegion, i)
            {
                label regionI = cellRegion[i];
                regionSize[regionI]++;
            }
            Pstream::listCombineReduce(regionSize, plusOp<label>());

            forAll(locationsInMesh, i)
            {
                // Get location and index of zone ("none" for cellZone -1)
                const point& insidePoint = locationsInMesh[i];

                // Find the region containing the keepPoint
                label keepRegionI = -1;

                label celli = findCell
                (
                    insidePoint,
                    mesh_,
                    meshCutter_
                );

                if (celli != -1)
                {
                    keepRegionI = cellRegion[celli];
                }
                reduce(keepRegionI, maxOp<label>());
/*
                // Find the region containing the insidePoint
                label keepRegionI = findRegion
                (
                    mesh_,
                    cellRegion,
                    mergeDistance_*vector(1,1,1),
                    insidePoint
                );
*/
                keepRegions[keepRegionI] = true;
            }

            while (true)
            {
                label nUpdated = 0;
                forAll(baffles, baffleI)
                {
                    label face0 = baffles[baffleI].first();
                    label face1 = baffles[baffleI].second();
                    if (face0 < 0 || face1 < 0)
                    {
                        continue;
                    }

                    label own0 = mesh_.faceOwner()[face0];
                    label own1 = mesh_.faceOwner()[face1];

                    label region0 = cellRegion[own0];
                    label region1 = cellRegion[own1];

                    if (keepRegions[region0] || keepRegions[region1])
                    {
                        if (region0 != region1)
                        {
                            if
                            (
                                !keepRegions[region0]
                                && regionSize[region0] > minZoneRegionSize
                            )
                            {
                                keepRegions[region0] = true;
                                nUpdated++;
                            }
                            if
                            (
                                !keepRegions[region1]
                                && regionSize[region1] > minZoneRegionSize
                            )
                            {
                                keepRegions[region1] = true;
                                nUpdated++;
                            }
                        }
                    }
                }
                if (returnReduce(nUpdated, sumOp<label>()) == 0)
                {
                    break;
                }
                else
                {
                    Pstream::listCombineReduce(keepRegions, orOp<bool>());
                }
            }

            forAll(mesh_.cells(), celli)
            {
                label regionI = cellRegion[celli];
                if (!keepRegions[regionI])
                {
                    regionCellsToRemove.append(celli);
                }
            }
            regionCellsToRemove.shrink();

            if (returnReduce(regionCellsToRemove.size(), sumOp<label>()) != 0)
            {
                Info<<"Removing "
                    << returnReduce(regionCellsToRemove.size(), sumOp<label>())
                    <<" cells in un-connected regions" <<endl;

                facePatch = -1;
                forAll(regionCellsToRemove, i)
                {
                    label celli = regionCellsToRemove[i];
                    const cell c = mesh_.cells()[celli];

                    label patchI = -1;
                    forAll(c, cfI)
                    {
                        label facei = c[cfI];
                        label patchFaceI = patches.whichPatch(facei);

                        if
                        (
                            patchFaceI != -1
                            && !patches[patchFaceI].coupled()
                        )
                        {
                            patchI = patchFaceI;
                            break;
                        }
                    }

                    forAll(c, cfI)
                    {
                        label facei = c[cfI];
                        facePatch[facei] = patchI;
                    }
                }

                // Pick up patches for exposed faces
                labelList exposedFaces
                    (cellRemover.getExposedFaces(regionCellsToRemove));
                labelList exposedPatches(exposedFaces.size());

                forAll(exposedFaces, i)
                {
                    label facei = exposedFaces[i];

                    if (facePatch[facei] != -1)
                    {
                        exposedPatches[i] = facePatch[facei];
                    }
                    else
                    {
                        exposedPatches[i] = defaultPatch;
                    }
                }

                autoPtr<mapPolyMesh> map = doRemoveCells
                (
                    regionCellsToRemove,
                    exposedFaces,
                    exposedPatches,
                    cellRemover
                );
                updateBaffles(map, baffles);
            }
            else
            {
                break;
            }
        }
    }
}


Foam::labelList Foam::meshRefinement::freeStandingBaffleFaces
(
    const labelList& faceToZone,
    const labelList& cellToZone,
    const labelList& neiCellZone
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();


    // We want to pick up the faces to orient. These faces come in
    // two variants:
    // - faces originating from stand-alone faceZones
    //   (these will most likely have no cellZone on either side so
    //    ownZone and neiZone both -1)
    // - sticky-up faces originating from a 'bulge' in a outside of
    //   a cellZone. These will have the same cellZone on either side.
    //   How to orient these is not really clearly defined so do them
    //   same as stand-alone faceZone faces for now. (Normally these will
    //   already have been removed by the 'allowFreeStandingZoneFaces=false'
    //   default setting)

    // Note that argument neiCellZone will have -1 on uncoupled boundaries.

    DynamicList<label> faceLabels(mesh_.nFaces()/100);

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (faceToZone[facei] != -1)
        {
            // Free standing baffle?
            label ownZone = cellToZone[faceOwner[facei]];
            label neiZone = cellToZone[faceNeighbour[facei]];
            if (ownZone == neiZone)
            {
                faceLabels.append(facei);
            }
        }
    }
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        forAll(pp, i)
        {
            label facei = pp.start()+i;
            if (faceToZone[facei] != -1)
            {
                // Free standing baffle?
                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                if (ownZone == neiZone)
                {
                    faceLabels.append(facei);
                }
            }
        }
    }
    return faceLabels.shrink();
}


void Foam::meshRefinement::calcPatchNumMasterFaces
(
    const PackedBoolList& isMasterFace,
    const indirectPrimitivePatch& patch,
    labelList& nMasterFacesPerEdge
) const
{
    // Number of (master)faces per edge
    nMasterFacesPerEdge.setSize(patch.nEdges());
    nMasterFacesPerEdge = 0;

    forAll(patch.addressing(), facei)
    {
        const label meshFaceI = patch.addressing()[facei];

        if (isMasterFace[meshFaceI])
        {
            const labelList& fEdges = patch.faceEdges()[facei];
            forAll(fEdges, fEdgeI)
            {
                nMasterFacesPerEdge[fEdges[fEdgeI]]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        patch.meshEdges(mesh_.edges(), mesh_.pointEdges()),
        nMasterFacesPerEdge,
        plusEqOp<label>(),
        label(0)
    );
}


Foam::label Foam::meshRefinement::markPatchZones
(
    const indirectPrimitivePatch& patch,
    const labelList& nMasterFacesPerEdge,
    labelList& faceToZone
) const
{
    List<patchEdgeFaceRegion> allEdgeInfo(patch.nEdges());
    List<patchEdgeFaceRegion> allFaceInfo(patch.size());


    // Protect all non-manifold edges
    {
//        label nProtected = 0;

        forAll(nMasterFacesPerEdge, edgeI)
        {
            if (nMasterFacesPerEdge[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = -2;
//                nProtected++;
            }
        }
        //Info<< "Protected from visiting "
        //    << returnReduce(nProtected, sumOp<label>())
        //    << " non-manifold edges" << nl << endl;
    }


    // Hand out zones

    DynamicList<label> changedEdges;
    DynamicList<patchEdgeFaceRegion> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchEdgeFaceRegion
    >::propagationTol();

    int dummyTrackData;

    const globalIndex globalFaces(patch.size());

    label facei = 0;

    label currentZoneI = 0;

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        for (; facei < allFaceInfo.size(); facei++)
        {
            if (!allFaceInfo[facei].valid(dummyTrackData))
            {
                globalSeed = globalFaces.toGlobal(facei);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        label procI = globalFaces.whichProcID(globalSeed);
        label seedFaceI = globalFaces.toLocal(procI, globalSeed);

        //Info<< "Seeding zone " << currentZoneI
        //    << " from processor " << procI << " face " << seedFaceI
        //    << endl;
        label nUpdatedEdges = changedEdges.size();
        if (procI == Pstream::myProcNo())
        {
            patchEdgeFaceRegion& faceInfo = allFaceInfo[seedFaceI];


            // Set face
            faceInfo = currentZoneI;

            // .. and seed its edges
            const labelList& fEdges = patch.faceEdges()[seedFaceI];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchEdgeFaceRegion& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFaceI,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }

        nUpdatedEdges = changedEdges.size() - nUpdatedEdges;
        if (returnReduce(nUpdatedEdges , sumOp<label>()) != 0)
        {
            // Walk
           PatchEdgeFaceWave
           <
               indirectPrimitivePatch,
               patchEdgeFaceRegion
           > calc
           (
               mesh_,
               patch,
               changedEdges,
               changedInfo,
               allEdgeInfo,
               allFaceInfo,
               returnReduce(patch.nEdges(), sumOp<label>())
           );
        }
        currentZoneI++;
    }


    faceToZone.setSize(patch.size());
    forAll(allFaceInfo, facei)
    {
        if (!allFaceInfo[facei].valid(dummyTrackData))
        {
            FatalErrorInFunction
                << "Problem: unvisited face " << facei
                << " at " << patch.faceCentres()[facei]
                << exit(FatalError);
        }
        faceToZone[facei] = allFaceInfo[facei].region();
    }

    return currentZoneI;
}


void Foam::meshRefinement::consistentOrientation
(
    const PackedBoolList& isMasterFace,
    const indirectPrimitivePatch& patch,
    const labelList& nMasterFacesPerEdge,
    const labelList& faceToZone,
    const Map<label>& zoneToOrientation,
    PackedBoolList& meshFlipMap
) const
{
    const polyBoundaryMesh& bm = mesh_.boundaryMesh();

    // Data on all edges and faces
    List<patchFaceOrientation> allEdgeInfo(patch.nEdges());
    List<patchFaceOrientation> allFaceInfo(patch.size());

    // Make sure we don't walk through
    // - slaves of coupled faces
    // - non-manifold edges
    {
//        label nProtected = 0;

        forAll(patch.addressing(), facei)
        {
            const label meshFaceI = patch.addressing()[facei];
            const label patchI = bm.whichPatch(meshFaceI);

            if
            (
                patchI != -1
             && bm[patchI].coupled()
             && !isMasterFace[meshFaceI]
            )
            {
                // Slave side. Mark so doesn't get visited.
                allFaceInfo[facei] = orientedSurface::NOFLIP;
//                nProtected++;
            }
        }
        //Info<< "Protected from visiting "
        //    << returnReduce(nProtected, sumOp<label>())
        //    << " slaves of coupled faces" << nl << endl;
    }
    {
        label nProtected = 0;

        forAll(nMasterFacesPerEdge, edgeI)
        {
            if (nMasterFacesPerEdge[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }

        Info<< "Protected from visiting "
            << returnReduce(nProtected, sumOp<label>())
            << " non-manifold edges" << nl << endl;
    }



    DynamicList<label> changedEdges;
    DynamicList<patchFaceOrientation> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchFaceOrientation
    >::propagationTol();

    int dummyTrackData;

    globalIndex globalFaces(patch.size());

    while (true)
    {
        // Pick an unset face
        label globalSeed = labelMax;
        forAll(allFaceInfo, facei)
        {
            if (allFaceInfo[facei] == orientedSurface::UNVISITED)
            {
                globalSeed = globalFaces.toGlobal(facei);
                break;
            }
        }

        reduce(globalSeed, minOp<label>());

        if (globalSeed == labelMax)
        {
            break;
        }

        label procI = globalFaces.whichProcID(globalSeed);
        label seedFaceI = globalFaces.toLocal(procI, globalSeed);

        //Info<< "Seeding from processor " << procI << " face " << seedFaceI
        //    << endl;

        label nUpdatedEdges = changedEdges.size();

        if (procI == Pstream::myProcNo())
        {
            // Determine orientation of seedFace

            patchFaceOrientation& faceInfo = allFaceInfo[seedFaceI];

            // Start off with correct orientation
            faceInfo = orientedSurface::NOFLIP;

            if (zoneToOrientation[faceToZone[seedFaceI]] < 0)
            {
                faceInfo.flip();
            }


            const labelList& fEdges = patch.faceEdges()[seedFaceI];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchFaceOrientation& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFaceI,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }

        nUpdatedEdges = changedEdges.size() - nUpdatedEdges;
        if (returnReduce(nUpdatedEdges , sumOp<label>()) != 0)
        {
            // Walk
            PatchEdgeFaceWave
            <
                indirectPrimitivePatch,
                patchFaceOrientation
            > calc
            (
                mesh_,
                patch,
                changedEdges,
                changedInfo,
                allEdgeInfo,
                allFaceInfo,
                returnReduce(patch.nEdges(), sumOp<label>())
            );
        }
    }


    // Push master zone info over to slave (since slave faces never visited)
    {
        labelList neiStatus
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            orientedSurface::UNVISITED
        );

        forAll(patch.addressing(), i)
        {
            const label meshFaceI = patch.addressing()[i];
            if (!mesh_.isInternalFace(meshFaceI))
            {
                neiStatus[meshFaceI-mesh_.nInternalFaces()] =
                    allFaceInfo[i].flipStatus();
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiStatus);

        forAll(patch.addressing(), i)
        {
            const label meshFaceI = patch.addressing()[i];
            const label patchI = bm.whichPatch(meshFaceI);

            if
            (
                patchI != -1
             && bm[patchI].coupled()
             && !isMasterFace[meshFaceI]
            )
            {
                // Slave side. Take flipped from neighbour
                label bFaceI = meshFaceI-mesh_.nInternalFaces();

                if (neiStatus[bFaceI] == orientedSurface::NOFLIP)
                {
                    allFaceInfo[i] = orientedSurface::FLIP;
                }
                else if (neiStatus[bFaceI] == orientedSurface::FLIP)
                {
                    allFaceInfo[i] = orientedSurface::NOFLIP;
                }
                else
                {
                    FatalErrorInFunction
                        << "Incorrect status for face " << meshFaceI
                        << abort(FatalError);
                }
            }
        }
    }


    // Convert to meshFlipMap and adapt faceZones

    meshFlipMap.setSize(mesh_.nFaces());
    meshFlipMap = false;

    forAll(allFaceInfo, facei)
    {
        label meshFaceI = patch.addressing()[facei];

        if (allFaceInfo[facei] == orientedSurface::NOFLIP)
        {
            meshFlipMap[meshFaceI] = false;
        }
        else if (allFaceInfo[facei] == orientedSurface::FLIP)
        {
            meshFlipMap[meshFaceI] = true;
        }
        else
        {
            FatalErrorInFunction
                << "Problem : unvisited face " << facei
                << " centre:" << mesh_.faceCentres()[meshFaceI]
                << abort(FatalError);
        }
    }
}


void Foam::meshRefinement::zonify
(
    // Get per face whether is it master (of a coupled set of faces)
    const PackedBoolList& isMasterFace,
    const labelList& cellToZone,
    const labelList& neiCellZone,
    const labelList& faceToZone,
    const PackedBoolList& meshFlipMap,
    polyTopoChange& meshMod
) const
{
    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        label faceZoneI = faceToZone[facei];

        if (faceZoneI != -1)
        {
            // Orient face zone to have slave cells in min cell zone.
            // Note: logic to use flipMap should be consistent with logic
            //       to pick up the freeStandingBaffleFaces!

            label ownZone = cellToZone[faceOwner[facei]];
            label neiZone = cellToZone[faceNeighbour[facei]];

            bool flip;

            if (ownZone == neiZone)
            {
                // free-standing face. Use geometrically derived orientation
                flip = meshFlipMap[facei];
            }
            else
            {
                flip =
                (
                    ownZone == -1
                 || (neiZone != -1 && ownZone > neiZone)
                );
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[facei],           // modified face
                    facei,                          // label of face
                    faceOwner[facei],               // owner
                    faceNeighbour[facei],           // neighbour
                    false,                          // face flip
                    -1,                             // patch for face
                    false,                          // remove from zone
                    faceZoneI,                      // zone for face
                    flip                            // face flip in zone
                )
            );
        }
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Set owner as no-flip
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        label facei = pp.start();

        forAll(pp, i)
        {
            label faceZoneI = faceToZone[facei];

            if (faceZoneI != -1)
            {
                label ownZone = cellToZone[faceOwner[facei]];
                label neiZone = neiCellZone[facei-mesh_.nInternalFaces()];

                bool flip;

                if (ownZone == neiZone)
                {
                    // free-standing face. Use geometrically derived orientation
                    flip = meshFlipMap[facei];
                }
                else
                {
                    flip =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                }

                meshMod.setAction
                (
                    polyModifyFace
                    (
                        mesh_.faces()[facei],           // modified face
                        facei,                          // label of face
                        faceOwner[facei],               // owner
                        -1,                             // neighbour
                        false,                          // face flip
                        patchI,                         // patch for face
                        false,                          // remove from zone
                        faceZoneI,                      // zone for face
                        flip                            // face flip in zone
                    )
                );
            }
            facei++;
        }
    }


    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(cellToZone, celli)
    {
        label zoneI = cellToZone[celli];

        if (zoneI >= 0)
        {
            meshMod.setAction
            (
                polyModifyCell
                (
                    celli,
                    false,          // removeFromZone
                    zoneI
                )
            );
        }
    }
}


void Foam::meshRefinement::resetFaceZoneFlipMap()
{
    labelList cellToZone(mesh_.nCells(), -1);

    const cellZoneMesh& cellZones = mesh_.cellZones();

    forAll(cellZones, cellZoneI)
    {
        const cellZone& cZone = mesh_.cellZones()[cellZoneI];
        forAll(cZone, czi)
        {
            label celli = cZone[czi];
            cellToZone[celli] = cellZoneI;
        }
    }

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

    const faceZoneMesh& faceZones = mesh_.faceZones();
    forAll(faceZones, faceZoneI)
    {
        faceZone& fZone = mesh_.faceZones()[faceZoneI];
        boolList zoneFlipMap = fZone.flipMap();
        labelUList addressing = fZone;

        forAll(fZone, fzi)
        {
            label facei = fZone[fzi];
            label own = mesh_.faceOwner()[facei];
            label ownZone = cellToZone[own];

            label neiZone = -1;
            if (facei < mesh_.nInternalFaces())
            {
                label nei = mesh_.faceNeighbour()[facei];
                neiZone = cellToZone[nei];
            }
            else
            {
                neiZone = neiCellZone[facei - mesh_.nInternalFaces()];
            }

            if (ownZone != neiZone)
            {
                bool flip =
                (
                    ownZone == -1
                 || (neiZone != -1 && ownZone > neiZone)
                );

                zoneFlipMap[fzi] = flip;
            }
        }

        fZone.resetAddressing(addressing, zoneFlipMap);
    }

    return;
}


void Foam::meshRefinement::allocateInterRegionFaceZone
(
    const label ownZone,
    const label neiZone,
    wordPairHashTable& zonesToFaceZone,
    HashTable<word, labelPair, labelPair::Hash<>>& zoneIDsToFaceZone
) const
{
    const cellZoneMesh& cellZones = mesh_.cellZones();

    if (ownZone != neiZone)
    {
        // Make sure lowest number cellZone is master. Non-cellZone
        // areas are slave
        bool swap =
        (
            ownZone == -1
         || (neiZone != -1 && ownZone > neiZone)
        );

        // Quick check whether we already have pair of zones
        labelPair key(ownZone, neiZone);
        if (swap)
        {
            Swap(key.first(), key.second());
        }

        HashTable<word, labelPair, labelPair::Hash<>>::
        const_iterator zoneFnd = zoneIDsToFaceZone.find
        (
            key
        );

        if (zoneFnd == zoneIDsToFaceZone.end())
        {
            // Not found. Allocate.
            const word ownZoneName =
            (
                ownZone != -1
              ? cellZones[ownZone].name()
              : "none"
            );
            const word neiZoneName =
            (
                neiZone != -1
              ? cellZones[neiZone].name()
              : "none"
            );

            // Get lowest zone first
            Pair<word> wordKey(ownZoneName, neiZoneName);
            if (swap)
            {
                Swap(wordKey.first(), wordKey.second());
            }

            word fzName = wordKey.first() + "_to_" + wordKey.second();

            zoneIDsToFaceZone.insert(key, fzName);
            zonesToFaceZone.insert(wordKey, fzName);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshRefinement::baffleAndSplitMesh
(
    const bool doHandleSnapProblems,
    const snapParameters& snapParams,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const refinementParameters& refineParams,
    const writer<scalar>& leakPathFormatter,
    const bool threaded /*= false*/
)
{
    const bool useTopologicalSnapDetection =
        refineParams.useTopologicalSnapDetection();
    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const pointField& locationsOutsideMesh =
        refineParams.locationsOutsideMesh();
    const bool extrudeExtraRemoval = refineParams.extrudeExtraRemoval();
    const bool split = refineParams.splitCells();

    // Introduce baffles
    // ~~~~~~~~~~~~~~~~~

    // Split the mesh along internal faces wherever there is a pierce between
    // two cell centres.

    Info<< "Introducing baffles for "
        << returnReduce(countHits(), sumOp<label>())
        << " faces that are intersected by the surface." << nl << endl;

    if
    (
        split && !removeEdgeConnectedCells
        && controller_.algorithm() == meshControl::STANDARD
        && controller_.mode() != meshControl::DRYRUN
    )
    {
        Info<<"Snapping and cutting cells"<<endl;
        setOldPoints(mesh_.points());
        snapAndCut();

        Info<< "Completed cutting of mesh cells in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }

    labelList ownPatch, nbrPatch;
    getBafflePatches
    (
        refineParams,
        globalToMasterPatch,
        ownPatch,
        nbrPatch
    );

    createBaffles(ownPatch, nbrPatch);

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }

    Info<< "Created baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After introducing baffles");

    if (debug&MESH)
    {
        Pout<< "Writing baffled mesh to time " << timeName()
            << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/"baffles"
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }


    // Introduce baffles to delete problem cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Create some additional baffles where we want surface cells removed.

    if
    (
        doHandleSnapProblems
        &&
        (
            controller_.algorithm() == meshControl::STANDARD
            || extrudeExtraRemoval
        )
    )
    {
        handleSnapProblems
        (
            snapParams,
            useTopologicalSnapDetection,
            removeEdgeConnectedCells,
            perpendicularAngle,
            motionDict,
            runTime,
            globalToMasterPatch,
            globalToSlavePatch,
            refineParams
        );

        // Removing additional cells might have created disconnected bits
        // so re-do the surface intersections
        {
            labelList ownPatch, nbrPatch;
            getBafflePatches
            (
                refineParams,
                globalToMasterPatch,

                ownPatch,
                nbrPatch
            );

            createBaffles(ownPatch, nbrPatch);
        }

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            checkData();
        }
    }


    // Select part of mesh
    // ~~~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Remove unreachable sections of mesh" << nl
        << "-----------------------------------" << nl
        << endl;

    if (debug)
    {
        runTime++;
    }

    splitMeshRegions
    (
        globalToMasterPatch,
        globalToSlavePatch,
        locationsInMesh,
        locationsOutsideMesh,
        leakPathFormatter,
        refineParams.fullLeakChecks()
    );

    if (debug)
    {
        // Debug:test all is still synced across proc patches
        checkData();
    }
    Info<< "Split mesh in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;

    printMeshInfo(debug, "After subsetting");

    if (debug&MESH)
    {
        runTime++;
        Pout<< "Writing subsetted mesh to time " << timeName()
            << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            runTime.path()/timeName()
        );
        Pout<< "Dumped debug data in = "
            << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


void Foam::meshRefinement::mergeFreeStandingBaffles
(
    const snapParameters& snapParams,
    const bool useTopologicalSnapDetection,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const scalar planarAngle,
    const dictionary& motionDict,
    Time& runTime,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const refinementParameters& refineParams,
    const writer<scalar>& leakPathFormatter,
    const bool threaded /*= false*/
)
{
    // Merge baffles
    // ~~~~~~~~~~~~~

    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const pointField& locationsOutsideMesh = refineParams.locationsOutsideMesh();

    Info<< nl
        << "Merge free-standing baffles" << nl
        << "---------------------------" << nl
        << endl;


    // List of pairs of freestanding baffle faces.
    List<labelPair> couples
    (
        freeStandingBaffles    // filter out freestanding baffles
        (
            localPointRegion::findDuplicateFacePairs(mesh_),
            planarAngle
        )
    );

    label nCouples = couples.size();
    reduce(nCouples, sumOp<label>());

    Info<< "Detected free-standing baffles : " << nCouples << endl;

    if (nCouples > 0)
    {
        // Actually merge baffles. Note: not exactly parallellized. Should
        // convert baffle faces into processor faces if they resulted
        // from them.
        mergeBaffles(couples, Map<label>(0), true);

        if
        (
            controller_.algorithm() == meshControl::STANDARD
            || controller_.algorithm() == meshControl::SHELL
        )
        {
            // Detect any problem cells resulting from merging of baffles
            // and delete them
            handleSnapProblems
            (
                snapParams,
                useTopologicalSnapDetection,
                removeEdgeConnectedCells,
                perpendicularAngle,
                motionDict,
                runTime,
                globalToMasterPatch,
                globalToSlavePatch,
                refineParams
             );
        }

        // Very occasionally removing a problem cell might create a disconnected
        // region so re-check

        Info<< nl
            << "Remove unreachable sections of mesh" << nl
            << "-----------------------------------" << nl
            << endl;

        if (debug)
        {
            runTime++;
        }

        splitMeshRegions
        (
            globalToMasterPatch,
            globalToSlavePatch,
            locationsInMesh,
            locationsOutsideMesh,
            leakPathFormatter,
            refineParams.fullLeakChecks()
        );


        if (debug)
        {
            // Debug:test all is still synced across proc patches
            checkData();
        }
    }
    Info<< "Merged free-standing baffles in = "
        << runTime.cpuTimeIncrement() << " s\n" << nl << endl;
}


void Foam::meshRefinement::splitMesh
(
    const refinementParameters& refineParams
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get per face whether is it master (of a coupled set of faces)
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));
    const pointField& locationsOutsideMesh =
        refineParams.locationsOutsideMesh();

    labelList unnamedSurfaces
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfaces_.surfZones())
    );

    labelList globalRegion1;
    labelList globalRegion2;
    labelList testFaces;
    if (checkForCracks())
    {
        testFaces = intersectedAndEdgeCellFaces();
    }
    else
    {
        testFaces = intersectedFaces();
    }

    getIntersections
    (
        mesh_.cellCentres(),
        unnamedSurfaces,  // surfacesToTest,
        testFaces,
        globalRegion1,
        globalRegion2,
        checkForCracks()
    );

    labelList bdyPts(mesh_.nPoints(), -1);
    forAll(mesh_.faces(), facei)
    {
        if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
        {
            face f = mesh_.faces()[facei];
            forAll(f,fp)
            {
                label pointi = f[fp];

                bdyPts[pointi] = max(bdyPts[pointi],globalRegion1[facei]);
            }
        }
    }
    syncTools::syncPointList
    (
        mesh_,
        bdyPts,
        maxEqOp<label>(),
        label(-1)
    );

    labelList expandedBdyFaces(mesh_.nFaces(), -1);
    forAll(mesh_.points(), pointi)
    {
        label mPt = bdyPts[pointi];
        if (mPt != -1)
        {
            const labelList& pCells = mesh_.pointCells()[pointi];
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                const cell& c = mesh_.cells()[celli];
                forAll(c, cFI)
                {
                    label facei = c[cFI];
                    expandedBdyFaces[facei] = max(expandedBdyFaces[facei], mPt);
                }
            }
        }
    }
    syncTools::syncFaceList(mesh_, expandedBdyFaces, maxEqOp<label>());

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces(), false);
    forAll(mesh_.faces(), facei)
    {
        if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());

    regionSplit cellRegion(mesh_, blockedFace);
    DynamicList<label> regionsToRemove(locationsOutsideMesh.size());

    labelList rCells(locationsOutsideMesh.size(), -1);

    // For all locationsInMesh find the cell
    forAll(locationsOutsideMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& outsidePoint = locationsOutsideMesh[i];

        // Find the region containing the keepPoint
        label removeRegionI = -1;

        rCells[i] = findCell
        (
            outsidePoint,
            mesh_,
            meshCutter_
        );

        if (rCells[i] != -1)
        {
            removeRegionI = cellRegion[rCells[i]];
        }
        reduce(removeRegionI, maxOp<label>());
        regionsToRemove.append(removeRegionI);
    }
    regionsToRemove.shrink();

    boolList regionIsRemoved(cellRegion.nRegions(), false);
    forAll(regionsToRemove, keepI)
    {
        regionIsRemoved[regionsToRemove[keepI]] = true;
    }

    Pstream::listCombineReduce(regionIsRemoved, orOp<bool>());

    boolList forRemoval(mesh_.nCells(), false);

    forAll(cellRegion, celli)
    {
        label regionI = cellRegion[celli];
        if (regionIsRemoved[regionI])
        {
            forRemoval[celli] = true;
        }
    }

    // Get coupled neighbour cellRegion
    boolList neiForRemoval;
    syncTools::swapBoundaryCellList(mesh_, forRemoval, neiForRemoval);
    forAll(mesh_.faces(), facei)
    {
        if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
        {
            label patchI = patches.whichPatch(facei);
            label own = mesh_.faceOwner()[facei];

            bool removeOwn = forRemoval[own];

            if (removeOwn)
            {
                bool removedNei = false;
                label nei = -1;
                if (patchI == -1)
                {
                    nei = mesh_.faceNeighbour()[facei];
                    removedNei = forRemoval[nei];
                }
                else if (patches[patchI].coupled())
                {
                    label bfacei = facei-mesh_.nInternalFaces();
                    removedNei = neiForRemoval[bfacei];
                }

                if (removedNei)
                {
                    forRemoval[own] = false;
                    if (patchI == -1)
                    {
                        forRemoval[nei] = false;
                    }
                }
            }
        }
    }

    syncTools::swapBoundaryCellList(mesh_, forRemoval, neiForRemoval);
    boolList externalFaces(mesh_.nFaces(), false);

    forAll(mesh_.faces(), facei)
    {
        label patchI = patches.whichPatch(facei);
        label own = mesh_.faceOwner()[facei];

        bool removeOwn = forRemoval[own];

        bool removeNei = true;

        if (patchI == -1)
        {
            label nei = mesh_.faceNeighbour()[facei];
            removeNei = forRemoval[nei];
        }
        else if (patches[patchI].coupled())
        {
            label bfacei = facei-mesh_.nInternalFaces();
            removeNei = neiForRemoval[bfacei];
        }

        if
        (
            (removeOwn && !removeNei)
            || (removeNei && !removeOwn)
        )
        {
            externalFaces[facei] = true;
        }
    }
    syncTools::syncFaceList(mesh_, externalFaces, orEqOp<bool>());


    //Check for pinched edges and points
    boolList pinchedEdges(mesh_.nEdges(), false);
    boolList pinchedPoints(mesh_.nPoints(), false);

    labelList nInterfaceEdgeFaces(mesh_.nEdges(), 0);
    forAll(mesh_.edges(), edgei)
    {
        const labelList& eFaces = mesh_.edgeFaces()[edgei];
        forAll(eFaces, eFI)
        {
            label facei = eFaces[eFI];
            if (isMasterFace[facei] && externalFaces[facei])
            {
                nInterfaceEdgeFaces[edgei]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        nInterfaceEdgeFaces,
        plusEqOp<label>(),
        label(0)
    );

    boolList  pinchedEdgePts(mesh_.nPoints(), false);
    forAll(mesh_.edges(), edgei)
    {
        if (nInterfaceEdgeFaces[edgei] == 4)
        {
            const edge e = mesh_.edges()[edgei];
            pinchedEdgePts[e[0]] = true;
            pinchedEdgePts[e[1]] = true;
            const labelList& eCells = mesh_.edgeCells()[edgei];
            forAll(eCells, eCI)
            {
                forRemoval[eCells[eCI]] = false;
            }
        }
    }


    syncTools::syncPointList
    (
        mesh_,
        pinchedEdgePts,
        orEqOp<bool>(),
        false
    );

    labelList nInterfacePointFaces(mesh_.nPoints(), 0);
    forAll(mesh_.points(), pointi)
    {
        const labelList& pFaces = mesh_.pointFaces()[pointi];
        forAll(pFaces, pFI)
        {
            label facei = pFaces[pFI];
            if (isMasterFace[facei] && externalFaces[facei])
            {
                nInterfacePointFaces[pointi]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        nInterfacePointFaces,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh_.points(), pointi)
    {
        if (nInterfacePointFaces[pointi] == 6 && !pinchedEdgePts[pointi])
        {
            const labelList& pCells = mesh_.pointCells()[pointi];
            forAll(pCells, pCI)
            {
                forRemoval[pCells[pCI]] = false;
            }
        }
    }

    syncTools::swapBoundaryCellList(mesh_, forRemoval, neiForRemoval);
    labelList nInternalFaces(mesh_.nCells(), 0);
    forAll(mesh_.cells(), celli)
    {
        if (forRemoval[celli])
        {
            nInternalFaces[celli] = -1;
            continue;
        }

        const cell& c = mesh_.cells()[celli];
        label nNeiKept = 0;
        forAll(c, cFI)
        {
            label facei = c[cFI];
            label patchI = patches.whichPatch(facei);
            if (patchI == -1)
            {
                label own = mesh_.faceOwner()[facei];
                label nei = -1;
                if (own == celli)
                {
                    nei  = mesh_.faceNeighbour()[facei];
                }
                else
                {
                    nei = own;
                }
                if (!forRemoval[nei])
                {
                    nNeiKept++;
                }
            }
            else if (patches[patchI].coupled())
            {
                label bfacei = facei-mesh_.nInternalFaces();
                if (!neiForRemoval[bfacei])
                {
                    nNeiKept++;
                }
            }
        }

        nInternalFaces[celli] = nNeiKept;
    }

    labelList neinInternalFaces;
    syncTools::swapBoundaryCellList(mesh_, nInternalFaces, neinInternalFaces);
    forAll(mesh_.cells(), celli)
    {
        if (nInternalFaces[celli] == -1)
        {
            continue;
        }

        if (nInternalFaces[celli] == 1)
        {
            forRemoval[celli] = true;
        }
        else if (nInternalFaces[celli] == 2)
        {
            const cell& c = mesh_.cells()[celli];
            bool removeCell = true;
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchI = patches.whichPatch(facei);
                if (patchI == -1)
                {
                    label own = mesh_.faceOwner()[facei];
                    label nei = -1;
                    if (own == celli)
                    {
                        nei  = mesh_.faceNeighbour()[facei];
                    }
                    else
                    {
                        nei = own;
                    }
                    if (nInternalFaces[nei] == 3 && !forRemoval[nei])
                    {
                        removeCell = false;
                    }
                }
                else if (patches[patchI].coupled())
                {
                    label bfacei = facei-mesh_.nInternalFaces();
                    if
                    (
                        neinInternalFaces[bfacei] == 3 &&
                        !neiForRemoval[bfacei]
                    )
                    {
                        removeCell = false;
                    }
                }
            }
            if (removeCell)
            {
                forRemoval[celli] = true;
            }
        }
    }

    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(mesh_.cells(), celli)
    {
        if (forRemoval[celli])
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label nCellsToKeep = mesh_.nCells() - cellsToRemove.size();
    reduce(nCellsToKeep, sumOp<label>());

    Info<< "Removing all cells containing points "
        << locationsOutsideMesh << endl
        << "Selected for keeping : " << nCellsToKeep << " cells." << endl;

    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatches(exposedFaces.size());
    label defaultPatch = 0;
    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];

        if (expandedBdyFaces[facei] != -1)
        {
            exposedPatches[i] = expandedBdyFaces[facei];
        }
        else
        {
            WarningInFunction
                << "For exposed face " << facei
                << " fc:" << mesh_.faceCentres()[facei]
                << " found no patch." << endl
                << "    Taking patch " << defaultPatch
                << " instead." << endl;
            exposedPatches[i] = defaultPatch;
        }
    }

    doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );

    return;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitMesh
(
    const label nBufferLayers,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,

    const refinementParameters& refineParams,
    const writer<scalar>& leakPathFormatter
)
{
    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const pointField& locationsOutsideMesh = refineParams.locationsOutsideMesh();

    // Determine patches to put intersections into
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find intersections with all unnamed(!) surfaces
    labelList ownPatch, nbrPatch;
    getBafflePatches
    (
        refineParams,
        globalToMasterPatch,
        ownPatch,
        nbrPatch
    );

    // Analyse regions. Reuse regionsplit
    boolList blockedFace(mesh_.nFaces(), false);

    forAll(ownPatch, facei)
    {
        if (ownPatch[facei] != -1 || nbrPatch[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());

    regionSplit cellRegion(mesh_, blockedFace);

    // Set unreachable cells to -1
    findRegions
    (
        mesh_,
        meshCutter_,
        mergeDistance_*vector(1,1,1),   // perturbVec
        locationsInMesh,
        locationsOutsideMesh,
        leakPathFormatter,
        refineParams.fullLeakChecks(),
        cellRegion.nRegions(),
        cellRegion,
        blockedFace
    );

    blockedFace.clear();

    // Walk out nBufferlayers from region boundary
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // (modifies cellRegion, ownPatch)
    // Takes over face patch onto points and then back to faces and cells
    // (so cell-face-point walk)

    const labelList& faceOwner = mesh_.faceOwner();
    const labelList& faceNeighbour = mesh_.faceNeighbour();

    // Patch for exposed faces for lack of anything sensible.
    label defaultPatch = 0;
    if (globalToMasterPatch.size())
    {
        defaultPatch = globalToMasterPatch[0];
    }

    for (label i = 0; i < nBufferLayers; i++)
    {
        // 1. From cells (via faces) to points

        labelList pointBaffle(mesh_.nPoints(), -1);

        // Get coupled neighbour cellRegion
        labelList neiCellRegion;
        syncTools::swapBoundaryCellList(mesh_, cellRegion, neiCellRegion);

        forAll(faceNeighbour, facei)
        {
            const face& f = mesh_.faces()[facei];

            label ownRegion = cellRegion[faceOwner[facei]];
            label neiRegion = cellRegion[faceNeighbour[facei]];

            if (ownRegion == -1 && neiRegion != -1)
            {
                // Note max(..) since possibly regionSplit might have split
                // off extra unreachable parts of mesh. Note: or can this only
                // happen for boundary faces?
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = max(defaultPatch, ownPatch[facei]);
                }
            }
            else if (ownRegion != -1 && neiRegion == -1)
            {
                label newPatchI = nbrPatch[facei];
                if (newPatchI == -1)
                {
                    newPatchI = max(defaultPatch, ownPatch[facei]);
                }
                forAll(f, fp)
                {
                    pointBaffle[f[fp]] = newPatchI;
                }
            }
        }

        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];
            bool coupledPatch(pp.coupled() ? true : false);
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label ownRegion = cellRegion[faceOwner[facei]];
                if (ownRegion != -1)
                {
                    label cellRegionNei = neiCellRegion[facei-mesh_.nInternalFaces()];
                    if (!coupledPatch || (coupledPatch && cellRegionNei == -1))
                    {
                        label newPatchI = max(defaultPatch, ownPatch[facei]);
                        const face& f = mesh_.faces()[facei];
                        forAll(f, fp)
                        {
                            pointBaffle[f[fp]] = newPatchI;
                        }
                    }
                }
            }
        }

        // Sync
        syncTools::syncPointList
        (
            mesh_,
            pointBaffle,
            maxEqOp<label>(),
            label(-1)           // null value
        );


        // 2. From points back to faces

        const labelListList& pointFaces = mesh_.pointFaces();

        forAll(pointFaces, pointI)
        {
            if (pointBaffle[pointI] != -1)
            {
                const labelList& pFaces = pointFaces[pointI];

                forAll(pFaces, pFaceI)
                {
                    label facei = pFaces[pFaceI];

                    if (ownPatch[facei] == -1)
                    {
                        ownPatch[facei] = pointBaffle[pointI];
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());


        // 3. From faces to cells (cellRegion) and back to faces (ownPatch)

        labelList newOwnPatch(ownPatch);

        forAll(ownPatch, facei)
        {
            if (ownPatch[facei] != -1)
            {
                label own = faceOwner[facei];

                if (cellRegion[own] == -1)
                {
                    cellRegion[own] = labelMax;

                    const cell& ownFaces = mesh_.cells()[own];
                    forAll(ownFaces, j)
                    {
                        if (ownPatch[ownFaces[j]] == -1)
                        {
                            newOwnPatch[ownFaces[j]] = ownPatch[facei];
                        }
                    }
                }
                if (mesh_.isInternalFace(facei))
                {
                    label nei = faceNeighbour[facei];

                    if (cellRegion[nei] == -1)
                    {
                        cellRegion[nei] = labelMax;

                        const cell& neiFaces = mesh_.cells()[nei];
                        forAll(neiFaces, j)
                        {
                            if (ownPatch[neiFaces[j]] == -1)
                            {
                                newOwnPatch[neiFaces[j]] = ownPatch[facei];
                            }
                        }
                    }
                }
            }
        }

        ownPatch.transfer(newOwnPatch);

        syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
    }


    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] == -1)
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label nCellsToKeep = mesh_.nCells() - cellsToRemove.size();
    reduce(nCellsToKeep, sumOp<label>());

    Info<< "Keeping all cells containing points " << locationsInMesh << endl
        << "Selected for keeping : " << nCellsToKeep << " cells." << endl;

    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatches(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];

        if (ownPatch[facei] != -1)
        {
            exposedPatches[i] = ownPatch[facei];
        }
        else
        {
            WarningInFunction
                << "For exposed face " << facei
                << " fc:" << mesh_.faceCentres()[facei]
                << " found no patch." << endl
                << "    Taking patch " << defaultPatch
                << " instead." << endl;
            exposedPatches[i] = defaultPatch;
        }
    }

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints
(
    const localPointRegion& regionSide
)
{
    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nNonManifPoints = returnReduce
    (
        regionSide.meshPointMap().size(),
        sumOp<label>()
    );

    Info<< "dupNonManifoldPoints : Found : " << nNonManifPoints
        << " non-manifold points (out of "
        << mesh_.globalData().nTotalPoints()
        << ')' << endl;


    autoPtr<mapPolyMesh> map;

    if (nNonManifPoints)
    {
        // Topo change engine
        duplicatePoints pointDuplicator(mesh_);

        // Insert changes into meshMod
        pointDuplicator.setRefinement(regionSide, meshMod);

        // Change the mesh (no inflation, parallel sync)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Is mapping only (no faces created, positions
        // stay same) so no need to recalculate intersections.
        updateMesh(map, labelList(0));
    }

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::dupNonManifoldPoints()
{
    // Analyse which points need to be duplicated
    localPointRegion regionSide(mesh_);

    return dupNonManifoldPoints(regionSide);
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::mergePoints
(
    const labelList& pointToDuplicate
)
{
    label nPointPairs = 0;
    forAll(pointToDuplicate, pointI)
    {
        label otherPointI = pointToDuplicate[pointI];
        if (otherPointI != -1)
        {
            nPointPairs++;
        }
    }

    autoPtr<mapPolyMesh> map;

    if (returnReduce(nPointPairs, sumOp<label>()))
    {
        Map<label> pointToMaster(2*nPointPairs);
        forAll(pointToDuplicate, pointI)
        {
            label otherPointI = pointToDuplicate[pointI];
            if (otherPointI != -1)
            {
                // Slave point
                pointToMaster.insert(pointI, otherPointI);
            }
        }

        // Topochange container
        polyTopoChange meshMod(mesh_);

        // Insert changes
        polyMeshAdder::mergePoints(mesh_, pointToMaster, meshMod);

        // Change the mesh (no inflation, parallel sync)
        map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // Update intersections. Is mapping only (no faces created, positions
        // stay same) so no need to recalculate intersections.
        updateMesh(map, labelList(0));
    }

    return map;
}


// Duplicate points on 'boundary' zones. Do not duplicate points on
// 'internal' or 'baffle' zone. Whether points are on normal patches does
// not matter
Foam::autoPtr<Foam::mapPolyMesh>
Foam::meshRefinement::dupNonManifoldBoundaryPoints()
{
    const labelList boundaryFaceZones
    (
        getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::BOUNDARY
            )
        )
    );
    labelList internalOrBaffleFaceZones;
    {
        List<surfaceZonesInfo::faceZoneType> fzTypes(2);
        fzTypes[0] = surfaceZonesInfo::INTERNAL;
        fzTypes[1] = surfaceZonesInfo::BAFFLE;
        internalOrBaffleFaceZones = getZones(fzTypes);
    }



    // 0 : point used by normal, unzoned boundary faces
    // 1 : point used by 'boundary' zone
    // 2 : point used by internal/baffle zone
    PackedList<2> pointStatus(mesh_.nPoints(), 0u);

    forAll(boundaryFaceZones, j)
    {
        const faceZone& fZone = mesh_.faceZones()[boundaryFaceZones[j]];
        forAll(fZone, i)
        {
            const face& f = mesh_.faces()[fZone[i]];
            forAll(f, fp)
            {
                pointStatus[f[fp]] = max(pointStatus[f[fp]], 1u);
            }
        }
    }
    forAll(internalOrBaffleFaceZones, j)
    {
        const faceZone& fZone = mesh_.faceZones()[internalOrBaffleFaceZones[j]];
        forAll(fZone, i)
        {
            const face& f = mesh_.faces()[fZone[i]];
            forAll(f, fp)
            {
                pointStatus[f[fp]] = max(pointStatus[f[fp]], 2u);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pointStatus,
        maxEqOp<unsigned int>(),    // combine op
        0u                          // null value
    );

    // Pick up points on boundary zones that are not on internal/baffle zones
    label n = 0;
    forAll(pointStatus, pointI)
    {
        if (pointStatus[pointI] == 1u)
        {
            n++;
        }
    }

    label globalNPoints = returnReduce(n, sumOp<label>());
    Info<< "Duplicating " << globalNPoints << " points on"
        << " faceZones of type "
        << surfaceZonesInfo::faceZoneTypeNames[surfaceZonesInfo::BOUNDARY]
        << endl;

    autoPtr<mapPolyMesh> map;

    if (globalNPoints)
    {
        labelList candidatePoints(n);
        n = 0;
        forAll(pointStatus, pointI)
        {
            if (pointStatus[pointI] == 1u)
            {
                candidatePoints[n++] = pointI;
            }
        }
        localPointRegion regionSide(mesh_, candidatePoints);
        map = dupNonManifoldPoints(regionSide);
    }
    return map;
}


// Zoning
void Foam::meshRefinement::zonify
(
    const refinementParameters& refineParams,
    const writer<scalar>& setFormatter,
    wordPairHashTable& zonesToFaceZone,
    bool keepAlreadyZonedCells,
    bool resetUnsetZones,
    bool removeUnsetCells
)
{
    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const wordList& zonesInMesh = refineParams.zonesInMesh();

    if (locationsInMesh.size() != zonesInMesh.size())
    {
        FatalErrorInFunction << "problem" << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();

    // Add any faceZones, cellZones originating from surface to the mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    labelList surfaceToCellZone;
    labelList surfaceToFaceZone;

    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));
    if (namedSurfaces.size())
    {
        Info<< "Setting cellZones according to named surfaces:" << endl;
        forAll(namedSurfaces, i)
        {
            label surfI = namedSurfaces[i];

            Info<< "Surface : " << surfaces_.names()[surfI] << nl
                << "    faceZone : " << surfZones[surfI].faceZoneName() << nl
                << "    cellZone : " << surfZones[surfI].cellZoneName() << endl;
        }
        Info<< endl;

        // Add zones to mesh
        surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
        surfaceToFaceZone = surfaceZonesInfo::addFaceZonesToMesh
        (
            surfZones,
            namedSurfaces,
            mesh_
        );
    }

    // Put the cells into the correct zone
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Zone per cell:
    // -2 : unset : not allowed!
    // -1 : not in any zone (zone 'none')
    // >=0: zoneID
    // namedSurfaceIndex:
    // -1  : face not intersecting a named surface
    // >=0 : index of named surface
    labelList cellToZone(mesh_.nCells(),-2);
    if (keepAlreadyZonedCells)
    {
        const cellZoneMesh& cellZones = mesh_.cellZones();

        forAll(cellZones, cellZoneI)
        {
            const cellZone& cZone = cellZones[cellZoneI];
            forAll(cZone, czi)
            {
                label celli = cZone[czi];
                cellToZone[celli] = cellZoneI;
            }
        }
    }
    labelList origCellZone(cellToZone);

    const label minZoneRegionSize = refineParams.minZoneRegionSize();

    labelList namedSurfaceIndex;
    labelList namedIntersections;
    PackedBoolList posOrientation;
    labelList testFaces;
    labelList globalRegion1;
    labelList globalRegion2;

    zonify
    (
        -1,                 // Set all cells with cellToZone -2 to -1
        refineParams,
        minZoneRegionSize,

        testFaces,
        globalRegion1,
        globalRegion2,

        cellToZone,
        namedSurfaceIndex,
        namedIntersections,
        posOrientation
    );

    if (refineParams.zoneLeakChecks())
    {
        checkZoneLeakPaths
        (
            refineParams,
            namedIntersections,
            globalRegion1,
            globalRegion2,
            cellToZone,
            setFormatter
        );
    }

    if (namedSurfaces.size() && refineParams.filterFreeStandingHoles())
    {
        filterFreeStandingHoles
        (
            refineParams,
            namedSurfaceIndex
        );
    }

    //Check for boundary zoned cells that might have been incorrectly allocated
    if (keepAlreadyZonedCells)
    {
        const cellZoneMesh& cellZones = mesh_.cellZones();
        boolList boundaryNamed(cellZones.size(), false);

        //Check for named surfaces
        const PtrList<surfaceZonesInfo>& surfZones = surfaces().surfZones();

        forAll(surfZones, surfI)
        {
            const word& faceZoneName = surfZones[surfI].faceZoneName();

            if (faceZoneName.size())
            {
                const surfaceZonesInfo::faceZoneType& faceType =
                    surfZones[surfI].faceType();

                if (faceType == surfaceZonesInfo::BOUNDARY)
                {
                    const word& cellZoneName = surfZones[surfI].cellZoneName();
                    if (cellZoneName.size())
                    {
                        label zoneI = cellZones.findZoneID(cellZoneName);
                        if (zoneI != -1)
                        {
                            boundaryNamed[zoneI] =  true;
                        }
                    }
                }
            }
        }

        for (int i = 0; i < 2; i++)
        {
            labelList neiCellZone;
            syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

            label nReset = 0;

            forAll(mesh_.cells(), celli)
            {
                label zonei = cellToZone[celli];
                if
                (
                    zonei > -1
                    && boundaryNamed[zonei]
                    && origCellZone[celli] != cellToZone[celli]
                )
                {
                    const cell c = mesh_.cells()[celli];
                    bool nonZonedNei = false;
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        label patchI = patches.whichPatch(facei);
                        if (patchI == -1)
                        {
                            label own = mesh_.faceOwner()[facei];
                            label nei
                            (
                                own == celli
                                ? mesh_.faceNeighbour()[facei] : own
                            );

                            if (cellToZone[nei] == -1)
                            {
                                nonZonedNei = true;
                                break;
                            }
                        }
                        else if (patches[patchI].coupled())
                        {
                            label bFaceI = facei-mesh_.nInternalFaces();
                            if (neiCellZone[bFaceI] == -1)
                            {
                               nonZonedNei = true;
                               break;
                            }
                        }
                    }
                    if (nonZonedNei)
                    {
                        cellToZone[celli] = -1;
                        nReset++;
                    }
                }
            }

            if (returnReduce(nReset, sumOp<label>()) == 0)
            {
                break;
            }
        }
    }

    DynamicList<label> cellsToRemove(mesh_.nCells()/1000);
    boolList removedCellFaces(mesh_.nFaces(), false);

    //Resetting unset zones. Used for dual mesh where cells can be produced
    // the wrong side of boundary. Set unset to max of neighbour cellZone values
    if (resetUnsetZones)
    {
        for (int i = 0; i < 2; i++)
        {
            labelList neiCellZone;
            syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
            label nReset = 0;

            forAll(mesh_.cells(), celli)
            {
                if (cellToZone[celli] == -2)
                {
                    const cell c = mesh_.cells()[celli];
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        label patchI = patches.whichPatch(facei);
                        if (patchI == -1)
                        {
                            label own = mesh_.faceOwner()[facei];
                            label nei
                            (
                                own == celli
                                ? mesh_.faceNeighbour()[facei] : own
                            );
                            cellToZone[celli] =
                                max(cellToZone[celli], cellToZone[nei]);
                            nReset++;
                        }
                        else if (patches[patchI].coupled())
                        {
                            label bFaceI = facei-mesh_.nInternalFaces();
                            cellToZone[celli] =
                                max(cellToZone[celli], neiCellZone[bFaceI]);
                            nReset++;
                        }
                    }
                }
            }

            if (returnReduce(nReset, sumOp<label>()) == 0)
            {
                break;
            }
        }

        // Check for boundary baffle faces in dual mesh and reset
        // if connected to named cell zone
        if (namedSurfaces.size())
        {
            for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
            {
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                label neiZone = cellToZone[mesh_.faceNeighbour()[facei]];

                if (namedSurfaceIndex[facei] == -1 && (ownZone != neiZone))
                {
                    // Give face the zone of min cell zone (but only if the
                    // cellZone originated from a closed, named surface)

                    label minZone;
                    if (ownZone == -1)
                    {
                        minZone = neiZone;
                    }
                    else if (neiZone == -1)
                    {
                        minZone = ownZone;
                    }
                    else
                    {
                        minZone = min(ownZone, neiZone);
                    }

                    // Make sure the cellZone originated from a closed surface
                    label geomSurfI = findIndex(surfaceToCellZone, minZone);

                    if (geomSurfI != -1)
                    {
                        namedSurfaceIndex[facei] = geomSurfI;
                    }
                }
            }

            labelList neiCellZone;
            syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

            const polyBoundaryMesh& patches = mesh_.boundaryMesh();

            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (pp.coupled())
                {
                    forAll(pp, i)
                    {
                        label facei = pp.start()+i;
                        label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                        label neiZone =
                            neiCellZone[facei-mesh_.nInternalFaces()];

                        if
                        (
                            namedSurfaceIndex[facei] == -1
                            && (ownZone != neiZone)
                        )
                        {
                            // Give face the min cell zone
                            label minZone;
                            if (ownZone == -1)
                            {
                                minZone = neiZone;
                            }
                            else if (neiZone == -1)
                            {
                                minZone = ownZone;
                            }
                            else
                            {
                                minZone = min(ownZone, neiZone);
                            }

                            // Make sure the cellZone originated
                            // from a closed surface
                            label geomSurfI =
                                findIndex(surfaceToCellZone, minZone);

                            if (geomSurfI != -1)
                            {
                                namedSurfaceIndex[facei] = geomSurfI;
                            }
                        }
                    }
                }
            }
            // Sync
            syncTools::syncFaceList(mesh_, namedSurfaceIndex, maxEqOp<label>());
        }
    }
    else
    {
        //Mark all unvisited cells to -1. This might indicate a problem
        // so it is probably better to use minZoneRegionSize to regroup zones
        // instead
        forAll(mesh_.cells(), celli)
        {
            if (cellToZone[celli] == -2)
            {
                cellToZone[celli] = -1;
                if (removeUnsetCells)
                {
                    cellsToRemove.append(celli);
                    const cell& c = mesh_.cells()[celli];
                    forAll(c,cFI)
                    {
                        removedCellFaces[c[cFI]] = true;
                    }
                }
            }
        }
    }

    if (removeUnsetCells)
    {
        syncTools::syncFaceList(mesh_, removedCellFaces, orEqOp<bool>());
    }

    // Convert namedSurfaceIndex (index of named surfaces) to
    // actual faceZone index

    //- Per face index of faceZone or -1
    labelList faceToZone(mesh_.nFaces(), -1);

    forAll(namedSurfaceIndex, facei)
    {
        label surfI = namedSurfaceIndex[facei];
        if (surfI != -1)
        {
            faceToZone[facei] = surfaceToFaceZone[surfI];
        }
    }

    // Allocate and assign faceZones from cellZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // 1. Detect inter-region face and allocate names

        HashTable<word, labelPair, labelPair::Hash<>> zoneIDsToFaceZone;

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            if (faceToZone[facei] == -1 && !removedCellFaces[facei])// Not named
            {
                // Face not yet in a faceZone. (it might already have been
                // done so by a 'named' surface). Check if inbetween different
                // cellZones
                allocateInterRegionFaceZone
                (
                    cellToZone[mesh_.faceOwner()[facei]],
                    cellToZone[mesh_.faceNeighbour()[facei]],
                    zonesToFaceZone,
                    zoneIDsToFaceZone
                );
            }
        }

        labelList neiCellZone;
        syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);

        forAll(neiCellZone, bFaceI)
        {
            label facei = bFaceI + mesh_.nInternalFaces();

            if (faceToZone[facei] == -1 && !removedCellFaces[facei])
            {
                allocateInterRegionFaceZone
                (
                    cellToZone[mesh_.faceOwner()[facei]],
                    neiCellZone[bFaceI],
                    zonesToFaceZone,
                    zoneIDsToFaceZone
                );
            }
        }


        // 2.Combine faceZoneNames allocated on different processors

        Pstream::mapCombineGather(zonesToFaceZone, eqOp<word>());
        Pstream::mapCombineScatter(zonesToFaceZone);


        // 3. Allocate faceZones from (now synchronised) faceZoneNames
        //    Note: the faceZoneNames contain the same data but in different
        //          order. We could sort the contents but instead just loop
        //          in sortedToc order.

        Info<< "Setting faceZones according to neighbouring cellZones:"
            << endl;

        // From cellZone indices to faceZone index
        HashTable<label, labelPair, labelPair::Hash<>> fZoneLookup
        (
            zonesToFaceZone.size()
        );

        const cellZoneMesh& cellZones = mesh_.cellZones();

        {
            List<Pair<word>> czs(zonesToFaceZone.sortedToc());

            forAll(czs, i)
            {
                const Pair<word>& cz = czs[i];
                const word& fzName = zonesToFaceZone[cz];

                Info<< indent<< "cellZones : "
                    << cz[0] << ' ' << cz[1] << nl
                    << "    faceZone : " << fzName << endl;

                label faceZoneI = surfaceZonesInfo::addFaceZone
                (
                    fzName,                 // name
                    labelList(0),           // addressing
                    boolList(0),            // flipMap
                    mesh_
                );

                label cz0 = cellZones.findZoneID(cz[0]);
                label cz1 = cellZones.findZoneID(cz[1]);

                fZoneLookup.insert(labelPair(cz0, cz1), faceZoneI);
            }
        }

        // 4. Set faceToZone with new faceZones

        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            if (faceToZone[facei] == -1 && !removedCellFaces[facei])
            {
                // Face not yet in a faceZone. (it might already have been
                // done so by a 'named' surface). Check if inbetween different
                // cellZones

                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                label neiZone = cellToZone[mesh_.faceNeighbour()[facei]];
                if (ownZone != neiZone)
                {
                    bool swap =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                    labelPair key(ownZone, neiZone);
                    if (swap)
                    {
                        Swap(key.first(), key.second());
                    }
                    faceToZone[facei] = fZoneLookup[key];
                }
            }
        }
        forAll(neiCellZone, bFaceI)
        {
            label facei = bFaceI + mesh_.nInternalFaces();
            if (faceToZone[facei] == -1 && !removedCellFaces[facei])
            {
                label ownZone = cellToZone[mesh_.faceOwner()[facei]];
                label neiZone = neiCellZone[bFaceI];
                if (ownZone != neiZone)
                {
                    bool swap =
                    (
                        ownZone == -1
                     || (neiZone != -1 && ownZone > neiZone)
                    );
                    labelPair key(ownZone, neiZone);
                    if (swap)
                    {
                        Swap(key.first(), key.second());
                    }
                    faceToZone[facei] = fZoneLookup[key];
                }
            }
        }
        Info<< endl;
    }

    // Get coupled neighbour cellZone. Set to -1 on non-coupled patches.
    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (!pp.coupled())
        {
            label bFaceI = pp.start()-mesh_.nInternalFaces();
            forAll(pp, i)
            {
                neiCellZone[bFaceI++] = -1;
            }
        }
    }

    // Get per face whether is it master (of a coupled set of faces)
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    // faceZones
    // ~~~~~~~~~
    // Faces on faceZones come in two variants:
    // - faces on the outside of a cellZone. They will be oriented to
    //   point out of the maximum cellZone.
    // - free-standing faces. These will be oriented according to the
    //   local surface normal. We do this in a two step algorithm:
    //      - do a consistent orientation
    //      - check number of faces with consistent orientation
    //      - if <0 flip the whole patch
    PackedBoolList meshFlipMap(mesh_.nFaces(), false);
    {
        // Collect all data on zone faces without cellZones on either side.
        const indirectPrimitivePatch patch
        (
            IndirectList<face>
            (
                mesh_.faces(),
                freeStandingBaffleFaces
                (
                    faceToZone,
                    cellToZone,
                    neiCellZone
                )
            ),
            mesh_.points()
        );

        label nFreeStanding = returnReduce(patch.size(), sumOp<label>());
        if (nFreeStanding > 0)
        {
            Info<< "Detected " << nFreeStanding << " free-standing zone faces"
                << endl;

            if (debug)
            {
                OBJstream str(mesh_.time().path()/"freeStanding.obj");
                str.write(patch.localFaces(), patch.localPoints(), false);
            }


            // Detect non-manifold edges
            labelList nMasterFacesPerEdge;
            calcPatchNumMasterFaces(isMasterFace, patch, nMasterFacesPerEdge);

            // Mark zones. Even a single original surface might create multiple
            // disconnected/non-manifold-connected zones
            labelList faceToConnectedZone;
            const label nZones = markPatchZones
            (
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone
            );

            Map<label> nPosOrientation(2*nZones);
            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                nPosOrientation.insert(zoneI, 0);
            }

            // Make orientations consistent in a topological way. This just
            // checks  the first face per zone for whether nPosOrientation
            // is negative (which it never is at this point)
            consistentOrientation
            (
                isMasterFace,
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone,
                nPosOrientation,

                meshFlipMap
            );

            // Count per region the number of orientations (taking the new
            // flipMap into account)
            forAll(patch.addressing(), facei)
            {
                label meshFaceI = patch.addressing()[facei];

                if (isMasterFace[meshFaceI])
                {
                    label n = 1;
                    if
                    (
                        bool(posOrientation[meshFaceI])
                     == meshFlipMap[meshFaceI]
                    )
                    {
                        n = -1;
                    }

                    nPosOrientation.find(faceToConnectedZone[facei])() += n;
                }
            }
            Pstream::mapCombineGather(nPosOrientation, plusEqOp<label>());
            Pstream::mapCombineScatter(nPosOrientation);


            Info<< "Split " << nFreeStanding << " free-standing zone faces"
                << " into " << nZones << " disconnected regions with size"
                << " (negative denotes wrong orientation) :"
                << endl;

            for (label zoneI = 0; zoneI < nZones; zoneI++)
            {
                Info<< "    " << zoneI << "\t" << nPosOrientation[zoneI]
                    << endl;
            }
            Info<< endl;


            // Re-apply with new counts (in nPosOrientation). This will cause
            // zones with a negative count to be flipped.
            consistentOrientation
            (
                isMasterFace,
                patch,
                nMasterFacesPerEdge,
                faceToConnectedZone,
                nPosOrientation,

                meshFlipMap
            );
        }
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);

    // Insert changes to put cells and faces into zone
    zonify
    (
        isMasterFace,
        cellToZone,
        neiCellZone,
        faceToZone,
        meshFlipMap,
        meshMod
    );

    // Change the mesh (no inflation, parallel sync)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Print some stats (note: zones are synchronised)
    if (mesh_.cellZones().size() > 0)
    {
        Info<< "CellZones:" << endl;
        forAll(mesh_.cellZones(), zoneI)
        {
            const cellZone& cz = mesh_.cellZones()[zoneI];
            Info<< "    " << cz.name()
                << "\tsize:" << returnReduce(cz.size(), sumOp<label>())
                << endl;
        }
        Info<< endl;
    }
    if (mesh_.faceZones().size() > 0)
    {
        Info<< "FaceZones:" << endl;
        forAll(mesh_.faceZones(), zoneI)
        {
            const faceZone& fz = mesh_.faceZones()[zoneI];
            Info<< "    " << fz.name()
                << "\tsize:" << returnReduce(fz.size(), sumOp<label>())
                << endl;
        }
        Info<< endl;
    }

    // None of the faces has changed, only the zones. Still...
    updateMesh(map, labelList());

    if (returnReduce(cellsToRemove.size(), sumOp<label>()) != 0)
    {
        // Remove unset cells
        removeCells cellRemover(mesh_);

        labelList exposedFaceID(mesh_.nFaces(), 0);
        forAll(cellsToRemove, celli)
        {
            const cell& c = mesh_.cells()[celli];
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchi = patches.whichPatch(facei);
                if (patchi != -1 && !patches[patchi].coupled())
                {
                    exposedFaceID[facei] = max(exposedFaceID[facei], patchi);
                }
            }
        }
        syncTools::syncFaceList(mesh_, exposedFaceID, maxEqOp<label>());

        // Pick up patches for exposed faces
        labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
        labelList exposedPatches(exposedFaces.size());

        forAll(exposedFaces, i)
        {
            label facei = exposedFaces[i];
            exposedPatches[i] = exposedFaceID[facei];
        }

        doRemoveCells
        (
            cellsToRemove,
            exposedFaces,
            exposedPatches,
            cellRemover
        );
    }

    return;
}


void Foam::meshRefinement::filterFreeStandingHoles
(
    const refinementParameters& refineParams,
    labelList& namedSurfaceIndex
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    //Check for named surfaces
    const PtrList<surfaceZonesInfo>& surfZones = surfaces().surfZones();

    //Mark free-standing blocked faces
    boolList blockedFace(mesh_.nFaces(),false);
    label nFreeStanding = 0;
    forAll(mesh_.faces(), facei)
    {
        label surfI = namedSurfaceIndex[facei];
        if
        (
           surfI != -1
           && surfZones[surfI].freeStanding()
           && surfZones[surfI].faceType() != surfaceZonesInfo::BOUNDARY
        )
        {
            blockedFace[facei] = true;
            nFreeStanding++;
        }
    }

    if (returnReduce(nFreeStanding, sumOp<label>()) > 0)
    {
        // Set region per cell based on walking
        regionSplit cellRegion(mesh_, blockedFace);

        const pointField& locationsInMesh = refineParams.locationsInMesh();
        DynamicList<label> regionsToKeep(cellRegion.nRegions());

        // For all locationsInMesh find the cell
        forAll(locationsInMesh, i)
        {
            // Get location and index of zone ("none" for cellZone -1)
            const point& insidePoint = locationsInMesh[i];

            // Find the region containing the keepPoint
            label keepRegionI = -1;

            label celli = findCell
            (
                insidePoint,
                mesh_,
                meshCutter_
            );

            if (celli != -1)
            {
                keepRegionI = cellRegion[celli];
            }
            reduce(keepRegionI, maxOp<label>());

            bool newRegion = true;
            forAll(regionsToKeep, regioni)
            {
                if (regionsToKeep[regioni] == keepRegionI)
                {
                    newRegion = false;
                    break;
                }
            }
            if (newRegion)
            {
                regionsToKeep.append(keepRegionI);
            }
        }
        regionsToKeep.shrink();

        labelHashSet regionsToKeepSet(regionsToKeep);
        if (cellRegion.nRegions() > regionsToKeep.size())
        {
            boolList markedRegions(mesh_.nCells(),true);
//            label nUnmarked = 0;
            forAll(mesh_.cells(), celli)
            {
                if (!regionsToKeepSet.found(cellRegion[celli]))
                {
                    markedRegions[celli] = false;
//                    nUnmarked++;
                }
            }

            boolList neiMarkedRegions;
            syncTools::swapBoundaryCellList
            (
               mesh_,
               markedRegions,
               neiMarkedRegions
            );
            DynamicList<label> outerShellFaces(nFreeStanding);
            forAll(mesh_.faces(), facei)
            {
                if (blockedFace[facei])
                {
                    label patchI = patches.whichPatch(facei);
                    label own = mesh_.faceOwner()[facei];
                    if (patchI == -1)
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        if (!markedRegions[own] && !markedRegions[nei])
                        {
                            namedSurfaceIndex[facei] = -1;
                        }
                        else if
                        (
                            (markedRegions[own] && !markedRegions[nei])
                            || (!markedRegions[own] && markedRegions[nei])
                        )
                        {
                            outerShellFaces.append(facei);
                        }
                    }
                    else if (patches[patchI].coupled())
                    {
                        label bFaceI = facei-mesh_.nInternalFaces();
                        if (!neiMarkedRegions[bFaceI] && !markedRegions[own])
                        {
                            namedSurfaceIndex[facei] = -1;
                        }
                        else if
                        (
                           (!neiMarkedRegions[bFaceI] && markedRegions[own])
                           || (neiMarkedRegions[bFaceI] && !markedRegions[own])
                        )
                        {
                            outerShellFaces.append(facei);
                        }
                    }
                }
            }

            //find closest snapped patches
            const labelList freeStandingSurfs =
               surfaceZonesInfo::getNonBoundaryFreeStandingSurfaces
               (
                   surfaces_.surfZones()
               );

            pointField outerPts(outerShellFaces.size(), vector::zero);
            forAll(outerShellFaces, i)
            {
                label facei = outerShellFaces[i];
                outerPts[i] = mesh_.faceCentres()[facei];
            }

            List<pointIndexHit> nearestInfo;
            labelList nearestSurface;
            surfaces_.findNearest
            (
                freeStandingSurfs,
                outerPts,
                scalarField(outerPts.size(), sqr(GREAT)), // sqr of attraction
                nearestSurface,
                nearestInfo
            );

            // Do normal testing per surface.
            vectorField nearestNormal(nearestInfo.size(),vector::zero);

            forAll(freeStandingSurfs, sI)
            {
                label surfI = freeStandingSurfs[sI];
                DynamicList<pointIndexHit> localHits;

                forAll(nearestSurface, i)
                {
                    if (nearestSurface[i] == surfI)
                    {
                        localHits.append(nearestInfo[i]);
                    }
                }

                label geomI = surfaces_.surfaces()[surfI];

                pointField localNormals;
                surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

                label localI = 0;
                forAll(nearestSurface, i)
                {
                    if (nearestSurface[i] == surfI)
                    {
                        nearestNormal[i] = localNormals[localI];
                        localI++;
                    }
                }
            }

            forAll(outerShellFaces, i)
            {
                if (nearestInfo[i].hit())
                {
                    label facei = outerShellFaces[i];
                    label own = mesh_.faceOwner()[facei];
                    vector outNorm = mesh_.faceAreas()[facei];
                    scalar outNormMag = mag(outNorm);
                    if (outNormMag > SMALL)
                    {
                        outNorm /= outNormMag;
                        if (markedRegions[own])
                        {
                            outNorm = -outNorm;
                        }
                        if ((outNorm & nearestNormal[i]) > SMALL)
                        {
                            namedSurfaceIndex[facei] = -1;
                        }
                    }
                }
            }
        }
    }

    return;
}


void Foam::meshRefinement::checkZoneFaces() const
{
    const faceZoneMesh& fZones = mesh_.faceZones();

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        if (isA<processorPolyPatch>(pp))
        {
            forAll(pp, i)
            {
                label facei = pp.start()+i;
                label zoneI = fZones.whichZone(facei);

                if (zoneI != -1)
                {
                    FatalErrorInFunction
                        << "face:" << facei << " on patch " << pp.name()
                        << " is in zone " << fZones[zoneI].name()
                        << exit(FatalError);
                }
            }
        }
    }
}

void Foam::meshRefinement::mergeSingleCouples()
{
    Info<< "Merging single couples" << endl;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    DynamicList<label> allPatchIDs(patches.size());

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            allPatchIDs.append(patchI);
        }
    }
    allPatchIDs.shrink();

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        makePatch
        (
            mesh_,
            allPatchIDs
         )
     );
    indirectPrimitivePatch& pp = ppPtr();

    // Identify coupled faces
    List<labelPair> couples(localPointRegion::findDuplicateFacePairs(mesh_));

    labelList coupled(mesh_.nFaces(), -1);
    forAll(couples, i)
    {
        coupled[couples[i].first()] =  i;
        coupled[couples[i].second()] =  i;
    }

    const labelList meshEdges(pp.meshEdges(mesh_.edges(), mesh_.pointEdges()));

    labelList nEdgeCoupledFaces(mesh_.nEdges(), 0);
    labelList nEdgePatchFaces(mesh_.nEdges(), 0);
    forAll(meshEdges, i)
    {
        label meshEdgeI = meshEdges[i];
        const labelList& edgeFaces = pp.edgeFaces()[i];

        nEdgePatchFaces[meshEdgeI] = edgeFaces.size();
        forAll(edgeFaces, j)
        {
            label meshFaceI = pp.addressing()[edgeFaces[j]];

            if (coupled[meshFaceI] != -1)
            {
                nEdgeCoupledFaces[meshEdgeI]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        nEdgeCoupledFaces,
        plusEqOp<label>(),
        label(0)                  // initial value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        nEdgePatchFaces,
        plusEqOp<label>(),
        label(0)                  // initial value
    );

    boolList foundCouple(couples.size(), false);
    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];

        if (coupled[meshFaceI] != -1)
        {
            const labelList& faceEdges = pp.faceEdges()[i];
            label n2 = 0, n4 = 0;
            forAll(faceEdges, j)
            {
                label meshEdgeI = meshEdges[faceEdges[j]];
                if
                (
                    nEdgeCoupledFaces[meshEdgeI] == 2
                    && nEdgePatchFaces[meshEdgeI] == 2
                )
                {
                    n2++;
                }
                if
                (
                    nEdgeCoupledFaces[meshEdgeI] == 4
                    && nEdgePatchFaces[meshEdgeI] == 4
                )
                {
                    n4++;
                }
            }

            if (n2 == 2 && n4 == 0)
            {
                foundCouple[coupled[meshFaceI]] = true;
            }
        }
    }

    const labelList allSurfaces(identity(surfaces_.surfaces().size()));

    List<pointIndexHit> nearestInfo;
    labelList nearestSurface;

    surfaces_.findNearest
    (
        allSurfaces,
        pp.faceCentres(),
        scalarField(pp.size(), sqr(GREAT)),    // sqr of attraction
        nearestSurface,
        nearestInfo
    );

    // Do normal testing per surface.
    vectorField nearestNormal(nearestInfo.size());

    forAll(allSurfaces, surfI)
    {
        DynamicList<pointIndexHit> localHits;

        forAll(nearestSurface, i)
        {
            if (nearestSurface[i] == surfI)
            {
                localHits.append(nearestInfo[i]);
            }
        }

        label geomI = surfaces_.surfaces()[surfI];

        pointField localNormals;
        surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

        label localI = 0;
        forAll(nearestSurface, i)
        {
            if (nearestSurface[i] == surfI)
            {
                nearestNormal[i] = localNormals[localI];
                localI++;
            }
        }
    }

    forAll(pp, i)
    {
        label meshFaceI = pp.addressing()[i];
        if (coupled[meshFaceI] != -1 && foundCouple[coupled[meshFaceI]])
        {
            if (nearestSurface[i] != -1)
            {
                vector hNorm = nearestNormal[i]
                    / (mag(nearestNormal[i])+SMALL);

                vector fNorm = mesh_.faceAreas()[meshFaceI]
                    / (mesh_.magFaceAreas()[meshFaceI] + SMALL);

                if (mag(hNorm & fNorm) > 0.342)
                {
                    foundCouple[coupled[meshFaceI]] = false;
                }
            }
        }
    }

    List<labelPair> couplesMerge(couples.size());
    label nMerge = 0;

    forAll(couples, i)
    {
        if (foundCouple[i])
        {
            couplesMerge[nMerge] = couples[i];
            nMerge++;
        }
    }
    couplesMerge.setSize(nMerge);

    mergeBaffles(couplesMerge, Map<label>(0), true);
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::removeDisconnectedRegions()
{
    Info<< "Removing disconnected regions" << endl;

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_);
    labelHashSet keepRegionSet(keepLargestRegions(cellRegion));

    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, celli)
    {
        if (!keepRegionSet.found(cellRegion[celli]))
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label nCellsToKeep = mesh_.nCells() - cellsToRemove.size();
    reduce(nCellsToKeep, sumOp<label>());

    Info<< "Keeping " << nCellsToKeep << " cells out of : "
        << returnReduce(mesh_.nCells(), sumOp<label>()) <<endl;

    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatches(exposedFaces.size(),0);

    return doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatches,
        cellRemover
    );
}


void Foam::meshRefinement::removeOtherSide
(
    const refinementParameters& refineParams,
    const labelList& globalToMasterPatch
)
{
    labelList testFaces = identity(mesh_.nFaces());//intersectedFaces();

    labelList globalRegion1;
    labelList globalRegion2;
    getIntersections
    (
        mesh_.cellCentres(),
        surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
        (
            surfaces_.surfZones()
        ),
        testFaces,
        globalRegion1,
        globalRegion2,
        false
    );

    boolList blockedFace(mesh_.nFaces(),false);
    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);

    const pointField& locationsInMesh = refineParams.locationsInMesh();
    DynamicList<label> regionsToKeep(cellRegion.nRegions());

    // For all locationsInMesh find the cell
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];

        // Find the region containing the keepPoint
        label keepRegionI = -1;

        label celli = findCell
        (
            insidePoint,
            mesh_,
            meshCutter_
        );

        if (celli != -1)
        {
            keepRegionI = cellRegion[celli];
        }
        reduce(keepRegionI, maxOp<label>());
/*
        // Find the region containing the keepPoint
        label keepRegionI = findRegion
        (
            mesh_,
            cellRegion,
            mergeDistance()*vector(1,1,1),
            insidePoint
        );
*/
        regionsToKeep.append(keepRegionI);
    }
    regionsToKeep.shrink();
    labelHashSet regionsToKeepSet(regionsToKeep);

    DynamicList<label> cellsToRemove(mesh_.nCells()/100);
    forAll(mesh_.cells(), celli)
    {
        if (!regionsToKeepSet.found(cellRegion[celli]))
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    Info<<"Removing "<< returnReduce(cellsToRemove.size(), sumOp<label>())
        <<" cells not connected to inside location" <<endl;

    labelList ownPatch(mesh_.nFaces(), -1);
    boolList marked(cellsToRemove.size(), false);
    label nUnset = 0;
    forAll(cellsToRemove, i)
    {
        label celli = cellsToRemove[i];

        cell c = mesh_.cells()[celli];

        label patchI = -1;

        forAll(c, cFI)
        {
            label facei = c[cFI];

            if (globalRegion1[facei] != -1)
            {
                patchI = globalToMasterPatch[globalRegion1[facei]];
                break;
            }
        }

        if (patchI != -1)
        {
            marked[i] = true;
            forAll(c, cFI)
            {
                label facei = c[cFI];
                ownPatch[facei] = patchI;
            }
        }
        else
        {
            nUnset++;
        }
    }
    syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());

    reduce(nUnset, sumOp<label>());
    if (nUnset != 0)
    {
        Info<<"Resetting " << nUnset
            <<" cells not connected to boundary by a face" <<endl;
        while (true)
        {
            label nSet = 0;
            forAll(cellsToRemove, i)
            {
                if (marked[i])
                {
                    continue;
                }
                label celli = cellsToRemove[i];
                cell c = mesh_.cells()[celli];

                label patchI = -1;

                forAll(c, cFI)
                {
                    label facei = c[cFI];

                    if (ownPatch[facei] != -1)
                    {
                        patchI = ownPatch[facei];
                        break;
                    }
                }

                if (patchI != -1)
                {
                    nSet++;
                    marked[i] = true;
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        ownPatch[facei] = patchI;
                    }
                }
            }
            if (returnReduce(nSet, sumOp<label>()) == 0)
            {
                forAll(cellsToRemove, i)
                {
                    if (!marked[i])
                    {
                        label celli = cellsToRemove[i];
                        cell c = mesh_.cells()[celli];
                        forAll(c, cFI)
                        {
                            label facei = c[cFI];
                            ownPatch[facei] = label(0);
                        }
                    }
                }
                syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
                break;
            }
            syncTools::syncFaceList(mesh_, ownPatch, maxEqOp<label>());
        }
    }

    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatchIDs(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];

        exposedPatchIDs[i] = ownPatch[facei];
    }

    doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        cellRemover
    );

    return;
}


void Foam::meshRefinement::checkRemovedAndRefine
(
    const refinementParameters& refineParams,
    const labelList& globalToMasterPatch,
    const label maxIter
)
{
    Info<<"Checking for holes at refinement interfaces and refine "<< nl << endl;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const pointField& locationsInMesh = refineParams.locationsInMesh();
    const labelList& cellLevel = meshCutter_.cellLevel();

    label iter = 0;

    for (;iter < maxIter; iter++)
    {
        labelList ownPatch, nbrPatch;
        getBafflePatches
        (
            refineParams,
            globalToMasterPatch,
            ownPatch,
            nbrPatch
        );

        boolList blockedFace(mesh_.nFaces(),false);
        forAll(mesh_.faces(), facei)
        {
            if (ownPatch[facei] != -1)
            {
                blockedFace[facei] = true;
            }
        }

        // Set region per cell based on walking
        regionSplit cellRegion(mesh_, blockedFace);

        DynamicList<label> regionsToKeep(cellRegion.nRegions());

        // For all locationsInMesh find the cell
        forAll(locationsInMesh, i)
        {
            // Get location and index of zone ("none" for cellZone -1)
            const point& insidePoint = locationsInMesh[i];

            // Find the region containing the keepPoint
            label keepRegionI = -1;

            label celli = findCell
            (
                insidePoint,
                mesh_,
                meshCutter_
            );

            if (celli != -1)
            {
                keepRegionI = cellRegion[celli];
            }
            reduce(keepRegionI, maxOp<label>());
            regionsToKeep.append(keepRegionI);
        }
        regionsToKeep.shrink();
        labelHashSet regionsToKeepSet(regionsToKeep);

        boolList removedCells(mesh_.nCells(), false);
        boolList markedCells(mesh_.nCells(), false);
        forAll(mesh_.cells(), celli)
        {
            if (!regionsToKeepSet.found(cellRegion[celli]))
            {
                removedCells[celli] = true;
            }
        }

        // Calculate neighbour cell level
        scalarField neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        boolList neiRemoved(mesh_.nFaces()-mesh_.nInternalFaces(), false);
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
         )
        {
            neiLevel[facei-mesh_.nInternalFaces()] =
                cellLevel[mesh_.faceOwner()[facei]];
            neiRemoved[facei-mesh_.nInternalFaces()] =
                removedCells[mesh_.faceOwner()[facei]];
        }
        syncTools::swapBoundaryFaceList(mesh_, neiLevel);
        syncTools::swapBoundaryFaceList(mesh_, neiRemoved);

        boolList refineLowerCell(mesh_.nFaces(), false);

        forAll(mesh_.cells(), celli)
        {
            if (removedCells[celli])
            {
                const labelList& cFaces = mesh_.cells()[celli];
                label oLevel = cellLevel[celli];

                DynamicList<label> higherFaces(cFaces.size());
                DynamicList<label> lowerFaces(cFaces.size());

                forAll(cFaces, cFI)
                {

                    label nLevel = -1;
                    label nRemoved = false;
                    label facei = cFaces[cFI];
                    if (mesh_.isInternalFace(facei))
                    {
                        label own = mesh_.faceOwner()[facei];
                        label nei = (own == celli)
                            ? mesh_.faceNeighbour()[facei] : own;
                        nLevel = cellLevel[nei];
                        nRemoved = removedCells[nei];
                    }
                    else
                    {
                        label patchI = patches.whichPatch(facei);
                        if (patches[patchI].coupled())
                        {
                            nLevel = neiLevel[facei-mesh_.nInternalFaces()];
                            nRemoved = neiRemoved[facei-mesh_.nInternalFaces()];
                        }
                        else
                        {
                            nRemoved = true;
                        }
                    }

                    if (!nRemoved && nLevel != -1 && oLevel < nLevel)
                    {
                        higherFaces.append(facei);
                    }
                    if (!nRemoved && oLevel == nLevel)
                    {
                        lowerFaces.append(facei);
                    }
                }

                if (lowerFaces.size() > 0 && higherFaces.size() > 0)
                {
                    bool oppositeFaces = false;

                    forAll(lowerFaces, lFI)
                    {
                        label lFaceI = lowerFaces[lFI];
                        forAll(higherFaces, hFI)
                        {
                            label hFaceI = higherFaces[hFI];
                            vector hfv = mesh_.faceAreas()[hFaceI];
                            hfv /= mag(hfv);

                            vector lfv = mesh_.faceAreas()[lFaceI];
                            lfv /= mag(lfv);

                            if (mag(hfv & lfv)> 0.707)
                            {
                                oppositeFaces = true;
                                break;
                            }
                        }
                        if (oppositeFaces)
                        {
                            break;
                        }
                    }
                    if (oppositeFaces)
                    {
                        markedCells[celli] = true;
                    }
                }

                lowerFaces.clear();
                higherFaces.clear();
                forAll(cFaces, cFI)
                {

                    label nLevel = -1;
                    label nRemoved = false;
                    label facei = cFaces[cFI];
                    if (mesh_.isInternalFace(facei))
                    {
                        label own = mesh_.faceOwner()[facei];
                        label nei = (own == celli)
                            ? mesh_.faceNeighbour()[facei] : own;
                        nLevel = cellLevel[nei];
                        nRemoved = removedCells[nei];
                    }
                    else
                    {
                        label patchI = patches.whichPatch(facei);
                        if (patches[patchI].coupled())
                        {
                            nLevel = neiLevel[facei-mesh_.nInternalFaces()];
                            nRemoved = neiRemoved[facei-mesh_.nInternalFaces()];
                        }
                        else
                        {
                            nRemoved = true;
                        }
                    }

                    if (!nRemoved && nLevel != -1 && oLevel == nLevel)
                    {
                        higherFaces.append(facei);
                    }
                    if (!nRemoved && oLevel > nLevel)
                    {
                        lowerFaces.append(facei);
                    }
                }

                if (lowerFaces.size() > 0 && higherFaces.size() > 0)
                {
                    label oppositeFaces = -1;

                    forAll(lowerFaces, lFI)
                    {
                        label lFaceI = lowerFaces[lFI];
                        forAll(higherFaces, hFI)
                        {
                            label hFaceI = higherFaces[hFI];
                            vector hfv = mesh_.faceAreas()[hFaceI];
                            hfv /= mag(hfv);

                            vector lfv = mesh_.faceAreas()[lFaceI];
                            lfv /= mag(lfv);

                            if (mag(hfv & lfv)> 0.707)
                            {
                                oppositeFaces = lFaceI;
                                break;
                            }
                        }
                        if (oppositeFaces != -1)
                        {
                            break;
                        }
                    }
                    if (oppositeFaces != -1)
                    {
                        refineLowerCell[oppositeFaces] = true;
                    }
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh_,
            refineLowerCell,
            orEqOp<bool>()
        );

        forAll(mesh_.faces(), facei)
        {
            if (refineLowerCell[facei])
            {
                if (mesh_.isInternalFace(facei))
                {
                    label own = mesh_.faceOwner()[facei];
                    label nei = mesh_.faceNeighbour()[facei];
                    label oLevel = cellLevel[own];
                    label nLevel = cellLevel[nei];
                    if (oLevel > nLevel)
                    {
                        markedCells[nei] = true;
                    }
                    else
                    {
                        markedCells[own] = true;
                    }
                }
                else
                {
                    label patchI = patches.whichPatch(facei);
                    if (patches[patchI].coupled())
                    {
                        label own = mesh_.faceOwner()[facei];
                        label nLevel = neiLevel[facei-mesh_.nInternalFaces()];
                        label oLevel = cellLevel[own];
                        if (oLevel < nLevel)
                        {
                            markedCells[own] = true;
                        }
                    }
                }
            }
        }

        cellSet candidateCellSet
        (
            mesh_,
            "candidateCells",
            mesh_.nCells()/1000
        );

        forAll(markedCells, celli)
        {
            if (markedCells[celli])
            {
                candidateCellSet.insert(celli);
            }
        }

        labelList cellsToRefine
        (
            meshCutter().consistentRefinement
            (
                candidateCellSet.toc(),
                true
            )
        );
        label nRefined = returnReduce(cellsToRefine.size(), sumOp<label>());

        if (nRefined != 0)
        {
            Info<<"Refine "<< nRefined <<" cells to prevent a void forming"<<endl;
            refine(cellsToRefine, refineParams.locationsInMesh());
        }
        else
        {
            break;
        }
    }
}


bool Foam::meshRefinement::checkRefinedAndRemove
(
    const labelList& cellsToRefine,
    const bool checkPointCells,
    const bool unnamedOnly
)
{
    labelList testFaces = identity(mesh_.nFaces());//intersectedFaces();

    labelList globalRegion1;
    labelList globalRegion2;
    labelList checkSurfaces;

    if (unnamedOnly)
    {
       checkSurfaces =
           surfaceZonesInfo::getUnnamedSurfaces(surfaces_.surfZones());
    }
    else
    {
       checkSurfaces =
           surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
           (surfaces_.surfZones());
    }


    getIntersections
    (
        mesh_.cellCentres(),
        checkSurfaces,
        testFaces,
        globalRegion1,
        globalRegion2,
        false
    );

    if (checkPointCells)
    {
        keepBoundaryPointCells
        (
            testFaces,
            globalRegion1,
            globalRegion2
        );
    }

    boolList blockedFace(mesh_.nFaces(),false);
    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);

    labelHashSet keepRegionSet(keepLargestRegions(cellRegion));
    boolList markedCells(mesh_.nCells(), false);
    forAll(mesh_.cells(), celli)
    {
        if (!keepRegionSet.found(cellRegion[celli]))
        {
            markedCells[celli] = true;
        }
    }

    //See if any gap cells present and prevent removal
    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh_.lookupObject<volScalarField>("gapCells");
        forAll(gapCells, celli)
        {
            if (gapCells[celli] > SMALL)
            {
                markedCells[celli] = false;
            }
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    DynamicList<label> cellsToRemove(cellsToRefine.size());

    forAll(cellsToRefine, i)
    {
        label celli = cellsToRefine[i];

        if (markedCells[celli])
        {
           cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    if (returnReduce(cellsToRemove.size(), sumOp<label>()) == 0)
    {
        return false;
    }
    else
    {
        Info<<"Removing "<< returnReduce(cellsToRemove.size(), sumOp<label>())
            <<" cells not connected to inside location" <<endl;

        labelList ownPatch(mesh_.nFaces(), -1);
        forAll(cellsToRemove, i)
        {
            label celli = cellsToRemove[i];
            cell c = mesh_.cells()[celli];
            label patchI = 0;
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label fp = patches.whichPatch(facei);
                if (fp != -1 && !patches[fp].coupled())
                {
                    patchI = fp;
                    break;
                }
            }

            forAll(c, cFI)
            {
                label facei = c[cFI];
                ownPatch[facei] = patchI;
            }
        }

        syncTools::syncFaceList
        (
            mesh_,
            ownPatch,
            maxEqOp<label>()
        );

        // Remove cells
        removeCells cellRemover(mesh_);

        // Pick up patches for exposed faces
        labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
        labelList exposedPatchIDs(exposedFaces.size());

        forAll(exposedFaces, i)
        {
            label facei = exposedFaces[i];
            exposedPatchIDs[i] = ownPatch[facei];
        }

        doRemoveCells
        (
            cellsToRemove,
            exposedFaces,
            exposedPatchIDs,
            cellRemover
        );
        return true;
    }
}


void Foam::meshRefinement::moveWrongSidedCells
(
    const refinementParameters& refineParams,
    bool wrongSidedCorner
)
{
    labelList testFaces = identity(mesh_.nFaces());//intersectedFaces();

    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc);

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    {
        labelList minLevel;
        calcCellCellRays
        (
            mesh_.cellCentres(),
            neiCc,
            labelList(neiCc.size(), -1),
            testFaces,
            start,
            end,
            minLevel
        );
    }

    labelList globalRegion1;
    labelList globalRegion2;

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;

    const labelList testSurfaces
    (
        wrongSidedCorner ?
        surfaces_.cornerCellSurfaces() :
        surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
        (
            surfaces_.surfZones()
        )
    );

    surfaces_.findNearestIntersection
    (
        testSurfaces,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
     );

    boolList blockedFace(mesh_.nFaces(),false);
    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        if (surface1[facei] != -1 || surface2[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);

    labelHashSet keepRegionSet(keepLargestRegions(cellRegion));
    labelList cellToZone(mesh_.nCells(), -1);

    boolList markedCells(mesh_.nCells(), false);
    label nMarked = 0;
    forAll(mesh_.cells(), celli)
    {
        if (!keepRegionSet.found(cellRegion[celli]))
        {
            markedCells[celli] = true;
            nMarked++;
        }
    }

    //See if any gap cells present and maintain them
    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh_.lookupObject<volScalarField>("gapCells");
        forAll(gapCells, celli)
        {
            if (gapCells[celli] > -1 && markedCells[celli])
            {
                markedCells[celli] = false;
                nMarked--;
            }
        }
    }

    Info<<"Number of wrong sided cells found : "
        << returnReduce(nMarked, sumOp<label>()) <<endl;

    pointField disp(mesh_.nPoints(), vector::zero);
    labelList nDisp(mesh_.nPoints(), 0);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList boundaryPoints(mesh_.nPoints(), -1);

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            const polyPatch& pp = patches[patchI];
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                const face f = mesh_.faces()[facei];
                forAll(f,fp)
                {
                    boundaryPoints[f[fp]] = patchI;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    forAll(boundaryPoints, pti)
    {
        if (boundaryPoints[pti] != -1)
        {
            const labelList& pFaces = mesh_.pointFaces()[pti];
            forAll(pFaces, pfi)
            {
                label facei = pFaces[pfi];
                blockedFace[facei] = false;
            }
        }
    }
    syncTools::syncFaceList(mesh_, blockedFace, andEqOp<bool>());

    nMarked = 0;
    forAll(mesh_.cells(), celli)
    {
        if (!markedCells[celli])
        {
            continue;
        }

        const labelList& c = mesh_.cells()[celli];
        vector transVec = vector::zero;
        label nHits = 0;

        const labelList& cPts = mesh_.cellPoints()[celli];

        bool allBoundary = true;
        forAll(cPts, cPtI)
        {
            if (boundaryPoints[cPts[cPtI]] == -1)
            {
                allBoundary = false;
                break;
            }
        }

        forAll(c, cFI)
        {
            label facei = c[cFI];
            label patchI = patches.whichPatch(facei);

            if
            (
                patchI == -1
                || (patchI != -1 && isA<processorPolyPatch>(patches[patchI]))
                || allBoundary
            )
            {
                if
                (
                    blockedFace[facei]
                    ||
                    (
                        allBoundary
                        && (surface1[facei] != -1 && surface2[facei] != -1)
                        &&
                        (
                            patchI != -1 &&
                            !isA<processorPolyPatch>(patches[patchI])
                        )
                    )
                )
                {
                    const point& cc = mesh_.cellCentres()[celli];
                    if (allBoundary)
                    {
                        point mPt =
                            0.5*(hit1[facei].hitPoint()+hit2[facei].hitPoint());
                        transVec += (mPt-cc);
                    }
                    else
                    {
                        vector vec1 = hit1[facei].hitPoint() - cc;
                        vector vec2 = hit2[facei].hitPoint() - cc;
                        if (mag(vec1) > mag(vec2))
                        {
                            transVec += vec1;
                        }
                        else
                        {
                            transVec += vec2;
                        }
                    }
                    nHits++;
                }
            }
        }

        if (mag(transVec) > SMALL)
        {
            transVec /= nHits;
            nMarked++;
            const labelList& cPts =  mesh_.cellPoints()[celli];
            forAll(cPts, cPtI)
            {
                label pointI = cPts[cPtI];
                disp[pointI] += transVec;
                nDisp[pointI]++;
            }
        }
    }

    DynamicList<label> projFaces(mesh_.nFaces()/10);
    DynamicList<scalar> projGapLength(mesh_.nFaces()/10);

    //Project gap internal and boundary faces
    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh_.lookupObject<volScalarField>("gapCells");

        boolList internalFaces(mesh_.nFaces(), false);
        label nInternalFaces = 0;

        scalarField neiGapCells;
        syncTools::swapBoundaryCellList(mesh_, gapCells, neiGapCells);

        forAll(mesh_.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);

            if (patchi != -1 && !patches[patchi].coupled())
            {
                continue;
            }

            label own = mesh_.faceOwner()[facei];
            scalar gapOwn = gapCells[own];
            scalar gapNei = 0;

            if (patchi == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                gapNei = gapCells[nei];
            }
            else if (patches[patchi].coupled())
            {
                gapNei = neiGapCells[facei-mesh_.nInternalFaces()];
            }

            if ((gapOwn > -1) && (gapOwn == gapNei))
            {
                internalFaces[facei] = true;
                nInternalFaces++;
            }
        }

        syncTools::syncFaceList(mesh_, internalFaces, orEqOp<bool>());

        DynamicList<label> testFaces(nInternalFaces);
        forAll(internalFaces, facei)
        {
            if (internalFaces[facei])
            {
                testFaces.append(facei);
            }
        }

        const indirectPrimitivePatch internalPatch
        (
            IndirectList<face>
            (
                mesh_.faces(),
                testFaces
            ),
            mesh_.points()
        );

        if (debug)
        {
            simpleVTKWriter
            (
                internalPatch.localFaces(),
                internalPatch.localPoints()
            ).write("moveWrongSidedGaps.vtk");
        }

        const labelList meshEdges
        (
            internalPatch.meshEdges(mesh_.edges(), mesh_.pointEdges())
        );
        const labelList& meshPoints = internalPatch.meshPoints();

        labelList nInternalEdgeFaces(internalPatch.nEdges());

        forAll(internalPatch.edges(), edgei)
        {
            const labelList& eFaces = internalPatch.edgeFaces()[edgei];
            nInternalEdgeFaces[edgei] = eFaces.size();
        }
        syncTools::syncEdgeList
        (
            mesh_,
            meshEdges,
            nInternalEdgeFaces,
            plusEqOp<label>(),
            label(0)
        );

        boolList internalPatchBoundaryPts(meshPoints.size(), false);

        forAll(nInternalEdgeFaces, edgei)
        {
            if (nInternalEdgeFaces[edgei] == 1)
            {
                edge e = internalPatch.edges()[edgei];
                internalPatchBoundaryPts[e[0]] = true;
                internalPatchBoundaryPts[e[1]] = true;
            }
        }

        boolList internalPatchBoundaryFaces(internalPatch.size(), false);
        forAll(meshPoints, ptI)
        {
            if (internalPatchBoundaryPts[ptI])
            {
                const labelList& pFaces =  internalPatch.pointFaces()[ptI];
                forAll(pFaces, pFI)
                {
                    internalPatchBoundaryFaces[pFaces[pFI]] = true;
                }
            }
        }

        pointField start(internalPatch.size());
        pointField end(internalPatch.size());
        nInternalFaces = 0;
        const scalar edge0Len = meshCutter_.level0EdgeLength();

        forAll(internalPatch, i)
        {
            label facei = internalPatch.addressing()[i];
            point fC = mesh_.faceCentres()[facei];
            point fA = mesh_.faceAreas()[facei];
            fA /= (mag(fA) + SMALL);
            label faceLevel = meshCutter_.faceLevel(facei);
            scalar projDist = 1.5*edge0Len/(1<<faceLevel);

            start[i] = fC - projDist*fA;
            end[i] = fC + projDist*fA;
        }

        labelList surface1;
        List<pointIndexHit> hit1;
        labelList region1;
        labelList surface2;
        List<pointIndexHit> hit2;
        labelList region2;

        labelList thinGapSurfaces = surfaces_.thinGapSurfaces();

        surfaces_.findNearestIntersection
        (
            thinGapSurfaces,
            start,
            end,

            surface1,
            hit1,
            region1,

            surface2,
            hit2,
            region2
        );

        pointField normal1(hit1.size(), vector(GREAT, GREAT, GREAT));
        pointField normal2(hit2.size(), vector(GREAT, GREAT, GREAT));

        forAll(thinGapSurfaces, sI)
        {
            label surfI = thinGapSurfaces[sI];
            label geomI = surfaces_.surfaces()[surfI];

            //normal for hit point1
            {
                DynamicList<pointIndexHit> localHits;

                forAll(hit1, i)
                {
                    if (surface1[i] == surfI)
                    {
                        localHits.append(hit1[i]);
                    }
                }
                localHits.shrink();

                pointField localNormals;
                surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

                label localI = 0;
                forAll(hit1, i)
                {
                    if (surface1[i] == surfI)
                    {
                        normal1[i] = localNormals[localI];
                        localI++;
                    }
                }
            }

            //normal for hit point2
            {
                DynamicList<pointIndexHit> localHits;
                forAll(hit2, i)
                {
                    if (surface2[i] == surfI)
                    {
                        localHits.append(hit2[i]);
                    }
                }
                localHits.shrink();

                pointField localNormals;

                surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

                label localI = 0;
                forAll(hit2, i)
                {
                    if (surface2[i] == surfI)
                    {
                        normal2[i] = localNormals[localI];
                        localI++;
                    }
                }
            }
        }

        //Store average normals at local points
        pointField aveCentre(meshPoints.size(), vector::zero);
        List<DynamicList<point>> aveNormal(meshPoints.size());
        scalarField aveDist(meshPoints.size(), scalar(0));
        labelList nHits(meshPoints.size(), label(0));

        //Store average normals at local edges
        pointField aveEdgeCentre(meshEdges.size(), vector::zero);
        List<DynamicList<point>> aveEdgeNormal(meshEdges.size());
        scalarField aveEdgeDist(meshEdges.size(), scalar(0));
        labelList nEdgeHits(meshEdges.size(), label(0));

        forAll(internalPatch, i)
        {
            if
            (
                hit1[i].hit() && hit2[i].hit()
                && !internalPatchBoundaryFaces[i]
            )
            {
                vector separationVec = hit2[i].hitPoint()-hit1[i].hitPoint();
                scalar separationDist = 0.5*
                (
                    mag(separationVec & normal1[i]) +
                    mag(separationVec & normal2[i])
                );

                if (separationDist < 1e-8)
                {
                    continue;
                }

                point aveHitPoint =
                    0.5*(hit2[i].hitPoint() + hit1[i].hitPoint());

                const face lf = internalPatch.localFaces()[i];

                vector nDir = normal1[i];
                if ((nDir & normal2[i]) > 0)
                {
                    nDir += normal2[i];
                }
                else
                {
                    nDir -= normal2[i];
                }

                nDir /= (mag(nDir) + SMALL);
                forAll(lf, fp)
                {
                    label localpointi = lf[fp];
                    aveCentre[localpointi] += aveHitPoint;
                    aveNormal[localpointi].append(nDir);
                    aveDist[localpointi] += separationDist;
                    nHits[localpointi]++;
                }

                const labelList& lfe = internalPatch.faceEdges()[i];
                forAll(lfe, fe)
                {
                    label localedgei = lfe[fe];
                    aveEdgeCentre[localedgei] += aveHitPoint;
                    aveEdgeNormal[localedgei].append(nDir);
                    aveEdgeDist[localedgei] += separationDist;
                    nEdgeHits[localedgei]++;
                }
            }
        }

        forAll(internalPatch, i)
        {
            if
            (
                ( hit1[i].hit() && hit2[i].hit() )
                || internalPatchBoundaryFaces[i]
            )
            {
                nMarked++;
                label facei = internalPatch.addressing()[i];
                vector separationVec = vector::zero;
                scalar separationDist = 0;

                if (!internalPatchBoundaryFaces[i])
                {
                    separationVec = hit2[i].hitPoint()-hit1[i].hitPoint();
                    separationDist = 0.5*
                    (
                        mag(separationVec & normal1[i]) +
                        mag(separationVec & normal2[i])
                    );

                    if (separationDist < 1e-8)
                    {
                        continue;
                    }
                }

                const face lf = internalPatch.localFaces()[i];
                const face f = mesh_.faces()[facei];

                point aveHitPoint = vector::zero;
                vector aveNorm = vector::zero;
                scalar aveSepDist = 0;

                point fC = mesh_.faceCentres()[facei];

                if (internalPatchBoundaryFaces[i])
                {
                    label nEdgesSet = 0;
                    const labelList& lfe = internalPatch.faceEdges()[i];

                    forAll(lfe,fe)
                    {
                        label localedgei = lfe[fe];
                        if (nEdgeHits[localedgei] > 0)
                        {
                            nEdgesSet++;
                        }
                    }

                    if (nEdgesSet != 0)
                    {
                        label nSet = 0;
                        forAll(lfe,fe)
                        {
                            label localedgei = lfe[fe];
                            if (nEdgeHits[localedgei] > 0)
                            {
                                nSet++;

                                forAll(aveEdgeNormal[localedgei], aNI)
                                {
                                    vector nVec = aveEdgeNormal[localedgei][aNI];
                                    if ((aveNorm & nVec) > 0)
                                    {
                                        aveNorm += nVec;
                                    }
                                    else
                                    {
                                        aveNorm -= nVec;
                                    }
                                    aveHitPoint += aveEdgeCentre[localedgei];
                                    aveSepDist += aveEdgeDist[localedgei];
                                }
                            }
                        }

                        aveHitPoint /= nSet;
                        aveNorm /= nSet;
                        aveSepDist /= nSet;

                        plane fPlane(aveHitPoint, aveNorm);
                        aveHitPoint = fPlane.nearestPoint(0.5*(fC+aveHitPoint));
                    }
                    else
                    {
                        label nSet = 0;
                        forAll(lf,fp)
                        {
                            label localpointi = lf[fp];
                            if (nHits[localpointi] > 0)
                            {
                                nSet++;
                                vector pNormal = vector::zero;
                                forAll(aveNormal[localpointi], aNI)
                                {
                                    vector nVec = aveNormal[localpointi][aNI];
                                    if (aNI == 0)
                                    {
                                        pNormal = nVec;
                                    }
                                    else
                                    {
                                        if ((pNormal & nVec) > 0)
                                        {
                                            pNormal += nVec;
                                        }
                                        else
                                        {
                                            pNormal -= nVec;
                                        }
                                    }
                                }
                                pNormal /= nHits[localpointi];

                                if (nSet == 1)
                                {
                                    aveNorm += pNormal;
                                }
                                else
                                {
                                    if ((aveNorm & pNormal) > 0)
                                    {
                                        aveNorm += pNormal;
                                    }
                                    else
                                    {
                                        aveNorm -= pNormal;
                                    }
                                }
                                aveHitPoint += aveCentre[localpointi]
                                    / nHits[localpointi];
                                aveSepDist += aveDist[localpointi]
                                    / nHits[localpointi];
                            }
                        }

                        if (nSet > 0)
                        {
                            aveHitPoint /= nSet;
                            aveNorm /= nSet;
                            aveSepDist /= nSet;

                            plane fPlane(aveHitPoint, aveNorm);
                            aveHitPoint = fPlane.nearestPoint
                            (
                                0.5*(fC+aveHitPoint)
                            );
                        }
                        else
                        {
                            continue;
                        }
                    }
                }
                else
                {
                    aveHitPoint =
                        0.5*(hit2[i].hitPoint() + hit1[i].hitPoint());
                    aveNorm = normal1[i];

                    if ((aveNorm & normal2[i]) > 0)
                    {
                        aveNorm += normal2[i];
                    }
                    else
                    {
                        aveNorm -= normal2[i];
                    }

                    aveNorm /= (mag(aveNorm) + SMALL);
                    aveSepDist = separationDist;
                }

                Pair<vector> n1n2;
                n1n2[0] = mesh_.faceAreas()[facei];
                n1n2[1] = aveNorm;
                if ((n1n2[1] & n1n2[0]) > 0)
                {
                    n1n2[1] = -n1n2[1];
                }

                n1n2[0] /= (mag(n1n2[0]) + VSMALL);
                n1n2[1] /= (mag(n1n2[1]) + VSMALL);

                // flatten the face before moving
                vectorField flatFace(f.size(), vector::zero);

                forAll(f, fp)
                {
                    label pointi = f[fp];
                    point pt = mesh_.points()[pointi];
                    flatFace[fp] = pt;
                    vector pos = pt - fC;
                    flatFace[fp] -= (pos & n1n2[0])*n1n2[0];
                }

                label nPoints = f.size();
                // calculate new face centre and area of flattened face
                vector newFaceCentre, newFaceArea;

                // If the face is a triangle, do a direct calculation for
                // efficiency and to avoid round-off error-related problems
                if (nPoints == 3)
                {
                    newFaceCentre = (1.0/3.0)*(flatFace[0] + flatFace[1]
                                  + flatFace[2]);
                    newFaceArea = 0.5*((flatFace[1] - flatFace[0])^(flatFace[2]
                                - flatFace[0]));
                }
                else
                {
                    vector sumN = vector::zero;
                    scalar sumA = 0.0;
                    vector sumAc = vector::zero;
                    point fCentre = flatFace[0];
                    for (label pi = 1; pi < nPoints; pi++)
                    {
                        fCentre += flatFace[pi];
                    }
                    fCentre /= nPoints;

                    for (label pi = 0; pi < nPoints; pi++)
                    {
                        const point& nextPoint = flatFace[(pi + 1) % nPoints];
                        vector c = flatFace[pi] + nextPoint + fCentre;
                        vector n = (nextPoint - flatFace[pi]) ^
                            (fCentre - flatFace[pi]);
                        scalar a = mag(n);
                        sumN += n;
                        sumA += a;
                        sumAc += a*c;
                    }

                    newFaceCentre = (1.0/3.0)*sumAc/(sumA + VSMALL);
                    newFaceArea = 0.5*sumN;
                }

                n1n2[0] = newFaceArea;
                n1n2[1] = aveNorm;
                if ((n1n2[1] & n1n2[0]) > 0)
                {
                    n1n2[1] = -n1n2[1];
                }

                n1n2[0] /= (mag(n1n2[0]) + VSMALL);
                n1n2[1] /= (mag(n1n2[1]) + VSMALL);

                tensor T = rotationTensor(n1n2[0], -n1n2[1]);

                pointField movedInternalPts(f.size());

                forAll(f, fp)
                {
                    label pointi = f[fp];
                    point startPt = mesh_.points()[pointi];

                    // rotate points
                    point p = flatFace[fp];
                    p -= newFaceCentre;
                    p = transform(T, p);
                    p += newFaceCentre;

                    p += (aveHitPoint - newFaceCentre);

                    movedInternalPts[fp] = p;
                    disp[pointi] += (p-startPt);
                    nDisp[pointi]++;
                }

                projFaces.append(facei);
                projGapLength.append(aveSepDist);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        disp,
        plusEqOp<vector>(),
        vector::zero        // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        nDisp,
        plusEqOp<label>(),
        label(0)        // null value
    );

    if (returnReduce(nMarked, sumOp<label>()) != 0)
    {
        forAll(mesh_.points(), pointI)
        {
            if (nDisp[pointI] > 0)
            {
                disp[pointI] /= nDisp[pointI];
            }
        }

        pointField newPoints(mesh_.points());
        newPoints += disp;
        mesh_.movePoints(newPoints);
    }

    // Now move boundary gap cells
    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        disp = vector::zero;
        nDisp = 0;
        forAll(projFaces, i)
        {
            label facei = projFaces[i];
            const face f = mesh_.faces()[facei];

            scalar aveSepDist = projGapLength[i];

            label own = mesh_.faceOwner()[facei];
            const cell& cOwn = mesh_.cells()[own];
            labelHashSet fIntSet(f);

            forAll(cOwn, cFI)
            {
                label ownfacei = cOwn[cFI];
                label patchi = patches.whichPatch(ownfacei);

                if (patchi != -1 && !patches[patchi].coupled())
                {
                    const face fOwn = mesh_.faces()[ownfacei];
                    vector sNorm = mesh_.faceAreas()[facei];
                    if ((sNorm & mesh_.faceAreas()[ownfacei]) < 0)
                    {
                        sNorm = -sNorm;
                    }

                    sNorm /= (mag(sNorm) + VSMALL);
                    vector fProj = 0.5*aveSepDist*sNorm;

                    forAll(fOwn,fp)
                    {
                        label pointi = fOwn[fp];
                        const labelList& pEdges =
                            mesh_.pointEdges()[pointi];

                        label intPt = -1;

                        forAll(pEdges, pEI)
                        {
                            const edge e = mesh_.edges()[pEdges[pEI]];
                            label otherPt(e[0] == pointi ? e[1] : e[0]);
                            if (fIntSet.found(otherPt))
                            {
                                intPt = otherPt;
                                break;
                            }
                        }

                        if (intPt != -1)
                        {
                            label fIndex = findIndex(f, intPt);
                            point startPt = mesh_.points()[pointi];
                            disp[pointi] +=
                            (
                                (fProj+mesh_.points()[f[fIndex]])
                                -startPt
                            );
                            nDisp[pointi]++;
                        }
                    }
                }
            }

            if (mesh_.isInternalFace(facei))
            {
                label nei = mesh_.faceNeighbour()[facei];
                const cell& cNei = mesh_.cells()[nei];

                forAll(cNei, cFI)
                {
                    label neifacei = cNei[cFI];
                    label patchi = patches.whichPatch(neifacei);

                    if (patchi != -1 && !patches[patchi].coupled())
                    {
                        const face fNei = mesh_.faces()[neifacei];
                        vector sNorm = mesh_.faceAreas()[facei];
                        if ((sNorm & mesh_.faceAreas()[neifacei]) < 0)
                        {
                            sNorm = -sNorm;
                        }
                        sNorm /= (mag(sNorm) + VSMALL);
                        vector fProj = 0.5*aveSepDist*sNorm;

                        forAll(fNei,fp)
                        {
                            label pointi = fNei[fp];

                            const labelList& pEdges =
                                mesh_.pointEdges()[pointi];

                            label intPt = -1;

                            forAll(pEdges, pEI)
                            {
                                const edge e = mesh_.edges()[pEdges[pEI]];
                                label otherPt(e[0] == pointi ? e[1] : e[0]);
                                if (fIntSet.found(otherPt))
                                {
                                    intPt = otherPt;
                                    break;
                                }
                            }

                            if (intPt != -1)
                            {
                                label fIndex = findIndex(f, intPt);
                                point startPt = mesh_.points()[pointi];
                                disp[pointi] += (fProj+mesh_.points()[f[fIndex]])
                                    -startPt;
                                nDisp[pointi]++;
                            }
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            disp,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh_,
            nDisp,
            plusEqOp<label>(),
            label(0)        // null value
        );

        forAll(mesh_.points(), pointI)
        {
            if (nDisp[pointI] > 0)
            {
                disp[pointI] /= nDisp[pointI];
            }
        }

        pointField newPoints(mesh_.points());
        newPoints += disp;
        mesh_.movePoints(newPoints);
    }

    // Smooth gap cells at interface
    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        disp = vector::zero;
        nDisp = 0;

        const volScalarField& gapCells =
            mesh_.lookupObject<volScalarField>("gapCells");

        scalarField neiGapCells;
        syncTools::swapBoundaryCellList(mesh_, gapCells, neiGapCells);

        //charactrise points
        labelList gapPtsType(mesh_.nPoints(), -1);

        forAll(mesh_.cells(), celli)
        {
            if (gapCells[celli] > -1)
            {
                const labelList& cPts = mesh_.cellPoints()[celli];
                forAll(cPts, i)
                {
                    gapPtsType[cPts[i]] = 0;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            gapPtsType,
            maxEqOp<label>(),
            label(-1)       // null value
        );

        forAll(mesh_.points(), pointi)
        {
            if (gapPtsType[pointi] == 0)
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, i)
                {
                    label celli = pCells[i];
                    if (gapCells[celli] < 0)
                    {
                        const labelList& cPts = mesh_.cellPoints()[celli];
                        forAll(cPts, j)
                        {
                            if (gapPtsType[cPts[j]] == -1)
                            {
                                gapPtsType[cPts[j]] = 1;
                            }
                        }
                    }
                }
            }
        }

        forAll(mesh_.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);

            if (patchi != -1 && !patches[patchi].coupled())
            {
                const face f = mesh_.faces()[facei];
                forAll(f,fp)
                {
                    if (gapPtsType[f[fp]] == 1)
                    {
                        gapPtsType[f[fp]] = 2;
                    }
                }
            }
        }

        forAll(mesh_.points(), pointi)
        {
            if (gapPtsType[pointi] == 1)
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, i)
                {
                    label celli = pCells[i];
                    disp[pointi] += mesh_.cellCentres()[celli];
                    nDisp[pointi]++;
                }
            }
            else if (gapPtsType[pointi] == 2)
            {
                const labelList& pFaces = mesh_.pointFaces()[pointi];
                forAll(pFaces, i)
                {
                    label facei = pFaces[i];
                    label patchi = patches.whichPatch(facei);

                    if (patchi != -1 && !patches[patchi].coupled())
                    {
                        disp[pointi] += mesh_.faceCentres()[facei];
                        nDisp[pointi]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            disp,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh_,
            nDisp,
            plusEqOp<label>(),
            label(0)        // null value
        );

        pointField newPoints(mesh_.points());

        forAll(mesh_.points(), pointi)
        {
            if (nDisp[pointi] > 0)
            {
                newPoints[pointi] = disp[pointi]/nDisp[pointi];
            }
        }
        mesh_.movePoints(newPoints);
    }
}


void Foam::meshRefinement::removeOtherSide
(
    const refinementParameters& refineParams,
    const labelList& globalToMasterPatch,
    const mapPolyMesh& map
)
{
    labelList testFaces = identity(mesh_.nFaces());//intersectedFaces();

    const labelList& unamedAndBoundarySurfaces =
        surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces
        (
            surfaces_.surfZones()
        );

    labelList globalRegion1;
    labelList globalRegion2;
    getIntersections
    (
        mesh_.cellCentres(),
        unamedAndBoundarySurfaces,
        testFaces,
        globalRegion1,
        globalRegion2,
        false
    );

    boolList blockedFace(mesh_.nFaces(),false);
    forAll(testFaces, i)
    {
        label facei = testFaces[i];

        if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
        {
            blockedFace[facei] = true;
        }
    }

    // Set region per cell based on walking
    regionSplit cellRegion(mesh_, blockedFace);

    const pointField& locationsInMesh = refineParams.locationsInMesh();
    DynamicList<label> regionsToKeep(cellRegion.nRegions());

    // For all locationsInMesh find the cell
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];

        // Find the region containing the keepPoint
        label keepRegionI = -1;

        label celli = findCell
        (
            insidePoint,
            mesh_,
            meshCutter_
        );

        if (celli != -1)
        {
            keepRegionI = cellRegion[celli];
        }
        reduce(keepRegionI, maxOp<label>());
        regionsToKeep.append(keepRegionI);
    }
    regionsToKeep.shrink();
    labelHashSet regionsToKeepSet(regionsToKeep);

    boolList markedCells(mesh_.nCells(), false);
    forAll(mesh_.cells(), celli)
    {
        if (!regionsToKeepSet.found(cellRegion[celli]))
        {
            markedCells[celli] = true;
        }
    }

    boolList origCellsRemove(map.nOldCells(), false);
    const labelList& oldCellMap = map.cellMap();
    forAll(mesh_.cells(), celli)
    {
        if (markedCells[celli])
        {
            label oldCellI = oldCellMap[celli];
            if (oldCellI != -1)
            {
                origCellsRemove[oldCellI] = true;
            }
        }
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    forAll(mesh_.cells(), celli)
    {
        label origCellI = oldCellMap[celli];
        if (origCellI != -1 && origCellsRemove[origCellI])
        {
            const labelList& cFaces = mesh_.cells()[celli];
            bool containsBoundaryFace = false;
            forAll(cFaces, cFI)
            {
                label facei = cFaces[cFI];
                label patchI = patches.whichPatch(facei);

                if (patchI != -1 && !patches[patchI].coupled())
                {
                    containsBoundaryFace = true;
                    break;
                }
            }
            if (containsBoundaryFace)
            {
                markedCells[celli] = true;
            }
        }
    }

    DynamicList<label> cellsToRemove(mesh_.nCells()/100);
    forAll(mesh_.cells(), celli)
    {
        if (markedCells[celli])
        {
            cellsToRemove.append(celli);
        }
    }


    cellsToRemove.shrink();

    Info<<"Removing "<< returnReduce(cellsToRemove.size(), sumOp<label>())
        <<" cells not connected to inside location" <<endl;

    labelList ownPatch(mesh_.nFaces(), -1);
    forAll(cellsToRemove, i)
    {
        label celli = cellsToRemove[i];
        cell c = mesh_.cells()[celli];
        label patchI = 0;

        forAll(c, cFI)
        {
            label facei = c[cFI];
            if (globalRegion1[facei] != -1)
            {
                patchI = globalToMasterPatch[globalRegion1[facei]];
                break;
            }
        }

        forAll(c, cFI)
        {
            label facei = c[cFI];
            ownPatch[facei] = patchI;
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        ownPatch,
        maxEqOp<label>()
    );

    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatchIDs(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];

        exposedPatchIDs[i] = ownPatch[facei];
    }

    doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        cellRemover
    );

    return;
}


label Foam::meshRefinement::removeHoleCells()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    const labelList& pointLevel = meshCutter_.pointLevel();
    const labelList& cellLevel = meshCutter_.cellLevel();

    labelList markedCells(mesh_.nCells(), -1);
    labelList boundaryFaces(mesh_.nEdges(), -1);
    labelList boundaryEdges(mesh_.nEdges(), -1);
    labelList boundaryPoints(mesh_.nPoints(), -1);

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            const polyPatch& pp = patches[patchI];
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                boundaryFaces[facei] = patchI;
                const labelList& fEdges = mesh_.faceEdges()[facei];
                forAll(fEdges, fEI)
                {
                    boundaryEdges[fEdges[fEI]] = patchI;
                }
                const face f = mesh_.faces()[facei];
                forAll(f,fp)
                {
                    boundaryPoints[f[fp]] = patchI;
                }
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        boundaryFaces,
        maxEqOp<label>()
    );

    syncTools::syncEdgeList
    (
        mesh_,
        boundaryEdges,
        maxEqOp<label>(),
        label(-1)
    );

    syncTools::syncPointList
    (
        mesh_,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    while (true)
    {
        label nMarked = 0;
        forAll(mesh_.cells(), celli)
        {
            if (markedCells[celli] != -1)
            {
                continue;
            }

            const labelList& cFaces = mesh_.cells()[celli];
            const labelList& cEdges = mesh_.cellEdges()[celli];
            const labelList& cPoints = mesh_.cellPoints()[celli];

            DynamicList<label> cellBoundEdges(cEdges.size());
            forAll(cEdges, cEI)
            {
                label edgeI = cEdges[cEI];
                if (boundaryEdges[edgeI] > -1)
                {
                    cellBoundEdges.append(edgeI);
                }
            }

            if (cellBoundEdges.size() > 3)
            {
                labelHashSet bEdges(cellBoundEdges);
                bool allInternal = true;
                forAll(cFaces, cFI)
                {
                    if (boundaryFaces[cFaces[cFI]] > -1)
                    {
                        allInternal = false;
                        break;
                    }
                }
                if (allInternal)
                {
                    forAll(cPoints, i)
                    {
                        label pointI = cPoints[i];
                        label startEdge = -1;
                        label nBE = 0;
                        label cLevel = cellLevel[celli];
                        if (pointLevel[pointI] <= cLevel)
                        {
                            const labelList& pEdges = mesh_.pointEdges()[pointI];
                            forAll(pEdges, pEI)
                            {
                                label edgeI = pEdges[pEI];
                                if (bEdges.found(edgeI))
                                {
                                    nBE++;
                                    startEdge = edgeI;
                                }
                            }
                        }

                        label curEdge = startEdge;
                        label curPt = pointI;
                        if (nBE == 1)
                        {
                            label nStringEdges = 0;
                            bool closedLoop = false;
                            while (true)
                            {
                                edge e = mesh_.edges()[curEdge];
                                label nextPt(e[0] == curPt ? e[1] : e[0]);
                                if (pointLevel[nextPt] <= cLevel)
                                {
                                    nStringEdges += 2;
                                }
                                else
                                {
                                    nStringEdges++;
                                }

                                if (nextPt == pointI)
                                {
                                    closedLoop = true;
                                    break;
                                }

                                const labelList& pEdges =
                                    mesh_.pointEdges()[nextPt];
                                nBE  = 0;
                                label nextEdge = -1;
                                forAll(pEdges, pEI)
                                {
                                    label edgeI = pEdges[pEI];
                                    if (bEdges.found(edgeI))
                                    {
                                        nBE++;
                                        if (edgeI != curEdge)
                                        {
                                            nextEdge = edgeI;
                                        }
                                    }
                                }

                                if (nBE != 2)
                                {
                                    break;
                                }
                                else
                                {
                                    curPt = nextPt;
                                    curEdge = nextEdge;
                                }
                            }

                            if (nStringEdges == 8 && ! closedLoop)
                            {
                                label patchI = boundaryEdges[cellBoundEdges[0]];
                                markedCells[celli] = patchI;
                                nMarked++;
                                forAll(cFaces, cFI)
                                {
                                    boundaryFaces[cFaces[cFI]] = patchI;
                                }

                                forAll(cEdges, cEI)
                                {
                                    boundaryFaces[cEdges[cEI]] = patchI;
                                }

                                forAll(cPoints, cPI)
                                {
                                    boundaryFaces[cPoints[cPI]] = patchI;
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }

        if (returnReduce(nMarked, sumOp<label>()) == 0)
        {
            break;
        }

        syncTools::syncFaceList
        (
            mesh_,
            boundaryFaces,
            maxEqOp<label>()
        );

        syncTools::syncEdgeList
        (
            mesh_,
            boundaryEdges,
            maxEqOp<label>(),
            label(-1)
        );

        syncTools::syncPointList
        (
            mesh_,
            boundaryPoints,
            maxEqOp<label>(),
            label(-1)
        );
    }

    DynamicList<label> cellsToRemove(mesh_.nCells()/100);
    labelList ownPatch(mesh_.nFaces(), 0);
    forAll(mesh_.cells(), celli)
    {
        label patchI = markedCells[celli];
        if (patchI != -1)
        {
            const cell& c = mesh_.cells()[celli];

            forAll(c, cFI)
            {
                ownPatch[c[cFI]] = patchI;
            }
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label totalRemove = returnReduce(cellsToRemove.size(), sumOp<label>());

    Info<<"Removing " << totalRemove <<" problem hole cells" <<endl;

    if (totalRemove == 0)
    {
        return 0;
    }

    syncTools::syncFaceList
    (
        mesh_,
        ownPatch,
        maxEqOp<label>()
    );

    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatchIDs(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        exposedPatchIDs[i] = ownPatch[facei];
    }

    doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        cellRemover
    );

    return cellsToRemove.size();
}


label Foam::meshRefinement::removeProblemDualisationCells
(
    const labelHashSet& noLayerPatchIDs,
    const bool additionalChecks
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList boundaryEdges(mesh_.nEdges(), -1);
    labelList boundaryPoints(mesh_.nPoints(), -1);
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            const polyPatch& pp = patches[patchI];
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                const labelList& fEdges = mesh_.faceEdges()[facei];
                forAll(fEdges, fEI)
                {
                    boundaryEdges[fEdges[fEI]] = patchI;
                }
                const face f = mesh_.faces()[facei];
                forAll(f,fp)
                {
                    boundaryPoints[f[fp]] = patchI;
                }
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        boundaryEdges,
        maxEqOp<label>(),
        label(-1)
    );

    syncTools::syncPointList
    (
        mesh_,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    labelList markedCells(mesh_.nCells(), -1);

    const labelList& owners = mesh_.faceOwner();
    const labelList& neighbours = mesh_.faceNeighbour();

    labelList nBoundaryEdges(mesh_.nPoints(), 0);
    labelList nPointCells(mesh_.nPoints(), 0);
    labelList nPointEdges(mesh_.nPoints(), 0);

    labelList nInternalEdges(mesh_.nCells(), 0);
    labelList neiInternalEdges(mesh_.nFaces()-mesh_.nInternalFaces());

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh_));

    forAll(mesh_.points(), pointI)
    {
        const labelList& pEdges = mesh_.pointEdges()[pointI];
        forAll(pEdges, pEI)
        {
            label meshEdgeI = pEdges[pEI];

            if (isMasterEdge[meshEdgeI])
            {
                nPointEdges[pointI]++;
            }
        }
    }
    syncTools::syncPointList
    (
        mesh_,
        nPointEdges,
        plusEqOp<label>(),
        label(0)
    );

    while (true)
    {
        label nMarked = 0;

        neiInternalEdges = 0;
        nBoundaryEdges = 0;
        nPointCells = 0;

        forAll(mesh_.points(), pointI)
        {
            const labelList& pEdges = mesh_.pointEdges()[pointI];
            forAll(pEdges, pEI)
            {
                label meshEdgeI = pEdges[pEI];

                if (isMasterEdge[meshEdgeI] && boundaryEdges[meshEdgeI] != -1)
                {
                    nBoundaryEdges[pointI]++;
                }
            }

            const labelList& pCells = mesh_.pointCells()[pointI];
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                if (markedCells[celli] == -1)
                {
                    nPointCells[pointI]++;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh_,
            nBoundaryEdges,
            plusEqOp<label>(),
            label(0)
        );
        syncTools::syncPointList
        (
            mesh_,
            nPointCells,
            plusEqOp<label>(),
            label(0)
        );

        forAll(mesh_.points(), pointI)
        {
            if (boundaryPoints[pointI] == -1)
            {
                continue;
            }

            if
            (
                nBoundaryEdges[pointI] == nPointEdges[pointI]
                && nPointCells[pointI] == 6
            )
            {
                const labelList& pCells = mesh_.pointCells()[pointI];
                label patchI = boundaryPoints[pointI];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];
                    if (markedCells[celli] == -1)
                    {
                        markedCells[celli] = patchI;
                        const labelList& c = mesh_.cells()[celli];
                        forAll(c, cFI)
                        {
                            label facei = c[cFI];
                            const labelList& fEdges = mesh_.faceEdges()[facei];
                            forAll(fEdges, fEI)
                            {
                                boundaryEdges[fEdges[fEI]] = patchI;
                            }
                            const face f = mesh_.faces()[facei];
                            forAll(f,fp)
                            {
                                boundaryPoints[f[fp]] = patchI;
                            }
                        }
                        nMarked++;
                    }
                }
            }
        }

        nInternalEdges = 0;

        forAll(mesh_.cells(), celli)
        {
            const labelList& cEdges = mesh_.cellEdges()[celli];
            forAll(cEdges, cEI)
            {
                label edgeI = cEdges[cEI];
                if (boundaryEdges[edgeI] == -1)
                {
                    nInternalEdges[celli]++;
                }
            }
        }

        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            label bFaceI = facei-mesh_.nInternalFaces();
            neiInternalEdges[bFaceI] = nInternalEdges[owners[facei]];
        }
        syncTools::swapBoundaryFaceList(mesh_, neiInternalEdges);


        forAll(mesh_.cells(), celli)
        {
            if (markedCells[celli] == -1)
            {
                const labelList& cPoints = mesh_.cellPoints()[celli];
//                label nBoundaryPts = 0;
                label patchI = 0;
                forAll(cPoints, cPI)
                {
                    if (boundaryPoints[cPoints[cPI]] != -1)
                    {
//                        nBoundaryPts++;
                        patchI = boundaryPoints[cPoints[cPI]];
                    }
                }

                if (nInternalEdges[celli] <= 1)
                {
                    const cell& c = mesh_.cells()[celli];
                    label nNei = 0;

                    forAll(c, cFI)
                    {
                        label meshFaceI = c[cFI];
                        label own = owners[meshFaceI];
                        label patchI = patches.whichPatch(meshFaceI);
                        if (patchI == -1)
                        {
                            label nei = neighbours[meshFaceI];

                            label otherSide(own == celli ? nei : own);
                            if (nInternalEdges[otherSide] <= 2)
                            {
                                nNei++;
                            }
                        }
                        else if (patches[patchI].coupled())
                        {
                            label bFaceI = meshFaceI-mesh_.nInternalFaces();

                            if (neiInternalEdges[bFaceI] <= 2)
                            {
                                nNei++;
                            }
                        }
                    }

                    if (nInternalEdges[celli] == 0 || nNei == 0)
                    {
                        markedCells[celli] = patchI;
                        const labelList& c = mesh_.cells()[celli];
                        forAll(c, cFI)
                        {
                            label facei = c[cFI];
                            const labelList& fEdges = mesh_.faceEdges()[facei];
                            forAll(fEdges, fEI)
                            {
                                boundaryEdges[fEdges[fEI]] = patchI;
                            }
                            const face f = mesh_.faces()[facei];
                            forAll(f,fp)
                            {
                                boundaryPoints[f[fp]] = patchI;
                            }
                        }
                        nMarked++;
                   }
                }
            }
        }

        //Check for internal points connected to two cells where
        // all other points are boundary
        if (additionalChecks)
        {
            forAll(mesh_.cells(), celli)
            {
                if (markedCells[celli] != -1)
                {
                    continue;
                }

                bool allBoundary = true;
                const labelList& cPoints = mesh_.cellPoints()[celli];
                label patchI = -1;
                forAll(cPoints, cptI)
                {
                    label pointI = cPoints[cptI];
                    if (boundaryPoints[pointI] == -1)
                    {
                        allBoundary = false;
                        break;
                    }
                    else
                    {
                        patchI = boundaryPoints[pointI];
                    }
                }

                if (patchI != -1 && allBoundary && nInternalEdges[celli] < 3)
                {
                    markedCells[celli] = patchI;
                    const labelList& c = mesh_.cells()[celli];
                    forAll(c, cFI)
                    {
                        label facei = c[cFI];
                        const labelList& fEdges = mesh_.faceEdges()[facei];
                        forAll(fEdges, fEI)
                        {
                            boundaryEdges[fEdges[fEI]] = patchI;
                        }
                        const face f = mesh_.faces()[facei];
                        forAll(f,fp)
                        {
                            boundaryPoints[f[fp]] = patchI;
                        }
                    }
                    nMarked++;
                }
            }

            labelList nSingleInternalPts(mesh_.nPoints(), 0);
            forAll(mesh_.points(), pointI)
            {
                const labelList& pCells = mesh_.pointCells()[pointI];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];
                    const labelList& cPoints = mesh_.cellPoints()[celli];
                    label nBoundaryPts = 0;
                    forAll(cPoints, cPI)
                    {
                        if (boundaryPoints[cPoints[cPI]] != -1)
                        {
                            nBoundaryPts++;
                        }
                    }

                    if (nBoundaryPts == cPoints.size() -1)
                    {
                        nSingleInternalPts[pointI]++;
                    }
                }
            }
            syncTools::syncPointList
            (
                mesh_,
                nSingleInternalPts,
                plusEqOp<label>(),
                label(0)
            );

            forAll(mesh_.points(), pointI)
            {
                if
                (
                    (nSingleInternalPts[pointI] == 2
                     && boundaryPoints[pointI] == -1)
                    ||
                    (nSingleInternalPts[pointI] == 4
                     && boundaryPoints[pointI] == -1)
                )
                {
                    const labelList& pCells = mesh_.pointCells()[pointI];

                    DynamicList<label> singleIntCells(pCells.size());
                    DynamicList<label> singleIntPatchID(pCells.size());

                    bool grownUpPatch = false;
                    forAll(pCells, pCI)
                    {
                        label celli = pCells[pCI];
                        const labelList& cPoints = mesh_.cellPoints()[celli];
                        label nBoundaryPts = 0;
                        label patchI = 0;
                        forAll(cPoints, cPI)
                        {
                            if (boundaryPoints[cPoints[cPI]] != -1)
                            {
                                nBoundaryPts++;
                                patchI = boundaryPoints[cPoints[cPI]];

                                if (noLayerPatchIDs.found(patchI))
                                {
                                    grownUpPatch = true;
                                }
                            }
                        }

                        if (nBoundaryPts == cPoints.size() -1)
                        {
                            singleIntCells.append(celli);
                            singleIntPatchID.append(patchI);
                        }
                    }
                    singleIntCells.shrink();
                    singleIntPatchID.shrink();

                    //Only remove cells if grown up patches detected
                    if (!grownUpPatch)
                    {
                        continue;
                    }

                    forAll(singleIntCells, pCI)
                    {
                        label celli = singleIntCells[pCI];
                        if (markedCells[celli] == -1)
                        {
                            label patchI = singleIntPatchID[pCI];
                            markedCells[celli] = patchI;
                            const labelList& c = mesh_.cells()[celli];
                            forAll(c, cFI)
                            {
                                label facei = c[cFI];
                                const labelList& fEdges =
                                    mesh_.faceEdges()[facei];
                                forAll(fEdges, fEI)
                                {
                                    boundaryEdges[fEdges[fEI]] = patchI;
                                }
                                const face f = mesh_.faces()[facei];
                                forAll(f,fp)
                                {
                                    boundaryPoints[f[fp]] = patchI;
                                }
                            }
                            nMarked++;
                        }
                    }
                }
            }
        }

        if (returnReduce(nMarked, sumOp<label>()) == 0)
        {
            break;
        }

        syncTools::syncEdgeList
        (
            mesh_,
            boundaryEdges,
            maxEqOp<label>(),
            label(-1)
         );

        syncTools::syncPointList
        (
            mesh_,
            boundaryPoints,
            maxEqOp<label>(),
            label(-1)
        );
    }

    //See if any gap cells present and prevent removal
    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh_.lookupObject<volScalarField>("gapCells");
        forAll(gapCells, celli)
        {
            if (gapCells[celli] > -1)
            {
                markedCells[celli] = -1;
            }
        }
    }

    DynamicList<label> cellsToRemove(mesh_.nCells()/100);
    labelList ownPatch(mesh_.nFaces(), 0);
    forAll(mesh_.cells(), celli)
    {
        label patchI = markedCells[celli];
        if (patchI != -1)
        {
            const cell& c = mesh_.cells()[celli];

            forAll(c, cFI)
            {
                ownPatch[c[cFI]] = patchI;
            }
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label totalRemove = returnReduce(cellsToRemove.size(), sumOp<label>());

    Info<<"Removing " << totalRemove <<" problem dualisation cells" <<endl;

    if (totalRemove == 0)
    {
        return 0;
    }

    syncTools::syncFaceList
    (
        mesh_,
        ownPatch,
        maxEqOp<label>()
    );

    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatchIDs(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        exposedPatchIDs[i] = ownPatch[facei];
    }

    doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        cellRemover
    );

    return cellsToRemove.size();
}


label Foam::meshRefinement::baffleHoles
(
    label closeHoleSize
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    bool dual(controller_.algorithm() == meshControl::DUAL);

    labelList boundaryEdges(mesh_.nEdges(), -1);
    labelList boundaryPoints(mesh_.nPoints(), -1);
    labelList baffleFaces(mesh_.nFaces(), -1);

    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled())
        {
            const polyPatch& pp = patches[patchI];
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                baffleFaces[facei] = patchI;
                const labelList& fEdges = mesh_.faceEdges()[facei];
                forAll(fEdges, fEI)
                {
                    boundaryEdges[fEdges[fEI]] = patchI;
                }
                const face f = mesh_.faces()[facei];
                forAll(f,fp)
                {
                    boundaryPoints[f[fp]] = patchI;
                }
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        boundaryEdges,
        maxEqOp<label>(),
        label(-1)
    );

    syncTools::syncPointList
    (
        mesh_,
        boundaryPoints,
        maxEqOp<label>(),
        label(-1)
    );

    // Get per face whether is it master (of a coupled set of faces)
    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    vector axisA = Zero;
    vector axisB = Zero;

    const labelList& pointLevel = meshCutter_.pointLevel();
    const labelList& cellLevel = meshCutter_.cellLevel();

    if (Pstream::master())
    {
        if (mesh_.nCells() > 0)
        {
            label celli = 0;
            label cLevel = cellLevel[celli];

            labelHashSet cFaces(mesh_.cells()[celli]);
            const labelList& cPoints = mesh_.cellPoints()[celli];

            label pointI = -1;
            forAll(cPoints, cPtI)
            {
                pointI = cPoints[cPtI];
                if (pointLevel[pointI] <= cLevel)
                {
                    break;
                }
            }

            const labelList& pFaces = mesh_.pointFaces()[pointI];
            label nFound = 0;

            forAll(pFaces, pFI)
            {
                label facei = pFaces[pFI];
                if (cFaces.found(facei))
                {
                    vector fA = mesh_.faceAreas()[facei];
                    fA /= (mag(fA) + SMALL);

                    if (nFound == 0)
                    {
                        axisA = fA;
                    }
                    else if (nFound == 1)
                    {
                        axisB = fA;
                        break;
                    }
                    nFound++;
                }
            }
        }
        else
        {
            axisA = vector(1, 0, 0);
            axisB = vector(0, 1, 0);
        }
    }

    reduce(
        std::tie(axisA, axisB),
        ParallelOp<maxMagSqrOp<vector>, maxMagSqrOp<vector>>{}
    );

    List<direction> faceDir(mesh_.nFaces(),direction(0));

    forAll(mesh_.faces(), facei)
    {
        vector fA = mesh_.faceAreas()[facei];
        fA /= (mag(fA) + SMALL);
        if (mag(axisA & fA) > 0.707)
        {
            faceDir[facei] = direction(0);
        }
        else if (mag(axisB & fA) > 0.707)
        {
            faceDir[facei] = direction(1);
        }
        else
        {
            faceDir[facei] = direction(2);
        }
    }


    //characterise number of face boundary edges:
    labelList faceType(mesh_.nFaces(), 0);

    labelList nDir0(mesh_.nPoints(), 0);
    labelList nDir1(mesh_.nPoints(), 0);
    labelList nDir2(mesh_.nPoints(), 0);

    List<labelVector> nft3dir(mesh_.nEdges(), Zero);
    List<labelVector> nnbrdir(mesh_.nEdges(), Zero);

    while (true)
    {
        label nMarked = 0;
        faceType = 0;
        forAll(mesh_.faces(), facei)
        {
            const labelList& fEdges = mesh_.faceEdges()[facei];

            //Don't handle split faces
            if (fEdges.size() != 4 || baffleFaces[facei] != -1)
            {
                continue;
            }

            label nBoundaryEdge = 0;

            forAll(fEdges, fEI)
            {
                label edgeI = fEdges[fEI];
                if (boundaryEdges[edgeI] != -1)
                {
                    nBoundaryEdge++;
                }
            }

            faceType[facei] = nBoundaryEdge;
        }

        nDir0 = 0;
        nDir1 = 0;
        nDir2 = 0;

        forAll(mesh_.points(), pointI)
        {
            if (boundaryPoints[pointI] == -1)
            {
                const labelList& pFaces = mesh_.pointFaces()[pointI];

                forAll(pFaces, pFI)
                {
                    label facei = pFaces[pFI];

                    if (isMasterFace[facei] && faceType[facei] == 2)
                    {
                        if (faceDir[facei] == 0)
                        {
                            nDir0[pointI]++;
                        }
                        else if (faceDir[facei] == 1)
                        {
                            nDir1[pointI]++;
                        }
                        else
                        {
                            nDir2[pointI]++;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            nDir0,
            plusEqOp<label>(),
            label(0)
        );

        syncTools::syncPointList
        (
            mesh_,
            nDir1,
            plusEqOp<label>(),
            label(0)
        );

        syncTools::syncPointList
        (
            mesh_,
            nDir2,
            plusEqOp<label>(),
            label(0)
        );

        //Remove 4 face holes
        if (closeHoleSize >= 4)
        {
            forAll(mesh_.points(), pointI)
            {
                DynamicList<label> fillDir(3);

                if (nDir0[pointI] == 4)
                {
                    fillDir.append(0);
                }

                if (nDir1[pointI] == 4)
                {
                    fillDir.append(1);
                }

                if (nDir2[pointI] == 4)
                {
                    fillDir.append(2);
                }

                if (fillDir.size())
                {
                    const labelList& pFaces = mesh_.pointFaces()[pointI];
                    labelHashSet fillDirSet(fillDir);

                    forAll(pFaces, pFI)
                    {
                        label facei = pFaces[pFI];
                        if
                        (
                            baffleFaces[facei] == -1
                            && fillDirSet.found(faceDir[facei])
                        )
                        {
                            label patchI = 0;
                            const face f = mesh_.faces()[facei];
                            forAll(f,fp)
                            {
                                if (boundaryPoints[f[fp]] != -1)
                                {
                                    patchI = boundaryPoints[f[fp]];
                                    break;
                                }
                            }

                            baffleFaces[facei] = patchI;

                            const labelList& fEdges = mesh_.faceEdges()[facei];
                            forAll(fEdges, fEI)
                            {
                                boundaryEdges[fEdges[fEI]] = patchI;
                            }

                            forAll(f,fp)
                            {
                                boundaryPoints[f[fp]] = patchI;
                            }
                            nMarked++;
                        }
                    }
                }
            }
        }

        nnbrdir = Zero;

        //Remove 6 face holes
        if (closeHoleSize >= 6)
        {
            forAll(mesh_.edges(), edgeI)
            {
                edge e = mesh_.edges()[edgeI];
                label v0 = e[0];
                label v1 = e[1];
                if (boundaryPoints[v0] == -1 && boundaryPoints[v1] == -1)
                {
                    if (nDir0[v0] == 2 && nDir0[v1] == 2)
                    {
                        nnbrdir[edgeI] += labelVector(1, 0, 0);
                    }
                    if (nDir1[v0] == 2 && nDir1[v1] == 2)
                    {
                        nnbrdir[edgeI] += labelVector(0, 1, 0);
                    }
                    if (nDir2[v0] == 2 && nDir2[v1] == 2)
                    {
                        nnbrdir[edgeI] += labelVector(0, 0, 1);
                    }
                }
            }

            syncTools::syncEdgeList
            (
                mesh_,
                nnbrdir,
                plusEqOp<labelVector>(),  // in-place add
                labelVector(Zero),            // initial value
                dummyTransform()
            );

            forAll(mesh_.edges(), edgeI)
            {
                DynamicList<label> fillDir(3);

                if (nnbrdir[edgeI].x() > 0)
                {
                    fillDir.append(0);
                }

                if (nnbrdir[edgeI].y() > 0)
                {
                    fillDir.append(1);
                }

                if (nnbrdir[edgeI].z() > 0)
                {
                    fillDir.append(2);
                }

                if (fillDir.size() > 0)
                {
                    const labelList& eFaces = mesh_.edgeFaces()[edgeI];
                    labelHashSet fillDirSet(fillDir);

                    forAll(eFaces, eFI)
                    {
                        label facei = eFaces[eFI];

                        if
                        (
                            baffleFaces[facei] == -1
                            && fillDirSet.found(faceDir[facei])
                        )
                        {
                            label patchI = 0;
                            const face f = mesh_.faces()[facei];
                            forAll(f,fp)
                            {
                                if (boundaryPoints[f[fp]] != -1)
                                {
                                    patchI = boundaryPoints[f[fp]];
                                    break;
                                }
                            }

                            baffleFaces[facei] = patchI;

                            const labelList& fEdges = mesh_.faceEdges()[facei];
                            forAll(fEdges, fEI)
                            {
                                boundaryEdges[fEdges[fEI]] = patchI;
                            }

                            forAll(f,fp)
                            {
                                boundaryPoints[f[fp]] = patchI;
                            }
                            nMarked++;
                        }
                    }
                }
            }
        }

        nft3dir = Zero;
        forAll(mesh_.edges(), edgeI)
        {
            if (boundaryEdges[edgeI] != -1)
            {
                continue;
            }

            const labelList& eFaces = mesh_.edgeFaces()[edgeI];
            forAll(eFaces, eFI)
            {
                label facei = eFaces[eFI];
                if (isMasterFace[facei] && faceType[facei] == 3)
                {
                    direction dir = faceDir[facei];
                    nft3dir[edgeI][dir]++;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            nft3dir,
            plusEqOp<labelVector>(),  // in-place add
            labelVector(Zero),            // initial value
            dummyTransform()
        );

        //Remove 2 and 3 face holes
        if (closeHoleSize >= 2)
        {
            forAll(mesh_.faces(), facei)
            {
                if
                (
                    baffleFaces[facei] == -1
                    && (faceType[facei] == 2 || faceType[facei] == 3)
                )
                {
                    const labelList& fEdges = mesh_.faceEdges()[facei];
                    labelVector countNbrs(Zero);

                    forAll(fEdges, fEI)
                    {
                        label edgeI = fEdges[fEI];
                        if (boundaryEdges[edgeI] == -1)
                        {
                            countNbrs += nft3dir[edgeI];
                        }
                    }

                    if
                    (
                        countNbrs[faceDir[facei]] == 2
                    )
                    {
                        label patchI = 0;
                        const face f = mesh_.faces()[facei];
                        forAll(f,fp)
                        {
                            if (boundaryPoints[f[fp]] != -1)
                            {
                                patchI = boundaryPoints[f[fp]];
                                break;
                            }
                        }

                        baffleFaces[facei] = patchI;

                        const labelList& fEdges = mesh_.faceEdges()[facei];
                        forAll(fEdges, fEI)
                        {
                            boundaryEdges[fEdges[fEI]] = patchI;
                        }

                        forAll(f,fp)
                        {
                            boundaryPoints[f[fp]] = patchI;
                        }
                        nMarked++;
                    }
                }
            }
        }

        //Remove single face holes
        if (closeHoleSize >= 1)
        {
            forAll(mesh_.faces(), facei)
            {
                if (baffleFaces[facei] != -1)
                {
                    continue;
                }

                const labelList& fEdges = mesh_.faceEdges()[facei];
                label nBoundary = 0;

                forAll(fEdges, fEI)
                {
                    if (boundaryEdges[fEdges[fEI]] != -1)
                    {
                        nBoundary++;
                    }
                }

                if
                (
                    nBoundary == fEdges.size() ||
                    (dual && closeHoleSize > 1 && nBoundary == fEdges.size()-1)
                )
                {
                    label patchI = 0;
                    const face f = mesh_.faces()[facei];
                    forAll(f,fp)
                    {
                        if (boundaryPoints[f[fp]] != -1)
                        {
                            patchI = boundaryPoints[f[fp]];
                            break;
                        }
                    }

                    baffleFaces[facei] = patchI;

                    const labelList& fEdges = mesh_.faceEdges()[facei];
                    forAll(fEdges, fEI)
                    {
                        boundaryEdges[fEdges[fEI]] = patchI;
                    }

                    forAll(f,fp)
                    {
                        boundaryPoints[f[fp]] = patchI;
                    }
                    nMarked++;
                }
            }
        }

        if (returnReduce(nMarked, sumOp<label>()) == 0)
        {
            break;
        }

        syncTools::syncEdgeList
        (
            mesh_,
            boundaryEdges,
            maxEqOp<label>(),
            label(-1)
         );

        syncTools::syncPointList
        (
            mesh_,
            boundaryPoints,
            maxEqOp<label>(),
            label(-1)
        );

    }


    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh_.lookupObject<volScalarField>("gapCells");

        scalarField neiGapCells;
        syncTools::swapBoundaryCellList(mesh_, gapCells, neiGapCells);
        forAll(mesh_.faces(), facei)
        {
            if (baffleFaces[facei] != -1)
            {
                label patchi = patches.whichPatch(facei);
                label own = mesh_.faceOwner()[facei];
                scalar gapOwn = gapCells[own];
                scalar gapNei = 0;

                if (patchi == -1)
                {
                    label nei = mesh_.faceNeighbour()[facei];
                    gapNei = gapCells[nei];
                }
                else if (patches[patchi].coupled())
                {
                    gapNei = neiGapCells[facei-mesh_.nInternalFaces()];
                }

                if (gapOwn > SMALL || gapNei > SMALL)
                {
                    baffleFaces[facei] = -1;
                }
            }
        }
    }

    // Topochange container
    polyTopoChange meshMod(mesh_);

    label nBaffles = 0;
    forAll(mesh_.faces(), facei)
    {
        label patchI = patches.whichPatch(facei);

        if (baffleFaces[facei] != -1)
        {
            if (patchI == -1 || patches[patchI].coupled())
            {
                // Create baffle or repatch face. Return label of
                // inserted baffle face.
                label bafflePatchI = baffleFaces[facei];

                createBaffle
                (
                    facei,
                    bafflePatchI,   // owner side patch
                    bafflePatchI,   // neighbour side patch
                    meshMod
                );
                nBaffles++;
            }
        }
    }

    reduce(nBaffles, sumOp<label>());

    Info<<"Closing "<<nBaffles<<" hole faces"<<endl;

    if (nBaffles > 0)
    {
        // Change the mesh (no inflation, parallel sync)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes.
            mesh_.clearOut();
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());

        // None of the faces has changed, only the zones. Still...
        updateMesh(map, labelList());
    }

    return nBaffles;
}


Foam::labelList Foam::meshRefinement::keepLargestRegions
(
    const regionSplit& cellRegion
)
{
    const labelIOList& masterRegionID = cellRegionID();

    labelList newRegionToMaster(cellRegion.nRegions(), -1);
    labelList cellsInRegion(cellRegion.nRegions(), 0);
    label nRegions = -1;

    forAll(cellRegion, celli)
    {
        label regionI = cellRegion[celli];
        newRegionToMaster[regionI] = masterRegionID[celli];
        cellsInRegion[regionI]++;
        nRegions = max(nRegions, masterRegionID[celli]+1);
    }

    Pstream::listCombineReduce(newRegionToMaster, maxOp<label>());
    Pstream::listCombineReduce(cellsInRegion, plusOp<label>());
    reduce(nRegions, maxOp<label>());

    labelList maxRegionSize(nRegions,-1);
    labelList keepRegions(nRegions,-1);
    forAll(newRegionToMaster, regionI)
    {
        label masterRegionI = newRegionToMaster[regionI];
        if (cellsInRegion[regionI] > maxRegionSize[masterRegionI])
        {
            maxRegionSize[masterRegionI] = cellsInRegion[regionI];
            keepRegions[masterRegionI] = regionI;
        }
    }

    if (debug)
    {
        const_cast<Time&>(mesh_.time())++;
        mesh_.write();

        volScalarField volCellToRegion
        (
            IOobject
            (
                "cellToRegion",
                timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(mesh_.cells(), celli)
        {
            volCellToRegion[celli] = cellRegion[celli];
        }
        volCellToRegion.write();
        Info<<"Regions to keep "<<keepRegions<<endl;
    }

    return  keepRegions;
}


void Foam::meshRefinement::updateMasterRegions()
{
    regionSplit cellRegion(mesh_);
    updateRegionID(cellRegion);
}


label Foam::meshRefinement::removeZoneHoles(const bool external)
{
    Info<<"Checking internal zones for holes"<<endl;

    //First merge external cells
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const cellZoneMesh& cellZones = mesh_.cellZones();
    const faceZoneMesh& faceZones = mesh_.faceZones();

    polyTopoChange meshMod(mesh_);
    labelList cellZoneID(mesh_.nCells(), -1);

    //Check for named surfaces
    const PtrList<surfaceZonesInfo>& surfZones = surfaces().surfZones();

    DynamicList<label> internalNamedCellZones(cellZones.size());
    DynamicList<label> internalNamedFaceZones(faceZones.size());
    forAll(surfZones, surfI)
    {
        const word& faceZoneName = surfZones[surfI].faceZoneName();

        if (faceZoneName.size())
        {
            const surfaceZonesInfo::faceZoneType& faceType =
                surfZones[surfI].faceType();
            if
            (
                !surfZones[surfI].freeStanding()
                && faceType != surfaceZonesInfo::BOUNDARY
            )
            {
                const word& cellZoneName = surfZones[surfI].cellZoneName();
                if (cellZoneName.size())
                {
                    label zoneI = cellZones.findZoneID(cellZoneName);
                    if (zoneI != -1)
                    {
                        internalNamedCellZones.append(zoneI);
                    }
                    zoneI = faceZones.findZoneID(faceZoneName);
                    if (zoneI != -1)
                    {
                        internalNamedFaceZones.append(zoneI);
                    }
                }
            }
        }
    }
    labelHashSet internalNamedCellZonesSet(internalNamedCellZones);
    labelHashSet internalNamedFaceZonesSet(internalNamedFaceZones);

    forAll(cellZones, cellZoneI)
    {
        if (internalNamedCellZonesSet.found(cellZoneI))
        {
            const cellZone& cZone = cellZones[cellZoneI];
            forAll(cZone, czi)
            {
                label celli = cZone[czi];
                cellZoneID[celli] = cellZoneI;
            }
        }
    }
    labelList neiCellZoneID;
    syncTools::swapBoundaryCellList(mesh_, cellZoneID, neiCellZoneID);

    //boundary faces (-1 patch face, > -1 zone face)
    labelList boundaryFaces(mesh_.nFaces(), -2);
    boolList boundaryEdges(mesh_.nEdges(), false);
    forAll(mesh_.faces(), facei)
    {
        label patchI = patches.whichPatch(facei);
        if (patchI != -1 && !patches[patchI].coupled())
        {
            boundaryFaces[facei] = -1;
            const labelList& fEdges = mesh_.faceEdges()[facei];
            forAll(fEdges, fEI)
            {
                boundaryEdges[fEdges[fEI]] = true;
            }
        }
    }

    forAll(faceZones, faceZoneI)
    {
        if (internalNamedFaceZonesSet.found(faceZoneI))
        {
            const faceZone& fZone = faceZones[faceZoneI];
            forAll(fZone, fzi)
            {
                label facei = fZone[fzi];
                label patchI = patches.whichPatch(facei);
                if (patchI == -1 || patches[patchI].coupled())
                {
                    boundaryFaces[ facei] = faceZoneI;
                    const labelList& fEdges = mesh_.faceEdges()[facei];
                    forAll(fEdges, fEI)
                    {
                        boundaryEdges[fEdges[fEI]] = true;
                    }
                }
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    label nChanged = 0;
    labelList newZoneID(mesh_.nCells(), -2);
    forAll(mesh_.cells(), celli)
    {
        if (external && cellZoneID[celli] != -1)
        {
            continue;
        }
        else if (!external && cellZoneID[celli] == -1)
        {
            continue;
        }

        const labelList& cEdges = mesh_.cellEdges()[celli];
        label nBoundEdges = 0;
        forAll(cEdges, cEI)
        {
            label edgeI = cEdges[cEI];
            if (boundaryEdges[edgeI])
            {
                nBoundEdges++;
            }
        }

        if (nBoundEdges > cEdges.size() - 3)
        {
            const labelList& cFaces = mesh_.cells()[celli];
            label zoneFace = -1;
            bool resetCell = false;
            forAll(cFaces, cFI)
            {
                label facei = cFaces[cFI];
                label faceZoneI = boundaryFaces[facei];
                if (faceZoneI > -1)
                {
                    if (zoneFace != -1 && zoneFace != faceZoneI)
                    {
                        resetCell = false;
                        break;
                    }
                    else
                    {
                        resetCell = true;
                        zoneFace = faceZoneI;
                        if (cellZoneID[celli] == -1)
                        {
                            label patchI = patches.whichPatch(facei);
                            if (patchI == -1)
                            {
                                label own = mesh_.faceOwner()[facei];
                                newZoneID[celli] = (own == celli) ?
                                    cellZoneID[mesh_.faceNeighbour()[facei]]
                                    : cellZoneID[own];
                            }
                            else if (patches[patchI].coupled())
                            {
                                newZoneID[celli] =
                                    neiCellZoneID[facei-mesh_.nInternalFaces()];
                            }
                        }
                    }
                }
            }

            if (resetCell)
            {
                if (newZoneID[celli] == -2)
                {
                    newZoneID[celli] = -1;
                }
                meshMod.modifyCell
                (
                    celli,
                    newZoneID[celli]
                );
                nChanged++;
            }
            else
            {
                newZoneID[celli] = -2;
            }
        }
    }

    reduce(nChanged, sumOp<label>());
    if (nChanged > 0)
    {
        forAll(cellZoneID, celli)
        {
            if (newZoneID[celli] != -2)
            {
                cellZoneID[celli] = newZoneID[celli];
            }
        }
        syncTools::swapBoundaryCellList(mesh_, cellZoneID, neiCellZoneID);

        labelList newFaceZone(mesh_.nFaces(), -2);
        forAll(newZoneID, celli)
        {
            if (newZoneID[celli] != -2)
            {
                const labelList& cFaces = mesh_.cells()[celli];
                label cellFaceZone = -1;
                forAll(cFaces, cFI)
                {
                    label facei = cFaces[cFI];
                    label faceZoneI = boundaryFaces[facei];
                    if (faceZoneI > -1)
                    {
                        cellFaceZone =  faceZoneI;
                        break;
                    }
                }

                forAll(cFaces, cFI)
                {
                    label facei = cFaces[cFI];
                    label patchI = patches.whichPatch(facei);

                    label own = mesh_.faceOwner()[facei];
                    label cellZoneOwn = cellZoneID[own];
                    label cellZoneNei = -1;
                    if (patchI == -1 || patches[patchI].coupled())
                    {
                        if (patchI == -1)
                        {
                            label nei = mesh_.faceNeighbour()[facei];
                            cellZoneNei = cellZoneID[nei];
                        }
                        else
                        {
                            cellZoneNei =
                                neiCellZoneID[facei-mesh_.nInternalFaces()];
                        }

                        if (cellZoneOwn != cellZoneNei)
                        {
                            newFaceZone[facei] = cellFaceZone;
                        }
                        else
                        {
                            newFaceZone[facei] = -1;
                        }
                    }
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh_,
            newFaceZone,
            maxEqOp<label>()
        );

        forAll(newFaceZone, facei)
        {
            label newZoneID = newFaceZone[facei];
            if (newZoneID != -2)
            {
                label own = mesh_.faceOwner()[facei];

                label patchI = patches.whichPatch(facei);
                label nei = -1;
                if (patchI == -1)
                {
                    nei = mesh_.faceNeighbour()[facei];
                }

                label ownZone = cellZoneID[own];
                label neiZone = -1;
                if (patchI == -1 || patches[patchI].coupled())
                {
                    if (patchI == -1)
                    {
                        neiZone = cellZoneID[nei];
                    }
                    else
                    {
                        neiZone = neiCellZoneID[facei-mesh_.nInternalFaces()];
                    }
                }

                bool flip =
                (
                    ownZone == -1
                 || (neiZone != -1 && ownZone > neiZone)
                );

                meshMod.modifyFace
                (
                    mesh_.faces()[facei], // modified face
                    facei,     // label of face being modified
                    own,            // owner
                    nei,             // neighbour
                    false,          // face flip
                    patchI,      // new patch for face
                    newZoneID, // zone for face
                    flip        // face flip in zone
                );
            }
        }

        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing might not do this)
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
            // Delete mesh volumes. No other way to do this?
            mesh_.clearOut();
        }
        updateMesh(map, labelList(0));

        checkCoupledFaceZones(mesh_);
    }

    return nChanged;
}


void Foam::meshRefinement::closeSingleCellGaps
(
    const labelList& closureSurfaces,
    const pointField& locationsInMesh,
    const labelList& intersectedFaces,
    const labelList& wrapLevels,
    labelList& wrappedFaces
) const
{
    Info<< "Closing single cell gaps "<<endl;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& cellLevel = meshCutter().cellLevel();

    boolList blockedFaces(mesh_.nFaces(), false);
    labelList markedCells(mesh_.nCells(), -1);
    labelList boundaryPts(mesh_.nPoints(), -1);
    // Mark cells connected to boundary or wrapped faces
    forAll(mesh_.faces(), facei)
    {
        bool markPointCells = false;
        label regioni
        (
            intersectedFaces[facei] != -1 ?
            intersectedFaces[facei] : wrappedFaces[facei]
        );

        if (regioni != -1)
        {
            blockedFaces[facei] = true;
            if (wrappedFaces[facei] != -1)
            {
                markPointCells = true;
            }
            else if
            (
                intersectedFaces[facei] != -1
                && wrapLevels[intersectedFaces[facei]] > -1
            )
            {
                markPointCells = true;
            }

            const face f = mesh_.faces()[facei];
            forAll(f, fp)
            {
                label meshPointI = f[fp];
                boundaryPts[meshPointI] = regioni;
                if (markPointCells)
                {
                    const labelList& pCells = mesh_.pointCells()[meshPointI];
                    forAll(pCells, pCI)
                    {
                        label celli = pCells[pCI];
                        markedCells[celli] = regioni;
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        boundaryPts,
        maxEqOp<label>(),
        label(-1)
    );

    label nTest = 0;
    forAll(markedCells, celli)
    {
        if (markedCells[celli] != -1)
        {
            nTest++;
        }
    }

    const scalar edge0Len = meshCutter().level0EdgeLength();
    const pointField& cellCentres =  mesh_.cellCentres();

    pointField testPts(nTest);
    scalarField testRadius(nTest);
    scalar weight = 0.866;

    nTest = 0;
    forAll(markedCells, celli)
    {
        if (markedCells[celli] != -1)
        {
            label clevel = cellLevel[celli];

            scalar len = weight*edge0Len / pow(2., clevel);
            testPts[nTest] = cellCentres[celli];
            testRadius[nTest] = magSqr(len);
            nTest++;
        }
    }
    testPts.setSize(nTest);
    testRadius.setSize(nTest);

    if (returnReduce(nTest, sumOp<label>()) == 0)
    {
        return;
    }

    List<pointIndexHit> nearestInfo;
    labelList nearestSurface;
    surfaces_.findNearest
    (
        closureSurfaces,
        testPts,
        testRadius,    // sqr of attraction
        nearestSurface,
        nearestInfo
    );

    // Force calculation of face decomposition (used in findCell)
    (void)mesh_.tetBasePtIs();

    nTest = 0;
    // Test if boundary touching cells contain geometry
    forAll(markedCells, celli)
    {
        if (markedCells[celli] != -1)
        {
            if (!nearestInfo[nTest].hit())
            {
                markedCells[celli] = -1;
            }
            nTest++;
        }
    }

    if (debug)
    {
        mesh_.write();

        DynamicList<label> closureCells;

        forAll(mesh_.cells(), celli)
        {
            if (markedCells[celli] != -1)
            {
                closureCells.append(celli);
            }
        }
        cellSet closeCells
        (
            mesh_,
            "closureCells",
            labelHashSet(closureCells)
        );

        closeCells.instance() = timeName();
        closeCells.write();
    }

    labelList neiMarkedCells;
    syncTools::swapBoundaryCellList(mesh_, markedCells, neiMarkedCells);
    forAll(markedCells, celli)
    {
        if (markedCells[celli] != -1)
        {
            const cell c = mesh_.cells()[celli];
            forAll(c, cFI)
            {
                label facei = c[cFI];
                blockedFaces[facei] = true;
            }
        }
    }

    syncTools::syncFaceList(mesh_, blockedFaces, orEqOp<bool>());
    regionSplit globalRegion(mesh_, blockedFaces);

    labelList neiGlobalRegion;
    syncTools::swapBoundaryCellList(mesh_, globalRegion, neiGlobalRegion);

    DynamicList<label> regionsToKeep(locationsInMesh.size());
    labelList keepCells(locationsInMesh.size(), -1);

    // For all locationsInMesh find the cell
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];

        // Find the region containing the keepPoint
        label keepRegionI = -1;

        keepCells[i] = findCell
        (
            insidePoint,
            mesh_,
            meshCutter_
        );

        if (keepCells[i] != -1)
        {
            keepRegionI = globalRegion[keepCells[i]];
        }
        reduce(keepRegionI, maxOp<label>());
        regionsToKeep.append(keepRegionI);
    }
    regionsToKeep.shrink();

    boolList regionIsKept(globalRegion.nRegions(), false);
    forAll(regionsToKeep, keepI)
    {
        regionIsKept[regionsToKeep[keepI]] = true;
    }

    Pstream::listCombineReduce(regionIsKept, orOp<bool>());

    labelList keepRegionPts(mesh_.nPoints(), -1);
    //Mark cells in keep region and extend to include neighbouring
    // non boundary point cells and face cells
    forAll(globalRegion, celli)
    {
        label regionI = globalRegion[celli];
        if (regionIsKept[regionI])
        {
            const labelList& cPts = mesh_.cellPoints()[celli];
            forAll(cPts, cPtI)
            {
                keepRegionPts[cPts[cPtI]] = regionI;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        keepRegionPts,
        maxEqOp<label>(),
        label(-1)
    );

    labelList globalRegionTmp = globalRegion;

    forAll(markedCells, celli)
    {
        if (markedCells[celli] != -1)
        {
            const labelList& cPts = mesh_.cellPoints()[celli];
            bool foundPt = false;
            forAll(cPts, cPtI)
            {
                label meshpointi = cPts[cPtI];
                if
                (
                    keepRegionPts[meshpointi] != -1
                    && boundaryPts[meshpointi] == -1
                )
                {
                    globalRegionTmp[celli] = keepRegionPts[meshpointi];
                    foundPt = true;
                    break;
                }
            }

            if (!foundPt)
            {
                const labelList& c = mesh_.cells()[celli];
                forAll(c, cFI)
                {
                    label facei = c[cFI];
                    label patchI = patches.whichPatch(facei);

                    if
                    (
                        intersectedFaces[facei] != -1 || wrappedFaces[facei] != -1
                        || (patchI != -1 && !patches[patchI].coupled())
                    )
                    {
                        continue;
                    }

                    label neiRegion = -1;
                    if (patchI == -1)
                    {
                        if (mesh_.faceOwner()[facei] == celli)
                        {
                            neiRegion =
                                globalRegion[mesh_.faceNeighbour()[facei]];
                        }
                        else
                        {
                            neiRegion = globalRegion[mesh_.faceOwner()[facei]];
                        }
                    }
                    else
                    {
                        neiRegion =
                            neiGlobalRegion[facei-mesh_.nInternalFaces()];
                    }

                    if (regionIsKept[neiRegion])
                    {
                        globalRegionTmp[celli] = neiRegion;
                        break;
                    }
                }
            }
        }
    }

    forAll(globalRegion, celli)
    {
        globalRegion[celli] = globalRegionTmp[celli];
    }

    syncTools::swapBoundaryCellList(mesh_, globalRegion, neiGlobalRegion);

    scalarField wrapArea(globalRegion.nRegions(), scalar(0));
    scalarField regionVolume(globalRegion.nRegions(), scalar(0));

    forAll(mesh_.cells(), celli)
    {
        label regioni = globalRegion[celli];
        regionVolume[regioni] += mesh_.cellVolumes()[celli];
        const labelList& cFaces = mesh_.cells()[celli];
        forAll(cFaces, cFI)
        {
            label facei = cFaces[cFI];
            label patchI = patches.whichPatch(facei);

            if
            (
                intersectedFaces[facei] != -1 ||
                (patchI != -1 && !patches[patchI].coupled())
            )
            {
                continue;
            }

            label neiRegion = -1;
            if (patchI == -1)
            {
                if (mesh_.faceOwner()[facei] == celli)
                {
                    neiRegion =
                        globalRegion[mesh_.faceNeighbour()[facei]];
                }
                else
                {
                    neiRegion = globalRegion[mesh_.faceOwner()[facei]];
                }
            }
            else
            {
                neiRegion = neiGlobalRegion[facei-mesh_.nInternalFaces()];
            }

            if (regioni != neiRegion)
            {
                wrapArea[regioni] = mesh_.magFaceAreas()[facei];
            }
        }
    }

    Pstream::listCombineReduce(wrapArea, plusOp<scalar>());
    Pstream::listCombineReduce(regionVolume, plusOp<scalar>());

    //Set region type (0: keep, 1: remove, 2: indeterminate)
    labelList regionType(mesh_.nCells(), 2);

    label nRegionsRemove = surfaces_.nRegionsWrapRemoval();
    if (nRegionsRemove > 0)
    {
        DynamicList<label> removalRegions(nRegionsRemove);
        SortableList<scalar> sortedVolumes(regionVolume);
        const labelList& sortedIndices = sortedVolumes.indices();

        forAllReverse(sortedIndices, i)
        {
            label regioni = sortedIndices[i];
            if (!regionIsKept[regioni] && wrapArea[regioni] > SMALL)
            {
                removalRegions.append(regioni);
                if (nRegionsRemove == removalRegions.size())
                {
                    break;
                }
            }
        }

        labelHashSet removalRegionsSet(removalRegions);
        forAll(mesh_.cells(), celli)
        {
            label regioni = globalRegion[celli];

            if (regionIsKept[regioni])
            {
                regionType[celli] = 0;
            }
            else if (removalRegionsSet.found(regioni))
            {
                regionType[celli] = 1;
            }
        }
    }
    else
    {
        scalarField areaToVol(globalRegion.nRegions(),GREAT);

        forAll(areaToVol, regioni)
        {
            if (regionVolume[regioni] > SMALL)
            {
                areaToVol[regioni] =
                    pow(wrapArea[regioni],0.5)/pow(regionVolume[regioni],0.333);
            }
        }

        forAll(mesh_.cells(), celli)
        {
            label regioni = globalRegion[celli];

            if (regionIsKept[regioni])
            {
                regionType[celli] = 0;
            }
            else if
            (
                regionVolume[regioni] > SMALL && wrapArea[regioni] > SMALL
                && areaToVol[regioni] < 0.05
            )
            {
                regionType[celli] = 1;
            }
        }
    }

    labelList neiRegionType;
    syncTools::swapBoundaryCellList(mesh_, regionType, neiRegionType);

    blockedFaces = false;
    forAll(mesh_.cells(), celli)
    {
        label ownRegion = regionType[celli];
        const labelList& cFaces = mesh_.cells()[celli];
        forAll(cFaces, cFI)
        {
            label facei = cFaces[cFI];
            label patchI = patches.whichPatch(facei);

            if
            (
                intersectedFaces[facei] != -1 ||
                (patchI != -1 && !patches[patchI].coupled())
            )
            {
                blockedFaces[facei] = true;
                continue;
            }

            label neiRegion = -1;
            if (patchI == -1)
            {
                if (mesh_.faceOwner()[facei] == celli)
                {
                    neiRegion =
                        regionType[mesh_.faceNeighbour()[facei]];
                }
                else
                {
                    neiRegion = regionType[mesh_.faceOwner()[facei]];
                }
            }
            else
            {
                neiRegion = neiRegionType[facei-mesh_.nInternalFaces()];
            }

            if (ownRegion != neiRegion)
            {
                blockedFaces[facei] = true;
            }
        }
    }
    syncTools::syncFaceList(mesh_, blockedFaces, orEqOp<bool>());
    regionSplit globalRegionFinal(mesh_, blockedFaces);

    scalarField keepArea(globalRegionFinal.nRegions(), scalar(0));
    scalarField removeArea(globalRegionFinal.nRegions(), scalar(0));
    forAll(mesh_.cells(), celli)
    {
        label regioni = globalRegionFinal[celli];
        label ownRegion = regionType[celli];

        if (ownRegion != 2)
        {
            continue;
        }
        const labelList& cFaces = mesh_.cells()[celli];
        forAll(cFaces, cFI)
        {
            label facei = cFaces[cFI];
            label patchI = patches.whichPatch(facei);

            if
            (
                intersectedFaces[facei] != -1 ||
                (patchI != -1 && !patches[patchI].coupled())
            )
            {
                continue;
            }

            label neiRegion = -1;
            if (patchI == -1)
            {
                if (mesh_.faceOwner()[facei] == celli)
                {
                    neiRegion =
                        regionType[mesh_.faceNeighbour()[facei]];
                }
                else
                {
                    neiRegion = regionType[mesh_.faceOwner()[facei]];
                }
            }
            else
            {
                neiRegion = neiRegionType[facei-mesh_.nInternalFaces()];
            }

            if (neiRegion == 0)
            {
                 keepArea[regioni] += mesh_.magFaceAreas()[facei];
            }
            else if (neiRegion == 1)
            {
                 removeArea[regioni] += mesh_.magFaceAreas()[facei];
            }
        }
    }

    forAll(mesh_.faces(), facei)
    {
        label patchI = patches.whichPatch(facei);

        if
        (
            intersectedFaces[facei] != -1 || wrappedFaces[facei] != -1
            || (patchI != -1 && !patches[patchI].coupled())
        )
        {
            continue;
        }

        label own = mesh_.faceOwner()[facei];
        label markedOwn = markedCells[own];
        label markedNei = -1;
        label ownType = regionType[own];
        label neiType = -1;
        if (patchI == -1)
        {
            label nei = mesh_.faceNeighbour()[facei];
            neiType = regionType[nei];
            markedNei = markedCells[nei];
        }
        else
        {
            neiType = neiRegionType[facei-mesh_.nInternalFaces()];
            markedNei = neiMarkedCells[facei-mesh_.nInternalFaces()];
        }

        if
        (
            ownType == 1 || neiType == 1
        )
        {
            wrappedFaces[facei] = max(0,max(markedOwn, markedNei));
        }
    }


    Pstream::listCombineReduce(removeArea, plusOp<scalar>());
    Pstream::listCombineReduce(keepArea, plusOp<scalar>());

    forAll(mesh_.cells(), celli)
    {
        label regioni = globalRegionFinal[celli];
        label ownRegion = regionType[celli];

        if (ownRegion == 2)
        {
            if (removeArea[regioni] > SMALL && keepArea[regioni] < SMALL)
            {
                regionType[celli] = 1;
            }
            else if (removeArea[regioni] < SMALL && keepArea[regioni] > SMALL)
            {
                regionType[celli] = 0;
            }
            else if (removeArea[regioni] > SMALL && keepArea[regioni] > SMALL)
            {
                if (removeArea[regioni] < keepArea[regioni])
                {
                    regionType[celli] = 0;
                }
                else
                {
                    regionType[celli] = 1;
                }
            }
            else
            {
                regionType[celli] = 1;
            }
        }
    }

    if (debug)
    {
        mesh_.write();
        volScalarField wrapRegionType
        (
            IOobject
            (
                "wrapRegionType",
                timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(mesh_.cells(), celli)
        {
            wrapRegionType[celli] = regionType[celli];
        }
        wrapRegionType.write();

        const_cast<Time&>(mesh_.time())++;
    }

    syncTools::swapBoundaryCellList(mesh_, regionType, neiRegionType);

    forAll(mesh_.faces(), facei)
    {
        label patchI = patches.whichPatch(facei);

        if
        (
            intersectedFaces[facei] != -1
            || (patchI != -1 && !patches[patchI].coupled())
        )
        {
            continue;
        }

        label own = mesh_.faceOwner()[facei];
        label markedOwn = markedCells[own];
        label markedNei = -1;
        label ownType = regionType[own];
        label neiType = -1;
        if (patchI == -1)
        {
            label nei = mesh_.faceNeighbour()[facei];
            neiType = regionType[nei];
            markedNei = markedCells[nei];
        }
        else
        {
            neiType = neiRegionType[facei-mesh_.nInternalFaces()];
            markedNei = neiMarkedCells[facei-mesh_.nInternalFaces()];
        }

        if
        (
            (ownType == 0 && neiType == 1)
            || (ownType == 1 && neiType == 0)
        )
        {
            wrappedFaces[facei] = max(0,max(markedOwn, markedNei));
        }
        else
        {
            wrappedFaces[facei] = -1;
        }
    }

    syncTools::syncFaceList(mesh_, wrappedFaces, maxEqOp<label>());

    return;
}


void Foam::meshRefinement::rebaffleThinGapCells
(
    const pointField& locationsInMesh,
    const pointField& cellCentres,
    labelList& globalRegion1,
    labelList& globalRegion2,
    labelList& testFaces
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList thinGapSurfaces = surfaces_.thinGapSurfaces();

    label defaultPatch = 0;

    if (thinGapSurfaces.size() == 0)
    {
        return;
    }

    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    if (debug)
    {
        DynamicList<label> boundaryFaces;

        forAll(mesh_.faces(), facei)
        {
            if (globalRegion1[facei] > -1 && globalRegion2[facei] > -1)
            {
                boundaryFaces.append(facei);
            }
        }

        const indirectPrimitivePatch boundaryPatch
        (
            IndirectList<face>
            (
                mesh_.faces(),
                boundaryFaces
            ),
            mesh_.points()
        );

        static label startGapIter = 0;

        simpleVTKWriter
        (
            boundaryPatch.localFaces(),
            boundaryPatch.localPoints()
        ).write("startingGapPatchFaces"+Foam::name(startGapIter)+".vtk");

        startGapIter++;
    }

    // Check for thin gap surfaces and keep owner/neighbour cells
    boolList gapCells;
    boolList protectedFaces(mesh_.nFaces(), false);
    forAll(globalRegion1, facei)
    {
        if (globalRegion1[facei] > -1)
        {
            protectedFaces[facei] = true;
        }
    }

    DynamicList<label> baffleFaces;
    markThinGapCells
    (
        cellCentres,
        gapCells,
        protectedFaces,
        baffleFaces
    );

    if (debug)
    {
        mesh_.write();
        DynamicList<label> gCells;

        forAll(mesh_.cells(), celli)
        {
            if (gapCells[celli])
            {
                gCells.append(celli);
            }
        }
        cellSet gapCellSet
        (
            mesh_,
            "gapCells",
            labelHashSet(gCells)
        );

        gapCellSet.instance() = timeName();
        gapCellSet.write();
    }

    boolList neiGapCells;
    syncTools::swapBoundaryCellList(mesh_, gapCells, neiGapCells);

    boolList gapCellExternalFaces(mesh_.nFaces(), false);

    label nGapCells = 0;
    forAll(mesh_.cells(), celli)
    {
        if (gapCells[celli])
        {
            nGapCells++;
            const cell& c = mesh_.cells()[celli];
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchI = patches.whichPatch(facei);
                bool neiGapCell = false;

                if (patchI == -1)
                {
                    label nei = mesh_.faceNeighbour()[facei];
                    if (nei == celli)
                    {
                        nei = mesh_.faceOwner()[facei];
                        neiGapCell = gapCells[nei];
                    }
                    else
                    {
                        neiGapCell = gapCells[nei];
                    }
                }
                else if (patches[patchI].coupled())
                {
                    neiGapCell =
                        neiGapCells[facei-mesh_.nInternalFaces()];
                }

                if (!neiGapCell)
                {
                    gapCellExternalFaces[facei] = true;
                }
            }
        }
    }
    syncTools::syncFaceList(mesh_, gapCellExternalFaces, orEqOp<bool>());

    Info<< "Selected " << returnReduce(nGapCells, sumOp<label>())
        << " gap cells" << endl;

    //Check for pinched edges and points
    boolList pinchedEdges(mesh_.nEdges(), false);
    boolList pinchedPoints(mesh_.nPoints(), false);

    labelList nGapEdgeFaces(mesh_.nEdges(), 0);
    forAll(mesh_.edges(), edgei)
    {
        const labelList& eFaces = mesh_.edgeFaces()[edgei];
        forAll(eFaces, eFI)
        {
            label facei = eFaces[eFI];
            if (isMasterFace[facei] && gapCellExternalFaces[facei])
            {
                nGapEdgeFaces[edgei]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        nGapEdgeFaces,
        plusEqOp<label>(),
        label(0)
    );

    boolList pinchedEdgePts(mesh_.nPoints(), false);

    forAll(mesh_.edges(), edgei)
    {
        if (nGapEdgeFaces[edgei] == 4)
        {
            const edge e = mesh_.edges()[edgei];
            pinchedEdges[edgei] = true;
            pinchedEdgePts[e[0]] = true;
            pinchedEdgePts[e[1]] = true;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        pinchedEdgePts,
        orEqOp<bool>(),
        false
    );

    labelList nGapPointFaces(mesh_.nPoints(), 0);
    forAll(mesh_.points(), pointi)
    {
        const labelList& pFaces = mesh_.pointFaces()[pointi];
        forAll(pFaces, pFI)
        {
            label facei = pFaces[pFI];
            if (isMasterFace[facei] && gapCellExternalFaces[facei])
            {
                nGapPointFaces[pointi]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        nGapPointFaces,
        plusEqOp<label>(),
        label(0)
    );

    forAll(mesh_.points(), pointi)
    {
        if (nGapPointFaces[pointi] == 6 && !pinchedEdgePts[pointi])
        {
            pinchedPoints[pointi] = true;
        }
    }

    //select regions to keep
    boolList blockedFaces(mesh_.nFaces(), false);
    forAll(mesh_.faces(), facei)
    {
        if (globalRegion1[facei] > -1 && globalRegion2[facei] > -1)
        {
            blockedFaces[facei] = true;
        }
    }

    syncTools::syncFaceList(mesh_, blockedFaces, orEqOp<bool>());
    regionSplit splitRegions(mesh_, blockedFaces);

    labelList neiSplitRegions;
    syncTools::swapBoundaryCellList(mesh_, splitRegions, neiSplitRegions);

    DynamicList<label> regionsToKeep(locationsInMesh.size());
    labelList selectedKeepCells(locationsInMesh.size(), -1);

    // For all locationsInMesh find the cell
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];

        // Find the region containing the keepPoint
        label keepRegionI = -1;

        selectedKeepCells[i] = findCell
        (
            insidePoint,
            mesh_,
            meshCutter_
        );

        if (selectedKeepCells[i] != -1)
        {
            keepRegionI = splitRegions[selectedKeepCells[i]];
        }
        reduce(keepRegionI, maxOp<label>());
        regionsToKeep.append(keepRegionI);
    }
    regionsToKeep.shrink();

    labelList regionIsKept(splitRegions.nRegions(), -1);
    forAll(regionsToKeep, keepI)
    {
        regionIsKept[regionsToKeep[keepI]] = keepI;
    }

    Pstream::listCombineReduce(regionIsKept, maxOp<label>());

    boolList keepCells(mesh_.nCells(), false);
    forAll(splitRegions, celli)
    {
        if (regionIsKept[splitRegions[celli]] > -1)
        {
            keepCells[celli] = true;
        }
    }

    boolList neiKeepCells;
    syncTools::swapBoundaryCellList(mesh_, keepCells, neiKeepCells);

    //Mark gap cells connect to keep cells by single edge
    labelList nEdgeGapCells(mesh_.nEdges(), 0);
    labelList nEdgeKeepCells(mesh_.nEdges(), 0);

    forAll(mesh_.edges(), edgei)
    {
        const labelList& eCells = mesh_.edgeCells()[edgei];
        forAll(eCells, eci)
        {
            label celli = eCells[eci];

            if (keepCells[celli])
            {
                nEdgeKeepCells[edgei]++;
            }

            if (gapCells[celli])
            {
                nEdgeGapCells[edgei]++;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        nEdgeGapCells,
        plusEqOp<label>(),
        label(0)
    );

    syncTools::syncEdgeList
    (
        mesh_,
        nEdgeKeepCells,
        plusEqOp<label>(),
        label(0)
    );


    //Type of gap cell: not in gap (-1), baffle (1), inside (2);
    labelList gapCellType(mesh_.nCells(), -1);

    forAll(baffleFaces, i)
    {
        label facei = baffleFaces[i];
        label patchI = patches.whichPatch(facei);
        gapCellType[mesh_.faceOwner()[facei]] = 1;
        if (patchI == -1)
        {
            gapCellType[mesh_.faceNeighbour()[facei]] = 1;
        }
    }

    forAll(gapCells, celli)
    {
        if (gapCells[celli] && gapCellType[celli] == -1)
        {
            gapCellType[celli] = 2;
        }
    }

    forAll(mesh_.edges(), edgei)
    {
        if
        (
            pinchedEdges[edgei]
            || (nEdgeGapCells[edgei] == 1 && nEdgeKeepCells[edgei] ==1)
        )
        {
            const labelList& eCells = mesh_.edgeCells()[edgei];
            forAll(eCells, eCI)
            {
                label celli = eCells[eCI];
                if (gapCellType[celli] == -1)
                {
                    gapCellType[celli] = 1;
                }
            }
        }
    }

    forAll(mesh_.points(), pointi)
    {
        if (pinchedPoints[pointi])
        {
            const labelList& pCells = mesh_.pointCells()[pointi];
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                if (gapCellType[celli] == -1)
                {
                    gapCellType[celli] = 1;
                }
            }
        }
    }

    //Mark additional gap cells if entry faces become pinched
    labelList gapFaceType(mesh_.nFaces(), -1);
    forAll(mesh_.cells(), celli)
    {
        label gCellType = gapCellType[celli];
        if (gCellType > -1)
        {
            const cell& c = mesh_.cells()[celli];
            forAll(c, cFI)
            {
                label facei = c[cFI];
                gapFaceType[facei] = max(gapFaceType[facei], gCellType);
            }
        }
    }

    labelList neiGapCellType;
    syncTools::swapBoundaryCellList(mesh_, gapCellType, neiGapCellType);

    syncTools::syncFaceList(mesh_, gapFaceType, maxEqOp<label>());
    DynamicList<label> interfaceFaces(mesh_.nFaces()/100);
    forAll(mesh_.faces(), facei)
    {
        if (gapFaceType[facei] > -1)
        {
            label patchI = patches.whichPatch(facei);
            label own = mesh_.faceOwner()[facei];
            bool oKeepCell = keepCells[own];
            bool nKeepCell = false;

            label oGapCellType = gapCellType[own];
            label nGapCellType = -1;

            if (patchI == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                nKeepCell = keepCells[nei];
                nGapCellType = gapCellType[nei];
            }
            else if (patches[patchI].coupled())
            {
                nKeepCell =
                    neiKeepCells[facei-mesh_.nInternalFaces()];
                nGapCellType =
                    neiGapCellType[facei-mesh_.nInternalFaces()];
            }

            if
            (
                (nGapCellType == -1 || oGapCellType == -1)
                &&
                (
                    (oKeepCell && nGapCellType != -1)
                    || (nKeepCell && oGapCellType != -1)
                )
            )
            {
                interfaceFaces.append(facei);
            }
        }
    }

    const indirectPrimitivePatch interfacePatch
    (
        IndirectList<face>
        (
            mesh_.faces(),
            interfaceFaces
        ),
        mesh_.points()
    );

    const labelList& interfacePts = interfacePatch.meshPoints();

    forAll(interfacePts, ptI)
    {
        if
        (
            interfacePatch.pointFaces()[ptI].size() == 2
            && interfacePatch.pointEdges()[ptI].size() == 4
        )
        {
            const labelList& pCells = mesh_.pointCells()[interfacePts[ptI]];
            forAll(pCells, pCI)
            {
                label celli = pCells[pCI];
                if (!keepCells[celli] && gapCellType[celli] == -1)
                {
                    gapCellType[celli] = 1;
                }
            }
        }
    }

    syncTools::swapBoundaryCellList(mesh_, gapCellType, neiGapCellType);

    if (createGapField_ && !mesh_.foundObject<volScalarField>("gapCells"))
    {
        Info<<"Creating gap field"<<endl;

        //Store the gap cells for later re-use
        volScalarField* gapCellPtr
        (
            new volScalarField
            (
                IOobject
                (
                    "gapCells",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                 ),
                mesh_,
                dimensionedScalar("gapCells", dimless, -1),
                zeroGradientFvPatchScalarField::typeName
             )
         );

        volScalarField& gapCellsField = *gapCellPtr;

        forAll(mesh_.cells(), celli)
        {
            if (gapCellType[celli] > -1)
            {
                gapCellsField[celli] = scalar(1);
            }
        }

        mesh_.objectRegistry::store(gapCellPtr);

        if (debug)
        {
            mesh_.write();
            gapCellsField.write();
        }
    }

    gapFaceType = -1;
    forAll(mesh_.cells(), celli)
    {
        label gCellType = gapCellType[celli];
        if (gCellType > -1)
        {
            const cell& c = mesh_.cells()[celli];
            forAll(c, cFI)
            {
                label facei = c[cFI];
                gapFaceType[facei] = max(gapFaceType[facei], gCellType);
            }
        }
    }
    syncTools::syncFaceList(mesh_, gapFaceType, maxEqOp<label>());

    DynamicList<label> addedFaces(mesh_.nFaces()/10);
    forAll(mesh_.faces(), facei)
    {
        if (gapFaceType[facei] > -1)
        {
            label patchI = patches.whichPatch(facei);
            label own = mesh_.faceOwner()[facei];
            bool oKeepCell = keepCells[own];
            bool nKeepCell = false;

            label oGapCellType = gapCellType[own];
            label nGapCellType = -1;

            if (patchI == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                nKeepCell = keepCells[nei];
                nGapCellType = gapCellType[nei];
            }
            else if (patches[patchI].coupled())
            {
                nKeepCell =
                    neiKeepCells[facei-mesh_.nInternalFaces()];
                nGapCellType =
                    neiGapCellType[facei-mesh_.nInternalFaces()];
            }

            if
            (
                nGapCellType != -1 && oGapCellType != -1
            )
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
            }
            else if
            (
                !protectedFaces[facei] &&
                (
                    (oKeepCell && nGapCellType != -1)
                    || (nKeepCell && oGapCellType != -1)
                )
            )
            {
                globalRegion1[facei] = -1;
                globalRegion2[facei] = -1;
            }
            else if (globalRegion1[facei] == -1)
            {
                globalRegion1[facei] = defaultPatch;
                globalRegion2[facei] = defaultPatch;
                addedFaces.append(facei);
            }
        }
    }

    //Ensure testFaces includes added patch faces
    if (addedFaces.size() > 0)
    {
        label sz = testFaces.size();
        testFaces.setSize(sz+addedFaces.size());
        forAll(addedFaces, i)
        {
            testFaces[sz+i] = addedFaces[i];
        }
    }

    if (debug)
    {
        DynamicList<label> boundaryFaces;

        forAll(mesh_.faces(), facei)
        {
            if (globalRegion1[facei] > -1 && globalRegion2[facei] > -1)
            {
                boundaryFaces.append(facei);
            }
        }

        const indirectPrimitivePatch boundaryPatch
        (
            IndirectList<face>
            (
                mesh_.faces(),
                boundaryFaces
            ),
            mesh_.points()
        );

        static label resetGapIter = 0;

        simpleVTKWriter
        (
            boundaryPatch.localFaces(),
            boundaryPatch.localPoints()
        ).write("resetGapPatchFaces"+Foam::name(resetGapIter)+".vtk");

        resetGapIter++;
    }

    return;
}


void Foam::meshRefinement::markThinGapCells
(
    const pointField& cellCentres,
    boolList& gapCells,
    boolList& protectedFaces,
    DynamicList<label>& baffleFaces
) const
{
    labelList thinGapSurfaces = surfaces_.thinGapSurfaces();

    if (thinGapSurfaces.size() == 0)
    {
        return;
    }
    else
    {
        Info<<"Checking for thin gap cells"<<endl;
    }

    scalar minGapSize = 1e-4;

    gapCells.setSize(mesh_.nCells(), false);
    baffleFaces.setCapacity(mesh_.nFaces()/10);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    labelList testFaces = intersectedAndNeighbouring();

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(cellCentres, neiLevel, neiCc);

    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());

    {
        labelList minLevel;
        calcCellCellRays
        (
            cellCentres,
            neiCc,
            labelList(neiCc.size(), -1),
            testFaces,
            start,
            end,
            minLevel
        );
    }

    // Do test for intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;

    surfaces_.findNearestIntersection
    (
        thinGapSurfaces,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
    );


    pointField normal1(hit1.size(), vector(GREAT, GREAT, GREAT));
    pointField normal2(hit2.size(), vector(GREAT, GREAT, GREAT));
    forAll(thinGapSurfaces, sI)
    {
        label surfI = thinGapSurfaces[sI];
        label geomI = surfaces_.surfaces()[surfI];

        //normal for hit point1
        {
            DynamicList<pointIndexHit> localHits;

            forAll(hit1, i)
            {
                if (surface1[i] == surfI)
                {
                    localHits.append(hit1[i]);
                }
            }
            localHits.shrink();

            pointField localNormals;
            surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

            label localI = 0;
            forAll(hit1, i)
            {
                if (surface1[i] == surfI)
                {
                    normal1[i] = localNormals[localI];
                    localI++;
                }
            }
        }

        //normal for hit point2
        {
            DynamicList<pointIndexHit> localHits;
            forAll(hit2, i)
            {
                if (surface2[i] == surfI)
                {
                    localHits.append(hit2[i]);
                }
            }
            localHits.shrink();

            pointField localNormals;

            surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

            label localI = 0;
            forAll(hit2, i)
            {
                if (surface2[i] == surfI)
                {
                    normal2[i] = localNormals[localI];
                    localI++;
                }
            }
        }
    }


    List<DynamicList<point>> cellNorm(mesh_.nCells());
    List<DynamicList<point>> cellHit(mesh_.nCells());

    DynamicList<label> retestFaces(testFaces.size());
    forAll(testFaces, i)
    {
        if (hit1[i].hit() && hit2[i].hit())
        {
            label facei = testFaces[i];

            protectedFaces[facei] = false;

            label patchI = patches.whichPatch(facei);
            label own = mesh_.faceOwner()[facei];

            cellNorm[own].append(normal1[i]);
            cellHit[own].append(hit1[i].hitPoint());
            if (patchI == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                cellNorm[nei].append(normal2[i]);
                cellHit[nei].append(hit2[i].hitPoint());
            }

            vector separationVec = hit2[i].hitPoint()-hit1[i].hitPoint();
            scalar separationDist = 0.5*
            (
                mag(separationVec & normal1[i]) +
                mag(separationVec & normal2[i])
            );

            if
            (
                separationDist < minGapSize
                || mag(normal1[i] & normal2[i]) < 0.707
            )
            {
                if (separationDist > minGapSize)
                {
                    retestFaces.append(facei);
                }
                continue;
            }

            baffleFaces.append(facei);

            gapCells[own] = true;

            if (patchI == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                gapCells[nei] = true;
            }
        }
    }

    //retest on selected faces
    {
        pointField reteststart(retestFaces.size());
        pointField retestend(retestFaces.size());

        {
            labelList minLevel;
            calcCellCellRays
            (
                cellCentres,
                neiCc,
                labelList(neiCc.size(), -1),
                retestFaces,
                reteststart,
                retestend,
                minLevel
            );
        }

        // Do test for intersections
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(thinGapSurfaces, i)
        {
            label surfi = thinGapSurfaces[i];
            label geomI = surfaces_.surfaces()[surfi];
            List<List<pointIndexHit>> hitInfo;
            surfaces_.geometry()[geomI].findLineAll
            (
                reteststart,
                retestend,
                hitInfo
            );

            forAll(hitInfo, hiti)
            {
                if (hitInfo[hiti].size() > 1)
                {
                    label facei = retestFaces[hiti];
                    pointField localNormals;
                    surfaces_.geometry()[geomI].getNormal
                    (
                        hitInfo[hiti],
                        localNormals
                     );

                    bool inGap = false;

                    for (label i = 0; i < localNormals.size(); i++)
                    {
                        for (label j = i+1; j < localNormals.size(); j++)
                        {
                            vector lnormal1 = localNormals[i];
                            vector lnormal2 = localNormals[j];

                            vector separationVec =
                                hitInfo[hiti][i].hitPoint()
                                -hitInfo[hiti][j].hitPoint();
                            scalar separationDist = 0.5*
                            (
                                mag(separationVec & lnormal1) +
                                mag(separationVec & lnormal2)
                            );

                            if
                            (
                                separationDist > minGapSize
                                && mag(lnormal1 & lnormal2) > 0.707
                            )
                            {
                                inGap = true;
                                break;
                            }
                        }
                        if (inGap)
                        {
                            break;
                        }
                    }
                    if (inGap)
                    {
                        label patchI = patches.whichPatch(facei);
                        label own = mesh_.faceOwner()[facei];

                        baffleFaces.append(facei);
                        gapCells[own] = true;

                        if (patchI == -1)
                        {
                            label nei = mesh_.faceNeighbour()[facei];
                            gapCells[nei] = true;
                        }
                    }
                }
            }
        }
    }

    forAll(cellNorm, celli)
    {
        const DynamicList<point>& cNorms = cellNorm[celli];
        const DynamicList<point>& cHits = cellHit[celli];

        if (cNorms.size() > 1)
        {
            for (label i = 0; i < cNorms.size(); i++)
            {
                for (label j = i+1; j < cNorms.size(); j++)
                {
                    vector hp1 = cHits[i];
                    vector hp2 = cHits[j];

                    vector np1 = cNorms[i];
                    vector np2 = cNorms[j];

                    vector separationVec = hp1-hp2;
                    scalar separationDist = 0.5*
                    (
                        mag(separationVec & np1) +
                        mag(separationVec & np2)
                    );

                    separationVec /= mag(separationVec) + SMALL;

                    if
                    (
                        separationDist < minGapSize
                        || mag(np1 & np2) < 0.707
                        || mag(separationVec & np1) < 0.5
                        || mag(separationVec & np2) < 0.5
                    )
                    {
                        continue;
                    }

                    gapCells[celli] = true;
                }
            }
        }
    }

    return;
}


void Foam::meshRefinement::createSingleCellGap()
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    labelList thinGapSurfaces = surfaces_.thinGapSurfaces();
    if (thinGapSurfaces.size() == 0)
    {
        return;
    }
    else
    {
        Info<<"Creating single cell gap"<<endl;
    }

    boolList gapCells(mesh_.nCells(), false);

    //lookup gap cell field
    {
        const volScalarField& gapCellsField =
            mesh_.lookupObject<volScalarField>("gapCells");

        forAll(gapCellsField, celli)
        {
            if (gapCellsField[celli] > -1)
            {
                gapCells[celli] = true;
            }
        }

        if (debug)
        {
            mesh_.write();
            gapCellsField.write();
        }
    }

    boolList blockedFaces(mesh_.nFaces(), false);

    //select gap cell regions
    boolList neiGapCells;
    syncTools::swapBoundaryCellList(mesh_, gapCells, neiGapCells);

    DynamicList<label> gapBoundaryFaces(mesh_.nFaces()/10);
    forAll(mesh_.faces(), facei)
    {
        label patchi = patches.whichPatch(facei);
        label own = mesh_.faceOwner()[facei];
        bool gapCellOwn = gapCells[own];

        if (patchi == -1 || patches[patchi].coupled())
        {
            bool gapCellNei = false;

            if (patchi == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                gapCellNei = gapCells[nei];
            }
            else
            {
                gapCellNei = neiGapCells[facei-mesh_.nInternalFaces()];
            }

            if
            (
                (gapCellOwn && !gapCellNei)
                || (!gapCellOwn && gapCellNei)
            )
            {
                blockedFaces[facei] = true;
            }
        }
        else if (gapCellOwn && !blockedFaces[facei])
        {
            gapBoundaryFaces.append(facei);
        }
    }

    syncTools::syncFaceList(mesh_, blockedFaces, orEqOp<bool>());
    regionSplit globalRegion(mesh_, blockedFaces);

    const indirectPrimitivePatch boundaryPatch
    (
        IndirectList<face>
        (
            mesh_.faces(),
            gapBoundaryFaces
         ),
        mesh_.points()
    );

    if (debug)
    {
        simpleVTKWriter
        (
            boundaryPatch.localFaces(),
            boundaryPatch.localPoints()
        ).write("gapBoundaryFaces.vtk");
    }

    labelList markedFaces(globalRegion.nRegions(), -1);
    labelList markedFacesProcID(globalRegion.nRegions(), -1);

    forAll(boundaryPatch, i)
    {
        label facei = boundaryPatch.addressing()[i];
        label own = mesh_.faceOwner()[facei];
        label regioni = globalRegion[own];

        if (gapCells[own] && markedFaces[regioni] == -1)
        {
            markedFaces[regioni] = i;
            markedFacesProcID[regioni] = Pstream::myProcNo();
        }
    }
    Pstream::listCombineReduce(markedFacesProcID, maxOp<label>());
    DynamicList<label> seedFaces(globalRegion.nRegions());

    forAll(boundaryPatch, i)
    {
        label facei = boundaryPatch.addressing()[i];
        label own = mesh_.faceOwner()[facei];
        label regioni = globalRegion[own];

        if
        (
            markedFaces[regioni] == i
            && markedFacesProcID[regioni] == Pstream::myProcNo()
        )
        {
            seedFaces.append(i);
        }
    }

    boolList visitedFaces(boundaryPatch.size(), false);
    boolList visitedEdges(boundaryPatch.nEdges(), false);

    forAll(seedFaces, i)
    {
        label facei = seedFaces[i];
        visitedFaces[facei] = true;
        const labelList& fEdges = boundaryPatch.faceEdges()[facei];
        forAll(fEdges, fEI)
        {
            visitedEdges[fEdges[fEI]] = true;
        }
    }


    labelList boundaryMeshEdges
    (
        boundaryPatch.meshEdges(mesh_.edges(), mesh_.pointEdges())
    );

    while (true)
    {
        syncTools::syncEdgeList
        (
            mesh_,
            boundaryMeshEdges,
            visitedEdges,
            orEqOp<bool>(),
            true
        );

        label nChanged = 0;
        forAll(boundaryPatch.edges(), edgei)
        {
            if (visitedEdges[edgei])
            {
                const edge bEdge = boundaryPatch.edges()[edgei];
                forAll(bEdge, eI)
                {
                     const labelList& pFaces =
                         boundaryPatch.pointFaces()[bEdge[eI]];
                     forAll(pFaces, pFI)
                     {
                         label facei = pFaces[pFI];
                         if (!visitedFaces[facei])
                         {
                             visitedFaces[facei] = true;
                             const labelList& fEdges =
                                 boundaryPatch.faceEdges()[facei];
                             forAll(fEdges, fEI)
                             {
                                 visitedEdges[fEdges[fEI]] = true;
                             }
                             nChanged++;
                         }
                     }
                }
            }
        }

        if (returnReduce(nChanged, sumOp<label>()) == 0)
        {
            break;
        }
    }

    DynamicList<label> sidedBoundaryFaces(mesh_.nFaces()/10);
    forAll(boundaryPatch, i)
    {
        if (visitedFaces[i])
        {
            label facei = boundaryPatch.addressing()[i];
            sidedBoundaryFaces.append(facei);
        }
    }

    indirectPrimitivePatch sidedPatch
    (
        IndirectList<face>
        (
            mesh_.faces(),
            sidedBoundaryFaces
        ),
        mesh_.points()
    );

    if (debug)
    {
        simpleVTKWriter
        (
            sidedPatch.localFaces(),
            sidedPatch.localPoints()
        ).write("singleSidedBoundaryFaces.vtk");
    }

    DynamicList<label> cellsToKeep(mesh_.nCells()/100);
    forAll(sidedPatch.meshPoints(), i)
    {
        label pointi = sidedPatch.meshPoints()[i];

        const labelList& pCells = mesh_.pointCells()[pointi];
        forAll(pCells, pCI)
        {
            label celli = pCells[pCI];
            if (gapCells[celli])
            {
                gapCells[celli] = false;
                cellsToKeep.append(celli);
            }
        }
    }

    DynamicList<label> cellsToRemove(mesh_.nCells()/100);
    forAll(gapCells,celli)
    {
        if (gapCells[celli])
        {
            cellsToRemove.append(celli);
        }
    }

    //walk other side

    forAll(sidedPatch, i)
    {
        label facei = sidedPatch.addressing()[i];
        blockedFaces[facei] = true;
    }

    gapCells = false;
    forAll(cellsToKeep, i)
    {
        label celli = cellsToKeep[i];
        gapCells[celli] = true;
    }

    syncTools::swapBoundaryCellList(mesh_, gapCells, neiGapCells);

    DynamicList<label> otherSide;
    forAll(mesh_.faces(), facei)
    {
        if (blockedFaces[facei])
        {
            continue;
        }

        label patchi = patches.whichPatch(facei);
        label own = mesh_.faceOwner()[facei];
        bool gapCellOwn = gapCells[own];

        if (patchi == -1 || patches[patchi].coupled())
        {
            bool gapCellNei = -1;

            if (patchi == -1)
            {
                label nei = mesh_.faceNeighbour()[facei];
                gapCellNei = gapCells[nei];
            }
            else
            {
                gapCellNei = neiGapCells[facei-mesh_.nInternalFaces()];
            }

            if
            (
                (gapCellOwn && !gapCellNei)
                || (!gapCellOwn && gapCellNei)
            )
            {
                otherSide.append(facei);
            }
        }
        else
        {
            if (gapCellOwn)
            {
                otherSide.append(facei);
            }
        }
    }

    sidedPatch.resetAddressing(otherSide);
    sidedPatch.clearOut();

    forAll(sidedPatch.meshPoints(), i)
    {
        label pointi = sidedPatch.meshPoints()[i];

        const labelList& pCells = mesh_.pointCells()[pointi];
        forAll(pCells, pCI)
        {
            label celli = pCells[pCI];
            if (gapCells[celli])
            {
                gapCells[celli] = false;
            }
        }
    }

    if (debug)
    {
        simpleVTKWriter
        (
            sidedPatch.localFaces(),
            sidedPatch.localPoints()
        ).write("otherSidedBoundaryFaces.vtk");
    }

    forAll(gapCells,celli)
    {
        if (gapCells[celli])
        {
            cellsToRemove.append(celli);
        }
    }

    labelList ownPatch(mesh_.nFaces(), 0);
    forAll(cellsToRemove, i)
    {
        label celli = cellsToRemove[i];

        cell c = mesh_.cells()[celli];

        label cellPatchID = -1;

        forAll(c, cFI)
        {
            label facei = c[cFI];
            label patchi = patches.whichPatch(facei);
            if (patchi != -1 && !patches[patchi].coupled())
            {
                cellPatchID = patchi;
                break;
            }
        }

        if (cellPatchID != -1)
        {
            forAll(c, cFI)
            {
                label facei = c[cFI];
                ownPatch[facei] = cellPatchID;
            }
        }
    }

    // Remove cells
    removeCells cellRemover(mesh_);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatchIDs(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];

        exposedPatchIDs[i] = ownPatch[facei];
    }

    doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        cellRemover
    );

    //reset gap baffle cells for subsequent cell removal
    {
        volScalarField& gapCellsField = const_cast<volScalarField&>
            (mesh_.lookupObject<volScalarField>("gapCells"));

        boolList gapCellPts(mesh_.nPoints(), false);
        forAll(mesh_.cells(), celli)
        {
            if (gapCellsField[celli] > -1)
            {
                const labelList& cPts = mesh_.cellPoints()[celli];
                forAll(cPts, cPtI)
                {
                    gapCellPts[cPts[cPtI]] = true;
                }
            }
        }

        // Sync
        syncTools::syncPointList
        (
            mesh_,
            gapCellPts,
            orEqOp<bool>(),
            false           // null value
        );

        boolList bdyPts(mesh_.nPoints(), false);
        forAll(mesh_.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            if (patchi != -1 && !patches[patchi].coupled())
            {
                const face f = mesh_.faces()[facei];
                forAll(f,fp)
                {
                    bdyPts[f[fp]] = true;
                }
            }
        }
        // Sync
        syncTools::syncPointList
        (
            mesh_,
            bdyPts,
            orEqOp<bool>(),
            false           // null value
        );

        boolList bdyCells(mesh_.nCells(), false);
        forAll(mesh_.cells(), celli)
        {
            const labelList& cPts = mesh_.cellPoints()[celli];
            forAll(cPts, cPtI)
            {
                if (bdyPts[cPts[cPtI]])
                {
                    bdyCells[celli] = true;
                    break;
                }
            }
        }

        boolList extPts(mesh_.nPoints(), false);
        forAll(mesh_.cells(), celli)
        {
            if (!bdyCells[celli])
            {
                const labelList& cPts = mesh_.cellPoints()[celli];
                forAll(cPts, cPtI)
                {
                    extPts[cPts[cPtI]] = true;
                }
            }
        }
        // Sync
        syncTools::syncPointList
        (
            mesh_,
            extPts,
            orEqOp<bool>(),
            false           // null value
        );

        forAll(extPts, pointi)
        {
            if (extPts[pointi])
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    bdyCells[pCells[pCI]] = false;
                }
            }
        }

        forAll(gapCellPts, pointi)
        {
            if (gapCellPts[pointi])
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    label celli = pCells[pCI];
                    if (bdyCells[celli])
                    {
                        gapCellsField[celli] = 1;
                    }
                }
            }
        }
    }
}


void Foam::meshRefinement::checkZoneLeakPaths
(
    const refinementParameters& refineParams,
    const labelList& namedSurfaceIndex,
    const labelList& globalRegion1,
    const labelList& globalRegion2,
    const labelList& cellToZone,
    const writer<scalar>& leakPathFormatter
) const
{
    pointField locationsInMesh;
    wordList zonesInMesh;

    if (refineParams.zonesToTestLocations().size())
    {
        locationsInMesh = refineParams.zonesToTestLocations();
        zonesInMesh = refineParams.zonesToTestNames();
    }
    else
    {
        locationsInMesh = refineParams.locationsInMesh();
        zonesInMesh = refineParams.zonesInMesh();
    }

    if (locationsInMesh.size() == 1)
    {
        return;
    }
    else
    {
        Info<<"Calculating zone leak paths "<<endl;
    }

    bool namedSet = false;
    if (namedSurfaceIndex.size() == mesh_.nFaces())
    {
        namedSet = true;
    }

    boolList blockedFace(mesh_.nFaces(), false);
    forAll(mesh_.faces(), facei)
    {
        if
        (
            globalRegion1[facei] != -1 || globalRegion2[facei] != -1
            || (namedSet && namedSurfaceIndex[facei] != -1)
        )
        {
            blockedFace[facei] = true;
        }
    }

    if (debug)
    {
        mesh_.write();
        faceSet intersectedFaces(mesh_, "intersectedFaces", mesh_.nFaces()/100);

        forAll(blockedFace, facei)
        {
            if (blockedFace[facei])
            {
                intersectedFaces.insert(facei);
            }
        }
        intersectedFaces.instance() = timeName();
        Pout<< "Dumping " << intersectedFaces.size()
            << " intersected faces to "
            << intersectedFaces.objectPath() << endl;
        intersectedFaces.write();
    }

    regionSplit cellRegion(mesh_, blockedFace);

    labelList insideRegions(locationsInMesh.size());
    forAll(insideRegions, i)
    {
        // Find the region containing the keepPoint
        label regioni = -1;

        label celli = findCell
        (
            locationsInMesh[i],
            mesh_,
            meshCutter_
        );

        if (celli != -1)
        {
            regioni = cellRegion[celli];
        }
        reduce(regioni, maxOp<label>());

        insideRegions[i] = regioni;
    }
    label nLeaks = 0;

    boolList alreadySet(locationsInMesh.size(), false);

    forAll(locationsInMesh, loci)
    {
        const word zoneName = zonesInMesh[loci];

        if (zoneName != "none")
        {
            alreadySet[loci] = true;
            forAll(locationsInMesh, locj)
            {
                if
                (
                    !alreadySet[locj]
                    && (insideRegions[loci] == insideRegions[locj])
                )
                {
                    word zoneTag = zoneName
                        + "_to_" + zonesInMesh[locj];
                    calculateLeakPath
                    (
                        "zoneLeakPath",
                        zoneTag,
                        mesh_,
                        blockedFace,
                        pointField(1,locationsInMesh[loci]),
                        pointField(1,locationsInMesh[locj]),
                        leakPathFormatter
                    );
                    nLeaks++;
                }
            }
        }
    }

    if (returnReduce(nLeaks, sumOp<label>()) > 0)
    {
        FatalErrorInFunction
            << "Zonal leak detected and exiting. Check output leak paths."
            << exit(FatalError);
    }

    return;
}

// ************************************************************************* //
