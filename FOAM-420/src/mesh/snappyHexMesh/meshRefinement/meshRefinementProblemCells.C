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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshRefinement/meshRefinement.H"
#include "fvMesh/fvMesh.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/Time/Time.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "sets/topoSets/pointSet.H"
#include "sets/topoSets/faceSet.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "sets/topoSets/cellSet.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "motionSmoother/polyMeshGeometry/polyMeshGeometry.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "global/unitConversion/unitConversion.H"
#include "snappyHexMeshDriver/snappySnapDriver.H"
#include "snappyHexMeshDriver/snapParameters/snapParameters.H"
#include "motionSmoother/motionSmoother.H"
#include "meshStructure/topoDistanceData.H"
#include "algorithms/FaceCellWave/FaceCellWave.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "snappyHexMeshDriver/layerParameters/layerParameters.H"
#include "edgeClassification/edgeClassification.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::meshRefinement::markBoundaryFace
(
    const label facei,
    boolList& isBoundaryFace,
    boolList& isBoundaryEdge,
    boolList& isBoundaryPoint
) const
{
    if (isBoundaryFace[facei])
    {
        return false;
    }
    else
    {
        isBoundaryFace[facei] = true;

        const labelList& fEdges = mesh_.faceEdges(facei);

        forAll(fEdges, fp)
        {
            isBoundaryEdge[fEdges[fp]] = true;
        }

        const face& f = mesh_.faces()[facei];

        forAll(f, fp)
        {
            isBoundaryPoint[f[fp]] = true;
        }

        return true;
    }
}


void Foam::meshRefinement::findNearest
(
    const labelList& meshFaces,
    List<pointIndexHit>& nearestInfo,
    labelList& nearestSurface,
    labelList& nearestRegion,
    vectorField& nearestNormal
) const
{
    pointField fc(meshFaces.size());
    forAll(meshFaces, i)
    {
        fc[i] = mesh_.faceCentres()[meshFaces[i]];
    }

    const labelList allSurfaces(identity(surfaces_.surfaces().size()));

    surfaces_.findNearest
    (
        allSurfaces,
        fc,
        scalarField(fc.size(), sqr(GREAT)),    // sqr of attraction
        nearestSurface,
        nearestInfo
    );

    // Do normal testing per surface.
    nearestNormal.setSize(nearestInfo.size());
    nearestRegion.setSize(nearestInfo.size());

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

        labelList localRegion;
        surfaces_.geometry()[geomI].getRegion(localHits, localRegion);

        label localI = 0;
        forAll(nearestSurface, i)
        {
            if (nearestSurface[i] == surfI)
            {
                nearestNormal[i] = localNormals[localI];
                nearestRegion[i] = localRegion[localI];
                localI++;
            }
        }
    }
}


Foam::Map<Foam::label> Foam::meshRefinement::findEdgeConnectedProblemCells
(
    const scalarField& perpendicularAngle,
    const labelList& globalToPatch
) const
{
    // Construct addressing engine from all patches added for meshing.
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh_,
            meshedPatches()
        )
    );
    const indirectPrimitivePatch& pp = ppPtr();


    // 1. Collect faces to test
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<label> candidateFaces(pp.size()/20);

    const labelListList& edgeFaces = pp.edgeFaces();

    const labelList& cellLevel = meshCutter_.cellLevel();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() == 2)
        {
            label face0 = pp.addressing()[eFaces[0]];
            label face1 = pp.addressing()[eFaces[1]];

            label cell0 = mesh_.faceOwner()[face0];
            label cell1 = mesh_.faceOwner()[face1];

            if (cellLevel[cell0] > cellLevel[cell1])
            {
                // cell0 smaller.
                const vector& n0 = pp.faceNormals()[eFaces[0]];
                const vector& n1 = pp.faceNormals()[eFaces[1]];

                if (mag(n0 & n1) < 0.1)
                {
                    candidateFaces.append(face0);
                }
            }
            else if (cellLevel[cell1] > cellLevel[cell0])
            {
                // cell1 smaller.
                const vector& n0 = pp.faceNormals()[eFaces[0]];
                const vector& n1 = pp.faceNormals()[eFaces[1]];

                if (mag(n0 & n1) < 0.1)
                {
                    candidateFaces.append(face1);
                }
            }
        }
    }
    candidateFaces.shrink();

    Info<< "Testing " << returnReduce(candidateFaces.size(), sumOp<label>())
        << " faces on edge-connected cells of differing level."
        << endl;

    if (debug&meshRefinement::MESH)
    {
        faceSet fSet
        (
            mesh_,
            "edgeConnectedFaces",
            labelHashSet(candidateFaces)
        );
        fSet.instance() = timeName();
        Pout<< "Writing " << fSet.size()
            << " with problematic topology to faceSet "
            << fSet.objectPath() << endl;
        fSet.write();
    }


    // 2. Find nearest surface on candidate faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    List<pointIndexHit> nearestInfo;
    labelList nearestSurface;
    labelList nearestRegion;
    vectorField nearestNormal;
    findNearest
    (
        candidateFaces,
        nearestInfo,
        nearestSurface,
        nearestRegion,
        nearestNormal
    );


    // 3. Test angle to surface
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    Map<label> candidateCells(candidateFaces.size());

    faceSet perpFaces(mesh_, "perpendicularFaces", labelHashSet(pp.size()/100));

    forAll(candidateFaces, i)
    {
        label facei = candidateFaces[i];

        vector n = mesh_.faceAreas()[facei];
        n /= mag(n);

        label region = surfaces_.globalRegion
        (
            nearestSurface[i],
            nearestRegion[i]
        );

        scalar angle = degToRad(perpendicularAngle[region]);

        if (angle >= 0)
        {
            if (mag(n & nearestNormal[i]) < Foam::sin(angle))
            {
                perpFaces.insert(facei);
                candidateCells.insert
                (
                    mesh_.faceOwner()[facei],
                    globalToPatch[region]
                );
            }
        }
    }

    if (debug&meshRefinement::MESH)
    {
        perpFaces.instance() = timeName();
        Pout<< "Writing " << perpFaces.size()
            << " faces that are perpendicular to the surface to set "
            << perpFaces.objectPath() << endl;
        perpFaces.write();
    }
    return candidateCells;
}


// Check if moving face to new points causes a 'collapsed' face.
// Uses new point position only for the face, not the neighbouring
// cell centres
bool Foam::meshRefinement::isCollapsedFace
(
    const pointField& points,
    const pointField& neiCc,
    const scalar minFaceArea,
    const scalar maxNonOrtho,
    const label facei
) const
{
    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(degToRad(maxNonOrtho));

    vector s = mesh_.faces()[facei].areaNormal(points);
    scalar magS = mag(s);

    // Check face area
    if (magS < minFaceArea)
    {
        return true;
    }

    // Check orthogonality
    const point& ownCc = mesh_.cellCentres()[mesh_.faceOwner()[facei]];

    if (mesh_.isInternalFace(facei))
    {
        label nei = mesh_.faceNeighbour()[facei];
        vector d = mesh_.cellCentres()[nei] - ownCc;

        scalar dDotS = (d & s)/(mag(d)*magS + VSMALL);

        if (dDotS < severeNonorthogonalityThreshold)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        label patchi = mesh_.boundaryMesh().whichPatch(facei);

        if (mesh_.boundaryMesh()[patchi].coupled())
        {
            vector d = neiCc[facei-mesh_.nInternalFaces()] - ownCc;

            scalar dDotS = (d & s)/(mag(d)*magS + VSMALL);

            if (dDotS < severeNonorthogonalityThreshold)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            // Collapsing normal boundary face does not cause problems with
            // non-orthogonality
            return false;
        }
    }
}


// Check if moving cell to new points causes it to collapse.
bool Foam::meshRefinement::isCollapsedCell
(
    const pointField& points,
    const scalar volFraction,
    const label celli
) const
{
    scalar vol = mesh_.cells()[celli].mag(points, mesh_.faces());

    if (vol/mesh_.cellVolumes()[celli] < volFraction)
    {
        return true;
    }
    else
    {
        return false;
    }
}


Foam::labelList Foam::meshRefinement::nearestPatch
(
    const labelList& adaptPatchIDs
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelList nearestAdaptPatch;

    if (adaptPatchIDs.size())
    {
        nearestAdaptPatch.setSize(mesh_.nFaces(), adaptPatchIDs[0]);


        // Count number of faces in adaptPatchIDs
        label nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            const polyPatch& pp = patches[adaptPatchIDs[i]];
            nFaces += pp.size();
        }

        // Field on cells and faces.
        List<topoDistanceData> cellData(mesh_.nCells());
        List<topoDistanceData> faceData(mesh_.nFaces());

        // Start of changes
        labelList patchFaces(nFaces);
        List<topoDistanceData> patchData(nFaces);
        nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            label patchi = adaptPatchIDs[i];
            const polyPatch& pp = patches[patchi];

            forAll(pp, i)
            {
                patchFaces[nFaces] = pp.start()+i;
                patchData[nFaces] = topoDistanceData(patchi, 0);
                nFaces++;
            }
        }

        // Propagate information inwards
        FaceCellWave<topoDistanceData> deltaCalc
        (
            mesh_,
            patchFaces,
            patchData,
            faceData,
            cellData,
            mesh_.globalData().nTotalCells()+1
        );

        // And extract

        bool haveWarned = false;
        forAll(faceData, facei)
        {
            if (!faceData[facei].valid(deltaCalc.data()))
            {
                if (!haveWarned)
                {
                    WarningInFunction
                        << "Did not visit some faces, e.g. face " << facei
                        << " at " << mesh_.faceCentres()[facei] << endl
                        << "Assigning  these cells to patch "
                        << adaptPatchIDs[0]
                        << endl;
                    haveWarned = true;
                }
            }
            else
            {
                nearestAdaptPatch[facei] = faceData[facei].data();
            }
        }
    }
    else
    {
        // Use patch 0
        nearestAdaptPatch.setSize(mesh_.nFaces(), 0);
    }

    return nearestAdaptPatch;
}


// Returns list with for every internal face -1 or the patch they should
// be baffled into. Gets run after createBaffles so all the unzoned surface
// intersections have already been turned into baffles. (Note: zoned surfaces
// are not baffled at this stage)
// Used to remove cells by baffling all their faces and have the
// splitMeshRegions chuck away non used regions.
Foam::labelList Foam::meshRefinement::markFacesOnProblemCells
(
    const dictionary& motionDict,
    const refinementParameters& refineParams,
    const bool removeEdgeConnectedCells,
    const scalarField& perpendicularAngle,
    const labelList& globalToPatch,
    const labelList& adaptPatchIDs,
    const labelList& checkSurfaces,
    const List<labelPair> removeFaces
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Mark all points and edges on baffle patches (so not on any inlets,
    // outlets etc.)
    boolList isBoundaryPoint(mesh_.nPoints(), false);
    boolList isBoundaryEdge(mesh_.nEdges(), false);
    boolList isBoundaryFace(mesh_.nFaces(), false);

    boolList removeFromChecks(mesh_.nFaces(), false);
    if (removeFaces.size())
    {
        forAll(removeFaces, i)
        {
            label face0 = removeFaces[i].first();
            label face1 = removeFaces[i].second();

            removeFromChecks[face0] = true;
            removeFromChecks[face1] = true;
        }
    }

    forAll(adaptPatchIDs, i)
    {
        const polyPatch& pp = patches[adaptPatchIDs[i]];

        label facei = pp.start();

        forAll(pp, j)
        {
            if (!removeFromChecks[facei])
            {
                markBoundaryFace
                (
                    facei,
                    isBoundaryFace,
                    isBoundaryEdge,
                    isBoundaryPoint
                 );
            }
            facei++;
        }
    }

    // Per internal face (boundary faces not used) the patch that the
    // baffle should get (or -1)
    labelList facePatch(mesh_.nFaces(), -1);

    // Swap neighbouring cell centres and cell level
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc);

    // Count of faces marked for baffling
    label nBaffleFaces = 0;
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    // Count of faces not baffled since would not cause a collapse
    label nPrevented = 0;

    if (removeEdgeConnectedCells && max(perpendicularAngle) >= 0)
    {
        Info<< "markFacesOnProblemCells :"
            << " Checking for edge-connected cells of highly differing sizes."
            << endl;

        // Pick up the cells that need to be removed and (a guess for)
        // the patch they should be patched with.
        Map<label> problemCells
        (
            findEdgeConnectedProblemCells
            (
                perpendicularAngle,
                globalToPatch
            )
        );

        // Baffle all faces of cells that need to be removed
        forAllConstIter(Map<label>, problemCells, iter)
        {
            const cell& cFaces = mesh_.cells()[iter.key()];

            forAll(cFaces, i)
            {
                label facei = cFaces[i];

                if (facePatch[facei] == -1 && mesh_.isInternalFace(facei))
                {
                    facePatch[facei] = getBafflePatch
                    (
                        adaptPatchIDs,
                        facePatch,
                        facei
                    );

                    nBaffleFaces++;

                    // Mark face as a 'boundary'
                    markBoundaryFace
                    (
                        facei,
                        isBoundaryFace,
                        isBoundaryEdge,
                        isBoundaryPoint
                    );
                }
                else if (facePatch[facei] == -1)
                {
                    label patchi = mesh_.boundaryMesh().whichPatch(facei);
                    if (mesh_.boundaryMesh()[patchi].coupled())
                    {
                        facePatch[facei] = getBafflePatch
                        (
                            adaptPatchIDs,
                            facePatch,
                            facei
                         );
                        nBaffleFaces++;

                        // Mark face as a 'boundary'
                        markBoundaryFace
                        (
                            facei,
                            isBoundaryFace,
                            isBoundaryEdge,
                            isBoundaryPoint
                         );
                    }
                }
            }
        }
        Info<< "markFacesOnProblemCells : Marked "
            << returnReduce(nBaffleFaces, sumOp<label>())
            << " additional internal faces to be converted into baffles"
            << " due to "
            << returnReduce(problemCells.size(), sumOp<label>())
            << " cells edge-connected to lower level cells." << endl;

        if (debug&meshRefinement::MESH)
        {
            cellSet problemCellSet
            (
                mesh_,
                "problemCells",
                labelHashSet(problemCells.toc())
            );
            problemCellSet.instance() = timeName();
            Pout<< "Writing " << problemCellSet.size()
                << " cells that are edge connected to coarser cell to set "
                << problemCellSet.objectPath() << endl;
            problemCellSet.write();
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

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>()
    );

    // See if checking for collapse
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Collapse checking parameters
    const scalar volFraction =
        motionDict.lookupOrDefault<scalar>("minVolCollapseRatio", -1);

    const bool checkCollapse = (volFraction >= 0);
    scalar minArea = -1;
    scalar maxNonOrtho = -1;

    // Find nearest (non-baffle) surface
    pointField newPoints,surfPointNormals;

    if (checkCollapse || refineParams.splitCells())
    {
        if (checkCollapse)
        {
            minArea = readScalar(motionDict.lookup("minArea"));
            maxNonOrtho = readScalar(motionDict.lookup("maxNonOrtho"));

            Info<< "markFacesOnProblemCells :"
                << " Deleting all-anchor surface cells only if"
                << " snapping them violates mesh quality constraints:" << nl
                << "    snapped/original cell volume < " << volFraction << nl
                << "    face area                    < " << minArea << nl
                << "    non-orthogonality            > " << maxNonOrtho << nl
                << endl;
        }

        // Construct addressing engine.
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh_,
                adaptPatchIDs
            )
        );
        const indirectPrimitivePatch& pp = ppPtr();
        const pointField& localPoints = pp.localPoints();
        const labelList& meshPoints = pp.meshPoints();

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces_.findNearest
        (
            checkSurfaces,
            localPoints,
            scalarField(localPoints.size(), sqr(GREAT)),    // sqr of attraction
            hitSurface,
            hitInfo
        );

        // Start off from current points
        newPoints = mesh_.points();
        surfPointNormals = vectorField(mesh_.points().size(), vector::max);

        forAll(checkSurfaces, sI)
        {
            label surfI = checkSurfaces[sI];
            DynamicList<pointIndexHit> localHits;
            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    localHits.append(hitInfo[i]);
                }
            }
            localHits.shrink();

            pointField localNormals;
            label geomI = surfaces_.surfaces()[surfI];
            surfaces_.geometry()[geomI].getNormal(localHits, localNormals);

            label localI = 0;
            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    surfPointNormals[meshPoints[i]] = localNormals[localI];
                    localI++;
                }
            }
        }
        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                newPoints[meshPoints[i]] = hitInfo[i].hitPoint();
            }
        }

        if (debug&meshRefinement::MESH)
        {
            const_cast<Time&>(mesh_.time())++;
            pointField oldPoints(mesh_.points());
            mesh_.movePoints(newPoints);
            Pout<< "Writing newPoints mesh to time " << timeName()
                << endl;
            write
            (
                debugType(debug),
                writeType(writeLevel() | WRITEMESH),
                mesh_.time().path()/"newPoints"
            );
            mesh_.movePoints(oldPoints);
        }
    }

    DynamicList<label> recheckCells(mesh_.nCells()/10000);

    // On-the-fly addressing storage.
    DynamicList<label> dynFEdges;
    DynamicList<label> dynCPoints;

    if (refineParams.splitCells())
    {
        const labelList& faceOwn = mesh_.faceOwner();
        boolList checkConn(mesh_.nCells(), false);

        label nCheckLoops = 0;
        bool checkAll = true;

        while (true)
        {
            label nMarked = 0;

            forAll(cellLevel, cellI)
            {
                const cell& cFaces = mesh_.cells()[cellI];

                // Count boundary faces.
                label nBfaces = 0;
                scalar cellPatchArea = 0.;
                scalar cellArea = 0.;
                forAll(cFaces, cFaceI)
                {
                    scalar fA = mesh_.magFaceAreas()[cFaces[cFaceI]];
                    cellArea += fA;

                    if (isBoundaryFace[cFaces[cFaceI]])
                    {
                        cellPatchArea += fA;
                        nBfaces++;
                    }
                }

                if
                (
                    cFaces.size() == nBfaces + 1 ||
                    (cellPatchArea/cellArea > 0.6 && checkAll)
                )
                {
                    // Block all faces of cell
                    forAll(cFaces, cf)
                    {
                        label faceI = cFaces[cf];

                        if
                        (
                            facePatch[faceI] == -1
                            && mesh_.isInternalFace(faceI)
                         )
                        {
                            facePatch[faceI] = getBafflePatch
                            (
                                adaptPatchIDs,
                                facePatch,
                                faceI
                            );
                            nBaffleFaces++;

                            // Mark face as a 'boundary'
                            if
                            (
                                markBoundaryFace
                                (
                                    faceI,
                                    isBoundaryFace,
                                    isBoundaryEdge,
                                    isBoundaryPoint
                                )
                            )
                            {
                                nMarked++;
                            }
                        }
                        else if (facePatch[faceI] == -1)
                        {
                            label patchI =
                                mesh_.boundaryMesh().whichPatch(faceI);
                            if (mesh_.boundaryMesh()[patchI].coupled())
                            {
                                facePatch[faceI] = getBafflePatch
                                (
                                    adaptPatchIDs,
                                    facePatch,
                                    faceI
                                 );
                                nBaffleFaces++;

                                // Mark face as a 'boundary'
                                if
                                (
                                    markBoundaryFace
                                    (
                                        faceI,
                                        isBoundaryFace,
                                        isBoundaryEdge,
                                        isBoundaryPoint
                                    )
                                )
                                {
                                    nMarked++;
                                }
                            }
                        }
                    }
                }
                const labelList& cPoints = mesh_.cellPoints(cellI, dynCPoints);
                label nBpts = 0;
                forAll(cPoints, i)
                {
                    label pointI = cPoints[i];

                    if (isBoundaryPoint[pointI])
                    {
                        nBpts++;
                    }
                }

                bool removeCell = false;
                if (cFaces.size() == 4)
                {
                    cell c = mesh_.cells()[cellI];
                    bool tet = true;
                    forAll(c, cFI)
                    {
                        face f = mesh_.faces()[c[cFI]];
                        if (f.size() != 3)
                        {
                            tet = false;
                            break;
                        }
                    }

                    if (tet)
                    {
                        labelHashSet cFaces(c);
                        labelList cPoints = c.labels(mesh_.faces());
                        DynamicList<label> oppositePoints(cPoints.size());

                        forAll(cPoints, cPtI)
                        {
                            labelList pFaces = mesh_.pointFaces()[cPoints[cPtI]];
                            label nFound = 0;
                            forAll(pFaces, pFI)
                            {
                                label faceI = pFaces[pFI];

                                if (cFaces.found(faceI) && isBoundaryFace[faceI])
                                {
                                    nFound++;
                                }
                            }

                            if (nFound == 1)
                            {
                                oppositePoints.append(cPoints[cPtI]);
                            }
                        }
                        oppositePoints.shrink();

                        if (oppositePoints.size() == 2)
                        {
                            label pt0 = oppositePoints[0];
                            label pt1 = oppositePoints[1];

                            if
                            (
                                surfPointNormals[pt0] != vector::max
                                && surfPointNormals[pt1] != vector::max
                            )
                            {
                                vector norm0 = surfPointNormals[pt0];
                                norm0 /= (mag(norm0) + SMALL);
                                vector norm1 = surfPointNormals[pt1];
                                norm1 /= (mag(norm1) + SMALL);

                                if (mag(norm0 & norm1) > 0.642)
                                {
                                    removeCell = true;
                                }
                            }
                            else
                            {
                                removeCell = true;
                            }
                        }
                    }
                }

                if (cPoints.size() == nBpts)
                {
                    if
                    (
                        !removeCell &&
                        !isCollapsedCell(newPoints, 0.25, cellI)
                     )
                    {
                            //Do nothing
                    }
                    else
                    {
                        // Block all faces of cell
                        forAll(cFaces, cf)
                        {
                            label faceI = cFaces[cf];

                            if
                            (
                                facePatch[faceI] == -1
                                && mesh_.isInternalFace(faceI)
                            )
                            {
                                facePatch[faceI] = getBafflePatch
                                (
                                    adaptPatchIDs,
                                    facePatch,
                                    faceI
                                );

                                nBaffleFaces++;

                                if
                                (
                                    // Mark face as a 'boundary'
                                    markBoundaryFace
                                    (
                                        faceI,
                                        isBoundaryFace,
                                        isBoundaryEdge,
                                        isBoundaryPoint
                                    )
                                )
                                {
                                    nMarked++;
                                }
                            }
                            else if (facePatch[faceI] == -1)
                            {
                                label patchI =
                                    mesh_.boundaryMesh().whichPatch(faceI);
                                if (mesh_.boundaryMesh()[patchI].coupled())
                                {
                                    facePatch[faceI] = getBafflePatch
                                    (
                                        adaptPatchIDs,
                                        facePatch,
                                        faceI
                                    );
                                    nBaffleFaces++;

                                    if
                                    (
                                        // Mark face as a 'boundary'
                                        markBoundaryFace
                                        (
                                            faceI,
                                            isBoundaryFace,
                                            isBoundaryEdge,
                                            isBoundaryPoint
                                        )
                                    )
                                    {
                                        nMarked++;
                                    }
                                }
                            }
                        }
                    }
                }
                else if (cPoints.size() == nBpts + 1)
                {
                    bool tet = false;
                    if (cFaces.size() == 4)
                    {
                        cell c = mesh_.cells()[cellI];
                        tet = true;
                        forAll(c, cFI)
                        {
                            face f = mesh_.faces()[c[cFI]];
                            if (f.size() != 3)
                            {
                                tet = false;
                                break;
                            }
                        }
                    }

                    if
                    (
                        (checkCollapse || tet)
                        && !isCollapsedCell(newPoints, 0.25, cellI)
                    )
                    {
                        //Do nothing
                    }
                    else if (tet)
                    {
                        // Block all faces of cell
                        forAll(cFaces, cf)
                        {
                            label faceI = cFaces[cf];

                            if
                            (
                                facePatch[faceI] == -1
                                && mesh_.isInternalFace(faceI)
                            )
                            {
                                facePatch[faceI] = getBafflePatch
                                (
                                    adaptPatchIDs,
                                    facePatch,
                                    faceI
                                );
                                nBaffleFaces++;

                                // Mark face as a 'boundary'
                                if
                                (
                                    markBoundaryFace
                                    (
                                        faceI,
                                        isBoundaryFace,
                                        isBoundaryEdge,
                                        isBoundaryPoint
                                    )
                                )
                                {
                                    nMarked++;
                                }
                            }
                            else if (facePatch[faceI] == -1)
                            {
                                label patchI =
                                    mesh_.boundaryMesh().whichPatch
                                    (faceI);
                                if (mesh_.boundaryMesh()[patchI].coupled())
                                {
                                    facePatch[faceI] = getBafflePatch
                                    (
                                        adaptPatchIDs,
                                        facePatch,
                                        faceI
                                    );
                                    nBaffleFaces++;

                                    // Mark face as a 'boundary'
                                    if
                                    (
                                        markBoundaryFace
                                        (
                                            faceI,
                                            isBoundaryFace,
                                            isBoundaryEdge,
                                            isBoundaryPoint
                                        )
                                    )
                                    {
                                        nMarked++;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (checkAll)
            {
                checkConn = false;

                forAll(cellLevel, cellI)
                {
                    const cell& c = mesh_.cells()[cellI];

                    label nBFaces = 0;
                    forAll(c, cFaceI)
                    {
                        if (isBoundaryFace[c[cFaceI]])
                        {
                            nBFaces++;
                        }
                    }


                    if ((c.size() - nBFaces)<= 2)
                    {
                        labelList cPoints = c.labels(mesh_.faces());
                        label nBpts = 0;
                        forAll(cPoints, i)
                        {
                            label pointI = cPoints[i];

                            if (isBoundaryPoint[pointI])
                            {
                                nBpts++;
                            }
                        }

                        if (nBpts == cPoints.size())
                        {
                            checkConn[cellI] = true;
                        }
                    }
                }

                // Calculate coupled blocked cells
                List<bool> neiCheckConn
                    (mesh_.nFaces()-mesh_.nInternalFaces(), false);
                for
                (
                    label facei = mesh_.nInternalFaces();
                    facei < mesh_.nFaces();
                    facei++
                 )
                {
                    neiCheckConn[facei-mesh_.nInternalFaces()] =
                        checkConn[faceOwn[facei]];
                }
                syncTools::swapBoundaryFaceList(mesh_, neiCheckConn);

                forAll(cellLevel, cellI)
                {
                    bool blockCell = false;
                    if (checkConn[cellI])
                    {
                        const cell& c = mesh_.cells()[cellI];

                        forAll(c, cFaceI)
                        {
                            label facei = c[cFaceI];
                            if (!isBoundaryFace[facei])
                            {
                                label own = mesh_.faceOwner()[facei];
                                if (mesh_.isInternalFace(facei))
                                {
                                    label nei = mesh_.faceNeighbour()[facei];
                                    if
                                    (
                                        own == cellI
                                        ? checkConn[nei] : checkConn[own]
                                    )
                                    {
                                        blockCell = true;
                                        break;
                                    }
                                }
                                else
                                {
                                    label patchi =
                                        mesh_.boundaryMesh().whichPatch(facei);

                                    if (mesh_.boundaryMesh()[patchi].coupled())
                                    {
                                        label bfacei = facei
                                            -mesh_.nInternalFaces();

                                        if (neiCheckConn[bfacei])
                                        {
                                            blockCell = true;
                                            break;
                                        }
                                    }
                                }
                            }
                        }

                        if (blockCell)
                        {
                            // Block all faces of cell
                            forAll(c, cf)
                            {
                                label faceI = c[cf];

                                if
                                (
                                    facePatch[faceI] == -1
                                    && mesh_.isInternalFace(faceI)
                                 )
                                {
                                    facePatch[faceI] = getBafflePatch
                                    (
                                        adaptPatchIDs,
                                        facePatch,
                                        faceI
                                     );
                                    nBaffleFaces++;

                                    // Mark face as a 'boundary'
                                    if
                                    (
                                        markBoundaryFace
                                        (
                                            faceI,
                                            isBoundaryFace,
                                            isBoundaryEdge,
                                            isBoundaryPoint
                                         )
                                     )
                                    {
                                        nMarked++;
                                    }
                                }
                                else if (facePatch[faceI] == -1)
                                {
                                    label patchI =
                                        mesh_.boundaryMesh().whichPatch
                                        (faceI);
                                    if (mesh_.boundaryMesh()[patchI].coupled())
                                    {
                                        facePatch[faceI] = getBafflePatch
                                        (
                                            adaptPatchIDs,
                                            facePatch,
                                            faceI
                                         );
                                        nBaffleFaces++;

                                        // Mark face as a 'boundary'
                                        if
                                        (
                                            markBoundaryFace
                                            (
                                                faceI,
                                                isBoundaryFace,
                                                isBoundaryEdge,
                                                isBoundaryPoint
                                             )
                                         )
                                        {
                                            nMarked++;
                                        }
                                    }
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
            else
            {
                if (nCheckLoops > 10)
                {
                    checkAll = false;
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

                syncTools::syncFaceList
                (
                    mesh_,
                    isBoundaryFace,
                    orEqOp<bool>()
                );
                nCheckLoops++;
            }
        }
    }
    else
    {
        // For each cell count the number of anchor points that are on
        // the boundary:
        // 8 : check the number of (baffle) boundary faces. If 3 or more block
        //     off the cell since the cell would get squeezed down to a diamond
        //     (probably; if the 3 or more faces are unrefined (only use the
        //      anchor points))
        // 7 : store. Used to check later on whether there are points with
        //     3 or more of these cells. (note that on a flat surface a boundary
        //     point will only have 4 cells connected to it)

        // Does cell have exactly 7 of its 8 anchor points on the boundary?
        PackedBoolList hasSevenBoundaryAnchorPoints(mesh_.nCells());
        // If so what is the remaining non-boundary anchor point?
        labelHashSet nonBoundaryAnchors(mesh_.nCells()/10000);

        forAll(cellLevel, celli)
        {
            const labelList& cPoints = mesh_.cellPoints(celli, dynCPoints);

            // Get number of anchor points (pointLevel <= cellLevel)

            label nBoundaryAnchors = 0;
//            label nNonAnchorBoundary = 0;
            label nonBoundaryAnchor = -1;

            forAll(cPoints, i)
            {
                label pointi = cPoints[i];

                if (pointLevel[pointi] <= cellLevel[celli])
                {
                    // Anchor point
                    if (isBoundaryPoint[pointi])
                    {
                        nBoundaryAnchors++;
                    }
                    else
                    {
                        // Anchor point which is not on the surface
                        nonBoundaryAnchor = pointi;
                    }
                }
                else if (isBoundaryPoint[pointi])
                {
//                    nNonAnchorBoundary++;
                }
            }

            if (nBoundaryAnchors == 8)
            {
                const cell& cFaces = mesh_.cells()[celli];
                /*
                // Count boundary faces.
                label nBfaces = 0;

                forAll(cFaces, cFacei)
                {
                    if (isBoundaryFace[cFaces[cFacei]])
                    {
                        nBfaces++;
                    }
                }
                */

                // If nBfaces > 1 make all non-boundary non-baffle faces baffles.
                // We assume that this situation is where there is a single
                // cell sticking out which would get flattened.

                // Eugene: delete cell no matter what.
                //if (nBfaces > 1)
                {
                    if
                    (
                        checkCollapse
                        && !isCollapsedCell(newPoints, volFraction, celli)
                    )
                    {
                        recheckCells.append(celli);
                        nPrevented++;
                        //Pout<< "Preventing baffling/removal of 8 anchor point"
                        // << " cell "
                        // << cellI << " at " << mesh_.cellCentres()[cellI]
                        // << " since new volume "
                        // << mesh_.cells()[celli].mag(newPoints, mesh_.faces())
                        // << " old volume " << mesh_.cellVolumes()[celli]
                        // << endl;
                    }
                    else
                    {
                        // Block all faces of cell
                        forAll(cFaces, cf)
                        {
                            label facei = cFaces[cf];

                            if
                            (
                                facePatch[facei] == -1
                                && mesh_.isInternalFace(facei)
                            )
                            {
                                facePatch[facei] = getBafflePatch
                                (
                                    adaptPatchIDs,
                                    facePatch,
                                    facei
                                );
                                nBaffleFaces++;

                                // Mark face as a 'boundary'
                                markBoundaryFace
                                (
                                    facei,
                                    isBoundaryFace,
                                    isBoundaryEdge,
                                    isBoundaryPoint
                                );
                            }
                            else if (facePatch[facei] == -1)
                            {
                                label patchI =
                                    mesh_.boundaryMesh().whichPatch(facei);
                                if (mesh_.boundaryMesh()[patchI].coupled())
                                {
                                    facePatch[facei] = getBafflePatch
                                    (
                                        adaptPatchIDs,
                                        facePatch,
                                        facei
                                    );
                                    nBaffleFaces++;

                                    // Mark face as a 'boundary'
                                    markBoundaryFace
                                    (
                                        facei,
                                        isBoundaryFace,
                                        isBoundaryEdge,
                                        isBoundaryPoint
                                    );
                                }
                            }
                        }
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
                        if
                        (
                            false //disable check as want to baffle
                            && !isCollapsedCell(newPoints, volFraction, celli)
                         )
                        {
                            recheckCells.append(celli);
                            nPrevented++;
                            //Pout<< "Preventing baffling of 7 anchor cell "
                            //    << celli
                            //    << " at " << mesh_.cellCentres()[celli]
                            //    << " since new volume "
                            //    << mesh_.cells()[celli].mag
                            //        (newPoints, mesh_.faces())
                            //    << " old volume " << mesh_.cellVolumes()[celli]
                            //    << endl;
                        }
                        else
                        {
                            const cell& cFaces = mesh_.cells()[celli];

                            forAll(cFaces, cf)
                            {
                                label facei = cFaces[cf];

                                if
                                (
                                    facePatch[facei] == -1
                                    && mesh_.isInternalFace(facei)
                                 )
                                {
                                    facePatch[facei] = getBafflePatch
                                    (
                                        adaptPatchIDs,
                                        facePatch,
                                        facei
                                     );
                                    nBaffleFaces++;

                                    // Mark face as a 'boundary'
                                    markBoundaryFace
                                    (
                                        facei,
                                        isBoundaryFace,
                                        isBoundaryEdge,
                                        isBoundaryPoint
                                     );
                                }
                                else if (facePatch[facei] == -1)
                                {
                                    label patchi = mesh_.boundaryMesh().whichPatch(facei);
                                    if (mesh_.boundaryMesh()[patchi].coupled())
                                    {
                                        facePatch[facei] = getBafflePatch
                                        (
                                            adaptPatchIDs,
                                            facePatch,
                                            facei
                                         );
                                        nBaffleFaces++;

                                        // Mark face as a 'boundary'
                                        markBoundaryFace
                                        (
                                            facei,
                                            isBoundaryFace,
                                            isBoundaryEdge,
                                            isBoundaryPoint
                                         );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    //Check cc-new cc and if it has moved through
    // a surface then baffle
    recheckCells.shrink();

    pointField newCC(recheckCells.size());
    pointField origCC(recheckCells.size());

    forAll(recheckCells, i)
    {
        label cellI = recheckCells[i];
        newCC[i] = mesh_.cells()[cellI].centre(newPoints, mesh_.faces());
        origCC[i] = mesh_.cellCentres()[cellI];
    }

    labelList surface1;
    labelList surface2;
    {
        List<pointIndexHit> hit1;
        labelList region1;
        List<pointIndexHit> hit2;
        labelList region2;
        surfaces_.findNearestIntersection
        (
            checkSurfaces,
            newCC,
            origCC,
            surface1,
            hit1,
            region1,
            surface2,
            hit2,
            region2
         );
    }

    forAll(recheckCells, i)
    {
        label cellI = recheckCells[i];
        if (surface1[i] != -1)
        {
            const cell& cFaces = mesh_.cells()[cellI];

            forAll(cFaces, cf)
            {
                label faceI = cFaces[cf];

                if
                (
                    facePatch[faceI] == -1
                    && mesh_.isInternalFace(faceI)
                 )
                {
                    facePatch[faceI] = getBafflePatch
                    (
                        adaptPatchIDs,
                        facePatch,
                        faceI
                     );
                    nBaffleFaces++;

                        // Mark face as a 'boundary'
                    markBoundaryFace
                    (
                        faceI,
                        isBoundaryFace,
                        isBoundaryEdge,
                        isBoundaryPoint
                     );
                }
                else if (facePatch[faceI] == -1)
                {
                    label patchI = mesh_.boundaryMesh().whichPatch(faceI);
                    if (mesh_.boundaryMesh()[patchI].coupled())
                    {
                        facePatch[faceI] = getBafflePatch
                        (
                            adaptPatchIDs,
                            facePatch,
                            faceI
                         );
                        nBaffleFaces++;

                        // Mark face as a 'boundary'
                        markBoundaryFace
                        (
                            faceI,
                            isBoundaryFace,
                            isBoundaryEdge,
                            isBoundaryPoint
                         );
                    }
                }
            }
        }
    }

    // Sync all. (note that pointdata and facedata not used anymore but sync
    // anyway)

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

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>()
    );

    const bool checkFaceCollapse = false;

    // optional baffling of faces with points all boundary ones
    const Switch baffleAllPtsBoundary =
        motionDict.lookupOrDefault<Switch>("baffleAllPointsBoundary", false);

    // Find faces with all edges on the boundary and make them baffles
    for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
    {
        if (facePatch[facei] == -1)
        {
            const labelList& fEdges = mesh_.faceEdges(facei, dynFEdges);
            label nFaceBoundaryEdges = 0;

            forAll(fEdges, fe)
            {
                if (isBoundaryEdge[fEdges[fe]])
                {
                    nFaceBoundaryEdges++;
                }
            }

            if (nFaceBoundaryEdges == fEdges.size())
            {
                if
                (
                    checkFaceCollapse
                && !isCollapsedFace
                    (
                        newPoints,
                        neiCc,
                        minArea,
                        maxNonOrtho,
                        facei
                    )
                )
                {
                    nPrevented++;
                    //Pout<< "Preventing baffling (to avoid collapse) of face "
                    //    << facei
                    //    << " with all boundary edges "
                    //    << " at " << mesh_.faceCentres()[facei]
                    //    << endl;
                }
                else
                {
                    facePatch[facei] = getBafflePatch
                                       (
                                           adaptPatchIDs,
                                           facePatch,
                                           facei
                                       );

                    nBaffleFaces++;

                    // Do NOT update boundary data since this would grow blocked
                    // faces across gaps.
                }
            }
            else if (baffleAllPtsBoundary)
            {
                const face f = mesh_.faces()[facei];
                label nFaceBoundaryPts = 0;
                forAll(f, fp)
                {
                    if (isBoundaryPoint[f[fp]])
                    {
                        nFaceBoundaryPts++;
                    }
                }

                if (nFaceBoundaryPts == f.size())
                {
                    facePatch[facei] = getBafflePatch
                    (
                        adaptPatchIDs,
                        facePatch,
                        facei
                    );
                    nBaffleFaces++;
                }
            }
        }
    }

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            label facei = pp.start();

            forAll(pp, i)
            {
                if (facePatch[facei] == -1)
                {
                    const labelList& fEdges = mesh_.faceEdges(facei, dynFEdges);
                    label nFaceBoundaryEdges = 0;

                    forAll(fEdges, fe)
                    {
                        if (isBoundaryEdge[fEdges[fe]])
                        {
                            nFaceBoundaryEdges++;
                        }
                    }

                    if (nFaceBoundaryEdges == fEdges.size())
                    {
                        if
                        (
                            checkFaceCollapse
                        && !isCollapsedFace
                            (
                                newPoints,
                                neiCc,
                                minArea,
                                maxNonOrtho,
                                facei
                            )
                        )
                        {
                            nPrevented++;
                            //Pout<< "Preventing baffling of coupled face "
                            //    << facei
                            //    << " with all boundary edges "
                            //    << " at " << mesh_.faceCentres()[facei]
                            //    << endl;
                        }
                        else
                        {
                            facePatch[facei] = getBafflePatch
                            (
                                adaptPatchIDs,
                                facePatch,
                                facei
                            );

                            if (isMasterFace[facei])
                            {
                                nBaffleFaces++;
                            }

                            // Do NOT update boundary data since this would grow
                            // blocked faces across gaps.
                        }
                    }
                }

                facei++;
            }
        }
    }


    // Because of isCollapsedFace one side can decide not to baffle whereas
    // the other side does so sync. Baffling is prefered over not baffling.
    if (checkCollapse)  // Or always?
    {
        syncTools::syncFaceList
        (
            mesh_,
            facePatch,
            maxEqOp<label>()
        );
    }

    Info<< "markFacesOnProblemCells : marked "
        << returnReduce(nBaffleFaces, sumOp<label>())
        << " additional internal faces to be converted into baffles."
        << endl;

    if (checkCollapse)
    {
        Info<< "markFacesOnProblemCells : prevented "
            << returnReduce(nPrevented, sumOp<label>())
            << " internal faces from getting converted into baffles."
            << endl;
    }

    //Check for internal faceZones
    bool foundZone = false;
    forAll(mesh_.faceZones(), fzonei)
    {
        const faceZone& fZone = mesh_.faceZones()[fzonei];
        forAll(fZone, fzi)
        {
            label facei = fZone[fzi];
            label patchi = mesh_.boundaryMesh().whichPatch(facei);
            if (patchi == -1 || patches[patchi].coupled())
            {
                foundZone = true;
                break;
            }
        }
        if (foundZone)
        {
            break;
        }
    }
    if (returnReduce(foundZone, orOp<bool>()))
    {
        foundZone = false;
        forAll(mesh_.cells(), celli)
        {
            const cell& c = mesh_.cells()[celli];
            DynamicList<label> zonedFaces(c.size());
            label nBoundFaces = 0;
            forAll(c,ci)
            {
                label facei = c[ci];
                label zoneID = mesh_.faceZones().whichZone(facei);
                label patchi = patches.whichPatch(facei);
                if (patchi != -1 && !mesh_.boundaryMesh()[patchi].coupled())
                {
                    nBoundFaces++;
                }
                else if (zoneID != -1)
                {
                    zonedFaces.append(facei);
                }
            }

            if (zonedFaces.size() > 0 && nBoundFaces > 0)
            {
                forAll(zonedFaces, j)
                {
                    label facei = zonedFaces[j];

                    if (!removeFromChecks[facei])
                    {
                        foundZone = true;
                        markBoundaryFace
                        (
                            facei,
                            isBoundaryFace,
                            isBoundaryEdge,
                            isBoundaryPoint
                         );
                    }
                }
            }
        }

        if (returnReduce(foundZone, orOp<bool>()))
        {
            syncTools::syncPointList
            (
                mesh_,
                isBoundaryPoint,
                orEqOp<bool>(),
                false               // null value
             );

            forAll(cellLevel, celli)
            {
                const labelList& cPoints = mesh_.cellPoints(celli, dynCPoints);
                label nBoundaryAnchors = 0;
                forAll(cPoints, i)
                {
                    label pointi = cPoints[i];

                    if (pointLevel[pointi] <= cellLevel[celli])
                    {
                        // Anchor point
                        if (isBoundaryPoint[pointi])
                        {
                            nBoundaryAnchors++;
                        }
                    }
                }

                if (nBoundaryAnchors == 8)
                {
                    const cell& cFaces = mesh_.cells()[celli];
                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];
                        if (facePatch[facei] == -1)
                        {
                            facePatch[facei] = getBafflePatch
                            (
                                adaptPatchIDs,
                                facePatch,
                                facei
                             );
                        }
                    }
                }
            }
            syncTools::syncFaceList
            (
                mesh_,
                facePatch,
                maxEqOp<label>()
            );
        }
    }

    // set all remaining faces not set with nearest patch value

    label testI = 0;
    forAll(facePatch, facei)
    {
        if (facePatch[facei] == -2)
        {
            testI++;
        }
    }

    pointField fc(testI);

    // Repeat (most of) snappySnapDriver::doSnap
    if (returnReduce(fc.size(), sumOp<label>()))
    {
        labelList faceMap(testI);
        testI = 0;
        forAll(facePatch, facei)
        {
            if (facePatch[facei] == -2)
            {
                fc[testI] = mesh_.faceCentres()[facei];
                faceMap[testI] = facei;
                testI++;
            }
        }

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces_.findNearest
        (
            checkSurfaces,
            fc,
            scalarField(testI, sqr(GREAT)),    // sqr of attraction
            hitSurface,
            hitInfo
         );


        labelList pRegions(hitInfo.size());
        forAll(checkSurfaces, sI)
        {
            DynamicList<pointIndexHit> localHits;
            label surfI = checkSurfaces[sI];

            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    localHits.append(hitInfo[i]);
                }
            }
            labelList localRegions;
            label geomI = surfaces_.surfaces()[surfI];
            surfaces_.geometry()[geomI].getRegion(localHits, localRegions);

            label localI = 0;
            forAll(hitSurface, i)
            {
                if (hitSurface[i] == surfI)
                {
                    pRegions[i] = localRegions[localI];
                    localI++;
                }
            }
        }
        forAll(hitInfo, i)
        {
            label meshFaceI = faceMap[i];

            if (hitInfo[i].hit())
            {
                label global =
                    surfaces_.globalRegion(hitSurface[i], pRegions[i]);
                facePatch[meshFaceI] = globalToPatch[global];
            }
            else
            {
                //assign to first available patch
                facePatch[meshFaceI] = 0;
            }
        }
    }

    return facePatch;
}


// Mark faces to be baffled to prevent snapping problems. Does
// test to find nearest surface and checks which faces would get squashed.
Foam::labelList Foam::meshRefinement::markFacesOnProblemCellsGeometric
(
    const snapParameters& snapParams,
    const dictionary& motionDict,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch
) const
{
    pointField oldPoints(mesh_.points());

    // Repeat (most of) autoSnapDriver::doSnap
    {
        labelList adaptPatchIDs(meshedPatches());

        // Construct addressing engine.
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh_,
                adaptPatchIDs
            )
        );
        indirectPrimitivePatch& pp = ppPtr();

        // Distance to attract to nearest feature on surface
        const scalarField snapDist
        (
            snappySnapDriver::calcSnapDistance(*this, snapParams, pp)
        );


        // Construct iterative mesh mover.
        Info<< "Constructing mesh displacer ..." << endl;
        Info<< "Using mesh parameters " << motionDict << nl << endl;

        const pointMesh& pMesh = pointMesh::New(mesh_);

        motionSmoother meshMover
        (
            mesh_,
            pp,
            adaptPatchIDs,
            meshRefinement::makeDisplacementField(pMesh, adaptPatchIDs)(),
            motionDict
        );


        // Check initial mesh
        Info<< "Checking initial mesh ..." << endl;
        labelHashSet wrongFaces(mesh_.nFaces()/100);
        label nInitErrors = motionSmoother::checkMesh(false, mesh_, motionDict, wrongFaces);

        Info<< "Detected " << nInitErrors << " illegal faces"
            << " (concave, zero area or negative cell pyramid volume)"
            << endl;


        Info<< "Checked initial mesh in = "
            << mesh_.time().cpuTimeIncrement() << " s\n" << nl << endl;

        // Pre-smooth patch vertices (so before determining nearest)
        snappySnapDriver::preSmoothPatch
        (
            *this,
            snapParams,
            nInitErrors,
            List<labelPair>(0),
            meshMover
        );

        pointField nearestPoint;
        vectorField nearestNormal;
        const vectorField disp
        (
            snappySnapDriver::calcNearestSurface
            (
                snapParams.strictRegionSnap(),
                *this,
                globalToMasterPatch,
                globalToSlavePatch,
                snapDist,   // attraction
                pp,
                nearestPoint,
                nearestNormal
            )
        );

        const labelList& meshPoints = pp.meshPoints();

        pointField newPoints(mesh_.points());
        forAll(meshPoints, i)
        {
            newPoints[meshPoints[i]] += disp[i];
        }

        syncTools::syncPointList
        (
            mesh_,
            newPoints,
            minMagSqrEqOp<point>(),     // combine op
            vector(GREAT, GREAT, GREAT) // null value (note: cannot use VGREAT)
        );

        mesh_.movePoints(newPoints);
    }


    // Per face the nearest adaptPatch
    const labelList nearestAdaptPatch(nearestPatch(meshedPatches()));

    // Per face (internal or coupled!) the patch that the
    // baffle should get (or -1).
    labelList facePatch(mesh_.nFaces(), -1);
    // Count of baffled faces
    label nBaffleFaces = 0;

    {
        faceSet wrongFaces(mesh_, "wrongFaces", 100);
        {
            //motionSmoother::checkMesh(false, mesh_, motionDict, wrongFaces);

            // Just check the errors from squashing
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            const labelList allFaces(identity(mesh_.nFaces()));
            label nWrongFaces = 0;

            //const scalar minV(readScalar(motionDict.lookup("minVol", true)));
            //if (minV > -GREAT)
            //{
            //    polyMeshGeometry::checkFacePyramids
            //    (
            //        false,
            //        minV,
            //        mesh_,
            //        mesh_.cellCentres(),
            //        mesh_.points(),
            //        allFaces,
            //        List<labelPair>(0),
            //        &wrongFaces
            //    );
            //
            //    label nNewWrongFaces = returnReduce
            //    (
            //        wrongFaces.size(),
            //        sumOp<label>()
            //    );
            //
            //    Info<< "    faces with pyramid volume < "
            //        << setw(5) << minV
            //        << " m^3                  : "
            //        << nNewWrongFaces-nWrongFaces << endl;
            //
            //    nWrongFaces = nNewWrongFaces;
            //}

            scalar minArea(readScalar(motionDict.lookup("minArea")));
            if (minArea > -SMALL)
            {
                polyMeshGeometry::checkFaceArea
                (
                    false,
                    minArea,
                    mesh_,
                    mesh_.faceAreas(),
                    allFaces,
                    &wrongFaces
                );

                label nNewWrongFaces = returnReduce
                (
                    wrongFaces.size(),
                    sumOp<label>()
                );

                Info<< "    faces with area < "
                    << setw(5) << minArea
                    << " m^2                            : "
                    << nNewWrongFaces-nWrongFaces << endl;

                nWrongFaces = nNewWrongFaces;
            }

            scalar minDet(readScalar(motionDict.lookup("minDeterminant")));
            if (minDet > -1)
            {
                polyMeshGeometry::checkCellDeterminant
                (
                    false,
                    minDet,
                    mesh_,
                    mesh_.faceAreas(),
                    allFaces,
                    polyMeshGeometry::affectedCells(mesh_, allFaces),
                    &wrongFaces
                );

                label nNewWrongFaces = returnReduce
                (
                    wrongFaces.size(),
                    sumOp<label>()
                );

                Info<< "    faces on cells with determinant < "
                    << setw(5) << minDet << "                : "
                    << nNewWrongFaces-nWrongFaces << endl;

                nWrongFaces = nNewWrongFaces;
            }
        }


        forAllConstIter(faceSet, wrongFaces, iter)
        {
            label patchi = mesh_.boundaryMesh().whichPatch(iter.key());

            if (patchi == -1 || mesh_.boundaryMesh()[patchi].coupled())
            {
                facePatch[iter.key()] = nearestAdaptPatch[iter.key()];
                nBaffleFaces++;

                //Pout<< "    " << iter.key()
                //    //<< " on patch " << mesh_.boundaryMesh()[patchi].name()
                //    << " is destined for patch " << facePatch[iter.key()]
                //    << endl;
            }
        }
    }


    // Restore points.
    mesh_.movePoints(oldPoints);


    Info<< "markFacesOnProblemCellsGeometric : marked "
        << returnReduce(nBaffleFaces, sumOp<label>())
        << " additional internal and coupled faces"
        << " to be converted into baffles." << endl;

    syncTools::syncFaceList
    (
        mesh_,
        facePatch,
        maxEqOp<label>()
    );

    return facePatch;
}


void Foam::meshRefinement::rezoneProblemFaces
(
    const dictionary& motionDict,
    const labelList& surfaceToCellZone,
    labelList& namedSurfaceIndex,
    labelList& cellToZone
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    // Mark all points and edges on baffle patches (so not on any inlets,
    // outlets etc.)
    boolList isBoundaryPoint(mesh_.nPoints(), false);
    boolList isBoundaryEdge(mesh_.nEdges(), false);
    boolList isBoundaryFace(mesh_.nFaces(), false);

    // Find nearest (non-baffle) surface
    pointField newPoints;
    boolList marked(mesh_.nCells(), false);

    const PtrList<surfaceZonesInfo>& surfZones = surfaces_.surfZones();
    labelList namedSurfaces(surfaceZonesInfo::getNamedSurfaces(surfZones));

    forAll(namedSurfaceIndex, faceI)
    {
        if (namedSurfaceIndex[faceI] >= 0)
        {
            markBoundaryFace
            (
                faceI,
                isBoundaryFace,
                isBoundaryEdge,
                isBoundaryPoint
             );
        }
    }


    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false              // null value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false              // null value
    );

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>()
    );


    // See if checking for collapse
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Collapse checking parameters
    scalar volFraction = 0.2;

    const bool checkCollapse = (volFraction >= 0);
    scalar minArea = -1;
    scalar maxNonOrtho = -1;

    if (checkCollapse)
    {
        minArea = readScalar(motionDict.lookup("minArea"));
        maxNonOrtho = readScalar(motionDict.lookup("maxNonOrtho"));

        Info<< "markFacesOnProblemCells :"
            << " Deleting all-anchor surface cells only if"
            << "snapping them violates mesh quality constraints:" << nl
            << "    snapped/original cell volume < " << volFraction << nl
            << "    face area                    < " << minArea << nl
            << "    non-orthogonality            > " << maxNonOrtho << nl
            << endl;


        pointField localPoints(mesh_.nPoints());
        labelList meshPoints(mesh_.nPoints());
        label nBoundaryPoints = 0;
        forAll(mesh_.points(), pointI)
        {
            if (isBoundaryPoint[pointI])
            {
                localPoints[nBoundaryPoints] = mesh_.points()[pointI];
                meshPoints[nBoundaryPoints] = pointI;
                nBoundaryPoints++;
            }
        }
        localPoints.setSize(nBoundaryPoints);
        meshPoints.setSize(nBoundaryPoints);

        List<pointIndexHit> hitInfo;
        labelList hitSurface;
        surfaces_.findNearest
        (
            namedSurfaces,
            localPoints,
            scalarField(localPoints.size(), sqr(GREAT)),    // sqr of attraction
            hitSurface,
            hitInfo
        );

        // Start of from current points
        newPoints = mesh_.points();

        forAll(hitInfo, i)
        {
            if (hitInfo[i].hit())
            {
                newPoints[meshPoints[i]] = hitInfo[i].hitPoint();
            }
        }
    }

    DynamicList<label> recheckCells(mesh_.nCells()/10000);

    {
        // For each cell count the number of anchor points that are on
        // the boundary:
        // 8 : check the number of (baffle) boundary faces. If 3 or more block
        //     off the cell since the cell would get squeezed down to a diamond
        //     (probably; if the 3 or more faces are unrefined (only use the
        //      anchor points))
        // 7 : store. Used to check later on whether there are points with
        //     3 or more of these cells. (note that on a flat surface a boundary
        //     point will only have 4 cells connected to it)

        // Does cell have exactly 7 of its 8 anchor points on the boundary?
        PackedBoolList hasSevenBoundaryAnchorPoints(mesh_.nCells());
        // If so what is the remaining non-boundary anchor point?
        labelHashSet nonBoundaryAnchors(mesh_.nCells()/10000);

        // On-the-fly addressing storage.
        DynamicList<label> dynFEdges;
        DynamicList<label> dynCPoints;

        forAll(cellLevel, cellI)
        {
            const labelList& cPoints = mesh_.cellPoints(cellI, dynCPoints);

            // Get number of anchor points (pointLevel <= cellLevel)

            label nBoundaryAnchors = 0;
//            label nNonAnchorBoundary = 0;
            label nonBoundaryAnchor = -1;

            forAll(cPoints, i)
            {
                label pointI = cPoints[i];

                if (pointLevel[pointI] <= cellLevel[cellI])
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
//                    nNonAnchorBoundary++;
                }
            }

            if (nBoundaryAnchors == 8)
            {
                /*
                const cell& cFaces = mesh_.cells()[cellI];

                // Count boundary faces.
                label nBfaces = 0;

                forAll(cFaces, cFaceI)
                {
                    if (isBoundaryFace[cFaces[cFaceI]])
                    {
                        nBfaces++;
                    }
                }
                */

                // If nBfaces > 1 make all non-boundary non-baffle faces baffles.
                // We assume that this situation is where there is a single
                // cell sticking out which would get flattened.

                // Eugene: delete cell no matter what.
                //if (nBfaces > 1)
                {
                    if
                    (
                        checkCollapse
                        && !isCollapsedCell(newPoints, volFraction, cellI)
                     )
                    {
                        recheckCells.append(cellI);
                    }
                    else
                    {
                        // Mark cells
                        marked[cellI] = true;
                    }
                }
            }
            else if (nBoundaryAnchors == 7)
            {
                // Mark the cell. Store the (single!) non-boundary anchor point.
                hasSevenBoundaryAnchorPoints.set(cellI, 1u);
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
            label(0)            // null value
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
                    label cellI = pCells[i];

                    if (hasSevenBoundaryAnchorPoints.get(cellI) == 1u)
                    {
                        if
                        (
                            false //disable check as want to baffle
                            && !isCollapsedCell(newPoints, volFraction, cellI)
                         )
                        {
                            recheckCells.append(cellI);
                        }
                        else
                        {
                            marked[cellI] = true;
                        }
                    }
                }
            }
        }
    }

    const labelList& faceOwn = mesh_.faceOwner();
    boolList toBaffle(mesh_.nCells(), false);
    boolList visited(mesh_.nCells(), false);
    DynamicList<label> okCells(mesh_.nCells());
    DynamicList<label> resetCells(mesh_.nCells());

    //loop setting new zone faces
    while (true)
    {
        //reset processor changed cells
        forAll(resetCells, i)
        {
            label cellI = resetCells[i];
            if (marked[cellI] && !visited[cellI])
            {
                DynamicList<label> walkedCells(mesh_.nCells()/10000);
                bool ownerSide = true;
                walkedCells.append(cellI);
                visited[cellI] = true;
                markConnected
                (
                    cellI,
                    marked,
                    isBoundaryFace,
                    toBaffle,
                    visited,
                    walkedCells,
                    ownerSide
                );
            }
        }

        //reset processor changed cells
        forAll(okCells, i)
        {
            label cellI = okCells[i];
            if (marked[cellI] && !visited[cellI])
            {
                DynamicList<label> walkedCells(mesh_.nCells()/10000);
                bool ownerSide = true;
                walkedCells.append(cellI);
                visited[cellI] = true;
                markConnected
                (
                    cellI,
                    marked,
                    isBoundaryFace,
                    toBaffle,
                    visited,
                    walkedCells,
                    ownerSide
                );

                walkedCells.shrink();

                if (ownerSide)
                {
                    forAll(walkedCells, i)
                    {
                        toBaffle[walkedCells[i]] = true;
                    }
                }
            }
        }

        //reset all unmarked cells
        forAll(mesh_.cells(), cellI)
        {
            if (marked[cellI] && !visited[cellI])
            {
                DynamicList<label> walkedCells(mesh_.nCells()/10000);
                bool ownerSide = true;
                walkedCells.append(cellI);
                visited[cellI] = true;
                markConnected
                (
                    cellI,
                    marked,
                    isBoundaryFace,
                    toBaffle,
                    visited,
                    walkedCells,
                    ownerSide
                );

                walkedCells.shrink();

                if (ownerSide)
                {
                    forAll(walkedCells, i)
                    {
                        toBaffle[walkedCells[i]] = true;
                    }
                }
            }
        }
        // Calculate coupled blocked cells
        List<bool> neiBlocked(mesh_.nFaces()-mesh_.nInternalFaces(), false);
        for
        (
            label faceI = mesh_.nInternalFaces();
            faceI < mesh_.nFaces();
            faceI++
        )
        {
            neiBlocked[faceI-mesh_.nInternalFaces()] = toBaffle[faceOwn[faceI]];
        }
        syncTools::swapBoundaryFaceList(mesh_, neiBlocked);
        resetCells.clear();
        okCells.clear();

        label nChanged = 0;
        forAll(mesh_.boundaryMesh(), patchI)
        {
            const polyPatch& patch = mesh_.boundaryMesh()[patchI];
            if (isA<processorPolyPatch>(patch))
            {
                const processorPolyPatch& pp =
                    static_cast<const processorPolyPatch&>(patch);
                if (pp.neighbProcNo() > pp.myProcNo())
                {
                    forAll(pp, i)
                    {
                        label faceI = pp.start() + i;
                        if (isBoundaryFace[faceI])
                        {
                            label own = faceOwn[faceI];
                            if
                            (
                                marked[own]
                                && neiBlocked[faceI-mesh_.nInternalFaces()]
                            )
                            {
                                if (toBaffle[own])
                                {
                                    nChanged++;
                                }
                                resetCells.append(own);
                            }
                        }
                        else
                        {
                            label own = faceOwn[faceI];
                            if
                            (
                                marked[own]
                                && neiBlocked[faceI-mesh_.nInternalFaces()]
                            )
                            {
                                if (!toBaffle[own])
                                {
                                    nChanged++;
                                }
                                okCells.append(own);
                            }
                        }
                    }
                }
            }
        }
        resetCells.shrink();
        okCells.shrink();

        if (nChanged)
        {
            visited =  false;
            toBaffle = false;
        }
        if (!returnReduce(nChanged, sumOp<label>()))
        {
            break;
        }
    }

    label nBaffleCells = 0;
    forAll(mesh_.cells(), cellI)
    {
        if (toBaffle[cellI])
        {
            nBaffleCells++;
            const cell& cellFaces = mesh_.cells()[cellI];
            label faceZoneI = -1;
            label cellZoneI = -1;
            forAll(cellFaces, i)
            {
                label faceI = cellFaces[i];
                if (namedSurfaceIndex[faceI] >= 0)
                {
                    faceZoneI = namedSurfaceIndex[faceI];
                    cellZoneI = surfaceToCellZone[namedSurfaceIndex[faceI]];
                    break;
                }
            }
            forAll(cellFaces, i)
            {
                label faceI = cellFaces[i];
                if (namedSurfaceIndex[faceI] >= 0)
                {
                    namedSurfaceIndex[faceI] = -1;
                }
                else
                {
                    namedSurfaceIndex[faceI] = faceZoneI;
                }
            }
            if (cellToZone[cellI] >= 0)
            {
                cellToZone[cellI] = -1;
            }
            else
            {
                cellToZone[cellI]  = cellZoneI;
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        namedSurfaceIndex,
        maxEqOp<label>()
     );

    Info<< "rezoneProblemFaces : marked "
        << returnReduce(nBaffleCells, sumOp<label>())
        << " to inverse baffle."
        << endl;

    return;
}


void Foam::meshRefinement::markConnected
(
    const label cellI,
    const boolList& marked,
    const boolList& isBoundaryFace,
    const boolList& toBaffle,
    boolList& visited,
    DynamicList<label>& walkedCells,
    bool& ownerSide
) const
{
    const cell& cellFaces = mesh_.cells()[cellI];

    forAll(cellFaces, i)
    {
        label faceI = cellFaces[i];

        label own = mesh_.faceOwner()[faceI];

        if (mesh_.isInternalFace(faceI))
        {
            label nei = mesh_.faceNeighbour()[faceI];

            label neiCellI = (own == cellI ? nei : own);
            if (!isBoundaryFace[faceI])
            {
                if (!visited[neiCellI] && marked[neiCellI])
                {
                    visited[neiCellI] = true;
                    walkedCells.append(neiCellI);
                    markConnected
                    (
                        neiCellI,
                        marked,
                        isBoundaryFace,
                        toBaffle,
                        visited,
                        walkedCells,
                        ownerSide
                    );
                }
            }
            else if (toBaffle[neiCellI] && marked[neiCellI])
            {
                ownerSide = false;
            }
        }
    }
}


label Foam::meshRefinement::removeProblemCells()
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    boolList isBoundaryPt(mesh_.nPoints(),false);
    boolList isBoundaryFace(mesh_.nFaces(),false);
    forAll(mesh_.faces(), faceI)
    {
        label patchI = patches.whichPatch(faceI);
        if (patchI != -1 && !patches[patchI].coupled())
        {
            const face& f = mesh_.faces()[faceI];
            isBoundaryFace[faceI] = true;
            forAll(f, fp)
            {
                isBoundaryPt[f[fp]] = true;
            }
        }
    }

    boolList markedForCheck(mesh_.nCells(),false);
    boolList markedForRemoval(mesh_.nCells(),false);
    labelList nMarkedCells(mesh_.nPoints(), 0);
    labelList nBoundaryFaces(mesh_.nPoints(), 0);
    while (true)
    {
        nMarkedCells = 0;
        nBoundaryFaces = 0;
        markedForCheck = false;

        label nRemoved = 0;
        forAll(markedForRemoval, cellI)
        {
            if (markedForRemoval[cellI])
            {
                const labelList& c = mesh_.cells()[cellI];
                forAll(c, cFI)
                {
                    label faceI = c[cFI];
                    const face& f = mesh_.faces()[faceI];
                    isBoundaryFace[faceI] = true;
                    forAll(f, fp)
                    {
                        isBoundaryPt[f[fp]] = true;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            isBoundaryPt,
            orEqOp<bool>(),
            false               // null value
        );

        syncTools::syncFaceList
        (
            mesh_,
            isBoundaryFace,
            orEqOp<bool>()
        );

        forAll(mesh_.cells(), cellI)
        {
            const labelList& c = mesh_.cells()[cellI];
            const labelList& cPts = mesh_.cellPoints()[cellI];
            label clevel = cellLevel[cellI];

            scalar totalBoundaryArea = 0;
            scalar totalCellArea = 0;
            label nBoundFaces = 0;
            forAll(c, cFI)
            {
                label faceI = c[cFI];
                scalar fA = mesh_.magFaceAreas()[faceI];
                totalCellArea += fA;

                if (isBoundaryFace[faceI])
                {
                    nBoundFaces++;
                    totalBoundaryArea += fA;
                }
            }

            label nBoundAnchors = 0;
            forAll(cPts, cPtI)
            {
                label pointI = cPts[cPtI];
                if (isBoundaryPt[pointI] && pointLevel[pointI] <= clevel)
                {
                    nBoundAnchors++;
                }
            }
            scalar bRatio = (totalBoundaryArea/totalCellArea);
            if (nBoundAnchors == 8 && bRatio > 0.65)
            {
                if (!markedForRemoval[cellI])
                {
                    markedForRemoval[cellI] = true;
                    nRemoved++;
                }
            }
            else if
            (
                (nBoundAnchors > 6 && bRatio > 0.49)
                || (c.size() == 6 && nBoundAnchors == 6 && nBoundFaces == 1)
            )
            {
                markedForCheck[cellI] = true;
            }
        }

        forAll(mesh_.points(), pointI)
        {
            const labelList& pCells = mesh_.pointCells()[pointI];
            forAll(pCells, pCI)
            {
                if (markedForCheck[pCells[pCI]])
                {
                    nMarkedCells[pointI]++;
                }
            }
            const labelList& pFaces = mesh_.pointFaces()[pointI];
            forAll(pFaces, pFI)
            {
                if (isBoundaryFace[pFaces[pFI]])
                {
                    nBoundaryFaces[pointI]++;
                }
            }

        }

        syncTools::syncPointList
        (
            mesh_,
            nMarkedCells,
            plusEqOp<label>(),
            label(0)               // null value
         );

        forAll(mesh_.points(), pointI)
        {
            if (nMarkedCells[pointI] > 3 && nBoundaryFaces[pointI] == 0)
            {
                const labelList& pCells = mesh_.pointCells()[pointI];
                forAll(pCells, pCI)
                {
                    label cellI = pCells[pCI];
                    if (markedForCheck[cellI] && !markedForRemoval[cellI])
                    {
                        markedForRemoval[cellI] = true;
                        nRemoved++;
                    }
                }
            }
        }
        if (returnReduce(nRemoved, sumOp<label>()) == 0)
        {
            break;
        }

    }

    //See if any gap cells present and prevent removal
    if (mesh_.foundObject<volScalarField>("gapCells"))
    {
        const volScalarField& gapCells =
            mesh_.lookupObject<volScalarField>("gapCells");
        forAll(gapCells, cellI)
        {
            if (gapCells[cellI] > SMALL)
            {
                markedForRemoval[cellI] = false;
            }
        }
    }

    DynamicList<label> cellsToRemove(mesh_.nCells()/100);
    forAll(mesh_.cells(), cellI)
    {
        if (markedForRemoval[cellI])
        {
            cellsToRemove.append(cellI);
        }
    }

    cellsToRemove.shrink();

    Info<<"Removing "<< returnReduce(cellsToRemove.size(), sumOp<label>())
        <<" cells with seven boundary anchors" <<endl;

    if (returnReduce(cellsToRemove.size(), sumOp<label>()) != 0)
    {
        labelList ownPatch(mesh_.nFaces(), -1);
        forAll(cellsToRemove, i)
        {
            label cellI = cellsToRemove[i];
            cell c = mesh_.cells()[cellI];
            label patchI = 0;
            forAll(c, cFI)
            {
                label faceI = c[cFI];
                label fp = patches.whichPatch(faceI);
                if (fp != -1 && !patches[fp].coupled())
                {
                    patchI = fp;
                    break;
                }
            }

            forAll(c, cFI)
            {
                label faceI = c[cFI];
                ownPatch[faceI] = patchI;
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
            label faceI = exposedFaces[i];
            exposedPatchIDs[i] = ownPatch[faceI];
        }

        doRemoveCells
        (
            cellsToRemove,
            exposedFaces,
            exposedPatchIDs,
            cellRemover
         );
    }

    return cellsToRemove.size();
}


label Foam::meshRefinement::removeSelectedCellSets
(
    const dictionary& dict
)
{
    PackedList<1> reportCells(mesh_.nCells(), 0);

    if (dict.found("cellRemovalSets"))
    {
        PtrList<entry> regions(dict.lookup("cellRemovalSets"));

        forAll(regions, regioni)
        {
            const entry& region = regions[regioni];

            autoPtr<topoSetSource> cellSelector =
                topoSetSource::New(region.keyword(), mesh_, region.dict());

            cellSet selectedCellSet
            (
                mesh_,
                "cellSet",
                mesh_.nCells()/10+1  // Reasonable size estimate.
            );

            cellSelector->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
            );

            Switch invert = region.dict().lookupOrDefault("invert", false);
            if (invert)
            {
                selectedCellSet.invert(mesh_.nCells());
            }

            forAllConstIter(labelHashSet, selectedCellSet, iter)
            {
                reportCells.set(iter.key(), 1);
            }
        }
    }

    DynamicList<label> cellsToRemove(mesh_.nCells()/100);
    forAll(mesh_.cells(), cellI)
    {
        if (reportCells.get(cellI) == 1)
        {
            cellsToRemove.append(cellI);
        }
    }
    cellsToRemove.shrink();

    label nCellstoRemove = cellsToRemove.size();
    reduce(nCellstoRemove, sumOp<label>());

    if (nCellstoRemove == 0)
    {
        return 0;
    }
    else
    {
        Info<< "Removing " << nCellstoRemove << " cells"
             << " using removal sets interface" << endl;

        word exposedPatchName =
            dict.lookupOrDefault<word>("patchName", "exposedPatch");

        dictionary patchInfo;
        patchInfo.set("type", wallPolyPatch::typeName);

        label patchI = addMeshedPatch
        (
            exposedPatchName,
            patchInfo
        );

        // Remove cells
        removeCells cellRemover(mesh_);

        // Pick up patches for exposed faces
        labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
        labelList exposedPatchIDs(exposedFaces.size(),patchI);

        doRemoveCells
        (
            cellsToRemove,
            exposedFaces,
            exposedPatchIDs,
            cellRemover
        );

        return  nCellstoRemove;
    }
}


void Foam::meshRefinement::removeExtrusionProblemCells
(
    const dictionary& meshDict,
    const refinementParameters& refineParams
)
{
    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");
    const dictionary& layerDict = meshDict.subDict("addLayersControls");

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& owners = mesh_.faceOwner();
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    //Need to ensure grown up patches don't have cells with all boundary pts
    layerParameters layerParams(layerDict,mesh_.boundaryMesh(),true);
    const labelList& numLayers = layerParams.numLayers();
    boolList grownUpPatches(patches.size(),false);
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && numLayers[patchI] == -1)
        {
            grownUpPatches[patchI] = true;
        }
    }

    bool removeBoundaryAnyGrown = refineParams.removeBoundaryAnyGrownUp();
    if (!removeBoundaryAnyGrown)
    {
        const dictionary& snapDict = meshDict.subDict("snapControls");
        removeBoundaryAnyGrown =
           snapDict.lookupOrDefault<bool>("preFaceMergeExtrude", false);
    }

    label holeSize = refineParams.holeSize();

    bool firstPass = true;
    bool secondPass = false;

    label nSubIter = 0;
    while (true)
    {
        nSubIter++;
        label nBaffleHoles = baffleHoles(holeSize);
        dupNonManifoldPoints();

        labelList boundaryPts(mesh_.nPoints(), -1);
        labelList boundaryFaces(mesh_.nFaces(), -1);
        labelList boundaryEdges(mesh_.nEdges(), -1);

        label nBdyFaces = 0;
        forAll(mesh_.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            if (patchi != -1 && !patches[patchi].coupled())
            {
                nBdyFaces++;
                boundaryFaces[facei] = patchi;
                const face& f = mesh_.faces()[facei];
                forAll(f, fp)
                {
                    boundaryPts[f[fp]] = patchi;
                }

                const labelList& fEdges = mesh_.faceEdges()[facei];
                forAll(fEdges, fe)
                {
                    boundaryEdges[fEdges[fe]] = patchi;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            boundaryPts,
            maxEqOp<label>(),
            label(-1)               // null value
        );

        syncTools::syncEdgeList
        (
            mesh_,
            boundaryEdges,
            maxEqOp<label>(),
            label(-1)               // null value
        );


        labelList addressing(nBdyFaces);
        nBdyFaces = 0;
        forAll(mesh_.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            if (patchi != -1 && !patches[patchi].coupled())
            {
                addressing[nBdyFaces++] = facei;
            }
        }
        autoPtr<indirectPrimitivePatch> ppPtr
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>(mesh_.faces(), addressing),
                mesh_.points()
             )
        );
        const indirectPrimitivePatch& pp = ppPtr();
        const labelList meshEdges
        (
            pp.meshEdges(mesh_.edges(), mesh_.pointEdges())
        );

        boolList excludedFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh_,
            mesh_.points(),
            pp,
            meshEdges,
            excludedFaces,
            0.707,//convex
            0.707,//concave
            -0.9848 //baffle
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        labelList nPointCells(mesh_.nPoints(), 0);
        labelList nPointFaces(mesh_.nPoints(), 0);
        labelList nEdgeFaces(mesh_.nEdges(), 0);

        PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));
        PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh_));

        forAll(mesh_.points(), pointi)
        {
            if (boundaryPts[pointi] != -1)
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    nPointCells[pointi]++;
                }

                const labelList& pFaces = mesh_.pointFaces()[pointi];
                forAll(pFaces, pFI)
                {
                    label facei = pFaces[pFI];
                    if
                    (
                        isMasterFace.get(facei) == 1
                        && boundaryFaces[facei] != -1
                    )
                    {
                        nPointFaces[pointi]++;
                    }
                }
            }
        }

        labelList edgeType(mesh_.nEdges(), -1);
        forAll(meshEdges, edgei)
        {
            label meshedgei = meshEdges[edgei];
            if
            (
                eType[edgei].first() == edgeClassification::CONCAVE
            )
            {
                edgeType[meshedgei] = 1;
            }
            else if
            (
                eType[edgei].first() == edgeClassification::CONVEX
            )
            {
                edgeType[meshedgei] = 3;
            }
            else if
            (
                eType[edgei].first() == edgeClassification::BAFFLE
            )
            {
                edgeType[meshedgei] = 4;
            }
            else
            {
                edgeType[meshedgei] = 2;
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            edgeType,
            maxEqOp<label>(),
            label(0)               // null value
        );

        forAll(mesh_.edges(), edgei)
        {
            if (boundaryEdges[edgei] != -1)
            {
                const labelList& eFaces = mesh_.edgeFaces()[edgei];
                forAll(eFaces, eFI)
                {
                    label facei = eFaces[eFI];
                    if
                    (
                        isMasterFace.get(facei) == 1
                        && boundaryFaces[facei] != -1
                    )
                    {
                        nEdgeFaces[edgei]++;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh_,
            nPointFaces,
            plusEqOp<label>(),
            label(0)               // null value
        );

        syncTools::syncPointList
        (
            mesh_,
            nPointCells,
            plusEqOp<label>(),
            label(0)               // null value
        );

        syncTools::syncEdgeList
        (
            mesh_,
            nEdgeFaces,
            plusEqOp<label>(),
            label(0)               // null value
        );

        labelList nConvexEdges(mesh_.nPoints(), 0);
        forAll(mesh_.points(), pointi)
        {
            const labelList& pEdges = mesh_.pointEdges()[pointi];
            forAll(pEdges, pei)
            {
                label edgei = pEdges[pei];
                if
                (
                    isMasterEdge.get(edgei) == 1
                    && boundaryEdges[edgei] != -1
                    && edgeType[edgei] == 3
                )
                {
                    nConvexEdges[pointi]++;
                }
            }
        }
        syncTools::syncPointList
        (
            mesh_,
            nConvexEdges,
            plusEqOp<label>(),
            label(0)               // null value
        );

        labelList allBoundaryCells(mesh_.nCells(), -1);
        boolList markedForRemoval(mesh_.nCells(), false);

        forAll(mesh_.cells(), celli)
        {
            const labelList& cPts = mesh_.cellPoints()[celli];
            label nBoundaryPts = 0;

            forAll(cPts, cptI)
            {
                if (boundaryPts[cPts[cptI]] != -1)
                {
                    nBoundaryPts++;
                }
            }
            if (nBoundaryPts == cPts.size())
            {
                const cell& c = mesh_.cells()[celli];
                label clevel = cellLevel[celli];
                label nAnchorBoundFaces = 0;

                bool grownUpPatch = false;
                bool allGrownUpPatches = true;

                forAll(c,cFI)
                {
                    label facei = c[cFI];
                    label patchi = patches.whichPatch(facei);
                    if (patchi != -1 && !patches[patchi].coupled())
                    {
                        const face& f = mesh_.faces()[facei];
                        bool refineFace = false;
                        forAll(f,fp)
                        {
                            if (pointLevel[f[fp]] > clevel)
                            {
                                refineFace = true;
                                break;
                            }
                        }
                        if (refineFace)
                        {
                            nAnchorBoundFaces++;
                        }
                        else
                        {
                            nAnchorBoundFaces += 4;
                        }

                        if (grownUpPatches[patchi])
                        {
                            grownUpPatch = true;
                        }
                        else
                        {
                            allGrownUpPatches = false;
                        }
                    }
                }

                if
                (
                    (allGrownUpPatches && grownUpPatch)
                    || (removeBoundaryAnyGrown && grownUpPatch)
                )
                {
                    markedForRemoval[celli] = true;
                }

                if (nAnchorBoundFaces >= 16)
                {
                    allBoundaryCells[celli] = label(2);
                }
                else if (nAnchorBoundFaces >= 12)
                {
                    allBoundaryCells[celli] = label(1);
                }
                else
                {
                    allBoundaryCells[celli] = label(0);
                }
            }
        }

        labelList neiAllBoundaryCells(mesh_.nFaces()-mesh_.nInternalFaces());
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            neiAllBoundaryCells[facei-mesh_.nInternalFaces()] =
                allBoundaryCells[owners[facei]];
        }
        syncTools::swapBoundaryFaceList(mesh_, neiAllBoundaryCells);

        boolList singleBoundaryNei(mesh_.nCells(), false);
        //Remove single connected all boundary cells
        forAll(mesh_.cells(), celli)
        {
            if (allBoundaryCells[celli] < 1)
            {
                continue;
            }

            const cell& c = mesh_.cells()[celli];
            label nBoundNei = 0;
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchi = patches.whichPatch(facei);
                if (patchi == -1)
                {
                    label nei = owners[facei];
                    if (nei == celli)
                    {
                        nei = mesh_.faceNeighbour()[facei];
                    }
                    if (allBoundaryCells[nei] > -1)
                    {
                        nBoundNei++;
                    }
                }
                else if (patches[patchi].coupled())
                {
                    if (neiAllBoundaryCells[facei-mesh_.nInternalFaces()] > -1)
                    {
                        nBoundNei++;
                    }
                }
            }

            const labelList& cPts = mesh_.cellPoints()[celli];
            bool singleGap = false;
            if (cPts.size() == 8)
            {
                const labelList& cEdges = mesh_.cellEdges()[celli];
                label nConcave = 0;
                label nConvex = 0;
                label nInternal = 0;
                forAll(cEdges, cEI)
                {
                    label edgei = cEdges[cEI];
                    if (edgeType[edgei] == 3)
                    {
                        nConvex++;
                    }
                    else if (edgeType[edgei] == 1)
                    {
                        nConcave++;
                    }
                    else if (edgeType[edgei] == -1)
                    {
                        nInternal++;
                    }
                }

                if (nConcave == 5 && nConvex == 6 && nInternal == 1)
                {
                    singleGap = true;
                }
            }

            if (nBoundNei == 1)
            {
                singleBoundaryNei[celli] = true;
            }

            if ((allBoundaryCells[celli] == 2 && nBoundNei == 0) || singleGap)
            {
                markedForRemoval[celli] = true;
            }
        }

        boolList neiSingleBoundaryNei(mesh_.nFaces()-mesh_.nInternalFaces());
        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            neiSingleBoundaryNei[facei-mesh_.nInternalFaces()] =
                singleBoundaryNei[owners[facei]];
        }
        syncTools::swapBoundaryFaceList(mesh_, neiSingleBoundaryNei);

        forAll(mesh_.cells(), celli)
        {
            if (!singleBoundaryNei[celli])
            {
                continue;
            }

            const cell& c = mesh_.cells()[celli];
            label nBoundNei = 0;
            forAll(c, cFI)
            {
                label facei = c[cFI];
                label patchi = patches.whichPatch(facei);
                if (patchi == -1)
                {
                    label nei = owners[facei];
                    if (nei == celli)
                    {
                        nei = mesh_.faceNeighbour()[facei];
                    }
                    if (singleBoundaryNei[nei])
                    {
                        nBoundNei++;
                    }
                }
                else if (patches[patchi].coupled())
                {
                    if (neiSingleBoundaryNei[facei-mesh_.nInternalFaces()])
                    {
                        nBoundNei++;
                    }
                }
            }

            if (nBoundNei == 1)
            {
                markedForRemoval[celli] = true;
            }
        }

        //Remove corner cells
        forAll(mesh_.cells(), celli)
        {
            const labelList& cPts = mesh_.cellPoints()[celli];
            const labelList& cEdges = mesh_.cellEdges()[celli];
            if (cPts.size() == 8 && cEdges.size() == 12)
            {
                label nCornerEdges = 0;
                label nConvexEdges = 0;
                label nInternal = 0;
                bool simpleCorner = true;
                forAll(cEdges, cEI)
                {
                    label edgei = cEdges[cEI];
                    if (edgeType[edgei] == 1)
                    {
                        edge e = mesh_.edges()[edgei];

                        nCornerEdges++;
                        if (nPointFaces[e[0]] > 5 || nPointFaces[e[1]] > 5)
                        {
                            simpleCorner = false;
                            break;
                        }
                    }
                    else if (edgeType[edgei] == 3)
                    {
                        nConvexEdges++;
                    }
                    else if (edgeType[edgei] == -1)
                    {
                        nInternal++;
                    }
                    else
                    {
                        break;
                    }
                }
                if
                (
                    simpleCorner && nCornerEdges == 3
                    && nConvexEdges == 6 && nInternal == 3
                )
                {
                    markedForRemoval[celli] = true;
                }
            }
        }

        forAll(mesh_.points(), pointi)
        {
            if
            (
                ( nPointCells[pointi] == 2 && nPointFaces[pointi] == 6 )
                ||
                (
                    nPointCells[pointi] == 6 && nPointFaces[pointi] == 6
                    && nConvexEdges[pointi] == 6
                )
            )
            {
                const labelList& pCells = mesh_.pointCells()[pointi];
                forAll(pCells, pCI)
                {
                    markedForRemoval[pCells[pCI]] = true;
                }
            }
        }

        //find processor non-manifold edges
        {
            labelList maxProcID(mesh_.nPoints(), Pstream::myProcNo());
            labelList minProcID(mesh_.nPoints(), Pstream::myProcNo());

            syncTools::syncPointList
            (
                mesh_,
                maxProcID,
                maxEqOp<label>(),
                label(0)               // null value
             );

            syncTools::syncPointList
            (
                mesh_,
                minProcID,
                minEqOp<label>(),
                label(0)               // null value
             );

            List<pointField> pECC(mesh_.nPoints(), pointField(0));
            forAll(mesh_.points(), pointi)
            {
                const labelList& pEdges = mesh_.pointEdges()[pointi];
                pointField checkPts(pEdges.size());
                label nPts = 0;
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    if (edgeType[edgei] == 1)
                    {
                        const edge& e = mesh_.edges()[edgei];
                        point eC = e.centre(mesh_.points());
                        checkPts[nPts++] = eC;
                    }
                }
                checkPts.setSize(nPts);
                pECC[pointi] = checkPts;
            }

            syncTools::syncPointList
            (
                mesh_,
                pECC,
                pointFieldCombine(),
                pointField() // initial value
            );

            forAll(mesh_.edges(), edgei)
            {
                if (edgeType[edgei] == 1)
                {
                    const edge& e = mesh_.edges()[edgei];
                    point eC = e.centre(mesh_.points());
                    label processorManifold = true;
                    forAll(e, ei)
                    {
                        label pointi = e[ei];
                        const pointField& pECentres = pECC[pointi];
                        label nEdgesFound = 0;
                        forAll(pECentres, peci)
                        {
                            if (mag(eC-pECentres[peci]) < SMALL)
                            {
                                nEdgesFound++;
                            }
                        }
                        if
                        (
                            nEdgesFound != 2
                            || (minProcID[pointi] == maxProcID[pointi])
                        )
                        {
                            processorManifold = false;
                            break;
                        }
                    }
                    if (processorManifold)
                    {
                        const labelList& eCells = mesh_.edgeCells()[edgei];
                        forAll(eCells, eCI)
                        {
                            markedForRemoval[eCells[eCI]] = true;
                        }
                    }
                }
            }
        }

        //Baffle based on edge topology
        forAll(mesh_.edges(), edgei)
        {
            bool baffleCells = false;
            edge e = mesh_.edges()[edgei];

            if (nEdgeFaces[edgei] == 4)
            {
                baffleCells = true;
            }
            else if
            (
                edgeType[edgei] == 1
                && (nPointCells[e[0]] == 6 || nPointCells[e[1]] == 6)
            )
            {
                baffleCells = true;
            }
            else if
            (
                edgeType[edgei] == 1
                &&
                (
                    (
                        nPointCells[e[0]] == 1 && nPointCells[e[1]] == 5
                        && nPointFaces[e[1]] > 5
                     )
                    ||
                    (
                        nPointCells[e[1]] == 1 && nPointCells[e[0]] == 5
                        && nPointFaces[e[0]] > 5
                     )
                 )
            )
            {
                baffleCells = true;
            }
            else if
            (
                nPointCells[e[0]] == 6 && nPointCells[e[1]] == 6
                && nPointFaces[e[0]] == 6 && nPointFaces[e[1]] == 6
            )
            {
                baffleCells = true;
            }
            else if
            (
                (nPointCells[e[0]] ==  6 && nPointCells[e[1]] == 5
                && nPointFaces[e[0]] == 6 && nPointFaces[e[1]] == 7)
                || (nPointCells[e[1]] ==  6 && nPointCells[e[0]] == 5
                && nPointFaces[e[1]] == 6 && nPointFaces[e[0]] == 7)
            )
            {
                baffleCells = true;
            }

            if
            (
                baffleCells
            )
            {
                const labelList& eCells = mesh_.edgeCells()[edgei];
                forAll(eCells, eCI)
                {
                    label celli = eCells[eCI];
                    markedForRemoval[celli] = true;
                }
            }
        }

        if ((firstPass || secondPass) && nSubIter < 5)
        {
            forAll(mesh_.faces(), facei)
            {
                if (boundaryFaces[facei] != -1)
                {
                    const labelList& fEdges = mesh_.faceEdges()[facei];
                    label nBaffleEdges = 0;
                    label nCornerEdges = 0;
                    label nConvexEdges = 0;

                    forAll(fEdges, fEI)
                    {
                        label edgei = fEdges[fEI];
                        if (edgeType[edgei] == 1)
                        {
                            nCornerEdges++;
                        }
                        else if (edgeType[edgei] == 3)
                        {
                            nConvexEdges++;
                        }
                        else if (edgeType[edgei] == 4)
                        {
                            nBaffleEdges++;
                        }
                    }

                    if (firstPass)
                    {
                        if
                        (
                            (nCornerEdges == 2 && nBaffleEdges == 2)
                            || (nCornerEdges > 2 && nBaffleEdges == 1)
                        )
                        {
                            markedForRemoval[mesh_.faceOwner()[facei]] = true;
                        }
                        else if
                        (
                            nCornerEdges == 1 && nConvexEdges == 0
                            && nBaffleEdges > 0
                        )
                        {
                            bool keepCell = false;
                            if (nBaffleEdges == 1)
                            {
                                const face& f = mesh_.faces()[facei];

                                forAll(f,fp)
                                {
                                    label currPt = f[fp];
                                    label nextIndex = f.fcIndex(fp);
                                    label nextPt = f[nextIndex];
                                    label meshEdgeI = meshTools::findEdge
                                    (
                                        mesh_.edges(),
                                        mesh_.pointEdges()[currPt],
                                        currPt,
                                        nextPt
                                    );

                                    if (edgeType[meshEdgeI] == 1)
                                    {
                                        label prevPt = f[f.rcIndex(fp)];
                                        label prevEdgeI = meshTools::findEdge
                                        (
                                            mesh_.edges(),
                                            mesh_.pointEdges()[prevPt],
                                            prevPt,
                                            currPt
                                        );
                                        label extendedPt =
                                            f[f.fcIndex(nextIndex)];
                                        label nextEdgeI = meshTools::findEdge
                                        (
                                            mesh_.edges(),
                                            mesh_.pointEdges()[nextPt],
                                            nextPt,
                                            extendedPt
                                        );
                                        if (prevEdgeI != -1 && nextEdgeI != -1)
                                        {
                                            if
                                            (
                                                edgeType[prevEdgeI] == 4
                                                || edgeType[nextEdgeI] == 4
                                            )
                                            {
                                                keepCell = true;
                                            }
                                        }
                                        break;
                                    }
                                }
                            }
                            if (!keepCell)
                            {
                                label own = mesh_.faceOwner()[facei];
                                markedForRemoval[own] = true;
                            }
                        }
                    }
                    else
                    {
                        if ((nCornerEdges+nBaffleEdges) == fEdges.size())
                        {
                            markedForRemoval[mesh_.faceOwner()[facei]] = true;
                        }
                    }
                }
            }
        }
        else
        {
            forAll(mesh_.faces(), facei)
            {
                if (boundaryFaces[facei] != -1)
                {
                    const labelList& fEdges = mesh_.faceEdges()[facei];
                    forAll(fEdges, fEI)
                    {
                        if (edgeType[fEdges[fEI]] == 4)
                        {
                            markedForRemoval[mesh_.faceOwner()[facei]] = true;
                            break;
                        }
                    }
                }
            }
        }

        label nMarkedForRemoval = 0;
        forAll(mesh_.cells(), celli)
        {
            if (markedForRemoval[celli])
            {
                nMarkedForRemoval++;
            }
        }

        Info<<"Removing "<< returnReduce(nMarkedForRemoval, sumOp<label>())
            <<" cells with problem extrusions" <<endl;

        if (nBaffleHoles + returnReduce(nMarkedForRemoval, sumOp<label>()) != 0)
        {
            boolList blockedFace(mesh_.nFaces(),false);
            forAll(mesh_.cells(), celli)
            {
                if (markedForRemoval[celli])
                {
                    const cell& c = mesh_.cells()[celli];
                    forAll(c, cfi)
                    {
                        blockedFace[c[cfi]] = true;
                    }
                }
            }

            syncTools::syncFaceList(mesh_, blockedFace, orEqOp<bool>());

            // Set region per cell based on walking
            regionSplit cellRegion(mesh_, blockedFace);

            labelHashSet keepRegionSet(keepLargestRegions(cellRegion));
            forAll(mesh_.cells(), celli)
            {
                if (!keepRegionSet.found(cellRegion[celli]))
                {
                    markedForRemoval[celli] = true;
                }
            }

            DynamicList<label> cellsToRemove(mesh_.nCells()/100);
            forAll(mesh_.cells(), celli)
            {
                if (markedForRemoval[celli])
                {
                    cellsToRemove.append(celli);
                }
            }
            cellsToRemove.shrink();

            labelList ownPatch(mesh_.nFaces(), -1);
            forAll(cellsToRemove, i)
            {
                label celli = cellsToRemove[i];
                const labelList& cPts = mesh_.cellPoints()[celli];
                label patchI = 0;
                forAll(cPts, cPI)
                {
                    label pointI = cPts[cPI];
                    if (boundaryPts[pointI] != -1)
                    {
                        patchI = boundaryPts[pointI];
                        break;
                    }
                }

                const cell& c = mesh_.cells()[celli];
                forAll(c, cFI)
                {
                    label faceI = c[cFI];
                    ownPatch[faceI] = patchI;
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
                label faceI = exposedFaces[i];
                exposedPatchIDs[i] = ownPatch[faceI];
            }

            doRemoveCells
            (
                cellsToRemove,
                exposedFaces,
                exposedPatchIDs,
                cellRemover
            );

            if (!firstPass && !secondPass)
            {
                if (refineParams.mergePreExtrude())
                {
                    mergePatchFacesUndo
                    (
                        Foam::cos(degToRad(45.0)),
                        Foam::cos(degToRad(45.0)),
                        -1, //max area ratio (-1 to disable)
                        meshedPatches(),
                        motionDict,
                        false, //update intersections
                        false //maintain patches
                     );
                    mergeEdgesUndo(-1., motionDict, false);
                }
            }
        }
        else if (firstPass)
        {
            if (refineParams.mergePreExtrude())
            {
                mergePatchFacesUndo
                (
                    Foam::cos(degToRad(45.0)),
                    Foam::cos(degToRad(45.0)),
                    -1, //max area ratio (-1 to disable)
                    meshedPatches(),
                    motionDict,
                    false, //update intersections
                    false //maintain patches
                 );

                mergeEdgesUndo(-1., motionDict, false);
            }
            firstPass = false;
            secondPass = true;
            nSubIter = 0;
        }
        else if (secondPass)
        {
            if (refineParams.mergePreExtrude())
            {
                mergePatchFacesUndo
                (
                    Foam::cos(degToRad(45.0)),
                    Foam::cos(degToRad(45.0)),
                    -1, //max area ratio (-1 to disable)
                    meshedPatches(),
                    motionDict,
                    false, //update intersections
                    false //maintain patches
                 );

                mergeEdgesUndo(-1., motionDict, false);
            }
            secondPass = false;
            nSubIter = 0;
        }
        else
        {
            break;
        }
    }

    return;
}


// ************************************************************************* //
