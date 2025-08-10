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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "boundaryLayerRefinement/boundaryLayerRefinement.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChanger/polyTopoChanger.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::boundaryLayerRefinement::boundaryLayerRefinement
(
    const indirectPrimitivePatch& pp,
    meshRefinement& meshRefiner
)
:
    pp_(pp),
    meshRefiner_(meshRefiner)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boundaryLayerRefinement::addTopFaces
(
    const fvMesh& mesh,
    const labelList& pointLevel,
    const labelList& newEdgePts,
    const label& cLevel,
    const face& f,
    const label& anchorPt,
    const label& newMidPt,
    DynamicList<label>& faceVerts
)
{
    faceVerts.append(anchorPt);
    label anchorFp = findIndex(f,anchorPt);
    label startFp = anchorFp;
    while (true)
    {
        label nextFp = f.fcIndex(startFp);
        label startPt = f[startFp];
        label nextPt = f[nextFp];
        label startLev = pointLevel[startPt];
        label nextLev = pointLevel[nextPt];
        label edgei =  meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[startPt],
            startPt,
            nextPt
        );
        label newEdgePt = newEdgePts[edgei];
        if (newEdgePt != -1)
        {
            faceVerts.append(newEdgePt);
            if (startLev > cLevel || nextLev > cLevel)
            {
                faceVerts.append(nextPt);
            }
            break;
        }

        if (nextLev <= cLevel)
        {
            break;
        }
        else if (nextLev <= cLevel+1)
        {
            faceVerts.append(nextPt);
            break;
        }
        else
        {
            faceVerts.append(nextPt);
        }
        startFp = nextFp;
    }
    faceVerts.append(newMidPt);

    startFp = anchorFp;
    DynamicList<label> backWalk(f.size());
    while (true)
    {
        label prevFp = f.rcIndex(startFp);
        label startPt = f[startFp];
        label prevPt = f[prevFp];
        label startLev = pointLevel[startPt];
        label prevLev = pointLevel[prevPt];
        label edgei =  meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[startPt],
            startPt,
            prevPt
        );

        label newEdgePt = newEdgePts[edgei];
        if (newEdgePt != -1)
        {
            if (startLev > cLevel || prevLev > cLevel)
            {
                faceVerts.append(prevPt);
            }
            backWalk.append(newEdgePt);
            break;
        }

        if (pointLevel[prevPt] <= cLevel)
        {
            break;
        }
        else if (pointLevel[prevPt] <= cLevel+1)
        {
            backWalk.append(prevPt);
            break;
        }
        else
        {
            backWalk.append(prevPt);
        }
        startFp = prevFp;
    }
    label backWalkSz = backWalk.size();
    for (label bwi = backWalkSz-1; bwi >= 0; bwi--)
    {
        faceVerts.append(backWalk[bwi]);
    }

    return;
}


void Foam::boundaryLayerRefinement::addAnchorToMid
(
    const fvMesh& mesh,
    const labelList& pointLevel,
    const labelList& newEdgePts,
    const label& cLevel,
    const face& f,
    const label& anchorPt,
    DynamicList<label>& faceVerts
)
{
    label pointi = f[anchorPt];
    faceVerts.append(pointi);
    label startFp = anchorPt;
    while (true)
    {
        label nextFp = f.fcIndex(startFp);
        label nextMeshPointI = f[nextFp];
        label edgei =  meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[f[startFp]],
            f[startFp],
            f[nextFp]
        );
        label newEdgePt = newEdgePts[edgei];

        if (newEdgePt != -1)
        {
            faceVerts.append(newEdgePt);
            break;
        }

        if (pointLevel[nextMeshPointI] <= cLevel)
        {
            break;
        }
        else if (pointLevel[nextMeshPointI] <= cLevel+1)
        {
            faceVerts.append(nextMeshPointI);
            break;
        }
        else
        {
            faceVerts.append(nextMeshPointI);
        }
        startFp = nextFp;
    }

    return;
}


label Foam::boundaryLayerRefinement::findMidPoint
(
    const fvMesh& mesh,
    const labelList& pointLevel,
    const labelList& newEdgePts,
    const label& cLevel,
    const face& f,
    const label& anchorStart,
    const label& anchorEnd
)
{
    label startFp = anchorStart;
    label midPt = -1;
    label maxLevel = labelMin;
    while (true)
    {
        label nextFp = f.fcIndex(startFp);
        label nextMeshPointI = f[nextFp];
        label startPt = f[startFp];
        label nextPt = f[nextFp];
        label edgei =  meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[startPt],
            startPt,
            nextPt
        );
        label newEdgePt = newEdgePts[edgei];
        label startLev = pointLevel[startPt];
        label nextLev = pointLevel[nextPt];

        if
        (
            newEdgePt != -1
            && startLev <= cLevel && nextLev <= cLevel
        )
        {
            midPt = newEdgePt;
            break;
        }
        label nextLevel = pointLevel[nextMeshPointI];

        if (nextLevel > maxLevel)
        {
            midPt = nextMeshPointI;
            maxLevel = nextLevel;
        }

        if (nextFp == anchorEnd)
        {
            break;
        }

        startFp = nextFp;
    }

    return midPt;
}


void Foam::boundaryLayerRefinement::setRefinement()
{
    Info<<"Refining boundary layer cells"<<endl;

    fvMesh& mesh = meshRefiner_.mesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    const labelList meshEdges(pp_.meshEdges(mesh.edges(), mesh.pointEdges()));

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbors = mesh.faceNeighbour();

    volScalarField& layerCells = const_cast<volScalarField&>
        (mesh.lookupObject<volScalarField>("layerStacks"));

    labelList neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        label own = owners[facei];
        neiLayerCells[facei-mesh.nInternalFaces()] =
            layerCells[own];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLayerCells);

    // New point/cell level. Copy of pointLevel for existing points.
    DynamicList<label> newCellLevel(cellLevel.size());
    forAll(cellLevel, celli)
    {
        newCellLevel.append(cellLevel[celli]);
    }
    DynamicList<label> newPointLevel(pointLevel.size());
    forAll(pointLevel, pointi)
    {
        newPointLevel.append(pointLevel[pointi]);
    }

    polyTopoChange meshMod(mesh);

    //New bottom and top centre point
    List<label> newTopBottomCentrePoint
    (
        mesh.nFaces(),
        -1
    );
    List<label> markedTopLevel
    (
        mesh.nFaces(),
        -1
    );

    //PP faces anchor cells
    List<Map<label>> newAnchorCells
    (
        pp_.size(),
        Map<label>(16)
    );

    //top face map
    List<label> topface
    (
        pp_.size(),
        -1
    );

    //Map from upper to lower anchors
    List<label> upperToLowerAnchorMap(mesh.nPoints(), -1);

    //Add edge mid points and any new faces
    labelList newEdgePts(mesh.nEdges(),-1);

    forAll(pp_, i)
    {
        label facei = pp_.addressing()[i];
        label celli = owners[facei];

        if (layerCells[celli] < 0)
        {
            continue;
        }

        const cell& c = mesh.cells()[celli];
        label cLevel = cellLevel[celli];

        //Add bottom split edge points
        labelList fEdges = mesh.faceEdges()[facei];
        forAll(fEdges, fei)
        {
            label edgei = fEdges[fei];
            if (newEdgePts[edgei] != -1)
            {
                continue;
            }
            const edge e = mesh.edges()[edgei];
            label mshPt0 = e[0];
            label mshPt1 = e[1];
            label pLevel0 = pointLevel[mshPt0];
            label pLevel1 = pointLevel[mshPt1];
            if (pLevel0 <= cLevel && pLevel1 <= cLevel)
            {
                point midPoint = 0.5*
                (
                    mesh.points()[mshPt0]
                    + mesh.points()[mshPt1]
                );

                label newpointi = meshMod.addPoint
                (
                    midPoint,         // point
                    -1,         // master point
                    -1,      // zone for point
                    true        // supports a cell
                );
                newEdgePts[edgei] = newpointi;
                newPointLevel(newpointi) = cLevel+1;
            }
        }


        //Identify top face
        forAll(c, cFI)
        {
            label cfacei = c[cFI];
            if (cfacei == facei)
            {
                continue;
            }

            label patchi = patches.whichPatch(cfacei);
            label neiLayerID = -1;

            if (patchi == -1 || patches[patchi].coupled())
            {
                if (patchi == -1)
                {
                    label nei = neighbors[cfacei];
                    neiLayerID  = (owners[cfacei] == celli)
                        ? layerCells[nei] : layerCells[owners[cfacei]];
                }
                else if (patches[patchi].coupled())
                {
                    neiLayerID =
                        neiLayerCells[cfacei-mesh.nInternalFaces()];
                }
            }
            else
            {
                continue;
            }
            if (neiLayerID == -1)
            {
                topface[i] = cfacei;
                break;
            }
        }
        label topfacei = topface[i];

        //Add top split edge points
        fEdges = mesh.faceEdges()[topfacei];
        forAll(fEdges, fei)
        {
            label edgei = fEdges[fei];
            if (newEdgePts[edgei] != -1)
            {
                continue;
            }
            const edge e = mesh.edges()[edgei];
            label mshPt0 = e[0];
            label mshPt1 = e[1];
            label pLevel0 = pointLevel[mshPt0];
            label pLevel1 = pointLevel[mshPt1];
            if (pLevel0 <= cLevel && pLevel1 <= cLevel)
            {
                point midPoint = 0.5*
                (
                    mesh.points()[mshPt0]
                    + mesh.points()[mshPt1]
                );

                label newpointi = meshMod.addPoint
                (
                    midPoint,         // point
                    -1,         // master point
                    -1,      // zone for point
                    true        // supports a cell
                );
                newEdgePts[edgei] = newpointi;
                newPointLevel(newpointi) = cLevel+1;
            }
        }
    }

    boolList updatedFaces(mesh.nFaces(), false);
    forAll(pp_, i)
    {
        label facei = pp_.addressing()[i];
        label celli = owners[facei];
        if (layerCells[celli] < 0)
        {
            continue;
        }

        const cell& c = mesh.cells()[celli];
        forAll(c,cFI)
        {
            updatedFaces[c[cFI]] = true;
        }

        face f = mesh.faces()[facei];

        DynamicList<label> anchorPts(f.size());
        label cLevel = cellLevel[celli];
        forAll(f,fp)
        {
            label pointi = f[fp];
            label pLevel = pointLevel[pointi];
            if (pLevel <= cLevel)
            {
                anchorPts.append(pointi);
            }
        }

        face af(anchorPts);
        if (af.size() < 3)
        {
            continue;
        }

        label topfacei = topface[i];
        if (topfacei == -1)
        {
            continue;
        }

        //Create map between bottom and top points
        labelHashSet topFaceSet(mesh.faces()[topfacei]);
        forAll(f,fp)
        {
            label pointi = f[fp];
            label pLevel = pointLevel[pointi];
            if (pLevel <= cLevel)
            {
                const labelList pEdges = mesh.pointEdges()[pointi];
                forAll(pEdges, pEI)
                {
                    label edgei = pEdges[pEI];
                    const edge& e = mesh.edges()[edgei];
                    label otherPt(pointi == e[0] ? e[1] : e[0]);
                    if (topFaceSet.found(otherPt))
                    {
                        upperToLowerAnchorMap[otherPt] = pointi;
                        break;
                    }
                }
            }
        }

        //Add boundary central point and remove face
        const point& fc =  mesh.faceCentres()[facei];
        label newpointi = meshMod.addPoint
        (
            fc,         // point
            -1,         // master point
            -1,      // zone for point
            true        // supports a cell
        );
        meshMod.removeFace(facei, -1);
        newTopBottomCentrePoint[facei] = newpointi;
        newPointLevel(newpointi) = cLevel+1;

        //Add top face central point and remove face
        const point& tfc =  mesh.faceCentres()[topfacei];
        newpointi = meshMod.addPoint
        (
            tfc,         // point
            -1,         // master point
            -1,      // zone for point
            true        // supports a cell
        );
        newTopBottomCentrePoint[topfacei] = newpointi;
        markedTopLevel[topfacei] = cLevel;
        meshMod.removeFace(topfacei, -1);
        newPointLevel(newpointi) = cLevel+1;

        //Add anchor cells

        label ownZoneI = mesh.cellZones().whichZone(celli);
        forAll(af, afi)
        {
            label pointi = af[afi];

            label newCellI = celli;
            if (afi == 0)
            {
                // Update cell level
                newCellLevel[celli] = cLevel+1;
            }
            else
            {
                newCellI = meshMod.addCell
                (
                    -1,             // master point
                    -1,             // master edge
                    -1,             // master face
                    celli,            //master
                    ownZoneI        // zone for cell
                 );
                newCellLevel(newCellI) = cLevel+1;
            }

            newAnchorCells[i].insert
            (
                pointi,
                newCellI
            );
        }

        //Add bottom faces
        {
            face f = mesh.faces()[facei];
            forAll(f,fp)
            {
                label pointi = f[fp];
                label pLevel = pointLevel[pointi];
                if (pLevel <= cLevel)
                {
                    DynamicList<label> faceVerts(f.size());
                    addTopFaces
                    (
                        mesh,
                        pointLevel,
                        newEdgePts,
                        cLevel,
                        f,
                        pointi,
                        newTopBottomCentrePoint[facei],
                        faceVerts
                    );
                    label newCellI = newAnchorCells[i][pointi];
                    face newFace(faceVerts);
                    label patchi =
                        patches.whichPatch(facei);
                    label zonei =
                        mesh.faceZones().whichZone(facei);
                    bool zoneFlip = false;
                    if (zonei >= 0)
                    {
                        const faceZone& fZone =
                            mesh.faceZones()[zonei];
                        zoneFlip =
                            fZone.flipMap()[fZone.whichFace(facei)];
                    }
                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            newFace,  // face
                            newCellI,    // owner
                            -1,       // neighbour
                            -1,       // master point
                            -1,       // master edge
                            facei,    // master face for addition
                            false,    // flux flip
                            patchi,   // patch for face
                            zonei,   // zone for face
                            zoneFlip  // face zone flip
                         )
                     );

                }
            }
        }

        //Add top faces
        {
            face f = mesh.faces()[topfacei];
            forAll(f,fp)
            {
                label pointi = f[fp];
                label pLevel = pointLevel[pointi];
                if (pLevel <= cLevel)
                {
                    DynamicList<label> faceVerts(f.size());
                    addTopFaces
                    (
                        mesh,
                        pointLevel,
                        newEdgePts,
                        cLevel,
                        f,
                        pointi,
                        newTopBottomCentrePoint[topfacei],
                        faceVerts
                     );

                    label baseAnchor = upperToLowerAnchorMap[pointi];
                    label newCellI = newAnchorCells[i][baseAnchor];
                    face newFace(faceVerts);
                    label patchi =
                        patches.whichPatch(topfacei);
                    label zonei =
                        mesh.faceZones().whichZone(topfacei);
                    bool zoneFlip = false;
                    if (zonei >= 0)
                    {
                        const faceZone& fZone =
                            mesh.faceZones()[zonei];
                        zoneFlip =
                            fZone.flipMap()[fZone.whichFace(topfacei)];
                    }

                    label own = newCellI;
                    label nei = -1;
                    if (patchi == -1)
                    {
                        nei = newCellI;
                        if (owners[topfacei] == celli)
                        {
                            own = neighbors[topfacei];
                            zoneFlip = !zoneFlip;
                            newFace.flip();
                        }
                        else
                        {
                            own = owners[topfacei];
                        }
                    }
                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            newFace,  // face
                            own,    // owner
                            nei,       // neighbour
                            -1,       // master point
                            -1,       // master edge
                            topfacei,    // master face for addition
                            false,    // flux flip
                            patchi,   // patch for face
                            zonei,   // zone for face
                            zoneFlip  // face zone flip
                         )
                     );

                }
            }
        }

        //Add new internal faces
        {
            const face& f = mesh.faces()[facei];
            label startAnchor = -1;
            forAll(f,fp)
            {
                label pointi = f[fp];
                label pLevel = pointLevel[pointi];
                if (pLevel <= cLevel)
                {
                    startAnchor = fp;
                    break;
                }
            }

            label topfacei = topface[i];
            face tf = mesh.faces()[topfacei];
            label toppatchi =
                patches.whichPatch(topfacei);
            //flip face so ordered same as bottom face
            if (owners[topfacei] == celli || toppatchi != -1)
            {
                tf.flip();
            }

            label nextFp = startAnchor;
            label prevFp = startAnchor;
            while (true)
            {
                nextFp = f.fcIndex(nextFp);
                label pointi = f[nextFp];
                label pLevel = pointLevel[pointi];

                if (pLevel <= cLevel)
                {
                    label startPt = f[prevFp];
                    label endPt = f[nextFp];

                    label newOwn = newAnchorCells[i][startPt];
                    label newNei = newAnchorCells[i][endPt];

                    DynamicList<label> newPts(f.size());

                    newPts.append(newTopBottomCentrePoint[facei]);

                    label midPt = findMidPoint
                    (
                        mesh,
                        pointLevel,
                        newEdgePts,
                        cLevel,
                        f,
                        prevFp,
                        nextFp
                    );
                    newPts.append(midPt);

                    label startUpper = -1;
                    label endUpper = -1;
                    forAll(tf,fp)
                    {
                        label upointi = tf[fp];
                        if (upperToLowerAnchorMap[upointi] == startPt)
                        {
                            startUpper = fp;
                        }
                        else if (upperToLowerAnchorMap[upointi] == endPt)
                        {
                            endUpper = fp;
                        }
                    }

                    midPt = findMidPoint
                    (
                        mesh,
                        pointLevel,
                        newEdgePts,
                        cLevel,
                        tf,
                        startUpper,
                        endUpper
                    );
                    newPts.append(midPt);
                    newPts.append(newTopBottomCentrePoint[topfacei]);

                    face newFace(newPts);
                    if (newOwn > newNei)
                    {
                       label swpOwn = newOwn;
                       newOwn = newNei;
                       newNei = swpOwn;
                       newFace.flip();
                    }

                    //Add new face
                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            newFace,  // face
                            newOwn,    // owner
                            newNei,       // neighbour
                            -1,       // master point
                            -1,       // master edge
                            -1,    // master face for addition
                            false,    // flux flip
                            -1,   // patch for face
                            -1,   // zone for face
                            false  // face zone flip
                         )
                     );

                    if (nextFp == startAnchor)
                    {
                        break;
                    }
                    else
                    {
                        prevFp = nextFp;
                    }
                }
            }
        }
    }
    boolList splitEdges(mesh.nEdges(), false);
    forAll(splitEdges, edgei)
    {
        if (newEdgePts[edgei] != -1)
        {
            splitEdges[edgei] = true;
        }
    }
    syncTools::syncEdgeList
    (
        mesh,
        splitEdges,
        orEqOp<bool>(),
        false // initial value
    );

    //Add processor split edge points
    forAll(splitEdges, edgei)
    {
        if (splitEdges[edgei] && newEdgePts[edgei] == -1)
        {
            const edge e = mesh.edges()[edgei];

            label pLevel0 = pointLevel[e[0]];
            label pLevel1 = pointLevel[e[1]];

            point midPoint = 0.5*
            (
                mesh.points()[e[0]]
                + mesh.points()[e[1]]
            );

            label newpointi = meshMod.addPoint
            (
                midPoint,         // point
                -1,         // master point
                -1,      // zone for point
                true        // supports a cell
            );
            newEdgePts[edgei] = newpointi;
            newPointLevel(newpointi) = max(pLevel0,pLevel1)+1;
        }
    }

    //Update side faces
    boolList processorSplitEdges(mesh.nEdges(), false);
    forAll(pp_.edges(), edgei)
    {
        const labelList& eFaces = pp_.edgeFaces()[edgei];
        label i = eFaces[0];
        label facei = pp_.addressing()[i];
        label celli = owners[facei];
        if (layerCells[celli] < 0)
        {
            continue;
        }
        label meshedgei = meshEdges[edgei];
        const edge& e = mesh.edges()[meshedgei];
        label pLevel0 = pointLevel[e[0]];
        label pLevel1 = pointLevel[e[1]];
        label cLevel = cellLevel[celli];

        const labelList& meshEdgeFaces = mesh.edgeFaces()[meshedgei];
        labelHashSet cFaceSet(mesh.cells()[celli]);
        label sidefacei = -1;
        forAll(meshEdgeFaces, mEFI)
        {
            label otherFaceI = meshEdgeFaces[mEFI];
            if (otherFaceI != facei && cFaceSet.found(otherFaceI))
            {
                sidefacei = otherFaceI;
                break;
            }
        }
        label patchi =
            patches.whichPatch(sidefacei);
        label zonei =
            mesh.faceZones().whichZone(sidefacei);
        bool zoneFlip = false;
        if (zonei >= 0)
        {
            const faceZone& fZone =
                mesh.faceZones()[zonei];
            zoneFlip =
                fZone.flipMap()[fZone.whichFace(sidefacei)];
        }

        bool splitEdge = false;
        label anchorPt = -1;
        if (pLevel0 <= cLevel && pLevel1 <= cLevel)
        {
            if (patchi != -1 && patches[patchi].coupled())
            {
                processorSplitEdges[meshedgei] = true;
            }
            splitEdge = true;
        }
        else if (newEdgePts[meshedgei] != -1)
        {
            anchorPt = (pLevel0 <= cLevel ? e[0] : e[1]);
            splitEdge = true;
        }

        if (splitEdge)
        {
            //Remove old side face
            meshMod.removeFace(sidefacei, -1);

            const face& f = mesh.faces()[facei];
            label topfacei = topface[i];
            face tf = mesh.faces()[topfacei];
            label toppatchi =
                patches.whichPatch(topfacei);
            bool flipUpper = false;
            //flip face so ordered same as bottom face
            if (owners[topfacei] == celli || toppatchi != -1)
            {
                flipUpper = true;
            }

            forAll(e, ei)
            {
                label pointi = e[ei];
                label otherPt(ei == 0 ? e[1] : e[0]);
                label own = celli;
                label newOwn = -1;
                label nei = -1;
                label newNei = -1;

                if (anchorPt != -1)
                {
                    newOwn = newAnchorCells[i][anchorPt];
                    if (patchi == -1)
                    {
                        label j = eFaces[1];
                        newNei = newAnchorCells[j][pointi];
                    }
                }
                else
                {
                    newOwn = newAnchorCells[i][pointi];
                    if (patchi == -1)
                    {
                        nei =
                        (
                            owners[sidefacei] == own
                            ? neighbors[sidefacei]
                            : owners[sidefacei]
                        );

                        if (eFaces.size() == 2)
                        {
                            label j = eFaces[1];
                            label nbrLev = cellLevel[nei];
                            label nbrAnchor = pointi;
                            if (pLevel0 > nbrLev)
                            {
                               nbrAnchor = e[1];
                            }
                            else if (pLevel1 > nbrLev)
                            {
                               nbrAnchor = e[0];
                            }
                            newNei = newAnchorCells[j][nbrAnchor];
                        }
                        else
                        {
                            newNei = nei;
                        }
                    }
                }
                bool flipLower = false;
                label lowerFp = findIndex(f,pointi);
                if (f[f.fcIndex(lowerFp)] != otherPt)
                {
                    flipLower = true;
                }

                face lf(flipLower ? f.reverseFace() : f);
                face uf = tf;
                if
                (
                    (flipLower && !flipUpper)
                    || (!flipLower && flipUpper)
                )
                {
                    uf.flip();
                }
                lowerFp = findIndex(lf,pointi);

                label topFp = -1;
                label topStart = -1;
                forAll(uf,fp)
                {
                    label upointi = uf[fp];
                    if (upperToLowerAnchorMap[upointi] == pointi)
                    {
                        topStart = upointi;
                        break;
                    }
                }
                topFp = findIndex(uf,topStart);
                DynamicList<label> faceVerts(f.size());
                addAnchorToMid
                (
                    mesh,
                    pointLevel,
                    newEdgePts,
                    cLevel,
                    lf,
                    lowerFp,
                    faceVerts
                );
                DynamicList<label> topPts(f.size());
                addAnchorToMid
                (
                    mesh,
                    pointLevel,
                    newEdgePts,
                    cLevel,
                    uf,
                    topFp,
                    topPts
                );

                label backWalkSz = topPts.size();
                for (label bwi = backWalkSz-1; bwi >= 0; bwi--)
                {
                    faceVerts.append(topPts[bwi]);
                }
                bool flipFace = false;
                if (flipLower)
                {
                    //out pointing
                    if (patchi == -1)
                    {
                        if (newOwn > newNei)
                        {
                            label tmpNewOwn = newOwn;
                            newOwn = newNei;
                            newNei = tmpNewOwn;
                            flipFace = true;
                        }
                    }
                }
                else
                {
                    //in pointing
                    if (patchi == -1)
                    {
                        if (newOwn > newNei)
                        {
                            label tmpNewOwn = newOwn;
                            newOwn = newNei;
                            newNei = tmpNewOwn;
                        }
                        else
                        {
                            flipFace = true;
                        }
                    }
                    else
                    {
                        flipFace = true;
                    }
                }

                face newFace(faceVerts);
                if (flipFace)
                {
                    newFace.flip();
                }
                meshMod.setAction
                (
                    polyAddFace
                    (
                        newFace,  // face
                        newOwn,   // owner
                        newNei,   // neighbour
                        -1,       // master point
                        -1,       // master edge
                        sidefacei,    // master face for addition
                        false,    // flux flip
                        patchi,   // patch for face
                        zonei,   // zone for face
                        zoneFlip  // face zone flip
                     )
                 );
            }
        }
        else
        {
            //faces already split
            label anchorpointi(pLevel0 <= cLevel ? e[0] : e[1]);
            label own = celli;
            label newOwn = newAnchorCells[i][anchorpointi];
            label newNei = -1;
            face f = mesh.faces()[sidefacei];
            bool flipFace = false;
            if (patchi == -1)
            {
                newNei = newAnchorCells[i][anchorpointi];
                if (owners[sidefacei] == own)
                {
                    newOwn = neighbors[sidefacei];
                    if (newOwn > newNei)
                    {
                       label swpOwn = newOwn;
                       newOwn = newNei;
                       newNei = swpOwn;
                    }
                    else
                    {
                       flipFace = true;
                    }
                }
                else
                {
                    newOwn = owners[sidefacei];
                    if (newOwn > newNei)
                    {
                       label swpOwn = newOwn;
                       newOwn = newNei;
                       newNei = swpOwn;
                       flipFace = true;
                    }
                }
                if (flipFace)
                {
                    f.flip();
                }
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    f,              // modified face
                    sidefacei,      // label of face
                    newOwn,         // owner
                    newNei,            // neighbour
                    false,          // face flip
                    patchi,         // patch for face
                    false,          // remove from zone
                    zonei,         // zone for face
                    zoneFlip        // face flip in zone
                 )
             );
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        processorSplitEdges,
        orEqOp<bool>(),
        false // initial value
    );

    //Add processor split side faces
    forAll(mesh.edges(), edgei)
    {
        if (processorSplitEdges[edgei])
        {
            const labelList& eFaces = mesh.edgeFaces()[edgei];
            label sidefacei = -1;
            label sidepatchi = -1;
            label zonei =  -1;
            bool zoneFlip = false;
            forAll(eFaces, eFI)
            {
                label facei = eFaces[eFI];
                if (updatedFaces[facei])
                {
                    continue;
                }
                label patchi =
                    patches.whichPatch(facei);
                if (patchi != -1)
                {
                    if (patches[patchi].coupled())
                    {
                        sidefacei = facei;
                        sidepatchi = patchi;
                        zonei =
                            mesh.faceZones().whichZone(facei);
                        if (zonei >= 0)
                        {
                            const faceZone& fZone =
                                mesh.faceZones()[zonei];
                            zoneFlip =
                                fZone.flipMap()[fZone.whichFace(facei)];
                        }
                    }
                }
            }

            if (sidefacei != -1)
            {
                //Remove old processor face
                meshMod.removeFace(sidefacei, -1);
                updatedFaces[sidefacei] = true;

                label own = owners[sidefacei];
                label cLevel = cellLevel[own];
                label nei = -1;
                const edge e = mesh.edges()[edgei];
                forAll(e, ei)
                {
                    label pointi = e[ei];
                    label otherPt(ei == 0 ? e[1] : e[0]);
                    face f = mesh.faces()[sidefacei];
                    label anchorFp = findIndex(f,pointi);
                    label nextFp = f.fcIndex(anchorFp);
                    label upperPt = -1;
                    bool flipFace = false;
                    if (f[nextFp] != otherPt)
                    {
                        upperPt = f[nextFp];
                    }
                    else
                    {
                        flipFace = true;
                        upperPt = f[f.rcIndex(anchorFp)];
                    }

                    label upperFp = -1;
                    if (flipFace)
                    {
                        f.flip();
                    }
                    upperFp = findIndex(f,upperPt);

                    DynamicList<label> faceVerts(f.size());
                    addAnchorToMid
                    (
                        mesh,
                        pointLevel,
                        newEdgePts,
                        cLevel,
                        f,
                        upperFp,
                        faceVerts
                    );

                    f.flip();
                    label lowerFp = findIndex(f,pointi);
                    DynamicList<label> topPts(f.size());
                    addAnchorToMid
                    (
                        mesh,
                        pointLevel,
                        newEdgePts,
                        cLevel,
                        f,
                        lowerFp,
                        topPts
                    );
                    label backWalkSz = topPts.size();
                    for (label bwi = backWalkSz-1; bwi >= 0; bwi--)
                    {
                        faceVerts.append(topPts[bwi]);
                    }

                    face newFace(faceVerts);
                    if (flipFace)
                    {
                        newFace.flip();
                    }
                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            newFace,  // face
                            own,   // owner
                            nei,   // neighbour
                            -1,       // master point
                            -1,       // master edge
                            sidefacei,    // master face for addition
                            false,    // flux flip
                            sidepatchi,   // patch for face
                            zonei,   // zone for face
                            zoneFlip  // face zone flip
                         )
                    );
                }
            }
        }
    }

    //Modify other side processor top faces
    {
        syncTools::syncFaceList
        (
            mesh,
            markedTopLevel,
            maxEqOp<label>()     // combine operator
        );
        forAll(mesh.faces(), facei)
        {
            if
            (
                newTopBottomCentrePoint[facei] == -1
                && markedTopLevel[facei] != -1
            )
            {
                const point fC = mesh.faceCentres()[facei];
                label newpointi = meshMod.addPoint
                (
                    fC,         // point
                    -1,         // master point
                    -1,      // zone for point
                    true        // supports a cell
                );
                newTopBottomCentrePoint[facei] = newpointi;
                meshMod.removeFace(facei, -1);
                updatedFaces[facei] = true;
                newPointLevel(newpointi) =
                    markedTopLevel[facei]+1;

                const face& f = mesh.faces()[facei];
                label own = owners[facei];
                label nei = -1;
                label patchi =
                    patches.whichPatch(facei);
                label zonei =
                    mesh.faceZones().whichZone(facei);
                bool zoneFlip = false;
                if (zonei >= 0)
                {
                    const faceZone& fZone =
                        mesh.faceZones()[zonei];
                    zoneFlip =
                        fZone.flipMap()[fZone.whichFace(facei)];
                }

                label cLevel = markedTopLevel[facei];
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    label pLevel = pointLevel[pointi];
                    if (pLevel <= cLevel)
                    {
                        DynamicList<label> faceVerts(f.size());
                        addTopFaces
                        (
                            mesh,
                            pointLevel,
                            newEdgePts,
                            cLevel,
                            f,
                            pointi,
                            newTopBottomCentrePoint[facei],
                            faceVerts
                        );

                        face newFace(faceVerts);
                        meshMod.setAction
                        (
                            polyAddFace
                            (
                                newFace,  // face
                                own,    // owner
                                nei,       // neighbour
                                -1,       // master point
                                -1,       // master edge
                                facei,    // master face for addition
                                false,    // flux flip
                                patchi,   // patch for face
                                zonei,   // zone for face
                                zoneFlip  // face zone flip
                             )
                         );
                    }
                }
            }
        }
    }

    //Modify non split cell connected faces
    forAll(mesh.faces(), facei)
    {
        if (!updatedFaces[facei])
        {
            const labelList& fEdges = mesh.faceEdges()[facei];
            label nAddedPts = 0;
            forAll(fEdges, fEI)
            {
                label edgei = fEdges[fEI];
                if (newEdgePts[edgei] != -1)
                {
                    nAddedPts++;
                }
            }

            if (nAddedPts > 0)
            {
                const face& f = mesh.faces()[facei];
                label sz = f.size();
                DynamicList<label> newFacePts(sz);
                forAll(f,fp)
                {
                    label pointi = f[fp];
                    newFacePts.append(pointi);
                    label nextFp = f.fcIndex(fp);
                    label edgei =  meshTools::findEdge
                    (
                        mesh.edges(),
                        mesh.pointEdges()[pointi],
                        pointi,
                        f[nextFp]
                     );
                    label newEdgePt = newEdgePts[edgei];
                    if (newEdgePt > -1)
                    {
                        newFacePts.append(newEdgePt);
                    }
                }

                label patchi =
                    patches.whichPatch(facei);
                label zonei =
                    mesh.faceZones().whichZone(facei);
                bool zoneFlip = false;
                if (zonei >= 0)
                {
                    const faceZone& fZone =
                        mesh.faceZones()[zonei];
                    zoneFlip =
                        fZone.flipMap()[fZone.whichFace(facei)];
                }

                label own = owners[facei];
                label nei = -1;
                if (patchi == -1)
                {
                    nei = neighbors[facei];
                }

                face newFace(newFacePts);
                meshMod.setAction
                (
                    polyModifyFace
                    (
                        newFace,              // modified face
                        facei,      // label of face
                        own,         // owner
                        nei,            // neighbour
                        false,          // face flip
                        patchi,         // patch for face
                        false,          // remove from zone
                        zonei,         // zone for face
                        zoneFlip        // face flip in zone
                     )
                );
            }
        }
    }

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh
    (
        mesh,
        false,
        true
    );

    // Update fields (problem when face added to zero sized patch)
    mesh.updateMesh(map);

    // Move mesh if in inflation mode
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        mesh.clearOut();
    }

    meshRefiner_.updateMesh(map,labelList(0));
    meshRefiner_.meshCutter().updateLevels(newPointLevel,newCellLevel);

    //Update layer cell field
    scalar maxLayerCount(gMax(layerCells));
    label nNewCells = map().cellMap().size() - map().nOldCells();

    globalIndex globalNewCells(nNewCells);
    label start  = maxLayerCount
        + globalNewCells.offset(Pstream::myProcNo());
    label count = 1;
    forAll(map().cellMap(), cellI)
    {
        if (cellI >= map().nOldCells())
        {
            label oldCellI = map().cellMap()[cellI];
            if (oldCellI != -1 && layerCells[oldCellI] == -1)
            {
                continue;
            }

            layerCells[cellI] = start + count;
            count++;
        }
    }

    return;
}

// ************************************************************************* //
