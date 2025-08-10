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
    (c) 2014-2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "externalDisplacementMeshMover/medialAxisMeshMover.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/pointPatchFields/basic/value/valuePointPatchFields.H"
#include "algorithms/PointEdgeWave/PointEdgeWave.H"
#include "meshRefinement/meshRefinement.H"
#include "global/unitConversion/unitConversion.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "surfaceFormats/obj/OBJstream.H"
#include "algorithms/PointEdgeWave/PointData.H"
#include "externalDisplacementMeshMover/zeroFixedValue/zeroFixedValuePointPatchFields.H"
#include "fields/pointPatchFields/derived/slip/slipPointPatchField.H"
#include "fields/pointPatchFields/derived/fixedNormalSlip/fixedNormalSlipPointPatchField.H"
#include "regionSplit/regionSplit.H"
#include "edgeClassification/edgeClassification.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(medialAxisMeshMover, 0);

    addToRunTimeSelectionTable
    (
        externalDisplacementMeshMover,
        medialAxisMeshMover,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::medialAxisMeshMover::updateMeshGeometric
(
    const polyMesh& mesh,
    const pointField& newPoints,
    pointField& newFaceCentres,
    vectorField& newFaceAreas,
    pointField& newCellCentres,
    scalarField& newCellVolumes,
    vectorField& cEst,
    labelField& nCellFaces
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const faceList& fs = mesh.faces();

    //Calculate updated mesh properties
    newCellCentres = vector::zero;
    newCellVolumes = 0;
    newFaceCentres = vector::zero;
    newFaceAreas  = vector::zero;
    cEst = vector::zero;
    nCellFaces = 0;

    forAll(fs, facei)
    {
        const labelList& f = fs[facei];
        label nPoints = f.size();

        if (nPoints == 3)
        {
            newFaceCentres[facei] = (1.0/3.0)
                *(newPoints[f[0]] + newPoints[f[1]] + newPoints[f[2]]);
            newFaceAreas[facei] =
                0.5*((newPoints[f[1]] - newPoints[f[0]])
                     ^(newPoints[f[2]] - newPoints[f[0]]));
        }
        else
        {
            vector sumN = vector::zero;
            scalar sumA = 0.0;
            vector sumAc = vector::zero;

            point fCentre = newPoints[f[0]];
            for (label pi = 1; pi < nPoints; pi++)
            {
                fCentre += newPoints[f[pi]];
            }
            fCentre /= nPoints;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = newPoints[f[(pi + 1) % nPoints]];

                vector c = newPoints[f[pi]] + nextPoint + fCentre;
                vector n = (nextPoint - newPoints[f[pi]])
                    ^(fCentre - newPoints[f[pi]]);
                scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            newFaceCentres[facei] = (1.0/3.0)*sumAc/(sumA + VSMALL);
            newFaceAreas[facei] = 0.5*sumN;
        }
    }

    forAll(own, facei)
    {
        cEst[own[facei]] += newFaceCentres[facei];
        nCellFaces[own[facei]] += 1;
    }

    forAll(nei, facei)
    {
        cEst[nei[facei]] += newFaceCentres[facei];
        nCellFaces[nei[facei]] += 1;
    }

    forAll(cEst, celli)
    {
        cEst[celli] /= nCellFaces[celli];
    }

    forAll(own, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol =max
        (
            newFaceAreas[facei] & (newFaceCentres[facei] - cEst[own[facei]]),
            VSMALL
        );

        // Calculate face-pyramid centre
        vector pc =
            (3.0/4.0)*newFaceCentres[facei] + (1.0/4.0)*cEst[own[facei]];

        // Accumulate volume-weighted face-pyramid centre
        newCellCentres[own[facei]] += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        newCellVolumes[own[facei]] += pyr3Vol;
    }

    forAll(nei, facei)
    {
        // Calculate 3*face-pyramid volume
        scalar pyr3Vol = max
        (
            newFaceAreas[facei] & (cEst[nei[facei]] - newFaceCentres[facei]),
            VSMALL
        );

        // Calculate face-pyramid centre
        vector pc =
            (3.0/4.0)*newFaceCentres[facei] + (1.0/4.0)*cEst[nei[facei]];

        // Accumulate volume-weighted face-pyramid centre
        newCellCentres[nei[facei]] += pyr3Vol*pc;

        // Accumulate face-pyramid volume
        newCellVolumes[nei[facei]] += pyr3Vol;
    }

    newCellCentres /= newCellVolumes;
    newCellVolumes *= (1.0/3.0);
}


void Foam::medialAxisMeshMover::smoothDisplacement
(
    const dictionary& motionDict,
    const labelList& planarPatches,
    const scalar& edge0Len,
    const polyMesh& mesh,

    pointField newPoints,
    vectorField& displacement
)
{
    label nDispSmooth = motionDict.lookupOrDefault<label>("nVolSmoothIter", 8);

    const labelIOList& cellLevel =
        mesh.lookupObject<labelIOList>("cellLevel");

    const labelIOList& pointLevel =
        mesh.lookupObject<labelIOList>("pointLevel");

    if (nDispSmooth == 0)
    {
        return;
    }

    Info<< "snappyHexMeshDriver: Smoothing displaced mesh interior points ..."
        << endl;

    bool flattenPoints =
        motionDict.lookupOrDefault<Switch>("smoothAlignedEdges", false);

    scalar smoothTime = mesh.time().elapsedCpuTime();

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            planarPatches
        )
    );
    const indirectPrimitivePatch& pp = ppPtr();

    // Precalculate meshEdge per pp edge
    labelList meshEdges(pp.nEdges());

    forAll(meshEdges, patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];

        label v0 = pp.meshPoints()[e[0]];
        label v1 = pp.meshPoints()[e[1]];
        meshEdges[patchEdgeI] = meshTools::findEdge
        (
            mesh.edges(),
            mesh.pointEdges()[v0],
            v0,
            v1
        );
    }

    PackedList<1> isFixedPoint(mesh.nPoints(), 0);
    PackedList<1> isFixedPatchPoint(mesh.nPoints(), 0);

    labelList checkFaces(identity(mesh.nFaces()));
    // Construct table of patches to include.
    // Surface points that are to remain fixed.
    {
        boolList stationaryPoints(mesh.nPoints(), false);
        boolList stationaryPatchPoints(mesh.nPoints(), true);

        forAll(pp.meshPoints(), pointI)
        {
            const label meshPointI = pp.meshPoints()[pointI];
            stationaryPatchPoints[meshPointI] = false;
        }

        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
        for (label patchI = 0; patchI < bMesh.size(); patchI++)
        {
            if (!isA<processorPolyPatch>(bMesh[patchI]))
            {
                forAll(bMesh[patchI].meshPoints(),pointI)
                {
                    label meshPointI = bMesh[patchI].meshPoints()[pointI];
                    stationaryPoints[meshPointI] = true;
                }
            }
        }

        // mark all zone points as stationary
        forAll(mesh.faceZones(), zoneI)
        {
            const faceZone& fz = mesh.faceZones()[zoneI];
            forAll(fz, faceI)
            {
                const face& f = mesh.faces()[fz[faceI]];
                bool internalFace = false;
                if (mesh.isInternalFace(fz[faceI]))
                {
                    internalFace = true;
                }

                forAll(f, fp)
                {
                    label meshPointI = f[fp];
                    stationaryPoints[meshPointI] = true;
                    if (internalFace)
                    {
                       stationaryPatchPoints[meshPointI] = true;
                    }
                }
            }
        }

        if (mesh.foundObject<volScalarField>("protectedCells"))
        {
            volScalarField& protectedCells =
                const_cast<volScalarField&>
                (
                    mesh.lookupObject<volScalarField>("protectedCells")
                );

            forAll(mesh.cells(), cellI)
            {
                if (protectedCells[cellI] != -1)
                {
                    const labelList& cPoints = mesh.cellPoints()[cellI];
                    forAll(cPoints, ptI)
                    {
                        label meshPointI = cPoints[ptI];
                        stationaryPoints[meshPointI] = true;
                        stationaryPatchPoints[meshPointI] = true;
                    }
                }
            }
        }

        if (mesh.foundObject<volScalarField>("gapCells"))
        {
            const volScalarField& gapCells =
                mesh.lookupObject<volScalarField>("gapCells");

            forAll(gapCells, celli)
            {
                if (gapCells[celli] > -1)
                {
                    const labelList& cPoints = mesh.cellPoints()[celli];
                    forAll(cPoints, ptI)
                    {
                        label meshPointI = cPoints[ptI];
                        stationaryPoints[meshPointI] = true;
                        stationaryPatchPoints[meshPointI] = true;
                    }
                }
            }
        }

        // prevent smoothing on high aspect ratio cells
        scalar maxAR = motionDict.lookupOrDefault<scalar>("smoothMaxAR", -1);
        if (maxAR > 0)
        {
            //Mark high aspect ratio cells
            forAll(mesh.cells(), cellI)
            {
                const cell c = mesh.cells()[cellI];
                vector sA = vector::zero;
                forAll(c, cFI)
                {
                    vector fA = mesh.faceAreas()[c[cFI]];
                    sA += cmptMag(fA);
                }
                scalar maxSA = cmptMax(sA);
                scalar minSA = cmptMin(sA);

                scalar ratio = maxSA/(minSA+SMALL);

                if (ratio > maxAR)
                {
                    labelList cmeshPts = c.labels(mesh.faces());
                    forAll(cmeshPts, ptI)
                    {
                        label meshPointI = cmeshPts[ptI];
                        stationaryPoints[meshPointI] = true;
                        stationaryPatchPoints[meshPointI] = true;
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            stationaryPoints,
            orEqOp<bool>(),
            false
        );

        forAll(stationaryPoints, pointI)
        {
            if (stationaryPoints[pointI])
            {
                isFixedPoint.set(pointI, 1);
            }
        }

        labelList nExternalEdge(mesh.nEdges(), 0);
        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(pp.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            const labelList& eFaces = edgeFaces[edgeI];
            nExternalEdge[meshEdgeI] = eFaces.size();
        }

        syncTools::syncEdgeList
        (
            mesh,
            nExternalEdge,
            plusEqOp<label>(),
            label(0)               // null value
         );

        forAll(mesh.edges(), edgeI)
        {
            if (nExternalEdge[edgeI] == 1)
            {
                const edge& e = mesh.edges()[edgeI];
                stationaryPatchPoints[e[0]] = true;
                stationaryPatchPoints[e[1]] = true;
            }
        }

        // Normal component of normals of connected faces.
        vectorField edgeNormal(mesh.nEdges(), vector(GREAT, GREAT, GREAT));

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            forAll(eFaces, i)
            {
                nomalsCombine()
                (
                    edgeNormal[meshEdgeI],
                    pp.faceNormals()[eFaces[i]]
                 );
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            edgeNormal,
            nomalsCombine(),
            vector(GREAT, GREAT, GREAT)   // null value
         );

        forAll(edgeFaces, edgeI)
        {
            const labelList& eFaces = pp.edgeFaces()[edgeI];

            label meshEdgeI = meshEdges[edgeI];

            const vector& n = edgeNormal[meshEdgeI];

            if (n != vector(GREAT, GREAT, GREAT))
            {
                scalar cos = n & pp.faceNormals()[eFaces[0]];
                const edge& e = pp.edges()[edgeI];

                if (cos < 0.9396)
                {
                    stationaryPatchPoints[pp.meshPoints()[e[0]]] = true;
                    stationaryPatchPoints[pp.meshPoints()[e[1]]] = true;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            stationaryPatchPoints,
            orEqOp<bool>(),
            false
         );

        forAll(stationaryPatchPoints, pointI)
        {
            if (stationaryPatchPoints[pointI])
            {
                isFixedPatchPoint.set(pointI, 1);
            }
        }
    }

    //extend fixed points by one layer for projection method
    PackedList<1> isFixedPointExtended(mesh.nPoints(), 0);
    {
        boolList stationaryPoints(mesh.nPoints(), false);

        forAll(mesh.points(), pointI)
        {
            if
            (
                !(isFixedPoint.get(pointI) == 0
                  || isFixedPatchPoint.get(pointI) == 0)
             )
            {
                stationaryPoints[pointI] = true;
                labelList pCells = mesh.pointCells()[pointI];

                forAll(pCells, pCI)
                {
                    cell c = mesh.cells()[pCells[pCI]];
                    labelList clabels = c.labels(mesh.faces());
                    forAll(clabels, cli)
                    {
                        stationaryPoints[clabels[cli]] = true;
                    }
                }
            }
        }

        forAll(mesh.cells(), cellI)
        {
            label cLevel = cellLevel[cellI];
            const cell c = mesh.cells()[cellI];

            labelList cPoints = c.labels(mesh.faces());

            edgeList cEdges = c.edges(mesh.faces());
            labelHashSet cellMeshEdges(cEdges.size());

            forAll(cEdges, cEI)
            {
                edge e = cEdges[cEI];
                label meshEdgeI = meshTools::findEdge
                (
                    mesh.edges(),
                    mesh.pointEdges()[e[0]],
                    e[0],
                    e[1]
                );
                cellMeshEdges.insert(meshEdgeI);
            }

            forAll(cPoints, cPtI)
            {
                label meshPointI = cPoints[cPtI];

                if
                (
                    isFixedPoint.get(meshPointI) == 0
                    || isFixedPatchPoint.get(meshPointI) == 0
                 )
                {
                    if (pointLevel[meshPointI] > cLevel)
                    {
                        label nFound = 0;
                        label nFound2 = 0;
                        const labelList pEdges = mesh.pointEdges()[meshPointI];
                        forAll(pEdges, pEI)
                        {
                            if (cellMeshEdges.found(pEdges[pEI]))
                            {
                                edge e = mesh.edges()[pEdges[pEI]];
                                label otherPt =
                                    (e[0] ==  meshPointI ? e[1] : e[0]);
                                if (pointLevel[otherPt] <= cLevel)
                                {
                                    nFound++;
                                }
                                if (pointLevel[otherPt] > cLevel)
                                {
                                    nFound2++;
                                }
                            }
                        }

                        if (nFound == 2 || nFound2 == 4)
                        {
                            stationaryPoints[meshPointI] = false;
                        }
                    }
                    else if (pointLevel[meshPointI] < cLevel)
                    {
                        stationaryPoints[meshPointI] = false;
                    }
                }
            }
        }

        forAll(mesh.points(), pointI)
        {
            if
            (
                !(isFixedPoint.get(pointI) == 0
                  || isFixedPatchPoint.get(pointI) == 0)
            )
            {
                labelList pCells = mesh.pointCells()[pointI];

                forAll(pCells, pCI)
                {
                    label cellI = pCells[pCI];
                    cell c = mesh.cells()[cellI];

                    labelList clabels = c.labels(mesh.faces());
                    label cLevel = cellLevel[cellI];

                    label nAnchors = 0;
                    forAll(clabels, cPtI)
                    {
                        label meshPointI = clabels[cPtI];

                        if (pointLevel[meshPointI] <= cLevel)
                        {
                            nAnchors++;
                        }
                    }

                    if (nAnchors < 8)
                    {
                        forAll(clabels, cli)
                        {
                            stationaryPoints[clabels[cli]] = true;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            stationaryPoints,
            orEqOp<bool>(),
            false
        );

        forAll(stationaryPoints, pointI)
        {
            if (stationaryPoints[pointI])
            {
                isFixedPointExtended.set(pointI, 1);
            }
        }
    }

    vectorField volCentre(mesh.nPoints(), vector::zero);
    scalarField volSum(mesh.nPoints(), 0.);

    PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
    pointField pDisp(mesh.nPoints(), vector::zero);
    labelList nSet(mesh.nPoints(), 0);

    pointField newCellCentres(mesh.nCells(), vector::zero);
    scalarField newCellVolumes(mesh.nCells(), 0);
    pointField newFaceCentres(mesh.nFaces(), vector::zero);
    vectorField newFaceAreas(mesh.nFaces(), vector::zero);
    vectorField cEst(mesh.nCells(), vector::zero);
    labelField nCellFaces(mesh.nCells(), 0);

    // Perform fixed number of smoothing iterations
    for (label i = 0; i < nDispSmooth; i++)
    {
        updateMeshGeometric
        (
            mesh,
            newPoints,
            newFaceCentres,
            newFaceAreas,
            newCellCentres,
            newCellVolumes,
            cEst,
            nCellFaces
        );

        if (i < nDispSmooth/2)
        {
            volCentre = vector::zero;
            volSum =  0.;

            forAll(mesh.points(),pointI)
            {
                if (isFixedPoint.get(pointI) == 0)
                {
                    const labelList& pCells = mesh.pointCells()[pointI];
                    forAll(pCells, cellI)
                    {
                        label  level = cellLevel[pCells[cellI]];
                        scalar len = edge0Len / pow(2., level);

                        scalar weight = 0.;
                        if (isFixedPointExtended.get(pointI) == 0)
                        {
                            weight = sqrt(newCellVolumes[pCells[cellI]])
                                / len;
                        }
                        else
                        {
                            weight = newCellVolumes[pCells[cellI]]
                                / len;
                        }

                        volCentre[pointI] +=  newCellCentres[pCells[cellI]]
                            * weight;
                        volSum[pointI] += weight;
                    }
                }
            }

            forAll(pp.meshPoints(), pointI)
            {
                const label meshPointI = pp.meshPoints()[pointI];

                if (isFixedPatchPoint.get(meshPointI) == 0)
                {
                    const labelList& pFaces = pp.pointFaces()[pointI];
                    forAll(pFaces, faceI)
                    {
                        const label pFaceI = pFaces[faceI];
                        const label meshFaceI = pp.addressing()[pFaceI];
                        const label own = mesh.faceOwner()[meshFaceI];
                        label  level = cellLevel[own];
                        scalar len = edge0Len / pow(2., level);
                        scalar weight = sqrt(mag(newFaceAreas[meshFaceI]))
                            / len;

                        volCentre[meshPointI] +=  newFaceCentres[meshFaceI]
                            * weight;
                        volSum[meshPointI] += weight;
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                volCentre,
                plusEqOp<vector>(),
                vector::zero        // null value
             );

            syncTools::syncPointList
            (
                mesh,
                volSum,
                plusEqOp<scalar>(),
                scalar(0.)      // null value
             );
            point weightedAverage(vector::zero);

            forAll(mesh.points(), pointI)
            {
                if
                (
                    isFixedPoint.get(pointI) == 0
                    || isFixedPatchPoint.get(pointI) == 0
                 )
                {
                    weightedAverage =
                        volCentre[pointI]/ (volSum[pointI] + SMALL);

                    newPoints[pointI] =
                        0.5*(newPoints[pointI] + weightedAverage);
                }
            }
        }

        if (flattenPoints)
        {
            Field<pointField> nbrPts(mesh.nPoints(), pointField());
            pDisp = vector::zero;
            nSet = 0;

            forAll(nbrPts, pointI)
            {
                labelList pEdges = mesh.pointEdges()[pointI];
                DynamicList<point> nP(pEdges.size());

                forAll(pEdges, pEI)
                {
                    if (isMasterEdge[pEdges[pEI]])
                    {
                        edge e = mesh.edges()[pEdges[pEI]];
                        label otherPt0 = (e[0] ==  pointI ? e[1] : e[0]);
                        nP.append(newPoints[otherPt0]);
                    }
                }
                nbrPts[pointI] = nP.shrink();
            }

            syncTools::syncPointList
            (
                mesh,
                nbrPts,
                meshRefinement::pointFieldCombine(),
                pointField()      // null value
             );

            forAll(mesh.points(), pointI)
            {
                if (isFixedPointExtended.get(pointI) == 0)
                {
                    pointField nbrs = nbrPts[pointI];

                    boolList edgeSet(nbrs.size(), false);
                    forAll(nbrs, nbrsI)
                    {
                        if (edgeSet[nbrsI])
                        {
                            continue;
                        }
                        vector eVec =  nbrs[nbrsI] - newPoints[pointI];
                        scalar eMag = mag(eVec);

                        if (eMag < SMALL)
                        {
                            continue;
                        }
                        eVec /= (eMag + SMALL);

                        label alignedEdge = -1;
                        scalar maxAlign = -GREAT;

                        forAll(nbrs, nbrsJ)
                        {
                            if (nbrsI != nbrsJ)
                            {
                                vector eVecN =  nbrs[nbrsJ] - newPoints[pointI];
                                scalar eMagN = mag(eVecN);

                                if (eMagN < SMALL)
                                {
                                    continue;
                                }
                                eVecN /= (eMagN + SMALL);
                                if (mag(eVec & eVecN) > maxAlign)
                                {
                                    alignedEdge = nbrsJ;
                                    maxAlign = mag(eVec & eVecN);
                                }
                            }
                        }

                        if (maxAlign > 0.866)
                        {
                            edgeSet[nbrsI] = true;
                            edgeSet[alignedEdge] = true;
                            linePointRef line(nbrs[nbrsI], nbrs[alignedEdge]);
                            pointHit hit = line.nearestDist(newPoints[pointI]);

                            if (hit.hit())
                            {
                                pDisp[pointI] += hit.hitPoint()
                                    - newPoints[pointI];
                                nSet[pointI]++;
                            }
                        }
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                pDisp,
                plusEqOp<point>(),         // combine op
                vector::zero     // null value
             );

            syncTools::syncPointList
            (
                mesh,
                nSet,
                plusEqOp<label>(),         // combine op
                label(0)     // null value
             );


            forAll(pDisp, pointI)
            {
                if (nSet[pointI] > 0)
                {
                    pDisp[pointI] /= nSet[pointI];
                    newPoints[pointI] += pDisp[pointI];
                }
            }
        }
    }

    displacement = newPoints-mesh.points();

    Info<< "Finished smoothing displaced mesh in = "
        << mesh.time().elapsedCpuTime() - smoothTime << " s." << endl;
}


Foam::labelList Foam::medialAxisMeshMover::getSlipBCs
(
    const pointVectorField& fld
)
{
    DynamicList<label> slipPatchIDs;
    forAll(fld.boundaryField(), patchI)
    {
        const pointPatchField<vector>& patchFld =
            fld.boundaryField()[patchI];

        if
        (
            isA<slipPointPatchField<vector>>(patchFld)
            || isA<fixedNormalSlipPointPatchField<vector>>(patchFld)
        )
        {
            slipPatchIDs.append(patchI);
        }
    }
    return labelList(slipPatchIDs, true);
}


// Tries and find a medial axis point. Done by comparing vectors to nearest
// wall point for both vertices of edge.
bool Foam::medialAxisMeshMover::isMaxEdge
(
    const List<pointData>& pointWallDist,
    const label edgeI,
    const scalar minCos
) const
{
    const pointField& points = mesh().points();

    // Do not mark edges with one side on moving wall.

    const edge& e = mesh().edges()[edgeI];

    vector v0(points[e[0]] - pointWallDist[e[0]].origin());
    scalar magV0(mag(v0));

    if (magV0 < SMALL)
    {
        return false;
    }

    vector v1(points[e[1]] - pointWallDist[e[1]].origin());
    scalar magV1(mag(v1));

    if (magV1 < SMALL)
    {
        return false;
    }


    //- Detect based on vector to nearest point differing for both endpoints
    //v0 /= magV0;
    //v1 /= magV1;
    //
    //// Test angle.
    //if ((v0 & v1) < minCos)
    //{
    //    return true;
    //}
    //else
    //{
    //    return false;
    //}

    //- Detect based on extrusion vector differing for both endpoints
    //  the idea is that e.g. a sawtooth wall can still be extruded
    //  successfully as long as it is done all to the same direction.
    if ((pointWallDist[e[0]].v() & pointWallDist[e[1]].v()) < minCos)
    {
        vector origVec =
            pointWallDist[e[1]].origin() - pointWallDist[e[0]].origin();
        scalar magOrigVec = mag(origVec);

        if (magOrigVec > SMALL)
        {
            origVec /= magOrigVec;
            if
            (
                (origVec & pointWallDist[e[1]].v()) > 0
                && (origVec & pointWallDist[e[0]].v()) < 0
             )
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
            return false;
        }
    }
    else
    {
        return false;
    }
}


void Foam::medialAxisMeshMover::update(const dictionary& coeffDict)
{
    Info<< typeName
        << " : Calculating distance to Medial Axis ..." << endl;

    const pointField& points = mesh().points();

    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelList& meshPoints = pp.meshPoints();


    // Read a few parameters
    // ~~~~~~~~~~~~~~~~~~~~~

    //- Smooth surface normals
    const label nSmoothSurfaceNormals = coeffDict.lookupOrDefault<label>
    (
        "nSmoothSurfaceNormals",
        6
    );

    //- When is medial axis
    word angleKey = "minMedialAxisAngle";
    if (!coeffDict.found(angleKey))
    {
        // Backwards compatibility
        angleKey = "minMedianAxisAngle";
    }

    scalar minMedialAxisAngleCos = Foam::cos
    (
        degToRad
        (
            coeffDict.lookupOrDefault<scalar>
            (
                angleKey,
                90
            )
        )
    );

    //- Smooth internal normals
    const label nSmoothNormals = coeffDict.lookupOrDefault<label>
    (
        "nSmoothNormals",
        3
    );

    //- Number of edges walking out
    const label nMedialAxisIter = coeffDict.lookupOrDefault<label>
    (
        "nMedialAxisIter",
        mesh().globalData().nTotalPoints()
    );

    const scalar maxProjectionDist = coeffDict.lookupOrDefault<scalar>
    (
        "maxProjectionDistance",
        GREAT
     );

    const scalar maxCellDistortion = coeffDict.lookupOrDefault<scalar>
    (
        "maxCellDistortion",
        50
    );

    const Switch growZoneLayers = coeffDict.lookupOrDefault<Switch>
    (
        "growZoneLayers",
        false
    );

    const label medialRatioExp = coeffDict.lookupOrDefault<label>
    (
        "medialRatioExp",
        1
    );

    const Switch smoothBaffleEdgeNormals = coeffDict.lookupOrDefault<Switch>
    (
        "smoothBaffleEdgeNormals",
        false
    );

    // Predetermine mesh edges
    // ~~~~~~~~~~~~~~~~~~~~~~~

    // Precalulate (mesh) master point/edge (only relevant for shared pts/edges)
    const PackedBoolList isMeshMasterPoint(syncTools::getMasterPoints(mesh()));
    const PackedBoolList isMeshMasterEdge(syncTools::getMasterEdges(mesh()));
    // Precalculate meshEdge per pp edge
    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

    // Precalulate (patch) master point/edge
    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh(),
            meshPoints
        )
    );
    const PackedBoolList isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh(),
            meshEdges
        )
    );

    // Determine pointNormal
    // ~~~~~~~~~~~~~~~~~~~~~
    const vectorField& faceNormals = pp.faceNormals();
    pointField pointNormals(pp.nPoints(), vector::zero);
    {
        scalarField nPointFaces(pp.nPoints(), 0.);

        forAll(faceNormals, faceI)
        {
            const face& f = pp.localFaces()[faceI];

            forAll(f, fp)
            {
                scalar fA = mesh().magFaceAreas()[pp.addressing()[faceI]];
                pointNormals[f[fp]] += faceNormals[faceI]*fA;
                nPointFaces[f[fp]] += fA;
            }
        }

        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            pointNormals,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            nPointFaces,
            plusEqOp<scalar>(),
            scalar(0.)          // null value
        );

        forAll(pointNormals, i)
        {
            pointNormals[i] /= nPointFaces[i];
        }
    }
    pointNormals /= (mag(pointNormals) + SMALL);

    // Do not smooth on points marked as external
    boolList grownUpPoints(mesh().nPoints(),false);
    boolList isExternalEdge(mesh().nEdges(),false);

    const polyBoundaryMesh& patches = mesh().boundaryMesh();
    labelHashSet grownPatches(slipPatchIDs_);

    autoPtr<indirectPrimitivePatch> grownPPPtr
    (
        meshRefinement::makePatch
        (
            mesh(),
            slipPatchIDs_
        )
    );
    indirectPrimitivePatch& grownPP = grownPPPtr();

    boolList grownUp(mesh().nPoints(), false);
    boolList grownUpFaces(mesh().nFaces(), false);
    forAll(grownPP.meshPoints(), i)
    {
        grownUp[grownPP.meshPoints()[i]] = true;
    }

    forAll(grownPP, i)
    {
        grownUpFaces[grownPP.addressing()[i]] = true;
    }

    syncTools::syncPointList
    (
        mesh(),
        grownUp,
        orEqOp<bool>(),
        false           // null value
    );


    // Precalculate grown up patches  meshEdge per pp edge
    labelList grownUpMeshEdges(grownPP.nEdges());

    forAll(grownUpMeshEdges, patchEdgeI)
    {
        const edge& e = grownPP.edges()[patchEdgeI];

        label v0 = grownPP.meshPoints()[e[0]];
        label v1 = grownPP.meshPoints()[e[1]];
        grownUpMeshEdges[patchEdgeI] = meshTools::findEdge
        (
            mesh().edges(),
            mesh().pointEdges()[v0],
            v0,
            v1
        );
    }

    // unique edge centre to face centre vector.
    vectorField edgeNormal(mesh().nEdges(), vector(GREAT, GREAT, GREAT));

    const labelListList& edgeFaces = grownPP.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];
        label meshEdgeI = grownUpMeshEdges[edgeI];

        forAll(eFaces, i)
        {
            vector dir =
                mesh().faceAreas()[grownPP.addressing()[eFaces[i]]];
            dir /= (mag(dir) + SMALL);

            uniqueEdgeNormal()
            (
                edgeNormal[meshEdgeI],
                dir
            );
        }
    }

    syncTools::syncEdgeList
    (
        mesh(),
        edgeNormal,
        uniqueEdgeNormal(),
        vector(GREAT, GREAT, GREAT)  // null value
    );

    boolList cornerPts(mesh().nPoints(), false);
    {
        labelList nExternalEdge(mesh().nEdges(), 0);

        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(pp.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            const labelList& eFaces = edgeFaces[edgeI];
            nExternalEdge[meshEdgeI] = eFaces.size();
        }

        syncTools::syncEdgeList
        (
            mesh(),
            nExternalEdge,
            plusEqOp<label>(),
            label(0)              // null value
        );

        vectorField grownUpNormals(mesh().nPoints(), vector(GREAT, GREAT, GREAT));

        forAll(mesh().edges(), edgeI)
        {
            if (nExternalEdge[edgeI] == 1 || nExternalEdge[edgeI] == 2)
            {
                const labelList& meshEdgeFaces = mesh().edgeFaces()[edgeI];

                forAll(meshEdgeFaces, k)
                {
                    label faceI = meshEdgeFaces[k];
                    if (!mesh().isInternalFace(faceI))
                    {
                        label patchI = patches.whichPatch(faceI);

                        if (grownPatches.found(patchI))
                        {
                            isExternalEdge[edgeI] = true;
                            const edge& e = mesh().edges()[edgeI];
                            point norm  = mesh().faceAreas()[faceI];
                            norm /= (mag(norm) + SMALL);

                            grownUpNormals[e[0]] = norm;
                            grownUpNormals[e[1]] = norm;
                            grownUpPoints[e[0]] = true;
                            grownUpPoints[e[1]] = true;
                            break;
                        }
                    }
                }
            }
        }

        syncTools::syncPointList
        (
            mesh(),
            grownUpPoints,
            orEqOp<bool>(),
            false              // null value
        );

        syncTools::syncPointList
        (
            mesh(),
            grownUpNormals,
            minMagEqOp(),
            vector(GREAT, GREAT, GREAT)
        );

        PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh()));
        Field<pointField> nbrBoundaryPts(pp.meshPoints().size(), pointField());
        forAll(pp.meshPoints(), pointi)
        {
            label meshPointI = pp.meshPoints()[pointi];

            if (grownUpPoints[meshPointI])
            {
                labelList pEdges = mesh().pointEdges()[meshPointI];
                DynamicList<point> nP(pEdges.size());
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (isMasterEdge[edgeI] && isExternalEdge[edgeI])
                    {
                        edge e = mesh().edges()[pEdges[pEI]];
                        label otherPt0 = (e[0] ==  meshPointI ? e[1] : e[0]);
                        nP.append(mesh().points()[otherPt0]);
                    }
                }
                nbrBoundaryPts[pointi] = nP.shrink();
            }
        }

        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            nbrBoundaryPts,
            meshRefinement::pointFieldCombine(),
            pointField()      // null value
        );

        forAll(nbrBoundaryPts, pointi)
        {
            if (nbrBoundaryPts[pointi].size() == 2)
            {
                label meshPointI = pp.meshPoints()[pointi];
                point pt = mesh().points()[meshPointI];
                vector eV0 = nbrBoundaryPts[pointi][0] - pt;
                eV0 /= mag(eV0) + SMALL;
                vector eV1 = nbrBoundaryPts[pointi][1] - pt;
                eV1 /= mag(eV1) + SMALL;

                if ((eV0 & eV1) > -0.939)
                {
                    cornerPts[meshPointI] = true;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh(),
            cornerPts,
            orEqOp<bool>(),
            false              // null value
         );

        forAll(pointNormals, pointI)
        {
            label meshPointI = pp.meshPoints()[pointI];
            if (cornerPts[meshPointI])
            {
                const labelList& pEdges = mesh().pointEdges()[meshPointI];
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (!isExternalEdge[edgeI])
                    {
                        edge e = mesh().edges()[edgeI];
                        label otherPt = (e[0] ==  meshPointI ? e[1] : e[0]);

                        if (grownUp[otherPt])
                        {
                            const labelList& meshEdgeFaces = mesh().edgeFaces()[edgeI];
                            forAll(meshEdgeFaces, eFI)
                            {
                                label facei = meshEdgeFaces[eFI];
                                if (grownUpFaces[facei])
                                {
                                    vector dir =
                                        mesh().faceAreas()[facei];
                                    dir /= (mag(dir) + SMALL);

                                    if (mag(edgeNormal[edgeI] &  dir) < 0.939)
                                    {
                                        vector eVec = mesh().points()[meshPointI]
                                            -mesh().points()[otherPt];
                                        eVec /= (mag(eVec) + SMALL);
                                        pointNormals[pointI] = eVec;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else if (grownUpPoints[meshPointI])
            {
                point pNorm = grownUpNormals[meshPointI];
                //project point normal orthogonally onto grown up plane
                pointNormals[pointI] = pointNormals[pointI]
                    - (pointNormals[pointI] & pNorm) * pNorm;
                pointNormals[pointI] /= (mag(pointNormals[pointI]) + SMALL);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        pointNormals,
        minMagEqOp(),
        vector(GREAT, GREAT, GREAT)          // null value
    );

    // Find corner points to stop smoothing of patch normals at these locations
    vectorField pGrownUpNormal(meshPoints.size(), vector(GREAT, GREAT, GREAT));

    forAll(pp.localPoints(), pI)
    {
        label meshPointI = pp.meshPoints()[pI];

        if (grownUpPoints[meshPointI])
        {
            forAll(mesh().pointFaces()[meshPointI], fI)
            {
                label meshFaceI = mesh().pointFaces()[meshPointI][fI];
                if (!mesh().isInternalFace(meshFaceI))
                {
                    label patchI = patches.whichPatch(meshFaceI);

                    if (grownPatches.found(patchI))
                    {
                        vector fA = mesh().faceAreas()[meshFaceI];
                        fA /= (mag(fA) + SMALL);

                        nomalsCombine()
                        (
                            pGrownUpNormal[pI],
                            fA
                        );
                    }
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        pGrownUpNormal,
        nomalsCombine(),
        vector(GREAT, GREAT, GREAT)          // null value
    );

    boolList staticPoints = grownUpPoints;
    boolList staticEdges = isExternalEdge;

    if (!smoothBaffleEdgeNormals)
    {
        //Extact baffle edges
        const boolList excludeFaces(pp.size(), false);
        edgeClassification eClass
        (
            mesh(),
            mesh().points(),
            pp,
            meshEdges,
            excludeFaces,
            -0.5,
            -0.5
        );
        const List<Tuple2<edgeClassification::edgeType,scalar>>&
            eType = eClass.edgeTypes();

        pointField baffleDir(pp.nPoints(), vector::zero);
        scalarField nPointFaces(pp.nPoints(), 0.);

        forAll(eType, edgeI)
        {
            if (eType[edgeI].first() == edgeClassification::CONVEX)
            {
                label meshEdgeI = meshEdges[edgeI];
                edge e = mesh().edges()[meshEdgeI];

                staticPoints[e[0]] = true;
                staticPoints[e[1]] = true;
                staticEdges[meshEdgeI] = true;
                const labelList& ppEdgeFaces = pp.edgeFaces()[edgeI];
                forAll(ppEdgeFaces, ppEFI)
                {
                    label meshFaceI = pp.addressing()[ppEdgeFaces[ppEFI]];

                    face f = mesh().faces()[meshFaceI];
                    vector eVec = e.vec(mesh().points());
                    label fp = findIndex(f, e[0]);

                    vector eNorm = (eVec ^ mesh().faceAreas()[meshFaceI]);
                    eNorm /= (mag(eNorm) + SMALL);
                    scalar eMag = mag(eVec);

                    if (f[f.fcIndex(fp)] != e[1])
                    {
                        eNorm = -eNorm;
                    }

                    edge ppEdge = pp.edges()[edgeI];

                    baffleDir[ppEdge[0]] += eNorm*eMag;
                    baffleDir[ppEdge[1]] += eNorm*eMag;
                    nPointFaces[ppEdge[0]] += eMag;
                    nPointFaces[ppEdge[1]] += eMag;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh(),
            staticEdges,
            orEqOp<bool>(),
            false              // null value
        );

        syncTools::syncPointList
        (
            mesh(),
            staticPoints,
            orEqOp<bool>(),
            false              // null value
        );

        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            baffleDir,
            plusEqOp<vector>(),
            vector::zero        // null value
        );

        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            nPointFaces,
            plusEqOp<scalar>(),
            scalar(0.)          // null value
        );

        forAll(baffleDir, i)
        {
            if (nPointFaces[i] > 0 && mag(baffleDir[i]) > SMALL)
            {
                vector dir = baffleDir[i] / nPointFaces[i];
                dir /= mag(dir);
                pointNormals[i] = -dir;
            }
        }
    }

    // Smooth patch normal vectors
    fieldSmoother_.smoothPatchNormals
    (
        nSmoothSurfaceNormals,
        isPatchMasterPoint,
        isPatchMasterEdge,
        staticEdges,//isExternalEdge,
        staticPoints,//grownUpPoints,
        cornerPts,
        meshEdges,
        pp,
        pointNormals
    );

    if (debug)
    {
        simpleVTKWriter grownVTK
        (
            pp.localFaces(),
            pp.localPoints()
        );

        grownVTK.addPointData("pointNormals", pointNormals);
        grownVTK.write("grownPatchData.vtk");
    }

    // Calculate distance to pp points
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DynamicList<label> fixedPoints(mesh().nPoints());
    forAll(meshPoints, i)
    {
        fixedPoints.append(meshPoints[i]);
    }

    // Distance to wall
    List<pointData> pointWallDist(mesh().nPoints());

    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;

    boolList emptyRegions(mesh().nPoints(), false);

    regionSplit cellRegion(mesh());

    // 1. Calculate distance to points where displacement is specified.
    {
        // Seed data.
        DynamicList<pointData> wallInfo(meshPoints.size());
        DynamicList<label> wallPoints(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label meshPointI = meshPoints[patchPointI];
            wallPoints.append(meshPointI);
            wallInfo.append
            (
                pointData
                (
                    points[meshPointI],
                    0.0,
                    mag(pointDisplacement_[meshPointI]),
                    pointNormals[patchPointI]
                )
            );
        }

        if (cellRegion.nRegions() > 1)
        {
            Info<<"Multiregion mesh:"
                <<" Checking seeding for wave calculation"<<endl;

            const polyBoundaryMesh& patches = mesh().boundaryMesh();
            labelList numberSet(cellRegion.nRegions(),0);
            forAll(pp.addressing(), faceI)
            {
                label meshFaceI = pp.addressing()[faceI];
                label own = mesh().faceOwner()[meshFaceI];
                label regionI = cellRegion[own];
                numberSet[regionI]++;
            }
            Pstream::listCombineReduce(numberSet, plusOp<label>());

            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (!isA<processorPolyPatch>(pp))
                {
                    forAll(pp, faceI)
                    {
                        label meshFaceI = pp.start() + faceI;
                        label own = mesh().faceOwner()[meshFaceI];
                        label regionI = cellRegion[own];
                        if (numberSet[regionI] == 0)
                        {
                            face f = mesh().faces()[meshFaceI];
                            forAll(f, fI)
                            {
                                label meshPointI = f[fI];
                                wallPoints.append(meshPointI);
                                wallInfo.append
                                (
                                    pointData
                                    (
                                        points[meshPointI],
                                        0.0,
                                        0.0,
                                        vector::one
                                     )
                                );
                                fixedPoints.append(meshPointI);
                                emptyRegions[meshPointI] =  true;
                            }
                        }
                    }
                }
            }
        }

        wallInfo.shrink();
        wallPoints.shrink();

        List<pointData> edgeWallDist(mesh().nEdges());
        // Do all calculations
        PointEdgeWave<pointData> wallDistCalc
        (
            mesh(),
            wallPoints,
            wallInfo,

            pointWallDist,
            edgeWallDist,
            0,
            dummyTrackData
        );
        wallDistCalc.iterate(nMedialAxisIter);

        label nUnvisit = returnReduce
        (
            wallDistCalc.getUnsetPoints(),
            sumOp<label>()
        );

        if (nUnvisit > 0)
        {
            WarningInFunction
                << "Walking did not visit all points." << nl
                << "    Did not visit " << nUnvisit
                << " out of " << mesh().globalData().nTotalPoints()
                << " points. This is not necessarily a problem" << nl
                << "    and might be due to faceZones splitting of part"
                << " of the domain." << nl << endl;
        }
    }

    List<pointEdgePoint> pointMedialDist(mesh().nPoints());

    // 2. Find points with max distance and transport information back to
    //    wall.
    {
        List<pointEdgePoint> edgeMedialDist(mesh().nEdges());

        // Seed point data.
        DynamicList<pointEdgePoint> maxInfo(meshPoints.size());
        DynamicList<label> maxPoints(meshPoints.size());

        //If maximum projection distance set then set to medial points
        //if distance is greater than maxProjectionDist
        if (maxProjectionDist < GREAT)
        {
            forAll(pointWallDist, pointI)
            {
                if (sqrt(pointWallDist[pointI].distSqr()) > maxProjectionDist)
                {
                    maxPoints.append(pointI);
                    maxInfo.append
                    (
                        pointEdgePoint
                        (
                            points[pointI],
                            0.0
                        )
                    );
                }
            }
        }

        //set maximum points if displacement to cell size becomes too large
        {
            const labelIOList& pointLevel =
                mesh().lookupObject<labelIOList>("pointLevel");
            const scalar edge0Len = readScalar
            (
                coeffDict.lookup("edge0Len")
             );

            forAll(pointWallDist, pointI)
            {
                if (pointLevel[pointI] < 0)
                {
                    continue;
                }
                scalar thickness = pointWallDist[pointI].s();
                scalar dx = edge0Len/(1<<pointLevel[pointI]);

                if (thickness/dx > maxCellDistortion)
                {
                    maxPoints.append(pointI);
                    maxInfo.append
                    (
                        pointEdgePoint
                        (
                            points[pointI],
                            0.0
                        )
                    );
                }
            }
        }

        // 1. Medial axis points

        const edgeList& edges = mesh().edges();

        forAll(edges, edgeI)
        {
            const edge& e = edges[edgeI];

            if
            (
                !pointWallDist[e[0]].valid(dummyTrackData)
             || !pointWallDist[e[1]].valid(dummyTrackData)
            )
            {
                // Unvisited point. See above about nUnvisit warning
            }
            else if (isMaxEdge(pointWallDist, edgeI, minMedialAxisAngleCos))
            {
                // Both end points of edge have very different nearest wall
                // point. Mark both points as medial axis points.

                // Approximate medial axis location on edge.
                //const point medialAxisPt = e.centre(points);
                vector eVec = e.vec(points);
                scalar eMag = mag(eVec);
                if (eMag > VSMALL)
                {
                    eVec /= eMag;

                    // Calculate distance along edge
                    const point& p0 = points[e[0]];
                    const point& p1 = points[e[1]];
                    scalar dist0 = (p0-pointWallDist[e[0]].origin()) & eVec;
                    scalar dist1 = (pointWallDist[e[1]].origin()-p1) & eVec;
                    scalar s = 0.5*(dist1+eMag+dist0);

                    point medialAxisPt;
                    if (s <= dist0)
                    {
                        medialAxisPt = p0;
                    }
                    else if (s >= dist0+eMag)
                    {
                        medialAxisPt = p1;
                    }
                    else
                    {
                        medialAxisPt = p0+(s-dist0)*eVec;
                    }

                    forAll(e, ep)
                    {
                        label pointI = e[ep];
                        maxPoints.append(pointI);
                        maxInfo.append
                        (
                            pointEdgePoint
                            (
                                medialAxisPt,   //points[pointI],
                                magSqr(points[pointI]-medialAxisPt)//0.0
                            )
                        );
                    }
                }
            }
        }

        if (!growZoneLayers)
        {
            // 2. Zoned points
            forAll(mesh().faceZones(), zoneI)
            {
                const faceZone& fz = mesh().faceZones()[zoneI];
                forAll(fz, faceI)
                {
                    if (mesh().isInternalFace(fz[faceI]))
                    {
                        const face& f = mesh().faces()[fz[faceI]];
                        forAll(f, fp)
                        {
                            label meshPointI = f[fp];
                            if (!emptyRegions[meshPointI])
                            {
                                maxPoints.append(meshPointI);
                                maxInfo.append
                                (
                                    pointEdgePoint
                                    (
                                        points[meshPointI],
                                        0.0
                                    )
                                );
                            }
                        }
                    }
                }
            }
        }

        // 3. Seed non-adapt patches
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        labelHashSet adaptPatches(adaptPatchIDs_);
        labelHashSet grownPatches(slipPatchIDs_);
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if
            (
                !pp.coupled()
             && !isA<emptyPolyPatch>(pp)
             && !adaptPatches.found(patchI)
             && !grownPatches.found(patchI)
            )
            {
                const labelList& meshPoints = pp.meshPoints();

                forAll(meshPoints, i)
                {
                    label pointI = meshPoints[i];

                    if (!emptyRegions[pointI])
                    {
                        maxPoints.append(pointI);
                        maxInfo.append
                        (
                            pointEdgePoint
                            (
                                points[pointI],
                                0.0
                            )
                        );
                    }
                }
            }
        }

        if (slipPatchIDs_.size())
        {
            boolList baffleEdge(mesh().nEdges(), false);

            forAll(edgeFaces, edgeI)
            {
                const labelList& eFaces = edgeFaces[edgeI];
                label meshEdgeI = grownUpMeshEdges[edgeI];

                if (edgeNormal[meshEdgeI] != vector(GREAT, GREAT, GREAT))
                {
                    forAll(eFaces, i)
                    {
                        vector dir =
                            mesh().faceAreas()[grownPP.addressing()[eFaces[i]]];
                        dir /= (mag(dir) + SMALL);

                        if (mag(edgeNormal[meshEdgeI] &  dir) < 0.9848)
                        {
                            const edge e = mesh().edges()[meshEdgeI];

                            if
                            (
                                pointWallDist[e[0]].valid(dummyTrackData)
                                && pointWallDist[e[1]].valid(dummyTrackData)
                            )
                            {
                                vector disp0 = pointWallDist[e[0]].v();
                                disp0 /= (mag(disp0) + SMALL);
                                vector disp1 = pointWallDist[e[1]].v();
                                disp1 /= (mag(disp1) + SMALL);

                                vector eVec = e.vec(mesh().points());
                                eVec /= (mag(eVec) + SMALL);

                                if
                                (
                                    mag(eVec & disp0) < 0.984
                                    || mag(eVec & disp1) < 0.984
                                )
                                {
                                    baffleEdge[meshEdgeI] = true;
                                    break;
                                }
                            }
                            else
                            {
                                baffleEdge[meshEdgeI] = true;
                                break;
                            }
                        }
                    }
                }
            }

            syncTools::syncEdgeList
            (
                mesh(),
                baffleEdge,
                orEqOp<bool>(),
                false           // null value
             );

            forAll(mesh().edges(), edgeI)
            {
                if (baffleEdge[edgeI])
                {
                    const edge& e = mesh().edges()[edgeI];
                    label pointI = e[0];
                    if (!emptyRegions[pointI] && !cornerPts[pointI])
                    {
                        maxPoints.append(pointI);
                        maxInfo.append
                        (
                            pointEdgePoint
                            (
                                points[pointI],
                                0.0
                            )
                        );
                    }

                    pointI = e[1];
                    if (!emptyRegions[pointI] && !cornerPts[pointI])
                    {
                        maxPoints.append(pointI);
                        maxInfo.append
                        (
                            pointEdgePoint
                            (
                                points[pointI],
                                0.0
                            )
                        );
                    }
                }
            }
        }

        //if have not found medial point then set to maximum wall distance
        {
            boolList medialSetPts(mesh().nPoints(), false);
            forAll(maxPoints, ptI)
            {
                medialSetPts[maxPoints[ptI]] = true;
            }

            labelList numberRegionMedialPoints(cellRegion.nRegions(),0);
            List<scalar> maxRegionWallDist(cellRegion.nRegions(),0.);

            labelList pointRegions(mesh().nPoints(), -1);
            forAll(mesh().cells(), cellI)
            {
                label regionI = cellRegion[cellI];
                const cell c = mesh().cells()[cellI];

                labelList mPoints = c.labels(mesh().faces());

                forAll(mPoints, ptI)
                {
                    label pointI = mPoints[ptI];
                    pointRegions[pointI] = regionI;
                }
            }

            forAll(mesh().points(), pointI)
            {
                label regionI = pointRegions[pointI];
                if (medialSetPts[pointI])
                {
                    numberRegionMedialPoints[regionI]++;
                }
                maxRegionWallDist[regionI] = max
                (maxRegionWallDist[regionI],pointWallDist[pointI].distSqr());
            }

            Pstream::listCombineReduce(numberRegionMedialPoints, plusOp<label>());
            Pstream::listCombineReduce(maxRegionWallDist, maxOp<scalar>());

            bool setMedialPt = false;

            forAll(numberRegionMedialPoints, regionI)
            {
                if (numberRegionMedialPoints[regionI] == 0)
                {
                    Info<<"No medial point. Setting to maximum wall distance"
                        <<endl;
                    setMedialPt = true;
                }
            }

            if (setMedialPt)
            {
                forAll(mesh().points(), pointI)
                {
                    label regionI = pointRegions[pointI];

                    if (numberRegionMedialPoints[regionI] == 0)
                    {
                        if
                        (
                            pointWallDist[pointI].distSqr()
                            >= 0.99*maxRegionWallDist[regionI]
                            && !emptyRegions[pointI]
                         )
                        {
                            maxPoints.append(pointI);
                            maxInfo.append
                            (
                                pointEdgePoint
                                (
                                    points[pointI],
                                    0.0
                                 )
                            );
                        }
                    }
                }
            }
        }

        maxInfo.shrink();
        maxPoints.shrink();

        // Do all calculations
        PointEdgeWave<pointEdgePoint> medialDistCalc
        (
            mesh(),
            maxPoints,
            maxInfo,

            pointMedialDist,
            edgeMedialDist,
            0,
            dummyTrackData
        );
        medialDistCalc.iterate(2*nMedialAxisIter);

        // Extract medial axis distance as pointScalarField
        forAll(pointMedialDist, pointI)
        {
            medialDist_[pointI] = Foam::sqrt(pointMedialDist[pointI].distSqr());
            medialVec_[pointI] = pointMedialDist[pointI].origin();
        }
    }

    // Extract transported surface normals as pointVectorField
    forAll(dispVec_, i)
    {
        if (!pointWallDist[i].valid(dummyTrackData))
        {
            dispVec_[i] = vector(1, 0, 0);
        }
        else
        {
            dispVec_[i] = pointWallDist[i].v();
        }
    }

    label grownSz = grownPP.size();
    // 1. Calculate distance to grown up patch points
    if (returnReduce(grownSz, sumOp<label>()) > 0)
    {
        // Distance to gown edges
        List<PointData<vector>> pointExtDist(mesh().nPoints());
        // Seed point data.
        DynamicList<PointData<vector>> extInfo(meshPoints.size());
        DynamicList<label> extPoints(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label pointI = meshPoints[patchPointI];
            if (grownUpPoints[pointI] && !emptyRegions[pointI])
            {

                extPoints.append(pointI);
                extInfo.append
                (
                    PointData<vector>
                    (
                        points[pointI],
                        0.0,
                        pointNormals[patchPointI]     // surface normals
                    )
                );
            }
        }
        extPoints.shrink();
        extInfo.shrink();

        if (returnReduce(extPoints.size(), sumOp<label>()) > 0)
        {
            // Do all calculations
            List<PointData<vector>> edgeExtDist(mesh().nEdges());
            PointEdgeWave<PointData<vector>> wallDistCalc
            (
                mesh(),
                extPoints,
                extInfo,
                pointExtDist,
                edgeExtDist,
                0,
                dummyTrackData
            );
            wallDistCalc.iterate(2*nMedialAxisIter);

            forAll(mesh().points(), meshPointI)
            {
                if (grownUp[meshPointI] && !emptyRegions[meshPointI])
                {
                    dispVec_[meshPointI] = pointExtDist[meshPointI].data();
                    fixedPoints.append(meshPointI);
                }
            }
        }
    }

    fixedPoints.shrink();

    // Smooth normal vectors. Do not change normals on pp.meshPoints
    fieldSmoother_.smoothNormals
    (
        nSmoothNormals,
        isMeshMasterPoint,
        isMeshMasterEdge,
        fixedPoints,
        dispVec_
    );
    // Calculate ratio point medial distance to point wall distance
    forAll(medialRatio_, pointI)
    {
        if (!pointWallDist[pointI].valid(dummyTrackData))
        {
            medialRatio_[pointI] = 0.0;
        }
        else
        {
            scalar wDist2 = pointWallDist[pointI].distSqr();
            scalar mDist = medialDist_[pointI];

            if (wDist2 < sqr(SMALL) && mDist < SMALL)
            {
                medialRatio_[pointI] = 0.0;
            }
            else
            {
                medialRatio_[pointI] = mDist / (Foam::sqrt(wDist2) + mDist);
                if (medialRatioExp > 1)
                {
                    medialRatio_[pointI] =
                        pow(medialRatio_[pointI],medialRatioExp);
                }

                medialRatio_[pointI] =
                (
                    cos
                    (
                        (1. - medialRatio_[pointI])
                        * Foam::constant::mathematical::pi
                     ) + 1.
                ) / 2.;
            }
        }
    }

    if (debug)
    {
        Info<< typeName
            << " : Writing medial axis fields:" << nl
            << incrIndent
            << "ratio of medial distance to wall distance : "
            << medialRatio_.name() << nl
            << "distance to nearest medial axis           : "
            << medialDist_.name() << nl
            << "nearest medial axis location              : "
            << medialVec_.name() << nl
            << "normal at nearest wall                    : "
            << dispVec_.name() << nl
            << decrIndent << nl
            << endl;

        mesh().write();
        dispVec_.write();
        medialRatio_.write();
        medialDist_.write();
        medialVec_.write();
    }
}


bool Foam::medialAxisMeshMover::unmarkExtrusion
(
    const label patchPointI,
    pointField& patchDisp,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus
)
{
    if (extrudeStatus[patchPointI] == snappyLayerDriver::EXTRUDE)
    {
        extrudeStatus[patchPointI] = snappyLayerDriver::NOEXTRUDE;
        patchDisp[patchPointI] = Zero;
        return true;
    }
    else if (extrudeStatus[patchPointI] == snappyLayerDriver::EXTRUDEREMOVE)
    {
        extrudeStatus[patchPointI] = snappyLayerDriver::NOEXTRUDE;
        patchDisp[patchPointI] = Zero;
        return true;
    }
    else
    {
        return false;
    }
}


void Foam::medialAxisMeshMover::syncPatchDisplacement
(
    const scalarField& minThickness,
    pointField& patchDisp,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelList& meshPoints = pp.meshPoints();

//    label nChangedTotal = 0;

    while (true)
    {
        label nChanged = 0;

        // Sync displacement (by taking min)
        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            patchDisp,
            minMagSqrEqOp<vector>(),
            point::rootMax           // null value
        );

        // Unmark if displacement too small
        forAll(patchDisp, i)
        {
            if (mag(patchDisp[i]) < minThickness[i])
            {
                if (unmarkExtrusion(i, patchDisp, extrudeStatus))
                {
                    nChanged++;
                }
            }
        }

        //labelList syncPatchNLayers(patchNLayers);
        //
        //syncTools::syncPointList
        //(
        //    mesh(),
        //    meshPoints,
        //    syncPatchNLayers,
        //    minEqOp<label>(),
        //    labelMax            // null value
        //);
        //
        //// Reset if differs
        //// 1. take max
        //forAll(syncPatchNLayers, i)
        //{
        //    if (syncPatchNLayers[i] != patchNLayers[i])
        //    {
        //        if
        //        (
        //            unmarkExtrusion
        //            (
        //                i,
        //                patchDisp,
        //                patchNLayers,
        //                extrudeStatus
        //            )
        //        )
        //        {
        //            nChanged++;
        //        }
        //    }
        //}
        //
        //syncTools::syncPointList
        //(
        //    mesh(),
        //    meshPoints,
        //    syncPatchNLayers,
        //    maxEqOp<label>(),
        //    labelMin            // null value
        //);
        //
        //// Reset if differs
        //// 2. take min
        //forAll(syncPatchNLayers, i)
        //{
        //    if (syncPatchNLayers[i] != patchNLayers[i])
        //    {
        //        if
        //        (
        //            unmarkExtrusion
        //            (
        //                i,
        //                patchDisp,
        //                patchNLayers,
        //                extrudeStatus
        //            )
        //        )
        //        {
        //            nChanged++;
        //        }
        //    }
        //}

//        nChangedTotal += nChanged;

        if (!returnReduce(nChanged, sumOp<label>()))
        {
            break;
        }
    }

    //Info<< "Prevented extrusion on "
    //    << returnReduce(nChangedTotal, sumOp<label>())
    //    << " coupled patch points during syncPatchDisplacement." << endl;
}


// Stop layer growth where mesh wraps around edge with a
// large feature angle
void Foam::medialAxisMeshMover::
handleFeatureAngleLayerTerminations
(
    const scalar minCos,
    const labelList& meshEdges,
    const PackedList<1>& isConvexEdgePoint,
    const PackedList<1>& isConcaveEdgePoint,
    const vectorField& patchEdgeNormals,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp,
    label& nPointCounter
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();

    // Mark faces that have all points extruded
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nExtrudedFacePoints(pp.size(), 0);

    forAll(pp.localFaces(), faceI)
    {
        const face& f = pp.localFaces()[faceI];

        forAll(f, fp)
        {
            if (extrudeStatus[f[fp]] != snappyLayerDriver::NOEXTRUDE)
            {
                nExtrudedFacePoints[faceI]++;
            }
        }
    }

    //label nOldPointCounter = nPointCounter;

    // Detect situation where two featureedge-neighbouring faces are partly or
    // not extruded and the edge itself is extruded. In this case unmark the
    // edge for extrusion.

    forAll(pp.edgeFaces(), edgeI)
    {
        const labelList& eFaces = pp.edgeFaces()[edgeI];

        const edge& e = pp.edges()[edgeI];
        label v0 = e[0];
        label v1 = e[1];

        if
        (
            isConcaveEdgePoint.get(pp.meshPoints()[v0]) == 1
         && isConcaveEdgePoint.get(pp.meshPoints()[v1]) == 1
         && (isConvexEdgePoint.get(pp.meshPoints()[v0]) == 0
             || isConvexEdgePoint.get(pp.meshPoints()[v1]) == 0)
        )
        {
            if
            (
                extrudeStatus[v0] != snappyLayerDriver::NOEXTRUDE
             || extrudeStatus[v1] != snappyLayerDriver::NOEXTRUDE
            )
            {
                forAll(eFaces, faceI)
                {
                    if (nExtrudedFacePoints[eFaces[faceI]] == 1)
                    {
                        if
                        (
                            unmarkExtrusion
                            (
                                v0,
                                patchDisp,
                                extrudeStatus
                            )
                        )
                        {
                            nPointCounter++;
                        }
                        if
                        (
                            unmarkExtrusion
                            (
                                v1,
                                patchDisp,
                                extrudeStatus
                            )
                        )
                        {
                            nPointCounter++;
                        }
                    }
                }
            }
            else
            {
                continue;
            }
        }
        else if
        (
            isConcaveEdgePoint.get(pp.meshPoints()[v0]) == 0
         && isConcaveEdgePoint.get(pp.meshPoints()[v1]) == 0
         && (isConvexEdgePoint.get(pp.meshPoints()[v0]) == 1
             || isConvexEdgePoint.get(pp.meshPoints()[v1]) == 1)
        )
        {
            if
            (
                extrudeStatus[v0] != snappyLayerDriver::NOEXTRUDE
             || extrudeStatus[v1] != snappyLayerDriver::NOEXTRUDE
            )
            {
                forAll(eFaces, faceI)
                {
                    if (nExtrudedFacePoints[eFaces[faceI]] == 1)
                    {
                        if
                        (
                            unmarkExtrusion
                            (
                                v0,
                                patchDisp,
                                extrudeStatus
                            )
                        )
                        {
                            nPointCounter++;
                        }
                        if
                        (
                            unmarkExtrusion
                            (
                                v1,
                                patchDisp,
                                extrudeStatus
                            )
                        )
                        {
                            nPointCounter++;
                        }
                    }
                }
            }
            else
            {
                continue;
            }
        }


        if
        (
            extrudeStatus[v0] != snappyLayerDriver::NOEXTRUDE
         || extrudeStatus[v1] != snappyLayerDriver::NOEXTRUDE
        )
        {
            bool extrudeAll = true;
            bool triFace = false;
            forAll(eFaces, faceI)
            {
                if
                (
                    nExtrudedFacePoints[eFaces[faceI]] !=
                    pp.localFaces()[eFaces[faceI]].size()
                )
                {
                    extrudeAll = false;
                }
                if (pp.localFaces()[eFaces[faceI]].size() == 3)
                {
                    triFace = true;
                }
            }
            if (!extrudeAll)
            {
                const vector& n0 = pp.faceNormals()[eFaces[0]];
                scalar dotProd = (patchEdgeNormals[edgeI] & n0);

                if
                (
                    (triFace && dotProd < minCos)
                    || (!triFace && dotProd < 0.173)
                )
                {
                    if
                    (
                        unmarkExtrusion
                        (
                            v0,
                            patchDisp,
                            extrudeStatus
                        )
                    )
                    {
                        nPointCounter++;
                    }
                    if
                    (
                        unmarkExtrusion
                        (
                            v1,
                            patchDisp,
                            extrudeStatus
                        )
                    )
                    {
                        nPointCounter++;
                    }
                }
            }
        }
    }
}


// Find isolated islands (points, edges and faces and layer terminations)
// in the layer mesh and stop any layer growth at these points.
void Foam::medialAxisMeshMover::findIsolatedRegions
(
    const scalar minCosLayerTermination,
    const bool detectExtrusionIsland,
    const PackedBoolList& isPatchMasterEdge,
    const PackedList<1>& isConvexEdgePoint,
    const PackedList<1>& isConcaveEdgePoint,
    const labelList& meshEdges,
    const scalarField& minThickness,
    const vectorField& patchEdgeNormals,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp
) const
{
    const indirectPrimitivePatch& pp = adaptPatchPtr_();
    const labelListList& pointFaces = pp.pointFaces();
    const labelList& meshPoints = pp.meshPoints();

    Info<< typeName << " : Removing isolated regions ..." << nl
        << indent << "- if partially extruded faces make angle < "
        << Foam::radToDeg(Foam::acos(minCosLayerTermination)) <<  nl;
    if (detectExtrusionIsland)
    {
        Info<< indent << "- if exclusively surrounded by non-extruded points"
            << nl;
    }
    else
    {
        Info<< indent << "- if exclusively surrounded by non-extruded faces"
            << nl;
    }

    // Keep count of number of points unextruded
    label nPointCounter = 0;


    autoPtr<OBJstream> str;
    if (debug)
    {
        str.reset
        (
            new OBJstream
            (
                mesh().time().path()
              / "islandExcludePoints_"
              + mesh().time().timeName()
              + ".obj"
            )
        );
        Info<< typeName
            << " : Writing points surrounded by non-extruded points to "
            << str().name() << endl;
    }

    while (true)
    {
        // Stop layer growth where mesh wraps around edge with a
        // large feature angle
        if (minCosLayerTermination > -1)
        {
            handleFeatureAngleLayerTerminations
            (
                minCosLayerTermination,
                meshEdges,

                isConvexEdgePoint,
                isConcaveEdgePoint,
                patchEdgeNormals,

                extrudeStatus,
                patchDisp,
                nPointCounter
            );

            syncPatchDisplacement(minThickness, patchDisp, extrudeStatus);
        }

        // Detect either:
        // - point where all surrounding points are not extruded
        //   (detectExtrusionIsland)
        // or
        // - point where all the faces surrounding it are not fully
        //   extruded

        boolList keptPoints(pp.nPoints(), false);

        if (detectExtrusionIsland)
        {
            // Do not extrude from point where all neighbouring
            // points are not grown
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            labelList islandPoint(pp.size(), -1);
            forAll(pp, faceI)
            {
                const face& f = pp.localFaces()[faceI];

                forAll(f, fp)
                {
                    if (extrudeStatus[f[fp]] != snappyLayerDriver::NOEXTRUDE)
                    {
                        if (islandPoint[faceI] == -1)
                        {
                            // First point to extrude
                            islandPoint[faceI] = f[fp];
                        }
                        else if (islandPoint[faceI] != -2)
                        {
                            // Second or more point to extrude
                            islandPoint[faceI] = -2;
                        }
                    }
                }
            }

            // islandPoint:
            //  -1 : no point extruded on face
            //  -2 : >= 2 points extruded on face
            //  >=0: label of point extruded

            // Check all surrounding faces that I am the islandPoint
            forAll(pointFaces, patchPointI)
            {
                if (extrudeStatus[patchPointI] != snappyLayerDriver::NOEXTRUDE)
                {
                    const labelList& pFaces = pointFaces[patchPointI];

                    forAll(pFaces, i)
                    {
                        label faceI = pFaces[i];
                        if (islandPoint[faceI] != patchPointI)
                        {
                            keptPoints[patchPointI] = true;
                            break;
                        }
                    }
                }
            }
        }
        else
        {
            // Do not extrude from point where all neighbouring
            // faces are not grown
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            boolList extrudedFaces(pp.size(), true);
            forAll(pp.localFaces(), faceI)
            {
                const face& f = pp.localFaces()[faceI];
                forAll(f, fp)
                {
                    if (extrudeStatus[f[fp]] == snappyLayerDriver::NOEXTRUDE)
                    {
                        extrudedFaces[faceI] = false;
                        break;
                    }
                }
            }

            const labelListList& pointFaces = pp.pointFaces();

            forAll(keptPoints, patchPointI)
            {
                const labelList& pFaces = pointFaces[patchPointI];

                forAll(pFaces, i)
                {
                    label faceI = pFaces[i];
                    if (extrudedFaces[faceI])
                    {
                        keptPoints[patchPointI] = true;
                        break;
                    }
                }
            }
        }


        syncTools::syncPointList
        (
            mesh(),
            meshPoints,
            keptPoints,
            orEqOp<bool>(),
            false               // null value
        );

        label nChanged = 0;

        forAll(keptPoints, patchPointI)
        {
            if (!keptPoints[patchPointI])
            {
                if (unmarkExtrusion(patchPointI, patchDisp, extrudeStatus))
                {
                    nPointCounter++;
                    nChanged++;
                }
            }
        }

        //At layer terminations remove extrusion where truncated
        //edge alighned with neighbouring edge
        forAll(pp.localFaces(), faceI)
        {
            const face& f = pp.localFaces()[faceI];
            forAll(f, fp)
            {
                label next = f.fcIndex(fp);
                vector eVec = pp.localPoints()[f[next]]
                    -  pp.localPoints()[f[fp]];
                eVec /= (mag(eVec) + SMALL);

                if
                (
                    extrudeStatus[f[fp]] == snappyLayerDriver::NOEXTRUDE
                    && extrudeStatus[f[next]] == snappyLayerDriver::EXTRUDE
                )
                {
                    label next2 = f.fcIndex(next);
                    vector eVecNext = pp.localPoints()[f[next2]]
                        -  pp.localPoints()[f[next]];
                    eVecNext /= (mag(eVecNext) + SMALL);

                    if ((eVec & eVecNext) > 0.9848)
                    {
                        if
                        (
                            unmarkExtrusion
                            (
                                f[next],
                                patchDisp,
                                extrudeStatus
                            )
                         )
                        {
                            nPointCounter++;
                            nChanged++;
                        }
                    }
                }
                else if
                (
                    extrudeStatus[f[next]] == snappyLayerDriver::NOEXTRUDE
                    && extrudeStatus[f[fp]] == snappyLayerDriver::EXTRUDE
                )
                {
                    label prev = f.rcIndex(fp);
                    vector eVecPrev = pp.localPoints()[f[fp]]
                        -  pp.localPoints()[f[prev]];
                    eVecPrev /= (mag(eVecPrev) + SMALL);

                    if ((eVec & eVecPrev) > 0.9848)
                    {
                        if
                        (
                            unmarkExtrusion
                            (
                                f[fp],
                                patchDisp,
                                extrudeStatus
                            )
                         )
                        {
                            nPointCounter++;
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

    const edgeList& edges = pp.edges();


    // Count number of mesh edges using a point
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList isolatedPoint(pp.nPoints(),0);

    forAll(edges, edgeI)
    {
        if (isPatchMasterEdge[edgeI])
        {
            const edge& e = edges[edgeI];

            label v0 = e[0];
            label v1 = e[1];

            if (extrudeStatus[v1] != snappyLayerDriver::NOEXTRUDE)
            {
                isolatedPoint[v0] += 1;
            }
            if (extrudeStatus[v0] != snappyLayerDriver::NOEXTRUDE)
            {
                isolatedPoint[v1] += 1;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        isolatedPoint,
        plusEqOp<label>(),
        label(0)        // null value
    );

    // stop layer growth on isolated faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    forAll(pp, faceI)
    {
        const face& f = pp.localFaces()[faceI];
        bool failed = false;
        forAll(f, fp)
        {
            if (isolatedPoint[f[fp]] > 2)
            {
                failed = true;
                break;
            }
        }
        bool allPointsExtruded = true;
        if (!failed)
        {
            forAll(f, fp)
            {
                if (extrudeStatus[f[fp]] == snappyLayerDriver::NOEXTRUDE)
                {
                    allPointsExtruded = false;
                    break;
                }
            }

            if (allPointsExtruded)
            {
                forAll(f, fp)
                {
                    if
                    (
                        unmarkExtrusion
                        (
                            f[fp],
                            patchDisp,
                            extrudeStatus
                        )
                    )
                    {
                        nPointCounter++;
                        if (str.valid())
                        {
                            str().write(pp.points()[meshPoints[f[fp]]]);
                        }
                    }
                }
            }
        }
    }

    reduce(nPointCounter, sumOp<label>());
    Info<< typeName
        << " : Number of isolated points extrusion stopped : "<< nPointCounter
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::medialAxisMeshMover::medialAxisMeshMover
(
    const dictionary& dict,
    const List<labelPair>& baffles,
    pointVectorField& pointDisplacement
)
:
    externalDisplacementMeshMover(dict, baffles, pointDisplacement),
    slipPatchIDs_(getSlipBCs(pointDisplacement)),
    adaptPatchIDs_(getFixedValueBCs(pointDisplacement)),
    adaptPatchPtr_(getPatch(mesh(), adaptPatchIDs_)),
    scale_
    (
        IOobject
        (
            "scale",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("scale", dimless, 1.0)
    ),
    oldPoints_(mesh().points()),
    meshMover_
    (
        const_cast<polyMesh&>(mesh()),
        const_cast<pointMesh&>(pMesh()),
        adaptPatchPtr_(),
        pointDisplacement,
        scale_,
        oldPoints_,
        adaptPatchIDs_,
        dict
    ),
    fieldSmoother_(mesh()),
    dispVec_
    (
        IOobject
        (
            "dispVec",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pMesh(),
        dimensionedVector("dispVec", dimLength, Zero)
    ),
    medialRatio_
    (
        IOobject
        (
            "medialRatio",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pMesh(),
        dimensionedScalar("medialRatio", dimless, 0.0)
    ),
    medialDist_
    (
        IOobject
        (
            "pointMedialDist",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pMesh(),
        dimensionedScalar("pointMedialDist", dimLength, 0.0)
    ),
    medialVec_
    (
        IOobject
        (
            "medialVec",
            pointDisplacement.time().timeName(),
            pointDisplacement.db(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pMesh(),
        dimensionedVector("medialVec", dimLength, Zero)
    )
{
    update(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::medialAxisMeshMover::~medialAxisMeshMover()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::medialAxisMeshMover::calculateDisplacement
(
    const dictionary& coeffDict,
    const scalarField& minThickness,
    const PackedList<1>& isConvexEdgePoint,
    const PackedList<1>& isConcaveEdgePoint,
    List<snappyLayerDriver::extrudeMode>& extrudeStatus,
    pointField& patchDisp
)
{
    Info<< typeName << " : Smoothing using Medial Axis ..." << endl;

    const indirectPrimitivePatch& pp = adaptPatchPtr_;
    const labelList& meshPoints = pp.meshPoints();


    // Read settings
    // ~~~~~~~~~~~~~

    //- (lambda-mu) smoothing of internal displacement
    const label nSmoothDisplacement = coeffDict.lookupOrDefault<label>
    (
        "nSmoothDisplacement",
        0
    );

    //- Layer thickness too big
    const scalar maxThicknessToMedialRatio = coeffDict.lookupOrDefault<scalar>
    (
        "maxThicknessToMedialRatio",
        0.3
    );

    scalar minCosLayerTermination = 0.173;
    //- Feature angle when to stop adding layers
    if (coeffDict.found("featureAngleGrowAround"))
    {
        scalar angle = readScalar(coeffDict.lookup("featureAngleGrowAround"));
        minCosLayerTermination = Foam::cos
        (
            Foam::cos(degToRad(0.5*angle))
        );
    }

    //- Smoothing wanted patch thickness
    const label nSmoothPatchThickness = coeffDict.lookupOrDefault<label>
    (
        "nSmoothPatchThickness",
        10
    );

    //- Number of edges walking out
    const label nMedialAxisIter = coeffDict.lookupOrDefault<label>
    (
        "nMedialAxisIter",
        mesh().globalData().nTotalPoints()
    );

    //- Use strict extrusionIsland detection
    const Switch detectExtrusionIsland = coeffDict.lookupOrDefault<Switch>
    (
        "detectExtrusionIsland",
        false
    );


    // Precalulate master points/edge (only relevant for shared points/edges)
    const PackedBoolList isMeshMasterPoint(syncTools::getMasterPoints(mesh()));
    const PackedBoolList isMeshMasterEdge(syncTools::getMasterEdges(mesh()));

    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh().edges(),
            mesh().pointEdges()
        )
    );

    // Precalulate (patch) master point/edge
    const PackedBoolList isPatchMasterPoint
    (
        meshRefinement::getMasterPoints
        (
            mesh(),
            meshPoints
        )
    );

    const PackedBoolList isPatchMasterEdge
    (
        meshRefinement::getMasterEdges
        (
            mesh(),
            meshEdges
        )
    );

    // Normal component of normals of connected faces.
    // Edge normal for all pp. edges.
    vectorField patchEdgeNormals(pp.nEdges(), vector::zero);

    vectorField edgeNormal(mesh().nEdges(), vector(GREAT, GREAT, GREAT));

    const labelListList& edgeFaces = pp.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
         const labelList& eFaces = pp.edgeFaces()[edgeI];

         label meshEdgeI = meshEdges[edgeI];

         forAll(eFaces, i)
         {
             nomalsCombine()
             (
                 edgeNormal[meshEdgeI],
                 pp.faceNormals()[eFaces[i]]
             );
         }
    }

    syncTools::syncEdgeList
    (
        mesh(),
        edgeNormal,
        nomalsCombine(),
        vector(GREAT, GREAT, GREAT)          // null value
     );


    forAll(edgeFaces, edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        patchEdgeNormals[edgeI] = edgeNormal[meshEdgeI];
    }


    scalarField thickness(mag(patchDisp));

    forAll(thickness, patchPointI)
    {
        if (extrudeStatus[patchPointI] == snappyLayerDriver::NOEXTRUDE)
        {
            thickness[patchPointI] = 0.0;
        }
    }
    label numThicknessRatioExclude = 0;

    // reduce thickness where thickness/medial axis distance large
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<OBJstream> str;
    if (debug)
    {
        str.reset
        (
            new OBJstream
            (
                mesh().time().path()
              / "thicknessRatioExcludePoints_"
              + mesh().time().timeName()
              + ".obj"
            )
        );
        Info<< typeName
            << " : Writing points with too large an extrusion distance to "
            << str().name() << endl;
    }

    autoPtr<OBJstream> medialVecStr;
    if (debug)
    {
        medialVecStr.reset
        (
            new OBJstream
            (
                mesh().time().path()
              / "thicknessRatioExcludeMedialVec_"
              + mesh().time().timeName()
              + ".obj"
            )
        );
        Info<< typeName
            << " : Writing medial axis vectors on points with too large"
            << " an extrusion distance to " << medialVecStr().name() << endl;
    }

    //now check primitive edges for collisions on projection
    boolList marked(pp.localPoints().size(), false);
    scalarField reduceFactor(pp.localPoints().size(), GREAT);
    forAll(pp.edges(), patchEdgeI)
    {
        const edge& e = pp.edges()[patchEdgeI];
        label v0 = e[0];
        label v1 = e[1];

        vector n0 = patchDisp[v0] / (mag(patchDisp[v0]) + VSMALL);
        vector n1 = patchDisp[v1] / (mag(patchDisp[v1]) + VSMALL);

        vector eVec = e.vec(pp.localPoints());
        scalar eLen = mag(eVec);
        eVec /= (mag(eVec) + SMALL);
        scalar eLengthRatio = 0.5*eLen;

        scalar pDist0 = (thickness[v0] * n0) & eVec;
        scalar pDist1 = (thickness[v1] * n1) & eVec;

        scalar pDistTotal = 0;
        pDistTotal += pDist0;
        pDistTotal -= pDist1;

        if (pDistTotal > eLengthRatio)
        {
            scalar redRatio = eLengthRatio/pDistTotal;
            if (extrudeStatus[v0] != snappyLayerDriver::NOEXTRUDE)
            {
                reduceFactor[v0] = min(reduceFactor[v0],redRatio);
            }
            if (extrudeStatus[v1] != snappyLayerDriver::NOEXTRUDE)
            {
                reduceFactor[v1] = min(reduceFactor[v1],redRatio);
            }
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        reduceFactor,
        minEqOp<scalar>(),
        GREAT        // null value
    );

    forAll(reduceFactor, pointI)
    {
        if (reduceFactor[pointI] < GREAT)
        {
            vector n = patchDisp[pointI] / (mag(patchDisp[pointI]) + VSMALL);
            thickness[pointI] = reduceFactor[pointI]*thickness[pointI];
            patchDisp[pointI] = thickness[pointI]*n;
            numThicknessRatioExclude++;
            marked[pointI] = true;
        }
    }

    forAll(meshPoints, patchPointI)
    {
        if
        (
            !marked[patchPointI]
            && extrudeStatus[patchPointI] != snappyLayerDriver::NOEXTRUDE
        )
        {
            label pointI = meshPoints[patchPointI];

            //- Option 1: look only at extrusion thickness v.s. distance
            //  to nearest (medial axis or static) point.
            scalar mDist = medialDist_[pointI];
            scalar thicknessRatio = thickness[patchPointI]/(mDist+VSMALL);

            //- Option 2: Look at component in the direction
            //  of nearest (medial axis or static) point.
            vector n =
                patchDisp[patchPointI]
              / (mag(patchDisp[patchPointI]) + VSMALL);
            vector mVec = medialVec_[pointI]-mesh().points()[pointI];
            mVec /= mag(mVec)+VSMALL;
            thicknessRatio *= (n&mVec);

            if (thicknessRatio > maxThicknessToMedialRatio)
            {
                // Truncate thickness.
                if (debug)
                {
                    Pout<< "truncating displacement at "
                        << mesh().points()[pointI]
                        << " from " << thickness[patchPointI]
                        << " to "
                        <<  0.5
                           *(
                                minThickness[patchPointI]
                               +thickness[patchPointI]
                            )
                        << " medial direction:" << mVec
                        << " extrusion direction:" << n
                        << " with thicknessRatio:" << thicknessRatio
                        << endl;
                }

                thickness[patchPointI] =
                    0.5*(minThickness[patchPointI]+thickness[patchPointI]);

                patchDisp[patchPointI] = thickness[patchPointI]*n;

                numThicknessRatioExclude++;

                if (str.valid())
                {
                    const point& pt = mesh().points()[pointI];
                    str().write(linePointRef(pt, pt+patchDisp[patchPointI]));
                }
                if (medialVecStr.valid())
                {
                    const point& pt = mesh().points()[pointI];
                    medialVecStr().write
                    (
                        linePointRef
                        (
                            pt,
                            medialVec_[pointI]
                        )
                    );
                }
            }
        }
    }

    reduce(numThicknessRatioExclude, sumOp<label>());
    Info<< typeName << " : Reducing layer thickness at "
        << numThicknessRatioExclude
        << " nodes where thickness to medial axis distance is large " << endl;

    // find points where layer growth isolated to a lone point, edge or face
    findIsolatedRegions
    (
        minCosLayerTermination,
        detectExtrusionIsland,

        isPatchMasterEdge,
        isConvexEdgePoint,
        isConcaveEdgePoint,
        meshEdges,

        minThickness,
        patchEdgeNormals,

        extrudeStatus,
        patchDisp
    );

    // Update thickess for changed extrusion
    forAll(thickness, patchPointI)
    {
        if (extrudeStatus[patchPointI] == snappyLayerDriver::NOEXTRUDE)
        {
            thickness[patchPointI] = 0.0;
        }
    }

    syncTools::syncPointList
    (
        mesh(),
        meshPoints,
        thickness,
        minEqOp<scalar>(),  // combine op
        GREAT              // null value
    );

    // Smooth layer thickness on moving patch. Since some locations will have
    // disabled the extrusion this might cause big jumps in wanted displacement
    // for neighbouring patch points. So smooth the wanted displacement
    // before actually trying to move the mesh.
    fieldSmoother_.minSmoothField
    (
        nSmoothPatchThickness,
        isPatchMasterPoint,
        isPatchMasterEdge,
        meshEdges,
        pp,
        minThickness,
        thickness
    );

    // Dummy additional info for PointEdgeWave
    int dummyTrackData = 0;

    List<PointData<scalar>> pointWallDist(mesh().nPoints());

    const pointField& points = mesh().points();

    // 1. Calculate distance to points where displacement is specified.
    // This wave is used to transport layer thickness
    {
        // Dummy additional info for PointEdgeWave
        List<PointData<scalar>> edgeWallDist(mesh().nEdges());
        int dummyTrackData = 0;


        // Seed data.
        DynamicList<PointData<scalar>> wallInfo(meshPoints.size());
        DynamicList<label> wallPoints(meshPoints.size());

        forAll(meshPoints, patchPointI)
        {
            label meshPointI = meshPoints[patchPointI];
            wallPoints.append(meshPointI);
            wallInfo.append
            (
                PointData<scalar>
                (
                    points[meshPointI],
                    0.0,
                    thickness[patchPointI] // transport layer thickness
                )
            );
        }

        regionSplit cellRegion(mesh());
        if (cellRegion.nRegions() > 1)
        {
            Info<<"Multiregion mesh:"
                <<" Checking seeding for wave calculation"<<endl;

            const polyBoundaryMesh& patches = mesh().boundaryMesh();
            labelList numberSet(cellRegion.nRegions(),0);
            forAll(pp.addressing(), faceI)
            {
                label meshFaceI = pp.addressing()[faceI];
                label own = mesh().faceOwner()[meshFaceI];
                label regionI = cellRegion[own];
                numberSet[regionI]++;
            }
            Pstream::listCombineReduce(numberSet, plusOp<label>());

            forAll(patches, patchI)
            {
                const polyPatch& pp = patches[patchI];

                if (!isA<processorPolyPatch>(pp))
                {
                    forAll(pp, faceI)
                    {
                        label meshFaceI = pp.start() + faceI;
                        label own = mesh().faceOwner()[meshFaceI];
                        label regionI = cellRegion[own];
                        if (numberSet[regionI] == 0)
                        {
                            face f = mesh().faces()[meshFaceI];
                            forAll(f, fI)
                            {
                                label meshPointI = f[fI];
                                wallPoints.append(meshPointI);
                                wallInfo.append
                                (
                                    PointData<scalar>
                                    (
                                        points[meshPointI],
                                        0.0,
                                        0.0 // transport layer thickness
                                    )
                                );
                            }
                        }
                    }
                }
            }
        }
        wallInfo.shrink();
        wallPoints.shrink();

        // Do all calculations
        PointEdgeWave<PointData<scalar>> wallDistCalc
        (
            mesh(),
            wallPoints,
            wallInfo,
            pointWallDist,
            edgeWallDist,
            0,   // max iterations
            dummyTrackData
        );
        wallDistCalc.iterate(nMedialAxisIter);
    }


    // Calculate scaled displacement vector
    pointField& displacement = pointDisplacement_;

    forAll(displacement, pointI)
    {
        if (!pointWallDist[pointI].valid(dummyTrackData))
        {
            displacement[pointI] = Zero;
        }
        else
        {
            // 1) displacement on nearest wall point, scaled by medialRatio
            //    (wall distance / medial distance)
            // 2) pointWallDist[pointI].data() is layer thickness transported
            //    from closest wall point.
            // 3) shrink in opposite direction of addedPoints
            displacement[pointI] =
                -medialRatio_[pointI]
                *pointWallDist[pointI].data()
                *dispVec_[pointI];
        }
    }

    // Smear displacement away from fixed values (medialRatio=0 or 1)
    if (nSmoothDisplacement > 0)
    {
        PackedBoolList isToBeSmoothed(displacement.size(), false);

        forAll(displacement, i)
        {
            isToBeSmoothed[i] =
                medialRatio_[i] > SMALL && medialRatio_[i] < 1-SMALL;
        }

        fieldSmoother_.smoothLambdaMuDisplacement
        (
            nSmoothDisplacement,
            isMeshMasterPoint,
            isMeshMasterEdge,
            isToBeSmoothed,
            displacement
        );
    }
    else
    {
        const scalar edge0Len = readScalar
        (
            coeffDict.lookup("edge0Len")
        );

        // Smooth the displaced mesh
        smoothDisplacement
        (
            coeffDict,
            slipPatchIDs_,
            edge0Len,

            mesh(),
            meshMover_.curPoints()(),
            displacement
        );
    }
}


bool Foam::medialAxisMeshMover::shrinkMesh
(
    const dictionary& meshQualityDict,
    const label nAllowableErrors,
    labelList& checkFaces
)
{
    //- Number of attempts shrinking the mesh
    const label nSnap = meshQualityDict.lookupOrDefault<label>("nRelaxIter", 5);

    // Make sure displacement boundary conditions is uptodate with
    // internal field
    meshMover_.setDisplacementPatchFields();

    Info<< typeName << " : Moving mesh ..." << endl;
    scalar oldErrorReduction = -1;

    bool meshOk = false;

    for (label iter = 0; iter < 2*nSnap ; iter++)
    {
        Info<< typeName
            << " : Iteration " << iter << endl;
        if (iter == nSnap)
        {
            Info<< typeName
                << " : Displacement scaling for error reduction set to 0."
                << endl;
            oldErrorReduction = meshMover_.setErrorReduction(0.0);
        }

        if
        (
            meshMover_.scaleMesh
            (
                checkFaces,
                baffles_,
                meshMover_.paramDict(),
                meshQualityDict,
                true,
                nAllowableErrors
            )
        )
        {
            Info<< typeName << " : Successfully moved mesh" << endl;
            meshOk = true;
            break;
        }
    }

    if (oldErrorReduction >= 0)
    {
        meshMover_.setErrorReduction(oldErrorReduction);
    }

    Info<< typeName << " : Finished moving mesh ..." << endl;

    return meshOk;
}


bool Foam::medialAxisMeshMover::move
(
    const dictionary& moveDict,
    const label nAllowableErrors,
    const PackedList<1>& isConvexEdgePoint,
    const PackedList<1>& isConcaveEdgePoint,
    labelList& checkFaces
)
{
    // Read a few settings
    // ~~~~~~~~~~~~~~~~~~~

    //- Name of field specifying min thickness
    const word minThicknessName = word(moveDict.lookup("minThicknessName"));


    // Extract out patch-wise displacement
    const indirectPrimitivePatch& pp = adaptPatchPtr_();

    scalarField zeroMinThickness;
    if (minThicknessName == "none")
    {
        zeroMinThickness = scalarField(pp.nPoints(), 0.0);
    }
    const scalarField& minThickness =
    (
        (minThicknessName == "none")
      ? zeroMinThickness
      : mesh().lookupObject<scalarField>(minThicknessName)
    );

    pointField patchDisp
    (
        pointDisplacement_.internalField(),
        pp.meshPoints()
    );

    List<snappyLayerDriver::extrudeMode> extrudeStatus
    (
        pp.nPoints(),
        snappyLayerDriver::EXTRUDE
    );
    forAll(extrudeStatus, pointI)
    {
        if (mag(patchDisp[pointI]) <= minThickness[pointI]+SMALL)
        {
            extrudeStatus[pointI] = snappyLayerDriver::NOEXTRUDE;
        }
    }

    // Solve displacement
    calculateDisplacement
    (
        moveDict,
        minThickness,
        isConvexEdgePoint,
        isConcaveEdgePoint,
        extrudeStatus,
        patchDisp
    );

    //- Move mesh according to calculated displacement
    return shrinkMesh
    (
        moveDict,           // meshQualityDict,
        nAllowableErrors,   // nAllowableErrors
        checkFaces
    );
}


void Foam::medialAxisMeshMover::movePoints(const pointField& p)
{
    externalDisplacementMeshMover::movePoints(p);

    // Update motionSmoother for new geometry (moves adaptPatchPtr_)
    meshMover_.movePoints();

    // Assume corrent mesh location is correct (reset oldPoints, scale)
    meshMover_.correct();
}


// ************************************************************************* //
