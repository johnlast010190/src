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
    (c) 2012-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/polyMeshCheck/polyMeshTools.H"
#include "global/unitConversion/unitConversion.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polyMesh::checkFaceOrthogonality
(
    const vectorField& fAreas,
    const vectorField& cellCtrs,
    const bool report,
    const bool detailedReport,
    labelHashSet* setPtr,
    const scalar nonOrthThreshold
) const
{
    if (debug)
    {
        InfoInFunction << "Checking mesh non-orthogonality" << endl;
    }

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // Calculate orthogonality for all internal and coupled boundary faces
    // (1 for uncoupled boundary faces)
    tmp<scalarField> tortho = polyMeshTools::faceOrthogonality
    (
        *this,
        fAreas,
        cellCtrs
    );
    const scalarField& ortho = tortho.ref();

    scalar threshold
    (
        nonOrthThreshold > 0 ?
        nonOrthThreshold : primitiveMesh::nonOrthThreshold_
    );

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold =
        ::cos(degToRad(threshold));


    scalar minDDotS = GREAT;
    scalar sumDDotS = 0.0;
    label nSummed = 0;
    label severeNonOrth = 0;
    label errorNonOrth = 0;


    // Statistics only for internal and masters of coupled faces
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(*this));

    forAll(ortho, facei)
    {
        if (ortho[facei] < severeNonorthogonalityThreshold)
        {
            if (ortho[facei] > SMALL)
            {
                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                severeNonOrth++;
            }
            else
            {
                // Error : non-ortho too large
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
                if (detailedReport && errorNonOrth == 0)
                {
                    // Non-orthogonality greater than 90 deg
                    WarningInFunction
                        << "Severe non-orthogonality for face "
                        << facei
                        << " between cells " << own[facei]
                        << " and " << nei[facei]
                        << ": Angle = "
                        << radToDeg(::acos(min(1.0, max(-1.0, ortho[facei]))))
                        << " deg." << endl;
                }

                errorNonOrth++;
            }
        }

        if (isMasterFace[facei])
        {
            minDDotS = min(minDDotS, ortho[facei]);
            sumDDotS += radToDeg(::acos(min(1.0, max(-1.0, ortho[facei]))));
            nSummed++;
        }
    }

    reduce(minDDotS, minOp<scalar>());
    reduce(sumDDotS, sumOp<scalar>());
    reduce(nSummed, sumOp<label>());
    reduce(severeNonOrth, sumOp<label>());
    reduce(errorNonOrth, sumOp<label>());

    if (debug || report)
    {
        if (nSummed > 0)
        {
            if (debug || report)
            {
                Info<< "    Mesh non-orthogonality Max: "
                    << radToDeg(::acos(min(1.0, max(-1.0, minDDotS))))
                    << " average: "
                    << sumDDotS/nSummed
                    << endl;
            }
        }

        if (severeNonOrth > 0)
        {
            Info<< "   *Number of severely non-orthogonal (> "
                << threshold << " degrees) faces: "
                << severeNonOrth << "." << endl;
        }
    }

    if (errorNonOrth > 0)
    {
        if (debug || report)
        {
            Info<< " ***Number of non-orthogonality errors: "
                << errorNonOrth << "." << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Non-orthogonality check OK." << endl;
        }

        return false;
    }
}


bool Foam::polyMesh::checkFaceSkewness
(
    const pointField& points,
    const vectorField& fCtrs,
    const vectorField& fAreas,
    const vectorField& cellCtrs,
    const bool report,
    const bool detailedReport,
    labelHashSet* setPtr,
    const scalar skewnessThreshold
) const
{
    if (debug)
    {
        InfoInFunction << "Checking face skewness" << endl;
    }

    const labelList& own = faceOwner();
    const labelList& nei = faceNeighbour();

    // Warn if the skew correction vector is more than skewWarning times
    // larger than the face area vector

    tmp<scalarField> tskew = polyMeshTools::faceSkewness
    (
        *this,
        points,
        fCtrs,
        fAreas,
        cellCtrs
    );
    const scalarField& skew = tskew.ref();

    scalar maxSkew = max(skew);
    label nWarnSkew = 0;

    // Statistics only for all faces except slave coupled faces
    PackedBoolList isMasterFace(syncTools::getMasterFaces(*this));

    scalar threshold
    (
        skewnessThreshold > 0 ?
        skewnessThreshold : skewThreshold_
    );

    forAll(skew, facei)
    {
        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor mesh.
        if (skew[facei] > threshold)
        {
            if (setPtr)
            {
                setPtr->insert(facei);
            }
            if (detailedReport && nWarnSkew == 0)
            {
                // Non-orthogonality greater than 90 deg
                if (isInternalFace(facei))
                {
                    WarningInFunction
                        << "Severe skewness " << skew[facei]
                        << " for face " << facei
                        << " between cells " << own[facei]
                        << " and " << nei[facei];
                }
                else
                {
                    WarningInFunction
                        << "Severe skewness " << skew[facei]
                        << " for boundary face " << facei
                        << " on cell " << own[facei];
                }
            }

            if (isMasterFace[facei])
            {
                nWarnSkew++;
            }
        }
    }

    reduce(maxSkew, maxOp<scalar>());
    reduce(nWarnSkew, sumOp<label>());

    if (nWarnSkew > 0)
    {
        if (debug || report)
        {
            Info<< " *Max skewness = " << maxSkew
                << ", " << nWarnSkew << " highly skew faces detected"
                   " which may impair the quality of the results"
                << endl;
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    Max skewness = " << maxSkew << " OK." << endl;
        }

        return false;
    }
}


bool Foam::polyMesh::checkEdgeAlignment
(
    const pointField& p,
    const bool report,
    const Vector<label>& directions,
    labelHashSet* setPtr
) const
{
    // Check 1D/2Dness of edges. Gets passed the non-empty directions and
    // checks all edges in the mesh whether they:
    // - have no component in a non-empty direction or
    // - are only in a singe non-empty direction.
    // Empty direction info is passed in as a vector of labels (synchronised)
    // which are 1 if the direction is non-empty, 0 if it is.

    if (debug)
    {
        InfoInFunction << "Checking edge alignment" << endl;
    }

    label nDirs = 0;
    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        if (directions[cmpt] == 1)
        {
            nDirs++;
        }
        else if (directions[cmpt] != 0)
        {
            FatalErrorInFunction
                << "directions should contain 0 or 1 but is now " << directions
                << exit(FatalError);
        }
    }

    if (nDirs == vector::nComponents)
    {
        return false;
    }


    const faceList& fcs = faces();

    EdgeMap<label> edgesInError;

    forAll(fcs, facei)
    {
        const face& f = fcs[facei];

        forAll(f, fp)
        {
            label p0 = f[fp];
            label p1 = f.nextLabel(fp);
            if (p0 < p1)
            {
                vector d(p[p1]-p[p0]);
                scalar magD = mag(d);

                if (magD > ROOTVSMALL)
                {
                    d /= magD;

                    // Check how many empty directions are used by the edge.
                    label nEmptyDirs = 0;
                    label nNonEmptyDirs = 0;
                    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                    {
                        if (mag(d[cmpt]) > 1e-6)
                        {
                            if (directions[cmpt] == 0)
                            {
                                nEmptyDirs++;
                            }
                            else
                            {
                                nNonEmptyDirs++;
                            }
                        }
                    }

                    if (nEmptyDirs == 0)
                    {
                        // Purely in ok directions.
                    }
                    else if (nEmptyDirs == 1)
                    {
                        // Ok if purely in empty directions.
                        if (nNonEmptyDirs > 0)
                        {
                            edgesInError.insert(edge(p0, p1), facei);
                        }
                    }
                    else if (nEmptyDirs > 1)
                    {
                        // Always an error
                        edgesInError.insert(edge(p0, p1), facei);
                    }
                }
            }
        }
    }

    label nErrorEdges = returnReduce(edgesInError.size(), sumOp<label>());

    if (nErrorEdges > 0)
    {
        if (debug || report)
        {
            Info<< " ***Number of edges not aligned with or perpendicular to "
                << "non-empty directions: " << nErrorEdges << endl;
        }

        if (setPtr)
        {
            setPtr->resize(2*edgesInError.size());
            forAllConstIter(EdgeMap<label>, edgesInError, iter)
            {
                setPtr->insert(iter.key()[0]);
                setPtr->insert(iter.key()[1]);
            }
        }

        return true;
    }
    else
    {
        if (debug || report)
        {
            Info<< "    All edges aligned with or perpendicular to "
                << "non-empty directions." << endl;
        }
        return false;
    }
}


bool Foam::polyMesh::checkCellDeterminant
(
    const vectorField& faceAreas,
    const bool report,
    labelHashSet* setPtr,
    const Vector<label>& meshD
) const
{
    const scalar warnDet = 1e-3;

    if (debug)
    {
        InfoInFunction << "Checking for under-determined cells" << endl;
    }

    tmp<scalarField> tcellDeterminant = primitiveMeshTools::cellDeterminant
    (
        *this,
        meshD,
        faceAreas,
        syncTools::getInternalOrCoupledFaces(*this)
    );
    scalarField& cellDeterminant = tcellDeterminant.ref();


    label nErrorCells = 0;
    scalar minDet = min(cellDeterminant);
    scalar sumDet = sum(cellDeterminant);

    forAll(cellDeterminant, celli)
    {
        if (cellDeterminant[celli] < warnDet)
        {
            if (setPtr)
            {
                setPtr->insert(celli);
            }

            nErrorCells++;
        }
    }

    label nSummed = cellDeterminant.size();
    reduce(
        std::tie(nErrorCells, minDet, sumDet, nSummed),
        ParallelOp<sumOp<label>, minOp<scalar>, sumOp<scalar>, sumOp<label>>{}
    );

    if (debug || report)
    {
        if (nSummed > 0)
        {
            Info<< "    Cell determinant (wellposedness) : minimum: " << minDet
                << " average: " << sumDet/nSummed
                << endl;
        }
    }

    if (nErrorCells > 0)
    {
        if (debug || report)
        {
            Info<< " ***Cells with small determinant (< "
                << warnDet << ") found, number of cells: "
                << nErrorCells << endl;
        }

        return true;
    }

    if (debug || report)
    {
        Info<< "    Cell determinant check OK." << endl;
    }

    return false;
}


bool Foam::polyMesh::checkFaceWeight
(
    const vectorField& fCtrs,
    const vectorField& fAreas,
    const vectorField& cellCtrs,
    const bool report,
    const scalar minWeight,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        InfoInFunction << "Checking for low face interpolation weights" << endl;
    }

    tmp<scalarField> tfaceWght = polyMeshTools::faceWeights
    (
        *this,
        fCtrs,
        fAreas,
        cellCtrs
    );
    scalarField& faceWght = tfaceWght.ref();


    label nErrorFaces = 0;
    scalar minDet = GREAT;
    scalar sumDet = 0.0;
    label nSummed = 0;

    // Statistics only for internal and masters of coupled faces
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(*this));

    forAll(faceWght, facei)
    {
        if (faceWght[facei] < minWeight)
        {
            // Note: insert both sides of coupled faces
            if (setPtr)
            {
                setPtr->insert(facei);
            }

            nErrorFaces++;
        }

        // Note: statistics only on master of coupled faces
        if (isMasterFace[facei])
        {
            minDet = min(minDet, faceWght[facei]);
            sumDet += faceWght[facei];
            nSummed++;
        }
    }

    reduce(
        std::tie(nErrorFaces, minDet, sumDet, nSummed),
        ParallelOp<sumOp<label>, minOp<scalar>, sumOp<scalar>, sumOp<label>>{}
    );

    if (debug || report)
    {
        if (nSummed > 0)
        {
            Info<< "    Face interpolation weight : minimum: " << minDet
                << " average: " << sumDet/nSummed
                << endl;
        }
    }

    if (nErrorFaces > 0)
    {
        if (debug || report)
        {
            Info<< " ***Faces with small interpolation weight (< " << minWeight
                << ") found, number of faces: "
                << nErrorFaces << endl;
        }

        return true;
    }

    if (debug || report)
    {
        Info<< "    Face interpolation weight check OK." << endl;
    }

    return false;
}


bool Foam::polyMesh::checkVolRatio
(
    const scalarField& cellVols,
    const bool report,
    const scalar minRatio,
    labelHashSet* setPtr
) const
{
    if (debug)
    {
        InfoInFunction << "Checking for volume ratio < " << minRatio << endl;
    }

    tmp<scalarField> tvolRatio = polyMeshTools::volRatio(*this, cellVols);
    scalarField& volRatio = tvolRatio.ref();


    label nErrorFaces = 0;
    scalar minDet = GREAT;
    scalar sumDet = 0.0;
    label nSummed = 0;

    // Statistics only for internal and masters of coupled faces
    PackedBoolList isMasterFace(syncTools::getInternalOrMasterFaces(*this));

    forAll(volRatio, facei)
    {
        if (volRatio[facei] < minRatio)
        {
            // Note: insert both sides of coupled faces
            if (setPtr)
            {
                setPtr->insert(facei);
            }

            nErrorFaces++;
        }

        // Note: statistics only on master of coupled faces
        if (isMasterFace[facei])
        {
            minDet = min(minDet, volRatio[facei]);
            sumDet += volRatio[facei];
            nSummed++;
        }
    }

    reduce(
        std::tie(nErrorFaces, minDet, sumDet, nSummed),
        ParallelOp<sumOp<label>, minOp<scalar>, sumOp<scalar>, sumOp<label>>{}
    );

    if (debug || report)
    {
        if (nSummed > 0)
        {
            Info<< "    Face volume ratio : minimum: " << minDet
                << " average: " << sumDet/nSummed
                << endl;
        }
    }

    if (nErrorFaces > 0)
    {
        if (debug || report)
        {
            Info<< " ***Faces with small volume ratio (< " << minRatio
                << ") found, number of faces: "
                << nErrorFaces << endl;
        }

        return true;
    }

    if (debug || report)
    {
        Info<< "    Face volume ratio check OK." << endl;
    }

    return false;
}


bool Foam::polyMesh::checkFaceOrthogonality
(
    const bool report,
    labelHashSet* setPtr,
    const scalar nonOrthThreshold
) const
{
    return checkFaceOrthogonality
    (
        faceAreas(),
        cellCentres(),
        report,
        false,  // detailedReport
        setPtr,
        nonOrthThreshold
    );
}


bool Foam::polyMesh::checkFaceSkewness
(
    const bool report,
    labelHashSet* setPtr,
    const scalar skewnessThreshold
) const
{
    return checkFaceSkewness
    (
        points(),
        faceCentres(),
        faceAreas(),
        cellCentres(),
        report,
        false,  // detailedReport
        setPtr,
        skewnessThreshold
    );
}


bool Foam::polyMesh::checkEdgeAlignment
(
    const bool report,
    const Vector<label>& directions,
    labelHashSet* setPtr
) const
{
    return checkEdgeAlignment
    (
        points(),
        report,
        directions,
        setPtr
    );
}


bool Foam::polyMesh::checkCellDeterminant
(
    const bool report,
    labelHashSet* setPtr
) const
{
    return checkCellDeterminant
    (
        faceAreas(),
        report,
        setPtr,
        geometricD()
    );
}


bool Foam::polyMesh::checkFaceWeight
(
    const bool report,
    const scalar minWeight,
    labelHashSet* setPtr
) const
{
    return checkFaceWeight
    (
        faceCentres(),
        faceAreas(),
        cellCentres(),
        report,
        minWeight,
        setPtr
    );
}


bool Foam::polyMesh::checkVolRatio
(
    const bool report,
    const scalar minRatio,
    labelHashSet* setPtr
) const
{
    return checkVolRatio(cellVolumes(), report, minRatio, setPtr);
}


bool Foam::polyMesh::checkMeshMotion
(
    const pointField& newPoints,
    const bool report,
    const bool detailedReport
) const
{
    if (debug || report)
    {
        Pout<< "bool polyMesh::checkMeshMotion("
            << "const pointField&, const bool, const bool) const: "
            << "checking mesh motion" << endl;
    }

    vectorField fCtrs(nFaces());
    vectorField fAreas(nFaces());
    scalarField magfAreas(nFaces());

    makeFaceCentresAndAreas(newPoints, fCtrs, fAreas, magfAreas);

    // Check cell volumes and calculate new cell centres
    vectorField cellCtrs(nCells());
    scalarField cellVols(nCells());

    makeCellCentresAndVols(fCtrs, fAreas, cellCtrs, cellVols);

    // Check cell volumes
    bool error = checkCellVolumes
    (
        cellVols,       // vols
        report,         // report
        detailedReport, // detailedReport
        nullptr            // setPtr
    );


    // Check face areas
    bool areaError = checkFaceAreas
    (
        fAreas,
        report,         // report
        detailedReport, // detailedReport,
        nullptr            // setPtr
    );
    error = error || areaError;


    // Check pyramid volumes
    bool pyrVolError = checkFacePyramids
    (
        cellCtrs,
        fCtrs,
        fAreas,
        report,         // report,
        detailedReport, // detailedReport,
        -SMALL,         // minPyrVol
        nullptr            // setPtr
    );
    error = error || pyrVolError;


    // Check face non-orthogonality
    bool nonOrthoError = checkFaceOrthogonality
    (
        fAreas,
        cellCtrs,
        report,         // report
        detailedReport, // detailedReport
        nullptr            // setPtr
    );
    error = error || nonOrthoError;


    if (!error && (debug || report))
    {
        Pout<< "Mesh motion check OK." << endl;
    }

    return error;
}


void Foam::polyMesh::checkZones() const
{
    //Calculate zone statistics
    if (cellZones().size() || faceZones().size() || pointZones().size())
    {
        Info<< endl;
        if (cellZones().size())
        {
            //Calculate cellZone statistics
            Info<< "Checking cellZones..." << nl
                << nl;

            label maxZoneNameLen = 0;
            forAll(cellZones(), zoneID)
            {
                word zoneName = cellZones().names()[zoneID];
                maxZoneNameLen = max(maxZoneNameLen, label(zoneName.size()));
            }

            Info<< setf(ios_base::left) << setw(maxZoneNameLen) << "cellZone"
                 << setw(0) << "   Cells  Volume   Centroid" << endl;

            forAll(cellZones(), zoneID)
            {
                const labelList& zoneCells = cellZones()[zoneID];

                scalar zoneVolume = 0.0;
                vector zoneCentre = vector::zero;
                forAll(zoneCells, i)
                {
                    label cellI = zoneCells[i];
                    zoneVolume += cellVolumes()[cellI];
                    zoneCentre += cellCentres()[cellI]*cellVolumes()[cellI];
                }

                label numZoneCells = zoneCells.size();
                reduce(
                    std::tie(zoneVolume, zoneCentre, numZoneCells),
                    ParallelOp<sumOp<scalar>, sumOp<vector>, sumOp<label>>{}
                );
                if (zoneVolume > SMALL)
                {
                    zoneCentre /= zoneVolume;
                }

                Info<< setf(ios_base::left) << setw(maxZoneNameLen)
                    << cellZones().names()[zoneID]
                    << " " << setw(8)
                    << numZoneCells
                    << " " << setw(6) << zoneVolume
                    << " " << setw(0) << zoneCentre
                    << endl;
            }
            Info<<nl;
        }

        if (faceZones().size())
        {
            //Calculate faceZone statistics
            Info<< "Checking faceZones..." << nl
                << nl;

            label maxZoneNameLen = 0;
            forAll(faceZones(), zoneID)
            {
                word zoneName = faceZones().names()[zoneID];
                maxZoneNameLen = max(maxZoneNameLen, label(zoneName.size()));
            }

            Info<< setf(ios_base::left) << setw(maxZoneNameLen) << "faceZone"
                 << setw(0) << "   Faces   Area   Centroid" << endl;

            forAll(faceZones(), zoneID)
            {
                const labelList& zoneFaces = faceZones()[zoneID];

                scalar zoneArea = 0.0;
                vector zoneCentre = vector::zero;
                forAll(zoneFaces, i)
                {
                    label faceI = zoneFaces[i];
                    scalar faceArea = magFaceAreas()[faceI];
                    zoneArea += mag(faceArea);
                    zoneCentre += faceCentres()[faceI]*faceArea;
                }

                label numZoneFaces = zoneFaces.size();
                reduce(
                    std::tie(zoneArea, zoneCentre, numZoneFaces),
                    ParallelOp<sumOp<scalar>, sumOp<vector>, sumOp<label>>{}
                );

                if (zoneArea > SMALL)
                {
                    zoneCentre /= zoneArea;
                }

                Info<< setf(ios_base::left) << setw(maxZoneNameLen)
                    << faceZones().names()[zoneID]
                    << " " << setw(8)
                    << numZoneFaces
                    << " " << setw(6) << zoneArea
                    << " " << setw(0) << zoneCentre
                    << endl;
            }
            Info<<nl;
        }

        if (pointZones().size())
        {
            //Calculate pointZone statistics
            Info<< "Checking pointZones..." << nl
                << nl;

            label maxZoneNameLen = 0;
            forAll(pointZones(), zoneID)
            {
                word zoneName = pointZones().names()[zoneID];
                maxZoneNameLen = max(maxZoneNameLen, label(zoneName.size()));
            }

            Info<< setf(ios_base::left) << setw(maxZoneNameLen) << "pointZone"
                 << setw(0) << " Points" << endl;

            forAll(pointZones(), zoneID)
            {
                const labelList& zonePts = pointZones()[zoneID];

                Info<< setf(ios_base::left) << setw(maxZoneNameLen)
                    << pointZones().names()[zoneID]
                    << " " << setw(8)
                    << returnReduce(zonePts.size(), sumOp<label>())
                    << endl;
            }
            Info<<nl;
        }
    }
}

// ************************************************************************* //
