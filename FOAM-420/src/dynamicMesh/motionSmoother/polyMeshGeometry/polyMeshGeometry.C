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
    (c) 2010-2012 Esi Ltd.
    (c) 2015-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSmoother/polyMeshGeometry/polyMeshGeometry.H"
#include "meshes/polyMesh/polyMeshTetDecomposition/polyMeshTetDecomposition.H"
#include "meshes/primitiveShapes/tetrahedron/tetrahedron.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "global/unitConversion/unitConversion.H"
#include "meshes/primitiveMesh/primitiveMeshCheck/primitiveMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyMeshGeometry, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


bool Foam::polyMeshGeometry::isPyrCell(const polyMesh& mesh, const label cellI)
{
    const labelList& cFaces = mesh.cells()[cellI];

    if (cFaces.size() == 5)
    {
        label nQuad = 0;
        label nTri = 0;

        forAll(cFaces, cFI)
        {
            label facei = cFaces[cFI];

            const labelList& fEdges = mesh.faceEdges()[facei];
            label nEdges = 0;
            forAll(fEdges, fEI)
            {
                label edgeI = fEdges[fEI];
                edge e = mesh.edges()[edgeI];
                if (e.mag(mesh.points()) > SMALL)
                {
                    nEdges++;
                }
            }
            if (nEdges == 4)
            {
                nQuad++;
            }
            else if (nEdges == 3)
            {
                nTri++;
            }
        }

        if (nQuad == 1 && nTri == 4)
        {
            return true;
        }
    }

    return false;
}


void Foam::polyMeshGeometry::updateFaceCentresAndAreas
(
    const pointField& p,
    const labelList& changedFaces
)
{
    const faceList& fs = mesh_.faces();

    forAll(changedFaces, i)
    {
        label facei = changedFaces[i];

        const labelList& f = fs[facei];
        label nPoints = f.size();

        // If the face is a triangle, do a direct calculation for efficiency
        // and to avoid round-off error-related problems
        if (nPoints == 3)
        {
            faceCentres_[facei] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
            faceAreas_[facei] = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
        }
        else
        {
            vector sumN = Zero;
            scalar sumA = 0.0;
            vector sumAc = Zero;

            point fCentre = p[f[0]];
            for (label pi = 1; pi < nPoints; pi++)
            {
                fCentre += p[f[pi]];
            }

            fCentre /= nPoints;

            for (label pi = 0; pi < nPoints; pi++)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];

                vector c = p[f[pi]] + nextPoint + fCentre;
                vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
                scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            faceCentres_[facei] = (1.0/3.0)*sumAc/(sumA + VSMALL);
            faceAreas_[facei] = 0.5*sumN;
        }
    }
}


void Foam::polyMeshGeometry::updateCellCentresAndVols
(
    const labelList& changedCells,
    const labelList& changedFaces
)
{
    const labelList& own = mesh().faceOwner();
    const cellList& cells = mesh().cells();

    // Clear the fields for accumulation
    UIndirectList<vector>(cellCentres_, changedCells) = Zero;
    UIndirectList<scalar>(cellVolumes_, changedCells) = 0.0;


    // Re-calculate the changed cell centres and volumes
    forAll(changedCells, changedCellI)
    {
        const label cellI(changedCells[changedCellI]);
        const labelList& cFaces(cells[cellI]);

        // Estimate the cell centre and bounding box using the face centres
        vector cEst(Zero);
        boundBox bb(boundBox::invertedBox);

        forAll(cFaces, cFaceI)
        {
            const point& fc = faceCentres_[cFaces[cFaceI]];
            cEst += fc;
            bb.add(fc);
        }
        cEst /= cFaces.size();


        // Sum up the face-pyramid contributions
        forAll(cFaces, cFaceI)
        {
            const label facei(cFaces[cFaceI]);

            // Calculate 3* the face-pyramid volume
            scalar pyr3Vol = faceAreas_[facei] & (faceCentres_[facei] - cEst);

            if (own[facei] != cellI)
            {
                pyr3Vol = -pyr3Vol;
            }

            // Accumulate face-pyramid volume
            cellVolumes_[cellI] += pyr3Vol;

            // Calculate the face-pyramid centre
            const vector pCtr = (3.0/4.0)*faceCentres_[facei] + (1.0/4.0)*cEst;

            // Accumulate volume-weighted face-pyramid centre
            cellCentres_[cellI] += pyr3Vol*pCtr;
        }

        // Average the accumulated quantities

        if (mag(cellVolumes_[cellI]) > VSMALL)
        {
            point cc = cellCentres_[cellI] / cellVolumes_[cellI];

            // Do additional check for collapsed cells since some volumes
            // (e.g. 1e-33) do not trigger above but do return completely
            // wrong cell centre
            if (bb.contains(cc))
            {
                cellCentres_[cellI] = cc;
            }
            else
            {
                cellCentres_[cellI] = cEst;
            }
        }
        else
        {
            cellCentres_[cellI] = cEst;
        }

        cellVolumes_[cellI] *= (1.0/3.0);
    }
}


Foam::labelList Foam::polyMeshGeometry::affectedCells
(
    const polyMesh& mesh,
    const labelList& changedFaces
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    labelHashSet affectedCells(2*changedFaces.size());

    forAll(changedFaces, i)
    {
        label facei = changedFaces[i];

        affectedCells.insert(own[facei]);

        if (mesh.isInternalFace(facei))
        {
            affectedCells.insert(nei[facei]);
        }
    }
    return affectedCells.toc();
}


Foam::scalar Foam::polyMeshGeometry::checkNonOrtho
(
    const polyMesh& mesh,
    const bool verbose,
    const scalar severeNonorthogonalityThreshold,
    const label facei,
    const vector& s,    // face area vector
    const vector& d,    // cc-cc vector

    label& severeNonOrth,
    label& errorNonOrth,
    labelHashSet* setPtr
)
{
    scalar dDotS = (d & s)/(mag(d)*mag(s) + VSMALL);

    if (dDotS < severeNonorthogonalityThreshold)
    {
        label nei = -1;

        if (mesh.isInternalFace(facei))
        {
            nei = mesh.faceNeighbour()[facei];
        }

        if (dDotS > SMALL)
        {
            if (verbose)
            {
                // Severe non-orthogonality but mesh still OK
                Pout<< "Severe non-orthogonality for face " << facei
                    << " between cells " << mesh.faceOwner()[facei]
                    << " and " << nei
                    << ": Angle = "
                    << radToDeg(::acos(dDotS))
                    << " deg." << endl;
            }

            severeNonOrth++;
        }
        else
        {
            // Non-orthogonality greater than 90 deg
            if (verbose)
            {
                WarningInFunction
                    << "Severe non-orthogonality detected for face "
                    << facei
                    << " between cells " << mesh.faceOwner()[facei]
                    << " and " << nei
                    << ": Angle = "
                    << radToDeg(::acos(dDotS))
                    << " deg." << endl;
            }

            errorNonOrth++;
        }

        if (setPtr)
        {
            setPtr->insert(facei);
        }
    }
    return dDotS;
}


Foam::scalar Foam::polyMeshGeometry::calcSkewness
(
    const polyMesh& mesh,
    const point& ownCc,
    const point& neiCc,
    const point& fc,
    const vector& fa,
    const label& facei
)
{
    scalar dOwn = mag(fc - ownCc);
    scalar dNei = mag(fc - neiCc);

    point faceIntersection =
        ownCc*dNei/(dOwn+dNei+VSMALL)
      + neiCc*dOwn/(dOwn+dNei+VSMALL);

    scalar fSkew = mag(fc - faceIntersection)/
        (
            max(mag(neiCc-ownCc), sqrt(mag(fa)))
            + VSMALL
         );

    if (fSkew < 0.1)
    {
        return fSkew;
    }
    else
    {
         vector fN = fa / (mag(fa) + VSMALL);
        vector dir = faceIntersection - fc;

        faceIntersection -= (dir & fN) * fN;

        dir = fc - faceIntersection;
        dir /= (mag(dir) + SMALL);

        const face& f = mesh.faces()[facei];

        scalar maxVal = -GREAT;
        scalar minVal = GREAT;

        forAll(f, fI)
        {
            point pt = mesh.points()[f[fI]];
            scalar dp = pt & dir;
            maxVal = max(maxVal, dp);
            minVal = min(minVal, dp);
        }
        scalar charDist = maxVal - minVal;

        fSkew = max(mag(fc - faceIntersection)/(charDist + VSMALL),fSkew);

        return fSkew;
    }
}

// Create the neighbour pyramid - it will have positive volume
bool Foam::polyMeshGeometry::checkFaceTet
(
    const polyMesh& mesh,
    const bool verbose,
    const scalar minTetQuality,
    const pointField& p,
    const label facei,
    const point& fc,    // face centre
    const point& cc,    // cell centre

    labelHashSet* setPtr
)
{
    const face& f = mesh.faces()[facei];

    forAll(f, fp)
    {
        scalar tetQual = tetPointRef
        (
            p[f[fp]],
            p[f.nextLabel(fp)],
            fc,
            cc
        ).quality();

        if (tetQual < minTetQuality)
        {
            if (verbose)
            {
                Pout<< "bool polyMeshGeometry::checkFaceTets("
                    << "const bool, const scalar, const pointField&"
                    << ", const pointField&"
                    << ", const labelList&, labelHashSet*) : "
                    << "face " << facei
                    << " has a triangle that points the wrong way."
                     << endl
                    << "Tet quality: " << tetQual
                    << " Face " << facei
                    << endl;
            }

            if (setPtr)
            {
                setPtr->insert(facei);
            }
            return true;
        }
    }
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::polyMeshGeometry::polyMeshGeometry(const polyMesh& mesh)
:
    mesh_(mesh)
{
    correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMeshGeometry::correct()
{
    faceAreas_ = mesh_.faceAreas();
    faceCentres_ = mesh_.faceCentres();
    cellCentres_ = mesh_.cellCentres();
    cellVolumes_ = mesh_.cellVolumes();
}


void Foam::polyMeshGeometry::correct
(
    const pointField& p,
    const labelList& changedFaces
)
{
    // Update face quantities
    updateFaceCentresAndAreas(p, changedFaces);
    // Update cell quantities from face quantities
    updateCellCentresAndVols(affectedCells(mesh_, changedFaces), changedFaces);
}


Foam::label Foam::polyMeshGeometry::findStatIndex
(
    const PtrList<meshStatistics>* meshStats,
    const word& metricName
)
{
    if (meshStats != nullptr)
    {
        const PtrList<meshStatistics>& ms = *meshStats;
        forAll(ms, i)
        {
            if (ms[i].name() == metricName)
            {
                return i;
            }
        }
    }
    return -1;
}


bool Foam::polyMeshGeometry::checkFaceDotProduct
(
    const bool verbose,
    const scalar orthWarn,
    const scalar maxFaceCentreNonOrtho,
    const polyMesh& mesh,
    const vectorField& faceCentres,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr,
    PtrList<meshStatistics>* meshStats
)
{
    // for all internal and coupled faces check theat the d dot S product
    // is positive
    const word metric = "nonOrthogonality";

    const label fieldIndex = findStatIndex(meshStats, metric);

    DynamicList<scalar> qmetrics(checkFaces.size());
    autoPtr<PackedBoolList> isMasterFacePtr;
    if (fieldIndex != -1)
    {
        isMasterFacePtr.reset
        (
            new PackedBoolList(syncTools::getMasterFaces(mesh))
        );
        (*meshStats)[fieldIndex].field() = -GREAT;
    }

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Severe nonorthogonality threshold
    const scalar severeNonorthogonalityThreshold = ::cos(degToRad(orthWarn));
    const scalar severeFaceCentreOrthogonalityThreshold =
        ::cos(degToRad(maxFaceCentreNonOrtho));

    bool additionalOrthoCheck =
        (maxFaceCentreNonOrtho < 180.0-SMALL ? true : false );

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces() - mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    scalar minDDotS = GREAT;

    scalar sumDDotS = 0;
    label nDDotS = 0;

    label severeNonOrth = 0;

    label errorNonOrth = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const point& ownCc = cellCentres[own[facei]];

        if (mesh.isInternalFace(facei))
        {
            scalar dDotS = checkNonOrtho
            (
                mesh,
                verbose,
                severeNonorthogonalityThreshold,
                facei,
                faceAreas[facei],
                cellCentres[nei[facei]] - ownCc,

                severeNonOrth,
                errorNonOrth,
                setPtr
            );

            if (additionalOrthoCheck)
            {
                dDotS =
                    min
                    (
                        dDotS,
                        checkNonOrtho
                        (
                            mesh,
                            verbose,
                            severeFaceCentreOrthogonalityThreshold,
                            facei,
                            faceAreas[facei],
                            faceCentres[facei] - ownCc,

                            severeNonOrth,
                            errorNonOrth,
                            setPtr
                         )
                     );

                dDotS =
                    min
                    (
                        dDotS,
                        checkNonOrtho
                        (
                            mesh,
                            verbose,
                            severeFaceCentreOrthogonalityThreshold,
                            facei,
                            faceAreas[facei],
                            cellCentres[nei[facei]] - faceCentres[facei],

                            severeNonOrth,
                            errorNonOrth,
                            setPtr
                         )
                     );

            }

            if (dDotS < minDDotS)
            {
                minDDotS = dDotS;
            }
            scalar orthoAngle =
                radToDeg(::acos(min(scalar(1.0), max(scalar(-1.0), dDotS))));
            sumDDotS += orthoAngle;
            nDDotS++;
            if (fieldIndex != -1)
            {
                qmetrics.append(orthoAngle);
                (*meshStats)[fieldIndex].field()[own[facei]] =
                    max
                    (
                        (*meshStats)[fieldIndex].field()[own[facei]],
                        orthoAngle
                    );
                (*meshStats)[fieldIndex].field()[nei[facei]] =
                    max
                    (
                        (*meshStats)[fieldIndex].field()[nei[facei]],
                        orthoAngle
                    );
            }
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                scalar dDotS = checkNonOrtho
                (
                    mesh,
                    verbose,
                    severeNonorthogonalityThreshold,
                    facei,
                    faceAreas[facei],
                    neiCc[facei-mesh.nInternalFaces()] - ownCc,

                    severeNonOrth,
                    errorNonOrth,
                    setPtr
                );

                if (additionalOrthoCheck)
                {
                    dDotS = min
                    (
                        dDotS,
                        checkNonOrtho
                        (
                            mesh,
                            verbose,
                            severeFaceCentreOrthogonalityThreshold,
                            facei,
                            faceAreas[facei],
                            faceCentres[facei] - ownCc,
                            severeNonOrth,
                            errorNonOrth,
                            setPtr
                         )
                     );
                }

                if (dDotS < minDDotS)
                {
                    minDDotS = dDotS;
                }

                scalar orthoAngle =
                    radToDeg(::acos(min(scalar(1.0),max(scalar(-1.0), dDotS))));
                sumDDotS += orthoAngle;
                nDDotS++;
                if (fieldIndex != -1)
                {
                    if (isMasterFacePtr().get(facei))
                    {
                        qmetrics.append(orthoAngle);
                    }
                    (*meshStats)[fieldIndex].field()[own[facei]] =
                        max
                        (
                            (*meshStats)[fieldIndex].field()[own[facei]],
                            orthoAngle
                        );
                }
            }
            else if (additionalOrthoCheck)
            {
                scalar dDotS = checkNonOrtho
                (
                    mesh,
                    verbose,
                    severeFaceCentreOrthogonalityThreshold,
                    facei,
                    faceAreas[facei],
                    faceCentres[facei] - ownCc,
                    severeNonOrth,
                    errorNonOrth,
                    setPtr
                );

                if (dDotS < minDDotS)
                {
                    minDDotS = dDotS;
                }
                scalar orthoAngle =
                    radToDeg(::acos(min(scalar(1.0),max(scalar(-1.0), dDotS))));
                sumDDotS += orthoAngle;
                nDDotS++;
                if (fieldIndex != -1)
                {
                    qmetrics.append(orthoAngle);
                    (*meshStats)[fieldIndex].field()[own[facei]] =
                        max
                        (
                            (*meshStats)[fieldIndex].field()[own[facei]],
                            orthoAngle
                        );
                }
            }
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];

        scalar dDotS = checkNonOrtho
        (
            mesh,
            verbose,
            severeNonorthogonalityThreshold,
            face0,
            faceAreas[face0],
            cellCentres[own[face1]] - ownCc,

            severeNonOrth,
            errorNonOrth,
                setPtr
         );

        if (additionalOrthoCheck)
        {
            dDotS = min
            (
                dDotS,
                checkNonOrtho
                (
                    mesh,
                    verbose,
                    severeFaceCentreOrthogonalityThreshold,
                    face0,
                    faceAreas[face0],
                    cellCentres[own[face1]] - faceCentres[face0],

                    severeNonOrth,
                    errorNonOrth,
                    setPtr
                )
            );

            dDotS = min
            (
                dDotS,
                checkNonOrtho
                (
                    mesh,
                    verbose,
                    severeFaceCentreOrthogonalityThreshold,
                    face0,
                    faceAreas[face0],
                    faceCentres[face0] - ownCc,

                    severeNonOrth,
                    errorNonOrth,
                    setPtr
                )
            );
        }


        if (setPtr && dDotS < severeNonorthogonalityThreshold)
        {
            setPtr->insert(face1);
        }

        if (dDotS < minDDotS)
        {
            minDDotS = dDotS;
        }
        scalar orthoAngle =
            radToDeg(::acos(min(scalar(1.0),max(scalar(-1.0), dDotS))));
        sumDDotS += orthoAngle;
        nDDotS++;
    }

    reduce(
        std::tie(minDDotS, sumDDotS, nDDotS, severeNonOrth, errorNonOrth),
        ParallelOp<minOp<scalar>, sumOp<scalar>, sumOp<label>, sumOp<label>, sumOp<label>>{}
    );

    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].calcStats(qmetrics);
    }

    // Only report if there are some internal faces
    if (nDDotS > 0)
    {
        if (verbose && minDDotS < severeNonorthogonalityThreshold)
        {
            Info<< "Number of non-orthogonality errors: " << errorNonOrth
                << ". Number of severely non-orthogonal faces: "
                << severeNonOrth  << "." << endl;
        }
    }

    if (verbose)
    {
        if (nDDotS > 0)
        {
            Info<< "Mesh non-orthogonality Max: "
                << radToDeg(::acos(minDDotS))
                << " average: " << sumDDotS/nDDotS
                << endl;
        }
    }

    if (fieldIndex != -1)
    {
        const label combinedIndex = findStatIndex(meshStats, "meshQuality");
        if (combinedIndex != -1)
        {
            const volScalarField& mf = (*meshStats)[fieldIndex].field();
            volScalarField& cf = (*meshStats)[combinedIndex].field();
            forAll(cf, celli)
            {
                cf[celli] += min(mf[celli]/orthWarn,scalar(1.0));
            }
        }
    }

    if (errorNonOrth > 0)
    {
        if (verbose)
        {
            SeriousErrorInFunction
                << "Error in non-orthogonality detected" << endl;
        }

        return true;
    }
    else
    {
        if (verbose)
        {
            Info<< "Non-orthogonality check OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkFacePyramids
(
    const bool verbose,
    const scalar minPyrVol,
    const polyMesh& mesh,
    const pointField& faceCentres,
    const pointField& cellCentres,
    const vectorField& faceAreas,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr,
    PtrList<meshStatistics>* meshStats
)
{
    // check whether face area vector points to the cell with higher label
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();

    const faceList& f = mesh.faces();

    const word metric = "pyramids";

    const label fieldIndex = findStatIndex(meshStats, metric);

    DynamicList<scalar> qmetrics(checkFaces.size());
    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].field() = 1.0;
    }

    label nErrorPyrs = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        // Create the owner pyramid - it will have negative volume
        scalar pyrVol = primitiveMeshTools::pyramidVol
        (
            cellCentres[own[facei]],
            faceCentres[facei],
            faceAreas[facei]
        );

        if (fieldIndex != -1)
        {
            qmetrics.append(-pyrVol);
            (*meshStats)[fieldIndex].field()[own[facei]] =
                min
                (
                    (*meshStats)[fieldIndex].field()[own[facei]],
                    -pyrVol
                );
        }

        if (pyrVol > -minPyrVol)
        {
            if (verbose)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << facei << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVol
                    << " Face " << f[facei] << " area: " << f[facei].mag(p)
                    << " Owner cell: " << own[facei] << endl
                    << "Owner cell vertex labels: "
                    << mesh.cells()[own[facei]].labels(f)
                    << endl;
            }


            if (setPtr)
            {
                setPtr->insert(facei);
            }

            nErrorPyrs++;
        }

        if (mesh.isInternalFace(facei))
        {
            // Create the neighbour pyramid - it will have positive volume
            scalar pyrVol = primitiveMeshTools::pyramidVol
            (
                cellCentres[nei[facei]],
                faceCentres[facei],
                faceAreas[facei]
            );

            if (fieldIndex != -1)
            {
                qmetrics.append(pyrVol);
                (*meshStats)[fieldIndex].field()[nei[facei]] =
                    min
                    (
                        (*meshStats)[fieldIndex].field()[nei[facei]],
                        pyrVol
                    );
            }

            if (pyrVol < minPyrVol)
            {
                if (verbose)
                {
                    Pout<< "bool polyMeshGeometry::checkFacePyramids("
                        << "const bool, const scalar, const pointField&"
                        << ", const labelList&, labelHashSet*): "
                        << "face " << facei << " points the wrong way. " << endl
                        << "Pyramid volume: " << -pyrVol
                        << " Face " << f[facei] << " area: " << f[facei].mag(p)
                        << " Neighbour cell: " << nei[facei] << endl
                        << "Neighbour cell vertex labels: "
                        << mesh.cells()[nei[facei]].labels(f)
                        << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nErrorPyrs++;
            }
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];

        // Create the owner pyramid - it will have negative volume
        scalar pyrVolOwn = primitiveMeshTools::pyramidVol
        (
             ownCc,
             faceCentres[face0],
             faceAreas[face0]
        );

        if (pyrVolOwn > -minPyrVol)
        {
            if (verbose)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << face0 << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVolOwn
                    << " Face " << f[face0] << " area: " << f[face0].mag(p)
                    << " Owner cell: " << own[face0] << endl
                    << "Owner cell vertex labels: "
                    << mesh.cells()[own[face0]].labels(f)
                    << endl;
            }


            if (setPtr)
            {
                setPtr->insert(face0);
            }

            nErrorPyrs++;
        }

        // Create the neighbour pyramid - it will have positive volume
        scalar pyrVolNbr = primitiveMeshTools::pyramidVol
        (
             cellCentres[own[face1]],
             faceCentres[face0],
             faceAreas[face0]
        );

        if (pyrVolNbr < minPyrVol)
        {
            if (verbose)
            {
                Pout<< "bool polyMeshGeometry::checkFacePyramids("
                    << "const bool, const scalar, const pointField&"
                    << ", const labelList&, labelHashSet*): "
                    << "face " << face0 << " points the wrong way. " << endl
                    << "Pyramid volume: " << -pyrVolNbr
                    << " Face " << f[face0] << " area: " << f[face0].mag(p)
                    << " Neighbour cell: " << own[face1] << endl
                    << "Neighbour cell vertex labels: "
                    << mesh.cells()[own[face1]].labels(f)
                    << endl;
            }

            if (setPtr)
            {
                setPtr->insert(face1);
            }

            nErrorPyrs++;
        }
    }

    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].calcStats(qmetrics);
    }

    if (fieldIndex != -1)
    {
        const label combinedIndex = findStatIndex(meshStats, "meshQuality");
        if (combinedIndex != -1)
        {
            const volScalarField& mf = (*meshStats)[fieldIndex].field();
            volScalarField& cf = (*meshStats)[combinedIndex].field();
            forAll(cf, celli)
            {
                scalar cellQ
                (
                    mf[celli] < minPyrVol ? scalar(1.0) : scalar(0.0)
                );
                cf[celli] += cellQ;
            }
        }
    }

    reduce(nErrorPyrs, sumOp<label>());
    if (nErrorPyrs > 0)
    {
        if (verbose)
        {
            SeriousErrorInFunction
                << "Error in face pyramids: faces pointing the wrong way."
                << endl;
        }

        return true;
    }
    else
    {
        if (verbose)
        {
            Info<< "Face pyramids OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceTets
(
    const bool verbose,
    const scalar minTetQuality,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
)
{
    // check whether decomposing each cell into tets results in
    // positive volume, non-flat tets
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces() - mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei - mesh.nInternalFaces()] = cellCentres[own[facei]];
    }

    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    label nErrorTets = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        // Create the owner pyramid - note: exchange cell and face centre
        // to get positive volume.
        bool tetError = checkFaceTet
        (
            mesh,
            verbose,
            minTetQuality,
            p,
            facei,
            cellCentres[own[facei]],    // face centre
            faceCentres[facei],         // cell centre
            setPtr
        );

        if (tetError)
        {
            nErrorTets++;
        }

        if (mesh.isInternalFace(facei))
        {
            // Create the neighbour tets - they will have positive volume
            bool tetError = checkFaceTet
            (
                mesh,
                verbose,
                minTetQuality,
                p,
                facei,
                faceCentres[facei],         // face centre
                cellCentres[nei[facei]],    // cell centre
                setPtr
            );

            if (tetError)
            {
                nErrorTets++;
            }

            if
            (
                polyMeshTetDecomposition::findSharedBasePoint
                (
                    mesh,
                    facei,
                    minTetQuality,
                    verbose
                ) == -1
            )
            {
                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nErrorTets++;
            }
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                if
                (
                    polyMeshTetDecomposition::findSharedBasePoint
                    (
                        mesh,
                        facei,
                        neiCc[facei - mesh.nInternalFaces()],
                        minTetQuality,
                        verbose
                    ) == -1
                )
                {
                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }

                    nErrorTets++;
                }
            }
            else
            {
                if
                (
                    polyMeshTetDecomposition::findBasePoint
                    (
                        mesh,
                        facei,
                        minTetQuality,
                        verbose
                    ) == -1
                )
                {
                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }

                    nErrorTets++;
                }
            }
        }
    }
    reduce(nErrorTets, sumOp<label>());

    if (nErrorTets > 0)
    {
        if (verbose)
        {
            SeriousErrorIn
            (
                "polyMeshGeometry::checkFaceTets("
                "const bool, const scalar, const pointField&, const pointField&"
                ", const labelList&, labelHashSet*)"
            )   << "Error in face decomposition: negative tets."
                << endl;
        }

        return true;
    }
    else
    {
        if (verbose)
        {
            Info<< "Face tets OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceWarpage
(
    const bool verbose,
    const scalar internalWarpage,
    const scalar boundaryWarpage,
    const polyMesh& mesh,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr,
    PtrList<meshStatistics>* meshStats
)
{
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const word metric = "warpage";

    const label fieldIndex = findStatIndex(meshStats, metric);
    DynamicList<scalar> qmetrics(checkFaces.size());
    autoPtr<PackedBoolList> isMasterFacePtr;
    if (fieldIndex != -1)
    {
        isMasterFacePtr.reset
        (
            new PackedBoolList(syncTools::getMasterFaces(mesh))
        );
        (*meshStats)[fieldIndex].field() = 0.;
    }

    scalar maxWarp = 0;

    label nWarnWarp = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];


        bool boundaryFace(true);

        if
        (
            mesh.isInternalFace(facei)
            || patches[patches.whichPatch(facei)].coupled()
        )
        {
            boundaryFace = false;
        }

        const point& fc = faceCentres[facei];
        scalar area = mag(faceAreas[facei]);

        if (area > SMALL)
        {
            const face& f = mesh.faces()[facei];
            vector fn = faceAreas[facei] / area;
            scalarField vProj(f.size());

            forAll(f, fp)
            {
                vector n = p[f[fp]] - fc;
                vProj[fp] = (n & fn);
            }
            // Get normal 'span' of face
            scalar minVal = min(vProj);
            scalar maxVal = max(vProj);

            scalar warp = (maxVal - minVal)/sqrt(area);
            maxWarp = max(maxWarp, warp);

            if (fieldIndex != -1)
            {
                if (isMasterFacePtr().get(facei))
                {
                    qmetrics.append(warp);
                }
                (*meshStats)[fieldIndex].field()[own[facei]] =
                    max
                    (
                        (*meshStats)[fieldIndex].field()[own[facei]],
                        warp
                     );
                if (mesh.isInternalFace(facei))
                {
                    (*meshStats)[fieldIndex].field()[nei[facei]] =
                        max
                        (
                            (*meshStats)[fieldIndex].field()[nei[facei]],
                            warp
                        );
                }
            }

            if
            (
                (!boundaryFace && warp > internalWarpage
                 && internalWarpage > 0)
                || (boundaryFace && warp > boundaryWarpage
                    && boundaryWarpage > 0)
            )
            {
                if (verbose)
                {
                    Pout<< "Severe warpage for face " << facei
                        << " warpage = " << warp << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnWarp++;
            }
        }
    }

    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].calcStats(qmetrics);
    }

    reduce(
        std::tie(maxWarp, nWarnWarp),
        ParallelOp<maxOp<scalar>, sumOp<label>>{}
    );

    if (nWarnWarp > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << "Large face warpage detected.  Max warpage = "
                << maxWarp
                << "\nThis may impair the quality of the result." << nl
                << nWarnWarp << " highly warped faces detected."
                << endl;
        }

        return true;
    }
    else
    {
        if (verbose)
        {
            Info<< "Max warpage = " << maxWarp
                << " Face warpage OK.\n" << endl;
        }

        return false;
    }

}

bool Foam::polyMeshGeometry::checkFaceSkewness
(
    const bool verbose,
    const scalar internalSkew,
    const scalar boundarySkew,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr,
    PtrList<meshStatistics>* meshStats
)
{
    // Warn if the skew correction vector is more than skew times
    // larger than the face area vector

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc;
    syncTools::swapBoundaryCellPositions(mesh, cellCentres, neiCc);
    const word metric = "skewness";
    DynamicList<scalar> qmetrics(checkFaces.size());
    const label fieldIndex = findStatIndex(meshStats, metric);
    autoPtr<PackedBoolList> isMasterFacePtr;
    if (fieldIndex != -1)
    {
        isMasterFacePtr.reset
        (
            new PackedBoolList(syncTools::getMasterFaces(mesh))
        );
        (*meshStats)[fieldIndex].field() = 0.;
    }

    scalar maxSkew = 0;

    label nWarnSkew = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        if (mesh.isInternalFace(facei))
        {
            scalar skewness = calcSkewness
            (
                mesh,
                cellCentres[own[facei]],
                cellCentres[nei[facei]],
                faceCentres[facei],
                faceAreas[facei],
                facei
            );

            if (fieldIndex != -1)
            {
                qmetrics.append(skewness);
                (*meshStats)[fieldIndex].field()[own[facei]] =
                    max
                    (
                        (*meshStats)[fieldIndex].field()[own[facei]],
                        skewness
                    );
                (*meshStats)[fieldIndex].field()[nei[facei]] =
                    max
                    (
                        (*meshStats)[fieldIndex].field()[nei[facei]],
                        skewness
                    );
            }

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > internalSkew)
            {
                if (verbose)
                {
                    Pout<< "Severe skewness for face " << facei
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnSkew++;
            }

            maxSkew = max(maxSkew, skewness);
        }
        else if (patches[patches.whichPatch(facei)].coupled())
        {
            scalar skewness = calcSkewness
            (
                mesh,
                cellCentres[own[facei]],
                neiCc[facei-mesh.nInternalFaces()],
                faceCentres[facei],
                faceAreas[facei],
                facei
            );

            if (fieldIndex != -1)
            {
                if (isMasterFacePtr().get(facei))
                {
                    qmetrics.append(skewness);
                }
                (*meshStats)[fieldIndex].field()[own[facei]] =
                    max
                    (
                        (*meshStats)[fieldIndex].field()[own[facei]],
                        skewness
                    );
            }

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > internalSkew)
            {
                if (verbose)
                {
                    Pout<< "Severe skewness for coupled face " << facei
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnSkew++;
            }

            maxSkew = max(maxSkew, skewness);
        }
        else
        {
            // Boundary faces: consider them to have only skewness error.
            // (i.e. treat as if mirror cell on other side)

            vector faceNormal = faceAreas[facei];
            faceNormal /= mag(faceNormal) + ROOTVSMALL;

            vector fC = faceCentres[facei];

            vector dOwn = fC - cellCentres[own[facei]];

            vector dWall = faceNormal*(faceNormal & dOwn);

            point faceIntersection = cellCentres[own[facei]] + dWall;

            scalar skewness = mag(fC - faceIntersection)
                /(2*mag(dWall) + ROOTVSMALL);

            if (skewness > 0.1)
            {
                vector dir = faceIntersection - fC;

                faceIntersection -= (dir & faceNormal) * faceNormal;

                dir = fC - faceIntersection;
                dir /= (mag(dir) + SMALL);

                const face& f = mesh.faces()[facei];

                scalar maxVal = -GREAT;
                scalar minVal = GREAT;

                forAll(f, fI)
                {
                    point pt = mesh.points()[f[fI]];
                    scalar dp = pt & dir;
                    maxVal = max(maxVal, dp);
                    minVal = min(minVal, dp);
                }
                scalar charDist = maxVal - minVal;

                skewness = max
                (
                    mag(fC - faceIntersection)/(charDist + VSMALL),skewness
                );
            }

            if (fieldIndex != -1)
            {
                qmetrics.append(skewness);
                const label patchi = patches.whichPatch(facei);
                const label index = facei - patches[patchi].start();
                if (!(*meshStats)[fieldIndex].field().boundaryFieldRef()[patchi].empty())
                {
                    (*meshStats)[fieldIndex].field().boundaryFieldRef()[patchi][index]
                        = skewness;
                }
                (*meshStats)[fieldIndex].field()[own[facei]] =
                    max
                    (
                        (*meshStats)[fieldIndex].field()[own[facei]],
                        skewness
                    );
            }

            // Check if the skewness vector is greater than the PN vector.
            // This does not cause trouble but is a good indication of a poor
            // mesh.
            if (skewness > boundarySkew)
            {
                if (verbose)
                {
                    Pout<< "Severe skewness for boundary face " << facei
                        << " skewness = " << skewness << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnSkew++;
            }

            maxSkew = max(maxSkew, skewness);
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];

        scalar skewness = calcSkewness
        (
            mesh,
            ownCc,
            cellCentres[own[face1]],
            faceCentres[face0],
            faceAreas[face0],
            face0
        );

        // Check if the skewness vector is greater than the PN vector.
        // This does not cause trouble but is a good indication of a poor
        // mesh.
        if (skewness > internalSkew)
        {
            if (verbose)
            {
                Pout<< "Severe skewness for face " << face0
                    << " skewness = " << skewness << endl;
            }

            if (setPtr)
            {
                setPtr->insert(face0);
                setPtr->insert(face1);
            }

            nWarnSkew++;
        }

        maxSkew = max(maxSkew, skewness);
    }
    reduce(
        std::tie(maxSkew, nWarnSkew),
        ParallelOp<maxOp<scalar>, sumOp<label>>{}
    );

    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].calcStats(qmetrics);
    }

    if (fieldIndex != -1)
    {
        const label combinedIndex = findStatIndex(meshStats, "meshQuality");
        if (combinedIndex != -1)
        {
            const volScalarField& mf = (*meshStats)[fieldIndex].field();
            volScalarField& cf = (*meshStats)[combinedIndex].field();
            forAll(cf, celli)
            {
                cf[celli] += min(mf[celli]/internalSkew,scalar(1.0));
            }
        }
    }

    if (nWarnSkew > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << 100*maxSkew
                << " percent.\nThis may impair the quality of the result." << nl
                << nWarnSkew << " highly skew faces detected."
                << endl;
        }

        return true;
    }
    else
    {
        if (verbose)
        {
            Info<< "Max skewness = " << 100*maxSkew
                << " percent.  Face skewness OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceWeights
(
    const bool verbose,
    const scalar warnWeight,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceCentres,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr,
    PtrList<meshStatistics>* meshStats
)
{
    // Warn if the delta factor (0..1) is too large.

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces()-mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    const word metric = "weights";

    const label fieldIndex = findStatIndex(meshStats, metric);
    DynamicList<scalar> qmetrics(checkFaces.size());
    autoPtr<PackedBoolList> isMasterFacePtr;
    if (fieldIndex != -1)
    {
        isMasterFacePtr.reset
        (
            new PackedBoolList(syncTools::getMasterFaces(mesh))
        );
        (*meshStats)[fieldIndex].field() = 0.5;
    }

    scalar minWeight = GREAT;

    label nWarnWeight = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const point& fc = faceCentres[facei];
        const vector& fa = faceAreas[facei];

        scalar dOwn = mag(fa & (fc-cellCentres[own[facei]]));

        if (mesh.isInternalFace(facei))
        {
            scalar dNei = mag(fa & (cellCentres[nei[facei]]-fc));
            scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

            if (fieldIndex != -1)
            {
                qmetrics.append(weight);
                (*meshStats)[fieldIndex].field()[own[facei]] =
                    min
                    (
                        (*meshStats)[fieldIndex].field()[own[facei]],
                        weight
                    );
                (*meshStats)[fieldIndex].field()[nei[facei]] =
                    min
                    (
                        (*meshStats)[fieldIndex].field()[nei[facei]],
                        weight
                    );
            }

            if (weight < warnWeight)
            {
                if (verbose)
                {
                    Pout<< "Small weighting factor for face " << facei
                        << " weight = " << weight << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnWeight++;
            }

            minWeight = min(minWeight, weight);
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                scalar dNei = mag(fa & (neiCc[facei-mesh.nInternalFaces()]-fc));
                scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

                if (fieldIndex != -1)
                {
                    if (isMasterFacePtr().get(facei))
                    {
                        qmetrics.append(weight);
                    }
                    (*meshStats)[fieldIndex].field()[own[facei]] =
                        min
                        (
                            (*meshStats)[fieldIndex].field()[own[facei]],
                            weight
                        );
                }

                if (weight < warnWeight)
                {
                    if (verbose)
                    {
                        Pout<< "Small weighting factor for face " << facei
                            << " weight = " << weight << endl;
                    }

                    if (setPtr)
                    {
                        setPtr->insert(facei);
                    }

                    nWarnWeight++;
                }

                minWeight = min(minWeight, weight);
            }
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        const point& ownCc = cellCentres[own[face0]];
        const point& fc = faceCentres[face0];
        const vector& fa = faceAreas[face0];

        scalar dOwn = mag(fa & (fc-ownCc));
        scalar dNei = mag(fa & (cellCentres[own[face1]]-fc));
        scalar weight = min(dNei,dOwn)/(dNei+dOwn+VSMALL);

        if (weight < warnWeight)
        {
            if (verbose)
            {
                Pout<< "Small weighting factor for face " << face0
                    << " weight = " << weight << endl;
            }

            if (setPtr)
            {
                setPtr->insert(face0);
                setPtr->insert(face1);
            }

            nWarnWeight++;
        }
        minWeight = min(minWeight, weight);
    }

    reduce(
        std::tie(minWeight, nWarnWeight),
        ParallelOp<minOp<scalar>, sumOp<label>>{}
    );

    if (fieldIndex != -1)
    {
        const label combinedIndex = findStatIndex(meshStats, "meshQuality");
        if (combinedIndex != -1)
        {
            const volScalarField& mf = (*meshStats)[fieldIndex].field();
            volScalarField& cf = (*meshStats)[combinedIndex].field();
            forAll(cf, celli)
            {
                if (mf[celli] <= warnWeight)
                {
                    cf[celli] += scalar(1.0);
                }
                else
                {
                    cf[celli] += 1.0-max(mf[celli]/0.5,scalar(1.0));
                }
            }
        }
    }

    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].calcStats(qmetrics);
    }

    if (minWeight < warnWeight)
    {
        if (verbose)
        {
            WarningInFunction
                << minWeight << '.' << nl
                << nWarnWeight << " faces with small weights detected."
                << endl;
        }

        return true;
    }
    else
    {
        if (verbose)
        {
            Info<< "Min weight = " << minWeight
                << ".  Weights OK.\n" << endl;
        }

        return false;
    }
}


bool Foam::polyMeshGeometry::checkVolRatio
(
    const bool verbose,
    const scalar warnRatio,
    const polyMesh& mesh,
    const scalarField& cellVolumes,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr,
    PtrList<meshStatistics>* meshStats
)
{
    // Warn if the volume ratio between neighbouring cells is too large

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell vol
    scalarField neiVols(mesh.nFaces()-mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiVols[facei-mesh.nInternalFaces()] = cellVolumes[own[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiVols);

    const word metric = "volumeRatio";
    DynamicList<scalar> qmetrics(checkFaces.size());
    const label fieldIndex = findStatIndex(meshStats, metric);
    autoPtr<PackedBoolList> isMasterFacePtr;
    if (fieldIndex != -1)
    {
        isMasterFacePtr.reset
        (
            new PackedBoolList(syncTools::getMasterFaces(mesh))
        );
        (*meshStats)[fieldIndex].field() = 1.;
    }

    scalar minRatio = GREAT;

    label nWarnRatio = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        scalar ownVol = mag(cellVolumes[own[facei]]);

        scalar neiVol = -GREAT;

        if (mesh.isInternalFace(facei))
        {
            neiVol = mag(cellVolumes[nei[facei]]);
        }
        else
        {
            label patchi = patches.whichPatch(facei);

            if (patches[patchi].coupled())
            {
                neiVol = mag(neiVols[facei-mesh.nInternalFaces()]);
            }
        }

        if (neiVol >= 0)
        {
            scalar ratio = min(ownVol, neiVol) / (max(ownVol, neiVol) + VSMALL);

            if (fieldIndex != -1)
            {
                if (isMasterFacePtr().get(facei))
                {
                    qmetrics.append(ratio);
                }
                (*meshStats)[fieldIndex].field()[own[facei]] =
                    min
                    (
                        (*meshStats)[fieldIndex].field()[own[facei]],
                        ratio
                    );
                if (mesh.isInternalFace(facei))
                {
                    (*meshStats)[fieldIndex].field()[nei[facei]] =
                        min
                        (
                            (*meshStats)[fieldIndex].field()[nei[facei]],
                            ratio
                        );
                }
            }

            if (ratio < warnRatio)
            {
                if (verbose)
                {
                    Pout<< "Small ratio for face " << facei
                        << " ratio = " << ratio << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(facei);
                }

                nWarnRatio++;
            }

            minRatio = min(minRatio, ratio);
        }
    }

    forAll(baffles, i)
    {
        label face0 = baffles[i].first();
        label face1 = baffles[i].second();

        scalar ownVol = mag(cellVolumes[own[face0]]);

        scalar neiVol = mag(cellVolumes[own[face1]]);

        if (neiVol >= 0)
        {
            scalar ratio = min(ownVol, neiVol) / (max(ownVol, neiVol) + VSMALL);

            if (ratio < warnRatio)
            {
                if (verbose)
                {
                    Pout<< "Small ratio for face " << face0
                        << " ratio = " << ratio << endl;
                }

                if (setPtr)
                {
                    setPtr->insert(face0);
                    setPtr->insert(face1);
                }

                nWarnRatio++;
            }

            minRatio = min(minRatio, ratio);
        }
    }

    reduce(
        std::tie(minRatio, nWarnRatio),
        ParallelOp<minOp<scalar>, sumOp<label>>{}
    );

    if (fieldIndex != -1)
    {
        const label combinedIndex = findStatIndex(meshStats, "meshQuality");
        if (combinedIndex != -1)
        {
            const volScalarField& mf = (*meshStats)[fieldIndex].field();
            volScalarField& cf = (*meshStats)[combinedIndex].field();
            forAll(cf, celli)
            {
                if (mf[celli] <= warnRatio)
                {
                    cf[celli] += scalar(1.0);
                }
                else
                {
                    cf[celli] += 1.0-mf[celli];
                }
            }
        }
    }

    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].calcStats(qmetrics);
    }

    if (minRatio < warnRatio)
    {
        if (verbose)
        {
            WarningInFunction
                << minRatio << '.' << nl
                << nWarnRatio << " faces with small ratios detected."
                << endl;
        }

        return true;
    }
    else
    {
        if (verbose)
        {
            Info<< "Min ratio = " << minRatio
                << ".  Ratios OK.\n" << endl;
        }

        return false;
    }
}


// Check convexity of angles in a face. Allow a slight non-convexity.
// E.g. maxDeg = 10 allows for angles < 190 (or 10 degrees concavity)
// (if truly concave and points not visible from face centre the face-pyramid
//  check in checkMesh will fail)
bool Foam::polyMeshGeometry::checkFaceAngles
(
    const bool verbose,
    const scalar maxDeg,
    const scalar maxPyrDeg,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (maxDeg < -SMALL || maxDeg > 180+SMALL)
    {
        FatalErrorInFunction
            << "maxDeg should be [0..180] but is now " << maxDeg
            << abort(FatalError);
    }

    bool checkPyramids = false;
    scalar maxPyrSin = -1;
    if (maxPyrDeg > 0 && maxPyrDeg < 180.0-SMALL)
    {
        checkPyramids = true;
        maxPyrSin = Foam::sin(degToRad(maxPyrDeg));
    }

    const scalar maxSin = Foam::sin(degToRad(maxDeg));

    const faceList& fcs = mesh.faces();

    scalar maxEdgeSin = 0.0;

    label nConcave = 0;

    label errorfacei = -1;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const face& f = fcs[facei];

        vector faceNormal = faceAreas[facei];
        faceNormal /= mag(faceNormal) + VSMALL;

        // Get edge from f[0] to f[size-1];
        vector ePrev(p[f.first()] - p[f.last()]);
        scalar magEPrev = mag(ePrev);
        ePrev /= magEPrev + VSMALL;

        forAll(f, fp0)
        {
            // Get vertex after fp
            label fp1 = f.fcIndex(fp0);

            // Normalized vector between two consecutive points
            vector e10(p[f[fp1]] - p[f[fp0]]);
            scalar magE10 = mag(e10);
            e10 /= magE10 + VSMALL;

            if (magEPrev > SMALL && magE10 > SMALL)
            {
                vector edgeNormal = ePrev ^ e10;
                scalar magEdgeNormal = mag(edgeNormal);

                if (magEdgeNormal < maxSin)
                {
                    // Edges (almost) aligned -> face is ok.

                    // If required perform additional pyramid check
                    if (checkPyramids && f.size() != 3)
                    {
                        label own = mesh.faceOwner()[facei];

                        bool pyrCell(isPyrCell(mesh,own));
                        if (!pyrCell && mesh.isInternalFace(facei))
                        {
                            label nei = mesh.faceNeighbour()[facei];
                            pyrCell = isPyrCell(mesh,nei);
                        }

                        if (pyrCell && magEdgeNormal >= maxPyrSin)
                        {
                            // Check normal
                            edgeNormal /= magEdgeNormal;

                            if ((edgeNormal & faceNormal) < SMALL)
                            {
                                if (facei != errorfacei)
                                {
                                    // Count only one error per face.
                                    errorfacei = facei;
                                    nConcave++;
                                }

                                if (setPtr)
                                {
                                    setPtr->insert(facei);
                                }

                                maxEdgeSin = max(maxEdgeSin, magEdgeNormal);
                            }
                        }
                    }
                }
                else
                {
                    // Check normal
                    edgeNormal /= magEdgeNormal;

                    if ((edgeNormal & faceNormal) < SMALL)
                    {
                        if (facei != errorfacei)
                        {
                            // Count only one error per face.
                            errorfacei = facei;
                            nConcave++;
                        }

                        if (setPtr)
                        {
                            setPtr->insert(facei);
                        }

                        maxEdgeSin = max(maxEdgeSin, magEdgeNormal);
                    }
                }
            }

            ePrev = e10;
            magEPrev = magE10;
        }
    }

    reduce(
        std::tie(nConcave, maxEdgeSin),
        ParallelOp<sumOp<label>, maxOp<scalar>>{}
    );

    if (verbose)
    {
        if (maxEdgeSin > SMALL)
        {
            scalar maxConcaveDegr =
                radToDeg(Foam::asin(Foam::min(1.0, maxEdgeSin)));

            Info<< "There are " << nConcave
                << " faces with concave angles between consecutive"
                << " edges. Max concave angle = " << maxConcaveDegr
                << " degrees.\n" << endl;
        }
        else
        {
            Info<< "All angles in faces are convex or less than "  << maxDeg
                << " degrees concave.\n" << endl;
        }
    }

    if (nConcave > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nConcave  << " face points with severe concave angle (> "
                << min(maxPyrDeg,maxDeg) << " deg) found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// Check twist of faces. Is calculated as the difference between normals of
// individual triangles and the cell-cell centre edge.
bool Foam::polyMeshGeometry::checkFaceTwist
(
    const bool verbose,
    const scalar minTwist,
    const polyMesh& mesh,
    const vectorField& cellCentres,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minTwist < -1-SMALL || minTwist > 1+SMALL)
    {
        FatalErrorInFunction
            << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }


    const faceList& fcs = mesh.faces();


    label nWarped = 0;

//    forAll(checkFaces, i)
//    {
//        label facei = checkFaces[i];
//
//        const face& f = fcs[facei];
//
//        scalar magArea = mag(faceAreas[facei]);
//
//        if (f.size() > 3 && magArea > VSMALL)
//        {
//            const vector nf = faceAreas[facei] / magArea;
//
//            const point& fc = faceCentres[facei];
//
//            forAll(f, fpI)
//            {
//                vector triArea
//                (
//                    triPointRef
//                    (
//                        p[f[fpI]],
//                        p[f.nextLabel(fpI)],
//                        fc
//                    ).areaNormal()
//                );
//
//                scalar magTri = mag(triArea);
//
//                if (magTri > VSMALL && ((nf & triArea/magTri) < minTwist))
//                {
//                    nWarped++;
//
//                    if (setPtr)
//                    {
//                        setPtr->insert(facei);
//                    }
//
//                    break;
//                }
//            }
//        }
//    }


    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Calculate coupled cell centre
    pointField neiCc(mesh.nFaces()-mesh.nInternalFaces());

    for (label facei = mesh.nInternalFaces(); facei < mesh.nFaces(); facei++)
    {
        neiCc[facei-mesh.nInternalFaces()] = cellCentres[own[facei]];
    }
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            vector nf(Zero);

            if (mesh.isInternalFace(facei))
            {
                nf = cellCentres[nei[facei]] - cellCentres[own[facei]];
                nf /= mag(nf) + VSMALL;
            }
            else if (patches[patches.whichPatch(facei)].coupled())
            {
                nf =
                    neiCc[facei-mesh.nInternalFaces()]
                  - cellCentres[own[facei]];
                nf /= mag(nf) + VSMALL;
            }
            else
            {
                nf = faceCentres[facei] - cellCentres[own[facei]];
                nf /= mag(nf) + VSMALL;
            }

            if (nf != vector::zero)
            {
                const point& fc = faceCentres[facei];

                forAll(f, fpI)
                {
                    vector triArea
                    (
                        triPointRef
                        (
                            p[f[fpI]],
                            p[f.nextLabel(fpI)],
                            fc
                        ).areaNormal()
                    );

                    scalar magTri = mag(triArea);
                    scalar eLen = mag(p[f[fpI]] - p[f.nextLabel(fpI)]);

                    if
                    (
                        magTri > VSMALL && eLen > SMALL
                        && ((nf & triArea/magTri) < minTwist)
                    )
                    {
                        nWarped++;

                        if (setPtr)
                        {
                            setPtr->insert(facei);
                        }

                        break;
                    }
                }
            }
        }
    }

    reduce(nWarped, sumOp<label>());

    if (verbose)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with cosine of the angle"
                << " between triangle normal and face normal less than "
                << minTwist << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the cosine of the angle"
                << " between triangle normal and face normal less than "
                << minTwist << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nWarped  << " faces with severe warpage "
                << "(cosine of the angle between triangle normal and "
                << "face normal < " << minTwist << ") found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// Like checkFaceTwist but compares normals of consecutive triangles.
bool Foam::polyMeshGeometry::checkTriangleTwist
(
    const bool verbose,
    const scalar minTwist,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minTwist < -1-SMALL || minTwist > 1+SMALL)
    {
        FatalErrorInFunction
            << "minTwist should be [-1..1] but is now " << minTwist
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            const point& fc = faceCentres[facei];

            // Find starting triangle (at startFp) with non-zero area
            label startFp = -1;
            vector prevN;

            forAll(f, fp)
            {
                prevN = triPointRef
                (
                    p[f[fp]],
                    p[f.nextLabel(fp)],
                    fc
                ).areaNormal();

                scalar magTri = mag(prevN);

                if (magTri > VSMALL)
                {
                    startFp = fp;
                    prevN /= magTri;
                    break;
                }
            }

            if (startFp != -1)
            {
                label fp = startFp;
                bool ignoreNext = false;

                do
                {
                    fp = f.fcIndex(fp);

                    vector triN
                    (
                        triPointRef
                        (
                            p[f[fp]],
                            p[f.nextLabel(fp)],
                            fc
                        ).areaNormal()
                    );
                    scalar magTri = mag(triN);
                    scalar eLenNext = mag(p[f[fp]]-p[f.nextLabel(fp)]);
                    scalar eLenPrev = mag(p[f[fp]]-p[f.prevLabel(fp)]);

                    if (eLenNext > SMALL && eLenPrev > SMALL)
                    {
                        if (!ignoreNext)
                        {
                            if (magTri > VSMALL)
                            {
                                triN /= magTri;
                                if ((prevN & triN) < minTwist)
                                {
                                    nWarped++;
                                    if (setPtr)
                                    {
                                        setPtr->insert(facei);
                                    }

                                    break;
                                }

                                prevN = triN;
                            }
                            else if (minTwist > 0)
                            {
                                nWarped++;

                                if (setPtr)
                                {
                                    setPtr->insert(facei);
                                }

                                break;
                            }
                        }
                        else
                        {
                            if (magTri > VSMALL)
                            {
                                triN /= magTri;
                                prevN = triN;
                                ignoreNext = false;
                            }
                        }
                    }
                    else
                    {
                        ignoreNext = true;
                    }
                }
                while (fp != startFp);
            }
        }
    }

    reduce(nWarped, sumOp<label>());

    if (verbose)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with cosine of the angle"
                << " between consecutive triangle normals less than "
                << minTwist << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the cosine of the angle"
                << " between consecutive triangle normals is less than "
                << minTwist << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nWarped  << " faces with severe warpage "
                << "(cosine of the angle between consecutive triangle normals"
                << " < " << minTwist << ") found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceFlatness
(
    const bool verbose,
    const scalar minFlatness,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const vectorField& faceCentres,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    if (minFlatness < -SMALL || minFlatness > 1+SMALL)
    {
        FatalErrorInFunction
            << "minFlatness should be [0..1] but is now " << minFlatness
            << abort(FatalError);
    }

    const faceList& fcs = mesh.faces();

    label nWarped = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const face& f = fcs[facei];

        if (f.size() > 3)
        {
            const point& fc = faceCentres[facei];

            // Sum triangle areas
            scalar sumArea = 0.0;

            forAll(f, fp)
            {
                sumArea += triPointRef
                (
                    p[f[fp]],
                    p[f.nextLabel(fp)],
                    fc
                ).mag();
            }

            if (sumArea/mag(faceAreas[facei]) < minFlatness)
            {
                nWarped++;

                if (setPtr)
                {
                    setPtr->insert(facei);
                }
            }
        }
    }

    reduce(nWarped, sumOp<label>());

    if (verbose)
    {
        if (nWarped> 0)
        {
            Info<< "There are " << nWarped
                << " faces with area of invidual triangles"
                << " compared to overall area less than "
                << minFlatness << nl << endl;
        }
        else
        {
            Info<< "All faces are flat in that the area of invidual triangles"
                << " compared to overall area is less than "
                << minFlatness << nl << endl;
        }
    }

    if (nWarped > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nWarped  << " non-flat faces "
                << "(area of invidual triangles"
                << " compared to overall area"
                << " < " << minFlatness << ") found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceArea
(
    const bool verbose,
    const scalar minArea,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    label nZeroArea = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        if (mag(faceAreas[facei]) < minArea)
        {
            if (setPtr)
            {
                setPtr->insert(facei);
            }
            nZeroArea++;
        }
    }


    reduce(nZeroArea, sumOp<label>());

    if (verbose)
    {
        if (nZeroArea > 0)
        {
            Info<< "There are " << nZeroArea
                << " faces with area < " << minArea << '.' << nl << endl;
        }
        else
        {
            Info<< "All faces have area > " << minArea << '.' << nl << endl;
        }
    }

    if (nZeroArea > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nZeroArea  << " faces with area < " << minArea
                << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceEdgeLengths
(
    const bool verbose,
    const scalar minEdgeLength,
    const polyMesh& mesh,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    label nZeroEdgeLength = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];

        const labelList& fEdges = mesh.faceEdges()[facei];
        forAll(fEdges, fEI)
        {
            label edgeI = fEdges[fEI];
            edge e = mesh.edges()[edgeI];
            if (e.mag(mesh.points()) < minEdgeLength)
            {
                if (setPtr)
                {
                    setPtr->insert(facei);
                }
                nZeroEdgeLength++;
                break;
            }
        }
    }
    reduce(nZeroEdgeLength, sumOp<label>());

    if (verbose)
    {
        if (nZeroEdgeLength > 0)
        {
            Info<< "There are " << nZeroEdgeLength
                << " faces with minimum edge length < "
                << minEdgeLength << '.' << nl << endl;
        }
        else
        {
            Info<< "All faces have minimum edge length > "
                << minEdgeLength << '.' << nl << endl;
        }
    }

    if (nZeroEdgeLength > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nZeroEdgeLength  << " faces with area < " << minEdgeLength
                << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkCellDeterminant
(
    const bool verbose,
    const scalar warnDet,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr,
    PtrList<meshStatistics>* meshStats
)
{
    const word metric = "determinant";

    const label fieldIndex = findStatIndex(meshStats, metric);
    DynamicList<scalar> qmetrics(affectedCells.size());
    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].field() = 1.;
    }

    const cellList& cells = mesh.cells();

    scalar minDet = GREAT;
    scalar sumDet = 0.0;
    label nSumDet = 0;
    label nWarnDet = 0;

    forAll(affectedCells, i)
    {
        const cell& cFaces = cells[affectedCells[i]];

        tensor areaSum(Zero);
        scalar magAreaSum = 0;

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];

            scalar magArea = mag(faceAreas[facei]);

            magAreaSum += magArea;
            areaSum += faceAreas[facei]*(faceAreas[facei]/(magArea+VSMALL));
        }

        scalar scaledDet = det(areaSum/(magAreaSum+VSMALL))/0.037037037037037;
        if (fieldIndex != -1)
        {
            qmetrics.append(scaledDet);
            (*meshStats)[fieldIndex].field()[affectedCells[i]] = scaledDet;
        }

        minDet = min(minDet, scaledDet);
        sumDet += scaledDet;
        nSumDet++;

        if (scaledDet < warnDet)
        {
            if (setPtr)
            {
                // Insert all faces of the cell.
                forAll(cFaces, cFacei)
                {
                    label facei = cFaces[cFacei];
                    setPtr->insert(facei);
                }
            }
            nWarnDet++;
        }
    }

    reduce(
        std::tie(minDet, sumDet, nSumDet, nWarnDet),
        ParallelOp<minOp<scalar>, sumOp<scalar>, sumOp<label>, sumOp<label>>{}
    );

    if (fieldIndex != -1)
    {
        const label combinedIndex = findStatIndex(meshStats, "meshQuality");
        if (combinedIndex != -1)
        {
            const volScalarField& mf = (*meshStats)[fieldIndex].field();
            volScalarField& cf = (*meshStats)[combinedIndex].field();
            forAll(cf, celli)
            {
                if (mf[celli] <= warnDet)
                {
                    cf[celli] += scalar(1.0);
                }
                else
                {
                    cf[celli] += 1.0-mf[celli];
                }
            }
        }
    }

    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].calcStats(qmetrics);
    }

    if (verbose)
    {
        if (nSumDet > 0)
        {
            Info<< "Cell determinant (1 = uniform cube) : average = "
                << sumDet / nSumDet << "  min = " << minDet << endl;
        }

        if (nWarnDet > 0)
        {
            Info<< "There are " << nWarnDet
                << " cells with determinant < " << warnDet << '.' << nl
                << endl;
        }
        else
        {
            Info<< "All faces have determinant > " << warnDet << '.' << nl
                << endl;
        }
    }

    if (nWarnDet > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nWarnDet << " cells with determinant < " << warnDet
                << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::polyMeshGeometry::checkCellAspectRatio
(
    const bool verbose,
    const scalar warnAR,
    const polyMesh& mesh,
    const vectorField& faceAreas,
    const scalarField& cellVolumes,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr,
    PtrList<meshStatistics>* meshStats
)
{
    const cellList& cells = mesh.cells();
    const Vector<label>& meshD = mesh.geometricD();
    label nDims = 0;
    for (direction dir = 0; dir < vector::nComponents; dir++)
    {
        if (meshD[dir] == 1)
        {
            nDims++;
        }
    }

    const word metric = "aspectRatio";

    const label fieldIndex = findStatIndex(meshStats, metric);
    DynamicList<scalar> qmetrics(affectedCells.size());
    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].field() = 1.0;
    }

    scalar maxAR = -GREAT;
    scalar sumAR = 0.0;
    label nSumAR = 0;
    label nWarnAR = 0;

    vectorField sumMagClosed(mesh.nCells(), Zero);

    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    forAll(own, facei)
    {
        // Add to owner
        sumMagClosed[own[facei]] += cmptMag(faceAreas[facei]);
    }

    forAll(nei, facei)
    {
        // Subtract from neighbour
        sumMagClosed[nei[facei]] += cmptMag(faceAreas[facei]);
    }

    forAll(affectedCells, i)
    {
        label cellI = affectedCells[i];
       // Calculate the aspect ration as the maximum of Cartesian component
        // aspect ratio to the total area hydraulic area aspect ratio
        scalar minCmpt = VGREAT;
        scalar maxCmpt = -VGREAT;
        for (direction dir = 0; dir < vector::nComponents; dir++)
        {
            if (meshD[dir] == 1)
            {
                minCmpt = min(minCmpt, sumMagClosed[cellI][dir]);
                maxCmpt = max(maxCmpt, sumMagClosed[cellI][dir]);
            }
        }

        scalar metric = maxCmpt/(minCmpt + ROOTVSMALL);
        if (nDims == 3)
        {
            scalar v = max(ROOTVSMALL, cellVolumes[cellI]);

            metric = max
            (
                metric,
                1.0/6.0*cmptSum(sumMagClosed[cellI])/pow(v, 2.0/3.0)
            );
        }

        if (fieldIndex != -1)
        {
            qmetrics.append(metric);
            (*meshStats)[fieldIndex].field()[cellI] = metric;
        }

        maxAR = max(maxAR, metric);
        sumAR += maxAR;
        nSumAR++;

        if (metric > warnAR)
        {
            if (setPtr)
            {
                const cell& cFaces = cells[cellI];
                // Insert all faces of the cell.
                forAll(cFaces, cFaceI)
                {
                    label facei = cFaces[cFaceI];
                    setPtr->insert(facei);
                }
            }
            nWarnAR++;
        }

    }

    reduce(
        std::tie(maxAR, sumAR, nSumAR, nWarnAR),
        ParallelOp<maxOp<scalar>, sumOp<scalar>, sumOp<label>, sumOp<label>>{}
    );

    if (fieldIndex != -1)
    {
        (*meshStats)[fieldIndex].calcStats(qmetrics);
    }

    if (verbose)
    {
        if (nSumAR > 0)
        {
            Info<< "Cell Aspect Ratio : average = "
                << sumAR / nSumAR << "  max = " << maxAR << endl;
        }

        if (nWarnAR > 0)
        {
            Info<< "There are " << nWarnAR
                << " cells with Aspect Ratio > " << warnAR << '.' << nl
                << endl;
        }
        else
        {
            Info<< "All faces have Aspect Ratio < " << warnAR << '.' << nl
                << endl;
        }
    }

    if (nWarnAR > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nWarnAR << " cells with Aspect Ratio > " << warnAR
                << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkSnapVolume
(
    const bool verbose,
    const scalar minSnapRelativeVolume,
    const scalar minSnapRelativeTetVolume,
    const polyMesh& mesh,
    const scalarField& cellVolumes,
    const scalarField& initCellVolumes,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr
)
{
    const cellList& cells = mesh.cells();

    scalar minVRatio = GREAT;
    label nWarnRatio = 0;

    forAll(affectedCells, i)
    {
        label cellI = affectedCells[i];

        scalar vRatio = cellVolumes[cellI]/initCellVolumes[cellI];

        minVRatio = min(minVRatio, vRatio);

        const cell& cFaces = cells[cellI];

        if
        (
            vRatio < minSnapRelativeVolume
            || ((cFaces.size() == 4) && vRatio < minSnapRelativeTetVolume)
        )
        {
            if (setPtr)
            {
                // Insert all faces of the cell.
                forAll(cFaces, cFaceI)
                {
                    label facei = cFaces[cFaceI];
                    setPtr->insert(facei);
                }
            }
            nWarnRatio++;
        }
    }

    reduce(
        std::tie(minVRatio, nWarnRatio),
        ParallelOp<minOp<scalar>, sumOp<label>>{}
    );

    if (verbose)
    {
        if (minVRatio < GREAT)
        {
            Info<<"Minimum snapped volume ratio "<<minVRatio
                <<endl;
        }
        if (nWarnRatio > 0)
        {
            Info<< "There are " << nWarnRatio
                << " cells with an invalid snap volume change ratio. " << nl
                << endl;
        }
        else
        {
            Info<< "All cells have a valid snap volume change ratio." << nl
                << endl;
        }
    }

    if (nWarnRatio > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nWarnRatio << " cells with snap volume change ratio < "
                << minSnapRelativeVolume << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}

bool Foam::polyMeshGeometry::checkGaussGreenCentroid
(
    const bool verbose,
    const scalar maxGaussGreenCentroid,
    const polyMesh& mesh,
    const scalarField& initCellVolumes,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr
)
{
    const cellList& cells = mesh.cells();

    const vectorField& Sf(mesh.faceAreas());
    const vectorField& Cf(mesh.faceCentres());
    const vectorField& Cc(mesh.cellCentres());

    const labelUList& owner = mesh.faceOwner();
    const labelUList& neighbour = mesh.faceNeighbour();
    scalar onethird = 1.0/3.0;

    scalarField Vgg(mesh.nCells(),0.0);
    //calc cell volumes
    forAll(neighbour, facei)
    {
        Vgg[owner[facei]] += (Cf[facei] & Sf[facei]);
        Vgg[neighbour[facei]] -= (Cf[facei] & Sf[facei]);
    }

    forAll(mesh.boundaryMesh(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundaryMesh()[patchi].faceCells();

        label start = mesh.boundaryMesh()[patchi].start();
        forAll(mesh.boundaryMesh()[patchi], i)
        {
            label facei = start + i;
            Vgg[pFaceCells[i]] += (Cf[facei] & Sf[facei]);
        }
    }

    Vgg *= onethird;

    // calc cell centroids
    vectorField Cgg(mesh.nCells(), vector::zero);
    forAll(neighbour, facei)
    {
        Cgg[owner[facei]] += Cf[facei]*(Cf[facei] & Sf[facei]);
        Cgg[neighbour[facei]] -= Cf[facei]*(Cf[facei] & Sf[facei]);
    }

    forAll(mesh.boundaryMesh(), patchi)
    {
        const labelUList& pFaceCells =
            mesh.boundaryMesh()[patchi].faceCells();

        label start = mesh.boundaryMesh()[patchi].start();

        forAll(mesh.boundaryMesh()[patchi], i)
        {
            label facei = start + i;
            Cgg[pFaceCells[i]]
                += Cf[facei] * (Cf[facei] & Sf[facei]);
        }
    }

    scalar maxCCRatio = -GREAT;
    label nWarnRatio = 0;

    forAll(affectedCells, i)
    {
        label cellI = affectedCells[i];

        vector ggcc = Cgg[cellI] /(4*(mag(Vgg[cellI])+VSMALL));
        scalar charLength = pow(initCellVolumes[cellI],onethird);

        scalar ccRatio = mag(ggcc-Cc[cellI])/charLength;
        maxCCRatio = max(maxCCRatio, ccRatio);

        if (ccRatio > maxGaussGreenCentroid)
        {
            if (setPtr)
            {
                // Insert all faces of the cell.
                const cell& cFaces = cells[cellI];
                forAll(cFaces, cFaceI)
                {
                    label facei = cFaces[cFaceI];
                    setPtr->insert(facei);
                }
            }
            nWarnRatio++;
        }
    }

    reduce(
        std::tie(maxCCRatio, nWarnRatio),
        ParallelOp<maxOp<scalar>, sumOp<label>>{}
    );

    if (verbose)
    {
        if (maxCCRatio > -GREAT)
        {
            Info<<"Maximum snapped volume ratio "<<maxCCRatio
                <<endl;
        }
        if (nWarnRatio > 0)
        {
            Info<< "There are " << nWarnRatio
                << " cells with an an inavlid Gauss-Green Centroid. " << nl
                << endl;
        }
        else
        {
            Info<< "All cells have a valid Gauss-Green Centroid." << nl
                << endl;
        }
    }

    if (nWarnRatio > 0)
    {
        if (verbose)
        {
            WarningInFunction
                << nWarnRatio << " cells invalid Gauss-Green centroid < "
                << maxGaussGreenCentroid << " found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceFaceCells
(
    const bool verbose,
    const polyMesh& mesh,
    const labelList& checkFaces,
    labelHashSet* setPtr
)
{
    const cellList& cells = mesh.cells();
    const labelList& fOwn = mesh.faceOwner();
    const labelList& fNei = mesh.faceNeighbour();
    const labelListList& fEdges = mesh.faceEdges();
    const labelListList& eFaces = mesh.edgeFaces();
    labelList cellOwnerMarker(fOwn.size(), -1);
    labelList cellNeighMarker(fOwn.size(), -1);

    label nDegenerateFaceCells = 0;

    forAll(checkFaces, i)
    {
        label facei = checkFaces[i];
        label nAddOwnFaces = 0;

        label fownCellI = fOwn[facei];

        forAll(fEdges[facei], feI)
        {
            label cgEdge(fEdges[facei][feI]);

            forAll(eFaces[cgEdge], efI)
            {
                label cgFace(eFaces[cgEdge][efI]);

                if (cgFace != facei && cellOwnerMarker[cgFace] != facei)
                {
                    if (fOwn[cgFace] == fownCellI)
                    {
                        nAddOwnFaces++;
                        cellOwnerMarker[cgFace] = facei;
                    }
                    else if (cgFace < mesh.nInternalFaces())
                    {
                        if (fNei[cgFace] == fownCellI)
                        {
                            nAddOwnFaces++;
                            cellOwnerMarker[cgFace] = facei;
                        }
                    }
                }
            }
        }

        if (nAddOwnFaces < 3)
        {
            nDegenerateFaceCells++;
            if (setPtr)
            {
                // Insert all faces of the cell.
                const cell& cFaces = cells[fOwn[facei]];
                forAll(cFaces, cFaceI)
                {
                    label fI = cFaces[cFaceI];
                    setPtr->insert(fI);
                }
            }
        }


        if (facei < mesh.nInternalFaces())
        {
            label nAddNeiFaces = 0;

            label fneiCellI = fNei[facei];

            forAll(fEdges[facei], feI)
            {
                label cgEdge(fEdges[facei][feI]);

                forAll(eFaces[cgEdge], efI)
                {
                    label cgFace(eFaces[cgEdge][efI]);

                    if (cgFace != facei && cellNeighMarker[cgFace] != facei)
                    {
                        if (fOwn[cgFace] == fneiCellI)
                        {
                            nAddNeiFaces++;
                            cellNeighMarker[cgFace] = facei;
                        }
                        else if (cgFace < mesh.nInternalFaces())
                        {
                            if (fNei[cgFace] == fneiCellI)
                            {
                                nAddNeiFaces++;
                                cellNeighMarker[cgFace] = facei;
                            }
                        }
                    }
                }
            }

            if (nAddNeiFaces < 3)
            {
                nDegenerateFaceCells++;
                if (setPtr)
                {
                    // Insert all faces of the cell.
                    const cell& cFaces = cells[fneiCellI];
                    forAll(cFaces, cFaceI)
                    {
                        label fI = cFaces[cFaceI];
                        setPtr->insert(fI);
                    }
                }
            }
        }
    }

    reduce(nDegenerateFaceCells, sumOp<label>());

    if (verbose)
    {
        if (nDegenerateFaceCells > 0)
        {
            Info<< "There are " << nDegenerateFaceCells
                << " cells with face-face connectivity issues " << nl << endl;
        }
        else
        {
            Info<< "All cells have good face-face connectivity " << nl << endl;
        }
    }

    if (nDegenerateFaceCells > 0)
    {
        if (verbose)
        {
            WarningInFunction
                <<  nDegenerateFaceCells
                << " cells with face-face connectivity issues found.\n"
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::polyMeshGeometry::checkFaceDotProduct
(
    const bool verbose,
    const scalar orthWarn,
    const scalar maxFaceCentreNonOrtho,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceDotProduct
    (
        verbose,
        orthWarn,
        maxFaceCentreNonOrtho,
        mesh_,
        faceCentres_,
        cellCentres_,
        faceAreas_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFacePyramids
(
    const bool verbose,
    const scalar minPyrVol,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFacePyramids
    (
        verbose,
        minPyrVol,
        mesh_,
        faceCentres_,
        cellCentres_,
        faceAreas_,
        p,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceTets
(
    const bool verbose,
    const scalar minTetQuality,
    const pointField& p,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceTets
    (
        verbose,
        minTetQuality,
        mesh_,
        cellCentres_,
        faceCentres_,
        p,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceWarpage
(
    const bool verbose,
    const scalar internalWarpage,
    const scalar boundaryWarpage,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceWarpage
    (
        verbose,
        internalWarpage,
        boundaryWarpage,
        mesh_,
        faceCentres_,
        faceAreas_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceSkewness
(
    const bool verbose,
    const scalar internalSkew,
    const scalar boundarySkew,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceSkewness
    (
        verbose,
        internalSkew,
        boundarySkew,
        mesh_,
        cellCentres_,
        faceCentres_,
        faceAreas_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceWeights
(
    const bool verbose,
    const scalar warnWeight,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkFaceWeights
    (
        verbose,
        warnWeight,
        mesh_,
        cellCentres_,
        faceCentres_,
        faceAreas_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkVolRatio
(
    const bool verbose,
    const scalar warnRatio,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet* setPtr
) const
{
    return checkVolRatio
    (
        verbose,
        warnRatio,
        mesh_,
        cellVolumes_,
        checkFaces,
        baffles,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceAngles
(
    const bool verbose,
    const scalar maxDeg,
    const scalar maxPyrDeg,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceAngles
    (
        verbose,
        maxDeg,
        maxPyrDeg,
        mesh_,
        faceAreas_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceTwist
(
    const bool verbose,
    const scalar minTwist,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceTwist
    (
        verbose,
        minTwist,
        mesh_,
        cellCentres_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkTriangleTwist
(
    const bool verbose,
    const scalar minTwist,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkTriangleTwist
    (
        verbose,
        minTwist,
        mesh_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceFlatness
(
    const bool verbose,
    const scalar minFlatness,
    const pointField& p,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceFlatness
    (
        verbose,
        minFlatness,
        mesh_,
        faceAreas_,
        faceCentres_,
        p,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceArea
(
    const bool verbose,
    const scalar minArea,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceArea
    (
        verbose,
        minArea,
        mesh_,
        faceAreas_,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceEdgeLengths
(
    const bool verbose,
    const scalar minEdgeLength,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceEdgeLengths
    (
        verbose,
        minEdgeLength,
        mesh_,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkFaceFaceCells
(
    const bool verbose,
    const labelList& checkFaces,
    labelHashSet* setPtr
) const
{
    return checkFaceFaceCells
    (
        verbose,
        mesh_,
        checkFaces,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkCellDeterminant
(
    const bool verbose,
    const scalar warnDet,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr
) const
{
    return checkCellDeterminant
    (
        verbose,
        warnDet,
        mesh_,
        faceAreas_,
        checkFaces,
        affectedCells,
        setPtr
    );
}


bool Foam::polyMeshGeometry::checkCellAspectRatio
(
    const bool verbose,
    const scalar maxAspectRatio,
    const labelList& checkFaces,
    const labelList& affectedCells,
    labelHashSet* setPtr
) const
{
    return checkCellAspectRatio
    (
        verbose,
        maxAspectRatio,
        mesh_,
        faceAreas_,
        cellVolumes_,
        checkFaces,
        affectedCells,
        setPtr
    );
}


// ************************************************************************* //
