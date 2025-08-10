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
    (c) 2011-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "motionSmoother/motionSmootherAlgo.H"
#include "motionSmoother/polyMeshGeometry/polyMeshGeometry.H"
#include "db/IOstreams/IOstreams/IOmanip.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::motionSmootherAlgo::checkMesh
(
    const bool verbose,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    labelHashSet& wrongFaces,
    const bool report,
    PtrList<meshStatistics>* meshStats
)
{
    List<labelPair> emptyBaffles;
    return checkMesh
    (
        verbose,
        mesh,
        dict,
        checkFaces,
        emptyBaffles,
        wrongFaces,
        report,
        meshStats
    );
}

Foam::label Foam::motionSmootherAlgo::checkMesh
(
    const bool verbose,
    const polyMesh& mesh,
    const dictionary& dict,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& wrongFaces,
    const bool report,
    PtrList<meshStatistics>* meshStats
)
{
    const scalar maxNonOrtho
    (
        dict.lookupOrDefault<scalar>("maxNonOrtho", 180, true)
    );
    const scalar minVol
    (
        dict.lookupOrDefault<scalar>("minVol", -1E30, true)
    );
    const scalar minTetQuality
    (
        dict.lookupOrDefault<scalar>("minTetQuality", -1e30, true)
    );
    const scalar maxConcave
    (
        dict.lookupOrDefault<scalar>("maxConcave", 180, true)
    );
    const scalar maxPyrConcave
    (
        dict.lookupOrDefault<scalar>("maxPyrConcave", 180, true)
    );
    const scalar minArea
    (
        dict.lookupOrDefault<scalar>("minArea", -1e30, true)
    );
    const scalar minEdgeLength
    (
        dict.lookupOrDefault<scalar>("minEdgeLength", -1, true)
    );
    const scalar maxIntWarp
    (
        dict.lookupOrDefault<scalar>("maxInternalWarpage", -1, true)
    );
    const scalar maxBounWarp
    (
        dict.lookupOrDefault<scalar>("maxBoundaryWarpage", -1, true)
    );
    const scalar maxIntSkew
    (
        dict.lookupOrDefault<scalar>("maxInternalSkewness", -1, true)
    );
    const scalar maxBounSkew
    (
        dict.lookupOrDefault<scalar>("maxBoundarySkewness", -1, true)
    );
    const scalar minWeight
    (
        dict.lookupOrDefault<scalar>("minFaceWeight", -1, true)
    );
    const scalar minVolRatio
    (
        dict.lookupOrDefault<scalar>("minVolRatio", -1, true)
    );
    const scalar minTwist
    (
        dict.lookupOrDefault<scalar>("minTwist", -2, true)
    );
    const scalar minTriangleTwist
    (
        dict.lookupOrDefault<scalar>("minTriangleTwist", -2, true)
    );
    scalar minFaceFlatness = -1.0;
    dict.readIfPresent("minFaceFlatness", minFaceFlatness, true);
    const scalar minDet
    (
        dict.lookupOrDefault<scalar>("minDeterminant", -2, true)
    );
    const Switch faceFaceCells
    (
        dict.lookupOrDefault<Switch>("faceFaceCells", false)
    );
    const scalar minSnapRelativeVolume
    (
        dict.lookupOrDefault<scalar>("minSnapRelativeVolume", -1, true)
    );

    const scalar minSnapRelativeTetVolume
    (
        dict.lookupOrDefault<scalar>("minSnapRelativeTetVolume", -1, true)
    );

    const scalar maxGaussGreenCentroid
    (
        dict.lookupOrDefault<scalar>("maxGaussGreenCentroid", -1, true)
    );
    const scalar maxFaceCentreNonOrtho
    (
        dict.lookupOrDefault<scalar>("maxFaceCentreNonOrtho", 180)
    );

    const scalar maxCellAspectRatio
    (
        dict.lookupOrDefault<scalar>("maxCellAspectRatio", -1)
    );

    label nWrongFaces = 0;

    if (report)
    {
        Info<< "Checking faces in error :" << endl;
        //Pout.setf(ios_base::left);
    }

    label nBadDotProdFaces = 0;
    bool doDotProductCheck = maxNonOrtho < 180.0-SMALL;
    if (doDotProductCheck)
    {
        polyMeshGeometry::checkFaceDotProduct
        (
            verbose,
            maxNonOrtho,
            maxFaceCentreNonOrtho,
            mesh,
            mesh.faceCentres(),
            mesh.cellCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces,
            meshStats
        );

        nBadDotProdFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadDotProdFaces;
    }

    label nBadFacePyramidFaces = 0;
    bool doFacePyramidCheck = minVol > -GREAT;
    if (doFacePyramidCheck)
    {
        polyMeshGeometry::checkFacePyramids
        (
            verbose,
            minVol,
            mesh,
            mesh.faceCentres(),
            mesh.cellCentres(),
            mesh.faceAreas(),
            mesh.points(),
            checkFaces,
            baffles,
            &wrongFaces,
            meshStats
        );

        nBadFacePyramidFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFacePyramidFaces;
    }

    label nBadFaceTetFaces = 0;
    bool doFaceTetCheck = minTetQuality > -GREAT;
    if (doFaceTetCheck)
    {
        polyMeshGeometry::checkFaceTets
        (
            verbose,
            minTetQuality,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        nBadFaceTetFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceTetFaces;
    }

    label nBadFaceAngleFaces = 0;
    bool doFaceAngleCheck = maxConcave < 180.0-SMALL;
    if (doFaceAngleCheck)
    {
        polyMeshGeometry::checkFaceAngles
        (
            verbose,
            maxConcave,
            maxPyrConcave,
            mesh,
            mesh.faceAreas(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        nBadFaceAngleFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceAngleFaces;
    }

    label nBadFaceAreaFaces = 0;
    bool doFaceAreaCheck = minArea > -SMALL;
    if (minArea > -SMALL)
    {
        polyMeshGeometry::checkFaceArea
        (
            verbose,
            minArea,
            mesh,
            mesh.faceAreas(),
            checkFaces,
            &wrongFaces
        );

        nBadFaceAreaFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceAreaFaces;
    }

    label nBadEdgeLengthFaces = 0;
    bool doFaceEdgeLengthCheck = minEdgeLength > 0;
    if (minEdgeLength > 0)
    {
        polyMeshGeometry::checkFaceEdgeLengths
        (
            verbose,
            minEdgeLength,
            mesh,
            checkFaces,
            &wrongFaces
        );

        nBadEdgeLengthFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadEdgeLengthFaces;
    }

    label nBadFaceFaceCellsFaces = 0;
    if (faceFaceCells)
    {
        polyMeshGeometry::checkFaceFaceCells
        (
            verbose,
            mesh,
            checkFaces,
            &wrongFaces
        );

        nBadFaceFaceCellsFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceFaceCellsFaces;
    }

    label nBadFaceWarpageFaces = 0;
    bool doWarpageCheck = maxIntWarp > 0 || maxBounWarp > 0
                          || polyMeshGeometry::findStatIndex(meshStats,"warpage") != -1;
    if (doWarpageCheck)
    {
        polyMeshGeometry::checkFaceWarpage
        (
            verbose,
            maxIntWarp,
            maxBounWarp,
            mesh,
            mesh.faceCentres(),
            mesh.faceAreas(),
            mesh.points(),
            checkFaces,
            &wrongFaces,
            meshStats
        );

        nBadFaceWarpageFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceWarpageFaces;
    }

    label nBadFaceSkewnessFaces = 0;
    bool doSkewnessCheck = maxIntSkew > 0 || maxBounSkew > 0;
    if (doSkewnessCheck)
    {
        polyMeshGeometry::checkFaceSkewness
        (
            verbose,
            maxIntSkew,
            maxBounSkew,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces,
            meshStats
        );

        nBadFaceSkewnessFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceSkewnessFaces;
    }

    label nBadFaceWeightFaces = 0;
    bool doWeightCheck = minWeight >= 0 && minWeight < 1;
    if (doWeightCheck)
    {
        polyMeshGeometry::checkFaceWeights
        (
            verbose,
            minWeight,
            mesh,
            mesh.cellCentres(),
            mesh.faceCentres(),
            mesh.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces,
            meshStats
        );

        nBadFaceWeightFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceWeightFaces;
    }

    label nBadVolRatioFaces = 0;
    bool doVolRatioCheck = minVolRatio >= 0;
    if (minVolRatio >= 0)
    {
        polyMeshGeometry::checkVolRatio
        (
            verbose,
            minVolRatio,
            mesh,
            mesh.cellVolumes(),
            checkFaces,
            baffles,
            &wrongFaces,
            meshStats
        );

        nBadVolRatioFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadVolRatioFaces;
    }

    label nBadTwistFaces = 0;
    bool doTwistCheck = minTwist > -1;
    if (doTwistCheck)
    {
        //Pout<< "Checking face twist: dot product of face normal "
        //    << "with face triangle normals" << endl;
        polyMeshGeometry::checkFaceTwist
        (
            verbose,
            minTwist,
            mesh,
            mesh.cellCentres(),
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        nBadTwistFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadTwistFaces;
    }

    label nBadTriangleTwistFaces = 0;
    bool doTriangleTwistCheck = minTriangleTwist > -1;
    if (doTriangleTwistCheck)
    {
        //Pout<< "Checking triangle twist: dot product of consecutive triangle"
        //    << " normals resulting from face-centre decomposition" << endl;
        polyMeshGeometry::checkTriangleTwist
        (
            verbose,
            minTriangleTwist,
            mesh,
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        nBadTriangleTwistFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadTriangleTwistFaces;
    }

    label nBadFaceFlatnessFaces = 0;
    bool doFlatnessCheck = minFaceFlatness > -SMALL;
    if (doFlatnessCheck)
    {
        polyMeshGeometry::checkFaceFlatness
        (
            verbose,
            minFaceFlatness,
            mesh,
            mesh.faceAreas(),
            mesh.faceCentres(),
            mesh.points(),
            checkFaces,
            &wrongFaces
        );

        nBadFaceFlatnessFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceFlatnessFaces;
    }

    label nBadCellDeterminantFaces = 0;
    bool doCellDeterminantFaceCheck = minDet > -1;
    if (doCellDeterminantFaceCheck)
    {
        polyMeshGeometry::checkCellDeterminant
        (
            verbose,
            minDet,
            mesh,
            mesh.faceAreas(),
            checkFaces,
            polyMeshGeometry::affectedCells(mesh, checkFaces),
            &wrongFaces,
            meshStats
        );

        nBadCellDeterminantFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadCellDeterminantFaces;
    }

    label nBadSnapVolumeFaces = 0;
    bool doSnapVolumeCheck =    (minSnapRelativeVolume > SMALL  || minSnapRelativeTetVolume > SMALL)
                                && mesh.objectRegistry::foundObject<volScalarField>("preSnapVolume");
    if (doSnapVolumeCheck)
    {
        const volScalarField& initVol =
            mesh.lookupObject<volScalarField>("preSnapVolume");

        polyMeshGeometry::checkSnapVolume
        (
            verbose,
            minSnapRelativeVolume,
            minSnapRelativeTetVolume,
            mesh,
            mesh.cellVolumes(),
            initVol,
            checkFaces,
            polyMeshGeometry::affectedCells(mesh, checkFaces),
            &wrongFaces
        );

        nBadSnapVolumeFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadSnapVolumeFaces;
    }

    label nBadGaussGreenCentroidFaces = 0;
    bool doCentroidCheck = maxGaussGreenCentroid > SMALL
                           && mesh.objectRegistry::foundObject<volScalarField>("preSnapVolume");
    if (doCentroidCheck)
    {
        const volScalarField& initVol =
            mesh.lookupObject<volScalarField>("preSnapVolume");

        polyMeshGeometry::checkGaussGreenCentroid
        (
            verbose,
            maxGaussGreenCentroid,
            mesh,
            initVol,
            checkFaces,
            polyMeshGeometry::affectedCells(mesh, checkFaces),
            &wrongFaces
        );

        nBadGaussGreenCentroidFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadGaussGreenCentroidFaces;
    }

    label nBadCellAspectRatioFaces = 0;
    bool doAspectRatioCheck = maxCellAspectRatio > 0;
    if (doAspectRatioCheck)
    {
        polyMeshGeometry::checkCellAspectRatio
        (
            verbose,
            maxCellAspectRatio,
            mesh,
            mesh.faceAreas(),
            mesh.cellVolumes(),
            checkFaces,
            polyMeshGeometry::affectedCells(mesh, checkFaces),
            &wrongFaces,
            meshStats
        );

        nBadCellAspectRatioFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadCellAspectRatioFaces;
    }

    if (report) {
        // If the detailed report is needed, then reduce all the counters separately
        // so we can make a nice breakdown of what's wrong:
        reduce(
            std::tie(
                nBadDotProdFaces,
                nBadFacePyramidFaces,
                nBadFaceTetFaces,
                nBadFaceAngleFaces,
                nBadFaceAreaFaces,
                nBadEdgeLengthFaces,
                nBadFaceFaceCellsFaces,
                nBadFaceWarpageFaces,
                nBadFaceSkewnessFaces,
                nBadFaceWeightFaces,
                nBadVolRatioFaces,
                nBadTwistFaces,
                nBadTriangleTwistFaces,
                nBadFaceFlatnessFaces,
                nBadCellDeterminantFaces,
                nBadSnapVolumeFaces,
                nBadGaussGreenCentroidFaces,
                nBadCellAspectRatioFaces
            ),
            UniformParallelOp<sumOp<label>, 18>{}
        );

        if (doDotProductCheck) {
            Info<< "    non-orthogonality > "
                 << setw(3) << maxNonOrtho
                 << " degrees                        : "
                 << nBadDotProdFaces << endl;
        }
        if (doFacePyramidCheck) {
            Info<< "    faces with face pyramid volume < "
                 << setw(5) << minVol << "                 : "
                 << nBadFacePyramidFaces << endl;
        }
        if (doFaceTetCheck) {
            Info<< "    faces with face-decomposition tet quality < "
                 << setw(5) << minTetQuality << "      : "
                 << nBadFaceTetFaces << endl;
        }
        if (doFaceAngleCheck) {
            Info<< "    faces with concavity > "
                 << setw(3) << min(maxConcave, maxPyrConcave)
                 << " degrees                     : "
                 << nBadFaceAngleFaces << endl;
        }
        if (doFaceAreaCheck) {
            Info<< "    faces with area < "
                 << setw(5) << minArea
                 << " m^2                            : "
                 << nBadFaceAreaFaces << endl;
        }
        if (doFaceEdgeLengthCheck) {
            Info<< "    faces with edges < "
                 << setw(4) << minEdgeLength
                 << " m                             : "
                 << nBadEdgeLengthFaces << endl;
        }
        if (faceFaceCells) {
            Info<< "    face cells with faceFace connectivity issues  "
                 << "         : "
                 << nBadFaceFaceCellsFaces << endl;
        }
        if (doWarpageCheck) {
            Info<< "    faces with warpage > "
                 << setw(3) << maxIntWarp
                 << " (internal) or " << setw(3) << maxBounWarp
                 << " (boundary)  : " << nBadFaceWarpageFaces << endl;
        }
        if (doSkewnessCheck) {
            Info<< "    faces with skewness > "
                 << setw(3) << maxIntSkew
                 << " (internal) or " << setw(3) << maxBounSkew
                 << " (boundary) : " << nBadFaceSkewnessFaces << endl;
        }
        if (doWeightCheck) {
            Info<< "    faces with interpolation weights (0..1)  < "
                << setw(5) << minWeight
                << "       : "
                << nBadFaceWeightFaces << endl;
        }
        if (doVolRatioCheck) {
            Info<< "    faces with volume ratio of neighbour cells < "
                 << setw(5) << minVolRatio
                 << "     : "
                 << nBadVolRatioFaces << endl;
        }
        if (doTwistCheck) {
            Info<< "    faces with face twist < "
                 << setw(5) << minTwist
                 << "                          : "
                 << nBadTwistFaces << endl;
        }
        if (doTriangleTwistCheck) {
            Info<< "    faces with triangle twist < "
                 << setw(5) << minTriangleTwist
                 << "                      : "
                 << nBadTriangleTwistFaces << endl;
        }
        if (doFlatnessCheck) {
            Info<< "    faces with flatness < "
                 << setw(5) << minFaceFlatness
                 << "                      : "
                 << nBadFaceFlatnessFaces << endl;
        }
        if (doCellDeterminantFaceCheck) {
            Info<< "    faces on cells with determinant < "
                 << setw(5) << minDet << "               : "
                 << nBadCellDeterminantFaces << endl;
        }
        if (doSnapVolumeCheck) {
            Info<< "    faces on cells with snap volume ratio < "
                 << setw(3) << max(minSnapRelativeVolume,minSnapRelativeTetVolume)
                 << "            : "
                 << nBadSnapVolumeFaces << endl;
        }
        if (doCentroidCheck) {
            Info<< "    faces on cells with gauss-green centroid > "
                 << setw(4) << maxGaussGreenCentroid
                 << "        : "
                 << nBadGaussGreenCentroidFaces << endl;
        }
        if (doAspectRatioCheck) {
            Info<< "    faces on cells with aspect ratio > "
                 << setw(5) << maxCellAspectRatio << "               : "
                 << nBadCellAspectRatioFaces << endl;
        }
    }

    // If not reporting, this is going to just be a count of errors found by this node.
    // If we did reporting above, this will be all errors everywhere.
    nWrongFaces =
        nBadDotProdFaces +
        nBadFacePyramidFaces +
        nBadFaceTetFaces +
        nBadFaceAngleFaces +
        nBadFaceAreaFaces +
        nBadEdgeLengthFaces +
        nBadFaceFaceCellsFaces +
        nBadFaceWarpageFaces +
        nBadFaceSkewnessFaces +
        nBadFaceWeightFaces +
        nBadVolRatioFaces +
        nBadTwistFaces +
        nBadTriangleTwistFaces +
        nBadFaceFlatnessFaces +
        nBadCellDeterminantFaces +
        nBadSnapVolumeFaces +
        nBadGaussGreenCentroidFaces +
        nBadCellAspectRatioFaces;

    if (!report) {
        // If we didn't do the reporting, then we didn't distribute all the various
        // `n<SpecificProblem>Faces`, so the above calculation of `nWrongFaces` will
        // count only errors found on this node, so we must do a reduce of just that
        // one value to get the overall sum. This way, we only do the amount of reducing
        // we actually need: just this one integer if not reporting, but all of the
        // little counters if the breakdown is needed.
        reduce(nWrongFaces, sumOp<label>{});
    }

    //Pout.setf(ios_base::right);

    return nWrongFaces;
}


Foam::label Foam::motionSmootherAlgo::checkMesh
(
    const bool verbose,
    const polyMesh& mesh,
    const dictionary& dict,
    labelHashSet& wrongFaces,
    const bool report,
    PtrList<meshStatistics>* meshStats
)
{
    return checkMesh
    (
        verbose,
        mesh,
        dict,
        identity(mesh.nFaces()),
        wrongFaces,
        report,
        meshStats
    );
}

Foam::label Foam::motionSmootherAlgo::checkMesh
(
    const bool verbose,
    const dictionary& dict,
    const polyMeshGeometry& meshGeom,
    const pointField& points,
    const labelList& checkFaces,
    labelHashSet& wrongFaces,
    const bool report,
    PtrList<meshStatistics>* meshStats
)
{
    List<labelPair> emptyBaffles;

    return checkMesh
    (
        verbose,
        dict,
        meshGeom,
        points,
        checkFaces,
        emptyBaffles,
        wrongFaces,
        report,
        meshStats
     );
}


Foam::label Foam::motionSmootherAlgo::checkMesh
(
    const bool verbose,
    const dictionary& dict,
    const polyMeshGeometry& meshGeom,
    const pointField& points,
    const labelList& checkFaces,
    const List<labelPair>& baffles,
    labelHashSet& wrongFaces,
    const bool report,
    PtrList<meshStatistics>* meshStats
)
{
    const scalar maxNonOrtho
    (
        dict.lookupOrDefault<scalar>("maxNonOrtho", 180, true)
    );
    const scalar minVol
    (
        dict.lookupOrDefault<scalar>("minVol", -1E30, true)
    );
    const scalar minTetQuality
    (
        dict.lookupOrDefault<scalar>("minTetQuality", -1e30, true)
    );
    const scalar maxConcave
    (
        dict.lookupOrDefault<scalar>("maxConcave", 180, true)
    );
    const scalar maxPyrConcave
    (
        dict.lookupOrDefault<scalar>("maxPyrConcave", 180, true)
    );
    const scalar minArea
    (
        dict.lookupOrDefault<scalar>("minArea", -1e30, true)
    );
    const scalar minEdgeLength
    (
        dict.lookupOrDefault<scalar>("minEdgeLength", -1, true)
    );
    const scalar maxIntWarp
    (
        dict.lookupOrDefault<scalar>("maxInternalWarpage", -1, true)
    );
    const scalar maxBounWarp
    (
        dict.lookupOrDefault<scalar>("maxBoundaryWarpage", -1, true)
    );
    const scalar maxIntSkew
    (
        dict.lookupOrDefault<scalar>("maxInternalSkewness", -1, true)
    );
    const scalar maxBounSkew
    (
        dict.lookupOrDefault<scalar>("maxBoundarySkewness", -1, true)
    );
    const scalar minWeight
    (
        dict.lookupOrDefault<scalar>("minFaceWeight", -1, true)
    );
    const scalar minVolRatio
    (
        dict.lookupOrDefault<scalar>("minVolRatio", -1, true)
    );
    const scalar minTwist
    (
        dict.lookupOrDefault<scalar>("minTwist", -2, true)
    );
    const scalar minTriangleTwist
    (
        dict.lookupOrDefault<scalar>("minTriangleTwist", -2, true)
    );
    scalar minFaceFlatness = -1.0;
    dict.readIfPresent("minFaceFlatness", minFaceFlatness, true);
    const scalar minDet
    (
        dict.lookupOrDefault<scalar>("minDeterminant", -2, true)
    );
    const Switch faceFaceCells
    (
        dict.lookupOrDefault<Switch>("faceFaceCells", false)
    );
    const scalar maxFaceCentreNonOrtho
    (
        dict.lookupOrDefault<scalar>("maxFaceCentreNonOrtho", 180)
    );
    const scalar maxCellAspectRatio
    (
        dict.lookupOrDefault<scalar>("maxCellAspectRatio", -1)
    );

    label nWrongFaces = 0;

    if (report)
    {
        Info<< "Checking faces in error :" << endl;
        //Pout.setf(ios_base::left);
    }

    label nBadDotProdFaces = 0;
    bool doDotProdCheck = maxNonOrtho < 180.0-SMALL;
    if (doDotProdCheck)
    {
        meshGeom.checkFaceDotProduct
        (
            verbose,
            maxNonOrtho,
            maxFaceCentreNonOrtho,
            checkFaces,
            baffles,
            &wrongFaces
        );

        nBadDotProdFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadDotProdFaces;
    }

    label nBadFacePyramidFaces = 0;
    bool doFacePyramidCheck = minVol > -GREAT;
    if (doFacePyramidCheck)
    {
        meshGeom.checkFacePyramids
        (
            verbose,
            minVol,
            points,
            checkFaces,
            baffles,
            &wrongFaces
        );

        nBadFacePyramidFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFacePyramidFaces;
    }

    label nBadFaceTetFaces = 0;
    bool doTetQualityCheck = minTetQuality > -GREAT;
    if (doTetQualityCheck)
    {
        meshGeom.checkFaceTets
        (
            verbose,
            minTetQuality,
            points,
            checkFaces,
            baffles,
            &wrongFaces
        );

        nBadFaceTetFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceTetFaces;
    }

    label nBadFaceAngleFaces = 0;
    bool doFaceAngleCheck = maxConcave < 180.0-SMALL;
    if (doFaceAngleCheck)
    {
        meshGeom.checkFaceAngles
        (
            verbose,
            maxConcave,
            maxPyrConcave,
            points,
            checkFaces,
            &wrongFaces
        );

        nBadFaceAngleFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceAngleFaces;
    }

    label nBadFaceAreaFaces = 0;
    bool doFaceAreaCheck = minArea > -SMALL;
    if (doFaceAreaCheck)
    {
        meshGeom.checkFaceArea
        (
            verbose,
            minArea,
            checkFaces,
            &wrongFaces
        );

        nBadFaceAreaFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceAreaFaces;
    }

    label nBadEdgeLengthFaces = 0;
    bool doEdgeLengthCheck = minEdgeLength > 0;
    if (doEdgeLengthCheck)
    {
        meshGeom.checkFaceEdgeLengths
        (
            verbose,
            minEdgeLength,
            checkFaces,
            &wrongFaces
        );

        nBadEdgeLengthFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadEdgeLengthFaces;
    }

    label nBadFaceWarpageFaces = 0;
    bool doWarpageCheck = maxIntWarp > 0 || maxBounWarp > 0;
    if (doWarpageCheck)
    {
        meshGeom.checkFaceWarpage
        (
            verbose,
            maxIntWarp,
            maxBounWarp,
            meshGeom.mesh().points(),
            checkFaces,
            &wrongFaces
        );

        nBadFaceWarpageFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceWarpageFaces;
    }

    label nBadFaceFaceCellsFaces = 0;
    if (faceFaceCells)
    {
        meshGeom.checkFaceFaceCells
        (
            verbose,
            checkFaces,
            &wrongFaces
        );

        nBadFaceFaceCellsFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceFaceCellsFaces;
    }

    label nBadFaceSkewnessFaces = 0;
    bool doSkewnessCheck = maxIntSkew > 0 || maxBounSkew > 0;
    if (doSkewnessCheck)
    {
        polyMeshGeometry::checkFaceSkewness
        (
            verbose,
            maxIntSkew,
            maxBounSkew,
            meshGeom.mesh(),
            meshGeom.cellCentres(),
            meshGeom.faceCentres(),
            meshGeom.faceAreas(),
            checkFaces,
            baffles,
            &wrongFaces
        );

        nBadFaceSkewnessFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceSkewnessFaces;
    }

    label nBadFaceWeightFaces = 0;
    bool doWeightCheck = minWeight >= 0 && minWeight < 1;
    if (doWeightCheck)
    {
        meshGeom.checkFaceWeights
        (
            verbose,
            minWeight,
            checkFaces,
            baffles,
            &wrongFaces
        );

        nBadFaceWeightFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceWeightFaces;
    }

    label nBadVolRatioFaces = 0;
    bool doVolRatioCheck = minVolRatio >= 0;
    if (doVolRatioCheck)
    {
        meshGeom.checkVolRatio
        (
            verbose,
            minVolRatio,
            checkFaces,
            baffles,
            &wrongFaces
        );

        nBadVolRatioFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadVolRatioFaces;
    }

    label nBadTwistFaces = 0;
    bool doTwistCheck = minTwist > -1;
    if (doTwistCheck)
    {
        //Pout<< "Checking face twist: dot product of face normal "
        //    << "with face triangle normals" << endl;
        meshGeom.checkFaceTwist
        (
            verbose,
            minTwist,
            points,
            checkFaces,
            &wrongFaces
        );

        nBadTwistFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadTwistFaces;
    }

    label nBadTriangleTwistFaces = 0;
    bool doTriangleTwistCheck = minTriangleTwist > -1;
    if (doTriangleTwistCheck)
    {
        //Pout<< "Checking triangle twist: dot product of consecutive triangle"
        //    << " normals resulting from face-centre decomposition" << endl;
        meshGeom.checkTriangleTwist
        (
            verbose,
            minTriangleTwist,
            points,
            checkFaces,
            &wrongFaces
        );

        nBadTriangleTwistFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadTriangleTwistFaces;
    }

    label nBadFaceFlatnessFaces = 0;
    bool doFlatnessCheck = minFaceFlatness > -SMALL;
    if (doFlatnessCheck)
    {
        meshGeom.checkFaceFlatness
        (
            verbose,
            minFaceFlatness,
            points,
            checkFaces,
            &wrongFaces
        );

        nBadFaceFlatnessFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadFaceFlatnessFaces;
    }

    label nBadCellDeterminantFaces = 0;
    bool doDeterminantCheck = minDet > -1;
    if (doDeterminantCheck)
    {
        meshGeom.checkCellDeterminant
        (
            verbose,
            minDet,
            checkFaces,
            polyMeshGeometry::affectedCells(meshGeom.mesh(), checkFaces),
            &wrongFaces
        );

        nBadCellDeterminantFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadCellDeterminantFaces;
    }

    label nBadCellAspectRatioFaces = 0;
    bool doAspectRatioCheck = maxCellAspectRatio > 0;
    if (doAspectRatioCheck)
    {
        meshGeom.checkCellAspectRatio
        (
            verbose,
            maxCellAspectRatio,
            checkFaces,
            polyMeshGeometry::affectedCells(meshGeom.mesh(), checkFaces),
            &wrongFaces
        );

        nBadCellAspectRatioFaces = wrongFaces.size() - nWrongFaces;
        nWrongFaces += nBadCellAspectRatioFaces;
    }

    // TODO: Uh, this function seems to be the same as the other checkMesh function
    //       above, except this one skips snapVolumeFaces and GaussGreenCentroidFaces.
    //       ... Consider using a constexpr bool non-type-template-parameter to redefine
    //       as one function with zero runtime overhead. A kitten dies evey time you
    //       copy-paste a 500 line function.

    if (report) {
        // If the detailed report is needed, then reduce all the counters separately
        // so we can make a nice breakdown of what's wrong:
        reduce(
            std::tie(
                nBadDotProdFaces,
                nBadFacePyramidFaces,
                nBadFaceTetFaces,
                nBadFaceAngleFaces,
                nBadFaceAreaFaces,
                nBadEdgeLengthFaces,
                nBadFaceFaceCellsFaces,
                nBadFaceWarpageFaces,
                nBadFaceSkewnessFaces,
                nBadFaceWeightFaces,
                nBadVolRatioFaces,
                nBadTwistFaces,
                nBadTriangleTwistFaces,
                nBadFaceFlatnessFaces,
                nBadCellDeterminantFaces,
                nBadCellAspectRatioFaces
            ),
            UniformParallelOp<sumOp<label>, 16>{}
        );

        if (doDotProdCheck) {
            Info<< "    non-orthogonality > "
                 << setw(3) << maxNonOrtho
                 << " degrees                        : "
                 << nBadDotProdFaces << endl;
        }
        if (doFacePyramidCheck) {
            Info<< "    faces with face pyramid volume < "
                 << setw(5) << minVol << "                 : "
                 << nBadFacePyramidFaces << endl;
        }
        if (doTetQualityCheck) {
            Info<< "    faces with face-decomposition tet quality < "
                 << setw(5) << minTetQuality << "      : "
                 << nBadFaceTetFaces << endl;
        }
        if (doFaceAngleCheck) {
            Info<< "    faces with concavity > "
                 << setw(3) << min(maxConcave, maxPyrConcave)
                 << " degrees                     : "
                 << nBadFaceAngleFaces << endl;
        }
        if (doFaceAreaCheck) {
            Info<< "    faces with area < "
                 << setw(5) << minArea
                 << " m^2                            : "
                 << nBadFaceAreaFaces << endl;
        }
        if (doEdgeLengthCheck) {
            Info<< "    faces with edges < "
                 << setw(4) << minEdgeLength
                 << " m                             : "
                 << nBadEdgeLengthFaces << endl;
        }
        if (doWarpageCheck) {
            Info<< "    faces with warpage > "
                 << setw(3) << maxIntWarp
                 << " (internal) or " << setw(3) << maxBounWarp
                 << " (boundary) : " << nBadFaceWarpageFaces << endl;
        }
        if (faceFaceCells) {
            Info<< "    face cells with faceFace connectivity issues  "
                 << "         : "
                 << nBadFaceFaceCellsFaces << endl;
        }
        if (doSkewnessCheck) {
            Info<< "    faces with skewness > "
                 << setw(3) << maxIntSkew
                 << " (internal) or " << setw(3) << maxBounSkew
                 << " (boundary) : " << nBadFaceSkewnessFaces << endl;
        }
        if (doWeightCheck) {
            Info<< "    faces with interpolation weights (0..1)  < "
                 << setw(5) << minWeight
                 << "       : "
                 << nBadFaceWeightFaces << endl;
        }
        if (doVolRatioCheck) {
            Info<< "    faces with volume ratio of neighbour cells < "
                 << setw(5) << minVolRatio
                 << "     : "
                 << nBadVolRatioFaces << endl;
        }
        if (doTwistCheck) {
            Info<< "    faces with face twist < "
                 << setw(5) << minTwist
                 << "                          : "
                 << nBadTwistFaces << endl;
        }
         if (doTriangleTwistCheck) {
             Info<< "    faces with triangle twist < "
                  << setw(5) << minTriangleTwist
                  << "                      : "
                  << nBadTriangleTwistFaces << endl;
         }
         if (doFlatnessCheck) {
             Info<< "    faces with flatness < "
                  << setw(5) << minFaceFlatness
                  << "                      : "
                  << nBadFaceFlatnessFaces << endl;
         }
         if (doDeterminantCheck) {
             Info<< "    faces on cells with determinant < "
                  << setw(5) << minDet << "                : "
                  << nBadCellDeterminantFaces << endl;
         }
         if (doAspectRatioCheck) {
             Info<< "    faces on cells with aspect ratio > "
                  << setw(5) << maxCellAspectRatio << "               : "
                  << nBadCellAspectRatioFaces << endl;
         }
    }

    // If not reporting, this is going to just be a count of errors found by this node.
    // If we did reporting above, this will be all errors everywhere.
    nWrongFaces =
        nBadDotProdFaces +
        nBadFacePyramidFaces +
        nBadFaceTetFaces +
        nBadFaceAngleFaces +
        nBadFaceAreaFaces +
        nBadEdgeLengthFaces +
        nBadFaceFaceCellsFaces +
        nBadFaceWarpageFaces +
        nBadFaceSkewnessFaces +
        nBadFaceWeightFaces +
        nBadVolRatioFaces +
        nBadTwistFaces +
        nBadTriangleTwistFaces +
        nBadFaceFlatnessFaces +
        nBadCellDeterminantFaces +
        nBadCellAspectRatioFaces;

    if (!report) {
        // If we didn't do the reporting, then we didn't distribute all the various
        // `n<SpecificProblem>Faces`, so the above calculation of `nWrongFaces` will
        // count only errors found on this node, so we must do a reduce of just that
        // one value to get the overall sum. This way, we only do the amount of reducing
        // we actually need: just this one integer if not reporting, but all of the
        // little counters if the breakdown is needed.
        reduce(nWrongFaces, sumOp<label>{});
    }

    //Pout.setf(ios_base::right);

    return nWrongFaces;
}


// ************************************************************************* //
