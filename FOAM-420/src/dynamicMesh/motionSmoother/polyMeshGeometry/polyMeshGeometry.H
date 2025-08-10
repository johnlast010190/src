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

Class
    Foam::polyMeshGeometry

Description
    Updateable mesh geometry and checking routines.

    - non-ortho done across coupled faces.
    - faceWeight (delta factors) done across coupled faces.

SourceFiles
    polyMeshGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGeometry_H
#define polyMeshGeometry_H

#include "fields/GeometricFields/pointFields/pointFields.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "fields/volFields/volFields.H"
#include "primitives/Tuple2/Tuple2.H"
#include "motionSmoother/polyMeshGeometry/meshStatistics.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class polyMeshGeometry Declaration
\*---------------------------------------------------------------------------*/

class polyMeshGeometry
{
        //- Reference to polyMesh.
        const polyMesh& mesh_;

        //- Uptodate copy of face areas
        vectorField faceAreas_;

        //- Uptodate copy of face centres
        vectorField faceCentres_;

        //- Uptodate copy of cell centres
        vectorField cellCentres_;

        //- Uptodate copy of cell volumes
        scalarField cellVolumes_;


    // Private Member Functions

        //Check if cell is pyramid (check for collapsed edges)
        static bool isPyrCell(const polyMesh& mesh, const label cellI);

        //- Update face areas and centres on selected faces.
        void updateFaceCentresAndAreas
        (
            const pointField& p,
            const labelList& changedFaces
        );

        //- Update cell volumes and centres on selected cells. Requires
        //  cells and faces to be consistent set.
        void updateCellCentresAndVols
        (
            const labelList& changedCells,
            const labelList& changedFaces
        );

        //- Detect&report non-ortho error for single face.
        static scalar checkNonOrtho
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
        );

        //- Calculate skewness given two cell centres and one face centre.
        static scalar calcSkewness
        (
            const polyMesh& mesh,
            const point& ownCc,
            const point& neiCc,
            const point& fc,
            const vector& fa,
            const label& faceI
        );

        //- Detect&report incorrect face-triangle orientation
        static bool checkFaceTet
        (
            const polyMesh&,
            const bool verbose,
            const scalar minTetQuality,
            const pointField& p,
            const label facei,
            const point& fc,    // face centre
            const point& cc,    // cell centre
            labelHashSet* setPtr
        );


public:

    ClassName("polyMeshGeometry");

    // Constructors

        //- Construct from mesh
        polyMeshGeometry(const polyMesh&);


    // Member Functions

        // Access

            //- Find index of metric field
            static label findStatIndex
            (
                const PtrList<meshStatistics>* metricStats,
                const word& metricName
            );

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            const vectorField& faceAreas() const
            {
                return faceAreas_;
            }
            const vectorField& faceCentres() const
            {
                return faceCentres_;
            }
            const vectorField& cellCentres() const
            {
                return cellCentres_;
            }
            const scalarField& cellVolumes() const
            {
                return cellVolumes_;
            }

        // Edit

            //- Take over properties from mesh
            void correct();

            //- Recalculate on selected faces. Recalculates cell properties
            //  on owner and neighbour of these cells.
            void correct
            (
                const pointField& p,
                const labelList& changedFaces
            );

            //- Helper function: get affected cells from faces
            static labelList affectedCells
            (
                const polyMesh&,
                const labelList& changedFaces
            );


        // Checking of selected faces with supplied geometry (mesh only used for
        // topology). Coupled aware (coupled faces treated as internal ones)

            //- See primitiveMesh
            static bool checkFaceDotProduct
            (
                const bool verbose,
                const scalar orthWarn,
                const scalar maxFaceCentreNonOrtho,
                const polyMesh&,
                const vectorField& faceCentres,
                const vectorField& cellCentres,
                const vectorField& faceAreas,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet* setPtr,
                PtrList<meshStatistics>* meshStats = nullptr
            );

            //- See primitiveMesh
            static bool checkFacePyramids
            (
                const bool verbose,
                const scalar minPyrVol,
                const polyMesh&,
                const pointField& faceCentres,
                const pointField& cellCentres,
                const vectorField& faceAreas,
                const pointField& p,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet*,
                PtrList<meshStatistics>* meshStats = nullptr
            );

            //- See primitiveMesh
            static bool checkFaceTets
            (
                const bool verbose,
                const scalar minPyrVol,
                const polyMesh&,
                const vectorField& cellCentres,
                const vectorField& faceCentres,
                const pointField& p,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet*
            );

            //- See primitiveMesh
            static bool checkFaceWarpage
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
                PtrList<meshStatistics>* meshStats = nullptr
            );

            //- See primitiveMesh
            static bool checkFaceSkewness
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
                PtrList<meshStatistics>* meshStats = nullptr
            );

            //- Interpolation weights (0.5 for regular mesh)
            static bool checkFaceWeights
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
                PtrList<meshStatistics>* meshStats = nullptr
            );

            //- Cell volume ratio of neighbouring cells (1 for regular mesh)
            static bool checkVolRatio
            (
                const bool verbose,
                const scalar warnRatio,
                const polyMesh& mesh,
                const scalarField& cellVolumes,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet* setPtr,
                PtrList<meshStatistics>* meshStats = nullptr
            );

            //- See primitiveMesh
            static bool checkFaceAngles
            (
                const bool verbose,
                const scalar maxDeg,
                const scalar maxPyrDeg,
                const polyMesh& mesh,
                const vectorField& faceAreas,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            //- Triangle (from face-centre decomposition) normal v.s.
            //  average face normal
            static bool checkFaceTwist
            (
                const bool verbose,
                const scalar minTwist,
                const polyMesh&,
                const vectorField& cellCentres,
                const vectorField& faceAreas,
                const vectorField& faceCentres,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            //- Consecutive triangle (from face-centre decomposition) normals
            static bool checkTriangleTwist
            (
                const bool verbose,
                const scalar minTwist,
                const polyMesh&,
                const vectorField& faceAreas,
                const vectorField& faceCentres,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            //- Area of faces v.s. sum of triangle areas
            static bool checkFaceFlatness
            (
                const bool verbose,
                const scalar minFlatness,
                const polyMesh&,
                const vectorField& faceAreas,
                const vectorField& faceCentres,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            //- Small faces
            static bool checkFaceArea
            (
                const bool verbose,
                const scalar minArea,
                const polyMesh&,
                const vectorField& faceAreas,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            //- Small edges
            static bool checkFaceEdgeLengths
            (
                const bool verbose,
                const scalar minEdgeLength,
                const polyMesh&,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            //checking topological cell-faceFace usage
            static bool checkFaceFaceCells
            (
                const bool verbose,
                const polyMesh&,
                const labelList& checkFaces,
                labelHashSet* setPtr
            );

            //- Area of internal faces v.s. boundary faces
            static bool checkCellDeterminant
            (
                const bool verbose,
                const scalar minDet,
                const polyMesh&,
                const vectorField& faceAreas,
                const labelList& checkFaces,
                const labelList& affectedCells,
                labelHashSet* setPtr,
                PtrList<meshStatistics>* meshStats = nullptr
            );

            //- Specialised check using undistorted cell volumes
            static bool checkSnapVolume
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
             );

            //- specialised check of gauss-green centroid
            static bool checkGaussGreenCentroid
            (
                const bool verbose,
                const scalar maxGaussGreenCentroid,
                const polyMesh& mesh,
                const scalarField& initCellVolumes,
                const labelList& checkFaces,
                const labelList& affectedCells,
                labelHashSet* setPtr
             );

            //- specialised check of cell aspect ratio
            static bool checkCellAspectRatio
            (
                const bool verbose,
                const scalar warnAR,
                const polyMesh& mesh,
                const vectorField& faceAreas,
                const scalarField& cellVolumes,
                const labelList& checkFaces,
                const labelList& affectedCells,
                labelHashSet* setPtr,
                PtrList<meshStatistics>* meshStats = nullptr
            );

        // Checking of selected faces with local geometry. Uses above static
        // functions. Parallel aware.

            bool checkFaceDotProduct
            (
                const bool verbose,
                const scalar orthWarn,
                const scalar maxFaceCentreNonOrtho,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet* setPtr
            ) const;

            bool checkFacePyramids
            (
                const bool verbose,
                const scalar minPyrVol,
                const pointField& p,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet* setPtr
            ) const;

            bool checkFaceTets
            (
                const bool verbose,
                const scalar minTetQuality,
                const pointField& p,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet* setPtr
            ) const;

            bool checkFaceSkewness
            (
                const bool verbose,
                const scalar internalSkew,
                const scalar boundarySkew,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet* setPtr
            ) const;

            bool checkFaceWarpage
            (
                const bool verbose,
                const scalar internalWarpage,
                const scalar boundaryWarpage,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceWeights
            (
                const bool verbose,
                const scalar warnWeight,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet* setPtr
            ) const;

            bool checkVolRatio
            (
                const bool verbose,
                const scalar warnRatio,
                const labelList& checkFaces,
                const List<labelPair>& baffles,
                labelHashSet* setPtr
            ) const;

            bool checkFaceAngles
            (
                const bool verbose,
                const scalar maxDeg,
                const scalar maxPyrDeg,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceTwist
            (
                const bool verbose,
                const scalar minTwist,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkTriangleTwist
            (
                const bool verbose,
                const scalar minTwist,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceFlatness
            (
                const bool verbose,
                const scalar minFlatness,
                const pointField& p,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceArea
            (
                const bool verbose,
                const scalar minArea,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceEdgeLengths
            (
                const bool verbose,
                const scalar minEdgeLength,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkFaceFaceCells
            (
                const bool verbose,
                const labelList& checkFaces,
                labelHashSet* setPtr
            ) const;

            bool checkCellDeterminant
            (
                const bool verbose,
                const scalar warnDet,
                const labelList& checkFaces,
                const labelList& affectedCells,
                labelHashSet* setPtr
            ) const;

            bool checkCellAspectRatio
            (
                const bool verbose,
                const scalar warnAR,
                const labelList& checkFaces,
                const labelList& affectedCells,
                labelHashSet* setPtr
            ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
