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
    (c) 2015 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

Class
    Foam::medialAxisMeshMover

Description
    Mesh motion solver that uses a medial axis algorithm to work
    out a fraction between the (nearest point on a) moving surface
    and the (nearest point on a) fixed surface.

    This fraction is then used to scale the motion. Use
    - fixedValue on all moving patches
    - zeroFixedValue on stationary patches
    - slip on all slipping patches
    Note that the fixedValue boundary conditions might be changed by this
    solver to enforce consistency and a valid resulting mesh.

SourceFiles
    medialAxisMeshMover.C

\*---------------------------------------------------------------------------*/

#ifndef medialAxisMeshMover_H
#define medialAxisMeshMover_H

#include "externalDisplacementMeshMover/externalDisplacementMeshMover.H"
#include "motionSmoother/motionSmootherAlgo.H"
#include "snappyHexMeshDriver/snappyLayerDriver.H"
#include "externalDisplacementMeshMover/fieldSmoother/fieldSmoother.H"
#include "snappyHexMeshDriver/pointData/pointData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template <class DataType>
class PointData;

/*---------------------------------------------------------------------------*\
             Class medialAxisMeshMover Declaration
\*---------------------------------------------------------------------------*/

class medialAxisMeshMover
:
    public externalDisplacementMeshMover
{

    //- Combine operator class for equalizing displacements.
        class minMagEqOp
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (magSqr(y) < magSqr(x))
                {
                    x = y;
                }
            }
        };

    //- Combine operator class to combine normal with other normal.
        class nomalsCombine
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (y != vector(GREAT, GREAT, GREAT))
                {
                    if (x == vector(GREAT, GREAT, GREAT))
                    {
                        x = y;
                    }
                    else
                    {
                        x += y;
                        x /= mag(x) + SMALL;
                    }
                }
            }
        };

    //- Find a unique edge direction.
        struct uniqueEdgeNormal
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (x != vector(GREAT, GREAT, GREAT))
                {
                    return;
                }
                else if (y != vector(GREAT, GREAT, GREAT))
                {
                    x = y;
                }
            }
        };

    // Private data

        const labelList slipPatchIDs_;

        const labelList adaptPatchIDs_;

        autoPtr<indirectPrimitivePatch> adaptPatchPtr_;

        //- Scale factor for displacement
        pointScalarField scale_;

        //- Starting mesh position
        pointField oldPoints_;

        //- Mesh mover algorithm
        motionSmootherAlgo meshMover_;

        //- Field smoothing
        fieldSmoother fieldSmoother_;


    // Pre-calculated medial axis information

        //- Normal of nearest wall. Where this normal changes direction
        //  defines the medial axis
        pointVectorField dispVec_;

        //- Ratio of medial distance to wall distance.
        //  (1 at wall, 0 at medial axis)
        pointScalarField medialRatio_;

        //- Distance to nearest medial axis point
        pointScalarField medialDist_;

        //- Location on nearest medial axis point
        pointVectorField medialVec_;


    // Private Member Functions

        static labelList getSlipBCs(const pointVectorField&);

        //- Extract fixed-value patchfields
        //- Extract bc types. Replace fixedValue derivatives with fixedValue
        wordList getPatchFieldTypes(const pointVectorField& fld);


        // Calculation of medial axis information

            //- Is mesh edge on a cusp of displacement
            bool isMaxEdge
            (
                const List<pointData>& pointWallDist,
                const label edgeI,
                const scalar minCos
            ) const;

            //- Initialise medial axis. Uses pointDisplacement field only
            //  for boundary types, not values.
            void update(const dictionary&);


        // Calculation of mesh movement

            //- Unmark extrusion at given point
            static bool unmarkExtrusion
            (
                const label patchPointi,
                pointField& patchDisp,
                List<snappyLayerDriver::extrudeMode>& extrudeStatus
            );

            //- Synchronise extrusion
            void syncPatchDisplacement
            (
                const scalarField& minThickness,
                pointField& patchDisp,
                List<snappyLayerDriver::extrudeMode>& extrudeStatus
            ) const;

            //- Stop layer growth at feature edges
            void handleFeatureAngleLayerTerminations
            (
                const scalar minCos,
                const labelList& meshEdges,
                const PackedList<1>& isConvexEdgePoint,
                const PackedList<1>& isConcaveEdgePoint,
                const vectorField& patchEdgeNormals,
                List<snappyLayerDriver::extrudeMode>& extrudeStatus,
                pointField& patchDisp,
                label& nPointCounter
            ) const;

            //- Find isolated islands (points, edges and faces and layer
            //  terminations) in the layer mesh and stop any layer growth
            //  at these points
            void findIsolatedRegions
            (
                const scalar minCosLayerTermination,
                const bool detectExtrusionIsland,
                const PackedBoolList& isMasterEdge,
                const PackedList<1>& isConvexEdgePoint,
                const PackedList<1>& isConcaveEdgePoint,
                const labelList& meshEdges,
                const scalarField& minThickness,
                const vectorField& patchEdgeNormals,
                List<snappyLayerDriver::extrudeMode>& extrudeStatus,
                pointField& patchDisp
            ) const;


        //- Calculate desired displacement. Modifies per-patch displacement
        //  and calculates displacement for whole mesh
        //  (in pointDisplacement)
        void calculateDisplacement
        (
            const dictionary&,
            const scalarField& minThickness,
            const PackedList<1>& isConvexEdgePoint,
            const PackedList<1>& isConcaveEdgePoint,
            List<snappyLayerDriver::extrudeMode>& extrudeStatus,
            pointField& patchDisp
        );

        //- Move mesh according to calculated displacement
        bool shrinkMesh
        (
            const dictionary& meshQualityDict,
            const label nAllowableErrors,
            labelList& checkFaces
        );

        //- Disallow default bitwise copy construct
        medialAxisMeshMover(const medialAxisMeshMover&);

        //- Disallow default bitwise assignment
        void operator=(const medialAxisMeshMover&);


public:

    //- Runtime type information
    TypeName("displacementMedialAxis");


    // Constructors

        //- Construct from dictionary and displacement field
        medialAxisMeshMover
        (
            const dictionary& dict,
            const List<labelPair>& baffles,
            pointVectorField& pointDisplacement
        );


    // Destructor

        virtual ~medialAxisMeshMover();


    // Member Functions

        //- Move mesh using current pointDisplacement boundary values.
        //  Return true if succesful (errors on checkFaces less than
        //  allowable). Updates pointDisplacement.
        virtual bool move
        (
            const dictionary&,
            const label nAllowableErrors,
            const PackedList<1>& isConvexEdgePoint,
            const PackedList<1>& isConcaveEdgePoint,
            labelList& checkFaces
        );

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&)
        {
            NotImplemented;
        }

        //- Return faceCentre, faceArea, cellCentre and cellVolume
        //  fields from updated point field
        static void updateMeshGeometric
        (
            const polyMesh& mesh,
            const pointField& newPoints,
            pointField& newFaceCentres,
            vectorField& newFaceAreas,
            pointField& newCellCentres,
            scalarField& newCellVolumes,
            vectorField& cEst,
            labelField& nCellFaces
        );

        //- General mesh smoothing in the interior
        static void smoothDisplacement
        (
            const dictionary& motionDict,
            const labelList& planarPatches,
            const scalar& edge0Len,
            const polyMesh& mesh,

            pointField newPoints,
            vectorField& displacement
         );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
