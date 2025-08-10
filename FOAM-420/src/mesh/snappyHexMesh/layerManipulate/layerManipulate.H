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
    (c) 2014-1016 Esi Ltd.

Class
    Foam::layerManipulate

Description
    Perform manipulation to layers

SourceFiles
    layerManipulate.C

\*---------------------------------------------------------------------------*/

#ifndef layerManipulate_H
#define layerManipulate_H

#include "meshRefinement/meshRefinement.H"
#include "snappyHexMeshDriver/pointData/pointData.H"
#include "snappyHexMeshDriver/layerParameters/layerParameters.H"
#include "hessianMeshOptimization/layerCellHandling/layerPointStack/layerPointStack.H"
#include "meshControl/meshControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class layerManipulate Declaration
\*---------------------------------------------------------------------------*/

class layerManipulate
{
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

    class appendUniqueElements
    {
        public:

        void operator()(labelList& x, const labelList& y) const
        {
            if (y.size())
            {
                if (x.size())
                {
                    forAll(y, elemY)
                    {
                        label counter = 0;
                        forAll(x, elemX)
                        {
                            if (y[elemY] != x[elemX])
                            {
                                counter++;
                            }
                            if (counter == x.size())
                            {
                                x.append(y[elemY]);
                            }
                        }
                    }
                }
                else
                {
                    x = y;
                }
            }
        }
    };

    class mergePointStacks
    {
        public:

        void operator()(List<point>& x, const List<point>& y) const
        {
            label start = -1;
            label end = -1;
            if (x.size())
            {
                forAll(y, pI)
                {
                    if (x[0] == y[pI])
                    {
                        start = pI;
                    }

                    if (x[x.size()-1] == y[pI])
                    {
                        end = pI;
                    }
                }
            }
            else
            {
                x = y;
            }

            if (start != -1 && end != -1)
            {
                x = y;
            }
            else if (start > 0 && end == -1)
            {
                List<vector> newStackPoint(x.size()+start);
                for (int i = 0; i<start; i++)
                {
                    newStackPoint[i] = y[i];
                }
                for (int i = start; i<newStackPoint.size(); i++)
                {
                    int j = i - start;
                    newStackPoint[i] = x[j];
                }
                x = newStackPoint;
            }
            else if (start == -1 && end > -1)
            {
                List<vector> newStackPoint(x.size()+y.size()-end-1);
                for (int i = 0; i<x.size(); i++)
                {
                    newStackPoint[i] = x[i];
                }
                for (int i = x.size(); i<newStackPoint.size(); i++)
                {
                    int j = i - x.size() + end + 1;
                    newStackPoint[i] = y[j];
                }
                x = newStackPoint;
            }
        }
    };

private:

    // Private data

        //- Mesh+surface
        fvMesh& mesh_;

        //- Container for layer information
        const layerParameters& layerParams_;

        //- cellLevel
        const labelList& cellLevel_;

        //- pointLevel
        const labelList& pointLevel_;

        //- edge0 length
        const scalar edge0Length_;

public:

    //- Runtime type information
    TypeName("layerManipulate");

    // Constructors

        //- Construct from components
        layerManipulate
        (
            fvMesh& mesh,
            const layerParameters& layerParams,
            const labelList& cellLevel,
            const labelList& pointLevel,
            const scalar edge0Length
        );

    // Destructor

        virtual ~layerManipulate(){};

    // Member Functions

        autoPtr<indirectPrimitivePatch> makeFaceZonePatch();

        autoPtr<indirectPrimitivePatch> makeLayerPatch();

        autoPtr<indirectPrimitivePatch> makeGrownUpPatch
        (
            const bool onlySnap = false
        );

        static bool isMaxEdge
        (
            const fvMesh& mesh,
            const List<pointData>& pointWallDist,
            const label edgeI,
            const scalar minCos
        );

        //Update mesh geometric quantities
        void updateGeometry
        (
            const pointField& newPoints,
            pointField& newFaceCentres,
            vectorField& newFaceAreas,
            pointField& newCellCentres,
            scalarField& newCellVolumes
        );

        void deactivateSquishCells
        (
            const labelList& layerPointType,
            const labelList& boundaryPoints,
            const pointField& newFaceCentres,
            const vectorField& newFaceAreas,
            const pointField& newCellCentres,
            labelList& stationaryPts
        );

        void calculateLayerStacks
        (
            const indirectPrimitivePatch& pp,
            const boolList& boundaryFaces,
            const boolList& boundaryEdges,
            labelList& cellLayerNumber,
            labelList& layerFaceType,
            labelList& layerPointType,
            labelList& stackEdgeLayer
        );

        void writeLayerInfo();

        // Volume smooth layer interface and non-layer points
        void volumeSmoothPts
        (
            const meshControl& controller,
            const PackedBoolList& isMasterEdge,
            const boolList& hangingNodes,
            const List<labelList>& hangingNbrs,
            const labelList& layerPointType,
            const labelList& ftrPointOrigin,
            const labelList& boundaryPoints,
            const labelList& stationaryPts,
            const scalarField& layerWeights,
            const scalarField& nonLayerWeights,
            const boolList& stackEdges,
            const boolList& featureEdges,
            const scalarField& newCellVolumes,
            const vectorField& newFaceAreas,
            const pointField& newFaceCentres,
            const pointField& newCellCentres,
            const pointField& pointOrigin,
            const vectorField& pointSurfNormal,
            const pointField& tPts,
            pointField& newPoints,
            labelList& resetPts
        );

        // apply stretching to layer points
        void stretchLayerMesh
        (
            const meshControl& controller,
            const labelList& layerPointType,
            const boolList& boundaryEdges,
            const boolList& stackEdges,
            const labelList& boundaryPoints,
            const labelList& stationaryPts,
            const indirectPrimitivePatch& pp,
            const labelList& meshPoints,
            const List<Tuple2<word,scalar>>& layerMethod,
            const labelList& totalNumLayers,
            const labelList& layerCount,
            const vectorField& pointSurfNormal,
            const pointField& pointOrigin,
            const pointField& unsmoothedPointNormals,
            const pointField& outerPos,
            scalarField& initExpansion,
            pointField& newPoints,
            pointField& tPts
        );

        // reset and limit outer shell pts and project down to surface
        void reProjectOuter
        (
            const labelList& ftrPointOrigin,
            const labelList& layerPointType,
            const scalarField& maxLayerThickness,
            const boolList& boundaryEdges,
            const boolList& stackEdges,
            const vectorField& pointSurfNormal,
            const pointField& pointOrigin,
            const pointField& newPoints,
            pointField& outerPos
        );

        // optional resnap to grown up and zoned points
        void reSnap
        (
            const indirectPrimitivePatch& grownUpSnapPP,
            const indirectPrimitivePatch& fzonePP,
            const autoPtr<searchableSurfaces>& grownUpGeometryPtr,
            const autoPtr<searchableSurfaces>& grownUpZoneGeometryPtr,
            const labelList& boundaryPoints,
            const labelList& stationaryPts,
            const vectorField& newFaceAreas,
            const pointField& newFaceCentres,
            pointField& tPts,
            pointField& newPoints
        );

        // Try and correct concave outer shell points
        void correctConcaveCells
        (
            const labelList& layerPointType,
            const labelList& layerFaceType,
            const pointField& newFaceCentres,
            const pointField& newCellCentres,
            const vectorField& newFaceAreas,
            pointField& tPts,
            pointField& newPoints
        );

        // Reset medial points and rejected motion points
        void resetPoints
        (
            const label iter,
            const meshControl& controller,
            const labelList& layerPointType,
            const labelList& layerFaceType,
            const boolList& layerEdges,
            const scalarField& medialDist,
            const scalarField& medialLayerHeight,
            const pointField& pointOrigin,
            const pointField& pointSurfNormal,
            const pointField& tPts,
            const pointField& outerPos,
            scalarField& newCellVolumes,
            vectorField& newFaceAreas,
            pointField& newFaceCentres,
            pointField& newCellCentres,
            pointField& newPoints,
            labelList& resetPts,
            bool resetErrors
        );

        // Project outer layer shell to cell-centre to cell-centre mid point
        void projectOuterLayerFaces
        (
            const autoPtr<indirectPrimitivePatch>& outerShellPtr,
            const labelList& layerPointType,
            const labelList& layerFaceType,
            const labelList& ftrPointOrigin,
            const labelList& stationaryPts,
            const boolList& stackEdges,
            const PackedBoolList& isMasterFace,
            const vectorField& newFaceAreas,
            const pointField& newCellCentres,
            const pointField& pointOrigin,
            const vectorField& pointSurfNormal,
            const pointField& tPts,
            pointField& newPoints,
            labelList& resetPts
        );

        // Limit layer projection height
        void limitLayerHeight
        (
            const labelList& ftrPointOrigin,
            const labelList& layerPointType,
            const scalarField& maxLayerThickness,
            const vectorField& pointSurfNormal,
            const pointField& pointOrigin,
            pointField& pts
        );

        // Select layer growth methods and maximum layer thickness
        void setLayerMethod
        (
            const indirectPrimitivePatch& pp,
            const labelList& meshPoints,
            const labelList& ppPtLevel,
            const labelList& ppPtLevelMin,
            List<Tuple2<word,scalar>>& layerMethod,
            scalarField& ppMaxLayerThickness
        );

        // Adjust number of layers based on target expansion ratio
        // when using first cell height set as the free parameter
        void adjustNumLayers
        (
            const indirectPrimitivePatch& pp,
            const List<Tuple2<word,scalar>>& layerMethod,
            const scalarField& initExpansion,
            const pointField& pointOrigin,
            const pointField& outerPos,
            const labelList& totalNumLayers,
            labelList& adjustedNLayers
        );

        // Smooth layer cells and volume mesh
        void smoothLayerStack
        (
            const meshControl& controller,
            const label nSmoothIter,
            const dictionary& grownUpGeometryDict,
            const dictionary& grownUpZoneGeometryDict,
            const labelList& layerOffset,
            labelList& adjustedNLayers,
            const bool adjustFinalLayers = false
        );

        // Stretch layer points
        void stretchLayers
        (
            const meshControl& controller
        );

        void fitLayerPointStack();

        labelListList createLayerPointTagList();

        layerPointStack buildLayerStack
        (
            const label& pI,
            const labelListList& pointTagList,
            const label& pL,
            boolList& addedToStack,
            labelList& totalNuOfLayers
        );

        bool compareTags(const labelList& pt1, const labelList& pt2);

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
