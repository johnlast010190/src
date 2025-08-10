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
    (c) 2020 Esi Ltd.

Class
    Foam::autoSplitCells

Description
    Split boundary hex and split-hex cells into prisms based on
    interior angles and face warpage

SourceFiles
    autoSplitCells.C

\*---------------------------------------------------------------------------*/

#ifndef autoSplitCells_H
#define autoSplitCells_H

#include "meshRefinement/meshRefinement.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class autoSplitCells Declaration
\*---------------------------------------------------------------------------*/

class autoSplitCells
{
    // Private data

        //- Mesh+level information
        meshRefinement& meshRefiner_;

        const labelHashSet excludePatchSet_;

        //- Is merging on an extruded mesh type
        bool extruded_;

        //- Does extrude method perform a pre-merge step
        bool preMerged_;

        //- Revert if mesh contains errors
        bool revertErrors_;

    // Private Member Functions

        //- Is a valid patch to split
        bool isValidPatch(const label patchi);

        //- Check for edges collapses
        label checkEdges
        (
            const pointField& origPoints,
            labelHashSet& setPtr
        );

        //- Correct mesh points if errors introduced
        void correctPoints
        (
            vectorField& disp
        );

        //- Move feature points to improve split (feature face cutting)
        void moveFeaturePoints
        (
            const List<Tuple2<label,labelPair>>& cutFaces,
            const scalar featureAngle,
            const bool moveFeaturePts,
            const bool moveConvexPts
        );

        //- Move hanging points on split faces (warped face cutting)
        void moveSplitFaceHangingPts();

        //- Flatten refinement interface faces (wraped face cutting)
        void flattenRefinementFaces();

        //- Filter out faces which cannot produce valid cut cells
        void filterCutFaces
        (
            List<Tuple2<label,labelPair>>& cutFaces,
            const bool preventZoneCuts = false
        );

        //- Calculate and split hex cells into prisms
        List<labelPair> split
        (
            const List<Tuple2<label,labelPair>>& cutFaces,
            const bool checkCutAreas
        );

        //- Re-merge any split cells that introduce errors
        void remerge(const List<labelPair>& splitCells);

        //- Split and re-merge cells based on supplied cut faces
        void update
        (
            List<Tuple2<label,labelPair>>& cutFaces,
            bool checkCutAreas
        );

        //- Identify faces with warpage greater than supplied setting
        List<Tuple2<label,labelPair>> warpedFacesToSplit
        (
            const scalar warpage,
            const scalar fchwarpage,
            const List<wordList>& layerSpec,
            const labelList& patchNumLayers,
            const scalarField& patchFCH,
            const scalarField& patchExpansionRatio,
            const bool curvatureSplit,
            const bool mergedFaces,
            const bool fourAnchors
        );

        //- Identify faces on feature edges for splitting
        List<Tuple2<label,labelPair>> featureFacesToSplit
        (
            const scalar featureAngle
        );

public:

    //- Runtime type information
    ClassName("autoSplitCells");

    //- Cell splitting at warped faces
    void splitWarpedCells
    (
        const scalar warpage,
        const scalar fchwarpage,
        const List<wordList>& layerSpec,
        const labelList& patchNumLayers,
        const scalarField& patchFCH,
        const scalarField& patchExpansionRatio,
        const bool curvatureSplit,
        const bool convexSplit,
        const bool mergedFaces
    );

    //- Cell splitting at feature edge cells
    void splitFeatureCells
    (
        const scalar featureAngle,
        const bool moveFeaturePts,
        const bool moveConvexPts
    );

    //- Cell splitting using pre-calculated cuts
    void splitPreCalculated
    (
        List<Tuple2<label,labelPair>>& cutFaces
    );

    // Constructors
    autoSplitCells
    (
        meshRefinement& meshRefiner,
        labelList excludePatches,
        bool extruded,
        bool preMerged,
        bool revertErrors = true
    );

    //- Destructor
    ~autoSplitCells();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
