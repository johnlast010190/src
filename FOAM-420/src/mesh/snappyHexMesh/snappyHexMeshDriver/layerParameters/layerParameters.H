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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::layerParameters

Description
    Simple container to keep together layer specific information.

SourceFiles
    layerParameters.C

\*---------------------------------------------------------------------------*/

#ifndef layerParameters_H
#define layerParameters_H

#include "db/dictionary/dictionary.H"
#include "fields/Fields/scalarField/scalarField.H"
#include "primitives/ints/lists/labelList.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "primitives/bools/Switch/Switch.H"
#include "patchDist/wallPoint/wallPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyBoundaryMesh;
class refinementSurfaces;

/*---------------------------------------------------------------------------*\
                           Class layerParameters Declaration
\*---------------------------------------------------------------------------*/

class layerParameters
{
public:

        //- Whether to overwrite grown up patches
        enum grownUpType
        {
            GROWN,
            COLLAPSED,
            AUTO
        };

        static const NamedEnum<grownUpType, 3> grownUpTypeNames;

private:

    // Private classes

        // To combineReduce a labelList. Filters out duplicates.
        class uniqueEqOp
        {

        public:

            void operator()(labelList& x, const labelList& y) const
            {
                if (x.size() == 0)
                {
                    if (y.size() > 0)
                    {
                        x = y;
                    }
                }
                else
                {
                    forAll(y, yi)
                    {
                        if (findIndex(x, y[yi]) == -1)
                        {
                            label sz = x.size();
                            x.setSize(sz+1);
                            x[sz] = y[yi];
                        }
                    }
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
                        x /= (mag(x) + SMALL);
                    }
                }
            }
        };

    // Static data members

        //- Default angle for faces to be convcave
        static const scalar defaultConcaveAngle;


    // Private data

        const dictionary dict_;

        // Has the mesh been dualised
        const bool dualExtrude_;

        // Per patch (not region!) information

            //- How many layers to add.
            labelList numLayers_;

            //- Adjust the number of layers based on a target expansion
            // ratio range (extrude method)
            List<Tuple2<scalar, scalar>> targetExpansions_;

            //- manually set which patches to grow up .
            List<grownUpType> grownUp_;

            //- automatically detect which  patches to grow layers up
            labelList grownUpIDs_;

            //- Baffle patches to grow only on single side. Side determined
            // by surface normal direction
            labelList singleSidedIDs_;

            //- List of region entry's where layers cannot be grown
            PtrList<entry> excludedRegions_;

            scalar defaultExpansionRatio_;

            scalarField expansionRatio_;

            scalarField finalLayerThickness_;

            scalarField minThickness_;

            scalarField maxLayerThickness_;

            scalarField fch_;

            List<Switch> relativeSizes_;

            List<Switch> reSnap_;

            scalar mergePatchFacesAngle_;

            labelList minLayers_;

            scalar cosAngleTermination_;
            scalar concaveAngle_;

            label nGrow_;

            label maxLayerIter_;

            bool growConcaveEdge_;

            bool growConvexEdge_;

            bool truncateFromWall_;

            bool growUpPatches_;

            label terminationStrategy_;

            label nSmoothSurfaceNormals_;
            label nSmoothNormals_;
            label nSmoothThickness_;
            scalar maxFaceThicknessRatio_;

            scalar minMedianAxisAngleCos_;
            label nBufferCellsNoExtrude_;

            Switch additionalReporting_;

            word meshShrinker_;

            label nRelaxedIter_;

            label layerRecovery_;

            bool growZoneLayers_;

            bool preBalance_;

            bool noErrors_;

            bool splitBaffleTerminations_;

            bool writeVTK_;

            scalar projectGrownUp_;

            scalar maxCellDistortion_;

            scalar grownUpAngleTerminateCos_;

            scalar maxProjectionDist_;

            label nSmoothPerDualLayer_;

            scalar fanAngle_;

            scalar dualConcaveCollapse_;

            scalar layerInterfaceWeights_;

            scalar zoneLayersScaling_;

            Switch extrudeBlend_;

            scalar dualOrtho_;

            Switch dualReSnapZones_;

            List<Switch> fixedFCH_;

            scalar globalWarpedSplit_;

            scalar fchWarpedSplit_;

            Switch curvatureSplit_;

            Switch mergeSplitCells_;

            Switch fanTetSplit_;

            label maxMergePreIter_;
            label maxMergePostIter_;

            //When layers jump at patch interfaces
            // whether to increment from the lower layer count patch
            Switch incrementLower_;

            // whether to perform faster mesh gemotry update (extrude/dual)
            Switch fastGeomUpdate_;

            // whether to used fixed normal slip rather than slip
            Switch fixedNormalSlip_;

            // Estrude/Dual stop projecting when squish too high (-1 disable)
            scalar squishTol_;

            // whether to add additional extruded layers in two steps
            Switch twoStageExtrusion_;

            //  Limit minimum stretching value for extrude method
            scalar minStretch_;

            //  Number of layers with fch imposed (extrude method)
            label nFCHLayers_;

            List<Switch> cornerCollapse_;

            List<Switch> baffleCollapse_;

            List<wordList> layerSpec_;

    // Private Member Functions

        //- Set patch-wise layer growth method
        void setLayerGrowthMethod
        (
            const dictionary& dict,
            const polyBoundaryMesh& boundaryMesh
        );

        //- Set patch-wise patches to possibly collapse layers on
        void setCollapsePatches
        (
            const dictionary& dict,
            const polyBoundaryMesh& boundaryMesh
        );

        //- Disallow default bitwise copy construct
        layerParameters(const layerParameters&);

        //- Disallow default bitwise assignment
        void operator=(const layerParameters&);


public:

    // Constructors

        //- Construct from dictionary
        layerParameters
        (
            const dictionary& dict,
            const polyBoundaryMesh&,
            const bool dualExtrude = false
        );


    // Member Functions

        const dictionary& dict() const
        {
            return dict_;
        }

        // Per patch information

                //- Return layer growth type
                const List<wordList>& layerSpec() const
                {
                    return layerSpec_;
                }

                //- How many layers to add.
                const labelList& numLayers() const
                {
                    return numLayers_;
                }

                //- Target expansion range
                const List<Tuple2<scalar,scalar>>& targetExpansions() const
                {
                    return targetExpansions_;
                }

                //- Whether to grow up neighbouring patches
                const bool& growUpPatches()
                {
                    return growUpPatches_;
                }

                //- Read the patches to try growing layers up
                void setGrownUpIDs
                (
                    const polyBoundaryMesh& boundaryMesh
                );

                //- Which patches to grow up.
                const labelList& grownUpIDs() const
                {
                    return grownUpIDs_;
                }

                //- Which patches to grow on a single side.
                const labelList& singleSidedIDs() const
                {
                    return singleSidedIDs_;
                }

                const PtrList<entry>& excludedRegions() const
                {
                    return excludedRegions_;
                }

                //- First cell height.
                const scalarField& fch() const
                {
                    return fch_;
                }

                //- Whether to re-snap to grown up patches (dual mesh generator).
                const List<Switch>& reSnap() const
                {
                    return reSnap_;
                }

                //- Minimum number of layers to grow on a patch.
                const labelList& minLayers() const
                {
                    return minLayers_;
                }

                //- Maximum obtainable layer height.
                const scalarField& maxLayerThickness() const
                {
                    return maxLayerThickness_;
                }

                // Expansion factor for layer mesh
                const scalarField& expansionRatio() const
                {
                    return expansionRatio_;
                }

                //- Are size parameters relative to inner cell size or
                //  absolute distances.
                const List<Switch>& relativeSizes() const
                {
                    return relativeSizes_;
                }

                //- Wanted thickness of final added cell layer. If multiple
                //  layers is the thickness of the layer furthest away
                //  from the wall (i.e. nearest the original mesh)
                //  If relativeSize() this number is relative to undistorted
                //  size of the cell outside layer.
                const scalarField& finalLayerThickness() const
                {
                    return finalLayerThickness_;
                }

                //- Minimum thickness of cell layer. If for any reason layer
                //  cannot be above minThickness do not add layer.
                //  If relativeSize() this number is relative to undistorted
                //  size of the cell outside layer.
                const scalarField& minThickness() const
                {
                    return minThickness_;
                }

            scalar concaveAngle() const
            {
                return concaveAngle_;
            }

            //- If points get not extruded do nGrow layers of connected faces
            //  that are not grown. Is used to not do layers at all close to
            //  features.
            label nGrow() const
            {
                return nGrow_;
            }

            // Maximum number of layer generation iterations
            label maxLayerIter() const
            {
                return maxLayerIter_;
            }

            // Strategy for layer termination cell
            // (Default: 0, Merge cell with one above: 1)
            label terminationStrategy() const
            {
                return terminationStrategy_;
            }

            // Try growing layers around concave edges
            bool growConcaveEdge() const
            {
                return growConcaveEdge_;
            }

            // Try growing layers around convex edges
            bool growConvexEdge() const
            {
                return growConvexEdge_;
            }

            // Truncation method used at layer terminations
            bool truncateFromWall() const
            {
                return truncateFromWall_;
            }

            // Number of smoothing iterations of surface normals
            label nSmoothSurfaceNormals() const
            {
                return nSmoothSurfaceNormals_;
            }

            //- Number of smoothing iterations of interior mesh movement
            //  direction
            label nSmoothNormals() const
            {
                return nSmoothNormals_;
            }

            scalar mergePatchFacesAngle() const
            {
                return mergePatchFacesAngle_;
            }

            // Smooth layer thickness over surface patches
            label nSmoothThickness() const
            {
                return nSmoothThickness_;
            }

            // Stop layer growth on highly warped cells
            scalar maxFaceThicknessRatio() const
            {
                return maxFaceThicknessRatio_;
            }

            // Control to prevent layers growing around features
            scalar cosAngleTermination() const
            {
                return cosAngleTermination_;
            }

            //- Angle used to pick up medial axis points
            scalar grownUpAngleTerminateCos() const
            {
                return grownUpAngleTerminateCos_;
            }

            //- Maximum distance to project displacement
            scalar maxProjectionDist() const
            {
                return maxProjectionDist_;
            }

            //- Create buffer region for new layer terminations
            label nBufferCellsNoExtrude() const
            {
                return nBufferCellsNoExtrude_;
            }

            const Switch& additionalReporting() const
            {
                return additionalReporting_;
            }

            //- Type of mesh shrinker
            const word& meshShrinker() const
            {
                return meshShrinker_;
            }

            //- Number of iterations after which relaxed motion rules
            //  are to be used.
            label nRelaxedIter() const
            {
                return nRelaxedIter_;
            }

            //- Increment for regrowing collapsed layers
            label layerRecovery() const
            {
                return layerRecovery_;
            }

            //- Whether to grow layers at a zone interfaces. If true will
            // change zone interface conformity
            bool growZoneLayers() const
            {
                return growZoneLayers_;
            }

            //- Whether to rebalance mesh based on estimated cell size
            bool preBalance() const
            {
                return preBalance_;
            }

            //- Whether to exit layer method with errors
            bool noErrors() const
            {
                return noErrors_;
            }

            //- Whether split baffles at layer terminations
            bool splitBaffleTerminations() const
            {
                return splitBaffleTerminations_;
            }

            //- Whether VTK file showing number of layers grown
            bool writeVTK() const
            {
                return writeVTK_;
            }

            //- Whether to project points onto grown up patches
            // (assumes grown up patches are planar)
            scalar projectGrownUp() const
            {
                return projectGrownUp_;
            }

            //- Do not allow displacement of cell to be more than
            // this ratio of cell size
            scalar maxCellDistortion() const
            {
                return maxCellDistortion_;
            }

            //- Number of smoothing iterations for every added dual mesh layer
            label nSmoothPerDualLayer() const
            {
                return nSmoothPerDualLayer_;
            }

            //- Concave angle to fan layers around
            scalar fanAngle() const
            {
                return fanAngle_;
            }

            //- Concave angle to collapse layers around
            scalar dualConcaveCollapse() const
            {
                return dualConcaveCollapse_;
            }

            //- Weights for smoothing of layer and non-layer cells
            scalar layerInterfaceWeights() const
            {
                return layerInterfaceWeights_;
            }

            //- extend to which zones faces are allowed to move during smoothing
            scalar zoneLayersScaling() const
            {
                return zoneLayersScaling_;
            }

            //- fit layer stacks based on blended normal (extrusion method only)
            Switch extrudeBlend() const
            {
                return extrudeBlend_;
            }

            //- Disable smoothing when orthogonality greater than this angle
            // in degrees (180 degrees to disable)
            scalar dualOrtho() const
            {
                return dualOrtho_;
            }

            //- Re-snap to internal zones during layer addition
            Switch dualReSnapZones() const
            {
                return dualReSnapZones_;
            }

            //- Patches to collapse layers at baffle edges
            List<Switch> baffleCollapse() const
            {
                return baffleCollapse_;
            }

            //- Patches to collapse layers at complex corners
            List<Switch> cornerCollapse() const
            {
                return cornerCollapse_;
            }

            //- Set fixed fch for non-dual mesh generator
            List<Switch> fixedFCH() const
            {
                return fixedFCH_;
            }

            //- Split cells based on face warpage
            scalar globalWarpedSplit() const
            {
                return globalWarpedSplit_;
            }

            //- Split cells based on specified fch
            scalar fchWarpedSplit() const
            {
                return fchWarpedSplit_;
            }

            //- When splitting based on warpage whether split
            // based on curvature (true) or to prevent convex cells (false)
            Switch curvatureSplit() const
            {
                return curvatureSplit_;
            }

            //- Whether to remerge outer layer split cells
            Switch mergeSplitCells() const
            {
                return mergeSplitCells_;
            }

            //- Whether to split fan corner cells to form tets
            Switch fanTetSplit() const
            {
                return fanTetSplit_;
            }

            //- Whether to increment layers from fewer layers patch
            Switch incrementLower() const
            {
                return incrementLower_;
            }

            //- Whether to perform faster geometry update
            Switch fastGeomUpdate() const
            {
                return fastGeomUpdate_;
            }

            //- Whether to use fixed normal slip bc
            Switch fixedNormalSlip() const
            {
                return fixedNormalSlip_;
            }

            //- Whether to perform two stage added layer extrusion
            Switch twoStageExtrusion() const
            {
                return twoStageExtrusion_;
            }

            //- Limit minimum layer stretch
            scalar minStretch() const
            {
                return minStretch_;
            }

            //- Number of layers where fch imposed
            label nFCHLayers() const
            {
                return nFCHLayers_;
            }

            //- Whether to disable projection based on squish metric
            scalar squishTol() const
            {
                return squishTol_;
            }

            //- Whether to merge cells after layer addition
            label maxMergePreIter() const
            {
                return maxMergePreIter_;
            }

            //- Whether to merge cells post layer optimisation
            label maxMergePostIter() const
            {
                return maxMergePostIter_;
            }

            //- Whether to try collapsing layers at baffle edges
            // for dual or extrude method
            bool tryBaffleCollapse() const;

            //- Whether to try collapsing layers at corner points
            // for dual or extrude method
            bool tryCornerCollapse() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
