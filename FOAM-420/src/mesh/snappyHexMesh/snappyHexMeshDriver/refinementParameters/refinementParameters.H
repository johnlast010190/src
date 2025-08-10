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
    (c) 2015-2023 OpenCFD Ltd.
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::refinementParameters

Description
    Simple container to keep together refinement specific information.

SourceFiles
    refinementParameters.C

\*---------------------------------------------------------------------------*/

#ifndef refinementParameters_H
#define refinementParameters_H

#include "db/dictionary/dictionary.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "primitives/bools/Switch/Switch.H"
#include "snappyHexMeshDriver/refinementParameters/wordPairHashTable.H"
#include "refinementSurfaces/surfaceZonesInfo.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "cfdTools/general/include/fvCFD.H"
#include "polyTopoChange/hexRef8/hexRef8.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class refinementParameters Declaration
\*---------------------------------------------------------------------------*/

class refinementParameters
{
    // Private data

        //- Total number of cells
        const label maxGlobalCells_;

        //- Per processor max number of cells
        const label maxLocalCells_;

        //- When to stop refining
        const label minRefineCells_;

        //- Curvature
        scalar curvature_;

        //- Number of layers between different refinement levels
        const label nBufferLayers_;

        //- Planarity criterion
        scalar planarAngle_;


        // Selection of areas

            //- Areas not to keep
            const pointField locationsOutsideMesh_;

            //- Areas to keep
            mutable pointField locationsInMesh_;

            //- Region for location
            mutable wordList zonesInMesh_;

            //- Locations in zone to test for leakage
            pointField zonesToTestLocations_;

            //- Names of zones to test for leakage
            wordList zonesToTestNames_;

            //- Regions below this size removed from zone
            const label minZoneRegionSize_;

            //- Information on how to handle faces on faceZones
            dictionary faceZoneControls_;

        //- Use old topology based problem-cell removal (cells with 8 points
        //  on surface)
        Switch useTopologicalSnapDetection_;

        //- Enable re-zoning of problem cells
        Switch rezoneProblemCells_;

        //- Enable new zoning method
        Switch newZoning_;

        //- Enable merging of single couples
        Switch mergeSingleCouples_;

        //- Perform surface wrapping
        Switch wrap_;

        //- Perform balancing before refinement
        Switch balanceRefine_;

        //- Allowed load unbalance
        scalar maxLoadUnbalance_;

        //- Trigger cell count to start balancing
        label maxCellUnbalance_;

        //- Use estimated snap cell weighting for re-balancing mesh
        // before snapping (> 0 to activate)
        scalar snapWeights_;

        //- Define number of rays used in proximity check
        label nProxRays_;

        //- Define angle used for proximity checks
        scalar proximityAngle_;

        //- Whether to split cells at the end of refinement
        Switch splitCells_;

        //- Whether to perform crack detection on zoned surfaces
        Switch zoneCracks_;

        //- Whether to refine boundary problem cells
        Switch refineBoundaryProblemCells_;

        Switch handleSnapProblems_;

        Switch interfaceRefine_;

        //- Whether to remove boundary cells that cannot be split
        // prior to dualisation
        Switch removeUnsplittable_;

        //- Whether to perform additional dual refinement to improve
        // surface topology
        Switch additionalDualRefine_;

        //- If unnamed intersection between cell zones, keep the baffle face
        Switch interZoneBaffles_;

        //- If set perform leak check on all refined mesh removal stages
        // if locationsOutsideMesh is set
        Switch fullLeakChecks_;

        //- If set perform zonal leak checks
        // zones must be defined by names locationsInMesh
        Switch zoneLeakChecks_;

        //- Number of smoothing iterations of curvature field
        label curvatureSmooth_;

        //- Holes below this size are closed (extrude method)
        label holeSize_;

        //- Maximum number of iterations for resolving corner cells
        label maxCornerIter_;

        //- If resolving corner cells remove any associated baffle
        // faces connected to these regions
        Switch cornerRemoveBaffles_;

        //Allow removal of cells with all points on boundary
        // if only a single face is a grown up patch face (extrude/dual)
        Switch removeBoundaryAnyGrownUp_;

        //Remove from faceZone if faces are boundary
        Switch removeBoundaryFaceZones_;

        //Filter out feature lines below this length
        scalar minFeatureLength_;

        //Try to add outside cells at contact areas (extrude method)
        Switch addContactCells_;

        //Use added cell removal from standard method for extrude method
        Switch extrudeExtraRemoval_;

        //- Number added locations not set by locationsInMesh keyword
        mutable label nAddedLocations_;

        //Refine outside gap cells to neighbour level
        Switch refineOutsideGapCell_;

        //If named locations rezone trapped cellls between regions
        // with added minZoneRegionSize filter
        Switch namedLocationsRezone_;

        //Whether to merge boundary faces before extruding them
        Switch mergePreExtrude_;

        //Whether to try to remove holes for free standing faceType
        // baffle/internal using the surface orientation for sideness
        Switch filterFreeStandingHoles_;

        //- Wanted decomposition method for VDB grids
        word vdbDecompMethod_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        refinementParameters(const refinementParameters&);

        //- Disallow default bitwise assignment
        void operator=(const refinementParameters&);


public:

    // Constructors

        //- Construct from dictionary - new syntax
        refinementParameters(const dictionary& dict);


    // Member Functions

        // Access

            //- Total number of cells
            label maxGlobalCells() const
            {
                return maxGlobalCells_;
            }

            //- Per processor max number of cells
            label maxLocalCells() const
            {
                return maxLocalCells_;
            }

            //- When to stop refining
            label minRefineCells() const
            {
                return minRefineCells_;
            }

            //- Curvature
            scalar curvature() const
            {
                return curvature_;
            }

            //- Angle when two intersections are considered to be planar
            scalar planarAngle() const
            {
                return planarAngle_;
            }

            //- Number of layers between different refinement levels
            label nBufferLayers() const
            {
                return nBufferLayers_;
            }

            //- Number of gap refinement iterations
            label minZoneRegionSize() const
            {
                return minZoneRegionSize_;
            }

            //- Areas to keep
            const pointField& locationsInMesh() const
            {
                return locationsInMesh_;
            }

            //- Perform rezoning of problematic cells
            bool rezoneProblemCells() const
            {
                return rezoneProblemCells_;
            }

            //- Optional points which are checked to be outside the mesh
            const pointField& locationsOutsideMesh() const
            {
                return locationsOutsideMesh_;
            }

            //- Per area the zone name
            const wordList& zonesInMesh() const
            {
                return zonesInMesh_;
            }

            //- Areas to test
            const pointField& zonesToTestLocations() const
            {
                return zonesToTestLocations_;
            }

            //- per test area the zone name
            const wordList& zonesToTestNames() const
            {
                return zonesToTestNames_;
            }

            //- Use old topology based problem-cell removal
            //  (cells with 8 points on surface)
            bool useTopologicalSnapDetection() const
            {
                return useTopologicalSnapDetection_;
            }

            //- Set new method for creating zones
            bool newZoning() const
            {
                return newZoning_;
            }

            //- Set merging of single couples
            bool mergeSingleCouples() const
            {
                return mergeSingleCouples_;
            }

            //- Allow surface wrapping
            bool wrap() const
            {
                return wrap_;
            }

            //- Balancing before refinement
            bool balanceRefine() const
            {
                return balanceRefine_;
            }

            //- Allowed load unbalance
            scalar maxLoadUnbalance() const
            {
                return maxLoadUnbalance_;
            }

            //- Trigger cell count to start balancing
            label maxCellUnbalance() const
            {
                return maxCellUnbalance_;
            }

            //- Use weights to balance before snapping
            scalar snapWeights() const
            {
                return snapWeights_;
            }

            bool handleSnapProblems() const
            {
                return handleSnapProblems_;
            }

            //- Refine cell with opposite faces with different refinement level
            bool interfaceRefine() const
            {
                return interfaceRefine_;
            }

            //- Whether to remove unsplittable cells prior to dualisation
            bool removeUnsplittable() const
            {
                return removeUnsplittable_;
            }

            //- Whether to perform additional refinement on dual mesh
            bool additionalDualRefine() const
            {
                return additionalDualRefine_;
            }

            //- Whether to baffle un-named intersected faces between zones
            bool interZoneBaffles() const
            {
                return interZoneBaffles_;
            }

            //- Whether to perform final leak check
            bool fullLeakChecks() const
            {
                return fullLeakChecks_;
            }

            //- Whether to perform zone leak check
            bool zoneLeakChecks() const
            {
                return zoneLeakChecks_;
            }

            //- Number of smoothing iterations of curvature field
            label curvatureSmooth() const
            {
                return curvatureSmooth_;
            }

            //- Holes below this size are closed (extrude method)
            label holeSize() const
            {
                return holeSize_;
            }

            //- Maximum number iterations for resolving corner cells
            label maxCornerIter() const
            {
                return maxCornerIter_;
            }

            //- Whether to remove baffles in corner regions
            bool cornerRemoveBaffles() const
            {
                return cornerRemoveBaffles_;
            }

            //- Whether to remove cells with all boundary points
            // if contains only a single grown up patch face
            bool removeBoundaryAnyGrownUp() const
            {
                return removeBoundaryAnyGrownUp_;
            }

            //- Whether to remove faces from zone if type boundary
            bool removeBoundaryFaceZones() const
            {
                return removeBoundaryFaceZones_;
            }

            //- Filter feature lines below this length
            scalar minFeatureLength() const
            {
                return minFeatureLength_;
            }

            //- Whether to try to keep outside contact cells
            bool addContactCells() const
            {
                return addContactCells_;
            }

            //- Whether to remove additional cells for extrude meshing
            bool extrudeExtraRemoval() const
            {
                return extrudeExtraRemoval_;
            }

            //- Whether to refine outsider gap cells
            bool refineOutsideGapCell() const
            {
                return refineOutsideGapCell_;
            }

            //- Whether to rezone cells trapped between named locations
            bool namedLocationsRezone() const
            {
                return namedLocationsRezone_;
            }

            //- Whether to merge face pre extrusion
            bool mergePreExtrude() const
            {
                return mergePreExtrude_;
            }

            //- Whether to filter out free standing holes
            bool filterFreeStandingHoles() const
            {
                return  filterFreeStandingHoles_;
            }

            //- Number of proximity rays used
            label nProxRays() const
            {
                return nProxRays_;
            }

            //- Angle used for proximity checks
            scalar proximityAngle() const
            {
                return proximityAngle_;
            }

            //- Splitting cells at end of refinement
            bool splitCells() const
            {
                return splitCells_;
            }

            //- whether to perform crack detection on zoned surfaces
            bool zoneCracks() const
            {
                return zoneCracks_;
            }

            //- whether to refine boundary problem cells
            bool refineBoundaryProblemCells() const
            {
                return refineBoundaryProblemCells_;
            }

            //- decomposition method for VDB grids
            word vdbDecompMethod() const
            {
                return vdbDecompMethod_;
            }

        // Other

            //- Check if faceZoneControls has baffle checks
            bool isFaceZoneBaffleChecks
            (
                const word& fzName
            ) const;

            //- Check if faceZoneControls is boundary type
            bool isFaceZoneControlsBoundary
            (
                const word& fzName
            ) const;

            //- Get patchInfo and faceType for faceZone
            dictionary getZoneInfo
            (
                const word& fzName,
                surfaceZonesInfo::faceZoneType& faceType
            ) const;

            //- Add cellZones to mesh. Return indices of cellZones (or -1)
            labelList addCellZonesToMesh(polyMesh&) const;

            //- Checks that cells are in mesh. Returns cells (or -1) they
            //  are in.
            static labelList findCells
            (
                const bool checkInsideMesh,
                const polyMesh&,
                const pointField& locations
            );

            //- Extract indices of named locations (so excludes 'keepPoints')
            static labelList zonedLocations(const wordList& zonesInMesh);

            //- Extract indices of unnamed locations ('keepPoints')
            static labelList unzonedLocations(const wordList& zonesInMesh);

            //- Add additional locationsInMesh outside of keyword
            void addToLocationsInMesh
            (
                const searchableSurfaces& allGeometry,
                const dictionary& geometryDict
            );

            //- Mesh regionise and filter out duplicate locationsInMesh
            void filterLocations
            (
                const fvMesh& mesh,
                const hexRef8& meshCutter,
                const labelList& globalRegion1,
                const labelList& namedSurfaceIndex
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
