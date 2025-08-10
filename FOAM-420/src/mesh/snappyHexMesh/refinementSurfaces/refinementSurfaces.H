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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 ICON CFD
    (c) 2011-2014 OpenFOAM Foundation

Class
    Foam::refinementSurfaces

Description
    Container for data on surfaces used for surface-driven refinement.
    Contains all the data about the level of refinement needed per
    surface.

SourceFiles
    refinementSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef refinementSurfaces_H
#define refinementSurfaces_H

#include "triSurface/fields/triSurfaceGeoMesh.H"
#include "triSurface/fields/triSurfaceFields.H"
#include "primitives/Vector/lists/vectorList.H"
#include "meshes/primitiveShapes/objectHit/pointIndexHit.H"
#include "refinementSurfaces/surfaceZonesInfo.H"
#include "algorithms/indexedOctree/volumeType.H"
#include "meshes/primitiveShapes/point/pointList.H"
#include "primitives/Tuple2/Tuple2.H"
#include "meshes/boundBox/boundBox.H"
#include "regionSplit/regionSplit.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "primitives/bools/Switch/Switch.H"
#include "fvMesh/fvMesh.H"
#include "coordinate/systems/coordinateSystem.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;
class shellSurfaces;
class triSurfaceMesh;

/*---------------------------------------------------------------------------*\
                           Class refinementSurfaces Declaration
\*---------------------------------------------------------------------------*/

class refinementSurfaces
{
    // Private data

        //- Reference to all geometry.
        const searchableSurfaces& allGeometry_;

        //- Indices of surfaces that are refinement ones
        labelList surfaces_;

        //- Surface name (word)
        wordList names_;

        //- List of surface zone (face and cell zone) information
        PtrList<surfaceZonesInfo> surfZones_;

        //- From local region number to global region number
        labelList regionOffset_;

        //- Per 'interface' surface : Whether to move cell centres falling
        // on surface
        scalarField moveCentroids_;

        //- Per 'interface' surface : Whether to close gaps
        // on surface above the wrap level
        labelList wrapLevel_;

        //- Per 'interface' surface : Whether to use
        // for auto bounding box calculation
        boolList boundingSurface_;

        //- Per 'interface' surface : Whether to refine boundary up to
        // maximum level
        boolList refineBoundary_;

        //- Per 'interface' surface : Whether setup as single region
        boolList singleRegion_;

        //Remove this number of disconnected regions during wrap
        label nRegionsWrapRemoval_;

        //- Whether to orient unnamed surfaces
        Switch reorient_;

        //- From global region number to refinement level
        labelList minLevel_;

        //- From global region number to refinement level
        labelList maxLevel_;

        //- From global region number to small-gap level
        labelList gapLevel_;

        //- From global region number to small-gap level specification
        List<FixedList<label, 3>> extendedGapLevel_;

        //- From global region number to side of surface to detect
        List<volumeType> extendedGapMode_;

        //- From global region number to perpendicular angle
        scalarField perpendicularAngle_;

        //- From global region number to feature refinement angle
        scalarField featureRefineAngle_;

        //- From global region number to feature refinement only switch
        boolList featureRefineOnly_;

        //- From global region number to closure active
        boolList singleCellClosure_;

        //- From global region number to thin gap active
        boolList thinGap_;

        //- From global region number to where corner cell keeping active
        boolList cornerCellRegions_;

        //- From global region number to min baffle angle
        scalarField minBaffleAngle_;

        //- From global region number to the allowable increment on the
        //max surface refinement level applicable to proximity refinement
        labelList proxLevelIncr_;

        //- From global region number to proximity direction method
        // no checks (0), outward (-1), inward (1)
        labelList proxDir_;

        //- From global region number to curvature angle in Radians
        scalarField curvature_;

        //- From global region number to the maximum number of cells that can
        //span proximity gap
        labelList proxMaxCells_;

        //- From global region number to patchType
        PtrList<dictionary> patchInfo_;


    // Private Member Functions

        //- Given intersection results with geom detect local shell refinement
        //  level (possibly cached on triangles of geom)
        labelList findHigherLevel
        (
            const searchableSurface& geom,
            const shellSurfaces& shells,
            const List<pointIndexHit>& intersectionInfo,
            const labelList& surfaceLevel,
            const bool threaded = false
        ) const;

        //- Disallow default bitwise copy construct
        refinementSurfaces(const refinementSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const refinementSurfaces&);


public:

    //- Runtime type information
    ClassName("refinementSurfaces");

    // Constructors

        //- Construct from surfaces and dictionary
        refinementSurfaces
        (
            const searchableSurfaces& allGeometry,
            const dictionary& surfaceDict,
            const dictionary& refineDict,
            const label gapLevelIncrement
        );

        //- Construct from components
        refinementSurfaces
        (
            const searchableSurfaces& allGeometry,
            const labelList& surfaces,
            const wordList& names,
            const PtrList<surfaceZonesInfo>& surfZones,
            const labelList& regionOffset,
            const labelList& minLevel,
            const labelList& maxLevel,
            const labelList& gapLevel,
            const labelList& proxLevelIncr,
            const labelList& proxDir,
            const scalarField& curvature,
            const scalarField& perpendicularAngle,
            PtrList<dictionary>& patchInfo
        );


    // Member Functions

        // Access

            const searchableSurfaces& geometry() const
            {
                return allGeometry_;
            }

            const labelList& surfaces() const
            {
                return surfaces_;
            }

            //- Names of surfaces
            const wordList& names() const
            {
                return names_;
            }

            //- Per 'interface' surface : whether to move cell-centres sitting
            //  on surface
            const scalarField& moveCentroids() const
            {
                return moveCentroids_;
            }


            //- Per 'interface' surface : whether to close gaps
            //  above specified wrap level
            const labelList& wrapLevel() const
            {
                return wrapLevel_;
            }

            //- Per 'interface' surface : whether to refine to surface boundary
            const boolList& refineBoundary() const
            {
                return refineBoundary_;
            }

            //- Per 'interface' surface : whether single region
            const boolList& singleRegion() const
            {
                return singleRegion_;
            }

            //- Number of regions to remove during wrap checks
            const label& nRegionsWrapRemoval() const
            {
                return nRegionsWrapRemoval_;
            }

            //- Check whether unnamed surfaces contain a tri-surface mesh
            // or bounded primitive shape when running autoBlockMesh option
            bool checkForFiniteSurfaces() const;

            const PtrList<surfaceZonesInfo>& surfZones() const
            {
                return surfZones_;
            }

            //- From local region number to global region number
            const labelList& regionOffset() const
            {
                return regionOffset_;
            }

            //- From global region number to refinement level
            const labelList& minLevel() const
            {
                return minLevel_;
            }

            //- From global region number to refinement level
            const labelList& maxLevel() const
            {
                return maxLevel_;
            }

            //- From global region number to small gap refinement level
            const labelList& gapLevel() const
            {
                return gapLevel_;
            }

            //- All region refinement angles
            scalar minFeatureRefineAngle() const
            {
                scalar minFeatureAngle = 1.0;
                forAll(featureRefineAngle_, i)
                {
                    if (featureRefineAngle_[i] > SMALL)
                    {
                        minFeatureAngle =
                            min(minFeatureAngle , featureRefineAngle_[i]);
                    }
                }
                return minFeatureAngle;
            }

            //- From global region number to specification of gap and its
            //  refinement: 3 labels specifying
            //  - minimum wanted number of cells in the gap
            //  - minimum cell level when to start trying to detect gaps
            //  - maximum cell level to refine to (so do not detect gaps if
            //    cell >= maximum level)
            const List<FixedList<label, 3>>& extendedGapLevel() const
            {
                return extendedGapLevel_;
            }

            //- From global region number to side of surface to detect
            const List<volumeType>& extendedGapMode() const
            {
                return extendedGapMode_;
            }

            //- From global region number to perpendicular angle
            const scalarField& perpendicularAngle() const
            {
                return perpendicularAngle_;
            }

            //- From global region number to refinement level
            const labelList& proxLevelIncr() const
            {
                return proxLevelIncr_;
            }

            //- From global region number to proximity direction method
            const labelList& proxDir() const
            {
                return proxDir_;
            }

            //- From global region number to closure active
            const boolList& singleCellClosure() const
            {
                return singleCellClosure_;
            }

            //- From global region number to thin gap surfaces
            const boolList& thinGap() const
            {
                return thinGap_;
            }

            //- From global region number to thin gap surfaces
            const boolList& cornerCellRegions() const
            {
                return cornerCellRegions_;
            }

            //- From global region number to curvature angle in Radians
            const scalarField& curvature() const
            {
                return curvature_;
            }

            //- From global region number to maximum number gap cells
            const labelList& proxMaxCells() const
            {
                return proxMaxCells_;
            }


            //- From global region number to patch type
            const PtrList<dictionary>& patchInfo() const
            {
                return patchInfo_;
            }


        // Helper

            //- From surface and region on surface to global region
            label globalRegion(const label surfI, const label regionI) const
            {
                if (singleRegion()[surfI])
                {
                    return regionOffset_[surfI];
                }
                else
                {
                    return regionOffset_[surfI]+regionI;
                }
            }

            //- Min level for surface and region on surface
            label minLevel(const label surfI, const label regionI) const
            {
                return minLevel_[globalRegion(surfI, regionI)];
            }

            //- Max level for surface and region on surface
            label maxLevel(const label surfI, const label regionI) const
            {
                return maxLevel_[globalRegion(surfI, regionI)];
            }

            //- feature refinement angle for each region on surface
            scalar featureRefineAngle(const label surfI, const label regionI) const
            {
                return featureRefineAngle_[globalRegion(surfI, regionI)];
            }

            //- whether to refine and snap feature edges for each region on surface
            bool featureRefineOnly(const label surfI, const label regionI) const
            {
                return featureRefineOnly_[globalRegion(surfI, regionI)];
            }

            //- minimum baffle angle for each region on surface
            scalar minBaffleAngle(const label surfI, const label regionI) const
            {
                return minBaffleAngle_[globalRegion(surfI, regionI)];
            }

            //- Max level for surface and region on surface
            label maxProxLevel(const label surfI, const label regionI) const
            {
                return maxLevel_[globalRegion(surfI, regionI)]
                    + proxLevelIncr_[globalRegion(surfI,regionI)];
            }

            //- curvature refinement angle in Radians
            scalar curvature(const label surfI, const label regionI) const
            {
                return curvature_[globalRegion(surfI, regionI)];
            }

            //return surface to use in bounding box calculation
            labelList boundSurfaces() const;

            //return surface where cell-centroids can be moved
            labelList moveCentroidSurfaces() const;

            //return surface where single cell gaps can be closed
            labelList singleCellClosureSurfaces() const;

            //return surface where thin gap method selected
            labelList thinGapSurfaces() const;

            //return surface where corner cell keeping selected
            labelList cornerCellSurfaces() const;

            //return surface where wrap level specified
            labelList wrapLevelSurfaces() const;

            //return surfaces where curvature is active
            labelList curvatureSurfaces() const;

            //return surface where proximity is active
            labelList proximitySurfaces() const;

            label nRegions() const
            {
                return minLevel_.size();
            }

            //- Per surface the maximum extendedGapLevel over all its regions
            labelList maxGapLevel() const;

            //- Calculate curvature field with optional smoothing
            void setCurvatureFields(const label nSmoothIter);

            //- Calculate minLevelFields
            void setMinLevelFields
            (
                const autoPtr<shellSurfaces>& shellsPtr
            );

            //- Calculate bounding box of un-named surfaces
            boundBox getUnnamedBb();

            //- Calculate bounding box of all finite surfaces
            boundBox getFiniteSurfacesBb
            (
                const bool primitiveBounding
            );

            //- Calculate bounding box of all finite surfaces
            // in local coordinate direction
            void getFiniteSurfacesLocalBb
            (
                const coordinateSystem& coord,
                const bool primitiveBounding,
                point& minPt,
                point& maxPt
            ) const;

            //- Orientate named surfaces only
            void orient();

            //- Helper: orient (closed only) surfaces so keepPoint is outside.
            static void orientSurface
            (
                const point& keepPoint,
                triSurfaceMesh& surface
            );

            //- Whether to perform non-aligned baffle removal
            bool baffleRemoval() const;

        // Searching

            //- Find intersection of edge. Return -1 or first surface
            //  with higher (than currentLevel) minlevel.
            //  Return surface number and level.
            void findHigherIntersection
            (
                const shellSurfaces& shells,

                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                labelList& surfaces,
                labelList& surfaceLevel,
                const scalar edge0Len,

                const bool checkCurvature = false,
                const bool threaded = false
            ) const;

            //- Find all intersections of edge. Unsorted order.
            void findAllHigherIntersections
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                const labelList& globalRegionLevel, // level per surfregion
                List<vectorList>& surfaceNormal,
                labelListList& surfaceLevel
            ) const;

            void findAllIntersectedFaces
            (
                const pointField& start,
                const pointField& end,
                labelList& proxSurfaces,
                List<List<labelList>>& allHits
            ) const;

            //- Find all intersections of edge. Unsorted order.
            void findAllHigherIntersections
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level
                const labelList& globalRegionLevel, // level per surfregion

                List<pointList>& surfaceLocation,
                List<vectorList>& surfaceNormal,
                labelListList& surfaceLevel
            ) const;

            //- Find intersection nearest to the endpoints. surface1,2 are
            //  not indices into surfacesToTest but refinement surface indices.
            //  Returns surface, region on surface (so not global surface)
            //  and position on surface.
            void findNearestIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2
            ) const;

            //- findNearestIntersection but also get normals
            void findNearestIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                vectorField& normal1,

                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2,
                vectorField& normal2
            ) const;

            //- Find nearest (to start only) intersection of edge
            void findNearestIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                vectorField& normal
            ) const;

            //- Find nearest (to start only) intersection of edge
            void findNearestIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Used for debugging only: find intersection of edge.
            void findAnyIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            void findAnyIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            void findHighestIntersection
            (
                const pointField& start,
                const pointField& end,
                const labelList currentLevel,   // current cell refinement level

                labelList& maxLevel,
                labelList& hitSurface,
                List<pointIndexHit>& maxHit
            ) const;

            //- Find nearest point on surfaces.
            void findNearest
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find nearest point on surfaces. Return surface and region on
            //  surface (so not global surface)
            void findNearestRegion
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& hitSurface,
                labelList& hitRegion
            ) const;

            //- Find nearest point on surfaces. Return surface, region and
            //  normal on surface (so not global surface)
            void findNearestRegion
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& hitSurface,
                List<pointIndexHit>& hitInfo,
                labelList& hitRegion,
                vectorField& hitNormal
            ) const;

            //- Detect if a point is 'inside' (closed) surfaces.
            //  Returns -1 if not, returns first surface it is.
            void findInside
            (
                const labelList& surfacesToTest,
                const pointField& pt,
                labelList& insideSurfaces
            ) const;

            //- Find intersection nearest to the perturbed endpoints.
            //  surface1,2 are not indices into surfacesToTest but refinement
            //  surface indices. Returns surface, region on surface (so not
            //  global surface) and position on surface.
            void findNearestPerturbedIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,
                const scalar& tol,
                const bool& addedRays,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2
            ) const;

            //- Find intersection nearest to the perturbed endpoints
            //  and return surface normals.
            //  surface1,2 are not indices into surfacesToTest but refinement
            //  surface indices. Returns surface, region on surface (so not
            //  global surface) and position on surface.
            void findNearestPerturbedIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,
                const scalar& tol,
                const bool& addedRays,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                vectorField& normal1,

                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2,
                vectorField& normal2
            ) const;

            //- Find intersection of perturbed edges. Return -1 or first surface
            //  with higher (than currentLevel) minlevel.
            //  Return surface number and level.
            void findHigherPerturbedIntersection
            (
                const shellSurfaces& shells,

                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                labelList& surfaces,
                labelList& surfaceLevel,
                const scalar& tol,
                const bool& addedRays
            ) const;

            //- Find cell zones from named surfaces
            void findCellZone
            (
                const pointField& pts,
                const regionSplit& globalRegions,
                const labelList& surfaceToCellZone,
                labelList& cellToZone,
                label minRegionSize
            ) const;


            //- Find cell zones from named surfaces with perturbed points
            void findPerturbedCellZone
            (
                const labelList& testSurfaces,
                const pointField& pt,
                labelList& insideSurfaces
            ) const;

            //- Find surface normal and surface region for all hits
            void findAllProximityIntersections
            (
                const pointField& start,
                const pointField& end,

                List<List<pointIndexHit>>& allHits,
                List<List<vector>>& surfaceNormal,
                labelListList& surfaceRegion
             ) const;

            // Region wise searching

                //- Find nearest point on selected regions of surfaces.
                void findNearest
                (
                    const labelList& surfacesToTest,
                    const labelListList& regions,

                    const pointField& samples,
                    const scalarField& nearestDistSqr,

                    labelList& hitSurface,
                    List<pointIndexHit>& hitInfo
                ) const;

                //- Find nearest point on selected regions of surfaces.
                void findNearestRegion
                (
                    const labelList& surfacesToTest,
                    const labelListList& regions,

                    const pointField& samples,
                    const scalarField& nearestDistSqr,

                    labelList& hitSurface,
                    List<pointIndexHit>& hitInfo,
                    labelList& hitRegion,
                    vectorField& hitNormal
                ) const;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
