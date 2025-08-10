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
    (c) 2021 Esi Ltd.

Class
    Foam::VDBTree

\*---------------------------------------------------------------------------*/

#ifndef VDBTree_H
#define VDBTree_H

#include "meshes/treeBoundBox/treeBoundBox.H"
#include "meshes/primitiveShapes/objectHit/pointIndexHit.H"
#include "containers/Lists/FixedList/FixedList.H"
#include "db/IOstreams/IOstreams/Ostream.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "algorithms/indexedOctree/labelBits.H"
#include "containers/Lists/PackedList/PackedList.H"
#include "algorithms/indexedOctree/volumeType.H"

#include "include/OpenVDB.H"
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/RayIntersector.h>
#include <openvdb/math/Proximity.h> // for closestPointOnTriangleToPoint

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class VDBTree;
class triSurface;


/*---------------------------------------------------------------------------*\
                           Class VDBTree Declaration
\*---------------------------------------------------------------------------*/

class VDBTree
{
    typedef openvdb::tools::LinearSearchImpl<FloatGrid, /*iterations=*/2> SearchImplT;

    using RayIntersector = openvdb::tools::LevelSetRayIntersector<FloatGrid, SearchImplT>;
    //using RayIntersector = openvdb::tools::LevelSetRayIntersector<FloatGrid>;
    using Ray = RayIntersector::RayType;

    using FloatSampler = openvdb::tools::GridSampler<FloatGrid::ConstAccessor, openvdb::tools::PointSampler>; //nearest-neighbor (zero-order)
    using IndexSampler = openvdb::tools::GridSampler<IndexGrid::ConstAccessor, openvdb::tools::PointSampler>; //nearest-neighbor (zero-order)
    using Vec4ISampler = openvdb::tools::GridSampler<Vec4IGrid::ConstAccessor, openvdb::tools::PointSampler>; //nearest-neighbor (zero-order)
    using Vec3SSampler = openvdb::tools::GridSampler<Vec3SGrid::ConstAccessor, openvdb::tools::BoxSampler>; //nearest-neighbor (zero-order)
    //openvdb::tools::GridSampler<FloatGrid::ConstAccessor, openvdb::tools::BoxSampler>       //trilinear (first-order) interpolation
    //openvdb::tools::GridSampler<FloatGrid::ConstAccessor, openvdb::tools::QuadraticSampler> //triquadratic (second-order) interpolation

private:

    // Private data

        //- Underlying geometry
        const triSurface& surface_;

        //- Bounding box underlying geometry
        treeBoundBox bb_;

        // The offset to voxel-center points
        scalar offset_;

        // Background value of distance grid
        scalar backgroundSqr_;

        //- Grid Accessors
        FloatGrid::ConstAccessor distAcc_;
        IndexGrid::ConstAccessor idxAcc_;
        //Vec4IGrid::ConstAccessor idxAcc_;
        Vec3SGrid::ConstAccessor gradAcc_;

        //- Grid Samplers
        FloatSampler distSampler_;
        IndexSampler idxSampler_;
        //Vec4ISampler idxSampler_;
        Vec3SSampler gradSampler_;

        // Levelset ray intersector
        RayIntersector lsri_;

public:

    // Constructors

        //- Construct from components
        VDBTree
        (
            const triSurface& surface,
            const FloatGrid& distGrid,
            const IndexGrid& idxGrid,
            //const Vec4IGrid& idxGrid,
            const Vec3SGrid& gradGrid
        );

    // Member Functions

        // Access

            //- geometry bounding box
            const treeBoundBox& bb() const
            {
                return bb_;
            }

        // Queries

            pointIndexHit findNearest
            (
                const point& p,
                const scalar nearestDistSqr,
                point* levelsetPoint = nullptr,
                const scalar isoValue = 0,
                const bool pointOnLS = false
            ) const;

            //- Find nearest intersection of line between start and end.
            pointIndexHit findLine
            (
                const point& start,
                const point& end
            ) const;

            //- Find any intersection of line between start and end.
            pointIndexHit findLineAny
            (
                const point& start,
                const point& end
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            volumeType getVolumeType(const point&) const;

            //- Accelerated closest surface point queries for narrow band level sets
            // openVDB implementation, to be used only on Signed Distance Fields
            // without baffles
            // or Unsigned Distance Field with isoValue != 0
            /*inline*/ void closestPoints
            (
                const FloatGrid& distGrid,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>& info,
                const scalar isoValue = 0
            ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
