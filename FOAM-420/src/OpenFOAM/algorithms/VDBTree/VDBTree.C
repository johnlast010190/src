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

\*---------------------------------------------------------------------------*/

#include "VDBTree.H"
#include "meshes/primitiveShapes/line/linePointRef.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "memInfo/memInfo.H"
#include "triSurface/triSurface.H"

#include <tbb/tick_count.h>
#include <openvdb/util/Util.h>
#include <openvdb/tools/VolumeToSpheres.h> // for ClosestSurfacePoint

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::VDBTree::VDBTree
(
    const triSurface& surface,
    const FloatGrid& distGrid,
    const IndexGrid& idxGrid,
    //const Vec4IGrid& idxGrid,
    const Vec3SGrid& gradGrid
)
:
    surface_(surface),
    bb_(surface.localPoints()),
    offset_(-0.5*distGrid.transform().voxelSize()[0]),
    backgroundSqr_(sqr(scalar(distGrid.background()))),
    distAcc_(distGrid.getConstAccessor()),
    idxAcc_(idxGrid.getConstAccessor()),
    gradAcc_(gradGrid.getConstAccessor()),
    distSampler_(distAcc_, distGrid.transform()),
    idxSampler_(idxAcc_, idxGrid.transform()),
    gradSampler_(gradAcc_, gradGrid.transform()),
    lsri_(distGrid)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointIndexHit Foam::VDBTree::findNearest
(
    const point& p,
    const scalar nearestDistSqr,
    point* levelsetPoint, /*= nullptr*/
    const scalar isoValue, /*= 0*/
    const bool pointOnLS /*= false*/
) const
{
    assert(isoValue >= 0);

    //TODO consider when isoValue != 0
    if (nearestDistSqr > backgroundSqr_)
    {
        WarningInFunction
                << "nearestDistSqr (" << nearestDistSqr
                << ") larger than distance grid backgroundSqr (" << backgroundSqr_
                << "). Skipping findNearest for point " << p
                << ". Please specify nearDist < " << sqrt(backgroundSqr_)
                << endl;

        return pointIndexHit(false, Zero, -1);
    }

    // sample offset to voxel-center points
    const openvdb::Vec3d sample(p.x()+offset_, p.y()+offset_, p.z()+offset_);

    // Compute the value of the grid at a location in world space.
    FloatGrid::ValueType d = distSampler_.wsSample(sample);
    FloatGrid::ValueType dIso = d - isoValue;

    //TODO possible optimization OpenVDB grid have distSqr field (instead of dist)
    if (sqr(dIso) > nearestDistSqr)
    {
        return pointIndexHit(false, Zero, -1);
    }
    else
    {
        Vec3SGrid::ValueType grad = gradSampler_.wsSample(sample);

        openvdb::Vec3d cp = (dIso > 0 ? sample - grad*dIso : sample + grad*dIso) ;

        // offset back to world space
        cp -= offset_;

        if (levelsetPoint)
        {
            *levelsetPoint = point(cp.x(), cp.y(), cp.z());
        }
        if (isoValue > 0 || pointOnLS)
        {
            return pointIndexHit(true, *levelsetPoint, -1);
        }

        IndexGrid::ValueType nearestShapeI = idxSampler_.wsSample(cp);
        //Vec4IGrid::ValueType nearestShapes = idxSampler_.wsSample(cp);

        //std::cout<<"nearestShapes " << nearestShapeI <<std::endl;

        point nearestPoint(GREAT);
        label nearestShape = -1;

        //for (size_t i = 0; i < 4; i++)
        {
            //if (nearestShapes[i] == openvdb::Int32(openvdb::util::INVALID_IDX)) break;

            //label nearestShapeI = nearestShapes[i];

            nearestShapeI = min(max(0, nearestShapeI), surface_.size());

            const UList<point>& meshPoints = surface_.points();

            const face& f = surface_[nearestShapeI];

            //tbb::tick_count t0 = tbb::tick_count::now();
            //label nearType = -1, nearLabel = -1;

            //pointHit nearestPointH =
            //    triPointRef
            //    (
            //        meshPoints[f[0]],
            //        meshPoints[f[1]],
            //        meshPoints[f[2]]
            //    ).nearestPointClassify(p, nearType, nearLabel);

            //point nearestPoint = nearestPointH.rawPoint();

            //tbb::tick_count t1 = tbb::tick_count::now();

            openvdb::Vec3d uvw;
            openvdb::Vec3d a(meshPoints[f[0]].x(), meshPoints[f[0]].y(), meshPoints[f[0]].z());
            openvdb::Vec3d b(meshPoints[f[1]].x(), meshPoints[f[1]].y(), meshPoints[f[1]].z());
            openvdb::Vec3d c(meshPoints[f[2]].x(), meshPoints[f[2]].y(), meshPoints[f[2]].z());

            cp =
                openvdb::math::closestPointOnTriangleToPoint
                (
                    a, b, c, (sample-offset_), uvw
                );

            point nearestPointVDB(cp.x(), cp.y(), cp.z());

            //tbb::tick_count t2 = tbb::tick_count::now();

            //Info<<"nearestPoint   : " << nearestPoint    << " (in " << (t1-t0).seconds() << " sec)" << nl
            //    <<"nearestPointVDB: " << nearestPointVDB << " (in " << (t2-t1).seconds() << " sec)" << nl
            //    <<endl;
            if (mag(nearestPointVDB - p) < mag(nearestPoint - p))
            {
                nearestPoint = nearestPointVDB;
                nearestShape = nearestShapeI;
            }
        }
        return pointIndexHit(true, nearestPoint, nearestShape);
    }
} // findNearest


Foam::pointIndexHit Foam::VDBTree::findLine
(
    const point& start,
    const point& end
) const
{
    vector lineVec(end - start);
    const scalar lineLength = mag(lineVec);
    lineVec.normalise();

    Ray::Vec3T eye(start.x()+offset_, start.y()+offset_, start.z()+offset_);
    Ray::Vec3T dir(lineVec.x(), lineVec.y(), lineVec.z());

    Ray ray(eye, dir, SMALL, lineLength);
    //Ray ray(eye, dir);

    openvdb::Vec3d xyz(0);

    if (!lsri_.intersectsWS(ray, xyz))
    {
        return pointIndexHit(false, Zero, -1);
    }

    scalar distance = (xyz - eye).length();

    if (distance <= lineLength)
    {
        //std::cout<<"Ray :" << ray <<std::endl;

        //std::cout<<"xyz " << xyz <<std::endl;
        IndexGrid::ValueType nearestShapeI = idxSampler_.wsSample(xyz);
        //Vec4IGrid::ValueType nearestShapes = idxSampler_.wsSample(xyz);

        ////TODO check intersection on triangle ID
        //label nearestShapeI = nearestShapes[0];

        const point intersection(xyz.x()-offset_, xyz.y()-offset_, xyz.z()-offset_);

        return pointIndexHit(true, intersection, nearestShapeI);
    }
    else
    {
        return pointIndexHit(false, Zero, -1);
    }
} // findLine


Foam::pointIndexHit Foam::VDBTree::findLineAny
(
    const point& start,
    const point& end
) const
{
    return findLine
    (
        start,
        end
    );
}

Foam::volumeType Foam::VDBTree::getVolumeType
(
    const point& p
) const
{
    // Compute the value of the grid at a location in world space.
    const openvdb::Vec3d sample(p.x(), p.y(), p.z());

    FloatGrid::ValueType d = distSampler_.wsSample(sample);

    if (d > 0.0)
    {
        return volumeType::OUTSIDE;
    }
    else
    {
        return volumeType::INSIDE;
    }
} // getVolumeType


void Foam::VDBTree::closestPoints
(
    const FloatGrid& distGrid,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    List<pointIndexHit>& info,
    const scalar isoValue /*= 0*/
) const
{
    std::vector<openvdb::Vec3d> tmpPoints(samples.size());

    for (size_t n = 0, N = tmpPoints.size(); n < N; ++n)
    {
        tmpPoints[n][0] = samples[n].x() + offset_;
        tmpPoints[n][1] = samples[n].y() + offset_;
        tmpPoints[n][2] = samples[n].z() + offset_;
    }

    auto closestPoint =
        openvdb::tools::ClosestSurfacePoint<FloatGrid>::create(distGrid, isoValue);

    if (!closestPoint)
    {
        WarningInFunction
                << "Cannot create openvdb ClosestSurfacePoint"
                << endl;
        return;
    }

    std::vector<float> tmpDistances;

    closestPoint->searchAndReplace(tmpPoints, tmpDistances);

    for (size_t n = 0, N = tmpDistances.size(); n < N; ++n)
    {
        if (sqr(tmpDistances[n]) <= nearestDistSqr[n])
        {
            IndexGrid::ValueType nearestShape =
                idxSampler_.wsSample(tmpPoints[n]);

            point nearestPoint
            (
                tmpPoints[n].x() - offset_,
                tmpPoints[n].y() - offset_,
                tmpPoints[n].z() - offset_
            );

            info[n] = pointIndexHit(true, nearestPoint, nearestShape);
        }
        else
        {
            info[n] = pointIndexHit(false, Zero, -1);
        }
    }
} // closestPoints

// ************************************************************************* //
