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
    (c) 2015-2017 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "shellSurfaces/shellSurfaces.H"
#include "searchableSurfaces/searchableSurface/searchableSurface.H"
#include "meshes/boundBox/boundBox.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "meshes/primitiveShapes/objectHit/pointIndexHit.H"
#include "triSurface/orientedSurface/orientedSurface.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "algorithms/indexedOctree/volumeType.H"
#include "regionSplit/regionSplit.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


namespace Foam
{

template<>
const char*
NamedEnum<shellSurfaces::refineMode, 4>::
names[] =
{
    "inside",
    "outside",
    "distance",
    "direction"
};

const NamedEnum<shellSurfaces::refineMode, 4> shellSurfaces::refineModeNames_;

template<>
const char*
NamedEnum<shellSurfaces::refineType, 2>::
names[] =
{
    "iso",
    "aniso"
};

const NamedEnum<shellSurfaces::refineType, 2> shellSurfaces::refineTypeNames_;

} // End namespace Foam



// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary Foam::shellSurfaces::getCoordinateDict
(
    const dictionary& dict
)
{
    if (dict.found("coordinates"))
    {
        return dict.subDict("coordinates");
    }
    else
    {
        dictionary defaultDict("coordinates");

        defaultDict.add("origin",vector::zero,true);
        defaultDict.add(word("rotation"), dictionary(), false);
        dictionary& rotDict = defaultDict.subDict("rotation");
        rotDict.add("type","axesRotation",true);
        rotDict.add("e1",vector(1,0,0),true);
        rotDict.add("e2",vector(0,1,0),true);

        return defaultDict;
    }
}


void Foam::shellSurfaces::setAndCheckLevels
(
    const label shellI,
    const List<Tuple2<Tuple2<point,point>,label>>& distLevels
)
{
     // Extract information into separate distance and level
    distances_[shellI].setSize(distLevels.size());
    levels_[shellI].setSize(distLevels.size());

    forAll(distLevels, j)
    {
        point nPt = distLevels[j].first().first();
        point pPt = distLevels[j].first().second();
        label level = distLevels[j].second();

        if (cmptMin(nPt) == 0 || cmptMin(pPt) == 0)
        {
            const searchableSurface& shell = allGeometry_[shells_[shellI]];
            WarningInFunction
                << "Zero values for directional refinement values "
                << nPt << " "<< pPt << " with level " << level
                << " for shell " << shell.name()
                << " no volume refinement will be active for this level"
                << endl;
        }

        distances_[shellI][j] = Tuple2<point,point>
        (
            nPt,
            pPt
        );
        levels_[shellI][j] = labelVector(level);
    }
}


template<class T> void Foam::shellSurfaces::setAndCheckLevels
(
    const label shellI,
    const List<Tuple2<scalar,T>>& distLevels
)
{

    if (modes_[shellI] != DISTANCE && distLevels.size() != 1)
    {
        FatalErrorInFunction
            << "For refinement mode "
            << refineModeNames_[modes_[shellI]]
            << " specify only one distance+level."
            << " (its distance gets discarded)"
            << exit(FatalError);
    }
    // Extract information into separate distance and level
    distances_[shellI].setSize(distLevels.size());
    levels_[shellI].setSize(distLevels.size());

    forAll(distLevels, j)
    {
        distances_[shellI][j] = Tuple2<point,point>
        (
            distLevels[j].first(),
            distLevels[j].first()
        );
        levels_[shellI][j] = labelVector(distLevels[j].second());

        // Check in incremental order
        if (j > 0)
        {
            if
            (
                (distances_[shellI][j] <= distances_[shellI][j-1])
             || (levels_[shellI][j] > levels_[shellI][j-1])
            )
            {
                FatalErrorInFunction
                    << "For refinement mode "
                    << refineModeNames_[modes_[shellI]]
                    << " : Refinement should be specified in order"
                    << " of increasing distance"
                    << " (and decreasing refinement level)." << endl
                    << "Distance:" << distances_[shellI][j]
                    << " refinementLevel:" << levels_[shellI][j]
                    << exit(FatalError);
            }
        }
    }

    const searchableSurface& shell = allGeometry_[shells_[shellI]];

    if (modes_[shellI] == DISTANCE)
    {
        Info<< "Refinement level according to distance to "
            << shell.name() << endl;
        forAll(levels_[shellI], j)
        {
            Info<< "    level " << levels_[shellI][j]
                << " for all cells within " << distances_[shellI][j]
                << " metre." << endl;
        }
    }
    else
    {
        //Support has been added for non-manifold surfaces, so
        //following check is no longer required
/*
        if (!allGeometry_[shells_[shellI]].hasVolumeType())
        {
            FatalErrorInFunction
                << "Shell " << shell.name()
                << " does not support testing for "
                << refineModeNames_[modes_[shellI]] << endl
                << "Probably it is not closed."
                << exit(FatalError);
        }
*/
        if (modes_[shellI] == INSIDE)
        {
            Info<< "Refinement level " << levels_[shellI][0]
                << " for all cells inside " << shell.name() << endl;
        }
        else
        {
            Info<< "Refinement level " << levels_[shellI][0]
                << " for all cells outside " << shell.name() << endl;
        }
    }
}


void Foam::shellSurfaces::checkGapLevels
(
    const dictionary& shellDict,
    const label shellI,
    const List<FixedList<label, 3>>& levels
)
{
    const searchableSurface& shell = allGeometry_[shells_[shellI]];

    forAll(levels, regionI)
    {
        const FixedList<label, 3>& info = levels[regionI];

        if (info[2] > 0)
        {
            if (modes_[shellI] == DISTANCE)
            {
                FatalIOErrorInFunction(shellDict)
                    << "'gapLevel' specification cannot be used with mode "
                    << refineModeNames_[DISTANCE]
                    << " for shell " << shell.name()
                    << exit(FatalIOError);
            }
        }
    }

    // Hardcode for region 0
    if (levels[0][0] > 0)
    {
        Info<< "Refinement level up to " << levels[0][2]
            << " for all cells in gaps for shell "
            << shell.name() << endl;
    }
}



// Specifically orient triSurfaces using a calculated point outside.
// Done since quite often triSurfaces not of consistent orientation which
// is (currently) necessary for sideness calculation
void Foam::shellSurfaces::orient()
{
    // Determine outside point.
    boundBox overallBb = boundBox::invertedBox;

    bool hasSurface = false;

    forAll(shells_, shellI)
    {
        const searchableSurface& s = allGeometry_[shells_[shellI]];

        if (modes_[shellI] != DISTANCE && isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& shell = refCast<const triSurfaceMesh>(s);

            if (shell.triSurface::size())
            {
                hasSurface = true;
                // Assume surface is compact!
                overallBb.add(shell.points());
            }
        }
    }

    if (hasSurface)
    {
        const point outsidePt = overallBb.max() + overallBb.span();

        //Info<< "Using point " << outsidePt << " to orient shells" << endl;

        forAll(shells_, shellI)
        {
            const searchableSurface& s = allGeometry_[shells_[shellI]];

            if (modes_[shellI] != DISTANCE && isA<triSurfaceMesh>(s))
            {
                triSurfaceMesh& shell = const_cast<triSurfaceMesh&>
                (
                    refCast<const triSurfaceMesh>(s)
                );

                // Flip surface so outsidePt is outside.
                bool anyFlipped = orientedSurface::orient
                (
                    shell,
                    outsidePt,
                    true
                );

                if (anyFlipped)
                {
                    // orientedSurface will have done a clearOut of the surface.
                    // we could do a clearout of the triSurfaceMeshes::trees()
                    // but these aren't affected by orientation
                    // (except for cached
                    // sideness which should not be set at this point.
                    // !!Should check!)

                    Info<< "shellSurfaces : Flipped orientation of surface "
                        << s.name()
                        << " so point " << outsidePt << " is outside." << endl;
                }
            }
        }
    }
}


// Find maximum level of a shell.
void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const label shellI,
    labelList& maxLevel,
    const bool threaded /*= false*/
) const
{
    labelList minCmptLevel(levels_[shellI].size());
    forAll(levels_[shellI],levelI)
    {
        minCmptLevel[levelI] = cmptMin(levels_[shellI][levelI]);
    }

    if (modes_[shellI] == DIRECTION)
    {
         // Direction mode.
        const List<Tuple2<point,point>>& distances = distances_[shellI];

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        //TODO threaded
        // parallel_reduce
        forAll(maxLevel, pointi)
        {
            scalar maxSearchRadius = 0;
            forAllReverse(minCmptLevel, levelI)
            {
                if (minCmptLevel[levelI] > maxLevel[pointi])
                {
                    scalar maxLevelDist = max
                    (
                        magSqr(distances[levelI].first()),
                        magSqr(distances[levelI].second())
                    );

                    maxSearchRadius = max(maxSearchRadius,maxLevelDist);
                }
            }

            if (maxSearchRadius > 0)
            {
                candidates[candidateI] = pt[pointi];
                candidateMap[candidateI] = pointi;
                candidateDistSqr[candidateI] = maxSearchRadius;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo,
            threaded
        );

        const coordinateSystem& cSys = coord_[shellI]();

        // Update maxLevel
        if (threaded && nearInfo.size() > 2048)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    nearInfo.size() / tbb::this_task_arena::max_concurrency(),
                    1024
                );

            tbb::parallel_for
            (
                tbb::blocked_range<size_t>(0, nearInfo.size(), grainSize),
                [&](const tbb::blocked_range<size_t>& r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        if (nearInfo[i].hit())
                        {
                            label pointi = candidateMap[i];
                            forAllReverse(minCmptLevel, levelI)
                            {
                                if (minCmptLevel[levelI] > maxLevel[pointi])
                                {
                                    vector hitVec = candidates[i]-nearInfo[i].hitPoint();
                                    hitVec = cSys.localVector(hitVec);

                                    vector influence = vector::zero;

                                    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                                    {
                                        scalar denom =
                                        (
                                            hitVec[cmpt] > 0 ?
                                            distances[levelI].second()[cmpt] + SMALL
                                            : distances[levelI].first()[cmpt] + SMALL
                                         );

                                        influence[cmpt] = hitVec[cmpt] / denom;
                                    }

                                    if
                                    (
                                        mag(influence.x()) < 1.0
                                        && mag(influence.y()) < 1.0
                                        &&  mag(influence.z()) < 1.0
                                    )
                                    {
                                        maxLevel[pointi] = minCmptLevel[levelI];
                                        break;
                                    }
                                }
                            }
                        }
                    }
                },
                tbb::simple_partitioner()
            );
#endif
        }
        else
        {
            forAll(nearInfo, i)
            {
                if (nearInfo[i].hit())
                {
                    label pointi = candidateMap[i];
                    forAllReverse(minCmptLevel, levelI)
                    {
                        if (minCmptLevel[levelI] > maxLevel[pointi])
                        {
                            vector hitVec = candidates[i]-nearInfo[i].hitPoint();
                            hitVec = cSys.localVector(hitVec);

                            vector influence = vector::zero;

                            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                            {
                                scalar denom =
                                (
                                    hitVec[cmpt] > 0 ?
                                    distances[levelI].second()[cmpt] + SMALL
                                    : distances[levelI].first()[cmpt] + SMALL
                                 );

                                influence[cmpt] = hitVec[cmpt] / denom;
                            }

                            if
                            (
                                mag(influence.x()) < 1.0
                                && mag(influence.y()) < 1.0
                                &&  mag(influence.z()) < 1.0
                            )
                            {
                                maxLevel[pointi] = minCmptLevel[levelI];
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    else if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.

        const scalarField distances = shellIsoDistances(shellI);

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        //TODO threaded
        // parallel_reduce
        forAll(maxLevel, pointi)
        {
            forAllReverse(minCmptLevel, levelI)
            {
                if (minCmptLevel[levelI] > maxLevel[pointi])
                {
                    candidates[candidateI] = pt[pointi];
                    candidateMap[candidateI] = pointi;
                    candidateDistSqr[candidateI] = sqr(distances[levelI]);
                    candidateI++;
                    break;
                }
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo,
            threaded
        );

        // Update maxLevel
        if (threaded && nearInfo.size() > 2048)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    nearInfo.size() / tbb::this_task_arena::max_concurrency(),
                    1024
                );

            tbb::parallel_for
            (
                tbb::blocked_range<size_t>(0, nearInfo.size(), grainSize),
                [&](const tbb::blocked_range<size_t>& r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        if (nearInfo[i].hit())
                        {
                            // Check which level it actually is in.
                            label minDistI = findLower
                            (
                                distances,
                                mag(nearInfo[i].hitPoint()-candidates[i])
                            );

                            label pointi = candidateMap[i];

                            // pt is inbetween shell[minDistI] and shell[minDistI+1]
                            maxLevel[pointi] = minCmptLevel[minDistI+1];
                        }
                    }
                },
                tbb::simple_partitioner()
            );
#endif
        }
        else
        {
            forAll(nearInfo, i)
            {
                if (nearInfo[i].hit())
                {
                    // Check which level it actually is in.
                    label minDistI = findLower
                    (
                        distances,
                        mag(nearInfo[i].hitPoint()-candidates[i])
                    );

                    label pointi = candidateMap[i];

                    // pt is inbetween shell[minDistI] and shell[minDistI+1]
                    maxLevel[pointi] = minCmptLevel[minDistI+1];
                }
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        label candidateI = 0;

        //TODO threaded
        // parallel_reduce
        forAll(maxLevel, pointi)
        {
            if (minCmptLevel[0] > maxLevel[pointi])
            {
                candidates[candidateI] = pt[pointi];
                candidateMap[candidateI] = pointi;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<volumeType> volType;
        allGeometry_[shells_[shellI]].getVolumeType
        (
            candidates,
            volType,
            threaded
        );

        if (threaded && volType.size() > 2048)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    volType.size() / tbb::this_task_arena::max_concurrency(),
                    1024
                );

            tbb::parallel_for
            (
                tbb::blocked_range<size_t>(0, volType.size(), grainSize),
                [&](const tbb::blocked_range<size_t>& r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        label pointi = candidateMap[i];

                        if
                        (
                            (
                                modes_[shellI] == INSIDE
                             && volType[i] == volumeType::INSIDE
                            )
                         || (
                                modes_[shellI] == OUTSIDE
                             && volType[i] == volumeType::OUTSIDE
                            )
                        )
                        {
                            maxLevel[pointi] = minCmptLevel[0];
                        }
                    }
                },
                tbb::simple_partitioner()
            );
#endif
        }
        else
        {
            forAll(volType, i)
            {
                label pointi = candidateMap[i];

                if
                (
                    (
                        modes_[shellI] == INSIDE
                     && volType[i] == volumeType::INSIDE
                    )
                 || (
                        modes_[shellI] == OUTSIDE
                     && volType[i] == volumeType::OUTSIDE
                    )
                )
                {
                    maxLevel[pointi] = minCmptLevel[0];
                }
            }
        }
    }
} // findHigherLevel


// Find maximum level of a shell.
void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const label shellI,
    const labelList& regions,
    labelList& maxLevel,
    labelList& minCmptLevel,
    labelList& shellCheck,
    const refineType refType,
    const direction cmpt
) const
{

    const List<labelVector>& levels = levels_[shellI];
    if (modes_[shellI] == DIRECTION)
    {
         // Direction mode.
        const List<Tuple2<point,point>>& distances = distances_[shellI];

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(maxLevel, pointi)
        {
            scalar maxSearchRadius = 0;
            forAllReverse(levels, levelI)
            {
                if (levels[levelI][cmpt] > maxLevel[pointi])
                {
                    scalar maxLevelDist = max
                    (
                        magSqr(distances[levelI].first()),
                        magSqr(distances[levelI].second())
                    );

                    maxSearchRadius = max(maxSearchRadius,maxLevelDist);
                }
            }

            if (maxSearchRadius > 0)
            {
                candidates[candidateI] = pt[pointi];
                candidateMap[candidateI] = pointi;
                candidateDistSqr[candidateI] = maxSearchRadius;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        const coordinateSystem& cSys = coord_[shellI]();

        // Update maxLevel
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                label pointI = candidateMap[i];
                forAllReverse(levels, levelI)
                {
                    if (levels[levelI][cmpt] > maxLevel[pointI])
                    {
                        vector hitVec = candidates[i]-nearInfo[i].hitPoint();
                        hitVec = cSys.localVector(hitVec);

                        vector influence = vector::zero;

                        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                        {
                            scalar denom =
                            (
                                hitVec[cmpt] > 0 ?
                                distances[levelI].second()[cmpt] + SMALL
                                : distances[levelI].first()[cmpt] + SMALL
                             );

                            influence[cmpt] = hitVec[cmpt] / denom;
                        }

                        if
                        (
                            mag(influence.x()) < 1.0
                            && mag(influence.y()) < 1.0
                            &&  mag(influence.z()) < 1.0
                        )
                        {
                            maxLevel[pointI] = levels[levelI][cmpt];
                            minCmptLevel[pointI] = cmptMin(levels[levelI]);
                            shellCheck[pointI] = shellI;
                            break;
                        }
                    }
                }
            }
        }
    }
    else if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.

        const scalarField distances = shellIsoDistances(shellI);

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(maxLevel, pointI)
        {
            forAllReverse(levels, levelI)
            {
                if (levels[levelI][cmpt] > maxLevel[pointI])
                {
                    candidates[candidateI] = pt[pointI];
                    candidateMap[candidateI] = pointI;
                    candidateDistSqr[candidateI] = sqr(distances[levelI]);
                    candidateI++;
                    break;
                }
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        // Update maxLevel
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                // Check which level it actually is in.
                label minDistI = findLower
                (
                    distances,
                    mag(nearInfo[i].hitPoint()-candidates[i])
                );

                label pointI = candidateMap[i];

                // pt is inbetween shell[minDistI] and shell[minDistI+1]
                maxLevel[pointI] = levels[minDistI+1][cmpt];
                minCmptLevel[pointI] = cmptMin(levels[minDistI+1]);
                shellCheck[pointI] = shellI;
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateRegions(pt.size());
        labelList candidateMap(pt.size());
        label candidateI = 0;
        forAll(maxLevel, pointI)
        {
            if (levels[0][cmpt] > maxLevel[pointI])
            {
                candidates[candidateI] = pt[pointI];
                candidateRegions[candidateI] = regions[pointI];
                candidateMap[candidateI] = pointI;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateRegions.setSize(candidateI);
        candidateMap.setSize(candidateI);

        const searchableSurface& s = allGeometry_[shells_[shellI]];

        if (isA<triSurfaceMesh>(s))
        {
            boundBox overallBb
            (
                point(GREAT, GREAT, GREAT),
                point(-GREAT, -GREAT, -GREAT)
            );

            triSurfaceMesh& surf = const_cast<triSurfaceMesh&>
            (
                refCast<const triSurfaceMesh>(s)
            );
            tmp<pointField> tpoints(surf.points());
            const pointField& points = tpoints();
            boundBox surfBb(points, false);
            overallBb.min() = min(overallBb.min(), surfBb.min());
            overallBb.max() = max(overallBb.max(), surfBb.max());

            boolList visited(candidates.size(), false);
            // choose 3 end points to make sure initial insideness
            // check is robust
            point midBbPt = 0.5*(overallBb.max() + overallBb.min());
            scalar dr = overallBb.mag();

            pointField endPoints(3);
            endPoints[0] = midBbPt + point(dr, 0., 0.);
            endPoints[1] = midBbPt + point(0., dr, 0.);
            endPoints[2] = midBbPt + point(0., 0., dr);

            pointField start(1);
            pointField end(1);

            forAll(candidates, i)
            {
                if (!visited[i])
                {
                    start[0] = candidates[i];
                    label nIn = 0;
                    label nOut = 0;
                    forAll(endPoints, pointI)
                    {
                        end[0] = endPoints[pointI];
                        List<List<pointIndexHit>> hitInfo(1);
                        allGeometry_[shells_[shellI]].findLineAll
                        (
                            start,
                            end,
                            hitInfo
                        );

                        label nHits = hitInfo[0].size();
                        //Filter out coincident hit points
                        if (hitInfo[0].size() > 1)
                        {
                            for (int hitI = 1;hitI < hitInfo[0].size();hitI++)
                            {
                                scalar hitNbrDist = mag
                                (
                                    hitInfo[0][hitI].hitPoint()
                                    - hitInfo[0][hitI-1].hitPoint()
                                );
                                if (hitNbrDist < SMALL)
                                {
                                    nHits--;
                                }
                            }
                        }

                        if ((nHits % 2) == 0)
                        {
                            nOut++;
                        }
                        else
                        {
                            nIn++;
                        }
                    }
                    label regionI = candidateRegions[i];

                    bool outside = nOut > nIn;

                    if (outside)
                    {
                        if (modes_[shellI] == OUTSIDE)
                        {
                            visited[i] = true;
                            maxLevel[candidateMap[i]] = levels[0][cmpt];
                            minCmptLevel[candidateMap[i]] = cmptMin(levels[0]);
                            shellCheck[candidateMap[i]] = shellI;
                        }
                        for (label j=i+1; j < candidates.size(); j++)
                        {
                            if (candidateRegions[j] == regionI)
                            {
                                visited[j] = true;
                                if (modes_[shellI] == OUTSIDE)
                                {
                                    maxLevel[candidateMap[j]] = levels[0][cmpt];
                                    minCmptLevel[candidateMap[i]] =
                                        cmptMin(levels[0]);
                                    shellCheck[candidateMap[j]] = shellI;
                                }
                            }
                        }
                    }
                    else
                    {
                        if (modes_[shellI] == INSIDE)
                        {
                            visited[i] = true;
                            maxLevel[candidateMap[i]] = levels[0][cmpt];
                            minCmptLevel[candidateMap[i]] = cmptMin(levels[0]);
                            shellCheck[candidateMap[i]] = shellI;
                        }
                        for (label j=i+1; j < candidates.size(); j++)
                        {
                            if (candidateRegions[j] == regionI)
                            {
                                visited[j] = true;
                                if (modes_[shellI] == INSIDE)
                                {
                                    maxLevel[candidateMap[j]] = levels[0][cmpt];
                                    minCmptLevel[candidateMap[i]] =
                                        cmptMin(levels[0]);
                                    shellCheck[candidateMap[j]] = shellI;
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
        // Do the expensive nearest test only for the candidate points.
            List<volumeType> volType;
            allGeometry_[shells_[shellI]].getVolumeType(candidates, volType);

            forAll(volType, i)
            {
                label pointI = candidateMap[i];

                if
                (
                    (
                        modes_[shellI] == INSIDE
                     && volType[i] == volumeType::INSIDE
                    )
                 || (
                        modes_[shellI] == OUTSIDE
                     && volType[i] == volumeType::OUTSIDE
                     )
                )
                {
                    maxLevel[pointI] = levels[0][cmpt];
                    minCmptLevel[candidateMap[i]] = cmptMin(levels[0]);
                    shellCheck[pointI] = shellI;
                }
            }
        }
    }
}


void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    const label shellI,
    labelList& gapShell,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    //TBD: hardcoded for region 0 information
    const FixedList<label, 3>& info = extendedGapLevel_[shellI][0];
    volumeType mode = extendedGapMode_[shellI][0];

    if (info[2] == 0)
    {
        return;
    }


    // Collect all those points that have a current maxLevel less than the
    // shell.

    labelList candidateMap(pt.size());
    label candidateI = 0;

    forAll(ptLevel, pointI)
    {
        if (ptLevel[pointI] >= info[1] && ptLevel[pointI] < info[2])
        {
            candidateMap[candidateI++] = pointI;
        }
    }
    candidateMap.setSize(candidateI);

    // Do the expensive nearest test only for the candidate points.
    List<volumeType> volType;
    allGeometry_[shells_[shellI]].getVolumeType
    (
        pointField(pt, candidateMap),
        volType
    );

    forAll(volType, i)
    {
        label pointI = candidateMap[i];

        bool isInside = (volType[i] == volumeType::INSIDE);

        if
        (
            (
                (modes_[shellI] == INSIDE && isInside)
             || (modes_[shellI] == OUTSIDE && !isInside)
            )
         && info[2] > gapInfo[pointI][2]
        )
        {
            gapShell[pointI] = shellI;
            gapInfo[pointI] = info;
            gapMode[pointI] = mode;
        }
    }
}


void Foam::shellSurfaces::findLevel
(
    const pointField& pt,
    const label shellI,
    labelList& minLevel,
    labelList& shell,
    const refineType refType,
    const direction cmpt
) const
{

    if (types_[shellI] != refType)
    {
        return;
    }

    const List<labelVector>& levels = levels_[shellI];
    if (modes_[shellI] == DIRECTION)
    {
         // Direction mode.
        const List<Tuple2<point,point>>& distances = distances_[shellI];

        // Collect all those points that have a current maxLevel less than
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(minLevel, pointi)
        {
            scalar maxSearchRadius = 0;
            forAllReverse(levels, levelI)
            {
                if (levels[levelI][cmpt] <= minLevel[pointi])
                {
                    scalar maxLevelDist = max
                    (
                        magSqr(distances[levelI].first()),
                        magSqr(distances[levelI].second())
                    );

                    maxSearchRadius = max(maxSearchRadius,maxLevelDist);
                }
            }
            if (maxSearchRadius > 0)
            {
                candidates[candidateI] = pt[pointi];
                candidateMap[candidateI] = pointi;
                candidateDistSqr[candidateI] = maxSearchRadius;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        const coordinateSystem& cSys = coord_[shellI]();
        // Update maxLevel
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                label pointI = candidateMap[i];
                forAllReverse(levels, levelI)
                {
                    if (levels[levelI][cmpt] <= minLevel[pointI])
                    {
                        vector hitVec = candidates[i]-nearInfo[i].hitPoint();
                        hitVec = cSys.localVector(hitVec);

                        vector influence = vector::zero;

                        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
                        {
                            scalar denom =
                            (
                                hitVec[cmpt] > 0 ?
                                distances[levelI].second()[cmpt] + SMALL
                                : distances[levelI].first()[cmpt] + SMALL
                             );

                            influence[cmpt] = hitVec[cmpt] / denom;
                        }

                        if
                        (
                            mag(influence.x()) < 1.0
                            && mag(influence.y()) < 1.0
                            &&  mag(influence.z()) < 1.0
                        )
                        {
                            shell[pointI] = shellI;
                            minLevel[pointI] = levels[levelI][cmpt];
                            break;
                        }
                    }
                }
            }
        }
    }
    else if (modes_[shellI] == DISTANCE)
    {
        // Distance mode.

        const scalarField distances = shellIsoDistances(shellI);

        // Collect all those points that have a current level equal/greater
        // (any of) the shell. Also collect the furthest distance allowable
        // to any shell with a higher level.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        scalarField candidateDistSqr(pt.size());
        label candidateI = 0;

        forAll(shell, pointI)
        {
            if (shell[pointI] == -1)
            {
                forAllReverse(levels, levelI)
                {
                    if (levels[levelI][cmpt] <= minLevel[pointI])
                    {
                        candidates[candidateI] = pt[pointI];
                        candidateMap[candidateI] = pointI;
                        candidateDistSqr[candidateI] = sqr(distances[levelI]);
                        candidateI++;
                        break;
                    }
                }
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);
        candidateDistSqr.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<pointIndexHit> nearInfo;
        allGeometry_[shells_[shellI]].findNearest
        (
            candidates,
            candidateDistSqr,
            nearInfo
        );

        // Update maxLevel
        forAll(nearInfo, i)
        {
            if (nearInfo[i].hit())
            {
                // Check which level it actually is in.
                label minDistI = findLower
                (
                    distances,
                    mag(nearInfo[i].hitPoint()-candidates[i])
                );

                label pointI = candidateMap[i];

                // pt is inbetween shell[minDistI] and shell[minDistI+1]
                shell[pointI] = shellI;
                minLevel[pointI] = levels[minDistI+1][cmpt];
            }
        }
    }
    else
    {
        // Inside/outside mode

        // Collect all those points that have a current maxLevel less than the
        // shell.

        pointField candidates(pt.size());
        labelList candidateMap(pt.size());
        label candidateI = 0;

        forAll(shell, pointI)
        {
            if (shell[pointI] == -1 && levels[0][cmpt] <= minLevel[pointI])
            {
                candidates[candidateI] = pt[pointI];
                candidateMap[candidateI] = pointI;
                candidateI++;
            }
        }
        candidates.setSize(candidateI);
        candidateMap.setSize(candidateI);

        // Do the expensive nearest test only for the candidate points.
        List<volumeType> volType;
        allGeometry_[shells_[shellI]].getVolumeType(candidates, volType);

        forAll(volType, i)
        {
            if
            (
                (
                    modes_[shellI] == INSIDE
                 && volType[i] == volumeType::INSIDE
                )
             || (
                    modes_[shellI] == OUTSIDE
                 && volType[i] == volumeType::OUTSIDE
                )
            )
            {
                label pointI = candidateMap[i];
                shell[pointI] = shellI;
                minLevel[pointI] = levels[0][cmpt];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::shellSurfaces::shellSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& shellsDict,
    const meshControl& controller,
    bool globalAddedCheck
)
:
    allGeometry_(allGeometry)
{
    // Wilcard specification : loop over all surfaces and try to find a match.

    // Count number of shells.
    label shellI = 0;
    forAll(allGeometry.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        if (shellsDict.found(geomName))
        {
            shellI++;
        }
    }

    // Size lists
    shells_.setSize(shellI);
    modes_.setSize(shellI);
    types_.setSize(shellI);
    distances_.setSize(shellI);
    levels_.setSize(shellI);
    refineSurface_.setSize(shellI);
    additionalInsideCheck_.setSize(shellI);
    coord_.setSize(shellI);

    //Global refine surface boundaries
    extendedGapLevel_.setSize(shellI);
    extendedGapMode_.setSize(shellI);

    FixedList<label, 3> nullGapLevel;
    nullGapLevel[0] = 0;
    nullGapLevel[1] = 0;
    nullGapLevel[2] = 0;


    HashSet<word> unmatchedKeys(shellsDict.toc());
    shellI = 0;

    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        const entry* ePtr = shellsDict.lookupEntryPtr(geomName, false, true);

        if (ePtr)
        {
            const dictionary& dict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            shells_[shellI] = geomI;
            modes_[shellI] = refineModeNames_.read(dict.lookup("mode"));
            refineSurface_[shellI] = true;

            additionalInsideCheck_[shellI] = dict.lookupOrDefault
            (
                "additionalInsideCheck",globalAddedCheck
            );

            coord_[shellI] = coordinateSystem::New
            (
                "cartesian",
                getCoordinateDict(dict)
            );

            if (modes_[shellI] == DIRECTION)
            {
                types_[shellI] = ISO;
                const List<Tuple2<Tuple2<point,point>,label>> distLevels =
                    dict.lookup("levels");
                setAndCheckLevels(shellI, distLevels);
            }
            else
            {
                // Read pairs of distance+level
                if (dict.found("levels"))
                {
                    types_[shellI] = ISO;
                    const List<Tuple2<scalar,label>> distLevels =
                        dict.lookup("levels");
                    setAndCheckLevels(shellI, distLevels);
                }
                else if (dict.found("alevels"))
                {
                    if (controller.algorithm() == meshControl::STANDARD)
                    {
                        types_[shellI] = ANISO;
                        const List<Tuple2<scalar,labelVector>> distLevels =
                            dict.lookup("alevels");
                        setAndCheckLevels(shellI, distLevels);
                        refineSurface_[shellI] =
                           dict.lookupOrDefault("refineSurface",false);
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "Anisotropic shell refinement keyword alevels"
                            << " is not valid on region " << geomName << nl
                            << " when running extrude/dual method."
                            << " Replace all instances of alevels with levels"
                            << exit(FatalError);
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << " No level or alevel keyword for refinement region "
                        << geomName
                        << exit(FatalError);
                }
            }

            // Gap specification
            // ~~~~~~~~~~~~~~~~~


            // Shell-wide gap level specification
            const searchableSurface& surface = allGeometry_[geomI];
            const wordList& regionNames = surface.regions();

            FixedList<label, 3> gapSpec
            (
                dict.lookupOrDefault
                (
                    "gapLevel",
                    nullGapLevel
                )
            );
            extendedGapLevel_[shellI].setSize(regionNames.size());
            extendedGapLevel_[shellI] = gapSpec;

            volumeType gapModeSpec
            (
                volumeType::names
                [
                    dict.lookupOrDefault<word>
                    (
                        "gapMode",
                        volumeType::names[volumeType::MIXED]
                    )
                ]
            );
            extendedGapMode_[shellI].setSize(regionNames.size());
            extendedGapMode_[shellI] = gapModeSpec;

            // Override on a per-region basis?

            if (dict.found("regions"))
            {
                const dictionary& regionsDict = dict.subDict("regions");
                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );
                        FixedList<label, 3> gapSpec
                        (
                            regionDict.lookupOrDefault
                            (
                                "gapLevel",
                                nullGapLevel
                            )
                        );
                        extendedGapLevel_[shellI][regionI] = gapSpec;

                        volumeType gapModeSpec
                        (
                            volumeType::names
                            [
                                regionDict.lookupOrDefault<word>
                                (
                                    "gapMode",
                                    volumeType::names[volumeType::MIXED]
                                )
                            ]
                        );
                        extendedGapMode_[shellI][regionI] = gapModeSpec;
                    }
                }
            }


            checkGapLevels(dict, shellI, extendedGapLevel_[shellI]);

            shellI++;
        }
    }

    if (unmatchedKeys.size() > 0)
    {
        IOWarningInFunction
        (
            shellsDict
        )   << "Not all entries in refinementRegions dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }

    // Orient shell surfaces before any searching is done. Note that this
    // only needs to be done for inside or outside. Orienting surfaces
    // constructs lots of addressing which we want to avoid.
    orient();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Highest shell level
Foam::label Foam::shellSurfaces::maxInsideOutsideLevel() const
{
    label overallMax = 0;
    forAll(levels_, shellI)
    {
        if (modes_[shellI] != DISTANCE && modes_[shellI] != DIRECTION)
        {
            overallMax = max(overallMax, cmptMax(levels_[shellI][0]));
        }
    }
    return overallMax;
}

// Highest shell level
Foam::label Foam::shellSurfaces::maxIsoLevel() const
{
    label overallMax = 0;
    forAll(levels_, shellI)
    {
        if (types_[shellI] == ISO)
        {
            const List<labelVector>& levels = levels_[shellI];
            forAll(levels, levelI)
            {
                overallMax = max(overallMax, levels[levelI][0]);
            }
        }
    }
    return overallMax;
}

// Highest aniso-shell level
Foam::label Foam::shellSurfaces::maxAnisoLevel(direction cmpt) const
{
    label overallMax = 0;
    forAll(levels_, shellI)
    {
        if (types_[shellI] == ANISO)
        {
            const List<labelVector>& levels = levels_[shellI];
            forAll(levels, levelI)
            {
                overallMax = max(overallMax, levels[levelI][cmpt]);
            }
        }
    }
    return overallMax;
}


// Check whether anisotropic refinement
bool Foam::shellSurfaces::hasAnisotropicRefinement() const
{
    forAll(types_, shellI)
    {
        if (types_[shellI] == ANISO)
        {
            return true;
        }
    }

    return false;
}

Foam::labelList Foam::shellSurfaces::maxGapLevel() const
{
    labelList surfaceMax(extendedGapLevel_.size(), 0);

    forAll(extendedGapLevel_, shelli)
    {
        const List<FixedList<label, 3>>& levels = extendedGapLevel_[shelli];
        forAll(levels, i)
        {
            surfaceMax[shelli] = max(surfaceMax[shelli], levels[i][2]);
        }
    }
    return surfaceMax;
}


void Foam::shellSurfaces::findHigherLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& maxLevel,
    const bool threaded /*= false*/
) const
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;

    forAll(shells_, shelli)
    {
        findHigherLevel(pt, shelli, maxLevel, threaded);
    }
}


void Foam::shellSurfaces::findHigherLevel
(
    const fvMesh& mesh,
    const pointField& pt,
    const labelList& ptLevel,
    const labelList& ptMap,
    labelList& maxLevel,
    labelList& minCmptLevel,
    labelList& shellCheck,
    const refineType refType,
    const direction cmpt
) const
{
    // Maximum level of any shell. Start off with level of point.
    maxLevel = ptLevel;
    minCmptLevel.setSize(maxLevel.size(), -1);
    shellCheck.setSize(maxLevel.size(), -1);

    // Swap neighbouring cell centres
    pointField neiCc(mesh.nFaces()-mesh.nInternalFaces());
    const pointField& cellCentres = mesh.cellCentres();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        const labelUList& faceCells = pp.faceCells();
        const vectorField::subField faceCentres = pp.faceCentres();

        label bFaceI = pp.start()-mesh.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiCc[bFaceI] = cellCentres[faceCells[i]];
                bFaceI++;
            }
        }
        else
        {
            forAll(faceCells, i)
            {
                neiCc[bFaceI] = faceCentres[i];
                bFaceI++;
            }
        }
    }
    // Swap coupled boundaries. Apply separation to cc since is coordinate.
    syncTools::swapBoundaryFacePositions(mesh, neiCc);

    pointField start(mesh.nFaces());
    pointField end(mesh.nFaces());
    forAll(mesh.faces(), faceI)
    {
        label own = mesh.faceOwner()[faceI];

        if (mesh.isInternalFace(faceI))
        {
            label nei = mesh.faceNeighbour()[faceI];

            start[faceI] = cellCentres[own];
            end[faceI] = cellCentres[nei];
        }
        else
        {
            label bFaceI = faceI - mesh.nInternalFaces();

            start[faceI] = cellCentres[own];
            end[faceI] = neiCc[bFaceI];
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(ROOTSMALL*(end-start));
        start -= smallVec;
        end += smallVec;
    }

    forAll(shells_, shellI)
    {
        if (types_[shellI] != refType)
        {
            continue;
        }

        labelList testRegions(pt.size(), -1);

            //const searchableSurface& s = allGeometry_[shells_[shellI]];
            //if (isA<triSurfaceMesh>(s) && modes_[shellI] != DISTANCE)
        boolList blockedFaces(mesh.nFaces(), false);
        if (modes_[shellI] != DISTANCE)
        {
            List<pointIndexHit> hitInfo(start.size());
            allGeometry_[shells_[shellI]].findLineAny(start, end, hitInfo);
            // Determine global regions, separated by blockedFaces

            forAll(mesh.faces(), faceI)
            {
                if (hitInfo[faceI].hit())
                {
                    blockedFaces[faceI] = true;
                }
            }
            syncTools::syncFaceList(mesh, blockedFaces, orEqOp<bool>());

            regionSplit globalRegion(mesh, blockedFaces);
            forAll(pt, i)
            {
                testRegions[i] = globalRegion[ptMap[i]];
            }
        }

        findHigherLevel
        (
            pt,
            shellI,
            testRegions,
            maxLevel,
            minCmptLevel,
            shellCheck,
            refType,
            cmpt
        );

        if (additionalInsideCheck()[shellI] && modes_[shellI] != DISTANCE)
        {
            DynamicList<label> interfaceCells(mesh.nCells()/100);
            label nChecks = 0;
            forAll(pt, i)
            {
                if (shellCheck[i] != shellI)
                {
                    label cellI = ptMap[i];
                    cell c = mesh.cells()[cellI];
                    bool checkCell = false;
                    forAll(c, cFI)
                    {
                        label faceI = c[cFI];
                        if (blockedFaces[faceI])
                        {
                            checkCell = true;
                            break;
                        }
                    }

                    if (checkCell)
                    {
                        interfaceCells.append(cellI);
                        if (refType == ANISO)
                        {
                            nChecks += 2;
                        }
                        else
                        {
                            nChecks += 8;
                        }
                    }
                }
            }

            labelList cellID(nChecks);
            pointField start(nChecks);
            pointField end(nChecks);

            nChecks = 0;
            forAll(interfaceCells, i)
            {
                label cellI = interfaceCells[i];
                const labelList& cPoints = mesh.cellPoints()[cellI];
                treeBoundBox bbox(mesh.points(),cPoints);

                point cc = mesh.cellCentres()[cellI];

                if (refType == ANISO)
                {
                    vector dir = vector::zero;
                    dir[cmpt] = 1;
                    vector dC = 0.25*(bbox.max()[cmpt] - bbox.min()[cmpt])*dir;

                    start[nChecks] = cc;
                    end[nChecks] = cc+dC;
                    cellID[nChecks] = cellI;
                    nChecks++;

                    start[nChecks] = cc;
                    end[nChecks] = cc-dC;
                    cellID[nChecks] = cellI;
                    nChecks++;
                }
                else
                {
                    for (direction i = 0; i < 8; i++)
                    {
                        start[nChecks] = cc;
                        end[nChecks] = 0.5*(cc+bbox.corner(i));
                        cellID[nChecks] = cellI;
                        nChecks++;
                    }
                }
            }

            List<List<pointIndexHit>> hitInfo;
            allGeometry_[shells_[shellI]].findLineAll(start, end, hitInfo);

            forAll(start, startI)
            {
                if (hitInfo[startI].size() == 1)
                {
                    label i = cellID[startI];
                    maxLevel[i] = levels_[shellI][0][cmpt];
                    minCmptLevel[i] = cmptMin(levels_[shellI][0]);
                    shellCheck[i] = shellI;
                }
            }
        }
    }
}

void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& gapShell,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    gapShell.setSize(pt.size());
    gapShell = -1;

    FixedList<label, 3> nullGapLevel;
    nullGapLevel[0] = 0;
    nullGapLevel[1] = 0;
    nullGapLevel[2] = 0;

    gapInfo.setSize(pt.size());
    gapInfo = nullGapLevel;

    gapMode.setSize(pt.size());
    gapMode = volumeType::MIXED;

    forAll(shells_, shelli)
    {
        findHigherGapLevel(pt, ptLevel, shelli, gapShell, gapInfo, gapMode);
    }
}


void Foam::shellSurfaces::findHigherGapLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    List<FixedList<label, 3>>& gapInfo,
    List<volumeType>& gapMode
) const
{
    labelList gapShell;
    findHigherGapLevel(pt, ptLevel, gapShell, gapInfo, gapMode);
}


void Foam::shellSurfaces::findLevel
(
    const pointField& pt,
    const labelList& ptLevel,
    labelList& shell
) const
{
    shell.setSize(pt.size());
    shell = -1;

    labelList minLevel(ptLevel);

    forAll(shells_, shelli)
    {
        findLevel(pt, shelli, minLevel, shell);
    }
}


// ************************************************************************* //
