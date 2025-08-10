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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "refinementSurfaces/refinementSurfaces.H"
#include "triSurface/orientedSurface/orientedSurface.H"
#include "db/Time/Time.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "shellSurfaces/shellSurfaces.H"
#include "primitives/Pair/labelPair.H"
#include "searchableSurfaces/searchableSurfacesQueries/searchableSurfacesQueries.H"
#include "containers/Lists/ListOps/ListOps.H"
#include "containers/Lists/UPtrList/UPtrList.H"
#include "global/unitConversion/unitConversion.H"
#include "algorithms/indexedOctree/volumeType.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "distributedTriSurfaceMesh/distributedTriSurfaceMesh.H"
#include "searchableSurfaces/searchableSphere/searchableSphere.H"
#include "searchableSurfaces/searchableBox/searchableBox.H"
#include "searchableSurfaces/searchableCylinder/searchableCylinder.H"
#include "searchableSurfaces/searchableRing/searchableRing.H"
#include "searchableSurfaces/searchableCone/searchableCone.H"
#include "meshes/meshTools/simpleVTKWriter.H"

#ifdef FOAM_USE_TBB
  //#include <tbb/spin_mutex.h>
  #include <tbb/enumerable_thread_specific.h>
  #include <tbb/parallel_for.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(refinementSurfaces, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::refinementSurfaces::findHigherLevel
(
    const searchableSurface& geom,
    const shellSurfaces& shells,
    const List<pointIndexHit>& intersectionInfo,
    const labelList& surfaceLevel,      // current level
    const bool threaded /*= false*/
) const
{
    // See if a cached level field available
    labelList minLevelField;
    geom.getField(intersectionInfo, minLevelField, threaded);


    // Detect any uncached values and do proper search
    labelList localLevel(surfaceLevel);
    {
        // Check hits:
        // 1. cached value == -1 : store for re-testing
        // 2. cached value != -1 : use
        // 3. uncached : use region 0 value

        DynamicList<label> retestSet;
        label nHits = 0;

        //TODO threaded parallel_reduce
        forAll(intersectionInfo, i)
        {
            if (intersectionInfo[i].hit())
            {
                nHits++;

                // Check if minLevelField for this surface.
                if (minLevelField.size())
                {
                    if (minLevelField[i] == -1)
                    {
                        retestSet.append(i);
                    }
                    else
                    {
                        localLevel[i] = max(localLevel[i], minLevelField[i]);
                    }
                }
                else
                {
                    retestSet.append(i);
                }
            }
        }

        label nRetest = returnReduce(retestSet.size(), sumOp<label>());
        if (nRetest > 0)
        {
            reduce(nHits, sumOp<label>());

            //Info<< " Retesting " << nRetest
            //    << " out of " << nHits
            //    << " intersections on uncached elements on geometry "
            //    << geom.name() << endl;

            pointField samples(retestSet.size());
            if (threaded && retestSet.size() > 2048)
            {
#ifndef FOAM_USE_TBB
                WarningInFunction
                    << "FOAM not linked against TBB! Cannot run multithreaded."
                    << endl;
#else
                const auto grainSize =
                    std::max<size_t>
                    (
                        retestSet.size() / tbb::this_task_arena::max_concurrency(),
                        1024
                    );

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, retestSet.size(), grainSize),
                    [&](const tbb::blocked_range<size_t>& r)
                    {
                        for (size_t i = r.begin(); i < r.end(); ++i)
                        {
                            samples[i] = intersectionInfo[retestSet[i]].hitPoint();
                        }
                    },
                    tbb::simple_partitioner()
                );
#endif
            }
            else
            {
                forAll(retestSet, i)
                {
                    samples[i] = intersectionInfo[retestSet[i]].hitPoint();
                }
            }

            labelList shellLevel;
            //TODO threaded (parallel_reduce in shells.findHigherLevel)
            shells.findHigherLevel
            (
                samples,
                UIndirectList<label>(surfaceLevel, retestSet)(),
                shellLevel,
                threaded
            );

            if (threaded && retestSet.size() > 2048)
            {
#ifndef FOAM_USE_TBB
                WarningInFunction
                    << "FOAM not linked against TBB! Cannot run multithreaded."
                    << endl;
#else
                const auto grainSize =
                    std::max<size_t>
                    (
                        retestSet.size() / tbb::this_task_arena::max_concurrency(),
                        1024
                    );

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, retestSet.size(), grainSize),
                    [&](const tbb::blocked_range<size_t>& r)
                    {
                        for (size_t i = r.begin(); i < r.end(); ++i)
                        {
                            label sampleI = retestSet[i];
                            localLevel[sampleI] = max(localLevel[sampleI], shellLevel[i]);
                        }
                    },
                    tbb::simple_partitioner()
                );
#endif
            }
            else
            {
                forAll(retestSet, i)
                {
                    label sampleI = retestSet[i];
                    localLevel[sampleI] = max(localLevel[sampleI], shellLevel[i]);
                }
            }
        } // if (nRetest > 0)
    }

    return localLevel;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::refinementSurfaces::refinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& surfacesDict,
    const dictionary& refineDict,
    const label gapLevelIncrement
)
:
    allGeometry_(allGeometry),
    surfaces_(surfacesDict.size()),
    names_(surfacesDict.size()),
    surfZones_(surfacesDict.size()),
    regionOffset_(surfacesDict.size()),
    moveCentroids_(surfacesDict.size(), scalar(0)),
    wrapLevel_(surfacesDict.size(), label(-1)),
    boundingSurface_(surfacesDict.size(), false),
    refineBoundary_(surfacesDict.size(), false),
    singleRegion_(surfacesDict.size(), false),
    nRegionsWrapRemoval_
    (
        refineDict.lookupOrDefault<label>("nRegionsWrapRemove",0)
    ),
    reorient_(refineDict.lookupOrDefault<Switch>("reorient",false))
{
    // Wildcard specification : loop over all surface, all regions
    // and try to find a match.

    // Count number of surfaces.
    label surfI = 0;
    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        if (surfacesDict.found(geomName))
        {
            surfI++;
        }
    }

    // Size lists
    surfaces_.setSize(surfI);
    names_.setSize(surfI);
    surfZones_.setSize(surfI);
    regionOffset_.setSize(surfI);
    moveCentroids_.setSize(surfI,scalar(0));
    refineBoundary_.setSize(surfI,false);
    singleRegion_.setSize(surfI,false);
    wrapLevel_.setSize(surfI,label(-1));
    boundingSurface_.setSize(surfI,false);

    labelList globalMinLevel(surfI, 0);
    labelList globalMaxLevel(surfI, 0);
    labelList globalLevelIncr(surfI, 0);

    FixedList<label, 3> nullGapLevel;
    nullGapLevel[0] = 0;
    nullGapLevel[1] = 0;
    nullGapLevel[2] = 0;

    List<FixedList<label, 3>> globalGapLevel(surfI);
    List<volumeType> globalGapMode(surfI);

    scalarField globalAngle(surfI, -GREAT);

    boolList globalRefineOnly
    (
        surfI,
        refineDict.lookupOrDefault<Switch>("refineFeatureEdgesOnly", false)
    );

    boolList globalSingleCellClosure
    (
        surfI,
        refineDict.lookupOrDefault<Switch>("singleCellGapClosure", false)
    );

    boolList globalThinGap
    (
        surfI,
        refineDict.lookupOrDefault<Switch>("thinGap", false)
    );

    boolList globalCornerCells
    (
        surfI,
        refineDict.lookupOrDefault<Switch>("keepCornerCells", false)
    );

    scalar refineAngle = 0.866;
    if (refineDict.found("featureRefineAngle"))
    {
        scalar featRefineAngle =
            readScalar(refineDict.lookup("featureRefineAngle"));

        if (featRefineAngle < 0 || featRefineAngle > 180)
        {
            refineAngle = -GREAT;
        }
        else
        {
            refineAngle  =  Foam::cos(degToRad(featRefineAngle));
        }
    }

    scalarField globalFeatureAngle
    (
        surfI,
        refineAngle
    );

    scalar baffleAngle = GREAT;
    if (refineDict.found("minBaffleAngle"))
    {
        scalar angle =
            readScalar(refineDict.lookup("minBaffleAngle"));

        if (angle < 0 || angle > 180)
        {
            baffleAngle = GREAT;
        }
        else
        {
            baffleAngle  =  Foam::cos(degToRad(angle));
        }
    }

    scalarField globalBaffleAngle
    (
        surfI,
        baffleAngle
    );

    labelList globalProxLevelIncr
    (
        surfI,
        refineDict.lookupOrDefault<label>("proximityIncrement",-1)
    );

    labelList globalProxDir
    (
        surfI,
        refineDict.lookupOrDefault<label>("proxDir",0)
    );

    scalarField globalCurvature
    (
        surfI,
        degToRad(refineDict.lookupOrDefault<scalar>("curvature",-1))
    );

    labelList globalProxMaxCells
    (
        surfI,
        refineDict.lookupOrDefault<label>("maxCellsAcrossGap",2)
    );
    PtrList<dictionary> globalPatchInfo(surfI);

    List<Map<label>> regionMinLevel(surfI);
    List<Map<label>> regionMaxLevel(surfI);
    List<Map<label>> regionLevelIncr(surfI);
    List<Map<FixedList<label, 3>>> regionGapLevel(surfI);
    List<Map<volumeType>> regionGapMode(surfI);
    List<Map<scalar>> regionAngle(surfI);
    List<Map<scalar>> regionFeatureAngle(surfI);
    List<Map<bool>> regionRefineOnly(surfI);
    List<Map<bool>> regionSingleCellClosure(surfI);
    List<Map<bool>> regionThinGap(surfI);
    List<Map<bool>> regionCornerCells(surfI);
    List<Map<scalar>> regionBaffleAngle(surfI);
    List<Map<autoPtr<dictionary>>> regionPatchInfo(surfI);
    List<Map<label>> regionProxLevelIncr(surfI);
    List<Map<label>> regionProxDir(surfI);
    List<Map<scalar>> regionCurvature(surfI);
    List<Map<label>> regionProxMaxCells(surfI);

    //Global allow free standing zone faces
    bool globalFreeStanding =
        refineDict.lookupOrDefault<bool>("allowFreeStandingZoneFaces",false);

    //Global zone by walking method
    bool globalZoneWalk =
        refineDict.lookupOrDefault<bool>("zoneByWalk",false);

    //Global movement of centroids
    scalar globalMoveCentroids =
        refineDict.lookupOrDefault<scalar>("moveCentroidsTol",scalar(0));

    //Global whether to close surface at specified level
    label globalWrapLevel =
        refineDict.lookupOrDefault<label>("wrapLevel",label(-1));

    //Global whether to use surface in autoBlock calculation
    label globalBoundingSurface =
        refineDict.lookupOrDefault<bool>("autoBlockSurface",false);

    //Global refine surface boundaries
    bool globalRefineBoundary =
        refineDict.lookupOrDefault<bool>("refineSurfaceBoundary",false);

    HashSet<word> unmatchedKeys(surfacesDict.toc());

    surfI = 0;
    forAll(allGeometry_.names(), geomI)
    {
        const word& geomName = allGeometry_.names()[geomI];

        const entry* ePtr = surfacesDict.lookupEntryPtr(geomName, false, true);
        if (ePtr)
        {
            const dictionary& dict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            names_[surfI] = geomName;
            surfaces_[surfI] = geomI;

            const labelPair refLevel(dict.lookup("level"));
            globalMinLevel[surfI] = refLevel[0];
            globalMaxLevel[surfI] = refLevel[1];
            globalLevelIncr[surfI] = dict.lookupOrDefault
            (
                "gapLevelIncrement",
                gapLevelIncrement
            );

            if
            (
                globalMinLevel[surfI] < 0
             || globalMaxLevel[surfI] < globalMinLevel[surfI]
             || globalLevelIncr[surfI] < 0
            )
            {
                FatalIOErrorInFunction(dict)
                    << "Illegal level specification for surface "
                    << names_[surfI]
                    << " : minLevel:" << globalMinLevel[surfI]
                    << " maxLevel:" << globalMaxLevel[surfI]
                    << " levelIncrement:" << globalLevelIncr[surfI]
                    << exit(FatalIOError);
            }


            // Optional gapLevel specification

            globalGapLevel[surfI] = dict.lookupOrDefault
            (
                "gapLevel",
                nullGapLevel
            );
            globalGapMode[surfI] = volumeType::names
            [
                dict.lookupOrDefault<word>
                (
                    "gapMode",
                    volumeType::names[volumeType::MIXED]
                )
            ];
            if
            (
                globalGapMode[surfI] == volumeType::UNKNOWN
             || globalGapLevel[surfI][0] < 0
             || globalGapLevel[surfI][1] < 0
             || globalGapLevel[surfI][2] < 0
             || globalGapLevel[surfI][1] > globalGapLevel[surfI][2]
            )
            {
                FatalIOErrorInFunction(dict)
                    << "Illegal gapLevel specification for surface "
                    << names_[surfI]
                    << " : gapLevel:" << globalGapLevel[surfI]
                    << " gapMode:" << volumeType::names[globalGapMode[surfI]]
                    << exit(FatalIOError);
            }


            const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

            // Surface zones
            surfZones_.set
            (
                surfI,
                new surfaceZonesInfo
                (
                    surface,
                    dict,
                    globalFreeStanding,
                    globalZoneWalk
                )
            );

            // Global perpendicular angle
            if (dict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfI,
                    dict.subDict("patchInfo").clone()
                );
            }

            dict.readIfPresent("perpendicularAngle", globalAngle[surfI]);

            moveCentroids_[surfI] = dict.lookupOrDefault<scalar>
                ("moveCentroidsTol",globalMoveCentroids);

            wrapLevel_[surfI] =
                dict.lookupOrDefault<label>("wrapLevel",globalWrapLevel);

            boundingSurface_[surfI] =
                dict.lookupOrDefault<bool>
                (
                    "autoBlockSurface",
                    globalBoundingSurface
                );

            refineBoundary_[surfI] = dict.lookupOrDefault<bool>
                ("refineSurfaceBoundary",globalRefineBoundary);

            if (dict.found("featureRefineAngle"))
            {
                scalar featRefineAngle =
                    readScalar(dict.lookup("featureRefineAngle"));
                if (featRefineAngle < 0 || featRefineAngle > 180)
                {
                    globalFeatureAngle[surfI] = -GREAT;
                }
                else
                {
                    globalFeatureAngle[surfI]  =
                        Foam::cos(degToRad(featRefineAngle));
                }
            }

            if (dict.found("refineFeatureEdgesOnly"))
            {
                bool refineOnly =
                    dict.lookupOrDefault<Switch>
                    ("refineFeatureEdgesOnly", false);

                globalRefineOnly[surfI]  = refineOnly;
            }

            if (dict.found("singleCellGapClosure"))
            {
                bool gapClosure =
                    dict.lookupOrDefault<Switch>
                    ("singleCellGapClosure", false);

                globalSingleCellClosure[surfI]  = gapClosure;
            }

            if (dict.found("thinGap"))
            {
                bool thinGapMethod =
                    dict.lookupOrDefault<Switch>
                    ("thinGap", false);

                globalThinGap[surfI]  = thinGapMethod;
            }

            if (dict.found("keepCornerCells"))
            {
                bool cornerCellsMethod =
                    dict.lookupOrDefault<Switch>
                    ("keepCornerCells", false);

                globalCornerCells[surfI]  = cornerCellsMethod;
            }

            if (dict.found("minBaffleAngle"))
            {
                scalar angle =
                    readScalar(dict.lookup("minBaffleAngle"));
                if (angle < 0 || angle > 180)
                {
                    globalBaffleAngle[surfI] = GREAT;
                }
                else
                {
                    globalBaffleAngle[surfI]  = Foam::cos(degToRad(angle));
                }
            }

            if (dict.found("proximityIncrement"))
            {
                globalProxLevelIncr[surfI] = readLabel
                (
                    dict.lookup("proximityIncrement")
                );
            }

            if (dict.found("proxDir"))
            {
                globalProxDir[surfI] = readLabel
                (
                    dict.lookup("proxDir")
                );
            }

            if (dict.found("curvature"))
            {
                globalCurvature[surfI] = degToRad
                (
                    readScalar
                    (
                        dict.lookup("curvature")
                    )
                );
            }

            if (dict.found("maxCellsAcrossGap"))
            {
                globalProxMaxCells[surfI] = readLabel
                (
                    dict.lookup("maxCellsAcrossGap")
                );
            }

            singleRegion_[surfI] = geometry().singleRegion()[geomI];

            if (dict.found("regions") && !singleRegion_[surfI])
            {
                const dictionary& regionsDict = dict.subDict("regions");
                const wordList& regionNames = surface.regions();

                forAll(regionNames, regionI)
                {
                    if (regionsDict.found(regionNames[regionI]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regionI]
                        );

                        const labelPair refLevel(regionDict.lookup("level"));

                        regionMinLevel[surfI].insert(regionI, refLevel[0]);
                        regionMaxLevel[surfI].insert(regionI, refLevel[1]);
                        label levelIncr = regionDict.lookupOrDefault
                        (
                            "gapLevelIncrement",
                            gapLevelIncrement
                        );
                        regionLevelIncr[surfI].insert(regionI, levelIncr);



                        if
                        (
                            refLevel[0] < 0
                         || refLevel[1] < refLevel[0]
                         || levelIncr < 0
                        )
                        {
                            FatalIOErrorInFunction(dict)
                                << "Illegal level specification for surface "
                                << names_[surfI] << " region "
                                << regionNames[regionI]
                                << " : minLevel:" << refLevel[0]
                                << " maxLevel:" << refLevel[1]
                                << " levelIncrement:" << levelIncr
                                << exit(FatalIOError);
                        }

                        if (regionDict.found("featureRefineAngle"))
                        {
                            scalar featRefineAngle = readScalar
                                (regionDict.lookup("featureRefineAngle"));

                            if (featRefineAngle < 0 || featRefineAngle > 180)
                            {
                                featRefineAngle = -GREAT;
                            }
                            else
                            {
                                featRefineAngle  =
                                    Foam::cos(degToRad(featRefineAngle));
                            }
                            regionFeatureAngle[surfI].insert
                            (
                                regionI,
                                featRefineAngle
                            );
                        }

                        if (regionDict.found("refineFeatureEdgesOnly"))
                        {
                            bool refineOnly =
                                regionDict.lookupOrDefault<Switch>
                                ("refineFeatureEdgesOnly", false);

                            regionRefineOnly[surfI].insert
                            (
                                regionI,
                                refineOnly
                            );
                        }

                        if (regionDict.found("singleCellGapClosure"))
                        {
                            bool gapClosure =
                                regionDict.lookupOrDefault<Switch>
                                ("singleCellGapClosure", false);

                            regionSingleCellClosure[surfI].insert
                            (
                                regionI,
                                gapClosure
                            );
                        }

                        if (regionDict.found("thinGap"))
                        {
                            bool thinGapMethod =
                                regionDict.lookupOrDefault<Switch>
                                ("thinGap", false);

                            regionThinGap[surfI].insert
                            (
                                regionI,
                                thinGapMethod
                            );
                        }

                        if (regionDict.found("keepCornerCells"))
                        {
                            bool cornerCellsMethod =
                                regionDict.lookupOrDefault<Switch>
                                ("keepCornerCells", false);

                            regionCornerCells[surfI].insert
                            (
                                regionI,
                                cornerCellsMethod
                            );
                        }

                        if (regionDict.found("minBaffleAngle"))
                        {
                            scalar angle =
                                readScalar(dict.lookup("minBaffleAngle"));
                            if (angle < 0 || angle > 180)
                            {
                                angle = GREAT;
                            }
                            else
                            {
                                angle  = Foam::cos(degToRad(angle));
                            }
                            regionBaffleAngle[surfI].insert
                            (
                                regionI,
                                angle
                            );
                        }

                        // Optional gapLevel specification
                        FixedList<label, 3> gapSpec
                        (
                            regionDict.lookupOrDefault
                            (
                                "gapLevel",
                                nullGapLevel
                            )
                        );
                        regionGapLevel[surfI].insert(regionI, gapSpec);
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
                        regionGapMode[surfI].insert(regionI, gapModeSpec);
                        if
                        (
                            gapModeSpec == volumeType::UNKNOWN
                         || gapSpec[0] < 0
                         || gapSpec[1] < 0
                         || gapSpec[2] < 0
                         || gapSpec[1] > gapSpec[2]
                        )
                        {
                            FatalIOErrorInFunction(dict)
                                << "Illegal gapLevel specification for surface "
                                << names_[surfI]
                                << " : gapLevel:" << gapSpec
                                << " gapMode:" << volumeType::names[gapModeSpec]
                                << exit(FatalIOError);
                        }


                        if (regionDict.found("perpendicularAngle"))
                        {
                            regionAngle[surfI].insert
                            (
                                regionI,
                                readScalar
                                (
                                    regionDict.lookup("perpendicularAngle")
                                 )
                             );
                        }

                        if (regionDict.found("proximityIncrement"))
                        {
                            label levelIncr = readLabel
                            (
                                regionDict.lookup("proximityIncrement")
                            );
                            regionProxLevelIncr[surfI].insert
                                (regionI, levelIncr);
                        }

                        if (regionDict.found("proxDir"))
                        {
                            label proxDir = readLabel
                            (
                                regionDict.lookup("proxDir")
                            );
                            regionProxDir[surfI].insert
                                (regionI, proxDir);
                        }

                        if (regionDict.found("curvature"))
                        {
                            scalar curvature = degToRad
                            (
                                readScalar(regionDict.lookup("curvature"))
                            );
                            regionCurvature[surfI].insert
                                (regionI, curvature);
                        }

                        if (regionDict.found("maxCellsAcrossGap"))
                        {
                            label levelMaxCells = readLabel
                            (
                                regionDict.lookup("maxCellsAcrossGap")
                            );
                            regionProxMaxCells[surfI].insert
                                (regionI, levelMaxCells);
                        }

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfI].insert
                            (
                                regionI,
                                regionDict.subDict("patchInfo").clone()
                             );
                        }
                    }
                }
            }
            surfI++;
        }
    }

    if (unmatchedKeys.size() > 0)
    {
        IOWarningInFunction
        (
            surfacesDict
        )   << "Not all entries in refinementSurfaces dictionary were used."
            << " The following entries were not used : "
            << unmatchedKeys.sortedToc()
            << endl;
    }


    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces_, surfI)
    {
        regionOffset_[surfI] = nRegions;

        if (singleRegion_[surfI])
        {
            nRegions++;
        }
        else
        {
            nRegions += allGeometry_[surfaces_[surfI]].regions().size();
        }
    }

    // Rework surface specific information into information per global region
    minLevel_.setSize(nRegions);
    minLevel_ = 0;
    maxLevel_.setSize(nRegions);
    maxLevel_ = 0;
    gapLevel_.setSize(nRegions);
    gapLevel_ = -1;
    featureRefineAngle_.setSize(nRegions);
    featureRefineAngle_ = -GREAT;
    featureRefineOnly_.setSize(nRegions);
    featureRefineOnly_ = false;
    singleCellClosure_.setSize(nRegions);
    singleCellClosure_ = false;
    thinGap_.setSize(nRegions);
    thinGap_ = false;
    cornerCellRegions_.setSize(nRegions);
    cornerCellRegions_ = false;
    extendedGapLevel_.setSize(nRegions);
    extendedGapLevel_ = nullGapLevel;
    extendedGapMode_.setSize(nRegions);
    extendedGapMode_ = volumeType::UNKNOWN;
    perpendicularAngle_.setSize(nRegions);
    perpendicularAngle_ = -GREAT;
    patchInfo_.setSize(nRegions);
    proxLevelIncr_.setSize(nRegions);
    proxLevelIncr_ = -1;
    proxDir_.setSize(nRegions);
    proxDir_ = 0;
    curvature_.setSize(nRegions);
    curvature_ = -1;
    proxMaxCells_.setSize(nRegions);
    proxMaxCells_ = 0;
    minBaffleAngle_.setSize(nRegions);
    minBaffleAngle_ = GREAT;

    forAll(globalMinLevel, surfI)
    {
        label nRegions
        (
            singleRegion_[surfI] ?
            1 : allGeometry_[surfaces_[surfI]].regions().size()
        );

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegionI = regionOffset_[surfI] + i;
            minLevel_[globalRegionI] = globalMinLevel[surfI];
            maxLevel_[globalRegionI] = globalMaxLevel[surfI];
            gapLevel_[globalRegionI] = maxLevel_[globalRegionI]
              + globalLevelIncr[surfI];

            extendedGapLevel_[globalRegionI] = globalGapLevel[surfI];
            extendedGapMode_[globalRegionI] = globalGapMode[surfI];
            perpendicularAngle_[globalRegionI] = globalAngle[surfI];
            featureRefineAngle_[globalRegionI] = globalFeatureAngle[surfI];
            featureRefineOnly_[globalRegionI] = globalRefineOnly[surfI];
            singleCellClosure_[globalRegionI] = globalSingleCellClosure[surfI];
            thinGap_[globalRegionI] = globalThinGap[surfI];
            cornerCellRegions_[globalRegionI] = globalCornerCells[surfI];
            minBaffleAngle_[globalRegionI] = globalBaffleAngle[surfI];
            proxLevelIncr_[globalRegionI] = globalProxLevelIncr[surfI];
            proxDir_[globalRegionI] = globalProxDir[surfI];
            curvature_[globalRegionI] = globalCurvature[surfI];
            proxMaxCells_[globalRegionI] = globalProxMaxCells[surfI];

            if (globalPatchInfo.set(surfI))
            {
                patchInfo_.set
                (
                    globalRegionI,
                    globalPatchInfo[surfI].clone()
                );
            }
        }

        // Overwrite with region specific information
        forAllConstIter(Map<label>, regionMinLevel[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            minLevel_[globalRegionI] = iter();
            maxLevel_[globalRegionI] = regionMaxLevel[surfI][iter.key()];
            gapLevel_[globalRegionI] = maxLevel_[globalRegionI]
                + regionLevelIncr[surfI][iter.key()];
            extendedGapLevel_[globalRegionI] =
                regionGapLevel[surfI][iter.key()];
            extendedGapMode_[globalRegionI] =
                regionGapMode[surfI][iter.key()];
        }

        forAllConstIter(Map<scalar>, regionFeatureAngle[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            featureRefineAngle_[globalRegionI] =
                regionFeatureAngle[surfI][iter.key()];
        }

        forAllConstIter(Map<bool>, regionRefineOnly[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            featureRefineOnly_[globalRegionI] =
                regionRefineOnly[surfI][iter.key()];
        }

        forAllConstIter(Map<bool>, regionSingleCellClosure[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            singleCellClosure_[globalRegionI] =
                regionSingleCellClosure[surfI][iter.key()];
        }

        forAllConstIter(Map<bool>, regionThinGap[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            thinGap_[globalRegionI] = regionThinGap[surfI][iter.key()];
        }

        forAllConstIter(Map<bool>, regionCornerCells[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            cornerCellRegions_[globalRegionI] =
                regionCornerCells[surfI][iter.key()];
        }

        forAllConstIter(Map<scalar>, regionBaffleAngle[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            minBaffleAngle_[globalRegionI] =
                regionBaffleAngle[surfI][iter.key()];
        }

        forAllConstIter(Map<scalar>, regionAngle[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            perpendicularAngle_[globalRegionI] = regionAngle[surfI][iter.key()];
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfI];
        forAllConstIter(Map<autoPtr<dictionary>>, localInfo, iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            patchInfo_.set(globalRegionI, iter()().clone());
        }

        forAllConstIter(Map<label>, regionProxLevelIncr[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            proxLevelIncr_[globalRegionI] =
                regionProxLevelIncr[surfI][iter.key()];
        }

        forAllConstIter(Map<label>, regionProxDir[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            proxDir_[globalRegionI] = regionProxDir[surfI][iter.key()];
        }

        forAllConstIter(Map<scalar>, regionCurvature[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            curvature_[globalRegionI] = regionCurvature[surfI][iter.key()];
        }

        forAllConstIter(Map<label>, regionProxMaxCells[surfI], iter)
        {
            label globalRegionI = regionOffset_[surfI] + iter.key();

            proxMaxCells_[globalRegionI] =
                regionProxMaxCells[surfI][iter.key()];
        }

    }

    // Orient named surfaces before any searching is done. Note that this
    // only needs to be done for inside or outside. Orienting surfaces
    // constructs lots of addressing which we want to avoid.
    orient();
}


// check whether contains bounded surface
bool Foam::refinementSurfaces::checkForFiniteSurfaces() const
{
    forAll(surfaces_, surfI)
    {
        const searchableSurface& s = allGeometry_[surfaces_[surfI]];
        if
        (
            isA<triSurfaceMesh>(s) || isA<searchableSphere>(s)
            || isA<searchableBox>(s) || isA<searchableCylinder>(s)
            || isA<searchableRing>(s) || isA<searchableCone>(s)
        )
        {
            return true;
        }
    }

    return false;
}


// Get bounding box of un-named finite surfaces
Foam::boundBox Foam::refinementSurfaces::getFiniteSurfacesBb
(
    const bool primitiveBounding
)
{
    boundBox overallBb
    (
        point(GREAT, GREAT, GREAT),
        point(-GREAT, -GREAT, -GREAT)
    );

    const labelList unNamedSurfaces =
    (
        surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces(surfZones())
    );

    forAll(unNamedSurfaces, i)
    {
        label surfI = unNamedSurfaces[i];

        const searchableSurface& s = allGeometry_[surfaces_[surfI]];
        if
        (
            isA<triSurfaceMesh>(s)
            ||
            (
               primitiveBounding && ( isA<searchableSphere>(s)
               || isA<searchableBox>(s) || isA<searchableCylinder>(s)
               || isA<searchableRing>(s) || isA<searchableCone>(s))
            )
        )
        {
            const boundBox& surfBb = s.bounds();
            overallBb.min() = min(overallBb.min(), surfBb.min());
            overallBb.max() = max(overallBb.max(), surfBb.max());
        }
    }

    if (Pstream::parRun())
    {
        List<scalar> bbox(3);
        bbox[0] = overallBb.min().x();
        bbox[1] = overallBb.min().y();
        bbox[2] = overallBb.min().z();

        Pstream::listCombineReduce(bbox,minOp<scalar>());

        overallBb.min().x() = bbox[0];
        overallBb.min().y() = bbox[1];
        overallBb.min().z() = bbox[2];

        bbox[0] = overallBb.max().x();
        bbox[1] = overallBb.max().y();
        bbox[2] = overallBb.max().z();

        Pstream::listCombineReduce(bbox, maxOp<scalar>());

        overallBb.max().x() = bbox[0];
        overallBb.max().y() = bbox[1];
        overallBb.max().z() = bbox[2];
    }

   return overallBb;
}


// Get bounding box in local coordinate direction
void Foam::refinementSurfaces::getFiniteSurfacesLocalBb
(
    const coordinateSystem& coord,
    const bool primitiveBounding,
    point& minPt,
    point& maxPt
) const
{
    const labelList selectedSurfaces = boundSurfaces();

    const labelList unNamedSurfaces =
    (
        surfaceZonesInfo::getUnnamedAndBoundaryNamedSurfaces(surfZones())
    );

    const labelList checkSurfaces
    (
        selectedSurfaces.size() == 0
        ?  unNamedSurfaces
        :  selectedSurfaces
    );

    forAll(checkSurfaces, i)
    {
        label surfI = checkSurfaces[i];

        const searchableSurface& s = allGeometry_[surfaces_[surfI]];
        if
        (
            isA<triSurfaceMesh>(s)
        )
        {
            const triSurfaceMesh& triMesh = refCast<const triSurfaceMesh>(s);
            triSurface& tri = const_cast<triSurface&>
            (
                refCast<const triSurface>(triMesh)
            );
            const pointField& pts = tri.points();

            forAll(pts, pti)
            {
                point localPt = coord.invTransformPoint(pts[pti]);
                minPt = ::Foam::min(minPt, localPt);
                maxPt = ::Foam::max(maxPt, localPt);
            }
        }
        else if
        (
            primitiveBounding && ( isA<searchableSphere>(s)
         || isA<searchableBox>(s) || isA<searchableCylinder>(s)
         || isA<searchableRing>(s) || isA<searchableCone>(s))
        )
        {
            const boundBox& surfBb = s.bounds();
            const tmp<pointField> tpts = surfBb.points();
            const pointField& pts = tpts();
            forAll(pts, pti)
            {
                point localPt = coord.invTransformPoint(pts[pti]);
                minPt = ::Foam::min(minPt, localPt);
                maxPt = ::Foam::max(maxPt, localPt);
            }
        }
    }

    if (Pstream::parRun())
    {
        reduce(
            std::tie(minPt, maxPt),
            ParallelOp<minMagSqrOp<vector>, maxMagSqrOp<vector>>{}
        );
    }

    return;
}


// Get bounding box of un-named tri-surfaces
Foam::boundBox Foam::refinementSurfaces::getUnnamedBb()
{
    boundBox overallBb
    (
        point(GREAT, GREAT, GREAT),
        point(-GREAT, -GREAT, -GREAT)
    );

    const labelList unNamedSurfaces =
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfZones())
    );

    forAll(unNamedSurfaces, i)
    {
        label surfI = unNamedSurfaces[i];
        const searchableSurface& s = allGeometry_[surfaces_[surfI]];
        if (isA<triSurfaceMesh>(s))
        {
            triSurfaceMesh& surf = const_cast<triSurfaceMesh&>
            (
                refCast<const triSurfaceMesh>(s)
            );
            tmp<pointField> tpoints(surf.points());
            const pointField& points = tpoints();
            boundBox surfBb(points, false);
            overallBb.min() = min(overallBb.min(), surfBb.min());
            overallBb.max() = max(overallBb.max(), surfBb.max());
        }
    }

   return overallBb;
}


// Specifically orient zoned triSurfaces using a calculated point outside.
// Done since quite often triSurfaces not of consistent orientation which
// is (currently) necessary for sideness calculation
void Foam::refinementSurfaces::orient()
{
    // Determine outside point.
    boundBox overallBb
    (
        point(GREAT, GREAT, GREAT),
        point(-GREAT, -GREAT, -GREAT)
    );

    bool hasSurface = false;

    const labelList unNamedSurfaces =
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfZones())
    );

    forAll(unNamedSurfaces, i)
    {
        label surfI = unNamedSurfaces[i];
        const searchableSurface& s = allGeometry_[surfaces_[surfI]];

        if (isA<triSurfaceMesh>(s) && s.hasVolumeType())
        {
            triSurfaceMesh& surf = const_cast<triSurfaceMesh&>
            (
                refCast<const triSurfaceMesh>(s)
            );

            if (reorient_)
            {
                orientedSurface::orientConsistent(surf);
            }
            surf.clearOut();
        }
    }

    const labelList namedSurfaces =
    (
        surfaceZonesInfo::getNamedSurfaces(surfZones())
    );

    forAll(namedSurfaces, i)
    {
        label surfI = namedSurfaces[i];
        const searchableSurface& s = allGeometry_[surfaces_[surfI]];

        if (isA<triSurfaceMesh>(s))
        {
            const triSurfaceMesh& surf = refCast<const triSurfaceMesh>(s);

            if (surf.triSurface::size() > 0)
            {
                tmp<pointField> tpoints(surf.points());
                const pointField& points = tpoints();

                hasSurface = true;

                boundBox surfBb(points[0], points[0]);
                // Assume surface is compact!
                for (label i = 0; i < points.size(); i++)
                {
                    const point& pt = points[i];
                    surfBb.min() = min(surfBb.min(), pt);
                    surfBb.max() = max(surfBb.max(), pt);
                }

                overallBb.min() = min(overallBb.min(), surfBb.min());
                overallBb.max() = max(overallBb.max(), surfBb.max());
            }
        }
    }

    if (hasSurface)
    {
        const point outsidePt(2*overallBb.max() - overallBb.min());

        //Info<< "Using point " << outsidePt << " to orient surf" << endl;

        forAll(namedSurfaces, i)
        {
            label surfI = namedSurfaces[i];
            const searchableSurface& s = allGeometry_[surfaces_[surfI]];

            if (isA<triSurfaceMesh>(s) && s.hasVolumeType())
            {
                triSurfaceMesh& surf = const_cast<triSurfaceMesh&>
                (
                    refCast<const triSurfaceMesh>(s)
                );
                orientSurface(outsidePt, surf);
                surf.clearOut();
            }
        }
    }
}


// Orient surface (if they're closed) before any searching is done.
void Foam::refinementSurfaces::orientSurface
(
    const point& keepPoint,
    triSurfaceMesh& s
)
{
    // Flip surface so keepPoint is outside.
    orientedSurface::orientConsistent(s);
    bool anyFlipped = orientedSurface::orient(s, keepPoint, true);

    if (anyFlipped)
    {
        // orientedSurface will have done a clearOut of the surface.
        // we could do a clearout of the triSurfaceMeshes::trees()
        // but these aren't affected by orientation (except for cached
        // sideness which should not be set at this point. !!Should
        // check!)

        Info<< "orientSurfaces : Flipped orientation of surface "
            << s.searchableSurface::name()
            << " so point " << keepPoint << " is outside." << endl;
    }
}


Foam::refinementSurfaces::refinementSurfaces
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
)
:
    allGeometry_(allGeometry),
    surfaces_(surfaces),
    names_(names),
    surfZones_(surfZones),
    regionOffset_(regionOffset),
    minLevel_(minLevel),
    maxLevel_(maxLevel),
    gapLevel_(gapLevel),
    perpendicularAngle_(perpendicularAngle),
    proxLevelIncr_(proxLevelIncr),
    proxDir_(proxDir),
    curvature_(curvature),
    patchInfo_(patchInfo.size())
{
    forAll(patchInfo_, pI)
    {
        if (patchInfo.set(pI))
        {
            patchInfo_.set(pI, patchInfo.set(pI, nullptr));
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// // Pre-calculate the refinement level for every element
// void Foam::refinementSurfaces::wantedRefinementLevel
// (
//     const shellSurfaces& shells,
//     const label surfI,
//     const List<pointIndexHit>& info,    // Indices
//     const pointField& ctrs,             // Representative coordinate
//     labelList& minLevelField
// ) const
// {
//     const searchableSurface& geom = allGeometry_[surfaces_[surfI]];
//
//     // Get per element the region
//     labelList region;
//     geom.getRegion(info, region);
//
//     // Initialise fields to region wise minLevel
//     minLevelField.setSize(ctrs.size());
//     minLevelField = -1;
//
//     forAll(minLevelField, i)
//     {
//         if (info[i].hit())
//         {
//             minLevelField[i] = minLevel(surfI, region[i]);
//         }
//     }
//
//     // Find out if triangle inside shell with higher level
//     // What level does shell want to refine fc to?
//     labelList shellLevel;
//     shells.findHigherLevel(ctrs, minLevelField, shellLevel);
//
//     forAll(minLevelField, i)
//     {
//         minLevelField[i] = max(minLevelField[i], shellLevel[i]);
//     }
// }


Foam::labelList Foam::refinementSurfaces::maxGapLevel() const
{
    labelList surfaceMax(surfaces_.size(), 0);

    forAll(surfaces_, surfI)
    {
        const wordList& regionNames = allGeometry_[surfaces_[surfI]].regions();

        forAll(regionNames, regionI)
        {
            label globalI = globalRegion(surfI, regionI);
            const FixedList<label, 3>& gapInfo = extendedGapLevel_[globalI];
            surfaceMax[surfI] = max(surfaceMax[surfI], gapInfo[2]);
        }
    }
    return surfaceMax;
}


// Precalculate the curvature for every point of a triangulated searchable
// surface.
void Foam::refinementSurfaces::setCurvatureFields
(
    const label nSmoothIter
)
{
    labelList curveSurfaces(curvatureSurfaces());
    if (curveSurfaces.empty())
    {
        return;
    }

    Info<<"Calculating curvature on selected surfaces"<<endl;

    scalarField curvTol(curvature().size(), -1);
    forAll(curveSurfaces, i)
    {
        label surfI = curveSurfaces[i];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        scalarField triCurv(0);

        if (isA<distributedTriSurfaceMesh>(geom))
        {
            const triSurfaceMesh& triMesh =
                refCast<const triSurfaceMesh>(geom);
            triSurface& s = const_cast<triSurface&>
            (
                refCast<const triSurface>(triMesh)
            );
            triCurv = s.calcCurvature(geom.name(),nSmoothIter);
        }
        else if (isA<triSurface>(geom))
        {
            if (Pstream::master())
            {
                const triSurfaceMesh& triMesh =
                    refCast<const triSurfaceMesh>(geom);
                triSurface& s = const_cast<triSurface&>
                (
                    refCast<const triSurface>(triMesh)
                );
                triCurv = s.calcCurvature(geom.name(),nSmoothIter);
            }
        }

        if (returnReduce(triCurv.size(), sumOp<label>()) > 0)
        {
            if (!isA<distributedTriSurfaceMesh>(geom))
            {
                Pstream::scatter(triCurv);
            }
            const_cast<searchableSurface&>
                (geom).setCurvaturePointField(triCurv);
        }
    }
}


// Precalculate the refinement level for every element of the searchable
// surface.
void Foam::refinementSurfaces::setMinLevelFields
(
    const autoPtr<shellSurfaces>& shellsPtr
)
{
    Info<<"Caching minimum level field"<<endl;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Cache the refinement level (max of surface level and shell level)
        // on a per-element basis. Only makes sense if there are lots of
        // elements. Possibly should have 'enough' elements to have fine
        // enough resolution but for now just make sure we don't catch e.g.
        // searchableBox (size=6)
        if (geom.globalSize() > 10 || isA<triSurface>(geom))
        {
            // Representative local coordinates and bounding sphere
            pointField ctrs;
            scalarField radiusSqr;
            geom.boundingSpheres(ctrs, radiusSqr);

            labelList minLevelField(ctrs.size(), 0);
            {
                // Get per element the region
                labelList region;

                if
                (
                    isA<triSurfaceMesh>(geom)
                    && !isA<distributedTriSurfaceMesh>(geom)
                )
                {
                    const triSurface& ts = refCast<const triSurface>(geom);
                    forAll(minLevelField, i)
                    {
                        minLevelField[i] = minLevel(surfI, ts[i].region());
                    }
                }
                else
                {
                     // Get the element index in a roundabout way. Problem is
                     // e.g. distributed surface where local indices differ
                     // from global ones (needed for getRegion call)
                     List<pointIndexHit> info;
                     geom.findNearest(ctrs, radiusSqr, info);
                     geom.getRegion(info, region);
                     // From the region get the surface-wise refinement level
                     forAll(minLevelField, i)
                     {
                         if (info[i].hit()) //Note: should not be necessary
                         {
                             minLevelField[i] = minLevel(surfI, region[i]);
                         }
                     }
                }
            }

            if (shellsPtr.valid())
            {
                const shellSurfaces& shells = shellsPtr();
                // Find out if triangle inside shell with higher level
                // What level does shell want to refine fc to?

                labelList shellLevel;
                //If running in parallel share the load
                if (Pstream::parRun() && !isA<distributedTriSurfaceMesh>(geom))
                {
                    label binSz = ctrs.size();

                    labelList localToGlobal = identity(binSz);
                    labelList localMinLevelField(binSz);
                    pointField localCtrs(binSz);

                    binSz = 0;
                    label start = Pstream::myProcNo();
                    forAll(ctrs, facei)
                    {
                        if (facei == start)
                        {
                            localToGlobal[binSz] = facei;
                            localMinLevelField[binSz] = minLevelField[facei];
                            localCtrs[binSz] = ctrs[facei];
                            binSz++;
                            start += Pstream::nProcs();
                        }
                    }

                    localToGlobal.setSize(binSz);
                    localMinLevelField.setSize(binSz);
                    localCtrs.setSize(binSz);

                    labelList localShellLevel;
                    shells.findHigherLevel
                    (
                        localCtrs,
                        localMinLevelField,
                        localShellLevel
                    );

                    //Rework local level into global level
                    shellLevel.setSize(ctrs.size());
                    shellLevel = -1;
                    forAll(localToGlobal, locali)
                    {
                        label globali = localToGlobal[locali];
                        shellLevel[globali] = localShellLevel[locali];
                    }
                    Pstream::listCombineReduce(shellLevel, maxOp<label>());
                }
                else
                {
                    shells.findHigherLevel(ctrs, minLevelField, shellLevel);
                }

                // In case of triangulated surfaces only cache value if triangle
                // centre and vertices are in same shell
                if (isA<triSurface>(geom))
                {
                    label nUncached = 0;

                    // Check if points differing from ctr level

                    const triSurface& ts = refCast<const triSurface>(geom);
                    const pointField& points = ts.points();

                    // Determine minimum expected level to avoid having to
                    // test lots of points
                    labelList minPointLevel(points.size(), labelMax);
                    forAll(shellLevel, triI)
                    {
                        const labelledTri& t = ts[triI];
                        label level = minLevelField[triI];//shellLevel[triI];
                        forAll(t, tI)
                        {
                            minPointLevel[t[tI]] =
                                min(minPointLevel[t[tI]], level);
                        }
                    }

                    labelList pointLevel;
                    if (Pstream::parRun() && !isA<distributedTriSurfaceMesh>(geom))
                    {
                        label binSz = points.size();
                        labelList localToGlobalPts = identity(binSz);
                        labelList localMinPointLevel(binSz);
                        pointField localPts(binSz);
                        binSz = 0;
                        label start = Pstream::myProcNo();
                        forAll(points, pointi)
                        {
                            if (pointi == start)
                            {
                                localToGlobalPts[binSz] = pointi;
                                localMinPointLevel[binSz] =
                                    minPointLevel[pointi];
                                localPts[binSz] = points[pointi];
                                binSz++;
                                start += Pstream::nProcs();
                            }
                        }

                        localToGlobalPts.setSize(binSz);
                        localMinPointLevel.setSize(binSz);
                        localPts.setSize(binSz);

                        labelList localPointLevel;
                        shells.findHigherLevel
                        (
                            localPts,
                            localMinPointLevel,
                            localPointLevel
                        );

                        pointLevel.setSize(points.size());
                        pointLevel = -1;
                        forAll(localToGlobalPts, locali)
                        {
                            label globali = localToGlobalPts[locali];
                            pointLevel[globali] = localPointLevel[locali];
                        }
                        Pstream::listCombineReduce(pointLevel, maxOp<label>());
                    }
                    else
                    {
                        // See if inside any shells with higher refinement level
                        shells.findHigherLevel
                        (
                            points,
                            minPointLevel,
                            pointLevel
                        );
                    }

                    // See if triangle centre values differ from triangle points
                    forAll(shellLevel, triI)
                    {
                        const labelledTri& t = ts[triI];
                        label fLevel = shellLevel[triI];
                        if
                        (
                            (pointLevel[t[0]] < fLevel)
                            || (pointLevel[t[1]] < fLevel)
                            || (pointLevel[t[2]] < fLevel)
                        )
                        {
                            // Mark as uncached
                            shellLevel[triI] = -1;
                            nUncached++;
                        }
                    }

                    if (isA<distributedTriSurfaceMesh>(geom))
                    {
                        reduce(nUncached, sumOp<label>());
                    }

                    Info<< "For geometry " << geom.name() << " detected "
                        << nUncached << " uncached triangles out of "
                        << geom.globalSize() << endl;
                }

                // Combine overall level field with current shell level.
                // Make sure to preserve -1 (from triSurfaceMeshes with
                // triangles partly inside/outside)
                forAll(minLevelField, i)
                {
                    if (min(minLevelField[i], shellLevel[i]) < 0)
                    {
                        minLevelField[i] = -1;
                    }
                    else
                    {
                        minLevelField[i] = max(minLevelField[i], shellLevel[i]);
                    }
                }
            }
            // Store minLevelField on surface
            const_cast<searchableSurface&>(geom).setField(minLevelField);
        }
    }
}


// Find intersections of edge. Return -1 or first surface with higher minLevel
// number.
void Foam::refinementSurfaces::findHigherIntersection
(
    const shellSurfaces& shells,

    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    labelList& surfaces,
    labelList& surfaceLevel,
    scalar edge0Len,
    bool checkCurvature, /*= false*/
    const bool threaded /*= false*/
) const
{
    surfaces.setSize(start.size());
    surfaces = -1;
    surfaceLevel.setSize(start.size());
    surfaceLevel = -1;

    if (surfaces_.empty())
    {
        return;
    }

    if (surfaces_.size() == 1)
    {
        // Optimisation: single segmented surface. No need to duplicate
        // point storage.

        label surfI = 0;

        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        List<pointIndexHit> intersectionInfo(start.size());
        geom.findLineAny(start, end, intersectionInfo, threaded);


        // Surface-based refinement level
        labelList surfaceOnlyLevel(start.size(), -1);
        {
            // Get per intersection the region
            labelList region;
            geom.getRegion(intersectionInfo, region, threaded);

            if (threaded && intersectionInfo.size() > 2048)
            {
#ifndef FOAM_USE_TBB
                WarningInFunction
                    << "FOAM not linked against TBB! Cannot run multithreaded."
                    << endl;
#else
                const auto grainSize =
                    std::max<size_t>
                    (
                        intersectionInfo.size() / tbb::this_task_arena::max_concurrency(),
                        1024
                    );

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, intersectionInfo.size(), grainSize),
                    [&](const tbb::blocked_range<size_t>& r)
                    {
                        for (size_t i = r.begin(); i < r.end(); ++i)
                        {
                            if (intersectionInfo[i].hit())
                            {
                                surfaceOnlyLevel[i] = minLevel(surfI, region[i]);
                            }
                        }
                    },
                    tbb::simple_partitioner()
                );
#endif
            }
            else
            {
                forAll(intersectionInfo, i)
                {
                    if (intersectionInfo[i].hit())
                    {
                        surfaceOnlyLevel[i] = minLevel(surfI, region[i]);
                    }
                }
            }

            if (checkCurvature) //TODO threaded
            {
                // See if a cached curvature field available
                scalarField curvatureField;
                geom.getCurvatureField(intersectionInfo, curvatureField);
                if (curvatureField.size() > 0)
                {
                    forAll(intersectionInfo, i)
                    {
                        if (intersectionInfo[i].hit())
                        {
                            scalar dTheta = curvatureField[i]
                                * edge0Len/(1<<currentLevel[i]);
                            scalar regionCurvature =
                                curvature(surfI, region[i]);
                            if
                            (
                                regionCurvature > 0
                                && dTheta > regionCurvature
                            )
                            {
                                label maxCurvatureLevel = min
                                (
                                    currentLevel[i]+1,
                                    maxLevel(surfI, region[i])
                                );

                                surfaceOnlyLevel[i] = max
                                (
                                    surfaceOnlyLevel[i],
                                    maxCurvatureLevel
                                );
                            }
                        }
                    }
                }
            } //if checkCurvature
        }

        // Get shell refinement level if higher
        const labelList localLevel
        (
            findHigherLevel
            (
                geom,
                shells,
                intersectionInfo,
                surfaceOnlyLevel, // starting level
                threaded
            )
        );

        // Combine localLevel with current level
        if (threaded && localLevel.size() > 2048)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    localLevel.size() / tbb::this_task_arena::max_concurrency(),
                    1024
                );

            tbb::parallel_for
            (
                tbb::blocked_range<size_t>(0, localLevel.size(), grainSize),
                [&](const tbb::blocked_range<size_t>& r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        if (localLevel[i] > currentLevel[i])
                        {
                            surfaces[i] = surfI;
                            surfaceLevel[i] = localLevel[i];
                        }
                    }
                },
                tbb::simple_partitioner()
            );
#endif
        }
        else
        {
            forAll(localLevel, i)
            {
                if (localLevel[i] > currentLevel[i])
                {
                    surfaces[i] = surfI;    // index of surface
                    surfaceLevel[i] = localLevel[i];
                }
            }
        }

        return;
    }



    // Work arrays
    pointField p0(start);
    pointField p1(end);
    labelList intersectionToPoint(identity(start.size()));
    List<pointIndexHit> intersectionInfo(start.size());

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        geom.findLineAny(p0, p1, intersectionInfo, threaded);


        // Surface-based refinement level
        labelList surfaceOnlyLevel(intersectionInfo.size(), -1);
        {
            // Get per intersection the region
            labelList region;
            geom.getRegion(intersectionInfo, region, threaded);

            if (threaded && intersectionInfo.size() > 2048)
            {
#ifndef FOAM_USE_TBB
                WarningInFunction
                    << "FOAM not linked against TBB! Cannot run multithreaded."
                    << endl;
#else
                const auto grainSize =
                    std::max<size_t>
                    (
                        intersectionInfo.size() / tbb::this_task_arena::max_concurrency(),
                        1024
                    );

                tbb::parallel_for
                (
                    tbb::blocked_range<size_t>(0, intersectionInfo.size(), grainSize),
                    [&](const tbb::blocked_range<size_t>& r)
                    {
                        for (size_t i = r.begin(); i < r.end(); ++i)
                        {
                            if (intersectionInfo[i].hit())
                            {
                                surfaceOnlyLevel[i] = minLevel(surfI, region[i]);
                            }
                        }
                    },
                    tbb::simple_partitioner()
                );
#endif
            }
            else
            {
                forAll(intersectionInfo, i)
                {
                    if (intersectionInfo[i].hit())
                    {
                        surfaceOnlyLevel[i] = minLevel(surfI, region[i]);
                    }
                }
            }

            if (checkCurvature) //TODO threaded
            {
                // See if a cached curvature field available
                scalarField curvatureField;
                geom.getCurvatureField(intersectionInfo, curvatureField);
                if (curvatureField.size() > 0)
                {
                    forAll(intersectionInfo, i)
                    {
                        if (intersectionInfo[i].hit())
                        {
                            scalar dTheta = curvatureField[i]
                                * edge0Len/(1<<currentLevel[i]);
                            scalar regionCurvature =
                                curvature(surfI, region[i]);
                            if
                            (
                                regionCurvature > 0
                                && dTheta > regionCurvature
                            )
                            {
                                label maxCurvatureLevel = min
                                (
                                    currentLevel[i]+1,
                                    maxLevel(surfI, region[i])
                                );

                                surfaceOnlyLevel[i] = max
                                (
                                    surfaceOnlyLevel[i],
                                    maxCurvatureLevel
                                );
                            }
                        }
                    }
                }
            } // if (checkCurvature)
        }


        // Get shell refinement level if higher
        const labelList localLevel
        (
            findHigherLevel
            (
                geom,
                shells,
                intersectionInfo,
                surfaceOnlyLevel,
                threaded
            )
        );


        // Combine localLevel with current level
        label missI = 0;

        if (false && threaded && localLevel.size() > 2048) //not thread-safe
        //if (threaded && localLevel.size() > 2048)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    localLevel.size() / tbb::this_task_arena::max_concurrency(),
                    1024
                );

            //tbb::spin_mutex mut;
            ////std::atomic<label> missIa(0);

            tbb::parallel_for
            (
                tbb::blocked_range<size_t>(0, localLevel.size(), grainSize),
                [&](const tbb::blocked_range<size_t>& r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        label pointI = intersectionToPoint[i];

                        if (localLevel[i] > currentLevel[pointI])
                        {
                            // Mark point for refinement
                            surfaces[pointI] = surfI;
                            surfaceLevel[pointI] = localLevel[i];
                        }
                        //else
                        //{
                        //    //label Id = missIa.fetch_add(1);
                        //    //p0[Id] = start[pointI];
                        //    //p1[Id] = end[pointI];
                        //    //intersectionToPoint[Id] = pointI;

                        //    tbb::spin_mutex::scoped_lock lock(mut);
                        //    p0[missI] = start[pointI];
                        //    p1[missI] = end[pointI];
                        //    intersectionToPoint[missI] = pointI;
                        //    intersectionInfo[missI] = intersectionInfo[pointI]; //GGG
                        //    missI++;
                        //}
                    }
                },
                tbb::simple_partitioner()
            );

            //using ETSp = tbb::enumerable_thread_specific<std::vector<point>>;
            //using ETSl = tbb::enumerable_thread_specific<std::vector<label>>;
            //using ETSh = tbb::enumerable_thread_specific<std::vector<pointIndexHit>>;

            //ETSp pp0;//{grainSize};
            //ETSp pp1;//{grainSize};
            //ETSl pintersectionToPoint;//{grainSize};
            //ETSh pintersectionInfo;//{grainSize};

            //tbb::parallel_for
            //(
            //    tbb::blocked_range<size_t>(0, localLevel.size(), grainSize),
            //    [&](const tbb::blocked_range<size_t>& r)
            //    {
            //        ETSp::reference myP0 = pp0.local();
            //        ETSp::reference myP1 = pp1.local();
            //        ETSl::reference myIntersectionToPoint = pintersectionToPoint.local();
            //        ETSh::reference myIntersectionInfo = pintersectionInfo.local();

            //        for (size_t i = r.begin(); i < r.end(); ++i)
            //        {
            //            label pointI = intersectionToPoint[i];

            //            if (localLevel[i] > currentLevel[pointI])
            //            {
            //                // Mark point for refinement
            //                surfaces[pointI] = surfI;
            //                surfaceLevel[pointI] = localLevel[i];
            //            }
            //            else
            //            {
            //                myP0.emplace_back(start[pointI]);
            //                myP1.emplace_back(end[pointI]);
            //                myIntersectionToPoint.emplace_back(pointI);
            //                //myIntersectionInfo.emplace_back(intersectionInfo[pointI]);
            //            }
            //        }
            //    },
            //    tbb::simple_partitioner()
            //);

            //auto it0 = pp1.begin();
            //auto it1 = pintersectionToPoint.begin();
            //auto it2 = pintersectionInfo.begin();
            //for (auto it = pp0.begin(); it != pp0.end(); ++it)
            //{
            //    for (label j = 0; j < (*it).size(); ++j)
            //    {
            //        p0[missI] = (*it)[j];
            //        p1[missI] = (*it0)[j];
            //        intersectionToPoint[missI] = (*it1)[j];
            //        intersectionInfo[missI] = (*it2)[j];
            //        ++missI;
            //    }
            //    ++it0;
            //    ++it1;
            //    ++it2;
            //}
#endif
        }
        else
        {
            forAll(localLevel, i)
            {
                label pointI = intersectionToPoint[i];

                if (localLevel[i] > currentLevel[pointI])
                {
                    // Mark point for refinement
                    surfaces[pointI] = surfI;
                    surfaceLevel[pointI] = localLevel[i];
                }
                else
                {
                    p0[missI] = start[pointI];
                    p1[missI] = end[pointI];
                    intersectionToPoint[missI] = pointI;
                    ++missI;
                }
            }
        }

        // All done? Note that this decision should be synchronised
        if (returnReduce(missI, sumOp<label>()) == 0)
        {
            break;
        }

        // Trim misses
        p0.setSize(missI);
        p1.setSize(missI);
        intersectionToPoint.setSize(missI);
        intersectionInfo.setSize(missI);
    } //forAll surfaces
} //findHigherIntersection


void Foam::refinementSurfaces::findAllHigherIntersections
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    const labelList& globalRegionLevel,
    List<List<vector>>& surfaceNormal,

    labelListList& surfaceLevel
) const
{
    surfaceLevel.setSize(start.size());
    surfaceNormal.setSize(start.size());

    if (surfaces_.empty())
    {
        return;
    }

    // Work arrays
    List<List<pointIndexHit>> hitInfo;
    labelList pRegions;
    vectorField pNormals;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        surface.findLineAll(start, end, hitInfo);

        // Repack hits for surface into flat list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // To avoid overhead of calling getRegion for every point

        label n = 0;
        forAll(hitInfo, pointI)
        {
            n += hitInfo[pointI].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointI)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointI];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointI;
                n++;
            }
        }

        labelList surfRegion(n);
        vectorField surfNormal(n);
        surface.getRegion(surfInfo, surfRegion);
        surface.getNormal(surfInfo, surfNormal);

        surfInfo.clear();


        // Extract back into pointwise
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(surfRegion, i)
        {
            label region = globalRegion(surfI, surfRegion[i]);
            label pointI = pointMap[i];

            if (globalRegionLevel[region] > currentLevel[pointI])
            {
                // Append to pointI info
                label sz = surfaceNormal[pointI].size();
                surfaceNormal[pointI].setSize(sz+1);
                surfaceNormal[pointI][sz] = surfNormal[i];

                surfaceLevel[pointI].setSize(sz+1);
                surfaceLevel[pointI][sz] = globalRegionLevel[region];
            }
        }
    }
}


void Foam::refinementSurfaces::findAllHigherIntersections
(
    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    const labelList& globalRegionLevel,

    List<pointList>& surfaceLocation,
    List<vectorList>& surfaceNormal,
    labelListList& surfaceLevel
) const
{
    surfaceLevel.setSize(start.size());
    surfaceNormal.setSize(start.size());
    surfaceLocation.setSize(start.size());

    if (surfaces_.empty())
    {
        return;
    }

    // Work arrays
    List<List<pointIndexHit>> hitInfo;
    labelList pRegions;
    vectorField pNormals;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        surface.findLineAll(start, end, hitInfo);

        // Repack hits for surface into flat list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // To avoid overhead of calling getRegion for every point

        label n = 0;
        forAll(hitInfo, pointI)
        {
            n += hitInfo[pointI].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointI)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointI];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointI;
                n++;
            }
        }

        labelList surfRegion(n);
        vectorField surfNormal(n);
        surface.getRegion(surfInfo, surfRegion);
        surface.getNormal(surfInfo, surfNormal);

        // Extract back into pointwise
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(surfRegion, i)
        {
            label region = globalRegion(surfI, surfRegion[i]);
            label pointI = pointMap[i];

            if (globalRegionLevel[region] > currentLevel[pointI])
            {
                // Append to pointI info
                label sz = surfaceNormal[pointI].size();
                surfaceLocation[pointI].setSize(sz+1);
                surfaceLocation[pointI][sz] = surfInfo[i].hitPoint();

                surfaceNormal[pointI].setSize(sz+1);
                surfaceNormal[pointI][sz] = surfNormal[i];

                surfaceLevel[pointI].setSize(sz+1);
                surfaceLevel[pointI][sz] = globalRegionLevel[region];
            }
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
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
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        surface.findLine
        (
            start,
            nearest,
            nearestInfo
        );
        surface.getRegion
        (
            nearestInfo,
            region
        );

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit1[pointI] = nearestInfo[pointI];
                surface1[pointI] = surfI;
                region1[pointI] = region[pointI];
                nearest[pointI] = hit1[pointI].hitPoint();
            }
        }
    }


    // 2. intersection from end to last intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = hit1;
    region2 = region1;

    // Set current end of segment to test.
    forAll(nearest, pointI)
    {
        if (hit1[pointI].hit())
        {
            nearest[pointI] = hit1[pointI].hitPoint();
        }
        else
        {
            // Disable testing by setting to end.
            nearest[pointI] = end[pointI];
        }
    }

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        // See if any intersection between end and current nearest
        surface.findLine
        (
            end,
            nearest,
            nearestInfo
        );
        surface.getRegion
        (
            nearestInfo,
            region
        );

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit2[pointI] = nearestInfo[pointI];
                surface2[pointI] = surfI;
                region2[pointI] = region[pointI];
                nearest[pointI] = hit2[pointI].hitPoint();
            }
        }
    }


    // Make sure that if hit1 has hit something, hit2 will have at least the
    // same point (due to tolerances it might miss its end point)
    forAll(hit1, pointI)
    {
        if (hit1[pointI].hit() && !hit2[pointI].hit())
        {
            hit2[pointI] = hit1[pointI];
            surface2[pointI] = surface1[pointI];
            region2[pointI] = region1[pointI];
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
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
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());
    region1 = -1;
    normal1.setSize(start.size());
    normal1 = Zero;

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit1[pointI] = nearestInfo[pointI];
                surface1[pointI] = surfI;
                region1[pointI] = region[pointI];
                normal1[pointI] = normal[pointI];
                nearest[pointI] = hit1[pointI].hitPoint();
            }
        }
    }


    // 2. intersection from end to last intersection
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the nearest intersection from end to start. Note that we initialize
    // to the first intersection (if any).
    surface2 = surface1;
    hit2 = hit1;
    region2 = region1;
    normal2 = normal1;

    // Set current end of segment to test.
    forAll(nearest, pointI)
    {
        if (hit1[pointI].hit())
        {
            nearest[pointI] = hit1[pointI].hitPoint();
        }
        else
        {
            // Disable testing by setting to end.
            nearest[pointI] = end[pointI];
        }
    }

    forAll(surfacesToTest, testI)
    {
        label surfI = surfacesToTest[testI];
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between end and current nearest
        geom.findLine(end, nearest, nearestInfo);
        geom.getRegion(nearestInfo, region);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                hit2[pointI] = nearestInfo[pointI];
                surface2[pointI] = surfI;
                region2[pointI] = region[pointI];
                normal2[pointI] = normal[pointI];
                nearest[pointI] = hit2[pointI].hitPoint();
            }
        }
    }


    // Make sure that if hit1 has hit something, hit2 will have at least the
    // same point (due to tolerances it might miss its end point)
    forAll(hit1, pointI)
    {
        if (hit1[pointI].hit() && !hit2[pointI].hit())
        {
            hit2[pointI] = hit1[pointI];
            surface2[pointI] = surface1[pointI];
            region2[pointI] = region1[pointI];
            normal2[pointI] = normal1[pointI];
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    vectorField& normal1
) const
{
    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    normal1.setSize(start.size());
    normal1 = Zero;

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                surface1[pointI] = surfI;
                normal1[pointI] = normal[pointI];
                nearest[pointI] = nearestInfo[pointI].hitPoint();
            }
        }
    }
}


void Foam::refinementSurfaces::findNearestIntersection
(
    const pointField& start,
    const pointField& end,

    labelList& surface1,
    List<pointIndexHit>& hitInfo1,
    vectorField& normal1
) const
{
    // 1. intersection from start to end
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hitInfo1.setSize(start.size());
    hitInfo1 = pointIndexHit();
    normal1.setSize(start.size());
    normal1 = Zero;

    // Current end of segment to test.
    pointField nearest(end);
    // Work array
    List<pointIndexHit> nearestInfo(start.size());
    labelList region;
    vectorField normal;

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // See if any intersection between start and current nearest
        geom.findLine(start, nearest, nearestInfo);
        geom.getNormal(nearestInfo, normal);

        forAll(nearestInfo, pointI)
        {
            if (nearestInfo[pointI].hit())
            {
                surface1[pointI] = surfI;
                hitInfo1[pointI] = nearestInfo[pointI];
                normal1[pointI] = normal[pointI];
                nearest[pointI] = nearestInfo[pointI].hitPoint();
            }
        }
    }
}


void Foam::refinementSurfaces::findAnyIntersection
(
    const pointField& start,
    const pointField& end,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        surfaces_,
        start,
        end,
        hitSurface,
        hitInfo
    );
}


void Foam::refinementSurfaces::findNearest
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    labelList& hitRegion
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    List<pointIndexHit> hitInfo;
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const pointField& samples,
    const scalarField& nearestDistSqr,
    labelList& hitSurface,
    List<pointIndexHit>& hitInfo,
    labelList& hitRegion,
    vectorField& hitNormal
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;
    hitNormal.setSize(hitSurface.size());
    hitNormal = Zero;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        // Region
        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }

        // Normal
        vectorField localNormal;
        allGeometry_[surfaces_[surfI]].getNormal(localHits, localNormal);

        forAll(localIndices, i)
        {
            hitNormal[localIndices[i]] = localNormal[i];
        }
    }
}


void Foam::refinementSurfaces::findInside
(
    const labelList& testSurfaces,
    const pointField& pt,
    labelList& insideSurfaces
) const
{
    insideSurfaces.setSize(pt.size());
    insideSurfaces = -1;

    forAll(testSurfaces, i)
    {
        label surfI = testSurfaces[i];

        const surfaceZonesInfo::areaSelectionAlgo selectionMethod =
            surfZones_[surfI].zoneInside();

        if
        (
            selectionMethod != surfaceZonesInfo::INSIDE
            && selectionMethod != surfaceZonesInfo::OUTSIDE
        )
        {
            FatalErrorInFunction
                << "Trying to use surface "
                << allGeometry_[surfaces_[surfI]].name()
                << " which has non-geometric inside selection method "
                << surfaceZonesInfo::areaSelectionAlgoNames[selectionMethod]
                << exit(FatalError);
        }

        if (allGeometry_[surfaces_[surfI]].hasVolumeType())
        {
            List<volumeType> volType;
            allGeometry_[surfaces_[surfI]].getVolumeType(pt, volType);

            forAll(volType, pointI)
            {
                if (insideSurfaces[pointI] == -1)
                {
                    if
                    (
                        (
                            volType[pointI] == volumeType::INSIDE
                         && selectionMethod == surfaceZonesInfo::INSIDE
                        )
                     || (
                            volType[pointI] == volumeType::OUTSIDE
                         && selectionMethod == surfaceZonesInfo::OUTSIDE
                        )
                    )
                    {
                        insideSurfaces[pointI] = surfI;
                    }
                }
            }
        }
    }
}


void Foam::refinementSurfaces::findNearest
(
    const labelList& surfacesToTest,
    const labelListList& regions,

    const pointField& samples,
    const scalarField& nearestDistSqr,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        regions,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }
}


void Foam::refinementSurfaces::findNearestRegion
(
    const labelList& surfacesToTest,
    const labelListList& regions,

    const pointField& samples,
    const scalarField& nearestDistSqr,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo,
    labelList& hitRegion,
    vectorField& hitNormal
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    // Do the tests. Note that findNearest returns index in geometries.
    searchableSurfacesQueries::findNearest
    (
        allGeometry_,
        geometries,
        regions,
        samples,
        nearestDistSqr,
        hitSurface,
        hitInfo
    );

    // Rework the hitSurface to be surface (i.e. index into surfaces_)
    forAll(hitSurface, pointI)
    {
        if (hitSurface[pointI] != -1)
        {
            hitSurface[pointI] = surfacesToTest[hitSurface[pointI]];
        }
    }

    // Collect the region
    hitRegion.setSize(hitSurface.size());
    hitRegion = -1;
    hitNormal.setSize(hitSurface.size());
    hitNormal = Zero;

    forAll(surfacesToTest, i)
    {
        label surfI = surfacesToTest[i];

        // Collect hits for surfI
        const labelList localIndices(findIndices(hitSurface, surfI));

        List<pointIndexHit> localHits
        (
            UIndirectList<pointIndexHit>
            (
                hitInfo,
                localIndices
            )
        );

        // Region
        labelList localRegion;
        allGeometry_[surfaces_[surfI]].getRegion(localHits, localRegion);

        forAll(localIndices, i)
        {
            hitRegion[localIndices[i]] = localRegion[i];
        }

        // Normal
        vectorField localNormal;
        allGeometry_[surfaces_[surfI]].getNormal(localHits, localNormal);

        forAll(localIndices, i)
        {
            hitNormal[localIndices[i]] = localNormal[i];
        }
    }
}


void Foam::refinementSurfaces::findNearestPerturbedIntersection
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
) const
{
    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());

    surface2.setSize(start.size());
    surface2 = -1;
    hit2.setSize(start.size());
    region2.setSize(start.size());

    labelList surface1p;
    List<pointIndexHit> hit1p;
    labelList region1p;

    labelList surface2p;
    List<pointIndexHit> hit2p;
    labelList region2p;

    findNearestIntersection
    (
        surfacesToTest,
        start,
        end,

        surface1p,
        hit1p,
        region1p,

        surface2p,
        hit2p,
        region2p
    );

    label nMisses = 0;
    labelList missedFaces(start.size());
    forAll(hit1p, i)
    {
        if (hit1p[i].hit() && hit2p[i].hit())
        {
            surface1[i] = surface1p[i];
            hit1[i] = hit1p[i];
            region1[i] = region1p[i];

            surface2[i] = surface2p[i];
            hit2[i] = hit2p[i];
            region2[i] = region2p[i];
        }
        else
        {
            missedFaces[nMisses++] = i;
        }
    }
    missedFaces.setSize(nMisses);

    label nOuter(addedRays ? 2 : 1);

    pointField startPerturbed(nMisses);
    pointField endPerturbed(nMisses);
    for (int outer = 0; outer < nOuter; outer++)
    {
        for (int i = 0; i < 2; i++)
        {
            for (direction j = 0; j < vector::nComponents; j++)
            {
                vector dir = vector(0., 0., 0.);
                if (i == 0)
                {
                    dir[j] = tol;
                }
                else
                {
                    dir[j] = -tol;
                }

                forAll(startPerturbed, hitI)
                {
                    label missedFaceI = missedFaces[hitI];

                    if (addedRays)
                    {
                        if (outer == 0)
                        {
                            startPerturbed[hitI] = start[missedFaceI];
                            endPerturbed[hitI] = end[missedFaceI] + dir;
                        }
                        else if (outer == 1)
                        {
                            startPerturbed[hitI] = start[missedFaceI] +dir;
                            endPerturbed[hitI] = end[missedFaceI];
                        }
                    }
                    else
                    {
                        startPerturbed[hitI] = start[missedFaceI] + dir;
                        endPerturbed[hitI] = end[missedFaceI] + dir;
                    }
                }

                // Initialize arguments
                surface1p.setSize(nMisses);
                hit1p.setSize(nMisses);
                region1p.setSize(nMisses);

                surface2p.setSize(nMisses);
                hit2p.setSize(nMisses);
                region2p.setSize(nMisses);
                surface1p = -1;
                surface2p = -1;
                hit1p = pointIndexHit();
                hit2p = pointIndexHit();

                findNearestIntersection
                (
                    surfacesToTest,
                    startPerturbed,
                    endPerturbed,

                    surface1p,
                    hit1p,
                    region1p,

                    surface2p,
                    hit2p,
                    region2p
                );

                forAll(hit1p, hitI)
                {
                    label missedFaceI = missedFaces[hitI];
                    if
                    (
                        !hit1[missedFaceI].hit() && !hit2[missedFaceI].hit()
                        && hit1p[hitI].hit() && hit2p[hitI].hit()
                     )
                    {
                        surface1[missedFaceI] = surface1p[hitI];
                        hit1[missedFaceI] = hit1p[hitI];
                        region1[missedFaceI] = region1p[hitI];

                        surface2[missedFaceI] = surface2p[hitI];
                        hit2[missedFaceI] = hit2p[hitI];
                        region2[missedFaceI] = region2p[hitI];
                    }
                }
            }
        }
    }
}


void Foam::refinementSurfaces::findNearestPerturbedIntersection
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
) const
{
    // Initialize arguments
    surface1.setSize(start.size());
    surface1 = -1;
    hit1.setSize(start.size());
    region1.setSize(start.size());
    normal1.setSize(start.size());

    surface2.setSize(start.size());
    surface2 = -1;
    hit2.setSize(start.size());
    region2.setSize(start.size());
    normal2.setSize(start.size());

    labelList surface1p;
    List<pointIndexHit> hit1p;
    labelList region1p;
    vectorField normal1p;

    labelList surface2p;
    List<pointIndexHit> hit2p;
    labelList region2p;
    vectorField normal2p;

    findNearestIntersection
    (
        surfacesToTest,
        start,
        end,

        surface1p,
        hit1p,
        region1p,
        normal1p,

        surface2p,
        hit2p,
        region2p,
        normal2p
    );

    label nMisses = 0;
    labelList missedFaces(start.size());
    forAll(hit1p, i)
    {
        if (hit1p[i].hit() && hit2p[i].hit())
        {
            surface1[i] = surface1p[i];
            hit1[i] = hit1p[i];
            region1[i] = region1p[i];
            normal1[i] = normal1p[i];

            surface2[i] = surface2p[i];
            hit2[i] = hit2p[i];
            region2[i] = region2p[i];
            normal2[i] = normal2p[i];
        }
        else
        {
            missedFaces[nMisses++] = i;
        }
    }

    missedFaces.setSize(nMisses);
    pointField startPerturbed(nMisses);
    pointField endPerturbed(nMisses);

    label nOuter(addedRays ? 2 : 1);

    for (int outer = 0; outer < nOuter; outer++)
    {
        for (int i = 0; i < 2; i++)
        {
            for (direction j = 0; j < vector::nComponents; j++)
            {
                vector dir = vector(0., 0., 0.);
                if (i == 0)
                {
                    dir[j] = tol;
                }
                else
                {
                    dir[j] = -tol;
                }

                forAll(startPerturbed, hitI)
                {
                    label missedFaceI = missedFaces[hitI];
                    if (addedRays)
                    {
                        if (outer == 0)
                        {
                            startPerturbed[hitI] = start[missedFaceI];
                            endPerturbed[hitI] = end[missedFaceI] + dir;
                        }
                        else if (outer == 1)
                        {
                            startPerturbed[hitI] = start[missedFaceI] +dir;
                            endPerturbed[hitI] = end[missedFaceI];
                        }
                    }
                    else
                    {
                        startPerturbed[hitI] = start[missedFaceI] + dir;
                        endPerturbed[hitI] = end[missedFaceI] + dir;
                    }
                }

                // Initialize arguments
                surface1p.setSize(nMisses);
                hit1p.setSize(nMisses);
                region1p.setSize(nMisses);
                normal1p.setSize(nMisses);

                surface2p.setSize(nMisses);
                hit2p.setSize(nMisses);
                region2p.setSize(nMisses);
                normal2p.setSize(nMisses);

                surface1p = -1;
                surface2p = -1;
                hit1p = pointIndexHit();
                hit2p = pointIndexHit();

                findNearestIntersection
                (
                    surfacesToTest,
                    startPerturbed,
                    endPerturbed,

                    surface1p,
                    hit1p,
                    region1p,
                    normal1p,

                    surface2p,
                    hit2p,
                    region2p,
                    normal2p
                 );

                forAll(hit1p, hitI)
                {
                    label missedFaceI = missedFaces[hitI];
                    if
                    (
                        (!hit1[missedFaceI].hit() && !hit2[missedFaceI].hit())
                        && (hit1p[hitI].hit() && hit2p[hitI].hit())
                    )
                    {
                        surface1[missedFaceI] = surface1p[hitI];
                        hit1[missedFaceI] = hit1p[hitI];
                        region1[missedFaceI] = region1p[hitI];
                        normal1[missedFaceI] = normal1p[hitI];

                        surface2[missedFaceI] = surface2p[hitI];
                        hit2[missedFaceI] = hit2p[hitI];
                        region2[missedFaceI] = region2p[hitI];
                        normal2[missedFaceI] = normal2p[hitI];
                    }
                }
            }
        }
    }
}


// Find intersections of edge. Return -1 or first surface with higher minLevel
// number.
void Foam::refinementSurfaces::findHigherPerturbedIntersection
(
    const shellSurfaces& shells,

    const pointField& start,
    const pointField& end,
    const labelList& currentLevel,   // current cell refinement level

    labelList& surfaces,
    labelList& surfaceLevel,
    const scalar& tol,
    const bool& addedRays
) const
{
    surfaces.setSize(start.size());
    surfaces = -1;
    surfaceLevel.setSize(start.size());
    surfaceLevel = -1;

    if (surfaces_.empty())
    {
        return;
    }

    label nOuter(addedRays ? 2 : 1);

    if (surfaces_.size() == 1)
    {
        // Optimisation: single segmented surface. No need to duplicate
        // point storage.

        label surfI = 0;

        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];

        // Do intersection test
        List<pointIndexHit> intersectionInfo(start.size());

        labelList surfaceOnlyLevel(start.size(), -1);

        for (int outer = 0; outer < nOuter; outer++)
        {
            for (int i = 0; i < 2; i++)
            {
                for (direction j = 0; j < vector::nComponents; j++)
                {
                    vector dir = vector(0., 0., 0.);
                    if (i == 0)
                    {
                        dir[j] = tol;
                    }
                    else
                    {
                        dir[j] = -tol;
                    }

                    pointField startPerturbed = start;
                    pointField endPerturbed = end;

                    if (addedRays)
                    {
                        if (outer == 0)
                        {
                            startPerturbed = start;
                            endPerturbed = end + dir;
                        }
                        else if (outer == 1)
                        {
                            startPerturbed = start + dir;
                            endPerturbed = end;
                        }
                    }
                    else
                    {
                        startPerturbed = start + dir;
                        endPerturbed = end + dir;
                    }

                    geom.findLineAny(start, end, intersectionInfo);

                    // Surface-based refinement level
                    {
                        // Get per intersection the region
                        labelList region;
                        geom.getRegion(intersectionInfo, region);

                        forAll(intersectionInfo, i)
                        {
                            if (intersectionInfo[i].hit())
                            {
                                surfaceOnlyLevel[i] =
                                    minLevel(surfI, region[i]);
                            }
                        }
                    }
                }
            }
        }

        // Get shell refinement level if higher
        const labelList localLevel
        (
            findHigherLevel
            (
                geom,
                shells,
                intersectionInfo,
                surfaceOnlyLevel // starting level
            )
        );


        // Combine localLevel with current level
        forAll(localLevel, i)
        {
            if (localLevel[i] > currentLevel[i])
            {
                surfaces[i] = surfI;    // index of surface
                surfaceLevel[i] = localLevel[i];
            }
        }
        return;
    }

    // Work arrays
    pointField p0(start);
    pointField p1(end);
    labelList intersectionToPoint(identity(start.size()));
    List<pointIndexHit> intersectionInfo(start.size());

    forAll(surfaces_, surfI)
    {
        const searchableSurface& geom = allGeometry_[surfaces_[surfI]];
        labelList surfaceOnlyLevel(intersectionInfo.size(), -1);

        for (int outer = 0; outer < nOuter; outer++)
        {
            for (int i = 0; i < 2; i++)
            {
                for (direction j = 0; j < vector::nComponents; j++)
                {
                    vector dir = vector(0., 0., 0.);
                    if (i == 0)
                    {
                        dir[j] = tol;
                    }
                    else
                    {
                        dir[j] = -tol;
                    }

                    pointField startPerturbed = start;
                    pointField endPerturbed = end;

                    if (addedRays)
                    {
                        if (outer == 0)
                        {
                            startPerturbed = start;
                            endPerturbed = end + dir;
                        }
                        else if (outer == 1)
                        {
                            startPerturbed = start + dir;
                            endPerturbed = end;
                        }
                    }
                    else
                    {
                        startPerturbed = start + dir;
                        endPerturbed = end + dir;
                    }

                    // Do intersection test
                    geom.findLineAny(p0, p1, intersectionInfo);

                    // Surface-based refinement level
                    {
                        // Get per intersection the region
                        labelList region;
                        geom.getRegion(intersectionInfo, region);

                        forAll(intersectionInfo, i)
                        {
                            if (intersectionInfo[i].hit())
                            {
                                surfaceOnlyLevel[i] =
                                    minLevel(surfI, region[i]);
                            }
                        }
                    }
                }
            }
        }

        // Get shell refinement level if higher
        const labelList localLevel
        (
            findHigherLevel
            (
                geom,
                shells,
                intersectionInfo,
                surfaceOnlyLevel
            )
        );

        // Combine localLevel with current level
        label missI = 0;
        forAll(localLevel, i)
        {
            label pointI = intersectionToPoint[i];

            if (localLevel[i] > currentLevel[pointI])
            {
                // Mark point for refinement
                surfaces[pointI] = surfI;
                surfaceLevel[pointI] = localLevel[i];
            }
            else
            {
                p0[missI] = start[pointI];
                p1[missI] = end[pointI];
                intersectionToPoint[missI] = pointI;
                missI++;
            }
        }


        // All done? Note that this decision should be synchronised
        if (returnReduce(missI, sumOp<label>()) == 0)
        {
            break;
        }

        // Trim misses
        p0.setSize(missI);
        p1.setSize(missI);
        intersectionToPoint.setSize(missI);
        intersectionInfo.setSize(missI);
    }
}


void Foam::refinementSurfaces::findCellZone
(
    const pointField& pts,
    const regionSplit& globalRegions,
    const labelList& surfaceToCellZone,
    labelList& cellToZone,
    label minRegionSize
) const
{
    labelList nCellsPerRegion(globalRegions.nRegions(), 0);

    forAll(globalRegions, ptI)
    {
        label regionI = globalRegions[ptI];
        nCellsPerRegion[regionI]++;
    }

    Pstream::listCombineReduce(nCellsPerRegion, plusOp<label>());

    const labelList closedNamedSurfaces =
    (
        surfaceZonesInfo::getClosedNamedSurfaces(surfZones(),geometry(),surfaces())
    );

    forAll(closedNamedSurfaces, i)
    {
        label surfI = closedNamedSurfaces[i];
        label zoneI = surfaceToCellZone[surfI];
        if (zoneI == -1)
        {
            continue;
        }
        const searchableSurface& s = allGeometry_[surfaces_[surfI]];
        if (s.hasVolumeType())
        {

            const surfaceZonesInfo::areaSelectionAlgo selectionMethod =
                surfZones_[surfI].zoneInside();

            if (isA<triSurfaceMesh>(s))
            {
                boolList visited(pts.size(), false);

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

                // choose 3 end points to make sure initial insideness
                //check is robust
                point midBbPt = 0.5*(overallBb.max() + overallBb.min());
                scalar dr = overallBb.mag();

                pointField endPoints(3);
                endPoints[0] = midBbPt + point(dr, 0., 0.);
                endPoints[1] = midBbPt + point(0., dr, 0.);
                endPoints[2] = midBbPt + point(0., 0., dr);

                pointField start(1);
                pointField end(1);

                forAll(pts, i)
                {
                    label regionI = globalRegions[i];

                    if (!visited[i] && nCellsPerRegion[regionI] > minRegionSize)
                    {
                        visited[i] = true;
                        start[0] = pts[i];
                        label nIn = 0;
                        label nOut = 0;
                        forAll(endPoints, pointI)
                        {
                            end[0] = endPoints[pointI];
                            List<List<pointIndexHit>> hitInfo(1);
                            s.findLineAll(start, end, hitInfo);

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
                        bool foundZone = false;

                        if (nOut > nIn)
                        {
                            if (selectionMethod == surfaceZonesInfo::OUTSIDE)
                            {
                                if (cellToZone[i] < 0)
                                {
                                    cellToZone[i] = zoneI;
                                    foundZone = true;
                                }
                            }
                            for (label j=i+1; j < pts.size(); j++)
                            {
                                if (globalRegions[j] == regionI)
                                {
                                    visited[j] = true;
                                    if (foundZone)
                                    {
                                        cellToZone[j] = zoneI;
                                    }
                                }
                            }
                        }
                        else
                        {
                            if (selectionMethod == surfaceZonesInfo::INSIDE)
                            {
                                if (cellToZone[i] < 0)
                                {
                                    cellToZone[i] = zoneI;
                                    foundZone = true;
                                }
                            }
                            for (label j=i+1; j < pts.size(); j++)
                            {
                                if (globalRegions[j] == regionI)
                                {
                                    visited[j] = true;
                                    if (foundZone)
                                    {
                                        cellToZone[j] = zoneI;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                List<volumeType> volType;
                s.getVolumeType(pts, volType);
                forAll(volType, pointI)
                {
                    if
                    (
                        (
                            volType[pointI] == volumeType::INSIDE
                            && selectionMethod == surfaceZonesInfo::INSIDE
                         )
                        ||
                        (
                            volType[pointI] == volumeType::OUTSIDE
                            && selectionMethod == surfaceZonesInfo::OUTSIDE
                         )
                     )
                    {
                        if (cellToZone[pointI] < 0)
                        {
                            cellToZone[pointI] = zoneI;
                        }
                    }
                }
            }
        }
    }
}


void Foam::refinementSurfaces::findPerturbedCellZone
(
    const labelList& testSurfaces,
    const pointField& pt,
    labelList& insideSurfaces
) const
{
    insideSurfaces.setSize(pt.size());
    insideSurfaces = -1;

    pointField start(3*pt.size());
    pointField end(3*pt.size());

    forAll(testSurfaces, i)
    {
        label surfI = testSurfaces[i];
        const searchableSurface& s = allGeometry_[surfaces_[surfI]];
        if (s.hasVolumeType())
        {
            const surfaceZonesInfo::areaSelectionAlgo selectionMethod =
                surfZones_[surfI].zoneInside();

            if (isA<triSurfaceMesh>(s))
            {
                triSurfaceMesh& surf = const_cast<triSurfaceMesh&>
                (
                    refCast<const triSurfaceMesh>(s)
                );
                tmp<pointField> tpoints(surf.points());
                const pointField& points = tpoints();
                boundBox surfBb(points, false);

                boundBox overallBb
                (
                    point(GREAT, GREAT, GREAT),
                    point(-GREAT, -GREAT, -GREAT)
                );
                overallBb.min() = min(overallBb.min(), surfBb.min());
                overallBb.max() = max(overallBb.max(), surfBb.max());

                // choose 3 end points to make sure initial insideness
                //check is robust
                point midBbPt = 0.5*(overallBb.max() + overallBb.min());
                scalar dr = overallBb.mag();

                label k = 0;
                forAll(pt, pointI)
                {
                    start[k] = pt[pointI];
                    end[k++] = midBbPt + point(dr, 0., 0.);
                    start[k] = pt[pointI];
                    end[k++] = midBbPt + point(0., dr, 0.);
                    start[k] = pt[pointI];
                    end[k++] = midBbPt + point(0., 0., dr);
                }
                List<List<pointIndexHit>> hitInfo(start.size());
                s.findLineAll(start, end, hitInfo);

                k = 0;
                forAll(pt, pointI)
                {
                    label nIn = 0;
                    label nOut = 0;
                    for (int j = 0; j < 3; j++)
                    {
                        if (hitInfo[k].size() % 2 == 0)
                        {
                            nOut++;
                        }
                        else
                        {
                            nIn++;
                        }
                        k++;
                    }

                    if (insideSurfaces[pointI] == -1)
                    {
                        if (nOut > nIn)
                        {
                            if (selectionMethod == surfaceZonesInfo::OUTSIDE)
                            {
                                insideSurfaces[pointI] = surfI;
                            }
                        }
                        else
                        {
                            if (selectionMethod == surfaceZonesInfo::INSIDE)
                            {
                                insideSurfaces[pointI] = surfI;
                            }
                        }
                    }
                }
            }
            else
            {
                List<volumeType> volType;
                s.getVolumeType(pt, volType);

                forAll(volType, pointI)
                {
                    if (insideSurfaces[pointI] == -1)
                    {
                        if
                        (
                           (
                              selectionMethod == surfaceZonesInfo::INSIDE
                              && volType[pointI] == volumeType::INSIDE
                           )
                          ||
                           (
                              selectionMethod == surfaceZonesInfo::OUTSIDE
                              && volType[pointI] == volumeType::OUTSIDE
                           )
                         )
                         {
                            insideSurfaces[pointI] = surfI;
                         }
                    }
                }
            }
        }
    }
}


void Foam::refinementSurfaces::findAnyIntersection
(
    const labelList& surfacesToTest,
    const pointField& start,
    const pointField& end,

    labelList& hitSurface,
    List<pointIndexHit>& hitInfo
) const
{
    labelList geometries(UIndirectList<label>(surfaces_, surfacesToTest));

    searchableSurfacesQueries::findAnyIntersection
    (
        allGeometry_,
        geometries,
        start,
        end,
        hitSurface,
        hitInfo
    );
}



//// Find intersection with max of edge. Return -1 or the surface
//// with the highest maxLevel above currentLevel
void Foam::refinementSurfaces::findHighestIntersection
(
    const pointField& start,
    const pointField& end,
    const labelList currentLevel,   // current cell refinement level

    labelList& maxLevel,
    labelList& hitSurface,
    List<pointIndexHit>& maxHit
) const
{
    hitSurface.setSize(start.size());
    maxLevel.setSize(start.size());
    maxHit.setSize(start.size());

    // surface with highest maxlevel
    labelList maxSurface(start.size(), -1);
    // maxLevel of maxSurface
    maxLevel = currentLevel;

    List<List<pointIndexHit>> hitInfo;


    forAll(surfaces_, surfI)
    {
        allGeometry_[surfaces_[surfI]].findLineAll(start, end, hitInfo);

        label n = 0;
        forAll(hitInfo, pointI)
        {
            n += hitInfo[pointI].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointI)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointI];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointI;
                n++;
            }
        }

        labelList hitRegion(n);
        allGeometry_[surfaces_[surfI]].getRegion
        (
            surfInfo,
            hitRegion
        );

        n = 0;

        forAll(hitInfo, pointI)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointI];

            // Extract those hits that are on higher levelled surfaces.
            // Note: should move extraction of region, normal outside of loop
            // below for if getRegion/getNormal have high overhead.

            forAll(pHits, pHitI)
            {
                label region = globalRegion(surfI, hitRegion[n]);
                n++;

                if (maxLevel_[region] > maxLevel[pointI])
                {
                    // Append to pointI info
                    maxSurface[pointI] = surfI;
                    maxLevel[pointI] = maxLevel_[region];
                    hitSurface[pointI] = surfI;
                    maxHit[pointI] = pHits[pHitI];
                }
            }
        }
    }

    forAll(hitInfo, pointI)
    {
       if (maxSurface[pointI] == -1)
       {
        // maxLevel unchanged. No interesting surface hit.
          maxHit[pointI].setMiss();
       }
    }
}


//// Find all intersected surface faces (we're kind of assuming an STL here :-)
//// and return them organised by surface
void Foam::refinementSurfaces::findAllIntersectedFaces
(
    const pointField& start,
    const pointField& end,

    labelList& proxSurfaces,
    List<List<labelList>>& allHits
) const
{
    //for every segment, we create a Map from the surface index to all
    //the faces hit in that surface.  NB Most of the time this structure
    //will be empty or nearly empty, but we do need all the information
    //it preserves.
    proxSurfaces.setSize(surfaces_.size());
    label nProxSurfaces = 0;

    forAll(surfaces_, surfI)
    {
        bool proxSurf = false;

        label nRegions
        (
            singleRegion_[surfI] ?
            1 : allGeometry_[surfaces_[surfI]].regions().size()
        );
        for (int i = 0; i < nRegions; i++)
        {
            label globalRegionI = globalRegion(surfI, i);
            if (proxLevelIncr()[globalRegionI] >= 0)
            {
                proxSurf = true;
                break;
            }
        }

        const searchableSurface& s = allGeometry_[surfaces_[surfI]];
        if (isA<triSurfaceMesh>(s) && proxSurf)
        {
            proxSurfaces[nProxSurfaces] = surfI;
            nProxSurfaces++;
        }
    }
    proxSurfaces.resize(nProxSurfaces);

    allHits.setSize(start.size());
    forAll(allHits, i)
    {
        allHits[i].setSize(proxSurfaces.size(),labelList());
    }

    // Work arrays
    List<List<pointIndexHit>> hitInfo;
    forAll(proxSurfaces, i)
    {
        label surfI = proxSurfaces[i];

        allGeometry_[surfaces_[surfI]].findLineAll(start, end, hitInfo);

        forAll(hitInfo, pointI)
        {
            labelList &hitList = allHits[pointI][i];
            hitList.setSize(hitInfo[pointI].size());
                //...and populate with face numbers
            forAll(hitInfo[pointI],pHitI)
            {
                hitList[pHitI] = hitInfo[pointI][pHitI].index();
            }
        }
    }
}


bool Foam::refinementSurfaces::baffleRemoval() const
{

    forAll(surfaces_, surfI)
    {
        label nRegions
        (
            singleRegion_[surfI] ?
            1 : allGeometry_[surfaces_[surfI]].regions().size()
        );
        for (int i = 0; i < nRegions; i++)
        {
            if (minBaffleAngle(surfI,i) != GREAT)
            {
                return true;
            }
        }
    }

    return false;
}


Foam::labelList Foam::refinementSurfaces::singleCellClosureSurfaces() const
{
    DynamicList<label> singleCellSurfaces(surfaces_.size());

    forAll(surfaces_, surfI)
    {
        label nRegions
        (
            singleRegion_[surfI] ?
            1 : allGeometry_[surfaces_[surfI]].regions().size()
        );
        for (int i = 0; i < nRegions; i++)
        {
            label globalRegionI = globalRegion(surfI, i);
            if (singleCellClosure()[globalRegionI])
            {
                singleCellSurfaces.append(surfI);
                break;
            }
        }
    }

    return labelList(singleCellSurfaces, true);
}

Foam::labelList Foam::refinementSurfaces::thinGapSurfaces() const
{
    DynamicList<label> thinSurfaces(surfaces_.size());

    forAll(surfaces_, surfI)
    {
        label nRegions
        (
            singleRegion_[surfI] ?
            1 : allGeometry_[surfaces_[surfI]].regions().size()
        );
        for (int i = 0; i < nRegions; i++)
        {
            label globalRegionI = globalRegion(surfI, i);
            if (thinGap()[globalRegionI])
            {
                thinSurfaces.append(surfI);
                break;
            }
        }
    }

    return labelList(thinSurfaces, true);
}

Foam::labelList Foam::refinementSurfaces::cornerCellSurfaces() const
{
    DynamicList<label> cornerSurfaces(surfaces_.size());

    forAll(surfaces_, surfI)
    {
        label nRegions
        (
            singleRegion_[surfI] ?
            1 : allGeometry_[surfaces_[surfI]].regions().size()
        );
        for (int i = 0; i < nRegions; i++)
        {
            label globalRegionI = globalRegion(surfI, i);
            if (cornerCellRegions()[globalRegionI])
            {
                cornerSurfaces.append(surfI);
                break;
            }
        }
    }

    return labelList(cornerSurfaces, true);
}


Foam::labelList Foam::refinementSurfaces::wrapLevelSurfaces() const
{
    DynamicList<label> wrapLevelSurfaces(surfaces_.size());

    const labelList unNamedSurfaces =
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfZones())
    );

    forAll(unNamedSurfaces, i)
    {
        label surfI = unNamedSurfaces[i];
        if (wrapLevel_[surfI] > -1)
        {
            wrapLevelSurfaces.append(surfI);
        }
    }

    return labelList(wrapLevelSurfaces, true);
}


Foam::labelList Foam::refinementSurfaces::boundSurfaces() const
{
    DynamicList<label> bSurfaces(surfaces_.size());

    const labelList unNamedSurfaces =
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfZones())
    );

    forAll(unNamedSurfaces, i)
    {
        label surfI = unNamedSurfaces[i];
        if (boundingSurface_[surfI])
        {
            bSurfaces.append(surfI);
        }
    }

    return labelList(bSurfaces, true);
}


Foam::labelList Foam::refinementSurfaces::moveCentroidSurfaces() const
{
    DynamicList<label> centroidSurfaces(surfaces_.size());

    forAll(surfaces_, surfI)
    {
        if (moveCentroids_[surfI] > SMALL)
        {
            centroidSurfaces.append(surfI);
        }
    }

    return labelList(centroidSurfaces, true);
}


Foam::labelList Foam::refinementSurfaces::proximitySurfaces() const
{
    DynamicList<label> proxSurfaces(surfaces_.size());

    forAll(surfaces_, surfI)
    {
        bool proxSurf = false;

        label nRegions
        (
            singleRegion_[surfI] ?
            1 : allGeometry_[surfaces_[surfI]].regions().size()
        );
        for (int i = 0; i < nRegions; i++)
        {
            label globalRegionI = globalRegion(surfI, i);
            if (proxLevelIncr()[globalRegionI] >= 0)
            {
                proxSurf = true;
                break;
            }
        }

        if (proxSurf)
        {
            proxSurfaces.append(surfI);
        }
    }

    return labelList(proxSurfaces, true);
}


Foam::labelList Foam::refinementSurfaces::curvatureSurfaces() const
{
    DynamicList<label> curveSurfaces(surfaces_.size());

    forAll(surfaces_, surfI)
    {
        bool curveSurf = false;
        label nRegions
        (
            singleRegion_[surfI] ?
            1 : allGeometry_[surfaces_[surfI]].regions().size()
        );
        for (int i = 0; i < nRegions; i++)
        {
            label globalRegionI = globalRegion(surfI, i);
            if (curvature()[globalRegionI] > SMALL)
            {
                curveSurf = true;
                break;
            }
        }

        if (curveSurf)
        {
            curveSurfaces.append(surfI);
        }
    }

    return labelList(curveSurfaces, true);
}


void Foam::refinementSurfaces::findAllProximityIntersections
(
    const pointField& start,
    const pointField& end,

    List<List<pointIndexHit>>& allHits,
    List<List<vector>>& surfaceNormal,
    labelListList& surfaceRegion
) const
{
    allHits.setSize(start.size());
    surfaceNormal.setSize(start.size());
    surfaceRegion.setSize(start.size());

    labelList proxSurfaces(proximitySurfaces());

    if (proxSurfaces.empty())
    {
        return;
    }

    forAll(proxSurfaces, i)
    {
        label surfI = proxSurfaces[i];

        // Work arrays
        List<List<pointIndexHit>> hitInfo;
        labelList pRegions;
        vectorField pNormals;

        allGeometry_[surfaces_[surfI]].findLineAll(start, end, hitInfo);

        // Repack hits for surface into flat list
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // To avoid overhead of calling getRegion for every point

        label n = 0;
        forAll(hitInfo, pointI)
        {
            n += hitInfo[pointI].size();
        }

        List<pointIndexHit> surfInfo(n);
        labelList pointMap(n);
        n = 0;

        forAll(hitInfo, pointI)
        {
            const List<pointIndexHit>& pHits = hitInfo[pointI];

            forAll(pHits, i)
            {
                surfInfo[n] = pHits[i];
                pointMap[n] = pointI;
                n++;
            }
        }

        labelList surfRegion(n);
        vectorField surfNormal(n);
        allGeometry_[surfaces_[surfI]].getRegion(surfInfo, surfRegion);
        allGeometry_[surfaces_[surfI]].getNormal(surfInfo, surfNormal);

        // Extract back into pointwise
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        forAll(surfRegion, i)
        {
            label region = globalRegion(surfI, surfRegion[i]);
            label pointI = pointMap[i];

            // Append to pointI info
            label sz = surfaceNormal[pointI].size();
            surfaceNormal[pointI].setSize(sz+1);
            surfaceNormal[pointI][sz] = surfNormal[i];
            surfaceRegion[pointI].setSize(sz+1);
            surfaceRegion[pointI][sz] = region;
            allHits[pointI].setSize(sz+1);
            allHits[pointI][sz] = surfInfo[i];
        }

        surfInfo.clear();
    }
}

// ************************************************************************* //
