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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2010-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshRefinement/meshRefinement.H"
#include "volMesh/volMesh.H"
#include "fields/volFields/volFields.H"
#include "surfaceMesh/surfaceMesh.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/Time/Time.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "refinementFeatures/refinementFeatures.H"
#include "shellSurfaces/shellSurfaces.H"
#include "decompositionMethod/decompositionMethod.H"
#include "regionSplit/regionSplit.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "meshes/primitiveMesh/primitivePatch/indirectPrimitivePatch.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "regionSplit/localPointRegion.H"
#include "meshes/pointMesh/pointMesh.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "fields/pointPatchFields/derived/slip/slipPointPatchFields.H"
#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchFields.H"
#include "fields/pointPatchFields/basic/calculated/calculatedPointPatchFields.H"
#include "fields/pointPatchFields/constraint/cyclicSlip/cyclicSlipPointPatchFields.H"
#include "meshes/pointMesh/pointPatches/constraint/processor/processorPointPatch.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "meshTools/meshTools.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "geomDecomp/geomDecomp.H"
#include "primitives/random/Random/Random.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "meshes/treeBoundBox/treeBoundBox.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "fvMeshTools/fvMeshTools.H"
#include "motionSmoother/motionSmoother.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"
#include "sets/topoSets/cellSet.H"
#include "patchDist/wallPoint/wallPoint.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "polyTopoChange/polyTopoChange/removeObject/polyRemovePoint.H"
#include "polyTopoChange/hexRef8/hexRef8DataList.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "edgeClassification/edgeClassification.H"

// Leak path
#include "sampledSet/shortestPath/shortestPathSet.H"
#include "meshSearch/meshSearch.H"

#ifdef FOAM_USE_TBB
  #include <tbb/parallel_for.h>
  #include <tbb/parallel_reduce.h>
  #include "include/TBBTimer.H"
#endif

#include "renumberUtils.H"
#include "renumberMethod/renumberMethod.H"
#include "CuthillMcKeeRenumber/CuthillMcKeeRenumber.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshRefinement, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshRefinement::IOdebugType,
        5
    >::names[] =
    {
        "mesh",
        "intersections",
        "featureSeeds",
        "attraction",
        "layerInfo"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshRefinement::IOoutputType,
        1
    >::names[] =
    {
        "layerInfo"
    };

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshRefinement::IOwriteType,
        3
    >::names[] =
    {
        "mesh",
        "noRefinement",
        "scalarLevels"
    };

}

const Foam::NamedEnum<Foam::meshRefinement::IOdebugType, 5>
Foam::meshRefinement::IOdebugTypeNames;

const Foam::NamedEnum<Foam::meshRefinement::IOoutputType, 1>
Foam::meshRefinement::IOoutputTypeNames;

const Foam::NamedEnum<Foam::meshRefinement::IOwriteType, 3>
Foam::meshRefinement::IOwriteTypeNames;


Foam::meshRefinement::writeType Foam::meshRefinement::writeLevel_;

Foam::meshRefinement::outputType Foam::meshRefinement::outputLevel_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//- Calculate maximum cell level
label Foam::meshRefinement::maxCellLevel() const
{
    label maxLevel = 0;

    //find maximum feature level
    const labelListList& featureLevels = features()().levels();
    forAll(featureLevels, featI)
    {
        maxLevel = max(maxLevel, max(featureLevels[featI]));
    }

    //find maximum surface level
    labelList surfLevels = surfaces_.maxLevel();
    forAll(surfaces_.proxLevelIncr(), surfi)
    {
        label proxLevel = surfaces_.proxLevelIncr()[surfi];
        if (proxLevel > 0)
        {
            surfLevels[surfi] += proxLevel;
        }
    }

    const labelList& surfMaxGapLevel = surfaces_.maxGapLevel();

    maxLevel = max(maxLevel, max(surfLevels));
    maxLevel = max(maxLevel, max(surfMaxGapLevel));

    //find maximum shell level
    const labelList& shellMaxGapLevel = surfaces_.maxGapLevel();

    maxLevel = max(maxLevel, shells_.maxIsoLevel());
    maxLevel = max(maxLevel, shells_.maxAnisoLevel(0));
    maxLevel = max(maxLevel, shells_.maxAnisoLevel(1));
    maxLevel = max(maxLevel, shells_.maxAnisoLevel(2));

    maxLevel = max(maxLevel, max(shellMaxGapLevel));

    return maxLevel;
}


//- Estimate of crack tolerance
void Foam::meshRefinement::setCrackTol()
{
    const scalar edge0Len = meshCutter_.level0EdgeLength();
    labelList maxSurfLevels = surfaces_.maxLevel();
    forAll(surfaces_.proxLevelIncr(), surfi)
    {
        label proxLevel = surfaces_.proxLevelIncr()[surfi];
        if (proxLevel > 0)
        {
            maxSurfLevels[surfi] += proxLevel;
        }
    }

    const label maxLevel = max(max(maxSurfLevels), shells_.maxIsoLevel());

    crackTol_ = crackDetection_.second() * edge0Len/(1<<maxLevel);

    if (checkForCracks())
    {
        Info<<"Using crack detection to close cracks to tolerance: "
            <<crackTol_<<endl;
    }
}

void Foam::meshRefinement::calcNeighbourData
(
    const pointField& cellCentres,
    labelList& neiLevel,
    pointField& neiCc,
    const bool threaded /*= false*/
)  const
{
//#ifdef FOAM_USE_TBB
//    Timer t("meshRefinement::calcNeighbourData");
//#endif

    const labelList& cellLevel = meshCutter_.cellLevel();

    label nBoundaryFaces = mesh_.nFaces() - mesh_.nInternalFaces();

    if (neiLevel.size() != nBoundaryFaces || neiCc.size() != nBoundaryFaces)
    {
        FatalErrorInFunction
            << nBoundaryFaces << " neiLevel:" << neiLevel.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelHashSet addedPatchIDSet(meshedPatches());

    if (threaded)
    {
#ifndef FOAM_USE_TBB
        WarningInFunction
            << "FOAM not linked against TBB! Cannot run multithreaded."
            << endl;
#else
        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, patches.size()),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t patchi = r.begin(); patchi < r.end(); ++patchi)
                {
                    const polyPatch& pp = patches[patchi];

                    const labelUList& faceCells = pp.faceCells();
                    const vectorField::subField faceCentres = pp.faceCentres();
                    const vectorField::subField faceAreas = pp.faceAreas();

                    label bFacei = pp.start()-mesh_.nInternalFaces();
                    if (pp.coupled())
                    {
                        tbb::parallel_for
                        (
                            tbb::blocked_range<size_t>(0, faceCells.size()),
                            [&](const tbb::blocked_range<size_t>& rp)
                            {
                                for (size_t i = rp.begin(); i < rp.end(); ++i)
                                {
                                    neiLevel[bFacei+i] = cellLevel[faceCells[i]];
                                    neiCc[bFacei+i] = cellCentres[faceCells[i]];
                                }
                            }
                        );
                        //forAll(faceCells, i)
                        //{
                        //    neiLevel[bFacei+i] = cellLevel[faceCells[i]];
                        //    neiCc[bFacei+i] = cellCentres[faceCells[i]];
                        //}
                    }
                    else if (addedPatchIDSet.found(patchi))
                    {
                        // Face was introduced from cell-cell intersection. Try to
                        // reconstruct other side cell(centre). Three possibilities:
                        // - cells same size.
                        // - preserved cell smaller. Not handled.
                        // - preserved cell larger.
                        tbb::parallel_for
                        (
                            tbb::blocked_range<size_t>(0, faceCells.size()),
                            [&](const tbb::blocked_range<size_t>& rp)
                            {
                                for (size_t i = rp.begin(); i < rp.end(); ++i)
                                {
                                    // Extrapolate the face centre.
                                    vector fn = faceAreas[i];
                                    fn /= mag(fn)+VSMALL;

                                    label own = faceCells[i];
                                    label ownLevel = cellLevel[own];
                                    label faceLevel = meshCutter_.faceLevel(pp.start()+i);
                                    if (faceLevel < 0)
                                    {
                                        // Due to e.g. face merging no longer a consistent
                                        // refinementlevel of face. Assume same as cell
                                        faceLevel = ownLevel;
                                    }

                                    // Normal distance from face centre to cell centre
                                    scalar d = ((faceCentres[i] - cellCentres[own]) & fn);
                                    if (faceLevel > ownLevel)
                                    {
                                        // Other cell more refined. Adjust normal distance
                                        d *= 0.5;
                                    }
                                    neiLevel[bFacei+i] = faceLevel;
                                    // Calculate other cell centre by extrapolation
                                    neiCc[bFacei+i] = faceCentres[i] + d*fn;
                                }
                            }
                        );
                        //forAll(faceCells, i)
                        //{
                        //    // Extrapolate the face centre.
                        //    vector fn = faceAreas[i];
                        //    fn /= mag(fn)+VSMALL;

                        //    label own = faceCells[i];
                        //    label ownLevel = cellLevel[own];
                        //    label faceLevel = meshCutter_.faceLevel(pp.start()+i);
                        //    if (faceLevel < 0)
                        //    {
                        //        // Due to e.g. face merging no longer a consistent
                        //        // refinementlevel of face. Assume same as cell
                        //        faceLevel = ownLevel;
                        //    }

                        //    // Normal distance from face centre to cell centre
                        //    scalar d = ((faceCentres[i] - cellCentres[own]) & fn);
                        //    if (faceLevel > ownLevel)
                        //    {
                        //        // Other cell more refined. Adjust normal distance
                        //        d *= 0.5;
                        //    }
                        //    neiLevel[bFacei+i] = faceLevel;
                        //    // Calculate other cell centre by extrapolation
                        //    neiCc[bFacei+i] = faceCentres[i] + d*fn;
                        //}
                    }
                    else
                    {
                        tbb::parallel_for
                        (
                            tbb::blocked_range<size_t>(0, faceCells.size()),
                            [&](const tbb::blocked_range<size_t>& rp)
                            {
                                for (size_t i = rp.begin(); i < rp.end(); ++i)
                                {
                                    neiLevel[bFacei+i] = cellLevel[faceCells[i]];
                                    neiCc[bFacei+i] = faceCentres[i];
                                }
                            }
                        );
                        //forAll(faceCells, i)
                        //{
                        //    neiLevel[bFacei+i] = cellLevel[faceCells[i]];
                        //    neiCc[bFacei+i] = faceCentres[i];
                        //}
                    }
                } //for patchi in range
            } //parallel_for func
        );
#endif
    }
    else
    {
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];

            const labelUList& faceCells = pp.faceCells();
            const vectorField::subField faceCentres = pp.faceCentres();
            const vectorField::subField faceAreas = pp.faceAreas();

            label bFacei = pp.start()-mesh_.nInternalFaces();

            if (pp.coupled())
            {
                forAll(faceCells, i)
                {
                    neiLevel[bFacei] = cellLevel[faceCells[i]];
                    neiCc[bFacei] = cellCentres[faceCells[i]];
                    bFacei++;
                }
            }
            else if (addedPatchIDSet.found(patchi))
            {
                // Face was introduced from cell-cell intersection. Try to
                // reconstruct other side cell(centre). Three possibilities:
                // - cells same size.
                // - preserved cell smaller. Not handled.
                // - preserved cell larger.
                forAll(faceCells, i)
                {
                    // Extrapolate the face centre.
                    vector fn = faceAreas[i];
                    fn /= mag(fn)+VSMALL;

                    label own = faceCells[i];
                    label ownLevel = cellLevel[own];
                    label faceLevel = meshCutter_.faceLevel(pp.start()+i);
                    if (faceLevel < 0)
                    {
                        // Due to e.g. face merging no longer a consistent
                        // refinementlevel of face. Assume same as cell
                        faceLevel = ownLevel;
                    }

                    // Normal distance from face centre to cell centre
                    scalar d = ((faceCentres[i] - cellCentres[own]) & fn);
                    if (faceLevel > ownLevel)
                    {
                        // Other cell more refined. Adjust normal distance
                        d *= 0.5;
                    }
                    neiLevel[bFacei] = faceLevel;
                    // Calculate other cell centre by extrapolation
                    neiCc[bFacei] = faceCentres[i] + d*fn;
                    bFacei++;
                }
            }
            else
            {
                forAll(faceCells, i)
                {
                    neiLevel[bFacei] = cellLevel[faceCells[i]];
                    neiCc[bFacei] = faceCentres[i];
                    bFacei++;
                }
            }
        }
    }

    // Swap coupled boundaries. Apply separation to cc since is coordinate.
    syncTools::swapBoundaryFacePositions(mesh_, neiCc, threaded);
    syncTools::swapBoundaryFaceList(mesh_, neiLevel, threaded);
} //calcNeighbourData


void Foam::meshRefinement::calcPerturbedNeighbourData
(
    const pointField& cellCentres,
    const pointField& perturbCellCentres,
    labelList& neiLevel,
    pointField& neiCc
)  const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    label nBoundaryFaces = mesh_.nFaces() - mesh_.nInternalFaces();

    if (neiLevel.size() != nBoundaryFaces || neiCc.size() != nBoundaryFaces)
    {
        FatalErrorInFunction
            << nBoundaryFaces << " neiLevel:" << neiLevel.size()
            << abort(FatalError);
    }

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    labelHashSet addedPatchIDSet(meshedPatches());

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        const labelUList& faceCells = pp.faceCells();
        const vectorField::subField faceCentres = pp.faceCentres();
        const vectorField::subField faceAreas = pp.faceAreas();

        label bFacei = pp.start()-mesh_.nInternalFaces();

        if (pp.coupled())
        {
            forAll(faceCells, i)
            {
                neiLevel[bFacei] = cellLevel[faceCells[i]];
                neiCc[bFacei] = perturbCellCentres[faceCells[i]];
                bFacei++;
            }
        }
        else if (addedPatchIDSet.found(patchi))
        {
            // Face was introduced from cell-cell intersection. Try to
            // reconstruct other side cell(centre). Three possibilities:
            // - cells same size.
            // - preserved cell smaller. Not handled.
            // - preserved cell larger.
            forAll(faceCells, i)
            {
                // Extrapolate the face centre.
                vector fn = faceAreas[i];
                fn /= mag(fn)+VSMALL;

                label own = faceCells[i];
                label ownLevel = cellLevel[own];
                label faceLevel = meshCutter_.faceLevel(pp.start()+i);
                if (faceLevel < 0)
                {
                    // Due to e.g. face merging no longer a consistent
                    // refinementlevel of face. Assume same as cell
                    faceLevel = ownLevel;
                }

                // Normal distance from face centre to cell centre
                scalar d = ((faceCentres[i] - cellCentres[own]) & fn);
                if (faceLevel > ownLevel)
                {
                    // Other cell more refined. Adjust normal distance
                    d *= 0.5;
                }
                vector pertVec = perturbCellCentres[own] - cellCentres[own];

                neiLevel[bFacei] = faceLevel;
                // Calculate other cell centre by extrapolation
                neiCc[bFacei] = (faceCentres[i] + d*fn) + pertVec;
                bFacei++;
            }
        }
        else
        {
            forAll(faceCells, i)
            {
                neiLevel[bFacei] = cellLevel[faceCells[i]];
                neiCc[bFacei] = faceCentres[i];
                bFacei++;
            }
        }
    }

    // Swap coupled boundaries. Apply separation to cc since is coordinate.
    syncTools::swapBoundaryFacePositions(mesh_, neiCc);
    syncTools::swapBoundaryFaceList(mesh_, neiLevel);
}


void Foam::meshRefinement::calcCellCellRays
(
    const pointField& cellCentres,
    const pointField& neiCc,
    const labelList& neiLevel,
    const labelList& testFaces,
    pointField& start,
    pointField& end,
    labelList& minLevel,
    const bool threaded /*= false*/
) const
{
//#ifdef FOAM_USE_TBB
//    Timer timer("meshRefinement::calcCellCellRays");
//#endif

    const labelList& cellLevel = meshCutter_.cellLevel();

    start.setSize(testFaces.size());
    end.setSize(testFaces.size());
    minLevel.setSize(testFaces.size());

    if (threaded && testFaces.size() > 2048)
    {
#ifndef FOAM_USE_TBB
        WarningInFunction
            << "FOAM not linked against TBB! Cannot run multithreaded."
            << endl;
#else
        const auto grainSize =
            std::max<size_t>
            (
                testFaces.size() / tbb::this_task_arena::max_concurrency(),
                1024
            );

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, testFaces.size(), grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    label facei = testFaces[i];
                    label own = mesh_.faceOwner()[facei];

                    if (mesh_.isInternalFace(facei))
                    {
                        label nei = mesh_.faceNeighbour()[facei];

                        start[i] = cellCentres[own];
                        end[i] = cellCentres[nei];
                        minLevel[i] = min(cellLevel[own], cellLevel[nei]);
                    }
                    else
                    {
                        label bFacei = facei - mesh_.nInternalFaces();

                        start[i] = cellCentres[own];
                        end[i] = neiCc[bFacei];
                        minLevel[i] = min(cellLevel[own], neiLevel[bFacei]);
                    }
                }
            },
            tbb::simple_partitioner()
        );
#endif
    }
    else
    {
        forAll(testFaces, i)
        {
            label facei = testFaces[i];
            label own = mesh_.faceOwner()[facei];

            if (mesh_.isInternalFace(facei))
            {
                label nei = mesh_.faceNeighbour()[facei];

                start[i] = cellCentres[own];
                end[i] = cellCentres[nei];
                minLevel[i] = min(cellLevel[own], cellLevel[nei]);
            }
            else
            {
                label bFacei = facei - mesh_.nInternalFaces();

                start[i] = cellCentres[own];
                end[i] = neiCc[bFacei];
                minLevel[i] = min(cellLevel[own], neiLevel[bFacei]);
            }
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(ROOTSMALL*(end-start));
        start -= smallVec;
        end += smallVec;
    }
}


void Foam::meshRefinement::updateIntersections
(
    const labelList& changedFaces,
    const bool threaded /*= false*/
)
{
//#ifdef FOAM_USE_TBB
//    Timer timer("meshRefinement::updateIntersections");
//#endif

    // Stats on edges to test. Count proc faces only once.
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_, threaded));

    {
//#ifdef FOAM_USE_TBB
//        Timer timer("meshRefinement::updateIntersections::edgesToRetest");
//#endif

        label nMasterFaces = 0;

        if (threaded)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    mesh_.nFaces() / tbb::this_task_arena::max_concurrency(),
                    1024
                );
            const tbb::blocked_range<size_t> range(0, mesh_.nFaces(), grainSize);

            nMasterFaces =
                tbb::parallel_reduce
                (
                    range,
                    /*identity*/0,
                    [&](const tbb::blocked_range<size_t>& r, const label& init)
                    {
                        label res = init;
                        for (size_t facei = r.begin(); facei < r.end(); ++facei)
                        {
                            if (isMasterFace.get(facei) == 1)
                            {
                                res++;
                            }
                        }
                        return res;
                    },
                    /*reduction*/[](const label& x, const label& y)
                    {
                        return x+y;
                    },
                    tbb::simple_partitioner()
                );
#endif
        }
        else
        {
            forAll(isMasterFace, facei)
            {
                if (isMasterFace.get(facei) == 1)
                {
                    nMasterFaces++;
                }
            }
        }
        reduce(nMasterFaces, sumOp<label>());

        label nChangedFaces = 0;
        if (threaded && changedFaces.size() > 2048)
        {
#ifndef FOAM_USE_TBB
            WarningInFunction
                << "FOAM not linked against TBB! Cannot run multithreaded."
                << endl;
#else
            const auto grainSize =
                std::max<size_t>
                (
                    changedFaces.size() / tbb::this_task_arena::max_concurrency(),
                    1024
                );
            const tbb::blocked_range<size_t> range(0, changedFaces.size(), grainSize);

            nChangedFaces =
                tbb::parallel_reduce
                (
                    range,
                    /*identity*/0,
                    [&](const tbb::blocked_range<size_t>& r, const label& init)
                    {
                        label res = init;
                        for (size_t facei = r.begin(); facei < r.end(); ++facei)
                        {
                            if (isMasterFace.get(changedFaces[facei]) == 1)
                            {
                                res++;
                            }
                        }
                        return res;
                    },
                    /*reduction*/[](const label& x, const label& y)
                    {
                        return x+y;
                    },
                    tbb::simple_partitioner()
                );
#endif
        }
        else
        {
            forAll(changedFaces, i)
            {
                if (isMasterFace.get(changedFaces[i]) == 1)
                {
                    nChangedFaces++;
                }
            }
        }
        reduce(nChangedFaces, sumOp<label>());

        Info<< "Edge intersection testing:" << nl
            << "    Number of edges             : " << nMasterFaces << nl
            << "    Number of edges to retest   : " << nChangedFaces
            << endl;
    }


    // Get boundary face centre and level. Coupled aware.
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc, threaded);

    // Collect segments we want to test for
    pointField start(changedFaces.size());
    pointField end(changedFaces.size());
    {
        labelList minLevel;
        calcCellCellRays
        (
            mesh_.cellCentres(),
            neiCc,
            neiLevel,
            changedFaces,
            start,
            end,
            minLevel,
            threaded
        );
    }

    // Do tests in one go
    labelList surfaceHit;
    {
//#ifdef FOAM_USE_TBB
//        Timer timer("refinementSurfaces::findHigherIntersection");
//#endif
        labelList surfaceLevel;

        surfaces_.findHigherIntersection
        (
            shells_,
            start,
            end,
            labelList(start.size(), -1),    // accept any intersection
            surfaceHit,
            surfaceLevel,
            meshCutter_.level0EdgeLength(),
            /*checkCurvature*/false,
            threaded
         );
    }

    {
//#ifdef FOAM_USE_TBB
//    Timer t("meshRefinement::Keep just surface hit");
//#endif
    // Keep just surface hit
    if (threaded && surfaceHit.size() > 2048)
    {
#ifndef FOAM_USE_TBB
        WarningInFunction
            << "FOAM not linked against TBB! Cannot run multithreaded."
            << endl;
#else
        const auto grainSize =
            std::max<size_t>
            (
                surfaceHit.size() / tbb::this_task_arena::max_concurrency(),
                1024
            );

        tbb::parallel_for
        (
            tbb::blocked_range<size_t>(0, surfaceHit.size(), grainSize),
            [&](const tbb::blocked_range<size_t>& r)
            {
                for (size_t i = r.begin(); i < r.end(); ++i)
                {
                    surfaceIndex_[changedFaces[i]] = surfaceHit[i];
                }
            },
            tbb::simple_partitioner()
        );
#endif
    }
    else
    {
        forAll(surfaceHit, i)
        {
            surfaceIndex_[changedFaces[i]] = surfaceHit[i];
        }
    }
    } //Timer

    // Make sure both sides have same information. This should be
    // case in general since same vectors but just to make sure.
    syncTools::syncFaceList(mesh_, surfaceIndex_, maxEqOp<label>(), threaded);

    label nHits = countHits(threaded);
    label nTotHits = returnReduce(nHits, sumOp<label>());

    Info<< "    Number of intersected edges : " << nTotHits << endl;

    // Set files to same time as mesh
    setInstance(mesh_.facesInstance());
} // updateIntersections


// Block calculation (lower memory) to find intersections of edges
// (given by their two endpoints) with surfaces. Returns first intersection
//if there are more than one.
void Foam::meshRefinement::updateIntersectionsByBlock
(
    const labelList& changedFaces,
    const label numberBlocks
)
{
    const pointField& cellCentres = mesh_.cellCentres();

    // Stats on edges to test. Count proc faces only once.
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));

    {
        label nMasterFaces = 0;
        forAll(isMasterFace, faceI)
        {
            if (isMasterFace.get(faceI) == 1)
            {
                nMasterFaces++;
            }
        }
        reduce(nMasterFaces, sumOp<label>());

        label nChangedFaces = 0;
        forAll(changedFaces, i)
        {
            if (isMasterFace.get(changedFaces[i]) == 1)
            {
                nChangedFaces++;
            }
        }
        reduce(nChangedFaces, sumOp<label>());

        Info<< "Edge intersection testing:" << nl
            << "    Number of edges             : " << nMasterFaces << nl
            << "    Number of edges to retest   : " << nChangedFaces
            << endl;
    }

    // Get boundary face centre and level. Coupled aware.
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc);

    label blockSize = max(changedFaces.size() / numberBlocks, 1);

    label startIndex = 0;
    label endIndex = blockSize;

    for (int nb = 0; nb < numberBlocks + 1; nb++)
    {
        // Collect segments we want to test for
        label sz = endIndex - startIndex;
        pointField start(sz);
        pointField end(sz);

        for (int i = 0; i < sz; i++)
        {
            label faceI = changedFaces[startIndex+i];
            label own = mesh_.faceOwner()[faceI];

            start[i] = cellCentres[own];
            if (mesh_.isInternalFace(faceI))
            {
                end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
            }
            else
            {
                end[i] = neiCc[faceI-mesh_.nInternalFaces()];
            }
        }


        // Extend segments a bit
        {
            const vectorField smallVec(Foam::sqrt(SMALL)*(end-start));
            start -= smallVec;
            end += smallVec;
        }

        // Do tests in one go
        labelList surfaceHit;
        {
            labelList surfaceLevel;
            surfaces_.findHigherIntersection
            (
                shells_,
                start,
                end,
                labelList(start.size(), -1),    // accept any intersection
                surfaceHit,
                surfaceLevel,
                meshCutter_.level0EdgeLength()
             );
        }

        // Keep just surface hit
        forAll(surfaceHit, i)
        {
            surfaceIndex_[changedFaces[startIndex+i]] = surfaceHit[i];
        }


        startIndex = endIndex;
        endIndex = min(startIndex + blockSize, changedFaces.size());
    }

    label nHits = countHits();
    label nTotHits = returnReduce(nHits, sumOp<label>());

    Info<< "    Number of intersected edges : " << nTotHits << endl;

    // Set files to same time as mesh
    setInstance(mesh_.facesInstance());
}


void Foam::meshRefinement::testSyncPointList
(
    const string& msg,
    const polyMesh& mesh,
    const List<scalar>& fld
)
{
    if (fld.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << msg << endl
            << "fld size:" << fld.size() << " mesh points:" << mesh.nPoints()
            << abort(FatalError);
    }

    Pout<< "Checking field " << msg << endl;
    scalarField minFld(fld);
    syncTools::syncPointList
    (
        mesh,
        minFld,
        minEqOp<scalar>(),
        GREAT
    );
    scalarField maxFld(fld);
    syncTools::syncPointList
    (
        mesh,
        maxFld,
        maxEqOp<scalar>(),
        -GREAT
    );
    forAll(minFld, pointi)
    {
        const scalar& minVal = minFld[pointi];
        const scalar& maxVal = maxFld[pointi];
        if (mag(minVal-maxVal) > SMALL)
        {
            Pout<< msg << " at:" << mesh.points()[pointi] << nl
                << "    minFld:" << minVal << nl
                << "    maxFld:" << maxVal << nl
                << endl;
        }
    }
}


void Foam::meshRefinement::testSyncPointList
(
    const string& msg,
    const polyMesh& mesh,
    const List<point>& fld
)
{
    if (fld.size() != mesh.nPoints())
    {
        FatalErrorInFunction
            << msg << endl
            << "fld size:" << fld.size() << " mesh points:" << mesh.nPoints()
            << abort(FatalError);
    }

    Pout<< "Checking field " << msg << endl;
    pointField minFld(fld);
    syncTools::syncPointList
    (
        mesh,
        minFld,
        minMagSqrEqOp<point>(),
        point(GREAT, GREAT, GREAT)
    );
    pointField maxFld(fld);
    syncTools::syncPointList
    (
        mesh,
        maxFld,
        maxMagSqrEqOp<point>(),
        vector::zero
    );
    forAll(minFld, pointi)
    {
        const point& minVal = minFld[pointi];
        const point& maxVal = maxFld[pointi];
        if (mag(minVal-maxVal) > SMALL)
        {
            Pout<< msg << " at:" << mesh.points()[pointi] << nl
                << "    minFld:" << minVal << nl
                << "    maxFld:" << maxVal << nl
                << endl;
        }
    }
}


void Foam::meshRefinement::checkData()
{
    Pout<< "meshRefinement::checkData() : Checking refinement structure."
        << endl;
    meshCutter_.checkMesh();

    Pout<< "meshRefinement::checkData() : Checking refinement levels."
        << endl;
    meshCutter_.checkRefinementLevels(1, labelList(0));


    label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();

    Pout<< "meshRefinement::checkData() : Checking synchronization."
        << endl;

    // Check face centres
    {
        // Boundary face centres
        pointField::subList boundaryFc
        (
            mesh_.faceCentres(),
            nBnd,
            mesh_.nInternalFaces()
        );

        // Get neighbouring face centres
        pointField neiBoundaryFc(boundaryFc);
        syncTools::syncBoundaryFacePositions
        (
            mesh_,
            neiBoundaryFc,
            eqOp<point>()
        );

        // Compare
        testSyncBoundaryFaceList
        (
            mergeDistance_,
            "testing faceCentres : ",
            boundaryFc,
            neiBoundaryFc
        );
    }
    // Check meshRefinement
    {
        // Get boundary face centre and level. Coupled aware.
        labelList neiLevel(nBnd);
        pointField neiCc(nBnd);
        calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc);

        // Collect segments we want to test for
        pointField start(mesh_.nFaces());
        pointField end(mesh_.nFaces());

        forAll(start, facei)
        {
            start[facei] = mesh_.cellCentres()[mesh_.faceOwner()[facei]];

            if (mesh_.isInternalFace(facei))
            {
                end[facei] = mesh_.cellCentres()[mesh_.faceNeighbour()[facei]];
            }
            else
            {
                end[facei] = neiCc[facei-mesh_.nInternalFaces()];
            }
        }

        // Extend segments a bit
        {
            const vectorField smallVec(ROOTSMALL*(end-start));
            start -= smallVec;
            end += smallVec;
        }

        // Do tests in one go
        labelList surfaceHit;
        {
            labelList surfaceLevel;
            surfaces_.findHigherIntersection
            (
                shells_,
                start,
                end,
                labelList(start.size(), -1),    // accept any intersection
                surfaceHit,
                surfaceLevel,
                meshCutter_.level0EdgeLength()
            );
        }
        // Get the coupled hit
        labelList neiHit
        (
            SubList<label>
            (
                surfaceHit,
                nBnd,
                mesh_.nInternalFaces()
            )
        );
        syncTools::swapBoundaryFaceList(mesh_, neiHit);

        // Check
        forAll(surfaceHit, facei)
        {
            if (surfaceIndex_[facei] != surfaceHit[facei])
            {
                if (mesh_.isInternalFace(facei))
                {
                    WarningInFunction
                        << "Internal face:" << facei
                        << " fc:" << mesh_.faceCentres()[facei]
                        << " cached surfaceIndex_:" << surfaceIndex_[facei]
                        << " current:" << surfaceHit[facei]
                        << " ownCc:"
                        << mesh_.cellCentres()[mesh_.faceOwner()[facei]]
                        << " neiCc:"
                        << mesh_.cellCentres()[mesh_.faceNeighbour()[facei]]
                        << endl;
                }
                else if
                (
                    surfaceIndex_[facei]
                 != neiHit[facei-mesh_.nInternalFaces()]
                )
                {
                    WarningInFunction
                        << "Boundary face:" << facei
                        << " fc:" << mesh_.faceCentres()[facei]
                        << " cached surfaceIndex_:" << surfaceIndex_[facei]
                        << " current:" << surfaceHit[facei]
                        << " ownCc:"
                        << mesh_.cellCentres()[mesh_.faceOwner()[facei]]
                        << " neiCc:"
                        << neiCc[facei-mesh_.nInternalFaces()]
                        << " end:" << end[facei]
                        << " ownLevel:"
                        << meshCutter_.cellLevel()[mesh_.faceOwner()[facei]]
                        << " faceLevel:"
                        << meshCutter_.faceLevel(facei)
                        << endl;
                }
            }
        }
    }
    {
        labelList::subList boundarySurface
        (
            surfaceIndex_,
            mesh_.nFaces()-mesh_.nInternalFaces(),
            mesh_.nInternalFaces()
        );

        labelList neiBoundarySurface(boundarySurface);
        syncTools::swapBoundaryFaceList
        (
            mesh_,
            neiBoundarySurface
        );

        // Compare
        testSyncBoundaryFaceList
        (
            0,                              // tolerance
            "testing surfaceIndex() : ",
            boundarySurface,
            neiBoundarySurface
        );
    }


    // Find duplicate faces
    Pout<< "meshRefinement::checkData() : Counting duplicate faces."
        << endl;

    labelList duplicateFace
    (
        localPointRegion::findDuplicateFaces
        (
            mesh_,
            identity(mesh_.nFaces()-mesh_.nInternalFaces())
          + mesh_.nInternalFaces()
        )
    );

    // Count
    {
        label nDup = 0;

        forAll(duplicateFace, i)
        {
            if (duplicateFace[i] != -1)
            {
                nDup++;
            }
        }
        nDup /= 2;  // will have counted both faces of duplicate
        Pout<< "meshRefinement::checkData() : Found " << nDup
            << " duplicate pairs of faces." << endl;
    }
}


void Foam::meshRefinement::setInstance(const fileName& inst)
{
    meshCutter_.setInstance(inst);
    cellRegionID_.instance() = inst;
    surfaceIndex_.instance() = inst;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::doRemoveCells
(
    const labelList& cellsToRemove,
    const labelList& exposedFaces,
    const labelList& exposedPatchIDs,
    removeCells& cellRemover,
    const bool updateIntersections
)
{
    polyTopoChange meshMod(mesh_);

    // Arbitrary: put exposed faces into last patch.
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );

    // Change the mesh (no inflation)
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());
    setInstance(mesh_.facesInstance());

    // Update local mesh data
    cellRemover.updateMesh(map);

    // Update intersections. Recalculate intersections for exposed faces.
    labelList newExposedFaces = renumber
    (
        map().reverseFaceMap(),
        exposedFaces
    );

    //Pout<< "removeCells : updating intersections for "
    //    << newExposedFaces.size() << " newly exposed faces." << endl;

    updateMesh(map, newExposedFaces, updateIntersections);

    return map;
}


void Foam::meshRefinement::doSplitFaces
(
    const labelList& splitFaces,
    const labelPairList& splits,
    //const List<Pair<point>>& splitPoints,
    polyTopoChange& meshMod
) const
{
    forAll(splitFaces, i)
    {
        label facei = splitFaces[i];
        const face& f = mesh_.faces()[facei];

        // Split as start and end index in face
        const labelPair& split = splits[i];

        label nVerts = split[1]-split[0];
        if (nVerts < 0)
        {
            nVerts += f.size();
        }
        nVerts += 1;


        // Split into f0, f1
        face f0(nVerts);

        label fp = split[0];
        forAll(f0, i)
        {
            f0[i] = f[fp];
            fp = f.fcIndex(fp);
        }

        face f1(f.size()-f0.size()+2);
        fp = split[1];
        forAll(f1, i)
        {
            f1[i] = f[fp];
            fp = f.fcIndex(fp);
        }


        // Determine face properties
        label own = mesh_.faceOwner()[facei];
        label nei = -1;
        label patchi = -1;
        if (facei >= mesh_.nInternalFaces())
        {
            patchi = mesh_.boundaryMesh().whichPatch(facei);
        }
        else
        {
            nei = mesh_.faceNeighbour()[facei];
        }

        label zonei = mesh_.faceZones().whichZone(facei);
        bool zoneFlip = false;
        if (zonei != -1)
        {
            const faceZone& fz = mesh_.faceZones()[zonei];
            zoneFlip = fz.flipMap()[fz.whichFace(facei)];
        }


        if (debug)
        {
            Pout<< "face:" << facei << " verts:" << f
                << " split into f0:" << f0
                << " f1:" << f1 << endl;
        }

        // Change/add faces
        meshMod.modifyFace
        (
            f0,                         // modified face
            facei,                      // label of face
            own,                        // owner
            nei,                        // neighbour
            false,                      // face flip
            patchi,                     // patch for face
            zonei,                      // zone for face
            zoneFlip                    // face flip in zone
        );

        meshMod.addFace
        (
            f1,                         // modified face
            own,                        // owner
            nei,                        // neighbour
            -1,                         // master point
            -1,                         // master edge
            facei,                      // master face
            false,                      // face flip
            patchi,                     // patch for face
            zonei,                      // zone for face
            zoneFlip                    // face flip in zone
        );


        //// Move points
        //meshMod.modifyPoint
        //(
        //    f[split[0]],
        //    splitPoints[i][0],
        //    -1,
        //    true
        //);
        //meshMod.modifyPoint
        //(
        //    f[split[1]],
        //    splitPoints[i][1],
        //    -1,
        //    true
        //);
    }
}


Foam::label Foam::meshRefinement::splitFacesUndo
(
    const labelList& splitFaces,
    const labelPairList& splits,
    const dictionary& motionDict,

    labelList& duplicateFace,
    List<labelPair>& baffles
)
{
    label nSplit = returnReduce(splitFaces.size(), sumOp<label>());
    Info<< nl
        << "Splitting " << nSplit
        << " faces across diagonals" << "." << nl << endl;

    // To undo: store original faces
    faceList originalFaces(UIndirectList<face>(mesh_.faces(), splitFaces));
    labelPairList facePairs(splitFaces.size(), labelPair(-1, -1));


    {
        polyTopoChange meshMod(mesh_);
        meshMod.setCapacity
        (
            meshMod.points().size(),
            meshMod.faces().size()+splitFaces.size(),
            mesh_.nCells()
        );

        // Insert the mesh changes
        doSplitFaces(splitFaces, splits, meshMod);

        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing might not do this)
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());
        setInstance(mesh_.facesInstance());


        // Update local mesh data
        // ~~~~~~~~~~~~~~~~~~~~~~

        forAll(originalFaces, i)
        {
            inplaceRenumber(map().reversePointMap(), originalFaces[i]);
        }

        {
            Map<label> splitFaceToIndex(2*splitFaces.size());
            forAll(splitFaces, i)
            {
                splitFaceToIndex.insert(splitFaces[i], i);
            }

            forAll(map().faceMap(), facei)
            {
                label oldFacei = map().faceMap()[facei];
                Map<label>::iterator oldFaceFnd = splitFaceToIndex.find
                (
                    oldFacei
                );
                if (oldFaceFnd != splitFaceToIndex.end())
                {
                    labelPair& twoFaces = facePairs[oldFaceFnd()];
                    if (twoFaces[0] == -1)
                    {
                        twoFaces[0] = facei;
                    }
                    else if (twoFaces[1] == -1)
                    {
                        twoFaces[1] = facei;
                    }
                    else
                    {
                        FatalErrorInFunction
                            << "problem: twoFaces:" << twoFaces
                            << exit(FatalError);
                    }
                }
            }
        }


        // Update baffle data
        // ~~~~~~~~~~~~~~~~~~

        if (duplicateFace.size())
        {
            meshRefinement::updateList
            (
                map().faceMap(),
                label(-1),
                duplicateFace
            );
        }

        const labelList& oldToNewFaces = map().reverseFaceMap();
        forAll(baffles, i)
        {
            labelPair& baffle = baffles[i];
            baffle.first() = oldToNewFaces[baffle.first()];
            baffle.second() = oldToNewFaces[baffle.second()];

            if (baffle.first() == -1 || baffle.second() == -1)
            {
                FatalErrorInFunction
                    << "Removed baffle : faces:" << baffle
                    << exit(FatalError);
            }
        }


        // Update insersections
        // ~~~~~~~~~~~~~~~~~~~~

        {
            DynamicList<label> changedFaces(facePairs.size());
            forAll(facePairs, i)
            {
                changedFaces.append(facePairs[i].first());
                changedFaces.append(facePairs[i].second());
            }

            // Update intersections on changed faces
            updateMesh(map, growFaceCellFace(changedFaces));
        }
    }



    // Undo loop
    // ~~~~~~~~~
    // Maintains two bits of information:
    // facePairs     : two faces originating from the same face
    // originalFaces : original face in current vertices


    for (label iteration = 0; iteration < 100; iteration++)
    {
        Info<< nl
            << "Undo iteration " << iteration << nl
            << "----------------" << endl;


        // Check mesh for errors
        // ~~~~~~~~~~~~~~~~~~~~~

        faceSet errorFaces
        (
            mesh_,
            "errorFaces",
            mesh_.nFaces()-mesh_.nInternalFaces()
        );
        bool hasErrors = motionSmoother::checkMesh
        (
            false,  // report
            mesh_,
            motionDict,
            errorFaces
        );
        if (!hasErrors)
        {
            break;
        }

        // Extend faces
        {
            const labelList grownFaces(growFaceCellFace(errorFaces));
            errorFaces.clear();
            errorFaces.insert(grownFaces);
        }


        // Merge face pairs
        // ~~~~~~~~~~~~~~~~
        // (if one of the faces is in the errorFaces set)

        polyTopoChange meshMod(mesh_);

        // Indices (in facePairs) of merged faces
        labelHashSet mergedIndices(facePairs.size());
        forAll(facePairs, index)
        {
            const labelPair& twoFaces = facePairs[index];

            if
            (
                errorFaces.found(twoFaces.first())
             || errorFaces.found(twoFaces.second())
            )
            {
                const face& originalFace = originalFaces[index];


                // Determine face properties
                label own = mesh_.faceOwner()[twoFaces[0]];
                label nei = -1;
                label patchi = -1;
                if (twoFaces[0] >= mesh_.nInternalFaces())
                {
                    patchi = mesh_.boundaryMesh().whichPatch(twoFaces[0]);
                }
                else
                {
                    nei = mesh_.faceNeighbour()[twoFaces[0]];
                }

                label zonei = mesh_.faceZones().whichZone(twoFaces[0]);
                bool zoneFlip = false;
                if (zonei != -1)
                {
                    const faceZone& fz = mesh_.faceZones()[zonei];
                    zoneFlip = fz.flipMap()[fz.whichFace(twoFaces[0])];
                }

                // Modify first face
                meshMod.modifyFace
                (
                    originalFace,               // modified face
                    twoFaces[0],                // label of face
                    own,                        // owner
                    nei,                        // neighbour
                    false,                      // face flip
                    patchi,                     // patch for face
                    zonei,                      // zone for face
                    zoneFlip                    // face flip in zone
                );
                // Merge second face into first
                meshMod.removeFace(twoFaces[1], twoFaces[0]);

                mergedIndices.insert(index);
            }
        }

        label n = returnReduce(mergedIndices.size(), sumOp<label>());

        Info<< "Detected " << n
            << " split faces that will be restored to their original faces."
            << nl << endl;

        if (n == 0)
        {
            // Nothing to be restored
            break;
        }

        nSplit -= n;


        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing might not do this)
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }

        // Reset the instance for if in overwrite mode
        mesh_.setInstance(timeName());
        setInstance(mesh_.facesInstance());


        // Update local mesh data
        // ~~~~~~~~~~~~~~~~~~~~~~

        {
            const labelList& oldToNewFaces = map().reverseFaceMap();
            const labelList& oldToNewPoints = map().reversePointMap();

            // Compact out merged faces
            DynamicList<label> changedFaces(mergedIndices.size());

            label newIndex = 0;
            forAll(facePairs, index)
            {
                const labelPair& oldSplit = facePairs[index];
                label new0 = oldToNewFaces[oldSplit[0]];
                label new1 = oldToNewFaces[oldSplit[1]];

                if (!mergedIndices.found(index))
                {
                    // Faces still split
                    if (new0 < 0 || new1 < 0)
                    {
                        FatalErrorInFunction
                            << "Problem: oldFaces:" << oldSplit
                            << " newFaces:" << labelPair(new0, new1)
                            << exit(FatalError);
                    }

                    facePairs[newIndex] = labelPair(new0, new1);
                    originalFaces[newIndex] = renumber
                    (
                        oldToNewPoints,
                        originalFaces[index]
                    );
                    newIndex++;
                }
                else
                {
                    // Merged face. Only new0 kept.
                    if (new0 < 0 || new1 == -1)
                    {
                        FatalErrorInFunction
                            << "Problem: oldFaces:" << oldSplit
                            << " newFace:" << labelPair(new0, new1)
                            << exit(FatalError);
                    }
                    changedFaces.append(new0);
                }
            }

            facePairs.setSize(newIndex);
            originalFaces.setSize(newIndex);


            // Update intersections
            updateMesh(map, growFaceCellFace(changedFaces));
        }

        // Update baffle data
        // ~~~~~~~~~~~~~~~~~~
        {
            if (duplicateFace.size())
            {
                meshRefinement::updateList
                (
                    map().faceMap(),
                    label(-1),
                    duplicateFace
                );
            }

            const labelList& reverseFaceMap = map().reverseFaceMap();
            forAll(baffles, i)
            {
                labelPair& baffle = baffles[i];
                baffle.first() = reverseFaceMap[baffle.first()];
                baffle.second() = reverseFaceMap[baffle.second()];

                if (baffle.first() == -1 || baffle.second() == -1)
                {
                    FatalErrorInFunction
                        << "Removed baffle : faces:" << baffle
                        << exit(FatalError);
                }
            }
        }

    }

    return nSplit;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshRefinement::meshRefinement
(
    fvMesh& mesh,
    const scalar mergeDistance,
    const bool readLevel,
    const bool overwrite,
    const meshControl& controller,
    const Tuple2<bool, scalar> crackDetection,
    const bool addedRays,
    const refinementSurfaces& surfaces,
    const shellSurfaces& shells,
    const shellSurfaces& limitShells,
    autoPtr<refinementFeatures>& featuresPtr,
    const scalar manualLevel0,
    const bool useVDB /*= false*/
)
:
    mesh_(mesh),
    cellRegionID_
    (
        IOobject
        (
            "cellRegionID",
            mesh_.facesInstance(),
            polyMesh::meshSubDir,
            mesh_,
            (readLevel ? IOobject::READ_IF_PRESENT : IOobject::NO_READ),
            IOobject::NO_WRITE
        ),
        labelList(mesh_.nCells(), 0)
    ),
    mergeDistance_(mergeDistance),
    overwrite_(overwrite),
    controller_(controller),
    crackDetection_(crackDetection),
    addedRays_(addedRays),
    oldInstance_(mesh.pointsInstance()),
    surfaces_(surfaces),
    shells_(shells),
    limitShells_(limitShells),
    featuresPtr_(featuresPtr),
    meshCutter_
    (
        mesh,
        readLevel,
        false,               // do not try to read history.
        (useVDB ? false : controller.refine()), // check mesh levels if refining mesh
        manualLevel0
    ),
    surfaceIndex_
    (
        IOobject
        (
            "surfaceIndex",
            mesh_.facesInstance(),
            fvMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelList(mesh_.nFaces(), -1)
    ),
    wrapIndex_
    (
        IOobject
        (
            "wrapIndex",
            mesh_.facesInstance(),
            fvMesh::meshSubDir,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        labelList(mesh_.nFaces(), -1)
    ),
    userFaceData_(0),
    oldPointsPtr_(nullptr),
    wrapActive_(false),
    createGapField_(false)
{
    // set minimum edge length
    setCrackTol();

    if (!useVDB)
    {
        // recalculate intersections for all faces
        updateIntersections(identity(mesh_.nFaces()));
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::meshRefinement::countHits
(
    const bool threaded /*= false*/
) const
{
    // Stats on edges to test. Count proc faces only once.
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_, threaded));

    label nHits = 0;

    if (threaded)
    {
#ifndef FOAM_USE_TBB
        WarningInFunction
            << "FOAM not linked against TBB! Cannot run multithreaded."
            << endl;
#else
        const auto grainSize =
            std::max<size_t>
            (
                mesh_.nFaces() / tbb::this_task_arena::max_concurrency(),
                1024
            );
        const tbb::blocked_range<size_t> range(0, mesh_.nFaces(), grainSize);

        nHits =
            tbb::parallel_reduce
            (
                range,
                /*identity*/0,
                [&](const tbb::blocked_range<size_t>& r, const label& init)
                {
                    label res = init;
                    for (size_t facei = r.begin(); facei < r.end(); ++facei)
                    {
                        if
                        (
                            surfaceIndex_[facei] >= 0
                         && isMasterFace.get(facei) == 1
                        )
                        {
                            res++;
                        }
                    }
                    return res;
                },
                /*reduction*/[](const label& x, const label& y)
                {
                    return x+y;
                },
                tbb::simple_partitioner()
            );
#endif
    }
    else
    {
        forAll(surfaceIndex_, facei)
        {
            if (surfaceIndex_[facei] >= 0 && isMasterFace.get(facei) == 1)
            {
                nHits++;
            }
        }
    }
    return nHits;
}


void Foam::meshRefinement::calcSnapWeights
(
    const scalar& snapWeights,
    scalarField& weightField
) const
{
    Info<<"Using snap weights for rebalance with scaling factor : "
        << snapWeights << nl << endl;

    DynamicList<label> testFaces(mesh_.nFaces());
    forAll(mesh_.faces(), facei)
    {
        const label zoneID = mesh_.faceZones().whichZone(facei);
        const label patchi = mesh_.boundaryMesh().whichPatch(facei);
        if
        (
            patchi != -1 && !mesh_.boundaryMesh()[patchi].coupled()
        )
        {
            testFaces.append(facei);
        }
        else if (zoneID != -1)
        {
            const faceZone& fZone = mesh_.faceZones()[zoneID];
            bool flipMap = fZone.flipMap()[fZone.whichFace(facei)];
            if (patchi != -1)
            {
                if (!flipMap)
                {
                    testFaces.append(facei);
                }
            }
            else
            {
                testFaces.append(facei);
            }
        }
    }

    const scalar edge0Len = meshCutter().level0EdgeLength();
    const labelList& cellLevel = meshCutter().cellLevel();

    pointField testCentres(testFaces.size());
    scalarField testRadius(testFaces.size());

    forAll(testFaces, i)
    {
        label facei = testFaces[i];
        point fc = mesh_.faceCentres()[facei];
        testCentres[i] = fc;
        label own = mesh_.faceOwner()[facei];
        label  level = cellLevel[own];
        scalar len = edge0Len / pow(2., level);
        testRadius[i] = (len*len);
    }

    const labelList& refSurfaces = surfaces().surfaces();
    scalarField faceWeight(mesh_.nFaces(), 0);
    scalar totalWeights = 0;
    forAll(refSurfaces, i)
    {
        label surfI = refSurfaces[i];
        const searchableSurface& geom = surfaces().geometry()[surfI];
        if (isA<triSurface>(geom))
        {
            triSurfaceMesh& surf = const_cast<triSurfaceMesh&>
            (
                refCast<const triSurfaceMesh>(geom)
            );
            List<labelList> hitIndexes(testCentres.size());
            labelList nOverlapChecks(testCentres.size(),label(0));
            surf.findSphere(testCentres,testRadius,hitIndexes,nOverlapChecks);

            forAll(testFaces, testi)
            {
                label nHits = hitIndexes[testi].size();
                label nOverlaps = nOverlapChecks[testi];
                if (nHits > 1)
                {
                    label fWeight = 5*nHits + nOverlaps;
                    label facei = testFaces[testi];
                    faceWeight[facei] += scalar(fWeight);
                    totalWeights += scalar(fWeight);
                }
            }
        }
    }

    label numCells = returnReduce(mesh_.nCells(), sumOp<label>());
    reduce(totalWeights, sumOp<scalar>());

    if (totalWeights > 0)
    {
        scalar weightFactor = snapWeights
            * scalar(numCells) / scalar(totalWeights);
        forAll(mesh_.faces(), facei)
        {
            if (faceWeight[facei] > 0)
            {
                const label zoneID = mesh_.faceZones().whichZone(facei);
                const label patchi = mesh_.boundaryMesh().whichPatch(facei);
                label own = mesh_.faceOwner()[facei];

                if
                (
                    patchi != -1 && !mesh_.boundaryMesh()[patchi].coupled()
                )
                {
                    weightField[own] += weightFactor*faceWeight[facei];
                }
                else if (zoneID != -1)
                {
                    const faceZone& fZone = mesh_.faceZones()[zoneID];
                    bool flipMap = fZone.flipMap()[fZone.whichFace(facei)];
                    if (patchi != -1)
                    {
                        if (!flipMap)
                        {
                            weightField[own] += weightFactor*faceWeight[facei];
                        }
                    }
                    else
                    {
                        if (flipMap)
                        {
                            label nei = mesh_.faceNeighbour()[facei];
                            weightField[nei] += weightFactor*faceWeight[facei];
                        }
                        else
                        {
                            weightField[own] += weightFactor*faceWeight[facei];
                        }
                    }
                }
            }
        }
    }

    if (debug)
    {
        const_cast<Time&>(mesh_.time())++;
        mesh_.write();

        volScalarField wField
        (
            IOobject
            (
                "snapWeightField",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll(wField, celli)
        {
            wField[celli] = weightField[celli];
        }
        wField.write();
    }

    return;
}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::meshRefinement::balance
(
    const bool applyDictionaryConstraints,
    const bool keepZoneFaces,
    const bool keepBaffles,
    const scalarField& cellWeights,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    bool updateSurf
)
{
    autoPtr<mapDistributePolyMesh> map;

    if (Pstream::parRun())
    {
        // Wanted distribution
        labelList distribution;


        // Faces where owner and neighbour are not 'connected' so can
        // go to different processors.
        boolList blockedFace;
        label nUnblocked = 0;

        // Faces that move as block onto single processor
        PtrList<labelList> specifiedProcessorFaces;
        labelList specifiedProcessor;

        // Pairs of baffles
        List<labelPair> couples;

        // Constraints from decomposeParDict
        decomposer.setConstraints
        (
            mesh_,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            couples
        );

        if (keepZoneFaces || keepBaffles)
        {
            if (keepZoneFaces)
            {
                // Determine decomposition to keep/move surface zones
                // on one processor. The reason is that snapping will make these
                // into baffles, move and convert them back so if they were
                // proc boundaries after baffling&moving the points might be no
                // longer synchronised so recoupling will fail. To prevent this
                // keep owner&neighbour of such a surface zone on the same
                // processor.

                const faceZoneMesh& fZones = mesh_.faceZones();
                const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

                // Get faces whose owner and neighbour should stay together,
                // i.e. they are not 'blocked'.

                forAll(fZones, zonei)
                {
                    const faceZone& fZone = fZones[zonei];

                    forAll(fZone, i)
                    {
                        label facei = fZone[i];
                        if (blockedFace[facei])
                        {
                            if
                            (
                                mesh_.isInternalFace(facei)
                             || pbm[pbm.whichPatch(facei)].coupled()
                            )
                            {
                                blockedFace[facei] = false;
                                nUnblocked++;
                            }
                        }
                    }
                }


                // If the faceZones are not synchronised the blockedFace
                // might not be synchronised. If you are sure the faceZones
                // are synchronised remove below check.
                syncTools::syncFaceList
                (
                    mesh_,
                    blockedFace,
                    andEqOp<bool>()     // combine operator
                );
            }
            reduce(nUnblocked, sumOp<label>());

            if (keepZoneFaces)
            {
                Info<< "Found " << nUnblocked
                    << " zoned faces to keep together." << endl;
            }


            // Extend un-blockedFaces with any cyclics
            {
                boolList separatedCoupledFace(mesh_.nFaces(), false);
                selectSeparatedCoupledFaces(separatedCoupledFace);

                label nSeparated = 0;
                forAll(separatedCoupledFace, facei)
                {
                    if (separatedCoupledFace[facei])
                    {
                        if (blockedFace[facei])
                        {
                            blockedFace[facei] = false;
                            nSeparated++;
                        }
                    }
                }
                reduce(nSeparated, sumOp<label>());
                Info<< "Found " << nSeparated
                    << " separated coupled faces to keep together." << endl;

                nUnblocked += nSeparated;
            }


            if (keepBaffles)
            {
                label nBnd = mesh_.nFaces()-mesh_.nInternalFaces();

                labelList coupledFace(mesh_.nFaces(), -1);
                {
                    // Get boundary baffles that need to stay together
                    List<labelPair> allCouples
                    (
                        localPointRegion::findDuplicateFacePairs(mesh_)
                    );

                    // Merge with any couples from
                    // decompositionMethod::setConstraints
                    forAll(couples, i)
                    {
                        const labelPair& baffle = couples[i];
                        coupledFace[baffle.first()] = baffle.second();
                        coupledFace[baffle.second()] = baffle.first();
                    }
                    forAll(allCouples, i)
                    {
                        const labelPair& baffle = allCouples[i];
                        coupledFace[baffle.first()] = baffle.second();
                        coupledFace[baffle.second()] = baffle.first();
                    }
                }

                couples.setSize(nBnd);
                label nCpl = 0;
                forAll(coupledFace, facei)
                {
                    if (coupledFace[facei] != -1 && facei < coupledFace[facei])
                    {
                        couples[nCpl++] = labelPair(facei, coupledFace[facei]);
                    }
                }
                couples.setSize(nCpl);
            }
            label nCouples = returnReduce(couples.size(), sumOp<label>());

            if (keepBaffles)
            {
                Info<< "Found " << nCouples << " baffles to keep together."
                    << endl;
            }
        }


        // Make sure blockedFace not set on couples
        forAll(couples, i)
        {
            const labelPair& baffle = couples[i];
            blockedFace[baffle.first()] = false;
            blockedFace[baffle.second()] = false;
        }

        distribution = decomposer.decompose
        (
            mesh_,
            cellWeights,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            couples                 // explicit connections
        );

        if (applyDictionaryConstraints)
        {
            decomposer.applyConstraints
            (
                mesh_,
                blockedFace,
                specifiedProcessorFaces,
                specifiedProcessor,
                couples,
                distribution
            );
        }

        if (debug)
        {
            labelList nProcCells(distributor.countCells(distribution));
            Pout<< "Wanted distribution:" << nProcCells << endl;

            Pstream::listCombineReduce(nProcCells, plusOp<label>());

            Pout<< "Wanted resulting decomposition:" << endl;
            forAll(nProcCells, proci)
            {
                Pout<< "    " << proci << '\t' << nProcCells[proci] << endl;
            }
            Pout<< endl;
        }

        // Do actual sending/receiving of mesh
        map = distributor.distribute(distribution);

        // Update numbering of meshRefiner
        distribute(map,updateSurf);

        // Update manifold features if require tracking
        if (features().valid())
        {
            features()().manFeatures().distribute(map);
        }

        // Set correct instance (for if overwrite)
        mesh_.setInstance(timeName());
        setInstance(mesh_.facesInstance());

        if (debug && keepZoneFaces)
        {
            const faceZoneMesh& fZones = mesh_.faceZones();
            const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

            // Check that faceZone faces are indeed internal
            forAll(fZones, zonei)
            {
                const faceZone& fZone = fZones[zonei];

                forAll(fZone, i)
                {
                    label facei = fZone[i];
                    label patchi = pbm.whichPatch(facei);

                    if (patchi >= 0 && pbm[patchi].coupled())
                    {
                        WarningInFunction
                            << "Face at " << mesh_.faceCentres()[facei]
                            << " on zone " << fZone.name()
                            << " is on coupled patch " << pbm[patchi].name()
                            << endl;
                    }
                }
            }
        }

    }

    return map;
}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::meshRefinement::balance
(
    const boolList blockedFace,
    const scalarField& cellWeights,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    bool updateSurf
)
{
    autoPtr<mapDistributePolyMesh> map;

    if (Pstream::parRun())
    {
        // Faces that move as block onto single processor
        PtrList<labelList> specifiedProcessorFaces;
        labelList specifiedProcessor;

        labelList distribution = decomposer.decompose
        (
            mesh_,
            cellWeights,
            blockedFace,
            specifiedProcessorFaces,
            specifiedProcessor,
            List<labelPair>()  // explicit connections
        );

        labelList nProcCells(distributor.countCells(distribution));
        Pstream::listCombineReduce(nProcCells, plusOp<label>());

        Info<< "Calculated decomposition:" << endl;
        forAll(nProcCells, procI)
        {
            Info<< "    " << procI << '\t' << nProcCells[procI] << endl;
        }
        Info<< endl;

        // Do actual sending/receiving of mesh
        map = distributor.distribute(distribution);

        // Update numbering of meshRefiner
        distribute(map, updateSurf);

        // Set correct instance (for if overwrite)
        mesh_.setInstance(timeName());
        setInstance(mesh_.facesInstance());
    }

    return map;
}


void Foam::meshRefinement::redistributeToMany
(
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor
)
{
#ifdef FOAM_USE_TBB
    Timer timer("meshRefinement::redistributeToMany");
#endif
    const label nDestProcs = decomposer.nDomains();

    const cellZoneMesh& cellZones = mesh_.cellZones();

    labelList decomp(mesh_.nCells());

    bool balanceSurfaceCells = cellZones.findZoneID("innerGrid") >= 0;

    if (balanceSurfaceCells)
    {
        Info<< "\nBalancing surface cells" << endl;

        const pointField& cc = mesh_.cellCentres();

        //All mesh is already divided in inner and outer zones
        forAll(cellZones, i)
        {
            pointField zonePoints(cellZones[i].size());

            forAll(zonePoints, j)
            {
                zonePoints[j] = cc[cellZones[i][j]];
            }

            labelList zoneDist =
                decomposer.decompose
                (
                    zonePoints,
                    scalarField(zonePoints.size(), 1) //cellWeights
                );

            forAll(zoneDist, j)
            {
                decomp[cellZones[i][j]] = zoneDist[j];
            }
        }
    }
    else
    {
        decomp =
            decomposer.decompose
            (
                mesh_,
                scalarField() //cellWeights
            );
    }

    distributor.setNFinalProcs(nDestProcs);

    // Initialize empty list of meshes and fields
    // for each processor
    PtrList<Time>                  myProcTimes;
    PtrList<fvMesh>                myProcMeshes;
    PtrList<mapDistributePolyMesh> myProcDist;

    List<PtrList<volScalarField>>  myProcVolScalarFields;
    List<PtrList<volVectorField>>  myProcVolVectorFields;

    distributor.initMeshList
    (
        fvMesh::defaultRegion,
        nDestProcs,
        mesh_,

        myProcTimes,
        myProcMeshes,

        myProcVolScalarFields,
        myProcVolVectorFields
    );

    // Do all the distribution of mesh and fields
    autoPtr<mapDistributePolyMesh> rawMap =
        distributor.distribute
        (
            decomp,
            myProcMeshes,
            myProcDist
        );

    labelListList meshNoInProc
    (
        Pstream::nProcs(),
        labelList(myProcMeshes.size(), 0)
    );

    labelList nMeshesInProc(Pstream::nProcs());
    labelList meshToProc(identity(nDestProcs));
    meshNoInProc =
        fvMeshDistribute::calcMeshToProcMap
        (
            nDestProcs,
            nMeshesInProc,
            meshToProc
        );

    PtrList<mapPolyMesh> myProcMap(nMeshesInProc[Pstream::myProcNo()]);

    //TODO tbb::parallel_for
    forAll(myProcMeshes, i)
    {
        fvMesh& meshi = myProcMeshes[i];

        Time& runTimei = const_cast<Time&>(meshi.time());

        runTimei.setTime(mesh_.time());

        meshi.setInstance(runTimei.timeName());

        label meshNo = meshNoInProc[Pstream::myProcNo()][i];
        //if (fvMeshDistribute::debug)
        {
            Pout<< "Writing domain "
                << meshNo
                << " to " << meshi.time().path()/runTimei.timeName()
                << endl;
        }

        const cellZoneMesh& zonesI = meshi.cellZones();

        bool renumberZones = true;//false;//

        if (balanceSurfaceCells)
        {
            Info<< "\nRenumbering surface cells" << endl;
            label blockSize = 10000;

            autoPtr<decompositionMethod> decomposePtr;
            dictionary decomposeDict;
            decomposeDict.add("method", "kahip");

            bool oldParRun = UPstream::parRun();
            UPstream::parRun() = false;

            labelList cellToRegion(meshi.nCells());

            if (renumberZones)
            {
                labelList cellToZone(meshi.nCells());
                forAll(zonesI, i)
                {
                    forAll(zonesI[i], j)
                    {
                        cellToZone[zonesI[i][j]] = i;
                    }
                }

                label blockOffset = 0;
                //meshi is already divided in inner and outer zones
                forAll(zonesI, i)
                {
                    label nBlocks = zonesI[i].size()/blockSize;
                    Info<<"region: " << zonesI.names()[i]
                        <<" - nBlocks: " << nBlocks << endl;
                    decomposeDict.set("numberOfSubdomains", nBlocks);

                    decomposePtr.reset
                    (
                        decompositionMethod::New
                        (
                            decomposeDict
                        )
                    );

                    // Per region do a block decomposition.
                    fvMeshSubset subsetter(meshi);
                    subsetter.setLargeCellSubset(cellToZone, i);

                    const fvMesh& subMesh = subsetter.subMesh();

                    labelList subCellDist = decomposePtr().decompose
                    (
                        subMesh,
                        subMesh.cellCentres()
                    );

                    forAll(subCellDist, j)
                    {
                        cellToRegion[zonesI[i][j]] = subCellDist[j] + blockOffset;
                    }

                    blockOffset += nBlocks;
                }
            }
            else
            {
                label nBlocks = meshi.nCells()/blockSize;
                Info<< "nBlocks   = " << nBlocks << endl;
                decomposeDict.set("numberOfSubdomains", nBlocks);

                decomposePtr.reset
                (
                    decompositionMethod::New
                    (
                        decomposeDict
                    )
                );

                cellToRegion =
                    (
                        decomposePtr().decompose
                        (
                            meshi,
                            meshi.cellCentres()
                        )
                    );
            }

            //GGG //TODO update processor boundaries
            //// Restore state
            //UPstream::parRun() = oldParRun;

            autoPtr<renumberMethod> renumberPtr;
            dictionary renumberDict;
            renumberPtr.reset(new CuthillMcKeeRenumber(renumberDict));

            labelList cellOrder = regionRenumber(renumberPtr(), meshi, cellToRegion);

            // Determine new to old face order with new cell numbering
            labelList faceOrder = getRegionFaceOrder
            (
                meshi,
                cellOrder,
                cellToRegion
            );

            // Change the mesh.
            autoPtr<mapPolyMesh> map = reorderMesh(meshi, cellOrder, faceOrder);

            bool orderPoints = true;
            if (orderPoints)
            {
                polyTopoChange meshMod(meshi);
                autoPtr<mapPolyMesh> pointOrderMap = meshMod.changeMesh
                (
                    meshi,
                    false,      // inflate
                    true,       // syncParallel
                    false,      // orderCells
                    orderPoints // orderPoints
                );

                // Combine point reordering into map.
                const_cast<labelList&>(map().pointMap()) = UIndirectList<label>
                (
                    map().pointMap(),
                    pointOrderMap().pointMap()
                )();

                inplaceRenumber
                (
                    pointOrderMap().reversePointMap(),
                    const_cast<labelList&>(map().reversePointMap())
                 );
            }

            // Update fields
            meshi.updateMesh(map);

            //TODO update processorBoundaries
            //GGG
            // Restore state
            UPstream::parRun() = oldParRun;
            //GGG

            // Move mesh (since morphing might not do this)
            if (map().hasMotionPoints())
            {
                meshi.movePoints(map().preMotionPoints());
            }

            myProcMap.set(i, map);
        }

        meshi.write();

        //if (debug)
        {
            //count boundary faces per proc and cells in innerGrid
            const fvBoundaryMesh& patches = meshi.boundary();

            label boundaryFaces = 0;

            forAll(patches, patchi)
            {
                const polyPatch& pp = patches[patchi].patch();

                if
                (
                    pp.coupled()
                 || pp.name() == "defaultFaces"
                 || pp.name() == "allBoundaries"
                )
                {
                    Info<< "Skipping patch: " << pp.name() << endl;
                    continue;
                }

                boundaryFaces += pp.size();
            }

            label innerCells = 0;

            if (balanceSurfaceCells)
            {
                innerCells += zonesI[0].size();
            }

            Pout<< "Domain " << meshNo
                << " - boundaryFaces: " << boundaryFaces
                << " - nCells: " << meshi.nCells()
                << " - innerCells: " << innerCells
                << endl;

            tmp<volScalarField> tfld
            (
                new volScalarField
                (
                    IOobject
                    (
                        "newCellID",
                        meshi.time().timeName(),
                        meshi,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        false
                    ),
                    meshi,
                    dimensionedScalar("zero", dimless, 0),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            volScalarField& fld = tfld.ref();

            labelList elems(identity(meshi.nCells()));

            forAll(fld, celli)
            {
                fld[celli] = elems[celli];
            }

            tfld().write();

            tmp<volScalarField> procID
            (
                new volScalarField
                (
                    IOobject
                    (
                        "newProcID",
                        meshi.time().timeName(),
                        meshi,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        false
                    ),
                    meshi,
                    dimensionedScalar("newProcID", dimless, meshNo),
                    zeroGradientFvPatchScalarField::typeName
                )
            );
            procID().write();
        }
    } //forAll myProcMeshes

    // wait for all the processors to finish writing
    Pstream::barrier();

    printMeshData(myProcMeshes);

    //redistribute hexRef8 data
    bool oldParRun = UPstream::parRun();
    UPstream::parRun() = false;

    hexRef8DataList refData
    (
        myProcMeshes,
        meshToProc,
        meshNoInProc
    );

    UPstream::parRun() = oldParRun;

    //TODO refinementHistory not currently distributed
    refData.distribute(myProcDist);
    if (balanceSurfaceCells)
    {
        refData.reorder(myProcMap);
    }

    refData.write();
}


void Foam::meshRefinement::printMeshData(PtrList<fvMesh>& meshes)
{
    labelListList localSizes(Pstream::nProcs(), labelList(meshes.size()));
    labelListList localBoundFaces(Pstream::nProcs(), labelList(meshes.size()));

    labelListList myProcPatchNeiProcNo(meshes.size());
    labelListList myProcPatchSize(meshes.size());

    forAll(meshes, i)
    {
        fvMesh& mesh = meshes[i];
        localSizes[Pstream::myProcNo()][i] = mesh.nCells();

        localBoundFaces[Pstream::myProcNo()][i] =
            mesh.nFaces() - mesh.nInternalFaces();

        labelList processorNeighbours(mesh.boundaryMesh().size());
        labelList processorSize(mesh.boundaryMesh().size());

        label nNeighbours = 0;

        forAll(mesh.boundaryMesh(), patchi)
        {
            if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchi]))
            {
                const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchi]
                    );
                processorNeighbours[nNeighbours] = ppp.neighbProcNo();
                processorSize[nNeighbours++] = ppp.size();
            }
        }
        processorNeighbours.setSize(nNeighbours);
        processorSize.setSize(nNeighbours);

        myProcPatchNeiProcNo[i] = processorNeighbours;
        myProcPatchSize[i] = processorSize;
    }

    Pstream::gatherList(localSizes);
    Pstream::gatherList(localBoundFaces);

    labelListListList patchNeiProcNos(Pstream::nProcs());
    labelListListList patchSizes     (Pstream::nProcs());
    patchNeiProcNos[Pstream::myProcNo()] = myProcPatchNeiProcNo;
    patchSizes     [Pstream::myProcNo()] = myProcPatchSize;
    Pstream::gatherList(patchNeiProcNos);
    Pstream::gatherList(patchSizes);

    label nFinalProcs = meshes.size();
    reduce(nFinalProcs, sumOp<label>());

    if (!Pstream::master()) return;

    labelList nMeshesInProc(Pstream::nProcs());
    labelList meshToProc(identity(nFinalProcs));
    labelListList meshNoInProc =
        fvMeshDistribute::calcMeshToProcMap
        (
            nFinalProcs,
            nMeshesInProc,
            meshToProc
        );

    label         globalCells = 0;
    labelList     localSize          (nFinalProcs);
    labelList     globalBoundaryFaces(nFinalProcs);
    labelListList patchNeiProcNo     (nFinalProcs);
    labelListList patchSize          (nFinalProcs);
    forAll(localSizes, proci)
    {
        forAll(localSizes[proci], j)
        {
            label meshi = meshNoInProc[proci][j];

            globalCells += localSizes[proci][j];

            localSize[meshi] = localSizes[proci][j];

            globalBoundaryFaces[meshi] = localBoundFaces[proci][j];

            patchNeiProcNo[meshi] = patchNeiProcNos[proci][j];

            patchSize[meshi] = patchSizes[proci][j];
        }
    }

    FOAM_ASSERT (globalCells >= 0) {
        FatalErrorInFunction
            << "Overflow : sum of sizes " << globalCells
            << " exceeds capability of label (" << labelMax
            << "). Please recompile with larger datatype for label."
            << exit(FatalError);
    }

    // Print stats

    label maxProcCells = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;
    label maxProcFaces = 0;

    for (label procI = 0; procI < nFinalProcs; procI++)
    {
        Info<< endl
            << "Processor " << procI << nl
            << "    Number of cells = " << localSize[procI]
            << endl;

        label nProcFaces = 0;

        const labelList& nei = patchNeiProcNo[procI];

        forAll(patchNeiProcNo[procI], i)
        {
            Info<< "    Number of faces shared with processor "
                << patchNeiProcNo[procI][i] << " = " << patchSize[procI][i]
                << endl;

            nProcFaces += patchSize[procI][i];
        }

        Info<< "    Number of processor patches = " << nei.size() << nl
            << "    Number of processor faces = " << nProcFaces << nl
            << "    Number of boundary faces = "
            << globalBoundaryFaces[procI]-nProcFaces
            << endl;

        maxProcCells = max(maxProcCells, localSize[procI]);
        totProcFaces += nProcFaces;
        totProcPatches += nei.size();
        maxProcPatches = max(maxProcPatches, nei.size());
        maxProcFaces = max(maxProcFaces, nProcFaces);
    }

    // Stats

    // In case of all faces on one processor, use max to avoid division by 0.
    scalar avgProcCells   = max(1, scalar(globalCells)/nFinalProcs);
    scalar avgProcPatches = max(1, scalar(totProcPatches)/nFinalProcs);
    scalar avgProcFaces   = max(1, scalar(totProcFaces)/nFinalProcs);

    if (avgProcCells > 1)
    {
        Info<< nl
            << "Number of processor faces = " << totProcFaces/2 << nl
            << "Max number of cells = " << maxProcCells
            << " (" << 100.0*(maxProcCells-avgProcCells)/avgProcCells
            << "% above average " << avgProcCells << ")" << nl
            << "Max number of processor patches = " << maxProcPatches
            << " (" << 100.0*(maxProcPatches-avgProcPatches)/avgProcPatches
            << "% above average " << avgProcPatches << ")" << nl
            << "Max number of faces between processors = " << maxProcFaces
            << " (" << 100.0*(maxProcFaces-avgProcFaces)/avgProcFaces
            << "% above average " << avgProcFaces << ")" << nl
            << endl;
    }
}


// Helper function to get intersected faces and edge neighbours
Foam::labelList Foam::meshRefinement::intersectedAndEdgeCellFaces() const
{
    boolList markedEdges(mesh_.nEdges(), false);
    boolList markedFaces(mesh_.nFaces(), false);

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            markedFaces[faceI] = true;
            const labelList& fEdges = mesh_.faceEdges()[faceI];
            forAll(fEdges, fEI)
            {
                markedEdges[fEdges[fEI]] = true;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        markedEdges,
        orEqOp<bool>(),
        false
    );

    forAll(markedEdges, edgeI)
    {
        if (markedEdges[edgeI])
        {
            const labelList& eCells = mesh_.edgeCells()[edgeI];
            forAll(eCells, eCI)
            {
                cell c = mesh_.cells()[eCells[eCI]];
                forAll(c, cFI)
                {
                    markedFaces[c[cFI]] = true;
                }
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        markedFaces,
        orEqOp<bool>()
    );

    label nBoundaryFaces = 0;
    forAll(markedFaces, faceI)
    {
        if (markedFaces[faceI])
        {
            nBoundaryFaces++;
        }
    }

    labelList surfaceFaces(nBoundaryFaces);
    nBoundaryFaces = 0;
    forAll(markedFaces, faceI)
    {
        if (markedFaces[faceI])
        {
            surfaceFaces[nBoundaryFaces++] = faceI;
        }
    }

    return surfaceFaces;
}


// Helper function to get intersected faces and edge neighbours
Foam::labelList Foam::meshRefinement::intersectedAndNeighbouring() const
{
    boolList markedEdges(mesh_.nEdges(), false);
    boolList markedFaces(mesh_.nFaces(), false);

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            markedFaces[faceI] = true;
            const labelList& fEdges = mesh_.faceEdges()[faceI];
            forAll(fEdges, fEI)
            {
                markedEdges[fEdges[fEI]] = true;
            }
        }
    }

    syncTools::syncEdgeList
    (
        mesh_,
        markedEdges,
        orEqOp<bool>(),
        false
    );

    // Identify cells with open faces connected to more
    // than one marked edge and baffle all cell faces
    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] == -1)
        {
            const labelList& fEdges = mesh_.faceEdges()[faceI];

            label nMarkedEdges = 0;
            forAll(fEdges, fEI)
            {
                if (markedEdges[fEdges[fEI]])
                {
                    nMarkedEdges++;
                }
            }
            if (nMarkedEdges > 1)
            {
                label own = mesh_.faceOwner()[faceI];
                cell cOwn = mesh_.cells()[own];
                forAll(cOwn, cFI)
                {
                    markedFaces[cOwn[cFI]] = true;
                }

                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    cell cNei = mesh_.cells()[nei];
                    forAll(cNei, cFI)
                    {
                        markedFaces[cNei[cFI]] = true;
                    }
                }
            }
        }
    }
/*
    forAll(markedEdges, edgeI)
    {
        if (markedEdges[edgeI])
        {
            const labelList& eCells = mesh_.edgeCells()[edgeI];
            forAll(eCells, eCI)
            {
                cell c = mesh_.cells()[eCells[eCI]];
                forAll(c, cFI)
                {
                    markedFaces[c[cFI]] = true;
                }
            }
        }
    }
*/
    syncTools::syncFaceList
    (
        mesh_,
        markedFaces,
        orEqOp<bool>()
    );

    label nBoundaryFaces = 0;
    forAll(markedFaces, faceI)
    {
        if (markedFaces[faceI])
        {
            nBoundaryFaces++;
        }
    }

    labelList surfaceFaces(nBoundaryFaces);
    nBoundaryFaces = 0;
    forAll(markedFaces, faceI)
    {
        if (markedFaces[faceI])
        {
            surfaceFaces[nBoundaryFaces++] = faceI;
        }
    }

    return surfaceFaces;
}


Foam::labelList Foam::meshRefinement::intersectedFaces() const
{
    label nBoundaryFaces = 0;

    forAll(surfaceIndex_, facei)
    {
        if (surfaceIndex_[facei] != -1)
        {
            nBoundaryFaces++;
        }
    }

    labelList surfaceFaces(nBoundaryFaces);
    nBoundaryFaces = 0;

    forAll(surfaceIndex_, facei)
    {
        if (surfaceIndex_[facei] != -1)
        {
            surfaceFaces[nBoundaryFaces++] = facei;
        }
    }
    return surfaceFaces;
}

Foam::labelList Foam::meshRefinement::intersectedPoints() const
{
    const faceList& faces = mesh_.faces();

    // Mark all points on faces that will become baffles
    labelHashSet boundaryPoints(mesh_.nPoints()/100);

    forAll(surfaceIndex_, facei)
    {
        if (surfaceIndex_[facei] != -1)
        {
            const face& f = faces[facei];

            forAll(f, fp)
            {
                boundaryPoints.insert(f[fp]);
            }
        }
    }

    // Mark existing boundary points
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if (!pp.coupled())
        {
            label startFace = pp.start();
            forAll(pp, i)
            {
                label faceI = startFace + i;
                const face& f = faces[faceI];

                forAll(f, fp)
                {
                    boundaryPoints.insert(f[fp]);
                }
            }
        }
    }

    pointSet bPoints(mesh_, "usedPoints", boundaryPoints);
    bPoints.sync(mesh_);

    return bPoints.toc();
}


// Helper function to get cells from intersected faces
Foam::labelList Foam::meshRefinement::intersectedCells() const
{
    // Mark all points on faces that will become baffles
    labelHashSet boundaryCells(mesh_.nCells()/100);

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            label own = mesh_.faceOwner()[faceI];
            boundaryCells.insert(own);

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];
                boundaryCells.insert(nei);
            }
        }
    }
    cellSet bCells(mesh_, "usedCells", boundaryCells);

    return bCells.toc();
}


void Foam::meshRefinement::calcPatchAddressing
(
    const polyMesh &mesh,
    const labelList &patchIDs,
    labelList &addressing
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    // Count faces.
    label nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];
        nFaces += pp.size();
    }

    // Collect faces.
    addressing.setSize(nFaces);
    nFaces = 0;

    forAll(patchIDs, i)
    {
        const polyPatch& pp = patches[patchIDs[i]];
        label meshFacei = pp.start();

        forAll(pp, i)
        {
            addressing[nFaces++] = meshFacei++;
        }
    }
}


//- Create patch from set of patches
Foam::autoPtr<Foam::indirectPrimitivePatch> Foam::meshRefinement::makePatch
(
    const polyMesh& mesh,
    const labelList& patchIDs
)
{
    labelList addressing;
    calcPatchAddressing(mesh, patchIDs, addressing);

    return autoPtr<indirectPrimitivePatch>
    (
        new indirectPrimitivePatch
        (
            IndirectList<face>(mesh.faces(), addressing),
            mesh.points()
        )
    );
}


Foam::tmp<Foam::pointVectorField> Foam::meshRefinement::makeDisplacementField
(
    const pointMesh& pMesh,
    const labelList& adaptPatchIDs
)
{
    const polyMesh& mesh = pMesh();

    // Construct displacement field.
    const pointBoundaryMesh& pointPatches = pMesh.boundary();

    wordList patchFieldTypes
    (
        pointPatches.size(),
        slipPointPatchVectorField::typeName
    );

    forAll(adaptPatchIDs, i)
    {
        patchFieldTypes[adaptPatchIDs[i]] =
            fixedValuePointPatchVectorField::typeName;
    }

    forAll(pointPatches, patchi)
    {
        if (isA<processorPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = calculatedPointPatchVectorField::typeName;
        }
        else if (isA<cyclicPointPatch>(pointPatches[patchi]))
        {
            patchFieldTypes[patchi] = cyclicSlipPointPatchVectorField::typeName;
        }
    }

    // Note: time().timeName() instead of meshRefinement::timeName() since
    // postprocessable field.
    tmp<pointVectorField> tfld
    (
        new pointVectorField
        (
            IOobject
            (
                "pointDisplacement",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector("displacement", dimLength, Zero),
            patchFieldTypes
        )
    );

    return tfld;
}


void Foam::meshRefinement::checkCoupledFaceZones(const polyMesh& mesh)
{
    const faceZoneMesh& fZones = mesh.faceZones();

    // Check any zones are present anywhere and in same order

    {
        List<wordList> zoneNames(Pstream::nProcs());
        zoneNames[Pstream::myProcNo()] = fZones.names();
        Pstream::allGatherList(zoneNames);
        // All have same data now. Check.
        forAll(zoneNames, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                if (zoneNames[proci] != zoneNames[Pstream::myProcNo()])
                {
                    FatalErrorInFunction
                        << "faceZones are not synchronised on processors." << nl
                        << "Processor " << proci << " has faceZones "
                        << zoneNames[proci] << nl
                        << "Processor " << Pstream::myProcNo()
                        << " has faceZones "
                        << zoneNames[Pstream::myProcNo()] << nl
                        << exit(FatalError);
                }
            }
        }
    }

    // Check that coupled faces are present on both sides.

    labelList faceToZone(mesh.nFaces()-mesh.nInternalFaces(), -1);

    forAll(fZones, zonei)
    {
        const faceZone& fZone = fZones[zonei];

        forAll(fZone, i)
        {
            label bFacei = fZone[i]-mesh.nInternalFaces();

            if (bFacei >= 0)
            {
                if (faceToZone[bFacei] == -1)
                {
                    faceToZone[bFacei] = zonei;
                }
                else if (faceToZone[bFacei] == zonei)
                {
                    FatalErrorInFunction
                        << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is twice in zone!"
                        << abort(FatalError);
                }
                else
                {
                    FatalErrorInFunction
                        << "Face " << fZone[i] << " in zone "
                        << fZone.name()
                        << " is also in zone "
                        << fZones[faceToZone[bFacei]].name()
                        << abort(FatalError);
                }
            }
        }
    }

    labelList neiFaceToZone(faceToZone);
    syncTools::swapBoundaryFaceList(mesh, neiFaceToZone);

    forAll(faceToZone, i)
    {
        if (faceToZone[i] != neiFaceToZone[i])
        {
            FatalErrorInFunction
                << "Face " << mesh.nInternalFaces()+i
                << " is in zone " << faceToZone[i]
                << ", its coupled face is in zone " << neiFaceToZone[i]
                << abort(FatalError);
        }
    }
}


void Foam::meshRefinement::calculateEdgeWeights
(
    const polyMesh& mesh,
    const PackedBoolList& isMasterEdge,
    const boolList& isBoundaryEdge,
    const boolList& isBoundaryPoint,
    const labelList& meshEdges,
    const labelList& meshPoints,
    const edgeList& edges,
    scalarField& edgeWeights,
    scalarField& invSumWeight
)
{
    const pointField& pts = mesh.points();

    // Calculate edgeWeights and inverse sum of edge weights
    edgeWeights.setSize(isMasterEdge.size());
    invSumWeight.setSize(meshPoints.size());

    forAll(edges, edgei)
    {
        const edge& e = edges[edgei];
        scalar eMag = max
        (
            SMALL,
            mag
            (
                pts[meshPoints[e[1]]]
              - pts[meshPoints[e[0]]]
            )
        );
        edgeWeights[edgei] = 1.0/eMag;
    }

    // Sum per point all edge weights
    weightedSum
    (
        mesh,
        isMasterEdge,
        isBoundaryEdge,
        isBoundaryPoint,
        meshEdges,
        meshPoints,
        edges,
        edgeWeights,
        scalarField(meshPoints.size(), 1.0),  // data
        invSumWeight
    );

    // Inplace invert
    forAll(invSumWeight, pointi)
    {
        scalar w = invSumWeight[pointi];

        if (w > 0.0)
        {
            invSumWeight[pointi] = 1.0/w;
        }
    }
}


Foam::label Foam::meshRefinement::appendPatch
(
    fvMesh& mesh,
    const word& patchName,
    const dictionary& patchDict
)
{
    // Clear local fields and e.g. polyMesh parallelInfo.
    mesh.clearOut();

    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    label patchi = polyPatches.size();

    // Add polyPatch at the end
    polyPatches.setSize(patchi+1);
    polyPatches.set
    (
        patchi,
        polyPatch::New
        (
            patchName,
            patchDict,
            patchi,
            polyPatches
        )
    );

    fvPatches.setSize(patchi+1);
    fvPatches.set
    (
        patchi,
        fvPatch::New
        (
            polyPatches[patchi],  // point to newly added polyPatch
            mesh.boundary()
        )
    );

    addPatchFields<volScalarField>
    (
        mesh,
        calculatedFvPatchField<scalar>::typeName
    );
    addPatchFields<volVectorField>
    (
        mesh,
        calculatedFvPatchField<vector>::typeName
    );
    addPatchFields<volSphericalTensorField>
    (
        mesh,
        calculatedFvPatchField<sphericalTensor>::typeName
    );
    addPatchFields<volSymmTensorField>
    (
        mesh,
        calculatedFvPatchField<symmTensor>::typeName
    );
    addPatchFields<volTensorField>
    (
        mesh,
        calculatedFvPatchField<tensor>::typeName
    );

    // Surface fields

    addPatchFields<surfaceScalarField>
    (
        mesh,
        calculatedFvPatchField<scalar>::typeName
    );
    addPatchFields<surfaceVectorField>
    (
        mesh,
        calculatedFvPatchField<vector>::typeName
    );
    addPatchFields<surfaceSphericalTensorField>
    (
        mesh,
        calculatedFvPatchField<sphericalTensor>::typeName
    );
    addPatchFields<surfaceSymmTensorField>
    (
        mesh,
        calculatedFvPatchField<symmTensor>::typeName
    );
    addPatchFields<surfaceTensorField>
    (
        mesh,
        calculatedFvPatchField<tensor>::typeName
    );
    return patchi;
}


Foam::label Foam::meshRefinement::addPatch
(
    fvMesh& mesh,
    const word& patchName,
    const dictionary& patchInfo
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches = const_cast<fvBoundaryMesh&>(mesh.boundary());

    const label patchi = polyPatches.findPatchID(patchName);
    if (patchi != -1)
    {
        // Already there
        return patchi;
    }


    label insertPatchi = polyPatches.size();
    label startFacei = mesh.nFaces();

    forAll(polyPatches, patchi)
    {
        const polyPatch& pp = polyPatches[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            insertPatchi = patchi;
            startFacei = pp.start();
            break;
        }
    }

    dictionary patchDict(patchInfo);
    patchDict.set("nFaces", 0);
    patchDict.set("startFace", startFacei);

    // Below is all quite a hack. Feel free to change once there is a better
    // mechanism to insert and reorder patches.

    label addedPatchi = appendPatch(mesh, patchName, patchDict);

    // Create reordering list
    // patches before insert position stay as is
    labelList oldToNew(addedPatchi+1);
    for (label i = 0; i < insertPatchi; i++)
    {
        oldToNew[i] = i;
    }
    // patches after insert position move one up
    for (label i = insertPatchi; i < addedPatchi; i++)
    {
        oldToNew[i] = i+1;
    }
    // appended patch gets moved to insert position
    oldToNew[addedPatchi] = insertPatchi;

    // Shuffle into place
    polyPatches.reorder(oldToNew, true);
    fvPatches.reorder(oldToNew);

    reorderPatchFields<volScalarField>(mesh, oldToNew);
    reorderPatchFields<volVectorField>(mesh, oldToNew);
    reorderPatchFields<volSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<volSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<volTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceScalarField>(mesh, oldToNew);
    reorderPatchFields<surfaceVectorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSphericalTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceSymmTensorField>(mesh, oldToNew);
    reorderPatchFields<surfaceTensorField>(mesh, oldToNew);

    return insertPatchi;
}


Foam::label Foam::meshRefinement::addMeshedPatch
(
    const word& name,
    const dictionary& patchInfo
)
{
    label meshedi = findIndex(meshedPatches_, name);

    if (meshedi != -1)
    {
        // Already there. Get corresponding polypatch
        return mesh_.boundaryMesh().findPatchID(name);
    }
    else
    {
        polyBoundaryMesh& polyPatches =
            const_cast<polyBoundaryMesh&>(mesh_.boundaryMesh());

        if
        (
            (
                controller_.algorithm() == meshControl::STANDARD
                || controller_.algorithm() == meshControl::SHELL
            )
            && polyPatches.findPatchID(name) != -1
        )
        {
            WarningInFunction
                << "Trying to add a geometry patch : " << name
                << " that is already existing as a blockMesh patch"
                << endl;
        }
        const label patchi = addPatch(mesh_, name, patchInfo);

//        dictionary patchDict(patchInfo);
//        patchDict.set("nFaces", 0);
//        patchDict.set("startFace", 0);
//        autoPtr<polyPatch> ppPtr
//        (
//            polyPatch::New
//            (
//                name,
//                patchDict,
//                0,
//                mesh_.boundaryMesh()
//            )
//        );
//        label patchi = fvMeshTools::addPatch
//        (
//            mesh_,
//            ppPtr(),
//            dictionary(),       // optional field values
//            calculatedFvPatchField<scalar>::typeName,
//            true
//        );

        // Store
        meshedPatches_.append(name);

        // Clear patch based addressing
        faceToCoupledPatch_.clear();

        return patchi;
    }
}


void Foam::meshRefinement::removeZeroSizedPatches(fvMesh& mesh)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    label patchI = 0;

    labelList oldToNew(patches.size());
    forAll(patches, i)
    {
        const polyPatch& pp = patches[i];
        if (!isA<processorPolyPatch>(pp))
        {
            label sz = returnReduce(pp.size(), sumOp<label>());
            if (sz != 0)
            {
                oldToNew[i] = patchI;
                patchI++;
            }
        }
        else
        {
            if (pp.size())
            {
                oldToNew[i] = patchI;
                patchI++;
            }
        }
    }

    label nKeep = patchI;

    forAll(patches, i)
    {
        const polyPatch& pp = patches[i];
        if (!isA<processorPolyPatch>(pp))
        {
            label sz = returnReduce(pp.size(), sumOp<label>());
            if (sz == 0)
            {
                oldToNew[i] = patchI;
                patchI++;
            }
        }
        else
        {
            if (pp.size() == 0)
            {
                oldToNew[i] = patchI;
                patchI++;
            }
        }
    }

    fvMeshTools::reorderPatches(mesh, oldToNew, nKeep, true);

    return;
}


Foam::label Foam::meshRefinement::addFaceZone
(
    const word& fzName,
    const word& masterPatch,
    const word& slavePatch,
    const surfaceZonesInfo::faceZoneType& fzType
)
{
    label zonei = surfaceZonesInfo::addFaceZone
    (
        fzName,   //name
        labelList(0),   //addressing
        boolList(0),    //flipmap
        mesh_
    );

    faceZoneToMasterPatch_.insert(fzName, masterPatch);
    faceZoneToSlavePatch_.insert(fzName, slavePatch);
    faceZoneToType_.insert(fzName, fzType);

    return zonei;
}


bool Foam::meshRefinement::getFaceZoneInfo
(
    const word& fzName,
    label& masterPatchID,
    label& slavePatchID,
    surfaceZonesInfo::faceZoneType& fzType
) const
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    if (!faceZoneToMasterPatch_.found(fzName))
    {
        return false;
    }
    else
    {
        const word& masterName = faceZoneToMasterPatch_[fzName];
        masterPatchID = pbm.findPatchID(masterName);

        const word& slaveName = faceZoneToSlavePatch_[fzName];
        slavePatchID = pbm.findPatchID(slaveName);

        fzType = faceZoneToType_[fzName];

        return true;
    }
}


void Foam::meshRefinement::selectSeparatedCoupledFaces(boolList& selected) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchi)
    {
        // Check all coupled. Avoid using .coupled() so we also pick up AMi.
        if (isA<coupledPolyPatch>(patches[patchi]))
        {
            const coupledPolyPatch& cpp = refCast<const coupledPolyPatch>
            (
                patches[patchi]
            );

            if (cpp.transform().transformsPosition())
            {
                forAll(cpp, i)
                {
                    selected[cpp.start()+i] = true;
                }
            }
        }
    }
}


Foam::labelList Foam::meshRefinement::meshedPatches() const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    DynamicList<label> patchIDs(meshedPatches_.size());
    forAll(meshedPatches_, i)
    {
        label patchI = patches.findPatchID(meshedPatches_[i]);

        if (patchI == -1)
        {
            WarningInFunction
                << "Problem : did not find patch " << meshedPatches_[i]
                << endl << "Valid patches are " << patches.names()
                << endl;    //abort(FatalError);
            continue;
        }
        if (!polyPatch::constraintType(patches[patchI].type()))
        {
            patchIDs.append(patchI);
        }
        else if (controller_.dual())
        {
            WarningInFunction
                << "Block mesh type " << patches[patchI].type()
                << " for patch " << patches[patchI].name()
                << " is not supported by dual mesh generator"
                << endl;
        }
    }

    return labelList(patchIDs, true);
}


Foam::labelList Foam::meshRefinement::unmeshedPatches() const
{
    DynamicList<label> patchIDs(mesh_.boundaryMesh().size());
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    labelHashSet mPatches(meshedPatches());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if
        (
            !isA<processorPolyPatch>(pp)
            && !mPatches.found(patchI)
        )
        {
            patchIDs.append(patchI);
        }
    }

    return labelList(patchIDs, true);
}


Foam::label Foam::meshRefinement::findRegion
(
    const polyMesh& mesh,
    const labelList& cellToRegion,
    const vector& perturbVec,
    const point& p
)
{
    label regioni = -1;

    // Force calculation of base points (needs to be synchronised)
    (void)mesh.tetBasePtIs();

    label celli = mesh.findCell(p);
    //if (celli != -1)
    //{
    //    Pout<< "findRegion : Found point:" << p << " in cell " << celli
    //        << " at:" << mesh.cellCentres()[celli] << endl;
    //}
    //else
    //{
    //    Pout<< "findRegion : Found point:" << p << " in cell " << celli
    //        << endl;
    //}

    if (celli != -1)
    {
        regioni = cellToRegion[celli];
    }
    reduce(regioni, maxOp<label>());

    if (regioni == -1)
    {
        // See if we can perturb a bit
        celli = mesh.findCell(p+perturbVec);
        if (celli != -1)
        {
            regioni = cellToRegion[celli];
        }
        reduce(regioni, maxOp<label>());
    }
    return regioni;
}
//XXXXXXXX
//Foam::labelList Foam::meshRefinement::findRegion
//(
//    const polyMesh& mesh,
//    const labelList& cellToRegion,
//    const vector& perturbVec,
//    const pointField& pts
//)
//{
//    labelList regions(pts.size(), -1);
//
//    forAll(pts, i)
//    {
//        label celli = mesh.findCell(pts[i]);
//        if (celli != -1)
//        {
//            regions[i] = cellToRegion[celli];
//        }
//        reduce(regions[i], maxOp<label>());
//
//        if (regions[i] == -1)
//        {
//            // See if we can perturb a bit
//            celli = mesh.findCell(pts[i]+perturbVec);
//            if (celli != -1)
//            {
//                regions[i] = cellToRegion[celli];
//            }
//            reduce(regions[i], maxOp<label>());
//        }
//    }
//
//    forAll(regions, i)
//    {
//        if (regions[i] == -1)
//        {
//            FatalErrorInFunction
//                << "Point " << pts[i]
//                << " is not inside the mesh." << nl
//                << "Bounding box of the mesh:" << mesh.bounds()
//                //<< "All points " << pts
//                //<< " with corresponding regions " << regions
//                << exit(FatalError);
//        }
//    }
//
//    return regions;
//}
//XXXXXXXX

Foam::label Foam::meshRefinement::findCell
(
    const point& keepPoint,
    const polyMesh& mesh,
    const hexRef8& meshCutter
)
{
    label cellI = mesh.findCell(keepPoint, polyMesh::FACE_PLANES);

    label gCellI = cellI;
    reduce(gCellI, maxOp<label>());
    if (gCellI == -1)
    {
        Info<<"mesh.findCell not working for keepPoint: "<< keepPoint
            <<" trying alternative method"
            <<endl;

        label nearestCellI = mesh.findNearestCell(keepPoint);
        List<scalar> projDist(Pstream::nProcs(), GREAT);

        if (nearestCellI != -1)
        {
            const scalar edge0Len = meshCutter.level0EdgeLength();
            const vectorField& cf = mesh.faceCentres();
            const vectorField& Sf = mesh.faceAreas();

            const labelList& f = mesh.cells()[nearestCellI];
            projDist[Pstream::myProcNo()] = -GREAT;
            scalar& dist = projDist[Pstream::myProcNo()];

            forAll(f, facei)
            {
                label nFace = f[facei];

                vector proj = keepPoint - cf[nFace];
                vector normal = Sf[nFace];
                if (mesh.faceOwner()[nFace] != nearestCellI)
                {
                    normal = -normal;
                }
                dist = max(dist, (normal & proj));
            }
            dist /= (edge0Len/(1<<meshCutter.cellLevel()[nearestCellI]));
        }
        reduce(projDist, minOp<List<scalar>>());
        label minProcID = -1;
        scalar minDist = GREAT;

        forAll(projDist, i)
        {
            if (projDist[i] < minDist)
            {
                minDist = projDist[i];
                minProcID = i;
            }
        }
        if (Pstream::myProcNo() == minProcID && minDist < 0.05)
        {
            cellI = nearestCellI;
        }
    }

    gCellI = cellI;
    reduce(gCellI, maxOp<label>());
    if (gCellI == -1)
    {
        WarningInFunction
            << "Point " << keepPoint
            << " is not inside the mesh." << nl
            << "Bounding box of the mesh:" << mesh.bounds()
            << endl;
    }

    return cellI;
}


// Modify cellRegion to be consistent with locationsInMesh.
// - all regions not in locationsInMesh are set to -1
// - check that all regions inside locationsOutsideMesh are set to -1
void Foam::meshRefinement::findRegions
(
    const polyMesh& mesh,
    const hexRef8& meshCutter,
    const vector& perturbVec,
    const pointField& locationsInMesh,
    const pointField& locationsOutsideMesh,
    const writer<scalar>& leakPathFormatter,
    const bool& calcLeakPath,
    const label nRegions,
    labelList& cellRegion,
    const boolList& blockedFace
)
{
    PackedBoolList insideCell(mesh.nCells());


    // Mark all cells reachable from locationsInMesh
    labelList insideRegions(locationsInMesh.size());
    forAll(insideRegions, i)
    {
        // Find the region containing the keepPoint
        label regioni = -1;

        label cellI = findCell
        (
            locationsInMesh[i],
            mesh,
            meshCutter
        );

        if (cellI != -1)
        {
            regioni = cellRegion[cellI];
        }
        reduce(regioni, maxOp<label>());

/*
        // Find the region containing the point
        label regioni = findRegion
        (
            mesh,
            cellRegion,
            perturbVec,
            locationsInMesh[i]
        );
*/
        insideRegions[i] = regioni;

        // Mark all cells in the region as being inside
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == regioni)
            {
                insideCell[celli] = true;
            }
        }
    }

    // Check that all the locations outside the
    // mesh do not conflict with those inside
    DynamicList<labelPair> leakLocations(locationsOutsideMesh.size());
    forAll(locationsOutsideMesh, i)
    {
        // Find the region containing the point
        label regioni = findRegion
        (
            mesh,
            cellRegion,
            perturbVec,
            locationsOutsideMesh[i]
        );

        if (calcLeakPath && regioni != -1)
        {
            // Do a quick check for locationsOutsideMesh overlapping with
            // inside ones.
            label index = findIndex(insideRegions, regioni);
            word outTag = "outsideLocation" + name(i);
            if (index != -1)
            {
                calculateLeakPath
                (
                    "leakPath",
                    outTag,
                    mesh,
                    blockedFace,
                    pointField(1,locationsInMesh[index]),
                    pointField(1,locationsOutsideMesh[i]),
                    leakPathFormatter
                );
                leakLocations.append(labelPair(i,index));
            }
        }
        else if (regioni != -1)
        {
            label index = findIndex(insideRegions, regioni);
            WarningInFunction
                << "Location in mesh " << locationsInMesh[index]
                << " is inside same mesh region " << regioni
                << " as one of the locations outside mesh "
                << locationsOutsideMesh
                << endl;
        }
    }

    if (leakLocations.size())
    {
        fileName outputDir =
            mesh.time().rootPath()
            / mesh.time().globalCaseName()
            / "postProcessing"
            / mesh.pointsInstance();

        Info<<"Detected leak path from locationsInMesh "
            <<" to locationsOutsideMesh for the following locations "
            << endl;

        forAll(leakLocations, leaki)
        {
            label outLeak  = leakLocations[leaki].first();
            label inLeak  = leakLocations[leaki].second();
            Info<<"Leak detected from location : "
                << locationsOutsideMesh[outLeak]
                <<" to : " <<locationsInMesh[inLeak] <<endl;
        }

        Info<<"Leaks paths written to "<<outputDir<<endl;

        FatalErrorInFunction
            << "Exiting after leak detected. "
            << exit(FatalError);
    }

    // Now update cellRegion to -1 for unreachable cells
    forAll(insideCell, celli)
    {
        if (!insideCell[celli])
        {
            cellRegion[celli] = -1;
        }
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::splitMeshRegions
(
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const pointField& locationsInMesh,
    const pointField& locationsOutsideMesh,
    const writer<scalar>& leakPathFormatter,
    const bool& calcLeakPath
)
{
    // Force calculation of face decomposition (used in findCell)
    (void)mesh_.tetBasePtIs();

    // Determine connected regions. regionSplit is the labelList with the
    // region per cell.

    boolList blockedFace(mesh_.nFaces(), false);
    selectSeparatedCoupledFaces(blockedFace);

    regionSplit cellRegion(mesh_, blockedFace);

    findRegions
    (
        mesh_,
        meshCutter_,
        mergeDistance_*vector(1,1,1),   // perturbVec
        locationsInMesh,
        locationsOutsideMesh,
        leakPathFormatter,
        calcLeakPath,
        cellRegion.nRegions(),
        cellRegion,
        blockedFace
    );

    // Subset
    // ~~~~~~

    // Get cells to remove
    DynamicList<label> cellsToRemove(mesh_.nCells());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] == -1)
        {
            cellsToRemove.append(celli);
        }
    }
    cellsToRemove.shrink();

    label nTotCellsToRemove = returnReduce
    (
        cellsToRemove.size(),
        sumOp<label>()
    );


    autoPtr<mapPolyMesh> mapPtr;
    if (nTotCellsToRemove > 0)
    {
        label nCellsToKeep =
            mesh_.globalData().nTotalCells()
          - nTotCellsToRemove;

        Info<< "Keeping all cells containing points " << locationsInMesh << endl
            << "Selected for keeping : "
            << nCellsToKeep
            << " cells." << endl;


        // Remove cells
        removeCells cellRemover(mesh_);

        labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
        labelList exposedPatch;

        label nExposedFaces = returnReduce(exposedFaces.size(), sumOp<label>());
        if (nExposedFaces)
        {
            // FatalErrorInFunction
            //    << "Removing non-reachable cells should only expose"
            //    << " boundary faces" << nl
            //    << "ExposedFaces:" << exposedFaces << abort(FatalError);

            // Patch for exposed faces for lack of anything sensible.
            label defaultPatch = 0;
            if (globalToMasterPatch.size())
            {
                defaultPatch = globalToMasterPatch[0];
            }

            WarningInFunction
                << "Removing non-reachable cells exposes "
                << nExposedFaces << " internal or coupled faces." << endl
                << "    These get put into patch " << defaultPatch << endl;
            exposedPatch.setSize(exposedFaces.size(), defaultPatch);
        }

        mapPtr = doRemoveCells
        (
            cellsToRemove,
            exposedFaces,
            exposedPatch,
            cellRemover
        );
    }
    return mapPtr;
}


void Foam::meshRefinement::distribute
(
   const mapDistributePolyMesh& map,
   const bool& updateSurf
)
{
    // mesh_ already distributed; distribute my member data
    // refinement
    meshCutter_.distribute(map);

    // surfaceIndex is face data.
    map.distributeFaceData(surfaceIndex_);

    // wrapIndex is face data.
    map.distributeFaceData(wrapIndex_);

    // faceToPatch (baffles that were on coupled faces) is not maintained
    // (since baffling also disconnects points)
    faceToCoupledPatch_.clear();

    // Update cellRegionID
    map.distributeCellData(cellRegionID_);

    // maintainedFaces are indices of faces.
    forAll(userFaceData_, i)
    {
        map.distributeFaceData(userFaceData_[i].second());
    }

    // Redistribute surface and any fields on it.
    if (updateSurf)
    {
        Random rndGen(653213);

        {
            // Get local mesh bounding box. Single box for now.
            List<treeBoundBox> meshBb(1);
            treeBoundBox& bb = meshBb[0];

            point minBb = point::max;
            point maxBb = -point::max;

            if (mesh_.points().empty())
            {
                minBb = point::zero;
                maxBb = point::zero;
            }

            forAll(mesh_.points(), pointI)
            {
                minBb = Foam::min(minBb, mesh_.points()[pointI]);
                maxBb = Foam::max(maxBb, mesh_.points()[pointI]);
            }

            const pointField& cellCentres = mesh_.cellCentres();
            const pointField& faceCentres = mesh_.faceCentres();
            for
            (
                label faceI = mesh_.nInternalFaces();
                faceI < mesh_.nFaces();
                faceI++
            )
            {
                label own = mesh_.faceOwner()[faceI];
                const point& fc = faceCentres[faceI];
                const point& cc = cellCentres[own];
                point external = 2*fc - cc;
                minBb = Foam::min(minBb, external);
                maxBb = Foam::max(maxBb, external);
            }
            bb = treeBoundBox(minBb, maxBb);
            bb = bb.extend(rndGen, 1E-4);

            // Distribute all geometry (so refinementSurfaces and
            // shellSurfaces)
            searchableSurfaces& geometry =
                const_cast<searchableSurfaces&>(surfaces_.geometry());

            forAll(geometry, i)
            {
                autoPtr<mapDistribute> faceMap;
                autoPtr<mapDistribute> pointMap;
                geometry[i].distribute
                (
                    meshBb,
                    false,          // do not keep outside triangles
                    faceMap,
                    pointMap
                );

                if (faceMap.valid())
                {
                    // (ab)use the instance() to signal current modification time
                    geometry[i].instance() = geometry[i].time().timeName();
                }

                faceMap.clear();
                pointMap.clear();
            }
        }
/*
        // re-distribute feature lines
        {
            // compact BB for feature lines
            List<treeBoundBox> meshBb(1);
            meshBb[0] =
                treeBoundBox
                (
                    boundBox(mesh_.points(), false)
                 ).extend(rndGen, 1E-3);

            forAll(featuresPtr_(), featI)
            {
                edgeMesh& eMesh = featuresPtr_()[featI];;

                eMesh.distribute
                (
                    meshBb,
                    mergeDistance_
                 );
            }
            //clear feature line search trees
            featuresPtr_().clearTrees();
        }
*/
    }
}


void Foam::meshRefinement::updateMesh
(
    const mapPolyMesh& map,
    const labelList& changedFaces,
    const bool updateIntersections
)
{
    Map<label> dummyMap(0);

    updateMesh
    (
        map,
        changedFaces,
        dummyMap,
        dummyMap,
        dummyMap,
        updateIntersections
    );
}


void Foam::meshRefinement::storeData
(
    const labelList& pointsToStore,
    const labelList& facesToStore,
    const labelList& cellsToStore
)
{
    // For now only meshCutter has storable/retrievable data.
    meshCutter_.storeData
    (
        pointsToStore,
        facesToStore,
        cellsToStore
    );
}

void Foam::meshRefinement::checkMeshSize()
{
    //Check size of mesh for 32 bit label issue
    if constexpr (sizeof(Foam::label) == 4)
    {
        label globalNumCells = mesh_.nCells();
        label globalNumEdges = mesh_.nEdges();
        label globalNumPoints = mesh_.nPoints();
        label globalNumFaces = mesh_.nFaces();

        reduce(
            std::tie(globalNumCells, globalNumEdges, globalNumPoints, globalNumFaces),
            UniformParallelOp<sumOp<label>, 4>{}
        );

        label warnSz = 1500000000;
        if
        (
            globalNumCells > warnSz
            || globalNumEdges > warnSz
            || globalNumPoints > warnSz
            || globalNumFaces > warnSz
        )
        {
            WarningInFunction
                << "Mesh size getting close to 32 bit label size limit"
                << " Number Cells : " << globalNumCells
                << " Number Edges : " << globalNumEdges
                << " Number Points : " << globalNumPoints
                << " Number Faces : " << globalNumFaces
                << " Max label size 2e9. Try switching to 64 bit label size."
                << endl;
        }
    }
}


void Foam::meshRefinement::updateMesh
(
    const mapPolyMesh& map,
    const labelList& changedFaces,
    const Map<label>& pointsToRestore,
    const Map<label>& facesToRestore,
    const Map<label>& cellsToRestore,
    const bool update
)
{
    checkMeshSize();

    // For now only meshCutter has storable/retrievable data.

    // Update numbering of cells/vertices.
    meshCutter_.updateMesh
    (
        map,
        pointsToRestore,
        facesToRestore,
        cellsToRestore
    );

    // Update cell region ID
    {
        // Map data
        const labelList& cellMap = map.cellMap();

        labelList newCellRegion(cellMap.size());
        forAll(cellMap, newCelli)
        {
            label oldCelli = cellMap[newCelli];
            if (oldCelli == -1)
            {
                newCellRegion[newCelli] = -1;
            }
            else
            {
                newCellRegion[newCelli] = cellRegionID_[oldCelli];
            }
        }
        cellRegionID_.transfer(newCellRegion);
    }

    if (oldPointsPtr_)
    {
        // Map data
        const labelList& pointMap = map.pointMap();

        pointField& oldPts = oldPoints();
        pointField newPoints(pointMap.size());

        forAll(pointMap, newPointI)
        {
            label oldPointI = pointMap[newPointI];

            if (oldPointI == -1)
            {
                newPoints[newPointI] = vector::zero;
            }
            else
            {
                newPoints[newPointI] = oldPts[oldPointI];
            }
        }
        oldPts.transfer(newPoints);
    }

    // Update surfaceIndex
    updateList(map.faceMap(), label(-1), surfaceIndex_);

    // Update wrap index
    updateList(map.faceMap(), label(-1), wrapIndex_);

    // Update faceToCoupledPatch_
    {
        Map<label> newFaceToPatch(faceToCoupledPatch_.size());
        forAllConstIter(Map<label>, faceToCoupledPatch_, iter)
        {
            label newFacei = map.reverseFaceMap()[iter.key()];

            if (newFacei >= 0)
            {
                newFaceToPatch.insert(newFacei, iter());
            }
        }
        faceToCoupledPatch_.transfer(newFaceToPatch);
    }

    if (update && returnReduce(changedFaces.size(), sumOp<label>()))
    {
        // Update cached intersection information
        updateIntersections(changedFaces);
    }

    // Update maintained faces
    forAll(userFaceData_, i)
    {
        labelList& data = userFaceData_[i].second();

        if (userFaceData_[i].first() == KEEPALL)
        {
            // extend list with face-from-face data
            updateList(map.faceMap(), label(-1), data);
        }
        else if (userFaceData_[i].first() == MASTERONLY)
        {
            // keep master only
            labelList newFaceData(map.faceMap().size(), -1);

            forAll(newFaceData, facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0 && map.reverseFaceMap()[oldFacei] == facei)
                {
                    newFaceData[facei] = data[oldFacei];
                }
            }
            data.transfer(newFaceData);
        }
        else
        {
            // remove any face that has been refined i.e. referenced more than
            // once.

            // 1. Determine all old faces that get referenced more than once.
            // These get marked with -1 in reverseFaceMap

            labelList reverseFaceMap(map.reverseFaceMap());

            forAll(map.faceMap(), facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0)
                {
                    if (reverseFaceMap[oldFacei] != facei)
                    {
                        // facei is slave face. Mark old face.
                        reverseFaceMap[oldFacei] = -1;
                    }
                }
            }

            // 2. Map only faces with intact reverseFaceMap
            labelList newFaceData(map.faceMap().size(), -1);
            forAll(newFaceData, facei)
            {
                label oldFacei = map.faceMap()[facei];

                if (oldFacei >= 0)
                {
                    if (reverseFaceMap[oldFacei] == facei)
                    {
                        newFaceData[facei] = data[oldFacei];
                    }
                }
            }
            data.transfer(newFaceData);
        }
    }

    //Updated manifold features if required for tracking
    if (features().valid())
    {
        features()().manFeatures().update(map);
    }
}


bool Foam::meshRefinement::write(bool writeSurf) const
{
    bool writeOk = mesh_.write()  && cellRegionID_.write();

    if (debug)
    {
        writeOk = meshCutter_.write()
            && surfaceIndex_.write();
    }

    if (writeSurf)
    {
       // Make sure that any distributed surfaces (so ones which probably have
       // been changed) get written as well.
       // Note: should ideally have some 'modified' flag to say whether it
       // has been changed or not.
       searchableSurfaces& geometry =
          const_cast<searchableSurfaces&>(surfaces_.geometry());

       forAll(geometry, i)
       {
          searchableSurface& s = geometry[i];

          // Check if instance() of surface is not constant or system.
          // Is good hint that surface is distributed.
          if
          (
             s.instance() != s.time().system()
          && s.instance() != s.time().caseSystem()
          && s.instance() != s.time().constant()
          && s.instance() != s.time().caseConstant()
          )
          {
             // Make sure it gets written to current time, not constant.
             s.instance() = s.time().timeName();
             writeOk = writeOk && s.write();
          }
       }
    }

    return writeOk;
}


Foam::PackedBoolList Foam::meshRefinement::getMasterPoints
(
    const polyMesh& mesh,
    const labelList& meshPoints
)
{
    const globalIndex globalPoints(meshPoints.size());

    labelList myPoints(meshPoints.size());
    forAll(meshPoints, pointi)
    {
        myPoints[pointi] = globalPoints.toGlobal(pointi);
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        myPoints,
        minEqOp<label>(),
        labelMax
    );


    PackedBoolList isPatchMasterPoint(meshPoints.size());
    forAll(meshPoints, pointi)
    {
        if (myPoints[pointi] == globalPoints.toGlobal(pointi))
        {
            isPatchMasterPoint[pointi] = true;
        }
    }

    return isPatchMasterPoint;
}


Foam::PackedBoolList Foam::meshRefinement::getMasterEdges
(
    const polyMesh& mesh,
    const labelList& meshEdges
)
{
    const globalIndex globalEdges(meshEdges.size());

    labelList myEdges(meshEdges.size());
    forAll(meshEdges, edgei)
    {
        myEdges[edgei] = globalEdges.toGlobal(edgei);
    }

    syncTools::syncEdgeList
    (
        mesh,
        meshEdges,
        myEdges,
        minEqOp<label>(),
        labelMax
    );


    PackedBoolList isMasterEdge(meshEdges.size());
    forAll(meshEdges, edgei)
    {
        if (myEdges[edgei] == globalEdges.toGlobal(edgei))
        {
            isMasterEdge[edgei] = true;
        }
    }

    return isMasterEdge;
}


void Foam::meshRefinement::printMeshInfo(const bool debug, const string& msg)
const
{
    const globalMeshData& pData = mesh_.globalData();

    if (debug)
    {
        Pout<< msg.c_str()
            << " : cells(local):" << mesh_.nCells()
            << "  faces(local):" << mesh_.nFaces()
            << "  points(local):" << mesh_.nPoints()
            << endl;
    }

    {
        PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));
        label nMasterFaces = 0;
        forAll(isMasterFace, i)
        {
            if (isMasterFace[i])
            {
                nMasterFaces++;
            }
        }

        PackedBoolList isMeshMasterPoint(syncTools::getMasterPoints(mesh_));
        label nMasterPoints = 0;
        forAll(isMeshMasterPoint, i)
        {
            if (isMeshMasterPoint[i])
            {
                nMasterPoints++;
            }
        }

        Info<< msg.c_str()
            << " : cells:" << pData.nTotalCells()
            << "  faces:" << returnReduce(nMasterFaces, sumOp<label>())
            << "  points:" << returnReduce(nMasterPoints, sumOp<label>())
            << endl;
    }


    //if (debug)
    {
        const labelList& cellLevel = meshCutter_.cellLevel();

        label nLevels = gMax(cellLevel);
        if (nLevels > -1)
        {
            labelList nCells(nLevels+1, 0);

            forAll(cellLevel, celli)
            {
                nCells[cellLevel[celli]]++;
            }

            Pstream::listCombineReduce(nCells, plusOp<label>());

            Info<< "Cells per refinement level:" << endl;
            forAll(nCells, leveli)
            {
                Info<< "    " << leveli << '\t' << nCells[leveli]
                    << endl;
            }
        }
    }
}

//- Return either time().constant() or oldInstance
Foam::word Foam::meshRefinement::timeName() const
{
    if (overwrite_ && mesh_.time().timeIndex() == 0)
    {
        return oldInstance_;
    }
    else
    {
        return mesh_.time().timeName();
    }
}


void Foam::meshRefinement::dumpRefinementLevel() const
{
    // Note: use time().timeName(), not meshRefinement::timeName()
    // so as to dump the fields to 0, not to constant.
    {
        volScalarField volRefLevel
        (
            IOobject
            (
                "cellLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(volRefLevel, celli)
        {
            volRefLevel[celli] = cellLevel[celli];
        }

        volRefLevel.write();
    }

    // Dump pointLevel
    {
        const pointMesh& pMesh = pointMesh::New(mesh_);

        pointScalarField pointRefLevel
        (
            IOobject
            (
                "pointLevel",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pMesh,
            dimensionedScalar("zero", dimless, 0)
        );

        const labelList& pointLevel = meshCutter_.pointLevel();

        forAll(pointRefLevel, pointi)
        {
            pointRefLevel[pointi] = pointLevel[pointi];
        }

        pointRefLevel.write();
    }
}


// Dump cell centres
void Foam::meshRefinement::dumpIntersections(const fileName& prefix) const
{
    {
        OFstream str(prefix + "_edges.obj");
        label verti = 0;
        Pout<< "meshRefinement::dumpIntersections :"
            << " Writing cellcentre-cellcentre intersections to file "
            << str.name() << endl;


        // Redo all intersections
        // ~~~~~~~~~~~~~~~~~~~~~~

        // Get boundary face centre and level. Coupled aware.
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
        calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc);

        labelList intersectionFaces(intersectedFaces());

        // Collect segments we want to test for
        pointField start(intersectionFaces.size());
        pointField end(intersectionFaces.size());
        {
            labelList minLevel;
            calcCellCellRays
            (
                mesh_.cellCentres(),
                neiCc,
                labelList(neiCc.size(), -1),
                intersectionFaces,
                start,
                end,
                minLevel
            );
        }


        // Do tests in one go
        labelList surfaceHit;
        List<pointIndexHit> surfaceHitInfo;
        surfaces_.findAnyIntersection
        (
            start,
            end,
            surfaceHit,
            surfaceHitInfo
        );

        forAll(intersectionFaces, i)
        {
            if (surfaceHit[i] != -1)
            {
                meshTools::writeOBJ(str, start[i]);
                verti++;
                meshTools::writeOBJ(str, surfaceHitInfo[i].hitPoint());
                verti++;
                meshTools::writeOBJ(str, end[i]);
                verti++;
                str << "l " << verti-2 << ' ' << verti-1 << nl
                    << "l " << verti-1 << ' ' << verti << nl;
            }
        }
    }

    Pout<< endl;
}


void Foam::meshRefinement::write
(
    const debugType debugFlags,
    const writeType writeFlags,
    const fileName& prefix,
    bool writeSurf
) const
{
    if (writeFlags & WRITEMESH)
    {
        write(writeSurf);
    }

    if (writeFlags && !(writeFlags & NOWRITEREFINEMENT))
    {
        meshCutter_.write();
        surfaceIndex_.write();
    }

    if (writeFlags & WRITELEVELS)
    {
        dumpRefinementLevel();
    }

    if (debugFlags & OBJINTERSECTIONS && prefix.size())
    {
        dumpIntersections(prefix);
    }
}


void Foam::meshRefinement::removeFiles(const polyMesh& mesh)
{
    IOobject io
    (
        "dummy",
        mesh.facesInstance(),
        mesh.meshSubDir,
        mesh
    );
    fileName setsDir(io.path());

    if (topoSet::debug) DebugVar(setsDir);

    // Remove local files
    if (exists(setsDir/"surfaceIndex"))
    {
        rm(setsDir/"surfaceIndex");
    }

    // Remove other files
    hexRef8::removeFiles(mesh);
}

Foam::meshRefinement::writeType Foam::meshRefinement::writeLevel()
{
    return writeLevel_;
}


void Foam::meshRefinement::writeLevel(const writeType flags)
{
    writeLevel_ = flags;
}


Foam::meshRefinement::outputType Foam::meshRefinement::outputLevel()
{
    return outputLevel_;
}


void Foam::meshRefinement::outputLevel(const outputType flags)
{
    outputLevel_ = flags;
}


void Foam::meshRefinement::surfaceToPatch
(
    const dictionary& surfaceToPatchDict
)
{
    Info<<"Repatching faces based on input surfaces"<<endl;

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // all surface geometry
    const dictionary& geometryDict = surfaceToPatchDict.subDict("geometry");

    autoPtr<searchableSurfaces> geometryPtr
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                      // dummy name
                mesh_.time().constant(),     // directory
                "triSurface",               // instance
                mesh_.time(),                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            geometryDict,
            true,
            List<fileName>(0),
            labelMax
        )
    );

    const searchableSurfaces& geom = geometryPtr();

    labelList faceIndex(mesh_.nBoundaryFaces());
    pointField start(mesh_.nBoundaryFaces());
    pointField end(mesh_.nBoundaryFaces());

    label nChecked = 0;

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        if (!pp.coupled())
        {
            label startFace = pp.start();
            forAll(pp, i)
            {
                label facei = startFace + i;
                label own = mesh_.faceOwner()[facei];
                point cc = mesh_.cellCentres()[own];
                point fc = mesh_.faceCentres()[facei];
                point fA = mesh_.faceAreas()[facei];
                scalar area = mag(fA);
                fA /= (area + SMALL);

                scalar fW = 0;
                scalar cW = 0;
                {
                    face f = mesh_.faces()[facei];
                    scalarField fProj(f.size());
                    forAll(f, fp)
                    {
                        vector n = mesh_.points()[f[fp]] - fc;
                        fProj[fp] = (n & fA);
                    }

                    // Get normal 'span' of face
                    scalar minVal = min(fProj);
                    scalar maxVal = max(fProj);
                    fW = (maxVal - minVal);
                }

                {
                    const labelList& cPts = mesh_.cellPoints()[own];
                    scalarField cProj(cPts.size());
                    forAll(cPts, cPtI)
                    {
                        vector n = mesh_.points()[cPts[cPtI]] - cc;
                        cProj[cPtI] = (n & fA);
                    }

                    // Get normal 'span' of face
                    scalar minVal = min(cProj);
                    scalar maxVal = max(cProj);
                    cW = (maxVal - minVal);
                }

                scalar fcToCC = mag((cc - fc)&fA);
                scalar projDistFC = max(fW, fcToCC);
                scalar projDistCC = max(0.5*cW,fcToCC);

                start[nChecked] = fc - projDistCC*fA;
                end[nChecked] = fc + projDistFC*fA;
                faceIndex[nChecked] = facei;
                nChecked++;
            }
        }
    }
    faceIndex.setSize(nChecked);
    start.setSize(nChecked);
    end.setSize(nChecked);

    List<pointIndexHit> info;
    labelList nearestSurfaces;
    geom.findAnyIntersection
    (
        start,
        end,
        nearestSurfaces,
        info
    );

    labelList surfaceRegion(info.size(), -1);
    vectorField surfaceNormal(info.size(), vector::zero);
    forAll(geom, geomi)
    {
        labelList regionHit(info.size(),-1);
        DynamicList<pointIndexHit> localHits;
        forAll(info, hiti)
        {
            if (nearestSurfaces[hiti] == geomi)
            {
                localHits.append(info[hiti]);
            }
        }
        labelList localRegions;
        geom[geomi].getRegion(localHits, localRegions);
        vectorField localNormals;
        geom[geomi].getNormal(localHits, localNormals);
        label nSurfaceHits = 0;
        forAll(info, hiti)
        {
            if (nearestSurfaces[hiti] == geomi)
            {
                surfaceRegion[hiti] = localRegions[nSurfaceHits];
                surfaceNormal[hiti] = localNormals[nSurfaceHits];
                nSurfaceHits++;
            }
        }
    }

    const List<wordList>& allNames = geom.regionNames();
    List<labelList> regionToPatch(allNames.size());
    forAll(geom, geomi)
    {
        const wordList& regNames = allNames[geomi];
        forAll(regNames, i)
        {
            regionToPatch[i].setSize(regNames.size());
            dictionary patchInfo;
            patchInfo.set("type", wallPolyPatch::typeName);

            regionToPatch[geomi][i] = addMeshedPatch
            (
                regNames[i],
                patchInfo
             );
        }
    }

    // Topology changes container
    polyTopoChange meshMod(mesh_);

    label nUpdated = 0;
    forAll(geom, geomi)
    {
        labelList regionHit(info.size(),-1);
        DynamicList<pointIndexHit> localHits;
        forAll(info, hiti)
        {
            if (nearestSurfaces[hiti] == geomi)
            {
                label facei = faceIndex[hiti];

                label patchi = patches.whichPatch(facei);
                label regioni = surfaceRegion[hiti];
                label closestPatchID = regionToPatch[geomi][regioni];
                vector sNorm = surfaceNormal[hiti];
                vector fNorm = mesh_.faceAreas()[facei];
                scalar fArea = mag(fNorm);
                fNorm /= (fArea + SMALL);
                scalar dProd(fNorm & sNorm);

                if (patchi != closestPatchID && mag(dProd) > 0.707)
                {
                    label own = mesh_.faceOwner()[facei];
                    label zoneID = mesh_.faceZones().whichZone(facei);
                    bool zoneFlip = false;
                    if (zoneID >= 0)
                    {
                        const faceZone& fZone = mesh_.faceZones()[zoneID];
                        zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
                    }

                    meshMod.modifyFace
                    (
                        mesh_.faces()[facei], // modified face
                        facei,          // label of face being modified
                        own,            // owner
                        -1,             // neighbour
                        false,          // face flip
                        closestPatchID,   // new patch for face
                        zoneID,         // zone for face
                        zoneFlip        // face flip in zone
                     );
                    nUpdated++;
                }
            }
        }
    }

    if (returnReduce(nUpdated, sumOp<label>()) != 0)
    {
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

        // Update fields
        mesh_.updateMesh(map);

        // Move mesh (since morphing might not do this)
        if (map().hasMotionPoints())
        {
            mesh_.movePoints(map().preMotionPoints());
        }
        else
        {
                // Delete mesh volumes. No other way to do this?
            mesh_.clearOut();
        }
        updateMesh(map, labelList(0));
    }
}


void Foam::meshRefinement::repatchRegion
(
    const dictionary& repatchDict
)
{
    Info<<"Repatching region faces"<<endl;
    forAllConstIter(dictionary, repatchDict, iter)
    {
        const word& key = iter().keyword();

        dictionary patchInfo;
        patchInfo.set("type", wallPolyPatch::typeName);
        addPatch
        (
            mesh_,
            key,
            patchInfo
        );
    }

    const regionSplit cellRegion(mesh_);
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    polyTopoChange meshMod(mesh_);

    boolList flaggedCells(mesh_.nCells(), false);
    boolList flaggedFaces(mesh_.nFaces(), false);

    forAllConstIter(dictionary, repatchDict, iter)
    {
        const word& key = iter().keyword();
        const label newPatchI = patches.findPatchID(key);

        const dictionary& rdict = repatchDict.subDict(key);
        const word& zoneName =  rdict.lookup("zone");

        labelHashSet excludePatches = patches.patchSet
        (
            wordReList(rdict.lookup("excludePatches")), false, true
        );

        cellZoneMesh& cellZones = mesh_.cellZones();
        label zoneI = cellZones.findZoneID(zoneName);

        if (zoneI == -1)
        {
            Info<<"Zone "<<zoneName<<" not found in mesh"
                <<" using location to zone"<<endl;

            if (rdict.found("location"))
            {
                const point& location = rdict.lookup("location");
                // Find the region containing the insidePoint
                label keepRegionI = -1;

                label cellI = mesh_.findCell(location, polyMesh::FACE_PLANES);

                if (cellI != -1)
                {
                    keepRegionI = cellRegion[cellI];
                }
                reduce(keepRegionI, maxOp<label>());

                if (keepRegionI == -1)
                {
                    WarningInFunction
                        << "Point " << location
                        << " cannot be found inside the mesh." << nl
                        << "Bounding box of the mesh:" << mesh_.bounds()
                        << endl;

                    continue;
                }
                else
                {
                    Info<<"Found region "<<keepRegionI<<" of "
                        <<cellRegion.nRegions()<<endl;
                }

                cellZoneMesh& cellZones = mesh_.cellZones();

               //create new zone
                zoneI = cellZones.size();
                cellZones.setSize(zoneI+1);
                cellZones.set
                (
                    zoneI,
                    new cellZone
                    (
                        zoneName,   //name
                        labelList(0),           //addressing
                        zoneI,                  //index
                        cellZones               //cellZoneMesh
                     )
                 );

                // Set all cells with this region
                forAll(cellRegion, cellI)
                {
                    bool set = false;
                    if (cellRegion[cellI] == keepRegionI)
                    {
                        if (!flaggedCells[cellI])
                        {
                            flaggedCells[cellI] = true;
                            set = true;

                            meshMod.modifyCell
                            (
                                cellI,
                                zoneI
                             );
                        }
                    }

                    if (set)
                    {
                        const cell& c = mesh_.cells()[cellI];
                        forAll(c, i)
                        {
                            label faceI = c[i];
                            label patchI = patches.whichPatch(faceI);
                            if
                            (
                                patchI != -1 && !flaggedFaces[faceI]
                                && !patches[patchI].coupled()
                                && !excludePatches.found(patchI)
                             )
                            {
                                flaggedFaces[faceI] = true;

                                label own = mesh_.faceOwner()[faceI];
                                label zoneID =
                                    mesh_.faceZones().whichZone(faceI);
                                bool zoneFlip = false;
                                if (zoneID >= 0)
                                {
                                    const faceZone& fZone =
                                        mesh_.faceZones()[zoneID];
                                    zoneFlip =
                                        fZone.flipMap()[fZone.whichFace(faceI)];
                                }

                                meshMod.modifyFace
                                (
                                    mesh_.faces()[faceI], // modified face
                                    faceI,     // label of face being modified
                                    own,            // owner
                                    -1,             // neighbour
                                    false,          // face flip
                                    newPatchI,      // new patch for face
                                    zoneID,         // zone for face
                                    zoneFlip        // face flip in zone
                                 );
                            }
                        }
                    }
                }
            }
            else
            {
                WarningInFunction
                    << "New zone: " << zoneName
                    << " has no location keyword set"
                    << endl;
                continue;
            }
        }
        else
        {
            // Set all cells with this region
            forAll(cellZones[zoneI], i)
            {
                label cellI = cellZones[zoneI][i];

                bool set = false;
                if (!flaggedCells[cellI])
                {
                    flaggedCells[cellI] = true;
                    set = true;
                }

                if (set)
                {
                    const cell& c = mesh_.cells()[cellI];
                    forAll(c, i)
                    {
                        label faceI = c[i];
                        label patchI = patches.whichPatch(faceI);
                        if
                        (
                            patchI != -1 && !flaggedFaces[faceI]
                            && !patches[patchI].coupled()
                            && !excludePatches.found(patchI)
                        )
                        {
                            flaggedFaces[faceI] = true;

                            label own = mesh_.faceOwner()[faceI];
                            label zoneID = mesh_.faceZones().whichZone(faceI);
                            bool zoneFlip = false;
                            if (zoneID >= 0)
                            {
                                const faceZone& fZone =
                                    mesh_.faceZones()[zoneID];
                                zoneFlip =
                                    fZone.flipMap()[fZone.whichFace(faceI)];
                            }

                            meshMod.modifyFace
                            (
                                mesh_.faces()[faceI], // modified face
                                faceI,          // label of face being modified
                                own,            // owner
                                -1,             // neighbour
                                false,          // face flip
                                newPatchI,      // new patch for face
                                zoneID,         // zone for face
                                zoneFlip        // face flip in zone
                             );
                        }
                    }
                }
            }
        }
    }


    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh_.clearOut();
    }
    updateMesh(map, labelList(0));

    return;
}


Foam::edgeMesh Foam::meshRefinement::calcBoundaryEdgeMesh
(
    const polyMesh& mesh,
    const labelHashSet& patchSet
)
{
    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
           mesh,
           patchSet.toc()
        )
    );
    const indirectPrimitivePatch& pp = ppPtr();

    if (returnReduce(pp.size(), sumOp<label>()) == 0)
    {
        WarningInFunction
           << "Patches : "<< patchSet.toc()
           << " contain zero faces." <<endl;
        return edgeMesh(pointField(0), edgeList(0));
    }

    const labelList meshEdges
    (
        pp.meshEdges
        (
            mesh.edges(),
            mesh.pointEdges()
        )
    );

    boolList excludedFaces(pp.size(), false);
    edgeClassification eClass
    (
        mesh,
        mesh.points(),
        pp,
        meshEdges,
        excludedFaces,
        0.707,
        0.707
    );

    const List<Tuple2<edgeClassification::edgeType,scalar>>&
        eType = eClass.edgeTypes();
    DynamicList<edge> localEdges(eType.size());
    DynamicList<label> edgeMeshPts(eType.size());
    const labelList& meshPoints = pp.meshPoints();
    labelList addedPoints(meshPoints.size(), -1);

    label nAddedPts = 0;
    forAll(eType, edgei)
    {
        if (eType[edgei].first() == edgeClassification::BOUNDARY)
        {
            edge le(-1,-1);
            const edge& e = pp.edges()[edgei];
            forAll(e,ei)
            {
                label ePt = e[ei];
                if (addedPoints[ePt] == -1)
                {
                    addedPoints[ePt] = nAddedPts;
                    le[ei] = nAddedPts;
                    label meshpointi = meshPoints[ePt];
                    edgeMeshPts.append(meshpointi);
                    nAddedPts++;
                }
                else
                {
                    le[ei] = addedPoints[ePt];
                }
            }
            localEdges.append(le);
        }
    }

    if (returnReduce(nAddedPts, sumOp<label>()) == 0)
    {
        WarningInFunction
            << "Patches : "<< patchSet.toc()
            << " are manifold and contain no boundary edges." <<endl;
        return edgeMesh(pointField(0), edgeList(0));
    }

    // Find correspondence to master points
    labelList pointToGlobal;
    labelList uniqueMeshPoints;
    autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
    (
        edgeMeshPts,
        pointToGlobal,
        uniqueMeshPoints
    );
    // Gather all unique points on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()] = pointField
    (
        mesh.points(),
        uniqueMeshPoints
    );
    Pstream::gatherList(gatheredPoints);

    // Gather all faces
    List<edgeList> gatheredEdges(Pstream::nProcs());
    gatheredEdges[Pstream::myProcNo()] = localEdges;
    forAll(gatheredEdges[Pstream::myProcNo()],i)
    {
        inplaceRenumber
        (
            pointToGlobal,
            gatheredEdges[Pstream::myProcNo()][i]
        );
    }
    Pstream::gatherList(gatheredEdges);

    // On master combine all points, faces, zones
    if (Pstream::master())
    {
        pointField allPoints = ListListOps::combine<pointField>
        (
            gatheredPoints,
            accessOp<pointField>()
        );
        gatheredPoints.clear();

        edgeList allEdges = ListListOps::combine<edgeList>
        (
            gatheredEdges,
            accessOp<edgeList>()
        );
        gatheredEdges.clear();

        return edgeMesh(allPoints,allEdges);
    }
    else
    {
        return edgeMesh(pointField(0), edgeList(0));
    }
}


void Foam::meshRefinement::calcSourceTargetMatch
(
    const polyMesh& mesh,
    const wordReList& sourcePatches,
    const wordReList& targetPatches,
    const bool writeVTK
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();
    labelHashSet sourceSet = pbm.patchSet
    (
        sourcePatches,
        false,
        true
    );
    edgeMesh sourceEMesh = calcBoundaryEdgeMesh(mesh,sourceSet);

    if (returnReduce(sourceEMesh.edges().size(), sumOp<label>()) == 0)
    {
        WarningInFunction
            << "Zero sized boundary edge for patches : "
            << sourceSet.toc() << endl;
        return;
    }

    labelHashSet targetSet = pbm.patchSet
    (
        targetPatches,
        false,
        true
    );
    edgeMesh targetEMesh = calcBoundaryEdgeMesh(mesh,targetSet);

    if (returnReduce(targetEMesh.edges().size(), sumOp<label>()) == 0)
    {
        WarningInFunction
            << "Zero sized boundary edge for patches : "
            << targetSet.toc() << endl;
        return;
    }

    if (Pstream::master())
    {
        const pointField& points = targetEMesh.points();
        const edgeList& edges = targetEMesh.edges();

        // Calculate bb of all points
        treeBoundBox bb(points);

        // Random number generator. Bit dodgy since not exactly random ;-)
        Random rndGen(65431);

        // Slightly extended bb. Slightly off-centred just so on symmetric
        // geometry there are less face/edge aligned items.
        bb = bb.extend(rndGen, 1e-4);
        bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        indexedOctree<treeDataEdge> targetEdgeTree
        (
            treeDataEdge
            (
                false,                  // do not cache bb
                edges,
                points,
                identity(edges.size())
            ),
            bb,     // overall search domain
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        );

        scalarField targetDist(sourceEMesh.points().size(),scalar(-1));
        scalar aveError = scalar(0);
        scalar maxError = -GREAT;
        point maxLocation = vector::zero;

        label nHits = 0;

        forAll(sourceEMesh.points(),pointi)
        {
           const point& sample = sourceEMesh.points()[pointi];
           pointIndexHit info = targetEdgeTree.findNearest(sample, GREAT);
           if (info.hit())
           {
               scalar hitDist = mag(info.hitPoint()-sample);
               targetDist[pointi] = hitDist;
               if (hitDist > maxError)
               {
                   maxError = hitDist;
                   maxLocation = sample;
               }
               aveError += hitDist;
               nHits++;
           }
        }

        if (nHits > 0)
        {
            Info<< nl
                 << "Calculating patch boundary matching..." << nl
                 << "Source patches : " << sourcePatches << nl
                 << "Target patches : " << targetPatches << nl
                 << "Average error : "<< (aveError/nHits) << nl
                 << "Maximum error : "<< maxError
                 << " at location : " << maxLocation << nl <<endl;
        }

        if (writeVTK)
        {
            simpleVTKWriter sourceEMeshVTK
            (
                sourceEMesh.edges(),
                sourceEMesh.points()
            );
            sourceEMeshVTK.addPointData("error",targetDist);
            word sourceVTKName =
               "sourceEdgeMesh" + sourcePatches[0] + ".vtk";
            sourceEMeshVTK.write(sourceVTKName,false);

            simpleVTKWriter targetEMeshVTK
            (
                targetEMesh.edges(),
                targetEMesh.points()
            );
            word targetVTKName =
               "targetEdgeMesh" + targetPatches[0] + ".vtk";
            targetEMeshVTK.write(targetVTKName,false);
        }
    }

    return;
}


void Foam::meshRefinement::sourceTargetChecking
(
    const polyMesh& mesh,
    const dictionary& patchMatchingDict
)
{
    Info<< nl;
    Info<< "Checking patch boundary matching" <<endl;
    Info<< "--------------------------------" <<endl;

    forAllConstIter(dictionary, patchMatchingDict, iter)
    {
        const word& key = iter().keyword();
        const dictionary& dict = patchMatchingDict.subDict(key);

        if (!dict.found("sourcePatches"))
        {
            WarningInFunction
               << "No sourcePatches for dictionary " << dict << endl;
            continue;
        }

        if (!dict.found("targetPatches"))
        {
            WarningInFunction
               << "No targetPatches for dictionary " << dict << endl;
            continue;
        }

        const wordReList sourcePatches  = dict.lookup("sourcePatches");
        const wordReList targetPatches  = dict.lookup("targetPatches");
        bool writeVTK = dict.lookupOrDefault<Switch>
        (
            "writeVTK",
            true
        );

        calcSourceTargetMatch
        (
            mesh,
            sourcePatches,
            targetPatches,
            writeVTK
        );
    }
}


void Foam::meshRefinement::triangulateProblemBoundaryFaces
(
    const dictionary& triangulateDict
)
{
    polyTopoChange meshMod(mesh_);

    scalar maxWarpage =
        triangulateDict.lookupOrDefault<scalar>("maxWarpage",GREAT);

    label nWarped = 0;
    if (maxWarpage != GREAT)
    {
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();

        const wordList includedPatches =
            triangulateDict.lookupOrDefault("includedPatches", wordList());
        DynamicList<label> includedPatchesIDs(includedPatches.size());

        forAll(includedPatches, i)
        {
            word includeName = includedPatches[i];
            const label patchI = patches.findPatchID(includeName);
            if (patchI != -1)
            {
                includedPatchesIDs.append(patchI);
            }
        }
        includedPatchesIDs.shrink();
        if (includedPatchesIDs.size())
        {
            Info<<"Triangulating problem faces"<<endl;

            const scalar edge0Len = meshCutter_.level0EdgeLength();

            autoPtr<indirectPrimitivePatch> ppPtr
            (
                meshRefinement::makePatch
                (
                    mesh_,
                    includedPatchesIDs
                )
            );
            const indirectPrimitivePatch& pp = ppPtr();

            const pointField& pts = mesh_.points();
            const pointField& fAreas = mesh_.faceAreas();
            const pointField& fCentres = mesh_.faceCentres();

            forAll(pp.localFaces(), fI)
            {
                label meshFaceI = pp.addressing()[fI];
                const face& f = mesh_.faces()[meshFaceI];
                if (f.size() == 3)
                {
                    continue;
                }

                label own = mesh_.faceOwner()[meshFaceI];
                vector fn = fAreas[meshFaceI] / (mag(fAreas[meshFaceI]) + SMALL);
                scalarField vProj(f.size());

                point fc = fCentres[meshFaceI];
                forAll(f, fp)
                {
                    vector n = pts[f[fp]] - fc;
                    vProj[fp] = (n & fn);
                }
                // Get normal 'span' of face
                scalar minVal = min(vProj);
                scalar maxVal = max(vProj);

                scalar edgeLen = edge0Len/(1<<meshCutter_.cellLevel()[own]);

                scalar warp = (maxVal - minVal)/edgeLen;

                if (warp >= maxWarpage)
                {
                    nWarped++;
                    label patchI = patches.whichPatch(meshFaceI);

                    label zoneID =
                        mesh_.faceZones().whichZone(meshFaceI);
                    bool zoneFlip = false;
                    if (zoneID >= 0)
                    {
                        const faceZone& fZone =
                            mesh_.faceZones()[zoneID];
                        zoneFlip =
                            fZone.flipMap()[fZone.whichFace(meshFaceI)];
                    }

                    label zoneI = mesh_.pointZones().whichZone(f[0]);
                    label newPointI = meshMod.addPoint
                    (
                        fc,         // point
                        -1,         // master point
                        zoneI,      // zone for point
                        true        // supports a cell
                    );

                    forAll(f, fp)
                    {
                        face newFace(3);
                        label next = f.fcIndex(fp);
                        newFace[0] = f[fp];
                        newFace[1] = f[next];
                        newFace[2] = newPointI;

                        if (fp == 0)
                        {
                            meshMod.modifyFace
                            (
                                newFace, // modified face
                                meshFaceI,     // label of face being modified
                                own,            // owner
                                -1,             // neighbour
                                false,          // face flip
                                patchI,      // new patch for face
                                zoneID,         // zone for face
                                zoneFlip        // face flip in zone
                             );
                        }
                        else
                        {
                            meshMod.addFace
                            (
                                newFace,
                                own,
                                -1,
                                -1,
                                -1,
                                meshFaceI,
                                false,
                                patchI,
                                zoneID,
                                zoneFlip
                            );
                        }
                    }
                }

            }
            Info<<"Triangulating "<<returnReduce(nWarped, sumOp<label>())
                <<" boundary faces"<<endl;
        }
    }


    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, true);

    // Update fields
    mesh_.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh_.clearOut();
    }
    updateMesh(map, labelList(0));

    return;
}


void Foam::meshRefinement::removeUnusedPoints(polyMesh& mesh)
{
    PackedBoolList markedPts(mesh.nPoints());

    forAll(mesh.faces(), faceI)
    {
        const face& f = mesh.faces()[faceI];
        forAll(f,fp)
        {
            markedPts[f[fp]] = true;
        }
    }

    polyTopoChange meshMod(mesh);
    label nUnused = 0;
    forAll(markedPts, ptI)
    {
        if (!markedPts[ptI])
        {
            meshMod.setAction(polyRemovePoint(ptI));
            nUnused++;
        }
    }

    if (returnReduce(nUnused, sumOp<label>()) != 0)
    {
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);
    }
}


void Foam::meshRefinement::calculateLeakPath
(
    const word leakType,
    const word leakPathName,
    const polyMesh& mesh,
    const boolList& blockedFace,
    const pointField& locationsInMesh,
    const pointField& locationsOutsideMesh,
    const writer<scalar>& leakPathFormatter
)
{
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    fileName outputDir;
    if (Pstream::master())
    {
        outputDir =
            mesh.time().rootPath()
            / mesh.time().globalCaseName()
            / "postProcessing"
            / mesh.pointsInstance();
        outputDir.clean();
        mkDir(outputDir);
    }

    // Write the leak path

    meshSearch searchEngine(mesh);
    shortestPathSet leakPath
    (
        leakType,
        mesh,
        searchEngine,
        coordSet::coordFormatNames[coordSet::coordFormat::DISTANCE],
        false,  //true,
        50,     // tbd. Number of iterations
        pbm.groupPatchIDs()["wall"],
        locationsInMesh,
        locationsOutsideMesh,
        blockedFace
    );

    // Split leak path according to segment. Note: segment index
    // is global (= index in locationsInsideMesh)
    List<pointList> segmentPoints;
    List<scalarList> segmentDist;
    {
        label nSegments = 0;
        if (leakPath.segments().size())
        {
            nSegments = max(leakPath.segments())+1;
        }
        reduce(nSegments, maxOp<label>());

        labelList nElemsPerSegment(nSegments, 0);
        for (label segmenti : leakPath.segments())
        {
            nElemsPerSegment[segmenti]++;
        }
        segmentPoints.setSize(nElemsPerSegment.size());
        segmentDist.setSize(nElemsPerSegment.size());
        forAll(nElemsPerSegment, i)
        {
            segmentPoints[i].setSize(nElemsPerSegment[i]);
            segmentDist[i].setSize(nElemsPerSegment[i]);
        }
        nElemsPerSegment = 0;

        forAll(leakPath, elemi)
        {
            label segmenti = leakPath.segments()[elemi];
            pointList& points = segmentPoints[segmenti];
            scalarList& dist = segmentDist[segmenti];
            label& n = nElemsPerSegment[segmenti];

            points[n] = leakPath[elemi];
            dist[n] = leakPath.curveDist()[elemi];
            n++;
        }
    }

    PtrList<coordSet> allLeakPaths(segmentPoints.size());
    forAll(allLeakPaths, segmenti)
    {
        // Collect data from all processors
        List<pointList> gatheredPts(Pstream::nProcs());
        gatheredPts[Pstream::myProcNo()] =
            std::move(segmentPoints[segmenti]);
        Pstream::gatherList(gatheredPts);

        List<scalarList> gatheredDist(Pstream::nProcs());
        gatheredDist[Pstream::myProcNo()] =
            std::move(segmentDist[segmenti]);
        Pstream::gatherList(gatheredDist);

        // Combine processor lists into one big list.
        pointList allPts
        (
            ListListOps::combine<pointList>
            (
                gatheredPts, accessOp<pointList>()
            )
        );
        scalarList allDist
        (
            ListListOps::combine<scalarList>
            (
                gatheredDist, accessOp<scalarList>()
            )
        );

        // Sort according to curveDist
        labelList indexSet;
        Foam::sortedOrder(allDist, indexSet);

        allLeakPaths.set
        (
            segmenti,
            new coordSet
            (
                leakPath.name(),
                leakPath.axis(),
                pointList(allPts, indexSet),
                scalarList(allPts.size(), scalar(segmenti))
            )
        );
    }

    fileName fName;
    if (Pstream::master())
    {
        List<List<scalarField>> allLeakData(1);
        List<scalarField>& varData = allLeakData[0];
        varData.setSize(allLeakPaths.size());
        forAll(allLeakPaths, segmenti)
        {
            varData[segmenti] = allLeakPaths[segmenti].curveDist();
        }

        const wordList valueSetNames(1, leakPathName);

        fName = outputDir / leakPathFormatter.getFileName
        (
            allLeakPaths[0],
            valueSetNames
        );

        // Note scope to force writing to finish before
        // FatalError exit
        OFstream ofs(fName);
        if (ofs.opened())
        {
            leakPathFormatter.write
            (
                true,               // write tracks
                allLeakPaths,
                valueSetNames,
                allLeakData,
                ofs
            );
        }
    }

    Pstream::scatter(fName);

    return;
}


Foam::labelList Foam::meshRefinement::selectInnerFaces
() const
{
    const cellZoneMesh& cellZones = mesh_.cellZones();
    const label innerID = cellZones.findZoneID("innerGrid");

    if (innerID >= 0)
    {
        Info<< "Selecting faces in cellZone: innerGrid" << endl;

        labelHashSet zoneFaces(2*cellZones[innerID].size());

        const cellList& cells = mesh_.cells();

        forAll(cellZones[innerID], cellI)
        {
            const cell& c = cells[cellZones[innerID][cellI]];
            forAll(c, cf)
            {
                zoneFaces.insert(c[cf]);
            }
        }

        labelList checkFaces(zoneFaces.size());

        label n = 0;
        forAllConstIter(labelHashSet, zoneFaces, iter)
        {
            checkFaces[n++] = iter.key();
        }

        return checkFaces;
    }
    else
    {
        return labelList(identity(mesh_.nFaces()));
    }
} // selectInnerFaces

// ************************************************************************* //
