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
    (c) 2010-2021 Esi Ltd.
    (c) 2010-2012 ICON CFD
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2015-2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "snappyHexMeshDriver/snappyRefineDriver.H"
#include "meshRefinement/meshRefinement.H"
#include "fvMesh/fvMesh.H"
#include "db/Time/Time.H"
#include "sets/topoSets/cellSet.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "snappyHexMeshDriver/refinementParameters/refinementParameters.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "refinementFeatures/refinementFeatures.H"
#include "shellSurfaces/shellSurfaces.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "global/unitConversion/unitConversion.H"
#include "snappyHexMeshDriver/snapParameters/snapParameters.H"
#include "snappyHexMeshDriver/layerParameters/layerParameters.H"
#include "regionSplit/localPointRegion.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "primitives/Vector/labelVector/labelVector.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "triSurface/orientedSurface/orientedSurface.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "anisoRefiner/anisoRefiner.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "snappyHexMeshDriver/snappyVoxelMeshDriver.H"
#if !defined( WIN32 ) && !defined( WIN64 )
#include "global/profiling/profiling.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(snappyRefineDriver, 0);

} // End namespace Foam


void Foam::snappyRefineDriver::addEdgeToRegion
(
    const edge& e,
    const label region,
    const label maxLevel,
    EdgeMap<labelPair>& levelPerEdge
) const
{
    EdgeMap<labelPair>::iterator eFnd = levelPerEdge.find(e);
    if (eFnd != levelPerEdge.end())
    {
        if (maxLevel > eFnd()[1])
        {
            eFnd()[0] = region;
            eFnd()[1] = maxLevel;
        }
    }
    else
    {
        labelPair levelData(region, maxLevel);
        levelPerEdge.insert(e, levelData);
    }
}

Foam::labelList Foam::snappyRefineDriver::consistentRefinement
(
    const label nBuffer,
    const labelList& candidateCells
) const
{
    // Problem choosing starting faces for bufferlayers (bFaces)
    //  - we can't use the current intersected boundary faces
    //    (intersectedFaces) since this grows indefinitely
    //  - if we use 0 faces we don't satisfy bufferLayers from the
    //    surface.
    //  - possibly we want to have bFaces only the initial set of faces
    //    and maintain the list while doing the refinement.
    labelList bFaces
    (
        findIndices(meshRefiner_.userFaceData()[0].second(), 0)
     );

    //Info<< "Collected boundary faces : "
    //    << returnReduce(bFaces.size(), sumOp<label>()) << endl;

    labelList cellsToRefine;

    if (nBuffer <= 2)
    {
        return meshRefiner_.meshCutter().consistentSlowRefinement
        (
            nBuffer,
            candidateCells,                     // cells to refine
            bFaces,                             // faces for nBufferLayers
            1,                                  // point difference
            meshRefiner_.intersectedPoints()    // points to check
         );
    }
    else
    {
        return meshRefiner_.meshCutter().consistentSlowRefinement2
        (
            nBuffer,
            candidateCells,                 // cells to refine
            bFaces                          // faces for nBufferLayers
         );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::snappyRefineDriver::snappyRefineDriver
(
    meshRefinement& meshRefiner,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const meshControl& controller,
    const writer<scalar>& setFormatter
)
:
    meshRefiner_(meshRefiner),
    decomposer_(decomposer),
    distributor_(distributor),
    globalToMasterPatch_(globalToMasterPatch),
    globalToSlavePatch_(globalToSlavePatch),
    controller_(controller),
    setFormatter_(setFormatter)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::snappyRefineDriver::initialRefine
(
    const refinementParameters& refineParams,
    const label maxIter,
    const label minRefine
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(edge, "snappyHexMesh::refine::init");
    #endif
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    label nFtrEdges = meshRefiner_.features()().size();
    reduce(nFtrEdges, sumOp<label>());

    for (; iter < maxIter; iter++)
    {
        Info<< nl
            << "Initial refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                iter,
                refineParams,

                true,               // initialRefinement
                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // proximity refinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false               // spreadGapSize
            )
        );

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );

        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for initial refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        if (nCellsToRefine <= minRefine)
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }

        if (debug > 0)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if
        (
            refineParams.balanceRefine() &&
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
             )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "Initial refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "Initial refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }

    //Clear up manifold features
    meshRefiner_.features()().manFeatures().clear();

    return iter;
}


Foam::label Foam::snappyRefineDriver::featureEdgeRefine
(
    const refinementParameters& refineParams,
    const label maxIter,
    const label minRefine
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(edge, "snappyHexMesh::refine::edge");
    #endif
    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    label nFtrEdges = meshRefiner_.features()().size();
    reduce(nFtrEdges, sumOp<label>());

    if (nFtrEdges && maxIter > 0)
    {
        for (; iter < maxIter; iter++)
        {
            Info<< nl
                << "Feature refinement iteration " << iter << nl
                << "------------------------------" << nl
                << endl;

            labelList candidateCells
            (
                meshRefiner_.refineCandidates
                (
                    iter,
                    refineParams,

                    false,              // initialRefinement
                    true,               // featureRefinement
                    false,              // featureDistanceRefinement
                    false,              // internalRefinement
                    false,              // surfaceRefinement
                    false,              // curvatureRefinement
                    false,              // proximity refinement
                    false,              // smallFeatureRefinement
                    false,              // gapRefinement
                    false,              // bigGapRefinement
                    false               // spreadGapSize
                )
            );

            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );

            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;


            label nCellsToRefine = cellsToRefine.size();
            reduce(nCellsToRefine, sumOp<label>());

            Info<< "Selected for feature refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            if (nCellsToRefine <= minRefine)
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }

            if (debug > 0)
            {
                const_cast<Time&>(mesh.time())++;
            }

            if
            (
                refineParams.balanceRefine() &&
                returnReduce
                (
                    (mesh.nCells() >= refineParams.maxLocalCells()),
                    orOp<bool>()
                )
            )
            {
               meshRefiner_.balanceAndRefine
                (
                    "feature refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams
                );
            }
            else
            {
                meshRefiner_.refineAndBalance
                (
                    "feature refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams
                );
            }

        }
    }

    //Clear up manifold features
    meshRefiner_.features()().manFeatures().clear();

    return iter;
}


void Foam::snappyRefineDriver::refineAtBaffleEdges
(
    const refinementParameters& refineParams
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    const hexRef8& cutter = meshRefiner_.meshCutter();
    const labelList& cellLevel = cutter.cellLevel();
    const labelList& pointLevel = cutter.pointLevel();

//    label iter = 0;
    while (true)
    {
        labelList nEdgeCells(mesh.nEdges(), 0);
        labelList boundaryEdges(mesh.nEdges(), -1);
        labelList boundaryFaces(mesh.nFaces(), -1);

        forAll(mesh.faces(), facei)
        {
            label patchi = patches.whichPatch(facei);
            if (patchi != -1 && !patches[patchi].coupled())
            {
                boundaryFaces[facei] = patchi;
                const labelList& fEdges = mesh.faceEdges()[facei];
                forAll(fEdges, fe)
                {
                    boundaryEdges[fEdges[fe]] = patchi;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            boundaryEdges,
            maxEqOp<label>(),
            label(-1)               // null value
        );

        forAll(mesh.edges(), edgei)
        {
            if (boundaryEdges[edgei] != -1)
            {
                const labelList& eCells = mesh.edgeCells()[edgei];
                nEdgeCells[edgei] += eCells.size();
            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            nEdgeCells,
            plusEqOp<label>(),
            label(0)               // null value
        );

        labelList candidateCells;
        cellSet candidateCellSet(mesh, "candidateCells", mesh.nCells()/1000);

        forAll(mesh.edges(), edgei)
        {
            if (nEdgeCells[edgei] == 4)
            {
                const labelList& eFaces = mesh.edgeFaces()[edgei];
                forAll(eFaces, efi)
                {
                    label facei = eFaces[efi];
                    if (boundaryFaces[facei] != -1)
                    {
                        label own = mesh.faceOwner()[facei];
                        label ownLevel = cellLevel[own];
                        const edge e = mesh.edges()[edgei];
                        if
                        (
                            pointLevel[e[0]] > ownLevel
                            || pointLevel[e[1]] > ownLevel
                        )
                        {
                            candidateCellSet.insert(own);
                            break;
                        }
                    }
                }
            }
            if (nEdgeCells[edgei] == 3)
            {
                const labelList& eCells = mesh.edgeCells()[edgei];
                forAll(eCells, eci)
                {
                    label celli = eCells[eci];
                    const cell& c = mesh.cells()[celli];
                    labelHashSet cSet(c);
                    const labelList& eFaces = mesh.edgeFaces()[edgei];
                    DynamicList<label> ceFaces(eFaces.size());

                    forAll(eFaces, efi)
                    {
                        label facei = eFaces[efi];
                        if (cSet.found(facei))
                        {
                            ceFaces.append(facei);
                        }
                    }

                    if (ceFaces.size() == 2)
                    {
                        vector fN0 = mesh.faceAreas()[ceFaces[0]];
                        fN0 /= (mag(fN0) + SMALL);
                        vector fN1 = mesh.faceAreas()[ceFaces[1]];
                        fN1 /= (mag(fN1) + SMALL);
                        if ((fN0 & fN1) > 0.99)
                        {
                            candidateCellSet.insert(celli);
                            break;
                        }
                    }
                }
            }
        }

        candidateCells = candidateCellSet.toc();

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine. After a few iterations check if too
        // few cells
        if
        (
            nCellsToRefine == 0
        )
        {
            break;
        }
        else
        {
//            iter++;
        }

        meshRefiner_.refine(cellsToRefine,refineParams.locationsInMesh());
    }
}


Foam::label Foam::snappyRefineDriver::smallFeatureRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(feature, "snappyHexMesh::refine::smallFeature");
    #endif
    const fvMesh& mesh = meshRefiner_.mesh();


    label iter = 0;

    // See if any surface has an extendedGapLevel
    labelList surfaceMaxLevel(meshRefiner_.surfaces().maxGapLevel());
    labelList shellMaxLevel(meshRefiner_.shells().maxGapLevel());

    if (max(surfaceMaxLevel) == 0 && max(shellMaxLevel) == 0)
    {
        return iter;
    }

    for (; iter < maxIter; iter++)
    {
        Info<< nl
            << "Small surface feature refinement iteration " << iter << nl
            << "--------------------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                iter,
                refineParams,

                false,              // initialRefinement
                false,               // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // proximity refinement
                true,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false               // spreadGapSize
            )
        );

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if (nCellsToRefine == 0)
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "small feature refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "small feature refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::surfaceOnlyRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(surface, "snappyHexMesh::refine::surface");
    #endif
    const fvMesh& mesh = meshRefiner_.mesh();

    // Determine the maximum refinement level over all surfaces. This
    // determines the minumum number of surface refinement iterations.
    label overallMaxLevel = max
    (
        max(meshRefiner_.surfaces().maxLevel()),
        meshRefiner_.shells().maxInsideOutsideLevel()
    );

    // Mark current boundary faces with 0. Have meshRefiner maintain them.
    meshRefiner_.userFaceData().setSize(1);

    // mark list to remove any refined faces
    meshRefiner_.userFaceData()[0].first() = meshRefinement::REMOVE;
    meshRefiner_.userFaceData()[0].second() = createWithValues<labelList>
    (
        mesh.nFaces(),
        -1,
//        meshRefiner_.intersectedFaces(),
        meshRefiner_.intersectedAndNeighbouring(),
        0
    );

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Surface refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        labelList initialCells
        (
            meshRefiner_.refineCandidates
            (
                iter,
                refineParams,

                false,              // initialRefinement
                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                true,               // surfaceRefinement
                true,               // curvatureRefinement
                false,              // proximity refinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false               // spreadGapSize
            )
        );
        labelList candidateCells
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                initialCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = candidateCells.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if
        (
            refineParams.balanceRefine() &&
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                candidateCells,//cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                candidateCells,//cellsToRefine,
                refineParams
            );
        }
    }

    //Clear any fields created for curvature refinement
    const labelList& curvatureSurfaces =
        meshRefiner_.surfaces().curvatureSurfaces();
    forAll(curvatureSurfaces, i)
    {
        label surfI = curvatureSurfaces[i];
        const searchableSurface& geom =
            meshRefiner_.surfaces().geometry()[surfI];
        if (isA<triSurface>(geom))
        {
            const_cast<searchableSurface&>(geom).clearCurvaturePointField();
        }
    }

    return iter;
}


Foam::label Foam::snappyRefineDriver::gapOnlyRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    const fvMesh& mesh = meshRefiner_.mesh();

    // Determine the maximum refinement level over all surfaces. This
    // determines the minumum number of surface refinement iterations.

    label maxIncrement = 0;
    const labelList& maxLevel = meshRefiner_.surfaces().maxLevel();
    const labelList& gapLevel = meshRefiner_.surfaces().gapLevel();

    forAll(maxLevel, i)
    {
        maxIncrement = max(maxIncrement, gapLevel[i]-maxLevel[i]);
    }

    label iter = 0;

    if (maxIncrement == 0)
    {
        return iter;
    }

    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Gap refinement iteration " << iter << nl
            << "--------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                iter,
                refineParams,

                false,              // initialRefinement
                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,               // surfaceRefinement
                false,               // curvatureRefinement
                false,               // proximity refinement
                false,              // smallFeatureRefinement
                true,              // gapRefinement
                false,              // bigGapRefinement
                false               // spreadGapSize

            )
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromGap." << endl;
            cellSet c(mesh, "candidateCellsFromGap", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        // Grow by one layer to make sure we're covering the gap
        {
            boolList isCandidateCell(mesh.nCells(), false);
            forAll(candidateCells, i)
            {
                isCandidateCell[candidateCells[i]] = true;
            }

            for (label i=0; i<1; i++)
            {
                boolList newIsCandidateCell(isCandidateCell);

                // Internal faces
                for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
                {
                    label own = mesh.faceOwner()[facei];
                    label nei = mesh.faceNeighbour()[facei];

                    if (isCandidateCell[own] != isCandidateCell[nei])
                    {
                        newIsCandidateCell[own] = true;
                        newIsCandidateCell[nei] = true;
                    }
                }

                // Get coupled boundary condition values
                boolList neiIsCandidateCell;
                syncTools::swapBoundaryCellList
                (
                    mesh,
                    isCandidateCell,
                    neiIsCandidateCell
                );

                // Boundary faces
                for
                (
                    label facei = mesh.nInternalFaces();
                    facei < mesh.nFaces();
                    facei++
                )
                {
                    label own = mesh.faceOwner()[facei];
                    label bFacei = facei-mesh.nInternalFaces();

                    if (isCandidateCell[own] != neiIsCandidateCell[bFacei])
                    {
                        newIsCandidateCell[own] = true;
                    }
                }

                isCandidateCell.transfer(newIsCandidateCell);
            }

            label n = 0;
            forAll(isCandidateCell, celli)
            {
                if (isCandidateCell[celli])
                {
                    n++;
                }
            }
            candidateCells.setSize(n);
            n = 0;
            forAll(isCandidateCell, celli)
            {
                if (isCandidateCell[celli])
                {
                    candidateCells[n++] = celli;
                }
            }
        }


        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromGapPlusBuffer." << endl;
            cellSet c(mesh, "candidateCellsFromGapPlusBuffer", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }


        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= maxIncrement
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }
    return iter;
}


Foam::label Foam::snappyRefineDriver::proximityRefine
(
    const refinementParameters& refineParams,
    const label maxIter,
    const label minRefine,
    const bool oneMachine
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    const bool ldebug (false);

    const fvMesh& mesh = meshRefiner_.mesh();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Proximity refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;

        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~
        // Only look at surface intersections (minLevel and surface curvature),
        // do not do internal refinement (refinementShells)

        const PtrList<featureEdgeMesh> dummyFeatures;

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                iter,
                refineParams,

                false,              // initialRefinement
                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,               // surfaceRefinement
                false,               // curvatureRefinement
                true,               // proximity refinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false               // spreadGapSize

            )

        );

        if (ldebug)
        {
            Info<<"snappyRefineDriver::proximityRefine(): cells to refine:"
                <<endl;
            Info<< candidateCells<<endl;
        }

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );

        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum nessecary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine <= minRefine
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }

        if
        (
            refineParams.balanceRefine() &&
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "surface refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }
    return iter;
}

void Foam::snappyRefineDriver::removalInterfaceRefine
(
    const refinementParameters& refineParams
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return;
    }

    const fvMesh& mesh = meshRefiner_.mesh();

    labelList refineCells
    (
        meshRefiner_.markRemovalInterfaceRefine
        (
            refineParams
        )
    );

    labelList candidateCells
    (
        meshRefiner_.meshCutter().consistentRefinement
        (
            refineCells,
            true
        )
    );

    Info<< "Determined cells to refine in = "
        << mesh.time().cpuTimeIncrement() << " s" << endl;


    label nCellsToRefine = candidateCells.size();
    reduce(nCellsToRefine, sumOp<label>());

    Info<< "Selected for refinement : " << nCellsToRefine
        << " cells (out of " << mesh.globalData().nTotalCells()
        << ')' << endl;

    // Stop when no cells to refine or have done minimum necessary
    // iterations and not enough cells to refine.
    if
    (
        nCellsToRefine == 0
        || nCellsToRefine <= refineParams.minRefineCells()
    )
    {
        Info<< "Stopping refining since too few cells selected."
            << nl << endl;
        return;
    }


    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    if
    (
        refineParams.balanceRefine() &&
        returnReduce
        (
            (mesh.nCells() >= refineParams.maxLocalCells()),
            orOp<bool>()
        )
    )
    {
        meshRefiner_.balanceAndRefine
        (
            "removal refinement",
            decomposer_,
            distributor_,
            candidateCells,
            refineParams
        );
    }
    else
    {
        meshRefiner_.refineAndBalance
        (
            "removal refinement",
            decomposer_,
            distributor_,
            candidateCells,
            refineParams
         );
    }

    return;
}


Foam::label Foam::snappyRefineDriver::bigGapOnlyRefine
(
    const refinementParameters& refineParams,
    const bool spreadGapSize,
    const label maxIter
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    // See if any surface has an extendedGapLevel
    labelList surfaceMaxLevel(meshRefiner_.surfaces().maxGapLevel());
    labelList shellMaxLevel(meshRefiner_.shells().maxGapLevel());

    label overallMaxLevel(max(max(surfaceMaxLevel), max(shellMaxLevel)));

    if (overallMaxLevel == 0)
    {
        return iter;
    }


    for (; iter < maxIter; iter++)
    {
        Info<< nl
            << "Big gap refinement iteration " << iter << nl
            << "------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                iter,
                refineParams,

                false,              // initialRefinement
                false,              // featureRefinement
                false,              // featureDistanceRefinement
                false,              // internalRefinement
                false,               // surfaceRefinement
                false,               // curvatureRefinement
                true,               // proximity refinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                spreadGapSize       // spreadGapSize
            )
        );


        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromBigGap." << endl;
            cellSet c(mesh, "candidateCellsFromBigGap", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "big gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "big gap refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }
    return iter;
}


void Foam::snappyRefineDriver::removeWrappedCells
(
    const refinementParameters& refineParams,
    const dictionary& refineDict,
    const bool handleSnapProblems,
    const dictionary& motionDict
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return;
    }

    Info<< nl
        << "Removing mesh beyond surface and wrapped intersections" << nl
        << "------------------------------------------------------" << nl
        << endl;

    meshRefiner_.splitWrappedMesh
    (
        refineParams,
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineDict,
        motionDict,
        handleSnapProblems,             // detect&remove potential snap problem
        false,                          // perpendicular edge connected cells
        scalarField(0),                 // per region perpendicular angle
        setFormatter_
     );
}


// Refine cells to remove dangling cells
Foam::label Foam::snappyRefineDriver::danglingCellRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    const scalar oppositeCos = Foam::cos(Foam::degToRad(135));

    const bool interRefine = refineParams.interfaceRefine();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Dangling cells refinement iteration " << iter << nl
            << "-------------------------------------" << nl
            << endl;

        const fvMesh& mesh = meshRefiner_.mesh();
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        const hexRef8& cutter = meshRefiner_.meshCutter();

        const labelList& cellLevel = cutter.cellLevel();
        const labelList& pointLevel = cutter.pointLevel();

        const labelList& owners = mesh.faceOwner();
        const vectorField& fA = mesh.faceAreas();
        const labelList& surfaceIndex = meshRefiner_.surfaceIndex();

        labelList candidateCells;
        cellSet candidateCellSet(mesh, "candidateCells", mesh.nCells()/1000);
        labelList newCellLevel = cellLevel;

        labelList cellsToRefine;

        label maxRefined = labelMax;
        bool converged = true;
        while (true)
        {
            cellSet newCandidateCellSet
            (
                mesh,
                "newCandidateCells",
                mesh.nCells()/1000
            );
            candidateCells = candidateCellSet.toc();

            cellsToRefine = meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            );

            forAll(cellsToRefine, i)
            {
                label celli = cellsToRefine[i];
                newCellLevel[celli] = cellLevel[celli] + 1;
                if (candidateCellSet.insert(celli))
                {
                    newCandidateCellSet.insert(celli);
                }
            }

            // Calculate neighbour cell level
            scalarField neiLevel(mesh.nFaces()-mesh.nInternalFaces());
            for
            (
                label faceI = mesh.nInternalFaces();
                faceI < mesh.nFaces();
                faceI++
            )
            {
                neiLevel[faceI-mesh.nInternalFaces()] =
                    newCellLevel[owners[faceI]];
            }
            syncTools::swapBoundaryFaceList(mesh, neiLevel);

            PackedBoolList refinedInternalFace(mesh.nFaces());
            forAll(mesh.faces(), facei)
            {
                label own = owners[facei];

                if (mesh.isInternalFace(facei))
                {
                    label nei = mesh.faceNeighbour()[facei];

                    if (newCellLevel[own] != newCellLevel[nei])
                    {
                        refinedInternalFace.set(facei);
                    }
                }
                else
                {
                    label patchI = patches.whichPatch(facei);
                    if (patches[patchI].coupled())
                    {
                        if
                        (
                            newCellLevel[own]
                            != neiLevel[facei-mesh.nInternalFaces()]
                        )
                        {
                            refinedInternalFace.set(facei);
                        }
                    }
                }
            }

            PackedBoolList splitFaces(mesh.nFaces());
            forAll(mesh.cells(), celli)
            {
                label cLevel = cellLevel[celli];
                if (newCellLevel[celli] != cLevel)
                {
                    const labelList& cFaces = mesh.cells()[celli];
                    forAll(cFaces, cFI)
                    {
                        label facei = cFaces[cFI];
                        if (refinedInternalFace.get(facei))
                        {
                            const face& f = mesh.faces()[facei];
                            label nAnchors = 0;
                            forAll(f,fp)
                            {
                                label pointi = f[fp];
                                if (pointLevel[pointi] <= cLevel)
                                {
                                    nAnchors++;
                                }
                            }
                            if (nAnchors == 4)
                            {
                                splitFaces.set(facei);
                            }

                        }
                    }
                }
            }

            syncTools::syncFaceList
            (
                mesh,
                splitFaces,
                orEqOp<unsigned int>()
            );

            forAll(mesh.cells(), celli)
            {
                label cLevel = cellLevel[celli];
                if (newCellLevel[celli] != cLevel)
                {
                    continue;
                }
                const labelList& cFaces = mesh.cells()[celli];
                label nInterFaces = 0;
                forAll(cFaces, cFI)
                {
                    label facei = cFaces[cFI];
                    label bFacei = facei-mesh.nInternalFaces();
                    if (bFacei >= 0)
                    {
                        label patchi = patches.patchID()[bFacei];
                        if (!patches[patchi].coupled())
                        {
                            continue;
                        }
                    }

                    if (refinedInternalFace.get(facei))
                    {
                        if (splitFaces.get(facei))
                        {
                            nInterFaces += 4;
                        }
                        else
                        {
                            nInterFaces++;
                        }
                    }
                }

                if (nInterFaces > 15)
                {
                    if (candidateCellSet.insert(celli))
                    {
                        newCandidateCellSet.insert(celli);
                    }
                }
                else if (interRefine)
                {
                    // Detect opposite intersection
                    bool foundOpposite = false;

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        if
                        (
                            surfaceIndex[facei] != -1
                            && refinedInternalFace.get(facei)
                            &&
                            (
                                splitFaces.get(facei)
                                || cutter.faceLevel(facei) != cLevel
                            )
                        )
                        {
                            // Get outwards pointing normal
                            vector n = fA[facei]/mag(fA[facei]);
                            if (owners[facei] != celli)
                            {
                                n = -n;
                            }

                            // Check for any opposite intersection
                            forAll(cFaces, cFaceI2)
                            {
                                label face2i = cFaces[cFaceI2];

                                if
                                (
                                    face2i != facei
                                 && surfaceIndex[face2i] != -1
                                )
                                {
                                    // Get outwards pointing normal
                                    vector n2 = fA[face2i]/mag(fA[face2i]);
                                    if (owners[face2i] != celli)
                                    {
                                        n2 = -n2;
                                    }

                                    if ((n&n2) < oppositeCos)
                                    {
                                        label  nOpenFaces = 0;
                                        forAll(cFaces, cFaceI3)
                                        {
                                            label face3i = cFaces[cFaceI3];
                                            if (surfaceIndex[face3i] == -1)
                                            {
                                                nOpenFaces++;
                                            }
                                        }
                                        if (nOpenFaces > 1)
                                        {
                                            foundOpposite = true;
                                            break;
                                        }
                                    }
                                }
                            }

                            if (foundOpposite)
                            {
                                break;
                            }
                        }
                    }

                    if (foundOpposite)
                    {
                        if (candidateCellSet.insert(celli))
                        {
                            newCandidateCellSet.insert(celli);
                        }
                    }
                }
            }

            labelList newCandidateCells = newCandidateCellSet.toc();
            label newCellsToRefine = newCandidateCells.size();
            reduce(newCellsToRefine, sumOp<label>());
            if
            (
                newCellsToRefine == 0
            )
            {
                break;
            }
            else if (newCellsToRefine >= maxRefined)
            {
                converged = false;
                break;
            }
            maxRefined = newCellsToRefine;
        }

        if (converged)
        {
            cellsToRefine = candidateCellSet.toc();
        }
        else
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCellSet.toc(),
                true
            );
        }

        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;

        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

       Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        if (nCellsToRefine == 0)
        {
            break;
        }

        if
        (
            refineParams.balanceRefine() &&
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
             )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "Dangling refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "Dangling refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }

    return iter;
}


Foam::label Foam::snappyRefineDriver::singleLevelEdgeRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(dangling, "snappyHexMesh::refine::singleLevelEdgeRefine");
    #endif
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Single level edge jump refinement iteration " << iter << nl
            << "---------------------------------------------" << nl
            << endl;


        // Determine cells to refine
        // ~~~~~~~~~~~~~~~~~~~~~~~~~

        const cellList& cells = mesh.cells();

        labelList candidateCells;
        {
            cellSet candidateCellSet(mesh, "candidateCells", cells.size()/1000);

            labelList minEdgeLevel(mesh.nEdges(), labelMax);
            labelList maxEdgeLevel(mesh.nEdges(), -labelMax);
            forAll(mesh.edges(), edgeI)
            {
                const labelList& eCells = mesh.edgeCells()[edgeI];
                forAll(eCells, eCI)
                {
                    label celli = eCells[eCI];
                    label clevel = cellLevel[celli] ;
                    minEdgeLevel[edgeI] = min(clevel, minEdgeLevel[edgeI]);
                    maxEdgeLevel[edgeI] = max(clevel, maxEdgeLevel[edgeI]);
                }
            }

            syncTools::syncEdgeList
            (
                mesh,
                minEdgeLevel,
                minEqOp<label>(),
                labelMax               // null value
            );

            syncTools::syncEdgeList
            (
                mesh,
                maxEdgeLevel,
                maxEqOp<label>(),
                -labelMax               // null value
            );

            forAll(cells, celli)
            {
                label clevel = cellLevel[celli];
                const labelList& cEdges = mesh.cellEdges()[celli];
                forAll(cEdges, cEI)
                {
                    label edgei = cEdges[cEI];
                    label cdiff = maxEdgeLevel[edgei] - minEdgeLevel[edgei];
                    if (cdiff > 1 && clevel == minEdgeLevel[edgei])
                    {
                        candidateCellSet.insert(celli);
                        break;
                    }
                }
            }

            if (debug&meshRefinement::MESH)
            {
                Pout<< "Dumping " << candidateCellSet.size()
                    << " cells to cellSet candidateCellSet." << endl;
                candidateCellSet.instance() = meshRefiner_.timeName();
                candidateCellSet.write();
            }
            candidateCells = candidateCellSet.toc();
        }

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine. After a few iterations check if too
        // few cells
        if (nCellsToRefine == 0)
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "edge jump refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "edge jump refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }
    return iter;
}


void Foam::snappyRefineDriver::removeDualProblemCells
(
    const refinementParameters& refineParams
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    const labelListList& cEdgeList = mesh.cellEdges();
    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(interface, "snappyHexMesh::refine::transition");
    #endif

    const labelList meshedPatches = meshRefiner_.meshedPatches();

    autoPtr<indirectPrimitivePatch> ppPtr
    (
        meshRefinement::makePatch
        (
            mesh,
            meshedPatches
        )
    );
    const indirectPrimitivePatch& pp = ppPtr();

    boolList boundaryPts(mesh.nPoints(), false);
    forAll(pp.meshPoints(), ptI)
    {
        label meshPointI = pp.meshPoints()[ptI];

        boundaryPts[meshPointI] = true;
    }

    syncTools::syncPointList
    (
        mesh,
        boundaryPts,
        orEqOp<bool>(),
        false
    );

    boolList boundaryEdges(mesh.nEdges(), false);
    labelList meshEdges(pp.meshEdges(mesh.edges(), mesh.pointEdges()));
    forAll(meshEdges, edgeI)
    {
        label meshEdgeI = meshEdges[edgeI];
        boundaryEdges[meshEdgeI] = true;
    }

    syncTools::syncEdgeList
    (
        mesh,
        boundaryEdges,
        orEqOp<bool>(),
        false
    );

    // Edge status:
    //  >0 : label of number of cuts (1 or 2)
    //  -1 : no cuts
    labelList cutEdges(mesh.nEdges(), 0);

    forAll(mesh.edges(), edgeI)
    {
        if (!boundaryEdges[edgeI])
        {
            const edge& e = mesh.edges()[edgeI];
            if (boundaryPts[e[0]] || boundaryPts[e[1]])
            {
                if (boundaryPts[e[0]])
                {
                    cutEdges[edgeI]++;
                }

                if (boundaryPts[e[1]])
                {
                    cutEdges[edgeI]++;
                }
            }
            else
            {
                cutEdges[edgeI] = -1;
            }
        }
        else
        {
            cutEdges[edgeI] = -1;
        }
    }

    syncTools::syncEdgeList
    (
        mesh,
        cutEdges,
        maxEqOp<label>(),
        labelMin
    );

    label nChecks = 0;
    forAll(mesh.cells(), cellI)
    {
        const labelList& cellEdges = cEdgeList[cellI];

        forAll(cellEdges, cEI)
        {
            label edgeI = cellEdges[cEI];
            if (cutEdges[edgeI] != -1)
            {
                if (cutEdges[edgeI] == 1)
                {
                    nChecks++;
                }
                else
                {
                    nChecks += 2;
                }
            }
        }
    }

    labelList checkedCells(nChecks);
    pointField start(nChecks);
    pointField end(nChecks);

    nChecks = 0;
    forAll(mesh.cells(), cellI)
    {
        const labelList& cellEdges = cEdgeList[cellI];
        const point cc = mesh.cellCentres()[cellI];

        forAll(cellEdges, cEI)
        {
            label edgeI = cellEdges[cEI];
            if (cutEdges[edgeI] != -1)
            {
                edge e = mesh.edges()[edgeI];
                vector eVec = e.vec(mesh.points());

                if (cutEdges[edgeI] == 1)
                {
                    checkedCells[nChecks] = cellI;
                    start[nChecks] = cc;
                    vector midPoint = mesh.points()[e[0]]
                        + 0.5*eVec;

                    end[nChecks] = midPoint;
                    nChecks++;
                }
                else
                {
                    checkedCells[nChecks] = cellI;
                    start[nChecks] = cc;
                    vector cutPoint = mesh.points()[e[0]]
                        + 0.33*eVec;
                    end[nChecks] = cutPoint;
                    nChecks++;

                    checkedCells[nChecks] = cellI;
                    start[nChecks] = cc;
                    cutPoint = mesh.points()[e[0]]
                        + 0.66*eVec;
                    end[nChecks] = cutPoint;
                    nChecks++;
                }
            }
        }
    }

    labelList surface1;
    List<pointIndexHit> hit1;
    labelList region1;
    labelList surface2;
    List<pointIndexHit> hit2;
    labelList region2;

    const refinementSurfaces& surfaces = meshRefiner_.surfaces();

    // Surfaces that need to be baffled
    const labelList surfacesToBaffle
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfaces.surfZones())
    );

    surfaces.findNearestIntersection
    (
        surfacesToBaffle,
        start,
        end,

        surface1,
        hit1,
        region1,

        surface2,
        hit2,
        region2
     );

    // Get cells to remove
    label foundCell = -1;

    labelList ownPatch(mesh.nFaces(), 0);
    boolList blockedFace(mesh.nFaces(),false);
    boolList problemCells(mesh.nCells(), false);
    forAll(start, i)
    {
        label cellI = checkedCells[i];
        if (hit1[i].hit() && cellI != foundCell)
        {
            label patchI = globalToMasterPatch_
            [
                surfaces.globalRegion(surface1[i], region1[i])
            ];

            const cell& c = mesh.cells()[cellI];
            forAll(c, cFI)
            {
                label faceI = c[cFI];
                ownPatch[faceI] = patchI;
                blockedFace[faceI] = true;
            }

            problemCells[cellI] = true;
            foundCell = cellI;
        }
    }

    // Set region per cell based on walking
    regionSplit cellRegion(mesh, blockedFace);


    const pointField& locationsInMesh = refineParams.locationsInMesh();
    DynamicList<label> regionsToKeep(cellRegion.nRegions());

    // For all locationsInMesh find the cell
    forAll(locationsInMesh, i)
    {
        // Get location and index of zone ("none" for cellZone -1)
        const point& insidePoint = locationsInMesh[i];

        // Find the region containing the keepPoint
        label keepRegionI = -1;

        label cellI = meshRefiner_.findCell
        (
            insidePoint,
            mesh,
            meshRefiner_.meshCutter()
        );

        if (cellI != -1)
        {
            keepRegionI = cellRegion[cellI];
        }
        reduce(keepRegionI, maxOp<label>());
/*
        // Find the region containing the keepPoint
        label keepRegionI = meshRefiner_.findRegion
        (
            mesh,
            cellRegion,
            meshRefiner_.mergeDistance()*vector(1,1,1),
            insidePoint
        );
*/
        regionsToKeep.append(keepRegionI);
    }
    regionsToKeep.shrink();
    labelHashSet regionsToKeepSet(regionsToKeep);

    DynamicList<label> cellsToRemove(mesh.nCells()/100);
    forAll(mesh.cells(), cellI)
    {
        if (problemCells[cellI] || !regionsToKeepSet.found(cellRegion[cellI]))
        {
            cellsToRemove.append(cellI);
        }
    }
    cellsToRemove.shrink();

    // Remove cells
    removeCells cellRemover(mesh);

    // Pick up patches for exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
    labelList exposedPatchIDs(exposedFaces.size());

    forAll(exposedFaces, i)
    {
        label faceI = exposedFaces[i];

        exposedPatchIDs[i] = ownPatch[faceI];
    }

    meshRefiner_.doRemoveCells
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        cellRemover
    );
}


// Refine cells to aid interface layer addition for dualised mesh pre-processing
Foam::label Foam::snappyRefineDriver::interfaceProblemCellRefine
(
    const refinementParameters& refineParams,
    const dictionary& layerDict,
    const label maxIter
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();

    layerParameters layerParams(layerDict, patches, true);
    const labelList& numLayers = layerParams.numLayers();

    DynamicList<label> noLayers(patches.size());
    forAll(patches, patchI)
    {
        if (!patches[patchI].coupled() && numLayers[patchI] == -1)
        {
            noLayers.append(patchI);
        }
    }
    noLayers.shrink();
    labelHashSet noLayerPatchIDs(noLayers);

    const bool removeUnsplittable = refineParams.removeUnsplittable();

    //Whether to perform additional refinement - improves surface topology
    // but can results in refinement propogating further
    const bool addedRefine = refineParams.additionalDualRefine();

    //Perform initial removal of problem cells
    meshRefiner_.baffleHoles();
    meshRefiner_.checkRefinedAndRemove(identity(mesh.nCells()),false,false);


    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Interface cells refinement iteration " << iter << nl
            << "---------------------------------------" << nl
            << endl;

        labelList candidateCells;
        cellSet candidateCellSet(mesh, "candidateCells", mesh.nCells()/1000);

        PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));
        PackedBoolList isMasterEdge(syncTools::getMasterEdges(mesh));
        labelList newCellLevel = cellLevel;

        boolList isBoundaryEdge(mesh.nEdges(),false);
        boolList isBoundaryPt(mesh.nPoints(),false);
        boolList isBoundaryFace(mesh.nFaces(),false);
        forAll(mesh.faces(), faceI)
        {
            label patchI = patches.whichPatch(faceI);
            if (patchI != -1 && !patches[patchI].coupled())
            {
                isBoundaryFace[faceI] = true;
                const labelList& fEdges = mesh.faceEdges()[faceI];

                forAll(fEdges, fEI)
                {
                    label edgeI = fEdges[fEI];
                    isBoundaryEdge[edgeI] = true;
                }
                const face& f = mesh.faces()[faceI];
                forAll(f, fp)
                {
                    isBoundaryPt[f[fp]] = true;
                }

            }
        }

        syncTools::syncEdgeList
        (
            mesh,
            isBoundaryEdge,
            orEqOp<bool>(),
            false               // null value
        );

        syncTools::syncPointList
        (
            mesh,
            isBoundaryPt,
            orEqOp<bool>(),
            false               // null value
        );

        labelList cellsToRefine;

        while (true)
        {
            label typeA(0),typeB(0),typeC(0),typeD(0),typeE(0),typeF(0);
            label typeG(0),typeH(0),typeI(0),typeJ(0),typeK(0),typeL(0);
            label typeM(0);

            cellSet newCandidateCellSet
            (
                mesh,
                "newCandidateCells",
                mesh.nCells()/1000
            );
            candidateCells = candidateCellSet.toc();

            cellsToRefine = meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            );

            forAll(cellsToRefine, i)
            {
                label cellI = cellsToRefine[i];
                newCellLevel[cellI] = cellLevel[cellI] + 1;
                if (candidateCellSet.insert(cellI))
                {
                    newCandidateCellSet.insert(cellI);
                }
            }

            // Calculate neighbour cell level
            scalarField neiLevel(mesh.nFaces()-mesh.nInternalFaces());
            for
            (
                label faceI = mesh.nInternalFaces();
                faceI < mesh.nFaces();
                faceI++
            )
            {
                neiLevel[faceI-mesh.nInternalFaces()] =
                    newCellLevel[mesh.faceOwner()[faceI]];
            }
            syncTools::swapBoundaryFaceList(mesh, neiLevel);

            PackedBoolList refinedInternalFace(mesh.nFaces());
            PackedBoolList refinedInterfacePoint(mesh.nPoints());
            forAll(mesh.faces(), faceI)
            {
                label own = mesh.faceOwner()[faceI];

                if (mesh.isInternalFace(faceI))
                {
                    label nei = mesh.faceNeighbour()[faceI];

                    if (newCellLevel[own] != newCellLevel[nei])
                    {
                        refinedInternalFace.set(faceI, 1u);
                    }
                }
                else
                {
                    label patchI = patches.whichPatch(faceI);
                    if (patches[patchI].coupled())
                    {
                        if
                        (
                            newCellLevel[own]
                            != neiLevel[faceI-mesh.nInternalFaces()]
                        )
                        {
                            refinedInternalFace.set(faceI, 1u);
                        }
                    }
                }
                if (refinedInternalFace.get(faceI) == 1)
                {
                    face f = mesh.faces()[faceI];
                    forAll(f,fp)
                    {
                        refinedInterfacePoint.set(f[fp], 1u);
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                refinedInterfacePoint,
                orEqOp<unsigned int>(),
                0
            );
/*
            labelList maxPointCellLevel(mesh.nPoints(), labelMin);
            forAll(mesh.points(), pointI)
            {
                const labelList& pCells = mesh.pointCells()[pointI];

                forAll(pCells, pCI)
                {
                    label cellI = pCells[pCI];
                    label cLevel = newCellLevel[cellI];

                    if (cLevel > maxPointCellLevel[pointI])
                    {
                        maxPointCellLevel[pointI] = cLevel;
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                maxPointCellLevel,
                maxEqOp<label>(),
                labelMin              // null value
            );
*/
            labelList maxCellPointLevel(mesh.nCells(), labelMin);
            forAll(mesh.cells(), cellI)
            {
                const labelList& cPts = mesh.cellPoints()[cellI];
                forAll(cPts, cPtI)
                {
                    maxCellPointLevel[cellI] = max
                    (
                        maxCellPointLevel[cellI],
                        pointLevel[cPts[cPtI]]
                    );
                }
            }

            labelList maxPointCellLevel(mesh.nPoints(), labelMin);
            forAll(mesh.points(), pointI)
            {
                const labelList& pCells = mesh.pointCells()[pointI];

                forAll(pCells, pCI)
                {
                    label cellI = pCells[pCI];
                    label cLevel = maxCellPointLevel[cellI];

                    if (cLevel > maxPointCellLevel[pointI])
                    {
                        maxPointCellLevel[pointI] = cLevel;
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                maxPointCellLevel,
                maxEqOp<label>(),
                labelMin              // null value
            );

            //Refine boundary interface cells if cannot be
            //grown up (i.e not perpendicular)
            labelList nRefinementFaces(mesh.nEdges(), 0);
            labelList minEdgeLevel(mesh.nEdges(), labelMax);

            forAll(mesh.edges(), edgeI)
            {
                const labelList& eCells = mesh.edgeCells()[edgeI];
                forAll(eCells, eCI)
                {
                    label cellI = eCells[eCI];
                    minEdgeLevel[edgeI] =
                        min(newCellLevel[cellI], minEdgeLevel[edgeI]);
                }
            }

            syncTools::syncEdgeList
            (
                mesh,
                minEdgeLevel,
                minEqOp<label>(),
                labelMax               // null value
            );

            labelList nHigherLevels(mesh.nEdges(), 0);
            labelList nLowerLevels(mesh.nEdges(), 0);
            forAll(mesh.edges(), edgeI)
            {
                const labelList& eCells = mesh.edgeCells()[edgeI];

                forAll(eCells, eCI)
                {
                    label lev = newCellLevel[eCells[eCI]];
                    if (lev > minEdgeLevel[edgeI])
                    {
                        nHigherLevels[edgeI]++;
                    }
                    else
                    {
                        nLowerLevels[edgeI]++;
                    }
                }

                const labelList& eFaces = mesh.edgeFaces()[edgeI];
                forAll(eFaces, eFI)
                {
                    label faceI = eFaces[eFI];
                    if
                    (
                        isMasterFace.get(faceI) == 1
                        && refinedInternalFace.get(faceI) == 1
                    )
                    {
                        nRefinementFaces[edgeI]++;
                    }
                }
            }

            syncTools::syncEdgeList
            (
                mesh,
                nHigherLevels,
                plusEqOp<label>(),
                label(0)               // null value
            );

            syncTools::syncEdgeList
            (
                mesh,
                nLowerLevels,
                plusEqOp<label>(),
                label(0)               // null value
            );

            syncTools::syncEdgeList
            (
                mesh,
                nRefinementFaces,
                plusEqOp<label>(),
                label(0)          // null value
            );

            labelList edgeType(mesh.nEdges(), -1);
            forAll(mesh.edges(), edgeI)
            {
                if ((nHigherLevels[edgeI] + nLowerLevels[edgeI]) > 4)
                {
                    edge e = mesh.edges()[edgeI];
                    Pout<<"More than 4 edge cells detected. This might indicate "
                        <<"a sync problem at edge :"
                        <<edgeI<<" "<<e.centre(mesh.points())<<endl;
                }

                if (nRefinementFaces[edgeI] == 1 || nRefinementFaces[edgeI] == 2)
                {
                    if (nHigherLevels[edgeI] == 3)
                    {
                        if
                        (
                            isBoundaryEdge[edgeI] && nLowerLevels[edgeI] == 1
                           && nRefinementFaces[edgeI] == 1
                        )
                        {
                            edgeType[edgeI] = 0;
                        }
                        else
                        {
                            edgeType[edgeI] = 1;
                        }
                    }
                    else if (nHigherLevels[edgeI] == 2)
                    {
                        if (isBoundaryEdge[edgeI] && nRefinementFaces[edgeI] == 2)
                        {
                            if (nLowerLevels[edgeI] == 1)
                            {
                                edgeType[edgeI] = 1;
                            }
                            else
                            {
                                edgeType[edgeI] = 0;
                            }
                        }
                        else
                        {
                            edgeType[edgeI] = 0;
                        }
                    }
                    else if (nHigherLevels[edgeI] == 1)
                    {
                        if (isBoundaryEdge[edgeI])
                        {
                            edgeType[edgeI] = 0;
                        }
                        else
                        {
                            edgeType[edgeI] = 2;
                        }
                    }
                }
            }

            syncTools::syncEdgeList
            (
                mesh,
                edgeType,
                maxEqOp<label>(),
                labelMin               // null value
            );

            boolList modEdges(mesh.nEdges(), false);
            boolList modPts(mesh.nPoints(), false);

            forAll(mesh.edges(), edgeI)
            {
                if (nRefinementFaces[edgeI] > 0 && isBoundaryEdge[edgeI])
                {
                    const labelList eCells = mesh.edgeCells()[edgeI];

                    forAll(eCells, eCI)
                    {
                        label cellI = eCells[eCI];
                        if (newCellLevel[cellI] == minEdgeLevel[edgeI])
                        {
                            labelHashSet cFaceSet(mesh.cells()[cellI]);

                            const labelList& eFaces = mesh.edgeFaces()[edgeI];
                            label nBoundary = 0;
                            label nRefined = 0;

                            vector refVec(vector::zero);
                            vector bdyVec(vector::zero);
                            forAll(eFaces, eFI)
                            {
                                label faceI = eFaces[eFI];

                                if (cFaceSet.found(faceI))
                                {
                                    label patchI = patches.whichPatch(faceI);

                                    if (refinedInternalFace.get(faceI) == 1)
                                    {
                                        refVec = mesh.faceAreas()[faceI];
                                        refVec /= mag(refVec);
                                        if (patchI == -1)
                                        {
                                            if (cellI != mesh.faceOwner()[faceI])
                                            {
                                                refVec = -refVec;
                                            }
                                        }
                                        nRefined++;
                                    }

                                    if
                                    (
                                        patchI != -1
                                        && !patches[patchI].coupled()
                                    )
                                    {
                                        bdyVec = mesh.faceAreas()[faceI];
                                        bdyVec /= mag(bdyVec);
                                        nBoundary++;
                                    }
                                }
                            }
                            if (nBoundary == 0 && nRefined == 1)
                            {
                                modEdges[edgeI] = true;
                                if (candidateCellSet.insert(cellI))
                                {
                                    typeA++;
                                    newCandidateCellSet.insert(cellI);
                                    break;
                                }
                            }
                            else if (nBoundary == 1 && nRefined == 1)
                            {
                                if ((refVec&bdyVec) > 0.8)
                                {
                                    modEdges[edgeI] = true;
                                    if (candidateCellSet.insert(cellI))
                                    {
                                        typeB++;
                                        newCandidateCellSet.insert(cellI);
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //Refine single cells connected to refined cells
            forAll(mesh.cells(), cellI)
            {
                label cLevel = cellLevel[cellI];
                if (newCellLevel[cellI] == cLevel)
                {
                    const labelList& cPts = mesh.cellPoints()[cellI];
                    const cell c =  mesh.cells()[cellI];
                    // refine cells where all anchor pts are interface pts
                    if (cPts.size() >= 18 && c.size() >= 12)
                    {
                        bool allAnchorsInterface = true;

                        forAll(cPts, cPtI)
                        {
                            label pointI = cPts[cPtI];
                            if (pointLevel[pointI] <= cLevel)
                            {
                                if (refinedInterfacePoint.get(pointI) != 1)
                                {
                                    allAnchorsInterface = false;
                                    break;
                                }
                            }
                        }

                        if
                        (
                            allAnchorsInterface
                            && candidateCellSet.insert(cellI)
                        )
                        {
                            typeC++;
                            newCandidateCellSet.insert(cellI);
                        }

                    }

                    // refine isolated cells where all anchor pts are boundary
                    // and contains 4 refined neighbours
                    if (cPts.size() == 13 && c.size() == 9)
                    {
                        bool allBoundary = true;
                        forAll(cPts, cPtI)
                        {
                            label pointI = cPts[cPtI];
                            if (pointLevel[pointI] <= cLevel)
                            {
                                if (!isBoundaryPt[pointI])
                                {
                                    allBoundary = false;
                                    break;
                                }
                            }
                        }
                        if (allBoundary)
                        {
                            label nBF = 0;
                            bool interfaceFace =  false;
                            forAll(c, cFI)
                            {
                                label faceI = c[cFI];
                                if (refinedInternalFace.get(faceI) == 1)
                                {
                                    interfaceFace = true;
                                }
                                else
                                {
                                    label patchI = patches.whichPatch(faceI);
                                    if
                                    (
                                        patchI != -1
                                        && !patches[patchI].coupled()
                                    )
                                    {
                                        nBF++;
                                    }
                                }
                            }
                            if (interfaceFace && nBF >= 4)
                            {
                                if (candidateCellSet.insert(cellI))
                                {
                                    typeD++;
                                    newCandidateCellSet.insert(cellI);
                                }
                            }
                        }
                    }
                }
            }

            forAll(patches, patchI)
            {
                if (!patches[patchI].coupled())
                {
                    const polyPatch& pp = patches[patchI];
                    label start = pp.start();

                    forAll(pp, i)
                    {
                        label faceI = start + i;
                        if (mesh.faces()[faceI].size() == 8)
                        {
                            label own = mesh.faceOwner()[faceI];
                            if (candidateCellSet.insert(own))
                            {
                                typeE++;
                                newCandidateCellSet.insert(own);
                            }
                        }
                    }
                }
            }

            //check max to min point cell level not greater than one
            forAll(mesh.points(), pointI)
            {
                const labelList& pCells = mesh.pointCells()[pointI];
                forAll(pCells, pCI)
                {
                    label cellI = pCells[pCI];
                    if ((maxPointCellLevel[pointI]-newCellLevel[cellI])>1)
                    {
                        modPts[pointI] = true;
                        if (candidateCellSet.insert(cellI))
                        {
                            typeF++;
                            newCandidateCellSet.insert(cellI);
                        }
                    }
                }
            }

            forAll(nRefinementFaces, edgeI)
            {
                if (nRefinementFaces[edgeI] > 2 && !modEdges[edgeI])
                {
                    const labelList& eCells = mesh.edgeCells()[edgeI];
                    forAll(eCells, eCI)
                    {
                        label cellI = eCells[eCI];

                        if (minEdgeLevel[edgeI] == newCellLevel[cellI])
                        {
                            modEdges[edgeI] = true;
                            if (candidateCellSet.insert(cellI))
                            {
                                typeG++;
                                newCandidateCellSet.insert(cellI);
                            }
                        }
                    }
                }
            }

            //Refine convex corner top cell when boundary and
            // cannot be grown up
            labelList nEdgeBoundaryFaces(mesh.nEdges(), 0);
            forAll(mesh.edges(), edgeI)
            {
                const labelList& eFaces = mesh.edgeFaces()[edgeI];
                forAll(eFaces, eFI)
                {
                    label faceI = eFaces[eFI];
                    if
                    (
                        isBoundaryFace[faceI]
                    )
                    {
                        nEdgeBoundaryFaces[edgeI]++;
                    }
                }
            }

            syncTools::syncEdgeList
            (
                mesh,
                nEdgeBoundaryFaces,
                plusEqOp<label>(),
                label(0)          // null value
            );

            forAll(mesh.edges(), edgeI)
            {
                if (isBoundaryEdge[edgeI])
                {
                    if
                    (
                        nRefinementFaces[edgeI] == 2
                        && nHigherLevels[edgeI] == 2 && nLowerLevels[edgeI] == 1
                    )
                    {
                        const labelList& eCells = mesh.edgeCells()[edgeI];
                        forAll(eCells, eCI)
                        {
                            label cellI = eCells[eCI];

                            if (minEdgeLevel[edgeI] == newCellLevel[cellI])
                            {
                                if (candidateCellSet.insert(cellI))
                                {
                                    newCandidateCellSet.insert(cellI);
                                }
                            }
                        }
                    }
                }
            }

            labelList nPtRefinedInterfaces(mesh.nPoints(), 0);
            forAll(mesh.points(), pointI)
            {
                const labelList& pFaces = mesh.pointFaces()[pointI];
                forAll(pFaces, pFI)
                {
                    label faceI = pFaces[pFI];
                    if
                    (
                        isMasterFace.get(faceI) == 1
                        && refinedInternalFace.get(faceI) == 1
                    )
                    {
                        nPtRefinedInterfaces[pointI]++;
                    }
                }
            }

            syncTools::syncPointList
            (
                mesh,
                nPtRefinedInterfaces,
                plusEqOp<label>(),
                label(0)               // null value
             );

            forAll(mesh.points(), pointI)
            {
                if (nPtRefinedInterfaces[pointI] == 3 && isBoundaryPt[pointI])
                {
                    const labelList& pFaces = mesh.pointFaces()[pointI];
                    forAll(pFaces, pFI)
                    {
                        label faceI = pFaces[pFI];
                        if (!refinedInternalFace.get(faceI))
                        {
                            continue;
                        }

                        face f = mesh.faces()[faceI];
                        label startPt = findIndex(f, pointI);
                        label nextPt = f[f.fcIndex(startPt)];
                        label prevPt = f[f.rcIndex(startPt)];

                        label nextEdgeI = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[pointI],
                            pointI,
                            nextPt
                        );

                        label prevEdgeI = meshTools::findEdge
                        (
                            mesh.edges(),
                            mesh.pointEdges()[pointI],
                            pointI,
                            prevPt
                        );

                        if (edgeType[nextEdgeI] == 2 && edgeType[prevEdgeI] == 2)
                        {
                            label coarseCellI = -1;
                            label own = mesh.faceOwner()[faceI];
                            if (mesh.isInternalFace(faceI))
                            {
                                label nei = mesh.faceNeighbour()[faceI];
                                if (newCellLevel[own] > newCellLevel[nei])
                                {
                                    coarseCellI = nei;
                                }
                                else
                                {
                                    coarseCellI = own;
                                }
                            }
                            else
                            {
                                label patchI = patches.whichPatch(faceI);
                                if (patches[patchI].coupled())
                                {
                                    if
                                    (
                                        newCellLevel[own]
                                        < neiLevel[faceI-mesh.nInternalFaces()]
                                    )
                                    {
                                        coarseCellI = own;
                                    }
                                }
                            }

                            if (coarseCellI != -1)
                            {
                                labelHashSet cellEdges
                                (
                                    mesh.cellEdges()[coarseCellI]
                                );
                                const labelList& pEdges =
                                    mesh.pointEdges()[pointI];

                                forAll(pEdges, pEI)
                                {
                                    label edgeI = pEdges[pEI];
                                    if (edgeI != nextEdgeI && edgeI != prevEdgeI)
                                    {
                                        if (cellEdges.found(edgeI))
                                        {
                                            if (nLowerLevels[edgeI] == 3)
                                            {
                                                continue;
                                            }

                                            if (nEdgeBoundaryFaces[edgeI] < 4)
                                            {
                                                if
                                                (
                                                    candidateCellSet.insert
                                                    (coarseCellI)
                                                )
                                                {
                                                    typeH++;
                                                    newCandidateCellSet.insert
                                                        (coarseCellI);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            labelList nConcave(mesh.nPoints(), 0);
            labelList nConvex(mesh.nPoints(), 0);
            forAll(mesh.points(), pointI)
            {
                const labelList& pEdges = mesh.pointEdges()[pointI];
                forAll(pEdges, pEI)
                {
                    label edgeI = pEdges[pEI];
                    if (isMasterEdge.get(edgeI) == 1)
                    {
                        if (edgeType[edgeI] == 1)
                        {
                            nConcave[pointI]++;
                        }
                        else if (edgeType[edgeI] == 2)
                        {
                            nConvex[pointI]++;
                        }
                    }
                }
            }
            syncTools::syncPointList
            (
                mesh,
                nConcave,
                plusEqOp<label>(),
                label(0)               // null value
             );
            syncTools::syncPointList
            (
                mesh,
                nConvex,
                plusEqOp<label>(),
                label(0)               // null value
            );

            forAll(mesh.points(), pointI)
            {
                if
                (
                    refinedInterfacePoint.get(pointI) == 1 && !modPts[pointI]
                    && nConvex[pointI] == 1 && nConcave[pointI] == 0
                    && isBoundaryPt[pointI]
                )
                {
                    const labelList& pCells = mesh.pointCells()[pointI];

                    forAll(pCells, pCI)
                    {
                        label cellI = pCells[pCI];
                        if (newCellLevel[cellI]<maxPointCellLevel[pointI])
                        {
                            labelHashSet cellEdges(mesh.cellEdges()[cellI]);
                            label nBoundaryEdges = 0;
                            label nConvexEdges = 0;
                            const labelList& pEdges = mesh.pointEdges()[pointI];

                            forAll(pEdges, pEI)
                            {
                                label edgeI = pEdges[pEI];
                                if (cellEdges.found(edgeI))
                                {
                                    if (isBoundaryEdge[edgeI])
                                    {
                                        nBoundaryEdges++;
                                    }
                                    if (edgeType[edgeI] == 2)
                                    {
                                        nConvexEdges++;
                                    }
                                }
                            }

                            if (nConvexEdges == 1 && nBoundaryEdges == 2)
                            {
                                const cell cFaces = mesh.cells()[cellI];
                                labelHashSet cFacesSet(cFaces);
                                const labelList& pFaces =
                                    mesh.pointFaces()[pointI];
                                label nInternalFaces = 0;

                                forAll(pFaces, pFI)
                                {
                                    label faceI = pFaces[pFI];
                                    if (cFacesSet.found(faceI))
                                    {
                                        label patchI =
                                            patches.whichPatch(faceI);
                                        if
                                        (
                                            patchI == -1
                                            || patches[patchI].coupled()
                                        )
                                        {
                                            nInternalFaces++;
                                        }
                                    }
                                }

                                if (nInternalFaces == 3)
                                {
                                    if (candidateCellSet.insert(cellI))
                                    {
                                        typeI++;
                                        newCandidateCellSet.insert(cellI);
                                    }
                                }
                            }
                        }
                    }
                }

                if
                (
                    nConvex[pointI] == 3 && nConcave[pointI] == 0
                    && isBoundaryPt[pointI]
                )
                {
                    const labelList& pCells = mesh.pointCells()[pointI];
                    label maxPtLevel = maxPointCellLevel[pointI];

                    forAll(pCells, pCI)
                    {
                        label cellI = pCells[pCI];
                        if (newCellLevel[cellI] < maxPtLevel)
                        {
                            if (candidateCellSet.insert(cellI))
                            {
                                typeJ++;
                                newCandidateCellSet.insert(cellI);
                            }
                        }
                    }
                }

                if
                (
                    refinedInterfacePoint.get(pointI) == 1 && !modPts[pointI]
                    && nConvex[pointI] == 0 && nConcave[pointI] == 1
                )
                {
                    const labelList& pCells = mesh.pointCells()[pointI];
                    forAll(pCells, pCI)
                    {
                        label cellI = pCells[pCI];
                        if (newCellLevel[cellI]<maxPointCellLevel[pointI])
                        {
                            labelHashSet cellEdges(mesh.cellEdges()[cellI]);
                            label nBoundaryEdges = 0;
                            const labelList& pEdges = mesh.pointEdges()[pointI];

                            forAll(pEdges, pEI)
                            {
                                label edgeI = pEdges[pEI];
                                if
                                (
                                    isBoundaryEdge[edgeI]
                                    && cellEdges.found(edgeI)
                                )
                                {
                                    nBoundaryEdges++;
                                }
                            }

                            if (nBoundaryEdges == 3)
                            {
                                labelHashSet cellFaces(mesh.cells()[cellI]);
                                const labelList& pFaces =
                                    mesh.pointFaces()[pointI];
                                bool bFace = false;

                                forAll(pFaces, pFI)
                                {
                                    label faceI = pFaces[pFI];
                                    if
                                    (
                                        cellFaces.found(faceI)
                                        && isBoundaryFace[faceI]
                                    )
                                    {
                                        bFace = true;
                                        break;
                                    }
                                }

                                if (!bFace && candidateCellSet.insert(cellI))
                                {
                                    typeK++;
                                    newCandidateCellSet.insert(cellI);
                                }
                            }
                            else if (addedRefine && nBoundaryEdges == 2)
                            {
                                const labelList& cPts =
                                    mesh.cellPoints()[cellI];
                                label cLevel = cellLevel[cellI];
                                label nCellBoundaryPts = 0;
                                forAll(cPts, cptI)
                                {
                                    label cellMeshPt = cPts[cptI];
                                    if
                                    (
                                        isBoundaryPt[cellMeshPt]
                                        && pointLevel[cellMeshPt] <= cLevel
                                    )
                                    {
                                        nCellBoundaryPts++;
                                    }
                                }
                                if (nCellBoundaryPts > 6)
                                {

                                    if (candidateCellSet.insert(cellI))
                                    {
                                        typeL++;
                                        newCandidateCellSet.insert(cellI);
                                    }

                                }
                            }
                        }
                    }
                }
            }

            //Check for topologically valid connections for adding layers
            //and refine out non-valid ones
            forAll(mesh.points(), pointI)
            {
                if (refinedInterfacePoint.get(pointI) == 1 && !modPts[pointI])
                {
                    if
                    (
                        (nConvex[pointI] == 0 && nConcave[pointI] == 0)
                        || (nConvex[pointI] == 0 && nConcave[pointI] == 1)
                        || (nConvex[pointI] == 0 && nConcave[pointI] == 2)
                        || (nConvex[pointI] == 0 && nConcave[pointI] == 3)
                        || (nConvex[pointI] == 1 && nConcave[pointI] == 0)
                        || (nConvex[pointI] == 1 && nConcave[pointI] == 1)
                        || (nConvex[pointI] == 1 && nConcave[pointI] == 2)
                        || (nConvex[pointI] == 2 && nConcave[pointI] == 0)
                        || (nConvex[pointI] == 2 && nConcave[pointI] == 1)
                        || (nConvex[pointI] == 3 && nConcave[pointI] == 0)
                        || (nConvex[pointI] == 3 && nConcave[pointI] == 3)
                     )
                    {
                        continue;
                    }

                    label maxPtLevel = maxPointCellLevel[pointI];
                    const labelList& pCells = mesh.pointCells()[pointI];
                    forAll(pCells, pCI)
                    {
                        label cellI = pCells[pCI];
                        if (newCellLevel[cellI] < maxPtLevel)
                        {
                            if (candidateCellSet.insert(cellI))
                            {
                                typeM++;
                                newCandidateCellSet.insert(cellI);
                            }
                        }
                    }
                }
            }

            if (debug)
            {
                Info<<"Type A: "<<returnReduce(typeA, sumOp<label>())<<endl;
                Info<<"Type B: "<<returnReduce(typeB, sumOp<label>())<<endl;
                Info<<"Type C: "<<returnReduce(typeC, sumOp<label>())<<endl;
                Info<<"Type D: "<<returnReduce(typeD, sumOp<label>())<<endl;
                Info<<"Type E: "<<returnReduce(typeE, sumOp<label>())<<endl;
                Info<<"Type F: "<<returnReduce(typeF, sumOp<label>())<<endl;
                Info<<"Type G: "<<returnReduce(typeG, sumOp<label>())<<endl;
                Info<<"Type H: "<<returnReduce(typeH, sumOp<label>())<<endl;
                Info<<"Type I: "<<returnReduce(typeI, sumOp<label>())<<endl;
                Info<<"Type J: "<<returnReduce(typeJ, sumOp<label>())<<endl;
                Info<<"Type K: "<<returnReduce(typeK, sumOp<label>())<<endl;
                Info<<"Type L: "<<returnReduce(typeL, sumOp<label>())<<endl;
                Info<<"Type M: "<<returnReduce(typeM, sumOp<label>())<<endl;
            }

            labelList newCandidateCells = newCandidateCellSet.toc();

            if (returnReduce(newCandidateCells.size(), sumOp<label>()) == 0)
            {
                break;
            }
        }

        cellsToRefine = candidateCellSet.toc();

        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;

        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
            mesh.write();
            cellSet refinedCellSet
            (
                mesh,
                "refinedCells",
                labelHashSet(cellsToRefine)
            );
            refinedCellSet.write();
        }

        if (!meshRefiner_.checkRefinedAndRemove(cellsToRefine,false,false))
        {
            //Refine cells
            autoPtr<mapPolyMesh> map = meshRefiner_.refine
            (
                cellsToRefine,
                refineParams.locationsInMesh()
            );

            label nChanges = 0;

            //Baffle holes
            nChanges += meshRefiner_.baffleHoles();

            //Remove cells that might cause snapping issues
            if (removeUnsplittable)
            {
                nChanges += meshRefiner_.removeProblemCells();
            }

            //Remove cells that might cause dualisation issues
            nChanges +=
                meshRefiner_.removeProblemDualisationCells
                (
                    noLayerPatchIDs,
                    removeUnsplittable
                );

            nChanges += meshRefiner_.removeHoleCells();

            reduce(nChanges, sumOp<label>());

            // Stop when no cells to refine. After a few iterations check if too
            // few cells

            if (nCellsToRefine == 0 && nChanges == 0)
            {
                Info<< "Stopping refining interface cells"
                    << nl << endl;
                break;
            }
        }
    }

    return iter;
}


// Detect cells with opposing intersected faces of differing refinement
// level and refine them.
Foam::label Foam::snappyRefineDriver::refinementInterfaceRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    const fvMesh& mesh = meshRefiner_.mesh();

    label iter = 0;

    if (refineParams.interfaceRefine())
    {
        for (;iter < maxIter; iter++)
        {
            Info<< nl
                << "Refinement transition refinement iteration " << iter << nl
                << "--------------------------------------------" << nl
                << endl;

            const labelList& surfaceIndex = meshRefiner_.surfaceIndex();
            const hexRef8& cutter = meshRefiner_.meshCutter();
            const vectorField& fA = mesh.faceAreas();
            const labelList& faceOwner = mesh.faceOwner();


            // Determine cells to refine
            // ~~~~~~~~~~~~~~~~~~~~~~~~~

            const cellList& cells = mesh.cells();

            labelList candidateCells;
            {
                // Pass1: pick up cells with differing face level

                cellSet transitionCells
                (
                    mesh,
                    "transitionCells",
                    cells.size()/100
                );

                forAll(cells, celli)
                {
                    const cell& cFaces = cells[celli];
                    label cLevel = cutter.cellLevel()[celli];

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        if (surfaceIndex[facei] != -1)
                        {
                            label fLevel = cutter.faceLevel(facei);
                            if (fLevel != cLevel)
                            {
                                transitionCells.insert(celli);
                            }
                        }
                    }
                }


                cellSet candidateCellSet
                (
                    mesh,
                    "candidateCells",
                    cells.size()/1000
                );

                // Pass2: check for oppositeness

                //forAllConstIter(cellSet, transitionCells, iter)
                //{
                //    label celli = iter.key();
                //    const cell& cFaces = cells[celli];
                //    const point& cc = cellCentres[celli];
                //    const scalar rCVol = pow(cellVolumes[celli], -5.0/3.0);
                //
                //    // Determine principal axes of cell
                //    symmTensor R(Zero);
                //
                //    forAll(cFaces, i)
                //    {
                //        label facei = cFaces[i];
                //
                //        const point& fc = faceCentres[facei];
                //
                //        // Calculate face-pyramid volume
                //        scalar pyrVol = 1.0/3.0 * fA[facei] & (fc-cc);
                //
                //        if (faceOwner[facei] != celli)
                //        {
                //            pyrVol = -pyrVol;
                //        }
                //
                //        // Calculate face-pyramid centre
                //        vector pc = (3.0/4.0)*fc + (1.0/4.0)*cc;
                //
                //        R += pyrVol*sqr(pc-cc)*rCVol;
                //    }
                //
                //    //- MEJ: Problem: truncation errors cause complex evs
                //    vector lambdas(eigenValues(R));
                //    const tensor axes(eigenVectors(R, lambdas));
                //
                //
                //    // Check if this cell has
                //    // - opposing sides intersected
                //    // - which are of different refinement level
                //    // - plus the inbetween face
                //
                //    labelVector plusFaceLevel(labelVector(-1, -1, -1));
                //    labelVector minFaceLevel(labelVector(-1, -1, -1));
                //
                //    forAll(cFaces, cFacei)
                //    {
                //        label facei = cFaces[cFacei];
                //
                //        if (surfaceIndex[facei] != -1)
                //        {
                //            label fLevel = cutter.faceLevel(facei);
                //
                //            // Get outwards pointing normal
                //            vector n = fA[facei]/mag(fA[facei]);
                //            if (faceOwner[facei] != celli)
                //            {
                //                n = -n;
                //            }
                //
                //            // What is major direction and sign
                //            direction cmpt = vector::X;
                //            scalar maxComp = (n&axes.x());
                //
                //            scalar yComp = (n&axes.y());
                //            scalar zComp = (n&axes.z());
                //
                //            if (mag(yComp) > mag(maxComp))
                //            {
                //                maxComp = yComp;
                //                cmpt = vector::Y;
                //            }
                //
                //            if (mag(zComp) > mag(maxComp))
                //            {
                //                maxComp = zComp;
                //                cmpt = vector::Z;
                //            }
                //
                //            if (maxComp > 0)
                //            {
                //                plusFaceLevel[cmpt] = max
                //                (
                //                    plusFaceLevel[cmpt],
                //                    fLevel
                //                );
                //            }
                //            else
                //            {
                //                minFaceLevel[cmpt] = max
                //                (
                //                    minFaceLevel[cmpt],
                //                    fLevel
                //                );
                //            }
                //        }
                //    }
                //
                //    // Check if we picked up any opposite differing level
                //    for (direction dir = 0; dir < vector::nComponents; dir++)
                //    {
                //        if
                //        (
                //            plusFaceLevel[dir] != -1
                //         && minFaceLevel[dir] != -1
                //         && plusFaceLevel[dir] != minFaceLevel[dir]
                //        )
                //        {
                //            candidateCellSet.insert(celli);
                //        }
                //    }
                //}

                const scalar oppositeCos = Foam::cos(Foam::degToRad(135));

                forAllConstIter(cellSet, transitionCells, iter)
                {
                    label celli = iter.key();
                    const cell& cFaces = cells[celli];
                    label cLevel = cutter.cellLevel()[celli];

                    // Detect opposite intersection
                    bool foundOpposite = false;

                    forAll(cFaces, cFacei)
                    {
                        label facei = cFaces[cFacei];

                        if
                        (
                            surfaceIndex[facei] != -1
                         && cutter.faceLevel(facei) > cLevel
                        )
                        {
                            // Get outwards pointing normal
                            vector n = fA[facei]/mag(fA[facei]);
                            if (faceOwner[facei] != celli)
                            {
                                n = -n;
                            }

                            // Check for any opposite intersection
                            forAll(cFaces, cFaceI2)
                            {
                                label face2i = cFaces[cFaceI2];

                                if
                                (
                                    face2i != facei
                                 && surfaceIndex[face2i] != -1
                                )
                                {
                                    // Get outwards pointing normal
                                    vector n2 = fA[face2i]/mag(fA[face2i]);
                                    if (faceOwner[face2i] != celli)
                                    {
                                        n2 = -n2;
                                    }


                                    if ((n&n2) < oppositeCos)
                                    {
                                        foundOpposite = true;
                                        break;
                                    }
                                }
                            }

                            if (foundOpposite)
                            {
                                break;
                            }
                        }
                    }


                    if (foundOpposite)
                    {
                        candidateCellSet.insert(celli);
                    }
                }

                if (debug&meshRefinement::MESH)
                {
                    Pout<< "Dumping " << candidateCellSet.size()
                        << " cells to cellSet candidateCellSet." << endl;
                    candidateCellSet.instance() = meshRefiner_.timeName();
                    candidateCellSet.write();
                }
                candidateCells = candidateCellSet.toc();
            }



            labelList cellsToRefine
            (
                meshRefiner_.meshCutter().consistentRefinement
                (
                    candidateCells,
                    true
                )
            );
            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;


            label nCellsToRefine = cellsToRefine.size();
            reduce(nCellsToRefine, sumOp<label>());

            Info<< "Selected for refinement : " << nCellsToRefine
                << " cells (out of " << mesh.globalData().nTotalCells()
                << ')' << endl;

            // Stop when no cells to refine. After a few iterations check if too
            // few cells
            if
            (
                nCellsToRefine == 0
             || (
                    iter >= 1
                 && nCellsToRefine <= refineParams.minRefineCells()
                )
            )
            {
                Info<< "Stopping refining since too few cells selected."
                    << nl << endl;
                break;
            }


            if (debug)
            {
                const_cast<Time&>(mesh.time())++;
            }


            if
            (
                returnReduce
                (
                    (mesh.nCells() >= refineParams.maxLocalCells()),
                    orOp<bool>()
                )
            )
            {
                meshRefiner_.balanceAndRefine
                (
                    "interface cell refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams
                );
            }
            else
            {
                meshRefiner_.refineAndBalance
                (
                    "interface cell refinement iteration " + name(iter),
                    decomposer_,
                    distributor_,
                    cellsToRefine,
                    refineParams
                );
            }
        }
    }
    return iter;
}


void Foam::snappyRefineDriver::removeInsideCells
(
    const refinementParameters& refineParams,
    const label nBufferLayers
)
{
    Info<< nl
        << "Removing mesh beyond surface intersections" << nl
        << "------------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (debug)
    {
       const_cast<Time&>(mesh.time())++;
    }

    meshRefiner_.splitMesh
    (
        nBufferLayers,                  // nBufferLayers
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams,
        setFormatter_
    );

    if (debug&meshRefinement::MESH)
    {
        const_cast<Time&>(mesh.time())++;

        Pout<< "Writing subsetted mesh to time "
            << meshRefiner_.timeName() << endl;
        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner_.timeName()
        );
        Pout<< "Dumped mesh in = "
            << mesh.time().cpuTimeIncrement() << " s\n" << nl << endl;
    }
}


void Foam::snappyRefineDriver::removeOutsideCells
(
    const refinementParameters& refineParams
)
{
    Info<< nl
        << "Removing mesh beyond surface intersections" << nl
        << "------------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (debug)
    {
       const_cast<Time&>(mesh.time())++;
    }

    meshRefiner_.splitMesh(refineParams);
}


Foam::label Foam::snappyRefineDriver::shellRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    if (controller_.mode() == meshControl::DRYRUN)
    {
        return 0;
    }

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(shell, "snappyHexMesh::refine::shell");
    #endif

    const fvMesh& mesh = meshRefiner_.mesh();

    // Mark current boundary faces with 0. Have meshRefiner maintain them.
    meshRefiner_.userFaceData().setSize(1);

    // mark list to remove any refined faces
    meshRefiner_.userFaceData()[0].first() = meshRefinement::REMOVE;
    meshRefiner_.userFaceData()[0].second() = createWithValues<labelList>
    (
        mesh.nFaces(),
        -1,
        meshRefiner_.intersectedFaces(),
        0
    );

    // Determine the maximum refinement level over all volume refinement
    // regions. This determines the minumum number of shell refinement
    // iterations.
    label overallMaxShellLevel = meshRefiner_.shells().maxIsoLevel();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Shell refinement iteration " << iter << nl
            << "----------------------------" << nl
            << endl;

        labelList candidateCells
        (
            meshRefiner_.refineCandidates
            (
                iter,
                refineParams,

                false,              // initialRefinement
                false,              // featureRefinement
                true,               // featureDistanceRefinement
                true,               // internalRefinement
                false,              // surfaceRefinement
                false,              // curvatureRefinement
                false,              // proximity refinement
                false,              // smallFeatureRefinement
                false,              // gapRefinement
                false,              // bigGapRefinement
                false               // spreadGapSize
            )
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Dumping " << candidateCells.size()
                << " cells to cellSet candidateCellsFromShells." << endl;

            cellSet c(mesh, "candidateCellsFromShells", candidateCells);
            c.instance() = meshRefiner_.timeName();
            c.write();
        }

        // Problem choosing starting faces for bufferlayers (bFaces)
        //  - we can't use the current intersected boundary faces
        //    (intersectedFaces) since this grows indefinitely
        //  - if we use 0 faces we don't satisfy bufferLayers from the
        //    surface.
        //  - possibly we want to have bFaces only the initial set of faces
        //    and maintain the list while doing the refinement.
        labelList bFaces
        (
            findIndices(meshRefiner_.userFaceData()[0].second(), 0)
        );

        //Info<< "Collected boundary faces : "
        //    << returnReduce(bFaces.size(), sumOp<label>()) << endl;

        labelList cellsToRefine;

        if (refineParams.nBufferLayers() <= 2)
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentSlowRefinement
            (
                refineParams.nBufferLayers(),
                candidateCells,                     // cells to refine
                bFaces,                             // faces for nBufferLayers
                1,                                  // point difference
                meshRefiner_.intersectedPoints()    // points to check
            );
        }
        else
        {
            cellsToRefine = meshRefiner_.meshCutter().consistentSlowRefinement2
            (
                refineParams.nBufferLayers(),
                candidateCells,                 // cells to refine
                bFaces                          // faces for nBufferLayers
            );
        }

        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for internal refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0
         || (
                iter >= overallMaxShellLevel
             && nCellsToRefine <= refineParams.minRefineCells()
            )
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "shell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "coarse cell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }

    return iter;
}


Foam::label Foam::snappyRefineDriver::refineBoundaryCells
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();

    label iter;
    for (iter = 0; iter < maxIter; iter++)
    {
        Info<< nl
            << "Problem boundary cell refinement iteration " << iter << nl
            << "--------------------------------------------" << nl
            << endl;

        labelList maxPointCellLevel(mesh.nPoints(), labelMin);
        labelList minPointCellLevel(mesh.nPoints(), labelMax);

        forAll(mesh.points(), pointI)
        {
            const labelList& pCells = mesh.pointCells()[pointI];

            forAll(pCells, pCI)
            {
                label cellI = pCells[pCI];
                label cLevel = cellLevel[cellI];

                if (cLevel > maxPointCellLevel[pointI])
                {
                    maxPointCellLevel[pointI] = cLevel;
                }
                if (cLevel < minPointCellLevel[pointI])
                {
                    minPointCellLevel[pointI] = cLevel;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            maxPointCellLevel,
            maxEqOp<label>(),
            labelMin              // null value
        );

        syncTools::syncPointList
        (
            mesh,
            minPointCellLevel,
            minEqOp<label>(),
            labelMax              // null value
        );


        labelList nLowerCells(mesh.nPoints(), 0);
        labelList nHigherCells(mesh.nPoints(), 0);
        labelList nInterfaceFaces(mesh.nPoints(), 0);
        labelList nPointBoundFaces(mesh.nPoints(), 0);

        forAll(mesh.points(), pointI)
        {
            const labelList& pCells = mesh.pointCells()[pointI];

            label maxPLevel = maxPointCellLevel[pointI];
            label minPLevel = minPointCellLevel[pointI];

            forAll(pCells, pCI)
            {
                label cellI = pCells[pCI];
                label cLevel = cellLevel[cellI];

                if (cLevel == maxPLevel)
                {
                    nHigherCells[pointI]++;
                }

                if (cLevel == minPLevel)
                {
                    nLowerCells[pointI]++;
                }
            }
        }

        syncTools::syncPointList
        (
            mesh,
            nLowerCells,
            plusEqOp<label>(),
            label(0)      // null value
        );

        syncTools::syncPointList
        (
            mesh,
            nHigherCells,
            plusEqOp<label>(),
            label(0)      // null value
        );

        boolList markedFace(mesh.nFaces(), false);
//        label nMarked = 0;
        forAll(mesh.points(), pointI)
        {
            if (nHigherCells[pointI] == 1 && nLowerCells[pointI] == 4)
            {
                const labelList& pCells = mesh.pointCells()[pointI];
                label higherCell = -1;
                label maxPLevel = maxPointCellLevel[pointI];
                forAll(pCells, i)
                {
                    label cellI = pCells[i];
                    if (cellLevel[cellI] == maxPLevel)
                    {
                        higherCell = cellI;
                        break;
                    }
                }

                if (higherCell != -1)
                {
                    const labelList& pFaces = mesh.pointFaces()[pointI];
                    labelHashSet cFaces(mesh.cells()[higherCell]);

                    forAll(pFaces, pFI)
                    {
                        label faceI = pFaces[pFI];

                        label patchI = patches.whichPatch(faceI);
                        if
                        (
                            cFaces.found(faceI)
                            && (patchI == -1 || patches[patchI].coupled())
                        )
                        {
                            markedFace[faceI] = true;
//                            nMarked++;
                        }
                    }
                }
            }
        }

        // Calculate cell level
        scalarField neiLevel(mesh.nFaces()-mesh.nInternalFaces());

        for
        (
            label faceI = mesh.nInternalFaces();
            faceI < mesh.nFaces();
            faceI++
         )
        {
            neiLevel[faceI-mesh.nInternalFaces()] =
                cellLevel[mesh.faceOwner()[faceI]];
        }
        syncTools::swapBoundaryFaceList(mesh, neiLevel);

        syncTools::syncFaceList
        (
            mesh,
            markedFace,
            orEqOp<bool>()
         );

        labelList candidateCells;
        cellSet candidateCellSet(mesh, "candidateCells", mesh.nCells()/1000);

        forAll(mesh.faces(), faceI)
        {
            if (markedFace[faceI])
            {
                label own = mesh.faceOwner()[faceI];

                if (mesh.isInternalFace(faceI))
                {
                    label nbr = mesh.faceNeighbour()[faceI];

                    if (cellLevel[own] < cellLevel[nbr])
                    {
                        candidateCellSet.insert(own);
                    }
                    else
                    {
                        candidateCellSet.insert(nbr);
                    }
                }
                else
                {
                    label patchI = patches.whichPatch(faceI);
                    if (patches[patchI].coupled())
                    {
                        if
                        (
                            cellLevel[own]
                            < neiLevel[faceI-mesh.nInternalFaces()]
                         )
                        {
                            candidateCellSet.insert(own);
                        }
                    }
                }
            }
        }

        candidateCells = candidateCellSet.toc();

        labelList cellsToRefine
        (
            meshRefiner_.meshCutter().consistentRefinement
            (
                candidateCells,
                true
            )
        );
        Info<< "Determined cells to refine in = "
            << mesh.time().cpuTimeIncrement() << " s" << endl;


        label nCellsToRefine = cellsToRefine.size();
        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine. After a few iterations check if too
        // few cells
        if
        (
            nCellsToRefine == 0
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }


        if
        (
            returnReduce
            (
                (mesh.nCells() >= refineParams.maxLocalCells()),
                orOp<bool>()
            )
        )
        {
            meshRefiner_.balanceAndRefine
            (
                "coarse cell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
        else
        {
            meshRefiner_.refineAndBalance
            (
                "coarse cell refinement iteration " + name(iter),
                decomposer_,
                distributor_,
                cellsToRefine,
                refineParams
            );
        }
    }

    return iter;
}


void Foam::snappyRefineDriver::baffleAndSplitMesh
(
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool handleSnapProblems,
    const dictionary& motionDict,
    const bool threaded /*= false*/
)
{
    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(split, "snappyHexMesh::refine::splitting");
    #endif

    Info<< nl
        << "Splitting mesh at surface intersections" << nl
        << "---------------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    if (debug)
    {
       const_cast<Time&>(mesh.time())++;
    }

    // Introduce baffles at surface intersections. Note:
    // meshRefiment::surfaceIndex() will
    // be like boundary face from now on so not coupled anymore.
    meshRefiner_.baffleAndSplitMesh
    (
        handleSnapProblems,             // detect&remove potential snap problem

        // Snap problem cell detection
        snapParams,
        false,                          // perpendicular edge connected cells
        scalarField(0),                 // per region perpendicular angle

        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams,
        setFormatter_,
        threaded
    );

    if (!handleSnapProblems) // merge free standing baffles?
    {
        meshRefiner_.mergeFreeStandingBaffles
        (
            snapParams,
            refineParams.useTopologicalSnapDetection(),
            false,                  // perpendicular edge connected cells
            scalarField(0),         // per region perpendicular angle
            refineParams.planarAngle(),
            motionDict,
            const_cast<Time&>(mesh.time()),
            globalToMasterPatch_,
            globalToSlavePatch_,
            refineParams,
            setFormatter_,
            threaded
        );
    }
}


void Foam::snappyRefineDriver::removeBoundaryZoneProblemCells
(
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const dictionary& motionDict,
    List<labelPair>& baffles
)
{
    const fvMesh& mesh = meshRefiner_.mesh();

    meshRefiner_.handleZoneSnapProblems
    (
        snapParams,
        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams,
        baffles
    );
}


// Zoning method used after dual mesh construction
void Foam::snappyRefineDriver::zonifyDual
(
    const refinementParameters& refineParams,
    const bool updateIntersections
)
{
    bool foundZone  = false;
    const wordList& zonesInMesh = refineParams.zonesInMesh();
    forAll(zonesInMesh, locationI)
    {
        if (zonesInMesh[locationI] != "none")
        {
            foundZone = true;
            break;
        }
    }

    const labelList nonBoundaryNamedSurfaces =
        surfaceZonesInfo::getNonBoundaryNamedSurfaces
        (meshRefiner_.surfaces().surfZones());

    if (foundZone || nonBoundaryNamedSurfaces.size())
    {
        const fvMesh& mesh = meshRefiner_.mesh();

        //Need to update list of intersected faces after dual mesh constructed.
        // Currently setting to list of all faces
        if (updateIntersections)
        {
            meshRefiner_.surfaceIndex() = identity(mesh.nFaces());
        }

        // Mesh is at its finest. Do optional zoning (cellZones and faceZones)
        wordPairHashTable zonesToFaceZone;


        Info<< nl
            << "Introducing zones for interfaces" << nl
            << "--------------------------------" << nl
            << endl;


        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        meshRefiner_.zonify
        (
            refineParams,
            setFormatter_,
            zonesToFaceZone,
            true,// whether to use already set zones
            true, // whether to reset unset zones
            false
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing zoned mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner_.timeName()
            );
        }

        // Check that all faces are synced
        meshRefinement::checkCoupledFaceZones(mesh);

        // Create pairs of patches for faceZones
        {
            HashTable<Pair<word>> faceZoneToPatches(zonesToFaceZone.size());

            //    Note: zonesToFaceZone contains the same data on different
            //          processors but in different order. We could sort the
            //          contents but instead just loop in sortedToc order.
            List<Pair<word>> czs(zonesToFaceZone.sortedToc());

            forAll(czs, i)
            {
                const Pair<word>& czNames = czs[i];
                const word& fzName = zonesToFaceZone[czNames];
                const word& masterName = fzName;
                const word slaveName = czNames.second() + "_to_"
                    + czNames.first();
                Pair<word> patches(masterName, slaveName);
                faceZoneToPatches.insert(fzName, patches);
            }
            addFaceZones(meshRefiner_, refineParams, faceZoneToPatches);
        }

        for (label iter = 0; iter < 4; iter++)
        {
            label nCellsSwitched = meshRefiner_.removeZoneHoles(true);
            nCellsSwitched += meshRefiner_.removeZoneHoles(false);
            if (nCellsSwitched == 0)
            {
                break;
            }
        }
    }
}


void Foam::snappyRefineDriver::zonify
(
    const refinementParameters& refineParams,
    wordPairHashTable& zonesToFaceZone
)
{
    // Mesh is at its finest. Do zoning
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This puts all faces with intersection across a zoneable surface
    // into that surface's faceZone. All cells inside faceZone get given the
    // same cellZone.

    const labelList namedSurfaces =
        surfaceZonesInfo::getNamedSurfaces(meshRefiner_.surfaces().surfZones());

    if
    (
        namedSurfaces.size()
     || refineParams.zonesInMesh().size()
    )
    {
        Info<< nl
            << "Introducing zones for interfaces" << nl
            << "--------------------------------" << nl
            << endl;

        const fvMesh& mesh = meshRefiner_.mesh();

        if (debug)
        {
            const_cast<Time&>(mesh.time())++;
        }

        meshRefiner_.zonify
        (
            refineParams,
            setFormatter_,
            zonesToFaceZone,
            false,// whether to use already set zones
            false, // whether to reset unset zones
            (controller_.algorithm() == meshControl::EXTRUDE) //remove unset
        );

        if (debug&meshRefinement::MESH)
        {
            Pout<< "Writing zoned mesh to time "
                << meshRefiner_.timeName() << endl;
            meshRefiner_.write
            (
                meshRefinement::debugType(debug),
                meshRefinement::writeType
                (
                    meshRefinement::writeLevel()
                  | meshRefinement::WRITEMESH
                ),
                mesh.time().path()/meshRefiner_.timeName()
            );
        }

        // Check that all faces are synced
        meshRefinement::checkCoupledFaceZones(mesh);
    }
}


void Foam::snappyRefineDriver::splitAndMergeBaffles
(
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool handleSnapProblems,
    const dictionary& motionDict
)
{
    Info<< nl
        << "Handling cells with snap problems" << nl
        << "---------------------------------" << nl
        << endl;

    const fvMesh& mesh = meshRefiner_.mesh();

    // Introduce baffles and split mesh
    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    const scalarField& perpAngle = meshRefiner_.surfaces().perpendicularAngle();

    meshRefiner_.baffleAndSplitMesh
    (
        handleSnapProblems,

        // Snap problem cell detection
        snapParams,
        handleSnapProblems,                 // remove perp edge connected cells
        perpAngle,                          // perp angle

        motionDict,
        const_cast<Time&>(mesh.time()),
        globalToMasterPatch_,
        globalToSlavePatch_,
        refineParams,
        setFormatter_
    );

    if (controller_.algorithm() != meshControl::EXTRUDE)
    {
        // Merge free-standing baffles always
        meshRefiner_.mergeFreeStandingBaffles
        (
            snapParams,
            refineParams.useTopologicalSnapDetection(),
            handleSnapProblems,
            perpAngle,
            refineParams.planarAngle(),
            motionDict,
            const_cast<Time&>(mesh.time()),
            globalToMasterPatch_,
            globalToSlavePatch_,
            refineParams,
            setFormatter_
        );
    }

    if (debug)
    {
        const_cast<Time&>(mesh.time())++;
    }

    if (refineParams.mergeSingleCouples())
    {
        meshRefiner_.mergeSingleCouples();
    }

    // Duplicate points on baffles that are on more than one cell
    // region. This will help snapping pull them to separate surfaces.
    meshRefiner_.dupNonManifoldPoints();

    // Merge all baffles that are still remaining after duplicating points.
    List<labelPair> couples(localPointRegion::findDuplicateFacePairs(mesh));

    label nCouples = returnReduce(couples.size(), sumOp<label>());

    Info<< "Detected unsplittable baffles : " << nCouples << endl;

    if (nCouples > 0)
    {
        // Actually merge baffles. Note: not exactly parallellized. Should
        // convert baffle faces into processor faces if they resulted
        // from them.
        meshRefiner_.mergeBaffles(couples, Map<label>(0), true);

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            meshRefiner_.checkData();
        }

        // Remove any now dangling parts
        meshRefiner_.splitMeshRegions
        (
            globalToMasterPatch_,
            globalToSlavePatch_,
            refineParams.locationsInMesh(),
            refineParams.locationsOutsideMesh(),
            setFormatter_,
            refineParams.fullLeakChecks()
        );

        if (debug)
        {
            // Debug:test all is still synced across proc patches
            meshRefiner_.checkData();
        }

        Info<< "Merged free-standing baffles in = "
            << mesh.time().cpuTimeIncrement() << " s." << endl;
    }

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing handleProblemCells mesh to time "
            << meshRefiner_.timeName() << endl;
        meshRefiner_.write
        (
            meshRefinement::debugType(debug),
            meshRefinement::writeType
            (
                meshRefinement::writeLevel()
              | meshRefinement::WRITEMESH
            ),
            mesh.time().path()/meshRefiner_.timeName()
        );
    }
}


void Foam::snappyRefineDriver::addFaceZones
(
    meshRefinement& meshRefiner,
    const refinementParameters& refineParams,
    const HashTable<Pair<word>>& faceZoneToPatches
)
{
    if (faceZoneToPatches.size())
    {
        Info<< nl
            << "Adding patches for face zones" << nl
            << "-----------------------------" << nl
            << endl;

        Info<< setf(ios_base::left)
            << setw(6) << "Patch"
            << setw(20) << "Type"
            << setw(30) << "Name"
            << setw(30) << "FaceZone"
            << setw(10) << "FaceType"
            << nl
            << setw(6) << "-----"
            << setw(20) << "----"
            << setw(30) << "----"
            << setw(30) << "--------"
            << setw(10) << "--------"
            << endl;

        const polyMesh& mesh = meshRefiner.mesh();

        // Add patches for added inter-region faceZones
        forAllConstIter(HashTable<Pair<word>>, faceZoneToPatches, iter)
        {
            const word& fzName = iter.key();
            const Pair<word>& patchNames = iter();

            // Get any user-defined faceZone data
            surfaceZonesInfo::faceZoneType fzType;
            dictionary patchInfo = refineParams.getZoneInfo(fzName, fzType);

            const word& masterName = fzName;
            //const word slaveName = fzName + "_slave";
            //const word slaveName = czNames.second()+"_to_"+czNames.first();
            const word& slaveName = patchNames.second();

            label mpi = meshRefiner.addMeshedPatch(masterName, patchInfo);

            Info<< setf(ios_base::left)
                << setw(6) << mpi
                << setw(20) << mesh.boundaryMesh()[mpi].type()
                << setw(30) << masterName
                << setw(30) << fzName
                << setw(10) << surfaceZonesInfo::faceZoneTypeNames[fzType]
                << nl;


            label sli = meshRefiner.addMeshedPatch(slaveName, patchInfo);

            Info<< setf(ios_base::left)
                << setw(6) << sli
                << setw(20) << mesh.boundaryMesh()[sli].type()
                << setw(30) << slaveName
                << setw(30) << fzName
                << setw(10) << surfaceZonesInfo::faceZoneTypeNames[fzType]
                << nl;

            meshRefiner.addFaceZone(fzName, masterName, slaveName, fzType);
        }

        Info<< endl;
    }
}


void Foam::snappyRefineDriver::mergePatchFaces
(
    const bool geometricMerge,
    const dictionary& motionDict
)
{
    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(merge, "snappyHexMesh::refine::merge");
    #endif

    Info<< nl
        << "Merge refined boundary faces" << nl
        << "----------------------------" << nl
        << endl;

    if (geometricMerge)
    {
        meshRefiner_.mergePatchFacesUndo
        (
            Foam::cos(degToRad(45.0)),
            Foam::cos(degToRad(45.0)),
            -1, //max area ratio (-1 to disable)
            meshRefiner_.meshedPatches(),
            motionDict,
            true, //update intersections
            false //maintain patches
        );
    }
    else
    {
        // Still merge refined boundary faces if all four are on same patch
        meshRefiner_.mergePatchFaces
        (
            Foam::cos(degToRad(45.0)),
            Foam::cos(degToRad(45.0)),
            4,          // only merge faces split into 4
            meshRefiner_.meshedPatches()
        );
    }

    if (debug)
    {
        meshRefiner_.checkData();
    }

    meshRefiner_.mergeEdgesUndo(-1., motionDict, true);

    if (debug)
    {
        meshRefiner_.checkData();
    }
}


Foam::label Foam::snappyRefineDriver::anisoRefine
(
    const refinementParameters& refineParams,
    const label maxIter
)
{
    const fvMesh& mesh = meshRefiner_.mesh();
    const labelList& cellLevel = meshRefiner_.meshCutter().cellLevel();
    const labelList& pointLevel = meshRefiner_.meshCutter().pointLevel();

    // Mark current boundary faces with 0. Have meshRefiner maintain them.
    meshRefiner_.userFaceData().setSize(1);

    // mark list to remove any refined faces
    meshRefiner_.userFaceData()[0].first() = meshRefinement::REMOVE;
    meshRefiner_.userFaceData()[0].second() = createWithValues<labelList>
    (
        mesh.nFaces(),
        -1,
        meshRefiner_.intersectedFaces(),
    //        meshRefiner_.intersectedAndBoundaryFaces(),
        0
    );

    anisoRefiner refiner
    (
        meshRefiner_,
        cellLevel,
        pointLevel
    );

    label totalIter(0);

    for (label iter = 0; iter < maxIter; iter++)
    {
        label nCellsToRefine = 0;
        for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
        {
            Info<< nl
                << "Shell aniso refinement direction "<<cmpt<<" iteration "
                << iter << nl
                << "----------------------------------------------" << nl
                << endl;

            nCellsToRefine += refiner.setRefinement(cmpt);

            Info<< "Determined cells to refine in = "
                << mesh.time().cpuTimeIncrement() << " s" << endl;

            totalIter++;
        }

        reduce(nCellsToRefine, sumOp<label>());

        Info<< "Selected for internal refinement : " << nCellsToRefine
            << " cells (out of " << mesh.globalData().nTotalCells()
            << ')' << endl;

        // Stop when no cells to refine or have done minimum necessary
        // iterations and not enough cells to refine.
        if
        (
            nCellsToRefine == 0 ||
            nCellsToRefine <= refineParams.minRefineCells()
        )
        {
            Info<< "Stopping refining since too few cells selected."
                << nl << endl;
            break;
        }

        if (Pstream::nProcs() > 1)
        {
            scalar nIdealCells =
                mesh.globalData().nTotalCells()
                / Pstream::nProcs();

            scalar unbalance = returnReduce
            (
                scalar(mag(1.0-mesh.nCells()/nIdealCells)),
                maxOp<scalar>()
             );

            if (unbalance <= refineParams.maxLoadUnbalance())
            {
                Info<< "Skipping balancing since max unbalance " << unbalance
                    << " is less than allowable "
                    << refineParams.maxLoadUnbalance()
                    << endl;
            }
            else
            {
                autoPtr<mapDistributePolyMesh> distMap = meshRefiner_.balance
                (
                    false,
                    true,
                    false,
                    scalarField(mesh.nCells(), 1), // dummy weights
                    decomposer_,
                    distributor_
                );
                refiner.distribute(distMap);
            }
        }
    }

    meshRefiner_.userFaceData().clear();

    return totalIter;
}


void Foam::snappyRefineDriver::resetBoundaryFaceZones()
{
    fvMesh& mesh = meshRefiner_.mesh();

    label nFaceZones = mesh.faceZones().size();
    boolList freeStanding(nFaceZones, true);

    if (nFaceZones > 0)
    {
        const PtrList<surfaceZonesInfo>& surfZones =
            meshRefiner_.surfaces().surfZones();
        const labelList zonedSurfaces = surfaceZonesInfo::getNamedSurfaces
        (
            meshRefiner_.surfaces().surfZones()
        );

        //Per face zone whether free-standing
        forAll(zonedSurfaces, i)
        {
            label surfi = zonedSurfaces[i];
            word faceZoneName = surfZones[surfi].faceZoneName();
            label zonei = mesh.faceZones().findZoneID(faceZoneName);
            if (zonei != -1 && !surfZones[surfi].freeStanding())
            {
                freeStanding[zonei] = false;
            }
        }
    }

    polyTopoChange meshMod(mesh);
    label nReset = 0;

    forAll(mesh.faces(), facei)
    {
        label patchi = mesh.boundaryMesh().whichPatch(facei);
        label zonei = mesh.faceZones().whichZone(facei);
        if
        (
            zonei != -1 && !freeStanding[zonei]
            && patchi != -1 && !mesh.boundaryMesh()[patchi].coupled()
        )
        {
            label own = mesh.faceOwner()[facei];
            meshMod.modifyFace
            (
                mesh.faces()[facei], // modified face
                facei,          // label of face being modified
                own,            // owner
                -1,             // neighbour
                false,          // face flip
                patchi,   // new patch for face
                -1,         // zone for face
                false       // face flip in zone
            );
            nReset++;
        }
    }

    if (returnReduce(nReset, sumOp<label>()) != 0)
    {
        // Change the mesh (no inflation)
        autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh, false, true);

        // Update fields (problem when face added to zero sized patch)
        mesh.updateMesh(map);

        // Move mesh if in inflation mode
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }
        else
        {
            mesh.clearOut();
        }
        meshRefiner_.updateMesh(map,labelList(0));
    }

    return;
}


void Foam::snappyRefineDriver::doThinShellRefine
(
    const dictionary& meshDict,
    const refinementParameters& refineParams
)
{
    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(refine, "snappyHexMesh::refine");
    #endif

    Info<< nl
        << "Thin Shell Refinement phase" << nl
        << "---------------------------" << nl
        << endl;

    fvMesh& mesh = meshRefiner_.mesh();

    // Refine around feature edges
    label maxIter(meshRefiner_.maxCellLevel());

    const label minCellsToRefine(0);

    // Refine cells containing surface bounding box smaller than cell size
    initialRefine
    (
        refineParams,
        maxIter,    // maxIter
        minCellsToRefine // min cells to refine
    );

    // Refine around feature edges
    featureEdgeRefine
    (
        refineParams,
        maxIter,    // maxIter
        minCellsToRefine // min cells to refine
    );

    if
    (
        max(meshRefiner_.surfaces().maxGapLevel()) > 0
     || max(meshRefiner_.shells().maxGapLevel()) > 0
    )
    {
        // In case we use automatic gap level refinement do some pre-refinement
        // (fast) since it is so slow.

        // Refine based on surface
        surfaceOnlyRefine
        (
            refineParams,
            maxIter     // maxIter
        );

        // Refine cells that contain a gap
        smallFeatureRefine
        (
            refineParams,
            maxIter     // maxIter
        );
    }

    // Refine based on surface
    surfaceOnlyRefine
    (
        refineParams,
        50
    );

    label minCellsToRefineProx
    (
        returnReduce(mesh.nCells(),sumOp<label>())
    );
    minCellsToRefineProx /= 1000000;

    // Refine based on proximity of portions of the STL surface
    proximityRefine
    (
        refineParams,
        maxIter,
        minCellsToRefineProx
    );

    gapOnlyRefine
    (
        refineParams,
        maxIter     // maxIter
    );

    // Refine consistently across narrow gaps (a form of shell refinement)
    bigGapOnlyRefine
    (
        refineParams,
        true,   // spreadGapSize
        maxIter     // maxIter
    );

    // Remove cells (a certain distance) beyond surface intersections
    removeOutsideCells
    (
        refineParams
    );

    mergePatchFaces(controller_.topoChanges(), motionDict);

    // Duplicate points on baffles that are on more than one cell
    // region. This will help snapping pull them to separate surfaces.
    meshRefiner_.dupNonManifoldPoints();
}


void Foam::snappyRefineDriver::doRefine
(
    const dictionary& meshDict,
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool prepareForSnapping
)
{
    bool dryRun = false;
    if (controller_.mode() == meshControl::DRYRUN)
    {
        dryRun = true;
    }
    // refinement parameters
    const dictionary& refineDict = meshDict.subDict("castellatedMeshControls");

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    #if !defined( WIN32 ) && !defined( WIN64 )
    addProfiling(refine, "snappyHexMesh::refine");
    #endif

    Info<< nl
        << "Refinement phase" << nl
        << "----------------" << nl
        << endl;

    fvMesh& mesh = meshRefiner_.mesh();

    // Check that all the keep points are inside the mesh.
    if (dryRun)
    {
        refineParams.findCells(true, mesh, refineParams.locationsOutsideMesh());
    }

    // Refine around feature edges
    label maxIter(meshRefiner_.maxCellLevel());
    if (controller_.mode() != meshControl::FAST)
    {
        maxIter *= 2;
    }
    else if (dryRun)
    {
        maxIter = 0;
    }

    const label minCellsToRefine(0);

    // Refine cells containing surface bounding box smaller than cell size
    initialRefine
    (
        refineParams,
        maxIter,    // maxIter
        minCellsToRefine // min cells to refine
    );

    //If wrapLevel keyword set refine surface to minimum wrap level
    if (meshRefiner_.surfaces().wrapLevelSurfaces().size() > 0)
    {
        // Refine around feature edges
        featureEdgeRefine
        (
            refineParams,
            maxIter,    // maxIter
            minCellsToRefine // min cells to refine
        );

        // Refine based on surface
        surfaceOnlyRefine
        (
            refineParams,
            (dryRun ? 0 : 50)
        );

        //Activate wrapper checking
        meshRefiner_.setWrapActive(true);

        //Mark initial wrap faces
        meshRefiner_.markWrapFaces(refineParams.locationsInMesh());
    }

    // Refine around feature edges
    featureEdgeRefine
    (
        refineParams,
        maxIter,    // maxIter
        minCellsToRefine // min cells to refine
    );

    if
    (
        max(meshRefiner_.surfaces().maxGapLevel()) > 0
     || max(meshRefiner_.shells().maxGapLevel()) > 0
    )
    {
        // In case we use automatic gap level refinement do some pre-refinement
        // (fast) since it is so slow.

        // Refine based on surface
        surfaceOnlyRefine
        (
            refineParams,
            maxIter     // maxIter
        );

        // Refine cells that contain a gap
        smallFeatureRefine
        (
            refineParams,
            maxIter     // maxIter
        );
    }

    // Refine based on surface
    surfaceOnlyRefine
    (
        refineParams,
        (dryRun ? 0 : 50)
    );

    label minCellsToRefineProx
    (
        returnReduce(mesh.nCells(),sumOp<label>())
    );
    minCellsToRefineProx /= 1000000;

    // Refine based on proximity of portions of the STL surface
    proximityRefine
    (
        refineParams,
        maxIter,
        minCellsToRefineProx
    );

    if (refineParams.refineOutsideGapCell())
    {
        // Refine cell identified for removal if neighbouring kept cell
        // at higher level and in a gap
        removalInterfaceRefine(refineParams);
    }

    // Remove cells (a certain distance) beyond surface intersections
    removeInsideCells
    (
        refineParams,
        1       // nBufferLayers
    );

    gapOnlyRefine
    (
        refineParams,
        maxIter     // maxIter
    );

    // Refine consistently across narrow gaps (a form of shell refinement)
    bigGapOnlyRefine
    (
        refineParams,
        true,   // spreadGapSize
        maxIter     // maxIter
    );

    // Internal mesh refinement
    shellRefine
    (
        refineParams,
        maxIter    // maxIter
    );

    if (controller_.algorithm() != meshControl::DUAL)
    {
        //refine cells with 4 or greater refined faces
        // or all anchor points are interface
        danglingCellRefine
        (
            refineParams,
            meshRefiner_.maxCellLevel()+1
        );
    }

    if (controller_.algorithm() == meshControl::EXTRUDE)
    {
        singleLevelEdgeRefine(refineParams,maxIter);
    }

    //Try and prevent voids at refinement interfaces
    if (controller_.algorithm() == meshControl::DUAL)
    {
        meshRefiner_.checkRemovedAndRefine
        (
            refineParams,
            globalToMasterPatch_,
            5
        );
    }

    //mesh at finest create gap field if active
    meshRefiner_.setGapFieldCreation(true);

    // Introduce baffles at surface intersections. Remove sections unreachable
    // from keepPoint.
    baffleAndSplitMesh
    (
        refineParams,
        snapParams,
        prepareForSnapping,
        motionDict
    );

    // Mesh is at its finest. Do optional zoning (cellZones and faceZones)
    wordPairHashTable zonesToFaceZone;
    zonify(refineParams, zonesToFaceZone);

    // Create pairs of patches for faceZones
    {
        HashTable<Pair<word>> faceZoneToPatches(zonesToFaceZone.size());

        //    Note: zonesToFaceZone contains the same data on different
        //          processors but in different order. We could sort the
        //          contents but instead just loop in sortedToc order.
        List<Pair<word>> czs(zonesToFaceZone.sortedToc());

        forAll(czs, i)
        {
            const Pair<word>& czNames = czs[i];
            const word& fzName = zonesToFaceZone[czNames];

            const word& masterName = fzName;
            const word slaveName = czNames.second() + "_to_" + czNames.first();
            Pair<word> patches(masterName, slaveName);
            faceZoneToPatches.insert(fzName, patches);
        }
        addFaceZones(meshRefiner_, refineParams, faceZoneToPatches);
    }

    // Pull baffles apart
    splitAndMergeBaffles
    (
        refineParams,
        snapParams,
        prepareForSnapping,
        motionDict
    );

    if
    (
        controller_.algorithm() == meshControl::STANDARD
        && refineParams.splitCells() && !dryRun
    )
    {
        //move points back to pre-cut state
        mesh.movePoints(meshRefiner_.oldPoints());
        meshRefiner_.clearOldPoints();
        mesh.clearMeshPhi();
    }

    if (refineParams.wrap())
    {
        if (Pstream::parRun())
        {
            meshRefiner_.balance
            (
                false,
                true,
                false,
                scalarField(mesh.nCells(), 1), // dummy weights
                decomposer_,
                distributor_
             );
        }

        removeWrappedCells
        (
            refineParams,
            refineDict,
            true,
            motionDict
         );
    }

    if (refineParams.refineBoundaryProblemCells())
    {
        refineBoundaryCells
        (
            refineParams,
            maxIter
        );
    }

    const labelList zonesToBaffle
    (
        meshRefiner_.getZones
        (
            List<surfaceZonesInfo::faceZoneType>
            (
                1,
                surfaceZonesInfo::BOUNDARY
            )
        )
    );

    if (controller_.algorithm() == meshControl::DUAL)
    {
        meshRefiner_.createSingleCellGap();

        const dictionary& layerDict = meshDict.subDict("addLayersControls");

        //Perform initial removal of problem cells
        meshRefiner_.baffleHoles();
        meshRefiner_.checkRefinedAndRemove(identity(mesh.nCells()),false,true);

        meshRefiner_.balance
        (
            false,
            true,                           // keepZoneFaces
            false,                          // keepBaffles
            scalarField(mesh.nCells(), 1),  // cellWeights
            decomposer_,
            distributor_
        );

        // Create baffles from boundary named surfaces
        if (zonesToBaffle.size())
        {
            List<labelPair> baffles;
            labelList originatingFaceZone;

            meshRefiner_.createZoneBaffles
            (
                zonesToBaffle,
                baffles,
                originatingFaceZone,
                true //update intersections
            );

            meshRefiner_.dupNonManifoldPoints();
        }

        //The mesh should now be correctly zoned reset master regions
        meshRefiner_.updateMasterRegions();

        interfaceProblemCellRefine
        (
            refineParams,
            layerDict,
            1000  // maxIter
        );
    }
    else if (controller_.algorithm() == meshControl::EXTRUDE)
    {
        // Do final balancing. Keep zoned faces on one processor since the
        // snap phase will convert them to baffles and this only works for
        // internal faces.
        meshRefiner_.balance
        (
            false,
            true,    // keepZoneFaces
            false,   // keepBaffles
            scalarField(mesh.nCells(), 1), // cellWeights
            decomposer_,
            distributor_
        );

        if (zonesToBaffle.size())
        {
            List<labelPair> baffles(0);

            // Create baffles from boundary named surfaces
            {
                labelList originatingFaceZone;
                meshRefiner_.createZoneBaffles
                (
                    zonesToBaffle,
                    baffles,
                    originatingFaceZone
                );
            }

            meshRefiner_.dupNonManifoldPoints();
        }
        meshRefiner_.updateMasterRegions();

        //Ensure cell level does not change are a baffle edge
        refineAtBaffleEdges(refineParams);

        //Perform initial removal of problem cells
        meshRefiner_.removeExtrusionProblemCells
        (
            meshDict,
            refineParams
        );
    }

    // Do something about cells with refined faces on the boundary
    if
    (
        prepareForSnapping
        &&
        (
            controller_.algorithm() != meshControl::EXTRUDE
            || refineParams.mergePreExtrude()
        )
    )
    {
        mergePatchFaces(controller_.topoChanges(), motionDict);
    }

    // Allow an-isotropic refinement
    if (meshRefiner_.shells().hasAnisotropicRefinement())
    {
        if
        (
            controller_.algorithm() == meshControl::STANDARD
        )
        {
            if (!dryRun)
            {
               anisoRefine(refineParams,maxIter);
               meshRefiner_.removeOtherSide(refineParams,globalToMasterPatch_);
            }
        }
        else
        {
           WarningInFunction
              << "The Dual and Extrude mesh algorithm does not support"
              << " anisotropic volumetric refinement. Disabling." << endl;
        }
    }

    //Optional removal of selected cells at the end of refinement
    if (refineDict.isDict("cellRemoval"))
    {
        const dictionary cellRemovalDict = refineDict.subDict
        (
            "cellRemoval"
        );
        meshRefiner_.removeSelectedCellSets(cellRemovalDict);
    }

    if (refineParams.removeBoundaryFaceZones())
    {
        //Make sure no patch faces are in face zone list
        resetBoundaryFaceZones();
    }

    if
    (
        !dryRun && Pstream::parRun()
        && controller_.algorithm() == meshControl::STANDARD
    )
    {
        Info<< nl
            << "Doing final balancing" << nl
            << "---------------------" << nl
            << endl;

        scalarField weightField(mesh.nCells(), 1);

        scalar snapWeights = refineParams.snapWeights();
        if (snapWeights > 0)
        {
            meshRefiner_.calcSnapWeights
            (
                snapWeights,
                weightField
            );
        }

        // Do final balancing. Keep zoned faces on one processor since the
        // snap phase will convert them to baffles and this only works for
        // internal faces.
        meshRefiner_.balance
        (
            false,
            true,                           // keepZoneFaces
            false,                          // keepBaffles
            weightField,
            decomposer_,
            distributor_
        );

        if (debug)
        {
            meshRefiner_.checkZoneFaces();
        }
    }

    List<labelPair> baffles(0);
    if (controller_.algorithm() == meshControl::STANDARD && zonesToBaffle.size())
    {
        // Create baffles from boundary named surfaces
        {
            labelList originatingFaceZone;
            meshRefiner_.createZoneBaffles
            (
                zonesToBaffle,
                baffles,
                originatingFaceZone
             );
        }
    }

    // Duplicate points on baffles that are on more than one cell
    // region. This will help snapping pull them to separate surfaces.
    meshRefiner_.dupNonManifoldPoints();

    if (controller_.algorithm() == meshControl::STANDARD && zonesToBaffle.size())
    {
        removeBoundaryZoneProblemCells
        (
            refineParams,
            snapParams,
            motionDict,
            baffles
        );
    }

    // Estimate cell sizes
    if (dryRun)
    {
        snappyVoxelMeshDriver voxelDriver
        (
            meshRefiner_,
            globalToMasterPatch_,
            globalToSlavePatch_
        );
        voxelDriver.doRefine(refineParams);
    }

    if
    (
        !refineParams.fullLeakChecks()
        && refineParams.locationsOutsideMesh().size()
    )
    {
        const boolList blockedFace(mesh.nFaces(), false);
        regionSplit cellRegion(mesh, blockedFace);
        meshRefiner_.findRegions
        (
            mesh,
            meshRefiner_.meshCutter(),
            meshRefiner_.mergeDistance()*vector(1,1,1),   // perturbVec
            refineParams.locationsInMesh(),
            refineParams.locationsOutsideMesh(),
            setFormatter_,
            true,
            cellRegion.nRegions(),
            cellRegion,
            blockedFace
        );
    }

    //Clear up as should no longer be required
    meshRefiner_.userFaceData().clear();

} //doRefine


#ifdef FOAM_USE_OPENVDB
void Foam::snappyRefineDriver::doVDBRefine
(
    const dictionary& meshDict,
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool prepareForSnapping /*= true*/
)
{

    Info<< nl
        << "VDB Refinement phase" << nl
        << "--------------------" << nl
        << endl;

    const Tuple2<scalar, label> blockData(meshDict.lookup("blockData"));

    //Calculate length scale and cell level
    scalar length = blockData.first();
    label level = blockData.second();

    //Convert length to level 0
    scalar level0Edge = length * Foam::pow(2, level);

    label maxLevel = meshRefiner_.maxCellLevel();
    scalar voxelSize = level0Edge / Foam::pow(2, maxLevel);

    Info<< "minVoxelSize: " << voxelSize
        << " at cellLevel " << maxLevel
        << "; level0Edge: " << level0Edge
        << nl << endl;

    foamVDB vdb(voxelSize, maxLevel);

    vdb.voxelize
    (
        meshRefiner_,
        refineParams,
        meshDict.subDict("VDBdomain")
    );
} //do VDBRefine


void Foam::snappyRefineDriver::doVDBzones
(
    const dictionary& meshDict,
    const refinementParameters& refineParams,
    const snapParameters& snapParams,
    const bool prepareForSnapping /*= true*/
)
{
    fvMesh& mesh = meshRefiner_.mesh();

    //update surfaceIndex
    meshRefiner_.updateRegionID
    (
        labelList(mesh.nCells(), 0)
    );
    meshRefiner_.surfaceIndex() = labelList(mesh.nFaces(), -1);
    meshRefiner_.wrapIndex() = labelList(mesh.nFaces(), -1);

    // recalculate intersections for all faces
    // updates meshRefiner_.surfaceIndex(), e.g. for every mesh face
    // index of intersecting surface (or -1 if not intersected)
    const Switch threaded(meshDict.lookupOrDefault("threaded", false));
    meshRefiner_.updateIntersections(identity(mesh.nFaces()), threaded);

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    // Introduce baffles at surface intersections. Remove sections unreachable
    // from keepPoint.
    baffleAndSplitMesh
    (
        refineParams,
        snapParams,
        prepareForSnapping,
        motionDict,
        threaded
    );

    // Mesh is at its finest. Do optional zoning (cellZones and faceZones)
    wordPairHashTable zonesToFaceZone;
    zonify(refineParams, zonesToFaceZone);

    // Create pairs of patches for faceZones
    {
        HashTable<Pair<word>> faceZoneToPatches(zonesToFaceZone.size());

        //    Note: zonesToFaceZone contains the same data on different
        //          processors but in different order. We could sort the
        //          contents but instead just loop in sortedToc order.
        List<Pair<word>> czs(zonesToFaceZone.sortedToc());

        forAll(czs, i)
        {
            const Pair<word>& czNames = czs[i];
            const word& fzName = zonesToFaceZone[czNames];

            const word& masterName = fzName;
            const word slaveName = czNames.second() + "_to_" + czNames.first();
            Pair<word> patches(masterName, slaveName);
            faceZoneToPatches.insert(fzName, patches);
        }
        addFaceZones(meshRefiner_, refineParams, faceZoneToPatches);
    }
} //doVDBzones
#endif

// ************************************************************************* //
