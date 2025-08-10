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
    (c) 2015 ICON CFD
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

//include this first because of VDB include conflicts down the line
#include "searchableSurfaces/triSurfaceMesh/triSurfaceMesh.H"
#include "meshRefinement/meshRefinement.H"
#include "trackedParticle/trackedParticle.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "db/Time/Time.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "refinementFeatures/refinementFeatures.H"
#include "shellSurfaces/shellSurfaces.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/pointSet.H"
#include "decompositionMethod/decompositionMethod.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"
#include "Cloud/Cloud.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#include "surfaceFormats/obj/OBJstream.H"
#include "algorithms/indexedOctree/treeDataCell.H"
#include "searchableSurfaces/searchableSphere/searchableSphere.H"
#include "searchableSurfaces/searchableBox/searchableBox.H"
#include "searchableSurfaces/searchableCylinder/searchableCylinder.H"
#include "searchableSurfaces/searchableRing/searchableRing.H"
#include "searchableSurfaces/searchableCone/searchableCone.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    //- To compare normals
    static bool less(const vector& x, const vector& y)
    {
        for (direction i = 0; i < vector::nComponents; i++)
        {
            if (x[i] < y[i])
            {
                return true;
            }
            else if (x[i] > y[i])
            {
                return false;
            }
        }
        // All components the same
        return false;
    }


    //- To compare normals
    class normalLess
    {
        const vectorList& values_;

    public:

        normalLess(const vectorList& values)
        :
            values_(values)
        {}

        bool operator()(const label a, const label b) const
        {
            return less(values_[a], values_[b]);
        }
    };


    //- Template specialization for pTraits<labelList> so we can have fields
    template<>
    class pTraits<labelList>
    {

    public:

        //- Component type
        typedef labelList cmptType;
    };

    //- Template specialization for pTraits<labelList> so we can have fields
    template<>
    class pTraits<vectorList>
    {

    public:

        //- Component type
        typedef vectorList cmptType;
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Get faces (on the new mesh) that have in some way been affected by the
// mesh change. Picks up all faces but those that are between two
// unrefined faces. (Note that of an unchanged face the edge still might be
// split but that does not change any face centre or cell centre.
Foam::labelList Foam::meshRefinement::getChangedFaces
(
    const mapPolyMesh& map,
    const labelList& oldCellsToRefine
)
{
    const polyMesh& mesh = map.mesh();

    labelList changedFaces;

    // For reporting: number of masterFaces changed
    label nMasterChanged = 0;
    PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh));

    {
        // Mark any face on a cell which has been added or changed
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // Note that refining a face changes the face centre (for a warped face)
        // which changes the cell centre. This again changes the cellcentre-
        // cellcentre edge across all faces of the cell.
        // Note: this does not happen for unwarped faces but unfortunately
        // we don't have this information.

        const labelList& faceOwner = mesh.faceOwner();
        const labelList& faceNeighbour = mesh.faceNeighbour();
        const cellList& cells = mesh.cells();
        const label nInternalFaces = mesh.nInternalFaces();

        // Mark refined cells on old mesh
        PackedBoolList oldRefineCell(map.nOldCells());

        forAll(oldCellsToRefine, i)
        {
            oldRefineCell.set(oldCellsToRefine[i], 1u);
        }

        // Mark refined faces
        PackedBoolList refinedInternalFace(nInternalFaces);

        // 1. Internal faces

        for (label faceI = 0; faceI < nInternalFaces; faceI++)
        {
            label oldOwn = map.cellMap()[faceOwner[faceI]];
            label oldNei = map.cellMap()[faceNeighbour[faceI]];

            if
            (
                oldOwn >= 0
             && oldRefineCell.get(oldOwn) == 0u
             && oldNei >= 0
             && oldRefineCell.get(oldNei) == 0u
            )
            {
                // Unaffected face since both neighbours were not refined.
            }
            else
            {
                refinedInternalFace.set(faceI, 1u);
            }
        }


        // 2. Boundary faces

        boolList refinedBoundaryFace(mesh.nFaces()-nInternalFaces, false);

        forAll(mesh.boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh.boundaryMesh()[patchI];

            label faceI = pp.start();

            forAll(pp, i)
            {
                label oldOwn = map.cellMap()[faceOwner[faceI]];

                if (oldOwn >= 0 && oldRefineCell.get(oldOwn) == 0u)
                {
                    // owner did exist and wasn't refined.
                }
                else
                {
                    refinedBoundaryFace[faceI-nInternalFaces] = true;
                }
                faceI++;
            }
        }

        // Synchronise refined face status
        syncTools::syncBoundaryFaceList
        (
            mesh,
            refinedBoundaryFace,
            orEqOp<bool>()
        );


        // 3. Mark all faces affected by refinement. Refined faces are in
        //    - refinedInternalFace
        //    - refinedBoundaryFace
        boolList changedFace(mesh.nFaces(), false);

        forAll(refinedInternalFace, faceI)
        {
            if (refinedInternalFace.get(faceI) == 1u)
            {
                const cell& ownFaces = cells[faceOwner[faceI]];
                forAll(ownFaces, ownI)
                {
                    changedFace[ownFaces[ownI]] = true;
                }
                const cell& neiFaces = cells[faceNeighbour[faceI]];
                forAll(neiFaces, neiI)
                {
                    changedFace[neiFaces[neiI]] = true;
                }
            }
        }

        forAll(refinedBoundaryFace, i)
        {
            if (refinedBoundaryFace[i])
            {
                const cell& ownFaces = cells[faceOwner[i+nInternalFaces]];
                forAll(ownFaces, ownI)
                {
                    changedFace[ownFaces[ownI]] = true;
                }
            }
        }

        syncTools::syncFaceList
        (
            mesh,
            changedFace,
            orEqOp<bool>()
        );


        // Now we have in changedFace marked all affected faces. Pack.
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        changedFaces = findIndices(changedFace, true);

        // Count changed master faces.
        nMasterChanged = 0;

        forAll(changedFace, faceI)
        {
            if (changedFace[faceI] && isMasterFace[faceI])
            {
                nMasterChanged++;
            }
        }

    }

    if (debug&meshRefinement::MESH)
    {
        Pout<< "getChangedFaces : Detected "
            << " local:" << changedFaces.size()
            << " global:" << returnReduce(nMasterChanged, sumOp<label>())
            << " changed faces out of " << mesh.globalData().nTotalFaces()
            << endl;

        faceSet changedFacesSet
        (
            mesh,
            "changedFaces",
            labelHashSet(changedFaces)
        );
        Pout<< "getChangedFaces : Writing " << changedFaces.size()
            << " changed faces to faceSet " << changedFacesSet.name()
            << endl;
        changedFacesSet.write();
    }

    return changedFaces;
}


// Mark cell for refinement (if not already marked). Return false if
// refinelimit hit. Keeps running count (in nRefine) of cells marked for
// refinement
bool Foam::meshRefinement::markForRefine
(
    const label markValue,
    const label nAllowRefine,

    label& cellValue,
    label& nRefine
)
{
    if (cellValue == -1)
    {
        cellValue = markValue;
        nRefine++;
    }

    return nRefine <= nAllowRefine;
}


void Foam::meshRefinement::markFeatureCellLevel
(
    const refinementParameters& refineParams,
    labelList& maxFeatureLevel
) const
{
    // We want to refine all cells containing a feature edge.
    // - don't want to search over whole mesh
    // - don't want to build octree for whole mesh
    // - so use tracking to follow the feature edges
    //
    // 1. find non-manifold points on feature edges (i.e. start of feature edge
    //    or 'knot')
    // 2. seed particle starting at keepPoint going to this non-manifold point
    // 3. track particles to their non-manifold point
    //
    // 4. track particles across their connected feature edges, marking all
    //    visited cells with their level (through trackingData)
    // 5. do 4 until all edges have been visited.

    const point& keepPoint = refineParams.locationsInMesh()[0];
    const scalar& minFeatureLength = refineParams.minFeatureLength();

    // Find all start cells of features. Is done by tracking from keepPoint.
    Cloud<trackedParticle> cloud
    (
        mesh_,
        "startPointCloud",
        IDLList<trackedParticle>()
    );

    label cellI = -1;
    label tetFaceI = -1;
    label tetPtI = -1;

    refinementFeatures& features = featuresPtr_();
    manifoldFeatures& mFeatures = features.manFeatures();
    const PtrList<Tuple2<labelList, label>>& featureMeshes =
        mFeatures.features(100,minFeatureLength);

    labelList& seedCells = mFeatures.seeds();
    label maxSeedCell = gMax(seedCells);

    vector startPt = vector::zero;
    if (maxSeedCell == -1)
    {
        //Find average feature start point
        forAll(featureMeshes, featI)
        {
            const label origFeatI =
                featureMeshes[featI].second();
            const label fPointI =
                featureMeshes[featI].first()[0];
            point basePoint =
                features[origFeatI].points()[fPointI];

            startPt += basePoint;
        }

        label nFtrs = featureMeshes.size();
        reduce(
            std::tie(nFtrs, startPt),
            ParallelOp<sumOp<label>, sumOp<point>>{}
        );

        if (nFtrs > 0)
        {
            startPt = startPt / nFtrs;
        }

        // Force construction of search tree even if processor holds no
        // cells
        (void)mesh_.cellTree();
        if (mesh_.nCells())
        {
            mesh_.findCellFacePt(startPt, cellI, tetFaceI, tetPtI, false);
        }

        if (returnReduce(cellI, maxOp<label>()) == -1)
        {
            startPt = keepPoint;
            if (mesh_.nCells())
            {
                mesh_.findCellFacePt(keepPoint, cellI, tetFaceI, tetPtI);
            }
        }
    }

    //Basic reporting on number of feature lines
    Info<<"Number of feature lines: "<<featureMeshes.size()<<endl;
    label nFtrPts = 0;
    forAll(featureMeshes, featI)
    {
        nFtrPts += featureMeshes[featI].first().size();
    }
    Info<<"Number of feature points: "<<nFtrPts<<endl;

    forAll(featureMeshes, featI)
    {
        label seedCellI = -1;
        if (maxSeedCell == -1)
        {
            seedCellI = cellI;
        }
        else
        {
            seedCellI = seedCells[featI];
            if (seedCellI != -1)
            {
                startPt = mesh_.cellCentres()[seedCellI];
                mesh_.findCellFacePt(startPt, seedCellI, tetFaceI, tetPtI, false);
            }
        }

        if (seedCellI != -1)
        {
            const label origFeatI =
                featureMeshes[featI].second();
            const label fPointI =
                featureMeshes[featI].first()[0];
            point basePoint =
                features[origFeatI].points()[fPointI];
            label baseLevel =
                features.levels()[origFeatI][0];

            //Create particle.
            cloud.addParticle
            (
                new trackedParticle
                (
                    mesh_,
                    startPt,
                    seedCellI,
                    tetFaceI,
                    tetPtI,
                    basePoint, // endpos
                    baseLevel, // level
                    featI,     // featureMesh
                    0          // stating feature point
                    -1         // feature edge
                )
            );
        }
    }

    seedCells = -1;

    // Database to pass into trackedParticle::move
    trackedParticle::trackingData td(cloud, maxFeatureLevel);

    // Track all particles to their end position (= starting feature point)
    cloud.move(td, GREAT);

    boolList trackingFtr(featureMeshes.size());
    boolList trackingFtrPrev(featureMeshes.size(),true);

    while (true)
    {
        trackingFtr = false;
        label nParticles = 0;
        // Make particle follow edge.
        forAllIter(Cloud<trackedParticle>, cloud, iter)
        {
            trackedParticle& tp = iter();

            label pointI = tp.j();
            label featI = tp.i();

            if (tp.onFtr() && !trackingFtrPrev[featI])
            {
                nParticles++;
            }
            else
            {
                if
                (
                    pointI == featureMeshes[featI].first().size()-1
                )
                {
                    cloud.deleteParticle(tp);
                    continue;
                }
                else
                {
                    // Keep particle
                    nParticles++;
                }
            }

            const label origFeatI = featureMeshes[featI].second();

            if (tp.onFtr() && tp.j() > 0 && !trackingFtrPrev[featI])
            {
                label nextPt = featureMeshes[featI].first()[tp.j()-1];
                tp.end() = features[origFeatI].points()[nextPt];
            }
            else
            {
                tp.j()++;
                label nextPt = featureMeshes[featI].first()[tp.j()];
                tp.end() = features[origFeatI].points()[nextPt];
            }

            if (tp.onFtr())
            {
                if (tp.j() == 0 || (tp.j() > 0 && !trackingFtrPrev[featI]))
                {
                    seedCells[featI] = tp.cell();
                }
                trackingFtr[featI]= true;
            }
        }

        reduce(nParticles, sumOp<label>());
        if (nParticles == 0)
        {
            break;
        }

        Pstream::listCombineReduce(trackingFtr, orOp<bool>());

        trackingFtrPrev = trackingFtr;

        // Track all particles to their end position.
        cloud.move(td, GREAT);
    }
}


// Calculates list of cells to refine based on surface bounding box size
// with respect to the local cell size.
Foam::label Foam::meshRefinement::markInitialRefinement
(
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const scalar edge0Len = meshCutter_.level0EdgeLength();
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    const searchableSurfaces& geometry = surfaces_.geometry();
    const labelList& surfaces = surfaces_.surfaces();
    const labelList& shells = shells_.shells();

    scalarField minLevel(geometry.size(), -1);
    scalarField bbMag(geometry.size(), -1);
    pointField bbCentre(geometry.size(), vector::zero);

    scalar scaledRefLength = 2*edge0Len;

    forAll(surfaces, surfI)
    {
        label geomI = surfaces[surfI];

        const searchableSurface& s = geometry[geomI];
        if
        (
            isA<triSurfaceMesh>(s) || isA<searchableSphere>(s)
            || isA<searchableBox>(s) || isA<searchableCylinder>(s)
            || isA<searchableRing>(s) || isA<searchableCone>(s)
        )
        {
            const boundBox& bb = s.bounds();
            bbCentre[geomI] = bb.midpoint();
            bbMag[geomI] = bb.mag();

            if (bbMag[geomI] < scaledRefLength)
            {
                forAll(s.regions(), regionI)
                {
                    label mLevel = surfaces_.minLevel(surfI,regionI);
                    minLevel[geomI] = max(minLevel[geomI], mLevel);
                }
            }
        }
    }

    forAll(shells, surfI)
    {
        label geomI = shells[surfI];

        const searchableSurface& s = geometry[geomI];
        if
        (
            isA<triSurfaceMesh>(s) || isA<searchableSphere>(s)
            || isA<searchableBox>(s) || isA<searchableCylinder>(s)
            || isA<searchableRing>(s) || isA<searchableCone>(s)
        )
        {
            const boundBox& bb = s.bounds();
            bbCentre[geomI] = bb.midpoint();
            bbMag[geomI] = bb.mag();

            if (bbMag[geomI] < scaledRefLength)
            {
                const List<labelVector>& dLevels = shells_.levels()[surfI];
                forAll(dLevels, dLI)
                {
                    minLevel[geomI] =
                        max(minLevel[geomI], cmptMax(dLevels[dLI]));
                }
            }
        }
    }

    if (max(minLevel) == -1)
    {
        // No surface smaller than a factor of largest edge length found
        return 0;
    }

    label oldNRefine = nRefine;

    forAll(mesh_.cells(), cellI)
    {
        point cc = cellCentres[cellI];
        label cLevel = cellLevel[cellI];

        scalar spacing = scaledRefLength / pow(2., cLevel);
        label markedMinLevel = -1;

        forAll(geometry, geomI)
        {
            if (bbMag[geomI] > SMALL)
            {
                if (cLevel >= minLevel[geomI])
                {
                    continue;
                }

                scalar ccToBoxCC = mag(cc - bbCentre[geomI]);
                if (ccToBoxCC < spacing && bbMag[geomI] < spacing)
                {
                    markedMinLevel = minLevel[geomI];
                    break;
                }
            }
        }

        if (markedMinLevel > -1)
        {
            bool reachedLimit = !markForRefine
            (
                markedMinLevel,    // mark with any positive value
                nAllowRefine,
                refineCell[cellI],
                nRefine
            );

            if (reachedLimit)
            {
                if (debug)
                {
                    Pout<< "Stopped refining seed cells"
                        << " since reaching my cell limit of "
                        << mesh_.nCells()+7*nRefine << endl;
                }
                break;
            }
        }
    }

    return returnReduce(nRefine-oldNRefine,  sumOp<label>());
}


// Calculates list of cells to refine based on intersection with feature edge.
Foam::label Foam::meshRefinement::markFeatureRefinement
(
    const refinementParameters& refineParams,
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    // Largest refinement level of any feature passed through
    labelList maxFeatureLevel(mesh_.nCells(), -1);
    markFeatureCellLevel(refineParams, maxFeatureLevel);

    // See which cells to refine. maxFeatureLevel will hold highest level
    // of any feature edge that passed through.

    const labelList& cellLevel = meshCutter_.cellLevel();

    label oldNRefine = nRefine;

    label minWrapLevel = min(surfaces_.wrapLevel());

    forAll(maxFeatureLevel, cellI)
    {
        if (maxFeatureLevel[cellI] > cellLevel[cellI])
        {
            if
            (
                minWrapLevel > -1
                &&!wrapActive_
                && cellLevel[cellI] >= minWrapLevel
            )
            {
                continue;
            }

            // Mark
            if
            (
               !markForRefine
                (
                    0,                      // surface (n/a)
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                )
            )
            {
                // Reached limit
                break;
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine,  sumOp<label>());
}


// Mark cells for distance-to-feature based refinement.
Foam::label Foam::meshRefinement::markInternalDistanceToFeatureRefinement
(
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    // Detect if there are any distance shells
    if (featuresPtr_().maxDistance() <= 0.0)
    {
        return 0;
    }

    label oldNRefine = nRefine;

    // Collect cells to test
    pointField testCc(cellLevel.size()-nRefine);
    labelList testLevels(cellLevel.size()-nRefine);
    label testI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            testCc[testI] = cellCentres[cellI];
            testLevels[testI] = cellLevel[cellI];
            testI++;
        }
    }

    // Do test to see whether cells are inside/outside shell with higher level
    labelList maxLevel;
    featuresPtr_().findHigherLevel(testCc, testLevels, maxLevel);

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.
    testI = 0;
    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            if (maxLevel[testI] > testLevels[testI])
            {
                bool reachedLimit = !markForRefine
                (
                    maxLevel[testI],    // mark with any positive value
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );

                if (reachedLimit)
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining internal cells"
                            << " since reaching my cell limit of "
                            << mesh_.nCells()+7*nRefine << endl;
                    }
                    break;
                }
            }
            testI++;
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Mark cells for non-surface intersection based refinement.
Foam::label Foam::meshRefinement::markInternalRefinement
(
    const label nAllowRefine,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // Collect cells to test
    pointField testCc(cellLevel.size()-nRefine);
    labelList testLevels(cellLevel.size()-nRefine);
    labelList testMap(cellLevel.size()-nRefine);
    label testI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            testCc[testI] = cellCentres[cellI];
            testLevels[testI] = cellLevel[cellI];
            testMap[testI] = cellI;
            testI++;
        }
    }

    // Do test to see whether cells are inside/outside shell with higher level
    labelList maxLevel;
    labelList minCmptLevel;
    labelList shellCheck;
    shells_.findHigherLevel
    (
        mesh_,
        testCc,
        testLevels,
        testMap,
        maxLevel,
        minCmptLevel,
        shellCheck
    );

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.
    testI = 0;
    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] == -1)
        {
            if (maxLevel[testI] > testLevels[testI])
            {
                bool reachedLimit = !markForRefine
                (
                    maxLevel[testI],    // mark with any positive value
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );

                if (reachedLimit)
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining internal cells"
                            << " since reaching my cell limit of "
                            << mesh_.nCells()+7*nRefine << endl;
                    }
                    break;
                }
            }
            testI++;
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


Foam::label Foam::meshRefinement::unmarkInternalRefinement
(
    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // Collect cells to test
    pointField testCc(nRefine);
    labelList testLevels(nRefine);
    label testI = 0;

    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] >= 0)
        {
            testCc[testI] = cellCentres[cellI];
            testLevels[testI] = cellLevel[cellI];
            testI++;
        }
    }

    // Do test to see whether cells are inside/outside shell with higher level
    labelList levelShell;
    limitShells_.findLevel(testCc, testLevels, levelShell);

    // Mark for refinement. Note that we didn't store the original cellID so
    // now just reloop in same order.
    testI = 0;
    forAll(cellLevel, cellI)
    {
        if (refineCell[cellI] >= 0)
        {
            if (levelShell[testI] != -1)
            {
                refineCell[cellI] = -1;
                nRefine--;
            }
            testI++;
        }
    }

    return returnReduce(oldNRefine-nRefine, sumOp<label>());
}


// Collect faces that are intersected and whose neighbours aren't yet marked
// for refinement.
Foam::labelList Foam::meshRefinement::getRefineCandidateFaces
(
    const labelList& refineCell
) const
{
    labelList testFaces(mesh_.nFaces());

    label nTest = 0;

    forAll(surfaceIndex_, faceI)
    {
        if (surfaceIndex_[faceI] != -1)
        {
            label own = mesh_.faceOwner()[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                if (refineCell[own] == -1 || refineCell[nei] == -1)
                {
                    testFaces[nTest++] = faceI;
                }
            }
            else
            {
                if (refineCell[own] == -1)
                {
                    testFaces[nTest++] = faceI;
                }
            }
        }
    }
    testFaces.setSize(nTest);

    return testFaces;
}


// Mark cells for surface intersection based refinement.
Foam::label Foam::meshRefinement::markSurfaceRefinement
(
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    label oldNRefine = nRefine;

    // Use cached surfaceIndex_ to detect if any intersection. If so
    // re-intersect to determine level wanted.

    // Collect candidate faces
    // ~~~~~~~~~~~~~~~~~~~~~~~

    labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    // ~~~~~~~~~~~~~~~~

    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    calcCellCellRays
    (
        mesh_.cellCentres(),
        neiCc,
        neiLevel,
        testFaces,
        start,
        end,
        minLevel
    );


    // Do test for higher intersections
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList surfaceHit;
    labelList surfaceMinLevel;
    surfaces_.findHigherIntersection
    (
        shells_,
        start,
        end,
        minLevel,

        surfaceHit,
        surfaceMinLevel,
        meshCutter_.level0EdgeLength(),
        true
    );
    // Mark cells for refinement
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    label minWrapLevel = min(surfaces_.wrapLevel());

    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        label surfI = surfaceHit[i];

        if (surfI != -1)
        {
            // Found intersection with surface with higher wanted
            // refinement. Check if the level field on that surface
            // specifies an even higher level. Note:this is weird. Should
            // do the check with the surfaceMinLevel whilst intersecting the
            // surfaces?

            label own = mesh_.faceOwner()[faceI];

            if (surfaceMinLevel[i] > cellLevel[own])
            {
                if
                (
                    minWrapLevel > -1
                    && !wrapActive_
                    && cellLevel[own] >= minWrapLevel
                )
                {
                    continue;
                }

                // Owner needs refining
                if
                (
                   !markForRefine
                    (
                        surfI,
                        nAllowRefine,
                        refineCell[own],
                        nRefine
                    )
                )
                {
                    break;
                }
            }

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];
                if (surfaceMinLevel[i] > cellLevel[nei])
                {
                    if
                    (
                        minWrapLevel > -1
                        && !wrapActive_
                        && cellLevel[nei] >= minWrapLevel
                    )
                    {
                        continue;
                    }

                    // Neighbour needs refining
                    if
                    (
                       !markForRefine
                        (
                            surfI,
                            nAllowRefine,
                            refineCell[nei],
                            nRefine
                        )
                    )
                    {
                        break;
                    }
                }
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Count number of matches of first argument in second argument
Foam::label Foam::meshRefinement::countMatches
(
    const List<point>& normals1,
    const List<point>& normals2,
    const scalar tol
) const
{
    label nMatches = 0;

    forAll(normals1, i)
    {
        const vector& n1 = normals1[i];

        forAll(normals2, j)
        {
            const vector& n2 = normals2[j];

            if (magSqr(n1-n2) < tol)
            {
                nMatches++;
                break;
            }
        }
    }

    return nMatches;
}


Foam::labelList Foam::meshRefinement::markRemovalInterfaceRefine
(
    const refinementParameters& refineParams
) const
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const labelList& cellLevel = meshCutter().cellLevel();
    const labelList& pointLevel = meshCutter().pointLevel();


    labelList cellToZone(mesh_.nCells(),-2);
    labelList namedSurfaceIndex;  //filtered named intersections
    labelList namedIntersections; //all named intersections
    PackedBoolList posOrientation;

    labelList testFaces;
    labelList globalRegion1;
    labelList globalRegion2;

    zonify
    (
        -2,                 // zone to put unreached cells into
        refineParams,
        0,

        testFaces,
        globalRegion1,
        globalRegion2,

        cellToZone,
        namedSurfaceIndex,
        namedIntersections,
        posOrientation
    );

    boolList isBoundaryPoint(mesh_.nPoints(), false);
    forAll(globalRegion1, facei)
    {
        if (globalRegion1[facei] != -1 || globalRegion2[facei] != -1)
        {
            const face& f = mesh_.faces()[facei];
            forAll(f,fp)
            {
                isBoundaryPoint[f[fp]] = true;
            }
        }
    }

    labelList cellBoundaryAnchors(mesh_.nCells(), 0);
    forAll(mesh_.cells(), celli)
    {
        label cZone = cellToZone[celli];

        if (cZone == -2)
        {
            continue;
        }
        const labelList& cPts = mesh_.cellPoints()[celli];
        label nBoundaryAnchors = 0;
        forAll(cPts, i)
        {
            label pointi = cPts[i];
            if (pointLevel[pointi] <= cellLevel[celli])
            {
                // Anchor point
                if (isBoundaryPoint[pointi])
                {
                    nBoundaryAnchors++;
                }
            }
        }

        cellBoundaryAnchors[celli] = nBoundaryAnchors;
    }

    labelList nPointGapBdyCell(mesh_.nPoints(), 0);
    labelList nPointInteriorCells(mesh_.nPoints(), 0);
    forAll(nPointGapBdyCell, pointi)
    {
        const labelList& pCells = mesh_.pointCells()[pointi];
        forAll(pCells, pCI)
        {
            label celli = pCells[pCI];

            if (cellBoundaryAnchors[celli] > 5)
            {
                nPointGapBdyCell[pointi]++;
            }
            if (cellToZone[celli] != -2)
            {
                nPointInteriorCells[pointi]++;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        nPointGapBdyCell,
        plusEqOp<label>(),
        label(0)              // null value
    );

    syncTools::syncPointList
    (
        mesh_,
        nPointInteriorCells,
        plusEqOp<label>(),
        label(0)              // null value
    );

    boolList gapCells(mesh_.nCells(), false);
    forAll(mesh_.cells(), celli)
    {
        if (cellBoundaryAnchors[celli] == 8)
        {
            gapCells[celli] = true;
        }
        else
        {
            const labelList& cPts = mesh_.cellPoints()[celli];
            bool mark = false;
            forAll(cPts, cPtI)
            {
                label pointi = cPts[cPtI];
                if
                (
                    nPointGapBdyCell[pointi] > 3
                    && nPointGapBdyCell[pointi] == nPointInteriorCells[pointi]
                )
                {
                    mark = true;
                    break;
                 }
            }
            if (mark)
            {
                gapCells[celli] = true;
            }
        }
    }

    boolList neiGapCells;
    syncTools::swapBoundaryCellList(mesh_, gapCells, neiGapCells);

    labelList neiCellZone;
    syncTools::swapBoundaryCellList(mesh_, cellToZone, neiCellZone);
    labelList neiCellLevel;
    syncTools::swapBoundaryCellList(mesh_, cellLevel, neiCellLevel);

    DynamicList<label> refineCells(mesh_.nCells()/10);

    forAll(cellToZone, celli)
    {
        label cZone = cellToZone[celli];

        if (cZone == -2 || (cZone !=-2 && gapCells[celli]))
        {
            const cell& c = mesh_.cells()[celli];
            label cLevel = cellLevel[celli];

            forAll(c, cfI)
            {
                label facei = c[cfI];
                label patchi = patches.whichPatch(facei);
                if (patchi == -1 || patches[patchi].coupled())
                {
                    label neiZone = -1;
                    label neiLevel = -1;
                    bool neiGap  = false;

                    if (patchi == -1)
                    {
                        label nei = mesh_.faceNeighbour()[facei];
                        if (celli == nei)
                        {
                            label own = mesh_.faceOwner()[facei];
                            neiZone = cellToZone[own];
                            neiLevel = cellLevel[own];
                            neiGap = gapCells[own];
                        }
                        else
                        {
                            neiZone  = cellToZone[nei];
                            neiLevel = cellLevel[nei];
                            neiGap = gapCells[nei];
                        }
                    }
                    else
                    {
                        neiZone = neiCellZone[facei-mesh_.nInternalFaces()];
                        neiLevel = neiCellLevel[facei-mesh_.nInternalFaces()];
                        neiGap = neiGapCells[facei-mesh_.nInternalFaces()];
                    }

                    if
                    (
                        neiLevel > cLevel &&
                        (
                            (cZone == -2 && neiZone != -2 && neiGap)
                            || (cZone != -2 && neiZone == -2)
                        )
                    )
                    {
                        refineCells.append(celli);
                        break;
                    }
                }
            }
        }
    }

    return labelList(refineCells, true);
}


// Mark cells for surface curvature based refinement.
Foam::label Foam::meshRefinement::markSurfaceCurvatureRefinement
(
    const scalar curvature,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const pointField& cellCentres = mesh_.cellCentres();

    label oldNRefine = nRefine;

    // 1. local test: any cell on more than one surface gets refined
    // (if its current level is < max of the surface max level)

    // 2. 'global' test: any cell on only one surface with a neighbour
    // on a different surface gets refined (if its current level etc.)


    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));


    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined.
    labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    // Note: uses isMasterFace otherwise could be call to calcCellCellRays
    forAll(testFaces, i)
    {
        label faceI = testFaces[i];

        label own = mesh_.faceOwner()[faceI];

        if (mesh_.isInternalFace(faceI))
        {
            label nei = mesh_.faceNeighbour()[faceI];

            start[i] = cellCentres[own];
            end[i] = cellCentres[nei];
            minLevel[i] = min(cellLevel[own], cellLevel[nei]);
        }
        else
        {
            label bFaceI = faceI - mesh_.nInternalFaces();

            start[i] = cellCentres[own];
            end[i] = neiCc[bFaceI];

            if (!isMasterFace[faceI])
            {
                Swap(start[i], end[i]);
            }

            minLevel[i] = min(cellLevel[own], neiLevel[bFaceI]);
        }
    }

    // Extend segments a bit
    {
        const vectorField smallVec(ROOTSMALL*(end-start));
        start -= smallVec;
        end += smallVec;
    }


    // Test for all intersections (with surfaces of higher max level than
    // minLevel) and cache per cell the interesting inter
    labelListList cellSurfLevels(mesh_.nCells());
    List<vectorList> cellSurfNormals(mesh_.nCells());

    {
        // Per segment the normals of the surfaces
        List<vectorList> surfaceNormal;
        // Per segment the list of levels of the surfaces
        labelListList surfaceLevel;

        surfaces_.findAllHigherIntersections
        (
            start,
            end,
            minLevel,           // max level of surface has to be bigger
                                // than min level of neighbouring cells

            surfaces_.maxLevel(),

            surfaceNormal,
            surfaceLevel
        );


        // Sort the data according to surface location. This will guarantee
        // that on coupled faces both sides visit the intersections in
        // the same order so will decide the same
        labelList visitOrder;
        forAll(surfaceNormal, pointI)
        {
            vectorList& pNormals = surfaceNormal[pointI];
            labelList& pLevel = surfaceLevel[pointI];

            sortedOrder(pNormals, visitOrder, normalLess(pNormals));

            pNormals = List<point>(pNormals, visitOrder);
            pLevel = UIndirectList<label>(pLevel, visitOrder);
        }

        // Clear out unnecessary data
        start.clear();
        end.clear();
        minLevel.clear();

        // Convert face-wise data to cell.
        forAll(testFaces, i)
        {
            label faceI = testFaces[i];
            label own = mesh_.faceOwner()[faceI];

            const vectorList& fNormals = surfaceNormal[i];
            const labelList& fLevels = surfaceLevel[i];

            forAll(fNormals, hitI)
            {
                if (fLevels[hitI] > cellLevel[own])
                {
                    cellSurfLevels[own].append(fLevels[hitI]);
                    cellSurfNormals[own].append(fNormals[hitI]);
                }

                if (mesh_.isInternalFace(faceI))
                {
                    label nei = mesh_.faceNeighbour()[faceI];
                    if (fLevels[hitI] > cellLevel[nei])
                    {
                        cellSurfLevels[nei].append(fLevels[hitI]);
                        cellSurfNormals[nei].append(fNormals[hitI]);
                    }
                }
            }
        }
    }



    // Bit of statistics
    if (debug)
    {
        label nSet = 0;
        label nNormals = 9;
        forAll(cellSurfNormals, cellI)
        {
            const vectorList& normals = cellSurfNormals[cellI];
            if (normals.size())
            {
                nSet++;
                nNormals += normals.size();
            }
        }
        reduce(
            std::tie(nSet, nNormals),
            UniformParallelOp<sumOp<label>, 2>{}
        );
        Info<< "markSurfaceCurvatureRefinement :"
            << " cells:" << mesh_.globalData().nTotalCells()
            << " of which with normals:" << nSet
            << " ; total normals stored:" << nNormals
            << endl;
    }



    bool reachedLimit = false;


    // 1. Check normals per cell
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    for
    (
        label cellI = 0;
        !reachedLimit && cellI < cellSurfNormals.size();
        cellI++
    )
    {
        const vectorList& normals = cellSurfNormals[cellI];
        const labelList& levels = cellSurfLevels[cellI];

        // n^2 comparison of all normals in a cell
        for (label i = 0; !reachedLimit && i < normals.size(); i++)
        {
            for (label j = i+1; !reachedLimit && j < normals.size(); j++)
            {
                if ((normals[i] & normals[j]) < curvature)
                {
                    label maxLevel = max(levels[i], levels[j]);

                    if (cellLevel[cellI] < maxLevel)
                    {
                        if
                        (
                            !markForRefine
                            (
                                maxLevel,
                                nAllowRefine,
                                refineCell[cellI],
                                nRefine
                            )
                        )
                        {
                            if (debug)
                            {
                                Pout<< "Stopped refining since reaching my cell"
                                    << " limit of " << mesh_.nCells()+7*nRefine
                                    << endl;
                            }
                            reachedLimit = true;
                            break;
                        }
                    }
                }
            }
        }
    }



    // 2. Find out a measure of surface curvature
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Look at normals between neighbouring surfaces
    // Loop over all faces. Could only be checkFaces, except if they're coupled

    // Internal faces
    for
    (
        label faceI = 0;
        !reachedLimit && faceI < mesh_.nInternalFaces();
        faceI++
    )
    {
        label own = mesh_.faceOwner()[faceI];
        label nei = mesh_.faceNeighbour()[faceI];

        const vectorList& ownNormals = cellSurfNormals[own];
        const labelList& ownLevels = cellSurfLevels[own];
        const vectorList& neiNormals = cellSurfNormals[nei];
        const labelList& neiLevels = cellSurfLevels[nei];


        // Special case: owner normals set is a subset of the neighbour
        // normals. Do not do curvature refinement since other cell's normals
        // do not add any information. Avoids weird corner extensions of
        // refinement regions.
        bool ownIsSubset =
            countMatches(ownNormals, neiNormals)
         == ownNormals.size();

        bool neiIsSubset =
            countMatches(neiNormals, ownNormals)
         == neiNormals.size();


        if (!ownIsSubset && !neiIsSubset)
        {
            // n^2 comparison of between ownNormals and neiNormals
            for (label i = 0; !reachedLimit &&  i < ownNormals.size(); i++)
            {
                for (label j = 0; !reachedLimit && j < neiNormals.size(); j++)
                {
                    // Have valid data on both sides. Check curvature.
                    if ((ownNormals[i] & neiNormals[j]) < curvature)
                    {
                        // See which side to refine.
                        if (cellLevel[own] < ownLevels[i])
                        {
                            if
                            (
                                !markForRefine
                                (
                                    ownLevels[i],
                                    nAllowRefine,
                                    refineCell[own],
                                    nRefine
                                )
                            )
                            {
                                if (debug)
                                {
                                    Pout<< "Stopped refining since reaching"
                                        << " my cell limit of "
                                        << mesh_.nCells()+7*nRefine << endl;
                                }
                                reachedLimit = true;
                                break;
                            }
                        }
                        if (cellLevel[nei] < neiLevels[j])
                        {
                            if
                            (
                                !markForRefine
                                (
                                    neiLevels[j],
                                    nAllowRefine,
                                    refineCell[nei],
                                    nRefine
                                )
                            )
                            {
                                if (debug)
                                {
                                    Pout<< "Stopped refining since reaching"
                                        << " my cell limit of "
                                        << mesh_.nCells()+7*nRefine << endl;
                                }
                                reachedLimit = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }


    // Send over surface normal to neighbour cell.
    List<vectorList> neiSurfaceNormals;
    syncTools::swapBoundaryCellList(mesh_, cellSurfNormals, neiSurfaceNormals);

    // Boundary faces
    for
    (
        label faceI = mesh_.nInternalFaces();
        !reachedLimit && faceI < mesh_.nFaces();
        faceI++
    )
    {
        label own = mesh_.faceOwner()[faceI];
        label bFaceI = faceI - mesh_.nInternalFaces();

        const vectorList& ownNormals = cellSurfNormals[own];
        const labelList& ownLevels = cellSurfLevels[own];
        const vectorList& neiNormals = neiSurfaceNormals[bFaceI];

        // Special case: owner normals set is a subset of the neighbour
        // normals. Do not do curvature refinement since other cell's normals
        // do not add any information. Avoids weird corner extensions of
        // refinement regions.
        bool ownIsSubset =
            countMatches(ownNormals, neiNormals)
         == ownNormals.size();

        bool neiIsSubset =
            countMatches(neiNormals, ownNormals)
         == neiNormals.size();


        if (!ownIsSubset && !neiIsSubset)
        {
            // n^2 comparison of between ownNormals and neiNormals
            for (label i = 0; !reachedLimit &&  i < ownNormals.size(); i++)
            {
                for (label j = 0; !reachedLimit && j < neiNormals.size(); j++)
                {
                    // Have valid data on both sides. Check curvature.
                    if ((ownNormals[i] & neiNormals[j]) < curvature)
                    {
                        if (cellLevel[own] < ownLevels[i])
                        {
                            if
                            (
                                !markForRefine
                                (
                                    ownLevels[i],
                                    nAllowRefine,
                                    refineCell[own],
                                    nRefine
                                )
                            )
                            {
                                if (debug)
                                {
                                    Pout<< "Stopped refining since reaching"
                                        << " my cell limit of "
                                        << mesh_.nCells()+7*nRefine
                                        << endl;
                                }
                                reachedLimit = true;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }


    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


bool Foam::meshRefinement::isGap
(
    const scalar planarCos,
    const vector& point0,
    const vector& normal0,

    const vector& point1,
    const vector& normal1
) const
{
    //- Hits differ and angles are oppositeish and
    //  hits have a normal distance
    vector d = point1-point0;
    scalar magD = mag(d);

    if (magD > mergeDistance())
    {
        scalar cosAngle = (normal0 & normal1);

        vector avg = Zero;
        if (cosAngle < (-1+planarCos))
        {
            // Opposite normals
            avg = 0.5*(normal0-normal1);
        }
        else if (cosAngle > (1-planarCos))
        {
            avg = 0.5*(normal0+normal1);
        }

        if (avg != vector::zero)
        {
            avg /= mag(avg);

            // Check normal distance of intersection locations
            if (mag(avg&d) > mergeDistance())
            {
                return true;
            }
            else
            {
                return  false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}


// Mark small gaps
bool Foam::meshRefinement::isNormalGap
(
    const scalar planarCos,
    const vector& point0,
    const vector& normal0,

    const vector& point1,
    const vector& normal1
) const
{
    //- Hits differ and angles are oppositeish and
    //  hits have a normal distance
    vector d = point1-point0;
    scalar magD = mag(d);

    if (magD > mergeDistance())
    {
        scalar cosAngle = (normal0 & normal1);

        vector avg = Zero;
        if (cosAngle < (-1+planarCos))
        {
            // Opposite normals
            avg = 0.5*(normal0-normal1);
        }
        else if (cosAngle > (1-planarCos))
        {
            avg = 0.5*(normal0+normal1);
        }

        if (avg != vector::zero)
        {
            avg /= mag(avg);
            d /= magD;

            // Check average normal with respect to intersection locations
            if (mag(avg&d) > Foam::cos(degToRad(45)))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

bool Foam::meshRefinement::checkProximity
(
    const scalar planarCos,
    const label nAllowRefine,

    const label surfaceLevel,       // current intersection max level
    const vector& surfaceLocation,  // current intersection location
    const vector& surfaceNormal,    // current intersection normal

    const label cellI,

    label& cellMaxLevel,        // cached max surface level for this cell
    vector& cellMaxLocation,    // cached surface normal for this cell
    vector& cellMaxNormal,      // cached surface normal for this cell

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    // Test if surface applicable
    if (surfaceLevel > cellLevel[cellI])
    {
        if (cellMaxLevel == -1)
        {
            // First visit of cell. Store
            cellMaxLevel = surfaceLevel;
            cellMaxLocation = surfaceLocation;
            cellMaxNormal = surfaceNormal;
        }
        else
        {
            // Second or more visit.
            // Check if
            //  - different location
            //  - opposite surface

            bool closeSurfaces = isNormalGap
            (
                planarCos,
                cellMaxLocation,
                cellMaxNormal,
                surfaceLocation,
                surfaceNormal
            );

            // Set normal to that of highest surface. Not really necessary
            // over here but we reuse cellMax info when doing coupled faces.
            if (surfaceLevel > cellMaxLevel)
            {
                cellMaxLevel = surfaceLevel;
                cellMaxLocation = surfaceLocation;
                cellMaxNormal = surfaceNormal;
            }


            if (closeSurfaces)
            {
                //Pout<< "Found gap:" << nl
                //    << "    location:" << surfaceLocation
                //    << "\tnormal  :" << surfaceNormal << nl
                ///    << "    location:" << cellMaxLocation
                //    << "\tnormal  :" << cellMaxNormal << nl
                //    << "\tcos     :" << (surfaceNormal&cellMaxNormal) << nl
                //    << endl;

                return markForRefine
                (
                    surfaceLevel,   // mark with any non-neg number.
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                );
            }
        }
    }

    // Did not reach refinement limit.
    return true;
}


void Foam::meshRefinement::proximityIntersectionChecks
(
    const refinementParameters& refineParams,
    const boolList& checkCells,
    const label nRays,
    const label maxGapCells,
    boolList& possibles,
    boolList& refine
)
{
    const scalar edge0Len = meshCutter_.level0EdgeLength();
    const labelList& cellLevel = meshCutter_.cellLevel();

    const scalar proximityAngle = refineParams.proximityAngle();

    DynamicList<label> intersectionCells(mesh_.nCells());
    forAll(mesh_.cells(), cellI)
    {
        if (checkCells[cellI])
        {
            intersectionCells.append(cellI);
        }
    }
    intersectionCells.shrink();

    labelList nHitRays(mesh_.nCells(), 0);
    scalarField maxSpanDotProd(mesh_.nCells(),-GREAT);

    for (int j = 0; j < nRays; j++)
    {
        scalar theta = (j*Foam::constant::mathematical::pi)/nRays;
        for (int k = 0; k < (j == 0 ? 1 : nRays); k++)
        {
            scalar psi = (k*Foam::constant::mathematical::pi)/nRays;

            DynamicList<label>
                updatedIntersectionCells(intersectionCells.size());
            forAll(intersectionCells, i)
            {
                label cellI = intersectionCells[i];
                if (!refine[cellI])
                {
                    updatedIntersectionCells.append(cellI);
                }
            }

            updatedIntersectionCells.shrink();

            pointField start(updatedIntersectionCells.size());
            pointField end(updatedIntersectionCells.size());

            forAll(updatedIntersectionCells, i)
            {
                label cellI = updatedIntersectionCells[i];
                point cc = mesh_.cellCentres()[cellI];

                scalar spacing = maxGapCells * edge0Len
                    / pow(2., cellLevel[cellI]);

                start[i] = cc;
                end[i] = cc;

                scalar dx = spacing*Foam::sin(theta)*Foam::cos(psi);
                scalar dy = spacing*Foam::sin(theta)*Foam::sin(psi);
                scalar dz = spacing*Foam::cos(theta);

                start[i].x() += dx;
                start[i].y() += dy;
                start[i].z() += dz;

                end[i].x() -= dx;
                end[i].y() -= dy;
                end[i].z() -= dz;
            }

            List<List<pointIndexHit>> allHits;
            List<List<vector>> surfaceNormal;
            labelListList surfaceRegion;

            surfaces_.findAllProximityIntersections
            (
                start,
                end,
                allHits,
                surfaceNormal,
                surfaceRegion
            );

            forAll(updatedIntersectionCells, i)
            {
                label cellI = updatedIntersectionCells[i];
                point cc = mesh_.cellCentres()[cellI];

                if (allHits[i].size() > 1)
                {
                    vector rayVector = end[i] - start[i];
                    rayVector /= (mag(rayVector) + SMALL);

                    SortableList<scalar> rayDist(allHits[i].size());
                    forAll(allHits[i], l)
                    {
                        rayDist[l] =  (allHits[i][l].hitPoint() - start[i])
                            & rayVector;
                    }
                    rayDist.sort();
                    const labelList& indices = rayDist.indices();

                    for (int index = 0; index < indices.size()-1; index++)
                    {
                        label l = indices[index];
                        label m = indices[index+1];


                        label region0 = surfaceRegion[i][l];
                        label region1 = surfaceRegion[i][m];

                        label region0Incr = surfaces_.proxLevelIncr()[region0];
                        label region1Incr = surfaces_.proxLevelIncr()[region1];

                        if
                        (
                            region0Incr < 0 || region1Incr < 0
                            || (region0Incr + region1Incr) == 0
                        )
                        {
                            continue;
                        }

                        label prox0Dir = surfaces_.proxDir()[region0];
                        label prox1Dir = surfaces_.proxDir()[region1];

                        vector interVec = allHits[i][l].hitPoint()
                            - allHits[i][m].hitPoint();
                        vector sN0 = surfaceNormal[i][l] /
                            (mag(surfaceNormal[i][l]) + SMALL);
                        vector sN1 = surfaceNormal[i][m] /
                            (mag(surfaceNormal[i][m]) + SMALL);
                        scalar sN0iv(sN0 & interVec);
                        scalar sN1iv(sN1 & interVec);
                        if (mag(sN0 & sN1) > proximityAngle)
                        {
                            if (prox0Dir == -1 && prox0Dir == -1)
                            {
                                if (sN1iv > 0 && sN0iv < 0)
                                {
                                    continue;
                                }
                            }
                            else if (prox0Dir == 1 && prox1Dir == 1)
                            {
                                if (sN1iv < 0 && sN0iv > 0)
                                {
                                    continue;
                                }
                            }
                            else if (prox0Dir == 1 && prox1Dir == -1)
                            {
                                if (sN1iv > 0 && sN0iv > 0)
                                {
                                    continue;
                                }
                            }
                            else if (prox0Dir == -1 && prox1Dir == 1)
                            {
                                if (sN1iv < 0 && sN0iv < 0)
                                {
                                    continue;
                                }
                            }

                            scalar gap = 0.5*(mag(sN0iv)+mag(sN1iv));
                            label maxLevel0 = -1;
                            label maxGapCells0 = -1;
                            if (region0Incr > 0)
                            {
                                maxLevel0 = surfaces_.maxLevel()[region0]
                                    + region0Incr;
                                maxGapCells0 =
                                    surfaces_.proxMaxCells()[region0];
                            }

                            label maxLevel1 = -1;
                            label maxGapCells1 = -1;
                            if (region1Incr > 0)
                            {
                                maxLevel1 = surfaces_.maxLevel()[region1]
                                    + region1Incr;
                                maxGapCells1 =
                                    surfaces_.proxMaxCells()[region1];
                            }

                            label maxCellsAcrossGap =
                                max(maxGapCells0, maxGapCells1);

                            label maxProxLevel = max(maxLevel0, maxLevel1);
                            scalar minSpacing = edge0Len
                                / pow(2., maxProxLevel);

                            scalar spacing = edge0Len
                                / pow(2., cellLevel[cellI]);

                            if
                            (
                                minSpacing < gap
                                && maxCellsAcrossGap*spacing > gap
                                && cellLevel[cellI] < maxProxLevel
                            )
                            {
                                vector ccHitsVector0 =
                                    allHits[i][l].hitPoint() - cc;
                                vector ccHitsVector1 =
                                    allHits[i][m].hitPoint() - cc;
                                scalar hitDistFromCC0 =
                                    ccHitsVector0 & rayVector;
                                scalar hitDistFromCC1 =
                                    ccHitsVector1 & rayVector;

                                if
                                (
                                    (hitDistFromCC0*hitDistFromCC1) < 0.
                                    || mag(hitDistFromCC0) < spacing
                                    || mag(hitDistFromCC1) < spacing
                                )
                                {
                                    nHitRays[cellI]++;
                                    vector unitInterVec = interVec
                                        / (mag(interVec) + SMALL);
                                    maxSpanDotProd[cellI] = max
                                    (
                                        maxSpanDotProd[cellI],
                                        mag(unitInterVec & sN0)
                                    );
                                    maxSpanDotProd[cellI] = max
                                    (
                                        maxSpanDotProd[cellI],
                                        mag(unitInterVec & sN1)
                                    );
                                    if
                                    (
                                        nHitRays[cellI] > 1
                                        && maxSpanDotProd[cellI] > 0.707
                                    )
                                    {
                                        refine[cellI]  = true;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    forAll(nHitRays, cellI)
    {
        if (nHitRays[cellI] == 1)
        {
            possibles[cellI] = true;
        }
    }
}

Foam::label Foam::meshRefinement::markProximityRefinement
(
    const scalar planarCos,
    const label nAllowRefine,
    const labelList& neiLevel,
    const pointField& neiCc,

    labelList& refineCell,
    label& nRefine
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();

    label oldNRefine = nRefine;

    // 1. local test: any cell on more than one surface gets refined
    // (if its current level is < max of the surface max level)

    // 2. 'global' test: any cell on only one surface with a neighbour
    // on a different surface gets refined (if its current level etc.)


    // Collect candidate faces (i.e. intersecting any surface and
    // owner/neighbour not yet refined.
    labelList testFaces(getRefineCandidateFaces(refineCell));

    // Collect segments
    pointField start(testFaces.size());
    pointField end(testFaces.size());
    labelList minLevel(testFaces.size());

    calcCellCellRays
    (
        mesh_.cellCentres(),
        neiCc,
        neiLevel,
        testFaces,
        start,
        end,
        minLevel
    );


    // Test for all intersections (with surfaces of higher gap level than
    // minLevel) and cache per cell the max surface level and the local normal
    // on that surface.
    labelList cellMaxLevel(mesh_.nCells(), -1);
    vectorField cellMaxNormal(mesh_.nCells(), Zero);
    pointField cellMaxLocation(mesh_.nCells(), Zero);

    {
        // Per segment the normals of the surfaces
        List<vectorList> surfaceLocation;
        List<vectorList> surfaceNormal;
        // Per segment the list of levels of the surfaces
        labelListList surfaceLevel;

        surfaces_.findAllHigherIntersections
        (
            start,
            end,
            minLevel,           // gap level of surface has to be bigger
                                // than min level of neighbouring cells

            surfaces_.gapLevel(),

            surfaceLocation,
            surfaceNormal,
            surfaceLevel
        );
        // Clear out unnecessary data
        start.clear();
        end.clear();
        minLevel.clear();

        //// Extract per cell information on the surface with the highest max
        //OBJstream str
        //(
        //    mesh_.time().path()
        //  / "findAllHigherIntersections_"
        //  + mesh_.time().timeName()
        //  + ".obj"
        //);
        //// All intersections
        //OBJstream str2
        //(
        //    mesh_.time().path()
        //  / "findAllHigherIntersections2_"
        //  + mesh_.time().timeName()
        //  + ".obj"
        //);

        forAll(testFaces, i)
        {
            label faceI = testFaces[i];
            label own = mesh_.faceOwner()[faceI];

            const labelList& fLevels = surfaceLevel[i];
            const vectorList& fPoints = surfaceLocation[i];
            const vectorList& fNormals = surfaceNormal[i];

            forAll(fLevels, hitI)
            {
                checkProximity
                (
                    planarCos,
                    nAllowRefine,

                    fLevels[hitI],
                    fPoints[hitI],
                    fNormals[hitI],

                    own,
                    cellMaxLevel[own],
                    cellMaxLocation[own],
                    cellMaxNormal[own],

                    refineCell,
                    nRefine
                );
            }

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];

                forAll(fLevels, hitI)
                {
                    checkProximity
                    (
                        planarCos,
                        nAllowRefine,

                        fLevels[hitI],
                        fPoints[hitI],
                        fNormals[hitI],

                        nei,
                        cellMaxLevel[nei],
                        cellMaxLocation[nei],
                        cellMaxNormal[nei],

                        refineCell,
                        nRefine
                    );
                }
            }
        }
    }

    // 2. Find out a measure of surface curvature and region edges.
    // Send over surface region and surface normal to neighbour cell.

    labelList neiBndMaxLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiBndMaxLocation(mesh_.nFaces()-mesh_.nInternalFaces());
    vectorField neiBndMaxNormal(mesh_.nFaces()-mesh_.nInternalFaces());

    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        label bFaceI = faceI-mesh_.nInternalFaces();
        label own = mesh_.faceOwner()[faceI];

        neiBndMaxLevel[bFaceI] = cellMaxLevel[own];
        neiBndMaxLocation[bFaceI] = cellMaxLocation[own];
        neiBndMaxNormal[bFaceI] = cellMaxNormal[own];
    }
    syncTools::swapBoundaryFaceList(mesh_, neiBndMaxLevel);
    syncTools::swapBoundaryFaceList(mesh_, neiBndMaxLocation);
    syncTools::swapBoundaryFaceList(mesh_, neiBndMaxNormal);

    // Loop over all faces. Could only be checkFaces.. except if they're coupled

    // Internal faces
    for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label nei = mesh_.faceNeighbour()[faceI];

        if (cellMaxLevel[own] != -1 && cellMaxLevel[nei] != -1)
        {
            // Have valid data on both sides. Check planarCos.
            if
            (
                isNormalGap
                (
                    planarCos,
                    cellMaxLocation[own],
                    cellMaxNormal[own],
                    cellMaxLocation[nei],
                    cellMaxNormal[nei]
                )
            )
            {
                // See which side to refine
                if (cellLevel[own] < cellMaxLevel[own])
                {
                    if
                    (
                        !markForRefine
                        (
                            cellMaxLevel[own],
                            nAllowRefine,
                            refineCell[own],
                            nRefine
                        )
                    )
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining since reaching my cell"
                                << " limit of " << mesh_.nCells()+7*nRefine
                                << endl;
                        }
                        break;
                    }
                }

                if (cellLevel[nei] < cellMaxLevel[nei])
                {
                    if
                    (
                        !markForRefine
                        (
                            cellMaxLevel[nei],
                            nAllowRefine,
                            refineCell[nei],
                            nRefine
                        )
                    )
                    {
                        if (debug)
                        {
                            Pout<< "Stopped refining since reaching my cell"
                                << " limit of " << mesh_.nCells()+7*nRefine
                                << endl;
                        }
                        break;
                    }
                }
            }
        }
    }
    // Boundary faces
    for (label faceI = mesh_.nInternalFaces(); faceI < mesh_.nFaces(); faceI++)
    {
        label own = mesh_.faceOwner()[faceI];
        label bFaceI = faceI - mesh_.nInternalFaces();

        if (cellLevel[own] < cellMaxLevel[own] && neiBndMaxLevel[bFaceI] != -1)
        {
            // Have valid data on both sides. Check planarCos.
            if
            (
                isNormalGap
                (
                    planarCos,
                    cellMaxLocation[own],
                    cellMaxNormal[own],
                    neiBndMaxLocation[bFaceI],
                    neiBndMaxNormal[bFaceI]
                )
            )
            {
                if
                (
                    !markForRefine
                    (
                        cellMaxLevel[own],
                        nAllowRefine,
                        refineCell[own],
                        nRefine
                    )
                )
                {
                    if (debug)
                    {
                        Pout<< "Stopped refining since reaching my cell"
                            << " limit of " << mesh_.nCells()+7*nRefine
                            << endl;
                    }
                    break;
                }
            }
        }
    }

    if
    (
        returnReduce(nRefine, sumOp<label>())
      > returnReduce(nAllowRefine, sumOp<label>())
    )
    {
        Info<< "Reached refinement limit." << endl;
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}


// Mark cells for proximity intersection based refinement.
Foam::label Foam::meshRefinement::markProximityIncrementRefinement
(
    const refinementParameters& refineParams,
    const label nAllowRefine,

    label iter,
    labelList& refineCell,
    label& nRefine
)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();
    const scalar nRays = refineParams.nProxRays();
    const labelList proxSurfaces = surfaces_.proximitySurfaces();

    label maxIncr = -1;
    forAll(surfaces_.proxLevelIncr(),i)
    {
        maxIncr = max(maxIncr,surfaces_.proxLevelIncr()[i]);
    }
    //just in case the geometry gets partitioned one day
    reduce(maxIncr,maxOp<label>());

    if (maxIncr == 0)
    {
        return 0;
    }

    label maxGapCells = -1;
    forAll(surfaces_.proxMaxCells(),i)
    {
        maxGapCells = max(maxGapCells,surfaces_.proxMaxCells()[i]);
    }
    //just in case the geometry gets partitioned one day
    reduce(maxGapCells,maxOp<label>());

    label oldNRefine = nRefine;

    //Determine boundary cells
    labelHashSet boundaryCells(mesh_.nCells()/100);
    {
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
        calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc);

        labelList intersectionFaces(intersectedFaces());

        // Collect segments we want to test for
        pointField start(intersectionFaces.size());
        pointField end(intersectionFaces.size());
        const pointField& cellCentres = mesh_.cellCentres();

        forAll(intersectionFaces, i)
        {
            label faceI = intersectionFaces[i];
            start[i] = cellCentres[mesh_.faceOwner()[faceI]];

            if (mesh_.isInternalFace(faceI))
            {
                end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
            }
            else
            {
                end[i] = neiCc[faceI-mesh_.nInternalFaces()];
            }
        }

        // Do tests in one go
        labelList surfaceHit;
        List<pointIndexHit> surfaceHitInfo;
        surfaces_.findAnyIntersection
        (
            proxSurfaces,
            start,
            end,
            surfaceHit,
            surfaceHitInfo
         );

        labelList pRegions(surfaceHitInfo.size(),-1);

        forAll(proxSurfaces, sI)
        {
            DynamicList<pointIndexHit> localHits;
            label surfI = proxSurfaces[sI];

            forAll(surfaceHit, i)
            {
                if (surfaceHit[i] == sI)
                {
                    localHits.append(surfaceHitInfo[i]);
                }
            }
            labelList localRegions;
            label geomI = surfaces_.surfaces()[surfI];
            surfaces_.geometry()[geomI].getRegion(localHits, localRegions);

            label localI = 0;
            forAll(surfaceHit, i)
            {
                if (surfaceHit[i] == sI)
                {
                    pRegions[i] = localRegions[localI];
                    localI++;
                }
            }
        }

        forAll(surfaceHitInfo, i)
        {
            if (surfaceHitInfo[i].hit())
            {
                label global =
                    surfaces_.globalRegion(proxSurfaces[surfaceHit[i]], pRegions[i]);

                if (surfaces_.proxLevelIncr()[global] >= 0)
                {
                    label faceI = intersectionFaces[i];
                    label own = mesh_.faceOwner()[faceI];
                    boundaryCells.insert(own);

                    if (mesh_.isInternalFace(faceI))
                    {
                        label nei = mesh_.faceNeighbour()[faceI];
                        boundaryCells.insert(nei);
                    }
                }
            }
        }
    }

    labelList bCells = boundaryCells.toc();
    boolList checkCells(mesh_.nCells(), false);
    forAll(bCells, i)
    {
        checkCells[bCells[i]] = true;
    }

    //Expand checked cells from boundary intersected ones
    boolList updatedCheckCells = checkCells;

    label nBuffer = ((maxGapCells+1)/2) - 1;
     for (int i = 0; i < nBuffer; i++)
    {
        boolList neiCheckCells(mesh_.nFaces()-mesh_.nInternalFaces(),false);
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            const labelUList& faceCells = pp.faceCells();

            label bFaceI = pp.start()-mesh_.nInternalFaces();

            if (pp.coupled())
            {
                forAll(faceCells, i)
                {
                    neiCheckCells[bFaceI] = checkCells[faceCells[i]];
                    bFaceI++;
                }
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiCheckCells);

        forAll(mesh_.faces(), faceI)
        {
            label own = mesh_.faceOwner()[faceI];

            if (mesh_.isInternalFace(faceI))
            {
                label nei = mesh_.faceNeighbour()[faceI];
                if (checkCells[own])
                {
                    updatedCheckCells[nei] = true;
                }
                if (checkCells[nei])
                {
                    updatedCheckCells[own] = true;
                }
            }
            else
            {
                if (neiCheckCells[faceI-mesh_.nInternalFaces()])
                {
                    updatedCheckCells[own] = true;
                }
            }
        }
        checkCells = updatedCheckCells;
    }

     boolList possibles(mesh_.nCells(), false);
     boolList refine(mesh_.nCells(),false);

     proximityIntersectionChecks
     (
         refineParams,
         checkCells,
         nRays,
         maxGapCells,
         possibles,
         refine
     );

     //second more expensive check on subset of checkCells
     boolList dummy(mesh_.nCells(), false);
     proximityIntersectionChecks
     (
         refineParams,
         possibles,
         2*nRays+1,
         maxGapCells,
         dummy,
         refine
     );

    forAll(refine, cellI)
    {
        if (refine[cellI])
        {
            if
            (
                !markForRefine
                (
                    1,
                    nAllowRefine,
                    refineCell[cellI],
                    nRefine
                 )
             )
            {
                break;
            }
        }
    }

    return returnReduce(nRefine-oldNRefine, sumOp<label>());
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Calculate list of cells to refine. Gets for any edge (start - end)
// whether it intersects the surface. Does more accurate test and checks
// the wanted level on the surface intersected.
// Does approximate precalculation of how many cells can be refined before
// hitting overall limit maxGlobalCells.
Foam::labelList Foam::meshRefinement::refineCandidates
(
    const label iterationNum,
    const refinementParameters& refineParams,

    const bool initialRefinement,
    const bool featureRefinement,
    const bool featureDistanceRefinement,
    const bool internalRefinement,
    const bool surfaceRefinement,
    const bool curvatureRefinement,
    const bool proximityRefinement,
    const bool smallFeatureRefinement,
    const bool gapRefinement,
    const bool bigGapRefinement,
    const bool spreadGapSize
)
{
    const scalar curvature = refineParams.curvature();
    const scalar planarAngle = refineParams.planarAngle();
    const label maxGlobalCells = refineParams.maxGlobalCells();

    label totNCells = mesh_.globalData().nTotalCells();

    labelList cellsToRefine;

    if (totNCells >= maxGlobalCells)
    {
        Info<< "No cells marked for refinement since reached limit "
            << maxGlobalCells << '.' << endl;
    }
    else
    {
        // Every cell I refine adds 7 cells. Estimate the number of cells
        // I am allowed to refine.
        // Assume perfect distribution so can only refine as much the fraction
        // of the mesh I hold. This prediction step prevents us having to do
        // lots of reduces to keep count of the total number of cells selected
        // for refinement.

        //scalar fraction = scalar(mesh_.nCells())/totNCells;
        //label myMaxCells = label(maxGlobalCells*fraction);
        //label nAllowRefine = (myMaxCells - mesh_.nCells())/7;
        ////label nAllowRefine = (maxLocalCells - mesh_.nCells())/7;
        //
        //Pout<< "refineCandidates:" << nl
        //    << "    total cells:" << totNCells << nl
        //    << "    local cells:" << mesh_.nCells() << nl
        //    << "    local fraction:" << fraction << nl
        //    << "    local allowable cells:" << myMaxCells << nl
        //    << "    local allowable refinement:" << nAllowRefine << nl
        //    << endl;

        //- Disable refinement shortcut. nAllowRefine is per processor limit.
        label nAllowRefine = labelMax / Pstream::nProcs();

        // Marked for refinement (>= 0) or not (-1). Actual value is the
        // index of the surface it intersects.
        labelList refineCell(mesh_.nCells(), -1);
        label nRefine = 0;


        // Swap neighbouring cell centres and cell level
        labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
        pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
        calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc);

        const scalar planarCos = Foam::cos(degToRad(planarAngle));

        // Initial refinment where surface bounding box smaller than cell size
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (initialRefinement)
        {
            label nFeatures = markInitialRefinement
            (
                nAllowRefine,

                refineCell,
                nRefine
            );

            Info<< "Marked for refinement due to surface bounding box seeding  "
                << ": " << nFeatures << " cells."  << endl;
        }

        // Cells pierced by feature lines
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (featureRefinement)
        {
            label nFeatures = markFeatureRefinement
            (
                refineParams,
                nAllowRefine,

                refineCell,
                nRefine
            );

            Info<< "Marked for refinement due to explicit features             "
                << ": " << nFeatures << " cells."  << endl;
        }

        // Inside distance-to-feature shells
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (featureDistanceRefinement)
        {
            label nShell = markInternalDistanceToFeatureRefinement
            (
                nAllowRefine,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to distance to explicit features "
                ": " << nShell << " cells."  << endl;
        }

        // Inside refinement shells
        // ~~~~~~~~~~~~~~~~~~~~~~~~

        if (internalRefinement)
        {
            label nShell = markInternalRefinement
            (
                nAllowRefine,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to refinement shells             "
                << ": " << nShell << " cells."  << endl;
        }

        // Refinement based on intersection of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (surfaceRefinement)
        {
            label nSurf = markSurfaceRefinement
            (
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to surface intersection          "
                << ": " << nSurf << " cells."  << endl;

            // Refine intersected-cells only inside gaps. See
            // markInternalGapRefinement to refine all cells inside gaps.
            if
            (
                planarCos >= -1
             && planarCos <= 1
             && max(shells_.maxGapLevel()) > 0
            )
            {
                label nGapSurf = markSurfaceGapRefinement
                (
                    planarCos,
                    nAllowRefine,
                    neiLevel,
                    neiCc,

                    refineCell,
                    nRefine
                );
                Info<< "Marked for refinement due to surface intersection"
                    << " (at gaps)"
                    << ": " << nGapSurf << " cells."  << endl;
            }
        }

        // Refinement based on curvature of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if
        (
            curvatureRefinement
         && curvature != -GREAT
         && (surfaces_.minLevel() != surfaces_.maxLevel())
        )
        {
            label nCurv = markSurfaceCurvatureRefinement
            (
                curvature,
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to curvature/regions             "
                << ": " << nCurv << " cells."  << endl;
        }

        // Refinement based on intersection of surface
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (proximityRefinement)
        {
            label nProx = markProximityIncrementRefinement
            (
                refineParams,
                nAllowRefine,

                iterationNum,
                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to surface proximity : "
                << nProx << " cells."  << endl;
        }


        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if
        (
            smallFeatureRefinement
         && (planarCos >= -1 && planarCos <= 1)
         && max(shells_.maxGapLevel()) > 0
        )
        {
            label nGap = markSmallFeatureRefinement
            (
                planarCos,
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to close opposite surfaces       "
                << ": " << nGap << " cells."  << endl;
        }


        // Refinement based on gap (only neighbouring cells)
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if
        (
            gapRefinement
         && (planarCos >= -1 && planarCos <= 1)
         && (max(surfaces_.gapLevel()) > -1)
        )
        {
            Info<< "Specified gap level : " << max(surfaces_.gapLevel())
                << ", planar angle " << planarAngle << endl;

            label nGap = markProximityRefinement
            (
                planarCos,
                nAllowRefine,
                neiLevel,
                neiCc,

                refineCell,
                nRefine
            );
            Info<< "Marked for refinement due to close opposite surfaces       "
                << ": " << nGap << " cells."  << endl;
        }


        // Refinement based on gaps larger than cell size
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if
        (
            bigGapRefinement
         && (planarCos >= -1 && planarCos <= 1)
         && max(shells_.maxGapLevel()) > 0
        )
        {
            // Refine based on gap information provided by shell and nearest
            // surface
            labelList numGapCells;
            scalarField gapSize;
            label nGap = markInternalGapRefinement
            (
                planarCos,
                spreadGapSize,
                nAllowRefine,

                refineCell,
                nRefine,
                numGapCells,
                gapSize
            );
            Info<< "Marked for refinement due to opposite surfaces"
                << "             "
                << ": " << nGap << " cells."  << endl;
        }


        // Limit refinement
        // ~~~~~~~~~~~~~~~~

        {
            label nUnmarked = unmarkInternalRefinement(refineCell, nRefine);
            if (nUnmarked > 0)
            {
                Info<< "Unmarked for refinement due to limit shells"
                    << "                : " << nUnmarked << " cells."  << endl;
            }
        }

        // Pack cells-to-refine
        // ~~~~~~~~~~~~~~~~~~~~

        cellsToRefine.setSize(nRefine);
        nRefine = 0;

        forAll(refineCell, cellI)
        {
            if (refineCell[cellI] != -1)
            {
                cellsToRefine[nRefine++] = cellI;
            }
        }
    }

    return cellsToRefine;
}


void Foam::meshRefinement::markWrapFaces
(
    const pointField& locationsInMesh
)
{
    labelList wrapSurfaces(surfaces().wrapLevelSurfaces());
    if (wrapSurfaces.size() == 0)
    {
        return;
    }

    Info<<"Marking faces for wrapping"<<endl;

    if (debug)
    {
        mesh_.write();
        DynamicList<label> wrapFaces;
        forAll(mesh_.faces(), facei)
        {
            if (wrapIndex_[facei] != -1)
            {
                wrapFaces.append(facei);
            }
        }
        faceSet wrapSet
        (
            mesh_,
            "preWrapFaces",
            labelHashSet(wrapFaces)
        );
        wrapSet.instance() = timeName();
        wrapSet.write();
    }

    const pointField& cellCentres = mesh_.cellCentres();
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(cellCentres, neiLevel, neiCc);

    labelList intersectionFaces(intersectedFaces());

    // Collect segments we want to test for
    pointField start(intersectionFaces.size());
    pointField end(intersectionFaces.size());

    forAll(intersectionFaces, i)
    {
        label faceI = intersectionFaces[i];
        start[i] = cellCentres[mesh_.faceOwner()[faceI]];

        if (mesh_.isInternalFace(faceI))
        {
            end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
        }
        else
        {
            end[i] = neiCc[faceI-mesh_.nInternalFaces()];
        }
    }

    // Do tests in one go
    labelList surfaceHit;
    List<pointIndexHit> surfaceHitInfo;
    surfaces_.findAnyIntersection
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfaces_.surfZones()),
        start,
        end,
        surfaceHit,
        surfaceHitInfo
    );


    labelList unnamedHits(mesh_.nFaces(), -1);
    forAll(intersectionFaces, i)
    {
        if (surfaceHit[i] != -1)
        {
            label faceI = intersectionFaces[i];
            unnamedHits[faceI] = surfaceHit[i];
        }
    }

    closeSingleCellGaps
    (
        wrapSurfaces,
        locationsInMesh,
        unnamedHits,
        surfaces().wrapLevel(),
        wrapIndex_
    );

    if (debug)
    {
        Info<<"Debug: Writing wrapped and baffle faces"<<endl;

        DynamicList<label> wrapFaces;
        DynamicList<label> interFaces;

        forAll(mesh_.faces(), facei)
        {
            if (wrapIndex_[facei] != -1)
            {
                wrapFaces.append(facei);
            }
            if (surfaceIndex_[facei] != -1)
            {
                interFaces.append(facei);
            }
        }

        faceSet wrapSet
        (
            mesh_,
            "wrapFaces",
            labelHashSet(wrapFaces)
        );

        faceSet interSet
        (
            mesh_,
            "interFaces",
            labelHashSet(interFaces)
        );

        interSet.instance() = timeName();
        wrapSet.instance() = timeName();
        interSet.write();
        wrapSet.write();

        Time& runTime = const_cast<Time&>(mesh_.time());
        runTime++;
    }

    return;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::meshRefinement::refine
(
    const labelList& cellsToRefine,
    const pointField& locationsInMesh
)
{
    // Mesh changing engine.
    polyTopoChange meshMod(mesh_);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (no inflation), return map from old to new mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    // Update fields
    mesh_.updateMesh(map);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    // Reset the instance for if in overwrite mode
    mesh_.setInstance(timeName());

    // Update intersection info
    updateMesh(map, getChangedFaces(map, cellsToRefine));

    if (wrapActive_)
    {
        // Reset wrapped faces if normal surface intersection detected
        if (surfaces().wrapLevelSurfaces().size() > 0)
        {
            forAll(wrapIndex_, facei)
            {
                if (surfaceIndex_[facei] != -1 && wrapIndex_[facei] != -1)
                {
                    wrapIndex_[facei] = -1;
                }
            }
        }

        //Optional marking of wrapped faces
        markWrapFaces(locationsInMesh);
    }

    return map;
}


Foam::labelList Foam::meshRefinement::refineGaps()
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();
    const pointField& cellCentres = mesh_.cellCentres();
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Mark all points and edges on baffle patches (so not on any inlets,
    // outlets etc.)
    boolList isBoundaryPoint(mesh_.nPoints(), false);
    boolList isBoundaryEdge(mesh_.nEdges(), false);
    boolList isBoundaryFace(mesh_.nFaces(), false);

    // Fill boundary data. All elements on meshed patches get marked.
    // Get the labels of added patches.
    labelList neiLevel(mesh_.nFaces()-mesh_.nInternalFaces());
    pointField neiCc(mesh_.nFaces()-mesh_.nInternalFaces());
    calcNeighbourData(mesh_.cellCentres(), neiLevel, neiCc);

    labelList intersectionFaces(intersectedFaces());

    // Collect segments we want to test for
    pointField start(intersectionFaces.size());
    pointField end(intersectionFaces.size());

    forAll(intersectionFaces, i)
    {
        label faceI = intersectionFaces[i];
        start[i] = cellCentres[mesh_.faceOwner()[faceI]];

        if (mesh_.isInternalFace(faceI))
        {
            end[i] = cellCentres[mesh_.faceNeighbour()[faceI]];
        }
        else
        {
            end[i] = neiCc[faceI-mesh_.nInternalFaces()];
        }
    }

    // Do tests in one go
    labelList surfaceHit;
    List<pointIndexHit> surfaceHitInfo;
    surfaces_.findAnyIntersection
    (
        surfaceZonesInfo::getUnnamedSurfaces(surfaces_.surfZones()),
        start,
        end,
        surfaceHit,
        surfaceHitInfo
    );

    forAll(intersectionFaces, i)
    {
        if (surfaceHit[i] != -1)
        {
            label faceI = intersectionFaces[i];

            markBoundaryFace
            (
                faceI,
                isBoundaryFace,
                isBoundaryEdge,
                isBoundaryPoint
             );
        }

    }

   labelList adaptPatchIDs(unmeshedPatches());

    forAll(adaptPatchIDs, i)
    {
        label patchI = adaptPatchIDs[i];

        const polyPatch& pp = patches[patchI];

        label faceI = pp.start();

        forAll(pp, j)
        {
            markBoundaryFace
            (
                faceI,
                isBoundaryFace,
                isBoundaryEdge,
                isBoundaryPoint
            );

            faceI++;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        isBoundaryPoint,
        orEqOp<bool>(),
        false              // null value
    );

    syncTools::syncEdgeList
    (
        mesh_,
        isBoundaryEdge,
        orEqOp<bool>(),
        false              // null value
    );

    syncTools::syncFaceList
    (
        mesh_,
        isBoundaryFace,
        orEqOp<bool>()
    );

    // For each cell count the number of anchor points that are on
    // the boundary:
    DynamicList<label> markedCells(mesh_.nCells()/10000);
    DynamicList<label> dynCPoints;

    forAll(cellLevel, cellI)
    {
        const labelList& cPoints = mesh_.cellPoints(cellI, dynCPoints);

        // Get number of anchor points (pointLevel <= cellLevel)

        label nBoundaryAnchors = 0;
//        label nNonAnchorBoundary = 0;

        forAll(cPoints, i)
        {
            label pointI = cPoints[i];

            if (pointLevel[pointI] <= cellLevel[cellI])
            {
                // Anchor point
                if (isBoundaryPoint[pointI])
                {
                    nBoundaryAnchors++;
                }
            }
            else if (isBoundaryPoint[pointI])
            {
//                nNonAnchorBoundary++;
            }
        }

        if (nBoundaryAnchors == 8)
        {
            markedCells.append(cellI);
        }
    }

    return markedCells.shrink();
}


// Load balancing
Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::meshRefinement::balance
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    labelList& cellsToRefine,
    const refinementParameters& refineParams
)
{
    const scalar maxLoadUnbalance = refineParams.maxLoadUnbalance();
    const scalar maxCellUnbalance = refineParams.maxCellUnbalance();

    autoPtr<mapDistributePolyMesh> distMap;

    if (Pstream::nProcs() > 1)
    {
        // First check if we need to balance at all. Precalculate number of
        // cells after refinement and see what maximum difference is.
        const scalar nNewCells =
            scalar(mesh_.nCells() + 7*cellsToRefine.size());
        const scalar nNewCellsAll = returnReduce(nNewCells, sumOp<scalar>());
        const scalar nIdealNewCells = nNewCellsAll / Pstream::nProcs();
        const scalar unbalance = returnReduce
        (
            mag(scalar(1.0)-nNewCells/nIdealNewCells),
            maxOp<scalar>()
        );
        // Trigger the balancing to avoid too early balancing for better
        // scaling performance.
        const scalar nNewCellsOnly = scalar(7*cellsToRefine.size());

        const label maxNewCells =
            label(returnReduce(nNewCellsOnly, maxOp<scalar>()));

        const label maxDeltaCells =
            label(mag(returnReduce(nNewCells, maxOp<scalar>())-nIdealNewCells));

        // New trigger to avoid too early balancing
        // 1. Check if globally one proc exceeds the maxCellUnbalance value
        //    related to the new added cells at the refinement loop
        // 2. Check if globally one proc exceeds the maxCellUnbalance based on
        //    the average cell count a proc should have
        if
        (
            (maxNewCells <= maxCellUnbalance)
         && (maxDeltaCells <= maxCellUnbalance)
        )
        {
            Info<< "Skipping balancing since trigger value not reached:" << "\n"
                << "    Trigger cell count: " << maxCellUnbalance << nl
                << "    Max new cell count in proc: " << maxNewCells << nl
                << "    Max difference between new cells and balanced: "
                << maxDeltaCells << nl
                << "    Max load unbalance " << maxLoadUnbalance
                << nl <<endl;
        }
        else
        {
            if (unbalance <= maxLoadUnbalance)
            {
                Info<< "Skipping balancing since max unbalance " << unbalance
                    << " is less than allowable " << maxLoadUnbalance
                    << endl;
            }
            else
            {
                scalarField cellWeights(mesh_.nCells(), 1);
                forAll(cellsToRefine, i)
                {
                    cellWeights[cellsToRefine[i]] += 7;
                }

                distMap = balance
                (
                    false,
                    false,  //keepZoneFaces
                    false,  //keepBaffles
                    cellWeights,
                    decomposer,
                    distributor
                );

                // Update cells to refine
                distMap().distributeCellIndices(cellsToRefine);

                Info<< "Balanced mesh in = "
                    << mesh_.time().cpuTimeIncrement() << " s" << endl;

                if (debug&meshRefinement::MESH)
                {
                    Pout<< "Writing balanced " << msg
                        << " mesh to time " << timeName() << endl;
                    write
                    (
                        debugType(debug),
                        writeType(writeLevel() | WRITEMESH),
                        mesh_.time().path()/timeName()
                    );
                    Pout<< "Dumped debug data in = "
                        << mesh_.time().cpuTimeIncrement() << " s" << endl;

                    // test all is still synced across proc patches
                    checkData();
                }
            }
        }
    }

    return distMap;
}



// Do refinement of consistent set of cells followed by truncation and
// load balancing.
Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::meshRefinement::refineAndBalance
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& cellsToRefine,
    const refinementParameters& refineParams
)
{
    const pointField locationsInMesh = refineParams.locationsInMesh();

    // Do all refinement
    refine(cellsToRefine, locationsInMesh);

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing refined but unbalanced " << msg
            << " mesh to time " << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            mesh_.time().path()/timeName()
        );
        Pout<< "Dumped debug data in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // test all is still synced across proc patches
        checkData();
    }

    Info<< "Refined mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s" << endl;
    printMeshInfo(debug, "After refinement " + msg);


    // Load balancing
    // ~~~~~~~~~~~~~~

    labelList noCellsToRefine;

    auto distMap = balance
    (
        msg,
        decomposer,
        distributor,
        noCellsToRefine,    // mesh is already refined; no need to predict
        refineParams
    );

    return distMap;
}


// Do load balancing followed by refinement of consistent set of cells.
Foam::autoPtr<Foam::mapDistributePolyMesh>
Foam::meshRefinement::balanceAndRefine
(
    const string& msg,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const labelList& initCellsToRefine,
    const refinementParameters& refineParams
)
{
    const pointField locationsInMesh = refineParams.locationsInMesh();

    labelList cellsToRefine(initCellsToRefine);

    // Load balancing
    // ~~~~~~~~~~~~~~

    auto distMap = balance
    (
        msg,
        decomposer,
        distributor,
        cellsToRefine,
        refineParams
    );

    // Refinement
    // ~~~~~~~~~~

    refine(cellsToRefine,locationsInMesh);

    if (debug&meshRefinement::MESH)
    {
        Pout<< "Writing refined " << msg
            << " mesh to time " << timeName() << endl;
        write
        (
            debugType(debug),
            writeType(writeLevel() | WRITEMESH),
            mesh_.time().path()/timeName()
        );
        Pout<< "Dumped debug data in = "
            << mesh_.time().cpuTimeIncrement() << " s" << endl;

        // test all is still synced across proc patches
        checkData();
    }

    Info<< "Refined mesh in = "
        << mesh_.time().cpuTimeIncrement() << " s" << endl;

    printMeshInfo(debug, "After refinement " + msg);

    return distMap;
}


// ************************************************************************* //
