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
    (c) Creative Fields, Ltd.
    (c) 2020 Esi Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description

\*---------------------------------------------------------------------------*/

#include <stdexcept>

#include "include/demandDrivenData.H"
#include "utilities/smoothers/geometry/meshOptimizer/meshOptimizer.H"
#include "utilities/meshes/polyMeshGenAddressing/polyMeshGenAddressing.H"
#include "utilities/meshes/polyMeshGenChecks/polyMeshGenChecks.H"
#include "utilities/meshes/partTetMesh/partTetMesh.H"
#include "containers/HashTables/HashSet/HashSet.H"

#include "utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/tetMeshOptimisation.H"
#include "utilities/smoothers/geometry/meshOptimizer/boundaryLayerOptimisation/boundaryLayerOptimisation.H"
#include "utilities/boundaryLayers/refineBoundaryLayers/refineBoundaryLayers.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"

//#define DEBUGSmooth
#include "utilities/helperFunctions/helperFunctions.H"
# ifdef DEBUGSmooth
#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshOptimizer::untangleMeshFV
(
    const label maxNumGlobalIterations,
    const label maxNumIterations,
    const label maxNumSurfaceIterations,
    const label errorBufferLayers,
    const bool relaxedCheck,
    const bool relaxedBoundaryCheck,
    const bool checkWarped,
    const scalar minFaceArea
)
{
    Info<< "Starting untangling the mesh" << endl;

    pointField minErrorPts = mesh_.points();
    label globalMinCriticalFaces = labelMax;
    label globalMinBadFaces = labelMax;

    # ifdef DEBUGSmooth
    partTetMesh tm(mesh_);
    forAll(tm.tets(), tetI)
        if (tm.tets()[tetI].mag(tm.points()) < 0.0)
            Info<< "Tet " << tetI << " is inverted!" << endl;
    polyMeshGen tetPolyMesh(mesh_.returnTime());
    tm.createPolyMesh(tetPolyMesh);
    polyMeshGenModifier(tetPolyMesh).removeUnusedVertices();
    forAll(tm.smoothVertex(), pI)
        if (!tm.smoothVertex()[pI])
            Info<< "Point " << pI << " cannot be moved!" << endl;

    const VRWGraph& pTets = tm.pointTets();
    forAll(pTets, pointI)
    {
        const LongList<partTet>& tets = tm.tets();
        forAllRow(pTets, pointI, i)
            if (tets[pTets(pointI, i)].whichPosition(pointI) < 0)
                FatalError << "Wrong partTet" << abort(FatalError);

        partTetMeshSimplex simplex(tm, pointI);
    }

    boolList boundaryVertex(tetPolyMesh.points().size(), false);
    const labelList& neighbour = tetPolyMesh.neighbour();
    forAll(neighbour, faceI)
        if (neighbour[faceI] == -1)
        {
            const face& f = tetPolyMesh.faces()[faceI];

            forAll(f, pI)
                boundaryVertex[f[pI]] = true;
        }

    forAll(boundaryVertex, pI)
    {
        if (boundaryVertex[pI] && tm.smoothVertex()[pI])
            FatalErrorInFunction
                << "Boundary vertex should not be moved!"
                << abort(FatalError);
    }
    # endif

    label nBadFaces, nCriticalErrors, nGlobalIter(0), nIter;

    const faceListPMG& faces = mesh_.faces();

    if (errorBufferLayers > -1)
    {
        lockErrorFreeCells(errorBufferLayers,minFaceArea);
    }

    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
            lockedPoints.append(pointI);
    }

    labelHashSet badFaces;

    label startRelax = 4;
    bool unconverged = false;

    do
    {
        nIter = 0;

        label minNumBadFaces(labelMax), minIter(-1);
        do
        {
            if (!relaxedCheck && nGlobalIter < startRelax)
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFaces
                    (
                        mesh_,
                        checkWarped,
                        minFaceArea,
                        badFaces,
                        nCriticalErrors,
                        false,
                        &changedFace
                    );
            }
            else
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFacesRelaxed
                    (
                        mesh_,
                        minFaceArea,
                        badFaces,
                        false,
                        &changedFace
                    );
                nCriticalErrors = nBadFaces;
            }

            //Look for minimum critical errors if not error free
            if (nCriticalErrors <= globalMinCriticalFaces)
            {
                if (nCriticalErrors == globalMinCriticalFaces)
                {
                    if (nBadFaces < globalMinBadFaces)
                    {
                        globalMinBadFaces = nBadFaces;
                        minErrorPts = mesh_.points();
                    }
                }
                else
                {
                    globalMinBadFaces = nBadFaces;
                    globalMinCriticalFaces = nCriticalErrors;
                    minErrorPts = mesh_.points();
                }
            }
            else if
            (
                nCriticalErrors > max(label(1000),10*globalMinCriticalFaces)
            )
            {
                unconverged = true;
                break;
             }

            Info<< "Iteration " << nIter
                 << ". Number of bad faces is " << nBadFaces << nl << endl;

            //- perform optimisation
            if (nBadFaces == 0)
                break;

            if (nBadFaces < minNumBadFaces)
            {
                minNumBadFaces = nBadFaces;
                minIter = nIter;
            }

            //- create a tet mesh from the mesh and the labels of bad faces
            partTetMesh tetMesh
            (
                mesh_,
                lockedPoints,
                badFaces,
                (nGlobalIter / 2) + 1
            );

            //- construct tetMeshOptimisation and improve positions of
            //- points in the tet mesh
            tetMeshOptimisation tmo(tetMesh);

            tmo.optimiseUsingKnuppMetric();

            tmo.optimiseUsingMeshUntangler();

            //- Volume optimizer is slow so only active after first pass
            if (nGlobalIter > 0)
            {
                tmo.optimiseUsingVolumeOptimizer();
            }

            //- update points in the mesh from the coordinates in the tet mesh
            tetMesh.updateOrigMesh(&changedFace);

        } while
        (
            (nIter < minIter+3) && (++nIter < maxNumIterations)
            && (10*max(globalMinCriticalFaces,label(1)) >= nCriticalErrors)
        );

        if
        (
            unconverged
            || (nBadFaces == 0)
            || (++nGlobalIter >= maxNumGlobalIterations)
        )
        {
            break;
        }

        // move boundary vertices
        nIter = 0;

        while (nIter++ < maxNumSurfaceIterations)
        {
            if
            (
                !relaxedCheck && !relaxedBoundaryCheck
                && nGlobalIter < startRelax
            )
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFaces
                    (
                        mesh_,
                        checkWarped,
                        minFaceArea,
                        badFaces,
                        nCriticalErrors,
                        false,
                        &changedFace
                    );
            }
            else
            {
                nBadFaces =
                    polyMeshGenChecks::findBadFacesRelaxed
                    (
                        mesh_,
                        minFaceArea,
                        badFaces,
                        false,
                        &changedFace
                    );
                nCriticalErrors = nBadFaces;
            }

            //Look for minimum critical errors if not error free
            if (nCriticalErrors <= globalMinCriticalFaces)
            {
                if (nCriticalErrors == globalMinCriticalFaces)
                {
                    if (nBadFaces < globalMinBadFaces)
                    {
                        globalMinBadFaces = nBadFaces;
                        minErrorPts = mesh_.points();
                    }
                }
                else
                {
                    globalMinBadFaces = nBadFaces;
                    globalMinCriticalFaces = nCriticalErrors;
                    minErrorPts = mesh_.points();
                }
            }
            else if
            (
                nCriticalErrors > max(label(1000),10*globalMinCriticalFaces)
            )
            {
                unconverged = true;
                break;
            }

            Info<< "Iteration " << nIter
                 << ". Number of bad faces is " << nBadFaces << nl << endl;

            //- perform optimisation
            if (nBadFaces == 0)
            {
                break;
            }
            else if (enforceConstraints_)
            {
                const label subsetId =
                    mesh_.addPointSubset(badPointsSubsetName_);

                forAllConstIter(labelHashSet, badFaces, it)
                {
                    const face& f = faces[it.key()];
                    forAll(f, pI)
                        mesh_.addPointToSubset(subsetId, f[pI]);
                }

                WarningInFunction
                    << "Writing mesh with " << badPointsSubsetName_
                    << " subset. These points cannot be untangled"
                    << " without sacrificing geometry constraints. Exitting.."
                    << endl;

                returnReduce(label(1), sumOp<label>());
                throw std::logic_error
                (
                    "void meshOptimizer::untangleMeshFV()"
                    "Cannot untangle mesh!!"
                );
            }

            //- create tethrahedral mesh from the cells which shall be smoothed
            partTetMesh tetMesh(mesh_, lockedPoints, badFaces, 0);

            //- contruct tetMeshOptimisation
            tetMeshOptimisation tmo(tetMesh);

            if (nGlobalIter < 3)
            {
                //- the point stays in the plane determined by the point normal
                tmo.optimiseBoundaryVolumeOptimizer(3, true);
            }
            else if (nGlobalIter < 5)
            {
                //- move points without any constraints on the movement
                tmo.optimiseBoundarySurfaceLaplace();
            }
            else
            {
                //- move boundary points without any constraints
                tmo.optimiseBoundaryVolumeOptimizer(3, false);
            }

            tetMesh.updateOrigMesh(&changedFace);

        }

    } while (nBadFaces && !unconverged);

    if (nBadFaces != 0)
    {
        mesh_.points() = minErrorPts;
        label subsetId = mesh_.faceSubsetIndex("badFaces");
        if (subsetId >= 0)
            mesh_.removeFaceSubset(subsetId);
        subsetId = mesh_.addFaceSubset("badFaces");

        forAllConstIter(labelHashSet, badFaces, it)
            mesh_.addFaceToSubset(subsetId, it.key());
    }

    Info<< "Finished untangling the mesh" << endl;
}

void meshOptimizer::optimizeBoundaryLayer(const bool addBufferLayer)
{
    if (mesh_.returnTime().foundObject<IOdictionary>("meshDict"))
    {
        const dictionary& meshDict =
            mesh_.returnTime().lookupObject<IOdictionary>("meshDict");

        bool smoothLayer(false);

        if (meshDict.found("boundaryLayers"))
        {
            const dictionary& layersDict = meshDict.subDict("boundaryLayers");

            if (layersDict.found("optimiseLayer"))
                smoothLayer = readBool(layersDict.lookup("optimiseLayer"));
        }

        if (!smoothLayer)
            return;

        if (addBufferLayer)
        {
            //- create a buffer layer which will not be modified by the smoother
            refineBoundaryLayers refLayers(mesh_);

            refineBoundaryLayers::readSettings(meshDict, refLayers);

            refLayers.activateSpecialMode();

            refLayers.refineLayers();

            clearSurface();
            calculatePointLocations();
        }

        Info<< "Starting optimising boundary layer" << endl;

        const meshSurfaceEngine& mse = meshSurface();
        const labelList& faceOwner = mse.faceOwners();

        boundaryLayerOptimisation optimiser(mesh_, mse);

        boundaryLayerOptimisation::readSettings(meshDict, optimiser);

        optimiser.optimiseLayer();

        //- check if the bnd layer is tangled somewhere
        labelLongList bndLayerCells;
        const boolList& baseFace = optimiser.isBaseFace();

        # ifdef DEBUGSmooth
        const label blCellsId = mesh_.addCellSubset("blCells");
        # endif

        forAll(baseFace, bfI)
        {
            if (baseFace[bfI])
            {
                bndLayerCells.append(faceOwner[bfI]);

                # ifdef DEBUGSmooth
                mesh_.addCellToSubset(blCellsId, faceOwner[bfI]);
                # endif
            }
        }

        clearSurface();
        mesh_.clearAddressingData();

        //- lock boundary layer points, faces and cells
        lockCells(bndLayerCells);

        # ifdef DEBUGSmooth
        pointField origPoints(mesh_.points().size());
        forAll(origPoints, pI)
            origPoints[pI] = mesh_.points()[pI];
        # endif

        //- optimize mesh quality
        optimizeMeshFV(5, 1, 50, 0);

        //- untangle remaining faces and lock the boundary layer cells
        untangleMeshFV(2, 50, 0);

        # ifdef DEBUGSmooth
        forAll(vertexLocation_, pI)
        {
            if (vertexLocation_[pI] & LOCKED)
            {
                if (mag(origPoints[pI] - mesh_.points()[pI]) > SMALL)
                    FatalError << "Locked points were moved"
                               << abort(FatalError);
            }
        }
        # endif

        //- unlock bnd layer points
        removeUserConstraints();

        Info<< "Finished optimising boundary layer" << endl;
    }
}

void meshOptimizer::untangleBoundaryLayer()
{
    bool untangleLayer(true);
    if (mesh_.returnTime().foundObject<IOdictionary>("meshDict"))
    {
        const dictionary& meshDict =
            mesh_.returnTime().lookupObject<IOdictionary>("meshDict");

        if (meshDict.found("boundaryLayers"))
        {
            const dictionary& layersDict = meshDict.subDict("boundaryLayers");

            if (layersDict.found("untangleLayers"))
            {
                untangleLayer =
                    readBool(layersDict.lookup("untangleLayers"));
            }
        }
    }

    if (!untangleLayer)
    {
        labelHashSet badFaces;
        polyMeshGenChecks::checkFacePyramids(mesh_, false, VSMALL, &badFaces);

        const label nInvalidFaces =
            returnReduce(badFaces.size(), sumOp<label>());

        if (nInvalidFaces != 0)
        {
            const labelList& owner = mesh_.owner();
            const labelList& neighbour = mesh_.neighbour();

            const label badBlCellsId =
                mesh_.addCellSubset("invalidBoundaryLayerCells");

            forAllConstIter(labelHashSet, badFaces, it)
            {
                mesh_.addCellToSubset(badBlCellsId, owner[it.key()]);

                if (neighbour[it.key()] < 0)
                    continue;

                mesh_.addCellToSubset(badBlCellsId, neighbour[it.key()]);
            }

            returnReduce(label(1), sumOp<label>());

            throw std::logic_error
            (
                "void meshOptimizer::untangleBoundaryLayer()"
                "Found invalid faces in the boundary layer."
                " Cannot untangle mesh!!"
            );
        }
    }
    else
    {
        optimizeLowQualityFaces();
        removeUserConstraints();
        untangleMeshFV(2, 50, 1, -1, true);
    }
}

void meshOptimizer::optimizeLowQualityFaces(const label maxNumIterations)
{
    label nBadFaces, nIter(0);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
            lockedPoints.append(pointI);
    }

    do
    {
        labelHashSet lowQualityFaces;
        nBadFaces =
            polyMeshGenChecks::findLowQualityFaces
            (
                mesh_,
                lowQualityFaces,
                false,
                &changedFace
            );

        changedFace = false;
        forAllConstIter(labelHashSet, lowQualityFaces, it)
            changedFace[it.key()] = true;

        Info<< "Iteration " << nIter
            << ". Number of bad faces is " << nBadFaces << endl;

        //- perform optimisation
        if (nBadFaces == 0)
            break;

        partTetMesh tetMesh(mesh_, lockedPoints, lowQualityFaces, 2);

        //- construct tetMeshOptimisation and improve positions
        //- of points in the tet mesh
        tetMeshOptimisation tmo(tetMesh);

        tmo.optimiseUsingVolumeOptimizer();

        //- update points in the mesh from the new coordinates in the tet mesh
        tetMesh.updateOrigMesh(&changedFace);

    } while (++nIter < maxNumIterations);
}

void meshOptimizer::optimizeMeshNearBoundaries
(
    const label maxNumIterations,
    const label numLayersOfCells
)
{
    label nIter(0);

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
            lockedPoints.append(pointI);
    }

    partTetMesh tetMesh(mesh_, lockedPoints, numLayersOfCells);
    tetMeshOptimisation tmo(tetMesh);
    Info<< "Iteration:" << flush;
    do
    {
        tmo.optimiseUsingVolumeOptimizer(1);

        tetMesh.updateOrigMesh(&changedFace);

        Info<< "." << flush;

    } while (++nIter < maxNumIterations);

    Info<< endl;
}

void meshOptimizer::optimizeMeshFV
(
    const label numLaplaceIterations,
    const label maxNumGlobalIterations,
    const label maxNumIterations,
    const label maxNumSurfaceIterations,
    const label errorBufferLayers,
    const bool relaxedCheck,
    const bool relaxedBoundaryCheck,
    const bool checkWarped,
    const scalar minFaceArea
)
{
    Info<< "Starting smoothing the mesh" << endl;

    laplaceSmoother lps(mesh_, vertexLocation_);
    lps.optimizeLaplacianPC(numLaplaceIterations);

    untangleMeshFV
    (
        maxNumGlobalIterations,
        maxNumIterations,
        maxNumSurfaceIterations,
        errorBufferLayers,
        relaxedCheck,
        relaxedBoundaryCheck,
        checkWarped,
        minFaceArea
    );

    Info<< "Finished smoothing the mesh" << endl;
}

void meshOptimizer::optimizeMeshFVBestQuality
(
    const label maxNumIterations,
    const scalar threshold
)
{
    label nBadFaces, nIter(0);
    label minIter(-1);

    pointField minErrorPts = mesh_.points();

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
            lockedPoints.append(pointI);
    }

    label minNumBadFaces(labelMax);
    do
    {
        labelHashSet lowQualityFaces;
        nBadFaces =
            polyMeshGenChecks::findWorstQualityFaces
            (
                mesh_,
                lowQualityFaces,
                false,
                &changedFace,
                threshold
            );

        changedFace = false;
        forAllConstIter(labelHashSet, lowQualityFaces, it)
            changedFace[it.key()] = true;

        Info<< "Iteration " << nIter
            << ". Number of worst quality faces is " << nBadFaces << endl;

        //- perform optimisation
        if (nBadFaces == 0)
        {
            minNumBadFaces = 0;
            break;
        }

        if (nBadFaces < minNumBadFaces)
        {
            minErrorPts = mesh_.points();
            minNumBadFaces = nBadFaces;

            //- update the iteration number when the minimum is achieved
            minIter = nIter;
        }

        partTetMesh tetMesh(mesh_, lockedPoints, lowQualityFaces, 2);

        //- construct tetMeshOptimisation and improve positions
        //- of points in the tet mesh
        tetMeshOptimisation tmo(tetMesh);

        tmo.optimiseUsingVolumeOptimizer(20);

        //- update points in the mesh from the new coordinates in the tet mesh
        tetMesh.updateOrigMesh(&changedFace);

    } while ((nIter < minIter+3) && (++nIter < maxNumIterations));

    if (minNumBadFaces != 0)
    {
        mesh_.points() = minErrorPts;
    }
}

void meshOptimizer::optimizeMeshFVOrthogonality
(
    const label maxNumIterations,
    const scalar threshold
)
{
    label nBadFaces, nIter(0);
    label minIter(-1);

    pointField minErrorPts = mesh_.points();

    const faceListPMG& faces = mesh_.faces();
    boolList changedFace(faces.size(), true);

    //- check if any points in the tet mesh shall not move
    labelLongList lockedPoints;
    forAll(vertexLocation_, pointI)
    {
        if (vertexLocation_[pointI] & LOCKED)
            lockedPoints.append(pointI);
    }

    label minNumBadFaces(labelMax);
    do
    {
        labelHashSet lowQualityFaces;
        nBadFaces =
            polyMeshGenChecks::findHighOrthogonalityFaces
            (
                mesh_,
                lowQualityFaces,
                false,
                &changedFace,
                threshold
            );

        changedFace = false;
        forAllConstIter(labelHashSet, lowQualityFaces, it)
            changedFace[it.key()] = true;

        Info<< "Iteration " << nIter
            << ". Number of high orthogonlaity faces is " << nBadFaces << endl;

        //- perform optimisation
        if (nBadFaces == 0)
        {
            minNumBadFaces = 0;
            break;
        }

        if (nBadFaces < minNumBadFaces)
        {
            minErrorPts = mesh_.points();
            minNumBadFaces = nBadFaces;

            //- update the iteration number when the minimum is achieved
            minIter = nIter;
        }

        partTetMesh tetMesh(mesh_, lockedPoints, lowQualityFaces, 2);

        //- construct tetMeshOptimisation and improve positions
        //- of points in the tet mesh
        tetMeshOptimisation tmo(tetMesh);

        tmo.optimiseUsingVolumeOptimizer(20);

        //- update points in the mesh from the new coordinates in the tet mesh
        tetMesh.updateOrigMesh(&changedFace);

    } while ((nIter < minIter+3) && (++nIter < maxNumIterations));

    if (minNumBadFaces != 0)
    {
        mesh_.points() = minErrorPts;
    }
}


void meshOptimizer::lockErrorFreeCells
(
    const label nBufferLayers,
    const scalar minFaceArea
)
{
    const faceListPMG& faces = mesh_.faces();
    const cellListPMG& cells = mesh_.cells();
    const VRWGraph& pointCells = mesh_.addressingData().pointCells();
    boolList changedFace(faces.size(), true);

    label nCriticalErrors;
    labelHashSet badFaces;
    polyMeshGenChecks::findBadFaces
    (
        mesh_,
        true, //check warped faces
        minFaceArea,
        badFaces,
        nCriticalErrors,
        false,
        &changedFace
     );

    List<direction> useCell(cells.size(), direction(0));

    //- select cells containing at least one vertex of the bad faces
    forAll(faces, faceI)
    {
        if (badFaces.found(faceI))
        {
            const face& f = faces[faceI];

            forAll(f, pI)
            {
                forAllRow(pointCells, f[pI], pcI)
                    useCell[pointCells(f[pI], pcI)] = 1;
            }
        }
    }

    //- add additional layer of cells
    for (direction layerI=1;layerI<nBufferLayers;++layerI)
    {
        forAll(useCell, cI)
        {
            if (useCell[cI] == layerI)
            {
                const cell& c = cells[cI];

                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, pI)
                    {
                        forAllRow(pointCells, f[pI], pcI)
                        {
                            const label cLabel = pointCells(f[pI], pcI);
                            if (!useCell[cLabel])
                                useCell[cLabel] = layerI + 1;
                        }
                    }
                }
            }
        }

        if (Pstream::parRun())
        {
            const labelLongList& globalPointLabel =
                mesh_.addressingData().globalPointLabel();
            const VRWGraph& pProcs = mesh_.addressingData().pointAtProcs();
            const Map<label>& globalToLocal =
                mesh_.addressingData().globalToLocalPointAddressing();

            std::map<label, LongList<label>> eData;
            forAllConstIter(Map<label>, globalToLocal, iter)
            {
                const label pointI = iter();

                forAllRow(pProcs, pointI, procI)
                {
                    const label neiProc = pProcs(pointI, procI);
                    if (neiProc == Pstream::myProcNo())
                        continue;

                    if (eData.find(neiProc) == eData.end())
                    {
                        eData.insert
                        (
                            std::make_pair(neiProc, LongList<label>())
                        );
                    }

                    forAllRow(pointCells, pointI, pcI)
                    {
                        if (useCell[pointCells(pointI, pcI)] == layerI)
                        {
                            eData[neiProc].append(globalPointLabel[pointI]);
                            break;
                        }
                    }
                }
            }

            //- exchange data with other processors
            labelLongList receivedData;
            help::exchangeMap(eData, receivedData);

            forAll(receivedData, i)
            {
                const label pointI = globalToLocal[receivedData[i]];

                forAllRow(pointCells, pointI, pcI)
                {
                    const label cLabel = pointCells(pointI, pcI);
                    if (!useCell[cLabel])
                        useCell[cLabel] = layerI + 1;
                }
            }
        }
    }

    //Mark cells to be locked
    labelLongList lc;

    forAll(useCell, cI)
    {
        if (!useCell[cI])
        {
            lc.append(cI);
        }
    }
    lockCells(lc);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
