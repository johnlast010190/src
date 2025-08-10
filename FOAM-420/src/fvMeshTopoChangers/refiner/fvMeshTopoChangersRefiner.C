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
    (c) 2011-2014 OpenFOAM Foundation
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.
    (c) 2016-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "refiner/fvMeshTopoChangersRefiner.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "sets/topoSets/cellSet.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "fvMesh/fvPatches/constraint/cyclicAMI/cyclicAMIFvPatch.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "decompositionMethod/decompositionMethod.H"
#include "turbulenceModel.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshTopoChangers
{
    defineTypeNameAndDebug(refiner, 0);
    addToRunTimeSelectionTable(fvMeshTopoChanger, refiner, fvMesh);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshTopoChangers::refiner::findProtectedCells()
{
    if (protectedCell_.size() == 0 || protectedCell_.size() != mesh().nCells())
    {
        protectedCell_.clear();
        protectedCell_.setSize(mesh().nCells(), 0);
    }

    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    // Set cells that should not be refined.
    // This is currently any cell which does not have 8 anchor points or
    // uses any face which does not have 4 anchor points.
    // Note: do not use cellPoint addressing.

    // Count number of points <= cellLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelList nAnchors(mesh().nCells(), 0);

    label nProtected = 0;

    forAll(mesh().pointCells(), pointI)
    {
        const labelList& pCells = mesh().pointCells()[pointI];

        forAll(pCells, i)
        {
            label cellI = pCells[i];

            if (!protectedCell_.get(cellI))
            {
                if (pointLevel[pointI] <= cellLevel[cellI])
                {
                    nAnchors[cellI]++;

                    if (nAnchors[cellI] > 8)
                    {
                        protectedCell_.set(cellI, 1);
                        nProtected++;
                    }
                }
            }
        }
    }

    // Count number of points <= faceLevel
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Bit tricky since proc face might be one more refined than the owner
    // since the coupled one is refined.
    {
        labelList neiLevel(mesh().nFaces());

        for (label faceI = 0; faceI < mesh().nInternalFaces(); faceI++)
        {
            neiLevel[faceI] = cellLevel[mesh().faceNeighbour()[faceI]];
        }
        for
        (
            label faceI = mesh().nInternalFaces();
            faceI < mesh().nFaces();
            faceI++
        )
        {
            neiLevel[faceI] = cellLevel[mesh().faceOwner()[faceI]];
        }
        syncTools::swapFaceList(mesh(), neiLevel);

        boolList protectedFace(mesh().nFaces(), false);

        forAll(mesh().faceOwner(), faceI)
        {
            label faceLevel = max
            (
                cellLevel[mesh().faceOwner()[faceI]],
                neiLevel[faceI]
            );

            const face& f = mesh().faces()[faceI];

            label nAnchors = 0;

            forAll(f, fp)
            {
                if (pointLevel[f[fp]] <= faceLevel)
                {
                    nAnchors++;

                    if (nAnchors > 4)
                    {
                        protectedFace[faceI] = true;
                        break;
                    }
                }
            }
        }

        syncTools::syncFaceList(mesh(), protectedFace, orEqOp<bool>());

        for (label faceI = 0; faceI < mesh().nInternalFaces(); faceI++)
        {
            if (protectedFace[faceI])
            {
                protectedCell_.set(mesh().faceOwner()[faceI], 1);
                nProtected++;
                protectedCell_.set(mesh().faceNeighbour()[faceI], 1);
                nProtected++;
            }
        }
        for
        (
            label faceI = mesh().nInternalFaces();
            faceI < mesh().nFaces();
            faceI++
        )
        {
            if (protectedFace[faceI])
            {
                protectedCell_.set(mesh().faceOwner()[faceI], 1);
                nProtected++;
            }
        }

        // Also protect any cells that are less than hex
        forAll(mesh().cells(), cellI)
        {
            const cell& cFaces = mesh().cells()[cellI];

            if (cFaces.size() < 6)
            {
                if (protectedCell_.set(cellI, 1))
                {
                    nProtected++;
                }
            }
            else
            {
                forAll(cFaces, cFaceI)
                {
                    if (mesh().faces()[cFaces[cFaceI]].size() < 4)
                    {
                        if (protectedCell_.set(cellI, 1))
                        {
                            nProtected++;
                        }
                        break;
                    }
                }
            }
        }

        // Check cells for 8 corner points
        checkEightAnchorPoints(protectedCell_, nProtected);

        // Protect the wall layer cells
        if (protectLayers_)
        {
            // Count the number of protected wall cells
            label nWallCellsProtected = 0;

            const fvPatchList& patches = mesh().boundary();
            forAll(patches, patchi)
            {
                const fvPatch& currPatch = patches[patchi];
                if (isType<wallFvPatch>(currPatch))
                {
                    forAll(currPatch, facei)
                    {
                        label faceCelli = currPatch.faceCells()[facei];
                        if (!protectedCell_.get(faceCelli))
                        {
                            if (protectedCell_.set(faceCelli, 1))
                            {
                                nWallCellsProtected++;
                                nProtected++;
                            }
                        }
                    }
                }
            }

            label nWCP = returnReduce(nWallCellsProtected, sumOp<label>());

            if (gnWallCellsProtected_ != nWCP)
            {
                Info<< " Protected " << nWCP
                    << " wall cells from refinement" << endl;

                // Set new global count
                gnWallCellsProtected_ = nWCP;
            }
        }

        // Protect AMI patches automatically to prevent
        // AMI interpolation issues
        {
            // Count the number of protected AMI cells
            label nAMIPatchCellsProtected = 0;

            const fvPatchList& patches = mesh().boundary();
            forAll(patches, patchi)
            {
                const fvPatch& currPatch = patches[patchi];

                if (isType<cyclicAMIFvPatch>(currPatch))
                {
                    forAll(currPatch, facei)
                    {
                        label faceCelli = currPatch.faceCells()[facei];
                        if (!protectedCell_.get(faceCelli))
                        {
                            if (protectedCell_.set(faceCelli, 1))
                            {
                                nAMIPatchCellsProtected++;
                                nProtected++;
                            }
                        }
                    }
                }
            }

            label nAMIPCP =
                returnReduce(nAMIPatchCellsProtected, sumOp<label>());

            if (nAMIPCP)
            {
                Info<< " Protected "<< nAMIPCP
                     << " cells from refinement on AMI patches" << endl;
            }
        }

        if (protectedPatches_.size() > 0)
        {
            // Count the number of protected wall cells
            label nPatchCellsProtected = 0;

            const fvPatchList& patches = mesh().boundary();
            forAll(patches, patchi)
            {
                forAll(protectedPatches_,pI)
                {
                    const fvPatch& currPatch = patches[patchi];
                    if (currPatch.name() == protectedPatches_[pI])
                    {
                        forAll(currPatch, facei)
                        {
                            label faceCelli = currPatch.faceCells()[facei];
                            if (!protectedCell_.get(faceCelli))
                            {
                                if (protectedCell_.set(faceCelli, 1))
                                {
                                    nPatchCellsProtected++;
                                    nProtected++;
                                }
                            }
                        }

                    }
                }
            }

            label nPCP = returnReduce(nPatchCellsProtected, sumOp<label>());

            if (gnPatchCellsProtected_ != nPCP)
            {
                Info<< " Protected " << nPCP
                     << " cells from refinement on patches in"
                     << " protectPatches list " << protectedPatches_ << endl;

                // set new global count
                gnPatchCellsProtected_ = nPCP;
            }
        }

        // Protect cells in specific cellZones
        if (protectedCellZones_.size() > 0)
        {
            // Count the number of protected cell zone cells
            label nCellZoneCellsProtected = 0;

            forAll(protectedCellZones_, cZI)
            {
                const label cellZoneID = mesh().cellZones().
                    findZoneID(protectedCellZones_[cZI]);

                if (cellZoneID < 0)
                {
                    FatalError
                        << " Unable to find cell zone "
                        << protectedCellZones_[cZI] << endl
                        << " to protect cells in zone from refinement" << endl
                        << exit(FatalError);
                }

                const cellZone& cZone = mesh().cellZones()[cellZoneID];

                forAll(cZone, cI)
                {
                    if (!protectedCell_.get(cZone[cI]))
                    {
                        if (protectedCell_.set(cZone[cI], 1))
                        {
                            nCellZoneCellsProtected++;
                            nProtected++;
                        }
                    }
                }
            }

            label nCZCP =
                returnReduce(nCellZoneCellsProtected, sumOp<label>());

            if (gnCellZoneCellsProtected_ != nCZCP)
            {
                Info<< " Protected " << nCZCP
                     << " cells from refinement, defined in cellZones"
                     << " in protectCellZones list "
                     << protectedCellZones_ << endl;

                // Set new global count
                gnCellZoneCellsProtected_ = nCZCP;
            }
        }
    }

    // Add all cells from field protectedCells_
    if (mesh().foundObject<volScalarField>("protectedCells"))
    {
        forAll(protectedCells_(), cellI)
        {
            if (protectedCells_()[cellI] > 0)
            {
                if (protectedCell_.set(cellI, 1))
                {
                    nProtected++;
                }
            }
        }
    }

    if (returnReduce(nProtected, sumOp<label>()) == 0)
    {
        protectedCell_.clear();
    }
    else
    {
        cellSet protectedCells(mesh(), "protectedCells", nProtected);
        forAll(protectedCell_, cellI)
        {
            if (protectedCell_[cellI])
            {
                protectedCells.insert(cellI);
            }
        }

        Info<< " Detected " << returnReduce(nProtected, sumOp<label>())
            << " cells that are protected from refinement."
            << " Writing these to cellSet "
            << protectedCells.name()
            << "." << endl;

        if (mesh().time().write())
        {
            protectedCells.write();
        }
    }
}


Foam::label Foam::fvMeshTopoChangers::refiner::parentID(label p)
{
    label pParent = meshCutter().history().splitCells()[p].parent_;

    if (pParent < 0)
    {
        return p;
    }
    else
    {
        return parentID(pParent);
    }
}


bool Foam::fvMeshTopoChangers::refiner::walkNeighbours
(
    const label& patchMaster,
    const label& patchSlave,
    labelList& agglom,
    point& aveCC,
    label& nAMICells
)
{
    const labelList& owner = mesh().faceOwner();

    label nbr = -1;

    {
        const polyPatch& pp = mesh().boundaryMesh()[patchMaster];
        forAll(pp, i)
        {
            label own = owner[pp.start() + i];
            if (agglom[own] != -1)
            {
                nbr = agglom[own];
                break;
            }
        }
    }

    if (nbr == -1)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchSlave];
        forAll(pp, i)
        {
            label own = owner[pp.start() + i];
            if (agglom[own] != -1)
            {
                nbr = agglom[own];
                break;
            }
        }
    }

    bool reset = false;

    if (nbr != -1)
    {
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchMaster];
            forAll(pp, i)
            {
                label own = owner[pp.start() + i];
                if (agglom[own] == -1)
                {
                    agglom[own] = nbr;
                    aveCC += mesh().cellCentres()[own];
                    nAMICells++;
                    reset = true;
                }
            }
        }
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchSlave];
            forAll(pp, i)
            {
                label own = owner[pp.start() + i];
                if (agglom[own] == -1)
                {
                    agglom[own] = nbr;
                    aveCC += mesh().cellCentres()[own];
                    nAMICells++;
                    reset = true;
                }
            }
        }
    }

    return reset;
}


Foam::label Foam::fvMeshTopoChangers::refiner::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val)
        {
            n++;
        }
    }

    return n;
}


void Foam::fvMeshTopoChangers::refiner::calculateProtectedCells
(
    PackedBoolList& unrefineableCell
) const
{
    bool zeroProtected(protectedCell_.empty());
    if (returnReduce(zeroProtected, andOp<bool>()))
    {
        unrefineableCell.clear();
        return;
    }

    const labelList& cellLevel = meshCutter_.cellLevel();

    label nInvalidAnchors = 0;

    unrefineableCell = protectedCell_;

    // Check cells for 8 corner points
    checkEightAnchorPoints(unrefineableCell, nInvalidAnchors);
    Info<< "Detected " << returnReduce(nInvalidAnchors, sumOp<label>())
        << " cells that are protected due to invalid anchors." << endl;

    // Get neighbouring cell level
    labelList neiLevel(mesh().nFaces() - mesh().nInternalFaces());

    for
    (
        label facei = mesh().nInternalFaces();
        facei < mesh().nFaces();
        facei++
    )
    {
        neiLevel[facei - mesh().nInternalFaces()] =
            cellLevel[mesh().faceOwner()[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh(), neiLevel);

    while (true)
    {
        // Pick up faces on border of protected cells
        boolList seedFace(mesh().nFaces(), false);

        forAll(mesh().faceNeighbour(), facei)
        {
            label own = mesh().faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            label nei = mesh().faceNeighbour()[facei];
            bool neiProtected = unrefineableCell.get(nei);

            if (ownProtected && (cellLevel[nei] > cellLevel[own]))
            {
                seedFace[facei] = true;
            }
            else if (neiProtected && (cellLevel[own] > cellLevel[nei]))
            {
                seedFace[facei] = true;
            }
        }
        for
        (
            label facei = mesh().nInternalFaces();
            facei < mesh().nFaces();
            facei++
        )
        {
            label own = mesh().faceOwner()[facei];
            bool ownProtected = unrefineableCell.get(own);
            if
            (
                ownProtected
             && (neiLevel[facei - mesh().nInternalFaces()] > cellLevel[own])
            )
            {
                seedFace[facei] = true;
            }
        }

        syncTools::syncFaceList(mesh(), seedFace, orEqOp<bool>());

        // Extend unrefineableCell
        bool hasExtended = false;

        for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
        {
            if (seedFace[facei])
            {
                label own = mesh().faceOwner()[facei];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }

                label nei = mesh().faceNeighbour()[facei];
                if (unrefineableCell.get(nei) == 0)
                {
                    unrefineableCell.set(nei, 1);
                    hasExtended = true;
                }
            }
        }
        for
        (
            label facei = mesh().nInternalFaces();
            facei < mesh().nFaces();
            facei++
        )
        {
            if (seedFace[facei])
            {
                label own = mesh().faceOwner()[facei];
                if (unrefineableCell.get(own) == 0)
                {
                    unrefineableCell.set(own, 1);
                    hasExtended = true;
                }
            }
        }

        if (!returnReduce(hasExtended, orOp<bool>()))
        {
            break;
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::checkRefHistory
(
    PackedBoolList& hasRefHistory
) const
{
    // visibleCells[i] == -1 is a cell that has not been refined.
    // Only allow cells that have been refined to be unrefined.
    const labelList& visibleCells =
        meshCutter().history().visibleCells();

    forAll(hasRefHistory, cellI)
    {
        if (visibleCells[cellI] > -1)
        {
            hasRefHistory.set(cellI, true);
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::readDict()
{
    const dictionary& refineDict(dict());

    List<Pair<word>> fluxVelocities = List<Pair<word>>
    (
        refineDict.lookup("correctFluxes")
    );

    // Rework into hashtable
    correctFluxes_.resize(fluxVelocities.size());
    forAll(fluxVelocities, i)
    {
        correctFluxes_.insert(fluxVelocities[i][0], fluxVelocities[i][1]);
    }

    dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));

    protectLayers_ = refineDict.lookupOrDefault<Switch>("protectLayers", true);

    correctY_ = refineDict.lookupOrDefault<Switch>("correctY", true);

    if (refineDict.found("protectedPatches"))
        protectedPatches_ = wordList(refineDict.lookup("protectedPatches"));

    if (refineDict.found("protectedCellZones"))
        protectedCellZones_ =
            wordList(refineDict.lookup("protectedCellZones"));

    if (refineDict.found("refineStartTime"))
        refineStartTime_ = readScalar(refineDict.lookup("refineStartTime"));

    if (refineDict.found("refineStopTime"))
        refineStopTime_ = readScalar(refineDict.lookup("refineStopTime"));

    danglingRefine_ =
        refineDict.lookupOrDefault<Switch>("danglingRefine", false);

    if (refineDict.found("mapSurfaceFields"))
    {
        List<word> surfFlds = List<word>
        (
            refineDict.lookup("mapSurfaceFields")
        );

        // Rework into hashtable
        mapSurfaceFields_.resize(surfFlds.size());
        forAll(surfFlds, i)
        {
            mapSurfaceFields_.insert(surfFlds[i], surfFlds[i]);
        };
    }
}


// Refines cells, maps fields and recalculates (an approximate) flux
Foam::autoPtr<Foam::mapPolyMesh>
Foam::fvMeshTopoChangers::refiner::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine
    polyTopoChange meshMod(mesh());

    // Play refinement commands into mesh changer
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh(), true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh(), false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh().globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map
        for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
        {
            label oldFacei = map().faceMap()[facei];

            if (oldFacei >= mesh().nInternalFaces())
            {
                FatalErrorInFunction
                    << "New internal face:" << facei
                    << " fc:" << mesh().faceCentres()[facei]
                    << " originates from boundary oldFace:" << oldFacei
                    << abort(FatalError);
            }
        }
    }

    //    // Remove the stored tet base points
    //    tetBasePtIsPtr_.clear();
    //    // Remove the cell tree
    //    cellTreePtr_.clear();

    // Update fields
    mesh().updateMesh(map);

    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces
        // on the coarse cell that get split into four (or rather the
        // master face gets modified and three faces get added from the master)
        labelHashSet masterFaces(4*cellsToRefine.size());

        forAll(faceMap, facei)
        {
            label oldFacei = faceMap[facei];

            if (oldFacei >= 0)
            {
                label masterFacei = reverseFaceMap[oldFacei];

                if (masterFacei < 0)
                {
                    FatalErrorInFunction
                        << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << facei
                        << abort(FatalError);
                }
                else if (masterFacei != facei)
                {
                    masterFaces.insert(masterFacei);
                }
            }
        }
        if (debug)
        {
            Pout<< "Found " << masterFaces.size() << " split faces " << endl;
        }

        HashTable<surfaceScalarField*> fluxes
        (
            mesh().lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                if (debug)
                {
                    WarningInFunction
                        << "Cannot find surfaceScalarField " << iter.key()
                        << " in user-provided flux mapping table "
                        << correctFluxes_ << endl
                        << "    The flux mapping table is used to recreate the"
                        << " flux on newly created faces." << endl
                        << "    Either add the entry if it is a flux or use ("
                        << iter.key() << " none) to suppress this warning."
                        << endl;
                }
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (debug)
            {
                Pout<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    mesh().lookupObject<volVectorField>(UName)
                )
              & mesh().Sf()
            );

            // Recalculate new internal faces
            for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
            {
                label oldFacei = faceMap[facei];

                if (oldFacei == -1)
                {
                    // Inflated/appended
                    phi[facei] = phiU[facei];
                }
                else if (reverseFaceMap[oldFacei] != facei)
                {
                    // Face-from-masterface
                    phi[facei] = phiU[facei];
                }
            }

            // Recalculate new boundary faces
            surfaceScalarField::Boundary& phiBf =
                phi.boundaryFieldRef();
            forAll(phiBf, patchi)
            {
                fvsPatchScalarField& patchPhi = phiBf[patchi];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchi];

                label facei = patchPhi.patch().start();

                forAll(patchPhi, i)
                {
                    label oldFacei = faceMap[facei];

                    if (oldFacei == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFacei] != facei)
                    {
                        // Face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    facei++;
                }
            }

            // Update master faces
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label facei = iter.key();

                if (mesh().isInternalFace(facei))
                {
                    phi[facei] = phiU[facei];
                }
                else
                {
                    label patchi = mesh().boundaryMesh().whichPatch(facei);
                    label i = facei - mesh().boundaryMesh()[patchi].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchi];

                    fvsPatchScalarField& patchPhi = phiBf[patchi];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }

    // Map to newly generated (non-flux) surfaceFields
    const labelList& faceMap = map().faceMap();
    mapNewInternalFaces<scalar>(faceMap);
    mapNewInternalFaces<vector>(faceMap);

    // Update numbering of cells/vertices
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(mesh().nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            newProtectedCell.set(celli, protectedCell_.get(oldCelli));
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh>
Foam::fvMeshTopoChangers::refiner::unrefine
(
    const labelList& splitPoints
)
{
    polyTopoChange meshMod(mesh());

    // Play refinement commands into mesh changer
    meshCutter_.setUnrefinement(splitPoints, meshMod);

    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // For each face resulting of split of face into four, store the midpoint.
    Map<label> faceToSplitPoint(3*splitPoints.size());

    {
        forAll(splitPoints, i)
        {
            label pointi = splitPoints[i];

            const labelList& pEdges = mesh().pointEdges()[pointi];

            forAll(pEdges, j)
            {
                label otherPointi =
                    mesh().edges()[pEdges[j]].otherVertex(pointi);

                const labelList& pFaces = mesh().pointFaces()[otherPointi];

                forAll(pFaces, pFacei)
                {
                    faceToSplitPoint.insert(pFaces[pFacei], otherPointi);
                }
            }
        }
    }

    // Change mesh and generate map
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh(), true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh(), false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << mesh().globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    mesh().updateMesh(map);

    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified faces
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        HashTable<surfaceScalarField*> fluxes
        (
            mesh().lookupClass<surfaceScalarField>()
        );
        forAllIter(HashTable<surfaceScalarField*>, fluxes, iter)
        {
            if (!correctFluxes_.found(iter.key()))
            {
                if (debug)
                {
                    WarningInFunction
                        << "Cannot find surfaceScalarField " << iter.key()
                        << " in user-provided flux mapping table "
                        << correctFluxes_ << endl
                        << "    The flux mapping table is used to recreate the"
                        << " flux on newly created faces." << endl
                        << "    Either add the entry if it is a flux or use ("
                        << iter.key() << " none) to suppress this warning."
                        << endl;
                }
                continue;
            }

            const word& UName = correctFluxes_[iter.key()];

            if (UName == "none")
            {
                continue;
            }

            if (debug)
            {
                Info<< "Mapping flux " << iter.key()
                    << " using interpolated flux " << UName
                    << endl;
            }

            surfaceScalarField& phi = *iter();
            surfaceScalarField::Boundary& phiBf =
                phi.boundaryFieldRef();

            const surfaceScalarField phiU
            (
                fvc::interpolate
                (
                    mesh().lookupObject<volVectorField>(UName)
                )
              & mesh().Sf()
            );

            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                label oldFacei = iter.key();
                label oldPointi = iter();

                if (reversePointMap[oldPointi] < 0)
                {
                    // Midpoint was removed. See if face still exists.
                    label facei = reverseFaceMap[oldFacei];

                    if (facei >= 0)
                    {
                        if (mesh().isInternalFace(facei))
                        {
                            phi[facei] = phiU[facei];
                        }
                        else
                        {
                            label patchi =
                                mesh().boundaryMesh().whichPatch(facei);
                            label i =
                                facei - mesh().boundaryMesh()[patchi].start();

                            const fvsPatchScalarField& patchPhiU =
                                phiU.boundaryField()[patchi];
                            fvsPatchScalarField& patchPhi = phiBf[patchi];
                            patchPhi[i] = patchPhiU[i];
                        }
                    }
                }
            }
        }
    }

    // Update numbering of cells/vertices
    meshCutter_.updateMesh(map);

    // Update numbering of protectedCell_
    if (protectedCell_.size())
    {
        PackedBoolList newProtectedCell(mesh().nCells());

        forAll(newProtectedCell, celli)
        {
            label oldCelli = map().cellMap()[celli];
            if (oldCelli >= 0)
            {
                newProtectedCell.set(celli, protectedCell_.get(oldCelli));
            }
        }
        protectedCell_.transfer(newProtectedCell);
    }

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


Foam::scalarField
Foam::fvMeshTopoChangers::refiner::maxPointField(const scalarField& pFld) const
{
    scalarField vFld(mesh().nCells(), -GREAT);

    forAll(mesh().pointCells(), pointi)
    {
        const labelList& pCells = mesh().pointCells()[pointi];

        forAll(pCells, i)
        {
            vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointi]);
        }
    }

    return vFld;
}


Foam::scalarField
Foam::fvMeshTopoChangers::refiner::minCellField
(
    const volScalarField& vFld
) const
{
    scalarField pFld(mesh().nPoints(), GREAT);

    forAll(mesh().pointCells(), pointi)
    {
        const labelList& pCells = mesh().pointCells()[pointi];

        forAll(pCells, i)
        {
            pFld[pointi] = min(pFld[pointi], vFld[pCells[i]]);
        }
    }

    return pFld;
}


Foam::scalarField
Foam::fvMeshTopoChangers::refiner::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(mesh().nPoints());

    forAll(mesh().pointCells(), pointi)
    {
        const labelList& pCells = mesh().pointCells()[pointi];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointi] = sum/pCells.size();
    }

    return pFld;
}


Foam::scalarField Foam::fvMeshTopoChangers::refiner::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    scalarField c(fld.size(), -1);

    forAll(fld, i)
    {
        scalar err = min(fld[i] - minLevel, maxLevel - fld[i]);

        if (err >= 0)
        {
            c[i] = err;
        }
    }

    return c;
}


void Foam::fvMeshTopoChangers::refiner::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    scalarField cellError
    (
        maxPointField
        (
            error
            (
                cellToPoint(vFld),
                lowerRefineLevel,
                upperRefineLevel
            )
        )
    );

    // Mark cells that are candidates for refinement
    forAll(cellError, celli)
    {
        if (cellError[celli] > 0)
        {
            candidateCell.set(celli, 1);
        }
    }
}


Foam::labelList Foam::fvMeshTopoChangers::refiner::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const PackedBoolList& candidateCell
) const
{
    // Every refined cell causes 7 extra cells
    label nTotToRefine = (maxCells - mesh().globalData().nTotalCells())/7;

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Mark cells that cannot be refined since they would trigger refinement
    // of protected cells (since 2:1 cascade)
    PackedBoolList unrefineableCell;
    calculateProtectedCells(unrefineableCell);

    // Count current selection
    label nCandidates = returnReduce(count(candidateCell, 1), sumOp<label>());

    // Collect all cells
    DynamicList<label> candidates(mesh().nCells());

    if (nCandidates < nTotToRefine)
    {
        forAll(candidateCell, celli)
        {
            if
            (
                cellLevel[celli] < maxRefinement
             && candidateCell.get(celli)
             && (
                    unrefineableCell.empty()
                 || !unrefineableCell.get(celli)
                )
            )
            {
                candidates.append(celli);
            }
        }
    }
    else
    {
        // Sort by error? For now just truncate.
        for (label level = 0; level < maxRefinement; level++)
        {
            forAll(candidateCell, celli)
            {
                if
                (
                    cellLevel[celli] == level
                 && candidateCell.get(celli)
                 && (
                        unrefineableCell.empty()
                     || !unrefineableCell.get(celli)
                    )
                )
                {
                    candidates.append(celli);
                }
            }

            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }
    }

    // Guarantee 2:1 refinement after refinement
    labelList consistentSet
    (
        meshCutter_.consistentRefinement
        (
            candidates.shrink(),
            true               // Add to set to guarantee 2:1
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << mesh().globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


Foam::labelList Foam::fvMeshTopoChangers::refiner::selectUnrefinePoints
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // Find if cell has been refined before
    PackedBoolList hasRefHistory(mesh().nCells(), false);
    checkRefHistory(hasRefHistory);

    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getSplitPoints());

    DynamicList<label> newSplitPoints(splitPoints.size());

    forAll(splitPoints, i)
    {
        label pointi = splitPoints[i];

        if (pFld[pointi] < unrefineLevel)
        {
            // Check that all cells are not marked
            const labelList& pCells = mesh().pointCells()[pointi];

            bool hasMarked = false;

            bool hasHistory = false;

            forAll(pCells, pCelli)
            {
                if (hasRefHistory.get(pCells[pCelli]))
                {
                    // Could be unrefined
                    hasHistory = true;
                }

                if (markedCell.get(pCells[pCelli]))
                {
                    // Can't be unrefined
                    hasMarked = true;
                    break;
                }
            }

            if (!hasMarked && hasHistory)
            {
                newSplitPoints.append(pointi);
            }
        }
    }

    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement
        (
            newSplitPoints,
            false
        )
    );
    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}


void Foam::fvMeshTopoChangers::refiner::extendMarkedCells
(
    PackedBoolList& markedCell
) const
{
    // Mark faces using any marked cell
    boolList markedFace(mesh().nFaces(), false);

    forAll(markedCell, celli)
    {
        if (markedCell.get(celli))
        {
            const cell& cFaces = mesh().cells()[celli];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(mesh(), markedFace, orEqOp<bool>());

    // Update cells using any markedFace
    for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
    {
        if (markedFace[facei])
        {
            markedCell.set(mesh().faceOwner()[facei], 1);
            markedCell.set(mesh().faceNeighbour()[facei], 1);
        }
    }
    for
    (
        label facei = mesh().nInternalFaces();
        facei < mesh().nFaces();
        facei++
    )
    {
        if (markedFace[facei])
        {
            markedCell.set(mesh().faceOwner()[facei], 1);
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::checkEightAnchorPoints
(
    PackedBoolList& protectedCell,
    label& nProtected
) const
{
    const labelList& cellLevel = meshCutter_.cellLevel();
    const labelList& pointLevel = meshCutter_.pointLevel();

    labelList nAnchorPoints(mesh().nCells(), 0);

    forAll(pointLevel, pointi)
    {
        const labelList& pCells = mesh().pointCells(pointi);

        forAll(pCells, pCelli)
        {
            label celli = pCells[pCelli];

            if (pointLevel[pointi] <= cellLevel[celli])
            {
                // Check if cell has already 8 anchor points -> protect cell
                if (nAnchorPoints[celli] == 8)
                {
                    if (protectedCell.set(celli, true))
                    {
                        nProtected++;
                    }
                }

                if (!protectedCell[celli])
                {
                    nAnchorPoints[celli]++;
                }
            }
        }
    }

    forAll(protectedCell, celli)
    {
        if (!protectedCell[celli] && nAnchorPoints[celli] != 8)
        {
            protectedCell.set(celli, true);
            nProtected++;
        }
    }
}


template<class T>
void Foam::fvMeshTopoChangers::refiner::mapNewInternalFaces
(
    const labelList& faceMap
)
{
    typedef GeometricField<T, fvsPatchField, surfaceMesh> GeoField;
    HashTable<GeoField*> sFlds(mesh().objectRegistry::lookupClass<GeoField>());

    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();
    const dimensionedScalar deltaN = 1e-8/pow(average(mesh().V()), 1.0/3.0);

    forAllIter(typename HashTable<GeoField*>, sFlds, iter)
    {
        GeoField& sFld = *iter();
        if (mapSurfaceFields_.found(iter.key()))
        {
            if (debug)
            {
                Info<< "Correct new internal faces for field "
                    << sFld.name() << endl;
            }

            // Copy into Field<Type> to prevent looping over
            // boundary patches in averaging algorithm
            Field<T> tsFld(mesh().nFaces(), pTraits<T>::zero);
            forAll(sFld.internalField(), iFace)
            {
                tsFld[iFace] = sFld.internalField()[iFace];
            }
            forAll(sFld.boundaryField(), iPatch)
            {
                label globalIdx = mesh().boundaryMesh()[iPatch].start();
                forAll(sFld.boundaryField()[iPatch], faceI)
                {
                    tsFld[faceI+globalIdx] =
                        sFld.boundaryField()[iPatch][faceI];
                }
            }

            // Loop over all faces
            for (label facei = 0; facei < mesh().nInternalFaces(); facei++)
            {
                label oldFacei = faceMap[facei];

                // Map surface field on newly generated faces
                if (oldFacei == -1)
                {
                    // Find owner and neighbour cell
                    cell faceOwner = mesh().cells()[owner[facei]];
                    cell faceNeighbour = mesh().cells()[neighbour[facei]];

                    // Loop over all owner/neighbour cell faces
                    // and find already mapped ones
                    T tmpValue = pTraits<T>::zero;
                    scalar magFld = 0;
                    label counter = 0;

                    // Simple averaging of all neighbour master-faces
                    forAll(faceOwner, iFace)
                    {
                        if (faceMap[faceOwner[iFace]] != -1)
                        {
                            tmpValue += tsFld[faceOwner[iFace]];
                            magFld += mag(tsFld[faceOwner[iFace]]);
                            counter++;
                        }
                    }
                    forAll(faceNeighbour, iFace)
                    {
                        if (faceMap[faceNeighbour[iFace]] != -1)
                        {
                            tmpValue += tsFld[faceNeighbour[iFace]];
                            magFld += mag(tsFld[faceNeighbour[iFace]]);
                            counter++;
                        }
                    }
                    if (counter > 0)
                    {
                        if
                        (
                            GeometricField<T, fvsPatchField, surfaceMesh>::
                            typeName == "surfaceScalarField"
                        )
                        {
                            tmpValue /= counter;
                        }
                        else if
                        (
                            GeometricField<T, fvsPatchField, surfaceMesh>::
                            typeName == "surfaceVectorField"
                        )
                        {
                            magFld /= counter;
                            tmpValue *=
                                magFld/(mag(tmpValue) + deltaN.value());
                        }
                        else
                        {
                            FatalErrorInFunction
                                << "Mapping implementation only valid for"
                                << " scalar and vector fields! \n Field "
                                << sFld.name() << " is of type: "
                                << GeometricField
                                    <
                                        T,
                                        fvsPatchField,
                                        surfaceMesh
                                    >::typeName
                                << abort(FatalError);
                        }
                    }

                    sFld[facei] = tmpValue;
                }
            }
        }
    }
}


void Foam::fvMeshTopoChangers::refiner::danglingCellRefine
(
    labelList& cellsToRefine,
    const label nFaces
)
{
    // Determine cells to refine
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    const cellList& cells = mesh().cells();
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    labelList candidateCells;
    {
        cellSet candidateCellSet(mesh(), "candidateCells", cells.size()/1000);

        forAll(cells, cellI)
        {
            const cell& cFaces = cells[cellI];

            label nIntFaces = 0;
            forAll(cFaces, i)
            {
                label bFaceI = cFaces[i] - mesh().nInternalFaces();
                if (bFaceI < 0)
                {
                    nIntFaces++;
                }
                else
                {
                    label patchI = pbm.patchID()[bFaceI];
                    if (pbm[patchI].coupled())
                    {
                        nIntFaces++;
                    }
                }
            }

            if (nIntFaces == nFaces)
            {
                candidateCellSet.insert(cellI);
            }
        }

        candidateCells = candidateCellSet.toc();
    }

    cellsToRefine.append(candidateCells);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::refiner::refiner(fvMesh& mesh)
:
    fvMeshTopoChanger(mesh),
    meshCutter_(mesh, true, true),
    dumpLevel_(false),
    nRefinementIterations_(0),
    protectedCell_(mesh.nCells(), 0),
    protectLayers_(true),
    correctY_(true),
    protectedPatches_(0),
    protectedCellZones_(0),
    gnProtected_(0),
    gnWallCellsProtected_(0),
    gnPatchCellsProtected_(0),
    gnCellZoneCellsProtected_(0),
    refineStartTime_(-1),
    refineStopTime_(VGREAT),
    danglingRefine_(false),
    protectedCells_(),
    refineMeshOuterCorr_(true),
    isFirstIter_(true)
{
    // Read static part of dictionary
    readDict();

    if (refineStartTime_ > -1)
        Info<< " Refinement will start at time " << refineStartTime_ << endl;

    if (refineStopTime_ < VGREAT)
        Info<< " Refinement will stop after time " << refineStopTime_ << endl;

    IOobject protectedCellsHeader
    (
        "protectedCells",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (protectedCellsHeader.typeHeaderOk<volScalarField>(true))
    {
        protectedCells_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "protectedCells",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        // Alternatively check if file exists in constant
        IOobject protectedCellsHeaderConstant
        (
            "protectedCells",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        );

        if (protectedCellsHeaderConstant.typeHeaderOk<volScalarField>(true))
        {
            protectedCells_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        "protectedCells",
                        mesh.time().constant(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::AUTO_WRITE
                     ),
                    mesh
                 )
             );
        }
    }

    // Find the protectedCells_
    findProtectedCells();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshTopoChangers::refiner::~refiner()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fvMeshTopoChangers::refiner::update()
{
    // Re-read dictionary. Usually small so takes trivial amount of time
    // compared to actual refinement. Also very useful to be able to modify
    // on-the-fly.
    const dictionary& refineDict(dict());

    label refineInterval = readLabel(refineDict.lookup("refineInterval"));
    label rebalanceInterval =
        refineDict.lookupOrDefault<label>("rebalanceInterval", label(1));

    bool hasChanged = false;

    if (refineInterval == 0)
    {
        mesh().topoChanging(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorInFunction
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }

    if (rebalanceInterval < 1)
    {
        FatalErrorInFunction
            << "Illegal rebalanceInterval " << rebalanceInterval << nl
            << "The rebalanceInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }

    // Refine during refinement times.
    // Default start time is 0 and stop time is VGREAT.
    // Re-read if the start and stop times have changed.
    if (refineDict.found("refineStartTime"))
        refineStartTime_ = readScalar(refineDict.lookup("refineStartTime"));

    if (refineDict.found("refineStopTime"))
        refineStopTime_ = readScalar(refineDict.lookup("refineStopTime"));

    fvMesh& mesh(this->mesh());

    // Note: if PIMPLE first iteration we can perform AMR once, by setting
    //       refineMeshOuterCorrectors to false in PIMPLE dictionary.
    //       Otherwise, the use of AMR is controlled by the 'refineInterval'
    //       parameter in the dynamicMeshDict
    if (mesh.foundObject<pimpleControl>("solutionControl"))
    {
        const pimpleControl& pimple =
            mesh.lookupObject<pimpleControl>("solutionControl");
        isFirstIter_ = pimple.firstIter();
        refineMeshOuterCorr_ =
            pimple.dict().lookupOrDefault<Switch>
            (
                "refineMeshOuterCorrectors",
                true
            );
    }

    if
    (
        mesh.time().value() >= refineStartTime_
     && mesh.time().value() <= refineStopTime_
    )
    {
        // Note: cannot refine at time 0 since no V0 present,
        // since mesh not moved yet
        /*
        if
        (
            mesh.time().timeIndex() > 0
         && mesh.time().timeIndex() % refineInterval == 0
        )
        */
        if
        (
            mesh.time().timeIndex() > 0
         && mesh.time().timeIndex() % refineInterval == 0
         && (isFirstIter_ || refineMeshOuterCorr_)
        )
        {
            label maxCells = readLabel(refineDict.lookup("maxCells"));

            if (maxCells <= 0)
            {
                FatalErrorInFunction
                    << "Illegal maximum number of cells " << maxCells << nl
                    << "The maxCells setting in the dynamicMeshDict should"
                    << " be > 0." << nl
                    << exit(FatalError);
            }

            label maxRefinement =
                readLabel(refineDict.lookup("maxRefinement"));

            if (maxRefinement <= 0)
            {
                FatalErrorInFunction
                    << "Illegal maximum refinement level "
                    << maxRefinement << nl
                    << "The maxCells setting in the dynamicMeshDict should"
                    << " be > 0." << nl
                    << exit(FatalError);
            }

            const word fieldName(refineDict.lookup("field"));

            const volScalarField& vFld =
                mesh.lookupObject<volScalarField>(fieldName);

            const scalar lowerRefineLevel =
                readScalar(refineDict.lookup("lowerRefineLevel"));
            const scalar upperRefineLevel =
                readScalar(refineDict.lookup("upperRefineLevel"));
            const scalar unrefineLevel =
                readScalar(refineDict.lookup("unrefineLevel"));
            const label nBufferLayers =
                readLabel(refineDict.lookup("nBufferLayers"));
            const scalar maxLoadImbalance = refineDict.
                lookupOrDefault<scalar>("maxLoadImbalance", 0);

            // Cells marked for refinement or otherwise protected
            // from unrefinement
            PackedBoolList refineCell(mesh.nCells());

            if (mesh.globalData().nTotalCells() < maxCells)
            {
                // Determine candidates for refinement (looking at field only)
                selectRefineCandidates
                (
                    lowerRefineLevel,
                    upperRefineLevel,
                    vFld,
                    refineCell
                );

                // Select subset of candidates. Take into account max allowable
                // cells, refinement level and protected cells.
                labelList cellsToRefine
                (
                    selectRefineCells
                    (
                        maxCells,
                        maxRefinement,
                        refineCell
                    )
                );

                if (danglingRefine_)
                {
                    // Refine cells with almost all sides refined.
                    // Alter marked cells analogous to snappeRefineDriver.
                    // Refine any hexes with 5 or 6 faces refined to make
                    // smooth edges.
                    danglingCellRefine
                    (
                        cellsToRefine,
                        18      // 1 coarse face + 5 refined faces
                    );
                    danglingCellRefine
                    (
                        cellsToRefine,
                        21     // 1 coarse face + 5 refined faces
                    );
                    danglingCellRefine
                    (
                        cellsToRefine,
                        24     // 0 coarse faces + 6 refined faces
                    );
                }

                label nCellsToRefine = returnReduce
                (
                    cellsToRefine.size(), sumOp<label>()
                );

                if (nCellsToRefine > 0)
                {
                    // Refine/update mesh and map fields
                    autoPtr<mapPolyMesh> map = refine(cellsToRefine);

                    // Update refineCell. Note that some of the marked ones
                    // have not been refined due to constraints.
                    {
                        const labelList& cellMap = map().cellMap();
                        const labelList& reverseCellMap =
                            map().reverseCellMap();

                        PackedBoolList newRefineCell(cellMap.size());

                        forAll(cellMap, celli)
                        {
                            label oldCelli = cellMap[celli];

                            if (oldCelli < 0)
                            {
                                newRefineCell.set(celli, 1);
                            }
                            else if (reverseCellMap[oldCelli] != celli)
                            {
                                newRefineCell.set(celli, 1);
                            }
                            else
                            {
                                newRefineCell.set
                                (
                                    celli,
                                    refineCell.get(oldCelli)
                                );
                            }
                        }
                        refineCell.transfer(newRefineCell);
                    }

                    // Extend with a buffer layer to prevent neighbouring
                    // points being unrefined
                    for (label i = 0; i < nBufferLayers; i++)
                    {
                        extendMarkedCells(refineCell);
                    }

                    hasChanged = true;
                }
            }

            {
                // Select unrefineable points that are not marked in refineCell
                labelList pointsToUnrefine
                (
                    selectUnrefinePoints
                    (
                        unrefineLevel,
                        refineCell,
                        minCellField(vFld)
                    )
                );

                label nSplitPoints = returnReduce
                (
                    pointsToUnrefine.size(),
                    sumOp<label>()
                );

                if (nSplitPoints > 0)
                {
                    // Refine/update mesh
                    unrefine(pointsToUnrefine);

                    hasChanged = true;
                }
            }

            if ((nRefinementIterations_ % 10) == 0)
            {
                // Compact refinement history occassionally (how often?).
                // Unrefinement causes holes in the refinementHistory.
                const_cast<refinementHistory&>
                (
                    meshCutter().history()
                ).compact();
            }

            bool loadImbalanced = false;

            // Load balance after refine/unrefine
            if ((maxLoadImbalance > 0) && hasChanged)
            {
                if (Pstream::nProcs() > 1)
                {
                    scalar nIdealCells =
                        mesh.globalData().nTotalCells()/Pstream::nProcs();

                    scalar imbalance = returnReduce
                    (
                        scalar(mag(1.0 - mesh.nCells()/nIdealCells)),
                        maxOp<scalar>()
                    );

                    if (imbalance <= maxLoadImbalance)
                    {
                        Info<< "Skipping balancing since max imbalance "
                            << imbalance << " is less than allowable "
                            << maxLoadImbalance << endl;
                    }
                    else if
                    (
                        imbalance > maxLoadImbalance
                     && (nRefinementIterations_ % rebalanceInterval) != 0
                    )
                    {
                        Info<< "Skipping balancing since this is not a rebalance interval \n"
                            << " max imbalance " << imbalance
                            << " max allowable imbalance " << maxLoadImbalance << "\n"
                            << " current interval "<< nRefinementIterations_ % rebalanceInterval
                            << " of "<< rebalanceInterval <<" before rebalance."
                            << endl;
                    }
                    else
                    {
                        Info<< "Balancing mesh since max imbalance "
                            << imbalance << " is greater than allowable "
                            << maxLoadImbalance << endl;

                        // Create decomposeParDict for dynamic rebalancing
                        dictionary decompDict;

                        // By default if method absent, use ptscotch
                        if (!refineDict.found("decompositionMethod"))
                        {
                            decompDict.add("method", word("ptscotch"));
                            decompDict.add
                            (
                                "numberOfSubdomains",
                                label(Pstream::nProcs())
                            );
                        }
                        else
                        {
                            // Read in method and coeffs
                            word decompMethod
                            (
                                refineDict.lookup("decompositionMethod")
                            );

                            decompDict.add("method", decompMethod);

                            decompDict.add
                            (
                                "numberOfSubdomains",
                                label(Pstream::nProcs())
                            );

                            if (refineDict.found(decompMethod + "Coeffs"))
                            {
                                decompDict.add
                                (
                                    word(decompMethod + "Coeffs"),
                                    refineDict.subDict
                                    (
                                        decompMethod + "Coeffs"
                                    )
                                );
                            }
                        }

                        // Add refinementHistory constraint
                        decompDict.add(word("constraints"), dictionary());

                        decompDict.subDict("constraints").add
                        (
                            word("refinement"),
                            dictionary()
                        );
                        decompDict.subDict("constraints").subDict("refinement").add
                        (
                            "type",
                            "refinementHistory"
                        );
                        decompDict.subDict("constraints").subDict("refinement").add
                        (
                            "enabled",
                            "true"
                        );

                        autoPtr<decompositionMethod> decomposerPtr
                        (
                            decompositionMethod::New(decompDict)
                        );

                        // Always preserve AMI patches if found
                        bool preserveAMI = false;

                        const fvPatchList& patches = mesh.boundary();
                        forAll(patches, patchi)
                        {
                            const fvPatch& currPatch = patches[patchi];

                            if (isType<cyclicAMIFvPatch>(currPatch))
                            {
                                // AMI not supported at present due to MPI issues
                                FatalErrorInFunction
                                    << "Dynamic load balancing does not support "
                                    << "boundary type cyclicAMI on patch "
                                    << currPatch.name()
                                    << exit(FatalError);
                                preserveAMI = true;
                                preserveAMI = returnReduce
                                (
                                    preserveAMI,
                                    orOp<bool>()
                                );

                                break;
                            }
                        }

                        if (preserveAMI)
                        {
                            // Preserve AMI faces
                            Info<< "Selected decomposition method to preserve AMI interfaces"
                                << " over a single processor " << endl;

                            labelList agglom(mesh.nCells(), -1);
                            pointField agglomPts(mesh.nCells());
                            scalarField weights(mesh.nCells());

                            DynamicList<label> groupCount(mesh.boundaryMesh().size());
                            DynamicList<point> groupCC(mesh.boundaryMesh().size());
                            DynamicList<DynamicList<label>> amiGroup;

                            const labelList& owner = mesh.faceOwner();

                            label coarseID = 0;

                            const labelList& visibleCells =
                                meshCutter().history().visibleCells();

                            forAll(mesh.boundaryMesh(), patchI)
                            {
                                const polyPatch& startpp = mesh.boundaryMesh()[patchI];
                                point aveCC = vector::zero;    // average cell center
                                label nAMICells = 0;

                                label shadowPatchI;

                                if (isA<cyclicAMIPolyPatch>(startpp))
                                {
                                    const cyclicAMIPolyPatch& ami =
                                        dynamic_cast<const cyclicAMIPolyPatch&>(startpp);

                                    shadowPatchI = ami.nbrPatchID();

                                    bool found = false;

                                    {
                                        const polyPatch& startpp = mesh.boundaryMesh()[patchI];
                                        forAll(startpp, i)
                                        {
                                            label own = owner[startpp.start() + i];
                                            if (agglom[own] == -1)
                                            {
                                                agglom[own] = coarseID;
                                                aveCC += mesh.cellCentres()[own];
                                                nAMICells++;
                                                found = true;
                                            }
                                        }
                                    }
                                    {
                                        const polyPatch& startpp =
                                            mesh.boundaryMesh()[shadowPatchI];
                                        forAll(startpp, i)
                                        {
                                            label own = owner[startpp.start() + i];
                                            if (agglom[own] == -1)
                                            {
                                                agglom[own] = coarseID;
                                                aveCC += mesh.cellCentres()[own];
                                                nAMICells++;
                                                found = true;
                                            }
                                        }
                                    }

                                    if (found)
                                    {
                                        amiGroup.append(DynamicList<label>());
                                        DynamicList<label> &group = amiGroup[amiGroup.size()-1];
                                        group.append(patchI);
                                        group.append(shadowPatchI);
                                        bool foundNbr = true;
                                        while (foundNbr)
                                        {
                                            foundNbr = false;
                                            for (label j = patchI; j < mesh.boundaryMesh().size(); j++)
                                            {
                                                const polyPatch& pp = mesh.boundaryMesh()[j];

                                                label sj;

                                                if (isA<cyclicAMIPolyPatch>(pp))
                                                {
                                                    const cyclicAMIPolyPatch& nami =
                                                        dynamic_cast<const cyclicAMIPolyPatch&>(pp);

                                                    sj = nami.nbrPatchID();

                                                    if
                                                    (
                                                        walkNeighbours
                                                        (
                                                            j,
                                                            sj,
                                                            agglom,
                                                            aveCC,
                                                            nAMICells
                                                         )
                                                    )
                                                    {
                                                        group.append(j);
                                                        group.append(sj);
                                                        foundNbr = true;
                                                    }
                                                }
                                            }
                                        }

                                        if (nAMICells)
                                        {
                                            aveCC /= nAMICells;
                                        }
                                        groupCC.append(aveCC);
                                        groupCount.append(nAMICells);
                                        coarseID++;
                                        group.shrink();
                                     }
                                }
                            }

                            amiGroup.shrink();
                            groupCount.shrink();
                            groupCC.shrink();

                            scalar amiFactor =
                                refineDict.lookupOrDefault<scalar>
                                (
                                    "amiFactor",
                                    1000.
                                );

                            scalar maxAMICells = mesh.nCells()*amiFactor;

                            agglom = -1;
                            coarseID = 0;

                            forAll(groupCount, i)
                            {
                                const labelList& group = amiGroup[i];
                                if (groupCount[i] > maxAMICells)
                                {
                                    Info<<"Too many cells for selected group of AMI's"
                                        <<" the following AMI patches will not span a "
                                        <<"single processor: " << endl;
                                    forAll(group, groupI)
                                    {
                                        label patchI = group[groupI];
                                        const polyPatch& pp = mesh.boundaryMesh()[patchI];
                                        Info<< "Group patch: " << pp.name() << endl;

                                        forAll(pp, i)
                                        {
                                            label own = owner[pp.start() + i];
                                            if (agglom[own] == -1)
                                            {
                                                agglomPts[coarseID] = mesh.cellCentres()[own];
                                                weights[coarseID] = 1.0;
                                                agglom[own] = coarseID++;
                                            }
                                        }
                                    }
                                }
                                else
                                {
                                    forAll(group, groupI)
                                    {
                                        label patchI = group[groupI];
                                        const polyPatch& pp = mesh.boundaryMesh()[patchI];

                                        forAll(pp, i)
                                        {
                                            label own = owner[pp.start() + i];
                                            agglom[own] = coarseID;
                                        }
                                    }
                                    agglomPts[coarseID] = groupCC[i];
                                    weights[coarseID] = groupCount[i];
                                    coarseID++;
                                }
                            }

                            // Update all unset points
                            // TODO mimic non AMI case
                            forAll(mesh.cells(), cellI)
                            {
                                if (agglom[cellI] == -1)
                                {
                                    // Cell has been refined
                                    if (visibleCells[cellI] > -1)
                                    {
                                        // Work to keep parent and children
                                        // on the same processor
                                        label pID = parentID
                                        (
                                            meshCutter().history().parentIndex(cellI)
                                        );

                                        label parentCoarseID = agglom[pID];

                                        if (parentCoarseID == -1)
                                        {
                                            agglomPts[coarseID] = mesh.cellCentres()[pID];
                                            weights[coarseID] = 1.;
                                            agglom[pID] = coarseID++;
                                        }

                                        // Child cell has same coarseID
                                        agglom[cellI] = agglom[pID];

                                        // Add to weights
                                        weights[agglom[cellI]]++;

                                    }
                                    // Cell is not refined
                                    else
                                    {
                                        agglomPts[coarseID] = mesh.cellCentres()[cellI];
                                        weights[coarseID] = 1.0;
                                        agglom[cellI] = coarseID++;
                                    }
                                 }
                            }

                            agglomPts.setSize(coarseID);
                            weights.setSize(coarseID);

                            // Create decomposer
                            decompositionMethod& decomposer = decomposerPtr();

                            // Create distribution
                            labelList distribution = decomposer.decompose
                            (
                                mesh,
                                agglom,
                                agglomPts,
                                weights
                            );

                            // Do actual sending/receiving of mesh
                            fvMeshDistribute distributor
                                (
                                    mesh,
                                    globalMeshData::matchTol_*mesh.bounds().mag()
                                );

                            autoPtr<mapDistributePolyMesh> map
                                = distributor.distribute(distribution);

                            meshCutter_.distribute(map);

                            loadImbalanced = true;

                        }
                        else
                        {
                            const labelIOList& cellLevel = meshCutter().cellLevel();
                            const labelList& visibleCells =
                                meshCutter().history().visibleCells();
                            Map<label> coarseIDmap(mesh.nCells());
                            labelList uniqueIndex(mesh.nCells(), 0);

                            label nCoarse = 0;

                            forAll(mesh.cells(), cellI)
                            {
                                //if (cellLevel[cellI] > 0 || visibleCells[cellI] > -1)
                                // Has cell been refined?

                                if (visibleCells[cellI] > -1)
                                {
                                    label pI = meshCutter().history().parentIndex(cellI);
                                    if (pI > -1)
                                    {
                                        uniqueIndex[cellI] = mesh.nCells() + parentID(pI);
                                    }
                                    else
                                    {
                                        uniqueIndex[cellI] = cellI;
                                    }
                                }
                                else
                                {
                                    uniqueIndex[cellI] = cellI;
                                }

                                if (coarseIDmap.insert(uniqueIndex[cellI], nCoarse))
                                {
                                    ++nCoarse;
                                }
                            }

                            // Convert to local sequential indexing and calculate coarse
                            // points and weights
                            labelList localIndex(mesh.nCells(), 0);
                            pointField coarsePoints(nCoarse, vector::zero);
                            scalarField coarseWeights(nCoarse, 0.0);

                            forAll(uniqueIndex, cellI)
                            {
                                localIndex[cellI] = coarseIDmap[uniqueIndex[cellI]];

                                label w = (1 << (3*cellLevel[cellI]));

                                coarseWeights[localIndex[cellI]] += 1.0;
                                coarsePoints[localIndex[cellI]] += mesh.C()[cellI]/w;
                            }

                            decompositionMethod& decomposer = decomposerPtr();

                            labelList distribution = decomposer.decompose
                            (
                                mesh,
                                localIndex,
                                coarsePoints,
                                coarseWeights
                            );

                            // Do actual sending/receiving of mesh
                            fvMeshDistribute distributor
                                (
                                    mesh,
                                    globalMeshData::matchTol_*mesh.bounds().mag()
                                );

                            autoPtr<mapDistributePolyMesh> map
                                = distributor.distribute(distribution);

                            meshCutter_.distribute(map);

                            loadImbalanced = true;
                        }

                        // Correct values on all cyclic patches
                        correctBoundaries<scalar>();

                        correctBoundaries<vector>();
                        correctBoundaries<vector1>();
                        correctBoundaries<vector2>();
                        correctBoundaries<vector4>();
                        correctBoundaries<vector6>();
                        correctBoundaries<vector8>();

                        correctBoundaries<diagTensor1>();
                        correctBoundaries<diagTensor2>();
                        correctBoundaries<diagTensor4>();
                        correctBoundaries<diagTensor6>();
                        correctBoundaries<diagTensor8>();

                        correctBoundaries<sphericalTensor>();
                        correctBoundaries<sphericalTensor1>();
                        correctBoundaries<sphericalTensor2>();
                        correctBoundaries<sphericalTensor4>();
                        correctBoundaries<sphericalTensor6>();
                        correctBoundaries<sphericalTensor8>();
                        correctBoundaries<symmTensor>();
                        correctBoundaries<tensor>();
                        correctBoundaries<tensor>();
                        correctBoundaries<tensor1>();
                        correctBoundaries<tensor2>();
                        correctBoundaries<tensor4>();
                        correctBoundaries<tensor6>();
                        correctBoundaries<tensor8>();

                        // TODO merge all coupled types in here to allow for rebalancing
                        Info<< "Balanced mesh in = "
                            << mesh.time().cpuTimeIncrement() << " s" << endl;

                        // Recalculate the protectedCell_ set
                        findProtectedCells();
                    }
                }
            }

            hasChanged = returnReduce(hasChanged, orOp<bool>());

            if
            (
                (hasChanged && !protectLayers_)
             || loadImbalanced
            )
            {
                mesh.topoChanging(hasChanged);

                const turbulenceModel& TM = mesh.lookupObject
                        <turbulenceModel>("turbulenceProperties");

                if (TM.type() != "laminar" && correctY_)
                {
                    // Update nearWallDist once
                    nearWallDist& y = const_cast<nearWallDist&>(TM.y());
                    Info<<" Correcting near wall distance y " << endl;
                    y.correct();
                }
            }

            nRefinementIterations_++;
        }
    }

    hasChanged = returnReduce(hasChanged, orOp<bool>());

    mesh.topoChanging(hasChanged);

    return hasChanged;
}


template<class Type>
void Foam::fvMeshTopoChangers::refiner::correctBoundaries()
{
    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;

    HashTable<GeoField*> flds(mesh().objectRegistry::lookupClass<GeoField>());

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        // Mimic "evaluate" but only for coupled patches (processor or cyclic)
        // and only for blocking or nonBlocking comms (no scheduled comms)
        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll(fld.boundaryField(), patchi)
            {
                if (fld.boundaryField()[patchi].coupled())
                {
                    fld.boundaryFieldRef()[patchi].initEvaluate
                    (
                        Pstream::defaultCommsType
                    );
                }
            }

            // Block for any outstanding requests
            if
            (
                Pstream::parRun()
             && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
            )
            {
                Pstream::waitRequests(nReq);
            }

            forAll(fld.boundaryField(), patchi)
            {
                if (fld.boundaryField()[patchi].coupled())
                {
                    fld.boundaryFieldRef()[patchi].evaluate
                    (
                        Pstream::defaultCommsType
                    );
                }
            }
        }
        else
        {
            // Scheduled patch updates not supported
            FatalErrorInFunction
                << "Unsuported communications type "
                << Pstream::commsTypeNames[Pstream::defaultCommsType]
                << exit(FatalError);
        }
    }
}


bool Foam::fvMeshTopoChangers::refiner::write(const bool write) const
{
    // Force refinement data to go to the current time directory
    const_cast<hexRef8&>(meshCutter_).setInstance(mesh().time().timeName());

    bool writeOk = meshCutter_.write(write);

    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh(),
            dimensionedScalar("cellLevel", dimless, 0)
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(cellLevel, celli)
        {
            scalarCellLevel[celli] = cellLevel[celli];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }

    return writeOk;
}


// ************************************************************************* //
