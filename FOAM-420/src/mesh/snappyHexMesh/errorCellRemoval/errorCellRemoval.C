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
    (c) 2018-2020 Esi Ltd.
\*---------------------------------------------------------------------------*/

#include "errorCellRemoval/errorCellRemoval.H"
#include "regionSplit/regionSplit.H"
#include "polyTopoChange/removeCells/removeCells.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(errorCellRemoval, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * //

void Foam::errorCellRemoval::remove
(
    const bool updateIntersections
)
{
    Switch removeErrorCells = removalDict_.lookupOrDefault("removeCells", true);
    Switch writeVTK = removalDict_.lookupOrDefault("writeVTK", true);

    if (!writeVTK && !removeErrorCells)
    {
        return;
    }

    fvMesh& mesh = meshRefiner_.mesh();

    wordList removalChecks(0);

    if (removalDict_.found("checks"))
    {
        removalChecks = wordList(removalDict_.lookup("checks"));
    }
    else
    {
        removalChecks.setSize(3);
        removalChecks[0] = "volumes";
        removalChecks[1] = "facePyramids";
        removalChecks[2] = "weights";
    }

    scalar nonOrthoThreshold =
        removalDict_.lookupOrDefault("maxNonOrtho", 90.0);
    scalar maxAspectRatio =
        removalDict_.lookupOrDefault("maxAspectRatio", GREAT);
    scalar maxSkewness =
        removalDict_.lookupOrDefault("maxSkewness", GREAT);

    if (removalChecks.size() != 0)
    {
        Info<<"Removing/reporting of error cells : "
            << removalChecks <<endl;
    }
    else
    {
        return;
    }

    cellSet errorCells(mesh, "errorCells", mesh.nCells()/100+1);
    faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+1);

    List<direction> eType(mesh.nFaces(), NOERROR);

    forAll(removalChecks, checkI)
    {
        word check = removalChecks[checkI];

        if (check == "volumes")
        {
            labelHashSet oldErrorCells(errorCells);
            mesh.checkCellVolumes(true, &errorCells);
            if (writeVTK)
            {
                forAllConstIter(cellSet, errorCells, iter)
                {
                    label celli = iter.key();
                    if (!oldErrorCells.found(celli))
                    {
                        const cell& c = mesh.cells()[celli];
                        forAll(c, cfi)
                        {
                            eType[c[cfi]] |= VOL;
                        }
                    }
                }
            }
        }
        else if (check == "facePyramids")
        {
            labelHashSet oldErrorFaces(errorFaces);
            mesh.checkFacePyramids(true, -SMALL, &errorFaces);
            if (writeVTK)
            {
                forAllConstIter(faceSet, errorFaces, iter)
                {
                    label facei = iter.key();
                    if (!oldErrorFaces.found(facei))
                    {
                        eType[facei] |= PYR;
                    }
                }
            }
        }
        else if (check == "orthogonality")
        {
            labelHashSet oldErrorFaces(errorFaces);
            mesh.checkFaceOrthogonality(true, &errorFaces, nonOrthoThreshold);
            if (writeVTK)
            {
                forAllConstIter(faceSet, errorFaces, iter)
                {
                    label facei = iter.key();
                    if (!oldErrorFaces.found(facei))
                    {
                        eType[facei] |= ORTHO;
                    }
                }
            }
        }
        else if (check == "weights")
        {
            labelHashSet oldErrorFaces(errorFaces);
            const labelUList& owner = mesh.faceOwner();
            const labelUList& neighbour = mesh.faceNeighbour();
            const vectorField& faceCentres = mesh.faceCentres();
            const vectorField& cellCentres = mesh.cellCentres();
            const vectorField& faceAreas = mesh.faceAreas();
            const polyBoundaryMesh& patches = mesh.boundaryMesh();

            pointField neiCellCentres;
            syncTools::swapBoundaryCellList(mesh, cellCentres, neiCellCentres);

            forAll(mesh.faces(), facei)
            {
                label patchi = patches.whichPatch(facei);
                bool coupled = false;
                if (patchi != -1)
                {
                    if (patches[patchi].coupled())
                    {
                        coupled = true;
                    }
                    else
                    {
                        continue;
                    }
                }

                point neiCC;
                if (coupled)
                {
                    neiCC = neiCellCentres[facei-mesh.nInternalFaces()];
                }
                else
                {
                    neiCC = cellCentres[neighbour[facei]];
                }

                vector Sfi = faceAreas[facei];
                point fC = faceCentres[facei];
                point ownCC = cellCentres[owner[facei]];
                scalar SfdOwn = mag(Sfi & (fC - ownCC));
                scalar SfdNei = mag(Sfi & (neiCC - fC));

                if (mag(Sfi) < VSMALL || (SfdOwn + SfdNei) < VSMALL)
                {
                    errorFaces.insert(facei);
                }

            }

            if (writeVTK)
            {
                forAllConstIter(faceSet, errorFaces, iter)
                {
                    label facei = iter.key();
                    if (!oldErrorFaces.found(facei))
                    {
                        eType[facei] |= WEIGHTS;
                    }
                }
            }
        }
        else if (check == "aspectRatio")
        {
            const labelUList& owner = mesh.faceOwner();
            const labelUList& neighbour = mesh.faceNeighbour();
            const vectorField& faceAreas = mesh.faceAreas();
            const scalarField& cellVolumes = mesh.cellVolumes();
            labelHashSet oldErrorCells(errorCells);
            vectorField sumMagClosed(mesh.nCells(), Zero);

            forAll(owner, facei)
            {
                // Add to owner
                label own = owner[facei];
                sumMagClosed[own] += cmptMag(faceAreas[facei]);
            }

            forAll(neighbour, facei)
            {
                label nei = neighbour[facei];
                sumMagClosed[nei] += cmptMag(faceAreas[facei]);
            }

            forAll(mesh.cells(), celli)
            {
                scalar minCmpt = VGREAT;
                scalar maxCmpt = -VGREAT;
                for (direction dir = 0; dir < vector::nComponents; dir++)
                {
                    minCmpt = min(minCmpt, sumMagClosed[celli][dir]);
                    maxCmpt = max(maxCmpt, sumMagClosed[celli][dir]);
                }

                scalar aspectRatio = maxCmpt/(minCmpt + ROOTVSMALL);
                scalar v = max(ROOTVSMALL,  cellVolumes[celli]);

                aspectRatio = max
                (
                    aspectRatio,
                    1.0/6.0*cmptSum(sumMagClosed[celli])/pow(v, 2.0/3.0)
                );

                if (aspectRatio > maxAspectRatio)
                {
                    errorCells.insert(celli);
                }
            }

            if (writeVTK)
            {
                forAllConstIter(cellSet, errorCells, iter)
                {
                    label celli = iter.key();
                    if (!oldErrorCells.found(celli))
                    {
                        const cell& c = mesh.cells()[celli];
                        forAll(c, cfi)
                        {
                            eType[c[cfi]] |= AR;
                        }
                    }
                }
            }
        }
        else if (check == "skewness")
        {
            labelHashSet oldErrorFaces(errorFaces);
            mesh.checkFaceSkewness(true, &errorFaces, maxSkewness);
            if (writeVTK)
            {
                forAllConstIter(faceSet, errorFaces, iter)
                {
                    label facei = iter.key();
                    if (!oldErrorFaces.found(facei))
                    {
                        eType[facei] |= SKEW;
                    }
                }
            }
        }
        else
        {
            WarningInFunction
                << "Could not find removal method" << check
                << "Availble methods volumes, facePyramids "
                << "orthogonality, weights (zero), aspectRatio, "
                << "and skewness."<<endl;
        }
    }

    if (writeVTK)
    {
        DynamicList<label> outputFaces(mesh.nFaces()/10);
        forAll(mesh.faces(), facei)
        {
            if
            (
                eType[facei] & ORTHO
                || eType[facei] & PYR
                || eType[facei] & VOL
                || eType[facei] & WEIGHTS
                || eType[facei] & AR
                || eType[facei] & SKEW
             )
            {
                outputFaces.append(facei);
            }
        }

        autoPtr<indirectPrimitivePatch> errorPP
        (
            new indirectPrimitivePatch
            (
                IndirectList<face>(mesh.faces(), outputFaces),
                mesh.points()
             )
        );

        simpleVTKWriter errorVTK
        (
            errorPP().localFaces(),
            errorPP().localPoints()
        );

        labelList orthoFaceData(errorPP().size(), 0);
        labelList volFaceData(errorPP().size(), 0);
        labelList pyrFaceData(errorPP().size(), 0);
        labelList weightsFaceData(errorPP().size(), 0);
        labelList skewFaceData(errorPP().size(), 0);
        labelList arFaceData(errorPP().size(), 0);

        forAll(errorPP(), i)
        {
            label facei = errorPP().addressing()[i];
            if (eType[facei] & ORTHO)
            {
                orthoFaceData[i] = 1;
            }
            if (eType[facei] & VOL)
            {
                volFaceData[i] = 1;
            }
            if (eType[facei] & PYR)
            {
                pyrFaceData[i] = 1;
            }
            if (eType[facei] & WEIGHTS)
            {
                weightsFaceData[i] = 1;
            }
            if (eType[facei] & AR)
            {
                arFaceData[i] = 1;
            }
            if (eType[facei] & SKEW)
            {
                skewFaceData[i] = 1;
            }
        }

        errorVTK.addFaceData("ortho", orthoFaceData);
        errorVTK.addFaceData("pyr", pyrFaceData);
        errorVTK.addFaceData("vol", volFaceData);
        errorVTK.addFaceData("weights", weightsFaceData);
        errorVTK.addFaceData("aspectRatio", arFaceData);
        errorVTK.addFaceData("skewness", skewFaceData);

        errorVTK.write("errorFaces.vtk");
    }

    if (removeErrorCells)
    {
        boolList errors(mesh.nCells(), false);
        boolList blockedFace(mesh.nFaces(),false);
        forAllConstIter(labelHashSet, errorFaces, iter)
        {
            label own = mesh.faceOwner()[iter.key()];
            errors[own] = true;

            const labelList& cOwn =  mesh.cells()[own];
            forAll(cOwn, cFI)
            {
                blockedFace[cOwn[cFI]] = true;
            }

            if (iter.key() < mesh.nInternalFaces())
            {
                label nei = mesh.faceNeighbour()[iter.key()];
                errors[nei] = true;

                const labelList& cNei =  mesh.cells()[nei];
                forAll(cNei, cFI)
                {
                    blockedFace[cNei[cFI]] = true;
                }
            }
        }

        forAllConstIter(labelHashSet, errorCells, iter)
        {
            label cellI = iter.key();
            errors[cellI] = true;
            const labelList& c =  mesh.cells()[cellI];
            forAll(c, cFI)
            {
                blockedFace[c[cFI]] = true;
            }
        }

        syncTools::syncFaceList(mesh,blockedFace,orEqOp<bool>());

        // Set region per cell based on walking
        regionSplit cellRegion(mesh, blockedFace);
        labelHashSet keepRegionSet(meshRefiner_.keepLargestRegions(cellRegion));

        DynamicList<label> cellsToRemove(mesh.nCells()/100);
        forAll(mesh.cells(), cellI)
        {
            if (errors[cellI] || !keepRegionSet.found(cellRegion[cellI]))
            {
                cellsToRemove.append(cellI);
            }
        }
        cellsToRemove.shrink();

        label nCellsToRemove =
            returnReduce(cellsToRemove.size(), sumOp<label>());
        Info<<"Selected " << nCellsToRemove << " cells for removal. "<<endl;

        if (nCellsToRemove != 0)
        {
            dictionary patchInfo;
            patchInfo.set("type", wallPolyPatch::typeName);
            label removedID = meshRefiner_.addPatch
            (
                mesh,
                removalDict_.lookupOrDefault<word>("patch","oldInternalFaces"),
                patchInfo
             );

            // Remove cells
            removeCells cellRemover(mesh);

            // Pick up patches for exposed faces
            labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));
            labelList exposedPatchIDs(exposedFaces.size(),removedID);

            meshRefiner_.doRemoveCells
            (
                cellsToRemove,
                exposedFaces,
                exposedPatchIDs,
                cellRemover,
                updateIntersections
            );
        }
    }

    return;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::errorCellRemoval::errorCellRemoval
(
    meshRefinement& meshRefiner,
    const dictionary& dict
)
:
    meshRefiner_(meshRefiner),
    removalDict_(dict)
{}

// ************************************************************************* //
