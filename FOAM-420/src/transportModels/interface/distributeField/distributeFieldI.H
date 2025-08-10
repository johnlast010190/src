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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template <typename T>
inline void Foam::distributeField::distributeVolField
(
    GeometricField<T, fvPatchField, volMesh>& field,
    volScalarField& currField,
    const volScalarField& destField,
    int maxLoops,
    int minLoops,
    bool normalizeField
)const
{
    const fvMesh& mesh_ = field.mesh();

    // copy to contain current distribution
    volScalarField currDist(currField);
    bool finished = false;
    int iLoop = 0;

    while (((!finished) && (iLoop < maxLoops)) || iLoop < minLoops)
    {
        //- marker field for cells to be filled next
        volScalarField toBeDistributed(currField*0);

        //- pass field to neighbours if not yet set there
        forAll(currDist, cellI)
        {
            if (currDist[cellI])
            {
                labelList cellStencil = mesh_.cellCells()[cellI];
                forAll(cellStencil, i)
                {
                    label nei = cellStencil[i];
                    if (!currDist[nei])
                    {
                        toBeDistributed[nei] += 1;
                        field[nei] += field[cellI];
                    }
                }
            }
        }

        currDist.correctBoundaryConditions();
        field.correctBoundaryConditions();

        // boundary treatment
        const polyBoundaryMesh& patches = mesh_.boundaryMesh();
        forAll(patches, patchi)
        {
            const polyPatch& pp = patches[patchi];
            // parallel communication
            if (pp.coupled())
            {
                tmp<scalarField> tcurrDistBP
                (
                    currDist.boundaryField()[patchi].patchInternalField()
                );
                const scalarField& currDistBP = tcurrDistBP();
                tmp<scalarField> tcurrDistBN
                (
                    currDist.boundaryField()[patchi].patchNeighbourField()
                );
                const scalarField& currDistBN = tcurrDistBN();
                tmp<Field<T>> tfieldBN
                (
                     field.boundaryFieldRef()[patchi].patchNeighbourField()
                );
                const Field<T>& fieldBN = tfieldBN();

                // loop over all faces of proc-I coupled boundary and check patchNeighbourField
                label iFaceGlobal = pp.start();
                forAll(pp, facei)
                {
                    // if currDist in neighbour true and in interalField (=owner) false, mark owner cell toBeDistributed
                    if (currDistBN[facei] && (!currDistBP[facei]))
                    {
                        label own = mesh_.faceOwner()[iFaceGlobal];
                        toBeDistributed[own] += 1;
                        field[own] += fieldBN[facei];
                    }
                    iFaceGlobal++;
                }
            }
            else  // write cell value to non-coupled boundary
            {
                if (pp.type() != "empty")
                {
                    label iFaceGlobal = pp.start();

                    forAll(pp, facei)
                    {
                        label own = mesh_.faceOwner()[iFaceGlobal];
                        if (currDist[own] == true)
                        {
                            field.boundaryFieldRef()[patchi][facei] = field[own];
                        }
                        iFaceGlobal++;
                    }
                }
            }
        }

        // perform averaging of field and update current distribution according to newly set cells
        forAll(field, celli)
        {
            if (toBeDistributed[celli])
            {
                //- finish avarage calculation
                if (normalizeField)
                {
                    field[celli] /= mag(field[celli]) + VSMALL;
                }
                else
                {
                    field[celli] /= toBeDistributed[celli];
                }

                //re-set current distribution
                currDist[celli] = 1;
            }
        }

        //- check if further loop is required
        finished = true;
        forAll(currDist, celli)
        {
            // if one target distribution cell not yet distributed -> continue
            if (destField[celli] && (!currDist[celli]))
            {
                finished = false;
                break;
            }
        }
        reduce(finished, andOp<bool>());
        iLoop++;
    }

    field.correctBoundaryConditions();
}


inline void Foam::distributeField::calculateDistance
(
    boolList& currField,
    volScalarField& distance,
    volVectorField& nHatv,
    const boolList& destField,
    scalar maxDistance,
    int maxLoops,
    int minLoops
)
{
    const fvMesh& mesh_ = distance.mesh();

    // list of cells of which information is distributed in current loop
    boolList distribute(currField);

    bool finished = false;

    int iLoop = 0;
    while (((!finished) && (iLoop < maxLoops)) || iLoop < minLoops)
    {

        boolList   updated   (mesh_.nCells(), false);
        scalarList distancef (mesh_.nFaces(), 0.0);

        forAll(nHatv, iCell)
        {
            if (distribute[iCell] == true)
            {
                //- mark unset neighbour cells
                const labelList& cellCells = mesh_.cellCells()[iCell];
                forAll(cellCells, i)
                {
                    label neigh = cellCells[i];
                    if (currField[neigh] == false)
                    {
                        nHatv[neigh] += nHatv[iCell];
                        updated[neigh] = true;
                    }
                }

                //- calculate distance face-wise
                const labelList& faces = mesh_.cells()[iCell];
                vector xp = mesh_.C()[iCell];
                forAll(faces, i)
                {
                    label iFace = faces[i];
                    vector xf = mesh_.faceCentres()[iFace];
                    distancef[iFace] = distance[iCell] + (nHatv[iCell] & (xf - xp));
                }
            }
        }

        //- handling of coupled boundaries
        List<vector> nHatProc (mesh_.nFaces(), vector(0.0,0.0,0.0));
        boolList setProc   (mesh_.nFaces(), false);
        forAll(nHatv.boundaryField(), iPatch)
        {
            if (
                (nHatv.boundaryField()[iPatch].type() == "processor")
                || (nHatv.boundaryField()[iPatch].type() == "cyclic")
                || (nHatv.boundaryField()[iPatch].type() == "cyclicAMI")
            )
            {
                forAll(nHatv.boundaryField()[iPatch], iFace)
                {
                    label iFaceGlobal = iFace + mesh_.boundaryMesh()[iPatch].start();
                    label own = mesh_.faceOwner()[iFaceGlobal];
                    if (distribute[own] == true)
                    {
                        nHatProc[iFaceGlobal] = nHatv[own];
                        setProc[iFaceGlobal] = true;
                    }
                }
            }
            else
            {
                //- store value of nHats and distance on (real) boundaries
                forAll(nHatv.boundaryField()[iPatch], iFace)
                {
                    label iFaceGlobal = iFace + mesh_.boundaryMesh()[iPatch].start();
                    label own = mesh_.faceOwner()[iFaceGlobal];
                    if (distribute[own] == true)
                    {
                        nHatv.boundaryFieldRef()[iPatch][iFace] = nHatv[own];
                        distance.boundaryFieldRef()[iPatch][iFace] = distancef[iFaceGlobal];
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, nHatProc, plusEqOp<vector>());
        syncTools::syncFaceList(mesh_, distancef, plusEqOp<scalar>());
        syncTools::syncFaceList(mesh_, setProc, orEqOp<bool>());

        //- handling of coupled boundaries
        forAll(nHatv.boundaryField(), iPatch)
        {
            if (
                (nHatv.boundaryField()[iPatch].type() == "processor")
                || (nHatv.boundaryField()[iPatch].type() == "cyclic")
                || (nHatv.boundaryField()[iPatch].type() == "cyclicAMI")
            )
            {
                forAll(nHatv.boundaryField()[iPatch], iFace)
                {
                    label iFaceGlobal = iFace + mesh_.boundaryMesh()[iPatch].start();
                    label own = mesh_.faceOwner()[iFaceGlobal];
                    if ((currField[own] == false) && (setProc[iFaceGlobal]))
                    {
                        nHatv[own] += nHatProc[iFaceGlobal];
                        updated[own] = true;
                    }
                }
            }
        }

        forAll(nHatv, iCell)
        {
            //- if cell has already distributed its information it has
            //  not to be distributed again in next loop
            if (distribute[iCell] == true)
            {
                distribute[iCell] = false;
            }

            //- if cell value was updated, distribute it in the next loop
            if (updated[iCell] == true)
            {
                nHatv[iCell] /= mag(nHatv[iCell])+VSMALL;
                distribute[iCell] = true;
                currField[iCell] = true;

                //- calculate distance, get nearest face-interface distance and
                //  add distance to cell center
                const labelList& faces = mesh_.cells()[iCell];
                vector xp = mesh_.C()[iCell];
                forAll(faces, i)
                {
                    label iFace = faces[i];

                    scalar& dist = distance[iCell];
                    if (distancef[iFace] != 0.0) // NOT ALWAYS VALID!!!
                    {
                        vector xf = mesh_.faceCentres()[iFace];
                        scalar newDist = distancef[iFace] + (nHatv[iCell] & (xp - xf));
                        //- smaller distance is chosen
                        if ((mag(dist) < SMALL) || (mag(newDist) < mag(dist)))
                        {
                            dist = newDist;
                        }
                    }
                }

                //- check if distribution can be stopped here
                if (
                        (mag(distance[iCell]) > maxDistance)
                        && ((currField[iCell] == true) && (destField[iCell] == false))

                    )
                {
                    distribute[iCell] = false;
                }
            }
        }
        iLoop++;

        //- check if further loop is required
        finished = true;
        forAll(distribute, iCell)
        {
            if (distribute[iCell])
            {
                finished = false;
            }
        }
        reduce(finished, andOp<bool>());
    }
    distance.correctBoundaryConditions();
    nHatv.correctBoundaryConditions();
}

// ************************************************************************* //
