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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "patchDist/patchDistFuncs/patchDistFuncs.H"
#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<template<class> class MapType>
void Foam::patchDistFuncs::correctBoundaryFaceCells
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    MapType<labelPair>& nearestPatchAndFace
)
{
    // Size neighbours array for maximum possible (= size of largest patch)
    DynamicList<label> neighbours(maxPatchSize(mesh, patchIDs));

    // Correct all cells with face on wall
    const vectorField& cellCentres = mesh.cellCentres();

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Check cells with face on wall
        forAll(patch, patchFacei)
        {
            const label celli = patch.faceCells()[patchFacei];

            getPointNeighbours
            (
                patch,
                patchFacei,
                neighbours
            );

            label minPatchFacei = -1;

            wallDistCorrected[celli] =
                smallestDist
                (
                    cellCentres[celli],
                    patch,
                    neighbours,
                    minPatchFacei
                );

            // Store wallCell and its nearest neighbour
            nearestPatchAndFace.insert
            (
                celli,
                labelPair(patchi, minPatchFacei)
            );
        }
    }
}


template<template<class> class PatchField, template<class> class MapType>
void Foam::patchDistFuncs::correctBoundaryFaceFaceCells
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    FieldField<PatchField, scalar>& wallDistCorrected,
    MapType<labelPair>& nearestPatchAndFace
)
{
    // Size neighbours array for maximum possible (= size of largest patch)
    DynamicList<label> neighbours(maxPatchSize(mesh, patchIDs));

    // Correct all faces on a wall
    const vectorField& cellCentres = mesh.cellCentres();

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        // Check cells with face on wall
        forAll(patch, patchFacei)
        {
            const label celli = patch.faceCells()[patchFacei];

            getPointNeighbours
            (
                patch,
                patchFacei,
                neighbours
            );

            label minPatchFacei = -1;

            wallDistCorrected[patchi][patchFacei] =
                smallestDist
                (
                    cellCentres[celli],
                    patch,
                    neighbours,
                    minPatchFacei
                );

            // Store wallCell and its nearest neighbour
            nearestPatchAndFace.insert
            (
                celli,
                labelPair(patchi, minPatchFacei)
            );
        }
    }
}


template<template<class> class MapType>
void Foam::patchDistFuncs::correctBoundaryPointCells
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& wallDistCorrected,
    MapType<labelPair>& nearestPatchAndFace
)
{
    // Correct all (non-visited) cells with point on wall
    const vectorField& cellCentres = mesh.cellCentres();

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        const label patchi = iter.key();
        const polyPatch& patch = mesh.boundaryMesh()[patchi];

        const labelList& meshPoints = patch.meshPoints();
        const labelListList& pointFaces = patch.pointFaces();

        forAll(meshPoints, meshPointi)
        {
            const labelList& neighbours =
                mesh.pointCells(meshPoints[meshPointi]);

            forAll(neighbours, neighbourI)
            {
                const label celli = neighbours[neighbourI];

                if (!nearestPatchAndFace.found(celli))
                {
                    const labelList& wallFaces = pointFaces[meshPointi];

                    label minPatchFacei = -1;

                    wallDistCorrected[celli] =
                        smallestDist
                        (
                            cellCentres[celli],
                            patch,
                            wallFaces,
                            minPatchFacei
                        );

                    // Store wallCell and its nearest neighbour
                    nearestPatchAndFace.insert
                    (
                        celli,
                        labelPair(patchi, minPatchFacei)
                    );
                }
            }
        }
    }
}


// ************************************************************************* //
