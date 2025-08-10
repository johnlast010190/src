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

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
Foam::patchDistFuncs::NoMap<Foam::labelPair>
Foam::patchDistFuncs::NoMap<Foam::labelPair>::null =
Foam::patchDistFuncs::NoMap<Foam::labelPair>();


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::patchDistFuncs::smallestDist
(
    const point& p,
    const polyPatch& patch,
    const labelUList& wallFaces,
    label& minPatchFacei
)
{
    const pointField& points = patch.points();

    scalar minDist = GREAT;
    minPatchFacei = -1;

    forAll(wallFaces, wallFacei)
    {
        const label patchFacei = wallFaces[wallFacei];

        pointHit curHit = patch[patchFacei].nearestPoint(p, points);

        if (curHit.distance() < minDist)
        {
            minDist = curHit.distance();
            minPatchFacei = patchFacei;
        }
    }

    return minDist;
}


void Foam::patchDistFuncs::getPointNeighbours
(
    const polyPatch& patch,
    const label patchFacei,
    DynamicList<label>& neighbours
)
{
    neighbours.clear();

    // Add myself
    neighbours.append(patchFacei);

    // Add all face neighbours
    const labelList& faceNeighbours = patch.faceFaces()[patchFacei];
    forAll(faceNeighbours, faceNeighbourI)
    {
        // neighbours.append(faceNeighbours[faceNeighbourI]);
        label nbrFaceI = faceNeighbours[faceNeighbourI];

        if
        (
            findIndex
            (
                SubList<label>(neighbours, neighbours.size()),
                nbrFaceI
            ) == -1
        )
        {
            neighbours.append(nbrFaceI);
        }
    }

    // Add all point-only neighbours by linear searching in edge neighbours.
    // Assumes that point-only neighbours are not using multiple points on
    // face.
    const face& f = patch.localFaces()[patchFacei];
    forAll(f, fp)
    {
        const label pointi = f[fp];

        const labelList& pointNbs = patch.pointFaces()[pointi];

        forAll(pointNbs, nbI)
        {
            const label facei = pointNbs[nbI];

            // Check for facei in edge-neighbours part of neighbours
            if
            (
                findIndex
                (
                    SubList<label>(neighbours, neighbours.size()),
                    facei
                ) == -1
            )
            {
                neighbours.append(facei);
            }
        }
    }
}


Foam::label Foam::patchDistFuncs::maxPatchSize
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs
)
{
    label maxSize = 0;

    forAllConstIter(labelHashSet, patchIDs, iter)
    {
        maxSize = Foam::max(maxSize, mesh.boundaryMesh()[iter.key()].size());
    }

    return maxSize;
}


// ************************************************************************* //
