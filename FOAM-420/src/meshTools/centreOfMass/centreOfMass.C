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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "centreOfMass/centreOfMass.H"
#include "meshes/meshShapes/cell/pyramidPointFaceRef.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::point Foam::centreOfMass::centre
(
    const polyMesh& mesh,
    const labelHashSet& patchSet
)
{
    vector sumVc = Zero;
    scalar sumV = 0.0;
    scalar patchVol(0.);
    vector tempC = Zero;
    scalar tempA = 0.0;
    vector sumAc = Zero;
    scalar sumA = 0.0;

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
    const faceList& f = mesh.faces();
    const pointField& pts = mesh.points();

    //get area weigthed average mean faces center
    //objective: cEst needs to be in the symmetry plane for 2D cases
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        label patchi = iter.key();
        const word& patchName = bMesh[patchi].name();

        label patchI = bMesh.findPatchID(patchName);

        if (patchI != -1)
        {
            const polyPatch& patch = bMesh[patchI];
            label start = patch.start();

            forAll(patch, pfi)
            {
                label faceI = start + pfi;

                tempC = mesh.faceCentres()[faceI];
                tempA = mesh.magFaceAreas()[faceI];

                sumAc += tempC*tempA;
                sumA += tempA;
            }
        }
    }

    reduce(
        std::tie(sumAc, sumA),
        ParallelOp<sumOp<vector>, sumOp<scalar>>{}
    );

    point cEst = sumAc/(sumA + VSMALL);

    //get volume of pyramids
    forAllConstIter(labelHashSet, patchSet, iter)
    {
        label patchi = iter.key();
        const word& patchName = bMesh[patchi].name();

        label patchI = bMesh.findPatchID(patchName);

        if (patchI != -1)
        {
            const polyPatch& patch = bMesh[patchI];
            label start = patch.start();

            forAll(patch, pfi)
            {
                label faceI = start + pfi;
                pyramidPointFaceRef pPFRef(f[faceI], cEst);

                patchVol = pPFRef.mag(pts);
                vector patchCentre = pPFRef.centre(pts);

                sumVc += patchVol*patchCentre;
                sumV += patchVol;
            }
        }
    }

    reduce(
        std::tie(sumVc, sumV),
        ParallelOp<sumOp<vector>, sumOp<scalar>>{}
    );

    return sumVc/(sumV + VSMALL);
}


// ************************************************************************* //
