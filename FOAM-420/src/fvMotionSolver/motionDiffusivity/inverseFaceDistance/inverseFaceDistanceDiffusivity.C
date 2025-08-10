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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionDiffusivity/inverseFaceDistance/inverseFaceDistanceDiffusivity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "patchDist/wallPoint/wallPoint.H"
#include "algorithms/FaceCellWave/FaceCellWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseFaceDistanceDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        inverseFaceDistanceDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseFaceDistanceDiffusivity::inverseFaceDistanceDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    uniformDiffusivity(mesh, mdData),
    patchNames_(mdData)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseFaceDistanceDiffusivity::~inverseFaceDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inverseFaceDistanceDiffusivity::correct()
{
    const polyBoundaryMesh& bdry = mesh().boundaryMesh();

    labelHashSet patchSet(bdry.size());

    label nPatchFaces = 0;

    forAll(patchNames_, i)
    {
        const label pID = bdry.findPatchID(patchNames_[i]);

        if (pID > -1)
        {
            patchSet.insert(pID);
            nPatchFaces += bdry[pID].size();
        }
    }

    List<wallPoint> faceDist(nPatchFaces);
    labelList changedFaces(nPatchFaces);

    nPatchFaces = 0;

    forAllConstIter(labelHashSet, patchSet, iter)
    {
        const polyPatch& patch = bdry[iter.key()];

        const vectorField::subField fc(patch.faceCentres());

        forAll(fc, patchFacei)
        {
            changedFaces[nPatchFaces] = patch.start() + patchFacei;

            faceDist[nPatchFaces] = wallPoint(fc[patchFacei], 0);

            nPatchFaces++;
        }
    }
    faceDist.setSize(nPatchFaces);
    changedFaces.setSize(nPatchFaces);

    List<wallPoint> faceInfo(mesh().nFaces()), cellInfo(mesh().nCells());
    FaceCellWave<wallPoint> waveInfo
    (
        mesh(),
        changedFaces,
        faceDist,
        faceInfo,
        cellInfo,
        mesh().globalData().nTotalCells() + 1    // max iterations
    );

    for (label facei=0; facei<mesh().nInternalFaces(); facei++)
    {
        scalar dist = faceInfo[facei].distSqr();

        faceDiffusivity_[facei] = 1.0/sqrt(dist);
    }

    surfaceScalarField::Boundary& faceDiffusivityBf =
        faceDiffusivity_.boundaryFieldRef();

    forAll(faceDiffusivityBf, patchi)
    {
        fvsPatchScalarField& bfld = faceDiffusivityBf[patchi];

        const labelUList& faceCells = bfld.patch().faceCells();

        if (patchSet.found(patchi))
        {
            forAll(bfld, i)
            {
                scalar dist = cellInfo[faceCells[i]].distSqr();
                bfld[i] = 1.0/sqrt(dist);
            }
        }
        else
        {
            const label start = bfld.patch().start();

            forAll(bfld, i)
            {
                scalar dist = faceInfo[start+i].distSqr();
                bfld[i] = 1.0/sqrt(dist);
            }
        }
    }
}


// ************************************************************************* //
