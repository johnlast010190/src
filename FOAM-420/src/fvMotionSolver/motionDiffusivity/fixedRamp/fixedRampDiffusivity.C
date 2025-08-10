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
    (c) 2011 OpenFOAM Foundation
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionDiffusivity/fixedRamp/fixedRampDiffusivity.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "patchDist/wallPoint/wallPoint.H"
#include "algorithms/FaceCellWave/FaceCellWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fixedRampDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        fixedRampDiffusivity,
        Istream
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::fixedRampDiffusivity::fDiffusivity(scalar y)
{
    if (y < Lbuff_)
    {
        return DbuffRatio_;
    }
    else if (y > 2*Lbuff_)
    {
        return 1;
    }

    return cubicF_.value(y);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedRampDiffusivity::fixedRampDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    uniformDiffusivity(mesh, mdData),
    patchNames_(mdData),
    Lbuff_(readScalar(mdData)),
    DbuffRatio_(readScalar(mdData)),
    cubicF_()
{
    scalar x1 = Lbuff_;
    scalar x12 = sqr(x1);
    scalar x13 = x12*x1;
    scalar x2 = 2*x1;
    scalar x22 = sqr(x2);
    scalar x23 = x22*x2;

    scalar cf3Rcp = 3*x12*x2 - x23 -2*x13
        - 3.0/2.0 *(x12-x22)*(2*x1*x2-x22-x12)/(x1 - x2);

    cubicF_[3] = (DbuffRatio_ - 1) / cf3Rcp;

    cubicF_[2] = -3.0/2.0*cubicF_[3]*(x12-x22)/(x1-x2);

    cubicF_[1] = -2*cubicF_[2]*x2 - 3*cubicF_[3]*x22;

    cubicF_[0] = DbuffRatio_ - cubicF_[1]*x1 - cubicF_[2]*x12 - cubicF_[3]*x13;

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fixedRampDiffusivity::~fixedRampDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedRampDiffusivity::correct()
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

        forAll(fc, patchFaceI)
        {
            changedFaces[nPatchFaces] = patch.start() + patchFaceI;

            faceDist[nPatchFaces] = wallPoint(fc[patchFaceI], 0);

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

    for (label faceI=0; faceI<mesh().nInternalFaces(); faceI++)
    {
        scalar dist = sqrt(faceInfo[faceI].distSqr());

        faceDiffusivity_[faceI] = fDiffusivity(dist);
    }

    forAll(faceDiffusivity_.boundaryField(), patchI)
    {
        fvsPatchScalarField& bfld = faceDiffusivity_.boundaryFieldRef()[patchI];

        if (patchSet.found(patchI))
        {
            forAll(bfld, i)
            {
                bfld[i] = DbuffRatio_;
            }
        }
        else
        {
            const label start = bfld.patch().start();

            forAll(bfld, i)
            {
                scalar dist = sqrt(faceInfo[start+i].distSqr());
                bfld[i] = fDiffusivity(dist);
            }
        }
    }
}


// ************************************************************************* //
