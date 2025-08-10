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
    (c) 2018 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/movingGIBTools/constraintBoundaryPoints/constraintBoundaryPoints.H"
#include "meshes/polyMesh/polyPatches/constraint/empty/emptyPolyPatch.H"
#include "meshes/polyMesh/polyPatches/constraint/wedge/wedgePolyPatch.H"
#include "meshes/polyMesh/polyPatches/indirectPolyPatch/indirectPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void constraintBoundaryPoints::createDataStructure()
{
    markBoundaryPoints();
}


void constraintBoundaryPoints::markBoundaryPoints()
{
    forAll(pMesh_.boundaryMesh(), patchI)
    {
        const polyPatch& pp = pMesh_.boundaryMesh()[patchI];
        if
        (
            !isA<indirectPolyPatch>(pp) &&
            !isA<wedgePolyPatch>(pMesh_.boundaryMesh()[patchI]) &&
            !isA<emptyPolyPatch>(pMesh_.boundaryMesh()[patchI]) &&
            !pp.coupled()
        )
        {
            const labelList& bPoints = pp.meshPoints();
            forAll(bPoints, pI)
            {
                const label& gpI = bPoints[pI];
                boundaryPoints_[gpI] = true;
            }
        }
    }
}


void constraintBoundaryPoints::markGIBPoints
(
    const labelList& fl
)
{
    forAll(fl, fI)
    {
        label fII = fl[fI];
        const face& faceI = pMesh_.faces()[fII];
        forAll(faceI, pI)
        {
            gibPoints_[faceI[pI]] = true;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

constraintBoundaryPoints::constraintBoundaryPoints
(
    const polyMesh& pMesh,
    const pointField& basePoints,
    const vectorField& baseSf
)
:
    pMesh_(pMesh),
    basePoints_(basePoints),
    baseSf_(baseSf),
    boundaryPoints_(basePoints.size(), false),
    gibPoints_(basePoints.size(), false),
    gibBoundaryPoints_(basePoints.size(), false)
{
    createDataStructure();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void constraintBoundaryPoints::correctPoints
(
    pointField& newPoints,
    const labelList& fl
)
{
    markGIBPoints(fl);

    gibBoundaryPoints_ = false;
    forAll(gibBoundaryPoints_, pI)
    {
        if (boundaryPoints_[pI] && gibPoints_[pI])
        {
            gibBoundaryPoints_[pI] = true;
        }
    }

    const labelListList& pFaces = pMesh_.pointFaces();
    forAll(gibBoundaryPoints_, pI)
    {
        if (gibBoundaryPoints_[pI])
        {
            vector dis = newPoints[pI] - basePoints_[pI];

            const labelList& pFaceI = pFaces[pI];
            forAll(pFaceI, fI)
            {
                label faceI = pFaceI[fI];
                if (faceI >= pMesh_.nInternalFaces())
                {
                    label patchI = pMesh_.boundaryMesh().whichPatch(faceI);
                    const polyPatch& poly = pMesh_.boundaryMesh()[patchI];
                    if (isA<directPolyPatch>(poly))
                    {
                        if
                        (
                            !(
                                pMesh_.boundaryMesh()[patchI].coupled() ||
                                isA<wedgePolyPatch>(pMesh_.boundaryMesh()[patchI]) ||
                                isA<emptyPolyPatch>(pMesh_.boundaryMesh()[patchI])
                            )
                        )
                        {
                            vector n = baseSf_[faceI]/mag(baseSf_[faceI]);
                            dis -= (dis&n)*n;
                        }
                    }
                }
            }
            newPoints[pI] = basePoints_[pI] + dis;
        }
    }
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
