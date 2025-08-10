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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "pointPatchDist/pointPatchDist.H"
#include "pointPatchDist/externalPointEdgePoint.H"
#include "meshes/pointMesh/pointMesh.H"
#include "algorithms/PointEdgeWave/PointEdgeWave.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointPatchDist::pointPatchDist
(
    const pointMesh& pMesh,
    const labelHashSet& patchIDs,
    const pointField& points
)
:
    pointScalarField
    (
        IOobject
        (
            "pointDistance",
            pMesh.db().time().timeName(),
            pMesh.db()
        ),
        pMesh,
        dimensionedScalar("y", dimLength, GREAT)
    ),
    points_(points),
    patchIDs_(patchIDs),
    nUnset_(0)
{
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pointPatchDist::~pointPatchDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointPatchDist::correct()
{
    const pointBoundaryMesh& pbm = mesh().boundary();

    label nPoints = 0;

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchi = iter.key();
        nPoints += pbm[patchi].meshPoints().size();
    }

    externalPointEdgePoint::trackingData td(points_);

    // Set initial changed points to all the patch points(if patch present)
    List<externalPointEdgePoint> wallInfo(nPoints);
    labelList wallPoints(nPoints);
    nPoints = 0;

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        label patchi = iter.key();
        // Retrieve the patch now we have its index in patches.

        const labelList& mp = pbm[patchi].meshPoints();

        forAll(mp, ppI)
        {
            label meshPointi = mp[ppI];
            wallPoints[nPoints] = meshPointi;
            wallInfo[nPoints] = externalPointEdgePoint
            (
                td.points_[meshPointi],
                0.0
            );
            nPoints++;
        }
    }

    // Current info on points
    List<externalPointEdgePoint> allPointInfo(mesh()().nPoints());

    // Current info on edges
    List<externalPointEdgePoint> allEdgeInfo(mesh()().nEdges());

    PointEdgeWave
    <
        externalPointEdgePoint,
        externalPointEdgePoint::trackingData
    > wallCalc
    (
        mesh()(),
        wallPoints,
        wallInfo,

        allPointInfo,
        allEdgeInfo,
        mesh().globalData().nTotalPoints(), // max iterations
        td
    );

    pointScalarField& psf = *this;


    forAll(allPointInfo, pointi)
    {
        if (allPointInfo[pointi].valid(td))
        {
            psf[pointi] = Foam::sqrt(allPointInfo[pointi].distSqr());
        }
        else
        {
            nUnset_++;
        }
    }
}


// ************************************************************************* //
