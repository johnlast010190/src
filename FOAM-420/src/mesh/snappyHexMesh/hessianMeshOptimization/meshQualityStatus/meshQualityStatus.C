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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/meshQualityStatus/meshQualityStatus.H"

namespace Foam
{
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::meshQualityStatus::findWrongFaces(labelHashSet& wrongFaces)
{
    tmp<vectorField> cellCentres = calculateCorrectCellCentres();
    dictionary qualityControls = dict_.subDict("meshQualityControls");

    List<labelPair> emptyBaffles;

    Foam::polyMeshGeometry::checkFaceDotProduct
    (
        false,
        qualityControls.lookupOrDefault<scalar>("maxNonOrtho", 70, true),
        qualityControls.lookupOrDefault<scalar>("maxFaceCentreNonOrtho", 180),
        mesh_,
        metric_.fCtrs(),
        cellCentres,
        metric_.fAreas(),
        identity(mesh_.nFaces()),
        emptyBaffles,
        &wrongFaces
    );

    Foam::polyMeshGeometry::checkFaceWeights
    (
        false,
        qualityControls.lookupOrDefault<scalar>("minFaceWeight", -1, true),
        mesh_,
        cellCentres,
        metric_.fCtrs(),
        metric_.fAreas(),
        identity(mesh_.nFaces()),
        emptyBaffles,
        &wrongFaces
    );

     Foam::polyMeshGeometry::checkFaceTwist
     (
         false,
         qualityControls.lookupOrDefault<scalar>("minTwist", -2, true),
         mesh_,
         cellCentres,
         metric_.fAreas(),
         metric_.fCtrs(),
         metric_.getPoints(),
         identity(mesh_.nFaces()),
         &wrongFaces
     );

    Foam::polyMeshGeometry::checkFaceSkewness
    (
        false,
        qualityControls.lookupOrDefault<scalar>("maxInternalSkewness", 6),
        qualityControls.lookupOrDefault<scalar>("maxBoundarySkewness", 20),
        mesh_,
        cellCentres,
        metric_.fCtrs(),
        metric_.fAreas(),
        identity(mesh_.nFaces()),
        emptyBaffles,
        &wrongFaces
    );

    nWrongFaces_ = wrongFaces.size();
    reduce(nWrongFaces_, sumOp<label>());
    Info<<"----- Number of low quality faces :    "<<nWrongFaces_<<endl;

    isAWrongFace_ = false;

    forAllConstIter(labelHashSet, wrongFaces, i)
    {
        isAWrongFace_[i.key()] = true;
    }
}

void Foam::meshQualityStatus::calculateActivePointSet()
{
    activePoints_ = false;
    //Activate all points that belong to a wrong face
    forAll(mesh_.faces(), fI)
    {
        const labelList& facePoints = mesh_.faces()[fI];
        if (isAWrongFace_[fI])
        {
            forAll(facePoints, pI)
            {
                const label& pL = facePoints[pI];
                activePoints_[pL] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        activePoints_,
        maxEqOp<bool>(),
        false
    );

    //Increment the active points by a cell layer around each point
    //pointCells -> cellPoints

    label inc = dict_.lookupOrDefault("numberOfIncrementations", 10);
    Info<<"numberOfIncrementations.... "<<inc<<endl;
    for (int i=0;i<inc;i++)
    {
        incrementActivePointSet();
    }
}

void Foam::meshQualityStatus::incrementActivePointSet()
{
    boolList activeCells(mesh_.cells().size(), false);

    forAll(activePoints_, pI)
    {
        if (activePoints_[pI])
        {
            forAll(mesh_.pointCells()[pI], cI)
            {
                const label& cL = mesh_.pointCells()[pI][cI];
                activeCells[cL] = true;
            }
        }
    }

    forAll(mesh_.cells(), cI)
    {
        if (activeCells[cI])
        {
            forAll(mesh_.cellPoints()[cI], pI)
            {
                const label& pL = mesh_.cellPoints()[cI][pI];
                activePoints_[pL] = true;
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        activePoints_,
        maxEqOp<bool>(),
        false
    );
}

Foam::tmp<Foam::vectorField> Foam::meshQualityStatus::calculateCorrectCellCentres()
{
    tmp<vectorField> tCellCtrs(new vectorField(mesh_.cells().size()));
    vectorField& cellCtrs = tCellCtrs.ref();
    scalarField cellVols(mesh_.cells().size());
    const vectorField& fAreas = metric_.fAreas();
    const vectorField& fCtrs = metric_.fCtrs();

    forAll(mesh_.cells(), cI)
    {
        cellCtrs[cI] = vector::zero;
        cellVols[cI] = 0;
    }

    const label& nCells = mesh_.cells().size();

    vectorField cEst(nCells, vector::zero);
    labelField nCellFaces(nCells, 0);

    forAll(mesh_.cells(), cI)
    {
        const labelList& cFaces = mesh_.cells()[cI];
        forAll(cFaces, f)
        {
            const label& fL = cFaces[f];
            cEst[cI] += metric_.fCtrs()[fL];
            nCellFaces[cI]++;
        }
    }

    forAll(mesh_.cells(), celli)
    {
        cEst[celli] /= nCellFaces[celli];
    }

    const labelList& own = mesh_.faceOwner();

    forAll(mesh_.cells(), cI)
    {
        const labelList& cFaces = mesh_.cells()[cI];
        forAll(cFaces, f)
        {
            const label& fL = cFaces[f];
            scalar pyr3Vol;
            if (own[fL]==cI)
            {
                // Calculate 3*face-pyramid volume
                pyr3Vol =
                max((fAreas[fL]) & ((fCtrs[fL] - cEst[cI])), VSMALL);
            }
            else
            {
                   pyr3Vol =
                   max((fAreas[fL]) & ((cEst[cI] - fCtrs[fL])), VSMALL);
            }

               // Calculate face-pyramid centre
               vector pc = ((3.0/4.0)*fCtrs[fL] + (1.0/4.0)*cEst[cI]);

               // Accumulate volume-weighted face-pyramid centre
               cellCtrs[cI] += pyr3Vol*pc;

               // Accumulate face-pyramid volume
               cellVols[cI] += pyr3Vol;
        }
    }

    forAll(mesh_.cells(), cI)
    {
        cellCtrs[cI] /= cellVols[cI];
        cellVols[cI] *= (1.0/3.0);
    }
    return tCellCtrs;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::meshQualityStatus::meshQualityStatus
(
    const fvMesh& mesh,
    const meshMetric& metric,
    const dictionary& dict
)
:
    mesh_(mesh),
    metric_(metric),
    dict_(dict),
    nWrongFaces_(0),
    isAWrongFace_(mesh_.faces().size()),
    activePoints_(mesh_.points().size(), true)
{
}

meshQualityStatus::~meshQualityStatus()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::meshQualityStatus::checkMeshQuality()
{
    labelHashSet wrongFaces;
    findWrongFaces(wrongFaces);
    calculateActivePointSet();
}

const boolList& Foam::meshQualityStatus::getActivePoints() const
{
    return activePoints_;
}

const boolList& Foam::meshQualityStatus::getWrongFaces() const
{
    return isAWrongFace_;
}

const label& Foam::meshQualityStatus::numberOfWrongFaces() const
{
    return nWrongFaces_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} /* namespace Foam */
