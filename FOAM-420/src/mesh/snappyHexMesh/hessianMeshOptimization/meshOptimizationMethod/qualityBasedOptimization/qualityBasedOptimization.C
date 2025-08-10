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


#include "hessianMeshOptimization/meshOptimizationMethod/qualityBasedOptimization/qualityBasedOptimization.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(qualityBasedOptimization, 0);
    addToRunTimeSelectionTable
    (
        meshOptimizationMethod,
        qualityBasedOptimization,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void qualityBasedOptimization::calculateActiveSet()
{

    mqs_.checkMeshQuality();
    const boolList& activePoints = mqs_.getActivePoints();

    labelList pointSet = compactList(activePoints);

//Calculate active faces
    boolList isActiveFace(mesh_.faces().size(), false);
    forAll(mesh_.pointFaces(), pI)
    {
        if (activePoints[pI])
        {
            forAll(mesh_.pointFaces()[pI], fI)
            {
                const label& fL = mesh_.pointFaces()[pI][fI];
                isActiveFace[fL] = true;
            }
        }
    }
    labelList faceSet = compactList(isActiveFace);
    isActiveFace.clear();

//Calculate active cells
    boolList isActiveCell(mesh_.cells().size(), false);
    forAll(mesh_.pointCells(), pI)
    {
        if (activePoints[pI])
        {
            forAll(mesh_.pointCells()[pI], fI)
            {
                const label& fL = mesh_.pointCells()[pI][fI];
                isActiveCell[fL] = true;
            }
        }
    }
    labelList cellSet = compactList(isActiveCell);

//Calculate active faces for cell calculations
    boolList isActiveCellFace(mesh_.faces().size(), false);
    forAll(mesh_.cells(), cI)
    {
        if (isActiveCell[cI])
        {
            forAll(mesh_.cells()[cI], fI)
            {
                const label& fL = mesh_.cells()[cI][fI];
                isActiveCellFace[fL] = true;
            }
        }
    }

    syncTools::syncFaceList
    (
        mesh_,
        isActiveCellFace,
        maxEqOp<bool>()
    );

    labelList faceCellSet = compactList(isActiveCellFace);
    isActiveCellFace.clear();
    isActiveCell.clear();

//Calculate active edges
    boolList isActiveEdge(mesh_.edges().size(), false);
    forAll(mesh_.pointEdges(), pI)
    {
        if (activePoints[pI])
        {
            forAll(mesh_.pointEdges()[pI], fI)
            {
                const label& fL = mesh_.pointEdges()[pI][fI];
                isActiveEdge[fL] = true;
            }
        }
    }
    labelList edgeSet = compactList(isActiveEdge);
    isActiveEdge.clear();

    activeSet_.update(pointSet, faceSet, cellSet, faceCellSet, edgeSet);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
qualityBasedOptimization::qualityBasedOptimization
(
    const dictionary& dict,
    const fvMesh& mesh,
    const meshMetric& metric
)
:
    meshOptimizationMethod(dict,mesh, metric),
    mesh_(mesh),
    mqs_(mesh, metric, dict)
{
    calculateActiveSet();
}

qualityBasedOptimization::~qualityBasedOptimization()
{}

bool qualityBasedOptimization::updateActiveSet()
{
    calculateActiveSet();
    return true;
}

bool qualityBasedOptimization::converge()
{
    if (mqs_.numberOfWrongFaces()!=0)
    {
        return false;
    }

    return true;
}

void qualityBasedOptimization::merge(const labelList& pointSet)
{
    const labelList& localPointSet = activeSet_.activePoints();
    boolList isAnActivePoint(pointSet.size(), false);
    forAll(pointSet, pI)
    {
        forAll(localPointSet, pLI)
        {
            if (pointSet[pI]==localPointSet[pLI])
            {
                isAnActivePoint[pI] = true;
            }
        }
    }
    label counter = 0;
    forAll(pointSet, pI)
    {
        if (isAnActivePoint[pI])
        {
            counter++;
        }
    }
    labelList newPointSet(counter);

    counter = 0;
    forAll(pointSet, pI)
    {
        if (isAnActivePoint[pI])
        {
            newPointSet[counter++] = pointSet[pI];
        }
    }
    activeSet_.updatePointSet(newPointSet);
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

} /* namespace Foam */
