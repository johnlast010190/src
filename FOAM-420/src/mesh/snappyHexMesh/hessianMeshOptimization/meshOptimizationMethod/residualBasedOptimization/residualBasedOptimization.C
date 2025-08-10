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

#include "hessianMeshOptimization/meshOptimizationMethod/residualBasedOptimization/residualBasedOptimization.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(residualBasedOptimization, 0);
    addToRunTimeSelectionTable
    (
        meshOptimizationMethod,
        residualBasedOptimization,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void residualBasedOptimization::initialActiveSet()
{
    const labelList& pS = identity(mesh_.nPoints());
    const labelList& fS = identity(mesh_.nFaces());
    const labelList& cS = identity(mesh_.nCells());
    const labelList& eS = identity(mesh_.nEdges());
    activeSet_.update(pS, fS, cS, fS, eS);
}

void residualBasedOptimization::calculateActiveSet()
{
    boolList activePoints(points_.size(), false);

    activePoints[pL_] = true;
    pL_++;
    if (pL_ == points_.size())
    {
        pL_ = 0;
    }
//    scalar relTol = 0.00001;
//    scalar absTol = 0.00000001;
//
//    forAll(mesh_.points(), pI)
//    {
//        scalar disp = mag(points_[pI] - oldPoints_[pI]);
//
//        maxDisp_[pI] = max(maxDisp_[pI], disp);
//
//        if (maxDisp_[pI] > VSMALL)
//        {
//            scalar relDisp = disp/maxDisp_[pI];
//
//            if (relDisp < relTol || disp < absTol)
//            {
//                counter_[pI]++;
//            }
//            else
//            {
//                counter_[pI] = 0;
//            }
//
//            if (counter_[pI] == 3)
//            {
//                activePoints[pI] = false;
//            }
//        }
//    }

//    oldPoints_ = points_;

//    incrementActivePointsSet(activePoints);
//    incrementActivePointsSet(activePoints);

    labelList pointSet = compactList(activePoints);

    //calculate active faces
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
residualBasedOptimization::residualBasedOptimization
(
    const dictionary& dict,
    const fvMesh& mesh,
    const pointField& points
)
:
    meshOptimizationMethod(dict,mesh, points),
    mesh_(mesh),
    points_(points),
    oldPoints_(points),
    maxDisp_(mesh_.points().size(), 0.),
    counter_(mesh_.points().size(), 0),
    pL_(0)
{
//    initialActiveSet();
    calculateActiveSet();
}

residualBasedOptimization::~residualBasedOptimization()
{}

bool residualBasedOptimization::updateActiveSet()
{
    calculateActiveSet();
    return true;
}

bool residualBasedOptimization::converge()
{
//    if (activeSet_.activePoints().size() == 0 && activeSet_.updated())
//    {
//        return true;
//    }
//    else
//    {
//        return false;
//    }
    return false;
}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

} /* namespace Foam */
