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
    (c) 2013-2015 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "externalDisplacementMeshMover/displacementMeshMoverMotionSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "regionSplit/localPointRegion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(displacementMeshMoverMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        displacementMeshMoverMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::displacementMeshMoverMotionSolver::displacementMeshMoverMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    displacementMotionSolver(mesh, dict, typeName) // read pointDisplacement
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::displacementMeshMoverMotionSolver::
~displacementMeshMoverMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::externalDisplacementMeshMover&
Foam::displacementMeshMoverMotionSolver::meshMover() const
{
    if (!meshMoverPtr_.valid())
    {
        const word moverType(coeffDict().lookup("meshMover"));

        meshMoverPtr_ = externalDisplacementMeshMover::New
        (
            moverType,
            coeffDict().optionalSubDict(moverType + "Coeffs"),
            localPointRegion::findDuplicateFacePairs(mesh()),
            pointDisplacement_
        );
    }
    return meshMoverPtr_();
}


Foam::tmp<Foam::pointField>
Foam::displacementMeshMoverMotionSolver::curPoints() const
{
    // Return actual points. Cannot do a reference since complains about
    // assignment to self in polyMesh::movePoints
    return tmp<pointField>(new pointField(mesh().points()));
}


void Foam::displacementMeshMoverMotionSolver::solve()
{
    // The points have moved so before calculation update
    // the mesh and motionSolver accordingly
    movePoints(mesh().points());

    // Update any point motion bcs (e.g. timevarying)
    pointDisplacement().boundaryFieldRef().updateCoeffs();

    label nAllowableErrors = 0;
    labelList checkFaces(identity(mesh().nFaces()));

    // Points on a convex edge.
    PackedList<1> isConvexEdgePoint(mesh().nPoints(), 0);
    // Points on a concave edge.
    PackedList<1> isConcaveEdgePoint(mesh().nPoints(), 0);

    meshMover().move
    (
        coeffDict().optionalSubDict(meshMover().type() + "Coeffs"),
        nAllowableErrors,
        isConvexEdgePoint,
        isConcaveEdgePoint,
        checkFaces
    );

    // This will have updated the mesh and implicitly the pointDisplacement
    pointDisplacement().correctBoundaryConditions();
}


void Foam::displacementMeshMoverMotionSolver::movePoints(const pointField& p)
{
    displacementMotionSolver::movePoints(p);

    // Update meshMover for new geometry
    if (meshMoverPtr_.valid())
    {
        meshMover().movePoints(p);
    }
}


void Foam::displacementMeshMoverMotionSolver::updateMesh
(
    const mapPolyMesh& map
)
{
    displacementMotionSolver::updateMesh(map);

    // Update meshMover for new topology
    meshMoverPtr_.clear();
}


// ************************************************************************* //
