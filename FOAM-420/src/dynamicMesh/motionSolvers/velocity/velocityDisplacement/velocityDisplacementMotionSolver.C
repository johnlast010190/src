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
    (c) 2016 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "motionSolvers/velocity/velocityDisplacement/velocityDisplacementMotionSolver.H"
#include "motionSolvers/displacement/displacement/displacementMotionSolver.H"
#include "fields/pointPatchFields/basic/fixedValue/fixedValuePointPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityDisplacementMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        velocityDisplacementMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::wordList
Foam::velocityDisplacementMotionSolver::pointDisplacementBoundaryTypes() const
{
    const pointVectorField::Boundary& pmUbf(pointMotionU().boundaryField());

    wordList cmUbf = pmUbf.types();

    forAll(pmUbf, patchI)
    {
        if (isA<fixedValuePointPatchField<vector>>(pmUbf[patchI]))
        {
            cmUbf[patchI] = fixedValuePointPatchField<vector>::typeName;
        }
    }

    return cmUbf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityDisplacementMotionSolver::velocityDisplacementMotionSolver
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    velocityMotionSolver(mesh, dict, typeName),
    displacementMotionSolverPtr_()
{
    pointIOField points0(pointsMotionSolver::points0IO(mesh));

    pointVectorField pointDisplacement
    (
        IOobject
        (
            "pointVelocityDisplacement",
            mesh.time().timeName(),
            mesh
        ),
        pointMotionU().mesh(),
        dimLength,
        pointDisplacementBoundaryTypes()
    );

    pointDisplacement.primitiveFieldRef() = mesh.points() - points0;

    displacementMotionSolverPtr_.set
    (
        dynamic_cast<displacementMotionSolver*>
        (
            displacementMotionSolver::New
            (
                coeffDict().lookup("solver"),
                mesh,
                coeffDict(),
                pointDisplacement,
                points0
            ).ptr()
        )
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::velocityDisplacementMotionSolver::~velocityDisplacementMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::velocityDisplacementMotionSolver::curPoints() const
{
    return displacementMotionSolverPtr_->curPoints();
}


void Foam::velocityDisplacementMotionSolver::solve()
{
    movePoints(mesh().points());

    const scalar deltaT(mesh().time().deltaTValue());

    // Current and old point displacements
    pointVectorField& displacement
    (
        displacementMotionSolverPtr_->pointDisplacement()
    );
    const vectorField displacementOld
    (
        mesh().points() - displacementMotionSolverPtr_->points0()
    );

    // Update the velocity boundary conditions
    pointMotionU().correctBoundaryConditions();

    pointVectorField::Boundary& dispBf = displacement.boundaryFieldRef();

    // Update the displacement boundary conditions
    forAll(pointMotionU().boundaryField(), patchI)
    {
        const pointPatchVectorField& patchField
        (
            pointMotionU().boundaryField()[patchI]
        );

        dispBf[patchI].forceAssign(
            patchField.patchInternalField()*deltaT
          + patchField.patchInternalField(displacementOld)
        );
    }

    // Run the sub-solver
    displacementMotionSolverPtr_->solve();

    // Update the velocity
    pointMotionU().primitiveFieldRef() =
        (displacement.primitiveField() - displacementOld)/deltaT;
}


void Foam::velocityDisplacementMotionSolver::movePoints(const pointField& p)
{
    velocityMotionSolver::movePoints(p);

    displacementMotionSolverPtr_->movePoints(p);
}


void Foam::velocityDisplacementMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    velocityMotionSolver::updateMesh(mpm);

    displacementMotionSolverPtr_->updateMesh(mpm);
}


// ************************************************************************* //
