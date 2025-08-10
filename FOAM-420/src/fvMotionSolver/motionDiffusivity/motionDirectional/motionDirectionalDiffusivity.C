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
    (c) 2011-2012 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "motionDiffusivity/motionDirectional/motionDirectionalDiffusivity.H"
#include "interpolation/surfaceInterpolation/surfaceInterpolation/surfaceInterpolate.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "dynamicMotionSolverFvMesh/dynamicMotionSolverFvMesh.H"
#include "motionSolver/fvMeshMoversMotionSolver.H"
#include "motionSolvers/velocity/velocityMotionSolver.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionDirectionalDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        motionDirectionalDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionDirectionalDiffusivity::motionDirectionalDiffusivity
(
    const fvMesh& mesh,
    Istream& mdData
)
:
    uniformDiffusivity(mesh, mdData),
    diffusivityVector_(mdData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionDirectionalDiffusivity::~motionDirectionalDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::motionDirectionalDiffusivity::correct()
{
    static bool first = true;

    if (!first)
    {
        const volVectorField& cellMotionU =
            mesh().lookupObject<volVectorField>("cellMotionU");

        volVectorField D
        (
            IOobject
            (
                "D",
                mesh().time().timeName(),
                mesh()
            ),
            diffusivityVector_.y()*vector::one
          + (diffusivityVector_.x() - diffusivityVector_.y())*cellMotionU
           /(mag(cellMotionU) + dimensionedScalar("small", dimVelocity, SMALL)),
            zeroGradientFvPatchVectorField::typeName
        );
        D.correctBoundaryConditions();

        const surfaceVectorField n(mesh().Sf()/mesh().magSf());
        faceDiffusivity_.forceAssign((n & cmptMultiply(fvc::interpolate(D), n)));
    }
    else
    {
        first = false;

        const polyMesh& pMesh = mesh().boundaryMesh().mesh();

        // Get the motionSolver from the dynamic mesh or from the
        // mesh mover (if valid)
        const motionSolver& motion =
            mesh().hasChangers()
          ? refCast<const fvMeshMovers::motionSolver>(mesh().mover()).motion()
          : refCast<const dynamicMotionSolverFvMesh>(pMesh).motion();

        const velocityMotionSolver& mSolver =
            refCast<const velocityMotionSolver>(motion);

        const_cast<velocityMotionSolver&>(mSolver).solve();
        correct();
    }
}


// ************************************************************************* //
