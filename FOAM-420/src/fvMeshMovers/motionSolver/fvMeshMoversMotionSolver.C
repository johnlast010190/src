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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "motionSolver/fvMeshMoversMotionSolver.H"
#include "motionSolvers/motionSolver/motionSolver.H"
#include "fields/volFields/volFields.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fvMeshMovers
{
    defineTypeNameAndDebug(motionSolver, 0);
    addToRunTimeSelectionTable(fvMeshMover, motionSolver, fvMesh);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshMovers::motionSolver::motionSolver(fvMesh& mesh)
:
    fvMeshMover(mesh),
    motionPtr_(Foam::motionSolver::New(mesh, dict()))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMeshMovers::motionSolver::~motionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::motionSolver& Foam::fvMeshMovers::motionSolver::motion() const
{
    return motionPtr_();
}


bool Foam::fvMeshMovers::motionSolver::update()
{
    mesh().movePoints(motionPtr_->newPoints());

    if (mesh().foundObject<volVectorField>("U"))
    {
        mesh().lookupObjectRef<volVectorField>
        (
            "U"
        ).correctBoundaryConditions();
    }

    return true;
}


void Foam::fvMeshMovers::motionSolver::updateMesh(const mapPolyMesh& mpm)
{
    motionPtr_->updateMesh(mpm);
}


bool Foam::fvMeshMovers::motionSolver::write(const bool write) const
{
    if (write)
    {
        return motionPtr_->write();
    }
    else
    {
        return true;
    }
}


// ************************************************************************* //
