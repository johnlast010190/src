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
    (c) 2021-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMeshSolver.H"
#include "solverOption/SolverOption.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/CourantNumber/CourantNumber.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "referenceFrames/coordinateFrame.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(fvMeshSolver, 0);
    }
}

makeFvSolverOption(fvMeshSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fvMeshSolver::fvMeshSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::fvMeshSolver::initialise()
{
    // Ensure delta coeffs are calculated to prevent the calculation being
    // triggered in neighbouring regions from coupled BC,
    // which can cause synchronisation problems with thread teams
    mesh().deltaCoeffs();

    // Opt out if mesh is not dynamic
    return (isA<dynamicFvMesh>(mesh_) || mesh_.hasChangers());
}


void Foam::fv::fvMeshSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames = {"fvMesh"};
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


Foam::scalar Foam::fv::fvMeshSolver::getMaxTimeStep()
{
    const pimpleControl* pimplePtr =
        mesh_.lookupObjectPtr<pimpleControl>(solutionControl::typeName);
    if
    (
        pimplePtr
     && mesh_.changing()
     && pimplePtr->dict().lookupOrDefault("checkMeshCourantNo", false)
    )
    {
        scalar maxMeshCo =
            CourantNumber
            (
                mesh_,
                mesh_.time().deltaTValue(),
                geometricOneField(),
                mesh_.phi()
            );
        // Try region-specific maxMeshCo, otherwise revert to global
        scalar maxMaxCo(GREAT);
        if (pimplePtr->dict().found("maxMeshCo"))
        {
            pimplePtr->dict().lookup("maxMeshCo") >> maxMaxCo;
        }
        else
        {
            maxMaxCo =
                mesh_.time().controlDict().lookupOrDefault
                (
                    "maxMeshCo",
                    GREAT
                );
        }
        return maxMaxCo/maxMeshCo*mesh_.time().deltaTValue();
    }
    else
    {
        return GREAT;
    }
}


void Foam::fv::fvMeshSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    //- Reset update flags
    coordinateFrame::resetUpdates(mesh_.thisDb());

    //- Update reference frames
    coordinateFrame::updateStates(mesh_.thisDb());

    // Update dynamic mesh following old structure
    if (isA<dynamicFvMesh>(mesh_))
    {
        // Do any mesh changes
        const_cast<dynamicFvMesh&>
        (
            dynamicCast<const dynamicFvMesh>(mesh_)
        ).update();
    }

    // Update mesh with changers
    if (mesh_.hasChangers())
    {
        mesh_.update();
    }

    // Ensure delta coeffs are calculated to prevent the calculation being
    // triggered in neighbouring regions from coupled BC,
    // which can cause synchronisation problems with thread teams
    mesh().deltaCoeffs();
}


// ************************************************************************* //
