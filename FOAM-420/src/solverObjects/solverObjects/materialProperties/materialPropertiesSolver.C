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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "materialPropertiesSolver.H"
#include "fluidEnergySolver/fluidEnergySolver.H"
#include "general/referenceFields/referenceFields.H"
#include "foamSolve.H"


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(materialPropertiesSolver, 0);
}
}

makeFvSolverOption(materialPropertiesSolver);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::materialPropertiesSolver::materialPropertiesSolver
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    thermoPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::materialPropertiesSolver::initialise()
{
    // Don't create thermo; check for an existing one and unload this solver
    // object if not found. This allows it to work with solvers which don't use
    // a thermo object
    thermoPtr_ = basicThermo::lookupPtr(obr_);
    return thermoPtr_;
}


void Foam::fv::materialPropertiesSolver::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& derivedFields,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append("materialProperties");
    if (thermoPtr_->solvesTFromhe())
    {
        derivedFields.insert(solveNames[0], {thermoPtr_->T().name()});
    }
    else
    {
        derivedFields.insert(solveNames[0], {thermoPtr_->he().name()});
    }

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
    // Dependencies are added in the relevant solver objects
}


void Foam::fv::materialPropertiesSolver::correct
(
    const word& solveName,
    const word& regionName
)
{
    // Force the temperature to remain constant if there is no solver object
    // for energy
    basicThermo& thermo(*thermoPtr_);
    if
    (
        !this->solveController().fieldToSolveName().found
        (
            solveID(thermo.he().name(), regionName)
        )
    )
    {
        // Take the old temperature, make our he such as to produce
        // that T, and recalc
        thermo.he().forceAssign(thermo.he(thermo.p(), thermo.T().oldTime()));
        // The fixed-value patches (incl. calculated) won't be calculated in
        // correct; set them
        forAll(thermo.T().boundaryField(), patchi)
        {
            if (thermo.T().boundaryField()[patchi].fixesValue())
            {
                thermo.T().boundaryFieldRef()[patchi] =
                    thermo.T().oldTime().boundaryField()[patchi];
            }
        }
    }

    thermo.correct();
}


// ************************************************************************* //
