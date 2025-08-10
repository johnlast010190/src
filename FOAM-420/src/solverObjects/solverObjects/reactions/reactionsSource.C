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
    (c) 2019-2021 Esi Ltd.
    (c) 2013-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "reactionsSource.H"
#include "solverOption/SolverOption.H"
#include "transportModel/transportModel.H"
#include "basicThermo/basicThermo.H"
#include "solidThermo/solidThermo.H"
#include "fvMatrices/fvMatrices.H"
#include "finiteVolume/fvm/fvmSup.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(reactionsSource, 0);
}
}

makeFvSolverOption(reactionsSource);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::reactionsSource::reactionsSource
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict)
{
    Info<< "Creating reaction model\n" << endl;
    reaction_ =
        combustionModels::rhoCombustionModel::New
        (
            obr,
            combustionModels::rhoCombustionModel::combustionPropertiesName,
            phaseName_
        );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::reactionsSource::read(const dictionary& dict)
{}

bool Foam::fv::reactionsSource::initialise()
{
    return true;
}


void Foam::fv::reactionsSource::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    // Solves
    solveNames = {"reactions"};

    // Correctors
    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::reactionsSource::correct
(
    const word& solveName,
    const word& regionName
)
{
    reaction_->correct();
}


void Foam::fv::reactionsSource::getSourceGraph
(
    wordList& fields,
    HashTable<wordList>& sourceDependencies
)
{
    DynamicList<word> fn;
    fn.append(reaction_->thermo().he().name());
    auto& Y = reaction_->thermo().composition().Y();
    word inertSpecie =
        reaction_->thermo().lookupOrDefault("inertSpecie", word::null);
    forAll(Y, Yi)
    {
        if (reaction_->thermo().composition().species()[Yi] != inertSpecie)
        {
            fn.append(Y[Yi].name());
        }
    }
    fields.transfer(fn);

    forAll(fields, i)
    {
        sourceDependencies.insert(fields[i], {"reactions"});
    }
}


void Foam::fv::reactionsSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    bool superficial =
        reaction_->thermo().lookupOrDefault
        (
            "superficialReactionRate", false
        );
    scalar eps = 1;
    if (superficial)
    {
        eps = reaction_->thermo().lookupOrDefault("porosity", scalar(1));
    }

    if (fieldI > 0)
    {
        volScalarField& Yi = const_cast<volScalarField&>(eqn.psi());
        if (superficial)
        {
            eqn += 1/eps*reaction_->R(Yi);
        }
        else
        {
            eqn += reaction_->R(Yi);
        }
    }
    else
    {
        if (superficial)
        {
            eqn += reaction_->Qdot()/eps;
        }
        else
        {
            eqn += reaction_->Qdot();
        }
    }
}


// ************************************************************************* //
