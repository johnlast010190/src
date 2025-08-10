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
    (c) 2019-2022 Esi Ltd.
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "fixedTemperature.H"
#include "solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(fixedTemperature, 0);
}
}

makeFvSolverOption(fixedTemperature);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::fixedTemperature::fixedTemperature
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    thermoPtr_(nullptr)
{
    thermoPtr_ = &basicThermo::lookupOrCreate(mesh_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::fixedTemperature::initialise()
{
    thermoPtr_->T().oldTime();
    return true;
}


void Foam::fv::fixedTemperature::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& derivedFields,
    HashTable<wordList>& requiredDependencies,
    HashTable<wordList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    word heName = thermoPtr_->he().name();
    solveNames.append({heName});
    derivedFields.insert(heName, {IOobject::groupName("T", phaseName_)});

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::fixedTemperature::correct
(
    const word& solveName,
    const word& regionName
)
{
    basicThermo& thermo = *thermoPtr_;
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
    thermo.correct();
}

// ************************************************************************* //
