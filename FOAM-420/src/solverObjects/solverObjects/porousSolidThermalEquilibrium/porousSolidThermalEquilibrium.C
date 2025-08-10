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
    (c) 2019-2023 Esi Ltd.
    (c) 2016-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "porousSolidThermalEquilibrium.H"
#include "solverOption/SolverOption.H"
#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/diffusionNumber/diffusionNumber.H"
#include "cfdTools/general/porosityModel/porosityModel/IOporosityModelList.H"
#include "solidThermo/solidThermo.H"
#include "interpolation/surfaceInterpolation/schemes/midPoint/midPoint.H"
#include "finiteVolume/convectionSchemes/gaussConvectionScheme/gaussConvectionScheme.H"
#include "fvSolutionRegistry/fvSolutionRegistry.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(porousSolidThermalEquilibrium, 0);
}
}

makeFvSolverOption(porousSolidThermalEquilibrium);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::porousSolidThermalEquilibrium::porousSolidThermalEquilibrium
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    solverObject(name, obr, dict),
    reactionPtr_(nullptr),
    thermoPtr_(nullptr),
    fluidThermoPtr_(nullptr),
    fluidRegion_(dict.lookup("fluidRegion")),
    fluidDb_(mesh_.lookupObject<fvSolutionRegistry>(fluidRegion_)),
    solnControlPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::porousSolidThermalEquilibrium::read(const dictionary& dict)
{}

bool Foam::fv::porousSolidThermalEquilibrium::initialise()
{
    thermoPtr_ = &basicThermo::lookupOrCreate(obr_);
    fluidThermoPtr_ = &basicThermo::lookupOrCreate(fluidDb_);

    mesh_.schemes().setFluxRequired(fluidThermoPtr_->he().name());
    mesh_.schemes().setFluxRequired(fluidThermoPtr_->T().name());

    solnControlPtr_ = &solutionControl::lookupOrCreate(mesh_, obr_);

    return true;
}


void Foam::fv::porousSolidThermalEquilibrium::getSolveGraph
(
    wordList& solveNames,
    HashTable<wordList>& derivedSolves,
    HashTable<solveList>& requiredDependencies,
    HashTable<solveList>& optionalDependencies,
    HashTable<wordList>& correctorMembers
)
{
    solveNames.append(thermoPtr_->he().name());
    optionalDependencies.insert
    (
        thermoPtr_->he().name(),
        {solveID(fluidThermoPtr_->he().name(), fluidRegion_)}
    );

    correctorMembers.insert(solverObject::outerCorrectorName, solveNames);
}


void Foam::fv::porousSolidThermalEquilibrium::correct
(
    const word& solveName,
    const word& regionName
)
{
    basicThermo& thermo = *thermoPtr_;
    // Take the temperature from the fluid, make our he such as to produce
    // that T, and recalc
    thermo.he().forceAssign(thermo.he(thermo.p(), fluidThermoPtr_->T()));
    // The fixed-value patches (incl. calculated) won't be calculated in
    // correct; set them
    forAll(thermo.T().boundaryField(), patchi)
    {
        if (thermo.T().boundaryField()[patchi].fixesValue())
        {
            thermo.T().boundaryFieldRef()[patchi] =
                fluidThermoPtr_->T().boundaryField()[patchi];
        }
    }
    thermo.correct();
}


void Foam::fv::porousSolidThermalEquilibrium::getSourceGraph
(
    solveList& fields, SolveTable<solveList>& sourceDependencies
)
{
    // We act as a source for the fluid thermal equation, to add our
    // corrections to the effective density and diffusivity
    fields = {solveID(fluidThermoPtr_->he().name(), fluidRegion_)};

    // Forward any source terms, and hence their dependencies, from solid to
    // fluid equation
    SolveTable<solveList> dependencies;
    fvOptions().addSourceDependencies(dependencies);

    sourceDependencies.insert
    (
        fields[0],
        dependencies.lookup
        (
            solveID(thermoPtr_->he().name(), regionName_), solveList()
        )
    );
}


void Foam::fv::porousSolidThermalEquilibrium::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    // Adding terms to the fluid equation to correct the effective density
    // and diffusivity

    scalar eps;
    fluidThermoPtr_->lookup("porosity") >> eps;

    if (!isA<simpleControl>(*solnControlPtr_))
    {
        if (!porousRhoEffCorr_.valid())
        {
            porousRhoEffCorr_.set
            (
                new volScalarField
                (
                    "porousRhoEffCorr",
                    eqn.psi().db(),
                    mesh_,
                    dimDensity
                )
            );
        }

        // Being added to RHS, so signs reversed
        porousRhoEffCorr_() =
            (
                (1-eps)/eps
                *thermoPtr_->rho()*thermoPtr_->he()/fluidThermoPtr_->he()
            );
        eqn -=
            fvm::ddt
            (
                porousRhoEffCorr_(),
                fluidThermoPtr_->he(),
                "ddt("+rho.name()+","+fluidThermoPtr_->he().name()+")"
            );
    }

    // Add laplacian(kappa, T) to RHS

    word scheme = "laplacian(kappa";
    if (thermoPtr_->alpha().group() != word::null)
    {
        scheme += "." + thermoPtr_->alpha().group();
    }
    scheme += "," + thermoPtr_->T().name() + ")";

    tmp<volScalarField> kappaEff((1-eps)/eps*thermoPtr_->kappa());

    eqn +=
        changeVariable
        (
            fvm::laplacian(kappaEff, fluidThermoPtr_->T(), scheme),
            1/fluidThermoPtr_->Cpv(),
            fluidThermoPtr_->he()
        );

    // Any source terms that were meant for the solid equation get added
    // to the fluid equation multiplied by (1-eps)/eps
    eqn +=
        (1-eps)/eps
       *fvOptions()
        (
            thermoPtr_->rho(),
            fluidThermoPtr_->he(),
            thermoPtr_->he().name()
        );
}

// ************************************************************************* //
