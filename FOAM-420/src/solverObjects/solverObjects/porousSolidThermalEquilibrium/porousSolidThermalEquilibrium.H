/*---------------------------------------------------------------------------* \
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
    (c) 2016-2017 OpenCFD Ltd.

Class
    Foam::fv::porousSolidThermalEquilibrium

Group
    grpFvOptionsSolvers

Description
    Calculates the solid energy assuming thermal equilibrium in a porous
    matrix. Adds required source terms to fluid energy equation to simulate
    effective properties. Only formally correct where porosity is uniform.

SourceFiles
    porousSolidThermalEquilibrium.C

SeeAlso
    Foam::fv::fluidEnergySolver

\*---------------------------------------------------------------------------*/

#ifndef porousSolidThermalEquilibrium_H
#define porousSolidThermalEquilibrium_H

#include "solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "rhoCombustionModel/rhoCombustionModel/rhoCombustionModel.H"
#include "basicThermo/basicThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvSolutionRegistry;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class porousSolidThermalEquilibrium Declaration
\*---------------------------------------------------------------------------*/

class porousSolidThermalEquilibrium
:
    public solverObject
{

protected:

    // The fluid reaction
    combustionModels::rhoCombustionModel* reactionPtr_;

    // The thermos
    basicThermo* thermoPtr_;
    basicThermo* fluidThermoPtr_;

    word fluidRegion_;
    const fvSolutionRegistry& fluidDb_;

    // Correction to the effective density - stored for time derivative
    autoPtr<volScalarField> porousRhoEffCorr_;

    solutionControl* solnControlPtr_;

public:

//- Runtime type information
    TypeName("porousSolidThermalEquilibrium");

// Constructors

    // - Construct from components
    porousSolidThermalEquilibrium
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );

    porousSolidThermalEquilibrium(porousSolidThermalEquilibrium&) = delete;

    void operator=(porousSolidThermalEquilibrium&) = delete;

//- Destructor
    virtual ~porousSolidThermalEquilibrium()
    {}

// Member Functions

    // Solver functions

        virtual void read(const dictionary& dict) override;

        virtual bool initialise() override;

        virtual void getSolveGraph
        (
            wordList& solveNames,
            HashTable<wordList>& derivedSolves,
            HashTable<solveList>& requiredDependencies,
            HashTable<solveList>& optionalDependencies,
            HashTable<wordList>& correctorMembers
        ) override;

        virtual void correct(const word&, const word&) override;

    // Source term functions

        virtual void getSourceGraph
        (
            solveList& fields,
            SolveTable<solveList>& sourceDependencies
        ) override;

        //- Add contribution to compressible equation
        virtual void addSup
        (
            const volScalarField& rho,
            fvMatrix<scalar>& eqn,
            const label fieldI
        ) override;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
