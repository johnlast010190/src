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
    (c) 2019-2024 Esi Ltd.

Class
    Foam::fv::fluidEnergySolver

Group
    grpFvOptionsSolvers

Description

    solverObject for enthalpy or energy solution with optional porous solid
    (thermal equilibrium assumed)

SourceFiles
    fluidEnergySolver.C

\*---------------------------------------------------------------------------*/
#ifndef fluidEnergySolver_H
#define fluidEnergySolver_H

#include "solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "primitives/enums/Enum.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class basicThermo;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class fluidEnergySolver Declaration
\*---------------------------------------------------------------------------*/

class fluidEnergySolver
:
    public solverObject
{
protected:

    const solutionControl* solnControlPtr_;

    basicThermo* thermoPtr_;

    //- Does it have MRF?
    bool hasMRF_;

    //- Stabilisation mechanism for low Froude number flows
    Switch lowFroudeStabilisation_;

    //- Coefficient for stabilisation of low Froude number flows (0-1)
    scalar lowFroudeStabilisationCoeff_;

    enum formulationType
    {
        totalEnergyEnergy = 0,
        totalEnergyTemperature,
        boussinesqEnergy,
        boussinesqTemperature
    } formulation_;

    autoPtr<volScalarField> K_;

    //- create K field
    virtual void createK();

    //- Update K field
    const volScalarField& updateK(const volVectorField& U);


protected:

public:

    static const Enum<formulationType> formulationTypeNames_;

//- Runtime type information
    TypeName("fluidEnergySolver");

// Constructors

    // - Construct from components
    fluidEnergySolver
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );


//- Destructor
    virtual ~fluidEnergySolver()
    {}

// Member Functions

    virtual bool initialise() override;

    virtual void getSolveGraph
    (
        wordList& solveNames,
        HashTable<wordList>& requiredDependencies,
        HashTable<wordList>& optionalDependencies,
        HashTable<wordList>& correctorMembers
    ) override;

    virtual bool isFinalCorrector
    (
        const label corrector,
        const word& correctorName
    ) override;

    virtual tmp<fvScalarMatrix> assembleScalarMatrix
    (
        const word& solveName,
        bool& finalSolve,
        word& dictName
    ) override;

    virtual void correct(const word& solveName, const word& regionName)
    override;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
