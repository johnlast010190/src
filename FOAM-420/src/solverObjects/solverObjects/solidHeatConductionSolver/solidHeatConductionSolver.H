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

Class
    Foam::fv::solidHeatConductionSolver

Group
    grpFvOptionsSolvers

Description

    solverObject for enthalpy or energy solution

SourceFiles
    solidHeatConductionSolver.C

\*---------------------------------------------------------------------------*/
#ifndef solidHeatConductionSolver_H
#define solidHeatConductionSolver_H

#include "solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "solidThermo/solidThermo.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "referenceFrames/coordinateFrame.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "coordinate/systems/coordinateSystem.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class solidHeatConductionSolver Declaration
\*---------------------------------------------------------------------------*/

class solidHeatConductionSolver
:
    public solverObject
{
protected:
    solidThermo *thermoPtr_;
    const solutionControl* solnControlPtr_;

    coordinateFrame* coorFramePtr_;
    autoPtr<volSymmTensorField> aniAlpha_;
    autoPtr<volScalarField> betav_;
    autoPtr<volVectorField> gradT_;

    label nOuterCorrectors_;
    label outerCorr_;
    label nonOrthCorr_;

    // Maximum diffusion number
    scalar maxDi_;

protected:
    void calcDiffusionNumber();
    void updateAniAlpha();

public:

//- Runtime type information
    TypeName("solidHeatConductionSolver");

// Constructors

    // - Construct from components
    solidHeatConductionSolver
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );

//- Destructor
    virtual ~solidHeatConductionSolver()
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

    virtual scalar getMaxTimeStep() override;

    virtual bool isFinalCorrector
    (
        const label corrector,
        const word& correctorName
    ) override;

    virtual void beginIteration
    (
        const label corrector,
        const word& correctorName,
        const bool finalIter
    ) override;

    virtual tmp<fvScalarMatrix> assembleScalarMatrix
    (
        const word& fieldName,
        bool& finalSolve,
        word& solveName
    ) override;

    virtual void correct(const word&, const word&) override;

    virtual void endIteration
    (
        const label corrector,
        const word& correctorName,
        const bool finalIter
    ) override;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
