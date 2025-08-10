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
    (c) 2023-2024 Esi Ltd.

Class
    Foam::fv::MULESVolumeFractionSolver

Group
    grpFvOptionsSolvers

Description
    Solver object using MULES to transport the volume fraction (two- or multi-
    phase)

SourceFiles
    MULESVolumeFractionSolver.C

\*---------------------------------------------------------------------------*/
#ifndef MULESVolumeFractionSolver_H
#define MULESVolumeFractionSolver_H

#include "solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "solidThermo/solidThermo.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
#include "coordinate/systems/coordinateSystem.H"
#include "mixtures/speciesMassFractions/speciesMassFractions.H"
#include "chemistryModel/basicChemistryModel/basicChemistryModel.H"
#include "solidThermo/solidThermo.H"
#include "rhoCombustionModel/rhoCombustionModel/rhoCombustionModel.H"
#include "chemistryModel/rhoChemistryModel/rhoChemistryModel.H"
#include "interpolation/surfaceInterpolation/multivariateSchemes/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationScheme.H"
#include "mixtures/basicSpecieMixture/basicSpecieMixture.H"

#include "materialModels/materialTables/materialTables.H"
#include "materialModels/baseModels/materialModels.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class multiphaseThermo;
class solutionControl;

namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class MULESVolumeFractionSolver Declaration
\*---------------------------------------------------------------------------*/

class MULESVolumeFractionSolver
:
    public solverObject
{
protected:
    multiphaseThermo* thermoPtr_;

    //- Number of alpha-equation sub-cycles
    label nAlphaSubCycles_;

    //- Interface compression coefficient
    scalar cAlpha_;

    //- Number of alpha correctors
    label nAlphaCorr_;

    //- Implicit MULES
    bool MULESCorr_;

    //- Whether to apply previous iteration's corrector
    //  Apply the compression correction from the previous iteration
    //  Improves efficiency for steady-simulations but can only be applied
    //  once the alpha field is reasonably steady, i.e. fully developed
    bool alphaApplyPrevCorr_;

    //- Index of the phase which isn't solved for
    label passiveIndex_;

    //- Pointer to volume flux
    const surfaceScalarField* phivPtr_;

    //- Pointer to mass flux
    surfaceScalarField* phiPtr_;

    const volVectorField *UPtr_;

    //- Maximum interface Courant number
    scalar maxAlphaCo_;
    autoPtr<Function1<scalar>> maxAlphaCoDataPtr_;

    solutionControl* solnControlPtr_;

    //- Whether alpha scheme is Crank-Nicolson
    bool crankNicolson_;

    //- Crank-Nicolson off-centering coefficient
    scalar cnCoeff_;

    //- Phase fluxes. Only persistent if using Crank Nicolson
    PtrList<surfaceScalarField> alphaPhis_;

    //- Previous iteration correctors
    PtrList<surfaceScalarField> alphaPhiCorr0_;

    //- Local timestepping
    bool LTS_;

    //- Local deltaT if local timestepping
    autoPtr<volScalarField> rDeltaT_;

    void updateTimeSchemeInfo();

    void calcSuSp
    (
        const label phasei,
        autoPtr<volScalarField::Internal>& Su,
        autoPtr<volScalarField::Internal>& Sp
    );

    virtual void solveAlphas();

public:

//- Runtime type information
    TypeName("MULESVolumeFractionSolver");

// Constructors

    // - Construct from components
    MULESVolumeFractionSolver
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );

//- Destructor
    virtual ~MULESVolumeFractionSolver()
    {}

// Member Functions

    virtual void read(const dictionary& dict) override;

    virtual bool initialise() override;

    virtual void getSolveGraph
    (
        wordList& solveNames,
        HashTable<wordList>& derivedFields,
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

    virtual void correct(const word&, const word&) override;

    virtual void endIteration
    (
        const label corrector,
        const word& correctorName,
        const bool finalIter
    ) override;

    virtual void updateMesh(const mapPolyMesh& mpm) override;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
