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
    (c) 2019-2023 Esi Ltd.

Class
    Foam::fv::singlePhaseSegregatedFlowSolver

Group
    grpFvOptionsSolvers

Description

    solverOption for fluid pressure-based U-p solve

SourceFiles
    singlePhaseSegregatedFlowSolver.C

\*---------------------------------------------------------------------------*/
#ifndef singlePhaseSegregatedFlowSolver_H
#define singlePhaseSegregatedFlowSolver_H

#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "rhoThermo/rhoThermo.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "cfdTools/general/pressureControl/pressureControl.H"
#include "flowSolver/flowSolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class singlePhaseSegregatedFlowSolver Declaration
\*---------------------------------------------------------------------------*/

class singlePhaseSegregatedFlowSolver
:
    public flowSolver
{
protected:
    volScalarField* p_;
    autoPtr<volVectorField> U_;
    rhoThermo *thermoPtr_;
    autoPtr<volScalarField> rho_;
    autoPtr<surfaceScalarField> phi_;

    //- Volumetric flux - only used for transonic
    autoPtr<surfaceScalarField> phiv_;

    autoPtr<surfaceVectorField> rhoUf_;
    label pRefCell_;
    scalar pRefValue_;
    bool closedVolume_;
    autoPtr<dimensionedScalar> initialMass_;
    autoPtr<dimensionedScalar> compressibility_;
    autoPtr<dimensionedScalar> rhoMin_;
    autoPtr<dimensionedScalar> rhoMax_;
    compressible::turbulenceModel* turbulence_;
    Switch printContErr_;
    scalar cumulativeContErr_;
    scalar maxCo_;

    autoPtr<Function1<scalar>> maxCoDataPtr_;
    solutionControl* solnControlPtr_;
    bool useGradP_;
    autoPtr<pressureControl> pressureControl_;

    bool transient_;

    //- Is the case buoyant?
    bool buoyant_;

    //- Does it have different density definition for buoyancy term?
    bool distinctBuoyancy_;

    //- Density used for buoyancy term
    autoPtr<volScalarField> bRho_;

    // Indicates that phi was not read from disk at initialisation
    bool recreatePhi_;

    label corr_;
    label nonOrthCorr_;

    tmp<volVectorField> HbyA_;
    tmp<surfaceScalarField> phiHbyA_;
    tmp<volScalarField> rAU_;
    tmp<surfaceScalarField> rhorAUf_;
    tmp<surfaceScalarField> phig_;
    tmp<volScalarField> p0_;
    tmp<volScalarField> divRhoU_;
    tmp<volVectorField> rhoU_;

    // Use tmps for proper reference counting (allows cleanup from outside)
    tmp<fvVectorMatrix> UEqnLHS_;
    tmp<fvScalarMatrix> pDDtEqn_;
    tmp<fvScalarMatrix> pEqn_;

protected:
    void assembleUEqnLHS();
    void assemblepEqn();
    void solveContinuity();
    void printIncoContinuityErrors(const word&);
    void printCompContinuityErrors(const word&);

    //- Return pressure
    const volScalarField& p() const override
    {
        return *p_;
    }

    //- Return thermo model
    const fluidThermo& thermo() const override
    {
        #ifdef FULLDEBUG
            if (!thermoPtr_)
            {
                FatalErrorInFunction
                    << "thermoPtr_ not allocated"
                    << exit(FatalError);
            }
        #endif
        return *thermoPtr_;
    }

    //- Return solver rho
    const volScalarField& rho() const override
    {
        return rho_();
    }

    const volScalarField& buoyantRho() const override;

    const solutionControl& solnControl() const override
    {
        return *solnControlPtr_;
    }

    const surfaceScalarField& phiv() const override
    {
        return phiv_();
    }

    void correctPhi(bool defaultCorrectPhi);

public:

//- Runtime type information
    TypeName("singlePhaseSegregatedFlowSolver");

// Constructors

    // - Construct from components
    singlePhaseSegregatedFlowSolver
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );


//- Destructor
    virtual ~singlePhaseSegregatedFlowSolver()
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
        word& dictName
    ) override;
    virtual tmp<fvVectorMatrix> assembleVectorMatrix
    (
        const word& fieldName,
        bool& finalSolve,
        word& dictName
    ) override;

    virtual void correct(const word&, const word&) override;

    virtual void endIteration
    (
        const label corrector,
        const word& correctorName,
        const bool finalIter
    ) override;

    virtual bool movePoints() override;

    virtual void updateMesh(const mapPolyMesh& mpm) override;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
