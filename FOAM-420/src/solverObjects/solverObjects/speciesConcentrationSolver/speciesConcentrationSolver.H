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
    Foam::fv::speciesConcentrationSolver

Group
    grpFvOptionsSolvers

Description
    solverObject for fluid species concentration transport, using Fick's law
    for the diffusive component.

    Laminar species diffusivity defaults to using Reynolds analogy, with mass
    diffusivity (in units of [Mass]/[Length]/[Time]) equal to the effective
    viscosity (i.e. a Schmidt number of 1).

    Diffusivity of selected species can be overridden by specifying
    values of concentration diffusivity (as a Function1, with optional
    temperature dependence) in units of [Length]^2/[Time] in the
    materialProperties dictionary, e.g.

    speciesDiffusivity
    {
        CO2 constant 0.01;
    }

    Turbulent diffusivity is added as supplied by the turbulence model. Usually
    this dominates over the laminar diffusivity and is equal to mut/Sct where
    mut is the turbulent viscosity and Sct is the turbulent Schmidt number,
    which can be specified in the turbulence coefficients subdictionary and
    typically defaults to 0.7.

SourceFiles
    speciesConcentrationSolver.C

\*---------------------------------------------------------------------------*/
#ifndef speciesConcentrationSolver_H
#define speciesConcentrationSolver_H

#include "solverObject/solverObject.H"
#include "fields/UniformDimensionedFields/uniformDimensionedFields.H"
#include "solidThermo/solidThermo.H"
#include "cfdTools/general/solutionControl/solutionControl/solutionControl.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiationModels/radiationModel/radiationModel.H"
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
namespace fv
{

/*---------------------------------------------------------------------------*\
                    Class speciesConcentrationSolver Declaration
\*---------------------------------------------------------------------------*/

class speciesConcentrationSolver
:
    public solverObject
{
protected:
    const solutionControl* solnControlPtr_;

    //- If phasic - the phasic thermo
    rhoThermo* thermoPtr_;
    //- If phasic - the top-level thermo
    rhoThermo* globalThermoPtr_;
    const compressible::turbulenceModel* turbulencePtr_;

    basicSpecieMixture* compositionPtr_;
    label inertIndex_;

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable
        multivariateConvectionFields_;

    label nOuterCorrectors_;
    label outerCorr_;
    label nonOrthCorr_;

    // Optionally override species diffusivity with a manually specified value
    PtrList<Function1<scalar>> Di_;

public:

//- Runtime type information
    TypeName("speciesConcentrationSolver");

// Constructors

    // - Construct from components
    speciesConcentrationSolver
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict
    );

//- Destructor
    virtual ~speciesConcentrationSolver()
    {}

// Member Functions

    virtual void read(const dictionary& dict) override;

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

    // Source term contributions to other equations

        virtual void getSourceGraph
        (
            wordList& fields,
            HashTable<wordList>& sourceDependencies
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
