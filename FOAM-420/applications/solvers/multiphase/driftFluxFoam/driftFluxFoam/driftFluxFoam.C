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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

Application
    driftFluxFoam

Group
    grpMultiphaseSolvers

Description
    Solver for 2 incompressible fluids using the mixture approach with the
    drift-flux approximation for relative motion of the phases.

    Used for simulating the settling of the dispersed phase and other similar
    separation problems.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "fvMatrices/solvers/MULES/CMULES.H"
#include "algorithms/subCycle/subCycle.H"
#include "incompressibleTwoPhaseInteractingMixture.H"
#include "relativeVelocityModel/relativeVelocityModel.H"
#include "turbulenceModel.H"
#include "CompressibleTurbulenceModel/CompressibleTurbulenceModel.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.H"
#include "finiteVolume/snGradSchemes/uncorrectedSnGrad/uncorrectedSnGrad.H"

#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "interpolation/surfaceInterpolation/schemes/downwind/downwind.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "createFields.H"
    #include "cfdTools/general/include/initContinuityErrs.H"

    volScalarField& alpha2(mixture.alpha2());
    const dimensionedScalar& rho1 = mixture.rhod();
    const dimensionedScalar& rho2 = mixture.rhoc();
    relativeVelocityModel& UdmModel(UdmModelPtr());

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "cfdTools/general/include/readTimeControls.H"
        #include "cfdTools/incompressible/CourantNo.H"
        #include "cfdTools/general/include/setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"

            UdmModel.correct();

            #include "alphaEqnSubCycle.H"

            mixture.correct();

            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            fvOptions.correct();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
