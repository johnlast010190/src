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
    (c) 2011-2017 OpenFOAM Foundation

Application
    rhoReactingFoam

Group
    grpCombustionSolvers

Description
    Solver for combustion with chemical reactions using density based
    thermodynamics package.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "rhoCombustionModel/rhoCombustionModel/rhoCombustionModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "interpolation/surfaceInterpolation/multivariateSchemes/multivariateScheme/multivariateScheme.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/pressureControl/pressureControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "finiteVolume/fvc/fvcSmooth/fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "cfdTools/general/include/initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "cfdTools/compressible/compressibleCourantNo.H"
        #include "cfdTools/general/include/setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "cfdTools/general/include/readTimeControls.H"

        if (LTS)
        {
            #include "../reactingFoam/setRDeltaT.H"
        }
        else
        {
            #include "cfdTools/compressible/compressibleCourantNo.H"
            #include "cfdTools/general/include/setDeltaT.H"
        }
        #include "../reactingFoam/readSolverControls.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "cfdTools/compressible/rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "../reactingFoam/UEqn.H"
            #include "../reactingFoam/YEqn.H"
            #include "../reactingFoam/EEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
                    #include "../../../compressible/rhoPimpleFoam/rhoPimpleFoam/pcEqn.H"
                }
                else
                {
                    #include "../../../compressible/rhoPimpleFoam/rhoPimpleFoam/pEqn.H"
                }
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }

            fvOptions.correct();
        }

        rho = thermo.rho();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
