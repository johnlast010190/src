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
    rhoPorousSimpleFoam

Group
    grpCompressibleSolvers

Description
    Steady-state solver for turbulent flow of compressible fluids, with
    implicit or explicit porosity treatment and optional sources.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "fluidThermo/fluidThermo.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "radiation/radiationSolver.H"
#include "cfdTools/general/solutionControl/simpleControl/simpleControl.H"
#include "cfdTools/general/pressureControl/pressureControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/porosityModel/porosityModel/IOporosityModelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "createFields.H"
    #include "createZones.H"
    #include "cfdTools/general/include/initContinuityErrs.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "readSolverControls.H"

        // Pressure-velocity SIMPLE corrector
        if (runTime.timeIndex()%hSubCycles == 0)
        {
            // Pressure-velocity SIMPLE corrector
            {
                #include "UEqn.H"
                #include "EEqn.H"
                #include "pEqn.H"
            }

            turbulence->correct();
        }
        else
        {
            #include "EEqn.H"
        }

        fvOptions.correct();

        runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
