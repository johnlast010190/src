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
    (c) 2019 Esi Ltd.

Application
    foamBendingWave

Description
    Solver for bending waves traveling at a surface excited by an external
    force.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "include/faCFD.H"

#include "pressureVibratingWindow/pressureVibratingWindowFvPatchScalarField.H"
#include "fields/faPatchFields/derived/fixedValueZeroGradient/fixedValueZeroGradientFaPatchFields.H"

//#include "algorithms/subCycle/subCycle.H"
#include "db/Time/subCycleTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "createFaMesh.H"
    #include "readVibroacousticProperties.H"

    #include "createFaFields.H"
    #include "createFvFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting vibroacoustic time loop " << nl << endl;

    while (runTime.loop())
    {
        Info<< "Vibration subloops for time " << runTime.timeName() << endl;

        Foam::TimeState subTime = runTime.subCycle(nCycles);
        solverPerformance::debug = false;

        for (int subCycleI=0;subCycleI < nCycles; subCycleI++)
        {
            ++runTime;
            Info<< "  o-Time = " << runTime.value() << endl;

            #include "setExcitation.H"
            #include "wEqn.H"

            if (subCycleI == nCycles-1)
            {
                #include "mapFaToFv.H"
            }
        }

        Info<< "    ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

        solverPerformance::debug = true;
        runTime.endSubCycle();

        Info<< "Interior propagation loop for time " << runTime.timeName() << nl << endl;

        #include "pIntEqn.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
