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

Application
    multiphaseEulerFoam

Group
    grpMultiphaseSolvers

Description
    Solver for a system of many compressible fluid phases including
    heat-transfer.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "multiphaseSystem.H"
#include "phaseModel/phaseModel.H"
#include "interfacialModels/dragModels/dragModel/dragModel.H"
#include "interfacialModels/heatTransferModels/heatTransferModel/heatTransferModel.H"
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"
#include "fields/fvsPatchFields/basic/fixedValue/fixedValueFvsPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "createFields.H"
    #include "cfdTools/general/include/initContinuityErrs.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "correctPhi.H"
    #include "CourantNo.H"
    #include "cfdTools/general/include/setInitialDeltaT.H"

    scalar slamDampCoeff
    (
        fluid.lookupOrDefault<scalar>("slamDampCoeff", 1)
    );

    dimensionedScalar maxSlamVelocity
    (
        "maxSlamVelocity",
        dimVelocity,
        fluid.lookupOrDefault<scalar>("maxSlamVelocity", GREAT)
    );

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "cfdTools/general/include/readTimeControls.H"
        #include "CourantNo.H"
        #include "cfdTools/general/include/setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            turbulence->correct();
            fluid.solve();
            rho = fluid.rho();
            #include "zonePhaseVolumes.H"

            //#include "TEqns.H"
            #include "UEqns.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            #include "DDtU.H"

            fvOptions.correct();
        }

        runTime.write();

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
