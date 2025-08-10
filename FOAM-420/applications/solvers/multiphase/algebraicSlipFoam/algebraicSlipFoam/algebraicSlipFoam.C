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
    algebraicSlipFoam

Group
    grpMultiphaseSolvers

Description
    Solver for 2 incompressible fluids using the mixture approach with the
    algebraicSlip-flux approximation for relative motion of the phases.

    Used for simulating any dispersed phase separation problems.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "fvMatrices/solvers/MULES/MULES.H"
#include "algorithms/subCycle/subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "relativeVelocityModel/relativeVelocityModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "finiteVolume/laplacianSchemes/gaussLaplacianScheme/gaussLaplacianScheme.H"
#include "finiteVolume/snGradSchemes/uncorrectedSnGrad/uncorrectedSnGrad.H"

#include "finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/CrankNicolsonDdtScheme/CrankNicolsonDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"

    #include "include/createTime.H"
    #include "include/createMesh.H"

    pimpleControl pimple(mesh);

    #include "cfdTools/general/include/createTimeControls.H"
    #include "createFields.H"
    #include "cfdTools/general/include/createFvOptions.H"
    #include "cfdTools/general/include/initContinuityErrs.H"

    turbulence->validate();

    volScalarField divPhi
    (
        IOobject
        (
            "divPhi",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("divPhi", phiAlphaAdvection.dimensions()/dimVolume, 0.)
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "cfdTools/general/include/readTimeControls.H"
        #include "cfdTools/incompressible/CourantNo.H"
        #include "cfdTools/general/include/setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        UrelModel.correct();
        surfaceScalarField phir(fvc::flux(UrelModel.Urel()/rho)*rho2);

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "alphaControls.H"

            #include "alphaEqnSubCycle.H"

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

            divPhi = fvc::div(phiAlphaAdvection);
            Info<< "min/max(divPhiTot): " << min(divPhi) << " , " << max(divPhi) << endl;
        }

        Uair = U + (1.-alpha1)*rho2/rho * UrelModel.Urel();
        Uwater = U - alpha1*rho1/rho * UrelModel.Urel();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
