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
    dnsFoam

Group
    grpDNSSolvers

Description
    Direct numerical simulation solver for boxes of isotropic turbulence.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "Kmesh/Kmesh.H"
#include "processes/UOprocess/UOprocess.H"
#include "fft/fft.H"
#include "fft/calcEk.H"
#include "graph/graph.H"
#include "cfdTools/general/solutionControl/pisoControl/pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMeshNoClear.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "createFields.H"
    #include "cfdTools/general/include/initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< nl << "Starting time loop" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        force.primitiveFieldRef() = ReImSum
        (
            fft::reverseTransform
            (
                K/(mag(K) + 1.0e-6) ^ forceGen.newField(), K.nn()
            )
        );

        #include "globalProperties.H"

        fvVectorMatrix UEqn
        (
            fvm::ddt(U)
          + fvm::div(phi, U)
          - fvm::laplacian(nu, U)
         ==
            force
        );

        solve(UEqn == -fvc::grad(p));


        // --- PISO loop
        while (piso.correct())
        {
            volScalarField rAU(1.0/UEqn.A());
            surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
            volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
            surfaceScalarField phiHbyA
            (
                "phiHbyA",
                fvc::flux(HbyA)
              + rAUf*fvc::ddtCorr(U, phi)
            );

            // Update the pressure BCs to ensure flux consistency
            constrainPressure(p, U, phiHbyA, rAUf);

            fvScalarMatrix pEqn
            (
                fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
            );

            pEqn.solve(mesh.solution().solver(p.select(piso.finalInnerIter())));

            phi = phiHbyA - pEqn.flux();

            #include "cfdTools/incompressible/continuityErrs.H"

            U = HbyA - rAU*fvc::grad(p);
            U.correctBoundaryConditions();
        }

        runTime.write();

        if (runTime.writeTime())
        {
            calcEk(U, K).write
            (
                runTime.path()/"graphs"/runTime.timeName(),
                "Ek",
                runTime.graphFormat()
            );
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
