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
    shallowWaterFoam

Group
    grpIncompressibleSolvers

Description
    Transient solver for inviscid shallow-water equations with rotation.

    If the geometry is 3D then it is assumed to be one layers of cells and the
    component of the velocity normal to gravity is removed.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "db/functionObjects/functionObjectList/postProcess.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/solutionControl/createControl.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "\n Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            surfaceScalarField phiv("phiv", phi/fvc::interpolate(h));

            fvVectorMatrix hUEqn
            (
                fvm::ddt(hU)
              + fvm::div(phiv, hU)
            );

            hUEqn.relax();

            if (pimple.momentumPredictor())
            {
                if (rotating)
                {
                    solve(hUEqn + (F ^ hU) == -magg*h*fvc::grad(h + h0));
                }
                else
                {
                    solve(hUEqn == -magg*h*fvc::grad(h + h0));
                }

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                    hU.correctBoundaryConditions();
                }
            }

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                volScalarField rAU(1.0/hUEqn.A());
                surfaceScalarField ghrAUf(magg*fvc::interpolate(h*rAU));

                surfaceScalarField phih0(ghrAUf*mesh.magSf()*fvc::snGrad(h0));

                volVectorField HbyA("HbyA", hU);
                if (rotating)
                {
                    HbyA = rAU*(hUEqn.H() - (F ^ hU));
                }
                else
                {
                    HbyA = rAU*hUEqn.H();
                }

                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)
                  + fvc::interpolate(rAU)*fvc::ddtCorr(h, hU, phi)
                  - phih0
                );

                while (pimple.correctNonOrthogonal())
                {
                    fvScalarMatrix hEqn
                    (
                        fvm::ddt(h)
                      + fvc::div(phiHbyA)
                      - fvm::laplacian(ghrAUf, h)
                    );

                    hEqn.solve(mesh.solution().solver(h.select(pimple.finalInnerIter())));

                    if (pimple.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA + hEqn.flux();
                    }
                }

                hU = HbyA - rAU*h*magg*fvc::grad(h + h0);

                // Constrain the momentum to be in the geometry if 3D geometry
                if (mesh.nGeometricD() == 3)
                {
                    hU -= (gHat & hU)*gHat;
                }

                hU.correctBoundaryConditions();
            }
        }

        U.forceAssign(hU/h);
        hTotal.forceAssign(h + h0);

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
