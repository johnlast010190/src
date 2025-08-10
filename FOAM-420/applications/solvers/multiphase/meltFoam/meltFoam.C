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
    (c) 2011 OpenFOAM Foundation
    (c) 2010-2016 Esi Ltd.

Reference
    Solver from:
        @article{rosler_shell-and-tube_2011,
            title = {Shell-and-tube type latent heat thermal energy storage:
                    numerical analysis and comparison with experiments},
            volume = {47},
            issn = {0947-7411, 1432-1181},
            shorttitle = {Shell-and-tube type latent heat thermal energy storage},
            url = {http://link.springer.com/article/10.1007/s00231-011-0866-9},
            doi = {10.1007/s00231-011-0866-9},
            language = {en},
            number = {8},
            urldate = {2013-01-07},
            journal = {Heat and Mass Transfer},
            author = {Rösler, Fabian and Brüggemann, Dieter},
            month = aug,
            year = {2011},
            pages = {1027--1033}
        }

Application
    meltFoam

Description
    Solves a convection dominated solid/liquid phase change process
    using the enthalpy-porousity method with bouyancy in laminar
    systems.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "global/constants/mathematical/mathematicalConstants.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "cfdTools/general/include/readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "cfdTools/general/include/initContinuityErrs.H"
    #include "cfdTools/general/include/createTimeControls.H"
    #include "cfdTools/incompressible/CourantNo.H"
    #include "cfdTools/general/include/setInitialDeltaT.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "cfdTools/general/include/createTimeControls.H"
        #include "cfdTools/incompressible/CourantNo.H"
        #include "cfdTools/general/include/setDeltaT.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

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
