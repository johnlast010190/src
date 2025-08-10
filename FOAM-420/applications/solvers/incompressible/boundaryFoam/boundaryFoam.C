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
    (c) 2010-2016 Esi Ltd.

Application
    boundaryFoam

Group
    grpIncompressibleSolvers

Description
    Steady-state solver for incompressible, 1D turbulent flow, typically to
    generate boundary layer conditions at an inlet.

    Boundary layer code to calculate the U, k and epsilon distributions.
    Used to create inlet boundary conditions for experimental comparisons
    for which U and k have not been measured.
    Turbulence model is runtime selectable.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "graphField/makeGraph.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();

    #include "include/addRegionOption.H"

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    word regionName(fvMesh::defaultRegion);

    if (args.optionReadIfPresent("region", regionName))
    {
        Info<< "Create mesh " << regionName << " for time = "
            << runTime.timeName() << nl << endl;
    }

    #include "createMesh.H"
    #include "createFields.H"
    #include "interrogateWallPatches.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        fvVectorMatrix divR(turbulence->divDevReff(U));
        divR.source() = flowMask & divR.source();

        fvVectorMatrix UEqn
        (
            divR == gradP + fvOptions(U)
        );

        UEqn.relax();

        fvOptions.constrain(UEqn);

        UEqn.solve();

        fvOptions.correct(U);


        // Correct driving force for a constant volume flow rate
        dimensionedVector UbarStar = flowMask & U.weightedAverage(mesh.V());

        U += (Ubar - UbarStar);
        gradP += (Ubar - UbarStar)/(1.0/UEqn.A())().weightedAverage(mesh.V());

        laminarTransport.correct();
        turbulence->correct();

        fvOptions.correct();

        Info<< "Uncorrected Ubar = " << (flowDirection & UbarStar.value())
            << ", pressure gradient = " << (flowDirection & gradP.value())
            << endl;

        #include "evaluateNearWall.H"

        if (runTime.writeTime())
        {
            #include "makeGraphs.H"
        }

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
