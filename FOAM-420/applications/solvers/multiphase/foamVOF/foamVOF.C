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
    foamVOF

Group
    grpMultiphaseSolvers

Description
    Solver for multi-phase flows with two or more fluids

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "dynamicFvMesh/dynamicFvMesh.H"
#include "staticFvMesh/staticFvMesh.H"
#include "fvMatrices/solvers/MULES/CMULES.H"
#include "finiteVolume/ddtSchemes/EulerDdtScheme/EulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "finiteVolume/ddtSchemes/CrankNicolsonDdtScheme/CrankNicolsonDdtScheme.H"
#include "algorithms/subCycle/subCycle.H"
#include "immiscibleIncompressibleTwoPhaseMixture.H"
#include "barotropicCompressibilityModel/barotropicCompressibilityModel.H"
#include "incompressibleTwoPhaseMixture/incompressibleTwoPhaseMixture.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"
#include "cfdTools/general/solutionControl/pimpleControl/pimpleControl.H"
#include "cfdTools/general/fvOptions/fvOptions.H"
#include "cfdTools/general/CorrectPhi/CorrectPhi.H"
#include "finiteVolume/fvc/fvcSmooth/fvcSmooth.H"
#include "primitives/functions/Function1/Function1/Function1.H"
#include "matprop.H"
#include "vof.H"
#include "mulesAlgo.H"
#include "barocav.H"
#include "vofmodels.H"
#include "cavitate.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/checkAndCreateDynamicFvMesh.H"

    MatProp prop(runTime, mesh);

    autoPtr<FoamVOF> algo;
    label baroCavitate=runTime.controlDict().lookupOrDefault<label>("baroCavitate", 0);

    if (baroCavitate==0)
    {
        algo=new MulesAlgo(runTime, mesh,prop);
    }
    else
    {
        algo=new BaroCav(runTime, mesh,prop);
    }

    // #include "db/functionObjects/functionObjectList/postProcess.H"
    //#include "cfdTools/general/solutionControl/createControl.H"

    pimpleControl pimple(mesh);

    bool moveMeshOuterCorrectors
    (
      pimple.dict().lookupOrDefault<Switch>("moveMeshOuterCorrectors", false)
    );

    algo->moveMeshOuterCorrectors=moveMeshOuterCorrectors;
    algo->createTimeControls();

    algo->initContinuityErrs();

    algo->createFields(pimple);

    autoPtr<incompressible::turbulenceModel> turbulence;
    if (algo->baroCavitate()==0)
    {
        immiscibleIncompressibleTwoPhaseMixture &mixture=algo->mixture();

        turbulence=incompressible::turbulenceModel::New(algo->U(), algo->phi(), mixture);

    }
    else
    {
       incompressibleTwoPhaseMixture &mixture=algo->mixture2();
       turbulence=incompressible::turbulenceModel::New(algo->U(), algo->phi(), mixture);
    }

    // create physical models
    algo->createModels();

    fv::options& fvOptions(fv::options::New(mesh));
    turbulence->validate();

    if (algo->baroCavitate()==0)
    {
        algo->createAlphaFluxes();
        algo->correctPhi(pimple);

        if (!algo->LTS)
        {
            algo->readTimeControls();
            algo->CourantNo();
            algo->setInitialDeltaT();
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
      algo->updateFields
      (
          pimple,
          turbulence(),
          fvOptions
      );

    if (algo->solveT()>0)
    {
       Info<<"Solving heat transfer in VOF..."<<endl;
       algo->updateCp();
       algo->updateK();
       algo->TEqnSolve(turbulence());

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
