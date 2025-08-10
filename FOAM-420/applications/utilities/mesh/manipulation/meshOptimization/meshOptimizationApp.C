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
    (c) 2013 OpenFOAM Foundation

Application
    meshOptimization

Description
    A code for improving the mesh quality based on cell Sphericity.
\*---------------------------------------------------------------------------*/
//#include "global/argList/argList.H"
//#include "db/Time/Time.H"
//#include "fvMesh/fvMesh.H"
#include "hessianMeshOptimization/meshOptimization/hessianMeshOptimization.H"
#include "snappyHexMeshDriver/layerParameters/layerParameters.H"
#include "layerManipulate/layerManipulate.H"

using namespace Foam;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"

    Info<< "Starting Optimization Loop "<<endl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    IOdictionary dict
    (
        IOobject
        (
            "meshOptimizationDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );
    IOdictionary meshDict
    (
        IOobject
        (
            "foamHexMeshDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    hessianMeshOptimization meshOptim(mesh,dict);
    tmp<pointField> newPoints = meshOptim.newPoints();
    mesh.movePoints(newPoints);

    // layer addition parameters
    const dictionary& layerDict = meshDict.subDict("addLayersControls");
    layerParameters layerParams(layerDict, mesh.boundaryMesh());

    hexRef8 meshCutter
    (
        mesh,
        true,   //whether to try read level fields
        false   // do not try to read history.
    );

    layerManipulate layerManip
    (
        mesh,
        layerParams,
        meshCutter.cellLevel(),
        meshCutter.pointLevel(),
        meshCutter.level0EdgeLength()
    );
    layerManip.fitLayerPointStack();

    runTime++;
    mesh.write();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "\nEnd\n" << endl;
    return 0;

}


// ************************************************************************* //
