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
    (c) 2011-2013 OpenFOAM Foundation

Application
    engineCompRatio

Group
    grpPostProcessingUtilities

Description
    Calculate the geometric compression ratio.

    Note: if you have valves and/or extra volumes it will not work,
          since it calculates the volume at BDC and TCD.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "engineTime/engineTime.H"
#include "engineMesh/engineMesh/engineMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createEngineTime.H"
    #include "include/createEngineMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    scalar eps = 1.0e-10;
    scalar fullCycle = 360.0;

    scalar ca0 = -180.0;
    scalar ca1 = 0.0;

    while (runTime.theta() > ca0)
    {
        ca0 += fullCycle;
        ca1 += fullCycle;
    }

    while (mag(runTime.theta() - ca0) > eps)
    {
        scalar t0 = runTime.userTimeToTime(ca0 - runTime.theta());
        runTime.setDeltaT(t0);
        runTime++;
        Info<< "CA = " << runTime.theta() << endl;
        mesh.moveMesh();
    }

    scalar Vmax = sum(mesh.V().field());

    while (mag(runTime.theta() - ca1) > eps)
    {
        scalar t1 = runTime.userTimeToTime(ca1 - runTime.theta());
        runTime.setDeltaT(t1);
        runTime++;
        Info<< "CA = " << runTime.theta() << endl;
        mesh.moveMesh();
    }

    scalar Vmin = sum(mesh.V().field());

    Info<< "\nVmax = " << Vmax
        << ", Vmin = " << Vmin << nl
        << "Vmax/Vmin = " << Vmax/Vmin << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
