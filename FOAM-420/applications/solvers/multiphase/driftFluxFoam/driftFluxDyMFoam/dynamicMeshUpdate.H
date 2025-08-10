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
    (c) 2018-2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

if (pimple.firstIter() || moveMeshOuterCorrectors)
{
    scalar timeBeforeMeshUpdate = runTime.elapsedCpuTime();

    mesh.update();

    if (mesh.changing())
    {
        Info<< "Execution time for mesh.update() = "
            << runTime.elapsedCpuTime() - timeBeforeMeshUpdate
            << " s" << endl;

        // Do not apply previous time-step mesh compression flux
        // if the mesh topology changed
        if (mesh.topoChanging())
        {
            talphaPhiCorr0.clear();
        }

        gh = (g & mesh.C()) - ghRef;
        ghf = (g & mesh.Cf()) - ghRef;
    }

    if ((mesh.changing() && correctPhi) || mesh.topoChanging())
    {
        // Calculate absolute flux from the mapped surface velocity
        // Note: temporary fix until mapped Uf is assessed
        Uf = fvc::interpolate(U);

        // Calculate absolute flux from the mapped surface velocity
        phi = mesh.Sf() & Uf;

        #include "interDyMFoam/correctPhi.H"

        // Make the flux relative to the mesh motion
        fvc::makeRelative(phi, U);

        mixture.correct();
    }

    if (mesh.changing() && checkMeshCourantNo)
    {
        #include "include/meshCourantNo.H"
    }
}
