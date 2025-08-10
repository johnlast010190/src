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
    (c) 2016-2017 OpenCFD Ltd.

Application
    lumpedPointZones

Description
    Produces a VTK PolyData file \c lumpedPointZones.vtp in which the
    segmentation of the pressure integration zones can be visualized
    for diagnostic purposes. Does not use external coupling.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "db/Time/timeSelector.H"

#include "lumpedPointTools.H"
#include "lumpedPointIOMovement.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create lumpedPointZones.vtp to verify the segmentation of "
        "pressure integration zones used by lumpedPoint BC."
    );
    argList::noParallel();    // The VTP writer is not yet in parallel

    argList::noFunctionObjects();  // Never use function objects
    argList::addBoolOption
    (
        "verbose",
        "increased verbosity"
    );

    #include "include/addRegionOption.H"
    #include "include/setRootCase.H"

    // const bool verbose = args.optionFound("verbose");

    #include "include/createTime.H"

    runTime.setTime(instant(0, runTime.constant()), 0);

    #include "include/createNamedPolyMesh.H"

    autoPtr<lumpedPointIOMovement> movement = lumpedPointIOMovement::New
    (
        runTime
    );

    if (!movement.valid())
    {
        Info<< "no valid movement found" << endl;
        return 1;
    }

    const labelList patchLst = lumpedPointTools::lumpedPointPatchList(mesh);
    if (patchLst.empty())
    {
        Info<< "no patch list found" << endl;
        return 2;
    }

    pointIOField points0 = lumpedPointTools::points0Field(mesh);
    movement().setMapping(mesh, patchLst, points0);

    // Initial geometry, but with zone colouring
    movement().writeZonesVTP("lumpedPointZones.vtp", mesh, points0);

    // Initial positions/rotations
    movement().writeStateVTP("initialState.vtp");

    Info<< nl
        << "wrote 'lumpedPointZones.vtp'" << nl
        << "wrote 'initialState.vtp'" << nl
        << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
