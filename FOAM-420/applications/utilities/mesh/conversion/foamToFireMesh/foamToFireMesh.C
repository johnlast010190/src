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
    (c) 2016 OpenCFD Ltd.

Application
    foamToFireMesh

Description
    Reads an OpenFOAM mesh and writes an AVL/FIRE fpma format

Usage
    \b foamToFireMesh [OPTION]

    Options:
      - \par -ascii
        Write in ASCII format instead of binary

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 1 (ie, no scaling).

See also
    Foam::cellTable, Foam::meshWriter and Foam::fileFormats::FIREMeshWriter

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "fire/FIREMeshWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "read OpenFOAM mesh and write an AVL/FIRE fpma format"
    );
    argList::noParallel();
    timeSelector::addOptions();

    argList::addBoolOption
    (
        "ascii",
        "write in ASCII format instead of binary"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1 (none)"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    fileName exportName = meshWriter::defaultMeshName;
    if (args.optionFound("case"))
    {
        exportName += '-' + args.globalCaseName();
    }


    // write control options
    // ~~~~~~~~~~~~~~~~~~~~~
    fileFormats::FIREMeshWriter::binary = !args.optionFound("ascii");

    // default: rescale from [m] to [mm]
    scalar scaleFactor = 1;
    if (args.optionReadIfPresent("scale", scaleFactor))
    {
        if (scaleFactor <= 0)
        {
            scaleFactor = 1;
        }
    }

    #include "include/createPolyMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != polyMesh::UNCHANGED)
        {
            fileFormats::FIREMeshWriter writer(mesh, scaleFactor);

            fileName meshName(exportName);
            if (state != polyMesh::UNCHANGED)
            {
                meshName += '_' + runTime.timeName();
            }

            writer.write(meshName);
        }

        Info<< nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
