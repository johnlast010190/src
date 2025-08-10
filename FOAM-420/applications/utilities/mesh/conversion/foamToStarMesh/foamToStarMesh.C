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
    (c) 2016 OpenCFD Ltd.

Application
    foamToStarMesh

Group
    grpMeshConversionUtilities

Description
    Reads an OpenFOAM mesh and writes a pro-STAR (v4) bnd/cel/vrt format.

Usage
    \b foamToStarMesh [OPTION]

    Options:
      - \par -noBnd
        Suppress writing the \c .bnd file

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 1000 (scale \em [m] to \em [mm]).

Note
    The cellTable information available in the files
    \c constant/cellTable and \c constant/polyMesh/cellTableId
    will be used if available. Otherwise the cellZones are used when
    creating the cellTable information.

See also
    Foam::cellTable, Foam::meshWriter and Foam::fileFormats::STARCDMeshWriter

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/timeSelector.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "starcd/STARCDMeshWriter.H"
#include "db/IOobjects/IOdictionary/IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "read OpenFOAM mesh and write a pro-STAR (v4) bnd/cel/vrt format"
    );
    argList::noParallel();
    timeSelector::addOptions();

    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1000 ([m] to [mm])"
    );
    argList::addBoolOption
    (
        "noBnd",
        "suppress writing a boundary (.bnd) file"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    fileName exportName = meshWriter::defaultMeshName;
    if (args.optionFound("case"))
    {
        exportName += '-' + args.globalCaseName();
    }

    // Default rescale from [m] to [mm]
    const scalar scaleFactor = args.optionLookupOrDefault("scale", 1000.0);
    const bool  writeBndFile = !args.optionFound("noBnd");

    #include "include/createPolyMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != polyMesh::UNCHANGED)
        {
            fileFormats::STARCDMeshWriter writer
            (
                mesh,
                scaleFactor,
                writeBndFile
            );

            mkDir(runTime.rootPath()/runTime.caseName()/"meshExport"/"star");
            fileName meshName("meshExport/star"/exportName);
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
