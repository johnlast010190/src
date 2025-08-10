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
    (c) 2015 OpenCFD Ltd.
    (c) 2011-2015 OpenFOAM Foundation

Application
    surfaceFeatureConvert

Group
    grpSurfaceUtilities

Description
    Convert between edgeMesh formats.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"

#include "edgeMesh/edgeMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert between edgeMesh formats"
    );
    argList::noParallel();
    argList::validArgs.append("inputFile");
    argList::validArgs.append("outputFile");
    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1"
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    const fileName importName = args[1];
    const fileName exportName = args[2];

    // Disable inplace editing
    if (importName == exportName)
    {
        FatalErrorInFunction
            << "Output file " << exportName << " would overwrite input file."
            << exit(FatalError);
    }

    // Check that reading/writing is supported
    if
    (
        !edgeMesh::canReadType(importName.ext(), true)
     || !edgeMesh::canWriteType(exportName.ext(), true)
    )
    {
        return 1;
    }

    edgeMesh mesh(importName);

    Info<< "\nRead edgeMesh " << importName << nl;
    mesh.writeStats(Info);
    Info<< nl
        << "\nwriting " << exportName;

    scalar scaleFactor = 0;
    if (args.optionReadIfPresent("scale", scaleFactor) && scaleFactor > 0)
    {
        Info<< " with scaling " << scaleFactor << endl;
        mesh.scalePoints(scaleFactor);
    }
    else
    {
        Info<< " without scaling" << endl;
    }

    mesh.write(exportName);
    mesh.writeStats(Info);

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
