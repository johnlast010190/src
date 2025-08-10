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
    (c) 2017 OpenCFD Ltd.
    (c) 2012-2014 OpenFOAM Foundation

Application
    foamHelp

Group
    grpMiscUtilities

Description
    Top level wrapper utility around foam help utilities

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "helpTypes/helpType/helpType.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/addRegionOption.H"
    #include "addToolOption.H"

    // Intercept request for help
    if ((argc > 0) && (strcmp(argv[1], "-help") == 0))
    {
        #include "include/setRootCase.H"
    }

    if (argc < 2)
    {
        FatalError
            << "No help utility has been supplied" << nl
            << exit(FatalError);
    }

    word utilityName = argv[1];
    Foam::autoPtr<Foam::helpType> utility
    (
        helpType::New(utilityName)
    );

    utility().init();

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createNamedMesh.H"

    utility().execute(args, mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
