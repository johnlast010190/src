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
    (c) 2017 OpenFOAM Foundation

Application
    Test-decomposedBlockData

Description
    Convert decomposedBlockData into its components.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "db/IOobjects/decomposedBlockData/decomposedBlockData.H"
#include "db/IOstreams/Fstreams/OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::validArgs.append("file");
    #include "setRootCase.H"

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "Run in parallel" << exit(FatalError);
    }

    #include "createTime.H"

    const fileName file(args[1]);

    Info<< "Reading " << file << nl << endl;
    decomposedBlockData data
    (
        Pstream::worldComm,
        IOobject
        (
            file,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    data.rename(data.name() + "Data");
    fileName objPath(data.objectPath());
    mkDir(objPath.path());
    Info<< "Opening output file " << objPath << nl << endl;
    OFstream os
    (
        objPath,
        IOstream::BINARY,
        IOstream::currentVersion,
        runTime.writeCompression()
    );
    if (!os.good())
    {
        FatalErrorInFunction
            << "Failed opening " << objPath << exit(FatalError);
    }

    if (!data.writeData(os))
    {
        FatalErrorInFunction
            << "Failed writing " << objPath << exit(FatalError);
    }

    return 0;
}


// ************************************************************************* //
