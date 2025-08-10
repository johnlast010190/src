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
    wallFunctionTable

Group
    grpPreProcessingUtilities

Description
    Generates a table suitable for use by tabulated wall functions.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "tabulatedWallFunction/tabulatedWallFunction/tabulatedWallFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"

    IOdictionary dict
    (
        IOobject
        (
            "wallFunctionDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED
        )
    );

    autoPtr<tabulatedWallFunctions::tabulatedWallFunction>
        twf(tabulatedWallFunctions::tabulatedWallFunction::New(dict, mesh));

//    twf->writeData(Info);

    twf->write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
