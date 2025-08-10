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
    (c) 2019 Esi Ltd.

Application
    foamBendingWave

Description
    Solver for bending waves traveling at a surface excited by an external
    force.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "include/faCFD.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createMesh.H"
    #include "createFaMesh.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    //List<vectorField> obsProc(Pstream::nProcs());
    //obsProc[Pstream::myProcNo()] = aMesh.areaCentres();
    //Pstream::gatherList(obsProc);
    //Pstream::scatterList(obsProc);
    //vectorField obs = ListListOps::combine<vectorField>(obsProc,accessOp<vectorField>());
    //Info<< obs;

    if (aMesh.areaCentres().size() > 0)
    {
        fileName dirPath;

        if (Pstream::parRun())
        {
            dirPath = fileName("..")/".."/"CAA"/"processor"+std::to_string(Pstream::myProcNo());
        }
        else
        {
            dirPath = fileName("..")/"CAA";
        }

        vectorIOField obs
        (
            IOobject
            (
                "observers",
                runTime.system(),
                dirPath,
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            aMesh.areaCentres()
        );

        obs.write();
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
