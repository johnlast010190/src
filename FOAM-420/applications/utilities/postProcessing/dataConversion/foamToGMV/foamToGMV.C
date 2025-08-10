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
    (c) 2011-2017 OpenFOAM Foundation

Application
    foamToGMV

Group
    grpPostProcessingUtilities

Description
    Translates foam output to GMV readable files.

    A free post-processor with available binaries from
    http://www-xdiv.lanl.gov/XCM/gmv/

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/Time/instant/instantList.H"
#include "db/IOobjectList/IOobjectList.H"
#include "itoa.H"
#include "Cloud/Cloud.H"
#include "passiveParticle/passiveParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    const label nTypes = 4;
    const word fieldTypes[] =
    {
        "volScalarField",
        "volVectorField",
        "surfaceScalarField",
        cloud::prefix
    };

    #include "include/setRootCase.H"

    #include "include/createTime.H"
    #include "include/createMesh.H"

    #include "readConversionProperties.H"

    // get the available time-steps
    instantList TimeList = runTime.times();
    Info<< TimeList << endl;
    label nTimes = TimeList.size();

    for (label n=1; n < nTimes; n++)
    {
        if (TimeList[n].value() > startTime)
        {
            Info<< "Time = " << TimeList[n].value() << nl;

            // Set Time
            runTime.setTime(TimeList[n], n);
            word CurTime = runTime.timeName();

            IOobjectList objects(mesh, runTime.timeName());

            #include "moveMesh.H"

            // set the filename of the GMV file
            fileName gmvFileName = "plotGMV." + itoa(n);
            OFstream gmvFile(args.rootPath()/args.caseName()/gmvFileName);

            #include "gmvOutputHeader.H"
            #include "gmvOutput.H"
            #include "gmvOutputTail.H"
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
