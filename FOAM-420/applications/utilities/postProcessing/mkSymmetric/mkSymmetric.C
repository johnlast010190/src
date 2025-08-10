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
    (c) 2019-2020 Esi Ltd.

Description
    Make symmetric fields on patch surfaces based on non-symmetric solution fields
    using Smart mapping approach in Foam

\*---------------------------------------------------------------------------*/
#include "cfdTools/general/include/fvCFD.H"
#include "surfaceMap.H"

int main(int argc, char *argv[])
{

    #include "include/setRootCase.H"

    Foam::Info<< "Create time\n" << Foam::endl;

    Foam::Time runTime(Foam::Time::controlDictName, args);

    fvMesh mesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    surfaceMap surfMap
    (
        &mesh,
        &runTime
    );

    IOdictionary dict
    (
        IOobject
        (
            "foamMapDict",
            runTime.system(),
            "",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    surfMap.setInput(dict, args);

    surfMap.mapCase_=cwd();

    instantList sourceTimes = runTime.times();
    label sourceTimeIndex = runTime.timeIndex();

    if (surfMap.mapTimeName_ == "latestTime")
    {
        sourceTimeIndex = sourceTimes.size() - 1;
    }
    else
    {
        IStringStream is(surfMap.mapTimeName_);
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            readScalar(is)
        );
    }

    runTime.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);
    surfMap.mapTimeName_ = runTime.timeName();

    Info<< "\nSampling time " << runTime.timeName() << endl;

    surfMap.source_ =
        new readFields
        (
            &mesh,
            &runTime,
            &surfMap
        );

    surfMap.myid_       = Pstream::myProcNo();
    surfMap.masterid_   = Pstream::masterNo();
    surfMap.nprocs_     = Pstream::nProcs();
    surfMap.parInit(surfMap.nprocs_);

    surfMap.source().createFields(surfMap.mapTimeName_);
    surfMap.source().setInputs();
    surfMap.createMirrorFields(surfMap.mapTimeName_);

    surfMap.getParallelFields();
    surfMap.buildSearchTree();

    Info<<"Number of bnd maps:"<<surfMap.bndMaps_.size()<<endl;

    surfMap.getMirrorFields();

    surfMap.getSymFields();
    runTime.writeNow();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
