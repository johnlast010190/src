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
    (c) 2009-2019 Esi Ltd.

Description
    Sample a number of sub-blocks with uniform grid from a flow field
    of unstructured mesh.

\*---------------------------------------------------------------------------*/
#include "cfdTools/general/include/fvCFD.H"
#include "foamMap.H"
#include "blockSample.H"
#include "sampleParrun.H"

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "sourceCase",
        "dir",
        "Specify source case directory"
    );

    argList::addOption
    (
        "sourceTime",
        "scalar|latestTime",
        "Specify the time step/iteration of the source case to map from (default: latestTime)"
    );

    argList::addOption
    (
        "mapScalarFields",
        "p,nut",
        "Specify the scalar field names(separated by \",\") to be mapped from source to target"
    );

    argList::addOption
    (
        "mapVectorFields",
        "U",
        "Specify the scalar field names (separated by \",\") to be mapped from source to target"
    );

    argList::addOption
    (
        "nwdist",
        "200",
        "Specify the number of divisions in the wall distance range[0,1]"
    );

    argList::addOption
    (
        "fieldTypes",
        "p:pressure,U:velocity,nut:turbulence viscosity,k:turbulence energy",
        "Specify the property types of each maping field"
    );

    argList::addOption
    (
        "interpolation",
        "false",
        "Specify whether interpolation will be used in the mapping"
    );

    argList:: addNote("foamSample  samples a number of sub-blocks with uniform grid from a flow field\n\
with arbitrary unstructured mesh.");

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

    foamSample *gridSample;

    if (!Pstream::parRun())
    {
        gridSample =
            new foamSample
            (
                &runTime,
                &mesh
            );
    }
    else
    {
        gridSample =
            new samplePar
            (
                &runTime,
                &mesh
            );
    }

    gridSample->setInput(dict, args);
    gridSample->nproc_=Pstream::nProcs();
    gridSample->casedir_=cwd();

    if (Pstream::parRun())
    {
        gridSample->myid_     = Pstream::myProcNo();
        gridSample->masterid_ = Pstream::masterNo();
        gridSample->parInit(Pstream::nProcs());
    }

    gridSample->setOptions(args);

    instantList sourceTimes = runTime.times();
    label sourceTimeIndex = runTime.timeIndex();

    if (gridSample->mapTimeName_ == "latestTime")
    {
        sourceTimeIndex = sourceTimes.size() - 1;
    }
    else
    {
        IStringStream is(gridSample->mapTimeName_);
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            readScalar(is)
        );
    }

    runTime.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);
    gridSample->mapTimeName_ = runTime.timeName();

    Info<< "\nSampling time " << runTime.timeName() << endl;

    gridSample->source_ =
        new readFields
        (
            &mesh,
            &runTime,
            gridSample
        );

    gridSample->source().createFields(gridSample->mapTimeName_);
    gridSample->source().setInputs();

    if (!Pstream::parRun())
    {
      gridSample->buildSearchTrees();
      gridSample->source().storeFields();
      gridSample->getSampleData();//no normalization of sample data
    }
    else
    {
        const boundBox &bbox=gridSample->sampleBox();
        gridSample->source().storeFields(bbox);

        gridSample->wait();

        gridSample->combineFields();
        gridSample->getBodySurface();

        gridSample->wait();
        Info<<"construct knn...."<<endl;
        if (gridSample->master())
        {
            gridSample->constructKnn();
            gridSample->getSampleData();
        }

    }


    //label itime=atoi(gridSample->source().mapTime_.c_str())-atoi(runTime.timeName().c_str());
    //runTime+=itime;

    if (Pstream::parRun())
    {
        gridSample->wait();
    }

    if (gridSample->master() || !Pstream::parRun())
    {
        gridSample->writeNodeField();
        gridSample->writeVisualization();
        gridSample->saveToDatabase();
    }

    runTime.writeNow();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;
    delete gridSample;

    return 0;
}

// ************************************************************************* //
