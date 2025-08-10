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
    (c) 2019-2021 Esi Ltd.

Description
    Sample a block around a solid object with uniform grid to generate
    the mask of the object.

\*---------------------------------------------------------------------------*/
#include "cfdTools/general/include/fvCFD.H"
#include "foamMap.H"
#include "blockSample.H"
#include "sampleParrun.H"

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

    //gridSample->source().createFields(gridSample->mapTimeName_);
    gridSample->source().setInputs();

    if (!Pstream::parRun())
    {
        gridSample->buildSearchTrees();
        gridSample->getSampleData("geometry");
    }
    else
    {
        // const boundBox &bbox=gridSample->sampleBox();
        // gridSample->source().storeFields(bbox);

        gridSample->wait();

       // gridSample->combineFields();
        gridSample->getBodySurface();

        gridSample->wait();
        Info<<"construct knn...."<<endl;
        if (gridSample->master())
        {
            gridSample->constructKnn();
            gridSample->getSampleData("geometry");
        }

    }

    if (Pstream::parRun())
    {
        gridSample->wait();
    }


    if (!Pstream::parRun() || gridSample->master())
    {

        gridSample->saveObjectMasks(gridSample->caseType_);
    }

    Info<<"Geometry mask saved to file: "<<gridSample->databasePath()<<endl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    Info<< "End\n" << endl;

    delete gridSample;

	return 0;
}

// ************************************************************************* //
