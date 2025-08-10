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
    Find the best map of flow configuration from the database to the
    current case based on machine-learning method.

\*---------------------------------------------------------------------------*/
#include "cfdTools/general/include/fvCFD.H"
#include <unistd.h>
#include "foamMap.H"
#include "blockSample.H"
#include "trainGeoModel.H"

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


	trainGeoModel train
    (
		&runTime,
		&mesh
	);

	train.casedir_=cwd();

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

	train.setInput
	(
		dict,
        args
	);

	Info<< "Time now = " << runTime.timeName() << endl;
	word dbsfile=foamSample::databaseDir_+"/"+train.caseType_+"/objectMasks.dbs";

	bool exists=is_file(dbsfile.c_str());

	if (!exists)
	{
		FatalErrorInFunction<<"Error mapping flow field, database file: "
		<<dbsfile<<" not found or cannot open."
		<< exit(FatalError);
	}

	train.source_=new
	readFields
    (
		&mesh,
		&runTime,
		&train
	);

	train.source().setInputs();
	Info<<"Read database, this will take a while...."<<endl;
	train.readDatabase(dbsfile);
	std::map<word,gridField>::iterator it;
	it=train.sampleGrids_.begin();
	label m=it->second.nx_+1;
	label n=it->second.ny_+1;
	label k=it->second.nzpt();
	Info<<"nx,ny,nz:"<<m<<" "<<n<<" "<<k<<endl;
	train.build_model( m, n,  k);

	train.buildSearchTrees();
	train.getGeomSample();


	train.get_prediction();
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
