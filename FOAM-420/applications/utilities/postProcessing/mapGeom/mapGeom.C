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
#include "blockSample.H"
#include "triSurface/triSurface.H"
#include "array3d.H"
#include "foamMap.H"
#include "trainGeoModel.H"

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addNote("Map geometric object using surface mesh file in constant/triSurface");

    #include "include/setRootCase.H"

    Foam::Info<< "Create time\n" << Foam::endl;

	Foam::Time runTime(Foam::Time::controlDictName, args);


	IOdictionary dict
	(
		IOobject
		(
            "sampleGeomDict",
			runTime.system(),
			"",
			runTime,
			IOobject::MUST_READ,
			IOobject::NO_WRITE,
			false
		)
	);


	Info<< "Time now = " << runTime.timeName() << endl;

    word home="~/";
    fileName databaseDir(home / "foamAI"/"database"/"car3d");

    fileName dbsDir = dict.lookupOrDefault<fileName> ("databaseDir","default");
    if (dbsDir !="default")
    {
        databaseDir=dbsDir;
    }

    databaseDir= databaseDir.expand();
    Info<<"databaseDir:"<<databaseDir<<endl;

    fileName caseListFile(databaseDir / "caselist.dbs");
    fileName geomlistFile(databaseDir / "objectMask.dbs");
    bool exists=is_file(geomlistFile.c_str());

	if (!exists)
	{
		FatalErrorInFunction<<"Error mapping flow field, database file: "
		<<geomlistFile<<" not found or cannot open."
		<< exit(FatalError);
	}

    trainGeoModel train
    (
		&runTime,
		dict
	);


    Info<<"Read database, this will take a while...."<<endl;
	train.readDatabase(geomlistFile);

    vector ndiv0(127, 63, 63);

    vector ndiv=dict.lookupOrDefault<vector> ("ndivisions",ndiv0);

    label nx=label(ndiv[0]);
    label ny=label(ndiv[1]);
    label nz=label(ndiv[2]);

    Info<<"Grid divisions: "<<nx<<" "<<ny<<" "<<nz<<endl;

    train.build_model( nx, ny,  nz);

    train.get_prediction();
	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
