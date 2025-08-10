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


Description
    Sample a block around a solid object with uniform grid to generate
    the mask of the object, with level set values for each grid node.

\*---------------------------------------------------------------------------*/
#include "cfdTools/general/include/fvCFD.H"
#include "foamMap.H"
#include "blockSample.H"
#include "fieldSample.H"

int main(int argc, char *argv[])
{
	DynamicList<word> argNames;
	DynamicList<word> argValues;
	std::vector<word> keys;
	std::vector<word> values;

	argList::addOption
    (
        "blockDef",
        "blockDef",
        "Specify the block definition file name."
    );


	#include "include/setRootCase.H"

	if (argc>1)
    {
		for (label i=1;i<argc-1;i+=2)
		{
			argNames.append(argv[i]);
			argValues.append(argv[i+1]);
		}
	}

	word blockDef="blockDef";
	if (argNames.size()>=1)
	{
		for (label i=0;i<argNames.size();i++)
		{
			if (argNames[i]=="-blockDef")
			{
				blockDef=argValues[i];
				continue;
			}
		}
	}

    std::ifstream is(blockDef.c_str());
    if (!is)
    {
		Foam::Info<< "Error running level set utility, block definition file "
		<<blockDef<<" cannot open."<< Foam::endl;
		return -1;
	}



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

    const volScalarField& y = wallDist::New(mesh).y();
    scalar maxwdist=gMax(y)+1.0e-20;
    volScalarField y1( y/maxwdist );

    word s1;
	std::getline(is,s1,'\n');
	std::vector<word> vtmp;

	foamMap::stringtok(vtmp,s1," ");
	label nx=atoi(vtmp[0].c_str());
	label ny=atoi(vtmp[1].c_str());
	label nz=atoi(vtmp[2].c_str());
	scalar xmin=atof(vtmp[3].c_str());
	scalar ymin=atof(vtmp[4].c_str());
	scalar zmin=atof(vtmp[5].c_str());

	scalar xmax=atof(vtmp[6].c_str());
	scalar ymax=atof(vtmp[7].c_str());
	scalar zmax=atof(vtmp[8].c_str());
	while(!is.eof())
	{
		std::getline(is,s1,'\n');
		if (s1.length()<3) continue;
		foamMap::stringtok(vtmp,s1,"=");
		if (vtmp.size()!=2) continue;
		keys.push_back(vtmp[0]);
		values.push_back(vtmp[1]);
	}


    fieldSample grid
    (
    	nx,ny,nz,xmin,ymin,zmin,xmax,ymax,zmax
    );

	label sz=keys.size();
	for (label k=0;k<sz;k++)
	{
		word key=keys[k];
		word val=values[k];
		if (key=="bodySurfaceType")
		{
			grid.bodySurfaceType_=val;
			continue;
		}
		if (key=="byPhysicalType")
		{
			bool bypy=true;
			if (val=="false")
			{
				bypy=false;
			}
	        grid.byPhysicalType_=bypy;
			continue;
		}
		if (key=="bodySurfaceName")
		{
			word bname=val;
			grid.bodySurfaceName_=bname;
			continue;
		}
		if (key=="byName")
		{
			bool byname=true;
			if (val=="false" || val=="0")
			{
				byname=false;
			}
			grid.byName_=byname;
			continue;
		}
	}

    Info<< "Time now = " << runTime.timeName() << endl;
  	grid.construct_field_knn(&mesh);
  	grid.construct_body_knn(&mesh);
  	grid.sampleField (y1,&mesh,true);
	grid.write_samples("levelSet");

	 Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //
