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
    (c) 2021 Esi Ltd.

Application
    sampleGeom

Description
    Sample a geometry file defined by STL or obj format, and save it to database.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/include/fvCFD.H"
#include "foamMap.H"
#include "blockSample.H"
#include "triSurface/triSurface.H"
#include "array3d.H"

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addNote("Sample geometric object using surface mesh file in constant/triSurface");

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

    word home="~/";
    fileName databaseDir(home / "foamAI"/"database"/"car3d");

    fileName dbsDir(dict.lookupOrDefault<fileName>("databaseDir","default"));
    if (dbsDir !="default")
    {
        databaseDir=dbsDir;
    }

    databaseDir= databaseDir.expand();
    Info<<"databaseDir:"<<databaseDir<<endl;

    mkDir(databaseDir);


    fileName caseListFile(databaseDir / "caselist.dbs");
    fileName geomlistFile(databaseDir / "objectMask.dbs");

    word gfile=dict.lookupOrDefault<word> ("geofile","car.stl");

    word cdir=cwd();
    fileName bfile(cdir /"constant"/"triSurface"/ gfile);

    triSurface surfaces(bfile);

    const vectorField& cf=surfaces.Cf();
    Info<<"Number of surface elements : "<<cf.size()<<endl;

    point pmin=min(cf);
    point pmax=max(cf);

    point dp=pmax-pmin;
    Info<<"Bounding box for surface mesh: "<<pmin<<" "<<pmax<<endl;
    Info<<"Span of geometry: "<<dp<<endl;

    point p0(-0.5,-0.5, -0.1);
    point p1(4, 2.2, 2);

    boundBox b0(p0,p1);

    boundBox bbox=dict.lookupOrDefault<boundBox> ("boundBox",b0);

    Info<<"Boundbox for surface grid: "<<bbox<<endl;

    vector ndiv0(127, 63, 63);

    vector ndiv=dict.lookupOrDefault<vector> ("ndivisions",ndiv0);

    label nx=label(ndiv[0]);
    label ny=label(ndiv[1]);
    label nz=label(ndiv[2]);

    Info<<"Grid divisions: "<<nx<<" "<<ny<<" "<<nz<<endl;

    Array3d<label> mask;
    mask.setSize(nx,ny,nz);
    mask=0;

    pmin=bbox.min();
    pmax=bbox.max();

    dp=pmax-pmin;
    scalar dx=dp[0]/float(nx-1);
    scalar dy=dp[1]/float(ny-1);
    label nzpt=max(1,nz-1);
    scalar dz=dp[2]/float(nzpt);

    // update mask

    forAll(cf, ip)
    {
        point pt=cf[ip];
        if (!bbox.contains(pt))
        {
            continue;
        }

        label i=label((pt[0]-pmin[0])/dx);
        label j=label((pt[1]-pmin[1])/dy);
        label k=label((pt[2]-pmin[2])/dz);
        i=min(i,nx-1);
        j=min(j,ny-1);
        k=min(k,nz-1);
        mask(i,j,k)=1;
    }

    label cnt=0;
    for (label i=0; i<nx; i++)
    {
        for (label j=0; j<ny; j++)
        {
            for (label k=0; k<nz; k++)
            {
                cnt+=mask(i,j,k);
            }
        }
    }

    Info<<"Number of surface grids:"<<cnt<<endl;
    Info<<"Total grids:"<<nx*ny*nz<<endl;

    // save to database

    std::vector<word> cases;
    std::ifstream is(caseListFile.c_str());

    if (is)
    {
        word str;
        while (!is.eof())
        {
            std::getline(is,str,'\n');
            if (str.length()<2)
            {
                continue;
            }
            cases.push_back(str);
        }
        is.close();
    }

    bool exists=false;
    if (cases.size()>=1)
    {
        std::vector<word>::iterator it;
        it=std::find (cases.begin(), cases.end(), cdir);
        if (it!=cases.end())
        {
            exists=true;
        }
    }
    if (exists)
    {
        Info<<"Warning: Current case already in database."<<endl;
    }
    else
    {
        std::ofstream os1(caseListFile.c_str(),std::ios::app);
        os1<<cdir<<std::endl;
        os1.close();
        word type="car3d";
        word name="surfaceGrid";

        std::ofstream os(geomlistFile.c_str(),std::ios::app);
        word desc=type+","
             +cdir+","+name
             +","+std::to_string(nx)
             +","+std::to_string(ny)
             +","+std::to_string(nz)
             +","+std::to_string(pmin[0])
             +","+std::to_string(pmin[1])
             +","+std::to_string(pmin[2])
             +","+std::to_string(pmax[0])
             +","+std::to_string(pmax[1])
             +","+std::to_string(pmax[2]);
        os<<desc<<std::endl;
        word mcnts="mask";
        for (label i=0;i<nx;++i)
        {
            for (label j=0;j<ny;++j)
            {
                for (label k=0;k<nz;++k)
                {
                    label insolid=mask.getval(i,j,k);
                    if (insolid==1)
                    {
                        mcnts+=(","
                            +std::to_string(i)
                            +":"+std::to_string(j)
                            +":"+std::to_string(k));
                    }
                }
            }
        }
        os<<mcnts<<std::endl;
    }

    Info<<"Geometry mask saved to file: "<< geomlistFile <<endl;
    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	  << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	  << nl << endl;

    Info<< "End\n" << endl;

	return 0;
}

// ************************************************************************* //

