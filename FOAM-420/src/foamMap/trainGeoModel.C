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
    Train a machine-learning model to recognize the closest flow topology
    in the database to the one in consideration.


\*---------------------------------------------------------------------------*/

#include "trainGeoModel.H"
#include "npgrid.H"


Foam::trainGeoModel::trainGeoModel()
:
    foamSample()
{
	task_="train";
	K_=2;
}

Foam::trainGeoModel::trainGeoModel
(
    const Time* runTime,
    const dictionary& dict
)
:
    foamSample(runTime, nullptr)
{
	task_="train";
	K_=2;

    word gfile=dict.lookupOrDefault<word> ("geofile","car.stl");

    word cdir=cwd();
    casedir_ = cdir;
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

    gridField grid(bbox, nx-1,ny-1,nzpt);

    grid.setFieldSize();
    grid.insolids_ = mask;

    sampleGrids_["geomObject"] = grid;

    Info<<"Number of surface grids:"<<cnt<<endl;
    Info<<"Total grids:"<<nx*ny*nz<<endl;

}



Foam::trainGeoModel::trainGeoModel
(
	const Time* runTime,
	const fvMesh* mesh
)
:
    foamSample(runTime, mesh)
{
	task_="train";
	K_=2;
}

void Foam::trainGeoModel::readDatabase1
(
    const word& dbsfile
)
{
	word input;
	npGrid::read_text_file(input,dbsfile);
	std::vector<word> vtmp;
	stringtok(vtmp,input,"\n");
	vcaseItems_.resize(0);
	for (unsigned i=0; i<vtmp.size()-1;i+=2)
	{
		const word &s1=vtmp[i];
		const word &s2=vtmp[i+1];
		if (s1.length()<3 ||s2.size()<3) break;
		dataItem ditem;
		ditem.setDataInfo(s1);
		ditem.setData(s2);
		vcaseItems_.push_back(ditem);
	}
	Info<<"number of database items read:"<<vcaseItems_.size()<<endl;
}

void Foam::trainGeoModel::readDatabase
(
    const word& dbsfile
)
{
	std::ifstream is(dbsfile.c_str());
	if (!is)
	{
		FatalErrorInFunction<<"Error read database, database file: "
		<<dbsfile<<" not found or cannot open."
		<< exit(FatalError);
	}

	word line;
	word cnts;

	vcaseItems_.resize(0);
	while(!is.eof())
	{
		std::getline(is,line,'\n');
		if (line.length()<3)
		{
			break;
		}

		dataItem ditem;
		ditem.setDataInfo(line);

		std::getline(is,cnts,'\n');
		if (cnts.length()<3)
		{
			break;
		}

		ditem.setData(cnts);
		vcaseItems_.push_back(ditem);
	}
	Info<<"number of database items read:"<<vcaseItems_.size()<<endl;
}

void Foam::trainGeoModel::setOptions
(
	const DynamicList<word>& argNames,
	const DynamicList<word>& argValues
)
{
	if (argNames.size()==0)
    {
        return;
    }
	for (label i=0; i<argNames.size();i++)
	{
		word key=argNames[i];
		word val=argValues[i];

		if (key=="-databaseDir")
		{
			foamSample::databaseDir_=val;
			continue;
		}

		if (key=="-task")
		{
			task_=val;
			continue;
		}

	}
}

void Foam::trainGeoModel::getCaseLabels()
{
	caseNames_.resize(0);

	std::map<word,dataItem>::iterator it;
	for (it=caseItems_.begin();it!=caseItems_.end();it++)
	{
		caseNames_.push_back(it->first);
	}

	std::sort(caseNames_.begin(),caseNames_.end());
	for (unsigned i=0;i<caseNames_.size();i++)
	{
		 caseLabels_[caseNames_[i]]=i;
	}

	for (it=caseItems_.begin();it!=caseItems_.end();it++)
	{
		word cname=it->first;
		it->second.identity(caseLabels_[cname]);
	}
}


void Foam::trainGeoModel::readCaseList()
{
	word casefile=foamSample::databaseDir_+"/caselist.dbs";
	std::ifstream is(casefile.c_str());
	if (!is)
	{
		FatalErrorInFunction
			<< "Error read training case list, \n"<<casefile<<" cannot open for reading."
			<< exit(FatalError);
	}
	word str;
	while (!is.eof())
	{
		std::getline(is,str,'\n');
		if (str.length()<6)
		{
			continue;
		}
		dataItem ditem(str);
		caseItems_[ditem.caseName()]=ditem;
	}
	is.close();
	Info<<"number of cases read:"<<caseItems_.size()<<endl;
}

void Foam::trainGeoModel::readVector
(
	std::vector<scalar>& vect,
    std::ifstream& is
)
{
	word str;
	for (unsigned i=0;i<vect.size();i++)
	{
		is>>vect[i];
	}
	is>>str; //eat the end of line sign
}

void Foam::trainGeoModel::clearGeoData()
{
	for (unsigned i=0; i<vcaseItems_.size();i++)
	{
		vcaseItems_[i].clearGeoData();
	}
}

void Foam::trainGeoModel::build_model
(
	label m,
	label n,
	label k
)
{
	model_ =new KNN();
	farray2d trainData;
	label cnt=0;
	for (unsigned i=0; i<vcaseItems_.size();i++)
	{
		if
		(
			vcaseItems_[i].m()!=m
			||vcaseItems_[i].n()!=n
			||vcaseItems_[i].k()!=k
		)
		{
			continue;
		}
		trainData.push_back(vcaseItems_[i].geoData_);
		idmap_[cnt]=i;
		cnt++;
	}

	if (trainData.empty())
	{
		Info<<"Warning: no data available to build KNN search tree."<<endl;
	}
	else
	{
		model_().build(trainData);
	}

	Info<<"Number of training data:"<<trainData.size()<<endl;

//	clearGeoData();
}

void Foam::trainGeoModel::build_model()
{
	model_ =new KNN();
	label ncases=vcaseItems_.size();
	farray2d trainData;
//	label cnt=0;
	std::map<word,dataItem>::iterator it;
	for (label i=0;i<ncases;i++)
	{
		word cname=caseNames_[i];
		if (vcaseItems_[i].geoData_.size()>10)
		{
			trainData.push_back(vcaseItems_[i].geoData_);
//			cnt++;
		}
	}

	if (trainData.empty())
	{
		Info<<"Warning: no data available to build KNN search tree."<<endl;
	}
	else
	{
		model_().build(trainData);
	}

	clearGeoData();
}

void Foam::trainGeoModel::readGeometryFile()
{
	word geofile=foamSample::databaseDir_+"/geom.dbs";
	std::ifstream is(geofile.c_str());
	if (!is)
	{
		FatalErrorInFunction
			<< "Error read geometry representation file, \n"<<geofile<<" cannot open for reading."
			<< exit(FatalError);
	}
	caseNames_.resize(0);
	vcaseItems_.resize(0);
	word case_name;
	word sdimension;
	word sdata;
	Info<<"reading geometric representation file, this will take some time....."<<endl;
//	label cnt=0;
	while (!is.eof())
	{
		std::getline(is,case_name,'\n');
		if (case_name.size()<3) break;
		// Info<<"reading case: "<<case_name<<endl;
		std::getline(is,sdimension,'\n');
		if (sdimension.size()<3) break;
		std::getline(is,sdata,'\n');
		if (sdata.size()<3) break;


		std::vector<word> vtmp;
		foamMap::stringtok(vtmp,sdimension," ");
		if (vtmp.size()!=4)
		{
			FatalErrorInFunction<<"Error read geometry representation file, dimension data line "<<sdimension<<" incorrect"
			<< exit(FatalError);
		}
		label npoints=atoi(vtmp[0].c_str());
		label m=atoi(vtmp[1].c_str())+1;
		label n=atoi(vtmp[2].c_str())+1;
		label k=atoi(vtmp[3].c_str());
		if (k>1)
		{
			k+=1;
		}

		std::vector<scalar> vdata(npoints);

		std::vector<word> svdata;

		split(svdata,sdata,' ');
		label nsample=svdata.size();
		if (nsample!=npoints)
		{
			FatalErrorInFunction<<"Error read geometry representation file, incorrect data item for "<<case_name
			<< exit(FatalError);
		}

		for (label i=0;i<npoints;i++)
		{
			vdata[i]=atof(svdata[i].c_str());
		}

		dataItem ditem;
		ditem.setData(m,n,k,vdata);
		vcaseItems_.push_back(ditem);
		caseNames_.push_back(case_name);
//		cnt++;
	}
	is.close();
	Info<<"number of cases read:"<<vcaseItems_.size()<<endl;

}

void Foam::trainGeoModel::save_model()
{
	// std::string geomdl = foamSample::databaseDir_ + "/geoModel.mdl";
	fileName geomdl = foamSample::databaseDir_ + "/geoModel.mdl";
	Info<< "Save geometry recognition model into database...." << endl;
	model_().saveModel(geomdl);

}

Foam::scalar Foam::trainGeoModel::getDist
(
    label i,
    label j
)
{
	const std::vector<scalar> v1=vcaseItems_[i].geoData_;
	const std::vector<scalar> v2=vcaseItems_[j].geoData_;
	scalar sum=0;

	for (unsigned t=0;t<v1.size();t++)
	{
		sum+=fabs(v1[t]-v2[t]);
	}
	return sum;
}


void Foam::trainGeoModel::getGeomSample()
{
	std::map<word,gridField>::iterator it;

	for (it=sampleGrids_.begin();it!=sampleGrids_.end();it++)
	{
		foamSample::getGeomSample(it->second);
	}
}

void Foam::trainGeoModel::get_prediction()
{
	std::map<word,gridField>::iterator it;
	for (it=sampleGrids_.begin();it!=sampleGrids_.end();it++)
	{
		gridField &grid=it->second;
		word cnts="mask";
		for (label i=0;i<grid.nx_+1;++i)
		{
			for (label j=0;j<grid.ny_+1;++j)
			{
				for (label k=0;k<grid.nzpt();++k)
				{
					label insolid=grid.insolids_(i,j,k);
					cnts+=(","+std::to_string(insolid));
				}
			}
		}

		std::vector<word> vtmp0;
		stringtok(vtmp0,cnts,",");
		std::vector<scalar> vv0(vtmp0.size()-1);
		for (unsigned s=1;s<vtmp0.size();++s)
		{
			vv0[s-1]=atof(vtmp0[s].c_str());
		}

		std::vector<int> nb(2);
		std::vector<scalar> vdists(2);
		std::vector<scalar> vquery=vv0;

		model_().search
		(
			vquery,
			nb,
			vdists
		);
		label id=nb[0];
		label id1=nb[1];

		if (idmap_.size()==vcaseItems_.size())
		{
			id=idmap_[id];
			id1=idmap_[id1];
			Info<<"Using idmap..."<<id<<" "<<id1<<endl;
		}

		vcaseItems_[id].show();
		Info<<"distance:"<<vdists[0]<<endl<<endl;

		vcaseItems_[id1].show();
		Info<<"distance:"<<vdists[1]<<endl<<endl;
	}
}

void Foam::trainGeoModel::getSolidLabels
(
	gridField& grid,
	std::vector<label>& in_solid
)
{
	Info<<"total number of grid nodes:"<<grid.size()<<endl;
	in_solid.resize(grid.size());
	label cnt=0;
	for (label i=0;i<grid.nx_+1;++i)
	{
		for (label j=0;j<grid.ny_+1;++j)
		{
			for (label k=0;k<grid.nzpt();++k)
			{
				label idx=grid.ijk(i,j,k);
				in_solid[idx]=grid.insolids_(i,j,k);
				cnt+=in_solid[idx];
			}
		}
	}
	Info<<"insolids:"<<cnt<<endl;
	return;
}

void Foam::trainGeoModel::load_model()
{
	Info<<"loading pre-trained flow-topology-recognition model..."<<endl;
	word geomdl=foamSample::databaseDir_+"/geoModel.mdl";
	model_=new KNN(geomdl);
}

