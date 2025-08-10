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
    (c) 2009-2021 Esi Ltd.

Description
    Database item for AI training and prediction

\*---------------------------------------------------------------------------*/

#include "dataItem.H"

Foam::dataItem::dataItem ()
{
	caseName_="block1";
	caseType_="2dcar";
	Uref_.x()=30;
	Uref_.y()=0;
	Uref_.z()=0;
	rhoRef_=1.205;
	identity_=0;
}

Foam::scalar Foam::dataItem::getDist
(
	const std::vector<scalar>& vdata
) const
{
	scalar sum=0;
	for (unsigned i=0;i<vdata.size();i++)
	{
		sum+=fabs(vdata[i]-geoData_[i]);
	}
	return sum;

}

Foam::label Foam::dataItem::numSolids() const
{
	scalar sum=0;
	for (unsigned i=0;i<geoData_.size();i++)
	{
		sum+=geoData_[i];
	}
	return label(sum+0.1);
}

void Foam::dataItem::show() const
{
	Info<<"case name:"<<caseName_<<endl;
	Info<<"Case type:"<<caseType_<<endl;
	Info<<"Case path:"<<caseDir_<<endl;
	Info<<"Sample in x,y,z:"<<m_<<" "<<n_<<" "<<k_<<endl;
	Info<<"Bound box:"<<xmin_<<" "<<ymin_<<" "<<zmin_<<" "<<xmax_<<" "<<ymax_<<" "<<zmax_<<endl;
	Info<<"identity:"<<identity_<<endl;
	Info<<"Ref U:"<<Uref_<<endl;
	Info<<"Ref rho:"<<rhoRef_<<endl;
}

Foam::scalar Foam::dataItem::gdataSum() const
{
	scalar sum=0;
	for (unsigned i=0;i<geoData_.size();i++)
	{
		sum+=geoData_[i];
	}
	return sum;
}

void Foam::dataItem::setData
(
	const word &data
)
{
	std::vector<word> vtmp;

	split(vtmp,data,',');

	//mask,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
	if (vtmp.size()<2)
	{
		return;
	}

    geoData_.resize(m_*n_*k_,0);
    for (unsigned ik=1;ik<vtmp.size();ik++)
	{
		std::vector<word> vtmp1;
        split(vtmp1,vtmp[ik],':');
        label i=atoi(vtmp1[0].c_str());
        label j=atoi(vtmp1[1].c_str());
        label k=atoi(vtmp1[2].c_str());
        label ijk=i*n_*k_+j*k_+k;
        geoData_[ijk]=1;
	}
}

void Foam::dataItem::setDataInfo(const word& dat)
{
	std::vector<word> vtmp;
	foamMap::stringtok(vtmp,dat,",");
	//2dcar,/home/jianbo/2dcars/data_sample/27,block1,
	//256,256,1,
	// -1.500000,-0.317000,-0.007000,
	// 4.200000,1.800000,0.010000
	if (vtmp.size()!=12)
	{
		FatalErrorInFunction<<"Error set data item info, incorrect data format."
		<< exit(FatalError);
	}
	caseType_=vtmp[0];
	caseDir_=vtmp[1];
	caseName_=vtmp[2];
	m_=atoi(vtmp[3].c_str());
	n_=atoi(vtmp[4].c_str());
	k_=atoi(vtmp[5].c_str());
	xmin_=atof(vtmp[6].c_str());
	ymin_=atof(vtmp[7].c_str());
	zmin_=atof(vtmp[8].c_str());
	xmax_=atof(vtmp[9].c_str());
	ymax_=atof(vtmp[10].c_str());
	zmax_=atof(vtmp[11].c_str());
}

Foam::dataItem::dataItem(const word &dat)
{
	std::vector<word> vtmp;
	foamMap::stringtok(vtmp,dat,"=");
	caseType_="2dcar";
	identity_=0;

	if (vtmp.size()==2)
	{
		caseName_=vtmp[0];
		std::vector<word> vtmp1;
		foamMap::stringtok(vtmp1,vtmp[1]," ");
		//1.205 85.2512 8.41421 0 -1 -1 0 1 1 0.1
		rhoRef_=atof(vtmp1[0].c_str());
		scalar ux=atof(vtmp1[1].c_str());
		scalar uy=atof(vtmp1[2].c_str());
		scalar uz=atof(vtmp1[3].c_str());
		Uref_[0]=ux;
		Uref_[1]=uy;
		Uref_[2]=uz;
		xmin_=atof(vtmp1[4].c_str());
		ymin_=atof(vtmp1[5].c_str());
		zmin_=atof(vtmp1[6].c_str());
		xmax_=atof(vtmp1[7].c_str());
		ymax_=atof(vtmp1[8].c_str());
		zmax_=atof(vtmp1[9].c_str());
	}
}


void Foam::dataItem::setData
(
	label m,
	label n,
	label k,
	const std::vector<scalar> &data
)
{
	m_=m;
	n_=n;
	k_=k;
	geoData_=data;
}




