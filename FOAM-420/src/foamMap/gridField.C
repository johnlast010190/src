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
    (c) 2019 Esi Ltd.


Class
    Foam::gridField

Group
    grpfoamMap

Description
   Class for sampling flow fields from unstructured grid into uniform grid

SourceFile
    gridField.C

\*---------------------------------------------------------------------------*/

#include "gridField.H"
#include "foamMap.H"

using namespace Foam;

scalar gridField::rhoRef_;
vector gridField::Uref_;

label gridField::inlist
(
    const word &key,
    const std::vector<word> &vals
)
{
    label sz=vals.size();
    for (label i=0;i<sz;i++)
    {
        if (vals[i]==key)
        {
            return i;
        }
    }
    return -1;
}


void gridField::setDivision
(
    label nx,
    label ny,
    label nz
)
{
    //called after set field names
    nx_=nx;
    ny_=ny;
    nz_=nz;
}

void gridField::setDelta()
{
    dx_=(xmax()-xmin())/scalar(nx_);
    dy_=(ymax()-ymin())/scalar(ny_);
    dz_=(zmax()-zmin())/scalar(nz_);
}

void gridField::setSampleNames
(
    const word &scalars, //names are separated by ,
    const word &vectors
)
{
    foamMap::stringtok
    (
        scalarSamples_,
        scalars,
        ","
    );

    if (vectors!="none")
    {
        foamMap::stringtok
        (
            vectorSamples_,
            vectors,
            ","
        );
    }

}

void gridField::generateNodes()
{
    //generate node coordinates and store in gridNodes
    dx_=(xmax()-xmin())/scalar(nx_);
    dy_=(ymax()-ymin())/scalar(ny_);
    dz_=(zmax()-zmin())/scalar(nzpt());

    if (dx_ < SMALL || dy_ < SMALL || dz_ < SMALL)
    {
        FatalErrorInFunction
            << "One of more of the dx,dy,dz are too small,please check \nwhether bounding box is set properly."
            <<"dx="<<dx_<<",dy="<<dy_<<",dz="<<dz_
            << exit(FatalError);
    }

    if (nz_>1)
    {
        gridNodes_.setSize(nx_+1,ny_+1,nz_+1);
    }
    else
    {
        gridNodes_.setSize(nx_+1,ny_+1,1);
    }

    for (label i=0;i<nx_+1;++i)
    {
        scalar xi=xmin()+i*dx_;
        for (label j=0;j<ny_+1;++j)
        {
            scalar yj=ymin()+j*dy_;
            for (label k=0;k<nzpt();++k)
            {
                scalar zk=zmin()+k*dz_;
                point node(xi,yj,zk);
                gridNodes_(i,j,k)=node;
            }
        }
    }

    return;
}

void gridField::getBbox(const boundBox &bbox)
{
    xmin_=bbox.min().x();
    ymin_=bbox.min().y();
    zmin_=bbox.min().z();
    xmax_=bbox.max().x();
    ymax_=bbox.max().y();
    zmax_=bbox.max().z();

}

void gridField::getBbox
(
    const word  &val
)
{
    //boundBox=xmin,ymin,zmin,xmax,yax,zmax
    std::vector<word> vtmp;
    foamMap::stringtok(vtmp,val,",");
    if (vtmp.size()!=6)
    {
        FatalErrorInFunction
            << "Incorrect data:"<<val
            << exit(FatalError);
    }
    xmin_=atof(vtmp[0].c_str());
    ymin_=atof(vtmp[1].c_str());
    zmin_=atof(vtmp[2].c_str());
    xmax_=atof(vtmp[3].c_str());
    ymax_=atof(vtmp[4].c_str());
    zmax_=atof(vtmp[5].c_str());
}


gridField::gridField()
{
    nx_=20;
    ny_=20;
    nz_=1;
    dx_=0;
    dy_=0;
    dz_=0;
    interp_=false;
    interp2_=false;
    containSolidObject_=true;
}

gridField::gridField
(
    const boundBox  &bbox,
    label nx,
    label ny,
    label nz,
    bool contain_obj
)
{
    getBbox(bbox);
    setDivision
    (
        nx,
        ny,
        nz
    );
    setDelta();
    interp_=false;
    containSolidObject_=contain_obj;
}

gridField::gridField
(
    const word  &bbox,
    label nx,
    label ny,
    label nz,
    bool contain_obj
)
{
    getBbox(bbox);
    setDivision
    (
        nx,
        ny,
        nz
    );
    setDelta();
    interp_=false;
    containSolidObject_=contain_obj;
}

void gridField::setSampleNames
(
    const HashSet<word>& scalars,
    const HashSet<word>& vectors
)
{
    forAllConstIter(HashSet<word>, scalars, iter)
    {
        scalarSamples_.push_back(iter());
    }

    forAllConstIter(HashSet<word>, vectors, iter)
    {
        vectorSamples_.push_back(iter());
    }
    setFieldSize();
}

void gridField::setFieldSize()
{
    if (gridNodes_.size()==0)
    {
        gridNodes_.setSize(nx_+1,ny_+1,nzpt());
    }

    insolids_.setSize(nx_+1,ny_+1,nzpt());

    if (scalarSamples_.size()>=1)
    {
        gridFields_.setSize(nx_+1,ny_+1,nzpt());
    }

    if (vectorSamples_.size()>=1)
    {
        gridVectors_.setSize(nx_+1,ny_+1,nzpt());
    }

    dists_.setSize(nx_+1,ny_+1,nzpt());

    return;
}

