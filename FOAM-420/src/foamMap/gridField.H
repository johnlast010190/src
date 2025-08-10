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
    A class designed to create a symmetric scalar field from an unsymmetric field

SourceFile
    gridField.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "array3d.H"
#include "db/Time/Time.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"

namespace Foam
{
typedef DynamicList<scalar> nodeField;
typedef DynamicList<vector> nodeVectors;

class gridField
{
    protected:
        scalar xmin_;
        scalar xmax_;
        scalar ymin_;
        scalar ymax_;
        scalar zmin_;
        scalar zmax_;

    public:
        static scalar rhoRef_;
        static vector Uref_;
        bool containSolidObject_;

        label nx_;
        label ny_;
        label nz_;
        scalar dx_;
        scalar dy_;
        scalar dz_;
        word name_;
        bool interp_;
        bool interp2_;
        std::vector<word> scalarSamples_;
        std::vector<word> vectorSamples_;
        Array3d<label> insolids_;
        Array3d<nodeField> gridFields_;

        Array3d<scalar> dists_;
        Array3d<point> gridNodes_;
        Array3d<nodeVectors> gridVectors_;


        inline label scalarId(const word &scname);
        inline label vectorId(const word &vecname);

        void generateNodes();

        inline bool inblock
        (
            scalar x,
            scalar y,
            scalar z
        );

        static label inlist(const word &key,const std::vector<word> &vals);

        inline point gridPoint(label i,label j,label k); //coordinate of a grid point

        void setScalar
        (
            const word &scName,
            label i,
            label j,
            label k,
            scalar val
        )
        {
            label isc=scalarId(scName);
            gridFields_(i,j,k)[isc]=val;
        }


        void setVector
        (
            const word &vecName,
            label i,
            label j,
            label k,
            point val
        )
        {
            label iv=vectorId(vecName);
            gridVectors_(i,j,k)[iv]=val;
        }

        inline scalar p(label i,label j,label k,label ip);
        inline point U(label i,label j,label k,label iv);

        inline scalar getScalar
        (
            word scName,
            scalar x,
            scalar y,
            scalar z
        );

        inline point getVector
        (
            word vecName,
            scalar x,
            scalar y,
            scalar z
        );

        scalar xmin()
        {
            return xmin_;
        }
        scalar xmax()
        {
            return xmax_;
        }
        word name()
        {
            return name_;
        }

        label nzpt()
        {
            if (nz_==1)
            {
                return nz_;
            }
            else
            {
                return nz_+1;
            }
        }

        scalar ymin()
        {
            return ymin_;
        }
        scalar ymax()
        {
            return ymax_;
        }
        scalar zmin()
        {
            return zmin_;
        }

        scalar zmax()
        {
            return zmax_;
        }

        label nnodes()
        {
            if (nz_>1)
            {
                return (nx_+1)*(ny_+1)*(nz_+1);
            }
            else
            {
                return (nx_+1)*(ny_+1);
            }
        }

        label dimension()
        {
            if (nz_==1)
            {
                return 2;
            }
            else
            {
                return 3;
            }
        }

        scalar xbar(scalar x)
        {
            return (x-xmin())/(xmax()-xmin());
        }

        scalar ybar(scalar y)
        {
            return (y-ymin())/(ymax()-ymin());
        }

        scalar zbar(scalar z)
        {
            return (z-zmin())/(zmax()-zmin());
        }

        label cell_ijk(label i,label j,label k)
        {
            if (nz_==1)
            {
                return i+j*nx_;
            }
            else
            {
                return i+j*nx_+k*nx_*ny_;
            }
        }

        label ijk(label i,label j,label k)
        {
            if (nz_==1)
            {
                return i+j*(nx_+1);
            }
            else
            {
                return i+j*(nx_+1)+k*(nx_+1)*(ny_+1);
            }
        }

        void getBbox(const word  &val);
        void getBbox(const boundBox &val);

        void setDelta();
        void setSampleNames
        (
            const HashSet<word>& scalars,
            const HashSet<word>& vectors
        );

        void setSampleNames
        (
            const word &scalars, //names are separated by ,
            const word &vectors="none"
        );

        void setDivision
        (
            label nx,
            label ny,
            label nz=1
        );

        inline void getGridId
        (
            label &i,
            label &j,
            label &k,
            scalar x,
            scalar y,
            scalar z
        );

        label size()
        {
            return (nx_+1)*(ny_+1)*nzpt();
        }

        void setFieldSize();

        gridField();

        gridField
        (
            const word  &bbox, //bounding box xmin,ymin,zmin,xmax,ymax,zmax
            label nx,
            label ny,
            label nz=1,
            bool contain_obj=true
        );

        gridField
        (
            const boundBox  &bbox, //bounding box xmin,ymin,zmin,xmax,ymax,zmax
            label nx,
            label ny,
            label nz=1,
            bool contain_obj=true
        );

        virtual ~gridField(){;}
};

point gridField::gridPoint(label i,label j,label k)
{
    point pt(xmin()+dx_*i,ymin()+dy_*j,zmin()+dz_*k);
    return pt;
}

scalar gridField::p(label i,label j,label k,label ip)
{
    return gridFields_(i,j,k)[ip];

}

point gridField::U(label i,label j,label k,label iv)
{
    return gridVectors_(i,j,k)[iv];
}

void gridField::getGridId
(
    label &i,
    label &j,
    label &k,
    scalar x,
    scalar y,
    scalar z
)
{
    i=label((x-xmin())/dx_);
    j=label((y-ymin())/dy_);
    k=0;
    if (dimension()==3)
    {
        k=label((z-zmin())/dz_);
    }
}

label gridField::scalarId
(
    const word &scname
)
{
    for (unsigned isc=0;isc<scalarSamples_.size();isc++)
    {
        if (scalarSamples_[isc]==scname)
        {
            return isc;
        }
    }
    return -1;
}

point gridField::getVector
(
    word vecName,
    scalar x,
    scalar y,
    scalar z
)
{
    label iv=vectorId(vecName);
    if (iv<0)
    {
        FatalErrorInFunction<<"Error get vector from grid, "<<vecName<<" not found.\n"
        << exit(FatalError);
    }

    label i,j,k;
    getGridId(i,j,k,x,y,z);

    point v1=gridVectors_(i,j,k)[iv];
    if (dimension()==3)
    {
        return v1;
    }
    if (!interp2_)
    {
        return v1;
    }
    label i1=i+1;
    if (i1>nx_) i1=nx_;
    label j1=j+1;
    if (j1>ny_) j1=ny_;

    point v2=gridVectors_(i1,j,k)[iv];
    point v3=gridVectors_(i1,j1,k)[iv];
    point v4=gridVectors_(i,j1,k)[iv];

    point p0(x,y,z);
    point p1(xmin()+dx_*i,ymin()+dy_*j,z);
    point p2(xmin()+dx_*i1,ymin()+dy_*j,z);
    point p3(xmin()+dx_*i1,ymin()+dy_*j1,z);
    point p4(xmin()+dx_*i,ymin()+dy_*j1,z);
    point p01=p0-p1;
    point p02=p0-p2;
    point p03=p0-p3;
    point p04=p0-p4;
    scalar s01=mag(p01)+1.0e-20;
    scalar s02=mag(p02)+1.0e-20;
    scalar s03=mag(p03)+1.0e-20;
    scalar s04=mag(p04)+1.0e-20;
    point val=v1/s01+v2/s02+v3/s03+v4/s04;
    scalar ds=1/s01+1/s02+1/s03+1/s04;
    return val/ds;
}

scalar gridField::getScalar
(
    word scName,
    scalar x,
    scalar y,
    scalar z
)
{
    label isc=scalarId(scName);
    if (isc<0)
    {
        FatalErrorInFunction<<"Error get scalar from grid, "<<scName<<" not found.\n"
        << exit(FatalError);
    }
    label i,j,k;
    getGridId(i,j,k,x,y,z);
    scalar s1=gridFields_(i,j,k)[isc];
    if (dimension()==3)
    {
        return s1;
    }

    if (!interp2_)
    {
        return s1;
    }

    label i1=i+1;
    if (i1>nx_) i1=nx_;
    label j1=j+1;
    if (j1>ny_) j1=ny_;

    scalar s2=gridFields_(i1,j,k)[isc];
    scalar s3=gridFields_(i1,j1,k)[isc];
    scalar s4=gridFields_(i,j1,k)[isc];
    point p0(x,y,z);
    point p1(xmin()+dx_*i,ymin()+dy_*j,z);
    point p2(xmin()+dx_*i1,ymin()+dy_*j,z);
    point p3(xmin()+dx_*i1,ymin()+dy_*j1,z);
    point p4(xmin()+dx_*i,ymin()+dy_*j1,z);
    point p01=p0-p1;
    point p02=p0-p2;
    point p03=p0-p3;
    point p04=p0-p4;
    scalar s01=mag(p01)+1.0e-20;
    scalar s02=mag(p02)+1.0e-20;
    scalar s03=mag(p03)+1.0e-20;
    scalar s04=mag(p04)+1.0e-20;
    scalar val=s1/s01+s2/s02+s3/s03+s4/s04;
    scalar ds=1/s01+1/s02+1/s03+1/s04;
    return val/ds;
}


label gridField::vectorId
(
    const word &vecname
)
{
    for (unsigned iv=0;iv<vectorSamples_.size();iv++)
    {
        if (vectorSamples_[iv]==vecname)
        {
            return iv;
        }
    }
    return -1;
}

bool gridField::inblock
(
    scalar x,
    scalar y,
    scalar z
)
{
    if
    (
         x>=xmin() && x<=xmax()
       &&y>=ymin() && y<=ymax()
       &&z>=zmin() && z<=zmax()
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}

}//Foam
