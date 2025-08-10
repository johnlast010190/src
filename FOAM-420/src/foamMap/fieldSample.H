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
    Foam::fieldSample

Group
    grpfoamMap

Description
    A class designed to sample scalar and vector fields from unstructured mesh

SourceFile
    fieldSample .C

\*---------------------------------------------------------------------------*/



#pragma once

#include "gridField.H"
#include "readFields.H"
namespace Foam
{
class fieldSample : public gridField
{
    public:
        KNN bodyMap_;
        KNN internalMap_;
        word bodySurfaceType_;
        word bodySurfaceName_;
        bool byName_;
        bool byNameList_;
        bool byPhysicalType_;
        const volScalarField *field_;
        DynamicList<std::pair<label,label>> bodyFaces_;
        Array3d<scalar> gridField_;


        void write_samples
        (
            const word &fldname
        );

        void init();
        void sampleField
        (
            const volScalarField &fld,
            const fvMesh *mesh,
            bool setLevel=false //set distance inside solid
        );

        inline void getSolidField
        (
            label i,
            label j,
            label k,
            const std::vector<label> &nb,
            const std::vector<scalar> &dist,
            bool setLevel=false
        );


        fieldSample
        (
            label nx,
            label ny,
            label nz,
            scalar x0,
            scalar y0,
            scalar z0,
            scalar x1,
            scalar y1,
            scalar z1
        );

        inline void getGridField
        (
            label i,
            label j,
            label k,
            const std::vector<label> &nb,
            const std::vector<scalar> &dist
        );

        void construct_field_knn
        (
            const fvMesh *mesh
        );

        void construct_body_knn
        (
            const fvMesh *mesh
        );

        fieldSample();
        virtual ~fieldSample(){;}

};

//get the scalar field at a nodal point by kdTree search
void fieldSample::getGridField
(
    label i,
    label j,
    label k,
    const std::vector<label> &nb,
    const std::vector<scalar> &dist
)
{
    label ic=nb[0];
    scalar node=0;
    if (!interp_)
    {
        node=(*field_)[ic];
    }
    else
    {
        scalar sum=0;
        scalar sum1=1.0e-30;
        label sz=nb.size();
        for (label ik=0;ik<sz;ik++)
        {
            label ick=nb[ik];
            scalar tmp=1.0/(dist[ik]+1.0e-30);
            sum+=((*field_)[ick]*tmp);
            sum1+=tmp;
        }
        node=sum/sum1;
    }

    gridField_(i,j,k)=node;
    return;
}

void fieldSample::getSolidField
(
    label i,
    label j,
    label k,
    const std::vector<label> &nb,
    const std::vector<scalar> &dist,
    bool setLevel
)
{
    if (!setLevel)
    {
        gridField_(i,j,k)=0;
    }
    else
    {
        gridField_(i,j,k)=-dist[0];
    }
}

}//Foam
