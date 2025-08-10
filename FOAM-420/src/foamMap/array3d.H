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
Description:
    STL-based two-dimensional container for data storage and manipulation

\*---------------------------------------------------------------------------*/

#pragma once
#include <map>
#include <vector>
#include <string>
#include "global/argList/argList.H"
#include "array2d.H"

namespace Foam
{

template<class type>
class Array3d
{
    public:
        std::vector<Array2d<type>> data_;
        Array3d(){;}
        void resize(label n1)
        {
            data_.resize(n1);
        }

        void setSize(label nx, label ny, label nz)
        {
            data_.resize(nx);
            for (label i=0;i<nx;i++)
            {
                data_[i].resize(ny);
            }
            for (label i=0;i<nx;i++)
            {
                for (label j=0;j<ny;j++)
                {
                    data_[i][j].resize(nz);
                }
            }
        }

        Array2d<type> &operator[](label i)
        {
            return data_[i];
        }

        const Array2d<type> &operator[](label i)const
        {
            return data_[i];
        }

        void push_back(Array2d<type> &data)
        {
            data_.push_back(data);
        }

        label dim1()
        {
            return data_.size();
        }

        label dim2()
        {
            return data_[0].size();
        }

        label dim3()
        {
            return data_[0][0].size();
        }

        type &operator()(label i, label j, label k)
        {
            return data_[i][j][k];
        }

        void setval
        (
            label i,
            label j,
            label k,
            const type &val
        )
        {
            data_[i][j][k]=val;
        }

        type getval
        (
            label i,
            label j,
            label k
        )
        {
           return data_[i][j][k];
        }

        label size()
        {
            return data_.size();
        }

        type operator()(label i, label j, label k) const
        {
            return data_[i][j][k];
        }

        inline void operator=(const type &val);//assign a value to the whole array
        void operator+=(const type &val);
        void operator-=(const type &val);
        void operator*=(const type &val);
        void operator/=(const type &val);

        bool empty()
        {
            return data_.empty();
        }

        ~Array3d(){;}
};

template<class type>
void Array3d<type>::operator=(const type &val)
{
    for (unsigned i=0;i<data_.size();i++)
    {
        for (unsigned j=0;j<data_[i].size();j++)
        {
            for (unsigned k=0;k<data_[i][j].size();k++)
            {
                data_[i][j][k]=val;
            }
        }
    }
    return;
}


template<class type>
class Tri
{
    public:
        type a;
        type b;
        type c;
        Tri (const type &a0, const type &b0, const type &c0)
        {
            a=a0;b=b0;c=c0;
        }
        Tri(){;}

        ~Tri(){;}
};

}//Foam
