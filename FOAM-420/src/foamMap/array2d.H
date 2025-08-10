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
namespace Foam
{

template<class type>
class Array2d
{
    public:
        std::vector<std::vector<type>> data_;
        Array2d(){;}

        void setSize(label nx, label ny)
        {
            data_.resize(nx);
            for (label i=0;i<nx;i++)
            {
                data_[i].resize(ny);
            }
        }

        type &operator()(label i, label j)
        {
            return data_[i][j];
        }
        void resize(label m)
        {
            data_.resize(m);
        }

        unsigned size() {return data_.size();}

        std::vector<type> &operator[](label i)
        {
            return data_[i];
        }

        const std::vector<type> &operator[](label i)
        const
        {
            return data_[i];
        }

        label dim1()
        {
            return data_.size();
        }

        label dim2()
        {
            return data_[0].size();
        }

        void resize(label m,label n)
        {
            data_.resize(m);
            for (label i=0;i<m;i++)
            {
                data_[i].resize(n);
            }
        }

        Array2d(label m,label n)
        {
            resize(m,n);
        }
        void setval(label i, label j, const type &val)
        {
            data_[i][j]=val;
        }


         type operator()(label i, label j) const
        {
            return data_[i][j];
        }

//        void operator=(const type &val);//assign a value to the whole array
        void operator+=(const type &val);
        void operator-=(const type &val);
        void operator*=(const type &val);
        void operator/=(const type &val);

        bool empty() const
        {
            return data_.empty();
        }

        void operator=(type val)
        {
            for (unsigned i=0;i<data_.size();i++)
            {
                for (unsigned j=0;j<data_[i].size();j++)
                {
                    data_[i][j]=val;
                }
            }
        }
        ~Array2d(){;}
};


}//Foam
