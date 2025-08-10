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

#include "array3d.H"

using namespace Foam;

/*
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
*/

template<class type>
void Array3d<type>::operator+=(const type &val)
{
    for (unsigned i=0;i<data_.size();i++)
    {
        for (unsigned j=0;j<data_[i].size();j++)
        {
            for (unsigned k=0;k<data_[i][j].size();k++)
            {
                data_[i][j][k]+=val;
            }
        }
    }
    return;
}

template<class type>
void Array3d<type>::operator-=(const type &val)
{
    for (unsigned i=0;i<data_.size();i++)
    {
        for (unsigned j=0;j<data_[i].size();j++)
        {
            for (unsigned k=0;k<data_[i][j].size();k++)
            {
                data_[i][j][k]-=val;
            }
        }
    }
    return;
}

template<class type>
void Array3d<type>::operator*=(const type &val)
{
    for (unsigned i=0;i<data_.size();i++)
    {
        for (unsigned j=0;j<data_[i].size();j++)
        {
            for (unsigned k=0;k<data_[i][j].size();k++)
            {
                data_[i][j][k]*=val;
            }
        }
    }
    return;
}

template<class type>
void Array3d<type>::operator/=(const type &val)
{
    if (fabs(val)<1.0e-30)
    {
        return;
    }

    for (unsigned i=0;i<data_.size();i++)
    {
        for (unsigned j=0;j<data_[i].size();j++)
        {
            for (unsigned k=0;k<data_[i][j].size();k++)
            {
                data_[i][j][k]/=val;
            }
        }
    }
    return;
}













