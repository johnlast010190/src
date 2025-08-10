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
    (c) 2010 Ivor Clifford

Type
    diagTensor1

Description
    DiagTensorN of 1 scalars.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Assign to a SphericalTensorN
template<>
inline void diagTensor1::operator=(const sphericalTensor1& st)
{
    this->v_[0] = st[0];
}

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two diagonal tensors
inline diagTensor1
operator&(const diagTensor1& dt1, const diagTensor1& dt2)
{
    diagTensor1 result;

    result[0] = dt1[0]*dt2[0];

    return result;
}


//- Inner-product between a spherical tensor and a diagonal tensor
inline diagTensor1
operator&(const sphericalTensor1& st1, const diagTensor1& dt2)
{
    diagTensor1 result;

    result[0] = st1[0]*dt2[0];

    return result;
}

//- Inner-product between a diagonal tensor and a spherical tensor
inline diagTensor1
operator&(const diagTensor1& dt1, const sphericalTensor1& st2)
{
    diagTensor1 result;

    result[0] = dt1[0]*st2[0];

    return result;
}


//- Inner-product between a diagonal tensor and a vector
inline vector1
operator&(const diagTensor1& dt, const vector1& v)
{
    vector1 result;

    result[0] = dt[0]*v[0];

    return result;
}


//- Inner-product between a vector and a diagonal tensor
inline vector1
operator&(const vector1& v, const diagTensor1& dt)
{
    vector1 result;

    result[0] = v[0]*dt[0];

    return result;
}

//- Inner-product between a diagonal tensor and a vector
inline vector
operator&(const diagTensor1& dt, const vector& v)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and diagTensor1 "
        << "vector has 3 components while diagTensor1 has 1"
        << exit(FatalError);

    return result;
}


//- Inner-product between a vector and a diagonal tensor
inline vector
operator&(const vector& v, const diagTensor1& dt)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and diagTensor1 "
        << "vector has 3 components while diagTensor1 has 1"
        << exit(FatalError);

    return result;
}

//- Division of a scalar by a diagTensor1
inline diagTensor1 operator/(const scalar s, const diagTensor1& dt)
{
    diagTensor1 result;

    result[0] = s/dt[0];

    return result;
}


//- Return the inverse of a diagTensor1
inline diagTensor1 inv(const diagTensor1& dt)
{
    diagTensor1 result;

    result[0] = 1.0/dt[0];

    return result;
}


//- Return the component sum
inline scalar sum(const diagTensor1& dt)
{
    return dt[0];
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
