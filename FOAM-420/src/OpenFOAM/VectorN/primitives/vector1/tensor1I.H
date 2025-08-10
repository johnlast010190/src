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
    tensor1

Description
    TensorN of 1 scalars.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//- Return tensor transpose
template<>
inline tensor1 tensor1::T() const
{
    tensor1 transpose;

    transpose[0] = this->operator[](0);

    return transpose;
}

//- Assign to a sphericalTensor1
template<>
inline void tensor1::operator=(const sphericalTensor1& st)
{
    this->v_[0] = st[0];
}


//- Assign to a diagTensor1
template<>
inline void tensor1::operator=(const diagTensor1& dt)
{
    this->v_[0] = dt[0];
}

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- overide  VectorSpace functions that assume more than 1 component
inline scalar magSqr(const tensor1& v1)
{
    return magSqr(v1.v_[0]);
}

inline scalar cmptMax(const tensor1& v1)
{
    return v1.v_[0];
}

inline scalar cmptMin(const tensor1& v1)
{
    return v1.v_[0];
}

inline scalar cmptSum(const tensor1& v1)
{
    return v1.v_[0];
}

inline scalar cmptProduct(const tensor1& v1)
{
    return v1.v_[0];
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two tensors
inline tensor1
operator&(const tensor1& t1, const tensor1& t2)
{
    tensor1 result;

    result[0] = t1[0]*t2[0];


    return result;
}


//- Inner-product between a diagonal tensors and a tensor
inline tensor1
operator&(const diagTensor1& dt1, const tensor1& t2)
{
    tensor1 result;

    result[0] = dt1[0]*t2[0];

    return result;
}

//- Inner-product between a tensor and diagonal tensor
inline tensor1
operator&(const tensor1& t1, const diagTensor1& dt2)
{
    tensor1 result;

    result[0] = t1[0]*dt2[0];

    return result;
}


//- Inner-product between a spherical tensor and a tensor
inline tensor1
operator&(const sphericalTensor1& st1, const tensor1& t2)
{
    tensor1 result;

    result[0] = st1[0]*t2[0];

    return result;
}

//- Inner-product between a tensor and spherical tensor
inline tensor1
operator&(const tensor1& t1, const sphericalTensor1& st2)
{
    tensor1 result;

    result[0] = t1[0]*st2[0];

    return result;
}


//- Inner-product between a tensor and a vector
inline vector1
operator&(const tensor1& t, const vector1& v)
{
    vector1 result;

    result[0] = t[0]*v[0];


    return result;
}


//- Inner-product between a vector and a tensor
inline vector1
operator&(const vector1& v, const tensor1& t)
{
    vector1 result;

    result[0] = v[0]*t[0];

    return result;
}

//- Inner-product between a tensor and a vector
inline vector
operator&(const tensor1& t, const vector& v)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and tensor1 "
        << "vector has 3 components while tensor1 has 1"
        << exit(FatalError);

    return result;
}


//- Inner-product between a vector and a tensor
inline vector
operator&(const vector v, const tensor1& t)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and tensor1 "
        << "vector has 3 components while tensor1 has 1"
        << exit(FatalError);

    return result;
}


//- Outer-product between two vectors
inline tensor1
operator*(const vector1& v1, const vector1& v2)
{
    tensor1 result;

    result[0] = v1[0]*v2[0];

    return result;
}


//- Return the determinant of a tensor
inline scalar det(const tensor1& t)
{
    return
    (
        t[0]
    );
}

//- Return the inverse of a tensor given the determinant
inline tensor1 inv(const tensor1& t)
{
    tensor1 cofactor;

    cofactor[0] = 1/t[0];


    return cofactor;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
