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
    sphericalTensor2

Description
    SphericalTensorN of 2 scalars.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two spherical tensors
inline sphericalTensor2
operator&(const sphericalTensor2& st1, const sphericalTensor2& st2)
{
    return sphericalTensor2(st1[0]*st2[0]);
}


//- Inner-product between a spherical tensor and a vector
inline vector2
operator&(const sphericalTensor2& st, const vector2& v)
{
    vector2 result;

    result[0] = st[0]*v[0];
    result[1] = st[0]*v[1];

    return result;
}


//- Inner-product between a vector and a spherical tensor
inline vector2
operator&(const vector2& v, const sphericalTensor2& st)
{
    vector2 result;

    result[0] = v[0]*st[0];
    result[1] = v[1]*st[0];

    return result;
}

//- Inner-product between a spherical tensor and a vector
inline vector
operator&(const sphericalTensor2& st, const vector& v)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and sphericalTensor2 "
        << "vector has 3 components while sphericalTensor2 has 2"
        << exit(FatalError);

    return result;
}


//- Inner-product between a vector and a spherical tensor
inline vector
operator&(const vector& v, const sphericalTensor2& st)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and sphericalTensor2 "
        << "vector has 3 components while spericalTensor2 has 2"
        << exit(FatalError);

    return result;
}

//- Return the inverse of a tensor given the determinant
inline sphericalTensor2 inv(const sphericalTensor2& st)
{
    return sphericalTensor2(1/st[0]);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
