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
    sphericalTensor1

Description
    SphericalTensorN of 1 scalars.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Inner-product between two spherical tensors
inline sphericalTensor1
operator&(const sphericalTensor1& st1, const sphericalTensor1& st2)
{
    return sphericalTensor1(st1[0]*st2[0]);
}


//- Inner-product between a spherical tensor and a vector
inline vector1
operator&(const sphericalTensor1& st, const vector1& v)
{
    vector1 result;

    result[0] = st[0]*v[0];

    return result;
}


//- Inner-product between a vector and a spherical tensor
inline vector1
operator&(const vector1& v, const sphericalTensor1& st)
{
    vector1 result;

    result[0] = v[0]*st[0];

    return result;
}

//- Inner-product between a spherical tensor and a vector
inline vector
operator&(const sphericalTensor1& st, const vector& v)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and sphericalTensor1 "
        << "vector has 3 components while sphericalTensor1 has 1"
        << exit(FatalError);

    return result;
}

//- Inner-product between a vector and a spherical tensor
inline vector
operator&(const vector& v, const sphericalTensor1& st)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatible types of vector and sphericalTensor1 "
        << "vector has 3 components while sphericalTensor1 has 1"
        << exit(FatalError);

    return result;
}

//- Return the inverse of a tensor given the determinant
inline sphericalTensor1 inv(const sphericalTensor1& st)
{
    return sphericalTensor1(1/st[0]);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
