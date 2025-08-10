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
    tensor4

Description
    TensorN of 4 scalars.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


//- Return the determinant of a tensor
inline scalar det(const tensor4& t)
{
    return
    (
        //- G.Karpouzas - code optimization,
        //  Reducing operations by grouping them
        t[0]*
            (
                t[5]*(t[10]*t[15]-t[11]*t[14])
              + t[6]*(t[11]*t[13]-t[9]*t[15])
              + t[7]*(t[9]*t[14]-t[10]*t[13])
            )
       +t[1]*
            (
                t[4]*(t[11]*t[14]-t[10]*t[15])
              + t[6]*(t[8]*t[15]-t[11]*t[12])
              + t[7]*(t[10]*t[12]-t[8]*t[14])
            )
       +t[2]*
            (
                t[4]*(t[9]*t[15]-t[11]*t[13])
              + t[5]*(t[11]*t[12]-t[8]*t[15])
              + t[7]*(t[8]*t[13]-t[9]*t[12])
            )
       +t[3]*
            (
                t[4]*(t[10]*t[13]-t[9]*t[14])
              + t[5]*(t[8]*t[14]-t[10]*t[12])
              + t[6]*(t[9]*t[12]-t[8]*t[13])
            )
    );
}


//- Return the inverse of a tensor given the determinant
inline tensor4 inv(const tensor4& t)
{
    //- G.Karpouzas - code optimization,
    //  Reducing operations by grouping them
    //  Compute C_transpose --> A^-1 = C_transpose/det;
    tensor4 cofactor;

    cofactor[0] =
           t[5]*(t[10]*t[15] - t[11]*t[14])
         + t[6]*(t[11]*t[13] - t[9]*t[15])
         + t[7]*(t[9]*t[14] - t[10]*t[13]);

    cofactor[1] =
           t[1]*(t[11]*t[14] - t[10]*t[15])
         + t[2]*(t[9]*t[15] - t[11]*t[13])
         + t[3]*(t[10]*t[13] - t[9]*t[14]);

    cofactor[2] =
           t[1]*(t[6]*t[15]-t[7]*t[14])
         + t[2]*(t[7]*t[13]-t[5]*t[15])
         + t[3]*(t[5]*t[14]-t[6]*t[13]);

    cofactor[3] =
           t[1]*(t[7]*t[10]-t[6]*t[11])
         + t[2]*(t[5]*t[11]-t[7]*t[9])
         + t[3]*(t[6]*t[9]-t[5]*t[10]);

    cofactor[4] =
           t[4]*(t[11]*t[14]-t[10]*t[15])
         + t[6]*(t[8]*t[15]-t[11]*t[12])
         + t[7]*(t[10]*t[12]-t[8]*t[14]);

    cofactor[5] =
           t[0]*(t[10]*t[15]-t[11]*t[14])
         + t[2]*(t[11]*t[12]-t[8]*t[15])
         + t[3]*(t[8]*t[14]-t[10]*t[12]);

    cofactor[6] =
           t[0]*(t[7]*t[14]-t[6]*t[15])
         + t[2]*(t[4]*t[15]-t[7]*t[12])
         + t[3]*(t[6]*t[12]-t[4]*t[14]);

    cofactor[7] =
           t[0]*(t[6]*t[11]-t[7]*t[10])
         + t[2]*(t[7]*t[8]-t[4]*t[11])
         + t[3]*(t[4]*t[10]-t[6]*t[8]);

    cofactor[8] =
           t[4]*(t[9]*t[15]-t[11]*t[13])
         + t[5]*(t[11]*t[12]-t[8]*t[15])
         + t[7]*(t[8]*t[13]-t[9]*t[12]);

    cofactor[9] =
           t[0]*(t[11]*t[13]-t[9]*t[15])
         + t[1]*(t[8]*t[15]-t[11]*t[12])
         + t[3]*(t[9]*t[12]-t[8]*t[13]);

    cofactor[10]=
           t[0]*(t[5]*t[15]-t[7]*t[13])
         + t[1]*(t[7]*t[12]-t[4]*t[15])
         + t[3]*(t[4]*t[13]-t[5]*t[12]);

    cofactor[11]=
           t[0]*(t[7]*t[9]-t[5]*t[11])
         + t[1]*(t[4]*t[11]-t[7]*t[8])
         + t[3]*(t[5]*t[8]-t[4]*t[9]);

    cofactor[12]=
           t[4]*(t[10]*t[13]-t[9]*t[14])
         + t[5]*(t[8]*t[14]-t[10]*t[12])
         + t[6]*(t[9]*t[12]-t[8]*t[13]);

    cofactor[13]=
           t[0]*(t[9]*t[14]-t[10]*t[13])
         + t[1]*(t[10]*t[12]-t[8]*t[14])
         + t[2]*(t[8]*t[13]-t[9]*t[12]);

    cofactor[14]=
           t[0]*(t[6]*t[13]-t[5]*t[14])
         + t[1]*(t[4]*t[14]-t[6]*t[12])
         + t[2]*(t[5]*t[12]-t[4]*t[13]);

    cofactor[15]=
           t[0]*(t[5]*t[10]-t[6]*t[9])
         + t[1]*(t[6]*t[8]-t[4]*t[10])
         + t[2]*(t[4]*t[9]-t[5]*t[8]);

    return cofactor/det(t);
}

//- Inner-product between a tensor and a vector
inline vector
operator&(const tensor4& st, const vector& v)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatable types of vector and tensor4 "
        << "vector has 3 components while tensor4 has 4"
        << exit(FatalError);

    return result;
}

//- Inner-product between a spherical tensor and a vector
inline vector
operator&(const vector& v, const tensor4& st)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatable types of vector and tensor4 "
        << "vector has 3 components while tensor4 has 4"
        << exit(FatalError);

    return result;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
