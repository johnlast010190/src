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
    (c) 2017 Esi Ltd.

Type
    sphericalTensor8

Description
    SphericalTensorN of 8 scalars.

\*---------------------------------------------------------------------------*/


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //


//- Inner-product between a spherical tensor and a vector
inline vector
operator&(const sphericalTensor8& st, const vector& v)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatable types of vector and sphericalTensor8 "
        << "vector has 3 components while sphericalTensor8 has 8"
        << exit(FatalError);
    return result;
}

//- Inner-product between a spherical tensor and a vector
inline vector
operator&(const vector& v, const sphericalTensor8& st)
{
    vector result(0);

    FatalErrorInFunction
        << "Incompatable types of vector and sphericalTensor8 "
        << "vector has 3 components while sphericalTensor8 has 8"
        << exit(FatalError);

    return result;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
