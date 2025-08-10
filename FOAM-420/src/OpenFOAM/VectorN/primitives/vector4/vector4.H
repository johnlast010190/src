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
    Vector4

Description
    Vector4 obtained from generic VectorN

SourceFiles
    Vector4.C

\*---------------------------------------------------------------------------*/

#ifndef Vector4_H
#define Vector4_H

#include "primitives/Scalar/scalar/scalar.H"
#include "VectorN/primitives/VectorN.H"
#include "primitives/contiguous/contiguous.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef VectorN<scalar, 4> vector4;


//- Specify data associated with vector4 type is contiguous
template<>
constexpr bool contiguous<vector4>() {return true;}


//- Transform the spherical tensor
//- The components are assumed to be individual scalars
//- i.e. transform has no effect
template<>
inline vector4 transform
(
    const tensor& tt,
    const vector4& v
)
{
    vector4 vel4(v);

    vector vel(vector::zero);
    vel[0] = vel4[0];
    vel[1] = vel4[1];
    vel[2] = vel4[2];

    vel = tt & vel;

    vel4[0] = vel[0];
    vel4[1] = vel[1];
    vel4[2] = vel[2];

    return vel4;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// Include inline implementations
#include "VectorN/primitives/vector4/vector4I.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
