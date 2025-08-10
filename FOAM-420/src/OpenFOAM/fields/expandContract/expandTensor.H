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
    (c) 2011 H. Jasak

Description
    Global functions for expansion and contraction of tensor coefficient
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef expandTensor_H
#define expandTensor_H

#include "primitives/Tensor/Tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Return the diagonal of a tensor as a scalar
template <class Cmpt>
inline void contractScalar(Cmpt& result, const Tensor<Cmpt>& t)
{
    result = 1.0/3.0*(t.xx() + t.yy() + t.zz());
}


//- Return the diagonal of a tensor as a scalar
template <class Cmpt>
inline Cmpt contractScalar(const Tensor<Cmpt>& t)
{
    Cmpt result;
    contractScalar(result, t);
    return result;
}


//- Return the diagonal of a tensor as a vector
template <class Cmpt>
inline void contractLinear(Vector<Cmpt>& result, const Tensor<Cmpt>& t)
{
    result = Vector<Cmpt>(t.xx(), t.yy(), t.zz());
}


//- Return the diagonal of a tensor as a vector
template <class Cmpt>
inline Vector<Cmpt> contractLinear(const Tensor<Cmpt>& t)
{
    Vector<Cmpt> result;
    contractLinear(result, t);
    return result;
}


//- Return the vector given a scalar
template <class Cmpt>
inline void expandScalar(Vector<Cmpt>& result, const Cmpt& v)
{
    result = Vector<Cmpt>(v, v, v);
}


//- Return the tensor given a scalar
template <class Cmpt>
inline void expandScalar(Tensor<Cmpt>& result, const Cmpt& v)
{
    result = Tensor<Cmpt>
    (
        v,                    pTraits<Cmpt>::zero,     pTraits<Cmpt>::zero,
        pTraits<Cmpt>::zero,  v,                       pTraits<Cmpt>::zero,
        pTraits<Cmpt>::zero,  pTraits<Cmpt>::zero,     v
    );
}


//- Return the tensor given a diagonal vector
template <class Cmpt>
inline void expandLinear(Tensor<Cmpt>& result, const Vector<Cmpt>& v)
{
    result = Tensor<Cmpt>
    (
        v.x(),                pTraits<Cmpt>::zero,     pTraits<Cmpt>::zero,
        pTraits<Cmpt>::zero,  v.y(),                   pTraits<Cmpt>::zero,
        pTraits<Cmpt>::zero,  pTraits<Cmpt>::zero,     v.z()
    );
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
