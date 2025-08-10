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
    (c) 2004-2011 H. Jasak

Description
    Global functions for expansion and contraction of tensor field
    to diagonal type

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    ExpandTensorNField.C

\*---------------------------------------------------------------------------*/

#ifndef expandTensorNField_H
#define expandTensorNField_H

#include "fields/Fields/Field/Field.H"
#include "VectorN/Fields/VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UNARY_FUNCTION(typeF1, typeF2, FUNC)                                  \
void FUNC(Field<typeF1>& f1, const UList<typeF2>& f2);                        \
void FUNC(Field<typeF1>& f1, const tmp<Field<typeF2>>& tf2);


#define ExpandFieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
        vectorType, cmptType, args...)                                        \
                                                                              \
UNARY_FUNCTION(cmptType, tensorType, contractScalar)                          \
UNARY_FUNCTION(cmptType, diagTensorType, contractScalar)                      \
UNARY_FUNCTION(cmptType, sphericalTensorType, contractScalar)                 \
UNARY_FUNCTION(cmptType, vectorType, contractScalar)                          \
                                                                              \
UNARY_FUNCTION(vectorType, tensorType, contractLinear)                        \
UNARY_FUNCTION(vectorType, diagTensorType, contractLinear)                    \
UNARY_FUNCTION(vectorType, sphericalTensorType, contractLinear)               \
                                                                              \
UNARY_FUNCTION(vectorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(tensorType, cmptType, expandScalar)                            \
UNARY_FUNCTION(diagTensorType, cmptType, expandScalar)                        \
UNARY_FUNCTION(sphericalTensorType, cmptType, expandScalar)                   \
                                                                              \
UNARY_FUNCTION(tensorType, vectorType, expandLinear)                          \
UNARY_FUNCTION(diagTensorType, vectorType, expandLinear)                      \
UNARY_FUNCTION(sphericalTensorType, vectorType, expandLinear)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

forAllVectorTensorNTypes(ExpandFieldFunctions)


// template <class Cmpt, int length>
// void contractScalar
// (
//     Field<Cmpt>& res,
//     const UList<TensorN<Cmpt, length>>& f
// );
//
// template <class Cmpt, int length>
// void contractScalar
// (
//     Field<Cmpt>& res,
//     const UList<DiagTensorN<Cmpt, length>>& f
// );
//
// template <class Cmpt, int length>
// void contractScalar
// (
//     Field<Cmpt>& res,
//     const UList<SphericalTensorN<Cmpt, length>>& f
// );
//
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// template <class Cmpt, int length>
// void contractLinear
// (
//     Field<VectorN<Cmpt, length>>& res,
//     const UList<TensorN<Cmpt, length>>& f
// );
//
// template <class Cmpt, int length>
// void contractLinear
// (
//     Field<VectorN<Cmpt, length>>& res,
//     const UList<DiagTensorN<Cmpt, length>>& f
// );
//
// template <class Cmpt, int length>
// void contractLinear
// (
//     Field<VectorN<Cmpt, length>>& res,
//     const UList<SphericalTensorN<Cmpt, length>>& f
// );
//
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// template <class Cmpt, int length>
// void expandScalar
// (
//     Field<VectorN<Cmpt, length>>& res,
//     const UList<Cmpt>& f
// );
//
// template <class Cmpt, int length>
// void expandScalar
// (
//     Field<TensorN<Cmpt, length>>& res,
//     const UList<Cmpt>& f
// );
//
// template <class Cmpt, int length>
// void expandScalar
// (
//     Field<DiagTensorN<Cmpt, length>>& res,
//     const UList<Cmpt>& f
// );
//
// template <class Cmpt, int length>
// void expandScalar
// (
//     Field<SphericalTensorN<Cmpt, length>>& res,
//     const UList<Cmpt>& f
// );
//
// // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// template <class Cmpt, int length>
// void expandLinear
// (
//     Field<TensorN<Cmpt, length>>& res,
//     const UList<VectorN<Cmpt, length>>& f
// );
//
// template <class Cmpt, int length>
// void expandLinear
// (
//     Field<DiagTensorN<Cmpt, length>>& res,
//     const UList<VectorN<Cmpt, length>>& f
// );
//
// template <class Cmpt, int length>
// void expandLinear
// (
//     Field<SphericalTensorN<Cmpt, length>>& res,
//     const UList<VectorN<Cmpt, length>>& f
// );

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef UNARY_FUNCTION
#undef ExpandFieldFunctions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// #ifdef NoRepository
// #include "VectorN/expandContract/ExpandTensorNField.C"
// #endif

// ************************************************************************* //

#endif
