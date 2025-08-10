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

Description
    Include the header files for all the primitive types that Fields are
    instantiated for.

\*---------------------------------------------------------------------------*/

#ifndef VectorNFieldTypes_H
#define VectorNFieldTypes_H

#include "VectorN/primitives/vector1/vector1.H"
#include "VectorN/primitives/vector2/vector2.H"
#include "VectorN/primitives/vector4/vector4.H"
#include "VectorN/primitives/vector6/vector6.H"
#include "VectorN/primitives/vector8/vector8.H"

#include "VectorN/primitives/vector1/tensor1.H"
#include "VectorN/primitives/vector2/tensor2.H"
#include "VectorN/primitives/vector4/tensor4.H"
#include "VectorN/primitives/vector6/tensor6.H"
#include "VectorN/primitives/vector8/tensor8.H"

#include "VectorN/primitives/vector1/diagTensor1.H"
#include "VectorN/primitives/vector2/diagTensor2.H"
#include "VectorN/primitives/vector4/diagTensor4.H"
#include "VectorN/primitives/vector6/diagTensor6.H"
#include "VectorN/primitives/vector8/diagTensor8.H"

#include "VectorN/primitives/vector1/sphericalTensor1.H"
#include "VectorN/primitives/vector2/sphericalTensor2.H"
#include "VectorN/primitives/vector4/sphericalTensor4.H"
#include "VectorN/primitives/vector6/sphericalTensor6.H"
#include "VectorN/primitives/vector8/sphericalTensor8.H"
#include "include/macros.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define CAPITALIZE_vector1 Vector1
#define CAPITALIZE_vector2 Vector2
#define CAPITALIZE_vector4 Vector4
#define CAPITALIZE_vector6 Vector6
#define CAPITALIZE_vector8 Vector8

#define CAPITALIZE_tensor1 Tensor1
#define CAPITALIZE_tensor2 Tensor2
#define CAPITALIZE_tensor4 Tensor4
#define CAPITALIZE_tensor6 Tensor6
#define CAPITALIZE_tensor8 Tensor8

#define CAPITALIZE_diagTensor1 DiagTensor1
#define CAPITALIZE_diagTensor2 DiagTensor2
#define CAPITALIZE_diagTensor4 DiagTensor4
#define CAPITALIZE_diagTensor6 DiagTensor6
#define CAPITALIZE_diagTensor8 DiagTensor8

#define CAPITALIZE_sphericalTensor1 SphericalTensor1
#define CAPITALIZE_sphericalTensor2 SphericalTensor2
#define CAPITALIZE_sphericalTensor4 SphericalTensor4
#define CAPITALIZE_sphericalTensor6 SphericalTensor6
#define CAPITALIZE_sphericalTensor8 SphericalTensor8

#define FOR_ALL_FIELD_TYPES(Macro, ...)                                        \
    Macro(scalar, __VA_ARGS__)                                                 \
    Macro(vector, __VA_ARGS__)                                                 \
    Macro(sphericalTensor, __VA_ARGS__)                                        \
    Macro(symmTensor, __VA_ARGS__)                                             \
    Macro(tensor, __VA_ARGS__)

#define forAllVectorNTypes(m, args...)              \
    m(vector1, Vector1, args)                       \
    m(vector2, Vector2, args)                       \
    m(vector4, Vector4, args)                       \
    m(vector6, Vector6, args)                       \
    m(vector8, Vector8, args)

#define forAllTensorNTypes(m, args...)              \
    m(tensor1, Tensor1, args)                       \
    m(tensor2, Tensor2, args)                       \
    m(tensor4, Tensor4, args)                       \
    m(tensor6, Tensor6, args)                       \
    m(tensor8, Tensor8, args)

#define forAllDiagTensorNTypes(m, args...)          \
    m(diagTensor1, DiagTensor1, args)               \
    m(diagTensor2, DiagTensor2, args)               \
    m(diagTensor4, DiagTensor4, args)               \
    m(diagTensor6, DiagTensor6, args)               \
    m(diagTensor8, DiagTensor8, args)

#define forAllSphericalTensorNTypes(m, args...)     \
    m(sphericalTensor1, SphericalTensor1, args)     \
    m(sphericalTensor2, SphericalTensor2, args)     \
    m(sphericalTensor4, SphericalTensor4, args)     \
    m(sphericalTensor6, SphericalTensor6, args)     \
    m(sphericalTensor8, SphericalTensor8, args)

#define forAllVectorTensorNTypes(m, args...)                            \
    m(tensor1, diagTensor1, sphericalTensor1, vector1, scalar, args)    \
    m(tensor2, diagTensor2, sphericalTensor2, vector2, scalar, args)    \
    m(tensor4, diagTensor4, sphericalTensor4, vector4, scalar, args)    \
    m(tensor6, diagTensor6, sphericalTensor6, vector6, scalar, args)    \
    m(tensor8, diagTensor8, sphericalTensor8, vector8, scalar, args)

#define INSTANTIATE_CLASS(Type, ClassName)                                     \
    template class ClassName<Type>;

#define INSTANTIATE_CLASS2(type, Type, ClassName)                              \
    template class ClassName<type>;

#define INSTANTIATE_FOR_ALL_FIELD_TYPES(ClassName)                             \
    FOR_ALL_FIELD_TYPES(INSTANTIATE_CLASS, ClassName)

#define INSTANTIATE_FOR_ALL_VECTOR_TENSOR_N_TYPES(ClassName)                   \
    forAllVectorNTypes(INSTANTIATE_CLASS2, ClassName)                           \
    forAllTensorNTypes(INSTANTIATE_CLASS2, ClassName)                           \
    forAllDiagTensorNTypes(INSTANTIATE_CLASS2, ClassName)                       \
    forAllSphericalTensorNTypes(INSTANTIATE_CLASS2, ClassName)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
