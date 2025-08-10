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

Class
    tensorNFields

Description
    Specialisation of Field<T> for TensorN.

SourceFiles
    TensorNFields.C

\*---------------------------------------------------------------------------*/

#ifndef TensorNFields_H
#define TensorNFields_H

#include "VectorN/Fields/VectorTensorNFieldsFwd.H"
#include "fields/Fields/Field/Field.H"

#define TEMPLATE
#include "fields/Fields/Field/FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define TensorN_FieldFunctions(tensorType,diagTensorType,                    \
                               sphericalTensorType,vectorType,CmptType,      \
                               args...)                                      \
                                                                             \
UNARY_FUNCTION(tensorType, tensorType, inv)                                  \
UNARY_FUNCTION(diagTensorType, tensorType, diag)                             \
UNARY_FUNCTION(tensorType, tensorType, negSumDiag)                           \
                                                                             \
BINARY_OPERATOR(tensorType, CmptType, tensorType, /, divide)                 \
BINARY_TYPE_OPERATOR(tensorType, CmptType, tensorType, /, divide)            \
                                                                             \
BINARY_OPERATOR(vectorType, vectorType, tensorType, /, divide)               \
BINARY_TYPE_OPERATOR(vectorType, vectorType, tensorType, /, divide)          \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, tensorType, /, divide)               \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, /, divide)          \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, /, divide)           \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, /, divide)      \
                                                                             \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, /, divide)           \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, /, divide)      \
                                                                             \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, /, divide)      \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, /, divide) \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, /, divide)      \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, /, divide) \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, tensorType, +, add)                  \
BINARY_OPERATOR(tensorType, tensorType, tensorType, -, subtract)             \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, +, add)             \
BINARY_TYPE_OPERATOR(tensorType, tensorType, tensorType, -, subtract)        \
                                                                             \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, +, add)              \
BINARY_OPERATOR(tensorType, diagTensorType, tensorType, -, subtract)         \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, +, add)         \
BINARY_TYPE_OPERATOR(tensorType, diagTensorType, tensorType, -, subtract)    \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, +, add)              \
BINARY_OPERATOR(tensorType, tensorType, diagTensorType, -, subtract)         \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, +, add)         \
BINARY_TYPE_OPERATOR(tensorType, tensorType, diagTensorType, -, subtract)    \
                                                                             \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, +, add)         \
BINARY_OPERATOR(tensorType, sphericalTensorType, tensorType, -, subtract)    \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, +, add)    \
BINARY_TYPE_OPERATOR(tensorType, sphericalTensorType, tensorType, -, subtract) \
                                                                             \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, +, add)         \
BINARY_OPERATOR(tensorType, tensorType, sphericalTensorType, -, subtract)    \
                                                                             \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, +, add)    \
BINARY_TYPE_OPERATOR(tensorType, tensorType, sphericalTensorType, -, subtract)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forAllVectorTensorNTypes(TensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef TensorN_FieldFunctions

#include "fields/Fields/Field/undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
