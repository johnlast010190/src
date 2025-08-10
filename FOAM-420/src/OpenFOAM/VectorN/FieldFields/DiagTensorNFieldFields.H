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
    DiagTensorNFieldFields

Description
    Specialisation of FieldField<Field,T> for DiagTensorN

SourceFiles
    DiagTensorNFieldFields.C

\*---------------------------------------------------------------------------*/

#ifndef DiagTensorNFieldFields_H
#define DiagTensorNFieldFields_H

#include "VectorN/FieldFields/VectorTensorNFieldFieldsFwd.H"
#include "fields/FieldFields/FieldField/FieldField.H"

#define TEMPLATE template<template<class> class Field>
#include "fields/FieldFields/FieldField/FieldFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define DiagTensorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType, \
    vectorType, CmptType, args...)                                                  \
                                                                                    \
UNARY_FUNCTION(diagTensorType, diagTensorType, inv)                                 \
UNARY_FUNCTION(diagTensorType, diagTensorType, diag)                                \
UNARY_FUNCTION(vectorType, diagTensorType, contractLinear)                          \
UNARY_FUNCTION(CmptType, diagTensorType, contractScalar)                            \
                                                                                    \
BINARY_OPERATOR(diagTensorType, CmptType, diagTensorType, /, divide)                \
BINARY_TYPE_OPERATOR(diagTensorType, CmptType, diagTensorType, /, divide)           \
                                                                                    \
BINARY_OPERATOR(vectorType, vectorType, diagTensorType, /, divide)                  \
BINARY_TYPE_OPERATOR(vectorType, vectorType, diagTensorType, /, divide)             \
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /, divide)          \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, /, divide)     \
                                                                                    \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /, divide)     \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, /, divide)\
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /, divide)     \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, /, divide)\
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +, add)             \
BINARY_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -, subtract)        \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, +, add)        \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, diagTensorType, -, subtract)   \
                                                                                    \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +, add)        \
BINARY_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -, subtract)   \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, +, add)   \
BINARY_TYPE_OPERATOR(diagTensorType, sphericalTensorType, diagTensorType, -, subtract)  \
                                                                                    \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +, add)        \
BINARY_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -, subtract)   \
                                                                                    \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, +, add)   \
BINARY_TYPE_OPERATOR(diagTensorType, diagTensorType, sphericalTensorType, -, subtract)  \


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forAllVectorTensorNTypes(DiagTensorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef DiagTensorN_FieldFunctions

#include "fields/Fields/Field/undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "VectorN/FieldFields/DiagTensorNFieldFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
