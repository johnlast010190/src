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
    VectorNFieldFields

Description
    Specialisations of FieldField<Field,T> for VectorN types.

SourceFiles
    VectorNFieldFields.C

\*---------------------------------------------------------------------------*/

#ifndef VectorNFieldFields_H
#define VectorNFieldFields_H

#include "VectorN/FieldFields/VectorTensorNFieldFieldsFwd.H"
#include "fields/FieldFields/FieldField/FieldField.H"

#define TEMPLATE template<template<class> class Field>
#include "fields/FieldFields/FieldField/FieldFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define VectorN_FieldFunctions(tensorType, diagTensorType, sphericalTensorType,         \
    vectorType, CmptType, args...)                                                      \
                                                                                        \
UNARY_FUNCTION(CmptType, vectorType, cmptSum)                                           \
                                                                                        \
BINARY_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)                       \
BINARY_TYPE_FUNCTION(vectorType, vectorType, vectorType, cmptMultiply)                  \
                                                                                        \
BINARY_OPERATOR(vectorType, vectorType, CmptType, /,divide)                             \
BINARY_TYPE_OPERATOR(vectorType, vectorType, CmptType, /,divide)                        \
                                                                                        \
BINARY_OPERATOR(vectorType, vectorType, vectorType, +, add)                             \
BINARY_OPERATOR(vectorType, vectorType, vectorType, -, subtract)                        \
                                                                                        \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, +, add)                        \
BINARY_TYPE_OPERATOR(vectorType, vectorType, vectorType, -, subtract)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

forAllVectorTensorNTypes(VectorN_FieldFunctions)

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef VectorN_FieldFunctions

#include "fields/Fields/Field/undefFieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "VectorN/FieldFields/VectorNFieldFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
