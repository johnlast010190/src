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
    (c) 2011 OpenFOAM Foundation

Typedef
    Foam::vectorFieldIOField

Description
    vectorFieldField with IO.

\*---------------------------------------------------------------------------*/

#ifndef VectorTensorNFieldIOField_H
#define VectorTensorNFieldIOField_H

#include "VectorN/Fields/VectorTensorNFieldsFwd.H"
#include "db/IOobjects/CompactIOField/CompactIOField.H"

#define TEMPLATE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define VectorTensorN_FieldIOFieldTypes(tensorType,diagTensorType,                                                        \
                               sphericalTensorType,vectorType,CmptType,                                                 \
                               args...)                                                                                   \
                                                                                                                         \
    typedef IOField<tensorType##Field > tensorType##FieldIOField;                                                         \
    typedef CompactIOField<tensorType##Field, tensorType > tensorType##FieldCompactIOField;                             \
    defineTemplateTypeNameAndDebugWithName                                                                                \
    (                                                                                                                    \
        tensorType##FieldIOField,                                                                                        \
        #tensorType"FieldField",                                                                                        \
        0                                                                                                                \
    );                                                                                                                    \
                                                                                                                        \
    defineTemplateTypeNameAndDebugWithName                                                                                \
    (                                                                                                                    \
        tensorType##FieldCompactIOField,                                                                                \
       #tensorType"FieldCompactField",                                                                                    \
        0                                                                                                                \
    );                                                                                                                    \
    typedef IOField<diagTensorType##Field > diagTensorType##FieldIOField;                                                 \
    typedef CompactIOField<diagTensorType##Field, diagTensorType > diagTensorType##FieldCompactIOField;                 \
    defineTemplateTypeNameAndDebugWithName                                                                                \
    (                                                                                                                    \
        diagTensorType##FieldIOField,                                                                                    \
        #diagTensorType"FieldField",                                                                                    \
        0                                                                                                                \
    );                                                                                                                    \
                                                                                                                        \
    defineTemplateTypeNameAndDebugWithName                                                                                \
    (                                                                                                                    \
        diagTensorType##FieldCompactIOField,                                                                            \
        #diagTensorType"FieldCompactField",                                                                                \
        0                                                                                                                \
    );                                                                                                                    \
    typedef IOField<sphericalTensorType##Field > sphericalTensorType##FieldIOField;                                         \
    typedef CompactIOField<sphericalTensorType##Field, sphericalTensorType > sphericalTensorType##FieldCompactIOField;     \
    defineTemplateTypeNameAndDebugWithName                                                                                \
    (                                                                                                                    \
        sphericalTensorType##FieldIOField,                                                                                \
        #sphericalTensorType"FieldField",                                                                                \
        0                                                                                                                \
    );                                                                                                                    \
                                                                                                                        \
    defineTemplateTypeNameAndDebugWithName                                                                                \
    (                                                                                                                    \
        sphericalTensorType##FieldCompactIOField,                                                                        \
        #sphericalTensorType"FieldCompactField",                                                                        \
        0                                                                                                                \
    );                                                                                                                    \
    typedef IOField<vectorType##Field > vectorType##FieldIOField;                                                         \
    typedef CompactIOField<vectorType##Field, vectorType > vectorType##FieldCompactIOField;                             \
    defineTemplateTypeNameAndDebugWithName                                                                                \
    (                                                                                                                    \
        vectorType##FieldIOField,                                                                                        \
        #vectorType"FieldField",                                                                                        \
        0                                                                                                                \
    );                                                                                                                    \
                                                                                                                        \
    defineTemplateTypeNameAndDebugWithName                                                                                \
    (                                                                                                                    \
        vectorType##FieldCompactIOField,                                                                                \
        #vectorType"FieldCompactField",                                                                                    \
        0                                                                                                                \
    );                                                                                                                    \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
