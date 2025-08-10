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
    Foam::vectorIOField

Description
    vectorField with IO.

\*---------------------------------------------------------------------------*/

#ifndef VectorTensorNIOField_H
#define VectorTensorNIOField_H

#include "VectorN/Fields/VectorTensorNFieldsFwd.H"
#include "db/IOobjects/IOField/IOField.H"

#define TEMPLATE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define VectorTensorN_FieldIOTypes(tensorType,diagTensorType,                                                \
                               sphericalTensorType,vectorType,CmptType,                                     \
                               args...)                                                                       \
                                                                                                             \
    typedef IOField<tensorType> tensorType##IOField;                                                         \
    defineTemplateTypeNameAndDebugWithName(tensorType##IOField, #tensorType"Field", 0);                    \
                                                                                                            \
    typedef IOField<diagTensorType> diagTensorType##IOField;                                                 \
    defineTemplateTypeNameAndDebugWithName(diagTensorType##IOField, #diagTensorType"Field", 0);            \
                                                                                                            \
    typedef IOField<sphericalTensorType> sphericalTensorType##IOField;                                         \
    defineTemplateTypeNameAndDebugWithName(sphericalTensorType##IOField, #sphericalTensorType"Field", 0);    \
                                                                                                            \
    typedef IOField<vectorType> vectorType##IOField;                                                         \
    defineTemplateTypeNameAndDebugWithName(vectorType##IOField, #vectorType"Field", 0);                    \

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
