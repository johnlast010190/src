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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2018 OpenCFD Ltd.

InNamespace
    Foam

Description
    transformFieldField
    Spatial transformation functions for FieldField.

SourceFiles
    transformFieldField.C

\*---------------------------------------------------------------------------*/

#ifndef transformFieldField_H
#define transformFieldField_H

#include "primitives/transform/transform.H"
#include "fields/FieldFields/tensorFieldField/tensorFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
void transform
(
    FieldField<Field, Type>&,
    const FieldField<Field, tensor>&,
    const FieldField<Field, Type>&
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const FieldField<Field, tensor>&,
    const FieldField<Field, Type>&
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const FieldField<Field, tensor>&,
    const tmp<FieldField<Field, Type>>&
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tmp<FieldField<Field, tensor>>&,
    const FieldField<Field, Type>&
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tmp<FieldField<Field, tensor>>&,
    const tmp<FieldField<Field, Type>>&
);


template<template<class> class Field, class Type>
void transform
(
    FieldField<Field, Type>&,
    const tensor&,
    const FieldField<Field, Type>&
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tensor&,
    const FieldField<Field, Type>&
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tensor&,
    const tmp<FieldField<Field, Type>>&
);

// invTransform()

template<template<class> class Field, class Type>
void invTransform
(
    FieldField<Field, Type>& result,
    const tensor& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
void invTransform
(
    FieldField<Field, Type>& result,
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const FieldField<Field, tensor>& rot,
    const tmp<FieldField<Field, Type>>& tfld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const tmp<FieldField<Field, Type>>& tfld
);


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const tensor& rot,
    const FieldField<Field, Type>& fld
);

template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> invTransform
(
    const tensor& rot,
    const tmp<FieldField<Field, Type>>& tfld
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fields/FieldFields/transformFieldField/transformFieldField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
