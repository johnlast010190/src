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

Description
    Spatial transformation functions for FieldField.

\*---------------------------------------------------------------------------*/

#include "fields/FieldFields/transformFieldField/transformFieldField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

template<template<class> class Field, class Type>
void transform
(
    FieldField<Field, Type>& rtf,
    const FieldField<Field, tensor>& trf,
    const FieldField<Field, Type>& tf
)
{
    forAll(rtf, i)
    {
        transform(rtf[i], trf[i], tf[i]);
    }
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const FieldField<Field, tensor>& trf,
    const FieldField<Field, Type>& tf
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(tf)
    );
    transform(tranf(), trf, tf);
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const FieldField<Field, tensor>& trf,
    const tmp<FieldField<Field, Type>>& ttf
)
{
    tmp<FieldField<Field, Type>> tranf(ttf.ptr());
    transform(tranf(), trf, tranf());
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tmp<FieldField<Field, tensor>>& ttrf,
    const FieldField<Field, Type>& tf
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(tf)
    );
    transform(tranf(), ttrf(), tf);
    ttrf.clear();
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tmp<FieldField<Field, tensor>>& ttrf,
    const tmp<FieldField<Field, Type>>& ttf
)
{
    tmp<FieldField<Field, Type>> tranf(ttf.ptr());
    transform(tranf(), ttrf(), tranf());
    ttrf.clear();
    return tranf;
}


template<template<class> class Field, class Type>
void transform
(
    FieldField<Field, Type>& rtf,
    const tensor& t,
    const FieldField<Field, Type>& tf
)
{
    forAll(rtf, i)
    {
        transform(rtf[i], t, tf[i]);
    }
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tensor& t,
    const FieldField<Field, Type>& tf
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(tf)
    );
    transform(tranf(), t, tf);
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>> transform
(
    const tensor& t,
    const tmp<FieldField<Field, Type>>& ttf
)
{
    tmp<FieldField<Field, Type>> tranf(ttf.ptr());
    transform(tranf(), t, tranf());
    return tranf;
}

template<template<class> class Field, class Type>
void invTransform
(
    FieldField<Field, Type>& result,
    const tensor& rot,
    const FieldField<Field, Type>& fld
)
{
    forAll(result, i)
    {
        invTransform(result[i], rot, fld[i]);
    }
}


template<template<class> class Field, class Type>
void invTransform
(
    FieldField<Field, Type>& result,
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
)
{
    forAll(result, i)
    {
        invTransform(result[i], rot[i], fld[i]);
    }
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>>
invTransform
(
    const FieldField<Field, tensor>& rot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tranf
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    invTransform(tranf(), rot, fld);
    return tranf;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>>
invTransform
(
    const FieldField<Field, tensor>& rot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    invTransform(tresult(), rot, tresult());
    return tresult;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>>
invTransform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tresult
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    invTransform(tresult(), trot(), fld);
    trot.clear();
    return tresult;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>>
invTransform
(
    const tmp<FieldField<Field, tensor>>& trot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    invTransform(tresult(), trot(), tresult());
    trot.clear();
    return tresult;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>>
invTransform
(
    const tensor& rot,
    const FieldField<Field, Type>& fld
)
{
    tmp<FieldField<Field, Type>> tresult
    (
        FieldField<Field, Type>::NewCalculatedType(fld)
    );
    invTransform(tresult(), rot, fld);
    return tresult;
}


template<template<class> class Field, class Type>
tmp<FieldField<Field, Type>>
invTransform
(
    const tensor& rot,
    const tmp<FieldField<Field, Type>>& tfld
)
{
    tmp<FieldField<Field, Type>> tresult(tfld.ptr());
    invTransform(tresult(), rot, tresult());
    return tresult;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
