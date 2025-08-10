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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/Fields/transformField/transformField.H"
#include "fields/Fields/Field/FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

template<class Type>
List<Type> transform(const tensor& t, const UList<Type>& list)
{
    List<Type> newList(list.size());

    forAll(list, i)
    {
        newList[i] = transform(t, list[i]);
    }

    return newList;
}


template<class Type>
void transform
(
    Field<Type>& rtf,
    const tensorField& trf,
    const Field<Type>& tf
)
{
    if (trf.size() == 1)
    {
        return transform(rtf, trf[0], tf);
    }
    else
    {
        TFOR_ALL_F_OP_FUNC_F_F
        (
            Type, rtf, =, transform, tensor, trf, Type, tf
        )
    }
}


template<class Type>
tmp<Field<Type>> transform
(
    const tensorField& trf,
    const Field<Type>& tf
)
{
    tmp<Field<Type>> tranf(new Field<Type> (tf.size()));
    transform(tranf.ref(), trf, tf);
    return tranf;
}


template<class Type>
tmp<Field<Type>> transform
(
    const tensorField& trf,
    const tmp<Field<Type>>& ttf
)
{
    tmp<Field<Type>> tranf = New(ttf);
    transform(tranf.ref(), trf, ttf());
    ttf.clear();
    return tranf;
}


template<class Type>
tmp<Field<Type>> transform
(
    const tmp<tensorField>& ttrf,
    const Field<Type>& tf
)
{
    tmp<Field<Type>> tranf(new Field<Type> (tf.size()));
    transform(tranf.ref(), ttrf(), tf);
    ttrf.clear();
    return tranf;
}


template<class Type>
tmp<Field<Type>> transform
(
    const tmp<tensorField>& ttrf,
    const tmp<Field<Type>>& ttf
)
{
    tmp<Field<Type>> tranf = New(ttf);
    transform(tranf.ref(), ttrf(), ttf());
    ttf.clear();
    ttrf.clear();
    return tranf;
}


template<class Type>
void transform
(
    Field<Type>& rtf,
    const tensor& t,
    const Field<Type>& tf
)
{
    TFOR_ALL_F_OP_FUNC_S_F(Type, rtf, =, transform, tensor, t, Type, tf)
}


template<class Type>
tmp<Field<Type>> transform
(
    const tensor& t,
    const Field<Type>& tf
)
{
    tmp<Field<Type>> tranf(new Field<Type>(tf.size()));
    transform(tranf.ref(), t, tf);
    return tranf;
}


template<class Type>
tmp<Field<Type>> transform
(
    const tensor& t,
    const tmp<Field<Type>>& ttf
)
{
    tmp<Field<Type>> tranf = New(ttf);
    transform(tranf.ref(), t, ttf());
    ttf.clear();
    return tranf;
}

template<class Type>
void invTransform
(
    Field<Type>& result,
    const tensor& rot,
    const Field<Type>& fld
)
{
    TFOR_ALL_F_OP_FUNC_S_F
    (
        Type, result, =, invTransform, tensor, rot, Type, fld
    )
}


template<class Type>
void invTransform
(
    Field<Type>& result,
    const tensorField& rot,
    const Field<Type>& fld
)
{
    if (rot.size() == 1)
    {
        return invTransform(result, rot.first(), fld);
    }

    TFOR_ALL_F_OP_FUNC_F_F
    (
        Type, result, =, invTransform, tensor, rot, Type, fld
    )
}


template<class Type>
tmp<Foam::Field<Type>>
invTransform
(
    const tensorField& rot,
    const Field<Type>& fld
)
{
    auto tresult = tmp<Field<Type>>::New(fld.size());
    invTransform(tresult.ref(), rot, fld);
    return tresult;
}


template<class Type>
tmp<Foam::Field<Type>>
invTransform
(
    const tensorField& rot,
    const tmp<Field<Type>>& tfld
)
{
    tmp<Field<Type>> tresult = New(tfld);
    invTransform(tresult.ref(), rot, tfld());
    tfld.clear();
    return tresult;
}


template<class Type>
tmp<Foam::Field<Type>>
invTransform
(
    const tmp<tensorField>& trot,
    const Field<Type>& fld
)
{
    auto tresult = tmp<Field<Type>>::New(fld.size());
    invTransform(tresult.ref(), trot(), fld);
    trot.clear();
    return tresult;
}


template<class Type>
tmp<Foam::Field<Type>>
invTransform
(
    const tmp<tensorField>& trot,
    const tmp<Field<Type>>& tfld
)
{
    tmp<Field<Type>> tresult = New(tfld);
    invTransform(tresult.ref(), trot(), tfld());
    trot.clear();
    tfld.clear();
    return tresult;
}


template<class Type>
tmp<Foam::Field<Type>>
invTransform
(
    const tensor& rot,
    const Field<Type>& fld
)
{
    auto tresult = tmp<Field<Type>>::New(fld.size());
    invTransform(tresult.ref(), rot, fld);
    return tresult;
}


template<class Type>
tmp<Foam::Field<Type>>
invTransform
(
    const tensor& rot,
    const tmp<Field<Type>>& tfld
)
{
    tmp<Field<Type>> tresult = New(tfld);
    invTransform(tresult.ref(), rot, tfld());
    tfld.clear();
    return tresult;
}


template<class Type1, class Type2>
tmp<Field<Type1>> transformFieldMask(const Field<Type2>& f)
{
    return f;
}

template<class Type1, class Type2>
tmp<Field<Type1>> transformFieldMask(const tmp<Field<Type2>>& tf)
{
    return tmp<Field<Type1>>(tf.ptr());
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
