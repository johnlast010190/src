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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/Fields/symmTransformField/symmTransformField.H"
#include "fields/Fields/Field/FieldM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * global functions  * * * * * * * * * * * * * //

template<class Type>
void transform
(
    Field<Type>& rtf,
    const symmTensorField& trf,
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
            Type, rtf, =, transform, symmTensor, trf, Type, tf
        )
    }
}


template<class Type>
tmp<Field<Type>> transform
(
    const symmTensorField& trf,
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
    const symmTensorField& trf,
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
    const tmp<symmTensorField>& ttrf,
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
    const tmp<symmTensorField>& ttrf,
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
    const symmTensor& t,
    const Field<Type>& tf
)
{
    TFOR_ALL_F_OP_FUNC_S_F(Type, rtf, =, transform, tensor, t, Type, tf)
}


template<class Type>
tmp<Field<Type>> transform
(
    const symmTensor& t,
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
    const symmTensor& t,
    const tmp<Field<Type>>& ttf
)
{
    tmp<Field<Type>> tranf = New(ttf);
    transform(tranf.ref(), t, ttf());
    ttf.clear();
    return tranf;
}


template<class Type>
tmp<Field<Type>> symmetryValue(const Field<Type>& f, const vectorField& n)
{
    tmp<vectorField> nHat = n/stabilise(mag(n), VSMALL);
    return (f + transform(I - 2.0*sqr(nHat), f))/2.0;
}

template<class Type>
tmp<Field<Type>> symmetryValue(const tmp<Field<Type>>& f, const vectorField& n)
{
    return symmetryValue(f(), n);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
