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
    (c) 1991-2005 OpenCFD Ltd.
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/faePatchFields/basic/calculated/calculatedFaePatchField.H"
#include "fields/faPatchFields/faPatchField/faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
const word& faePatchField<Type>::calculatedType()
{
    return calculatedFaePatchField<Type>::typeName;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF
)
:
    faePatchField<Type>(p, iF)
{}


template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const calculatedFaePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faePatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const dictionary& dict
)
:
    faePatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const calculatedFaePatchField<Type>& ptf
)
:
    faePatchField<Type>(ptf)
{}


template<class Type>
calculatedFaePatchField<Type>::calculatedFaePatchField
(
    const calculatedFaePatchField<Type>& ptf,
    const DimensionedField<Type, faEdgeMesh>& iF
)
:
    faePatchField<Type>(ptf, iF)
{}


template<class Type>
tmp<faePatchField<Type>> faePatchField<Type>::NewCalculatedType
(
    const faPatch& p
)
{
    typename patchConstructorTable::iterator patchTypeCstrIter =
        patchConstructorTable_().find(p.type());

    if (patchTypeCstrIter != patchConstructorTable_().end())
    {
        return patchTypeCstrIter->second
        (
            p,
            DimensionedField<Type, faEdgeMesh>::null()
        );
    }
    else
    {
        return tmp<faePatchField<Type>>
        (
            new calculatedFaePatchField<Type>
            (
                p,
                DimensionedField<Type, faEdgeMesh>::null()
            )
        );
    }
}

template<class Type>
template<class Type2>
tmp<faePatchField<Type>> faePatchField<Type>::NewCalculatedType
(
    const faePatchField<Type2>& pf
)
{
    return NewCalculatedType(pf.patch());
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
