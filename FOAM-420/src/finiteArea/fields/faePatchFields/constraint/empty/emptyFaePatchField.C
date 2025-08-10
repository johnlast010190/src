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

#include "fields/faePatchFields/constraint/empty/emptyFaePatchField.H"
#include "fields/faPatchFields/faPatchField/faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
emptyFaePatchField<Type>::emptyFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF
)
:
    faePatchField<Type>(p, iF, Field<Type>(0))
{}


template<class Type>
emptyFaePatchField<Type>::emptyFaePatchField
(
    const emptyFaePatchField<Type>&,
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const faPatchFieldMapper&
)
:
    faePatchField<Type>(p, iF, Field<Type>(0))
{
    if (!isType<emptyFaPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
emptyFaePatchField<Type>::emptyFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const dictionary& dict
)
:
    faePatchField<Type>(p, iF, Field<Type>(0))
{
    if (!isType<emptyFaPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "patch " << this->patch().index() << " not empty type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
emptyFaePatchField<Type>::emptyFaePatchField
(
    const emptyFaePatchField<Type>& ptf
)
:
    faePatchField<Type>
    (
        ptf.patch(),
        ptf.internalField(),
        Field<Type>(0)
    )
{}


template<class Type>
emptyFaePatchField<Type>::emptyFaePatchField
(
    const emptyFaePatchField<Type>& ptf,
    const DimensionedField<Type, faEdgeMesh>& iF
)
:
    faePatchField<Type>(ptf.patch(), iF, Field<Type>(0))
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
