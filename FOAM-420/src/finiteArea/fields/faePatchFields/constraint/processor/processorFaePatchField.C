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

#include "fields/faePatchFields/constraint/processor/processorFaePatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
processorFaePatchField<Type>::processorFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF
)
:
    coupledFaePatchField<Type>(p, iF),
    procPatch_(refCast<const processorFaPatch>(p))
{}


template<class Type>
processorFaePatchField<Type>::processorFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const Field<Type>& f
)
:
    coupledFaePatchField<Type>(p, iF, f),
    procPatch_(refCast<const processorFaPatch>(p))
{}


// Construct by mapping given processorFaePatchField<Type>
template<class Type>
processorFaePatchField<Type>::processorFaePatchField
(
    const processorFaePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaePatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const processorFaPatch>(p))
{
    if (!isType<processorFaPatch>(this->patch()))
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
processorFaePatchField<Type>::processorFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, faEdgeMesh>& iF,
    const dictionary& dict
)
:
    coupledFaePatchField<Type>(p, iF, dict),
    procPatch_(refCast<const processorFaPatch>(p))
{
    if (!isType<processorFaPatch>(p))
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
processorFaePatchField<Type>::processorFaePatchField
(
    const processorFaePatchField<Type>& ptf
)
:
    coupledFaePatchField<Type>(ptf),
    procPatch_(refCast<const processorFaPatch>(ptf.patch()))
{}


template<class Type>
processorFaePatchField<Type>::processorFaePatchField
(
    const processorFaePatchField<Type>& ptf,
    const DimensionedField<Type, faEdgeMesh>& iF
)
:
    coupledFaePatchField<Type>(ptf, iF),
    procPatch_(refCast<const processorFaPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
processorFaePatchField<Type>::~processorFaePatchField()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
