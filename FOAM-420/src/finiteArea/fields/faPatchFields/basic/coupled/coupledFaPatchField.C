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
    (c) 2016 Esi Ltd.
    (c) 2016-2017 Wikki Ltd

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/basic/coupled/coupledFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& f
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF, f)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    lduInterfaceField(refCast<const lduInterface>(p)),
    faPatchField<Type>(p, iF, dict)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    faPatchField<Type>(ptf)
{}


template<class Type>
coupledFaPatchField<Type>::coupledFaPatchField
(
    const coupledFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    lduInterfaceField(refCast<const lduInterface>(ptf.patch())),
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::coupledFaPatchField<Type>::snGrad() const
{
    return
        (patchNeighbourField() - this->patchInternalField())
       *this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type>> coupledFaPatchField<Type>::snGrad
(
    const scalarField& deltaCoeffs
) const
{
    return
        deltaCoeffs
       *(this->patchNeighbourField() - this->patchInternalField());
}


template<class Type>
void coupledFaPatchField<Type>::initEvaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }
}


template<class Type>
void coupledFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    Field<Type>::operator=
    (
        this->patch().weights()*this->patchInternalField()
      + (1.0 - this->patch().weights())*patchNeighbourField()
    );
}


template<class Type>
tmp<Field<Type>> coupledFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*w;
}


template<class Type>
tmp<Field<Type>> coupledFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& w
) const
{
    return Type(pTraits<Type>::one)*(1.0 - w);
}

template<class Type>
tmp<Field<Type>> coupledFaPatchField<Type>::gradientInternalCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return -Type(pTraits<Type>::one)*deltaCoeffs;
}

template<class Type>
tmp<Field<Type>> coupledFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::coupledFaPatchField<Type>::gradientBoundaryCoeffs
(
    const scalarField& deltaCoeffs
) const
{
    return -this->gradientInternalCoeffs(deltaCoeffs);
}


template<class Type>
tmp<Field<Type>> coupledFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return -gradientInternalCoeffs();
}


template<class Type>
void coupledFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
