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

Description

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/derived/fixedValueOutflow/fixedValueOutflowFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF)
{}


template<class Type>
fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, Field<Type>("value", dict, p.size()))
{}


template<class Type>
fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const fixedValueOutflowFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


template<class Type>
fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const fixedValueOutflowFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const fixedValueOutflowFaPatchField<Type>& ptf
)
:
    faPatchField<Type>(ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>> fixedValueOutflowFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& weights
) const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(Type(pTraits<Type>::one)*weights)
    );
}


template<class Type>
tmp<Field<Type>> fixedValueOutflowFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& weights
) const
{
    return (1.0 - weights)*(*this);
}


template<class Type>
tmp<Field<Type>>
fixedValueOutflowFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type>>
fixedValueOutflowFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


// Write
template<class Type>
void fixedValueOutflowFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
