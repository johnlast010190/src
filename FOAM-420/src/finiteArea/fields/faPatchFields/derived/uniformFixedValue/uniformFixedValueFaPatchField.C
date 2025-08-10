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
    (c) 2016-2017 Wikki Ltd
    (c) 2018 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/derived/uniformFixedValue/uniformFixedValueFaPatchField.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::
uniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(p, iF),
    uniformValue_()
{}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::
uniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchField<Type>(p, iF),
    uniformValue_(Function1<Type>::New("uniformValue", dict))
{
    if (dict.found("value"))
    {
        this->forceAssign(Field<Type>("value", dict, p.size()));
    }
    else
    {
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::
uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchField<Type>(ptf, p, iF, mapper),
    uniformValue_(ptf.uniformValue_, false)
{
    // Evaluate since value not mapped
    this->evaluate();
}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::
uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf
)
:
    fixedValueFaPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_, false)
{}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::
uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_, false)
{
    // Evaluate the profile if defined
    if (ptf.uniformValue_.valid())
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedValueFaPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    faPatchField<Type>::forceAssign(uniformValue_->value(t));

    fixedValueFaPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedValueFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    uniformValue_->writeData(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
