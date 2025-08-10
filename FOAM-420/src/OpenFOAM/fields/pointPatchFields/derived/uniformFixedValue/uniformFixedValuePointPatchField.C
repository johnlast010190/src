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

\*---------------------------------------------------------------------------*/

#include "fields/pointPatchFields/derived/uniformFixedValue/uniformFixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(p, iF),
    uniformValue_()
{}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<Type>(p, iF, dict, false),
    uniformValue_(Function1<Type>::New("uniformValue", dict))
{
    if (dict.found("value"))
    {
        fixedValuePointPatchField<Type>::forceAssign
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        const scalar t = this->db().time().timeOutputValue();
        fixedValuePointPatchField<Type>::operator=(uniformValue_->value(t));
    }
}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<Type>(ptf, p, iF, mapper),
    uniformValue_(ptf.uniformValue_, false)
{
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();
    fixedValuePointPatchField<Type>::operator=(uniformValue_->value(t));
}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf),
    uniformValue_(ptf.uniformValue_, false)
{}


template<class Type>
Foam::uniformFixedValuePointPatchField<Type>::
uniformFixedValuePointPatchField
(
    const uniformFixedValuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF),
    uniformValue_(ptf.uniformValue_, false)
{
    // For safety re-evaluate
    const scalar t = this->db().time().timeOutputValue();
    fixedValuePointPatchField<Type>::forceAssign(uniformValue_->value(t));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    fixedValuePointPatchField<Type>::forceAssign(uniformValue_->value(t));

    fixedValuePointPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedValuePointPatchField<Type>::
write(Ostream& os) const
{
    // Note: write value
    fixedValuePointPatchField<Type>::write(os);
    uniformValue_->writeData(os);
}


// ************************************************************************* //
