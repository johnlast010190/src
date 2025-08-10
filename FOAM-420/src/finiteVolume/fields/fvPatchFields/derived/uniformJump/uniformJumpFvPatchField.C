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
    (c) 2012-2017 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/uniformJump/uniformJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpFvPatchField<Type>(p, iF),
    jumpTable_(new Function1<Type>("jumpTable"))
{}


template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const uniformJumpFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<Type>(ptf, p, iF, mapper),
    jumpTable_(ptf.jumpTable_, false)
{}


template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedJumpFvPatchField<Type>(p, iF),
    jumpTable_(new Function1<Type>("jumpTable"))
{
    if (this->cyclicPatch().owner())
    {
        jumpTable_ = Function1<Type>::New("jumpTable", dict);
    }

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate(Pstream::commsTypes::blocking);
    }
}


template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const uniformJumpFvPatchField<Type>& ptf
)
:
    fixedJumpFvPatchField<Type>(ptf),
    jumpTable_(ptf.jumpTable_, false)
{}


template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const uniformJumpFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpFvPatchField<Type>(ptf, iF),
    jumpTable_(ptf.jumpTable_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformJumpFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicPatch().owner())
    {
        this->jump_ = jumpTable_->value(this->db().time().value());
    }

    fixedJumpFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformJumpFvPatchField<Type>::write(Ostream& os) const
{
    fixedJumpFvPatchField<Type>::write(os);
    if (this->cyclicPatch().owner())
    {
        jumpTable_->writeData(os);
    }
}


// ************************************************************************* //
