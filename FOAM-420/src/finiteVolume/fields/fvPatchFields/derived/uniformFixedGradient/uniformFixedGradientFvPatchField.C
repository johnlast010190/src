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
    (c) 2017 OpenCFD Ltd.
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/uniformFixedGradient/uniformFixedGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    uniformGradient_()
{}

template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchField<Type>(p, iF),
    uniformGradient_(Function1<Type>::New("uniformGradient", dict))
{
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    this->evaluate();
}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchField<Type>(ptf, p, iF, mapper),
    uniformGradient_(ptf.uniformGradient_, false)
{}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf
)
:
    fixedGradientFvPatchField<Type>(ptf),
    uniformGradient_(ptf.uniformGradient_, false)
{}


template<class Type>
Foam::uniformFixedGradientFvPatchField<Type>::uniformFixedGradientFvPatchField
(
    const uniformFixedGradientFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedGradientFvPatchField<Type>(ptf, iF),
    uniformGradient_(ptf.uniformGradient_, false)
{
    // Evaluate the profile if defined
    if (ptf.uniformGradient_.valid())
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedGradientFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    this->gradient() = uniformGradient_->value(t);

    fixedGradientFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedGradientFvPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFvPatchField<Type>::write(os);
    uniformGradient_->writeData(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
