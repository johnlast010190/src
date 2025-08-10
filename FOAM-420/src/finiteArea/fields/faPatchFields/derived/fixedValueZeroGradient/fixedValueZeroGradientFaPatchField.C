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

#include "fields/faPatchFields/derived/fixedValueZeroGradient/fixedValueZeroGradientFaPatchField.H"
#include "faMatrices/faMatrix/faMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedValueZeroGradientFaPatchField<Type>::
fixedValueZeroGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(p, iF)
{}


template<class Type>
Foam::fixedValueZeroGradientFaPatchField<Type>::
fixedValueZeroGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchField<Type>(p, iF)
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
Foam::fixedValueZeroGradientFaPatchField<Type>::
fixedValueZeroGradientFaPatchField
(
    const fixedValueZeroGradientFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchField<Type>(ptf, p, iF, mapper)
{
    // Evaluate since value not mapped
    this->evaluate();
}


template<class Type>
Foam::fixedValueZeroGradientFaPatchField<Type>::
fixedValueZeroGradientFaPatchField
(
    const fixedValueZeroGradientFaPatchField<Type>& ptf
)
:
    fixedValueFaPatchField<Type>(ptf)
{}


template<class Type>
Foam::fixedValueZeroGradientFaPatchField<Type>::
fixedValueZeroGradientFaPatchField
(
    const fixedValueZeroGradientFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedValueZeroGradientFaPatchField<Type>::manipulateMatrix
(
    faMatrix<Type>& matrix
)
{
    if (this->manipulatedMatrix())
    {
        return;
    }

    Field<Type>& fld = *this;

    matrix.setValues(faPatchField<Type>::patch().edgeFaces(), fld);

    faPatchField<Type>::manipulateMatrix(matrix);
}

template<class Type>
void Foam::fixedValueZeroGradientFaPatchField<Type>::manipulateMatrix
(
    faMatrix<Type>& matrix,
    const Field<Type>& weights
)
{
    if (this->manipulatedMatrix())
    {
        return;
    }

    Field<Type>& fld = *this;
    //fld *= weights;
    matrix.setValues(faPatchField<Type>::patch().edgeFaces(), fld);

    faPatchField<Type>::manipulateMatrix(matrix);
}

template<class Type>
void Foam::fixedValueZeroGradientFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
