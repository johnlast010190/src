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

#include "fields/faPatchFields/basic/basicSymmetry/basicSymmetryFaPatchField.H"
#include "fields/Fields/symmTransformField/symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(p, iF)
{}


template<class Type>
basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const basicSymmetryFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    transformFaPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    transformFaPatchField<Type>(p, iF, dict)
{
    this->evaluate();
}


template<class Type>
basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const basicSymmetryFaPatchField<Type>& ptf
)
:
    transformFaPatchField<Type>(ptf)
{}


template<class Type>
basicSymmetryFaPatchField<Type>::basicSymmetryFaPatchField
(
    const basicSymmetryFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    transformFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return gradient at boundary
template<class Type>
tmp<Field<Type>> basicSymmetryFaPatchField<Type>::snGrad() const
{
    tmp<vectorField> nHat = this->patch().edgeNormals();

    const Field<Type> iF(this->patchInternalField());

    return
    (
        transform(I - 2.0*sqr(nHat), iF) - iF)
        *(this->patch().deltaCoeffs()/2.0);
}


// Evaluate the field on the patch
template<class Type>
void basicSymmetryFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    vectorField nHat ( this->patch().edgeNormals() );

    const Field<Type> iF(this->patchInternalField());

    Field<Type>::operator=
    (
        (iF + transform(I - 2.0*sqr(nHat), iF))/2.0
    );

    transformFaPatchField<Type>::evaluate();
}


// Return defining fields
template<class Type>
tmp<Field<Type>> basicSymmetryFaPatchField<Type>::snGradTransformDiag() const
{
    vectorField nHat ( this->patch().edgeNormals() );

    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
