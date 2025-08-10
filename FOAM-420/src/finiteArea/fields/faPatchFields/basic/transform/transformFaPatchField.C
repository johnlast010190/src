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

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/basic/transform/transformFaPatchField.H"
#include "db/IOstreams/IOstreams.H"
#include "fields/Fields/transformField/transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
transformFaPatchField<Type>::transformFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
   faPatchField<Type>(p, iF)
{}


template<class Type>
transformFaPatchField<Type>::transformFaPatchField
(
    const transformFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
transformFaPatchField<Type>::transformFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, dict)
{}


template<class Type>
transformFaPatchField<Type>::transformFaPatchField
(
    const transformFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return the matrix diagonal coefficients corresponding to the
// evaluation of the value of this patchField
template<class Type>
tmp<Field<Type>> transformFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return pTraits<Type>::one - snGradTransformDiag();
}


// Return the matrix source coefficients corresponding to the
// evaluation of the value of this patchField
template<class Type>
tmp<Field<Type>> transformFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
        *this
      - cmptMultiply
        (
            valueInternalCoeffs(this->patch().weights()),
            this->patchInternalField()
        );
}


// Return the matrix diagonal coefficients corresponding to the
// evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type>> transformFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -this->patch().deltaCoeffs()*snGradTransformDiag();
}


// Return the matrix source coefficients corresponding to the
// evaluation of the gradient of this patchField
template<class Type>
tmp<Field<Type>> transformFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        snGrad()
      - cmptMultiply(gradientInternalCoeffs(), this->patchInternalField());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void transformFaPatchField<Type>::operator=
(
    const faPatchField<Type>& ptf
)
{
    this->evaluate();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
