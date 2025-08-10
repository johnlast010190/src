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
    (c) 2018-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/faPatchFields/constraint/cyclic/cyclicFaPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{}


template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{
    if (!isA<cyclicFaPatch>(this->patch()))
    {
         FatalErrorInFunction
             << "' not constraint type '" << typeName << "'"
             << "\n    for patch " << p.name()
             << " of field " << this->internalField().name()
             << " in file " << this->internalField().objectPath()
             << exit(FatalIOError);

    }
}


template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    coupledFaPatchField<Type>(p, iF, dict),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{
     if (!isA<cyclicFaPatch>(p))
     {
         FatalIOErrorInFunction
         (
             dict
         )   << "    patch type '" << p.type()
             << "' not constraint type '" << typeName << "'"
             << "\n    for patch " << p.name()
             << " of field " << this->internalField().name()
             << " in file " << this->internalField().objectPath()
             << exit(FatalIOError);
     }

    this->evaluate(Pstream::commsTypes::blocking);
}


template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    coupledFaPatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
cyclicFaPatchField<Type>::cyclicFaPatchField
(
    const cyclicFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    coupledFaPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicFaPatchField<Type>::patchNeighbourField() const
{
    const Field<Type>& iField = this->primitiveField();
    const labelUList& faceCells = cyclicPatch_.faceCells();

    tmp<Field<Type>> tpnf(new Field<Type>(this->size()));
    Field<Type>& pnf = tpnf.ref();

    label sizeby2 = this->size()/2;

    if (doTransform())
    {
        for (label facei=0; facei<sizeby2; ++facei)
        {
            pnf[facei] = Foam::transform
            (
                forwardT()[0], iField[faceCells[facei + sizeby2]]
            );

            pnf[facei + sizeby2] = Foam::transform
            (
                reverseT()[0], iField[faceCells[facei]]
            );
        }
    }
    else
    {
        for (label facei=0; facei<sizeby2; ++facei)
        {
            pnf[facei] = iField[faceCells[facei + sizeby2]];
            pnf[facei + sizeby2] = iField[faceCells[facei]];
        }
    }

    return tpnf;
}


template<class Type>
void cyclicFaPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes commsType
) const
{
    scalarField pnf(this->size());

    label sizeby2 = this->size()/2;
    const labelUList& faceCells = cyclicPatch_.faceCells();

    for (label facei = 0; facei < sizeby2; ++facei)
    {
        pnf[facei] = psiInternal[faceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[faceCells[facei]];
    }

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    if (add)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}

template<class Type>
void cyclicFaPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes commsType
) const
{
    Field<Type> pnf(this->size());

    label sizeby2 = this->size()/2;
    const labelUList& faceCells = cyclicPatch_.faceCells();

    for (label facei = 0; facei < sizeby2; ++facei)
    {
        pnf[facei] = psiInternal[faceCells[facei + sizeby2]];
        pnf[facei + sizeby2] = psiInternal[faceCells[facei]];
    }

    // Multiply the field by coefficients and add into the result
    if (add)
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] += coeffs[elemI]*pnf[elemI];
        }
    }
    else
    {
        forAll(faceCells, elemI)
        {
            result[faceCells[elemI]] -= coeffs[elemI]*pnf[elemI];
        }
    }
}


template<class Type>
void cyclicFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
