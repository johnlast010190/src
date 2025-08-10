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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/constraint/cyclic/cyclicFvPatchField.H"
#include "fields/Fields/transformField/transformField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFvPatch>(p)),
    calculatedMode_(true)
{}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    coupledFvPatchField<Type>(p, iF, dict, false),
    cyclicPatch_(refCast<const cyclicFvPatch>(p)),
    calculatedMode_(true)
{
    if (!isA<cyclicFvPatch>(p))
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
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvPatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFvPatch>(p)),
    calculatedMode_(true)
{
    if (!isA<cyclicFvPatch>(this->patch()))
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
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf
)
:
    cyclicLduInterfaceField(),
    coupledFvPatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_),
    calculatedMode_(true)
{}


template<class Type>
Foam::cyclicFvPatchField<Type>::cyclicFvPatchField
(
    const cyclicFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    coupledFvPatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_),
    calculatedMode_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicFvPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes
) const
{
    if (calculatedMode_)
    {
        // Because the underlying patch is coupled, the boundary values were
        // obtained by linear interpolation from internal and neighbour values
        // (in the default implementation of evaluate()). We use this to deduce
        // the neighbour field values. This is necessary because a derived type
        // may have changed the implemention of patchNeighbourField below,
        // but will be substituted by this constraint type in field calculations
        return
            (*this - this->patch().weights()*this->patchInternalField())
           /(1-this->patch().weights());
    }
    else
    {
        const Field<Type>& iField = this->primitiveField();
        const labelUList& nbrFaceCells =
            cyclicPatch().nbrPatch().faceCells();

        tmp<Field<Type>> tpnf(new Field<Type>(this->size()));
        Field<Type>& pnf = tpnf.ref();

        forAll(pnf, facei)
        {
            pnf[facei] = transform().transform(iField[nbrFaceCells[facei]]);
        }

        return tpnf;
    }
}


template<class Type>
const Foam::cyclicFvPatchField<Type>&
Foam::cyclicFvPatchField<Type>::nbrPatchField() const
{
    const GeometricField<Type, fvPatchField, volMesh>& fld =
    static_cast<const GeometricField<Type, fvPatchField, volMesh>&>
    (
        this->primitiveField()
    );

    return refCast<const cyclicFvPatchField<Type>>
    (
        fld.boundaryField()[this->cyclicPatch().nbrPatchID()]
    );
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    calculatedMode_ = false;
    coupledFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    scalarField& result,
    const bool add,
    const scalarField& psiInternal,
    const scalarField& coeffs,
    const direction cmpt,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicPatch().nbrPatch().faceCells();

    scalarField pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf, cmpt);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    Field<Type>& result,
    const bool add,
    const Field<Type>& psiInternal,
    const scalarField& coeffs,
    const Pstream::commsTypes
) const
{
    const labelUList& nbrFaceCells =
        cyclicPatch().nbrPatch().faceCells();

    Field<Type> pnf(psiInternal, nbrFaceCells);

    // Transform according to the transformation tensors
    transformCoupleField(pnf);

    // Multiply the field by coefficients and add into the result
    this->addToInternalField(result, !add, coeffs, pnf);
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::updateInterfaceMatrix
(
    const Field<Type>& psiInternal,
    Field<Type>& result,
    const BlockLduMatrix<Type>&,
    const CoeffField<Type>& coeffs,
    const Pstream::commsTypes commsType
) const
{
    NotImplemented;
}


template<class Type>
void Foam::cyclicFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
}


// ************************************************************************* //
