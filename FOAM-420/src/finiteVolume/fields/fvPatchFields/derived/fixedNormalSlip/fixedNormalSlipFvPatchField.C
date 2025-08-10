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
    (c) 2017 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/fixedNormalSlip/fixedNormalSlipFvPatchField.H"
#include "fields/Fields/symmTransformField/symmTransformField.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(p, iF),
    referenceFrameFvPatch<Type>(p, iF),
    fixedValue_(p.size(), Zero),
    relax_(0.5)
{}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchField<Type>(ptf, p, iF, mapper),
    referenceFrameFvPatch<Type>
    (
        p,
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        mapper,
        iF
    ),
    fixedValue_(mapper(ptf.fixedValue_)),
    relax_(ptf.relax_)
{}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchField<Type>(p, iF),
    referenceFrameFvPatch<Type>(dict, p, iF),
    fixedValue_("fixedValue", dict, p.size()),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.5))
{
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    evaluate();
}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf
)
:
    transformFvPatchField<Type>(ptf),
    referenceFrameFvPatch<Type>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        ptf.internalField()
    ),
    fixedValue_(ptf.fixedValue_),
    relax_(ptf.relax_)
{}


template<class Type>
Foam::fixedNormalSlipFvPatchField<Type>::fixedNormalSlipFvPatchField
(
    const fixedNormalSlipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    transformFvPatchField<Type>(ptf, iF),
    referenceFrameFvPatch<Type>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        iF
    ),
    fixedValue_(ptf.fixedValue_),
    relax_(ptf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    transformFvPatchField<Type>::autoMap(m);
    referenceFrameFvPatch<Type>::autoMap(m);
    m(fixedValue_, fixedValue_);
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    transformFvPatchField<Type>::rmap(ptf, addr);

    const fixedNormalSlipFvPatchField<Type>& dmptf =
        refCast<const fixedNormalSlipFvPatchField<Type>>(ptf);

    fixedValue_.rmap(dmptf.fixedValue_, addr);
    if (this->inputValue().size())
    {
        const_cast<Field<Type>&>(this->inputValue()).rmap
        (
            dmptf.inputValue(),
            addr
        );
    }
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    transformFvPatchField<Type>::autoMapGIB(mapper);
    referenceFrameFvPatch<Type>::autoMapGIB(mapper);
    mapper.map(fixedValue_, pTraits<Type>::zero);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedNormalSlipFvPatchField<Type>::snGrad() const
{
    const vectorField nHat(this->patch().nf());
    const Field<Type> pif(this->patchInternalField());
    Field<Type> fixedValue(fixedValue_);
    this->makeVectorGlobal(fixedValue);

    return
    (
        (nHat*(nHat & (fixedValue + getFrameVelocity()))
      + transform(I - sqr(nHat), pif)) - pif
    )*this->patch().deltaCoeffs();
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const vectorField nHat(this->patch().nf());
    Field<Type> fixedValue(fixedValue_);
    this->makeVectorGlobal(fixedValue);

    Field<Type>::operator=
    (
        nHat*(nHat & (fixedValue + getFrameVelocity()))
      + transform(I - sqr(nHat), this->patchInternalField())
    );

    transformFvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedNormalSlipFvPatchField<Type>::valueDivInternalCoeffs
(
    const tmp<Field<scalar>>& w
) const
{
    return tmp<Field<Type>>
    (
        new Field<Type>(this->size(), Zero)
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedNormalSlipFvPatchField<Type>::valueDivBoundaryCoeffs
(
    const tmp<Field<scalar>>& w
) const
{
    const vectorField nHat(this->patch().nf());
    Field<Type> fixedValue(fixedValue_);
    this->makeVectorGlobal(fixedValue);

    return tmp<Field<Type>>
    (
        new Field<Type>(nHat*(nHat & (fixedValue + getFrameVelocity())))
    );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedNormalSlipFvPatchField<Type>::snGradTransformDiag() const
{
    const vectorField nHat(this->patch().nf());
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


template<class Type>
void Foam::fixedNormalSlipFvPatchField<Type>::write(Ostream& os) const
{
    transformFvPatchField<Type>::write(os);
    referenceFrameFvPatch<Type>::write(os);
    fixedValue_.writeEntry("fixedValue", os);
    this->writeEntry("value", os);
    this->template writeEntryIfDifferent<scalar>(os, "relax", 0.5, relax_);
}


// ************************************************************************* //
