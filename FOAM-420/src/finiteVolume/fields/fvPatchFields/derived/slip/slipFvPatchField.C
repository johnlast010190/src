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
    (c) 2011-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/slip/slipFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::slipFvPatchField<Type>::slipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    basicSymmetryFvPatchField<Type>(p, iF),
    referenceFrameFvPatch<Type>(p, iF),
    relax_(0.5)
{}


template<class Type>
Foam::slipFvPatchField<Type>::slipFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    basicSymmetryFvPatchField<Type>(p, iF, dict),
    referenceFrameFvPatch<Type>(dict, p, iF),
    relax_(dict.lookupOrDefault<scalar>("relax", 0.5))
{
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);
    evaluate();
}


template<class Type>
Foam::slipFvPatchField<Type>::slipFvPatchField
(
    const slipFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    basicSymmetryFvPatchField<Type>(ptf, p, iF, mapper),
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
    relax_(ptf.relax_)
{}


template<class Type>
Foam::slipFvPatchField<Type>::slipFvPatchField
(
    const slipFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    basicSymmetryFvPatchField<Type>(ptf, iF),
    referenceFrameFvPatch<Type>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        iF
    ),
    relax_(ptf.relax_)
{}


template<class Type>
Foam::slipFvPatchField<Type>::slipFvPatchField
(
    const slipFvPatchField<Type>& ptf
)
:
    basicSymmetryFvPatchField<Type>(ptf),
    referenceFrameFvPatch<Type>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        ptf.internalField()
    ),
    relax_(ptf.relax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::slipFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    transformFvPatchField<Type>::autoMap(m);
    referenceFrameFvPatch<Type>::autoMap(m);
}


template<class Type>
void Foam::slipFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    transformFvPatchField<Type>::rmap(ptf, addr);

    const slipFvPatchField<Type>& dmptf =
        refCast<const slipFvPatchField<Type>>(ptf);

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
void Foam::slipFvPatchField<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    transformFvPatchField<Type>::autoMapGIB(mapper);
    referenceFrameFvPatch<Type>::autoMapGIB(mapper);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slipFvPatchField<Type>::snGrad() const
{
    const vectorField nHat(this->patch().nf());

    const Field<Type> iF(this->patchInternalField());
    const Field<Type> frameVel(getFrameVelocity());

    return
        (transform(2.0*sqr(nHat), frameVel)
       + transform(I - 2.0*sqr(nHat), iF) - iF)
       *(this->patch().deltaCoeffs()/2.0);
}


template<class Type>
void Foam::slipFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    const vectorField nHat(this->patch().nf());
    const Field<Type> frameVel(getFrameVelocity());

    Field<Type>::operator=
    (
        transform(sqr(nHat), frameVel)
      + transform(I - sqr(nHat), this->patchInternalField())
    );

    transformFvPatchField<Type>::evaluate();
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slipFvPatchField<Type>::valueDivInternalCoeffs
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
Foam::slipFvPatchField<Type>::valueDivBoundaryCoeffs
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
void Foam::slipFvPatchField<Type>::write(Ostream& os) const
{
    transformFvPatchField<Type>::write(os);
    referenceFrameFvPatch<Type>::write(os);
    this->writeEntry("value", os);
    this->template writeEntryIfDifferent<scalar>(os, "relax", 0.5, relax_);
}

// ************************************************************************* //
