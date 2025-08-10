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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/basic/mixed/mixedFvPatchField.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    referenceFrameFvPatch<Type>(p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{}

template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool readValues
)
:
    fvPatchField<Type>(p, iF, dict, false),
    referenceFrameFvPatch<Type>(dict, p, iF),
    refValue_(p.size()),
    refGrad_(p.size()),
    valueFraction_(p.size())
{
    if (readValues)
    {
        refValue_ = Field<Type>("refValue", dict, p.size());
        refGrad_ = Field<Type>("refGradient", dict, p.size());
        valueFraction_ = scalarField("valueFraction", dict, p.size());

        evaluate();
    }
}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const mixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fvPatchField<Type>(ptf, p, iF, mapper),
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
    refValue_(mapper(ptf.refValue_)),
    refGrad_(mapper(ptf.refGrad_)),
    valueFraction_(mapper(ptf.valueFraction_))
{
    if (!mapper.suppressWarning() && notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const mixedFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    referenceFrameFvPatch<Type>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        ptf.internalField()
    ),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


template<class Type>
Foam::mixedFvPatchField<Type>::mixedFvPatchField
(
    const mixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    referenceFrameFvPatch<Type>
    (
        ptf.patch(),
        ptf.coorFramePtr_,
        ptf.inputValue(),
        ptf.inletFlux_,
        ptf.frameName_,
        iF
    ),
    refValue_(ptf.refValue_),
    refGrad_(ptf.refGrad_),
    valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mixedFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    referenceFrameFvPatch<Type>::autoMap(m);
    m(refValue_, refValue_);
    m(refGrad_, refGrad_);
    m(valueFraction_, valueFraction_);
}


template<class Type>
void Foam::mixedFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);
    referenceFrameFvPatch<Type>::rmap
    (
        dynamic_cast<const referenceFrameFvPatch<Type>&>(ptf),
        addr
    );

    const mixedFvPatchField<Type>& mptf =
        refCast<const mixedFvPatchField<Type>>(ptf);

    refValue_.rmap(mptf.refValue_, addr);
    refGrad_.rmap(mptf.refGrad_, addr);
    valueFraction_.rmap(mptf.valueFraction_, addr);
}

template<class Type>
void Foam::mixedFvPatchField<Type>::autoMapGIB
(
    const gibFvPatchFieldMapper& mapper
)
{
    fvPatchField<Type>::autoMapGIB(mapper);
    referenceFrameFvPatch<Type>::autoMapGIB(mapper);
    mapper.map(refValue_, pTraits<Type>::zero);
    mapper.map(refGrad_, pTraits<Type>::zero);
    mapper.map(valueFraction_, scalar(0));
}


template<class Type>
void Foam::mixedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    Field<Type>::operator=
    (
        valueFraction_*(refValue_ + getFrameVelocity())
      +
        (1.0 - valueFraction_)*
        (
            this->patchInternalField()
          + refGrad_/this->patch().deltaCoeffs()
        )
    );
    fvPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::snGrad() const
{
    return
        valueFraction_
       *(refValue_ + getFrameVelocity() - this->patchInternalField())
       *this->patch().deltaCoeffs()
      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return Type(pTraits<Type>::one)*(1.0 - valueFraction_);

}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
         valueFraction_*(refValue_ + getFrameVelocity())
       + (1.0 - valueFraction_)*refGrad_/this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*valueFraction_*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mixedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        valueFraction_*this->patch().deltaCoeffs()*(refValue_ + getFrameVelocity())
      + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
void Foam::mixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    referenceFrameFvPatch<Type>::write(os);
    refValue_.writeEntry("refValue", os);
    refGrad_.writeEntry("refGradient", os);
    valueFraction_.writeEntry("valueFraction", os);
    this->writeEntry("value", os);
}


// ************************************************************************* //
