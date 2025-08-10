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
    (c) 2012 OpenFOAM Foundation
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/patchMean/patchMeanFvPatchField.H"
#include "fields/volFields/volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::patchMeanFvPatchField<Type>::patchMeanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    sourcePatchRegExp_(),
    patchIDs_(0),
    offset_()
{}


template<class Type>
Foam::patchMeanFvPatchField<Type>::patchMeanFvPatchField
(
    const patchMeanFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    sourcePatchRegExp_(ptf.sourcePatchRegExp_),
    patchIDs_(ptf.patchIDs_),
    offset_
    (
        ptf.offset_.valid()
      ? ptf.offset_().clone().ptr()
      : nullptr
    )
{}


template<class Type>
Foam::patchMeanFvPatchField<Type>::patchMeanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    sourcePatchRegExp_(dict.lookup("patches")),
    patchIDs_(p.patch().boundaryMesh().patchSet(sourcePatchRegExp_)),
    offset_()
{
    if (dict.found("offset"))
    {
        offset_ = Function1<Type>::New("offset", dict);
    }
}


template<class Type>
Foam::patchMeanFvPatchField<Type>::patchMeanFvPatchField
(
    const patchMeanFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    sourcePatchRegExp_(ptf.sourcePatchRegExp_),
    patchIDs_(ptf.patchIDs_),
    offset_
    (
        ptf.offset_.valid()
      ? ptf.offset_().clone().ptr()
      : nullptr
    )
{}


template<class Type>
Foam::patchMeanFvPatchField<Type>::patchMeanFvPatchField
(
    const patchMeanFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    sourcePatchRegExp_(ptf.sourcePatchRegExp_),
    patchIDs_(ptf.patchIDs_),
    offset_
    (
        ptf.offset_.valid()
      ? ptf.offset_().clone().ptr()
      : nullptr
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::patchMeanFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    scalar totalArea = 0;
    Type meanValue = Zero;

    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;
    const volFieldType& cvf
    (
        dynamic_cast<const volFieldType&>(this->internalField())
    );

    forAllConstIter(labelHashSet, patchIDs_, iter)
    {
        const label patchi = iter.key();
        if (patchi == this->patch().index())
        {
            FatalErrorInFunction
                << "Cannot use patchMean on the same patch "
                << this->patch().name() << exit(FatalError);
        }

        totalArea += gSum(this->patch().boundaryMesh()[patchi].magSf());
        meanValue +=
            gSum
            (
                this->patch().boundaryMesh()[patchi].magSf()
               *cvf.boundaryField()[patchi]
            );
    }

    meanValue /= totalArea;

    if (offset_.valid())
    {
        meanValue += offset_->value(this->db().time().timeOutputValue());
    }

    // In general the boundary will already in the frame since
    // the boundaries that we are mapping from are already in the
    // frame. However, the frame setup could still be used to rotate
    // the velocity boundary field.
    Field<Type> ft(this->patch().size(), meanValue);
    this->frameFieldUpdate(ft);

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::patchMeanFvPatchField<Type>::write(Ostream& os) const
{
    fixedValueFvPatchField<Type>::write(os);
    os.writeEntry("patches", sourcePatchRegExp_);

    if (offset_.valid())
    {
        offset_->writeData(os);
    }
}


// ************************************************************************* //
