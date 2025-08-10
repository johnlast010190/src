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
    (c) 2018-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/interpolatedInletOutlet/interpolatedInletOutletFvPatchField.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
interpolatedInletOutletFvPatchField<Type>::
interpolatedInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    distribution_(new fieldProfile<Type>(p)),
    phiName_("phi")
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;
}


template<class Type>
interpolatedInletOutletFvPatchField<Type>::
interpolatedInletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict, false),
    distribution_(new fieldProfile<Type>(p, dict)),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
    this->refValue() = pTraits<Type>::zero;
    this->refGrad() = pTraits<Type>::zero;
    this->valueFraction() = 0.0;

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=
        (
            fvPatchField<Type>::patchInternalField()()
        );
    }
}


template<class Type>
interpolatedInletOutletFvPatchField<Type>::
interpolatedInletOutletFvPatchField
(
    const interpolatedInletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_)
{}


template<class Type>
interpolatedInletOutletFvPatchField<Type>::
interpolatedInletOutletFvPatchField
(
    const interpolatedInletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_)
{}


template<class Type>
interpolatedInletOutletFvPatchField<Type>::
interpolatedInletOutletFvPatchField
(
    const interpolatedInletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    distribution_(ptf.distribution_, false),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void interpolatedInletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (distribution_->update())
    {
        this->refValue() = distribution_->value();
    }

    if (phiName_ != "none")
    {
        const Field<scalar>& phip = this->patch().lookupPatchFieldInDb
        (
            this->db(),
            phiName_,
            reinterpret_cast<const surfaceScalarField*>(0),
            reinterpret_cast<const scalar*>(0)
        );

        this->valueFraction() = 1.0 - pos0(phip);
    }
    else
    {
        this->valueFraction() = 1.0;
    }

    mixedFvPatchField<Type>::updateCoeffs();
}

template<class Type>
tmp<Field<Type>> interpolatedInletOutletFvPatchField<Type>::gradient() const
{
    return distribution_->gradient();
}


template<class Type>
void interpolatedInletOutletFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);

    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }

    distribution_->write(os);

    referenceFrameFvPatch<Type>::write(os);
    this->writeEntry("value", os);
}


} // End namespace Foam

// ************************************************************************* //
