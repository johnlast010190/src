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

#include "fields/fvPatchFields/derived/inletOutlet/inletOutletFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::inletOutletFvPatchField<Type>::fluxValueFraction
(
    const word& phiName,
    const fvPatch& patch,
    scalar scaleFactor
)
{
    const Field<scalar>& phip =
        this->patch().template lookupPatchFieldInDb<surfaceScalarField,
        scalar>
        (
            this->db(),
            phiName
        );

    scalarField fluxThreshold
    (
        scaleFactor*patch.magSf()
        *gSum(mag(phip))/gSum(patch.magSf())
    );

    fluxThreshold = max(scaleFactor*SMALL*patch.magSf(), fluxThreshold);

    tmp<scalarField> vfracTmp(new scalarField(patch.size(), Zero));
    scalarField& vfrac = vfracTmp.ref();

    forAll(phip, fi)
    {
        if (phip[fi] < -fluxThreshold[fi])
        {
            vfrac[fi] = 1.0;
        }
        else if (phip[fi] > fluxThreshold[fi])
        {
            vfrac[fi] = 0.0;
        }
        else
        {
            vfrac[fi] = phip[fi]/(2*fluxThreshold[fi]) + 0.5;
        }
    }

    return vfracTmp;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::inletOutletFvPatchField<Type>::inletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phi"),
    zeroFluxFixedValue_(false),
    fluxTransition_(false)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::inletOutletFvPatchField<Type>::inletOutletFvPatchField
(
    const inletOutletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    zeroFluxFixedValue_(ptf.zeroFluxFixedValue_),
    fluxTransition_(ptf.fluxTransition_)
{}


template<class Type>
Foam::inletOutletFvPatchField<Type>::inletOutletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF, dict, false),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    zeroFluxFixedValue_(dict.lookupOrDefault<Switch>("zeroFluxFixedValue", false)),
    fluxTransition_(dict.lookupOrDefault<Switch>("fluxTransition", false))
{
    referenceFrameFvPatch<Type>::read(dict);
    this->patchType() = dict.lookupOrDefault<word>("patchType", word::null);

    this->refValue() = Field<Type>("inletValue", dict, p.size());

    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::inletOutletFvPatchField<Type>::inletOutletFvPatchField
(
    const inletOutletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_),
    zeroFluxFixedValue_(ptf.zeroFluxFixedValue_),
    fluxTransition_(ptf.fluxTransition_)
{}


template<class Type>
Foam::inletOutletFvPatchField<Type>::inletOutletFvPatchField
(
    const inletOutletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_),
    zeroFluxFixedValue_(ptf.zeroFluxFixedValue_),
    fluxTransition_(ptf.fluxTransition_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::inletOutletFvPatchField<Type>::fixesValue() const
{
    bool fixesValue(false);
    label thisSize(this->size());
    const scalarField& vFrac = this->valueFraction();

    for (label fI=0; fI<thisSize; fI++)
    {
        if (vFrac[fI]>0)
        {
            fixesValue = true;
            break;
        }
    }
    reduce(fixesValue, orOp<bool>());

    return fixesValue;
}


template<class Type>
void Foam::inletOutletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (fluxTransition_)
    {
        this->valueFraction()
            = fluxValueFraction(phiName_, this->patch(), 1e-6)();
    }
    else
    {
        const Field<scalar>& phip =
            this->patch().template lookupPatchFieldInDb<surfaceScalarField,
            scalar>
            (
                this->db(),
                phiName_
            );

        if (zeroFluxFixedValue_)
        {
            this->valueFraction() = 1.0 + neg(phip);
        }
        else
        {
            this->valueFraction() = 1.0 - pos0(phip);
        }
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::inletOutletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }
    if (zeroFluxFixedValue_)
    {
        os.writeEntry("zeroFluxFixedValue", zeroFluxFixedValue_);
    }
    if (fluxTransition_)
    {
        os.writeEntry("fluxTransition", fluxTransition_);
    }

    this->refValue().writeEntry("inletValue", os);
    referenceFrameFvPatch<Type>::write(os);
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::inletOutletFvPatchField<Type>::operator=
(
    const fvPatchField<Type>& ptf
)
{
    fvPatchField<Type>::operator=
    (
        this->valueFraction()*this->refValue()
        + (1 - this->valueFraction())*ptf
    );
}


// ************************************************************************* //
