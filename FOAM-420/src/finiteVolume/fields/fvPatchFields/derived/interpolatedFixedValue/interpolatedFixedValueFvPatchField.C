/*---------------------------------------------------------------------------*\
| Modified 2010-2012 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
                    |
     F ield         | FOAM: The Open Source CFD Toolbox
     O peration     |
     A nd           | Copyright (C) 2009 Icon CG Ltd.
     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of FOAM.

    FOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    FOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/interpolatedFixedValue/interpolatedFixedValueFvPatchField.H"
#include "fvMesh/wallDist/wallDist/wallDist.H"
#include "db/Time/Time.H"
#include "primitives/Scalar/scalar/scalar.H"
#include "finiteVolume/fvc/fvcGrad.H"


namespace Foam
{

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
interpolatedFixedValueFvPatchField<Type>::
interpolatedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    distribution_()
{}


template<class Type>
interpolatedFixedValueFvPatchField<Type>::
interpolatedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    distribution_(new fieldProfile<Type>(p, dict))
{
    this->updateCoeffs();
}


template<class Type>
interpolatedFixedValueFvPatchField<Type>::
interpolatedFixedValueFvPatchField
(
    const interpolatedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    distribution_(ptf.distribution_, false)
{}


template<class Type>
interpolatedFixedValueFvPatchField<Type>::
interpolatedFixedValueFvPatchField
(
    const interpolatedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    distribution_(ptf.distribution_, false)
{}


template<class Type>
interpolatedFixedValueFvPatchField<Type>::
interpolatedFixedValueFvPatchField
(
    const interpolatedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    distribution_(ptf.distribution_, false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void interpolatedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (distribution_->update())
    {
        Field<Type> ft(distribution_->value());
        this->frameFieldUpdate(ft);
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
tmp<Field<Type>> interpolatedFixedValueFvPatchField<Type>::gradient() const
{
    tmp<Field<Type>> gradient(distribution_->gradient());
    this->addFrameGradient(gradient.ref());
    return gradient;
}


template<class Type>
void interpolatedFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fixedValueFvPatchField<Type>::write(os);

    if (distribution_.valid())
    {
        distribution_->write(os);
    }
}


} // End namespace Foam

// ************************************************************************* //
