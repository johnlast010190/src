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
    (c) 2015-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/outletFixedValue/outletFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
outletFixedValueFvPatchField<Type>::outletFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    phiName_("phi")
{}


template<class Type>
outletFixedValueFvPatchField<Type>::outletFixedValueFvPatchField
(
    const outletFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
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
outletFixedValueFvPatchField<Type>::outletFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{}


template<class Type>
outletFixedValueFvPatchField<Type>::outletFixedValueFvPatchField
(
    const outletFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_)
{}


template<class Type>
outletFixedValueFvPatchField<Type>::outletFixedValueFvPatchField
(
    const outletFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type>> outletFixedValueFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    const Field<scalar>& phip =
        this->patch().template lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            this->db(),
            phiName_
        );

    tmp<Field<Type>> vic
    (
        new Field<Type>(this->size(), pTraits<Type>::one)
    );

    vic.ref() *= pos0(phip);

    return vic;
}


template<class Type>
tmp<Field<Type>> outletFixedValueFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    const Field<scalar>& phip =
        this->patch().template lookupPatchFieldInDb<surfaceScalarField, scalar>
        (
            this->db(),
            phiName_
        );

    tmp<Field<Type>> vbc
    (
        new Field<Type>(*this)
    );

    vbc.ref() *= neg(phip);

    return vbc;
}

template<class Type>
void outletFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }
    this->writeEntry("value", os);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
