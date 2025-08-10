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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvsPatchFields/constraint/cyclicAMI/cyclicAMIFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(p, iF),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{}


template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const cyclicAMIFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    coupledFvsPatchField<Type>(ptf, p, iF, mapper),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    coupledFvsPatchField<Type>(p, iF, dict),
    cyclicAMIPatch_(refCast<const cyclicAMIFvPatch>(p))
{
    if (!isA<cyclicAMIFvPatch>(p))
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "patch " << this->patch().index() << " not cyclicAMI type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const cyclicAMIFvsPatchField<Type>& ptf
)
:
    coupledFvsPatchField<Type>(ptf),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


template<class Type>
Foam::cyclicAMIFvsPatchField<Type>::cyclicAMIFvsPatchField
(
    const cyclicAMIFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    coupledFvsPatchField<Type>(ptf, iF),
    cyclicAMIPatch_(ptf.cyclicAMIPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::cyclicAMIFvsPatchField<Type>::coupled() const
{
    if
    (
        Pstream::parRun()
     || (
            this->cyclicAMIPatch_.size()
         && this->cyclicAMIPatch_.cyclicAMIPatch().nbrPatch().size()
        )
    )
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::cyclicAMIFvsPatchField<Type>::patchNeighbourField
(
    const Pstream::commsTypes commsType
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> geoField;
    const geoField& gf = refCast<const geoField>(this->internalField());

    const cyclicAMIFvPatch& cp =
        refCast<const cyclicAMIFvPatch>(this->patch());

    const Field<Type>& pnf = gf.boundaryField()[cp.nbrPatchID()];

    tmp<Field<Type>> tpnf;
    if (cp.applyLowWeightCorrection())
    {
        tpnf = cyclicAMIPatch_.interpolate(pnf, *this);
    }
    else
    {
        tpnf = cyclicAMIPatch_.interpolate(pnf);
    }

    cp.transform().transform(tpnf.ref(), tpnf());

    return tpnf;
}


// ************************************************************************* //
