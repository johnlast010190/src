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

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedFixedInternalValue/mappedFixedInternalValueFvPatchField.H"
#include "containers/Lists/UIndirectList/UIndirectList.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mappedFixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const mappedFixedInternalValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mappedFixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mappedFixedValueFvPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const mappedFixedInternalValueFvPatchField<Type>& ptf
)
:
    mappedFixedValueFvPatchField<Type>(ptf)
{}


template<class Type>
Foam::mappedFixedInternalValueFvPatchField<Type>::
mappedFixedInternalValueFvPatchField
(
    const mappedFixedInternalValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mappedFixedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFixedInternalValueFvPatchField<Type>::updateCoeffs()
{
    typedef GeometricField<Type, fvPatchField, volMesh> FieldType;

    if (this->updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    // Retrieve the neighbour values and assign to this patch boundary field
    mappedFixedValueFvPatchField<Type>::updateCoeffs();

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp =
        refCast<const mappedPatchBase>(this->patch().patch());
    const fvMesh& nbrMesh = refCast<const fvMesh>(mpp.sampleMesh());

    Field<Type> nbrIntFld;

    switch (mpp.mode())
    {
        case mappedPatchBase::NEARESICELL:
        {
            FatalErrorInFunction
                << "Cannot apply "
                << mappedPatchBase::sampleModeNames_
                   [
                       mappedPatchBase::NEARESICELL
                   ]
                << " mapping mode for patch " << this->patch().name()
                << exit(FatalError);

            break;
        }
        case mappedPatchBase::NEARESIPATCHFACE:
        case mappedPatchBase::NEARESIPATCHFACEAMI:
        {
            const label samplePatchi = mpp.samplePolyPatch().index();
            const fvPatchField<Type>& nbrPatchField =
                this->sampleField().boundaryField()[samplePatchi];
            nbrIntFld = nbrPatchField.patchInternalField();
            mpp.distribute(nbrIntFld);

            break;
        }
        case mappedPatchBase::NEARESIFACE:
        {
            Field<Type> allValues(nbrMesh.nFaces(), Zero);

            const FieldType& nbrField = this->sampleField();

            forAll(nbrField.boundaryField(), patchi)
            {
                const fvPatchField<Type>& pf = nbrField.boundaryField()[patchi];
                const Field<Type> pif(pf.patchInternalField());

                label faceStart = pf.patch().start();

                forAll(pf, facei)
                {
                    allValues[faceStart++] = pif[facei];
                }
            }

            mpp.distribute(allValues);
            nbrIntFld.transfer(allValues);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown sampling mode: " << mpp.mode()
                << abort(FatalError);
        }
    }

    // Restore tag
    UPstream::msgType() = oldTag;

    // Assign to (this) patch internal field its neighbour values
    Field<Type>& intFld = const_cast<Field<Type>&>(this->primitiveField());
    UIndirectList<Type>(intFld, this->patch().faceCells()) = nbrIntFld;
}


template<class Type>
void Foam::mappedFixedInternalValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    mappedFixedValueFvPatchField<Type>::write(os);
}


// ************************************************************************* //
