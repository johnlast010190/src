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
    (c) 2018 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/fixedValueZone/fixedValueZoneFvPatchField.H"
#include "fields/volFields/volFields.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fixedValueZoneFvPatchField<Type>::fixedValueZoneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    cellZoneValue_(pTraits<Type>::zero)
{}


template<class Type>
fixedValueZoneFvPatchField<Type>::fixedValueZoneFvPatchField
(
    const fixedValueZoneFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    cellZoneValue_(ptf.cellZoneValue_)
{}


template<class Type>
fixedValueZoneFvPatchField<Type>::fixedValueZoneFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    cellZoneValue_
    (
        dict.lookupOrDefault<Type>
        (
            "cellZoneValue", pTraits<Type>::zero
        )
    )
{}


template<class Type>
fixedValueZoneFvPatchField<Type>::fixedValueZoneFvPatchField
(
    const fixedValueZoneFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    cellZoneValue_(ptf.cellZoneValue_)
{}


template<class Type>
fixedValueZoneFvPatchField<Type>::fixedValueZoneFvPatchField
(
    const fixedValueZoneFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    cellZoneValue_(ptf.cellZoneValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fixedValueZoneFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Type>::forceAssign(Field<Type>(this->size(), cellZoneValue_));
    fvPatchField<Type>::evaluate();
}


template<class Type>
void Foam::fixedValueZoneFvPatchField<Type>::manipulateMatrix
(
    fvMatrix<Type>& matrix
)
{
    const polyPatch& pthis = this->patch().patch();
    const polyMesh& pMesh =  pthis.boundaryMesh().mesh();

    const indirectPolyPatch& gibPolyPatch =
        refCast<const indirectPolyPatch>(pthis);

    const word faceZoneName = pMesh.faceZones()[gibPolyPatch.zoneId()].name();

    const label& zoneId = pMesh.cellZones().findZoneID
    (
        "inactive_"+faceZoneName
    );

    matrix.setValues
    (
        pMesh.cellZones()[zoneId],
        Field<Type>(pMesh.cellZones()[zoneId].size(), cellZoneValue_)
    );
    fvPatchField<Type>::manipulateMatrix(matrix);
}


template<class Type>
void Foam::fixedValueZoneFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    this->writeEntry("value", os);
    os.writeEntry("cellZoneValue", cellZoneValue_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
