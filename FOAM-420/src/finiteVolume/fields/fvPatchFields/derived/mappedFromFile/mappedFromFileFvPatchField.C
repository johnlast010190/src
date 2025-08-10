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
    (c) 1991-2009 OpenCFD Ltd.
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fields/fvPatchFields/derived/mappedFromFile/mappedFromFileFvPatchField.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFieldMapper.H"
#include "fields/volFields/volFields.H"
#include "containers/Lists/DynamicList/DynamicList.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::mappedFromFileFvPatchField<Type>::mapField
(
    const pointField& targetPoints,
    const vectorField& targetNormals
)
{
    FatalErrorInFunction
    << "Function only defined for scalar types."
    << exit(FatalError);

    return
        tmp<Field<Type>>
        (
            Field<Type>(this->patch().size(), pTraits<Type>::zero)
        );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class Type>
Foam::mappedFromFileFvPatchField<Type>::mappedFromFileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    mappedData_(false)
{}


template<class Type>
Foam::mappedFromFileFvPatchField<Type>::mappedFromFileFvPatchField
(
    const mappedFromFileFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    fileData_(ptf.fileData_),
    mappedData_(false)
{}


template<class Type>
Foam::mappedFromFileFvPatchField<Type>::mappedFromFileFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    fileData_(dict),
    mappedData_(false)
{
    /* This was functional
    tmp<Field<Type>> mappedField
    (
        mapField
        (
            this->patch().Cf(),
            this->patch().Sf()
        )
    );

    if (mappedField.valid())
    {
        this->operator==(mappedField());
    }
    */
}


template<class Type>
Foam::mappedFromFileFvPatchField<Type>::mappedFromFileFvPatchField
(
    const mappedFromFileFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    fileData_(),
    mappedData_(false)
{}


template<class Type>
Foam::mappedFromFileFvPatchField<Type>::mappedFromFileFvPatchField
(
    const mappedFromFileFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    fileData_(),
    mappedData_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::mappedFromFileFvPatchField<Type>::updateCoeffs()
{

    if (this->updated())
    {
        return;
    }

    if (!mappedData_)
    {
        tmp<Field<Type>> mappedField
        (
            mapField
            (
                this->patch().Cf(),
                this->patch().Sf()
            )
        );

        if (mappedField.valid())
        {
            this->forceAssign(mappedField());
        }

        mappedData_ = true;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedFromFileFvPatchField<Type>::write(Ostream& os) const
{
    this->fvPatchField<Type>::write(os);
    os.writeEntry("file", fileData_.filename());
    os.writeEntry("field", fileData_.fieldName());
    os.writeEntry("mapBack", fileData_.mapBack());
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// ************************************************************************* //
