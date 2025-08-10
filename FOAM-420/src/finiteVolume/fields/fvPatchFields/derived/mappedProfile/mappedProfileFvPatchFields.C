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

#include "fields/fvPatchFields/derived/mappedProfile/mappedProfileFvPatchFields.H"
#include "volMesh/volMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOstreams/Fstreams/IFstream.H"
#include "fields/fvPatchFields/derived/mappedProfile/mappedProfileFvPatchField.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeTemplatePatchTypeField(scalar,mappedProfile);
makeTemplatePatchTypeField(vector,mappedProfile);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template <>
tmp<Field<scalar>> mappedProfileFvPatchField<scalar>::mapProfileToField
(
    const pointField& targetPoints,
    word fieldName,
    word regionName
)
{
    label inputIndex = -1;
    forAll(inputFileNameList_, ilI)
    {
        if (inputFileNameList_[ilI] == proFile_)
        {
            inputIndex = ilI;
        }
    }

    if (inputIndex == -1)
    {
        inputIndex = inputFileNameList_.size();
        inputFileNameList_.setSize(inputIndex+1);
        inputFileNameList_[inputIndex] = proFile_;

        IFstream is(proFile_);
        if (!is.good())
        {
            FatalIOErrorIn
            (
        "mappedProfileFvPatchField::mappedProfileFvPatchField(..)",
                is
             )
                << "Istream not OK for reading profile data "
                << exit(FatalIOError);
        }

        inputData_.setSize(inputFileNameList_.size());
        inputData_.set(inputIndex, new pointProfileDataList(is));
    }


    return tmp<Field<scalar>>
    (
        inputData_[inputIndex].mapData
        (
            targetPoints,
            fieldName,
            regionName
        )
    );
}

template <>
tmp<Field<vector>> mappedProfileFvPatchField<vector>::mapProfileToField
(
    const pointField& targetPoints,
    word fieldName,
    word regionName
)
{
    label inputIndex = -1;
    forAll(inputFileNameList_, ilI)
    {
        if (inputFileNameList_[ilI] == proFile_)
        {
            inputIndex = ilI;
        }
    }

    if (inputIndex == -1)
    {
        inputIndex = inputFileNameList_.size();
        inputFileNameList_.setSize(inputIndex+1);
        inputFileNameList_[inputIndex] = proFile_;

        IFstream is(proFile_);
        if (!is.good())
        {
            FatalIOErrorIn
            (
        "mappedProfileFvPatchField::mappedProfileFvPatchField(..)",
                is
             )
                << "Istream not OK for reading profile data "
                << exit(FatalIOError);
        }

        inputData_.setSize(inputFileNameList_.size());
        inputData_.set(inputIndex, new pointProfileDataList(is));
    }


    tmp<Field<scalar>> xcomponent
    (
        inputData_[inputIndex].mapData
        (
            targetPoints,
            "x-" + fieldName,
            regionName
        )
    );
    tmp<Field<scalar>> ycomponent
    (
        inputData_[inputIndex].mapData
        (
            targetPoints,
            "y-" + fieldName,
            regionName
        )
    );
    tmp<Field<scalar>> zcomponent
    (
        inputData_[inputIndex].mapData
        (
            targetPoints,
            "z-" + fieldName,
            regionName
        )
    );

    tmp<Field<vector>> vf(new vectorField(this->size(), vector::zero));

    vf->Field<vector>::replace(0, xcomponent());
    vf->Field<vector>::replace(1, ycomponent());
    vf->Field<vector>::replace(2, zcomponent());

    return vf;
}


template<>
void Foam::mappedProfileFvPatchField<scalar>::setTypeValue
(
    DynamicList<scalar>& field,
    const List<label>& locField,
    const List<string>& splitted
)
{
    scalar p0 = readScalar(IStringStream(splitted[locField[0]])());
    field.append(p0);
}


template<>
void Foam::mappedProfileFvPatchField<scalar>::transformTypeValue
(
    scalar& data,
    const vector& global
)
{
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
