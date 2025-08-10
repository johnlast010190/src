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
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "materialModel.H"
#include "materialModels/baseModels/baseModels.H"


template<class Type>
bool Foam::materialModel::isModelConst() const
{
    baseModels<typename Type::ModelType>* mod =
        materialTables_.sTable
        (
            phaseName_, specieName_
        ).lookup(Type::typeName, nullptr);
    return (mod == nullptr || mod->isConst());
}


template<class Type>
Foam::GeometricField<Type, fvPatchField, volMesh>*
Foam::materialModel::lookupPtr
(
    const word& fieldSubName,
    const word& phaseName
)
{
    word fName(phasePropertyName(fieldSubName, phaseName));
    if (phaseName == word::null && phaseName_ != word::null)
    {
        fName = phasePropertyName(fieldSubName, phaseName_);
    }

    typedef GeometricField<Type, fvPatchField, volMesh> GeoField;
    GeoField* fieldPtr = obr_.lookupObjectRefPtr<GeoField>(fName);

    if (fieldPtr == nullptr)
    {
        FatalErrorInFunction
            << "Field with name: " << fName
            << " not found." << nl
            << "Available fields: "
            << obr_.names<GeoField>()
            << exit(FatalError);
    }
    return fieldPtr;
}


template<class Type>
const Foam::referenceFields<Type>*
Foam::materialModel::constructOrReturnRefFieldPtr(const word& fieldName) const
{
    const dictionary& dict =
        materialsDict().optionalSubDict("referenceFields");

    typedef referenceFields<Type> refField;
    const word phaseFieldName = phasePropertyName(fieldName, phaseName_);
    if (!obr_.foundObject<refField>(phaseFieldName + "Ref"))
    {
        regIOobject::store(new refField(obr_, phaseFieldName, dict));
    }
    return obr_.lookupObjectPtr<refField>(phaseFieldName + "Ref");
}


// ************************************************************************* //
