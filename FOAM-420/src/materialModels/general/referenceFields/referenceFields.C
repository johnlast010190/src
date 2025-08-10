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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "referenceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::referenceFields<Type>::referenceFields
(
    const objectRegistry& obr,
    const word& fieldName,
    const dictionary& dict
)
:
    regIOobject(IOobject(fieldName + "Ref", obr.time().timeName(), obr)),
    fieldName_(fieldName),
    refFieldName_(fieldName + "Ref"),
    fieldRef_(dict.lookupOrDefault<dimensioned<Type>>(fieldName, dimensioned<Type>())),
    isConst_(dict.found(fieldName + "Const")),
    constVal_(dict.lookupOrDefault<dimensioned<Type>>(fieldName + "Const", dimensioned<Type>())),
    field_
    (
        obr.lookupObjectPtr<GeometricField<Type, fvPatchField, volMesh>>
        (
            fieldName
        )
    ),
    mesh_
    (
        isA<fvSolutionRegistry>(obr)
      ? dynamic_cast<const fvSolutionRegistry&>(obr).mesh()
      : dynamic_cast<const fvMesh&>(obr)
    ),
    boundary_(*this)
{
    if (field_ == nullptr && !isConst_)
    {
        FatalErrorInFunction
            << "Couldn't find field \"" << fieldName << "\" and "
            << fieldName << "Const is set to \"false\"." << nl
            << exit(FatalError);
    }
    else if (field_ == nullptr && isConst_)
    {
        // Type check => needs to have specified dimmensions
        dict.lookup<dimensioned<Type>>(fieldName + "Const");
    }
    else if (field_ != nullptr)
    {
        if (dict.found(fieldName_) && fieldRef_.dimensions() != field_->dimensions())
        {
            FatalErrorInFunction
                << "Dimensions of the field \"" << fieldName << "\""
                << "and specified reference are not matching." << nl
                << exit(FatalError);
        }
        // Reset dimensions to correct ones but the reference not
        else if (field_->dimensions() != fieldRef_.dimensions())
        {
            fieldRef_.dimensions().reset(field_->dimensions());
            fieldRef_.name() = refFieldName_;
        }
    }
    fieldRef_.name() = refFieldName_;
    constVal_.name() = fieldName + "Const";
    Info<< fieldRef_ << endl;
    if (isConst_)
    {
        Info<< constVal_ << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::referenceFields<Type>::primitiveField() const
{
    if (field_ != nullptr && !isConst_)
    {
        return field_->primitiveField() + fieldRef_.value();
    }
    return
        tmp<Field<Type>>
        (
            new Field<Type>(mesh_.nCells(), fieldRef_.value() + constVal_.value())
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::referenceFields<Type>::patchField(const label patchi) const
{
    if (field_ != nullptr && !isConst_)
    {
        return field_->boundaryField()[patchi] + fieldRef_.value();
    }
    return
        tmp<Field<Type>>
        (
            new Field<Type>
            (
                mesh_.boundary()[patchi].size(),
                fieldRef_.value() + constVal_.value()
            )
        );
}


template<class Type>
const Foam::dimensioned<Type>& Foam::referenceFields<Type>::offset() const
{
    return fieldRef_;
}


template<class Type>
const Foam::dimensioned<Type>& Foam::referenceFields<Type>::constant() const
{
    return constVal_;
}


template<class Type>
bool Foam::referenceFields<Type>::isConst() const
{
    return isConst_;
}


template<class Type>
void Foam::referenceFields<Type>::makeConst(const Type& value)
{
    isConst_ = true;
    constVal_.value() = value;
}


template<class Type>
void Foam::referenceFields<Type>::makeNonConst()
{
    isConst_ = false;
}


template<class Type>
bool Foam::referenceFields<Type>::updateScalarField
(
    const word& fieldName,
    const volScalarField& volField
)
{
    if (fieldName_ == fieldName && !isConst_)
    {
        field_ = &volField;
        return true;
    }
    return false;
}


template<class Type>
bool Foam::referenceFields<Type>::updateVectorField
(
    const word& fieldName,
    const volVectorField& volField
)
{
    if (fieldName_ == fieldName && !isConst_)
    {
        field_ = &volField;
        return true;
    }
    return false;
}


template<class Type>
bool Foam::referenceFields<Type>::updateTensorField
(
    const word& fieldName,
    const volTensorField& volField
)
{
    if (fieldName_ == fieldName && !isConst_)
    {
        field_ = &volField;
        return true;
    }
    return false;
}


template<class Type>
bool Foam::referenceFields<Type>::writeData(Ostream& os) const
{
    NotImplemented;
}


template<class Type>
bool Foam::referenceFields<Type>::write(const bool valid) const
{
    NotImplemented;
}


// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::referenceFields<Type>::operator()() const
{
    if (field_ != nullptr && !isConst_)
    {
        return
            tmp<geoField>(new geoField((*field_))) + fieldRef_;
    }

    return
        tmp<geoField>
        (
            new geoField
            (
                IOobject
                (
                    refFieldName_,
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                fieldRef_ + constVal_
            )
        );
}


// ************************************************************************* //
