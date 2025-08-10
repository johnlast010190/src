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
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "pointFieldDecomposer.H"
#include "fields/pointPatchFields/constraint/processor/processorPointPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::pointFieldDecomposer::decomposeField
(
    const GeometricField<Type, pointPatchField, pointMesh>& field
) const
{
    // Create and map the internal field values
    Field<Type> internalField(field.primitiveField(), pointAddressing_);

    // Create a list of pointers for the patchFields
    PtrList<pointPatchField<Type>> patchFields(procMesh_.boundary().size());

    // Create and map the patch field values
    forAll(procMesh_.boundary(), patchi)
    {
        if (patchi < completeMesh_.boundary().size())
        {
            patchFields.set
            (
                patchi,
                pointPatchField<Type>::New
                (
                    field.boundaryField()[patchi],
                    procMesh_.boundary()[patchi],
                    DimensionedField<Type, pointMesh>::null(),
                    patchFieldDecomposers_[patchi]
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchi,
                new processorPointPatchField<Type>
                (
                    procMesh_.boundary()[patchi],
                    DimensionedField<Type, pointMesh>::null()
                )
            );
        }
    }

    // Create the field for the processor
    return tmp<GeometricField<Type, pointPatchField, pointMesh>>
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                field.name(),
                procMesh_().time().timeName(),
                procMesh_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            procMesh_,
            field.dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class GeoField>
void Foam::pointFieldDecomposer::decomposeFields
(
    const PtrList<GeoField>& fields
) const
{
    forAll(fields, fieldi)
    {
        decomposeField(fields[fieldi])().write();
    }
}


// ************************************************************************* //
