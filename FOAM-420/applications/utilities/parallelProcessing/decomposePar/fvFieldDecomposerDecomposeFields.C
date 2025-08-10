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

#include "fvFieldDecomposer.H"
#include "fields/fvPatchFields/constraint/processor/processorFvPatchField.H"
#include "fields/fvsPatchFields/constraint/processor/processorFvsPatchField.H"
#include "fields/fvPatchFields/constraint/processorCyclic/processorCyclicFvPatchField.H"
#include "fields/fvsPatchFields/constraint/processorCyclic/processorCyclicFvsPatchField.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvFieldDecomposer::mapCellToFace
(
    const labelUList& owner,
    const labelUList& neighbour,
    const Field<Type>& field,
    const labelUList& addressing
)
{
    tmp<Field<Type>> tfld(new Field<Type>(addressing.size()));
    Field<Type>& fld = tfld.ref();

    forAll(addressing, i)
    {
        fld[i] =
            field
            [
                addressing[i] > 0
              ? neighbour[addressing[i] - 1]
              : owner[- addressing[i] + 1]
            ];
    }

    return tfld;
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvFieldDecomposer::mapFaceToFace
(
    const Field<Type>& field,
    const labelUList& addressing,
    const bool isFlux
)
{
    tmp<Field<Type>> tfld(new Field<Type>(addressing.size()));
    Field<Type>& fld = tfld.ref();

    forAll(addressing, i)
    {
        fld[i] =
            (isFlux && addressing[i] < 0 ? -1 : +1)
           *field[mag(addressing[i]) - 1];
    }

    return tfld;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    const bool allowUnknownPatchFields
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    // Create dummy patch fields
    PtrList<fvPatchField<Type>> patchFields(procMesh_.boundary().size());
    forAll(procMesh_.boundary(), procPatchi)
    {
        patchFields.set
        (
            procPatchi,
            fvPatchField<Type>::New
            (
                calculatedFvPatchField<Type>::typeName,
                procMesh_.boundary()[procPatchi],
                DimensionedField<Type, volMesh>::null()
            )
        );
    }

    // Create the processor field with the dummy patch fields
    tmp<VolFieldType> tresF
    (
        new VolFieldType
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            procMesh_,
            field.dimensions(),
            Field<Type>(field.primitiveField(), cellAddressing_),
            patchFields
        )
    );
    VolFieldType& resF = tresF.ref();
    resF.oriented() = field().oriented();

    // Change the patch fields to the correct type using a mapper constructor
    // (with reference to the now correct internal field)
    typename VolFieldType::Boundary& bf = resF.boundaryFieldRef();

    forAll(bf, procPatchi)
    {
        const fvPatch& procPatch = procMesh_.boundary()[procPatchi];

        // Determine the index of the corresponding complete patch
        const label completePatchi = completePatchID(procPatchi);

        if (completePatchi == procPatchi)
        {
            bf.set
            (
                procPatchi,
                fvPatchField<Type>::New
                (
                    field.boundaryField()[completePatchi],
                    procPatch,
                    resF(),
                    patchFieldDecomposers_[procPatchi]
                )
            );
        }
        else if (isA<processorCyclicFvPatch>(procPatch))
        {
            const label nbrCompletePatchi =
                refCast<const processorCyclicFvPatch>(procPatch)
               .referPatch().nbrPatchID();

            bf.set
            (
                procPatchi,
                new processorCyclicFvPatchField<Type>
                (
                    procPatch,
                    resF(),
                    mapCellToFace
                    (
                        labelUList(),
                        completeMesh_.lduAddr().patchAddr(nbrCompletePatchi),
                        field.primitiveField(),
                        faceAddressingBf_[procPatchi]
                    )
                )
            );
        }
        else if (isA<processorFvPatch>(procPatch))
        {
            bf.set
            (
                procPatchi,
                new processorFvPatchField<Type>
                (
                    procPatch,
                    resF(),
                    mapCellToFace
                    (
                        completeMesh_.owner(),
                        completeMesh_.neighbour(),
                        field.primitiveField(),
                        faceAddressingBf_[procPatchi]
                    )
                )
            );
        }
        else if (allowUnknownPatchFields)
        {
            bf.set
            (
                procPatchi,
                new emptyFvPatchField<Type>(procPatch, resF())
            );
        }
        else if (isA<indirectPolyPatch>(procPatch.patch()))
        {
            label completeIndID = completeMesh_.boundaryMesh().findPatchID
            (
                procPatch.patch().name()
            );

            bf.set
            (
                procPatchi,
                fvPatchField<Type>::New
                (
                    field.boundaryField()[completeIndID],
                    procPatch,
                    resF(),
                    indirectPatchFieldDecomposers_[procPatchi]
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unknown type." << abort(FatalError);
        }
    }

    return tresF;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    const SubList<label> faceAddressingIf
    (
        faceAddressing_,
        procMesh_.nInternalFaces()
    );

    // Create dummy patch fields
    PtrList<fvsPatchField<Type>> patchFields(procMesh_.boundary().size());
    forAll(procMesh_.boundary(), procPatchi)
    {
        patchFields.set
        (
            procPatchi,
            fvsPatchField<Type>::New
            (
                calculatedFvsPatchField<Type>::typeName,
                procMesh_.boundary()[procPatchi],
                DimensionedField<Type, surfaceMesh>::null()
            )
        );
    }

    // Create the processor field with the dummy patch fields
    tmp<SurfaceFieldType> tresF
    (
        new SurfaceFieldType
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            procMesh_,
            field.dimensions(),
            mapFaceToFace
            (
                field,
                faceAddressingIf,
                isFlux(field)
            ),
            patchFields
        )
    );
    SurfaceFieldType& resF = tresF.ref();
    resF.oriented() = field().oriented();

    // Change the patch fields to the correct type using a mapper constructor
    // (with reference to the now correct internal field)
    typename SurfaceFieldType::Boundary& bf = resF.boundaryFieldRef();

    forAll(procMesh_.boundary(), procPatchi)
    {
        const fvPatch& procPatch = procMesh_.boundary()[procPatchi];

        const label completePatchi = completePatchID(procPatchi);

        if (completePatchi == procPatchi)
        {
            bf.set
            (
                procPatchi,
                fvsPatchField<Type>::New
                (
                    field.boundaryField()[procPatchi],
                    procPatch,
                    resF(),
                    patchFieldDecomposers_[procPatchi]
                )
            );
        }
        else if (isA<processorCyclicFvPatch>(procPatch))
        {
            bf.set
            (
                procPatchi,
                new processorCyclicFvsPatchField<Type>
                (
                    procPatch,
                    resF(),
                    mapFaceToFace
                    (
                        field.boundaryField()[completePatchi],
                        faceAddressingBf_[procPatchi],
                        isFlux(field)
                    )
                )
            );
        }
        else if (isA<processorFvPatch>(procPatch))
        {
            bf.set
            (
                procPatchi,
                new processorFvsPatchField<Type>
                (
                    procPatch,
                    resF(),
                    mapFaceToFace
                    (
                        field.primitiveField(),
                        faceAddressingBf_[procPatchi],
                        isFlux(field)
                    )
                )
            );
        }
        else if (isA<indirectPolyPatch>(procPatch.patch()))
        {
            label completeIndID = completeMesh_.boundaryMesh().findPatchID
            (
                procPatch.patch().name()
            );

            bf.set
            (
                procPatchi,
                fvsPatchField<Type>::New
                (
                    field.boundaryField()[completeIndID],
                    procPatch,
                    resF(),
                    indirectPatchFieldDecomposers_[procPatchi]
                )
            );
        }
        else
        {
            FatalErrorInFunction
                << "Unknown type." << abort(FatalError);
        }
    }

    return tresF;
}


template<class GeoField>
void Foam::fvFieldDecomposer::decomposeFields
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
