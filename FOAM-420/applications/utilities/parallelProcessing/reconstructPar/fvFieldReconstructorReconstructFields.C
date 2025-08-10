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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvFieldReconstructor.H"
#include "db/Time/Time.H"
#include "containers/Lists/PtrList/PtrList.H"
#include "fields/fvPatchFields/fvPatchField/fvPatchFields.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fields/fvPatchFields/constraint/empty/emptyFvPatchField.H"
#include "fields/fvsPatchFields/constraint/empty/emptyFvsPatchField.H"
#include "fvMesh/fvPatches/constraint/processorCyclic/processorCyclicFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fvFieldReconstructor::rmapFaceToFace
(
    Field<Type>& toField,
    const Field<Type>& fromField,
    const labelUList& addressing,
    const bool isFlux
)
{
    forAll(addressing, i)
    {
        toField[mag(addressing[i]) - 1] =
            (isFlux && addressing[i] < 0 ? -1 : +1)*fromField[i];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructFvVolumeInternalField
(
    const IOobject& fieldIoObject,
    const PtrList<DimensionedField<Type, volMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(completeMesh_.nCells());

    forAll(procMeshes_, proci)
    {
        const DimensionedField<Type, volMesh>& procField = procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.field(),
            cellProcAddressing_[proci]
        );
    }

    tmp<DimensionedField<Type, volMesh>> tfield
    (
        new DimensionedField<Type, volMesh>
        (
            fieldIoObject,
            completeMesh_,
            procFields[0].dimensions(),
            internalField
        )
    );

    tfield.ref().oriented() = procFields[0].oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructFvVolumeInternalField
(
    const IOobject& fieldIoObject
) const
{
    PtrList<DimensionedField<Type, volMesh>>
        procFields(procMeshes_.size());

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new DimensionedField<Type, volMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().timeName(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMeshes_[proci]
            )
        );
    }

    return reconstructFvVolumeInternalField
    (
        IOobject
        (
            fieldIoObject.name(),
            completeMesh_.time().timeName(),
            completeMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        procFields
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject,
    const PtrList<GeometricField<Type, fvPatchField, volMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(completeMesh_.nCells());

    // Create the patch fields
    PtrList<fvPatchField<Type>> patchFields(completeMesh_.boundary().size());

    forAll(procFields, proci)
    {
        const GeometricField<Type, fvPatchField, volMesh>& procField =
            procFields[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.primitiveField(),
            cellProcAddressing_[proci]
        );

        // Set the boundary patch values in the reconstructed field
        forAll(procField.boundaryField(), procPatchi)
        {
            const fvPatch& procPatch =
                procMeshes_[proci].boundary()[procPatchi];

            // Determine the index of the corresponding complete patch
            const label completePatchi = completePatchID(proci, procPatchi);

            // Check if the boundary patch is not a processor patch
            // nor an indirect patch
            if (completePatchi == procPatchi)
            {
                if (!patchFields(completePatchi))
                {
                    patchFields.set
                    (
                        completePatchi,
                        fvPatchField<Type>::New
                        (
                            procField.boundaryField()[procPatchi],
                            completeMesh_.boundary()[completePatchi],
                            DimensionedField<Type, volMesh>::null(),
                            fvPatchFieldReconstructor
                            (
                                completeMesh_.boundary()
                                [
                                    completePatchi
                                ].size()
                            )
                        )
                    );
                }

                patchFields[completePatchi].rmap
                (
                    procField.boundaryField()[procPatchi],
                    faceProcAddressingBf_[proci][procPatchi] - 1
                );
            }
            else if (isA<processorCyclicFvPatch>(procPatch))
            {
                if (!patchFields(completePatchi))
                {
                    patchFields.set
                    (
                        completePatchi,
                        fvPatchField<Type>::New
                        (
                            completeMesh_.boundary()[completePatchi].type(),
                            completeMesh_.boundary()[completePatchi],
                            DimensionedField<Type, volMesh>::null()
                        )
                    );
                }

                patchFields[completePatchi].rmap
                (
                    procField.boundaryField()[procPatchi],
                    faceProcAddressingBf_[proci][procPatchi] - 1
                );
            }
            else if (isA<indirectPolyPatch>(procPatch.patch()))
            {
                // Determine the index of the complete indirect patch
                const label completeIndID =
                    completeMesh_.boundaryMesh().findPatchID
                    (
                        procPatch.patch().name()
                    );

                if (!patchFields(completeIndID))
                {
                    patchFields.set
                    (
                        completeIndID,
                        fvPatchField<Type>::New
                        (
                            procField.boundaryField()[procPatchi],
                            completeMesh_.boundary()[completeIndID],
                            DimensionedField<Type, volMesh>::null(),
                            fvPatchFieldReconstructor
                            (
                                completeMesh_.boundary()[completeIndID].size()
                            )
                        )
                    );
                }

                const indirectPolyPatch& procIndPP =
                    refCast<const indirectPolyPatch>
                    (
                            procField.mesh().boundary()[procPatchi].patch()
                    );

                const indirectPolyPatch& recIndPP =
                    refCast<const indirectPolyPatch>
                    (
                        completeMesh_.boundary()[completeIndID].patch()
                    );

                const labelList& pAddr = procIndPP.fAddr();
                const labelList& rAddr = recIndPP.fAddr();

                labelList reverseAddressing(procIndPP.size(), -1);
                forAll(reverseAddressing, fI)
                {
                    const label pfI = pAddr[fI];
                    const label recfI =
                        mag(faceProcAddressing_[proci][pfI]) - 1;

                    bool found = false;
                    forAll(rAddr, gfI)
                    {
                        if (rAddr[gfI] == recfI)
                        {
                            reverseAddressing[fI] = gfI;
                            found = true;
                        }
                        if (found) break;
                    }
                    if (!found) Pout<< "Mapping is not found" << endl;
                }

                patchFields[completeIndID].rmap
                (
                    procField.boundaryField()[procPatchi],
                    reverseAddressing
                );
            }
        }
    }

    // Construct and return the field
    tmp<GeometricField<Type, fvPatchField, volMesh>> tfield
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            fieldIoObject,
            completeMesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );

    tfield.ref().oriented() = procFields[0].oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::fvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject
) const
{
    PtrList<GeometricField<Type, fvPatchField, volMesh>>
        procFields(procMeshes_.size());

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeometricField<Type, fvPatchField, volMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().timeName(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMeshes_[proci]
            )
        );
    }

    return reconstructFvVolumeField
    (
        IOobject
        (
            fieldIoObject.name(),
            completeMesh_.time().timeName(),
            completeMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        procFields
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject,
    const PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>& procFields
) const
{
    // Create the internalField
    Field<Type> internalField(completeMesh_.nInternalFaces());

    // Create the patch fields
    PtrList<fvsPatchField<Type>> patchFields(completeMesh_.boundary().size());

    forAll(procMeshes_, proci)
    {
        const GeometricField<Type, fvsPatchField, surfaceMesh>& procField =
            procFields[proci];

        // Set the internal face values in the reconstructed field
        rmapFaceToFace
        (
            internalField,
            procField.primitiveField(),
            SubList<label>
            (
                faceProcAddressing_[proci],
                procMeshes_[proci].nInternalFaces()
            ),
            isFlux(procFields[proci])
        );

        // Set the boundary patch values in the reconstructed field
        forAll(procField.boundaryField(), procPatchi)
        {
            const fvPatch& procPatch =
                procMeshes_[proci].boundary()[procPatchi];

            // Determine the index of the corresponding complete patch
            const label completePatchi = completePatchID(proci, procPatchi);

            // Check if the boundary patch is not a processor patch
            // nor an indirect patch
            if (completePatchi == procPatchi)
            {
                if (!patchFields(completePatchi))
                {
                    patchFields.set
                    (
                        completePatchi,
                        fvsPatchField<Type>::New
                        (
                            procField.boundaryField()[procPatchi],
                            completeMesh_.boundary()[completePatchi],
                            DimensionedField<Type, surfaceMesh>::null(),
                            fvPatchFieldReconstructor
                            (
                                completeMesh_.boundary()[completePatchi].size()
                            )
                        )
                    );
                }

                patchFields[completePatchi].rmap
                (
                    procField.boundaryField()[procPatchi],
                    faceProcAddressingBf_[proci][procPatchi] - 1
                );
            }
            else if (isA<processorCyclicFvPatch>(procPatch))
            {
                if (!patchFields(completePatchi))
                {
                    patchFields.set
                    (
                        completePatchi,
                        fvsPatchField<Type>::New
                        (
                            completeMesh_.boundary()[completePatchi].type(),
                            completeMesh_.boundary()[completePatchi],
                            DimensionedField<Type, surfaceMesh>::null()
                        )
                    );
                }

                patchFields[completePatchi].rmap
                (
                    procField.boundaryField()[procPatchi],
                    faceProcAddressingBf_[proci][procPatchi] - 1
                );
            }
            else if (isA<processorFvPatch>(procPatch))
            {
                rmapFaceToFace
                (
                    internalField,
                    procField.boundaryField()[procPatchi],
                    faceProcAddressingBf_[proci][procPatchi],
                    isFlux(procFields[proci])
                );
            }
            else if (isA<indirectPolyPatch>(procPatch.patch()))
            {
                // Determine the index of the complete indirect patch
                const label completeIndID =
                    completeMesh_.boundaryMesh().findPatchID
                    (
                        procPatch.patch().name()
                    );

                if (!patchFields(completeIndID))
                {
                    patchFields.set
                    (
                        completeIndID,
                        fvsPatchField<Type>::New
                        (
                            procField.boundaryField()[procPatchi],
                            completeMesh_.boundary()[completeIndID],
                            DimensionedField<Type, surfaceMesh>::null(),
                            fvPatchFieldReconstructor
                            (
                                completeMesh_.boundary()[completeIndID].size()
                            )
                        )
                    );
                }

                const indirectPolyPatch& procIndPP =
                    refCast<const indirectPolyPatch>
                    (
                        procField.mesh().boundary()[procPatchi].patch()
                    );

                const indirectPolyPatch& recIndPP =
                    refCast<const indirectPolyPatch>
                    (
                        completeMesh_.boundary()[completeIndID].patch()
                    );

                const labelList& pAddr = procIndPP.fAddr();
                const labelList& rAddr = recIndPP.fAddr();

                labelList reverseAddressing(procIndPP.size(), -1);
                forAll(reverseAddressing, fI)
                {
                    const label pfI = pAddr[fI];
                    const label recfI =
                        mag(faceProcAddressing_[proci][pfI]) - 1;

                    bool found = false;
                    forAll(rAddr, gfI)
                    {
                        if (rAddr[gfI] == recfI)
                        {
                            reverseAddressing[fI] = gfI;
                            found = true;
                        }
                        if (found) break;
                    }
                    if (!found) Pout<< "Mapping is not found" << endl;
                }

                patchFields[completeIndID].rmap
                (
                    procField.boundaryField()[procPatchi],
                    reverseAddressing
                );
            }
        }
    }

    // Construct and return the field
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfield
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            fieldIoObject,
            completeMesh_,
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );

    tfield.ref().oriented() = procFields[0].oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject
) const
{
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh>>
        procFields(procMeshes_.size());

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeometricField<Type, fvsPatchField, surfaceMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci].time().timeName(),
                    procMeshes_[proci],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                procMeshes_[proci]
            )
        );
    }

    return reconstructFvSurfaceField
    (
        IOobject
        (
            fieldIoObject.name(),
            completeMesh_.time().timeName(),
            completeMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        procFields
    );
}


template<class Type>
void Foam::fvFieldReconstructor::reconstructFvVolumeInternalFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    const word& fieldClassName = DimensionedField<Type, volMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructFvVolumeInternalField<Type>(*fieldIter())().write();

                nReconstructed_++;
            }
        }
        Info<< endl;
    }
}


template<class Type>
void Foam::fvFieldReconstructor::reconstructFvVolumeFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    const word& fieldClassName =
        GeometricField<Type, fvPatchField, volMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructFvVolumeField<Type>(*fieldIter())().write();

                nReconstructed_++;
            }
        }
        Info<< endl;
    }
}


template<class Type>
void Foam::fvFieldReconstructor::reconstructFvSurfaceFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    const word& fieldClassName =
        GeometricField<Type, fvsPatchField, surfaceMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructFvSurfaceField<Type>(*fieldIter())().write();

                nReconstructed_++;
            }
        }
        Info<< endl;
    }
}


// ************************************************************************* //
