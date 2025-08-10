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

Application
    setFields

Group
    grpPreProcessingUtilities

Description
    Set values on a selected set of cells/patchfaces through a dictionary.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "sets/topoSetSource/topoSetSource.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "fields/volFields/volFields.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "containers/Lists/CompactListList/CompactListList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

template<class Type>
bool setCellFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedCells,
    Istream& fieldValueStream
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);

    // Check the current time directory
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.typeHeaderOk<fieldType>(true))
    {
        fieldHeader = IOobject
        (
            fieldName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        );
    }

    // Check field exists
    if (fieldHeader.typeHeaderOk<fieldType>(true))
    {
        Info<< "    Setting internal values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        fieldType field(fieldHeader, mesh);

        const Type value = pTraits<Type>(fieldValueStream);

        if (selectedCells.size() == field.size())
        {
            field.primitiveFieldRef() = value;
        }
        else
        {
            forAll(selectedCells, celli)
            {
                field[selectedCells[celli]] = value;
            }
        }

        typename GeometricField<Type, fvPatchField, volMesh>::
            Boundary& fieldBf = field.boundaryFieldRef();

        forAll(field.boundaryField(), patchi)
        {
            fieldBf[patchi] = fieldBf[patchi].patchInternalField();
        }

        if (!field.write())
        {
            FatalErrorInFunction
              << "Failed writing field " << fieldName << endl;
        }
    }
    else
    {
        WarningInFunction
          << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}


class setCellField
{

public:

    setCellField()
    {}

    autoPtr<setCellField> clone() const
    {
        return autoPtr<setCellField>(new setCellField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedCells_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedCells)
        :
            mesh_(mesh),
            selectedCells_(selectedCells)
        {}

        autoPtr<setCellField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setCellFieldType<scalar>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<vector>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<symmTensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                 || setCellFieldType<tensor>
                        (fieldType, mesh_, selectedCells_, fieldValues)
                )
            )
            {
                WarningInFunction
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setCellField>(new setCellField());
        }
    };
};


template<class Type>
bool setFaceFieldType
(
    const word& fieldTypeDesc,
    const fvMesh& mesh,
    const labelList& selectedFaces,
    Istream& fieldValueStream
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (fieldTypeDesc != fieldType::typeName + "Value")
    {
        return false;
    }

    word fieldName(fieldValueStream);

    // Check the current time directory
    IOobject fieldHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ
    );

    // Check the "constant" directory
    if (!fieldHeader.typeHeaderOk<fieldType>(true))
    {
        fieldHeader = IOobject
        (
            fieldName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ
        );
    }

    // Check field exists
    if (fieldHeader.typeHeaderOk<fieldType>(true))
    {
        Info<< "    Setting patchField values of "
            << fieldHeader.headerClassName()
            << " " << fieldName << endl;

        // Read the field
        fieldType field(fieldHeader, mesh);
        typename fieldType::Boundary& fieldBf = field.boundaryFieldRef();

        // Read the value
        const Type value = pTraits<Type>(fieldValueStream);

        // Determine the number of non-processor patches
        label nNonProcPatches = 0;
        forAll(fieldBf, patchi)
        {
            nNonProcPatches = patchi;
            if (isA<processorFvPatch>(mesh.boundary()[patchi]))
            {
                break;
            }
        }

        // Create a copy of the boundary field
        typename fieldType::Boundary fieldBfCopy
        (
            fieldType::Internal::null(),
            fieldBf
        );

        // Loop selected faces and set values in the copied boundary field
        bool haveWarnedInternal = false, haveWarnedProc = false;
        labelList nonProcPatchNChangedFaces(nNonProcPatches, 0);
        forAll(selectedFaces, i)
        {
            const label facei = selectedFaces[i];

            if (mesh.isInternalFace(facei))
            {
                if (!haveWarnedInternal)
                {
                    WarningInFunction
                        << "Ignoring internal face " << facei
                        << ". Suppressing further warnings." << endl;
                    haveWarnedInternal = true;
                }
            }
            else
            {
                const labelUList patches =
                    mesh.polyBFacePatches()[facei - mesh.nInternalFaces()];
                const labelUList patchFaces =
                    mesh.polyBFacePatchFaces()[facei - mesh.nInternalFaces()];

                forAll(patches, i)
                {
                    if (patches[i] >= nNonProcPatches)
                    {
                        if (!haveWarnedProc)
                        {
                            WarningInFunction
                                << "Ignoring face " << patchFaces[i]
                                << " of processor patch " << patches[i]
                                << ". Suppressing further warnings." << endl;
                            haveWarnedProc = true;
                        }
                    }
                    else
                    {
                        fieldBfCopy[patches[i]][patchFaces[i]] = value;
                        nonProcPatchNChangedFaces[patches[i]] ++;
                    }
                }
            }
        }

        Pstream::listCombineGather
        (
            nonProcPatchNChangedFaces,
            plusEqOp<label>()
        );
        Pstream::listCombineScatter
        (
            nonProcPatchNChangedFaces
        );

        // Reassign boundary values
        forAll(nonProcPatchNChangedFaces, patchi)
        {
            if (nonProcPatchNChangedFaces[patchi] > 0)
            {
                Info<< "    On patch "
                    << field.boundaryField()[patchi].patch().name()
                    << " set " << nonProcPatchNChangedFaces[patchi]
                    << " values" << endl;
                fieldBf[patchi].forceAssign(fieldBfCopy[patchi]);
            }
        }

        if (!field.write())
        {
            FatalErrorInFunction
                << "Failed writing field " << field.name() << exit(FatalError);
        }
    }
    else
    {
        WarningInFunction
          << "Field " << fieldName << " not found" << endl;

        // Consume value
        (void)pTraits<Type>(fieldValueStream);
    }

    return true;
}


class setFaceField
{

public:

    setFaceField()
    {}

    autoPtr<setFaceField> clone() const
    {
        return autoPtr<setFaceField>(new setFaceField());
    }

    class iNew
    {
        const fvMesh& mesh_;
        const labelList& selectedFaces_;

    public:

        iNew(const fvMesh& mesh, const labelList& selectedFaces)
        :
            mesh_(mesh),
            selectedFaces_(selectedFaces)
        {}

        autoPtr<setFaceField> operator()(Istream& fieldValues) const
        {
            word fieldType(fieldValues);

            if
            (
               !(
                    setFaceFieldType<scalar>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<vector>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<sphericalTensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<symmTensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                 || setFaceFieldType<tensor>
                        (fieldType, mesh_, selectedFaces_, fieldValues)
                )
            )
            {
                WarningInFunction
                    << "field type " << fieldType << " not currently supported"
                    << endl;
            }

            return autoPtr<setFaceField>(new setFaceField());
        }
    };
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "include/addDictOption.H"
    #include "include/addRegionOption.H"
    #include "include/setRootCase.H"
    #include "include/createTime.H"
    #include "include/createNamedMesh.H"

    const word dictName("setFieldsDict");
    #include "include/setSystemMeshDictionaryIO.H"

    Info<< "Reading " << dictName << "\n" << endl;

    IOdictionary setFieldsDict(dictIO);

    if (setFieldsDict.found("defaultFieldValues"))
    {
        Info<< "Setting field default values" << endl;
        PtrList<setCellField> defaultFieldValues
        (
            setFieldsDict.lookup("defaultFieldValues"),
            setCellField::iNew(mesh, labelList(mesh.nCells()))
        );
        Info<< endl;
    }

    Info<< "Setting field region values" << endl;

    PtrList<entry> regions(setFieldsDict.lookup("regions"));

    forAll(regions, regionI)
    {
        const entry& region = regions[regionI];

        autoPtr<topoSetSource> source =
            topoSetSource::New(region.keyword(), mesh, region.dict());

        if (source().setType() == topoSetSource::CELLSETSOURCE)
        {
            cellSet selectedCellSet
            (
                mesh,
                "cellSet",
                mesh.nCells()/10+1  // Reasonable size estimate.
            );

            source->applyToSet
            (
                topoSetSource::NEW,
                selectedCellSet
            );

            PtrList<setCellField> fieldValues
            (
                region.dict().lookup("fieldValues"),
                setCellField::iNew(mesh, selectedCellSet.toc())
            );
        }
        else if (source().setType() == topoSetSource::FACESETSOURCE)
        {
            faceSet selectedFaceSet
            (
                mesh,
                "faceSet",
                (mesh.nFaces()-mesh.nInternalFaces())/10+1
            );

            source->applyToSet
            (
                topoSetSource::NEW,
                selectedFaceSet
            );

            PtrList<setFaceField> fieldValues
            (
                region.dict().lookup("fieldValues"),
                setFaceField::iNew(mesh, selectedFaceSet.toc())
            );
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
