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
    (c) 2016-2021 Esi Ltd

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "cfdTools/general/pointBoundarySetup/pointBoundarySetup.H"
#include "cfdTools/general/surfaceBoundarySetup/surfaceBoundarySetup.H"
#include "surfaceMesh/surfaceMesh.H"
#include "volMesh/volMesh.H"

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class Type>
void Foam::stateFunctions::regionState::storeInternalField
(
    const word& name,
    const Field<Type>&
)
{
    //do nothing - should not reach here due to specialisations
    FatalErrorInFunction
        << "Illegal call to unspecialised storage function."
        << exit(FatalError);
}

template<class Type>
bool Foam::stateFunctions::regionState::setInternalValue(const word& name)
{
    //do nothing - should not reach here due to specialisations
    FatalErrorInFunction
        << "Illegal call to unspecialised assignment function."
        << exit(FatalError);

    return false;
}

template<class Type>
void Foam::stateFunctions::regionState::createVolField
(
    const fvMesh& mesh,
    const dictionary& fieldDict,
    const word fieldName
)
{
    IOdictionary cbcDict
    (
        IOobject
        (
            fieldName,
            mesh.time().caseSystem(),
            "boundaryConditions",
            localDb().registry(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        boundarySetup<Type>
        (
            mesh,
            fieldName,
            fieldDict,
            switches().strictPatchNameChecking(),
            switches().printBoundaries()
        ).boundaryDict()
    );

    dictionary fieldDef;

    fieldDef.add("dimensions", fieldDict.lookup("dimensions"));
    fieldDef.add("internalField", fieldDict.lookup("internalField"));

    fieldDef.set("boundaryField", dictionary());
    fieldDef.subDict("boundaryField").merge(cbcDict);

    if
    (
        switches().resetBoundaryFields()
      &&(collated_ || distributed_ || (Pstream::master() || !Pstream::parRun()))
    )
    {

        //if resetInternalFields = false, then coupled patch boundaries from
        // existing fields will be kept
        //this can cause issues the master processor is doing the CBC writing
        //(see regionState.C line 776)

        //delete coupled patch entries from CBC

        if (!switches().resetInternalFields())
        {
            forAll(meshPtr_->boundary(), pi)
            {
                if (meshPtr_->boundary()[pi].coupled())
                {
                    const fvPatch& cPatch(meshPtr_->boundary()[pi]);

                    cbcDict.remove(cPatch.name());
                }
            }
        }

        forAll(meshPtr_->boundary(), pi)
        {
            const fvPatch& cPatch(meshPtr_->boundary()[pi]);
            if (isA<processorFvPatch>(cPatch))
            {
                cbcDict.remove(cPatch.name());
            }
        }

        cbcDict.regIOobject::writeObject
        (
            mesh.time().writeFormat(),
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true
        );
    }


    mesh.objectRegistry::store
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                localDb().registry(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            fieldDef,
            true
        )
    );
}

template<class Type>
void Foam::stateFunctions::regionState::createSurfaceField
(
    const surfaceMesh& cmesh,
    const dictionary& fieldDict,
    const word fieldName
)
{
    //initial pointer to field
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> surfaceTypeField;

    autoPtr<surfaceTypeField> fieldPtr(nullptr);

    fieldPtr.reset
    (
        new surfaceTypeField
        (
            IOobject
            (
                fieldName,
                cmesh().time().timeName(),
                localDb().registry(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            cmesh(),
            dimensioned<Type>
            (
                fieldName,
                dimensionSet(fieldDict.lookup("dimensions")),
                pTraits<Type>::zero
            )
        )
    );

    fieldPtr->primitiveFieldRef()
        = Field<Type>("internalField", fieldDict, mesh().nFaces());

    fieldPtr->boundaryFieldRef().reset
    (
        typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary
        (
            cmesh().boundary(),
            fieldPtr->internalField(),
            surfaceBoundarySetup<Type>
            (
                cmesh,
                fieldPtr->internalField(),
                fieldDict,
                switches().strictPatchNameChecking()
            )
        )
    );

    fieldPtr->store(fieldPtr);
}

template<class Type>
void Foam::stateFunctions::regionState::createPointField
(
    const pointMesh& mesh,
    const dictionary& fieldDict,
    const word fieldName
)
{
    IOdictionary cbcDict
    (
        IOobject
        (
            fieldName,
            mesh.time().caseSystem(),
            "boundaryConditions",
            localDb().registry(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointBoundarySetup<Type>
        (
            mesh,
            fieldName,
            fieldDict,
            switches().strictPatchNameChecking()
        ).boundaryDict()
    );

    dictionary fieldDef;

    fieldDef.add("dimensions", fieldDict.lookup("dimensions"));
    fieldDef.add("internalField", fieldDict.lookup("internalField"));

    fieldDef.set("boundaryField", dictionary());
    fieldDef.subDict("boundaryField").merge(cbcDict);

    if
    (
        switches().resetBoundaryFields()
      &&(collated_ || distributed_ || (Pstream::master() || !Pstream::parRun()))
    )
    {

        //if resetInternalFields = false, then coupled patch boundaries from
        // existing fields will be kept
        //this can cause issues the master processor is doing the CBC writing
        //(see regionState.C line 776)

        //delete coupled patch entries from CBC

        if (!switches().resetInternalFields())
        {
            forAll(mesh.boundary(), pi)
            {
                if (mesh.boundary()[pi].coupled())
                {
                    const pointPatch& cPatch(mesh.boundary()[pi]);

                    cbcDict.remove(cPatch.name());
                }
            }
        }

        forAll(mesh.boundary(), pi)
        {
            const pointPatch& cPatch(mesh.boundary()[pi]);
            if (isA<processorPointPatch>(cPatch))
            {
                cbcDict.remove(cPatch.name());
            }
        }

        cbcDict.regIOobject::writeObject
        (
            mesh.time().writeFormat(),
            IOstream::currentVersion,
            IOstream::UNCOMPRESSED,
            true
        );
    }


    //mesh.objectRegistry::store
    mesh.thisDb().store
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                fieldName,
                mesh.time().timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            fieldDef,
            true
        )
    );
}


template<class Type>
void Foam::stateFunctions::regionState::createTypeField
(
    const dictionary& fd,
    const word fieldName
)
{
    meshType met
    (
        stateFunction::meshTypeNames_
        [fd.lookupOrDefault<word>("meshType", "volume")]
    );

    switch (met)
    {
        case meVolume :
        {
            createVolField<Type>
            (
                *meshPtr_, fd, fieldName
            );
        } break;

        case meSurface :
        {
            const surfaceMesh sMesh(mesh());
            createSurfaceField<Type>
            (
                sMesh, fd, fieldName
            );
        } break;

        case mePoint :
        {
            const pointMesh& pMesh = pointMesh::New(mesh());
            createPointField<Type>
            (
                pMesh, fd, fieldName
            );
        } break;


        default:
        {
            FatalErrorInFunction
                << "Invalid mesh type '"
                << fd.lookupOrDefault<word>("meshType", "volume")
                << "'" << exit(FatalError);

        }
    }
}

template<class Type, template<class> class PatchField, class Mesh, class Mesh2>
void Foam::stateFunctions::regionState::readInitialGeoField
(
    const objectRegistry& obr,
    const Mesh2& cmesh,
    const dictionary& fd,
    const word& fieldName
)
{
    //initial pointer to field
    autoPtr<GeometricField<Type, PatchField, Mesh>> fieldPtr(nullptr);

    //check if field exists
    IOobject fheader
    (
        fieldName,
        obr.time().timeName(),
        obr,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    if (fheader.typeHeaderOk<GeometricField<Type, PatchField, Mesh>>(true))
    {
        //- if reset boundary fields, transform them to calculated.
        //  By doing so, we will avoid comflicts with constraint types and
        //  change of boundary types where conflict with WF
        if (!switches().resetBoundaryFields())
        {
            fieldPtr.reset
            (
                new GeometricField<Type, PatchField, Mesh>
                (
                    fheader,
                    cmesh
                )
            );
        }
        else
        {
            //- it is dimless but it is overidden later
            fheader.readOpt() = IOobject::NO_READ;
            fieldPtr.reset
            (
                new GeometricField<Type, PatchField, Mesh>
                (
                    fheader,
                    cmesh,
                    dimless,
                    PatchField<Type>::calculatedType()
                )
            );
            fieldPtr->readInternalFields();
        }

        //store internal field
        storeInternalField<Type>(fieldName, fieldPtr->primitiveField());

        // store boundary field as dictionary
        OStringStream outbuf;
        fieldPtr->boundaryFieldRef().writeEntry("boundaryField", outbuf);
        IStringStream inbuf(outbuf.str());

        initialFieldBoundaryDictsPtr_->add
        (
            fieldName,
            dictionary(inbuf).subDict("boundaryField")
        );

        Info<< "   " << fieldName << ": stored" << endl;
    }
    else //no initial field found
    {
        Info<< tab << fieldName << ": not found" << endl;
    }

}


template<class Type>
void Foam::stateFunctions::regionState::readInitialTypeField
(
    const dictionary& fd,
    const word fieldName
)
{
    meshType met
    (
        stateFunction::meshTypeNames_
        [fd.lookupOrDefault<word>("meshType", "volume")]
    );

    const objectRegistry& obr(mesh());

    switch (met)
    {
        case meVolume :
        {
            readInitialGeoField<Type, fvPatchField, volMesh, fvMesh>
            (
                obr, mesh(), fd, fieldName
            );
        } break;

        case meSurface :
        {
            readInitialGeoField<Type, fvsPatchField, surfaceMesh, fvMesh>
            (
                obr, mesh(), fd, fieldName
            );
        } break;

        case mePoint :
        {
            const pointMesh& pMesh = pointMesh::New(mesh());
            readInitialGeoField<Type, pointPatchField, pointMesh, pointMesh>
            (
                obr, pMesh, fd, fieldName
            );
        } break;


        default:
        {
            FatalErrorInFunction
                << "Invalid mesh type '"
                << fd.lookupOrDefault<word>("meshType", "volume")
                << "'" << exit(FatalError);

        }
    }
}


template<template<class> class PatchType, class MeshType>
void Foam::stateFunctions::regionState::reInitBoundaryTypes(const fieldInit& fIS)
{
    reInitBoundaryType<scalar, PatchType, MeshType>(fIS);
    reInitBoundaryType<vector, PatchType, MeshType>(fIS);
    reInitBoundaryType<tensor, PatchType, MeshType>(fIS);
    reInitBoundaryType<symmTensor, PatchType, MeshType>(fIS);
    reInitBoundaryType<sphericalTensor, PatchType, MeshType>(fIS);
}


template<class Type, template<class> class PatchType, class MeshType>
void Foam::stateFunctions::regionState::reInitBoundaryType(const fieldInit& fIS)
{
    typedef GeometricField<Type, PatchType, MeshType> GeoField;
    if (localDb().registry().foundObject<GeoField>(fIS.name()))
    {
        GeoField& f(localDb().registry().lookupObjectRef<GeoField>(fIS.name()));

        OStringStream oss(IOstream::ASCII);
        f.boundaryFieldRef().writeEntries(oss);
        IStringStream iss(oss.str(), IOstream::ASCII);
        dictionary bdryDict(iss);
        f.boundaryFieldRef().readField(f.internalField(), bdryDict);
    }
}


// ************************************************************************* //
