/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/boundarySetup/boundarySetup.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/directFvPatchFieldMapper.H"

// * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * //

template<class T>
void Foam::addGIBBC
(
    const dictionary& bf,
    const word fieldName,
    const fvMesh& initMesh,
    const fvMesh& updatedMesh
)
{
    typedef GeometricField<T, fvPatchField, volMesh> GeometricFieldType;

    bool found = initMesh.foundObject<GeometricFieldType>(fieldName);
    if (found)
    {
        GeometricFieldType& field =
            const_cast<GeometricFieldType&>
            (
                initMesh.lookupObject<GeometricFieldType>(fieldName)
            );

        PtrList<fvPatchField<T>> nubf(updatedMesh.boundary().size());

        //-- create existing new fvPatchFields
        //   Use this ::New function because the initMesh has changed and patch_
        //   reference inside the BCs is not correct
        forAll(field.boundaryField(), pI)
        {
            labelList map(identity(updatedMesh.boundary()[pI].size()));
            nubf.set
            (
                pI,
                fvPatchField<T>::New
                (
                    field.boundaryField()[pI],
                    initMesh.boundary()[pI],
                    field.internalField(),
                    directFvPatchFieldMapper(map)
                )
            );

        }

        //-- create GIB BCs
        int nI = field.boundaryField().size();
        forAllConstIter(IDLList<entry>, bf, iter)
        {
            nubf.set
            (
                nI
                ,
                    fvPatchField<T>::New
                    (
                        updatedMesh.boundary()[nI],
                        field.internalField(),
                        iter().dict()
                    )
            );
            nI++;
        }


        //-- create GIB BCs
        typename GeometricField<T, fvPatchField, volMesh>::
            Boundary newBCs
            (
                updatedMesh.boundary(),
                field.internalField(),
                nubf
            );

        GeometricFieldType fieldNewMesh(
            IOobject
            (
                field.name(),
                updatedMesh.time().timeName(),
                updatedMesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            updatedMesh,
            field.dimensions(),
            field.internalField(),
            newBCs
        );
        fieldNewMesh.write();
    }
    else
    {
    }
}

// ************************************************************************* //
