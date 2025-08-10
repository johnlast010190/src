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
    (c) 2015-2017 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "readFields/readFields.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "surfFields/surfFields/surfFields.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::readFields::loadField(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, surfGeoMesh> SurfFieldType;

    if (foundObject<VolFieldType>(fieldName))
    {
        DebugInformation
            << "readFields : " << VolFieldType::typeName
            << " " << fieldName << " already in database"
            << endl;
    }
    else if (foundObject<SurfaceFieldType>(fieldName))
    {
        DebugInformation<< "readFields: " << SurfaceFieldType::typeName
            << " " << fieldName << " already exists in database"
            << " already in database" << endl;
    }
    else if (foundObject<SurfFieldType>(fieldName))
    {
        DebugInformation<< "readFields: " << SurfFieldType::typeName
            << " " << fieldName << " already exists in database"
            << " already in database" << endl;
    }
    else
    {
        IOobject fieldHeader
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (fieldHeader.typeHeaderOk<VolFieldType>(true, true, false))
        {
            // Store field on mesh database
            Log << "    Reading " << fieldName << endl;
            VolFieldType* vfPtr(new VolFieldType(fieldHeader, mesh_));
            mesh_.objectRegistry::store(vfPtr);
            return true;
        }
        else if (fieldHeader.typeHeaderOk<SurfaceFieldType>(true, true, false))
        {
            // Store field on mesh database
            Log << "    Reading " << fieldName << endl;
            SurfaceFieldType* sfPtr(new SurfaceFieldType(fieldHeader, mesh_));
            mesh_.objectRegistry::store(sfPtr);
            return true;
        }
        else if (fieldHeader.typeHeaderOk<SurfFieldType>(true))
        {
            if (isA<surfMesh>(obr()))
            {
                const surfMesh& s = dynamicCast<const surfMesh>(obr());

                // Store field on surfMesh database
                Log << "    Reading " << fieldName << endl;
                SurfFieldType* sfPtr(new SurfFieldType(fieldHeader, s));
                s.store(sfPtr);
                return true;
            }
            else
            {
                return false;
            }
        }
    }

    return false;
}


// ************************************************************************* //
