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
    (C) 2019 OpenCFD Ltd.
\*---------------------------------------------------------------------------*/

#include "polySurface/polySurface.H"
#include "polySurface/fields/polySurfaceFields.H"
#include "polySurface/fields/polySurfacePointFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoMeshType>
const Foam::regIOobject* Foam::polySurface::findFieldObject
(
    const word& fieldName
) const
{
    // Face Data first (main registry)

    const objectRegistry& obr = *this;

    const auto* ioptr = obr.cfindObject<regIOobject>(fieldName);

    if (ioptr)
    {
        return ioptr;
    }

    for (const regIOobject* ioObj : obr)
    {
        if (isA<objectRegistry>(ioObj))
        {
            const objectRegistry* subreg =
                static_cast<const objectRegistry*>(ioObj);

            if ((ioptr = subreg->cfindObject<regIOobject>(fieldName)))
            {
                return ioptr;
            }
        }
    }

    return ioptr;
}


template<class GeoMeshType>
const Foam::objectRegistry* Foam::polySurface::whichRegistry
(
    const word& fieldName
) const
{
    // Face Data first (main registry)

    const objectRegistry& obr = *this;

    if (obr.found(fieldName))
    {
        return this;
    }


    for (const regIOobject* ioObj : obr)
    {
        if (isA<objectRegistry>(ioObj))
        {
            const objectRegistry* subreg =
                static_cast<const objectRegistry*>(ioObj);

            if (subreg->found(fieldName))
            {
                return subreg;
            }
        }
    }

    return nullptr;
}


template<class Type, class GeoMeshType>
void Foam::polySurface::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    const Field<Type>& values
)
{
    // Force creates field database if needed.
    const objectRegistry& fieldDb = this->fieldData<GeoMeshType>();

    auto* dimfield =
        fieldDb.getObjectPtr<DimensionedField<Type, GeoMeshType>>(fieldName);

    if (dimfield)
    {
        dimfield->dimensions() = dims;
        dimfield->field() = values;
    }
    else
    {
        dimfield = new DimensionedField<Type, GeoMeshType>
        (
            IOobject
            (
                fieldName,
                fieldDb,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dims,
            values
        );

        dimfield->store();
    }
}


template<class Type, class GeoMeshType>
void Foam::polySurface::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    Field<Type>&& values
)
{
    // Force creates field database if needed.
    const objectRegistry& fieldDb = this->fieldData<GeoMeshType>();

    auto* dimfield =
        fieldDb.getObjectPtr<DimensionedField<Type, GeoMeshType>>(fieldName);

    if (dimfield)
    {
        dimfield->dimensions() = dims;
        dimfield->field() = std::move(values);
    }
    else
    {
        dimfield = new DimensionedField<Type, GeoMeshType>
        (
            IOobject
            (
                fieldName,
                fieldDb,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            dims,
            std::move(values)
        );

        dimfield->store();
    }
}


// ************************************************************************* //
