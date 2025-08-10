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
    (c) 2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class GeoFieldType>
bool Foam::functionObjects::components::calcFieldComponents()
{
    typedef typename GeoFieldType::value_type Type;

    const GeoFieldType& field(lookupObject<GeoFieldType>(fieldName_));

    resultNames_.setSize(Type::nComponents);

    bool stored = true;

    for (direction i=0; i<Type::nComponents; i++)
    {
        resultName_ = fieldName_ + word(Type::componentNames[i]);
        resultNames_[i] = resultName_;

        stored = stored && store(resultName_, field.component(i));
    }

    return stored;
}


template<class Type>
bool Foam::functionObjects::components::calcComponents()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    if (foundObject<VolFieldType>(fieldName_, false))
    {
        return calcFieldComponents<VolFieldType>();
    }
    else if (foundObject<SurfaceFieldType>(fieldName_, false))
    {
        return calcFieldComponents<SurfaceFieldType>();
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
