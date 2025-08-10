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
    (c) 2017 OpenCFD Ltd.

Class
    Foam::MapInternalField

Description
    Generic internal field mapper.  For "real" mapping, add template
    specialisations for mapping of internal fields depending on mesh
    type.

\*---------------------------------------------------------------------------*/

#ifndef MapGeometricFields_H
#define MapGeometricFields_H

#include "meshes/polyMesh/polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class MeshMapper, class GeoMesh>
class MapInternalField
{
public:

    MapInternalField()
    {}

    void operator()
    (
        DimensionedField<Type, GeoMesh>& field,
        const MeshMapper& mapper
    ) const;
};


//- Generic Geometric field mapper.
//  For "real" mapping, add template specialisations
//  for mapping of internal fields depending on mesh type.
template
<
    class Type,
    template<class> class PatchField,
    class MeshMapper,
    class GeoMesh
>
void MapGeometricFields
(
    const MeshMapper& mapper
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> FieldType;

    HashTable<const FieldType*> fields
    (
        mapper.thisDb().objectRegistry::template
            lookupClass<FieldType>()
    );

    // It is necessary to enforce that all old-time fields are stored
    // before the mapping is performed.  Otherwise, if the
    // old-time-level field is mapped before the field itself, sizes
    // will not match.

    forAllConstIter(typename HashTable<const FieldType*>, fields, fieldIter)
    {
        FieldType& field = const_cast<FieldType&>(*fieldIter());

        //Note: check can be removed once pointFields are actually stored on
        //      the pointMesh instead of now on the polyMesh!
        if (&field.mesh() == &mapper.mesh())
        {
            field.storeOldTimes();
        }
    }

    forAllConstIter(typename HashTable<const FieldType*>, fields, fieldIter)
    {
        FieldType& field = const_cast<FieldType&>(*fieldIter());

        if (&field.mesh() == &mapper.mesh())
        {
            if (polyMesh::debug)
            {
                Info<< "Mapping " << field.typeName << ' ' << field.name()
                    << endl;
            }

            // Map the internal field
            MapInternalField<Type, MeshMapper, GeoMesh>()(field.ref(), mapper);

            // Map the patch fields
            typename FieldType::Boundary& bfield = field.boundaryFieldRef();

            forAll(bfield, patchi)
            {
                // Cannot check sizes for patch fields because of
                // empty fields in FV and because point fields get their size
                // from the patch which has already been resized
                //

                bfield[patchi].autoMap(mapper.boundaryMap()[patchi]);
            }

            field.instance() = field.time().timeName();
        }
        else if (polyMesh::debug)
        {
            Info<< "Not mapping " << field.typeName << ' ' << field.name()
                << " since originating mesh differs from that of mapper."
                << endl;
        }
    }
}

template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void UpdateGIBFields
(
    polyMesh& mesh
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {

        typedef GeometricField<Type, PatchField, GeoMesh> FieldType;
        FieldType& field = const_cast<FieldType&>(*fieldIter());

        // Map the patch fields
        typename FieldType::Boundary& bfield = field.boundaryFieldRef();

        forAll(bfield, patchi)
        {
            // Cannot check sizes for patch fields because of
            // empty fields in FV and because point fields get their size
            // from the patch which has already been resized
            //
            bfield[patchi].updateGIB();
        }
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void StoreGIBFields
(
    polyMesh& mesh
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {

        typedef GeometricField<Type, PatchField, GeoMesh> FieldType;
        FieldType& field = const_cast<FieldType&>(*fieldIter());

        // Map the patch fields
        typename FieldType::Boundary& bfield = field.boundaryFieldRef();

        forAll(bfield, patchi)
        {
            bfield[patchi].storeGIB();
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
