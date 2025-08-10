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
    (c) 2015-2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/movingGIBTools/mapGIB/mapGIB.H"
#include "fields/fvPatchFields/fvPatchField/fieldMappers/gibFvPatchFieldMapper.H"
#include "fields/fvsPatchFields/fvsPatchField/fvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void mapGIB::StoreOldFieldsToPatch
(
    List<Tuple2<word, Field<Type>>>& ptrFields,
    const label& pI
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh_.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );
    ptrFields.resize(fields.size());
    label i = 0;
    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());
        if (field.name() != "meshPhi")
        {

            const List<Field<Type>>& bfield =
                field.boundaryField().oldTimeField();
            Type defaultValue = field.boundaryField()[pI].defaultGIBValue();

            ptrFields[i].first() = field.name();
            if (pI == mesh_.masterId())
            {
                ptrFields[i].second().setSize(sourcePatchm_->size());
            }
            else if (pI == mesh_.slaveId())
            {
                ptrFields[i].second().setSize(sourcePatchs_->size());
            }
            else
            {
                FatalErrorInFunction
                    << abort(FatalError);
            }

            forAll(ptrFields[i].second(), fI)
            {
                if (fI < bfield[pI].size())
                {
                    ptrFields[i].second()[fI] = bfield[pI][fI];
                }
                else
                {
                    ptrFields[i].second()[fI] = defaultValue;
                }
            }

            i++;
        }
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void mapGIB::MapGIBField
(
    List<Tuple2<word, Field<Type>>>& ptrFields,
    const label& pI,
    const AMIPatchToPatchInterpolation& ipToP
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh_.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

//    label i = 0;
    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        if (field.name() != "meshPhi")
        {
            typename GeometricField<Type, PatchField, GeoMesh>
            ::Boundary& bfield = field.boundaryFieldRef();

            Field<Type>& f = bfield[pI];

            forAll(ptrFields, j)
            {
                if (ptrFields[j].first() == field.name())
                {
                    const Field<Type>& oldFaceField = ptrFields[j].second();

                    tmp<Field<Type>> fnewt =
                        ipToP.interpolateToTarget(oldFaceField);
                    f = fnewt();

                    // Ask BC to map any stored fields
                    gibFvPatchFieldMapper mapper(ipToP);
                    bfield[pI].autoMapGIB(mapper);

//                    i++;
                }
            }
        }
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void mapGIB::Map1to1GIBField
(
    const List<Tuple2<word, Field<Type>>>& ptrmFields,
    const List<Tuple2<word, Field<Type>>>& ptrsFields
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh_.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

    const label& mId = mesh_.masterId();
    const label& sId = mesh_.slaveId();

//    label i = 0;
    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        if (field.name() != "meshPhi")
        {
            typename GeometricField<Type, PatchField, GeoMesh>
            ::Boundary& bfield = field.boundaryFieldRef();


            forAll(ptrmFields, j)
            {
                if (ptrmFields[j].first() == field.name())
                {
                    Field<Type>& fm = bfield[mId];
                    Field<Type>& fs = bfield[sId];

                    const Field<Type>& moldFaceField = ptrmFields[j].second();
                    const Field<Type>& soldFaceField = ptrsFields[j].second();
                    fm = moldFaceField;
                    fs = soldFaceField;
//                    i++;
                }
            }
        }
    }
}


template
<
    class Type,
    template<class> class PatchField,
    class GeoMesh
>
void mapGIB::MapNonOverlapFaces
(
    const List<Tuple2<word, Field<Type>>>& ptrmFields,
    const label mId,
    const labelList& dNonOverFaces,
    const labelList& nonOverFaces
)
{
    HashTable<const GeometricField<Type, PatchField, GeoMesh>*> fields
    (
        mesh_.thisDb().objectRegistry::template
            lookupClass<GeometricField<Type, PatchField, GeoMesh>>()
    );

//    label i = 0;
    for
    (
        typename HashTable<const GeometricField<Type, PatchField, GeoMesh>*>::
            iterator fieldIter = fields.begin();
        fieldIter != fields.end();
        ++fieldIter
    )
    {
        GeometricField<Type, PatchField, GeoMesh>& field =
            const_cast<GeometricField<Type, PatchField, GeoMesh>&>
            (*fieldIter());

        if (field.name() != "meshPhi")
        {
            typename GeometricField<Type, PatchField, GeoMesh>
            ::Boundary& bfield = field.boundaryFieldRef();


            forAll(ptrmFields, j)
            {
                if (ptrmFields[j].first() == field.name())
                {
                    Field<Type>& fm = bfield[mId];

                    const Field<Type>& moldFaceField = ptrmFields[j].second();

                    Type value = gAverage(moldFaceField);

                    forAll(nonOverFaces, fI)
                    {
                        label dfII = dNonOverFaces[fI];
                        label fII = nonOverFaces[fI];

                        //Pout<< field.name() << tab << dfII << tab << fII << tab
                        //     << endl;

                        if (dfII!=-1)
                        {
                            if (fII!=-1)
                            {
                                fm[dfII] = moldFaceField[fII];
                            }
                            else
                            {
                                fm[dfII] = value;
                            }
                        }
                    }
//                    i++;
                }
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// end namespace FOAM
}
// ************************************************************************* //
