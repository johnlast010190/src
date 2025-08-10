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
    (c) 2016-2018 OpenCFD Ltd.
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class GeoField>
Foam::tmp<GeoField> Foam::getField
(
    const IOobject* io,
    const typename GeoField::Mesh& mesh,
    const bool syncPar
)
{
    if (io)
    {
        return tmp<GeoField>::New(*io, mesh);
    }

    return tmp<GeoField>(nullptr);
}


template<class GeoField>
Foam::tmp<GeoField> Foam::getField
(
    const IOobject* io,
    const fvMeshSubsetProxy& proxy,
    const bool syncPar
)
{
    return
        proxy.interpolate
        (
            getField<GeoField>(io, proxy.baseMesh(), syncPar)
        );
}


template<class GeoField>
Foam::tmp<GeoField> Foam::getField
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const word& fieldName,
    const bool syncPar
)
{
    // Can do something with syncPar on failure ...

    return getField<GeoField>(objects.cfindObject(fieldName), mesh, syncPar);
}


template<class GeoField>
Foam::tmp<GeoField> Foam::getField
(
    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const word& fieldName,
    const bool syncPar
)
{
    // Can do something with syncPar on failure ...

    return getField<GeoField>(objects.cfindObject(fieldName), proxy, syncPar);
}


template<class GeoField>
Foam::PtrList<const GeoField> Foam::readFields
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const wordRes& selection
)
{
    const bool syncPar = true;

    // Available fields of type GeoField, sorted order
    wordList fieldNames;
    if (selection.empty())
    {
        fieldNames = objects.sortedNames();
    }
    else
    {
        const wordList all_names = objects.sortedNames();
        DynamicList<word> temp(all_names.size());
        for (auto name : all_names)
        {
            if (selection(name))
            {
                temp.append(name);
            }
        }
        fieldNames = temp;
    }

    // Construct the fields
    PtrList<const GeoField> fields(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        auto tfield =
            getField<GeoField>(mesh, objects, fieldName, syncPar);

        if (tfield.valid())
        {
            fields.set(nFields++, tfield.ptr());
        }
    }

    fields.resize(nFields);
    return fields;
}


template<class GeoField>
Foam::PtrList<const GeoField> Foam::readFields
(
    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const wordRes& selection
)
{
    const bool syncPar = true;

    wordList fieldNames;
    if (selection.empty())
    {
        fieldNames = objects.sortedNames();
    }
    else
    {
        const wordList all_names = objects.sortedNames();
        DynamicList<word> temp(all_names.size());
        for (auto name : all_names)
        {
            if (selection(name))
            {
                temp.append(name);
            }
        }
        fieldNames = temp;
    }

    // Construct the fields
    PtrList<const GeoField> fields(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        auto tfield =
            getField<GeoField>(proxy, objects, fieldName, syncPar);

        if (tfield.valid())
        {
            fields.set(nFields++, tfield.ptr());
        }
    }

    fields.resize(nFields);
    return fields;
}


// ************************************************************************* //
