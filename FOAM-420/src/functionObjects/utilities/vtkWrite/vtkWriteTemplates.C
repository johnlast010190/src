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
    (c) 2018 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField>
Foam::label Foam::functionObjects::vtkWrite::writeVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,
    const fvMeshSubset& proxy,
    const wordHashSet& acceptField
) const
{
    // TODO:  Investigate fvMeshProxy from 1906, and check that this work-around is okay
    // const fvMesh& baseMesh = proxy.baseMesh();
    const fvMesh& baseMesh = proxy.baseMesh();

    label count = 0;

    // TODO:  I assume what this is doing is going through the wordHashSet, picking out
    // the keys, and looping over the ones that are GeoFields.  It then iterates over
    // the words that sortedNames() spits out.
    // We don't have range-based for loops or the correct template for sortedNames (which
    // would have to be implemented in hashSet and objectRegistry, respectively)
    // for (const word& fieldName : baseMesh.sortedNames<GeoField>(acceptField))
    DynamicList<word> acceptFieldList;
    forAllConstIter(wordHashSet, acceptField, iter)
    {
        acceptFieldList.append(iter.key());
    }
    // auto test = baseMesh.sortedNames<GeoField>(acceptFieldList);
    forAllConstIter(DynamicList<word>, acceptFieldList, fieldName)
    {
        bool ok = false;
        const auto* fieldptr = &baseMesh.lookupObject<GeoField>(*fieldName);
        // const auto* fieldptr = baseMesh.lookupObject<GeoField>(fieldName);  // 1906

        if (!fieldptr)
        {
            continue;
        }

        // TODO:  Investigate fvMeshProxy from 1906, and check that this work-around is okay
        // auto tfield = fvMeshSubsetProxy::interpolate(proxy, *fieldptr);
        auto tfield = proxy.interpolate(*fieldptr);

        const auto& field = tfield();

        // Internal
        if (internalWriter.valid())
        {
            ok = true;
            internalWriter->write(field);
        }

        // Boundary
//        label writeri = 0;
        for (vtk::patchWriter& writer : patchWriters)
        {
            ok = true;
            writer.write(field);
//            ++writeri;
        }

        if (ok)
        {
            ++count;

            if (verbose_)
            {
                if (count == 1)
                {
                    Log << "    " << GeoField::typeName << '(';
                }
                else
                {
                    Log << ' ';
                }
                Log << fieldName;
            }
        }
    }

    if (verbose_ && count)
    {
        Log << ')' << endl;
    }

    return count;
}


template<class GeoField>
Foam::label Foam::functionObjects::vtkWrite::writeVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,
    UPtrList<vtk::patchWriter>& patchWriters,
    const UPtrList<PrimitivePatchInterpolation<primitivePatch>>& patchInterps,
    const fvMeshSubset& proxy,
    const wordHashSet& acceptField
) const
{
    // TODO:  Investigate fvMeshProxy from 1906, and check that this work-around is okay
    const fvMesh& baseMesh = proxy.baseMesh();
    // const fvMesh& baseMesh = proxy;

    label count = 0;

    // TODO:  I assume what this is doing is going through the wordHashSet, picking out
    // the keys, and looping over the ones that are GeoFields.  It then iterates over
    // the words that sortedNames() spits out.
    // We don't have range-based for loops or the correct template for sortedNames (which
    // would have to be implemented in hashSet and objectRegistry, respectively)
    // for (const word& fieldName : baseMesh.sortedNames<GeoField>(acceptField))
    DynamicList<word> acceptFieldList;
    forAllConstIter(wordHashSet, acceptField, iter)
    {
        acceptFieldList.append(iter.key());
    }
    // auto test = baseMesh.sortedNames<GeoField>(acceptFieldList);
    forAllConstIter(DynamicList<word>, acceptFieldList, fieldName)
    // for (const word& fieldName : baseMesh.sortedNames<GeoField>(acceptField))
    {
        bool ok = false;
        const auto* fieldptr = &baseMesh.lookupObject<GeoField>(*fieldName);
        // const auto* fieldptr = baseMesh.lookupObject<GeoField>(fieldName);  // 1906

        if (!fieldptr)
        {
            continue;
        }

        // TODO:  Investigate fvMeshProxy from 1906, and check that this work-around is okay
        // auto tfield = fvMeshSubsetProxy::interpolate(proxy, *fieldptr);
        auto tfield = proxy.interpolate(*fieldptr);
        const auto& field = tfield();

        // Internal
        if (internalWriter.valid() && pInterp.valid())
        {
            ok = true;
            internalWriter->write(field, *pInterp);
        }

        // Boundary
        label writeri = 0;
        for (vtk::patchWriter& writer : patchWriters)
        {
            if (writeri < patchInterps.size() && patchInterps.set(writeri))
            {
                ok = true;
                writer.write(field, patchInterps[writeri]);
            }
            ++writeri;
        }

        if (ok)
        {
            ++count;

            if (verbose_)
            {
                if (count == 1)
                {
                    Log << "    " << GeoField::typeName << "->point(";
                }
                else
                {
                    Log << ' ';
                }
                Log << fieldName;
            }
        }
    }

    if (verbose_ && count)
    {
        Log << ')' << endl;
    }

    return count;
}


// ************************************************************************* //
