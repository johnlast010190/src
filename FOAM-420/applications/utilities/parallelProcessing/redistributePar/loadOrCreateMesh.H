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
    (c) 2012 OpenFOAM Foundation

InNamespace
    Foam

Description
    Load or create (0 size) a mesh. Used in distributing meshes to a
    larger number of processors

SourceFiles
    loadOrCreateMesh.C

\*---------------------------------------------------------------------------*/

#ifndef loadOrCreateMesh_H
#define loadOrCreateMesh_H

#include "fvMesh/fvMesh.H"
#include "fvMeshSubset/fvMeshSubset.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "db/IOobjectList/IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Load (if it exists) or create zero cell mesh given an IOobject:
//      name     : regionName
//      instance : exact directory where to find mesh (i.e. does not
//                 do a findInstance
autoPtr<fvMesh> loadOrCreateMesh
(
    const IOobject& io,
    const bool removeDummyMesh=false
);


// Generic mesh-based field reading
template<class GeoField>
void readField
(
    const IOobject& io,
    const fvMesh& mesh,
    const label i,
    PtrList<GeoField>& fields
)
{
    fields.set(i, new GeoField(io, mesh));
}

// Definition of readField for GeometricFields only
template<class Type, template<class> class PatchField, class GeoMesh>
void readField
(
    const IOobject& io,
    const fvMesh& mesh,
    const label i,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields
)
{
    fields.set
    (
        i,
        new GeometricField<Type, PatchField, GeoMesh>(io, mesh, false, false)
    );
}

// Read vol or surface fields
template<class GeoField>
void readFields
(
    const boolList& haveMesh,
    const fvMesh& mesh,
    const autoPtr<fvMeshSubset>& subsetterPtr,
    IOobjectList& allObjects,
    PtrList<GeoField>& fields
)
{
    // Get my objects of type
    IOobjectList objects(allObjects.lookupClass(GeoField::typeName));

    // Check that we all have all objects
    wordList objectNames = objects.sortedNames();

    // Get master names
    wordList masterNames(objectNames);
    Pstream::scatter(masterNames);

    if (haveMesh[Pstream::myProcNo()] && objectNames != masterNames)
    {
        FatalErrorInFunction
            << "differing fields of type " << GeoField::typeName
            << " on processors." << endl
            << "Master has:" << masterNames << endl
            << Pstream::myProcNo() << " has:" << objectNames
            << abort(FatalError);
    }

    fields.setSize(masterNames.size());

    const bool oldParRun = Pstream::parRun();

    // Have master send all fields to processors that don't have a mesh
    if (Pstream::master())
    {
        forAll(masterNames, i)
        {
            const word& name = masterNames[i];
            IOobject& io = *objects[name];
            io.writeOpt() = IOobject::AUTO_WRITE;

            // do not sync BCs across processors
            Pstream::parRun() = false;

            // Load field (but not oldTime)
            readField(io, mesh, i, fields);

            Pstream::parRun() = oldParRun;

            // Create zero sized field and send
            if (subsetterPtr.valid())
            {
                tmp<GeoField> tsubfld = subsetterPtr().interpolate(fields[i]);

                // Send to all processors that don't have a mesh
                for (label procI = 1; procI < Pstream::nProcs(); procI++)
                {
                    if (!haveMesh[procI])
                    {
                        OPstream toProc
                        (
                            Pstream::commsTypes::blocking,
                            procI
                        );
                        toProc<< tsubfld();
                    }
                }
            }
        }
    }
    else if (!haveMesh[Pstream::myProcNo()])
    {
        // Don't have mesh (nor fields). Receive empty field from master.

        forAll(masterNames, i)
        {
            const word& name = masterNames[i];

            // Receive field
            IPstream fromMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            dictionary fieldDict(fromMaster);

            // do not sync BCs across processors
            Pstream::parRun() = false;

            fields.set
            (
                i,
                new GeoField
                (
                    IOobject
                    (
                        name,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    fieldDict
                )
            );

            Pstream::parRun() = oldParRun;

            //// Write it for next time round (since mesh gets written as well)
            //fields[i].write();
        }
    }
    else
    {
        // Have mesh so just try to load
        forAll(masterNames, i)
        {
            const word& name = masterNames[i];
            IOobject& io = *objects[name];
            io.writeOpt() = IOobject::AUTO_WRITE;

            // do not sync BCs across processors
            Pstream::parRun() = false;

            // Load field (but not oldTime)
            readField(io, mesh, i, fields);

            Pstream::parRun() = oldParRun;
        }
    }
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
