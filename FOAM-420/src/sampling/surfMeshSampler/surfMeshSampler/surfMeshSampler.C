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
    (c) 2016 OpenCFD Ltd.

\*---------------------------------------------------------------------------*/

#include "surfMeshSampler/surfMeshSampler/surfMeshSampler.H"
#include "MeshedSurface/MeshedSurfaces.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshSampler, 0);
    defineRunTimeSelectionTable(surfMeshSampler, word);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::surfMesh&
Foam::surfMeshSampler::getOrCreateSurfMesh() const
{
    if (!mesh().foundObject<surfMesh>(name()))
    {
        surfMesh* ptr = new surfMesh
        (
            IOobject
            (
                name(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            xferCopy(pointField()), // initially no points
            xferCopy(faceList()),   // initially no faces
            name()
        );
        ptr->setWriteOption(IOobject::NO_WRITE);

        mesh().objectRegistry::store(ptr);
    }

    return const_cast<surfMesh&>(surface());
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfMeshSampler>
Foam::surfMeshSampler::New
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    const word sampleType(dict.lookup("type"));

    const auto ctor = ctorTableLookup("sample type", wordConstructorTable_(), sampleType);
    return autoPtr<surfMeshSampler>(ctor(name, mesh, dict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSampler::surfMeshSampler
(
    const word& name,
    const polyMesh& mesh
)
:
    name_(name),
    mesh_(mesh)
{}


Foam::surfMeshSampler::surfMeshSampler
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    name_(name),
    mesh_(mesh)
{
    dict.readIfPresent("name", name_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMeshSampler::~surfMeshSampler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMeshSampler::create() const
{
    getOrCreateSurfMesh();
}


const Foam::surfMesh& Foam::surfMeshSampler::surface() const
{
    return mesh().lookupObject<surfMesh>(name());
}


// Demonstration of using separate tmp registry
// Foam::label Foam::surfMeshSampler::sample
// (
//     const objectRegistry& store,
//     const UList<word>& fields
// ) const
// {
//     label count = 0;
//     forAll(fields, fieldi)
//     {
//         const word& fieldName = fields[fieldi];
//         bool ok =
//         (
//             sampleAndStore(store, fieldName)
//          &&
//             (
//                 transferField<scalar>(store, fieldName)
//              || transferField<vector>(store, fieldName)
//              || transferField<sphericalTensor>(store, fieldName)
//              || transferField<symmTensor>(store, fieldName)
//              || transferField<tensor>(store, fieldName)
//             )
//         );
//
//         if (ok)
//         {
//             ++count;
//         }
//     }
//
//     return count;
// }


Foam::label Foam::surfMeshSampler::sample
(
    const UList<word>& fields
) const
{
    label count = 0;
    forAll(fields, fieldi)
    {
        if (sample(fields[fieldi]))
        {
            ++count;
        }
    }

    return count;
}


Foam::label Foam::surfMeshSampler::write(const wordReList& select) const
{
    label count =
    (
        writeFields<scalar>(select)
      + writeFields<vector>(select)
      + writeFields<sphericalTensor>(select)
      + writeFields<symmTensor>(select)
      + writeFields<tensor>(select)
    );

    if (count)
    {
        surface().write();
    }

    return count;
}



void Foam::surfMeshSampler::print(Ostream& os) const
{
    os << type();
}


// ************************************************************************* //
