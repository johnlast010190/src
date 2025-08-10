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

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyPatches/polyPatch/polyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/primitiveMesh/primitiveMesh.H"
#include "fields/Fields/Field/SubField.H"
#include "db/dictionary/entry/entry.H"
#include "db/dictionary/dictionary.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyPatch, 0);

    int polyPatch::disallowGenericPolyPatch
    (
        debug::debugSwitch("disallowGenericPolyPatch", 0)
    );

    defineRunTimeSelectionTable(polyPatch, word);
    defineRunTimeSelectionTable(polyPatch, indirect);
    defineRunTimeSelectionTable(polyPatch, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyPatch::polyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    patchIdentifier(name, index),
    start_(start)
{
    if
    (
        patchType != word::null
     && constraintType(patchType)
     && findIndex(inGroups(), patchType) == -1
    )
    {
        inGroups().append(patchType);
    }
}


Foam::polyPatch::polyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& physicalType,
    const wordList& inGroups
)
:
    patchIdentifier(name, index, physicalType, inGroups),
    start_(start)
{}


Foam::polyPatch::polyPatch
(
    const word& name,
    const label size,
    const label start,
    const label zoneI,
    const word zPType,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    patchIdentifier(name, index),
    start_(start)
{
    if
    (
        patchType != word::null
     && constraintType(patchType)
     && findIndex(inGroups(), patchType) == -1
    )
    {
        inGroups().append(patchType);
    }
}


Foam::polyPatch::polyPatch
(
    const word& name,
    const label start,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    patchIdentifier(name, dict, index),
    start_(start)
{
    if
    (
        patchType != word::null
     && constraintType(patchType)
     && findIndex(inGroups(), patchType) == -1
    )
    {
        inGroups().append(patchType);
    }
}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    patchIdentifier(pp),
    start_(pp.start())
{}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    patchIdentifier(pp, index),
    start_(newStart)
{}


Foam::polyPatch::polyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    patchIdentifier(pp, index),
    start_(newStart)
{}


Foam::polyPatch::polyPatch(const polyPatch& p)
:
    patchIdentifier(p),
    start_(p.start())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::polyPatch::constraintType(const word& pt)
{
    return pointPatchField<scalar>::pointPatchConstructorTable_().count(pt) > 0;
}


Foam::wordList Foam::polyPatch::constraintTypes()
{
    wordList cTypes(dictionaryConstructorTable_().size());

    label i = 0;

    for
    (
        dictionaryConstructorTable::iterator cstrIter =
            dictionaryConstructorTable_().begin();
        cstrIter != dictionaryConstructorTable_().end();
        ++cstrIter
    )
    {
        if (constraintType(cstrIter->first))
        {
            cTypes[i++] = cstrIter->first;
        }
    }

    cTypes.setSize(i);

    return cTypes;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyPatch& p)
{
    p.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
