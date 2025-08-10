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
#include "db/dictionary/dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& patchType,
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
{
    if (debug)
    {
        InfoInFunction << "Constructing polyPatch" << endl;
    }

    const auto ctor = ctorTableLookup("polyPatch type", wordConstructorTable_(), patchType);
    return autoPtr<polyPatch>
    (
        ctor
        (
            name,
            size,
            start,
            index,
            bm,
            patchType
        )
    );
}


Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& patchType,
    const word& name,
    const label size,
    const label zoneId,
    const word  zPType,
    const label index,
    const polyBoundaryMesh& bm
)
{
    if (debug)
    {
        Info<< "polyPatch::New(const word&, const word&, const label, "
               "const label, const label, const zoneBoundaryMesh&) : "
               "constructing polyPatch"
            << endl;
    }

    const auto ctor = ctorTableLookup("polyPatch type", indirectConstructorTable_(), patchType);
    return autoPtr<polyPatch>
    (
        ctor
        (
            name,
            size,
            zoneId,
            zPType,
            index,
            bm,
            patchType
        )
    );
}


Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
{
    if (debug)
    {
        InfoInFunction << "Constructing polyPatch" << endl;
    }

    word patchType(dict.lookup("type"));
    dict.readIfPresent("geometricType", patchType);

    return polyPatch::New(patchType, name, dict, index, bm);
}


Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
    const word& patchType,
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
{
    if (debug)
    {
        InfoInFunction << "Constructing polyPatch" << endl;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTable_().find(patchType);

    if (cstrIter == dictionaryConstructorTable_().end())
    {
        if (!disallowGenericPolyPatch)
        {
            cstrIter = dictionaryConstructorTable_().find("genericPatch");
        }

        if (cstrIter == dictionaryConstructorTable_().end())
        {
            FatalIOErrorInFunction(dict)
                << "Unknown polyPatch type "
                << patchType << " for patch " << name << nl << nl
                << exit(FatalIOError);
        }
    }

    return autoPtr<polyPatch>(cstrIter->second(name, dict, index, bm, patchType));
}


// ************************************************************************* //
