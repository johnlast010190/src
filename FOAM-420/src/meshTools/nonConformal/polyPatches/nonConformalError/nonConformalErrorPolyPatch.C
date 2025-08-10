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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nonConformal/polyPatches/nonConformalError/nonConformalErrorPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalErrorPolyPatch, 0);

    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalErrorPolyPatch,
        word
    );

    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalErrorPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalErrorPolyPatch::rename(const wordList& newNames)
{
    nonConformalPolyPatch::rename(newNames);
}


void Foam::nonConformalErrorPolyPatch::reorder(const labelUList& oldToNewIndex)
{
    directPolyPatch::reorder(oldToNewIndex);
    nonConformalPolyPatch::reorder(oldToNewIndex);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    directPolyPatch(name, size, start, index, bm, patchType),
    nonConformalPolyPatch(static_cast<const polyPatch&>(*this))
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& origPatchName
)
:
    directPolyPatch(name, size, start, index, bm, patchType),
    nonConformalPolyPatch
    (
        static_cast<const polyPatch&>(*this), origPatchName
    )
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    directPolyPatch(name, dict, index, bm, patchType),
    nonConformalPolyPatch(static_cast<const polyPatch&>(*this), dict)
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const nonConformalErrorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    directPolyPatch(pp, bm),
    nonConformalPolyPatch(static_cast<const nonConformalPolyPatch&>(pp))
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const nonConformalErrorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& origPatchName
)
:
    directPolyPatch(pp, bm, index, newSize, newStart),
    nonConformalPolyPatch
    (
        static_cast<const polyPatch&>(pp), origPatchName
    )
{}


Foam::nonConformalErrorPolyPatch::nonConformalErrorPolyPatch
(
    const nonConformalErrorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    directPolyPatch(pp, bm, index, mapAddressing, newStart),
    nonConformalPolyPatch(static_cast<const nonConformalPolyPatch&>(pp))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalErrorPolyPatch::~nonConformalErrorPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonConformalErrorPolyPatch::write(Ostream& os) const
{
    directPolyPatch::write(os);
    nonConformalPolyPatch::write(os);
}


// ************************************************************************* //
