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
    (c) 2011-2013 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "mappedPatches/mappedPolyPatch/mappedPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, mappedPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, mappedPolyPatch, dictionary);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPolyPatch::mappedPolyPatch
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
    mappedPatchBase(static_cast<const directPolyPatch&>(*this))
{
    //  mapped is not constraint type so add mapped group explicitly
    if (findIndex(inGroups(), typeName) == -1)
    {
        inGroups().append(typeName);
    }
}


Foam::mappedPolyPatch::mappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const mappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vectorField& offset,
    const polyBoundaryMesh& bm
)
:
    directPolyPatch(name, size, start, index, bm, typeName),
    mappedPatchBase
    (
        static_cast<const directPolyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::mappedPolyPatch::mappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const word& sampleRegion,
    const mappedPatchBase::sampleMode mode,
    const word& samplePatch,
    const vector& offset,
    const polyBoundaryMesh& bm
)
:
    directPolyPatch(name, size, start, index, bm, typeName),
    mappedPatchBase
    (
        static_cast<const directPolyPatch&>(*this),
        sampleRegion,
        mode,
        samplePatch,
        offset
    )
{}


Foam::mappedPolyPatch::mappedPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    directPolyPatch(name, dict, index, bm, patchType),
    mappedPatchBase(*this, dict)
{
    //  mapped is not constraint type so add mapped group explicitly
    if (findIndex(inGroups(), typeName) == -1)
    {
        inGroups().append(typeName);
    }
}


Foam::mappedPolyPatch::mappedPolyPatch
(
    const mappedPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    directPolyPatch(pp, bm),
    mappedPatchBase(*this, pp)
{}


Foam::mappedPolyPatch::mappedPolyPatch
(
    const mappedPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    directPolyPatch(pp, bm, index, newSize, newStart),
    mappedPatchBase(*this, pp)
{}


Foam::mappedPolyPatch::mappedPolyPatch
(
    const mappedPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    directPolyPatch(pp, bm, index, mapAddressing, newStart),
    mappedPatchBase(*this, pp, mapAddressing)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPolyPatch::~mappedPolyPatch()
{
    mappedPatchBase::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    directPolyPatch::initCalcGeometry(pBufs);
}


void Foam::mappedPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    directPolyPatch::calcGeometry(pBufs);
    mappedPatchBase::clearOut();
}


void Foam::mappedPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    directPolyPatch::initMovePoints(pBufs, p);
}


void Foam::mappedPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    directPolyPatch::movePoints(pBufs, p);
    mappedPatchBase::clearOut();
}


void Foam::mappedPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    directPolyPatch::initUpdateMesh(pBufs);
}


void Foam::mappedPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    directPolyPatch::updateMesh(pBufs);
    mappedPatchBase::clearOut();
}


void Foam::mappedPolyPatch::write(Ostream& os) const
{
    directPolyPatch::write(os);
    mappedPatchBase::write(os);
}


// ************************************************************************* //
