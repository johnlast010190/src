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
    (c) 2011-2014 OpenFOAM Foundation
    (c) 2015-2018 OpenCFD Ltd.
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyPatches/directPolyPatch/directPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "meshes/polyMesh/polyBoundaryMesh/polyBoundaryMesh.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/primitiveMesh/primitiveMesh.H"
#include "fields/Fields/Field/SubField.H"
#include "db/dictionary/entry/entry.H"
#include "db/dictionary/dictionary.H"
#include "fields/pointPatchFields/pointPatchField/pointPatchField.H"
#include "meshes/polyMesh/polyPatches/directPolyPatch/directPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, directPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, directPolyPatch, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::directPolyPatch::movePoints(const pointField& p)
{
    primitivePatch::clearGeom();
}


void Foam::directPolyPatch::movePoints(PstreamBuffers&, const pointField& p)
{
    primitivePatch::clearGeom();
}


void Foam::directPolyPatch::updateMesh(PstreamBuffers&)
{
    primitivePatch::clearGeom();
    clearAddressing();
}


void Foam::directPolyPatch::updateGIB()
{}


void Foam::directPolyPatch::clearGeom()
{
    primitivePatch::clearGeom();
}


void Foam::directPolyPatch::rename(const wordList& newNames)
{
    name_ = newNames[index_];
}


void Foam::directPolyPatch::reorder(const labelUList& oldToNewIndex)
{
    index_ = oldToNewIndex[index_];
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directPolyPatch::directPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, size, start, index, bm, patchType),
    primitivePatch
    (
        faceSubList(bm.mesh().faces(), size, start),
        bm.mesh().points()
    ),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr)
{}


Foam::directPolyPatch::directPolyPatch
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
    polyPatch(name, size, start, index, bm, physicalType, inGroups),
    primitivePatch
    (
        faceSubList(bm.mesh().faces(), size, start),
        bm.mesh().points()
    ),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr)
{}


Foam::directPolyPatch::directPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    polyPatch(name, readLabel(dict.lookup("startFace")), dict, index, bm, patchType),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            readLabel(dict.lookup("nFaces")),
            readLabel(dict.lookup("startFace"))
        ),
        bm.mesh().points()
    ),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr)
{}


Foam::label Foam::directPolyPatch::offset() const
{
    return start_ - boundaryMesh().start();
}


Foam::directPolyPatch::directPolyPatch
(
    const directPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            pp.size(),
            pp.start()
        ),
        bm.mesh().points()
    ),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr)
{}


Foam::directPolyPatch::directPolyPatch
(
    const directPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            newSize,
            newStart
        ),
        bm.mesh().points()
    ),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr)
{}


Foam::directPolyPatch::directPolyPatch
(
    const polyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            newSize,
            newStart
        ),
        bm.mesh().points()
    ),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr)
{
    this->start_ = newStart;
}

Foam::directPolyPatch::directPolyPatch
(
    const directPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    polyPatch(pp, bm, index, mapAddressing, newStart),
    primitivePatch
    (
        faceSubList
        (
            bm.mesh().faces(),
            mapAddressing.size(),
            newStart
        ),
        bm.mesh().points()
    ),
    boundaryMesh_(bm),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr)
{}


Foam::directPolyPatch::directPolyPatch(const directPolyPatch& p)
:
    polyPatch(p),
    primitivePatch(p),
    boundaryMesh_(p.boundaryMesh_),
    faceCellsPtr_(nullptr),
    mePtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directPolyPatch::~directPolyPatch()
{
    clearAddressing();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
const Foam::polyBoundaryMesh& Foam::directPolyPatch::boundaryMesh() const
{
    return boundaryMesh_;
}


const Foam::vectorField::subField Foam::directPolyPatch::faceCentres() const
{
    return patchSlice(boundaryMesh().mesh().faceCentres());
}


const Foam::vectorField::subField Foam::directPolyPatch::faceAreas() const
{
    return patchSlice(boundaryMesh().mesh().faceAreas());
}


const Foam::scalarField::subField Foam::directPolyPatch::magFaceAreas() const
{
    return patchSlice(boundaryMesh().mesh().magFaceAreas());
}


// Return the patch face neighbour cell centres
Foam::tmp<Foam::vectorField> Foam::directPolyPatch::faceCellCentres() const
{
    tmp<vectorField> tcc(new vectorField(size()));
    vectorField& cc = tcc.ref();

    // get reference to global cell centres
    const vectorField& gcc = boundaryMesh_.mesh().cellCentres();

    const labelUList& faceCells = this->faceCells();

    forAll(faceCells, facei)
    {
        cc[facei] = gcc[faceCells[facei]];
    }

    return tcc;
}


const Foam::labelUList& Foam::directPolyPatch::faceCells() const
{
    if (!faceCellsPtr_)
    {
        faceCellsPtr_ = new labelList::subList
        (
            patchSlice(boundaryMesh().mesh().faceOwner())
        );
    }

    return *faceCellsPtr_;
}


const Foam::labelList& Foam::directPolyPatch::meshEdges() const
{
    if (!mePtr_)
    {
        mePtr_ =
            new labelList
            (
                primitivePatch::meshEdges
                (
                    boundaryMesh().mesh().edges(),
                    boundaryMesh().mesh().pointEdges()
                )
            );
    }

    return *mePtr_;
}


void Foam::directPolyPatch::clearAddressing()
{
    primitivePatch::clearTopology();
    primitivePatch::clearPatchMeshAddr();
    deleteDemandDrivenData(faceCellsPtr_);
    deleteDemandDrivenData(mePtr_);
}


void Foam::directPolyPatch::write(Ostream& os) const
{
    os.writeEntry("type", this->type());
    patchIdentifier::write(os);
    os.writeEntry("nFaces", this->size());
    os.writeEntry("startFace", this->start());
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::directPolyPatch::operator=(const polyPatch& p)
{
    if (isA<directPolyPatch>(p))
    {
        const directPolyPatch& dpp =
            refCast<const directPolyPatch>(p);
        clearAddressing();
        patchIdentifier::operator=(dpp);
        primitivePatch::operator=(dpp);
        start_ = dpp.start();
    }
}

const Foam::face& Foam::directPolyPatch::operator[](const label i) const
{
    const Foam::face& f = primitivePatch::operator[](i);
    return f;
}
Foam::face& Foam::directPolyPatch::operator[](const label i)
{
    face& f = primitivePatch::operator[](i);
    return f;
}
// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const directPolyPatch& p)
{
    p.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
