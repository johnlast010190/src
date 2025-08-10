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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nonConformal/polyPatches/nonConformalCyclic/nonConformalCyclicPolyPatch.H"
#include "nonConformal/polyPatches/nonConformalError/nonConformalErrorPolyPatch.H"
#include "nonConformalBoundary/nonConformalBoundary.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCyclicPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, nonConformalCyclicPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        nonConformalCyclicPolyPatch,
        dictionary
    );
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalCyclicPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    cyclicPolyPatch::initCalcGeometry(pBufs);
    intersectionIsValid_ = false;
    raysIsValid_ = false;
}


void Foam::nonConformalCyclicPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (transformComplete())
    {
        static_cast<cyclicTransform&>(*this) =
            cyclicTransform
            (
                name(),
                origPatch().faceAreas(),
                *this,
                nbrPatchName(),
                nbrPatch(),
                matchTolerance()
            );
    }
    else
    {
        static_cast<cyclicTransform&>(*this) =
            cyclicTransform
            (
                name(),
                origPatch().faceCentres(),
                origPatch().faceAreas(),
                *this,
                nbrPatchName(),
                nbrPatch().origPatch().faceCentres(),
                nbrPatch().origPatch().faceAreas(),
                nbrPatch(),
                matchTolerance()
            );
    }
}


void Foam::nonConformalCyclicPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    cyclicPolyPatch::initMovePoints(pBufs, p);
    intersectionIsValid_ = false;
    raysIsValid_ = false;
}


void Foam::nonConformalCyclicPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    cyclicPolyPatch::initUpdateMesh(pBufs);
    intersectionIsValid_ = false;
    raysIsValid_ = false;
}


void Foam::nonConformalCyclicPolyPatch::clearGeom()
{
    cyclicPolyPatch::clearGeom();
    intersectionIsValid_ = false;
    raysIsValid_ = false;
}


void Foam::nonConformalCyclicPolyPatch::rename(const wordList& newNames)
{
    cyclicPolyPatch::rename(newNames);
    nonConformalCoupledPolyPatch::rename(newNames);
}


void Foam::nonConformalCyclicPolyPatch::reorder
(
    const labelUList& oldToNewIndex
)
{
    cyclicPolyPatch::reorder(oldToNewIndex);
    nonConformalCoupledPolyPatch::reorder(oldToNewIndex);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicPolyPatch(name, size, start, index, bm, patchType),
    nonConformalCoupledPolyPatch(static_cast<const polyPatch&>(*this)),
    intersectionIsValid_(false),
    intersection_(false),
    raysIsValid_(false),
    rays_(false)
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& nbrPatchName,
    const word& origPatchName,
    const cyclicTransform& transform
)
:
    cyclicPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        patchType,
        nbrPatchName,
        transform
    ),
    nonConformalCoupledPolyPatch
    (
        static_cast<const polyPatch&>(*this),
        origPatchName
    ),
    intersectionIsValid_(false),
    intersection_(false),
    raysIsValid_(false),
    rays_(false)
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicPolyPatch(name, dict, index, bm, patchType, true),
    nonConformalCoupledPolyPatch(static_cast<const polyPatch&>(*this), dict),
    intersectionIsValid_(false),
    intersection_(false),
    raysIsValid_(false),
    rays_(false)
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const nonConformalCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicPolyPatch(pp, bm),
    nonConformalCoupledPolyPatch
    (
        static_cast<const nonConformalCoupledPolyPatch&>(pp)
    ),
    intersectionIsValid_(false),
    intersection_(false),
    raysIsValid_(false),
    rays_(false)
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const nonConformalCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& origPatchName
)
:
    cyclicPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    nonConformalCoupledPolyPatch
    (
        static_cast<const polyPatch&>(*this),
        origPatchName
    ),
    intersectionIsValid_(false),
    intersection_(false),
    raysIsValid_(false),
    rays_(false)
{}


Foam::nonConformalCyclicPolyPatch::nonConformalCyclicPolyPatch
(
    const nonConformalCyclicPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicPolyPatch(pp, bm, index, mapAddressing, newStart),
    nonConformalCoupledPolyPatch
    (
        static_cast<const nonConformalCoupledPolyPatch&>(pp)
    ),
    intersectionIsValid_(false),
    intersection_(false),
    raysIsValid_(false),
    rays_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicPolyPatch::~nonConformalCyclicPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::patchToPatches::intersection&
Foam::nonConformalCyclicPolyPatch::intersection() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "The non-conformal cyclic intersection is only available to "
            << "the owner patch" << abort(FatalError);
    }

    if (!intersectionIsValid_)
    {
        const polyMesh& mesh = boundaryMesh().mesh();

        const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh);

        intersection_.update
        (
            refCast<const directPolyPatch>(origPatch()),
            ncb.patchPointNormals(origPatchID()),
            refCast<const directPolyPatch>(nbrPatch().origPatch()),
            transform()
        );

        intersectionIsValid_ = true;
    }

    return intersection_;
}


const Foam::patchToPatches::rays&
Foam::nonConformalCyclicPolyPatch::rays() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "The non-conformal cyclic rays are only available to "
            << "the owner patch" << abort(FatalError);
    }

    if (!raysIsValid_)
    {
        const polyMesh& mesh = boundaryMesh().mesh();

        const nonConformalBoundary& ncb = nonConformalBoundary::New(mesh);

        rays_.update
        (
            primitiveOldTimePatch
            (
                refCast<const directPolyPatch>(origPatch()),
                mesh.points(),
                mesh.oldPoints()
            ),
            ncb.patchPointNormals(origPatchID()),
            ncb.patchPointNormals0(origPatchID()),
            primitiveOldTimePatch
            (
                refCast<const directPolyPatch>(nbrPatch().origPatch()),
                mesh.points(),
                mesh.oldPoints()
            ),
            transform()
        );

        raysIsValid_ = true;
    }

    return rays_;
}


Foam::patchToPatch::procFace Foam::nonConformalCyclicPolyPatch::ray
(
    const scalar fraction,
    const label origFacei,
    const vector& p,
    const vector& n,
    point& nbrP
) const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    const nonConformalCyclicPolyPatch& ownerPatch =
        owner() ? *this : nbrPatch();

    auto ownerRaysMethod =
        owner()
      ? &patchToPatches::rays::srcToTgtRay
      : &patchToPatches::rays::tgtToSrcRay;

    return
        (ownerPatch.rays().*ownerRaysMethod)
        (
            primitiveOldTimePatch
            (
                refCast<const directPolyPatch>(nbrPatch().origPatch()),
                mesh.points(),
                mesh.oldPoints()
            ),
            fraction,
            origFacei,
            transform().invTransformPosition(p),
            transform().invTransform(n),
            nbrP
        );
}


void Foam::nonConformalCyclicPolyPatch::write(Ostream& os) const
{
    cyclicPolyPatch::write(os);
    nonConformalCoupledPolyPatch::write(os);
}


// ************************************************************************* //
