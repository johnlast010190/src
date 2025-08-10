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
    (c) 2011-2012 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "regionCoupled/patches/regionCoupledPolyPatch/regionCoupledPolyPatch.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(regionCoupledPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, regionCoupledPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, regionCoupledPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::regionCoupledPolyPatch::regionCoupledPolyPatch
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
    regionCoupledBase(static_cast<const directPolyPatch&>(*this))
{}


Foam::regionCoupledPolyPatch::regionCoupledPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    directPolyPatch(name, dict, index, bm, patchType),
    regionCoupledBase(*this, dict)
{}


Foam::regionCoupledPolyPatch::regionCoupledPolyPatch
(
    const regionCoupledPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    directPolyPatch(pp, bm),
    regionCoupledBase(*this, pp)
{}


Foam::regionCoupledPolyPatch::regionCoupledPolyPatch
(
    const regionCoupledPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    directPolyPatch(pp, bm, index, newSize, newStart),
    regionCoupledBase(*this, pp)
{}


Foam::regionCoupledPolyPatch::regionCoupledPolyPatch
(
    const regionCoupledPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    directPolyPatch(pp, bm, index, mapAddressing, newStart),
    regionCoupledBase(*this, pp)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionCoupledPolyPatch::~regionCoupledPolyPatch()
{
    regionCoupledBase::clearGeom();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionCoupledPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    directPolyPatch::initCalcGeometry(pBufs);
}


void Foam::regionCoupledPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    directPolyPatch::initMovePoints(pBufs, p);
}


void Foam::regionCoupledPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    directPolyPatch::movePoints(pBufs, p);
    regionCoupledBase::clearGeom();
}


void Foam::regionCoupledPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    directPolyPatch::initUpdateMesh(pBufs);
}


void Foam::regionCoupledPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    directPolyPatch::updateMesh(pBufs);
    regionCoupledBase::clearGeom();
}


void Foam::regionCoupledPolyPatch::write(Ostream& os) const
{
    directPolyPatch::write(os);
    regionCoupledBase::write(os);
}


// ************************************************************************* //
