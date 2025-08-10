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
    (c) 2011-2015 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPointPatch/cyclicAMIPointPatch.H"
#include "meshes/pointMesh/pointBoundaryMesh/pointBoundaryMesh.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPointPatch, 0);
    addToRunTimeSelectionTable
    (
        facePointPatch,
        cyclicAMIPointPatch,
        polyPatch
    );
    addNamedToRunTimeSelectionTable
    (
        facePointPatch,
        cyclicAMIPointPatch,
        polyPatch,
        cyclicPeriodicAMI
    );
    addNamedToRunTimeSelectionTable
    (
        facePointPatch,
        cyclicAMIPointPatch,
        polyPatch,
        discreteMixingPlane
    );
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::cyclicAMIPointPatch::initCalcGeometry(PstreamBuffers&)
{}


void Foam::cyclicAMIPointPatch::calcGeometry(PstreamBuffers&)
{}


void Foam::cyclicAMIPointPatch::initMovePoints
(
    PstreamBuffers&,
    const pointField&
)
{}


void Foam::cyclicAMIPointPatch::movePoints(PstreamBuffers&, const pointField&)
{}


void Foam::cyclicAMIPointPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::initUpdateMesh(pBufs);
//    cyclicAMIPointPatch::initCalcGeometry(pBufs);
}


void Foam::cyclicAMIPointPatch::updateMesh(PstreamBuffers& pBufs)
{
    facePointPatch::updateMesh(pBufs);
//    cyclicAMIPointPatch::calcGeometry(pBufs);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cyclicAMIPointPatch::cyclicAMIPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    cyclicAMIPolyPatch_(refCast<const cyclicAMIPolyPatch>(patch))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMIPointPatch::~cyclicAMIPointPatch()
{}


// ************************************************************************* //
