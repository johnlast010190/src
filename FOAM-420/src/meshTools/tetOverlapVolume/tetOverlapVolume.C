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
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2012-2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "tetOverlapVolume/tetOverlapVolume.H"
#include "meshes/primitiveShapes/tetrahedron/tetrahedron.H"
#include "meshes/primitiveShapes/tetrahedron/tetPoints.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/treeBoundBox/treeBoundBox.H"
#include "algorithms/indexedOctree/indexedOctree.H"
#include "algorithms/indexedOctree/treeDataCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tetOverlapVolume, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetOverlapVolume::tetOverlapVolume()
{}


// * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * * //

Foam::treeBoundBox Foam::tetOverlapVolume::pyrBb
(
    const pointField& points,
    const face& f,
    const point& fc
)
{
    treeBoundBox bb(fc);
    bb.add(points, f);

    return bb;
}


// * * * * * * * * * * * Public Member Functions  * * * * * * * * * * * * * //

bool Foam::tetOverlapVolume::cellCellOverlapMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,
    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB,
    const scalar threshold
) const
{
    hasOverlapOp overlapCheckOp(threshold);
    cellCellOverlapMinDecomp<hasOverlapOp>
    (
        meshA,
        cellAI,
        meshB,
        cellBI,
        cellBbB,
        overlapCheckOp
    );

    return overlapCheckOp.ok_;
}


Foam::scalar Foam::tetOverlapVolume::cellCellOverlapVolumeMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,

    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB
) const
{
    sumOverlapOp overlapSumOp;
    cellCellOverlapMinDecomp<sumOverlapOp>
    (
        meshA,
        cellAI,
        meshB,
        cellBI,
        cellBbB,
        overlapSumOp
    );

    return overlapSumOp.iop_.vol_;
}


Foam::Tuple2<Foam::scalar, Foam::point>
Foam::tetOverlapVolume::cellCellOverlapMomentMinDecomp
(
    const primitiveMesh& meshA,
    const label cellAI,

    const primitiveMesh& meshB,
    const label cellBI,
    const treeBoundBox& cellBbB
) const
{
    sumOverlapMomentOp overlapSumOp;
    cellCellOverlapMinDecomp<sumOverlapMomentOp>
    (
        meshA,
        cellAI,
        meshB,
        cellBI,
        cellBbB,
        overlapSumOp
    );

    return overlapSumOp.iop_.vol_;
}


Foam::labelList Foam::tetOverlapVolume::overlappingCells
(
    const polyMesh& fromMesh,
    const polyMesh& toMesh,
    const label iTo
) const
{
    const indexedOctree<treeDataCell>& treeA = fromMesh.cellTree();

    treeBoundBox bbB(toMesh.points(), toMesh.cellPoints()[iTo]);

    return treeA.findBox(bbB);
}


// ************************************************************************* //
