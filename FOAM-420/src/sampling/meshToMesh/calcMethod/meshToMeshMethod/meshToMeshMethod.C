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
    (c) 2013-2014 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshToMesh/calcMethod/meshToMeshMethod/meshToMeshMethod.H"
#include "tetOverlapVolume/tetOverlapVolume.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "db/Time/Time.H"
#include "meshes/treeBoundBox/treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshToMeshMethod, 0);
    defineRunTimeSelectionTable(meshToMeshMethod, components);
}

Foam::scalar Foam::meshToMeshMethod::tolerance_ = 1e-6;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::labelList Foam::meshToMeshMethod::maskCells() const
{
    boundBox intersectBb
    (
        max(src_.bounds().min(), tgt_.bounds().min()),
        min(src_.bounds().max(), tgt_.bounds().max())
    );

    intersectBb.inflate(0.01);

    const cellList& srcCells = src_.cells();
    const faceList& srcFaces = src_.faces();
    const pointField& srcPts = src_.points();

    DynamicList<label> cells(src_.nCells());
    forAll(srcCells, srcI)
    {
        boundBox cellBb(srcCells[srcI].points(srcFaces, srcPts), false);
        if (intersectBb.overlaps(cellBb))
        {
            cells.append(srcI);
        }
    }

    if (debug)
    {
        Pout<< "participating source mesh cells: " << cells.size() << endl;
    }

    return labelList(cells, true);
}


bool Foam::meshToMeshMethod::intersect
(
    const label srcCelli,
    const label tgtCelli
) const
{
    scalar threshold = tolerance_*src_.cellVolumes()[srcCelli];

    tetOverlapVolume overlapEngine;

    treeBoundBox bbTgtCell(tgt_.points(), tgt_.cellPoints()[tgtCelli]);

    return overlapEngine.cellCellOverlapMinDecomp
    (
        src_,
        srcCelli,
        tgt_,
        tgtCelli,
        bbTgtCell,
        threshold
    );
}


Foam::scalar Foam::meshToMeshMethod::interVol
(
    const label srcCelli,
    const label tgtCelli
) const
{
    tetOverlapVolume overlapEngine;

    treeBoundBox bbTgtCell(tgt_.points(), tgt_.cellPoints()[tgtCelli]);

    scalar vol = overlapEngine.cellCellOverlapVolumeMinDecomp
    (
        src_,
        srcCelli,
        tgt_,
        tgtCelli,
        bbTgtCell
    );

    return vol;
}


Foam::Tuple2<Foam::scalar, Foam::point>
Foam::meshToMeshMethod::interVolAndCentroid
(
    const label srcCellI,
    const label tgtCellI
)
{
    tetOverlapVolume overlapEngine;

    treeBoundBox bbTgtCell(tgt_.points(), tgt_.cellPoints()[tgtCellI]);

    Tuple2<scalar, point> volAndInertia =
    overlapEngine.cellCellOverlapMomentMinDecomp
    (
        src_,
        srcCellI,
        tgt_,
        tgtCellI,
        bbTgtCell
    );

    // Convert from inertia to centroid
    if (volAndInertia.first() <= ROOTVSMALL)
    {
        volAndInertia.first() = 0.0;
        volAndInertia.second() = Zero;
    }
    else
    {
        volAndInertia.second() /= volAndInertia.first();
    }

    return volAndInertia;
}

void Foam::meshToMeshMethod::appendNbrCells
(
    const label celli,
    const polyMesh& mesh,
    const DynamicList<label>& visitedCells,
    FIFOStack<label>& nbrCellIDs
) const
{
    const labelList& cellPts = mesh.cellPoints()[celli];

    forAll(cellPts, cPtI)
    {
        label pointi = cellPts[cPtI];
        const labelList& nbrCells = mesh.pointCells()[pointi];

        // filter out cells already visited from cell neighbours
        forAll(nbrCells, i)
        {
            label nbrCelli = nbrCells[i];

            if (findIndex(visitedCells, nbrCelli) == -1)
            {
                nbrCellIDs.push(nbrCelli);
            }
        }
    }
}


bool Foam::meshToMeshMethod::initialise
(
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght,
    labelListList& tgtToSrcAddr,
    scalarListList& tgtToSrcWght
) const
{
    srcToTgtAddr.setSize(src_.nCells());
    srcToTgtWght.setSize(src_.nCells());
    tgtToSrcAddr.setSize(tgt_.nCells());
    tgtToSrcWght.setSize(tgt_.nCells());

    if (!src_.nCells())
    {
        return false;
    }
    else if (!tgt_.nCells())
    {
        if (debug)
        {
            Pout<< "mesh interpolation: have " << src_.nCells() << " source "
                << " cells but no target cells" << endl;
        }

        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshToMeshMethod::meshToMeshMethod
(
    const polyMesh& src,
    const polyMesh& tgt,
    const dictionary& dict
)
:
    src_(src),
    tgt_(tgt),
    dict_(dict),
    V_(0.0)
{
    if (!src_.nCells() || !tgt_.nCells())
    {
        if (debug)
        {
            Pout<< "mesh interpolation: cells not on processor: Source cells = "
                << src_.nCells() << ", target cells = " << tgt_.nCells()
                << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMeshMethod::~meshToMeshMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshToMeshMethod::writeConnectivity
(
    const polyMesh& mesh1,
    const polyMesh& mesh2,
    const labelListList& mesh1ToMesh2Addr
) const
{
    Pout<< "Source size = " << mesh1.nCells() << endl;
    Pout<< "Target size = " << mesh2.nCells() << endl;

    word fName("addressing_" + mesh1.name() + "_to_" + mesh2.name());

    if (Pstream::parRun())
    {
        fName = fName +  "_proc" + Foam::name(Pstream::myProcNo());
    }

    OFstream os(src_.time().path()/fName + ".obj");

    label vertI = 0;
    forAll(mesh1ToMesh2Addr, i)
    {
        const labelList& addr = mesh1ToMesh2Addr[i];
        forAll(addr, j)
        {
            label celli = addr[j];
            const vector& c0 = mesh1.cellCentres()[i];

            const cell& c = mesh2.cells()[celli];
            const pointField pts(c.points(mesh2.faces(), mesh2.points()));
            forAll(pts, j)
            {
                const point& p = pts[j];
                os  << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
                vertI++;
                os  << "v " << c0.x() << ' ' << c0.y() << ' ' << c0.z()
                    << nl;
                vertI++;
                os  << "l " << vertI - 1 << ' ' << vertI << nl;
            }
        }
    }
}


// ************************************************************************* //
