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
    (c) 2015 OpenCFD Ltd.
    (c) 2015 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "meshToMesh/calcMethod/correctedCellVolumeWeight/correctedCellVolumeWeightMethod.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(correctedCellVolumeWeightMethod, 0);
    addToRunTimeSelectionTable
    (
        meshToMeshMethod,
        correctedCellVolumeWeightMethod,
        components
    );
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::correctedCellVolumeWeightMethod::calculateAddressing
(
    labelListList& srcToTgtCellAddr,
    scalarListList& srcToTgtCellWght,
    pointListList& srcToTgtCellVec,
    labelListList& tgtToSrcCellAddr,
    scalarListList& tgtToSrcCellWght,
    pointListList& tgtToSrcCellVec,
    const label srcSeedI,
    const label tgtSeedI,
    const labelList& srcCellIDs,
    boolList& mapFlag,
    label& startSeedI
)
{
    label srcCellI = srcSeedI;
    label tgtCellI = tgtSeedI;

    List<DynamicList<label>> srcToTgtAddr(src_.nCells());
    List<DynamicList<scalar>> srcToTgtWght(src_.nCells());
    List<DynamicList<point>> srcToTgtVec(src_.nCells());

    List<DynamicList<label>> tgtToSrcAddr(tgt_.nCells());
    List<DynamicList<scalar>> tgtToSrcWght(tgt_.nCells());
    List<DynamicList<point>> tgtToSrcVec(tgt_.nCells());

    // list of tgt cell neighbour cells
    FIFOStack<label> nbrTgtCells;

    // list of tgt cells currently visited for srcCellI to avoid multiple hits
    DynamicList<label> visitedTgtCells(10);

    // list to keep track of tgt cells used to seed src cells
    labelList seedCells(src_.nCells(), -1);
    seedCells[srcCellI] = tgtCellI;

    const scalarField& srcVol = src_.cellVolumes();
    const pointField& srcCc = src_.cellCentres();
    const pointField& tgtCc = tgt_.cellCentres();

    do
    {
        nbrTgtCells.clear();
        visitedTgtCells.clear();

        // append initial target cell and neighbours
        nbrTgtCells.push(tgtCellI);
        appendNbrCells(tgtCellI, tgt_, visitedTgtCells, nbrTgtCells);

        do
        {
            tgtCellI = nbrTgtCells.pop();
            visitedTgtCells.append(tgtCellI);

            Tuple2<scalar, point> vol = interVolAndCentroid
            (
                srcCellI,
                tgtCellI
            );

            // accumulate addressing and weights for valid intersection
            if (vol.first()/srcVol[srcCellI] > tolerance_)
            {
                // store src/tgt cell pair
                srcToTgtAddr[srcCellI].append(tgtCellI);
                srcToTgtWght[srcCellI].append(vol.first());
                srcToTgtVec[srcCellI].append(vol.second()-tgtCc[tgtCellI]);

                tgtToSrcAddr[tgtCellI].append(srcCellI);
                tgtToSrcWght[tgtCellI].append(vol.first());
                tgtToSrcVec[tgtCellI].append(vol.second()-srcCc[srcCellI]);

                appendNbrCells(tgtCellI, tgt_, visitedTgtCells, nbrTgtCells);

                // accumulate intersection volume
                V_ += vol.first();
            }
        }
        while (!nbrTgtCells.empty());

        mapFlag[srcCellI] = false;

        // find new source seed cell
        setNextCells
        (
            startSeedI,
            srcCellI,
            tgtCellI,
            srcCellIDs,
            mapFlag,
            visitedTgtCells,
            seedCells
        );
    }
    while (srcCellI != -1);

    // transfer addressing into persistent storage
    forAll(srcToTgtCellAddr, i)
    {
        srcToTgtCellAddr[i].transfer(srcToTgtAddr[i]);
        srcToTgtCellWght[i].transfer(srcToTgtWght[i]);
        srcToTgtCellVec[i].transfer(srcToTgtVec[i]);
    }

    forAll(tgtToSrcCellAddr, i)
    {
        tgtToSrcCellAddr[i].transfer(tgtToSrcAddr[i]);
        tgtToSrcCellWght[i].transfer(tgtToSrcWght[i]);
        tgtToSrcCellVec[i].transfer(tgtToSrcVec[i]);
    }


    if (debug%2)
    {
        // At this point the overlaps are still in volume so we could
        // get out the relative error
        forAll(srcToTgtCellAddr, cellI)
        {
            scalar srcVol = src_.cellVolumes()[cellI];
            scalar tgtVol = sum(srcToTgtCellWght[cellI]);

            if (mag(srcVol) > ROOTVSMALL && mag((tgtVol-srcVol)/srcVol) > 1e-6)
            {
                WarningInFunction
                    << "At cell " << cellI << " cc:"
                    << src_.cellCentres()[cellI]
                    << " vol:" << srcVol
                    << " total overlap volume:" << tgtVol
                    << endl;
            }
        }

        forAll(tgtToSrcCellAddr, cellI)
        {
            scalar tgtVol = tgt_.cellVolumes()[cellI];
            scalar srcVol = sum(tgtToSrcCellWght[cellI]);

            if (mag(tgtVol) > ROOTVSMALL && mag((srcVol-tgtVol)/tgtVol) > 1e-6)
            {
                WarningInFunction
                    << "At cell " << cellI << " cc:"
                    << tgt_.cellCentres()[cellI]
                    << " vol:" << tgtVol
                    << " total overlap volume:" << srcVol
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::correctedCellVolumeWeightMethod::correctedCellVolumeWeightMethod
(
    const polyMesh& src,
    const polyMesh& tgt,
    const dictionary& dict
)
:
    cellVolumeWeightMethod(src, tgt, dict)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

Foam::correctedCellVolumeWeightMethod::~correctedCellVolumeWeightMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::correctedCellVolumeWeightMethod::calculate
(
    labelListList&  srcToTgtAddr,
    scalarListList& srcToTgtWght,
    pointListList&  srcToTgtVec,

    labelListList&  tgtToSrcAddr,
    scalarListList& tgtToSrcWght,
    pointListList&  tgtToSrcVec
)
{
    bool ok = initialise
    (
        srcToTgtAddr,
        srcToTgtWght,
        tgtToSrcAddr,
        tgtToSrcWght
    );

    if (!ok)
    {
        return;
    }

    srcToTgtVec.setSize(srcToTgtAddr.size());
    tgtToSrcVec.setSize(tgtToSrcAddr.size());


    // (potentially) participating source mesh cells
    const labelList srcCellIDs(maskCells());

    // list to keep track of whether src cell can be mapped
    boolList mapFlag(src_.nCells(), false);
    UIndirectList<bool>(mapFlag, srcCellIDs) = true;

    // find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;

    bool startWalk =
        findInitialSeeds
        (
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (startWalk)
    {
        calculateAddressing
        (
            srcToTgtAddr,
            srcToTgtWght,
            srcToTgtVec,
            tgtToSrcAddr,
            tgtToSrcWght,
            tgtToSrcVec,
            srcSeedI,
            tgtSeedI,
            srcCellIDs,
            mapFlag,
            startSeedI
        );
    }
    else
    {
        // if meshes are collocated, after inflating the source mesh bounding
        // box tgt mesh cells may be transferred, but may still not overlap
        // with the source mesh
        return;
    }
}


// ************************************************************************* //
