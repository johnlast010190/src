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
    (c) 2017 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "dynamicGIBFvMesh/dynamicGIBFvMesh/dynamicGIBFvMesh.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "meshes/meshTools/simpleVTKWriter.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::dynamicGIBFvMesh::flipCells
(
    boolList& interPoints,
    bool regionType
) const
{
    interPoints = false;
    boolList markInterFace = boolList(this->faces().size(), false);

    const labelList& fll = fl();
    const labelList& reg = cRegion();

    //- mark interface points
    forAll(fll, fI)
    {
        const face& faceI = this->faces()[fll[fI]];
        markInterFace[fll[fI]] = true;
        forAll(faceI, pfI)
        {
            const label& pI = faceI[pfI];
            interPoints[pI] = true;
        }
    }
    syncTools::syncFaceList(*this, markInterFace, orEqOp<bool>());

    syncTools::syncPointList
    (
        *this,
        interPoints,
        plusEqOp<bool>(),
        true
    );

    //- mark all the cells in which all the points are interface points
    boolList flipCell = checkFlipCell(interPoints);

    forAll(flipCell, cI)
    {
        if (flipCell[cI])
        {
            if
            (
                (reg[cI]==0 && regionType==false)||
                (reg[cI]==1 && regionType==true)
            )
            {
//                Pout<< cI <<  tab << regionType << tab<< reg[cI] <<endl;
//                flipCell[cI] = false;
            }
        }
    }

    labelList cellFlipN = cellFlippingNumber(flipCell, markInterFace);

    if (!pyrPrismFlip_)
    {
        resetFlipCellsWithFewInterfaces(flipCell, cellFlipN);
    }

    labelList faceFlipType = faceFlippingType(flipCell, markInterFace);
/*
    findFlippingFaces
    (
        flipCell,
        markInterFace,
        interPoints,
        cellFlipN,
        faceFlipType
    );
    */

//    resetFlipHexCellsWith3Interfaces(flipCell, cellFlipN, faceFlipType);


//----------------------//
    boolList cellId = boolList(this->cells().size(), false);
    boolList faceId = boolList(this->faces().size(), false);
    forAll(reg, cI)
    {
        cellId[reg[cI]] = true;
    }
    forAll(faceId, fI)
    {
        if (fI<nInternalFaces())
        {
            faceId[fI] = !cellId[this->owner()[fI]]^cellId[this->neighbour()[fI]];
        }
        else
        {
            faceId[fI] = cellId[this->faceOwner()[fI]];
        }
    }
    syncTools::syncFaceList(*this, faceId, xxorEqOp());





//---------------------//


    //- snap - unsnap cell based on faces
    forAll(markInterFace, fI)
    {
        if (faceFlipType[fI] == 1)
        {
            markInterFace[fI] = true;
        }
        else if (faceFlipType[fI] == -1)
        {
            markInterFace[fI] = false;
        }
        else if (faceFlipType[fI] == 2)
        {
            if (faceId[fI])
            {
                markInterFace[fI] = false;
            }
            else
            {
                markInterFace[fI] = true;
            }
        }
        else if (faceFlipType[fI] == -2)
        {
            markInterFace[fI] = true;
        }
        else if (faceFlipType[fI] != 0)
        {
        }
    }

    syncTools::syncFaceList(*this, markInterFace, orEqOp<bool>());

    singleCellSuddenPopCellRemoval(markInterFace);

    interPoints = false;
    DynamicList<label> dfl(this->faces().size());
    forAll(this->faces(), fI)
    {
        if (markInterFace[fI])
        {
            if (fI<this->nInternalFaces())
            {
                forAll(this->faces()[fI], pfI)
                {
                    const label& pI = this->faces()[fI][pfI];
                    interPoints[pI] = true;
                }
                dfl.append(fI);
            }
            else
            {
                label patchI = this->boundaryMesh().whichPatch(fI);
                const polyPatch& pp = this->boundary()[patchI].patch();
                if (pp.coupled() || includeWalls())
                {
                    forAll(this->faces()[fI], pfI)
                    {
                        const label& pI = this->faces()[fI][pfI];
                        interPoints[pI] = true;
                    }
                    dfl.append(fI);
                }
            }
        }
    }

    syncTools::syncPointList
    (
        *this,
        interPoints,
        plusEqOp<bool>(),
        true
    );

    writeProblematicCells(interPoints);

    dfl.shrink();


    if (false)
    {
        const pointField& basePoints = *basePoints_;
        simpleVTKWriter
        (
            this->faces(),
            dfl,
            basePoints
        ).write
        (
            "fl0Flipped"+this->time().timeName()+".vtk"
        );
    }

    delete flPtr_;
    flPtr_ = new labelList (dfl);

    deleteDemandDrivenData(cRegionPtr_);
    deleteDemandDrivenData(fmPtr_);

}

void Foam::dynamicGIBFvMesh::singleCellSuddenPopCellRemoval
(
    boolList& markInterface
) const
{
    const boolList& markInterface0 = faceIndicator0();

    label nInternalFaces = this->nInternalFaces();
    forAll(this->cells(), cI)
    {
        const labelList& cellFaces(this->cells()[cI]);
        bool unsnap = true;
        forAll(cellFaces, cfI)
        {
            const label& fI = cellFaces[cfI];
            if
            (
                markInterface0[fI]
            ||  !markInterface[fI]
            ||  !(fI<nInternalFaces)
            )
            {
                unsnap = false;
            }
        }
        if (unsnap)
        {
            forAll(cellFaces, cfI)
            {
                const label& fI = cellFaces[cfI];
                markInterface[fI] = false;
            }
        }
    }
}


Foam::boolList Foam::dynamicGIBFvMesh::checkFlipCell
(
    const boolList& interPoints
) const
{
    const labelList& fll = fl();
    boolList flipCell = boolList(this->cells().size(), false);

    forAll(fll, fllI)
    {
        const label& fI = fll[fllI];
        bool flip = true;
        if (fI<this->nInternalFaces())
        {
            const label& on = this->owner()[fI];
            const label& nb = this->neighbour()[fI];
            const labelList& cPOn = this->cellPoints()[on];
            const labelList& cPNb = this->cellPoints()[nb];
            forAll(cPOn, cpI)
            {
                const label& pI = cPOn[cpI];
                if (!interPoints[pI])
                {
                    flip = false;
                }
            }
            flipCell[on] = flip;

            flip = true;
            forAll(cPNb, cpI)
            {
                const label& pI = cPNb[cpI];
                if (!interPoints[pI])
                {
                   flip = false;
                }
            }
            flipCell[nb] = flip;
        }
        else
        {
            flip = true;
            const label& on = this->faceOwner()[fI];
            const labelList& cPOn = this->cellPoints()[on];
            forAll(cPOn, cpI)
            {
                const label& pI = cPOn[cpI];
                if (!interPoints[pI])
                {
                    flip = false;
                }
            }
            flipCell[on] = flip;
        }
    }
    return flipCell;
}


Foam::labelList Foam::dynamicGIBFvMesh::faceFlippingType
(
    const boolList& flipCell,
    const boolList& markInterFace
) const
{
    labelList faceFlipType = labelList(this->faces().size(), 0);
    forAll(flipCell, cI)
    {
        if (flipCell[cI])
        {
            const labelList& cFaces = this->cells()[cI];
            forAll(cFaces, cFI)
            {
                const label& fI = cFaces[cFI];
                if (fI<this->nInternalFaces())
                {
                    if (!markInterFace[fI])
                    {
                        faceFlipType[fI] += 1;
                    }
                    else
                    {
                        faceFlipType[fI] -= 1;
                    }
                }
                else
                {
                    label patchI = this->boundaryMesh().whichPatch(fI);
                    if
                    (
                        !isA<wedgeFvPatch>(this->boundary()[patchI]) &&
                        !isA<emptyFvPatch>(this->boundary()[patchI]) &&
                        !isA<symmetryFvPatch>(this->boundary()[patchI])
                    )
                    {
                        if (!markInterFace[fI])
                        {
                            faceFlipType[fI] += 1;
                        }
                        else
                        {
                            faceFlipType[fI] -= 1;
                        }
                    }
                }
            }
        }
    }

    syncTools::syncFaceList
    (
        *this,
        faceFlipType,
        plusEqOp<label>()
    );
    return faceFlipType;
}


Foam::labelList Foam::dynamicGIBFvMesh::cellFlippingNumber
(
    const boolList& flipCell,
    const boolList& markInterFace
) const
{
    labelList cellFlippingN(labelList(this->cells().size(), 0));
    forAll(flipCell, cI)
    {
        if (flipCell[cI])
        {
            const labelList& cFaces = this->cells()[cI];
            forAll(cFaces, cFI)
            {
                const label& cFacesI = cFaces[cFI];
                if (markInterFace[cFacesI])
                {
                    cellFlippingN[cI] += 1;
                }
            }
        }
    }
    return cellFlippingN;
}


void Foam::dynamicGIBFvMesh::findFlippingFaces
(
    const boolList& flipCell,
    boolList& markInterFace,
    const boolList& interPoints,
    const labelList& cellFlipN,
    labelList& faceFlipType
) const
{
    forAll(faceFlipType, fI)
    {
        if (faceFlipType[fI] == -2)
        {
            if (fI<this->nInternalFaces())
            {
                const label& on = this->owner()[fI];
                const label& nb = this->neighbour()[fI];
                const label& cellFlipNon = cellFlipN[on];
                const label& cellFlipNnb = cellFlipN[nb];

                if
                (
                    (this->cellPoints()[on].size() == 8) &&
                    (this->cellPoints()[nb].size() == 8)
                )
                {
                    if ((cellFlipNnb == 3)&&(cellFlipNon == 4))
                    {
                        const labelList& cellFaces(this->cells()[on]);
                        Pout<< "A: " << on << tab << nb <<endl;
                        forAll(cellFaces, cfI)
                        {
                            const label& cellFacesI = cellFaces[cfI];
                            markInterFace[cellFacesI] =
                                !markInterFace[cellFacesI];
                            faceFlipType[cellFacesI] = 0;
                        }
                        const labelList& cellFacesnb(this->cells()[nb]);
                        forAll(cellFacesnb, cfI)
                        {
                            const label& cellFacesI = cellFacesnb[cfI];
                            faceFlipType[cellFacesI] = 0;
                        }
                    }


                    if ((cellFlipNnb == 4)&&(cellFlipNon == 3))
                    {
                        const labelList& cellFaces(this->cells()[nb]);
                        Pout<< "B: " << on << tab << nb <<endl;
                        forAll(cellFaces, cfI)
                        {
                            const label& cellFacesI = cellFaces[cfI];
                            markInterFace[cellFacesI] =
                                !markInterFace[cellFacesI];
                            faceFlipType[cellFacesI] = 0;
                        }
                        const labelList& cellFaceson(this->cells()[on]);
                        forAll(cellFaceson, cfI)
                        {
                            const label& cellFacesI = cellFaceson[cfI];
                            faceFlipType[cellFacesI] = 0;
                        }
                    }
                }
            }
        }
    }
}


void Foam::dynamicGIBFvMesh::resetFlipCellsWithFewInterfaces
(
    boolList& flipCell,
    const labelList& cellFlipN
) const
{
    forAll(flipCell, cI)
    {
        if (flipCell[cI])
        {
            if (this->cells()[cI].size() == 6)
            {
//                if (cellFlipN[cI] < 3)
                {
                    flipCell[cI] = false;
                }
            }
            else if (this->cells()[cI].size() == 5)
            {
//                if (cellFlipN[cI] < 2)
                {
                    flipCell[cI] = false;
                }
            }
        }
    }
}


void Foam::dynamicGIBFvMesh::resetFlipHexCellsWith3Interfaces
(
    boolList& flipCell,
    const labelList& cellFlipN,
    labelList& faceFlipType
) const
{
    forAll(flipCell, cI)
    {
        if (flipCell[cI])
        {
            if
            (
                (this->cells()[cI].size() == 6) &&
                (this->cellPoints()[cI].size() == 8)
            )
            {
                if (cellFlipN[cI] == 3)
                {
                    bool multipleflipCFould = false;
                    const labelList& cellFaces(this->cells()[cI]);
                    forAll(cellFaces, cfI)
                    {
                        const label& cellFacesI = cellFaces[cfI];
                        if
                        (
                            (faceFlipType[cellFacesI] == 2)
                            || (faceFlipType[cellFacesI] == -2)
                        )
                        {
                            multipleflipCFould = true;
                        }
                    }
                    if (!multipleflipCFould)
                    {
                        flipCell[cI] = false;
                        forAll(cellFaces, cfI)
                        {
                            const label& cellFacesI = cellFaces[cfI];
                            faceFlipType[cellFacesI] = 0;
                        }
                    }
                }
            }
            if (this->cells()[cI].size() == 5)
            {
                if (cellFlipN[cI] == 2)
                {
                    bool multipleflipCFould = false;
                    const labelList& cellFaces(this->cells()[cI]);
                    forAll(cellFaces, cfI)
                    {
                        const label& cellFacesI = cellFaces[cfI];
                        if
                        (
                            (faceFlipType[cellFacesI] == 2)
                            || (faceFlipType[cellFacesI] == -2)
                        )
                        {
                            multipleflipCFould = true;
                        }
                    }
                    if (!multipleflipCFould)
                    {
                        flipCell[cI] = false;
                        forAll(cellFaces, cfI)
                        {
                            const label& cellFacesI = cellFaces[cfI];
                            faceFlipType[cellFacesI] = 0;
                        }
                    }
                }
            }
        }
    }
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
