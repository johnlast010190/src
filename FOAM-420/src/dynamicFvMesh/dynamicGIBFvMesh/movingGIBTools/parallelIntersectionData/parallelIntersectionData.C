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

#include "dynamicGIBFvMesh/movingGIBTools/parallelIntersectionData/parallelIntersectionData.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::parallelIntersectionData::faceCellFaceItersections
(
    intersectionData& iD,
    LIFOStack<label>& lst,
    boolList& checkedFaces,
    boolList& checkedCells
)
{
    while
    (
        (lst.size()!=0) && (!iD.hitWall()) && (!iD.hitProc())
    )
    {
        label fI = lst.top();
        const point& fc = baseCf_[fI];
        const point& spI = iD.pS();
        const point& epI = iD.pE();

        const vector dis = epI - spI;

        pointHit faceInters = mesh_.faces()[fI].intersection
            (
                spI,
                dis,
                fc,
                basePoints_,
                intersection::FULL_RAY
            );
        if (faceInters.hit())
        {
            vector dis1 = faceInters.hitPoint() - spI;
            scalar magDis = mag(dis);
            scalar magDis1 = mag(dis1);
            if
            (
                ((dis1&dis)>0)&&
                (magDis1<magDis)
            )
            {
                if (!(fI<mesh_.nInternalFaces()))
                {
                    label patchI =
                        mesh_.boundaryMesh().whichPatch(fI);
                    if
                    (
                        !isA<wedgeFvPatch>(mesh_.boundary()[patchI]) &&
                        !isA<emptyFvPatch>(mesh_.boundary()[patchI]) &&
                        !isA<symmetryFvPatch>(mesh_.boundary()[patchI]) &&
                        !mesh_.boundary()[patchI].coupled()
                    )
                    {
                        iD.hitP() = faceInters.hitPoint();
                        iD.hitWall() = true;
                        iD.hitProc() = false;
                        iD.hitNone() = false;
                    }
                    else if (mesh_.boundary()[patchI].coupled())
                    {
                        const processorPolyPatch& cproPolyPatch =
                            dynamic_cast<const processorPolyPatch&>
                            (
                                mesh_.boundary()[patchI].patch()
                            );
                        processorPolyPatch& proPolyPatch =
                            const_cast<processorPolyPatch&>(cproPolyPatch);
                        label neiProcID = proPolyPatch.neighbProcNo();
                        label neiProcPatchID =
                            proPolyPatch.neighbProPatchID();

                        iD.hitProc() = true;
                        iD.hitNone() = false;
                        iD.hitWall() = false;

                        //iD.hitP() = faceInters.hitPoint();
                        iD.pAdd() = fI - proPolyPatch.start();
                        //iD.fromProc() = Pstream::myProcNo();
                        iD.toProc() = neiProcID;
                        iD.neiProcPatch() = neiProcPatchID;
                    }
                }
                else
                {
                    const label& on = mesh_.owner()[fI];
                    const label& nb = mesh_.neighbour()[fI];
                    label cellI = -1;
                    if ((checkedCells[on])&&(!checkedCells[nb]))
                    {
                        cellI = nb;
                    }
                    else if (!(checkedCells[on])&&(checkedCells[nb]))
                    {
                        cellI = on;
                    }

                    if (cellI != -1)
                    {
                        checkedCells[cellI] = true;
                        const labelList& cellFace = mesh_.cells()[cellI];
                        forAll(cellFace, fII)
                        {
                            const label& faceI = cellFace[fII];
                            if (!checkedFaces[faceI])
                            {
                                lst.push(faceI);
                                checkedFaces[faceI] = true;
                            }
                        }
                    }
                }
            }
        }
        lst.pop();
        if (lst.size()==0 && (!iD.hitProc()) && (!iD.hitWall()))
        {
            iD.hitNone() = true;
            iD.hitProc() = false;
            iD.hitWall() = false;
        }
    }
}



void Foam::parallelIntersectionData::itersectionChecking()
{
    while (doLoop())
    {
        DynamicList<intersectionData> interData(mesh_.points().size());
        DynamicList<label> addr(mesh_.points().size());
        sendPointsToProcessors(interData, addr);
        //Pout<< "A" << tab << interData.size() <<endl;

        forAll(interData, pI)
        {
            LIFOStack<label> lst;
            boolList checkedFaces(mesh_.faces().size(), false);
            boolList checkedCells(mesh_.cells().size(), false);

            const label& fI = interData[pI].pAdd();
            const label& patchI = interData[pI].neiProcPatch();
            const label startPatchFi = mesh_.boundary()[patchI].start();
            const label gfI = fI + startPatchFi;
            const label& cellI = mesh_.faceOwner()[gfI];
            const labelList& cellFace = mesh_.cells()[cellI];
            checkedCells[cellI] = true;
            checkedFaces[gfI] = true;
            forAll(cellFace, faceI)
            {
                const label& faceII = cellFace[faceI];
                if (!checkedFaces[faceII])
                {
                    checkedFaces[faceII] = true;
                    lst.push(faceII);
                }
            }

            faceCellFaceItersections
            (
                interData[pI],
                lst,
                checkedFaces,
                checkedCells
            );
        }
        //Pout<< "B" <<endl;
        receiveAndUpdateData(interData, addr);
    }
}


void Foam::parallelIntersectionData::sendPointsToProcessors
(
    DynamicList<intersectionData>& interData,
    DynamicList<label>& addr
)
{
    forAll(giDl_, procI)
    {
        List<intersectionData>& giDlI = giDl_[procI];
        forAll(giDlI, pI)
        {
            if
            (
                giDlI[pI].toProc() == Pstream::myProcNo() &&
                giDlI[pI].hitProc() == true
            )
            {
                giDlI[pI].hitProc() = false; // reset flag
                interData.append(giDlI[pI]);
                addr.append(pI);
            }
        }
    }
    interData.shrink();
    addr.shrink();
}


void Foam::parallelIntersectionData::receiveAndUpdateData
(
    const Foam::List<Foam::intersectionData>& interData,
    const Foam::List<label>& addr
)
{
   // Pout<< "After: " << interData <<endl;

    List<List<intersectionData>> ginterData(Pstream::nProcs());
    ginterData[Pstream::myProcNo()] = interData;
    Pstream::allGatherList(ginterData);

    List<List<label>> gaddr(Pstream::nProcs());
    gaddr[Pstream::myProcNo()] = addr;
    Pstream::allGatherList(gaddr);

    forAll(ginterData, procI)
    {
        const List<intersectionData>& ginterDataI = ginterData[procI];
        forAll(ginterDataI, pI)
        {
            const intersectionData& ginterDataII = ginterDataI[pI];
            if
            (
                ginterDataII.fromProc() == Pstream::myProcNo()
            )
            {
                const label& addI = gaddr[procI][pI];
                giDl_[Pstream::myProcNo()][addI] = ginterDataII;
            }
        }
    }
}


bool Foam::parallelIntersectionData::doLoop()
{
    int sumPoints = 0;
    forAll(giDl_[Pstream::myProcNo()], pI)
    {
        if (giDl_[Pstream::myProcNo()][pI].hitProc())
        {
            ++sumPoints;
        }
    }

    reduce(sumPoints, sumOp<int>());

    if (sumPoints>0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::parallelIntersectionData::parallelIntersectionData
(
    const dynamicGIBFvMesh& mesh,
    const pointField& basePoints,
    const vectorField& baseCf,
    const List<intersectionData>& iDl
)
:
    mesh_(mesh),
    basePoints_(basePoints),
    baseCf_(baseCf),
    iDl_(iDl),
    giDl_(Pstream::nProcs())
{
    giDl_[Pstream::myProcNo()] = iDl_;
    Pstream::allGatherList(giDl_);
   // Info<< giDl_ <<endl;
}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::parallelIntersectionData::modifyPointsPassingBoundary
(
    pointField& newPoints
)
{
    itersectionChecking();
    forAll(giDl_[Pstream::myProcNo()], pI)
    {
        if (giDl_[Pstream::myProcNo()][pI].hitWall())
        {
            label gpI = giDl_[Pstream::myProcNo()][pI].pp();
            newPoints[gpI] = giDl_[Pstream::myProcNo()][pI].hitP();
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
