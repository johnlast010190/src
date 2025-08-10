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

#include "dynamicGIBFvMesh/adjointGIBFvMesh/adjointGIBFvMesh.H"
#include "containers/LinkedLists/user/LIFOStack.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "fvMesh/fvPatches/constraint/empty/emptyFvPatch.H"
#include "fvMesh/fvPatches/constraint/symmetry/symmetryFvPatch.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "fvMesh/fvPatches/constraint/processor/processorFvPatch.H"
#include "dynamicGIBFvMesh/movingGIBTools/parallelIntersectionData/parallelIntersectionData.H"
#include "dynamicGIBFvMesh/movingGIBTools/intersectionData/intersectionData.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adjointGIBFvMesh::constraintPolyPointsOutMotion
(
    vectorField& newPoints
)
{
    const vectorField& baseCf = *baseCf_;
    const pointField& basePoints = *basePoints_;
    const labelListList& cBp = closeBoundaryPoints();

    const labelList& patchPoints = this->boundary()[masterGIB_].
        patch().meshPoints();

    const vectorField& pnf = this->boundary()[masterGIB_].
        patch().pointNormals();

    const boolList& conPoints = concavePoints();

    pointField newPointsTmp = newPoints;

    forAll(patchPoints, pI)
    {
        label gpI = patchPoints[pI];

        forAll(cBp[gpI], fI)
        {
            const label& faceI = cBp[gpI][fI];
            const point& startingPoint = this->points()[gpI];

            vector dis = newPoints[pI] - startingPoint;
            const point& fc = baseCf[faceI];

            pointHit faceInters = this->faces()[faceI].intersection
                (
                    startingPoint,
                    dis,
                    fc,
                    basePoints,
                    intersection::FULL_RAY
                );


            if (faceInters.hit())
            {
                const point& hitP =  faceInters.hitPoint();

                scalar magDis = mag(dis);

                vector vector2 = hitP - startingPoint;
                scalar magDis2 = mag(vector2);

                scalar dir = dis&vector2;


                if (dir>0)
                {
                    if ((magDis2<magDis))
                    {
                        if (!conPoints[gpI])
                        {
                            newPointsTmp[pI] = hitP;
                        }
                        else
                        {
                            newPointsTmp[pI] = basePoints[gpI];
                        }
                    }
                }
                else if (dir == 0)
                {
                    if (startingPoint == basePoints[gpI])
                    {
                        if ((pnf[pI]&dis) > 0)
                        {
                            newPointsTmp[pI] = basePoints[gpI];
                        }
                    }
                    else
                    {
                        newPointsTmp[pI] = hitP;
                    }
                }
                else
                {
                }
               /*
                if
                (
                    (basePoints[gpI] == debPoint1)||
                    (basePoints[gpI] == debPoint2)
                )
                {
                    Info<< "AAA" <<endl;
                    Info<< gpI <<endl;
                    Info<< startingPoint <<endl;
                    Info<< basePoints[gpI] <<endl;
                    Info<< newPointsTmp[pI] <<endl;
                    Info<< "dir: " << dir <<endl;
                    Info<< "magDis2: " << magDis2 <<endl;
                    Info<< "magDis: " << magDis <<endl;
                    Info<< "pnf: " << pnf[pI] <<endl;
                    Info<< "magDis2<magDis: " << (magDis2<magDis) <<endl;
                    Info<< "pnf[pI]&dis: " << (pnf[pI]&dis) <<endl;
                    Info<< "hitP: " << hitP <<endl;
                    Info<< "startingPoint == basePoints[gpI] "
                         << startingPoint <<tab << (basePoints[gpI]) <<tab <<
                    (startingPoint == basePoints[gpI]) <<endl;
                    Info<< "startingPoint - basePoints[gpI] " <<  (startingPoint
                    - basePoints[gpI]) <<endl;
                    Info<< endl;
                }
                */
            }
            else
            {
            }
        }
    }
    newPoints = newPointsTmp;
}



void Foam::adjointGIBFvMesh::nearBoundaryIntersectionsChecking
(
    vectorField& newPoints
)
{
    const vectorField& baseCf = *baseCf_;
    const pointField& basePoints = *basePoints_;
    const boolList& cBp = markedBoundaryPoints();

    const labelList& patchPoints = this->boundary()[masterGIB_].
        patch().meshPoints();

    DynamicList<intersectionData> interData(patchPoints.size());

    this->boundaryMesh().findNeighbProPatchIDs();

    vectorField newPointsCpy(newPoints);
    const pointField& pnf = this->boundary()[masterGIB_].patch().pointNormals();

    forAll(patchPoints, pI)
    {
        label gpI = patchPoints[pI];
        const point& startPoint = this->points()[gpI];
        const point& endPoint = newPointsCpy[pI];
        if (endPoint!=startPoint)
        {
            const vector dis = endPoint - startPoint;

            //- check points close to boundary
            if (cBp[gpI])
            {
                LIFOStack<label> lst;
                boolList checkedFaces(this->faces().size(), false);
                boolList checkedCells(this->cells().size(), false);

                const labelList& pCellsI = this->pointCells()[gpI];
                forAll(pCellsI, cI)
                {
                    const label& gcI = pCellsI[cI];
                    const labelList& cellFace =this->cells()[gcI];
                    forAll(cellFace, fI)
                    {
                        const label& gfI = cellFace[fI];
                        bool exists = false;
                        labelList lstcpy(lst);
                        forAll(lstcpy, fII)
                        {
                            if (lstcpy[fII] == gfI)
                            {
                                exists = true;
                            }
                        }
                        if (!exists)
                        {
                            lst.push(gfI);
                            checkedFaces[gfI] = true;
                        }
                    }
                    checkedCells[gcI] = true;
                }

                bool boundaryHit = false;

                while
                (
                    (lst.size()!=0) && (!boundaryHit)
                )
                {
                    label fII = lst.top();
                    const point& fc = baseCf[fII];
                    pointHit faceInters = this->faces()[fII].intersection
                        (
                            startPoint,
                            dis,
                            fc,
                            basePoints,
                            intersection::FULL_RAY
                        );
                    if (faceInters.hit())
                    {
                        vector dis1 = faceInters.hitPoint() - startPoint;
                        scalar magDis = mag(dis);
                        scalar magDis1 = mag(dis1);

                        if (magDis1>SMALL)
                        {
                            if
                            (
                                (((dis1&dis)>0))&&
                                (magDis1<magDis)
                            )
                            {
                                if (!(fII<this->nInternalFaces()))
                                {
                                    if (false)
                                    {
                                        Info<< "start: " << startPoint << endl;
                                        Info<< "end: " << endPoint << endl;
                                        Info<< "hitPoint: " << faceInters.hitPoint() << endl;
                                        Info<< "fII " << fII << tab << fc << endl;
                                        Info<< "dis1&dis " << (dis1&dis) << endl;
                                        Info<< "magDis1 magDis " << magDis1 << tab << magDis << endl;
                                        Info<< endl;
                                    }
                                    label patchI =
                                        this->boundaryMesh().whichPatch(fII);
                                    if
                                    (
                                        !isA<wedgeFvPatch>(this->boundary()[patchI]) &&
                                        !isA<emptyFvPatch>(this->boundary()[patchI]) &&
                                        !isA<symmetryFvPatch>(this->boundary()[patchI]) &&
                                        !this->boundary()[patchI].coupled()
                                    )
                                    {
                                        boundaryHit = true;
                                        newPointsCpy[pI] = faceInters.hitPoint();
                                    }
                                    else if (this->boundary()[patchI].coupled())
                                    {
                                        const processorPolyPatch& cproPolyPatch =
                                            dynamic_cast<const processorPolyPatch&>
                                            (
                                                this->boundary()[patchI].patch()
                                            );
                                        processorPolyPatch& proPolyPatch =
                                            const_cast<processorPolyPatch&>(cproPolyPatch);
                                        label neiProcID = proPolyPatch.neighbProcNo();
                                        label neiProcPatchID =
                                            proPolyPatch.neighbProPatchID();

                                        intersectionData iD
                                        (
                                            pI,
                                            patchI,
                                            fII-proPolyPatch.start(),
                                            Pstream::myProcNo(),
                                            neiProcID,
                                            neiProcPatchID,
                                            startPoint,
                                            endPoint,
                                            false,
                                            true,
                                            false,
                                            faceInters.hitPoint()
                                        );
                                        interData.append(iD);
                                    }
                                }
                                else
                                {
                                    const label& on = this->owner()[fII];
                                    const label& nb = this->neighbour()[fII];
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
                                        const labelList& cellFace =this->cells()[cellI];
                                        forAll(cellFace, fI)
                                        {
                                            const label& faceI = cellFace[fI];
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
                        else
                        {
                            if ((pnf[pI]&dis)>0)
                            {
                                boundaryHit = true;
                                newPointsCpy[pI] = startPoint;
                            }
                        }
                    }
                    lst.pop();
                }
            }
        }
    }

    interData.shrink();

    parallelIntersectionData pID
    (
        *this,
        basePoints,
        baseCf,
        interData
    );
    pID.modifyPointsPassingBoundary(newPointsCpy);

    syncProcBoundaryPoints(newPoints, newPointsCpy);
}


void Foam::adjointGIBFvMesh::checkConcaveBoundaryPatchPoints
(
    vectorField& newPoints
)
{
    const pointField& basePoints = *basePoints_;
    const labelList& patchPoints = this->boundary()[masterGIB_].
        patch().meshPoints();

    pointField startedpoints = basePoints;
    pointField endPoints = basePoints;
    forAll(patchPoints, pI)
    {
        label gpI = patchPoints[pI];
        endPoints[gpI] = newPoints[pI];
    }

    checkConcaveBoundaryPoints(startedpoints, endPoints);

    forAll(patchPoints, pI)
    {
        label gpI = patchPoints[pI];
        newPoints[pI] = endPoints[gpI];
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
