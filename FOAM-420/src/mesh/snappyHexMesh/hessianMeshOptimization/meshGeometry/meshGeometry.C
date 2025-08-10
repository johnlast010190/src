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
    (c) 2013 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/
#include "hessianMeshOptimization/meshGeometry/meshGeometry.H"
#include "hessianMeshOptimization/meshGeometry/faceTriangle.H"
// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::meshGeometry::updateAreas()
{
    const labelList& fs = mesh_.pointFaces()[pLabel_];
    const pointField& p = state_;
    int end = fs.size();
    for (int facei=0; facei<end; ++facei)
    {
        const label& fL = fs[facei];
        const face& fa = mesh_.faces()[fL];
        const labelList& f = fa;
        const label& nPoints = fa.size();

        if ((nPoints-3)== 0)
        {
            fCtrs_[fL] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
            fAreas_[fL] = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));

            const label& ownl = mesh_.faceOwner()[fL];
            const tensor& areaT = metric_.extProductTensor()[ownl];

            magfAreas_[fL][0] = mag(areaT&fAreas_[fL]);

            if (fL < mesh_.nInternalFaces())
            {
                const label& neil = mesh_.faceNeighbour()[fL];
                const tensor& areaT = metric_.extProductTensor()[neil];

                magfAreas_[fL][1] = mag(areaT&fAreas_[fL]);
            }
        }
        else
        {
            point fCentre = p[f[0]];
            for (int pi = 1; pi < nPoints; pi++)
            {
                fCentre += p[f[pi]];
            }

            fCentre /= nPoints;

            vector sumN(vector::zero);
            scalar sumA(0.0);
            scalar sumA_Own(0.0);
            scalar sumA_Nei(0.0);
            vector sumAc(vector::zero);

            for (int pi = 0; pi < nPoints; ++pi)
            {
                const point& nextPoint = p[f[(pi + 1) % nPoints]];
                vector c = p[f[pi]] + nextPoint + fCentre;
                vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);

                const label& ownl = mesh_.faceOwner()[fL];
                const tensor& areaT = metric_.extProductTensor()[ownl];

                sumA_Own += mag(areaT&n);

                if (fL < mesh_.nInternalFaces())
                {
                    const label& neil = mesh_.faceNeighbour()[fL];
                    const tensor& areaT = metric_.extProductTensor()[neil];

                    sumA_Nei += mag(areaT&n);
                }

                scalar a = mag(n);

                sumN += n;
                sumA += a;
                sumAc += a*c;
            }

            if (sumA < ROOTVSMALL)
            {
                fCtrs_[fL] = fCentre;
                fAreas_[fL] = vector::zero;
                magfAreas_[fL] = 0;
            }
            else
            {
                fCtrs_[fL] = (1.0/3.0)*sumAc/sumA;
                fAreas_[fL] = 0.5*sumN;
                magfAreas_[fL][0] = 0.5*sumA_Own;
                magfAreas_[fL][1] = 0.5*sumA_Nei;
            }
        }
    }
}

void Foam::meshGeometry::updateVolumes()
{
    const labelList& cells = mesh_.pointCells()[pLabel_];
    const tensorField& tT = metric_.tTensor();
    forAll(cells, cI)
    {
        const label& cL = cells[cI];
        cellCtrs_[cL] = vector::zero;
        cellVols_[cL] = 0;
    }

    const label& nCells = cells.size();

    vectorField cEst(nCells, vector::zero);
    labelField nCellFaces(nCells, 0);

    forAll(cells, cI)
    {
        const labelList& cFaces = mesh_.cells()[cells[cI]];
        forAll(cFaces, f)
        {
            const label& fL = cFaces[f];
            cEst[cI] += fCtrs_[fL];
            nCellFaces[cI]++;
        }
    }

    forAll(cells, celli)
    {
        cEst[celli] /= nCellFaces[celli];
    }

    const labelList& own = mesh_.faceOwner();

    forAll(cells, cI)
    {
        const label& cL = cells[cI];
        const labelList& cFaces = mesh_.cells()[cL];
        const tensor& areaT = metric_.extProductTensor()[cL];
        forAll(cFaces, f)
        {
            const label& fL = cFaces[f];
            scalar pyr3Vol;
            if (own[fL]==cL)
            {
                // Calculate 3*face-pyramid volume
                pyr3Vol =
                max((areaT&fAreas_[fL]) & (tT[cL]&(fCtrs_[fL] - cEst[cI])), VSMALL);
            }
            else
            {
                   pyr3Vol =
                   max((areaT&fAreas_[fL]) & (tT[cL]&(cEst[cI] - fCtrs_[fL])), VSMALL);
            }

               // Calculate face-pyramid centre
               vector pc = tT[cL]&((3.0/4.0)*fCtrs_[fL] + (1.0/4.0)*cEst[cI]);

               // Accumulate volume-weighted face-pyramid centre
               cellCtrs_[cL] += pyr3Vol*pc;

               // Accumulate face-pyramid volume
               cellVols_[cL] += pyr3Vol;
        }
    }

    forAll(cells, cI)
    {
        const label& cL = cells[cI];
        cellCtrs_[cL] /= cellVols_[cL];
        cellVols_[cL] *= (1.0/3.0);
    }
}

void Foam::meshGeometry::updateCellSurface()
{
    const labelList& cells = mesh_.pointCells()[pLabel_];
    forAll(cells, cI)
    {
        const label& cL = cells[cI];
        cellSurface_[cL] = 0.0;
    }
    forAll(cells, cI)
    {
        const label& cL = cells[cI];
        const labelList& cFaces = mesh_.cells()[cells[cI]];
        forAll(cFaces, f)
        {
            const label& fL = cFaces[f];
            if (mesh_.faceOwner()[fL]==cL)
            {
                cellSurface_[cL] += magfAreas_[fL][0];
            }
            else
            {
                cellSurface_[cL] += magfAreas_[fL][1];
            }
        }
    }
}

void Foam::meshGeometry::initializeCellSurface()
{
    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    cellSurface_ *= 0;
    forAll(mesh_.faces(), fI)
    {
        if (fI<mesh_.nInternalFaces())
        {
            const tensor& areaT = metric_.extProductTensor()[own[fI]];

            magfAreas_[fI][0] = mag(areaT&fAreas_[fI]);

            const tensor& areaTNei = metric_.extProductTensor()[nei[fI]];

            magfAreas_[fI][1] = mag(areaTNei&fAreas_[fI]);
        }
        else
        {
            const tensor& areaT = metric_.extProductTensor()[own[fI]];
            magfAreas_[fI][0] = mag(areaT&fAreas_[fI]);
        }
    }
    forAll(mesh_.faces(), fI)
    {
        if (fI<mesh_.nInternalFaces())
        {
            const label& neil = nei[fI];
            const label& ownl = own[fI];
            cellSurface_[neil] += magfAreas_[fI][0];
            cellSurface_[ownl] += magfAreas_[fI][1];
        }
        else
        {
            const label& ownl = own[fI];
            cellSurface_[ownl] += magfAreas_[fI][0];
        }
    }
}

void Foam::meshGeometry::updateGeometry(const label& pI)
{
    pLabel_ = pI;
    updateAreas();
    updateVolumes();
    updateCellSurface();
}

void Foam::meshGeometry::calculateDerivatives(const label& pI)
{
    pLabel_ = pI;
    calculateCellDerivatives();
}

void Foam::meshGeometry::calculateFaceCenterDerivatives()
{
    faceCentreDerivatives_.setSize(mesh_.faces().size(), tensor::zero);
    forAll(mesh_.faces(), fI)
    {
        const labelList& fPoints = mesh_.faces()[fI];
        const label& nPoints = fPoints.size();
        scalar invSize = 1./nPoints;
        faceCentreDerivatives_[fI].xx() = invSize;
        faceCentreDerivatives_[fI].yy() = invSize;
        faceCentreDerivatives_[fI].zz() = invSize;
    }
}

void Foam::meshGeometry::calculateCellDerivatives()
{
    const labelList& faces = mesh_.pointFaces()[pLabel_];
    const labelList& cells = mesh_.pointCells()[pLabel_];

    const labelList& own = mesh_.faceOwner();
    const labelList& nei = mesh_.faceNeighbour();
    const tensorField& tT = metric_.tTensor();
    const tensorField& areaT = metric_.extProductTensor();
    forAll(cellSurfaceDerivatives_, cI)
    {
        cellSurfaceDerivatives_[cI] *= 0;
        cellVolumeDerivatives_[cI] *= 0;
        hessian_[cI] *= 0;
    }

    hessian_.setSize(cells.size(), tensor::zero);
    cellSurfaceDerivatives_.resize(cells.size(), vector::zero);
    cellVolumeDerivatives_.resize(cells.size(), vector::zero);

    label owner = 1;
    label neigh = 0;
    label bound = -1;
    label pI;
    forAll(faces, I)
    {
        const label& fI = faces[I];
        const labelList& fPoints = mesh_.faces()[fI];

        forAll(fPoints, i)
        {
            if (fPoints[i]==pLabel_)
            {
                pI = i;
                break;
            }
        }

        if (fI<mesh_.nInternalFaces())
        {
            label ownl = own[fI];
            label neil = nei[fI];
            int ptrOwn = getTopologyIndex(ownl);
            int ptrNei = getTopologyIndex(neil);

            //TODO refactor the new surface derivative method
            faceTriangle tri(*this, fI, pI);
            ::pyramid pyrNei(*this, tri, neigh, metric_);
            ::pyramid pyrOwn(*this, tri, owner, metric_);

            label nPoints = mesh_.faces()[fI].size();
            label nextLabel = (pI + 1) % nPoints;
            label previousLabel = (nextLabel + nPoints - 2) % nPoints;

            point neP = tT[ownl]&state_[fPoints[nextLabel]];
            point prP = tT[ownl]&state_[fPoints[previousLabel]];
            vector transfAreas = areaT[ownl]&fAreas_[fI];

            vector contribution;

            for (int i = 0; i<2; i++)
            {
                label cL = 0;
                if (i==1)
                {
                    neP = tT[neil]&state_[fPoints[nextLabel]];
                    prP = tT[neil]&state_[fPoints[previousLabel]];
                    transfAreas = areaT[neil]&fAreas_[fI];
                    cL = 1;
                }

                contribution.x() = 0.5*(transfAreas.y()*(prP.z()-neP.z())
                                        + transfAreas.z()*(neP.y()-prP.y()));

                contribution.y() = 0.5*(transfAreas.x()*(neP.z()-prP.z())
                                        + transfAreas.z()*(prP.x()-neP.x()));

                contribution.z() = 0.5*(transfAreas.x()*(prP.y()-neP.y())
                                        + transfAreas.y()*(neP.x()-prP.x()));

                scalar magfAreasCL = magfAreas_[fI][cL] + VSMALL;

                contribution /= magfAreasCL;

                tensor t = tensor::zero;

                t.xx() = -contribution.x()*contribution.x()/magfAreasCL
                    + 0.25*((prP.z()-neP.z())*(prP.z()-neP.z())
                    +(neP.y()-prP.y())*(neP.y()-prP.y()))/magfAreasCL;

                t.yy() = -contribution.y()*contribution.y()/magfAreasCL
                    + 0.25*((prP.z()-neP.z())*(prP.z()-neP.z())
                    +(neP.x()-prP.x())*(neP.x()-prP.x()))/magfAreasCL;

                t.zz() = -contribution.z()*contribution.z()/magfAreasCL
                    + 0.25*((prP.y()-neP.y())*(prP.y()-neP.y())
                    +(neP.x()-prP.x())*(neP.x()-prP.x()))/magfAreasCL;

                t.xy() = -contribution.x()*contribution.y()/magfAreasCL
                    + 0.25*((neP.y()-prP.y())*(prP.x()-neP.x()))/magfAreasCL;

                t.xz() = -contribution.x()*contribution.z()/magfAreasCL
                    + 0.25*((neP.x()-prP.x())*(prP.z()-neP.z()))/magfAreasCL;

                t.yz() = -contribution.y()*contribution.z()/magfAreasCL
                    + 0.25*((neP.z()-prP.z())*(prP.y()-neP.y()))/magfAreasCL;

                if (i==0)
                {
                    cellSurfaceDerivatives_[ptrOwn] += contribution;
                    hessian_[ptrOwn] += t;
                }
                else
                {
                    cellSurfaceDerivatives_[ptrNei] += contribution;
                    hessian_[ptrNei] += t;
                }
            }

            contribution = pyrOwn.contribution();
            cellVolumeDerivatives_[ptrOwn] += contribution;

            contribution = pyrNei.contribution();
            cellVolumeDerivatives_[ptrNei] += contribution;

            while (tri.isANewPoint())
            {
                tri.next();
                contribution = pyrOwn.contribution();
                cellVolumeDerivatives_[ptrOwn] +=
                    contribution;

                contribution = pyrNei.contribution();
                cellVolumeDerivatives_[ptrNei] +=
                    contribution;
            }
        }
        else
        {
            label whichCell = own[fI];
            int ptrc = getTopologyIndex(whichCell);

            faceTriangle tri(*this, fI, pI);
            ::pyramid pyrOwn(*this, tri, bound, metric_);

            label nPoints = mesh_.faces()[fI].size();
            label nextLabel = (pI + 1) % nPoints;
            label previousLabel = (nextLabel + nPoints - 2) % nPoints;

            point neP = tT[whichCell]&state_[fPoints[nextLabel]];
            point prP = tT[whichCell]&state_[fPoints[previousLabel]];
            vector transfAreas = areaT[whichCell]&fAreas_[fI];

            vector contribution;

            contribution.x() = 0.5*(transfAreas.y()*(prP.z()-neP.z())
                    + transfAreas.z()*(neP.y()-prP.y()));

            contribution.y() = 0.5*(transfAreas.x()*(neP.z()-prP.z())
                    + transfAreas.z()*(prP.x()-neP.x()));

            contribution.z() = 0.5*(transfAreas.x()*(prP.y()-neP.y())
                    + transfAreas.y()*(neP.x()-prP.x()));

            scalar magfAreas0 = magfAreas_[fI][0] + VSMALL;

            contribution /= magfAreas0;

            cellSurfaceDerivatives_[ptrc] += contribution;

            hessian_[ptrc].xx() += -contribution.x()*contribution.x()/magfAreas0
                + 0.25*((prP.z()-neP.z())*(prP.z()-neP.z())
                +(neP.y()-prP.y())*(neP.y()-prP.y()))/magfAreas0;

            hessian_[ptrc].yy() += -contribution.y()*contribution.y()/magfAreas0
                + 0.25*((prP.z()-neP.z())*(prP.z()-neP.z())
                +(neP.x()-prP.x())*(neP.x()-prP.x()))/magfAreas0;

            hessian_[ptrc].zz() += -contribution.z()*contribution.z()/magfAreas0
                + 0.25*((prP.y()-neP.y())*(prP.y()-neP.y())
                +(neP.x()-prP.x())*(neP.x()-prP.x()))/magfAreas0;

            hessian_[ptrc].xy() += -contribution.x()*contribution.y()/magfAreas0
                + 0.25*((neP.y()-prP.y())*(prP.x()-neP.x()))/magfAreas0;

            hessian_[ptrc].xz() += -contribution.x()*contribution.z()/magfAreas0
                + 0.25*((neP.x()-prP.x())*(prP.z()-neP.z()))/magfAreas0;

            hessian_[ptrc].yz() += -contribution.y()*contribution.z()/magfAreas0
                + 0.25*((neP.z()-prP.z())*(prP.y()-neP.y()))/magfAreas0;

            contribution = pyrOwn.contribution();
            cellVolumeDerivatives_[ptrc] += contribution;
            while (tri.isANewPoint())
            {
                tri.next();
                contribution = pyrOwn.contribution();
                   cellVolumeDerivatives_[ptrc] += contribution;
            }
        }
    }

    forAll(cells, cI)
    {
        cellVolumeDerivatives_[cI] /= 6.;
    }
}

int Foam::meshGeometry::getTopologyIndex
(
    const label& cellLabel
)
{
   int index = -1;
   const label& nPointCells = mesh_.pointCells()[pLabel_].size();
   for (int i=0; i < nPointCells; i++)
   {
       if (mesh_.pointCells()[pLabel_][i] == cellLabel)
       {
           index = i;
       }
   }
   return index;
}

void Foam::meshGeometry::numberOfCellsPerPoint()
{
    forAll(mesh_.points(), pI)
    {
        const labelList& cPoints = mesh_.pointCells()[pI];
        pointCellSize_[pI] = cPoints.size();
    }

    syncTools::syncPointList
    (
        mesh_,
        pointCellSize_,
        plusEqOp<label>(),
        label(0)
    );
}

void Foam::meshGeometry::findBoundaryPoints()
{
    for (int i = mesh_.nInternalFaces(); i< mesh_.faces().size(); i++)
    {
        const face& f = mesh_.faces()[i];
        forAll(f, pI)
        {
            boundaryPoint_[f[pI]] = true;
        }
    }
}

void Foam::meshGeometry::findBoundaryProcessorPoints()
{
    boolList isAprocPoint(mesh_.points().size(), false);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const labelList& patchPoints = pp.meshPoints();

        if (pp.coupled())
        {
            forAll(patchPoints, pI)
            {
                const label& pL = patchPoints[pI];
                isAprocPoint[pL] = true;
            }
        }
    }

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const labelList& patchPoints = pp.meshPoints();
        if (!pp.coupled())
        {
            forAll(patchPoints, pI)
            {
                const label& pL = patchPoints[pI];
                if (isAprocPoint[pL])
                {
                    boundaryProcessorPoint_[pL] = true;
                }
            }
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        boundaryProcessorPoint_,
        maxEqOp<bool>(),
        false
    );

    label count = 0;
    forAll(mesh_.points(), pI)
    {
        if (isAprocPoint[pI])
        {
            count++;
        }
    }

    processorPoints_.setSize(count);
    count = 0;
    forAll(mesh_.points(), pI)
    {
        if (isAprocPoint[pI])
        {
            processorPoints_[count] = pI;
            count++;
        }
    }
}

void Foam::meshGeometry::splitProcessorPoints()
{
    const labelList& sharedLabels = mesh_.globalData().sharedPointLabels();
    boolList isASharedPoint(mesh_.points().size(), false);
    forAll(sharedLabels, pI)
    {
        isASharedPoint[sharedLabels[pI]] = true;
    }

    labelList colorTag(processorPoints_.size(), -1);
    bool colored = false;
    label color = -1;
    label counter = 0;
    while (!colored)
    {
        color++;
        forAll(processorPoints_, I)
        {
            bool valid = true;
            const label& gP = processorPoints_[I];
            if (isASharedPoint[gP] && colorTag[I]==-1)
            {
                const labelList& pointFaces = mesh_.pointFaces()[gP];
                forAll(pointFaces, fI)
                {
                    const label& fL = pointFaces[fI];
                    if (fL >= mesh_.nInternalFaces())
                    {
                        const face& fc = mesh_.faces()[fL];
                        forAll(fc, pI)
                        {
                            const label& meshP = fc[pI];
                            if (meshP!=gP)
                            {
                                label localLabel = -1;
                                forAll(processorPoints_, ppI)
                                {
                                    if (processorPoints_[ppI]==meshP)
                                    {
                                        localLabel = ppI;
                                    }
                                }
                                if (localLabel != -1 && colorTag[localLabel]==color)
                                {
                                    valid = false;
                                }
                            }
                        }
                    }
                }
                if (valid)
                {
                    colorTag[I] = color;
                    counter++;
                }
            }
        }
        if (counter == sharedLabels.size())
        {
            colored = true;
        }
    }

    syncTools::syncPointList
    (
        mesh_,
        processorPoints_,
        colorTag,
        eqOp<label>(),
        label(-1)
    );
    colored = false;
    color = -1;

    while (!colored)
    {
        color++;
        forAll(processorPoints_, I)
        {
            bool valid = true;
            const label& gP = processorPoints_[I];
            if (colorTag[I]==-1)
            {
                const labelList& pointFaces = mesh_.pointFaces()[gP];
                forAll(pointFaces, fI)
                {
                    const label& fL = pointFaces[fI];
                    if (fL >= mesh_.nInternalFaces())
                    {
                        const face& fc = mesh_.faces()[fL];
                        forAll(fc, pI)
                        {
                            const label& meshP = fc[pI];
                            if (meshP!=gP)
                            {
                                label localLabel = -1;
                                forAll(processorPoints_, ppI)
                                {
                                    if (processorPoints_[ppI]==meshP)
                                    {
                                        localLabel = ppI;
                                    }
                                }
                                if (localLabel != -1 && colorTag[localLabel]==color)
                                {
                                    valid = false;
                                }
                            }
                        }
                    }
                }
                if (valid)
                {
                    colorTag[I] = color;
                    counter++;
                }
            }
        }
        if (counter == processorPoints_.size())
        {
            colored = true;
        }
    }
    syncTools::syncPointList
    (
        mesh_,
        processorPoints_,
        colorTag,
        eqOp<label>(),
        label(-1)
    );

    label maxColor = 0;
    if (colorTag.size() > 0)
    {
        maxColor = max(colorTag)+1;
    }

    reduce(maxColor, maxOp<label>());
    Info<<"Max Number of Synchronizations "<<maxColor<<endl;
    procPointsLists_.setSize(maxColor);
    labelList procListCounters(maxColor, 0);
    forAll(processorPoints_, pI)
    {
        const label& color = colorTag[pI];
        procListCounters[color]++;
    }
    forAll(procPointsLists_, listI)
    {
        procPointsLists_[listI].setSize(procListCounters[listI]);
    }

    procListCounters = 0;
    forAll(processorPoints_, pI)
    {
        const label& cl = colorTag[pI];
        procPointsLists_[cl][procListCounters[cl]] = processorPoints_[pI];
        procListCounters[cl]++;
    }
}

void Foam::meshGeometry::calculateCellTypes()
{
    forAll(mesh_.cells(), cI)
    {
        label size = mesh_.cells()[cI].size();
        if (size == 4)
        {
//            cellTypeCoef_[cI] = 2*1.732051;
            cellTypeCoef_[cI] = 2;
        }
        else if (size == 5)
        {
//            cellTypeCoef_[cI] = 2*(1+1./3);
            cellTypeCoef_[cI] = 2;
        }
        else
        {
            cellTypeCoef_[cI] = 1;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::meshGeometry::meshGeometry
(
    const fvMesh& mesh,
    const pointField& state,
    const dictionary& dict
)
:
    mesh_(mesh),
    state_(state),
    fCtrs_(mesh.faceCentres()),
    fAreas_(mesh.faceAreas()),
    magfAreas_(mesh.faces().size()),
    cellCtrs_(mesh.cellCentres()),
    cellVols_(mesh.cellVolumes()),
    cellSurface_(mesh.cells().size(), 0),
    faceCentreDerivatives_(0),
    cellSurfaceDerivatives_(0),
    cellVolumeDerivatives_(0),
    hessian_(0),
    pointCellSize_(mesh.points().size(), 0),
    boundaryPoint_(mesh.points().size(), false),
    processorPoints_(0),
    boundaryProcessorPoint_(mesh.points().size(), false),
    procPointsLists_(0),
    cellTypeCoef_(mesh.cells().size()),
    metric_(mesh, dict),
    pLabel_(-1)
{
    forAll(mesh_.faces(), fI)
    {
        magfAreas_[fI].setSize(2, 0);
    }
    calculateCellTypes();
    initializeCellSurface();
    calculateFaceCenterDerivatives();
    numberOfCellsPerPoint();
    findBoundaryPoints();
    findBoundaryProcessorPoints();
    if (Pstream::parRun())
    {
        splitProcessorPoints();
    }
}

// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * * //
Foam::meshGeometry::~meshGeometry()
{}
// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
void Foam::meshGeometry::limitDisplacement(const label& pI, vector& ds)
{
    scalar maxDisp = findMaxDisplacement(pI);

    if (mag(ds)>maxDisp)
    {
        ds *= maxDisp/mag(ds);
    }
}

Foam::scalar Foam::meshGeometry::findMaxDisplacement(const label& pI)
{
    //Define the maxDisplacement as twice the cubic root of the mean volume
    const labelList& cells = mesh_.pointCells()[pI];
    scalar maxDisp = 0;
    scalar averageAR = 0;

    const scalarField& AR = metric_.aspectRatio();

    forAll(cells, cI)
    {
        const label& cL = cells[cI];
        maxDisp += cellVols_[cL];
        averageAR += AR[cL];
    }

    maxDisp /= pointCellSize_[pI];
    averageAR /= pointCellSize_[pI];

    maxDisp = 2*cbrt(maxDisp)/sqr(averageAR);

    return maxDisp;
}

void Foam::meshGeometry::processorPointsMaxDisplacement
(
    scalarField& maxDisp,
    const labelList& procPoints
)
{
    const scalarField& AR = metric_.aspectRatio();
    scalarField averageAR(maxDisp.size(), 0);
    forAll(procPoints, pI)
    {
        const label& pL = procPoints[pI];
        const labelList& cells = mesh_.pointCells()[pL];
        forAll(cells, cI)
        {
            const label& cL = cells[cI];
            maxDisp[pI] += cellVols_[cL];
            averageAR[pI] += AR[cL];
        }
        maxDisp[pI] /= pointCellSize_[pL];
        averageAR[pI] /= pointCellSize_[pL];
    }
    syncTools::syncPointList
    (
        mesh_,
        procPoints,
        maxDisp,
        plusEqOp<scalar>(),
        scalar(0)
    );
    syncTools::syncPointList
    (
        mesh_,
        procPoints,
        averageAR,
        plusEqOp<scalar>(),
        scalar(0)
    );
    forAll(procPoints, pI)
    {
        maxDisp[pI] = 2*cbrt(maxDisp[pI])/sqr(averageAR[pI]);
    }
}

const Foam::scalarField& Foam::meshGeometry::cellVols() const
{
    return cellVols_;
}

const Foam::scalarField& Foam::meshGeometry::cellSurface() const
{
    return cellSurface_;
}

const Foam::List<Foam::tensor>& Foam::meshGeometry::faceCentreDerivatives() const
{
    return faceCentreDerivatives_;
}

const Foam::List<Foam::vector>& Foam::meshGeometry::cellVolumeDerivatives()
const
{
    return cellVolumeDerivatives_;
}

const Foam::List<Foam::vector>& Foam::meshGeometry::cellSurfaceDerivatives()
const
{
    return cellSurfaceDerivatives_;
}

const Foam::List<Foam::tensor>& Foam::meshGeometry::surfaceHessian() const
{
    return hessian_;
}

const Foam::vectorField& Foam::meshGeometry::cellCtrs() const
{
    return cellCtrs_;
}

const Foam::vectorField& Foam::meshGeometry::fCtrs() const
{
    return fCtrs_;
}

const Foam::vectorField& Foam::meshGeometry::fAreas() const
{
    return fAreas_;
}

// ************************************************************************* //
