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

#include "refinementFeatures/manifoldFeatures.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(manifoldFeatures, 0);

} // End namespace Foam

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::manifoldFeatures::update
(
    const mapPolyMesh& map
)
{
    if (seeds_.size() > 0)
    {
        const labelList& rMap = map.reverseCellMap();

        forAll(seeds_, i)
        {
            label oldCellI = seeds_[i];
            if (oldCellI != -1)
            {
                label newCellI = rMap[oldCellI];
                seeds_[i] = newCellI;
            }
        }
    }
}

void Foam::manifoldFeatures::distribute
(
    const mapDistributePolyMesh& map
)
{
    label nManifold = seeds_.size();

    if (nManifold > 0)
    {
        labelList originProc(nManifold,-1);
        forAll(originProc, i)
        {
            if (seeds_[i] != -1)
            {
                originProc[i] = Pstream::myProcNo();
            }
        }
        Pstream::listCombineReduce(originProc, maxOp<label>());

        label nSeeds = 0;
        labelList seedMap(seeds_.size(),-1);
        forAll(originProc, i)
        {
            if (originProc[i] != -1)
            {
                originProc[nSeeds] = originProc[i];
                seedMap[nSeeds] = i;
                nSeeds++;
            }
        }
        seedMap.setSize(nSeeds);
        originProc.setSize(nSeeds);

        const labelListList& subCellMap = map.cellMap().subMap();

        label nOldCells = map.nOldCells();

        labelList destination(nOldCells, -1);
        labelList localSendIndex(nOldCells, -1);

        forAll(subCellMap, procI)
        {
            const labelList& newToOld = subCellMap[procI];

            forAll(newToOld, i)
            {
                label oldCellI = newToOld[i];
                destination[oldCellI] = procI;
                localSendIndex[oldCellI] = i;
            }
        }

        labelList destinationProc(nSeeds, -1);
        labelList sendingIndex(nSeeds, -1);

        forAll(destinationProc, i)
        {
            label oProc = originProc[i];
            label oldCellI = seeds_[seedMap[i]];

            if (oProc == Pstream::myProcNo())
            {
                destinationProc[i] = destination[oldCellI];
                sendingIndex[i] = localSendIndex[oldCellI];
            }
        }

        Pstream::listCombineReduce(destinationProc, maxOp<label>());
        Pstream::listCombineReduce(sendingIndex, maxOp<label>());

        const labelListList& subConstructMap = map.cellMap().constructMap();

        seeds_ = -1;
        forAll(destinationProc,i)
        {
            if (Pstream::myProcNo() == destinationProc[i])
            {
                const labelList& oldToNew = subConstructMap[originProc[i]];
                seeds_[seedMap[i]] = oldToNew[sendingIndex[i]];
            }
        }
    }
}


Foam::PtrList<Foam::Tuple2<Foam::labelList, Foam::label>>&
Foam::manifoldFeatures::features
(
    const label maxFtrSz,
    const scalar minFtrLength
)
{
    if (features_.empty())
    {
        calcManifoldFeatureMeshes(features_,maxFtrSz,minFtrLength);
        seeds_.setSize(features_.size(), -1);
    }

    return features_;
}


void Foam::manifoldFeatures::clear()
{
    if (!features_.empty())
    {
        features_.clear();
        seeds_.clear();
    }
}


Foam::scalar Foam::manifoldFeatures::calcFeatureLength
(
    const edgeMesh& eMesh,
    const labelList& meshPoints
)
{
    const pointField& pts = eMesh.points();

    scalar featureLength = 0;
    for (int i = 0; i < meshPoints.size()-1; i++)
    {
        point start = pts[meshPoints[i+1]];
        point end = pts[meshPoints[i]];
        featureLength += mag(end-start);
    }

    return featureLength;
}


//Calculate a set on manifold feature edges
void Foam::manifoldFeatures::calcManifoldFeatureMeshes
(
    PtrList<Tuple2<labelList,label>>& manifoldFeatures,
    const label maxFtrSz,
    const scalar minFtrLength
)
{
    label sz = 0;
    forAll(extFeatures_, featI)
    {
        const edgeMesh& eMesh = extFeatures_[featI];
        sz += eMesh.edges().size();
    }
    manifoldFeatures.setSize(sz);

    sz = 0;
    forAll(extFeatures_, featI)
    {
        if (checkRefinementOnly_ && refinementOnly_[featI])
        {
            continue;
        }

        const edgeMesh& eMesh = extFeatures_[featI];
        const pointField& ePoints = eMesh.points();
        const labelListList& pointEdges = eMesh.pointEdges();

        // Whether edge has been visited.
        PackedList<1> edgeVisited(eMesh.edges().size(), 0u);
        PackedList<1> visited(eMesh.edges().size(), 0u);

        forAll(pointEdges, pointI)
        {
            if (pointEdges[pointI].size() != 2)
            {
                forAll(pointEdges[pointI], i)
                {
                    label edgeI = pointEdges[pointI][i];

                    if (edgeVisited.set(edgeI, 1u))
                    {
                        PackedList<1> pointVisited(ePoints.size(), 0u);
                        DynamicList<label> meshPoints(ePoints.size());

                        // Unvisited edge. Make the particle go to
                        // the other point on the edge.
                        const edge& e = eMesh.edges()[edgeI];
                        label otherPointI = e.otherVertex(pointI);
                        meshPoints.append(pointI);
                        meshPoints.append(otherPointI);

                        pointVisited.set(pointI, 1u);
                        pointVisited.set(otherPointI, 1u);

                        manifoldEdgeWave
                        (
                            pointI,
                            otherPointI,
                            pointEdges,
                            eMesh,
                            edgeVisited,
                            pointVisited,
                            meshPoints,
                            maxFtrSz
                        );
                        meshPoints.shrink();

                        bool validFtr = true;
                        if
                        (
                            minFtrLength > -1 &&
                            calcFeatureLength(eMesh,meshPoints) < minFtrLength
                        )
                        {
                            validFtr = false;
                        }

                        if (validFtr)
                        {
                            manifoldFeatures.set
                            (
                                sz,
                                new Tuple2<labelList,label>
                                (
                                    meshPoints,
                                    featI
                                 )
                             );
                            sz++;
                        }
                    }
                }
            }
        }

        //start from a feature point
        forAll(pointEdges, pointI)
        {
            forAll(pointEdges[pointI], i)
            {
                label edgeI = pointEdges[pointI][i];

                bool foundFeat =  false;
                if (edgeVisited.get(edgeI) == 0u)
                {
                    const edge& e = eMesh.edges()[edgeI];
                    label otherPointI = e.otherVertex(pointI);

                    forAll(pointEdges[pointI], j)
                    {
                        label edgeJ = pointEdges[pointI][j];
                        if (edgeJ != edgeI)
                        {
                            const edge& ej = eMesh.edges()[edgeJ];
                            label otherPointJ = ej.otherVertex(pointI);

                            vector vecA = ePoints[otherPointI] - ePoints[pointI];
                            vecA /= (mag(vecA) + SMALL);

                            vector vecB = ePoints[pointI] - ePoints[otherPointJ];
                            vecB /= (mag(vecB) + SMALL);
                            if ((vecA & vecB) < 0.9659)
                            {
                                foundFeat = true;
                                break;
                            }
                        }
                    }
                }

                if (foundFeat && edgeVisited.set(edgeI, 1u))
                {
                    PackedList<1> pointVisited(ePoints.size(), 0u);
                    DynamicList<label> meshPoints(ePoints.size());
                    // Unvisited edge. Make the particle go to
                    // the other point on the edge.
                    const edge& e = eMesh.edges()[edgeI];
                    label otherPointI = e.otherVertex(pointI);
                    meshPoints.append(pointI);
                    meshPoints.append(otherPointI);

                    pointVisited.set(pointI, 1u);
                    pointVisited.set(otherPointI, 1u);

                    manifoldEdgeWave
                    (
                        pointI,
                        otherPointI,
                        pointEdges,
                        eMesh,
                        edgeVisited,
                        pointVisited,
                        meshPoints,
                        maxFtrSz
                    );
                    meshPoints.shrink();

                    bool validFtr = true;
                    if
                    (
                        minFtrLength > -1 &&
                        calcFeatureLength(eMesh,meshPoints) < minFtrLength
                    )
                    {
                        validFtr = false;
                    }

                    if (validFtr)
                    {
                        manifoldFeatures.set
                        (
                            sz,
                            new Tuple2<labelList,label>
                            (
                                meshPoints,
                                featI
                             )
                         );
                        sz++;
                    }
                }
            }
        }

        forAll(pointEdges, pointI)
        {
            forAll(pointEdges[pointI], i)
            {
                label edgeI = pointEdges[pointI][i];

                if (edgeVisited.set(edgeI, 1u))
                {
                    PackedList<1> pointVisited(ePoints.size(), 0u);
                    DynamicList<label> meshPoints(ePoints.size());
                    // Unvisited edge. Make the particle go to
                    // the other point on the edge.
                    const edge& e = eMesh.edges()[edgeI];
                    label otherPointI = e.otherVertex(pointI);
                    meshPoints.append(pointI);
                    meshPoints.append(otherPointI);

                    pointVisited.set(pointI, 1u);
                    pointVisited.set(otherPointI, 1u);

                    manifoldEdgeWave
                    (
                        pointI,
                        otherPointI,
                        pointEdges,
                        eMesh,
                        edgeVisited,
                        pointVisited,
                        meshPoints,
                        maxFtrSz
                    );
                    meshPoints.shrink();

                    bool validFtr = true;
                    if
                    (
                        minFtrLength > -1 &&
                        calcFeatureLength(eMesh,meshPoints) < minFtrLength
                    )
                    {
                        validFtr = false;
                    }

                    if (validFtr)
                    {
                        manifoldFeatures.set
                        (
                            sz,
                            new Tuple2<labelList,label>
                            (
                                meshPoints,
                                featI
                             )
                         );
                        sz++;
                    }
                }
            }
        }
    }

    manifoldFeatures.setSize(sz);

    return;
}
