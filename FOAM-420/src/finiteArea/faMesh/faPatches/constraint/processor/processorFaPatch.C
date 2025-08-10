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
    (c) 2016-2017 Wikki Ltd
    (c) 2018-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "faMesh/faPatches/constraint/processor/processorFaPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/IOstreams/Pstreams/IPstream.H"
#include "db/IOstreams/Pstreams/OPstream.H"
#include "fields/Fields/transformField/transformField.H"
#include "faMesh/faBoundaryMesh/faBoundaryMesh.H"
#include "faMesh/faMesh.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"
#include "db/IOstreams/Pstreams/PstreamBuffers.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(processorFaPatch, 0);
addToRunTimeSelectionTable(faPatch, processorFaPatch, dictionary);


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorFaPatch::~processorFaPatch()
{
    deleteDemandDrivenData(nbrPointsPtr_);
    deleteDemandDrivenData(nonGlobalPatchPointsPtr_);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void processorFaPatch::makeNonGlobalPatchPoints() const
{
    // If it is not running parallel or there are no global points
    // create a 1->1 map

    // Can not use faGlobalMeshData at this point yet

    if
    (
        !Pstream::parRun()
        || !boundaryMesh().mesh()().globalData().nGlobalPoints()
//         || !boundaryMesh().mesh().globalData().nGlobalPoints()
    )
    {
        nonGlobalPatchPointsPtr_ = new labelList(nPoints());
        labelList& ngpp = *nonGlobalPatchPointsPtr_;
        forAll(ngpp, i)
        {
            ngpp[i] = i;
        }
    }
    else
    {
        // Get reference to shared points
        const labelList& sharedPoints =
            boundaryMesh().mesh()().globalData().sharedPointLabels();

        nonGlobalPatchPointsPtr_ = new labelList(nPoints());
        labelList& ngpp = *nonGlobalPatchPointsPtr_;

        const labelList& faMeshPatchPoints = pointLabels();

        const labelList& meshPoints =
            boundaryMesh().mesh().patch().meshPoints();

        label noFiltPoints = 0;

        forAll(faMeshPatchPoints, pointI)
        {
            label curP = meshPoints[faMeshPatchPoints[pointI]];

            bool found = false;

            forAll(sharedPoints, sharedI)
            {
                if (sharedPoints[sharedI] == curP)
                {
                    found = true;
                    break;
                }
            }

            if (!found)
            {
                ngpp[noFiltPoints] = pointI;
                noFiltPoints++;
            }
        }

        ngpp.setSize(noFiltPoints);


//         // Get reference to shared points
//         const labelList& sharedPoints =
//             boundaryMesh().mesh().globalData().sharedPointLabels();

//         nonGlobalPatchPointsPtr_ = new labelList(nPoints());
//         labelList& ngpp = *nonGlobalPatchPointsPtr_;

//         const labelList& patchPoints = pointLabels();

//         label noFiltPoints = 0;

//         forAll(patchPoints, pointI)
//         {
//             label curP = patchPoints[pointI];

//             bool found = false;

//             forAll(sharedPoints, pI)
//             {
//                 if (sharedPoints[pI] == curP)
//                 {
//                     found = true;
//                     break;
//                 }
//             }

//             if (!found)
//             {
//                 ngpp[noFiltPoints] = pointI;
//                 noFiltPoints++;
//             }
//         }

//         ngpp.setSize(noFiltPoints);
    }
}


void processorFaPatch::initGeometry(PstreamBuffers& pBufs)
{
    faPatch::initGeometry(pBufs);
    if (Pstream::parRun())
    {
        UOPstream toNeighbProc(neighbProcNo(), pBufs);
//        OPstream toNeighbProc
//        (
//            Pstream::blocking,
//            neighbProcNo(),
//            3*(sizeof(label) + size()*sizeof(vector))
//        );
        toNeighbProc
            << edgeCentres()
            << edgeLengths()
            << edgeFaceCentres();
    }
}


void processorFaPatch::calcGeometry(PstreamBuffers& pBufs)
{
    faPatch::calcGeometry(pBufs);
    if (Pstream::parRun())
    {
        {
            UIPstream fromNeighbProc(neighbProcNo(), pBufs);
//            IPstream fromNeighbProc
//            (
//                Pstream::blocking,
//                neighbProcNo(),
//                3*(sizeof(label) + size()*sizeof(vector))
//            );
            fromNeighbProc
                >> nbrEdgeCentres_
                >> nbrEdgeLengths_
                >> nbrEdgeFaceCentres_;
        }

        const scalarField& magEl = magEdgeLengths();
        forAll(magEl, edgei)
        {
            label nEdge = nbrEdges()[edgei];
            scalar nmagEl = mag(nbrEdgeLengths_[nEdge]);
            scalar avEl = (magEl[edgei] + nmagEl)/2.0;

            if (mag(magEl[edgei] - nmagEl)/avEl > 1e-6)
            {
                FatalErrorInFunction
                    << "edge " << edgei
                    << " length does not match neighbour by "
                    << 100*mag(magEl[edgei] - nmagEl)/avEl
                    << "% -- possible edge ordering problem"
                    << exit(FatalError);
            }
        }

        calcTransformTensors
        (
            edgeCentres(),
            nbrEdgeCentres_,
            edgeNormals(),
            nbrEdgeLengths_/mag(nbrEdgeLengths_)
        );
    }
}


void processorFaPatch::initMovePoints(PstreamBuffers& pBufs, const pointField& p)
{
    faPatch::movePoints(pBufs, p);
    initGeometry(pBufs);
}


void processorFaPatch::movePoints(PstreamBuffers& pBufs, const pointField&)
{
    calcGeometry(pBufs);
}


void processorFaPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    // For completeness
    faPatch::initUpdateMesh(pBufs);
    deleteDemandDrivenData(nbrPointsPtr_);

    if (Pstream::parRun())
    {
        // Express all points as patch face and index in face.
        const edgeList::subList patchEdges =
            patchSlice(boundaryMesh().mesh().edges());

        labelList pointFace(nPoints());
        labelList pointIndex(nPoints());

        label fvProcPatch = ngbPolyPatchIndex();
        const polyMesh& M = boundaryMesh().mesh().mesh();
        const polyPatch& pp = M.boundaryMesh()[fvProcPatch];

        for (label patchPointI = 0; patchPointI < nPoints(); patchPointI++)
        {
            //Point in the indirect primitive patch of faMesh
            label pPoint = pointLabels()[patchPointI];
            //global fvMesh Point
            label gPoint = boundaryMesh().mesh().patch().meshPoints()[pPoint];
            //local point in the processor Patch
            label procP = pp.whichPoint(gPoint);

            label faceI = pp.pointFaces()[procP][0];
//            label edgeI = ptEdges[patchPointI][0];

//            patchEdge[patchPointI] = edgeI;
            pointFace[patchPointI] = faceI;

//            const edge& e = patchEdges[edgeI];
            const face& f = pp.localFaces()[faceI];

//            indexInEdge[patchPointI] =
//                findIndex
//                (
//                    e,
//                    pointLabels()[patchPointI]
//                );
            pointIndex[patchPointI] =
                findIndex
                (
                    f,
                    procP
                );
        }
        label nEdges = this->size();

        labelList edgeFace(nEdges);
        labelList edgeIndex(nEdges);

        labelList faceCells(boundaryMesh().mesh().faceLabels().size(), -1);

        forAll(faceCells, faceI)
        {
            label faceID = boundaryMesh().mesh().faceLabels()[faceI];

            faceCells[faceI] = M.faceOwner()[faceID];
        }

        labelList meshEdges =
            boundaryMesh().mesh().patch().meshEdges
            (
                M.edges(),
                M.cellEdges(),
                faceCells
            );

        for (label patchEdgeI = 0; patchEdgeI < nEdges; patchEdgeI++)
        {
            //Egde in the indirect primitive patch of faMesh
            edge dummyEdge =  patchEdges[patchEdgeI];

            label p0 = dummyEdge[0];
            label p1 = dummyEdge[1];
            //Point in the indirect primitive patch of faMesh
            //global fvMesh Point
            label gPoint0 = boundaryMesh().mesh().patch().meshPoints()[p0];
            label gPoint1 = boundaryMesh().mesh().patch().meshPoints()[p1];
            //local point in the processor Patch
            label procP0 = pp.whichPoint(gPoint0);
            label procP1 = pp.whichPoint(gPoint1);

            //local edge in the processor Patch
            dummyEdge[0] = procP0;
            dummyEdge[1] = procP1;

            label edgeLabel = pp.whichEdge(dummyEdge);

            label faceI = pp.edgeFaces()[edgeLabel][0];

            edgeFace[patchEdgeI] = faceI;

            const labelList& fEdges = pp.faceEdges()[faceI];

            edgeIndex[patchEdgeI] = findIndex(fEdges, edgeLabel);
        }

        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << edgeFace
            << edgeIndex
            << pointFace
            << pointIndex;
    }
}


void processorFaPatch::updateMesh(PstreamBuffers& pBufs)
{
    // For completeness
    faPatch::updateMesh(pBufs);

    if (Pstream::parRun())
    {
        labelList nbrEdgeFace;
        labelList nbrEdgeIndex;
        labelList nbrPointFace;
        labelList nbrPointIndex;

        {
            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> nbrEdgeFace
                >> nbrEdgeIndex
                >> nbrPointFace
                >> nbrPointIndex;
        }
        if (nbrPointFace.size() == nPoints())
        {
            // Convert neighbour edges and indices into face back into
            // my edges and points.
            nbrPointsPtr_ = new labelList(nPoints(), -1);
            labelList& nbrPoints = *nbrPointsPtr_;

            label fvProcPatch = ngbPolyPatchIndex();
            const polyMesh& M = boundaryMesh().mesh().mesh();
            const polyPatch& pp = M.boundaryMesh()[fvProcPatch];

            forAll(nbrPointFace, nbrPointI)
            {
                // Find edge and index in edge on this side.
//                const edge& e = patchEdges[nbrPatchEdge[nbrPointI]];
                const face& f = pp.localFaces()[nbrPointFace[nbrPointI]];
//                label index =  1 - nbrIndexInEdge[nbrPointI];
//                label patchPointI = findIndex(pointLabels(), e[index]);

                label index = (f.size() - nbrPointIndex[nbrPointI]) % f.size();
                //Point in the processor Patch
                label pPoint =  f[index];
                //Point on fvMesh
                label gPoint = pp.meshPoints()[pPoint];
                //Point on the indirect primitive Patch of faMesh
                label patchPoint = boundaryMesh().mesh().patch().whichPoint(gPoint);
                //Point of faProcessor Patch
                label patchPointI = findIndex(pointLabels(), patchPoint);

                if (nbrPoints[patchPointI] == -1)
                {
                    // First reference of point
                    nbrPoints[patchPointI] = nbrPointI;
                }
                else if (nbrPoints[patchPointI] >= 0)
                {
                    // Point already visited. Mark as duplicate.
                    Pout<<"duplicate" << pointEdges()[patchPointI]<<endl;
                    nbrPoints[patchPointI] = -2;
                }

            }
            forAll(nbrPoints, patchPointI)
            {
                if (nbrPoints[patchPointI] == -2)
                {
                    nbrPoints[patchPointI] = -1;
                }
            }
            // Convert edges.
            // ~~~~~~~~~~~~~~

            nbrEdgesPtr_ = new labelList(this->size(), -1);
            labelList& nbrEdges = *nbrEdgesPtr_;
            forAll(nbrEdgeFace, nbrEdgeI)
            {
                // Find face and index in face on this side.
                const labelList& f = pp.faceEdges()[nbrEdgeFace[nbrEdgeI]];

                label index = (f.size() - nbrEdgeIndex[nbrEdgeI] - 1) % f.size();
                //Edge in the processor Patch
                label pEdge = f[index];
                edge dummyEdge = pp.edges()[pEdge];
                label p0 = dummyEdge[0];
                label p1 = dummyEdge[1];
                //Point in the indirect primitive patch of faMesh
                //global fvMesh Point
                label gPoint0 = pp.meshPoints()[p0];
                label gPoint1 = pp.meshPoints()[p1];
                //local point in the processor Patch
                label procP0 = boundaryMesh().mesh().patch().whichPoint(gPoint0);
                label procP1 = boundaryMesh().mesh().patch().whichPoint(gPoint1);

                //local edge in the processor Patch
                dummyEdge[0] = procP0;
                dummyEdge[1] = procP1;

                label patchEdgeI = boundaryMesh().mesh().patch().whichEdge(dummyEdge);
                patchEdgeI = findIndex(*this, patchEdgeI);

                if (nbrEdges[patchEdgeI] == -1)
                {
                    // First reference of edge
                    nbrEdges[patchEdgeI] = nbrEdgeI;
                }
                else if (nbrEdges[patchEdgeI] >= 0)
                {
                    // Edge already visited. Mark as duplicate.
                    nbrEdges[patchEdgeI] = -2;
                }
            }
        }
        else
        {
            // Differing number of points. Probably patch includes
            // part of a cyclic.
//            nbrPointsPtr_ = nullptr;
        }
    }
}


const labelList& processorFaPatch::nbrPoints() const
{
    if (!nbrPointsPtr_)
    {
        // Was probably created from cyclic patch and hence the
        // number of edges or points might differ on both
        // sides of the processor patch since one side might have
        // it merged with another bit of geometry

        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << nl
            << "This can happen if the number of points  on both"
            << " sides of the two coupled patches differ." << nl
            << "This happens if the processorPatch was constructed from"
            << " part of a cyclic patch."
            << abort(FatalError);
    }

   return *nbrPointsPtr_;
}

const labelList& processorFaPatch::nbrEdges() const
{
    if (!nbrEdgesPtr_)
    {
        // Was probably created from cyclic patch and hence the
        // number of edges or points might differ on both
        // sides of the processor patch since one side might have
        // it merged with another bit of geometry

        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << nl
            << "This can happen if the number of points  on both"
            << " sides of the two coupled patches differ." << nl
            << "This happens if the processorPatch was constructed from"
            << " part of a cyclic patch."
            << abort(FatalError);
    }

   return *nbrEdgesPtr_;
}

// Make patch weighting factors
void processorFaPatch::makeWeights(scalarField& w) const
{
    if (Pstream::parRun())
    {
        // The face normals point in the opposite direction on the other side
        scalarField nbrEdgeCentresCn
        (
            (
                nbrEdgeLengths()
               /mag(nbrEdgeLengths())
            )
          & (
              nbrEdgeCentres()
            - nbrEdgeFaceCentres()
            )
        );
        w = nbrEdgeCentresCn/
            (
                (edgeNormals() & faPatch::delta())
              + nbrEdgeCentresCn
            );
    }
    else
    {
        w = 1.0;
    }
}


// Make patch edge - neighbour face distances
void processorFaPatch::makeDeltaCoeffs(scalarField& dc) const
{
    if (Pstream::parRun())
    {
        dc = (1.0 - weights())/(edgeNormals() & faPatch::delta());
    }
    else
    {
        dc = 1.0/(edgeNormals() & faPatch::delta());
    }
}


// Return delta (P to N) vectors across coupled patch
tmp<vectorField> processorFaPatch::delta() const
{
    if (Pstream::parRun())
    {
        // To the transformation if necessary
        if (parallel())
        {
            return
                coupledFaPatch::delta()
              - (
                    nbrEdgeCentres()
                  - nbrEdgeFaceCentres()
                );
        }
        else
        {
            return
                coupledFaPatch::delta()
              - Foam::transform
                (
                    forwardT(),
                    (
                        nbrEdgeCentres()
                      - nbrEdgeFaceCentres()
                    )
                );
        }
    }
    else
    {
        return coupledFaPatch::delta();
    }
}


const labelList& processorFaPatch::nonGlobalPatchPoints() const
{
    if (!nonGlobalPatchPointsPtr_)
    {
        makeNonGlobalPatchPoints();
    }

    return *nonGlobalPatchPointsPtr_;
}

tmp<labelField> processorFaPatch::interfaceInternalField
(
    const labelUList& internalData
) const
{
    return patchInternalField(internalData);
}


void processorFaPatch::initTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& interfaceData
) const
{
    send(commsType, interfaceData);
}


tmp<labelField> processorFaPatch::transfer
(
    const Pstream::commsTypes commsType,
    const labelUList&
) const
{
    return receive<label>(commsType, this->size());
}

void processorFaPatch::initInternalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList& iF
) const
{
    send(commsType, patchInternalField(iF)());
}

tmp<labelField> processorFaPatch::internalFieldTransfer
(
    const Pstream::commsTypes commsType,
    const labelUList&
) const
{
    return receive<label>(commsType, this->size());
}


// Write
void processorFaPatch::write(Ostream& os) const
{
    faPatch::write(os);
    os.writeEntry("myProcNo", myProcNo_);
    os.writeEntry("neighbProcNo", neighbProcNo_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
