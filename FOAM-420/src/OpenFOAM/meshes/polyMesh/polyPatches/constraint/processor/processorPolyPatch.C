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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2015 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "db/dictionary/dictionary.H"
#include "fields/Fields/Field/SubField.H"
#include "include/demandDrivenData.H"
#include "meshes/meshTools/matchPoints.H"
#include "db/IOstreams/Fstreams/OFstream.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"
#include "db/IOstreams/Pstreams/PstreamBuffers.H"
#include "containers/Circulators/ConstCirculator/ConstCirculator.H"
#include "meshes/polyMesh/globalMeshData/globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(processorPolyPatch, 0);
    addToRunTimeSelectionTable(polyPatch, processorPolyPatch, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::processorPolyPatch::processorPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const word& patchType
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    nbrFaceCentres_(),
    nbrFaceAreas_(),
    nbrFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const int myProcNo,
    const int neighbProcNo,
    const word& patchType
)
:
    coupledPolyPatch
    (
        newName(myProcNo, neighbProcNo),
        size,
        start,
        index,
        bm,
        patchType
    ),
    myProcNo_(myProcNo),
    neighbProcNo_(neighbProcNo),
    neighbProPatchID_(-1),
    nbrFaceCentres_(),
    nbrFaceAreas_(),
    nbrFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    myProcNo_(readLabel(dict.lookup("myProcNo"))),
    neighbProcNo_(readLabel(dict.lookup("neighbProcNo"))),
    neighbProPatchID_(-1),
    nbrFaceCentres_(),
    nbrFaceAreas_(),
    nbrFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbProPatchID_(pp.neighbProPatchID_),
    nbrFaceCentres_(),
    nbrFaceAreas_(),
    nbrFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbProPatchID_(pp.neighbProPatchID_),
    nbrFaceCentres_(),
    nbrFaceAreas_(),
    nbrFaceCellCentres_()
{}


Foam::processorPolyPatch::processorPolyPatch
(
    const processorPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    myProcNo_(pp.myProcNo_),
    neighbProcNo_(pp.neighbProcNo_),
    neighbProPatchID_(pp.neighbProPatchID_),
    nbrFaceCentres_(),
    nbrFaceAreas_(),
    nbrFaceCellCentres_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::processorPolyPatch::~processorPolyPatch()
{
    nbrPointsPtr_.clear();
    nbrEdgesPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::processorPolyPatch::newName
(
    const label myProcNo,
    const label neighbProcNo
)
{
    return
        "procBoundary"
      + Foam::name(myProcNo)
      + "to"
      + Foam::name(neighbProcNo);
}


void Foam::processorPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << faceCentres()
            << faceAreas()
            << faceCellCentres();
    }
}


void Foam::processorPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (Pstream::parRun())
    {
        {
            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> nbrFaceCentres_
                >> nbrFaceAreas_
                >> nbrFaceCellCentres_;
        }

        // My normals
        vectorField faceNormals(size());

        // Neighbour normals
        vectorField nbrFaceNormals(nbrFaceAreas_.size());

        // Face match tolerances
        const vectorField fc(faceCentres());
        scalarField tols = calcFaceTol(*this, points(), fc);

        // Calculate normals from areas and check
        forAll(faceNormals, facei)
        {
            scalar magSf = magFaceAreas()[facei];
            scalar nbrMagSf = mag(nbrFaceAreas_[facei]);
            scalar avSf = (magSf + nbrMagSf)/2.0;

            // For small face area calculation the results of the area
            // calculation have been found to only be accurate to ~1e-20
            if (magSf < SMALL || nbrMagSf < SMALL)
            {
                // Undetermined normal. Use dummy normal to force separation
                // check.
                faceNormals[facei] = point(1, 0, 0);
                nbrFaceNormals[facei] = -faceNormals[facei];
                tols[facei] = GREAT;
            }
            else if (mag(magSf - nbrMagSf) > matchTolerance()*sqr(tols[facei]))
            {
                const fileName patchOBJName
                (
                    boundaryMesh().mesh().time().path()/name() + "_faces.obj"
                );

                Pout<< "processorPolyPatch::calcGeometry : Writing my "
                    << size() << " faces to " << patchOBJName << endl;

                writeOBJ(patchOBJName, *this);

                const fileName centresOBJName
                (
                    boundaryMesh().mesh().time().path()/name()
                  + "_faceCentresConnections.obj"
                );

                Pout<< "processorPolyPatch::calcGeometry :"
                    << " Dumping lines between corresponding face centres to "
                    << centresOBJName.name() << endl;

                writeOBJ(centresOBJName, nbrFaceCentres_, faceCentres());

                FatalErrorInFunction
                    << "face " << facei << " area does not match neighbour by "
                    << 100*mag(magSf - nbrMagSf)/avSf
                    << "% -- possible face ordering problem." << endl
                    << "patch:" << name()
                    << " my area:" << magSf
                    << " neighbour area:" << nbrMagSf
                    << " matching tolerance:"
                    << matchTolerance()*sqr(tols[facei])
                    << endl
                    << "Mesh face:" << start()+facei
                    << " vertices:"
                    << UIndirectList<point>(points(), operator[](facei))()
                    << endl
                    << "If you are certain your matching is correct"
                    << " you can increase the 'matchTolerance' setting"
                    << " in the patch dictionary in the boundary file."
                    << endl
                    << "Rerun with processor debug flag set for"
                    << " more information." << exit(FatalError);
            }
            else
            {
                faceNormals[facei] = faceAreas()[facei]/magSf;
                nbrFaceNormals[facei] = nbrFaceAreas_[facei]/nbrMagSf;
            }
        }
    }
}


void Foam::processorPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    directPolyPatch::movePoints(pBufs, p);
    processorPolyPatch::initCalcGeometry(pBufs);
}


void Foam::processorPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField&
)
{
    processorPolyPatch::calcGeometry(pBufs);
}


void Foam::processorPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    directPolyPatch::initUpdateMesh(pBufs);

    if (Pstream::parRun())
    {
        // Express all points as patch face and index in face.
        labelList pointFace(nPoints());
        labelList pointIndex(nPoints());

        for (label patchPointi = 0; patchPointi < nPoints(); patchPointi++)
        {
            label facei = pointFaces()[patchPointi][0];

            pointFace[patchPointi] = facei;

            const face& f = localFaces()[facei];

            pointIndex[patchPointi] = findIndex(f, patchPointi);
        }

        // Express all edges as patch face and index in face.
        labelList edgeFace(nEdges());
        labelList edgeIndex(nEdges());

        for (label patchEdgeI = 0; patchEdgeI < nEdges(); patchEdgeI++)
        {
            label facei = edgeFaces()[patchEdgeI][0];

            edgeFace[patchEdgeI] = facei;

            const labelList& fEdges = faceEdges()[facei];

            edgeIndex[patchEdgeI] = findIndex(fEdges, patchEdgeI);
        }

        UOPstream toNeighbProc(neighbProcNo(), pBufs);

        toNeighbProc
            << pointFace
            << pointIndex
            << edgeFace
            << edgeIndex;
    }
}


void Foam::processorPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    // For completeness
    directPolyPatch::updateMesh(pBufs);

    nbrPointsPtr_.clear();
    nbrEdgesPtr_.clear();

    if (Pstream::parRun())
    {
        labelList nbrPointFace;
        labelList nbrPointIndex;
        labelList nbrEdgeFace;
        labelList nbrEdgeIndex;

        {
            // !Note: there is one situation where the opposite points and
            // do not exactly match and that is while we are distributing
            // meshes and multiple parts come together from different
            // processors. This can temporarily create the situation that
            // we have points which will be merged out later but we still
            // need the face connectivity to be correct.
            // So: cannot check here on same points and edges.

            UIPstream fromNeighbProc(neighbProcNo(), pBufs);

            fromNeighbProc
                >> nbrPointFace
                >> nbrPointIndex
                >> nbrEdgeFace
                >> nbrEdgeIndex;
        }

        // Convert neighbour faces and indices into face back into
        // my edges and points.

        // Convert points.
        // ~~~~~~~~~~~~~~~

        nbrPointsPtr_.reset(new labelList(nPoints(), -1));
        labelList& nbrPoints = nbrPointsPtr_();

        forAll(nbrPointFace, nbrPointi)
        {
            // Find face and index in face on this side.
            const face& f = localFaces()[nbrPointFace[nbrPointi]];

            label index = (f.size() - nbrPointIndex[nbrPointi]) % f.size();
            label patchPointi = f[index];

            if (nbrPoints[patchPointi] == -1)
            {
                // First reference of point
                nbrPoints[patchPointi] = nbrPointi;
            }
            else if (nbrPoints[patchPointi] >= 0)
            {
                // Point already visited. Mark as duplicate.
                nbrPoints[patchPointi] = -2;
            }
        }

        // Reset all duplicate entries to -1.
        forAll(nbrPoints, patchPointi)
        {
            if (nbrPoints[patchPointi] == -2)
            {
                nbrPoints[patchPointi] = -1;
            }
        }

        // Convert edges.
        // ~~~~~~~~~~~~~~

        nbrEdgesPtr_.reset(new labelList(nEdges(), -1));
        labelList& nbrEdges = nbrEdgesPtr_();

        forAll(nbrEdgeFace, nbrEdgeI)
        {
            // Find face and index in face on this side.
            const labelList& f = faceEdges()[nbrEdgeFace[nbrEdgeI]];
            label index = (f.size() - nbrEdgeIndex[nbrEdgeI] - 1) % f.size();
            label patchEdgeI = f[index];

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

        // Reset all duplicate entries to -1.
        forAll(nbrEdges, patchEdgeI)
        {
            if (nbrEdges[patchEdgeI] == -2)
            {
                nbrEdges[patchEdgeI] = -1;
            }
        }

        // Remove any addressing used for shared points/edges calculation
        // since mostly not needed.
        primitivePatch::clearOut();
    }
}


const Foam::labelList& Foam::processorPolyPatch::nbrPoints() const
{
    if (!nbrPointsPtr_.valid())
    {
        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << abort(FatalError);
    }
    return nbrPointsPtr_();
}


const Foam::labelList& Foam::processorPolyPatch::nbrEdges() const
{
    if (!nbrEdgesPtr_.valid())
    {
        FatalErrorInFunction
            << "No extended addressing calculated for patch " << name()
            << abort(FatalError);
    }
    return nbrEdgesPtr_();
}


Foam::label Foam::processorPolyPatch::getRotation
(
    const pointField& points,
    const face& f,
    const point& anchor,
    const scalar tol
)
{
    label anchorFp = -1;
    scalar minDistSqr = GREAT;

    forAll(f, fp)
    {
        scalar distSqr = magSqr(anchor - points[f[fp]]);

        if (distSqr < minDistSqr)
        {
            minDistSqr = distSqr;
            anchorFp = fp;
        }
    }

    if (anchorFp == -1 || Foam::sqrt(minDistSqr) > tol)
    {
        return -1;
    }
    else
    {
        // Check that anchor is unique.
        forAll(f, fp)
        {
            scalar distSqr = magSqr(anchor - points[f[fp]]);

            if (distSqr == minDistSqr && fp != anchorFp)
            {
                WarningInFunction
                    << "Cannot determine unique anchor point on face "
                    << UIndirectList<point>(points, f)
                    << endl
                    << "Both at index " << anchorFp << " and " << fp
                    << " the vertices have the same distance "
                    << Foam::sqrt(minDistSqr)
                    << " to the anchor " << anchor
                    << ". Continuing but results might be wrong."
                    << nl << endl;
            }
        }

        // Positive rotation
        return (f.size() - anchorFp) % f.size();
    }
}


bool Foam::processorPolyPatch::geometricMatch
(
    const primitivePatch& pp,
    const vectorField& masterCtrs,
    const vectorField& masterNormals,
    const vectorField& masterAnchors,
    const vectorField& masterFacePointAverages,
    labelList& faceMap,
    labelList& rotation
) const
{
    // Calculate typical distance from face centre
    scalarField tols
    (
        matchTolerance()*calcFaceTol(pp, pp.points(), pp.faceCentres())
    );

    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    bool change = false;

    if (owner())
    {
        // Do nothing (i.e. identical mapping, zero rotation).
        // See comment at top.
        forAll(faceMap, patchFacei)
        {
            faceMap[patchFacei] = patchFacei;
        }

        return change;
    }
    else
    {
        // Geometric match of face centre vectors
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // 1. Try existing ordering and transformation
        bool matchedAll = matchPoints
        (
            pp.faceCentres(),
            masterCtrs,
            pp.faceNormals(),
            masterNormals,
            tols,
            false,
            faceMap
        );

        // Fallback: try using face point average for matching
        if (!matchedAll)
        {
            const pointField& ppPoints = pp.points();

            pointField facePointAverages(pp.size(), Zero);
            forAll(pp, fI)
            {
                const labelList& facePoints = pp[fI];

                forAll(facePoints, pI)
                {
                    facePointAverages[fI] += ppPoints[facePoints[pI]];
                }

                facePointAverages[fI] /= facePoints.size();
            }

            scalarField tols2
            (
                matchTolerance()
                *calcFaceTol(pp, pp.points(), facePointAverages)
            );

            // Note that we do not use the faceNormals anymore for
            // comparison. Since we're
            // having problems with the face centres (e.g. due to extreme
            // aspect ratios) we will probably also have problems with
            // reliable normals calculation
            labelList faceMap2(faceMap.size(), -1);
            matchedAll = matchPoints
            (
                facePointAverages,
                masterFacePointAverages,
                tols2,
                true,
                faceMap2
            );

            forAll(faceMap, oldFacei)
            {
                if (faceMap[oldFacei] == -1)
                {
                    faceMap[oldFacei] = faceMap2[oldFacei];
                }
            }

            matchedAll = true;

            forAll(faceMap, oldFacei)
            {
                if (faceMap[oldFacei] == -1)
                {
                    matchedAll = false;
                }
            }
        }

        if (!matchedAll)
        {
            SeriousErrorInFunction
                << "in patch:" << name() << " : "
                << "Cannot match vectors to faces on both sides of patch"
                << endl
                << "    masterCtrs[0]:" << masterCtrs[0] << endl
                << "    ctrs[0]:" << pp.faceCentres()[0] << endl
                << "    Check your topology changes or maybe you have"
                << " multiple separated (from cyclics) processor patches"
                << endl
                << "    Continuing with incorrect face ordering from now on"
                << endl;

            return false;
        }

        // Set rotation.
        forAll(faceMap, oldFacei)
        {
            // The face f will be at newFacei (after morphing) and we want
            // its anchorPoint (= f[0]) to align with the anchorpoint for
            // the corresponding face on the other side.

            label newFacei = faceMap[oldFacei];

            const point& wantedAnchor = masterAnchors[newFacei];

            rotation[newFacei] = getRotation
            (
                pp.points(),
                pp[oldFacei],
                wantedAnchor,
                tols[oldFacei]
             );

            if (rotation[newFacei] == -1)
            {
                SeriousErrorInFunction
                    << "in patch " << name()
                    << " : "
                    << "Cannot find point on face " << pp[oldFacei]
                    << " with vertices "
                    << UIndirectList<point>(pp.points(), pp[oldFacei])()
                    << " that matches point " << wantedAnchor
                    << " when matching the halves of processor patch "
                    << name()
                    << "Continuing with incorrect face ordering from now on"
                    << endl;

                return false;
            }
        }

        forAll(faceMap, facei)
        {
            if (faceMap[facei] != facei || rotation[facei] != 0)
            {
                return true;
            }
        }

        return false;
    }
}


void Foam::processorPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    if (!Pstream::parRun())
    {
        return;
    }

    UOPstream toNeighbour(neighbProcNo(), pBufs);

    if (owner())
    {
        ownToNbrOrderData ownToNbr;
        autoPtr<ownToNbrDebugOrderData> ownToNbrDebugPtr
        (
            coupledPolyPatch::debug
          ? new ownToNbrDebugOrderData()
          : nullptr
        );

        coupledPolyPatch::initOrder
        (
            ownToNbr,
            ownToNbrDebugPtr,
            pp
        );

        const pointField& ppPoints = pp.points();
        pointField anchors(pp.size());
        // Return the first point
        forAll(pp, facei)
        {
            anchors[facei] = ppPoints[pp[facei][0]];
        }

        // Get the average of the points of each face. This is needed in
        // case the face centroid calculation is incorrect due to the face
        // having a very high aspect ratio.
        pointField facePointAverages(pp.size(), Zero);
        forAll(pp, fI)
        {
            const labelList& facePoints = pp[fI];

            forAll(facePoints, pI)
            {
                facePointAverages[fI] += ppPoints[facePoints[pI]];
            }

            facePointAverages[fI] /= facePoints.size();
        }

        toNeighbour << ownToNbr;
        toNeighbour << pp.faceCentres() << pp.faceNormals()
                    << anchors << facePointAverages << pp.nEdges();
        if (coupledPolyPatch::debug)
        {
            toNeighbour << ownToNbrDebugPtr();
        }
    }
    else
    {
        toNeighbour << pp.nEdges();
    }
}


bool Foam::processorPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    ownToNbrOrderData ownToNbr;
    autoPtr<ownToNbrDebugOrderData> ownToNbrDebugPtr
    (
        coupledPolyPatch::debug
      ? new ownToNbrDebugOrderData()
      : nullptr
    );

    vectorField masterCtrs;
    vectorField masterNormals;
    vectorField masterAnchors;
    vectorField masterFacePointAverages;
    label nNeiPPEdges = -1;

    UIPstream fromOwner(neighbProcNo(), pBufs);
    if (!owner())
    {
        fromOwner >> ownToNbr;
        ownToNbr.transform(transform());
        fromOwner >> masterCtrs >> masterNormals
                  >> masterAnchors >> masterFacePointAverages
                  >> nNeiPPEdges;
        if (coupledPolyPatch::debug)
        {
            fromOwner >> ownToNbrDebugPtr();
            ownToNbrDebugPtr->transform(transform());
        }
    }
    else
    {
        fromOwner >> nNeiPPEdges;
    }

    if (pp.nEdges() != nNeiPPEdges)
    {
        WarningInFunction << "Performing older geometric matching. "
            "If possible please supply testcase to Esi to debug." << endl;

        return
            geometricMatch
            (
                pp,
                masterCtrs,
                masterNormals,
                masterAnchors,
                masterFacePointAverages,
                faceMap,
                rotation
            );
    }
    else
    {
        return
            coupledPolyPatch::order
            (
                ownToNbr,
                ownToNbrDebugPtr,
                pp,
                faceMap,
                rotation
             );
    }
}


void Foam::processorPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);

    os.writeEntry("myProcNo", myProcNo_);
    os.writeEntry("neighbProcNo", neighbProcNo_);
}


// ************************************************************************* //
