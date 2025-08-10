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
    (c) 2016-2024 Esi Ltd.
    (c) 2016-2017 Wikki Ltd

\*---------------------------------------------------------------------------*/

#include "faMesh/faMesh.H"
#include "faMesh/faGlobalMeshData/faGlobalMeshData.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/primitiveMesh/primitiveMesh.H"
#include "include/demandDrivenData.H"
#include "containers/Lists/IndirectList/IndirectList.H"
#include "fields/areaFields/areaFields.H"
#include "fields/edgeFields/edgeFields.H"
#include "faMesh/faMeshLduAddressing.H"
#include "faMesh/faPatches/constraint/wedge/wedgeFaPatch.H"
#include "faMesh/faPatches/faPatch/faPatchData.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "faMesh/faPatches/constraint/symmetry/symmetryFaPatch.H"
#include "faMesh/faPatches/constraint/empty/emptyFaPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faMesh, 0);
}

Foam::word Foam::faMesh::meshSubDir = "faMesh";

const int Foam::faMesh::quadricsFit_ = 0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faMesh::setPrimitiveMeshData()
{
    DebugInFunction << "Setting primitive data" << endl;

    const indirectPrimitivePatch& bp = patch();

    // Set faMesh edges
    edges_.setSize(bp.nEdges());

    label edgeI = -1;


    label nIntEdges = bp.nInternalEdges();

    for (label curEdge = 0; curEdge < nIntEdges; ++curEdge)
    {
        edges_[++edgeI] = bp.edges()[curEdge];
    }

    forAll(boundary(), patchI)
    {
        const labelList& curP = boundary()[patchI];

        forAll(curP, eI)
        {
            edges_[++edgeI] = bp.edges()[curP[eI]];
        }
    }

    nEdges_ = edges_.size();
    nInternalEdges_ = nIntEdges;

    // Set edge owner and neighbour
    edgeOwner_.setSize(nEdges());
    edgeNeighbour_.setSize(nInternalEdges());

    edgeI = -1;

    const labelListList& bpef = bp.edgeFaces();

    for (label curEdge = 0; curEdge < nIntEdges; ++curEdge)
    {
        edgeOwner_[++edgeI] = bpef[curEdge][0];

        edgeNeighbour_[edgeI] = bpef[curEdge][1];
    }

    forAll(boundary(), patchI)
    {
        const labelList& curP = boundary()[patchI];

        forAll(curP, eI)
        {
            edgeOwner_[++edgeI] = bpef[curP[eI]][0];
        }
    }

    // Set number of faces
    nFaces_ = bp.size();

    // Set number of points
    nPoints_ = bp.nPoints();
}


void Foam::faMesh::clearGeomNotOldAreas() const
{
    DebugInFunction << "Clearing geometry" << endl;

    deleteDemandDrivenData(SPtr_);
    deleteDemandDrivenData(patchPtr_);
    deleteDemandDrivenData(patchStartsPtr_);
    deleteDemandDrivenData(LePtr_);
    deleteDemandDrivenData(magLePtr_);
    deleteDemandDrivenData(centresPtr_);
    deleteDemandDrivenData(edgeCentresPtr_);
    deleteDemandDrivenData(faceAreaNormalsPtr_);
    deleteDemandDrivenData(edgeAreaNormalsPtr_);
    deleteDemandDrivenData(pointAreaNormalsPtr_);
    deleteDemandDrivenData(faceCurvaturesPtr_);
    deleteDemandDrivenData(edgeTransformTensorsPtr_);
}

void Foam::faMesh::updateGeomNotOldAreas()
{
    bool haveS = (SPtr_ != nullptr);
    bool haveLe = (LePtr_ != nullptr);
    bool haveMagLe = (magLePtr_ != nullptr);
    bool haveC = (centresPtr_ != nullptr);
    bool haveEC = (edgeCentresPtr_ != nullptr);
    bool haveEdgeNormals = (edgeAreaNormalsPtr_ != nullptr);

    clearGeomNotOldAreas();

    // Now recreate the fields
    if (haveS)
    {
        (void)S();
    }
    if (haveLe)
    {
        (void)Le();
    }
    if (haveMagLe)
    {
        (void)magLe();
    }
    if (haveC)
    {
        (void)areaCentres();
    }
    if (haveEC)
    {
        (void)edgeCentres();
    }
    if (haveEdgeNormals)
    {
        (void)edgeAreaNormals();
    }
}

void Foam::faMesh::clearGeom() const
{
    DebugInFunction << "Clearing geometry" << endl;

    clearGeomNotOldAreas();
    deleteDemandDrivenData(S0Ptr_);
    deleteDemandDrivenData(S00Ptr_);
    deleteDemandDrivenData(correctPatchPointNormalsPtr_);
}


void Foam::faMesh::clearAddressing() const
{
    DebugInFunction << "Clearing addressing" << endl;

    deleteDemandDrivenData(lduPtr_);
}

void Foam::faMesh::clearMeshPhi()
{
    deleteDemandDrivenData(phiPtr_);
}

void Foam::faMesh::storeOldSurf(const scalarField& S)
{

    if (debug)
    {
        InfoInFunction
            << " Storing old time areas since from time " << curTimeIndex_
            << " and time now " << time().timeIndex()
            << " S:" << S.size()
            << endl;
    }

    if (S00Ptr_ && S0Ptr_)
    {
        *S00Ptr_ = *S0Ptr_;
    }

    if (S0Ptr_)
    {
        S0Ptr_->scalarField::operator=(S);
    }
    else
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimArea
        );
        scalarField& S0 = *S0Ptr_;
        // Note: S0 now sized with current areaMesh, not with (potentially
        //       different size) S.
        S0.setSize(S.size());
        S0 = S;
    }

    curTimeIndex_ = time().timeIndex();
}

void Foam::faMesh::clearOut() const
{
    clearGeom();
    clearAddressing();

    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);

    deleteDemandDrivenData(globalMeshDataPtr_);
}

Foam::tmp<Foam::scalarField> Foam::faMesh::sweptAreas
(
    const pointField& newPoints,
    const pointField& oldPoints
)
{
    const edgeList& eList = edges();

    tmp<scalarField> tsweptAreas(new scalarField(eList.size()));
    scalarField& sweptAreas = tsweptAreas.ref();
    const labelList& global = patch().meshPoints();
    if (mesh().calcSolverQties())
    {
           forAll(eList, edgeI)
           {
               const edge& e = eList[edgeI];
               label fOwner = edgeOwner()[edgeI];
               //Temporal workaround.
               vector vN = faceAreaNormals()[fOwner];
               sweptAreas[edgeI] =-
               (0.5*
               (
                   ((newPoints[global[e[1]]]-oldPoints[e[0]])
                 ^ (newPoints[global[e[0]]]-oldPoints[e[0]]))
               )&vN)
               -(0.5*
               (
                   ((oldPoints[e[1]]-oldPoints[e[0]])
                 ^ (newPoints[global[e[1]]]-oldPoints[e[0]]))
               )&vN);
           }
    }
    else
    {
        forAll(eList, edgeI)
        {
           sweptAreas[edgeI] = 0.0;
        }
    }

    return tsweptAreas;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMesh::faMesh(const polyMesh& pMesh)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>(pMesh),
    edgeInterpolation(*this),
    data(mesh()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            time().findInstance(meshDir(), "faceLabels"),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            time().findInstance(meshDir(), "faBoundary"),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        *this
    ),
    comm_(Pstream::worldComm),
    registry_(pMesh.thisDb()),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    phiPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    DebugInFunction << "Creating faMesh from IOobject" << endl;

    setPrimitiveMeshData();

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    if (isFile(pMesh.time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );

        S00();
    }

    // Check the existance of the mesh fluxes, read if present and set the
    // mesh to be moving
    if (isFile(time().timePath()/"faMeshPhi"))
    {
        phiPtr_ = new edgeScalarField
        (
            IOobject
            (
                "faMeshPhi",
                time().timeName(),
                mesh(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        // The mesh is now considered moving so the old-time cell volumes
        // will be required for the time derivatives so if they haven't been
        // read initialise to the current cell volumes
        if (!S0Ptr_)
        {
            S0Ptr_ = new DimensionedField<scalar, areaMesh>
            (
                IOobject
                (
                    "V0",
                    time().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                S()
            );
        }
    }
}


// Construct from components without boundary.
Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const labelList& faceLabels
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>(pMesh),
    edgeInterpolation(*this),
    data(mesh()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        faceLabels
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    comm_(Pstream::worldComm),
    registry_(pMesh.thisDb()),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    phiPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    DebugInFunction << "Creating faMesh from components" << endl;
}


// Construct from definition field
Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const fileName& defFile
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, Foam::UpdateableMeshObject, faMesh>(pMesh),
    edgeInterpolation(*this),
    data(mesh()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        List<label>(0)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    comm_(Pstream::worldComm),
    registry_(pMesh.thisDb()),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    phiPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    DebugInFunction << "Creating faMesh from definition file" << endl;

    // Reading faMeshDefinition dictionary
    IOdictionary faMeshDefinition
    (
        IOobject
        (
            defFile,
            mesh().time().constant(),
            meshSubDir,
            mesh(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    const wordList polyMeshPatches(faMeshDefinition.lookup("polyMeshPatches"));

    const dictionary& bndDict = faMeshDefinition.subDict("boundary");

    const wordList faPatchNames(bndDict.toc());

    List<faPatchData> faPatches(faPatchNames.size() + 1);

    const polyBoundaryMesh& pbm = pMesh.boundaryMesh();

    forAll(faPatchNames, patchI)
    {
        dictionary curPatchDict = bndDict.subDict(faPatchNames[patchI]);

        faPatches[patchI].name_ = faPatchNames[patchI];

        faPatches[patchI].type_ = word(curPatchDict.lookup("type"));

        faPatches[patchI].ownPolyPatchID_ =
            pbm.findPatchID(word(curPatchDict.lookup("ownerPolyPatch")));

        faPatches[patchI].ngbPolyPatchID_ =
            pbm.findPatchID(word(curPatchDict.lookup("neighbourPolyPatch")));
    }

    // Setting faceLabels list size
    label size = 0;

    labelList patchIDs(polyMeshPatches.size(), -1);

    forAll(polyMeshPatches, patchI)
    {
        patchIDs[patchI] = pbm.findPatchID(polyMeshPatches[patchI]);

        size += pbm[patchIDs[patchI]].size();
    }

    faceLabels_ = labelList(size, -1);


    // Filling of faceLabels list
    label faceI = -1;

    sort(patchIDs);

    forAll(polyMeshPatches, patchI)
    {
        label start = pbm[patchIDs[patchI]].start();
        label size  = pbm[patchIDs[patchI]].size();

        for (label i = 0; i < size; ++i)
        {
            faceLabels_[++faceI] = start + i;
        }
    }


    // Determination of faPatch ID for each boundary edge.
    // Result is in the bndEdgeFaPatchIDs list
    labelList faceCells(faceLabels_.size(), -1);

    forAll(faceCells, faceI)
    {
        label faceID = faceLabels_[faceI];

        faceCells[faceI] = mesh().faceOwner()[faceID];
    }

    labelList meshEdges =
        patch().meshEdges
        (
            mesh().edges(),
            mesh().cellEdges(),
            faceCells
        );

    const labelListList& edgeFaces = mesh().edgeFaces();

    const label nTotalEdges = patch().nEdges();
    const label nInternalEdges = patch().nInternalEdges();

    labelList bndEdgeFaPatchIDs(nTotalEdges - nInternalEdges, -1);

    for (label edgeI = nInternalEdges; edgeI < nTotalEdges; ++edgeI)
    {
        label curMeshEdge = meshEdges[edgeI];

        labelList curEdgePatchIDs(2, label(-1));

        label patchI = -1;

        forAll(edgeFaces[curMeshEdge], faceI)
        {
            label curFace = edgeFaces[curMeshEdge][faceI];

            label curPatchID = pbm.whichPatch(curFace);

            if (curPatchID != -1)
            {
                curEdgePatchIDs[++patchI] = curPatchID;
            }
        }

        for (label pI = 0; pI < faPatches.size() - 1; ++pI)
        {
            if
            (
                (
                    curEdgePatchIDs[0] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[1] == faPatches[pI].ngbPolyPatchID_
                )
             ||
                (
                    curEdgePatchIDs[1] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[0] == faPatches[pI].ngbPolyPatchID_
                )
            )
            {
                bndEdgeFaPatchIDs[edgeI - nInternalEdges] = pI;
                break;
            }
        }
    }


    // Set edgeLabels for each faPatch
    for (label pI = 0; pI < (faPatches.size() - 1); ++pI)
    {
        DynamicList<label> tmpList;

        forAll(bndEdgeFaPatchIDs, eI)
        {
            if (bndEdgeFaPatchIDs[eI] == pI)
            {
                tmpList.append(nInternalEdges + eI);
            }
        }

        faPatches[pI].edgeLabels_ = tmpList;
    }

    // Check for undefined edges
    DynamicList<label> tmpList;

    forAll(bndEdgeFaPatchIDs, eI)
    {
        if (bndEdgeFaPatchIDs[eI] == -1)
        {
            tmpList.append(nInternalEdges + eI);
        }
    }

    if (tmpList.size() > 0)
    {
        // Check for processor edges
        labelList allUndefEdges = tmpList;
        labelList ngbPolyPatch(allUndefEdges.size(), -1);
        forAll(ngbPolyPatch, edgeI)
        {
            label curEdge = allUndefEdges[edgeI];

            label curPMeshEdge = meshEdges[curEdge];

            forAll(edgeFaces[curPMeshEdge], faceI)
            {
                label curFace = edgeFaces[curPMeshEdge][faceI];

                if (findIndex(faceLabels_, curFace) == -1)
                {
                    label polyPatchID =
                        pMesh.boundaryMesh().whichPatch(curFace);

                    if (polyPatchID != -1)
                    {
                        ngbPolyPatch[edgeI] = polyPatchID;
                    }
                }
            }
        }

        // Count ngb processorPolyPatches
        labelHashSet processorPatchSet;
        forAll(ngbPolyPatch, edgeI)
        {
            if (ngbPolyPatch[edgeI] != -1)
            {
                if
                (
                    isA<processorPolyPatch>
                    (
                        pMesh.boundaryMesh()[ngbPolyPatch[edgeI]]
                    )
                )
                {
                    if (!processorPatchSet.found(ngbPolyPatch[edgeI]))
                    {
                        processorPatchSet.insert(ngbPolyPatch[edgeI]);
                    }
                }
            }
        }
        labelList processorPatches(processorPatchSet.toc());
        faPatches.setSize(faPatches.size() + processorPatches.size());

        for (label i = 0; i < processorPatches.size(); i++)
        {
            DynamicList<label> tmpLst;

            forAll(ngbPolyPatch, eI)
            {
                if (ngbPolyPatch[eI] == processorPatches[i])
                {
                    tmpLst.append(allUndefEdges[eI]);
                }
            }

            faPatches[faPatchNames.size() + i].edgeLabels_ = tmpLst;

            faPatches[faPatchNames.size() + i].name_ =
                pMesh.boundaryMesh()[processorPatches[i]].name();

            faPatches[faPatchNames.size() + i].type_ =
                processorFaPatch::typeName;

            faPatches[faPatchNames.size() + i].ngbPolyPatchID_ =
                processorPatches[i];
        }

        // Remaining undefined edges
        DynamicList<label> undefEdges;
        forAll(ngbPolyPatch, eI)
        {
            if (ngbPolyPatch[eI] == -1)
            {
                undefEdges.append(allUndefEdges[eI]);
            }
            else if
            (
               !isA<processorPolyPatch>
                (
                    pMesh.boundaryMesh()[ngbPolyPatch[eI]]
                )
            )
            {
                undefEdges.append(allUndefEdges[eI]);
            }
        }

        if (undefEdges.size() > 0)
        {
            label pI = faPatches.size()-1;

            faPatches[pI].name_ = "undefined";
            faPatches[pI].type_ = "patch";
            faPatches[pI].edgeLabels_ = undefEdges;
        }
        else
        {
            faPatches.setSize(faPatches.size() - 1);
        }
    }
    else
    {
        faPatches.setSize(faPatches.size() - 1);
    }


    // Reorder processorFaPatch using
    // ordering of ngb processorPolyPatch
    forAll(faPatches, patchI)
    {
        if (faPatches[patchI].type_ == processorFaPatch::typeName)
        {
            labelList ngbFaces(faPatches[patchI].edgeLabels_.size(), -1);

            forAll(ngbFaces, edgeI)
            {
                label curEdge = faPatches[patchI].edgeLabels_[edgeI];

                label curPMeshEdge = meshEdges[curEdge];

                forAll(edgeFaces[curPMeshEdge], faceI)
                {
                    label curFace = edgeFaces[curPMeshEdge][faceI];

                    label curPatchID =
                        pMesh.boundaryMesh().whichPatch(curFace);

                    if (curPatchID == faPatches[patchI].ngbPolyPatchID_)
                    {
                        ngbFaces[edgeI] = curFace;
                    }
                }
            }

            SortableList<label> sortedNgbFaces(ngbFaces);
            labelList reorderedEdgeLabels(ngbFaces.size(), -1);
            for (label i = 0; i < reorderedEdgeLabels.size(); i++)
            {
                reorderedEdgeLabels[i] =
                    faPatches[patchI].edgeLabels_
                    [
                        sortedNgbFaces.indices()[i]
                    ];
            }

            faPatches[patchI].edgeLabels_ = reorderedEdgeLabels;
        }
    }


    // Add good patches to faMesh
    DynamicList<faPatch*> faPatchLst;

    for (label pI = 0; pI < faPatches.size(); pI++)
    {
        faPatches[pI].dict_.add("type", faPatches[pI].type_);
        faPatches[pI].dict_.add("edgeLabels", faPatches[pI].edgeLabels_);
        faPatches[pI].dict_.add
        (
            "ngbPolyPatchIndex",
            faPatches[pI].ngbPolyPatchID_
        );

        if (faPatches[pI].type_ == processorFaPatch::typeName)
        {
            if (faPatches[pI].ngbPolyPatchID_ == -1)
            {
                FatalErrorInFunction
                    << "ngbPolyPatch is not defined for processorFaPatch: "
                    << faPatches[pI].name_
                    << abort(FatalError);
            }

            const processorPolyPatch& procPolyPatch =
                refCast<const processorPolyPatch>
                (
                    pMesh.boundaryMesh()[faPatches[pI].ngbPolyPatchID_]
                );

            faPatches[pI].dict_.add("myProcNo", procPolyPatch.myProcNo());
            faPatches[pI].dict_.add
            (
                "neighbProcNo",
                procPolyPatch.neighbProcNo()
            );
        }

        faPatchLst.append
        (
            faPatch::New
            (
                faPatches[pI].name_,
                faPatches[pI].dict_,
                pI,
                boundary()
            ).ptr()
        );
    }

    addFaPatches(List<faPatch*>(faPatchLst));

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();

    if (isFile(mesh().time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
}


// Construct from polyPatch
Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const label polyPatchID
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, UpdateableMeshObject, faMesh>
    (
        pMesh,
        mesh().boundaryMesh()[polyPatchID].name()
    ),
    edgeInterpolation(*this),
    data(pMesh.thisDb()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        labelList(pMesh.boundaryMesh()[polyPatchID].size(), -1)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    registry_
    (
        pMesh.thisDb().subRegistry
        (
            "faMesh-"+mesh().boundaryMesh()[polyPatchID].name(), true
        )
    ),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    phiPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    if (debug)
    {
        Info<< "faMesh::faMesh(...) : "
            << "Creating faMesh from polyPatch" << endl;
    }

    // Set faceLabels
    forAll(faceLabels_, faceI)
    {
        faceLabels_[faceI] =
            mesh().boundaryMesh()[polyPatchID].start() + faceI;
    }

    // Add one faPatch
    const indirectPrimitivePatch& bp = patch();

    const label nTotalEdges = bp.nEdges();

    const label nInternalEdges = bp.nInternalEdges();

    labelList edgeLabels(nTotalEdges - nInternalEdges, -1);

    forAll(edgeLabels, edgeI)
    {
        edgeLabels[edgeI] = nInternalEdges + edgeI;
    }

    dictionary patchDict;

    patchDict.add("type", "patch");
    patchDict.add("edgeLabels", edgeLabels);
    patchDict.add("ngbPolyPatchIndex", -1);

    List<faPatch*> faPatchLst(1);

    faPatchLst[0] =
        faPatch::New("default", patchDict, 0, boundary()).ptr();

    addFaPatches(faPatchLst);

    setPrimitiveMeshData();

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}

Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const word& patchName
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, UpdateableMeshObject, faMesh>(pMesh, patchName),
    edgeInterpolation(*this),
    data(static_cast<const objectRegistry&>(pMesh)),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        List<label>(0)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    registry_(pMesh.thisDb().subRegistry("faMesh-"+patchName, true)),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    phiPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    const label& polyPatchID = mesh().boundaryMesh().findPatchID(patchName);

    if (polyPatchID == -1)
    {
        FatalErrorInFunction
            << "there is no boundary patch with name"<<nl
            << patchName <<nl
            << "available boundary patches are" <<nl
            << mesh().boundaryMesh().names()
            << abort(FatalError);
    }
    faceLabels_.setSize(mesh().boundaryMesh()[polyPatchID].size(), -1);

    forAll(faceLabels_, faceI)
    {
        faceLabels_[faceI] =
            mesh().boundaryMesh()[polyPatchID].start() + faceI;
    }

    const label& nTotalEdges = patch().nEdges();
    const label& nInternalEdges = patch().nInternalEdges();

    const labelListList& edgeFaces = mesh().edgeFaces();

    labelList faceCells(faceLabels_.size(), -1);

    forAll(faceCells, faceI)
    {
        label faceID = faceLabels_[faceI];

        faceCells[faceI] = mesh().faceOwner()[faceID];
    }

    labelList meshEdges =
        patch().meshEdges
        (
            mesh().edges(),
            mesh().cellEdges(),
            faceCells
        );

    labelList bndEdgeFaPatchIDs(nTotalEdges - nInternalEdges, -1);
    labelList isAProcEdge(nTotalEdges - nInternalEdges, 0);

    labelList bndEdgeList(isAProcEdge.size(), -1);
    wordList boundaryPatchNames(nTotalEdges - nInternalEdges);
    wordList boundaryPatchTypes(nTotalEdges - nInternalEdges);

    for (label edgeI = nInternalEdges; edgeI < nTotalEdges; edgeI++)
    {
        const label& curMeshEdge = meshEdges[edgeI];
        bndEdgeList[edgeI-nInternalEdges] = curMeshEdge;

        forAll(edgeFaces[curMeshEdge], faceI)
        {
            const label& curFace = edgeFaces[curMeshEdge][faceI];
            const label& curPatchID = mesh().boundaryMesh().whichPatch(curFace);

            if (curPatchID!=-1)
            {
                if (curPatchID != polyPatchID)
                {
                    const word& pName = mesh().boundaryMesh()[curPatchID].name();
                    boundaryPatchNames[edgeI-nInternalEdges] = pName;

                    const word& pType = mesh().boundaryMesh()[curPatchID].type();
                    boundaryPatchTypes[edgeI-nInternalEdges] = pType;
                }

                const polyPatch& pp = mesh().boundaryMesh()[curPatchID];
                if (pp.coupled())
                {
                    bndEdgeFaPatchIDs[edgeI-nInternalEdges] = curPatchID;
                    isAProcEdge[edgeI-nInternalEdges] = 1;
                }
            }
        }
    }
    // correct hanging Edges Patches
    //Sync edges that belong to proc patch but in the finite area mesh are boundary
    //edges

    syncTools::syncEdgeList
    (
        mesh(),
        bndEdgeList,
        isAProcEdge,
        minEqOp<label>(),
        label(-1)
    );

//    label nProcHang = 0;

    DynamicList<label> hangingEdgesLocal;

    forAll(isAProcEdge, eI)
    {
        if (isAProcEdge[eI]==-1)
        {
//            nProcHang++;
            hangingEdgesLocal.append(nInternalEdges+eI);
           }
    }
    if (Pstream::parRun())
    {
        PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

           const polyPatch& aPatch = mesh().boundaryMesh()[polyPatchID];

        forAll(mesh().boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchI];
            if (pp.coupled())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(mesh().boundaryMesh()[patchI]);

                DynamicList<label> procEdges;
                forAll(hangingEdgesLocal, eI)
                {
                    const label& eL = hangingEdgesLocal[eI];
                    if (bndEdgeFaPatchIDs[eL - nInternalEdges]==patchI)
                    {
                         edge dummyEdge = aPatch.edges()[eL];

                         label p0 = dummyEdge[0];
                         label p1 = dummyEdge[1];
                         //Point in the indirect primitive patch of faMesh
                         //global fvMesh Point
                         label gPoint0 = aPatch.meshPoints()[p0];
                         label gPoint1 = aPatch.meshPoints()[p1];

                         //local point in the processor Patch
                         label procP0 = pp.whichPoint(gPoint0);
                         label procP1 = pp.whichPoint(gPoint1);

                         //local edge in the processor Patch
                         dummyEdge[0] = procP0;
                         dummyEdge[1] = procP1;

                         label edgeLabel = pp.whichEdge(dummyEdge);

                         procEdges.append(procPatch.nbrEdges()[edgeLabel]);
                    }
                }

                labelList sendE = procEdges;
                UOPstream toNeighbProc(procPatch.neighbProcNo(), pBufs);

                toNeighbProc
                    <<sendE;
            }
        }

        pBufs.finishedSends();
        labelListList nbrEdges(mesh().boundaryMesh().size());

        forAll(mesh().boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchI];

            if (pp.coupled())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(mesh().boundaryMesh()[patchI]);

                UIPstream fromNeighbProc(procPatch.neighbProcNo(), pBufs);

                fromNeighbProc
                    >> nbrEdges[patchI];
            }
        }

        pBufs.finishedSends();
        forAll(mesh().boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchI];

            if (pp.coupled())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(mesh().boundaryMesh()[patchI]);

                wordList pNames(nbrEdges[patchI].size());
                wordList pTypes(nbrEdges[patchI].size());
                labelList edgeLabels(pTypes.size());
                forAll(nbrEdges[patchI], eI)
                {
                    const label& eL = nbrEdges[patchI][eI];
                    const label& curMeshEdge = pp.meshEdges()[eL];
                    forAll(mesh().edgeFaces()[curMeshEdge], faceI)
                    {
                        const label& curFace = edgeFaces[curMeshEdge][faceI];
                        const label& curPatchID = mesh().boundaryMesh().whichPatch(curFace);

                        if (curPatchID!=-1 && curPatchID != patchI)
                        {
                            const word pName = mesh().boundaryMesh().names()[curPatchID];
                            pNames[eI] = pName;
                            const word pType = mesh().boundaryMesh().types()[curPatchID];
                            pTypes[eI] = pType;
                        }
                    }
                    edgeLabels[eI] = procPatch.nbrEdges()[eL];

                }

                UOPstream toNeighbProc(procPatch.neighbProcNo(), pBufs);
                toNeighbProc
                    << pNames
                    << pTypes
                    << edgeLabels;
            }
        }
        pBufs.finishedSends();

        forAll(mesh().boundaryMesh(), patchI)
        {
            const polyPatch& pp = mesh().boundaryMesh()[patchI];

            if (pp.coupled())
            {
                const processorPolyPatch& procPatch =
                    refCast<const processorPolyPatch>(mesh().boundaryMesh()[patchI]);

                wordList pNamesSent;
                wordList pTypesSent;
                labelList edgeLabelsSent;

                UIPstream fromNeighbProc(procPatch.neighbProcNo(), pBufs);

                fromNeighbProc
                    >> pNamesSent
                    >> pTypesSent
                    >> edgeLabelsSent;


                forAll(edgeLabelsSent, eI)
                {
                       const label& eL = edgeLabelsSent[eI];
                      edge dummyEdge = pp.edges()[eL];

                       label p0 = dummyEdge[0];
                       label p1 = dummyEdge[1];
                       //Point in the indirect primitive patch of faMesh
                       //global fvMesh Point
                       label gPoint0 = pp.meshPoints()[p0];
                       label gPoint1 = pp.meshPoints()[p1];

                       //local point in the processor Patch
                       label procP0 = aPatch.whichPoint(gPoint0);
                       label procP1 = aPatch.whichPoint(gPoint1);

                       //local edge in the processor Patch
                       dummyEdge[0] = procP0;
                       dummyEdge[1] = procP1;

                       label edgeLabel = aPatch.whichEdge(dummyEdge);

                       boundaryPatchNames[edgeLabel-nInternalEdges] = pNamesSent[eI];
                       boundaryPatchTypes[edgeLabel-nInternalEdges] = pTypesSent[eI];
                       bndEdgeFaPatchIDs[edgeLabel-nInternalEdges] = -1;
                }
               }
        }
        pBufs.finishedSends();
    }
    DynamicList<word> faPatchNames;
    DynamicList<word> faPatchTypes;

    if (boundaryPatchNames.size()>0)
    {
        faPatchNames.append(boundaryPatchNames[0]);
        faPatchTypes.append(boundaryPatchTypes[0]);
        for (int i=1; i < boundaryPatchNames.size(); i++)
        {
            label counter = 0;
            for (int j=0; j<faPatchNames.size(); j++)
            {
                if (boundaryPatchNames[i]!=faPatchNames[j])
                {
                    counter++;
                }
            }
            if (counter==faPatchNames.size())
            {
                   faPatchNames.append(boundaryPatchNames[i]);
                   faPatchTypes.append(boundaryPatchTypes[i]);
            }
        }
    }

    forAll(mesh().boundaryMesh(), patchI)
    {
        label counter = 0;
        for (int i=0; i < faPatchNames.size(); i++)
        {
            if (mesh().boundaryMesh().names()[patchI]!=faPatchNames[i])
            {
                counter++;
            }
            if (counter==faPatchNames.size())
            {
                faPatchNames.append(mesh().boundaryMesh().names()[patchI]);
                if
                (
                  mesh().boundaryMesh().types()[patchI]== "symmetry"
                )
                {
                    faPatchTypes.append
                    (mesh().boundaryMesh().types()[patchI]);
                }
                else if
                (
                  mesh().boundaryMesh().types()[patchI]== "processor"
                )
                {
                    faPatchTypes.append
                    (mesh().boundaryMesh().types()[patchI]);
                }
                else if
                (
                  mesh().boundaryMesh().types()[patchI]== "empty"
                )
                {
                    faPatchTypes.append
                    (mesh().boundaryMesh().types()[patchI]);
                }
                else
                {
                    faPatchTypes.append
                    ("patch");
                }
            }
        }
    }
    //If a processor does not have a patch in the FA mesh
    if (!faPatchNames.size())
    {
        forAll(mesh().boundaryMesh(), patchI)
        {
            faPatchNames.append(mesh().boundaryMesh().names()[patchI]);
            if
            (
                mesh().boundaryMesh().types()[patchI]== "symmetry"
            )
            {
                faPatchTypes.append
                (mesh().boundaryMesh().types()[patchI]);
            }
            else if
            (
                mesh().boundaryMesh().types()[patchI]== "processor"
            )
            {
                faPatchTypes.append
                (mesh().boundaryMesh().types()[patchI]);
            }
            else if
            (
                mesh().boundaryMesh().types()[patchI]== "empty"
            )
            {
                faPatchTypes.append
                (mesh().boundaryMesh().types()[patchI]);
            }
            else
            {
                faPatchTypes.append
                ("patch");
            }
        }
    }

    List<faPatchData> faPatches(faPatchNames.size());
    forAll(faPatchNames, patchI)
    {
        faPatches[patchI].name_ = faPatchNames[patchI];

        faPatches[patchI].type_ = faPatchTypes[patchI];
        if
        (
            faPatches[patchI].type_ == "wall"
         || faPatches[patchI].type_ == "mappedWall"
        )
        {
            faPatches[patchI].type_ = "patch";
        }

        faPatches[patchI].ownPolyPatchID_ =
            mesh().boundaryMesh().findPatchID(patchName);

        faPatches[patchI].ngbPolyPatchID_  =
            mesh().boundaryMesh().findPatchID
            (
                faPatchNames[patchI]
            );
    }

    for (label edgeI = nInternalEdges; edgeI < nTotalEdges; edgeI++)
    {
        label curMeshEdge = meshEdges[edgeI];

        labelList curEdgePatchIDs(2, label(-1));

        label patchI = -1;

        forAll(edgeFaces[curMeshEdge], faceI)
        {
            label curFace = edgeFaces[curMeshEdge][faceI];

            label curPatchID = mesh().boundaryMesh().whichPatch(curFace);

            if (curPatchID != -1)
            {
                curEdgePatchIDs[++patchI] = curPatchID;
            }
        }

        for (label pI = 0; pI < faPatches.size(); pI++)
        {
            if
            (
                (
                    curEdgePatchIDs[0] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[1] == faPatches[pI].ngbPolyPatchID_
                )
             ||
                (
                    curEdgePatchIDs[1] == faPatches[pI].ownPolyPatchID_
                 && curEdgePatchIDs[0] == faPatches[pI].ngbPolyPatchID_
                )
            )
            {
                bndEdgeFaPatchIDs[edgeI - nInternalEdges] = pI;
                break;
            }
        }
    }
    //Correct hanging edges
    forAll(hangingEdgesLocal, eI)
    {
        const label& eL = hangingEdgesLocal[eI];
        const word& pN = boundaryPatchNames[eL-nInternalEdges];
        for (label pI = 0; pI < faPatches.size(); pI++)
        {
            if (pN == faPatches[pI].name_)
            {
                bndEdgeFaPatchIDs[eL-nInternalEdges] = pI;
            }
        }
    }

    // Set edgeLabels for each faPatch
    for (label pI = 0; pI < faPatches.size(); pI++)
    {
        DynamicList<label> tmpList;

        forAll(bndEdgeFaPatchIDs, eI)
        {
            if (bndEdgeFaPatchIDs[eI] == pI)
            {
                tmpList.append(nInternalEdges + eI);
            }
        }
        faPatches[pI].edgeLabels_ = tmpList;
    }
    // Reorder processorFaPatch using
    // ordering of ngb processorPolyPatch

    forAll(faPatches, patchI)
    {
        if (faPatches[patchI].type_ == processorFaPatch::typeName)
        {
            labelList ngbFaces(faPatches[patchI].edgeLabels_.size(), -1);

            forAll(ngbFaces, edgeI)
            {
                label curEdge = faPatches[patchI].edgeLabels_[edgeI];

                label curPMeshEdge = meshEdges[curEdge];

                forAll(edgeFaces[curPMeshEdge], faceI)
                {
                    label curFace = edgeFaces[curPMeshEdge][faceI];

                    label curPatchID =
                        mesh().boundaryMesh().whichPatch(curFace);

                    if (curPatchID == faPatches[patchI].ngbPolyPatchID_)
                    {
                        ngbFaces[edgeI] = curFace;
                    }
                }
            }
            SortableList<label> sortedNgbFaces(ngbFaces);
            labelList reorderedEdgeLabels(ngbFaces.size(), -1);
            for (label i = 0; i < reorderedEdgeLabels.size(); i++)
            {
                reorderedEdgeLabels[i] =
                    faPatches[patchI].edgeLabels_
                    [
                        sortedNgbFaces.indices()[i]
                    ];
            }
            faPatches[patchI].edgeLabels_ = reorderedEdgeLabels;
        }
    }

    // Add good patches to faMesh
    DynamicList<faPatch*> faPatchLst;

    for (label pI = 0; pI < faPatches.size(); pI++)
    {
        faPatches[pI].dict_.add("type", faPatches[pI].type_);
        faPatches[pI].dict_.add("edgeLabels", faPatches[pI].edgeLabels_);
        faPatches[pI].dict_.add
        (
            "ngbPolyPatchIndex",
            faPatches[pI].ngbPolyPatchID_
        );

        if (faPatches[pI].type_ == processorFaPatch::typeName)
        {
            const processorPolyPatch& procPolyPatch =
                refCast<const processorPolyPatch>
                (
                    pMesh.boundaryMesh()[faPatches[pI].ngbPolyPatchID_]
                );

            faPatches[pI].dict_.add("myProcNo", procPolyPatch.myProcNo());
            faPatches[pI].dict_.add
            (
                "neighbProcNo",
                procPolyPatch.neighbProcNo()
            );
        }

        faPatchLst.append
        (
            faPatch::New
            (
                faPatches[pI].name_,
                faPatches[pI].dict_,
                pI,
                boundary()
            ).ptr()
        );
    }

    addFaPatches(List<faPatch*>(faPatchLst));


    if (isFile(mesh().time().timePath()/"S0"))
    {
        S0Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S0",
                time().timeName(),
                meshSubDir,
                mesh(),
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );
    }
}

Foam::faMesh::faMesh
(
    const polyMesh& pMesh,
    const wordReList& patchNames
)
:
    GeoMesh<polyMesh>(pMesh),
    MeshObject<polyMesh, UpdateableMeshObject, faMesh>(pMesh),
    edgeInterpolation(*this),
    data(pMesh.thisDb()),
    faceLabels_
    (
        IOobject
        (
            "faceLabels",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        List<label>(0)
    ),
    boundary_
    (
        IOobject
        (
            "faBoundary",
            mesh().facesInstance(),
            meshSubDir,
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        *this,
        0
    ),
    registry_(pMesh.thisDb()),
    patchPtr_(nullptr),
    lduPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    SPtr_(nullptr),
    S0Ptr_(nullptr),
    S00Ptr_(nullptr),
    patchStartsPtr_(nullptr),
    LePtr_(nullptr),
    magLePtr_(nullptr),
    centresPtr_(nullptr),
    edgeCentresPtr_(nullptr),
    faceAreaNormalsPtr_(nullptr),
    edgeAreaNormalsPtr_(nullptr),
    pointAreaNormalsPtr_(nullptr),
    faceCurvaturesPtr_(nullptr),
    edgeTransformTensorsPtr_(nullptr),
    correctPatchPointNormalsPtr_(nullptr),
    phiPtr_(nullptr),
    globalMeshDataPtr_(nullptr)
{
    labelList patchIDs(patchNames.size(), -1);

    label size = 0;

    forAll(patchNames, patchI)
    {
        patchIDs[patchI] =
            mesh().boundaryMesh().findPatchID(patchNames[patchI]);

        if (patchIDs[patchI] < 0)
        {
            FatalErrorInFunction
                << "Patch " << patchNames[patchI] << " does not exist"
                << exit(FatalError);
        }

        size += mesh().boundaryMesh()[patchIDs[patchI]].size();
    }

    faceLabels_ = labelList(size, -1);

    sort(patchIDs);

    label faceI = -1;

    forAll(patchNames, patchI)
    {
        label start = mesh().boundaryMesh()[patchIDs[patchI]].start();

        label size  = mesh().boundaryMesh()[patchIDs[patchI]].size();

        for (label i = 0; i < size; i++)
        {
            faceLabels_[++faceI] = start + i;
        }
    }

    const indirectPrimitivePatch& bP = patch();
    const label nTotalEdges = bP.nEdges();
    const label nInternalEdges = bP.nInternalEdges();
    const globalMeshData& pd = mesh().globalData();
    const labelList& coupledEdges = pd.coupledPatchMeshEdges();
    const labelListList& edgeFaces = mesh().edgeFaces();

    boolList isAFaMeshFace(mesh().faces().size(), false);
    forAll(faceLabels_, fI)
    {
        isAFaMeshFace[faceLabels_[fI]] = true;
    }

    boolList procEdges(mesh().nEdges(), false);
    forAll(coupledEdges, cEI)
    {
        label edgeI = coupledEdges[cEI];
        procEdges[edgeI] = true;
    }

    labelList faceCells(faceLabels_.size(), -1);

    forAll(faceCells, faceI)
    {
        label faceID = faceLabels_[faceI];

        faceCells[faceI] = mesh().faceOwner()[faceID];
    }

    labelList meshEdges =
        bP.meshEdges
        (
            mesh().edges(),
            mesh().cellEdges(),
            faceCells
        );

    labelList bndEdgeFaPatchIDs(nTotalEdges - nInternalEdges, -1);
    labelList isAProcEdge(nTotalEdges - nInternalEdges, 0);
    labelList bndEdgeList(isAProcEdge.size(), -1);

    enum patchType {processorPatch = 1, fixedValue = -2, symmetry = -3, empty = -4};

    for (label edgeI = nInternalEdges; edgeI < nTotalEdges; edgeI++)
    {
        const label& curMeshEdge = meshEdges[edgeI];
        bndEdgeList[edgeI-nInternalEdges] = curMeshEdge;

        forAll(edgeFaces[curMeshEdge], faceI)
        {
            label curFace = edgeFaces[curMeshEdge][faceI];

            label curPatchID = mesh().boundaryMesh().whichPatch(curFace);

            if (curPatchID != -1 && !isAFaMeshFace[curFace])
            {
                const polyPatch& pp = mesh().boundaryMesh()[curPatchID];
                if (pp.coupled())
                {
                    bndEdgeFaPatchIDs[edgeI-nInternalEdges] = curPatchID;
                    isAProcEdge[edgeI-nInternalEdges] = 1;
                }
                else if (mesh().boundaryMesh().types()[curPatchID]=="symmetryPlane"
                        || mesh().boundaryMesh().types()[curPatchID]=="symmetry")
                {
                    bndEdgeFaPatchIDs[edgeI-nInternalEdges] = symmetry;
                }
                else if (mesh().boundaryMesh().types()[curPatchID]=="empty")
                {
                    bndEdgeFaPatchIDs[edgeI-nInternalEdges] = empty;
                }
                else
                {
                    bndEdgeFaPatchIDs[edgeI-nInternalEdges] = fixedValue;
                }
            }
        }
    }
    //Sync edges that belong to proc patch but in the finite area mesh() are boundary
    //edges

    syncTools::syncEdgeList
    (
        mesh(),
        bndEdgeList,
        isAProcEdge,
        minEqOp<label>(),
        label(0)
    );

    forAll(bndEdgeList, edgeI)
    {
        label patch = bndEdgeFaPatchIDs[edgeI];
        //check which processor patches after the sync became boundary edges
        if
        (
            patch != fixedValue
         && patch != symmetry
         && patch != empty
         && mesh().boundaryMesh()[patch].coupled()
         && !isAProcEdge[edgeI])
        {
            bndEdgeFaPatchIDs[edgeI] = -2;
        }
    }

    //Find number of total faPathces
    labelList sL = bndEdgeFaPatchIDs;
       DynamicList<label> faPatchID;
    if (sL.size()>0)
    {
        sort(sL);
        faPatchID.append(sL[0]);
        for (int i=0; i < sL.size()-1; i++)
        {
            if (sL[i]<sL[i+1])
            {
                faPatchID.append(sL[i+1]);
            }
        }
    }
    else
    {
        faPatchID.append(-2);
    }
    List<faPatchData> faPatches(faPatchID.size());

    //build faPatchData

    forAll(faPatches, pI)
    {
        SLList<label> tmpList;

        forAll(bndEdgeFaPatchIDs, eI)
        {
            if (bndEdgeFaPatchIDs[eI]==faPatchID[pI])
            {
                tmpList.append(nInternalEdges+eI);
            }
        }
        faPatches[pI].edgeLabels_ = tmpList;
        if (faPatchID[pI] == fixedValue)
        {
            faPatches[pI].name_ = "constBounds";
            faPatches[pI].type_ = "patch";
        }
        else if (faPatchID[pI] == symmetry)
        {
            faPatches[pI].name_ = "symmetryPlane";
            faPatches[pI].type_ = "symmetry";
        }
        else if (faPatchID[pI] == empty)
        {
            faPatches[pI].name_ = "emptyPatches";
            faPatches[pI].type_ = "empty";
        }
        else
        {
               faPatches[pI].name_ = mesh().boundaryMesh()[faPatchID[pI]].name();
            faPatches[pI].type_ = processorFaPatch::typeName;
            faPatches[pI].ngbPolyPatchID_ = faPatchID[pI];
        }
    }

    // Reorder processorFaPatch using
    // ordering of ngb processorPolyPatch
    forAll(faPatches, patchI)
    {
        if (faPatches[patchI].type_ == processorFaPatch::typeName)
        {
            labelList ngbFaces(faPatches[patchI].edgeLabels_.size(), -1);

            forAll(ngbFaces, edgeI)
            {
                label curEdge = faPatches[patchI].edgeLabels_[edgeI];

                label curPMeshEdge = meshEdges[curEdge];

                forAll(edgeFaces[curPMeshEdge], faceI)
                {
                    label curFace = edgeFaces[curPMeshEdge][faceI];

                    label curPatchID =
                        mesh().boundaryMesh().whichPatch(curFace);

                    if (curPatchID == faPatches[patchI].ngbPolyPatchID_)
                    {
                        ngbFaces[edgeI] = curFace;
                    }
                }
            }
            SortableList<label> sortedNgbFaces(ngbFaces);
            labelList reorderedEdgeLabels(ngbFaces.size(), -1);
            for (label i = 0; i < reorderedEdgeLabels.size(); i++)
            {
                reorderedEdgeLabels[i] =
                    faPatches[patchI].edgeLabels_
                    [
                        sortedNgbFaces.indices()[i]
                    ];
            }
            faPatches[patchI].edgeLabels_ = reorderedEdgeLabels;
        }
    }

    List<faPatch*> faPatchLst(faPatchID.size());

    forAll(faPatchLst, pI)
    {
        faPatches[pI].dict_.add("type", faPatches[pI].type_);
        faPatches[pI].dict_.add("edgeLabels", faPatches[pI].edgeLabels_);
        if (faPatches[pI].type_ == processorFaPatch::typeName)
        {
            const processorPolyPatch& procPatch
                = refCast<const processorPolyPatch>(mesh().boundaryMesh()[faPatchID[pI]]);
            faPatches[pI].dict_.add("myProcNo", procPatch.myProcNo());
            faPatches[pI].dict_.add
            (
                "neighbProcNo",
                procPatch.neighbProcNo()
            );
            faPatches[pI].dict_.add
            (
                "ngbPolyPatchIndex",
                faPatchID[pI]
            );
        }
        else if (faPatches[pI].type_ == symmetryFaPatch::typeName)
        {
            faPatches[pI].dict_.add
            (
                "ngbPolyPatchIndex",
                100001
            );
        }
        else if (faPatches[pI].type_ == emptyFaPatch::typeName)
        {
            faPatches[pI].dict_.add
            (
                "ngbPolyPatchIndex",
                100002
            );
        }
        else
        {
            faPatches[pI].dict_.add
            (
                "ngbPolyPatchIndex",
                100000
            );
        }

        faPatchLst[pI] =
            faPatch::New
            (
                faPatches[pI].name_,
                faPatches[pI].dict_,
                pI,
                boundary()
            ).ptr();
    }

    //
    Info<<"add Patches"<<endl;
    addFaPatches(faPatchLst);

}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faMesh::~faMesh()
{
    clearOut();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::faMesh::meshDir() const
{
    return mesh().dbDir()/meshSubDir;
}


const Foam::Time& Foam::faMesh::time() const
{
    return mesh().time();
}


const Foam::fileName& Foam::faMesh::pointsInstance() const
{
    return mesh().pointsInstance();
}


const Foam::fileName& Foam::faMesh::facesInstance() const
{
    return mesh().facesInstance();
}


const Foam::indirectPrimitivePatch& Foam::faMesh::patch() const
{
    if (!patchPtr_)
    {
        patchPtr_ = new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh().faces(),
                faceLabels_
            ),
            mesh().points()
        );
    }

    return *patchPtr_;
}


Foam::indirectPrimitivePatch& Foam::faMesh::patch()
{
    if (!patchPtr_)
    {
        patchPtr_ = new indirectPrimitivePatch
        (
            IndirectList<face>
            (
                mesh().faces(),
                faceLabels_
            ),
            mesh().points()
        );
    }

    return *patchPtr_;
}


const Foam::pointField& Foam::faMesh::points() const
{
    return patch().localPoints();
}


const Foam::edgeList& Foam::faMesh::edges() const
{
    return edges_;
}


const Foam::faceList& Foam::faMesh::faces() const
{
    return patch().localFaces();
}


void Foam::faMesh::addFaPatches(const List<faPatch*>& p)
{
    if (debug)
    {
        Info<< "void faMesh::addFaPatches(const List<faPatch*>& p) : "
            << "Adding patches to faMesh" << endl;
    }

    if (boundary().size() > 0)
    {
        FatalErrorInFunction
            << "boundary already exists"
            << abort(FatalError);
    }

    boundary_.setSize(p.size());

    forAll(p, patchI)
    {
        boundary_.set(patchI, p[patchI]);
    }

    setPrimitiveMeshData();

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }
    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    //Have to trigger the geometry calculation simultaneously
    //In the finiteArea Not all the processors have processor patches!!!
    calcLe();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


const Foam::objectRegistry& Foam::faMesh::thisDb() const
{
    return registry_;
}


const Foam::faBoundaryMesh& Foam::faMesh::boundary() const
{
    return boundary_;
}


const Foam::labelList& Foam::faMesh::patchStarts() const
{
    if (!patchStartsPtr_)
    {
        calcPatchStarts();
    }

    return *patchStartsPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::Le() const
{
    if (!LePtr_)
    {
        calcLe();
    }

    return *LePtr_;
}


const Foam::edgeScalarField& Foam::faMesh::magLe() const
{
    if (!magLePtr_)
    {
        calcMagLe();
    }

    return *magLePtr_;
}


const Foam::areaVectorField& Foam::faMesh::areaCentres() const
{
    if (!centresPtr_)
    {
        calcAreaCentres();
    }

    return *centresPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::edgeCentres() const
{
    if (!edgeCentresPtr_)
    {
        calcEdgeCentres();
    }

    return *edgeCentresPtr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S() const
{
    if (!SPtr_)
    {
        calcS();
    }

    return *SPtr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S0() const
{
    if (!S0Ptr_)
    {
        FatalErrorInFunction
            << "S0 is not available"
            << abort(FatalError);
    }

    return *S0Ptr_;
}


const Foam::DimensionedField<Foam::scalar, Foam::areaMesh>&
Foam::faMesh::S00() const
{
    if (!S00Ptr_)
    {
        S00Ptr_ = new DimensionedField<scalar, areaMesh>
        (
            IOobject
            (
                "S00",
                time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            S0()
        );

        S0Ptr_->writeOpt() = IOobject::AUTO_WRITE;
    }

    return *S00Ptr_;
}

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::areaMesh>>
Foam::faMesh::Ssc() const
{
    if (moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar tFrac =
        (
            ts.value() - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (tFrac < (1 - SMALL))
        {
            return S0() + tFrac*(S() - S0());
        }
        else
        {
            return S();
        }
    }
    else
    {
        return S();
    }
}

Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::areaMesh>>
Foam::faMesh::Ssc0() const
{
    if (moving() && time().subCycling())
    {
        const TimeState& ts = time();
        const TimeState& ts0 = time().prevTimeState();

        scalar t0Frac =
        (
            (ts.value() - ts.deltaTValue())
          - (ts0.value() - ts0.deltaTValue())
        )/ts0.deltaTValue();

        if (t0Frac > SMALL)
        {
            return S0() + t0Frac*(S() - S0());
        }
        else
        {
            return S0();
        }
    }
    else
    {
        return S0();
    }
}



const Foam::areaVectorField& Foam::faMesh::faceAreaNormals() const
{
    if (!faceAreaNormalsPtr_)
    {
        calcFaceAreaNormals();
    }

    return *faceAreaNormalsPtr_;
}


const Foam::edgeVectorField& Foam::faMesh::edgeAreaNormals() const
{
    if (!edgeAreaNormalsPtr_)
    {
        calcEdgeAreaNormals();
    }

    return *edgeAreaNormalsPtr_;
}


const Foam::vectorField& Foam::faMesh::pointAreaNormals() const
{
    if (!pointAreaNormalsPtr_)
    {
        calcPointAreaNormals();

        if (quadricsFit_ > 0)
        {
            calcPointAreaNormalsByQuadricsFit();
        }
    }

    return *pointAreaNormalsPtr_;
}


const Foam::areaScalarField& Foam::faMesh::faceCurvatures() const
{
    if (!faceCurvaturesPtr_)
    {
        calcFaceCurvatures();
    }

    return *faceCurvaturesPtr_;
}


const Foam::FieldField<Foam::Field, Foam::tensor>&
Foam::faMesh::edgeTransformTensors() const
{
    if (!edgeTransformTensorsPtr_)
    {
        calcEdgeTransformTensors();
    }

    return *edgeTransformTensorsPtr_;
}


// Return global mesh data
const Foam::faGlobalMeshData& Foam::faMesh::globalData() const
{
    if (!globalMeshDataPtr_)
    {
        globalMeshDataPtr_ = new faGlobalMeshData(*this);
    }

    return *globalMeshDataPtr_;
}


const Foam::lduAddressing& Foam::faMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        calcLduAddressing();
    }

    return *lduPtr_;
}


bool Foam::faMesh::movePoints()
{
   // Grab point motion from polyMesh
    const vectorField& newPoints = mesh().points();

    // Grab old time areas if the time has been incremented
    if (curTimeIndex_ < time().timeIndex())
    {
        storeOldSurf(S());
    }
    if (!phiPtr_)
    {
        // Create mesh motion flux
        phiPtr_ = new edgeScalarField
        (
            IOobject
            (
                "faMeshPhi",
                this->time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimArea/dimTime
        );
    }
    else
    {
        // Grab old time mesh motion fluxes if the time has been incremented
        if (phiPtr_->timeIndex() != time().timeIndex())
        {
            phiPtr_->oldTime();
        }
    }
    edgeScalarField& phi = *phiPtr_;

    scalar rDeltaT = 1.0/time().deltaTValue();

    tmp<scalarField> tsweptAreas = sweptAreas(newPoints, points());
    scalarField& sweptAreas = tsweptAreas.ref();


    phi.primitiveFieldRef() = scalarField::subField(sweptAreas, nInternalEdges());
    phi.primitiveFieldRef() *= rDeltaT;

    const faPatchList& patches = boundary();

    forAll(boundary(), pI)
    {
        phi.boundaryFieldRef()[pI] = patches[pI].patchSlice(sweptAreas);
        phi.boundaryFieldRef()[pI] *= rDeltaT;
    }

    updateGeomNotOldAreas();

    // Move boundary points
    const_cast<faBoundaryMesh&>(boundary_).movePoints(newPoints);

    // Move interpolation
    const edgeInterpolation& cei = *this;
    const_cast<edgeInterpolation&>(cei).edgeInterpolation::movePoints();

    return true;
}

const Foam::edgeScalarField& Foam::faMesh::phi() const
{
    if (!phiPtr_)
    {
        FatalErrorInFunction
            << "mesh flux field does not exist, is the mesh actually moving?"
            << abort(FatalError);
    }

    // Set zero current time
    // mesh motion fluxes if the time has been incremented
    if (phiPtr_->timeIndex() != time().timeIndex())
    {
        (*phiPtr_) = dimensionedScalar("0", dimArea/dimTime, 0.0);
    }

    return *phiPtr_;
}


bool Foam::faMesh::correctPatchPointNormals(const label patchID) const
{
    if ((patchID < 0) || (patchID >= boundary().size()))
    {
        FatalErrorInFunction
            << "patchID is not in the valid range"
            << abort(FatalError);
    }

    if (correctPatchPointNormalsPtr_)
    {
        return (*correctPatchPointNormalsPtr_)[patchID];
    }

    return false;
}


//- Set patch point normals corrections
Foam::boolList& Foam::faMesh::correctPatchPointNormals() const
{
    if (!correctPatchPointNormalsPtr_)
    {
        correctPatchPointNormalsPtr_ =
            new boolList(boundary().size(), false);
    }

    return *correctPatchPointNormalsPtr_;
}


const Foam::faSchemes& Foam::faMesh::schemes() const
{
    if (!faSchemes_.valid())
    {
        faSchemes_ = new faSchemes(mesh());
    }

    return faSchemes_;
}


const Foam::faSolution& Foam::faMesh::solution() const
{
    if (!faSolution_.valid())
    {
        faSolution_ = new faSolution(mesh());
    }

    return faSolution_;
}


Foam::faSchemes& Foam::faMesh::schemes()
{
    if (!faSchemes_.valid())
    {
        faSchemes_ = new faSchemes(mesh());
    }

    return *faSchemes_;
}


Foam::faSolution& Foam::faMesh::solution()
{
    if (!faSolution_.valid())
    {
        faSolution_ = new faSolution(mesh());
    }

    return *faSolution_;
}


bool Foam::faMesh::write(const bool valid) const
{
    bool ok = true;
    if (phiPtr_)
    {
        ok = phiPtr_->write();
    }
    faceLabels_.write();
    boundary_.write();

    return ok;
}

const Foam::Vector<Foam::label>& Foam::faMesh::solutionD() const
{
    return mesh().solutionD();
}
// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::faMesh::operator!=(const faMesh& m) const
{
    return &m != this;
}


bool Foam::faMesh::operator==(const faMesh& m) const
{
    return &m == this;
}


// ************************************************************************* //
