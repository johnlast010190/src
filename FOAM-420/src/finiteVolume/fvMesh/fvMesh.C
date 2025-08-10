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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2016-2017 OpenCFD Ltd.
    (c) 2010-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/volFields/slicedVolFields.H"
#include "fields/surfaceFields/slicedSurfaceFields.H"
#include "fields/Fields/Field/SubField.H"
#include "include/demandDrivenData.H"
#include "fvMesh/fvMeshLduAddressing.H"
#include "meshes/polyMesh/mapPolyMesh/mapPolyMesh.H"
#include "interpolation/mapping/fvFieldMappers/MapFvFields.H"
#include "fvMesh/fvMeshMapper/fvMeshMapper.H"
#include "fields/cloud/mapClouds.H"
#include "meshes/MeshObject/MeshObject.H"
#include "fvMatrices/fvMatrix/fvMatrix.H"
#include "VectorN/finiteVolume/fields/volFields/volVectorNFields.H"
#include "containers/Lists/CompactListList/CompactListList.H"
#include "fvMesh/fvMeshStitchers/fvMeshStitcher/fvMeshStitcher.H"
#include "fvMesh/fvMeshTopoChangers/fvMeshTopoChanger/fvMeshTopoChanger.H"
#include "fvMesh/fvMeshMovers/fvMeshMover/fvMeshMover.H"
#include "nonConformal/nonConformalFuncs/nonConformalFuncs.H"
#include "fvMesh/fvPatches/derived/nonConformalOrig/nonConformalOrigFvPatch.H"
#include "fvMesh/fvPatches/constraint/nonConformal/nonConformalFvPatch.H"
#include "fields/fvsPatchFields/basic/nonConformalCalculated/nonConformalCalculatedFvsPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::clearGeomNotOldVol()
{
    if (debug)
    {
        InfoInFunction << "clearGeomNotOldVol" << endl;
    }

    meshObject::clearUpto
    <
        fvMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    meshObject::clearUpto
    <
        lduMesh,
        GeometricMeshObject,
        MoveableMeshObject
    >(*this);

    deleteDemandDrivenData(VPtr_);
    deleteDemandDrivenData(SfSlicePtr_);
    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfSlicePtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CSlicePtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfSlicePtr_);
    deleteDemandDrivenData(CfPtr_);
}


void Foam::fvMesh::updateGeomNotOldVol()
{
    bool haveV = (VPtr_ != nullptr);
    bool haveSf = (SfSlicePtr_ != nullptr || SfPtr_ != nullptr);
    bool haveMagSf = (magSfSlicePtr_ != nullptr || magSfPtr_ != nullptr);
    bool haveCP = (CSlicePtr_ != nullptr || CPtr_ != nullptr);
    bool haveCf = (CfSlicePtr_ != nullptr || CfPtr_ != nullptr);

    clearGeomNotOldVol();

    // Now recreate the fields
    if (haveV)
    {
        (void)V();
    }
    if (haveSf)
    {
        (void)Sf();
    }
    if (haveMagSf)
    {
        (void)magSf();
    }
    if (haveCP)
    {
        (void)C();
    }
    if (haveCf)
    {
        (void)Cf();
    }
}


void Foam::fvMesh::clearGeom()
{
    if (debug)
    {
        InfoInFunction << "Clearing geometric data" << endl;
    }

    clearGeomNotOldVol();

    deleteDemandDrivenData(phiPtr_);
    deleteDemandDrivenData(V0Ptr_);
    deleteDemandDrivenData(V00Ptr_);

    // Mesh motion flux cannot be deleted here because the old-time flux
    // needs to be saved.
}


void Foam::fvMesh::clearAddressing(const bool isMeshUpdate)
{
    if (debug)
    {
        InfoInFunction << "isMeshUpdate: " << isMeshUpdate << endl;
    }

    if (isMeshUpdate)
    {
        // Part of a mesh update. Keep meshObjects that have an updateMesh
        // callback
        meshObject::clearUpto
        <
            fvMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
        meshObject::clearUpto
        <
            lduMesh,
            TopologicalMeshObject,
            UpdateableMeshObject
        >
        (
            *this
        );
    }
    else
    {
        meshObject::clear<fvMesh, TopologicalMeshObject>(*this);
        meshObject::clear<lduMesh, TopologicalMeshObject>(*this);
    }

    deleteDemandDrivenData(lduPtr_);
    deleteDemandDrivenData(polyFacesBfPtr_);
    deleteDemandDrivenData(polyBFaceOffsetsPtr_);
    deleteDemandDrivenData(polyBFaceOffsetPatchesPtr_);
    deleteDemandDrivenData(polyBFaceOffsetPatchFacesPtr_);
    deleteDemandDrivenData(polyBFacePatchesPtr_);
    deleteDemandDrivenData(polyBFacePatchFacesPtr_);
}


void Foam::fvMesh::storeOldVol(const scalarField& V)
{
    if (curTimeIndex_ < time().timeIndex())
    {
        if (debug)
        {
            InfoInFunction
                << " Storing old time volumes since from time "
                << curTimeIndex_ << " and time now " << time().timeIndex()
                << " V:" << V.size() << endl;
        }

        if (V00Ptr_ && V0Ptr_)
        {
            // Copy V0 into V00 storage
            *V00Ptr_ = *V0Ptr_;
        }

        if (V0Ptr_)
        {
            // Copy V into V0 storage
            V0Ptr_->scalarField::operator=(V);
        }
        else
        {
            // Allocate V0 storage, fill with V
            V0Ptr_ = new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "V0",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                *this,
                dimVolume
            );
            scalarField& V0 = *V0Ptr_;
            // Note: V0 now sized with current mesh, not with (potentially
            //       different size) V.
            V0.setSize(V.size());
            V0 = V;
        }

        curTimeIndex_ = time().timeIndex();

        if (debug)
        {
            InfoInFunction
                << " Stored old time volumes V0:" << V0Ptr_->size()
                << endl;
            if (V00Ptr_)
            {
                InfoInFunction
                    << " Stored oldold time volumes V00:" << V00Ptr_->size()
                    << endl;
            }
        }
    }
}


Foam::wordList Foam::fvMesh::polyFacesPatchTypes() const
{
    wordList wantedPatchTypes
    (
        boundary().size(),
        calculatedFvsPatchLabelField::typeName
    );

    forAll(boundary(), patchi)
    {
        const fvPatch& fvp = boundary()[patchi];

        if (isA<nonConformalFvPatch>(fvp))
        {
            wantedPatchTypes[patchi] =
                nonConformalCalculatedFvsPatchLabelField::typeName;
        }
    }

    return wantedPatchTypes;
}


Foam::surfaceLabelField::Boundary& Foam::fvMesh::polyFacesBfRef()
{
    if (!polyFacesBfPtr_)
    {
        polyFacesBf();
    }

    return *polyFacesBfPtr_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMesh::fvMesh(const IOobject& io, const bool changers)
:
    polyMesh(io),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    stitcher_(fvMeshStitcher::New(*this, changers)),
    topoChanger_(nullptr),
    mover_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvMesh from IOobject" << endl;
    }

    if (changers && needsStitching())
    {
        // Construct non-conformal couples on-the-fly
        nonConformalFuncs::createNonConformalCouples(*this);
    }
    else
    {
        // Stitch or re-stitch if necessary
        stitcher_->connect(false, changers, true);
    }

    // Construct changers
    if (changers)
    {
        topoChanger_.set(fvMeshTopoChanger::New(*this).ptr());
        mover_.set(fvMeshMover::New(*this).ptr());
    }

    // Check the existence of the cell volumes and read if present
    // and set the storage of V00
    if (fileHandler().isFile(time().timePath()/"V0"))
    {
        V0Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "V0",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this
        );

        V00();
    }

    // Check the existence of the mesh fluxes and read if present
    if (fileHandler().isFile(time().timePath()/"meshPhi"))
    {
        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                true
            ),
            *this,
            true,
            false
        );
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<labelList>& allOwner,
    const Xfer<labelList>& allNeighbour,
    const bool syncPar
)
:
    polyMesh(io, points, faces, allOwner, allNeighbour, syncPar),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    mover_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvMesh from components" << endl;
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<cellList>& cells,
    const bool syncPar,
    const bool autoWrite
)
:
    polyMesh(io, points, faces, cells, syncPar, autoWrite),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    mover_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvMesh from components" << endl;
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const cellShapeList& shapes,
    const faceListList& boundaryFaces,
    const wordList& boundaryPatchNames,
    const PtrList<dictionary>& boundaryDicts,
    const word& defaultBoundaryPatchName,
    const word& defaultBoundaryPatchType,
    const bool syncPar
)
:
    polyMesh
    (
       io,
       points,
       shapes,
       boundaryFaces,
       boundaryPatchNames,
       boundaryDicts,
       defaultBoundaryPatchName,
       defaultBoundaryPatchType,
       syncPar
    ),
    surfaceInterpolation(*this),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this,boundaryMesh()),
    stitcher_(nullptr),
    topoChanger_(nullptr),
    mover_(nullptr),
    lduPtr_(nullptr),
    polyFacesBfPtr_(nullptr),
    polyBFaceOffsetsPtr_(nullptr),
    polyBFaceOffsetPatchesPtr_(nullptr),
    polyBFaceOffsetPatchFacesPtr_(nullptr),
    polyBFacePatchesPtr_(nullptr),
    polyBFacePatchFacesPtr_(nullptr),
    curTimeIndex_(time().timeIndex()),
    VPtr_(nullptr),
    V0Ptr_(nullptr),
    V00Ptr_(nullptr),
    SfSlicePtr_(nullptr),
    SfPtr_(nullptr),
    magSfSlicePtr_(nullptr),
    magSfPtr_(nullptr),
    CSlicePtr_(nullptr),
    CPtr_(nullptr),
    CfSlicePtr_(nullptr),
    CfPtr_(nullptr),
    phiPtr_(nullptr)
{
    if (debug)
    {
        InfoInFunction << "Constructing fvMesh from blockMesh" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMesh::~fvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMesh::addFvPatches
(
    const List<polyPatch*> & p,
    const bool validBoundary
)
{
    if (boundary().size())
    {
        FatalErrorInFunction
            << " boundary already exists"
            << abort(FatalError);
    }

    // first add polyPatches
    addPatches(p, validBoundary);
    boundary_.addPatches(boundaryMesh());
}


Foam::polyMesh::readUpdateState Foam::fvMesh::readUpdate()
{
    if (debug)
    {
        InfoInFunction << "Updating fvMesh.  ";
    }

    polyMesh::readUpdateState state = polyMesh::readUpdate();

    if (stitcher_.valid() && state != polyMesh::UNCHANGED)
    {
        stitcher_->disconnect(false, false);
    }

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        if (debug)
        {
            Info<< "Boundary and topological update" << endl;
        }

        boundary_.readUpdate(boundaryMesh());

        clearOut();

    }
    else if (state == polyMesh::TOPO_CHANGE)
    {
        if (debug)
        {
            Info<< "Topological update" << endl;
        }

        clearOut();
    }
    else if (state == polyMesh::POINTS_MOVED)
    {
        if (debug)
        {
            Info<< "Point motion update" << endl;
        }

        clearGeom();
    }
    else
    {
        if (debug)
        {
            Info<< "No update" << endl;
        }
    }

    if (stitcher_.valid() && state != polyMesh::UNCHANGED)
    {
        stitcher_->connect(false, false, true);
    }

    return state;
}


const Foam::lduAddressing& Foam::fvMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new fvMeshLduAddressing(*this);
    }

    return *lduPtr_;
}


Foam::SolverPerformance<Foam::scalar> Foam::fvMesh::solve
(
    fvMatrix<scalar>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::vector> Foam::fvMesh::solve
(
    fvMatrix<vector>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::sphericalTensor> Foam::fvMesh::solve
(
    fvMatrix<sphericalTensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::symmTensor> Foam::fvMesh::solve
(
    fvMatrix<symmTensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::SolverPerformance<Foam::tensor> Foam::fvMesh::solve
(
    fvMatrix<tensor>& m,
    const dictionary& dict
) const
{
    // Redirect to fvMatrix solver
    return m.solveSegregatedOrCoupled(dict);
}


Foam::IOobject Foam::fvMesh::polyFacesBfIO(const IOobject::readOption r) const
{
    return
        IOobject
        (
            "polyFaces",
            pointsInstance(),
            typeName,
            *this,
            r,
            IOobject::NO_WRITE,
            false
        );
}


const Foam::surfaceLabelField::Boundary& Foam::fvMesh::polyFacesBf() const
{
    if (!polyFacesBfPtr_)
    {
        polyFacesBfPtr_ =
            new surfaceLabelField::Boundary
            (
                boundary(),
                surfaceLabelField::null(),
                polyFacesPatchTypes(),
                boundaryMesh().types()
            );

        forAll(boundary(), patchi)
        {
            const polyPatch& pp = boundaryMesh()[patchi];
            (*polyFacesBfPtr_)[patchi] =
                labelList(identity(pp.size()) + pp.start());
        }
    }

    return *polyFacesBfPtr_;
}


const Foam::CompactListList<Foam::label>&
Foam::fvMesh::polyBFacePatches() const
{
    if (!polyBFacePatchesPtr_)
    {
        const label nPolyBFaces = nFaces() - nInternalFaces();

        // Count face-poly-bFaces to get the offsets
        polyBFaceOffsetsPtr_ = new labelList(nPolyBFaces + 1, 0);
        labelList& offsets = *polyBFaceOffsetsPtr_;
        forAll(boundary(), patchi)
        {
            if (!isA<indirectPolyPatch>(boundary()[patchi].patch()))
            {
                forAll(boundary()[patchi], patchFacei)
                {
                    const label polyBFacei =
                        (
                            polyFacesBfPtr_
                        ? (*polyFacesBfPtr_)[patchi][patchFacei]
                        : boundary()[patchi].start() + patchFacei
                        )
                    - nInternalFaces();
                    offsets[polyBFacei + 1] ++;
                }
            }
        }
        for (label polyBFacei = 0; polyBFacei < nPolyBFaces; ++ polyBFacei)
        {
            offsets[polyBFacei + 1] += offsets[polyBFacei];
        }

        // Set the poly-bFace patches and patch-faces, using the offsets as
        // counters
        polyBFaceOffsetPatchesPtr_ = new labelList(offsets.last());
        polyBFaceOffsetPatchFacesPtr_ = new labelList(offsets.last());
        labelUList& patches = *polyBFaceOffsetPatchesPtr_;
        labelUList& patchFaces = *polyBFaceOffsetPatchFacesPtr_;
        forAll(boundary(), patchi)
        {
            if (!isA<indirectPolyPatch>(boundary()[patchi].patch()))
            {
                forAll(boundary()[patchi], patchFacei)
                {
                    const label polyBFacei =
                        (
                            polyFacesBfPtr_
                            ? (*polyFacesBfPtr_)[patchi][patchFacei]
                            : boundary()[patchi].start() + patchFacei
                        )
                        - nInternalFaces();
                    patches[offsets[polyBFacei]] = patchi;
                    patchFaces[offsets[polyBFacei]] = patchFacei;
                    offsets[polyBFacei] ++;
                }
            }
        }

        // Restore the offsets by removing the count
        for
        (
            label polyBFacei = nPolyBFaces - 1;
            polyBFacei >= 0;
            -- polyBFacei
        )
        {
            offsets[polyBFacei + 1] = offsets[polyBFacei];
        }
        offsets[0] = 0;

        // List-lists
        polyBFacePatchesPtr_ =
            new CompactListList<label>(offsets, patches);
        polyBFacePatchFacesPtr_ =
            new CompactListList<label>(offsets, patchFaces);
    }

    return *polyBFacePatchesPtr_;
}


const Foam::CompactListList<Foam::label>&
Foam::fvMesh::polyBFacePatchFaces() const
{
    if (!polyBFacePatchFacesPtr_)
    {
        polyBFacePatches();
    }

    return *polyBFacePatchFacesPtr_;
}


const Foam::fvMeshStitcher& Foam::fvMesh::stitcher() const
{
    return stitcher_();
}


const Foam::fvMeshTopoChanger& Foam::fvMesh::topoChanger() const
{
    return topoChanger_();
}


const Foam::fvMeshMover& Foam::fvMesh::mover() const
{
    return mover_();
}


bool Foam::fvMesh::hasChangers() const
{
    return topoChanger_.valid() || mover_.valid();
}


bool Foam::fvMesh::needsStitching() const
{
    bool hasNonConformalCouples = false;

    forAll(boundary(), patchi)
    {
        if (isA<nonConformalOrigFvPatch>(boundary()[patchi]))
        {
            hasNonConformalCouples = true;
            break;
        }
    }

    IOobject nccDictIO
    (
        "nonConformalCouplesDict",
        this->time().system(),
        this->dbDir(),
        *this,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE,
        false
    );

    if (nccDictIO.typeHeaderOk<IOdictionary>())
    {
        IOdictionary nccDict(nccDictIO);
        if (nccDict.toc().size())
        {
            hasNonConformalCouples = true;
        }
        else if (hasNonConformalCouples)
        {
            FatalErrorInFunction
                << "Non-conformal original patches found in the mesh,"
                << " but empty 'nonConformalCouplesDict'. present"
                << exit(FatalError);
        }
    }
    else
    {
        if (hasNonConformalCouples)
        {
            FatalErrorInFunction
                << "Non-conformal original patches found in the mesh,"
                << " but cannot find file 'nonConformalCouplesDict'."
                << exit(FatalError);
        }
    }

    return hasNonConformalCouples && !stitcher_->stitches();
}


bool Foam::fvMesh::dynamic() const
{
    return topoChanger_->dynamic() || mover_->dynamic();
}


bool Foam::fvMesh::update()
{
    if (!conformal()) stitcher_->disconnect(true, true);

    bool updated = topoChanger_->update();

    updated = move() || updated;

    return updated;
}


bool Foam::fvMesh::move()
{
    if (!conformal()) stitcher_->disconnect(true, true);

    const bool moved = mover_->update();

    stitcher_->connect(true, true, false);

    return moved;
}


void Foam::fvMesh::clearMeshPhi()
{
    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);
}


void Foam::fvMesh::clearOut()
{
    clearGeom();

    surfaceInterpolation::clearOut();

    clearAddressing();

    polyMesh::clearOut();
}


void Foam::fvMesh::updateMesh(const mapPolyMesh& mpm)
{
    // Update polyMesh. This needs to keep volume existent!
    polyMesh::updateMesh(mpm);

    if (VPtr_ && mpm.hasOldCellVolumes())
    {
        // Grab old time volumes if the time has been incremented
        // This will update V0, V00
        storeOldVol(mpm.oldCellVolumes());

        // Few checks
        if (VPtr_ && (VPtr_->size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V:" << VPtr_->size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
        if (V0Ptr_ && (V0Ptr_->size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V0:" << V0Ptr_->size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
        if (V00Ptr_ && (V00Ptr_->size() != mpm.nOldCells()))
        {
            FatalErrorInFunction
                << "V0:" << V00Ptr_->size()
                << " not equal to the number of old cells "
                << mpm.nOldCells()
                << exit(FatalError);
        }
    }

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Map all fields
    mapFields(mpm);

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    meshObject::updateMesh<fvMesh>(*this, mpm);
    meshObject::updateMesh<lduMesh>(*this, mpm);

    if (mover_.valid())
    {
        mover_->updateMesh(mpm);

        // Reset the old-time cell volumes prior to mesh-motion
        if (V0Ptr_)
        {
            *V0Ptr_ = V();
        }
    }
}


void Foam::fvMesh::updateGIB()
{
    polyMesh::updateGIB();

    UpdateGIBFields<scalar, fvPatchField, volMesh>(*this);
    UpdateGIBFields<vector, fvPatchField, volMesh>(*this);
    UpdateGIBFields<sphericalTensor, fvPatchField, volMesh>(*this);
    UpdateGIBFields<symmTensor, fvPatchField, volMesh>(*this);
    UpdateGIBFields<tensor, fvPatchField, volMesh>(*this);

    UpdateGIBFields<scalar, fvsPatchField, surfaceMesh>(*this);
    UpdateGIBFields<vector, fvsPatchField, surfaceMesh>(*this);
    UpdateGIBFields<sphericalTensor, fvsPatchField, surfaceMesh>(*this);
    UpdateGIBFields<symmTensor, fvsPatchField, surfaceMesh>(*this);
    UpdateGIBFields<tensor, fvsPatchField, surfaceMesh>(*this);

    updateGeomNotOldVol();
    surfaceInterpolation::clearOut();
    clearAddressing(true);
    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);
}


void Foam::fvMesh::storeGIBFields()
{
    polyMesh::updateGIB();
    StoreGIBFields<scalar, fvPatchField, volMesh>(*this);
    StoreGIBFields<vector, fvPatchField, volMesh>(*this);
    StoreGIBFields<sphericalTensor, fvPatchField, volMesh>(*this);
    StoreGIBFields<symmTensor, fvPatchField, volMesh>(*this);
    StoreGIBFields<tensor, fvPatchField, volMesh>(*this);

    StoreGIBFields<scalar, fvsPatchField, surfaceMesh>(*this);
    StoreGIBFields<vector, fvsPatchField, surfaceMesh>(*this);
    StoreGIBFields<sphericalTensor, fvsPatchField, surfaceMesh>(*this);
    StoreGIBFields<symmTensor, fvsPatchField, surfaceMesh>(*this);
    StoreGIBFields<tensor, fvsPatchField, surfaceMesh>(*this);

    updateGeomNotOldVol();
}


Foam::tmp<Foam::scalarField> Foam::fvMesh::movePoints(const pointField& p)
{
    // Grab old time volumes if the time has been incremented
    // This will update V0, V00
    if (curTimeIndex_ < time().timeIndex())
    {
        storeOldVol(V());
    }

    if (calcSolverQties())
    {
        if (!phiPtr_)
        {
            // Create mesh motion flux
            phiPtr_ = new surfaceScalarField
            (
                IOobject
                (
                    "meshPhi",
                    this->time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    true
                ),
                *this,
                dimVolume/dimTime
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
    }

     // Move the polyMesh and set the mesh motion fluxes to the swept-volumes
    tmp<scalarField> tsweptVols = polyMesh::movePoints(p);

    if (calcSolverQties())
    {
        scalarField& sweptVols = tsweptVols.ref();
        surfaceScalarField& phi = *phiPtr_;

        scalar rDeltaT = 1.0/time().deltaTValue();
        phi.primitiveFieldRef() =
            scalarField::subField(sweptVols, nInternalFaces());
        phi.primitiveFieldRef() *= rDeltaT;

        const fvPatchList& patches = boundary();

        surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();

        forAll(patches, patchi)
        {
            phibf[patchi] = patches[patchi].patchSlice(sweptVols);
            phibf[patchi] *= rDeltaT;
        }
    }

    // Update or delete the local geometric properties as early as possible so
    // they can be used if necessary. These get recreated here instead of
    // demand driven since they might do parallel transfers which can conflict
    // with when they're actually being used.
    // Note that between above "polyMesh::movePoints(p)" and here nothing
    // should use the local geometric properties.
    updateGeomNotOldVol();

    // Update other local data
    boundary_.movePoints();
    surfaceInterpolation::movePoints();

    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);

    return tsweptVols;
}


Foam::tmp<Foam::scalarField> Foam::fvMesh::moveGIBPoints(const pointField& p)
{
    if (V0Ptr_)
    {
        // Copy V into V0 storage
        V0Ptr_->scalarField::operator=(V());
    }
    else
    {
        // Allocate V0 storage, fill with V
        V0Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "V0",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            *this,
            dimVolume
        );
        scalarField& V0 = *V0Ptr_;
        // Note: V0 now sized with current mesh, not with (potentially
        //       different size) V.
        V0.setSize(V().size());
        V0 = V();
    }

    if (!phiPtr_)
    {
        // Create mesh motion flux
       phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimVolume/dimTime
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

    surfaceScalarField& phi = *phiPtr_;

    // Move the polyMesh and set the mesh motion fluxes to the swept-volumes

    scalar rDeltaT = 1.0/time().deltaTValue();

    tmp<scalarField> tsweptVols = polyMesh::movePoints(p);
    scalarField& sweptVols = tsweptVols.ref();

    // Clear geometry
    clearGeomNotOldVol();

    phi.primitiveFieldRef() =
        scalarField::subField(sweptVols, nInternalFaces());
    phi.primitiveFieldRef() *= rDeltaT;

    boundary_.movePoints();

    surfaceInterpolation::clearOut();

    clearAddressing(true);

    meshObject::movePoints<fvMesh>(*this);
    meshObject::movePoints<lduMesh>(*this);
    return tsweptVols;
}


void Foam::fvMesh::conform(const surfaceScalarField& phi)
{
    // Clear the geometry fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    // Modify the mesh fluxes, if necessary
    if (notNull(phi) && phiPtr_)
    {
        for (label i = 0; i <= phi.nOldTimes(); i++)
        {
            setPhi().oldTime(i).forceAssign(phi.oldTime(i));
        }
    }
}


void Foam::fvMesh::unconform
(
    const surfaceLabelField::Boundary& polyFacesBf,
    const surfaceVectorField& Sf,
    const surfaceVectorField& Cf,
    const surfaceScalarField& phi,
    const bool sync
)
{
    // Clear the geometry fields
    clearGeomNotOldVol();

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    // Clear any non-updateable addressing
    clearAddressing(true);

    // Create non-sliced copies of geometry fields
    SfRef();
    magSfRef();
    CRef();
    CfRef();

    // Set the topology
    polyFacesBfRef().forceAssign(polyFacesBf);

    // Set the face geometry
    SfRef().forceAssign(Sf);
    magSfRef().forceAssign(mag(Sf));
    CRef().boundaryFieldRef().forceAssign(Cf.boundaryField());
    CfRef().forceAssign(Cf);

    // Communicate processor-coupled cell geometry. Cell-centre processor patch
    // fields must contain the (transformed) cell-centre locations on the other
    // side of the coupling. This is so that non-conformal patches can
    // construct weights and deltas without reference to the poly mesh
    // geometry.
    //
    // Note that the initEvaluate/evaluate communication does a transformation,
    // but it is wrong in this instance. A vector field gets transformed as if
    // it were a displacement, but the cell-centres need a positional
    // transformation. That's why there's the un-transform and re-transform bit
    // below just after the evaluate call.
    //
    // This transform handling is a bit of a hack. It would be nicer to have a
    // field attribute which identifies a field as needing a positional
    // transformation, and for it to apply automatically within the coupled
    // patch field. However, at the moment, the cell centres field is the only
    // vol-field containing an absolute position, so the hack is functionally
    // sufficient for now.
    if (sync && (Pstream::parRun() || !time().processorCase()))
    {
        volVectorField::Boundary& CBf = CRef().boundaryFieldRef();

        const label nReq = Pstream::nRequests();

        forAll(CBf, patchi)
        {
            if (isA<processorFvPatch>(CBf[patchi].patch()))
            {
                CBf[patchi].initEvaluate(Pstream::defaultCommsType);
            }
        }

        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(CBf, patchi)
        {
            if (isA<processorFvPatch>(CBf[patchi].patch()))
            {
                CBf[patchi].evaluate(Pstream::defaultCommsType);

                const transformer& t =
                    refCast<const processorFvPatch>(CBf[patchi].patch())
                   .transform();

                t.invTransform(CBf[patchi], CBf[patchi]);
                t.transformPosition(CBf[patchi], CBf[patchi]);
            }
        }
    }

    // Modify the mesh fluxes, if necessary
    if (notNull(phi) && phiPtr_)
    {
        for (label i = 0; i <= phi.nOldTimes(); i++)
        {
            setPhi().oldTime(i).forceAssign(phi.oldTime(i));
        }
    }
}


void Foam::fvMesh::mapFields(const mapPolyMesh& meshMap)
{
    if (debug)
    {
        InfoInFunction
            << " nOldCells:" << meshMap.nOldCells()
            << " nCells:" << nCells()
            << " nOldFaces:" << meshMap.nOldFaces()
            << " nFaces:" << nFaces()
            << endl;
    }


    // We require geometric properties valid for the old mesh
    if
    (
        meshMap.cellMap().size() != nCells()
     || meshMap.faceMap().size() != nFaces()
    )
    {
        FatalErrorInFunction
            << "mapPolyMesh does not correspond to the old mesh."
            << " nCells:" << nCells()
            << " cellMap:" << meshMap.cellMap().size()
            << " nOldCells:" << meshMap.nOldCells()
            << " nFaces:" << nFaces()
            << " faceMap:" << meshMap.faceMap().size()
            << " nOldFaces:" << meshMap.nOldFaces()
            << exit(FatalError);
    }

    // Create a mapper
    const fvMeshMapper mapper(*this, meshMap);

    // Map all the volFields in the objectRegistry
    #define mapVolFieldType(Type, nullArg)                                     \
        MapGeometricFields<Type, fvPatchField, fvMeshMapper, volMesh>(mapper);
    FOR_ALL_FIELD_TYPES(mapVolFieldType);

    //- Map vectorN-type volFields
    MapGeometricFields<vector1, fvPatchField, fvMeshMapper, volMesh>(mapper);
    MapGeometricFields<vector4, fvPatchField, fvMeshMapper, volMesh>(mapper);
    MapGeometricFields<tensor4, fvPatchField, fvMeshMapper, volMesh>(mapper);

    // Map all the surfaceFields in the objectRegistry
    #define mapSurfaceFieldType(Type, nullArg)                                 \
        MapGeometricFields<Type, fvsPatchField, fvMeshMapper, surfaceMesh>     \
        (mapper);
    FOR_ALL_FIELD_TYPES(mapSurfaceFieldType);

    // Map all the dimensionedFields in the objectRegistry
    #define mapVolInternalFieldType(Type, nullArg)                             \
        MapDimensionedFields<Type, fvMeshMapper, volMesh>(mapper);
    FOR_ALL_FIELD_TYPES(mapVolInternalFieldType);

    //- Map vectorN-type dimensionedFields
    MapDimensionedFields<vector1, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<vector4, fvMeshMapper, volMesh>(mapper);
    MapDimensionedFields<tensor4, fvMeshMapper, volMesh>(mapper);

    // Map all the clouds in the objectRegistry
    mapClouds(*this, meshMap);

    const labelList& cellMap = meshMap.cellMap();

    // Map the old volume. Just map to new cell labels.
    if (V0Ptr_)
    {
        scalarField& V0 = *V0Ptr_;

        scalarField savedV0(V0);
        V0.setSize(nCells());

        forAll(V0, i)
        {
            if (cellMap[i] > -1)
            {
                V0[i] = savedV0[cellMap[i]];
            }
            else
            {
                V0[i] = 0.0;
            }
        }

        // Inject volume of merged cells
        label nMerged = 0;
        forAll(meshMap.reverseCellMap(), oldCelli)
        {
            label index = meshMap.reverseCellMap()[oldCelli];

            if (index < -1)
            {
                label celli = -index-2;

                V0[celli] += savedV0[oldCelli];

                nMerged++;
            }
        }

        if (debug)
        {
            Info<< "Mapping old time volume V0. Merged "
                << nMerged << " out of " << nCells() << " cells" << endl;
        }
    }


    // Map the old-old volume. Just map to new cell labels.
    if (V00Ptr_)
    {
        scalarField& V00 = *V00Ptr_;

        scalarField savedV00(V00);
        V00.setSize(nCells());

        forAll(V00, i)
        {
            if (cellMap[i] > -1)
            {
                V00[i] = savedV00[cellMap[i]];
            }
            else
            {
                V00[i] = 0.0;
            }
        }

        // Inject volume of merged cells
        label nMerged = 0;
        forAll(meshMap.reverseCellMap(), oldCelli)
        {
            label index = meshMap.reverseCellMap()[oldCelli];

            if (index < -1)
            {
                label celli = -index-2;

                V00[celli] += savedV00[oldCelli];
                nMerged++;
            }
        }

        if (debug)
        {
            Info<< "Mapping old time volume V00. Merged "
                << nMerged << " out of " << nCells() << " cells" << endl;
        }
    }
}


void Foam::fvMesh::removeFvBoundary()
{
    if (debug)
    {
        InfoInFunction << "Removing boundary patches." << endl;
    }

    // Remove fvBoundaryMesh data first.
    boundary_.clear();
    boundary_.setSize(0);
    polyMesh::removeBoundary();

    clearOut();
}


Foam::tmp<Foam::vectorField> Foam::fvMesh::velocityCorrect
(
    const vectorField& pc
) const
{
    tmp<vectorField> zeroFieldt
        (
            new vectorField
            (
                pc.size(), vector::zero
            )
        );
    return zeroFieldt;
}


bool Foam::fvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp,
    const bool valid
) const
{
    bool ok = true;

    if (!conformal())
    {
        // Create a full surface field with the polyFacesBf boundary field then
        // overwrite all conformal faces with an index of -1 to save disk space

        surfaceLabelField polyFaces
        (
            polyFacesBfIO(IOobject::NO_READ),
            *this,
            dimless,
            labelField(nInternalFaces(), -1),
            *polyFacesBfPtr_
        );

        forAll(boundary(), patchi)
        {
            const fvPatch& fvp = boundary()[patchi];
            if (!isA<nonConformalFvPatch>(fvp))
            {
                polyFaces.boundaryFieldRef()[patchi] = -1;
            }
        }

        ok = ok & polyFaces.write(valid);
    }

    if (phiPtr_)
    {
        ok = ok && phiPtr_->write(valid);

        // NOTE: The old old time mesh phi might be necessary for certain
        // solver smooth restart using second order time schemes.
        //ok = phiPtr_->oldTime().write(valid);
    }

    // For second-order restarts we need to write V0
    if (V00Ptr_)
    {
        ok = ok && V0Ptr_->write(valid);
    }

    if (topoChanger_.valid())
    {
        topoChanger_->write(valid);
    }

    if (mover_.valid())
    {
        mover_->write(valid);
    }

    return ok && polyMesh::writeObject(fmt, ver, cmp, valid);
}


bool Foam::fvMesh::write(const bool valid) const
{
    return polyMesh::write(valid);
}


const Foam::fvSchemes& Foam::fvMesh::schemes() const
{
    if (!fvSchemes_.valid())
    {
        fvSchemes_ = new fvSchemes(*this);
    }

    return fvSchemes_();
}


const Foam::fvSolution& Foam::fvMesh::solution() const
{
    if (!fvSolution_.valid())
    {
        fvSolution_ = new fvSolution(*this);
    }

    return fvSolution_();
}


Foam::fvSchemes& Foam::fvMesh::schemes()
{
    if (!fvSchemes_.valid())
    {
        fvSchemes_ = new fvSchemes(*this);
    }

    return *fvSchemes_;
}


Foam::fvSolution& Foam::fvMesh::solution()
{
    if (!fvSolution_.valid())
    {
        fvSolution_ = new fvSolution(*this);
    }

    return *fvSolution_;
}


template<>
//typename Foam::pTraits<Foam::sphericalTensor>::labelType
Foam::pTraits<Foam::sphericalTensor>::labelType
Foam::fvMesh::validComponents<Foam::sphericalTensor>() const
{
    return Foam::pTraits<Foam::sphericalTensor>::labelType(1);
}


template<>
Foam::pTraits<Foam::sphericalTensor>::labelType
Foam::fvMesh::validComponents2<Foam::sphericalTensor>() const
{
    return Foam::pTraits<Foam::sphericalTensor>::labelType(1);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::fvMesh::operator!=(const fvMesh& bm) const
{
    return &bm != this;
}


bool Foam::fvMesh::operator==(const fvMesh& bm) const
{
    return &bm == this;
}


// ************************************************************************* //
