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
    (c) 2016-2020 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"
#include "fields/Fields/transformField/transformField.H"
#include "fields/Fields/Field/SubField.H"
#include "meshes/polyMesh/polyMesh.H"
#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "AMIInterpolation/faceAreaIntersect/faceAreaIntersect.H"
#include "primitives/ops/ops.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicAMIPolyPatch, dictionary);
}

const Foam::scalar Foam::cyclicAMIPolyPatch::tolerance_ = 1e-10;

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

void Foam::cyclicAMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{
    if (owner())
    {
        AMIPtr_.clear();

        const polyPatch& nbr = nbrPatch();
        pointField srcPoints(localPoints());
        pointField nbrPoints(nbr.localPoints());

        if (debug)
        {
            const Time& t = boundaryMesh().mesh().time();
            OFstream os(t.path()/name() + "_neighbourPatch-org.obj");
            meshTools::writeOBJ(os, nbrPatch().localFaces(), nbrPoints);
        }

        label patchSize0 = size();
        label nbrPatchSize0 = nbr.size();

        if (createAMIFaces_)
        {
            // AMI is created based on the original patch faces (non-extended patch)
            if (srcFaceIDs_.size())
            {
                patchSize0 = srcFaceIDs_.size();
            }
            if (tgtFaceIDs_.size())
            {
                nbrPatchSize0 = tgtFaceIDs_.size();
            }
        }

        // Transform neighbour patch to local system
        transform().transformPosition(nbrPoints, nbrPoints);

        primitivePatch nbrPatch0
        (
            SubList<face>(nbr.localFaces(), nbrPatchSize0),
            nbrPoints
        );

        primitivePatch patch0
        (
            SubList<face>(localFaces(), patchSize0),
            srcPoints
        );

        if (debug)
        {
            const Time& t = boundaryMesh().mesh().time();
            OFstream osN(t.path()/name() + "_neighbourPatch-trans.obj");
            meshTools::writeOBJ(osN, nbrPatch0.localFaces(), nbrPoints);

            OFstream osO(t.path()/name() + "_ownerPatch.obj");
            meshTools::writeOBJ(osO, this->localFaces(), localPoints());
        }

        {
            label srcSize = returnReduce(patch0.size(), sumOp<label>());
            label tgtSize = returnReduce(nbrPatch0.size(), sumOp<label>());

            Info<< "AMI: Creating addressing and weights between "
                << srcSize << " source faces (" << name() << ") and "
                << tgtSize << " target faces (" << nbrPatch().name() << ")"
                << endl;
        }


        // Construct/apply AMI interpolation to determine addressing and weights
        AMIPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                patch0,
                nbrPatch0,
                surfPtr(),
                faceAreaIntersect::tmMesh,
                AMIRequireMatch_,
                AMIMethod,
                AMILowWeightCorrection_,
                AMIReverse_,
                AMIPatchToPatchInterpolation::relTolDef_,
                AMIPatchToPatchInterpolation::tolDef_,
                AMIPatchToPatchInterpolation::cosMatchAngleDef_,
                AMIPatchToPatchInterpolation::maxAMIWeightScaleDef_,
                false
            )
        );

        if (true)
        {
            const scalarField& srcWSum = AMIPtr_->srcWeightsSum();
            label nFace = returnReduce(srcWSum.size(), sumOp<label>());
            label nLowWeight = 0;

            const scalar lowWeightThreshold = AMIPtr_->lowWeightCorrection();

            forAll(srcWSum, fI)
            {
                if (srcWSum[fI]<lowWeightThreshold) nLowWeight++;
            }

            if (nFace)
            {
                Info<< indent
                    << "AMI: Patch " << name()
                    << " sum(weights) min/max/average = "
                    << gMin(srcWSum) << ", "
                    << gMax(srcWSum) << ", "
                    << gAverage(srcWSum) << endl;
            }

            label nLow = returnReduce(nLowWeight, sumOp<label>());

            if (nLow)
            {
                Info<< indent
                    << "AMI: Patch " << name()
                    << " identified " << nLow
                    << " faces with weights less than " << lowWeightThreshold
                    << endl;
            }

            const scalarField& tgtWSum = AMIPtr_->tgtWeightsSum();
            nFace = returnReduce(tgtWSum.size(), sumOp<label>());
            nLowWeight = 0;

            forAll(tgtWSum, fI)
            {
                if (tgtWSum[fI]<lowWeightThreshold) nLowWeight++;
            }

            if (nFace)
            {
                Info<< indent
                    << "AMI: Patch " << nbrPatch().name()
                    << " sum(weights) min/max/average = "
                    << gMin(tgtWSum) << ", "
                    << gMax(tgtWSum) << ", "
                    << gAverage(tgtWSum) << endl;
            }

            nLow = returnReduce(nLowWeight, sumOp<label>());

            if (nLow)
            {
                Info<< indent
                    << "AMI: Patch " << nbrPatch().name()
                    << " identified " << nLow
                    << " faces with weights less than " << lowWeightThreshold
                    << endl;
            }
        }

        if (debug)
        {
            Pout<< "cyclicAMIPolyPatch : " << name()
                << " constructed AMI with " << nl
                << "    " << "srcAddress:" << AMIPtr_().srcAddress().size()
                << nl
                << "    " << "tgAddress :" << AMIPtr_().tgtAddress().size()
                << nl << endl;
        }
    }
}


void Foam::cyclicAMIPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    // The AMI is no longer valid. Leave it up to demand-driven calculation
    if (AMIPtr_.valid())
    {
        AMIPtr_->upToDate() = false;
        if (!createAMIFaces_)
        {
            AMIPtr_.clear();
        }
    }

    directPolyPatch::initCalcGeometry(pBufs);
}


void Foam::cyclicAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    const Time& t = this->boundaryMesh().mesh().time();

    if
    (
        !Pstream::parRun()
     && !t.processorCase()
     && !isDir(t.rootPath()/t.globalCaseName()/word("processor0"))
    )
    {
        static_cast<cyclicTransform&>(*this) =
            cyclicTransform
            (
                name(),
                this->primitivePatch::faceCentres(),
                this->primitivePatch::faceAreas(),
                *this,
                nbrPatchName(),
                nbrPatch().faceCentres(),
                nbrPatch().faceAreas(),
                nbrPatch(),
                matchTolerance()
            );
    }
    else
    {
        static_cast<cyclicTransform&>(*this) =
            cyclicTransform
            (
                name(),
                this->primitivePatch::faceAreas(),
                *this,
                nbrPatchName(),
                nbrPatch(),
                matchTolerance()
            );
    }

    // If a parallel pre-processing phase is executed for cyclicAMI patches
    // with incomplete transformation input, the missing information still
    // needs to be calculated from the patch geometry
    if (Pstream::parRun() && !transformComplete())
    {
        static_cast<cyclicTransform&>(*this) =
            cyclicTransform
            (
                name(),
                faceCentres(),
                faceAreas(),
                *this,
                nbrPatchName(),
                nbrPatch().faceCentres(),
                nbrPatch().faceAreas(),
                nbrPatch(),
                matchTolerance()
            );
    }

    // Calculate number of sectors for a rotational transformation.
    // Needed by cyclicPeriodicAMI's.
    if (transformType() == ROTATIONAL && transformComplete())
    {
        const scalar radAngle = degToRad(rotationAngle());

        if (radAngle != 0)
        {
            sectors_ = sign(radAngle)*constant::mathematical::twoPi/radAngle;
        }
        else
        {
            FatalErrorInFunction
                << "Patch " << name()
                << " computed a null angle between the periodics."
                << " Please, check setup or mesh generation."
                << exit(FatalError);
        }
    }
}


void Foam::cyclicAMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    // The AMI is no longer valid. Leave it up to demand-driven calculation

    directPolyPatch::initMovePoints(pBufs, p);

    // See below. Clear out any local geometry
    primitivePatch::clearGeom();

    // Note: processorPolyPatch::initMovePoints calls
    // processorPolyPatch::initCalcGeometry which will trigger calculation of
    // patch faceCentres() and cell volumes...

    if (createAMIFaces_)
    {
        // Note: AMI should have been updated in setTopology

        // faceAreas() and faceCentres() have been reset and will be
        // recalculated on-demand using the mesh points and no longer
        // correspond to the scaled areas!
        restoreScaledGeometry();

        // deltas need to be recalculated to use new face centres!
    }
    else
    {
        AMIPtr_.clear();
    }
}


void Foam::cyclicAMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    directPolyPatch::movePoints(pBufs, p);
}


void Foam::cyclicAMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    directPolyPatch::initUpdateMesh(pBufs);

    if (createAMIFaces_ && boundaryMesh().mesh().topoChanging() && owner())
    {
        setAMIFaces();
    }
}


void Foam::cyclicAMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    directPolyPatch::updateMesh(pBufs);
}


void Foam::cyclicAMIPolyPatch::clearGeom()
{
    AMIPtr_.clear();
    directPolyPatch::clearGeom();
}


void Foam::cyclicAMIPolyPatch::rename(const wordList& newNames)
{
    directPolyPatch::rename(newNames);
    nbrPatch().nbrPatchName_ = newNames[index()];
}


void Foam::cyclicAMIPolyPatch::reorder(const labelUList& oldToNewIndex)
{
    directPolyPatch::reorder(oldToNewIndex);
    if (nbrPatchID_ != -1)
    {
        nbrPatchID_ = oldToNewIndex[nbrPatchID_];
    }
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, size, start, index, bm, patchType),
    cyclicTransform(true),
    nbrPatchName_(word::null),
    nbrPatchID_(-1),
    AMIPtr_(nullptr),
    AMIMethod_(AMIPatchToPatchInterpolation::imFaceAreaWeight),
    AMIReverse_(false),
    AMIRequireMatch_(true),
    AMILowWeightCorrection_(0.01),
    surfPtr_(nullptr),
    surfDict_(fileName("surface")),
    createAMIFaces_(false),
    moveFaceCentres_(false),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_()
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    coupledPolyPatch(name, dict, index, bm, patchType),
    cyclicTransform(dict, true),
    nbrPatchName_(dict.lookupOrDefault<word>("neighbourPatch", "")),
    coupleGroup_(dict),
    nbrPatchID_(-1),
    AMIPtr_(nullptr),
    AMIMethod_
    (
        AMIPatchToPatchInterpolation::wordTointerpolationMethod
        (
            dict.lookupOrDefault
            (
                "method",
                AMIPatchToPatchInterpolation::interpolationMethodToWord
                (
                    AMIPatchToPatchInterpolation::imFaceAreaWeight
                )
            )
        )
    ),
    AMIReverse_(dict.lookupOrDefault<bool>("flipNormals", false)),
    AMIRequireMatch_(true),
    AMILowWeightCorrection_(dict.lookupOrDefault("lowWeightCorrection", 0.01)),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface")),
    createAMIFaces_(dict.lookupOrDefault("createAMIFaces", false)),
    moveFaceCentres_(false),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_()
{
    if (nbrPatchName_ == word::null && !coupleGroup_.valid())
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "No \"neighbourPatch\" or \"coupleGroup\" provided."
            << exit(FatalIOError);
    }

    if (nbrPatchName_ == name)
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible

    if (createAMIFaces_)
    {
        srcFaceIDs_.setSize(readLabel(dict.lookup("srcSize")));
        tgtFaceIDs_.setSize(readLabel(dict.lookup("tgtSize")));
        moveFaceCentres_ = dict.lookupOrDefault<bool>("moveFaceCentres", true);
    }
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    coupledPolyPatch(pp, bm),
    cyclicTransform(pp),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_),
    createAMIFaces_(pp.createAMIFaces_),
    moveFaceCentres_(pp.moveFaceCentres_),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_()
{
    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    coupledPolyPatch(pp, bm, index, newSize, newStart),
    cyclicTransform(pp),
    nbrPatchName_(nbrPatchName),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_),
    createAMIFaces_(pp.createAMIFaces_),
    moveFaceCentres_(pp.moveFaceCentres_),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_()
{
    if (nbrPatchName_ == name())
    {
        FatalErrorInFunction
            << "Neighbour patch name " << nbrPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Neighbour patch might not be valid yet so no transformation
    // calculation possible
}


Foam::cyclicAMIPolyPatch::cyclicAMIPolyPatch
(
    const cyclicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    coupledPolyPatch(pp, bm, index, mapAddressing, newStart),
    cyclicTransform(pp),
    nbrPatchName_(pp.nbrPatchName_),
    coupleGroup_(pp.coupleGroup_),
    nbrPatchID_(-1),
    AMIPtr_(nullptr),
    AMIMethod_(pp.AMIMethod_),
    AMIReverse_(pp.AMIReverse_),
    AMIRequireMatch_(pp.AMIRequireMatch_),
    AMILowWeightCorrection_(pp.AMILowWeightCorrection_),
    surfPtr_(nullptr),
    surfDict_(pp.surfDict_),
    createAMIFaces_(pp.createAMIFaces_),
    moveFaceCentres_(pp.moveFaceCentres_),
    srcFaceIDs_(),
    tgtFaceIDs_(),
    faceAreas0_(),
    faceCentres0_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicAMIPolyPatch::~cyclicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cyclicAMIPolyPatch::nbrPatchID() const
{
    if (nbrPatchID_ == -1)
    {
        nbrPatchID_ = this->boundaryMesh().findPatchID(nbrPatchName());

        if (nbrPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal neighbour patch name " << nbrPatchName()
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a cyclic AMI patch
        const cyclicAMIPolyPatch& nbrPatch =
            refCast<const cyclicAMIPolyPatch>
            (
                this->boundaryMesh()[nbrPatchID_]
            );

        if (nbrPatch.nbrPatchName() != name())
        {
            WarningInFunction
                << "Patch " << name()
                << " specifies neighbour patch " << nbrPatchName()
                << nl << " but that in return specifies "
                << nbrPatch.nbrPatchName() << endl;
        }
    }

    return nbrPatchID_;
}


bool Foam::cyclicAMIPolyPatch::owner() const
{
    return index() < nbrPatchID();
}


const Foam::cyclicAMIPolyPatch& Foam::cyclicAMIPolyPatch::nbrPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
    return refCast<const cyclicAMIPolyPatch>(pp);
}


const Foam::autoPtr<Foam::searchableSurface>&
Foam::cyclicAMIPolyPatch::surfPtr() const
{
    const word surfType(surfDict_.lookupOrDefault<word>("type", "none"));

    if (!surfPtr_.valid() && owner() && surfType != "none")
    {
        word surfName(surfDict_.lookupOrDefault("name", name()));

        const polyMesh& mesh = boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


const Foam::AMIPatchToPatchInterpolation& Foam::cyclicAMIPolyPatch::AMI() const
{
    if (!owner())
    {
        FatalErrorInFunction
            << "AMI interpolator only available to owner patch"
            << abort(FatalError);
    }

    if (!AMIPtr_.valid())
    {
        resetAMI(AMIMethod_);
    }

    return AMIPtr_();
}


bool Foam::cyclicAMIPolyPatch::applyLowWeightCorrection() const
{
    if (owner())
    {
        return AMI().applyLowWeightCorrection();
    }
    else
    {
        return nbrPatch().AMI().applyLowWeightCorrection();
    }
}


const Foam::scalarField& Foam::cyclicAMIPolyPatch::lowCorrAMIWeights() const
{
    if (owner())
    {
        return AMI().srcWeightsSum();
    }
    else
    {
        return nbrPatch().AMI().tgtWeightsSum();
    }
}


void Foam::cyclicAMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{}


bool Foam::cyclicAMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    faceMap.setSize(pp.size());
    faceMap = -1;

    rotation.setSize(pp.size());
    rotation = 0;

    return false;
}


Foam::label Foam::cyclicAMIPolyPatch::pointFace
(
    const label facei,
    const vector& n,
    point& p
) const
{
    point prt(transform().invTransformPosition(p));
    vector nrt(transform().invTransform(n));

    label nbrFacei = -1;

    if (owner())
    {
        nbrFacei = AMI().tgtPointFace
        (
            *this,
            nbrPatch(),
            nrt,
            facei,
            prt
        );
    }
    else
    {
        nbrFacei = nbrPatch().AMI().srcPointFace
        (
            nbrPatch(),
            *this,
            nrt,
            facei,
            prt
        );
    }

    if (nbrFacei >= 0)
    {
        p = prt;
    }

    return nbrFacei;
}


void Foam::cyclicAMIPolyPatch::write(Ostream& os) const
{
    coupledPolyPatch::write(os);

    if (!nbrPatchName_.empty())
    {
        os.writeEntry("neighbourPatch", nbrPatchName_);
    }

    coupleGroup_.write(os);

    cyclicTransform::write(os);

    if (AMIMethod_ != AMIPatchToPatchInterpolation::imFaceAreaWeight)
    {
        os.writeEntry
        (
            "method",
            AMIPatchToPatchInterpolation::interpolationMethodToWord(AMIMethod_)
        );
    }

    if (AMIReverse_)
    {
        os.writeEntry("flipNormals", AMIReverse_);
    }

    os.writeEntry("lowWeightCorrection", AMILowWeightCorrection_);

    if (!surfDict_.empty())
    {
        os.writeKeyword(surfDict_.dictName());
        os  << surfDict_;
    }

    if (createAMIFaces_)
    {
        os.writeEntry("createAMIFaces", createAMIFaces_);
        os.writeEntry("srcSize", srcFaceIDs_.size());
        os.writeEntry("tgtSize", tgtFaceIDs_.size());
        os.writeEntry("moveFaceCentres", moveFaceCentres_);
    }
}


// ************************************************************************* //
