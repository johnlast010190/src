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
    (c) 2013-2016 OpenFOAM Foundation
    (c) 2017 OpenCFD Ltd.
    (c) 2017-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/patches/cyclicACMI/cyclicACMIPolyPatch/cyclicACMIPolyPatch.H"
#include "fields/Fields/Field/SubField.H"
#include "db/Time/Time.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, dictionary);
}

const Foam::scalar Foam::cyclicACMIPolyPatch::tolerance_ = 1e-10;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::cyclicACMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod&
) const
{
    if (owner())
    {
        const polyPatch& nonOverlapPatch = this->nonOverlapPatch();

        if (debug)
        {
            Pout<< "cyclicACMIPolyPatch::resetAMI : recalculating weights"
                << " for " << name() << " and " << nonOverlapPatch.name()
                << endl;
        }

        if (boundaryMesh().mesh().hasCellCentres())
        {
            if (debug)
            {
                Pout<< "cyclicACMIPolyPatch::resetAMI : clearing cellCentres"
                    << " for " << name() << " and " << nonOverlapPatch.name()
                    << endl;
            }

            //WarningInFunction
            //    << "The mesh already has cellCentres calculated when"
            //    << " resetting ACMI " << name() << "." << endl
            //    << "This is a problem since ACMI adapts the face areas"
            //    << " (to close cells) so this has" << endl
            //    << "to be done before cell centre calculation." << endl
            //    << "This can happen if e.g. the cyclicACMI is after"
            //    << " any processor patches in the boundary." << endl;
            const_cast<polyMesh&>
            (
                boundaryMesh().mesh()
            ).primitiveMesh::clearGeom();
        }


        // Trigger re-building of faceAreas
        (void)boundaryMesh().mesh().faceAreas();


        // Calculate the AMI using partial face-area-weighted. This leaves
        // the weights as fractions of local areas (sum(weights) = 1 means
        // face is fully covered)
        cyclicAMIPolyPatch::resetAMI
        (
            AMIPatchToPatchInterpolation::imPartialFaceAreaWeight
        );

        AMIPatchToPatchInterpolation& AMI =
            const_cast<AMIPatchToPatchInterpolation&>(this->AMI());

        // Output some stats. AMIInterpolation will have already output the
        // average weights ("sum(weights) min/max/average = 1, 1, 1")
        {
            const scalarField& wghtsSum = AMI.srcWeightsSum();

            label nUncovered = 0;
            label nCovered = 0;
            forAll(wghtsSum, facei)
            {
                scalar sum = wghtsSum[facei];
                if (sum < tolerance_)
                {
                    nUncovered++;
                }
                else if (sum > scalar(1)-tolerance_)
                {
                    nCovered++;
                }
            }
            label nTotal = wghtsSum.size();
            reduce(
                std::tie(nUncovered, nCovered, nTotal),
                UniformParallelOp<sumOp<label>, 3>{}
            );

            Info<< "ACMI: Patch source uncovered/blended/covered = "
                << nUncovered << ", " << nTotal-nUncovered-nCovered
                << ", " << nCovered << endl;
        }
        {
            const scalarField& wghtsSum = AMI.tgtWeightsSum();

            label nUncovered = 0;
            label nCovered = 0;
            forAll(wghtsSum, facei)
            {
                scalar sum = wghtsSum[facei];
                if (sum < tolerance_)
                {
                    nUncovered++;
                }
                else if (sum > scalar(1)-tolerance_)
                {
                    nCovered++;
                }
            }
            label nTotal = wghtsSum.size();
            reduce(
                std::tie(nUncovered, nCovered, nTotal),
                UniformParallelOp<sumOp<label>, 3>{}
            );

            Info<< "ACMI: Patch target uncovered/blended/covered = "
                << nUncovered << ", " << nTotal-nUncovered-nCovered
                << ", " << nCovered << endl;
        }

        srcMask_ =
            min(scalar(1) - tolerance_, max(tolerance_, AMI.srcWeightsSum()));

        tgtMask_ =
            min(scalar(1) - tolerance_, max(tolerance_, AMI.tgtWeightsSum()));


        // Adapt owner side areas. Note that in uncoupled situations (e.g.
        // decomposePar) srcMask, tgtMask can be zero size.
        if (srcMask_.size())
        {
            vectorField::subField Sf = faceAreas();
            scalarField::subField magSf = magFaceAreas();
            vectorField::subField noSf = nonOverlapPatch.faceAreas();
            scalarField::subField noMagSf = nonOverlapPatch.magFaceAreas();

            forAll(Sf, facei)
            {
                Sf[facei] *= srcMask_[facei];
                magSf[facei] = mag(Sf[facei]);
                noSf[facei] *= 1.0 - srcMask_[facei];
                noMagSf[facei] = mag(noSf[facei]);
            }
        }
        // Adapt slave side areas
        if (tgtMask_.size())
        {
            const cyclicACMIPolyPatch& cp =
                refCast<const cyclicACMIPolyPatch>(this->nbrPatch());
            const polyPatch& pp = cp.nonOverlapPatch();

            vectorField::subField Sf = cp.faceAreas();
            scalarField::subField magSf = cp.magFaceAreas();
            vectorField::subField noSf = pp.faceAreas();
            scalarField::subField noMagSf = pp.magFaceAreas();

            forAll(Sf, facei)
            {
                Sf[facei] *= tgtMask_[facei];
                magSf[facei] = mag(Sf[facei]);
                noSf[facei] *= 1.0 - tgtMask_[facei];
                noMagSf[facei] = mag(noSf[facei]);
            }
        }

        // Re-normalise the weights since the effect of overlap is already
        // accounted for in the area.
        {
            scalarListList& srcWeights = AMI.srcWeights();
            scalarField& srcWeightsSum = AMI.srcWeightsSum();
            forAll(srcWeights, i)
            {
                scalarList& wghts = srcWeights[i];
                if (wghts.size())
                {
                    scalar& sum = srcWeightsSum[i];

                    forAll(wghts, j)
                    {
                        wghts[j] /= sum;
                    }
                    sum = 1.0;
                }
            }
        }
        {
            scalarListList& tgtWeights = AMI.tgtWeights();
            scalarField& tgtWeightsSum = AMI.tgtWeightsSum();
            forAll(tgtWeights, i)
            {
                scalarList& wghts = tgtWeights[i];
                if (wghts.size())
                {
                    scalar& sum = tgtWeightsSum[i];
                    forAll(wghts, j)
                    {
                        wghts[j] /= sum;
                    }
                    sum = 1.0;
                }
            }
        }

        // Set the updated flag
        updated_ = true;
    }
}


void Foam::cyclicACMIPolyPatch::initCalcGeometry(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initCalcGeometry : " << name() << endl;
    }

    cyclicAMIPolyPatch::initCalcGeometry(pBufs);

    // Initialise the AMI
    resetAMI();
}


void Foam::cyclicACMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::calcGeometry : " << name() << endl;
    }

    static_cast<cyclicTransform&>(*this) =
        cyclicTransform
        (
            name(),
            faceAreas(),
            *this,
            nbrPatchName(),
            nbrPatch(),
            matchTolerance()
        );
}


void Foam::cyclicACMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initMovePoints : " << name() << endl;
    }
    cyclicAMIPolyPatch::initMovePoints(pBufs, p);

    // Initialise the AMI
    resetAMI();
}


void Foam::cyclicACMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::movePoints : " << name() << endl;
    }
    cyclicAMIPolyPatch::movePoints(pBufs, p);
}


void Foam::cyclicACMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::initUpdateMesh : " << name() << endl;
    }
    cyclicAMIPolyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::updateMesh : " << name() << endl;
    }
    cyclicAMIPolyPatch::updateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::clearGeom()
{
    if (debug)
    {
        Pout<< "cyclicACMIPolyPatch::clearGeom : " << name() << endl;
    }
    cyclicAMIPolyPatch::clearGeom();
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::srcMask() const
{
    return srcMask_;
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::tgtMask() const
{
    return tgtMask_;
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch(name, size, start, index, bm, patchType),
    nonOverlapPatchName_(word::null),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    DeprecationWarningInFunction("cyclicACMI", "patch type", 40200)
        << "Please replace it with partially-overlapping non-conformal couples."
        << nl << endl;

    AMILowWeightCorrection_ = -1.0;
    AMIRequireMatch_ = false;

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch(name, dict, index, bm, patchType),
    nonOverlapPatchName_(dict.lookup("nonOverlapPatch")),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    DeprecationWarningInFunction("cyclicACMI", "patch type", 40200)
        << "Please replace it with partially-overlapping non-conformal couples."
        << nl << endl;

    //- correction of the low area correction default weights
    AMILowWeightCorrection_ = dict.lookupOrDefault("lowWeightCorrection", -1.0);

    AMIRequireMatch_ = false;

    if (nonOverlapPatchName_ == name)
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    AMIRequireMatch_ = false;

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& nonOverlapPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    nonOverlapPatchName_(nonOverlapPatchName),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    AMIRequireMatch_ = false;

    if (nonOverlapPatchName_ == name())
    {
        FatalErrorInFunction
            << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    updated_(false)
{
    AMIRequireMatch_ = false;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::~cyclicACMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::cyclicACMIPolyPatch& Foam::cyclicACMIPolyPatch::nbrPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[nbrPatchID()];
    return refCast<const cyclicACMIPolyPatch>(pp);
}


Foam::label Foam::cyclicACMIPolyPatch::nonOverlapPatchID() const
{
    if (nonOverlapPatchID_ == -1)
    {
        nonOverlapPatchID_ =
            this->boundaryMesh().findPatchID(nonOverlapPatchName_);

        if (nonOverlapPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal non-overlapping patch name " << nonOverlapPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        if (nonOverlapPatchID_ < index())
        {
            FatalErrorInFunction
                << "Boundary ordering error: " << type()
                << " patch must be defined prior to its non-overlapping patch"
                << nl
                << type() << " patch: " << name() << ", ID:" << index() << nl
                << "Non-overlap patch: " << nonOverlapPatchName_
                << ", ID:" << nonOverlapPatchID_ << nl
                << exit(FatalError);
        }

        const polyPatch& noPp = this->boundaryMesh()[nonOverlapPatchID_];

        bool ok = true;

        if (size() == noPp.size())
        {
            const scalarField magSf(magFaceAreas());
            const scalarField noMagSf(noPp.magFaceAreas());

            forAll(magSf, facei)
            {
                scalar ratio = mag(magSf[facei]/(noMagSf[facei] + ROOTVSMALL));

                if (ratio - 1 > tolerance_)
                {
                    ok = false;
                    break;
                }
            }
        }
        else
        {
            ok = false;
        }

        if (!ok)
        {
            FatalErrorInFunction
                << "Inconsistent ACMI patches " << name() << " and "
                << noPp.name() << ".  Patches should have identical topology"
                << exit(FatalError);
        }
    }

    return nonOverlapPatchID_;
}


void Foam::cyclicACMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    cyclicAMIPolyPatch::initOrder(pBufs, pp);
}


bool Foam::cyclicACMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    return cyclicAMIPolyPatch::order(pBufs, pp, faceMap, rotation);
}


void Foam::cyclicACMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    os.writeEntry("nonOverlapPatch", nonOverlapPatchName_);
}


// ************************************************************************* //
