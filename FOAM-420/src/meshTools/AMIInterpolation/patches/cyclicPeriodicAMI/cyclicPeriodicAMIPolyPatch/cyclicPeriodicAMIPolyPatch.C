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
    (c) 2015 OpenCFD Ltd.
    (c) 2018-2020 OpenFOAM Foundation
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/patches/cyclicPeriodicAMI/cyclicPeriodicAMIPolyPatch/cyclicPeriodicAMIPolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// For debugging
#include "surfaceFormats/obj/OBJstream.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicPeriodicAMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicPeriodicAMIPolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        cyclicPeriodicAMIPolyPatch,
        dictionary
    );
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::cyclicPeriodicAMIPolyPatch::computeSector() const
{
    if (nSectors_ == 0)
    {
        if (periodicPatchID()==-1)
        {
            //- Corresponds to 360 | periodicPatchName_ == word::null
            nSectors_ = 1;
        }
        else
        {
            const cyclicTransform& cycTransfPatch
            (
                dynamic_cast<const cyclicTransform&>
                (
                    boundaryMesh()[periodicPatchID()]
                )
            );

            const coupledPolyPatch& cpp
            (
                refCast<const coupledPolyPatch>
                (
                    boundaryMesh()[periodicPatchID()]
                )
            );

            if (cycTransfPatch.transformType() != ROTATIONAL)
            {
                FatalErrorInFunction
                    << "Periodic patch of cyclicPeriodicAMI patch "
                    << this->name() << " with name " << cpp.name()
                    << "has " << cycTransfPatch.transformType()
                    << " transformation type. " << endl
                    << "Currently only rotational is supported for "
                    << "cyclicPeriodicAMI patches."
                    << exit(FatalError);
            }

            scalar sec = cpp.sectors();
            scalar nSec = round(sec);

            if (mag(nSec-sec) > matchTolerance())
            {
                FatalErrorInFunction
                    << "Sector numbers are not computed correctly for patch name "
                    << this->name()
                    << ", sector = "
                    << sec
                    << this->name()
                    << endl
                    << "Check mesh generation or if sector numbers is close to"
                    << " integer, then increase match tolerance: "
                    << matchTolerance()
                    << exit(FatalError);
            }

            nSectors_ = nSec;

        }

        Info<< "Sectors not defined. "
             << "Computed automatically for patch: " << this->name()
             << ", sectors = "
             << nSectors_
             << endl;
    }
}


void Foam::cyclicPeriodicAMIPolyPatch::syncTransforms() const
{
    if (owner())
    {
        // At this point we guarantee that the transformations have been
        // updated. There is one particular case, where if the periodic patch
        // locally has zero faces then its transformation will not be set. We
        // need to synchronise the transforms over the zero-sized patches as
        // well.
        //
        // We can't put the logic into cyclicPolyPatch as
        // processorCyclicPolyPatch uses cyclicPolyPatch functionality.
        // Synchronisation will fail because processor-type patches do not exist
        // on all processors.
        //
        // The use in cyclicPeriodicAMI is special; we use the patch even
        // though we have no faces in it. Ideally, in future, the transformation
        // logic will be abstracted out, and we won't need a periodic patch
        // here. Until then, we can just fix the behaviour of the zero-sized
        // coupled patches here

        // Get the periodic patch
        const coupledPolyPatch& periodicPatch
        (
            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[periodicPatchID()]
            )
        );

        // If there are any zero-sized periodic patches
        if (returnReduce((size() && !periodicPatch.size()), orOp<bool>()))
        {
            // Note that zero-sized patches will have zero-sized fields for the
            // separation vector, forward and reverse transforms. These need
            // replacing with the transformations from other processors.

            // Parallel in this context refers to a parallel transformation,
            // rather than a rotational transformation.

            // Note that a cyclic with zero faces is considered parallel so
            // explicitly check for that.
            bool isParallel =
            (
                periodicPatch.size()
             && !periodicPatch.transform().transforms()
            );
            reduce(isParallel, orOp<bool>());

            if (isParallel)
            {
                // Sync a list of separation vectors
                List<vector> sep(Pstream::nProcs());
                sep[Pstream::myProcNo()] = periodicPatch.transform().t();
                Pstream::allGatherList(sep);

                // If locally we have zero faces pick the first one that has a
                // separation vector
                if (!periodicPatch.size())
                {
                    forAll(sep, procI)
                    {
                        if (sep[procI].size())
                        {
                            const_cast<vector&>
                            (
                                periodicPatch.transform().t()
                            ) = sep[procI];

                            break;
                        }
                    }
                }
            }
            else
            {
                // Sync a list of forward and reverse transforms
                List<tensor> forwardT(Pstream::nProcs());
                forwardT[Pstream::myProcNo()] = periodicPatch.transform().T();
                Pstream::allGatherList(forwardT);

                // If locally we have zero faces pick the first one that has a
                // transformation vector
                if (!periodicPatch.size())
                {
                    forAll(forwardT, procI)
                    {
                        if (forwardT[procI].size())
                        {
                            const_cast<tensor&>
                            (
                                periodicPatch.transform().T()
                            ) = forwardT[procI];

                            break;
                        }
                    }
                }
            }
        }
    }
}


Foam::transformer Foam::cyclicPeriodicAMIPolyPatch::getTransform
(
    const coupledPolyPatch& transformPatch
) const
{
    // Get the transform associated with the transform patch
    transformer t;
    {
        Tuple2<bool, transformer> bt
        (
            transformPatch.size(),
            transformPatch.transform()
        );

        reduce(bt, keepIfTrueOp<transformer>());

        if (!bt.first())
        {
            FatalErrorInFunction
                << "Transform patch " << transformPatch.name() << " for "
                << typeName << " patch " << name() << " has zero faces. It may "
                << "have been converted to a processor cyclic during "
                << "decomposition. Consider adding " << transformPatch.name()
                << " and it's neighbour to the list of preserved patches."
                << exit(FatalError);
        }

        t = bt.second();
    }

    return t;
}


void Foam::cyclicPeriodicAMIPolyPatch::writeOBJ
(
    const primitivePatch& p,
    OBJstream& str
) const
{
    // Collect faces and points
    pointField allPoints;
    faceList allFaces;
    labelList pointMergeMap;
    PatchTools::gatherAndMerge
    (
        -1.0,           // do not merge points
        p,
        allPoints,
        allFaces,
        pointMergeMap
    );

    if (Pstream::master())
    {
        // Write base geometry
        str.write(allFaces, allPoints);
    }
}


void Foam::cyclicPeriodicAMIPolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{

    computeSector();

    if (owner())
    {
        // Get the periodic patch
        const coupledPolyPatch& oPeriodicPatch
        (
            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[smallerSectorPatchID()]
            )
        );

        srcAMITransforms_.resize(0);
        tgtAMITransforms_.resize(0);

        const transformer oneTr = transformer::I;

        // Synchronise the transforms
        syncTransforms();

        // Create copies of both patches' points, transformed to the owner
        pointField thisPoints0(localPoints());
        pointField nbrPoints0(nbrPatch().localPoints());
        transform().transformPosition(nbrPoints0, nbrPoints0);

        srcAMITransforms_.append(oneTr);
        tgtAMITransforms_.append(oneTr);

        /*
        // Apply the stored number of periodic transforms
        for (label i = 0; i < nTransforms_; ++ i)
        {
            oPeriodicPatch.transformPosition(thisPoints0);

            transformer t = getTransform(oPeriodicPatch);
            tgtAMITransforms_[0] &= t;
        }
        for (label i = 0; i > nTransforms_; -- i)
        {
            oPeriodicPatch.transformPosition(nbrPoints0);
            transformer t = getTransform(oPeriodicPatch);
            srcAMITransforms_[0] &= t;
        }
        */

        autoPtr<OBJstream> ownStr;
        autoPtr<OBJstream> neiStr;
        if (false)
        {
            const Time& runTime = boundaryMesh().mesh().time();

            fileName dir(runTime.rootPath()/runTime.globalCaseName());
            fileName postfix("_" + runTime.timeName()+"_expanded.obj");

            ownStr.reset(new OBJstream(dir/name() + postfix));
            neiStr.reset(new OBJstream(dir/nbrPatch().name() + postfix));

            InfoInFunction
                << "patch:" << name()
                << " writing accumulated AMI to " << ownStr().name()
                << " and " << neiStr().name() << endl;
        }

        // Create another copy
        pointField thisPoints(thisPoints0);
        pointField nbrPoints(nbrPoints0);

        // Create patches for all the points

        // Source patch at initial location
        const primitivePatch thisPatch0
        (
            SubList<face>(localFaces(), size()),
            thisPoints0
        );
        // Source patch that gets moved
        primitivePatch thisPatch
        (
            SubList<face>(localFaces(), size()),
            thisPoints
        );
        // Target patch at initial location
        const primitivePatch nbrPatch0
        (
            SubList<face>(nbrPatch().localFaces(), nbrPatch().size()),
            nbrPoints0
        );
        // Target patch that gets moved
        primitivePatch nbrPatch
        (
            SubList<face>
            (
                this->nbrPatch().localFaces(),
                this->nbrPatch().size()
            ),
            nbrPoints
        );

        {
            label srcSize = returnReduce(thisPatch0.size(), sumOp<label>());
            label tgtSize = returnReduce(nbrPatch0.size(), sumOp<label>());

            Info<< "AMI: Creating addressing and weights between "
                << srcSize << " source faces and " << tgtSize << " target faces"
                << endl;
        }

        // Construct a new AMI interpolation between the initial patch locations
        AMIPtr_.reset
        (
            new AMIPatchToPatchInterpolation
            (
                thisPatch0,
                nbrPatch0,
                surfPtr(),
                faceAreaIntersect::tmMesh,
                false,
                AMIPatchToPatchInterpolation::imPartialFaceAreaWeight,
                AMILowWeightCorrection_,
                AMIReverse_,
                0.01,
                0,
                degToRad
                (
                    debug::floatOptimisationSwitch
                    (
                        "allowedOverlappingAngle",
                        45
                    )
                ),
                1.25,
                false
            )
        );

        // Number of geometry replications
        label iter(0);
        label nTransformsOld(nTransforms_);

        if (ownStr.valid())
        {
            writeOBJ(thisPatch0, ownStr());
        }
        if (neiStr.valid())
        {
            writeOBJ(nbrPatch0, neiStr());
        }


        // Weight sum averages
        scalar srcSum(gAverage(AMIPtr_->srcWeightsSum()));
        scalar tgtSum(gAverage(AMIPtr_->tgtWeightsSum()));

        // Direction (or rather side of AMI : this or nbr patch) of
        // geometry replication
        bool direction = nTransforms_ >= 0;

        // Increase in the source weight sum for the last iteration in the
        // opposite direction. If the current increase is less than this, the
        // direction (= side of AMI to transform) is reversed.
        // We switch the side to replicate instead of reversing the transform
        // since at the coupledPolyPatch level there is no
        // 'reverseTransformPosition' functionality.
        scalar srcSumDiff = 0;

        if (debug)
        {
            InfoInFunction
                << "patch:" << name()
                << " srcSum:" << srcSum
                << " tgtSum:" << tgtSum
                << " direction:" << direction
                << endl;
        }

        // Loop, replicating the geometry
        while
        (
            (iter < maxIter_)
         && (
                (1 - srcSum > matchTolerance())
             || (1 - tgtSum > matchTolerance())
            )
        )
        {
            if (!onNbSmallerSector())
            //if (direction)
            {
                oPeriodicPatch.transform().transformPosition
                (
                    thisPoints,
                    thisPoints
                );

                if (debug)
                {
                    InfoInFunction
                        << "patch:" << name()
                        << " moving this side from:"
                        << gAverage(thisPatch.points())
                        << " to:" << gAverage(thisPoints) << endl;
                }

                thisPatch.clearGeom();

                if (debug)
                {
                    InfoInFunction
                        << "patch:" << name()
                        << " appending weights with untransformed slave side"
                        << endl;
                }

                AMIPtr_->append(thisPatch, nbrPatch0, iter+1);

                transformer t = getTransform(oPeriodicPatch);

                transformer& tBef = tgtAMITransforms_[iter];
                srcAMITransforms_.append(oneTr);
                tgtAMITransforms_.append(t&tBef);

                if (ownStr.valid())
                {
                    writeOBJ(thisPatch, ownStr());
                }
            }
            else
            {
                oPeriodicPatch.transform().transformPosition
                (
                    nbrPoints,
                    nbrPoints
                );

                if (debug)
                {
                    InfoInFunction
                        << "patch:" << name()
                        << " moving neighbour side from:"
                        << gAverage(nbrPatch.points())
                        << " to:" << gAverage(nbrPoints) << endl;
                }

                nbrPatch.clearGeom();

                AMIPtr_->append(thisPatch0, nbrPatch, iter+1);

                transformer t = getTransform(oPeriodicPatch);

                transformer& tBef = srcAMITransforms_[iter];
                tgtAMITransforms_.append(oneTr);
                srcAMITransforms_.append(t&tBef);

                if (neiStr.valid())
                {
                    writeOBJ(nbrPatch, neiStr());
                }
            }

            const scalar srcSumNew = gAverage(AMIPtr_->srcWeightsSum());
            const scalar srcSumDiffNew = srcSumNew - srcSum;

            if (srcSumDiffNew < srcSumDiff || srcSumDiffNew < SMALL)
            {
                direction = !direction;

                srcSumDiff = srcSumDiffNew;
            }

            srcSum = srcSumNew;
            tgtSum = gAverage(AMIPtr_->tgtWeightsSum());

            nTransforms_ += direction ? +1 : -1;

            ++iter;

            if (debug)
            {
                InfoInFunction
                    << "patch:" << name()
                    << " iteration:" << iter
                    << " srcSum:" << srcSum
                    << " tgtSum:" << tgtSum
                    << " direction:" << direction
                    << endl;
            }
        }


        // Close debug streams
        if (ownStr.valid())
        {
            ownStr.clear();
        }
        if (neiStr.valid())
        {
            neiStr.clear();
        }



        // Average the number of transformstions
        nTransforms_ = (nTransforms_ + nTransformsOld)/2;

        // Check that the match is complete
        if (iter == maxIter_)
        {
            // The matching algorithm has exited without getting the
            // srcSum and tgtSum above 1. This can happen because
            // - of an incorrect setup
            // - or because of non-exact meshes and truncation errors
            //   (transformation, accumulation of cutting errors)
            // so for now this situation is flagged as a SeriousError instead of
            // a FatalError since the default matchTolerance is quite strict
            // (0.001) and can get triggered far into the simulation.
            SeriousErrorInFunction
                << "Patches " << name() << " and " << this->nbrPatch().name()
                << " do not couple to within a tolerance of "
                << matchTolerance()
                << " when transformed according to the periodic patch "
                << oPeriodicPatch.name() << "." << nl
                << "The current sum of weights are for owner " << name()
                << " : " << srcSum << " and for neighbour "
                << this->nbrPatch().name() << " : " << tgtSum << nl
                << "This is only acceptable during post-processing"
                << "; not during running. Improve your mesh or increase"
                << " the 'matchTolerance' setting in the patch specification."
                << endl;
        }

        // Check that both patches are covered based on tolerance
        if
        (
            mag(srcSum) < 1-matchTolerance()
         || mag(tgtSum) < 1-matchTolerance()
        )
        {
            // This condition is currently enforced until there is more
            // experience with the matching algorithm and numerics.
            // This check means that e.g. different numbers of stator and
            // rotor partitions are not allowed.
            // Again see the section above about tolerances.
            SeriousErrorInFunction
                << "Patches " << name() << " and " << this->nbrPatch().name()
                << " do not overlap an integer number of times when transformed"
                << " according to the periodic patch "
                << oPeriodicPatch.name() << "." << nl
                << "The current matchTolerance : " << matchTolerance()
                << ", sum of owner weights : " << srcSum
                << ", sum of neighbour weights : " << tgtSum
                 << "." << nl
                << "This is only acceptable during post-processing"
                << "; not during running. Improve your mesh or increase"
                << " the 'matchTolerance' setting in the patch specification."
                << endl;
        }

        // Normalise the weights. Disable printing since weights are
        // still areas.
        AMIPtr_->normaliseWeights(true, false);

        {
            // Print some statistics
            const label nFace = returnReduce
            (
                AMIPtr_->srcWeights().size(), sumOp<label>()
            );

            if (nFace)
            {
                scalarField srcWghtSum(size(), 0);
                forAll(*this, faceI)
                {
                    srcWghtSum[faceI] = sum(AMIPtr_->srcWeights()[faceI]);
                }

                Info<< indent
                    << "AMI: Patch " << name()
                    << " sum(weights) min/max/average = "
                    << gMin(srcWghtSum) << ", "
                    << gMax(srcWghtSum) << ", "
                    << gAverage(srcWghtSum) << endl;
            }
        }
        {
            // Print some statistics
            const label nFace = returnReduce
            (
                AMIPtr_->tgtWeights().size(), sumOp<label>()
            );

            if (nFace)
            {
                scalarField tgtWghtSum(AMIPtr_->tgtWeights().size(), 0);
                forAll(tgtWghtSum, faceI)
                {
                    tgtWghtSum[faceI] = sum(AMIPtr_->tgtWeights()[faceI]);
                }
                Info<< indent
                    << "AMI: Patch " << this->nbrPatch().name()
                    << " sum(weights) min/max/average = "
                    << gMin(tgtWghtSum) << ", "
                    << gMax(tgtWghtSum) << ", "
                    << gAverage(tgtWghtSum) << endl;
            }
        }

        if (false)
        {
            AMIPtr_->visualiseWeights
            (
                *this, this->nbrPatch(),
                false ,
                this->name()+"grgEnd"+this->boundaryMesh().time().timeName()
            );
        }
    }
    AMIPtr_->setTransformations(srcAMITransforms_, tgtAMITransforms_);
}


void Foam::cyclicPeriodicAMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
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


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
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
    periodicPatchName_(word::null),
    periodicPatchID_(-1),
    nTransforms_(0),
    nSectors_(0),
    maxIter_(36),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
    AMILowWeightCorrection_ = -1.0;
}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType
)
:
    cyclicAMIPolyPatch(name, dict, index, bm, patchType),
    periodicPatchName_(dict.lookupOrDefault<word>("periodicPatch", word::null)),
    periodicPatchID_(-1),
    nTransforms_(dict.lookupOrDefault<label>("nTransforms", 0)),
    nSectors_(dict.lookupOrDefault<label>("nSectors", 0)),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 36)),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
    AMILowWeightCorrection_ = dict.lookupOrDefault("lowWeightCorrection", 0.01);
}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
}


Foam::cyclicPeriodicAMIPolyPatch::cyclicPeriodicAMIPolyPatch
(
    const cyclicPeriodicAMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    nTransforms_(pp.nTransforms_),
    nSectors_(pp.nSectors_),
    maxIter_(pp.maxIter_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cyclicPeriodicAMIPolyPatch::~cyclicPeriodicAMIPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cyclicPeriodicAMIPolyPatch::nSectors() const
{
    if (nSectors_==0)
    {
        computeSector();
    }
    if (nSectors_<1)
    {
        FatalErrorInFunction
            << "Sector for patch " << this->name()
            << " is not computed correctly"
            << endl
            << "nSectors: "
            << nSectors_
            << exit(FatalError);
    }
    return nSectors_;
}


Foam::label Foam::cyclicPeriodicAMIPolyPatch::periodicPatchID() const
{
    if
    (
        periodicPatchName_ == word::null
      ||periodicPatchName_ == "none"
    )
    {
        //- corresponds in 360 sector - no periodic patch name defined

        periodicPatchID_ = -1;

        return periodicPatchID_;
    }

    if (periodicPatchID_ == -1)
    {
        periodicPatchID_ = this->boundaryMesh().findPatchID(periodicPatchName_);

        if (periodicPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal periodicPatch name " << periodicPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        // Check that it is a coupled patch
        refCast<const coupledPolyPatch>
        (
            this->boundaryMesh()[periodicPatchID_]
        );
    }

    return periodicPatchID_;
}


Foam::label Foam::cyclicPeriodicAMIPolyPatch::neiPeriodicPatchID() const
{
    const cyclicPeriodicAMIPolyPatch& cAMIpp =
        refCast<const cyclicPeriodicAMIPolyPatch>
        (
            nbrPatch()
        );

    return cAMIpp.periodicPatchID();
}


Foam::label Foam::cyclicPeriodicAMIPolyPatch::smallerSectorPatchID() const
{
    return onNbSmallerSector() ? neiPeriodicPatchID() : periodicPatchID();
}


bool Foam::cyclicPeriodicAMIPolyPatch::onNbSmallerSector() const
{
    const cyclicPeriodicAMIPolyPatch& cAMIpp =
        refCast<const cyclicPeriodicAMIPolyPatch>
        (
            nbrPatch()
        );
    const label onSectors = this->nSectors();
    const label nbSectors = cAMIpp.nSectors();

    if ((onSectors == 1) && (nbSectors == 1))
    {
        FatalErrorInFunction
            << "typeName " << " patches: "
            << this->name() << " and " << cAMIpp.name()
            << " completing 360 degree." << endl
            << "Change the types to cyclicAMI."
            << exit(FatalError);
    }

    return (onSectors<nbSectors);
}


void Foam::cyclicPeriodicAMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    if (periodicPatchName_ != word::null)
    {
        os.writeEntry("periodicPatch", periodicPatchName_);
    }

    if (nTransforms_ != 0)
    {
        os.writeEntry("nTransforms", nTransforms_);
    }

    if (nSectors_ != 0)
    {
        os.writeEntry("nSectors", nSectors_);
    }

    if (maxIter_ != 36)
    {
        os.writeEntry("maxIter", maxIter_);
    }
}


// ************************************************************************* //
