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
    (c) 2021-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "discreteMixingPlanePolyPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "AMIInterpolation/patches/cyclicPeriodicAMI/cyclicPeriodicAMIPolyPatch/cyclicPeriodicAMIPolyPatch.H"

// For debugging
#include "surfaceFormats/obj/OBJstream.H"
#include "meshes/primitiveMesh/PatchTools/PatchTools.H"
#include "db/Time/Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(discreteMixingPlanePolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, discreteMixingPlanePolyPatch, word);
    addToRunTimeSelectionTable
    (
        polyPatch,
        discreteMixingPlanePolyPatch,
        dictionary
    );
}

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::discreteMixingPlanePolyPatch::computeSector() const
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


void Foam::discreteMixingPlanePolyPatch::syncTransforms() const
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


Foam::transformer Foam::discreteMixingPlanePolyPatch::getTransform
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


void Foam::discreteMixingPlanePolyPatch::writeOBJ
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


void Foam::discreteMixingPlanePolyPatch::resetAMI
(
    const AMIPatchToPatchInterpolation::interpolationMethod& AMIMethod
) const
{

    computeSector();

    if (owner())
    {
        srcAMITransforms_.resize(0);
        tgtAMITransforms_.resize(0);

        const transformer oneTr = transformer::I;

        // Synchronise the transforms
        syncTransforms();

        // Create copies of both patches' points, transformed to the owner
        pointField thisPoints0(localPoints());
        pointField nbrPoints0(nbrPatch().localPoints());
        transform().transformPosition(nbrPoints0, nbrPoints0);

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

        // Number of geometry replications
        label iter(0);

        if (ownStr.valid())
        {
            writeOBJ(thisPatch0, ownStr());
        }
        if (neiStr.valid())
        {
            writeOBJ(nbrPatch0, neiStr());
        }

        if (debug)
        {
            Info<< "onSectors(): " << onSectors() << endl;
            Info<< "nbSectors(): " << nbSectors() << endl;
        }

        srcAMITransforms_.append(oneTr);
        tgtAMITransforms_.append(oneTr);

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

        // Weight sum averages
        scalar srcSum(gAverage(AMIPtr_->srcWeightsSum()));
        scalar tgtSum(gAverage(AMIPtr_->tgtWeightsSum()));

        if (debug)
        {
            InfoInFunction
                << "patch:" << name()
                << " srcSum:" << srcSum
                << " tgtSum:" << tgtSum
                << endl;
        }

        nbrPoints = nbrPoints0;
        for (label iNb = 2; iNb <= nbSectors(); iNb++)
        {
            const coupledPolyPatch& nbPeriodicPatch
            (
                refCast<const coupledPolyPatch>
                (
                    boundaryMesh()[neiPeriodicPatchID()]
                )
            );

            nbPeriodicPatch.transform().transformPosition
            (
                nbrPoints,
                nbrPoints
            );
            nbrPatch.clearGeom();
            AMIPtr_->append(thisPatch, nbrPatch, iter+1);

            transformer t = getTransform(nbPeriodicPatch);
            transformer& tBef = srcAMITransforms_[iter];
            tgtAMITransforms_.append(oneTr);
            srcAMITransforms_.append(t&tBef);

            ++iter;

            if (ownStr.valid())
            {
                writeOBJ(thisPatch, ownStr());
                Info<< "writeOBJon: " << thisPatch.size() << endl;
            }
            if (neiStr.valid())
            {
                writeOBJ(nbrPatch, neiStr());
            }

            srcSum = gAverage(AMIPtr_->srcWeightsSum());
            tgtSum = gAverage(AMIPtr_->tgtWeightsSum());

            if (debug)
            {
                InfoInFunction
                    << "patch:" << name()
                    << " iteration:" << iter
                    << " srcSum:" << srcSum
                    << " tgtSum:" << tgtSum
                    << endl;
            }
        }

        for (label iOn = 2; iOn <= onSectors(); iOn++)
        {
            const coupledPolyPatch& onPeriodicPatch
            (
                refCast<const coupledPolyPatch>
                (
                    boundaryMesh()[periodicPatchID()]
                )
            );

            {
                onPeriodicPatch.transform().transformPosition
                (
                    thisPoints,
                    thisPoints
                );
                thisPatch.clearGeom();

            }

            transformer sF =
                getTransform(onPeriodicPatch) & tgtAMITransforms_[iter];

            for (label iNb = 1; iNb <= nbSectors(); iNb++)
            {
                //- reset to original position
                if (iNb == 1)
                {
                    nbrPoints = nbrPoints0;
                    nbrPatch.clearGeom();

                    AMIPtr_->append(thisPatch, nbrPatch, iter+1);
                    srcAMITransforms_.append(oneTr);

                }
                else
                {
                    const coupledPolyPatch& nbPeriodicPatch
                    (
                        refCast<const coupledPolyPatch>
                        (
                            boundaryMesh()[neiPeriodicPatchID()]
                        )
                    );

                    nbPeriodicPatch.transform().transformPosition
                    (
                        nbrPoints,
                        nbrPoints
                    );
                    nbrPatch.clearGeom();
                    AMIPtr_->append(thisPatch, nbrPatch, iter+1);

                    transformer tF(oneTr);
                    transformer t = getTransform(nbPeriodicPatch);
                    transformer& tBef = srcAMITransforms_[iter];
                    tF = t&tBef;
                    srcAMITransforms_.append(tF);
                }
                tgtAMITransforms_.append(sF);

                ++iter;
                {
                    if (ownStr.valid())
                    {
                        writeOBJ(thisPatch, ownStr());
                    }
                    if (neiStr.valid())
                    {
                        writeOBJ(nbrPatch, neiStr());
                    }

                    srcSum = gAverage(AMIPtr_->srcWeightsSum());
                    tgtSum = gAverage(AMIPtr_->tgtWeightsSum());

                    if (debug)
                    {
                        InfoInFunction
                            << "patch:" << name()
                            << " iteration:" << iter
                            << " srcSum:" << srcSum
                            << " tgtSum:" << tgtSum
                            << endl;
                    }
                }
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
    if (debug)
    {
        Info<< "srcAMITransforms_ " <<  srcAMITransforms_ <<endl;
        Info<< "tgtAMITransforms_ " <<  tgtAMITransforms_ <<endl;
    }
    AMIPtr_->setTransformations(srcAMITransforms_, tgtAMITransforms_);
}


void Foam::discreteMixingPlanePolyPatch::calcGeometry(PstreamBuffers& pBufs)
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

Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
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
    nSectors_(0),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
    AMILowWeightCorrection_ = -1.0;
}


Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
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
    nSectors_(dict.lookupOrDefault<label>("nSectors", 0)),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
    AMILowWeightCorrection_ = dict.lookupOrDefault("lowWeightCorrection", 0.01);
}


Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
(
    const discreteMixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    nSectors_(pp.nSectors_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
}


Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
(
    const discreteMixingPlanePolyPatch& pp,
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
    nSectors_(pp.nSectors_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
}


Foam::discreteMixingPlanePolyPatch::discreteMixingPlanePolyPatch
(
    const discreteMixingPlanePolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    periodicPatchName_(pp.periodicPatchName_),
    periodicPatchID_(-1),
    nSectors_(pp.nSectors_),
    srcAMITransforms_(0),
    tgtAMITransforms_(0)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::discreteMixingPlanePolyPatch::~discreteMixingPlanePolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::discreteMixingPlanePolyPatch::nSectors() const
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


Foam::label Foam::discreteMixingPlanePolyPatch::periodicPatchID() const
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


Foam::label Foam::discreteMixingPlanePolyPatch::neiPeriodicPatchID() const
{
    const discreteMixingPlanePolyPatch& cAMIpp =
        refCast<const discreteMixingPlanePolyPatch>
        (
            nbrPatch()
        );

    return cAMIpp.periodicPatchID();
}


Foam::label Foam::discreteMixingPlanePolyPatch::onSectors() const
{
    return this->nSectors();
}


Foam::label Foam::discreteMixingPlanePolyPatch::nbSectors() const
{
    const discreteMixingPlanePolyPatch& cAMIpp =
        refCast<const discreteMixingPlanePolyPatch>
        (
            nbrPatch()
        );
    const label nbSectors = cAMIpp.nSectors();

    return nbSectors;
}


void Foam::discreteMixingPlanePolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    if (periodicPatchName_ != word::null)
    {
        os.writeEntry("periodicPatch", periodicPatchName_);
    }

    if (nSectors_ != 0)
    {
        os.writeEntry("nSectors", nSectors_);
    }
}


// ************************************************************************* //
