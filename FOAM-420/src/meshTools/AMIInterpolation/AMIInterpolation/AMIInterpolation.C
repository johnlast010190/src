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
    (c) 2015-2020 OpenCFD Ltd.
    (c) 2017-2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIInterpolation.H"
#include "AMIInterpolation/AMIInterpolation/AMIMethod/AMIMethod/AMIMethod.H"
#include "meshTools/meshTools.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistribute.H"
#include "primitives/ops/flipOp.H"
#include "meshes/meshTools/simpleVTKWriter.H"
#if defined(WIN64) || defined(WIN32)
#include "primitives/bools/Switch/Switch.H"
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
const Foam::scalar Foam::AMIInterpolation<SourcePatch, TargetPatch>::relTolDef_
    = debug::floatOptimisationSwitch("AMIRelTol", 0.01);


template<class SourcePatch, class TargetPatch>
const Foam::scalar Foam::AMIInterpolation<SourcePatch, TargetPatch>::
tolDef_ = 0;


template<class SourcePatch, class TargetPatch>
const Foam::scalar Foam::AMIInterpolation<SourcePatch, TargetPatch>::
cosMatchAngleDef_ =
    cos
    (
        degToRad
        (
            debug::floatOptimisationSwitch
            (
                "allowedOverlappingAngle",
                45
            )
        )
    );


template<class SourcePatch, class TargetPatch>
const Foam::scalar Foam::AMIInterpolation<SourcePatch, TargetPatch>::
maxAMIWeightScaleDef_ =
    debug::floatOptimisationSwitch
    (
        "maxAMIWeightScale",
        1.25
    );


template<class SourcePatch, class TargetPatch>
const bool Foam::AMIInterpolation<SourcePatch, TargetPatch>::outputDef_ = true;


template<class SourcePatch, class TargetPatch>
Foam::word
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolationMethodToWord
(
    const interpolationMethod& im
)
{
    word method = "unknown-interpolationMethod";

    switch (im)
    {
        case imDirect:
        {
            method = "directAMI";
            break;
        }
        case imMapNearest:
        {
            method = "mapNearestAMI";
            break;
        }
        case imFaceAreaWeight:
        {
            method = "faceAreaWeightAMI";
            break;
        }
        case imPartialFaceAreaWeight:
        {
            method = "partialFaceAreaWeightAMI";
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unhandled interpolationMethod enumeration " << method
                << abort(FatalError);
        }
    }

    return method;
}


template<class SourcePatch, class TargetPatch>
typename Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolationMethod
Foam::AMIInterpolation<SourcePatch, TargetPatch>::wordTointerpolationMethod
(
    const word& im
)
{
    interpolationMethod method = imDirect;

    wordList methods
    (
        IStringStream
        (
            "("
                "directAMI "
                "mapNearestAMI "
                "faceAreaWeightAMI "
                "partialFaceAreaWeightAMI"
            ")"
        )()
    );

    if (im == "directAMI")
    {
        method = imDirect;
    }
    else if (im == "mapNearestAMI")
    {
        method = imMapNearest;
    }
    else if (im == "faceAreaWeightAMI")
    {
        method = imFaceAreaWeight;
    }
    else if (im == "partialFaceAreaWeightAMI")
    {
        method = imPartialFaceAreaWeight;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid interpolationMethod " << im
            << ".  Valid methods are:" << methods
            << exit(FatalError);
    }

    return method;
}


template<class SourcePatch, class TargetPatch>
template<class pType>
Foam::tmp<Foam::scalarField>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::patchMagSf
(
    const pType& patch,
    const faceAreaIntersect::triangulationMode triMode
)
{
    tmp<scalarField> tResult(new scalarField(patch.size(), Zero));
    scalarField& result = tResult.ref();

    const pointField& patchPoints = patch.localPoints();

    triFaceList patchFaceTris;

    forAll(result, patchFacei)
    {
        faceAreaIntersect::triangulate
        (
            patch.localFaces()[patchFacei],
            patchPoints,
            triMode,
            patchFaceTris
        );

        forAll(patchFaceTris, i)
        {
            result[patchFacei] += patchFaceTris[i].tri(patchPoints).mag();
        }
    }

    return tResult;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::projectPointsToSurface
(
    const searchableSurface& surf,
    pointField& pts
) const
{
    if (debug)
    {
        Info<< "AMI: projecting points to surface" << endl;
    }

    List<pointIndexHit> nearInfo;

    surf.findNearest(pts, scalarField(pts.size(), GREAT), nearInfo);

    label nMiss = 0;
    forAll(nearInfo, i)
    {
        const pointIndexHit& pi = nearInfo[i];

        if (pi.hit())
        {
            pts[i] = pi.hitPoint();
        }
        else
        {
            pts[i] = pts[i];
            nMiss++;
        }
    }

    if (nMiss > 0)
    {
        FatalErrorInFunction
            << "Error projecting points to surface: "
            << nMiss << " faces could not be determined"
            << abort(FatalError);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::visualiseWeights
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const bool applyNormalisation,
    const word& ext
) const
{
    simpleVTKWriter sVTK
    (
        srcPatch.localFaces(),
        srcPatch.localPoints()
    );

    simpleVTKWriter tVTK
    (
        tgtPatch.localFaces(),
        tgtPatch.localPoints()
    );

    tmp<scalarField> swArea = patchMagSf(srcPatch, triMode_);
    tmp<scalarField> twArea = patchMagSf(tgtPatch, triMode_);

    scalarField swSum(srcPatch.localFaces().size(), 0);
    scalarField swSym(srcPatch.localFaces().size(), 0);

    forAll(srcWeights_, fI)
    {
        swSum[fI] = sum(srcWeights_[fI]);
        if (applyNormalisation)
        {
            swSum[fI] /= swArea()[fI];
        }
        swSym[fI] = srcSymWeights_[fI]/swArea()[fI];
    }

    scalarField twSum(tgtPatch.localFaces().size(), 0);
    scalarField twSym(tgtPatch.localFaces().size(), 0);

    forAll(tgtWeights_, fI)
    {
        twSum[fI] = sum(tgtWeights_[fI]);
        if (applyNormalisation)
        {
            twSum[fI] /= twArea()[fI];
        }
        twSym[fI] = tgtSymWeights_[fI]/twArea()[fI];
    }

    tVTK.addFaceData("twSum", twSum);
    tVTK.addFaceData("twSym", twSym);

    sVTK.addFaceData("swSum", swSum);
    sVTK.addFaceData("swSym", swSym);

    sVTK.write("s_"+ext+".vtk");
    tVTK.write("t_"+ext+"_.vtk");
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::normaliseWeights
(
    const scalarField& patchAreas,
    const word& patchName,
    const labelListList& addr,
    scalarListList& wght,
    scalarList& symWght,
    scalarField& wghtSum,
    const bool conformal,
    const bool output,
    const scalar lowWeightTol,
    const scalar maxScale
)
{
    // Normalise the weights
    wghtSum.setSize(wght.size(), 0.0);
    label nLowWeight = 0;

    forAll(wght, facei)
    {
        scalarList& w = wght[facei];

        if (w.size())
        {
            scalar denom = patchAreas[facei];

            scalar s = sum(w);
            scalar t = s/patchAreas[facei];

            if (conformal)
            {
                if (s>0)
                {
                    scalar scale = 1/t;
                    scale = min(scale, maxScale);
                    denom /= scale;
                    wghtSum[facei] = s/denom;
                }
            }
            else
            {
                wghtSum[facei] = t;
            }

            forAll(w, i)
            {
                w[i] /= denom;
            }

            if (t < lowWeightTol)
            {
                nLowWeight++;
            }
        }
    }

    if (output)
    {
        const label nFace = returnReduce(wght.size(), sumOp<label>());

        if (nFace)
        {
            Info<< indent
                << "AMI: Patch " << patchName
                << " sum(weights) min/max/average = "
                << gMin(wghtSum) << ", "
                << gMax(wghtSum) << ", "
                << gAverage(wghtSum) << endl;

            const label nLow = returnReduce(nLowWeight, sumOp<label>());

            if (nLow)
            {
                Info<< indent
                    << "AMI: Patch " << patchName
                    << " identified " << nLow
                    << " faces with weights less than " << lowWeightTol
                    << endl;
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::cropWeights
(
    const scalarField& patchAreas,
    const word& patchName,
    scalarListList& wght,
    scalarField& wghtSum
)
{
    // Normalise the weights
    wghtSum.setSize(wght.size(), 0.0);

    forAll(wghtSum, facei)
    {
        if (wghtSum[facei] > 1)
        {
            scalar scale = 1/wghtSum[facei];
            forAll(wght[facei], fI)
            {
                wght[facei][fI] *= scale;
            }
            wghtSum[facei] = 1;
        }
    }

    if (debug)
    {
        const label nFace = returnReduce(wght.size(), sumOp<label>());
        if (nFace)
        {
            Info<< indent
                << "AMI: Patch  " << patchName
                << " cropped sum(weights) min/max/average = "
                << gMin(wghtSum) << ", "
                << gMax(wghtSum) << ", "
                << gAverage(wghtSum) << endl;
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::setTransformations
(
    const List<transformer>& srcTr,
    const List<transformer>& tgtTr
)
{
    srcAMITransforms_ = srcTr;
    tgtAMITransforms_ = tgtTr;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::agglomerate
(
    const autoPtr<mapDistribute>& targetMapPtr,
    const scalarField& fineSrcMagSf,
    const labelListList& fineSrcAddress,
    const scalarListList& fineSrcWeights,
    const scalarList& fineSrcSymWeights,

    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing,
    const scalar maxScale,

    scalarField& srcMagSf,
    labelListList& srcAddress,
    scalarListList& srcWeights,
    scalarList& srcSymWeights,
    scalarField& srcWeightsSum,
    autoPtr<mapDistribute>& tgtMap
)
{
    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label targetCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    // Agglomerate face areas
    {
        srcMagSf.setSize(sourceRestrictAddressing.size(), 0.0);
        forAll(sourceRestrictAddressing, facei)
        {
            label coarseFacei = sourceRestrictAddressing[facei];
            srcMagSf[coarseFacei] += fineSrcMagSf[facei];
        }
    }


    // Agglomerate weights and indices
    if (targetMapPtr.valid())
    {
        const mapDistribute& map = targetMapPtr();

        // Get all restriction addressing.
        labelList allRestrict(targetRestrictAddressing);
        map.distribute(allRestrict);

        // So now we have agglomeration of the target side in
        // allRestrict:
        //  0..size-1 : local agglomeration (= targetRestrictAddressing
        //              (but potentially permutated))
        //  size..    : agglomeration data from other processors


        // The trickiness in this algorithm is finding out the compaction
        // of the remote data (i.e. allocation of the coarse 'slots'). We could
        // either send across the slot compaction maps or just make sure
        // that we allocate the slots in exactly the same order on both sending
        // and receiving side (e.g. if the submap is set up to send 4 items,
        // the constructMap is also set up to receive 4 items.


        // Short note about the various types of indices:
        // - face indices : indices into the geometry.
        // - coarse face indices : how the faces get agglomerated
        // - transferred data : how mapDistribute sends/receives data,
        // - slots : indices into data after distribution (e.g. stencil,
        //           srcAddress/tgtAddress). Note: for fully local addressing
        //           the slots are equal to face indices.
        // A mapDistribute has:
        // - a subMap : these are face indices
        // - a constructMap : these are from 'transferred-date' to slots

        labelListList tgtSubMap(Pstream::nProcs());

        // Local subMap is just identity
        {
            tgtSubMap[Pstream::myProcNo()] = identity(targetCoarseSize);
        }

        forAll(map.subMap(), proci)
        {
            if (proci != Pstream::myProcNo())
            {
                // Combine entries that point to the same coarse element.
                // The important bit is to loop over the data (and hand out
                // compact indices ) in 'transferred data' order. This
                // guarantees that we're doing exactly the
                // same on sending and receiving side - e.g. the fourth element
                // in the subMap is the fourth element received in the
                // constructMap

                const labelList& elems = map.subMap()[proci];
                const labelList& elemsMap =
                    map.constructMap()[Pstream::myProcNo()];
                labelList& newSubMap = tgtSubMap[proci];
                newSubMap.setSize(elems.size());

                labelList oldToNew(targetCoarseSize, -1);
                label newi = 0;

                forAll(elems, i)
                {
                    label fineElem = elemsMap[elems[i]];
                    label coarseElem = allRestrict[fineElem];
                    if (oldToNew[coarseElem] == -1)
                    {
                        oldToNew[coarseElem] = newi;
                        newSubMap[newi] = coarseElem;
                        newi++;
                    }
                }
                newSubMap.setSize(newi);
            }
        }

        // Reconstruct constructMap by combining entries. Note that order
        // of handing out indices should be the same as loop above to compact
        // the sending map

        labelListList tgtConstructMap(Pstream::nProcs());

        // Local constructMap is just identity
        {
            tgtConstructMap[Pstream::myProcNo()] =
                identity(targetCoarseSize);
        }

        labelList tgtCompactMap(map.constructSize());

        {
            // Note that in special cases (e.g. 'appending' two AMIs) the
            // local size after distributing can be longer than the number
            // of faces. I.e. it duplicates elements.
            // Since we don't know this size instead we loop over all
            // reachable elements (using the local constructMap)

            const labelList& elemsMap = map.constructMap()[Pstream::myProcNo()];
            forAll(elemsMap, i)
            {
                label fineElem = elemsMap[i];
                label coarseElem = allRestrict[fineElem];
                tgtCompactMap[fineElem] = coarseElem;
            }
        }

        label compacti = targetCoarseSize;

        // Compact data from other processors
        forAll(map.constructMap(), proci)
        {
            if (proci != Pstream::myProcNo())
            {
                // Combine entries that point to the same coarse element. All
                // elements now are remote data so we cannot use any local
                // data here - use allRestrict instead.
                const labelList& elems = map.constructMap()[proci];

                labelList& newConstructMap = tgtConstructMap[proci];
                newConstructMap.setSize(elems.size());

                if (elems.size())
                {
                    // Get the maximum target coarse size for this set of
                    // received data.
                    label remoteTargetCoarseSize = labelMin;
                    forAll(elems, i)
                    {
                        remoteTargetCoarseSize = max
                        (
                            remoteTargetCoarseSize,
                            allRestrict[elems[i]]
                        );
                    }
                    remoteTargetCoarseSize += 1;

                    // Combine locally data coming from proci
                    labelList oldToNew(remoteTargetCoarseSize, -1);
                    label newi = 0;

                    forAll(elems, i)
                    {
                        label fineElem = elems[i];
                        // fineElem now points to section from proci
                        label coarseElem = allRestrict[fineElem];
                        if (oldToNew[coarseElem] == -1)
                        {
                            oldToNew[coarseElem] = newi;
                            tgtCompactMap[fineElem] = compacti;
                            newConstructMap[newi] = compacti++;
                            newi++;
                        }
                        else
                        {
                            // Get compact index
                            label compacti = oldToNew[coarseElem];
                            tgtCompactMap[fineElem] = newConstructMap[compacti];
                        }
                    }
                    newConstructMap.setSize(newi);
                }
            }
        }

        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);
        srcSymWeights.setSize(sourceCoarseSize, 0.0);

        forAll(fineSrcAddress, facei)
        {
            // All the elements contributing to facei. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[facei];
            const scalarList& weights = fineSrcWeights[facei];
            const scalar fineArea = fineSrcMagSf[facei];

            label coarseFacei = sourceRestrictAddressing[facei];

            labelList& newElems = srcAddress[coarseFacei];
            scalarList& newWeights = srcWeights[coarseFacei];

            forAll(elems, i)
            {
                label elemi = elems[i];
                label coarseElemi = tgtCompactMap[elemi];

                label index = findIndex(newElems, coarseElemi);
                if (index == -1)
                {
                    newElems.append(coarseElemi);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
            srcSymWeights[coarseFacei] += fineArea*fineSrcSymWeights[facei];
        }

        tgtMap.reset
        (
            new mapDistribute
            (
                compacti,
                tgtSubMap.xfer(),
                tgtConstructMap.xfer()
            )
        );
    }
    else
    {
        srcAddress.setSize(sourceCoarseSize);
        srcWeights.setSize(sourceCoarseSize);
        srcSymWeights.setSize(sourceCoarseSize, 0.0);

        forAll(fineSrcAddress, facei)
        {
            // All the elements contributing to facei. Are slots in
            // mapDistribute'd data.
            const labelList& elems = fineSrcAddress[facei];
            const scalarList& weights = fineSrcWeights[facei];
            const scalar fineArea = fineSrcMagSf[facei];

            label coarseFacei = sourceRestrictAddressing[facei];

            labelList& newElems = srcAddress[coarseFacei];
            scalarList& newWeights = srcWeights[coarseFacei];

            forAll(elems, i)
            {
                label elemi = elems[i];
                label coarseElemi = targetRestrictAddressing[elemi];

                label index = findIndex(newElems, coarseElemi);
                if (index == -1)
                {
                    newElems.append(coarseElemi);
                    newWeights.append(fineArea*weights[i]);
                }
                else
                {
                    newWeights[index] += fineArea*weights[i];
                }
            }
            srcSymWeights[coarseFacei] += fineArea*fineSrcSymWeights[facei];
        }
    }

    // Weights normalisation
    normaliseWeights
    (
        srcMagSf,
        "source",
        srcAddress,
        srcWeights,
        srcSymWeights,
        srcWeightsSum,
        true,
        false,
        -1,
        maxScale
    );
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::constructFromSurface
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr
)
{
    if (upToDate_)
    {
        return;
    }

    if (surfPtr.valid())
    {
        // Create new patches for source and target
        pointField srcPoints = srcPatch.points();
        SourcePatch srcPatch0
        (
            SubList<face>
            (
                srcPatch,
                srcPatch.size(),
                0
            ),
            srcPoints
        );

        if (debug)
        {
            OFstream os("amiSrcPoints.obj");
            forAll(srcPoints, i)
            {
                meshTools::writeOBJ(os, srcPoints[i]);
            }
        }

        pointField tgtPoints = tgtPatch.points();
        TargetPatch tgtPatch0
        (
            SubList<face>
            (
                tgtPatch,
                tgtPatch.size(),
                0
            ),
            tgtPoints
        );

        if (debug)
        {
            OFstream os("amiTgtPoints.obj");
            forAll(tgtPoints, i)
            {
                meshTools::writeOBJ(os, tgtPoints[i]);
            }
        }


        // Map source and target patches onto projection surface
        projectPointsToSurface(surfPtr(), srcPoints);
        projectPointsToSurface(surfPtr(), tgtPoints);


        // Calculate AMI interpolation
        update(srcPatch0, tgtPatch0);
    }
    else
    {
        update(srcPatch, tgtPatch);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::
constructOverAndNonOverlappingAddressing()
{
    //- Initialization fields
    //  Target
    DynamicList<label> dNonOverlapTgtAddr(tgtAddress_.size());
    DynamicList<label> dOverlapTgtAddr(tgtAddress_.size());

    //  Source
    DynamicList<label> dNonOverlapSrcAddr(srcAddress_.size());
    DynamicList<label> dOverlapSrcAddr(srcAddress_.size());


    //- Separate overlap and non overlap addressing
    //  Target
    forAll(tgtAddress_, fI)
    {
        if (!tgtAddress_[fI].size())
        {
            dNonOverlapTgtAddr.append(fI);
        }
        else
        {
            dOverlapTgtAddr.append(fI);
        }
    }
    dNonOverlapTgtAddr.shrink();
    nonOverlapTgtAddrPtr_ = new labelList(dNonOverlapTgtAddr);

    dOverlapTgtAddr.shrink();
    overlapTgtAddrPtr_ = new labelList(dOverlapTgtAddr);

    //  Source
    forAll(srcAddress_, fI)
    {
        if (!srcAddress_[fI].size())
        {
            dNonOverlapSrcAddr.append(fI);
        }
        else
        {
            dOverlapSrcAddr.append(fI);
        }
    }
    dNonOverlapSrcAddr.shrink();
    nonOverlapSrcAddrPtr_ = new labelList(dNonOverlapSrcAddr);

    dOverlapSrcAddr.shrink();
    overlapSrcAddrPtr_ = new labelList(dOverlapSrcAddr);
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::
constructOverAndNonOverlappingAddressing
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
)
{
    //- Initialization fields
    //  Target
    DynamicList<label> dNonOverlapTgtAddr(tgtAddress_.size());
    DynamicList<label> dOverlapTgtAddr(tgtAddress_.size());

    //  Source
    DynamicList<label> dNonOverlapSrcAddr(srcAddress_.size());
    DynamicList<label> dOverlapSrcAddr(srcAddress_.size());


    if (!srcPatch.size() || !tgtPatch.size())
    {
        //Pout<< "AMI: Patches not on processor: Source faces = "
        //    << srcPatch.size() << ", target faces = " << tgtPatch.size()
        //    << endl;
    }

    if (!srcPatch.size())
    {
        forAll(tgtPatch, facei)
        {
            dNonOverlapTgtAddr.append(facei);
        }

        dNonOverlapTgtAddr.shrink();
        nonOverlapTgtAddrPtr_ = new labelList(dNonOverlapTgtAddr);

        if (debug)
        {
            WarningInFunction << "Target patch with size: " << tgtPatch.size()
                              << ", but without source patch size 0."
                              << endl;
        }
    }
    else if (!tgtPatch.size())
    {
        forAll(srcPatch, facei)
        {
            dNonOverlapSrcAddr.append(facei);
        }
        dNonOverlapSrcAddr.shrink();
        nonOverlapSrcAddrPtr_ = new labelList(dNonOverlapSrcAddr);

        if (debug)
        {
            WarningInFunction << "Source patch with size: " << srcPatch.size()
                              << ", but without target patch size 0."
                              << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Size mismatch. "
            << "one of the two must have zero size."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const interpolationMethod& method,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const scalar& relTol,
    const scalar& tol,
    const scalar& cosMatchAngle,
    const scalar& maxAMIWeightScale,
    const bool output
)
:
    methodName_(interpolationMethodToWord(method)),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    cosMatchAngle_(cosMatchAngle),
    maxScale_(maxAMIWeightScale),
    relTol_(relTol),
    tol_(tol),
    output_(output),

    srcAddress_(),
    srcWeights_(),
    srcSymWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcTransforms_(),

    tgtAddress_(),
    tgtWeights_(),
    tgtSymWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtTransforms_(),

    overlapSrcAddrPtr_(nullptr),
    nonOverlapSrcAddrPtr_(nullptr),
    overlapTgtAddrPtr_(nullptr),
    nonOverlapTgtAddrPtr_(nullptr),

    triMode_(triMode),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    srcAMITransforms_(),
    tgtAMITransforms_(),
    upToDate_(false)
{
    label srcSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtSize = returnReduce(tgtPatch.size(), sumOp<label>());

    if (output)
    {
        Info<< "AMI: Creating addressing and weights between "
            << srcSize << " source faces and " << tgtSize << " target faces"
            << endl;
    }

    update(srcPatch, tgtPatch);

    if (!srcPatch.size() || !tgtPatch.size())
    {
        constructOverAndNonOverlappingAddressing(srcPatch, tgtPatch);
    }
    else
    {
        constructOverAndNonOverlappingAddressing();
    }
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const word& methodName,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const scalar& relTol,
    const scalar& tol,
    const scalar& cosMatchAngle,
    const scalar& maxAMIWeightScale,
    const bool output
)
:
    methodName_(methodName),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    cosMatchAngle_(cosMatchAngle),
    maxScale_(maxAMIWeightScale),
    relTol_(relTol),
    tol_(tol),
    output_(output),

    srcAddress_(),
    srcWeights_(),
    srcSymWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcTransforms_(),

    tgtAddress_(),
    tgtWeights_(),
    tgtSymWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtTransforms_(),

    overlapSrcAddrPtr_(nullptr),
    nonOverlapSrcAddrPtr_(nullptr),
    overlapTgtAddrPtr_(nullptr),
    nonOverlapTgtAddrPtr_(nullptr),

    triMode_(triMode),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    srcAMITransforms_(),
    tgtAMITransforms_(),
    upToDate_(false)
{
    label srcSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtSize = returnReduce(tgtPatch.size(), sumOp<label>());

    if (output)
    {
        Info<< "AMI: Creating addressing and weights between "
            << srcSize << " source faces and " << tgtSize << " target faces"
            << endl;
    }

    update(srcPatch, tgtPatch);

    if (!srcPatch.size() || !tgtPatch.size())
    {
        constructOverAndNonOverlappingAddressing(srcPatch, tgtPatch);
    }
    else
    {
        constructOverAndNonOverlappingAddressing();
    }
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const interpolationMethod& method,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const scalar& relTol,
    const scalar& tol,
    const scalar& cosMatchAngle,
    const scalar& maxAMIWeightScale,
    const bool output
)
:
    methodName_(interpolationMethodToWord(method)),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    cosMatchAngle_(cosMatchAngle),
    maxScale_(maxAMIWeightScale),
    relTol_(relTol),
    tol_(tol),
    output_(output),

    srcAddress_(),
    srcWeights_(),
    srcSymWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcTransforms_(),

    tgtAddress_(),
    tgtWeights_(),
    tgtSymWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtTransforms_(),

    overlapSrcAddrPtr_(nullptr),
    nonOverlapSrcAddrPtr_(nullptr),
    overlapTgtAddrPtr_(nullptr),
    nonOverlapTgtAddrPtr_(nullptr),

    triMode_(triMode),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    srcAMITransforms_(),
    tgtAMITransforms_(),
    upToDate_(false)
{
    constructFromSurface(srcPatch, tgtPatch, surfPtr);

    constructOverAndNonOverlappingAddressing();
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr,
    const faceAreaIntersect::triangulationMode& triMode,
    const bool requireMatch,
    const word& methodName,
    const scalar lowWeightCorrection,
    const bool reverseTarget,
    const scalar& relTol,
    const scalar& tol,
    const scalar& cosMatchAngle,
    const scalar& maxAMIWeightScale,
    const bool output
)
:
    methodName_(methodName),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    singlePatchProc_(-999),
    lowWeightCorrection_(lowWeightCorrection),
    cosMatchAngle_(cosMatchAngle),
    maxScale_(maxAMIWeightScale),
    relTol_(relTol),
    tol_(tol),
    output_(output),

    srcAddress_(),
    srcWeights_(),
    srcSymWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcTransforms_(),

    tgtAddress_(),
    tgtWeights_(),
    tgtSymWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtTransforms_(),

    overlapSrcAddrPtr_(nullptr),
    nonOverlapSrcAddrPtr_(nullptr),
    overlapTgtAddrPtr_(nullptr),
    nonOverlapTgtAddrPtr_(nullptr),

    triMode_(triMode),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    srcAMITransforms_(),
    tgtAMITransforms_(),
    upToDate_(false)
{
    constructFromSurface(srcPatch, tgtPatch, surfPtr);

    constructOverAndNonOverlappingAddressing();
}


template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::AMIInterpolation
(
    const AMIInterpolation<SourcePatch, TargetPatch>& fineAMI,
    const labelList& sourceRestrictAddressing,
    const labelList& targetRestrictAddressing
)
:
    methodName_(fineAMI.methodName_),
    reverseTarget_(fineAMI.reverseTarget_),
    requireMatch_(fineAMI.requireMatch_),
    singlePatchProc_(fineAMI.singlePatchProc_),
    lowWeightCorrection_(-1.0),
    cosMatchAngle_(fineAMI.cosMatchAngle_),
    maxScale_(fineAMI.maxScale_),
    relTol_(fineAMI.relTol_),
    tol_(fineAMI.tol_),
    output_(false),

    srcAddress_(),
    srcWeights_(),
    srcSymWeights_(),
    srcWeightsSum_(),
    srcCentroids_(),
    srcTransforms_(),

    tgtAddress_(),
    tgtWeights_(),
    tgtSymWeights_(),
    tgtWeightsSum_(),
    tgtCentroids_(),
    tgtTransforms_(),

    overlapSrcAddrPtr_(nullptr),
    nonOverlapSrcAddrPtr_(nullptr),
    overlapTgtAddrPtr_(nullptr),
    nonOverlapTgtAddrPtr_(nullptr),

    triMode_(fineAMI.triMode_),
    srcMapPtr_(nullptr),
    tgtMapPtr_(nullptr),
    srcAMITransforms_(),
    tgtAMITransforms_(),
    upToDate_(false)
{
    label sourceCoarseSize =
    (
        sourceRestrictAddressing.size()
      ? max(sourceRestrictAddressing)+1
      : 0
    );

    label neighbourCoarseSize =
    (
        targetRestrictAddressing.size()
      ? max(targetRestrictAddressing)+1
      : 0
    );

    if (debug & 2)
    {
        Pout<< "AMI: Creating addressing and weights as agglomeration of AMI :"
            << " source:" << fineAMI.srcAddress().size()
            << " target:" << fineAMI.tgtAddress().size()
            << " coarse source size:" << sourceCoarseSize
            << " neighbour source size:" << neighbourCoarseSize
            << endl;
    }

    if
    (
        fineAMI.srcAddress().size() != sourceRestrictAddressing.size()
     || fineAMI.tgtAddress().size() != targetRestrictAddressing.size()
    )
    {
        FatalErrorInFunction
            << "Size mismatch." << nl
            << "Source patch size:" << fineAMI.srcAddress().size() << nl
            << "Source agglomeration size:"
            << sourceRestrictAddressing.size() << nl
            << "Target patch size:" << fineAMI.tgtAddress().size() << nl
            << "Target agglomeration size:"
            << targetRestrictAddressing.size()
            << exit(FatalError);
    }

    // Agglomerate addresses and weights

    agglomerate
    (
        fineAMI.tgtMapPtr_,
        fineAMI.srcMagSf(),
        fineAMI.srcAddress(),
        fineAMI.srcWeights(),
        fineAMI.srcSymWeights(),

        sourceRestrictAddressing,
        targetRestrictAddressing,
        maxScale_,

        srcMagSf_,
        srcAddress_,
        srcWeights_,
        srcSymWeights_,
        srcWeightsSum_,
        tgtMapPtr_
    );

    agglomerate
    (
        fineAMI.srcMapPtr_,
        fineAMI.tgtMagSf(),
        fineAMI.tgtAddress(),
        fineAMI.tgtWeights(),
        fineAMI.tgtSymWeights(),

        targetRestrictAddressing,
        sourceRestrictAddressing,
        maxScale_,

        tgtMagSf_,
        tgtAddress_,
        tgtWeights_,
        tgtSymWeights_,
        tgtWeightsSum_,
        srcMapPtr_
    );

    constructOverAndNonOverlappingAddressing();
}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::~AMIInterpolation()
{
   deleteDemandDrivenData(overlapSrcAddrPtr_);
   deleteDemandDrivenData(nonOverlapSrcAddrPtr_);
   deleteDemandDrivenData(overlapTgtAddrPtr_);
   deleteDemandDrivenData(nonOverlapTgtAddrPtr_);
   srcMapPtr_.clear();
   tgtMapPtr_.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::update
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch
)
{
    label srcTotalSize = returnReduce(srcPatch.size(), sumOp<label>());
    label tgtTotalSize = returnReduce(tgtPatch.size(), sumOp<label>());

    if (srcTotalSize == 0)
    {
        if (debug)
        {
            Info<< "AMI: no source faces present - no addressing constructed"
                << endl;
        }

        return;
    }

    if (output_)
    {
        Info<< indent
            << "AMI: Creating addressing and weights between "
            << srcTotalSize << " source faces and "
            << tgtTotalSize << " target faces"
            << endl;
    }

    // Calculate face areas

    srcMagSf_ = patchMagSf(srcPatch, triMode_);
    tgtMagSf_ = patchMagSf(tgtPatch, triMode_);
    srcNf_.setSize(srcPatch.size());
    forAll(srcMagSf_, facei)
    {
        srcNf_[facei] = srcPatch[facei].areaNormal(srcPatch.points());
    }
    tgtNf_.setSize(tgtPatch.size());
    forAll(tgtMagSf_, facei)
    {
        tgtNf_[facei] = tgtPatch[facei].areaNormal(tgtPatch.points());
    }

    // Calculate if patches present on multiple processors
    singlePatchProc_ = calcDistribution(srcPatch, tgtPatch);

    if (singlePatchProc_ == -1)
    {
        // Convert local addressing to global addressing
        globalIndex globalSrcFaces(srcPatch.size());
        globalIndex globalTgtFaces(tgtPatch.size());

        // Create processor map of overlapping faces. This map gets
        // (possibly remote) faces from the tgtPatch such that they (together)
        // cover all of the srcPatch
        autoPtr<mapDistribute> mapPtr = calcProcMap(srcPatch, tgtPatch);
        const mapDistribute& map = mapPtr();

        // Create new target patch that fully encompasses source patch

        // Faces and points
        faceList newTgtFaces;
        pointField newTgtPoints;

        // Original faces from tgtPatch (in globalIndexing since might be
        // remote)
        labelList tgtFaceIDs;
        distributeAndMergePatches
        (
            map,
            tgtPatch,
            globalTgtFaces,
            newTgtFaces,
            newTgtPoints,
            tgtFaceIDs
        );

        TargetPatch
            newTgtPatch
            (
                SubList<face>
                (
                    newTgtFaces,
                    newTgtFaces.size()
                ),
                newTgtPoints
            );

        scalarField newTgtMagSf(patchMagSf(newTgtPatch, triMode_));

        if (false)
        {
            {
                OFstream os
                (
                    "newTargetFaces" + Foam::name(Pstream::myProcNo()) + ".obj"
                );
                meshTools::writeOBJ
                (
                    os,
                    newTgtFaces,
                    newTgtPatch.points()
                );
            }
            {
                OFstream os
                (
                    "sourcePatch_" + Foam::name(Pstream::myProcNo()) + ".obj"
                );
                meshTools::writeOBJ
                (
                    os,
                    srcPatch.localFaces(),
                    srcPatch.localPoints()
                );
            }
        }

        // Calculate AMI interpolation
        autoPtr<AMIMethod<SourcePatch, TargetPatch>> AMIPtr
        (
            AMIMethod<SourcePatch, TargetPatch>::New
            (
                methodName_,
                srcPatch,
                newTgtPatch,
                srcMagSf_,
                newTgtMagSf,
                triMode_,
                cosMatchAngle_,
                reverseTarget_,
                requireMatch_ && (lowWeightCorrection_ < 0)
            )
        );

        AMIPtr->calculate
        (
            srcAddress_,
            srcWeights_,
            srcSymWeights_,
            srcCentroids_,
            tgtAddress_,
            tgtWeights_,
            tgtSymWeights_
        );

        // Now
        // ~~~
        //  srcAddress_ :   per srcPatch face a list of the newTgtPatch (not
        //                  tgtPatch) faces it overlaps
        //  tgtAddress_ :   per newTgtPatch (not tgtPatch) face a list of the
        //                  srcPatch faces it overlaps

        if (debug)
        {
            writeFaceConnectivity(srcPatch, newTgtPatch, srcAddress_);
        }


        // Rework newTgtPatch indices into globalIndices of tgtPatch
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        forAll(srcAddress_, i)
        {
            labelList& addressing = srcAddress_[i];
            forAll(addressing, addri)
            {
                addressing[addri] = tgtFaceIDs[addressing[addri]];
            }
        }

        forAll(tgtAddress_, i)
        {
            labelList& addressing = tgtAddress_[i];
            forAll(addressing, addri)
            {
                addressing[addri] = globalSrcFaces.toGlobal(addressing[addri]);
            }
        }

        // Send data back to originating procs. Note that contributions
        // from different processors get added (ListAppendEqOp)

        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            false,                      // has flip
            map.subMap(),
            false,                      // has flip
            tgtAddress_,
            ListAppendEqOp<label>(),
            flipOp(),                   // flip operation
            labelList()
        );

        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            false,
            map.subMap(),
            false,
            tgtWeights_,
            ListAppendEqOp<scalar>(),
            flipOp(),
            scalarList()
        );


        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            tgtPatch.size(),
            map.constructMap(),
            false,
            map.subMap(),
            false,
            tgtSymWeights_,
            //eqOp<scalar>(),
            plusEqOp<scalar>(),
            flipOp(),
            scalar(0)
        );


        //- Only for debug
        //  visualiseWeights(srcPatch, tgtPatch, true ,"before");

        // weights normalisation
        AMIPtr->normaliseWeights(output_, *this);

        //- Only for debug
        //  visualiseWeights(srcPatch, tgtPatch, false, "after");

        // cache maps and reset addresses
        List<Map<label>> cMap;
        srcMapPtr_.reset(new mapDistribute(globalSrcFaces, tgtAddress_, cMap));
        tgtMapPtr_.reset(new mapDistribute(globalTgtFaces, srcAddress_, cMap));
    }
    else
    {

        // Calculate AMI interpolation
        autoPtr<AMIMethod<SourcePatch, TargetPatch>> AMIPtr
        (
            AMIMethod<SourcePatch, TargetPatch>::New
            (
                methodName_,
                srcPatch,
                tgtPatch,
                srcMagSf_,
                tgtMagSf_,
                triMode_,
                cosMatchAngle_,
                reverseTarget_,
                requireMatch_ && (lowWeightCorrection_ < 0)
            )
        );

        AMIPtr->calculate
        (
            srcAddress_,
            srcWeights_,
            srcSymWeights_,
            srcCentroids_,
            tgtAddress_,
            tgtWeights_,
            tgtSymWeights_
        );

        if (debug)
        {
            writeFaceConnectivity(srcPatch, tgtPatch, srcAddress_);
        }

        // weights normalisation
        AMIPtr->normaliseWeights(output_, *this);

    }

    if (debug)
    {
        Info<< "AMIInterpolation : Constructed addressing and weights" << nl
            << "    triMode        :"
            << faceAreaIntersect::triangulationModeNames_[triMode_] << nl
            << "    singlePatchProc:" << singlePatchProc_ << nl
            << "    srcMagSf       :" << gSum(srcMagSf_) << nl
            << "    tgtMagSf       :" << gSum(tgtMagSf_) << nl
            << endl;
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::reset
(
    autoPtr<mapDistribute>&& srcToTgtMap,
    autoPtr<mapDistribute>&& tgtToSrcMap,
    labelListList&& srcAddress,
    scalarListList&& srcWeights,
    labelListList&& tgtAddress,
    scalarListList&& tgtWeights
)
{
    DebugInFunction<< endl;

    srcAddress_.transfer(srcAddress);
    srcWeights_.transfer(srcWeights);
    tgtAddress_.transfer(tgtAddress);
    tgtWeights_.transfer(tgtWeights);

    // Reset the sums of the weights
    srcWeightsSum_.setSize(srcWeights_.size());
    forAll(srcWeights_, facei)
    {
        srcWeightsSum_[facei] = sum(srcWeights_[facei]);
    }

    tgtWeightsSum_.setSize(tgtWeights_.size());
    forAll(tgtWeights_, facei)
    {
        tgtWeightsSum_[facei] = sum(tgtWeights_[facei]);
    }

    srcMapPtr_ = srcToTgtMap;
    tgtMapPtr_ = tgtToSrcMap;

// SOS
    upToDate_ = true;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::append
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const label tranfI
)
{
    // Create a new interpolation
    autoPtr<AMIInterpolation<SourcePatch, TargetPatch>> newPtr
    (
        new AMIInterpolation<SourcePatch, TargetPatch>
        (
            srcPatch,
            tgtPatch,
            triMode_,
            requireMatch_,
            methodName_,
            lowWeightCorrection_,
            reverseTarget_,
            relTol_,
            tol_,
            cosMatchAngle_,
            maxScale_,
            false
        )
    );

    // If parallel then combine the mapDistribution and re-index
    if (singlePatchProc_ == -1)
    {
        labelListList& srcSubMap = srcMapPtr_->subMap();
        labelListList& srcConstructMap = srcMapPtr_->constructMap();

        labelListList& tgtSubMap = tgtMapPtr_->subMap();
        labelListList& tgtConstructMap = tgtMapPtr_->constructMap();

        labelListList& newSrcSubMap = newPtr->srcMapPtr_->subMap();
        labelListList& newSrcConstructMap = newPtr->srcMapPtr_->constructMap();

        labelListList& newTgtSubMap = newPtr->tgtMapPtr_->subMap();
        labelListList& newTgtConstructMap = newPtr->tgtMapPtr_->constructMap();

        // Re-calculate the source indices
        {
            labelList mapMap(0), newMapMap(0);
            forAll(srcSubMap, proci)
            {
                mapMap.append
                (
                    identity(srcConstructMap[proci].size())
                  + mapMap.size() + newMapMap.size()
                );
                newMapMap.append
                (
                    identity(newSrcConstructMap[proci].size())
                  + mapMap.size() + newMapMap.size()
                );
            }

            forAll(srcSubMap, proci)
            {
                forAll(srcConstructMap[proci], srci)
                {
                    srcConstructMap[proci][srci] =
                        mapMap[srcConstructMap[proci][srci]];
                }
            }

            forAll(srcSubMap, proci)
            {
                forAll(newSrcConstructMap[proci], srci)
                {
                    newSrcConstructMap[proci][srci] =
                        newMapMap[newSrcConstructMap[proci][srci]];
                }
            }

            forAll(tgtAddress_, tgti)
            {
                forAll(tgtAddress_[tgti], tgtj)
                {
                    tgtAddress_[tgti][tgtj] =
                        mapMap[tgtAddress_[tgti][tgtj]];
                }
            }

            forAll(newPtr->tgtAddress_, tgti)
            {
                forAll(newPtr->tgtAddress_[tgti], tgtj)
                {
                    newPtr->tgtAddress_[tgti][tgtj] =
                        newMapMap[newPtr->tgtAddress_[tgti][tgtj]];
                }
            }
        }

        // Re-calculate the target indices
        {
            labelList mapMap(0), newMapMap(0);
            forAll(srcSubMap, proci)
            {
                mapMap.append
                (
                    identity(tgtConstructMap[proci].size())
                  + mapMap.size() + newMapMap.size()
                );
                newMapMap.append
                (
                    identity(newTgtConstructMap[proci].size())
                  + mapMap.size() + newMapMap.size()
                );
            }

            forAll(srcSubMap, proci)
            {
                forAll(tgtConstructMap[proci], tgti)
                {
                    tgtConstructMap[proci][tgti] =
                        mapMap[tgtConstructMap[proci][tgti]];
                }
            }

            forAll(srcSubMap, proci)
            {
                forAll(newTgtConstructMap[proci], tgti)
                {
                    newTgtConstructMap[proci][tgti] =
                        newMapMap[newTgtConstructMap[proci][tgti]];
                }
            }

            forAll(srcAddress_, srci)
            {
                forAll(srcAddress_[srci], srcj)
                {
                    srcAddress_[srci][srcj] =
                        mapMap[srcAddress_[srci][srcj]];
                }
            }

            forAll(newPtr->srcAddress_, srci)
            {
                forAll(newPtr->srcAddress_[srci], srcj)
                {
                    newPtr->srcAddress_[srci][srcj] =
                        newMapMap[newPtr->srcAddress_[srci][srcj]];
                }
            }
        }

        // Sum the construction sizes
        srcMapPtr_->constructSize() += newPtr->srcMapPtr_->constructSize();
        tgtMapPtr_->constructSize() += newPtr->tgtMapPtr_->constructSize();

        // Combine the maps
        forAll(srcSubMap, proci)
        {
            srcSubMap[proci].append(newSrcSubMap[proci]);
            srcConstructMap[proci].append(newSrcConstructMap[proci]);

            tgtSubMap[proci].append(newTgtSubMap[proci]);
            tgtConstructMap[proci].append(newTgtConstructMap[proci]);
        }
    }
    if
    (
        (srcTransforms_.size() == 0)
     && (srcWeights_.size() != 0)
    )
    {
        srcTransforms_.setSize(srcWeights_.size());
        forAll(srcWeights_, sfI)
        {
            srcTransforms_[sfI] = labelList(srcWeights_[sfI].size(), 0);
        }
    }

    if
    (
        (tgtTransforms_.size() == 0)
     && (tgtWeights_.size() != 0)
    )
    {
        tgtTransforms_.setSize(tgtWeights_.size());
        forAll(tgtWeights_, sfI)
        {
            tgtTransforms_[sfI].setSize(tgtWeights_[sfI].size());
            tgtTransforms_[sfI] = 0;
        }
    }

    // Combine new and current source data
    forAll(srcMagSf_, srcFacei)
    {
        srcAddress_[srcFacei].append(newPtr->srcAddress()[srcFacei]);
        srcWeights_[srcFacei].append(newPtr->srcWeights()[srcFacei]);
        srcWeightsSum_[srcFacei] += newPtr->srcWeightsSum()[srcFacei];
        srcSymWeights_[srcFacei] += newPtr->srcSymWeights()[srcFacei];

        srcTransforms_[srcFacei].append
        (
            labelList(newPtr->srcWeights()[srcFacei].size(), tranfI)
        );
    }

    // Combine new and current target data
    forAll(tgtMagSf_, tgtFacei)
    {
        tgtAddress_[tgtFacei].append(newPtr->tgtAddress()[tgtFacei]);
        tgtWeights_[tgtFacei].append(newPtr->tgtWeights()[tgtFacei]);
        tgtWeightsSum_[tgtFacei] += newPtr->tgtWeightsSum()[tgtFacei];
        tgtSymWeights_[tgtFacei] += newPtr->tgtSymWeights()[tgtFacei];

        tgtTransforms_[tgtFacei].append
        (
            labelList(newPtr->tgtWeights()[tgtFacei].size(), tranfI)
        );
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::normaliseWeights
(
    const bool conformal,
    const bool output
)
{
    normaliseWeights
    (
        srcMagSf_,
        "source",
        srcAddress_,
        srcWeights_,
        srcSymWeights_,
        srcWeightsSum_,
        conformal,
        output,
        lowWeightCorrection_,
        maxScale_
    );

    normaliseWeights
    (
        tgtMagSf_,
        "target",
        tgtAddress_,
        tgtWeights_,
        tgtSymWeights_,
        tgtWeightsSum_,
        conformal,
        output,
        lowWeightCorrection_,
        maxScale_
    );

}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::cropWeights()
{
    cropWeights
    (
        srcMagSf_,
        "source",
        srcWeights_,
        srcWeightsSum_
    );

    cropWeights
    (
        tgtMagSf_,
        "target",
        tgtWeights_,
        tgtWeightsSum_
    );
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    if (fld.size() != srcAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to source patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != tgtAddress_.size())
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    target patch   = " << tgtAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(tgtAddress_.size());

    bool transformed = (tgtAMITransforms_.size()>1);

    if (singlePatchProc_ == -1)
    {
        const mapDistribute& map = srcMapPtr_();

        List<Type> work(fld);
        map.distribute(work);

        forAll(result, facei)
        {
            const labelList& faces = tgtAddress_[facei];
            const scalarList& weights = tgtWeights_[facei];
            scalar w2 = 1 - tgtWeightsSum_[facei];

            forAll(faces, i)
            {

                if (transformed)
                {
                    Type transformedType = work[faces[i]];

                    const labelList& transIDs = tgtTransforms_[facei];
                    const label transID = transIDs[i];
                    const transformer& tgtTrans = tgtAMITransforms_[transID];
                    transformedType = tgtTrans.transformVector(work[faces[i]]);

                    const transformer srcTrans =
                        inv(srcAMITransforms_[transID]);

                    transformedType = srcTrans.transformVector(transformedType);

                    cop(result[facei], facei, transformedType, weights[i]);
                }
                else
                {
                    cop(result[facei], facei, work[faces[i]], weights[i]);
                }
            }
            if (lowWeightCorrection_ > 0)
            {
                cop(result[facei], facei, defaultValues[facei], w2);
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            const labelList& faces = tgtAddress_[facei];
            const scalarList& weights = tgtWeights_[facei];
            scalar w2 = 1-tgtWeightsSum_[facei];

            forAll(faces, i)
            {
                if (transformed)
                {
                    Type transformedType = fld[faces[i]];

                    const labelList& transIDs = tgtTransforms_[facei];
                    const label transID = transIDs[i];
                    const transformer& tgtTrans = tgtAMITransforms_[transID];
                    transformedType = tgtTrans.transformVector(transformedType);

                    const transformer srcTrans =
                        inv(srcAMITransforms_[transID]);

                    transformedType = srcTrans.transformVector(transformedType);

                    cop(result[facei], facei, transformedType, weights[i]);
                }
                else
                {
                    cop(result[facei], facei, fld[faces[i]], weights[i]);
                }
            }
            if (lowWeightCorrection_ > 0)
            {
                cop(result[facei], facei, defaultValues[facei], w2);
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const UList<Type>& fld,
    const CombineOp& cop,
    List<Type>& result,
    const UList<Type>& defaultValues
) const
{
    if (fld.size() != tgtAddress_.size())
    {
        FatalErrorInFunction
            << "Supplied field size is not equal to target patch size" << nl
            << "    source patch   = " << srcAddress_.size() << nl
            << "    target patch   = " << tgtAddress_.size() << nl
            << "    supplied field = " << fld.size()
            << abort(FatalError);
    }

    if (lowWeightCorrection_ > 0)
    {
        if (defaultValues.size() != srcAddress_.size())
        {
            FatalErrorInFunction
                << "Employing default values when sum of weights falls below "
                << lowWeightCorrection_
                << " but supplied default field size is not equal to target "
                << "patch size" << nl
                << "    default values = " << defaultValues.size() << nl
                << "    source patch   = " << srcAddress_.size() << nl
                << abort(FatalError);
        }
    }

    result.setSize(srcAddress_.size());

    bool transformed = (srcAMITransforms_.size()>1);

    if (singlePatchProc_ == -1)
    {
        const mapDistribute& map = tgtMapPtr_();

        List<Type> work(fld);
        map.distribute(work);


        forAll(result, facei)
        {
            const labelList& faces = srcAddress_[facei];
            const scalarList& weights = srcWeights_[facei];
            scalar w2 = 1-srcWeightsSum_[facei];
            forAll(faces, i)
            {
                Type transformedType = work[faces[i]];

                if (transformed)
                {
                    const labelList& transIDs = srcTransforms_[facei];
                    const label transID = transIDs[i];
                    const transformer& srcTrans = srcAMITransforms_[transID];
                    transformedType = srcTrans.transformVector(work[faces[i]]);

                    const transformer tgtTrans =
                        inv(tgtAMITransforms_[transID]);

                    transformedType = tgtTrans.transformVector(transformedType);

                    cop(result[facei], facei, transformedType, weights[i]);
                }
                else
                {
                    cop(result[facei], facei, work[faces[i]], weights[i]);
                }
            }
            if (lowWeightCorrection_ > 0)
            {
                cop(result[facei], facei, defaultValues[facei], w2);
            }
        }
    }
    else
    {
        forAll(result, facei)
        {
            const labelList& faces = srcAddress_[facei];
            const scalarList& weights = srcWeights_[facei];
            scalar w2 = 1-srcWeightsSum_[facei];

            forAll(faces, i)
            {
                if (transformed)
                {
                    Type transformedType = fld[faces[i]];

                    const labelList& transIDs = srcTransforms_[facei];
                    const label transID = transIDs[i];
                    const transformer& srcTrans = srcAMITransforms_[transID];
                    transformedType = srcTrans.transformVector(transformedType);

                    const transformer tgtTrans =
                        inv(tgtAMITransforms_[transID]);

                    transformedType = tgtTrans.transformVector(transformedType);

                    cop(result[facei], facei, transformedType, weights[i]);
                }
                else
                {
                    cop(result[facei], facei, fld[faces[i]], weights[i]);
                }
            }
            if (lowWeightCorrection_ > 0)
            {
                cop(result[facei], facei, defaultValues[facei], w2);
            }
        }
    }
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            srcAddress_.size(),
            Zero
        )
    );

    interpolateToSource
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), cop, defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    tmp<Field<Type>> tresult
    (
        new Field<Type>
        (
            tgtAddress_.size(),
            Zero
        )
    );

    interpolateToTarget
    (
        fld,
        multiplyWeightedOp<Type, CombineOp>(cop),
        tresult.ref(),
        defaultValues
    );

    return tresult;
}


template<class SourcePatch, class TargetPatch>
template<class Type, class CombineOp>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const CombineOp& cop,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), cop, defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(fld, plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToSource
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToSource(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const Field<Type>& fld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(fld, plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::AMIInterpolation<SourcePatch, TargetPatch>::interpolateToTarget
(
    const tmp<Field<Type>>& tFld,
    const UList<Type>& defaultValues
) const
{
    return interpolateToTarget(tFld(), plusEqOp<Type>(), defaultValues);
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcPointFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const vector& n,
    const label tgtFacei,
    point& tgtPoint
)
const
{
    const pointField& srcPoints = srcPatch.points();

    // Source face addresses that intersect target face tgtFacei
    const labelList& addr = tgtAddress_[tgtFacei];

    pointHit nearest;
    nearest.setDistance(GREAT);
    label nearestFacei = -1;

    forAll(addr, i)
    {
        const label srcFacei = addr[i];
        const face& f = srcPatch[srcFacei];

        pointHit ray = f.ray(tgtPoint, n, srcPoints);

        if (ray.hit())
        {
            // tgtPoint = ray.rawPoint();
            return srcFacei;
        }
        else if (ray.distance() < nearest.distance())
        {
            nearest = ray;
            nearestFacei = srcFacei;
        }
    }

    if (nearest.hit() || nearest.eligibleMiss())
    {
        // tgtPoint = nearest.rawPoint();
        return nearestFacei;
    }

    return -1;
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtPointFace
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const vector& n,
    const label srcFacei,
    point& srcPoint
)
const
{
    const pointField& tgtPoints = tgtPatch.points();

    pointHit nearest;
    nearest.setDistance(GREAT);
    label nearestFacei = -1;

    // Target face addresses that intersect source face srcFacei
    const labelList& addr = srcAddress_[srcFacei];

    forAll(addr, i)
    {
        const label tgtFacei = addr[i];
        const face& f = tgtPatch[tgtFacei];

        pointHit ray = f.ray(srcPoint, n, tgtPoints);

        if (ray.hit() || ray.eligibleMiss())
        {
            // srcPoint = ray.rawPoint();
            return tgtFacei;
        }
        else if (ray.distance() < nearest.distance())
        {
            nearest = ray;
            nearestFacei = tgtFacei;
        }
    }

    if (nearest.hit() || nearest.eligibleMiss())
    {
        // srcPoint = nearest.rawPoint();
        return nearestFacei;
    }

    return -1;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIInterpolation<SourcePatch, TargetPatch>::writeFaceConnectivity
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const labelListList& srcAddress
)
const
{
    OFstream os("faceConnectivity" + Foam::name(Pstream::myProcNo()) + ".obj");

    label pti = 1;

    forAll(srcAddress, i)
    {
        const labelList& addr = srcAddress[i];
        const point& srcPt = srcPatch.faceCentres()[i];

        forAll(addr, j)
        {
            label tgtPti = addr[j];
            const point& tgtPt = tgtPatch.faceCentres()[tgtPti];

            meshTools::writeOBJ(os, srcPt);
            meshTools::writeOBJ(os, tgtPt);

            os  << "l " << pti << " " << pti + 1 << endl;

            pti += 2;
        }
    }
}


// ************************************************************************* //
