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
    (c) 2020 OpenCFD Ltd.
    (c) 2021 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIMethod/partialFaceAreaWeightAMI/partialFaceAreaWeightAMI.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::setNextFaces
(
    label& startSeedI,
    label& srcFacei,
    label& tgtFacei,
    const boolList& mapFlag,
    labelList& seedFaces,
    const DynamicList<label>& visitedFaces,
    const bool errorOnNotFound
) const
{
    faceAreaWeightAMI<SourcePatch, TargetPatch>::setNextFaces
    (
        startSeedI,
        srcFacei,
        tgtFacei,
        mapFlag,
        seedFaces,
        visitedFaces,
        false // no error on not found
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::
partialFaceAreaWeightAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const scalar& cosMatchAngle,
    const bool reverseTarget,
    const bool requireMatch
)
:
    faceAreaWeightAMI<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        srcMagSf,
        tgtMagSf,
        triMode,
        cosMatchAngle,
        reverseTarget,
        requireMatch
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::
~partialFaceAreaWeightAMI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
bool Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::conformal() const
{
    return false;
}


template<class SourcePatch, class TargetPatch>
void Foam::partialFaceAreaWeightAMI<SourcePatch, TargetPatch>::calculate
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    scalarList& srcSymWeights,
    pointListList& srcCentroids,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    scalarList& tgtSymWeights,
    label srcFacei,
    label tgtFacei
)
{
    bool ok =
        this->initialise
        (
            srcAddress,
            srcWeights,
            srcSymWeights,
            tgtAddress,
            tgtWeights,
            tgtSymWeights,
            srcFacei,
            tgtFacei
        );

    if (!ok)
    {
        return;
    }

    // temporary storage for addressing and weights
    List<DynamicList<label>> srcAddr(this->srcPatch_.size());
    List<DynamicList<scalar>> srcWght(srcAddr.size());
    List<scalar> srcSymWght(srcAddr.size(), Zero);
    List<DynamicList<label>> tgtAddr(this->tgtPatch_.size());
    List<DynamicList<scalar>> tgtWght(tgtAddr.size());
    List<scalar> tgtSymWght(tgtAddr.size(), Zero);
    List<DynamicList<point>> srcCtr(srcAddr.size());

    faceAreaWeightAMI<SourcePatch, TargetPatch>::calcAddressing
    (
        srcAddr,
        srcWght,
        srcSymWght,
        srcCtr,
        tgtAddr,
        tgtWght,
        tgtSymWght,
        srcFacei,
        tgtFacei
    );

    // Check for badly covered faces
    if (this->restartUncoveredFaces())
    {
        this->restartUncoveredSourceFace
        (
            srcAddr,
            srcWght,
            srcSymWght,
            srcCtr,
            tgtAddr,
            tgtWght,
            tgtSymWght
        );
    }


    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i].transfer(srcWght[i]);
        srcSymWeights[i] = srcSymWght[i];
    }
    forAll(tgtAddr, i)
    {
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i].transfer(tgtWght[i]);
        tgtSymWeights[i] = tgtSymWght[i];
    }
//    this->writeAddr(this->srcPatch_, this->tgtPatch_, srcAddress ,"addre");
}


// ************************************************************************* //
