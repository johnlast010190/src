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
    (c) 2016, 2020 OpenCFD Ltd.
    (c) 2017, 2020 Esi

\*---------------------------------------------------------------------------*/


template<class SourcePatch, class TargetPatch>
inline bool Foam::AMIInterpolation<SourcePatch, TargetPatch>::upToDate() const
{
    return upToDate_;
}


template<class SourcePatch, class TargetPatch>
inline bool& Foam::AMIInterpolation<SourcePatch, TargetPatch>::upToDate()
{
    return upToDate_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::labelList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::overlapTargetFaces() const
{
    return (*overlapTgtAddrPtr_);
}

template<class SourcePatch, class TargetPatch>
inline const Foam::labelList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::nonOverlapTargetFaces() const
{
    return (*nonOverlapTgtAddrPtr_);
}

template<class SourcePatch, class TargetPatch>
inline const Foam::labelList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::overlapSourceFaces() const
{
    return (*overlapSrcAddrPtr_);
}

template<class SourcePatch, class TargetPatch>
inline const Foam::labelList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::nonOverlapSourceFaces() const
{
    return (*nonOverlapSrcAddrPtr_);
}

////////////////////////////////////



template<class SourcePatch, class TargetPatch>
inline Foam::label
Foam::AMIInterpolation<SourcePatch, TargetPatch>::singlePatchProc() const
{
    return singlePatchProc_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalar
Foam::AMIInterpolation<SourcePatch, TargetPatch>::lowWeightCorrection() const
{
    return lowWeightCorrection_;
}


template<class SourcePatch, class TargetPatch>
inline bool
Foam::AMIInterpolation<SourcePatch, TargetPatch>::
applyLowWeightCorrection() const
{
    return lowWeightCorrection_ > 0;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcMagSf() const
{
    return srcMagSf_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcMagSf()
{
    return srcMagSf_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::vectorField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcNf() const
{
    return srcNf_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::labelListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcAddress() const
{
    return srcAddress_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::labelListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcAddress()
{
    return srcAddress_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcWeights() const
{
    return srcWeights_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcWeights()
{
    return srcWeights_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcSymWeights() const
{
    return srcSymWeights_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcSymWeights()
{
    return srcSymWeights_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcWeightsSum() const
{
    return srcWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcWeightsSum()
{
    return srcWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::pointListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcCentroids() const
{
    return srcCentroids_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::pointListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcCentroids()
{
    return srcCentroids_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::mapDistribute&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::srcMap() const
{
    return srcMapPtr_();
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtMagSf() const
{
    return tgtMagSf_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtMagSf()
{
    return tgtMagSf_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::vectorField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtNf() const
{
    return tgtNf_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::labelListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtAddress() const
{
    return tgtAddress_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::labelListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtAddress()
{
    return tgtAddress_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtWeights() const
{
    return tgtWeights_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarListList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtWeights()
{
    return tgtWeights_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtSymWeights() const
{
    return tgtSymWeights_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarList&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtSymWeights()
{
    return tgtSymWeights_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtWeightsSum() const
{
    return tgtWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline Foam::scalarField&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtWeightsSum()
{
    return tgtWeightsSum_;
}


template<class SourcePatch, class TargetPatch>
inline const Foam::mapDistribute&
Foam::AMIInterpolation<SourcePatch, TargetPatch>::tgtMap() const
{
    return tgtMapPtr_();
}


// ************************************************************************* //
