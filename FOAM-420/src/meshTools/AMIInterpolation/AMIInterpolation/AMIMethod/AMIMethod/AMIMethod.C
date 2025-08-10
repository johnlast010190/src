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
    (c) 2013-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIMethod/AMIMethod/AMIMethod.H"
#include "meshTools/meshTools.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistribute.H"
#include "global/unitConversion/unitConversion.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::checkPatches() const
{
    if (debug && (!srcPatch_.size() || !tgtPatch_.size()))
    {
        Pout<< "AMI: Patches not on processor: Source faces = "
            << srcPatch_.size() << ", target faces = " << tgtPatch_.size()
            << endl;
    }


    if (conformal())
    {
        const scalar maxBoundsError = 0.05;

        // check bounds of source and target
        boundBox bbSrc(srcPatch_.points(), srcPatch_.meshPoints(), true);
        boundBox bbTgt(tgtPatch_.points(), tgtPatch_.meshPoints(), true);

        boundBox bbTgtInf(bbTgt);
        bbTgtInf.inflate(maxBoundsError);

        if (!bbTgtInf.contains(bbSrc))
        {
            WarningInFunction
                << "Source and target patch bounding boxes are not similar"
                << nl
                << "    source box span     : " << bbSrc.span() << nl
                << "    target box span     : " << bbTgt.span() << nl
                << "    source box          : " << bbSrc << nl
                << "    target box          : " << bbTgt << nl
                << "    inflated target box : " << bbTgtInf << endl;
        }
    }
}


template<class SourcePatch, class TargetPatch>
bool Foam::AMIMethod<SourcePatch, TargetPatch>::initialise
(
    labelListList& srcAddress,
    scalarListList& srcWeights,
    scalarList& srcSymWeights,
    labelListList& tgtAddress,
    scalarListList& tgtWeights,
    scalarList& tgtSymWeights,
    label& srcFacei,
    label& tgtFacei
)
{
    checkPatches();

    // set initial sizes for weights and addressing - must be done even if
    // returns false below
    srcAddress.setSize(srcPatch_.size());
    srcWeights.setSize(srcPatch_.size());
    tgtAddress.setSize(tgtPatch_.size());
    tgtWeights.setSize(tgtPatch_.size());

    srcSymWeights.setSize(srcPatch_.size(), 0);
    tgtSymWeights.setSize(tgtPatch_.size(), 0);

    // check that patch sizes are valid
    if (!srcPatch_.size())
    {
        return false;
    }
    else if (!tgtPatch_.size())
    {
        if (debug)
        {
            WarningInFunction
                << srcPatch_.size()
                << " source faces but no target faces" << endl;
        }

        return false;
    }

    // reset the octree
    resetTree();

    // find initial face match using brute force/octree search
    if ((srcFacei == -1) || (tgtFacei == -1))
    {
        srcFacei = 0;
        tgtFacei = 0;
        bool foundFace = false;
        forAll(srcPatch_, facei)
        {
            tgtFacei = findTargetFace(facei);
            if (tgtFacei >= 0)
            {
                srcFacei = facei;
                foundFace = true;
                break;
            }
        }

#if defined(WIN64) || defined(WIN32)
        if (debug)
        {
        #endif
            if (!foundFace)
            {
                if (requireMatch_)
                {
                    FatalErrorInFunction
                        << "Unable to find initial target face"
                        << abort(FatalError);
                }

                return false;
            }
#if defined(WIN64) || defined(WIN32)
        }
        #endif
    }

    if (debug)
    {
        Pout<< "AMI: initial target face = " << tgtFacei << endl;
    }

    return true;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::writeIntersectionOBJ
(
    const scalar area,
    const face& f1,
    const face& f2,
    const pointField& f1Points,
    const pointField& f2Points
) const
{
    static label count = 1;

    const pointField f1pts = f1.points(f1Points);
    const pointField f2pts = f2.points(f2Points);

    Pout<< "Face intersection area (" << count <<  "):" << nl
        << "    f1 face = " << f1 << nl
        << "    f1 pts  = " << f1pts << nl
        << "    f2 face = " << f2 << nl
        << "    f2 pts  = " << f2pts << nl
        << "    area    = " << area
        << endl;

    OFstream os("areas" + name(count) + ".obj");

    forAll(f1pts, i)
    {
        meshTools::writeOBJ(os, f1pts[i]);
    }
    os<< "l";
    forAll(f1pts, i)
    {
        os<< " " << i + 1;
    }
    os<< " 1" << endl;


    forAll(f2pts, i)
    {
        meshTools::writeOBJ(os, f2pts[i]);
    }
    os<< "l";
    forAll(f2pts, i)
    {
        os<< " " << f1pts.size() + i + 1;
    }
    os<< " " << f1pts.size() + 1 << endl;

    count++;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::writeAddr
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const labelListList& srcAddress,
    word fname
) const
{
    OFstream os(fname + Foam::name(Pstream::myProcNo()) + ".obj");

    label ptI = 1;

    forAll(srcAddress, i)
    {
        const labelList& addr = srcAddress[i];
        const point& srcPt = srcPatch.faceCentres()[i];
        forAll(addr, j)
        {
            label tgtPtI = addr[j];
            const point& tgtPt = tgtPatch.faceCentres()[tgtPtI];

            meshTools::writeOBJ(os, srcPt);
            meshTools::writeOBJ(os, tgtPt);

            os  << "l " << ptI << " " << ptI + 1 << endl;

            ptI += 2;
        }
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::resetTree()
{
    // Clear the old octree
    treePtr_.clear();

    treeBoundBox bb(tgtPatch_.points(), tgtPatch_.meshPoints());
    bb.inflate(0.01);

    if (!treePtr_.valid())
    {
        treePtr_.reset
        (
            new indexedOctree<treeType>
            (
                treeType
                (
                    false,
                    tgtPatch_,
                    indexedOctree<treeType>::perturbTol()
                ),
                bb,                         // overall search domain
                8,                          // maxLevel
                10,                         // leaf size
                3.0                         // duplicity
            )
        );
    }
}


template<class SourcePatch, class TargetPatch>
Foam::label Foam::AMIMethod<SourcePatch, TargetPatch>::findTargetFace
(
    const label srcFacei
) const
{
    label targetFacei = -1;

    const pointField& srcPts = srcPatch_.points();
    const face& srcFace = srcPatch_[srcFacei];
    const point srcPt = srcFace.centre(srcPts);
    const scalar srcFaceArea = srcMagSf_[srcFacei];

    pointIndexHit sample = treePtr_->findNearest(srcPt, 10.0*srcFaceArea);

    if (sample.hit())
    {
        targetFacei = sample.index();

        if (debug)
        {
            Pout<< "Source point = " << srcPt << ", Sample point = "
                << sample.hitPoint() << ", Sample index = " << sample.index()
                << endl;
        }
    }

    return targetFacei;
}


template<class SourcePatch, class TargetPatch>
void Foam::AMIMethod<SourcePatch, TargetPatch>::appendNbrFaces
(
    const label srcFaceI,
    const label tgtFaceI,
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const DynamicList<label>& visitedFaces,
    DynamicList<label>& faceIDs
) const
{
    const labelList& nbrFaces = tgtPatch.faceFaces()[tgtFaceI];

    const vector& n1 = srcPatch.faceNormals()[srcFaceI];

    // filter out faces already visited from face neighbours
    forAll(nbrFaces, i)
    {
        label nbrFacei = nbrFaces[i];
        bool valid = true;
        forAll(visitedFaces, j)
        {
            if (nbrFacei == visitedFaces[j])
            {
                valid = false;
                break;
            }
        }

        if (valid)
        {
            forAll(faceIDs, j)
            {
                if (nbrFacei == faceIDs[j])
                {
                    valid = false;
                    break;
                }
            }
        }

        // prevent addition of face if it is not on the same plane-ish
        if (valid)
        {
            const vector& n2 = tgtPatch.faceNormals()[nbrFacei];

            scalar cosI = n1 & n2;
            if (!reverseTarget_)
            {
                cosI *= -1;
            }

            if (cosI > Foam::cos(degToRad(89.0)))
            {
                faceIDs.append(nbrFacei);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIMethod<SourcePatch, TargetPatch>::AMIMethod
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
    srcPatch_(srcPatch),
    tgtPatch_(tgtPatch),
    reverseTarget_(reverseTarget),
    requireMatch_(requireMatch),
    srcMagSf_(srcMagSf),
    tgtMagSf_(tgtMagSf),
    srcNonOverlap_(),
    triMode_(triMode)
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::AMIMethod<SourcePatch, TargetPatch>::~AMIMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
bool Foam::AMIMethod<SourcePatch, TargetPatch>::conformal() const
{
    return true;
}


// ************************************************************************* //
