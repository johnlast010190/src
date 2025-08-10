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

\*---------------------------------------------------------------------------*/

#include "AMIInterpolation/AMIInterpolation/AMIMethod/directAMI/directAMI.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::appendToDirectSeeds
(
    labelList& mapFlag,
    labelList& srcTgtSeed,
    DynamicList<label>& srcSeeds,
    DynamicList<label>& nonOverlapFaces,
    label& srcFacei,
    label& tgtFacei
) const
{
    const labelList& srcNbr = this->srcPatch_.faceFaces()[srcFacei];
    const labelList& tgtNbr = this->tgtPatch_.faceFaces()[tgtFacei];

    const pointField& srcPoints = this->srcPatch_.points();
    const pointField& tgtPoints = this->tgtPatch_.points();

    const vectorField& srcCf = this->srcPatch_.faceCentres();

    forAll(srcNbr, i)
    {
        label srcI = srcNbr[i];

        if ((mapFlag[srcI] == 0) && (srcTgtSeed[srcI] == -1))
        {
            // first attempt: match by comparing face centres
            const face& srcF = this->srcPatch_[srcI];
            const point& srcC = srcCf[srcI];

            scalar tol = GREAT;
            forAll(srcF, fpI)
            {
                const point& p = srcPoints[srcF[fpI]];
                scalar d2 = magSqr(p - srcC);
                if (d2 < tol)
                {
                    tol = d2;
                }
            }
            tol = max(SMALL, 0.0001*sqrt(tol));

            bool found = false;
            forAll(tgtNbr, j)
            {
                label tgtI = tgtNbr[j];
                const face& tgtF = this->tgtPatch_[tgtI];
                const point tgtC = tgtF.centre(tgtPoints);

                if (mag(srcC - tgtC) < tol)
                {
                    // new match - append to lists
                    found = true;

                    srcTgtSeed[srcI] = tgtI;
                    srcSeeds.append(srcI);

                    break;
                }
            }

            // second attempt: match by shooting a ray into the tgt face
            if (!found)
            {
                const vector srcN = srcF.areaNormal(srcPoints);

                forAll(tgtNbr, j)
                {
                    label tgtI = tgtNbr[j];
                    const face& tgtF = this->tgtPatch_[tgtI];
                    pointHit ray = tgtF.ray(srcCf[srcI], srcN, tgtPoints);

                    if (ray.hit())
                    {
                        // new match - append to lists
                        found = true;

                        srcTgtSeed[srcI] = tgtI;
                        srcSeeds.append(srcI);

                        break;
                    }
                }
            }

            // no match available for source face srcI
            if (!found)
            {
                mapFlag[srcI] = -1;
                nonOverlapFaces.append(srcI);

                if (debug)
                {
                    Pout<< "source face not found: id=" << srcI
                        << " centre=" << srcCf[srcI]
                        << " normal=" << srcF.areaNormal(srcPoints)
                        << " points=" << srcF.points(srcPoints)
                        << endl;

                    Pout<< "target neighbours:" << nl;
                    forAll(tgtNbr, j)
                    {
                        label tgtI = tgtNbr[j];
                        const face& tgtF = this->tgtPatch_[tgtI];

                        Pout<< "face id: " << tgtI
                            << " centre=" << tgtF.centre(tgtPoints)
                            << " normal=" << tgtF.areaNormal(tgtPoints)
                            << " points=" << tgtF.points(tgtPoints)
                            << endl;
                    }
                }
            }
        }
    }

    if (srcSeeds.size())
    {
        srcFacei = srcSeeds.remove();
        tgtFacei = srcTgtSeed[srcFacei];
    }
    else
    {
        srcFacei = -1;
        tgtFacei = -1;
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::restartAdvancingFront
(
    labelList& mapFlag,
    DynamicList<label>& nonOverlapFaces,
    label& srcFacei,
    label& tgtFacei
) const
{
    forAll(mapFlag, facei)
    {
        if (mapFlag[facei] == 0)
        {
            tgtFacei = this->findTargetFace(facei);

            if (tgtFacei < 0)
            {
                mapFlag[facei] = -1;
                nonOverlapFaces.append(facei);
            }
            else
            {
                srcFacei = facei;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::directAMI<SourcePatch, TargetPatch>::directAMI
(
    const SourcePatch& srcPatch,
    const TargetPatch& tgtPatch,
    const scalarField& srcMagSf,
    const scalarField& tgtMagSf,
    const faceAreaIntersect::triangulationMode& triMode,
    const scalar& cosMatchAngle_,
    const bool reverseTarget,
    const bool requireMatch
)
:
    AMIMethod<SourcePatch, TargetPatch>
    (
        srcPatch,
        tgtPatch,
        srcMagSf,
        tgtMagSf,
        triMode,
        cosMatchAngle_,
        reverseTarget,
        requireMatch
    )
{}


// * * * * * * * * * * * * * * * * Destructor * * * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
Foam::directAMI<SourcePatch, TargetPatch>::~directAMI()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::calculate
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
    List<DynamicList<label>> tgtAddr(this->tgtPatch_.size());


    // construct weights and addressing
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // list of faces currently visited for srcFacei to avoid multiple hits
    DynamicList<label> srcSeeds(10);

    // list to keep track of tgt faces used to seed src faces
    labelList srcTgtSeed(srcAddr.size(), -1);
    srcTgtSeed[srcFacei] = tgtFacei;

    // list to keep track of whether src face can be mapped
    // 1 = mapped, 0 = untested, -1 = cannot map
    labelList mapFlag(srcAddr.size(), 0);

    label nTested = 0;
    DynamicList<label> nonOverlapFaces;
    do
    {
        srcAddr[srcFacei].append(tgtFacei);
        tgtAddr[tgtFacei].append(srcFacei);

        mapFlag[srcFacei] = 1;

        nTested++;

        // Do advancing front starting from srcFacei, tgtFacei
        appendToDirectSeeds
        (
            mapFlag,
            srcTgtSeed,
            srcSeeds,
            nonOverlapFaces,
            srcFacei,
            tgtFacei
        );

        if (srcFacei < 0 && nTested < this->srcPatch_.size())
        {
            restartAdvancingFront(mapFlag, nonOverlapFaces, srcFacei, tgtFacei);
        }

    } while (srcFacei >= 0);

    if (nonOverlapFaces.size() != 0)
    {
        Pout<< "    AMI: " << nonOverlapFaces.size()
            << " non-overlap faces identified"
            << endl;

        this->srcNonOverlap_.transfer(nonOverlapFaces);
    }

    // transfer data to persistent storage
    forAll(srcAddr, i)
    {
        scalar magSf = this->srcMagSf_[i];
        srcAddress[i].transfer(srcAddr[i]);
        srcWeights[i] = scalarList(1, magSf);
    }
    forAll(tgtAddr, i)
    {
        scalar magSf = this->tgtMagSf_[i];
        tgtAddress[i].transfer(tgtAddr[i]);
        tgtWeights[i] = scalarList(1, magSf);
    }
}


template<class SourcePatch, class TargetPatch>
void Foam::directAMI<SourcePatch, TargetPatch>::normaliseWeights
(
    const bool verbose,
    AMIInterpolation<SourcePatch, TargetPatch>& inter
)
{
    inter.normaliseWeights(this->conformal(), verbose);
}


// ************************************************************************* //
