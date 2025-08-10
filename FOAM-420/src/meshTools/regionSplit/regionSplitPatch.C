/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM® : Professional Open-source CFD
|   o   O   o    |  Version : 4.2.0
|    o     o     |  Copyright © 2015 ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
License
    This file is part of FOAMcore.
    FOAMcore is based on OpenFOAM® <http://www.openfoam.org/>.

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

Portions Copyright © 2011 OpenFOAM Foundation.

\*---------------------------------------------------------------------------*/

#include "regionSplit/regionSplitPatch.H"
#include "meshes/polyMesh/syncTools/syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(regionSplitPatch, 0);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::regionSplitPatch::fillSeedMask
(
    labelList& faceRegion,
    labelList& edgeRegion,
    const label seedFaceID,
    const label markValue
) const
{
    // Do seed face
    faceRegion[seedFaceID] = markValue;

    // Collect edges on seed face
    const labelList& fEdges = pp_.faceEdges()[seedFaceID];

    label nEdges = 0;

    labelList changedEdges(fEdges.size());

    forAll(fEdges, i)
    {
        label edgeI = fEdges[i];

        if (edgeRegion[edgeI] == -1)
        {
            edgeRegion[edgeI] = markValue;
            changedEdges[nEdges++] = edgeI;
        }
    }
    changedEdges.setSize(nEdges);


    // Loop over changed faces. EdgeFaceWave in small.

    while (changedEdges.size())
    {
        DynamicList<label> changedFaces(changedEdges.size());

        forAll(changedEdges, i)
        {
            label edgeI = changedEdges[i];

            const labelList& eFaces = pp_.edgeFaces()[edgeI];

            forAll(eFaces, eFI)
            {
                label faceI =  eFaces[eFI];

                if (faceRegion[faceI] == -1)
                {
                    faceRegion[faceI] = markValue;
                    changedFaces.append(faceI);
                }
            }
        }

        // Loop over changedFaces and collect edges
        DynamicList<label> newChangedEdges(changedFaces.size());

        forAll(changedFaces, i)
        {
            label faceI = changedFaces[i];

            const labelList& fEdges = pp_.faceEdges()[faceI];

            forAll(fEdges, fEI)
            {
                label edgeI = fEdges[fEI];

                if (edgeRegion[edgeI] == -1)
                {
                    faceRegion[faceI] = markValue;
                    newChangedEdges.append(edgeI);
                }
            }
        }

        changedEdges.transfer(newChangedEdges);
    }
}


Foam::label Foam::regionSplitPatch::calcLocalRegionSplit
(
    const boolList& blockedEdge,
    labelList& faceRegion
) const
{
    // Region per edge.
    // -1 unassigned
    // -2 blocked
    labelList edgeRegion(pp_.nEdges(), -1);

    if (blockedEdge.size())
    {
        forAll(blockedEdge, edgeI)
        {
            if (blockedEdge[edgeI])
            {
                edgeRegion[edgeI] = -2;
            }
        }
    }

    // Assign local regions
    // ~~~~~~~~~~~~~~~~~~~~

    // Start with region 0
    label nLocalRegions = 0;

    label unsetFaceI = 0;

    do
    {
        // Find first unset face

        for (; unsetFaceI < pp_.size(); unsetFaceI++)
        {
            if (faceRegion[unsetFaceI] == -1)
            {
                break;
            }
        }

        if (unsetFaceI >= pp_.size())
        {
            break;
        }

        fillSeedMask
        (
            faceRegion,
            edgeRegion,
            unsetFaceI,
            nLocalRegions
        );

        // Current unsetFace has now been handled. Go to next region.
        nLocalRegions++;
        unsetFaceI++;
    }
    while (true);

    return nLocalRegions;
}


Foam::autoPtr<Foam::globalIndex> Foam::regionSplitPatch::calcRegionSplit
(
    const bool doGlobalRegions,
    const boolList& blockedEdge,
    labelList& faceRegion
) const
{
    // 1. Do local analysis
    label nLocalRegions = calcLocalRegionSplit
    (
        blockedEdge,
        faceRegion
    );

    if (!doGlobalRegions)
    {
        return autoPtr<globalIndex>(new globalIndex(nLocalRegions));
    }

    // 2. Assign global regions
    // ~~~~~~~~~~~~~~~~~~~~~~~~
    // Offset local regions to create unique global regions.

    globalIndex globalRegions(nLocalRegions);


    // Convert regions to global ones
    forAll(faceRegion, faceI)
    {
        faceRegion[faceI] = globalRegions.toGlobal(faceRegion[faceI]);
    }


    // 3. Merge global regions
    // ~~~~~~~~~~~~~~~~~~~~~~~
    // Regions across non-blocked proc patches get merged.
    // This will set merged global regions to be the min of both.
    // (this will create gaps in the global region list so they will get
    // merged later on)

    const globalMeshData& pd = mesh().globalData();
    const labelList& coupledEdges = pd.coupledPatchMeshEdges();
    boolList procEdges(mesh().nEdges(), false);

    forAll(coupledEdges, cEI)
    {
        label edgeI = coupledEdges[cEI];
        procEdges[edgeI] = true;
    }

    // Precalculate meshEdge per pp edge
    const labelList meshEdges
    (
        pp_.meshEdges(mesh().edges(), mesh().pointEdges())
    );

    while (true)
    {
        if (debug)
        {
            Pout<< nl << "-- Starting Iteration --" << endl;
        }

        labelList minRegion(pp_.nEdges(), labelMax);
        forAll(pp_.edges(), edgeI)
        {
            const labelList& eFaces = pp_.edgeFaces()[edgeI];
            forAll(eFaces, eFI)
            {
                label faceI = eFaces[eFI];
                minRegion[edgeI] = min(minRegion[edgeI], faceRegion[faceI]);
            }
        }

        syncTools::syncEdgeList
        (
            mesh(),
            meshEdges,
            minRegion,
            minEqOp<label>(),
            labelMax
        );

        Map<label> globalToMerged(pp_.size());
        forAll(pp_.edges(), edgeI)
        {
            label meshEdgeI = meshEdges[edgeI];
            if (procEdges[meshEdgeI])
            {
                if (blockedEdge.size() && !blockedEdge[edgeI])
                {
                    const labelList& eFaces = pp_.edgeFaces()[edgeI];
                    forAll(eFaces, eFI)
                    {
                        label faceI = eFaces[eFI];
                        if (minRegion[edgeI] < faceRegion[faceI])
                        {
                            globalToMerged.insert
                            (
                                faceRegion[faceI],
                                minRegion[edgeI]
                             );
                        }
                    }
                }
            }
        }

        label nMerged = returnReduce(globalToMerged.size(), sumOp<label>());

        if (debug)
        {
            Pout<< "nMerged:" << nMerged << endl;
        }

        if (nMerged == 0)
        {
            break;
        }

        // Renumber the regions according to the globalToMerged
        forAll(faceRegion, faceI)
        {
            label regionI = faceRegion[faceI];
            Map<label>::const_iterator iter = globalToMerged.find(regionI);
            if (iter != globalToMerged.end())
            {
                 faceRegion[faceI] = iter();
            }
        }
    }


    // Now our faceRegion will have non-local elements in it. So compact
    // it.

    // 4a: count. Use a labelHashSet to count regions only once.
    label nCompact = 0;
    {
        labelHashSet localRegion(pp_.size());
        forAll(faceRegion, faceI)
        {
            if
            (
                globalRegions.isLocal(faceRegion[faceI])
             && localRegion.insert(faceRegion[faceI])
            )
            {
                nCompact++;
            }
        }
    }

    if (debug)
    {
        Pout<< "Compacted from " << nLocalRegions
            << " down to " << nCompact << " local regions." << endl;
    }


    // Global numbering for compacted local regions
    autoPtr<globalIndex> globalCompactPtr(new globalIndex(nCompact));
    const globalIndex& globalCompact = globalCompactPtr();

    // 4b: renumber
    // Renumber into compact indices. Note that since we've already made
    // all regions global we now need a Map to store the compacting information
    // instead of a labelList - otherwise we could have used a straight
    // labelList.

    // Local compaction map
    Map<label> globalToCompact(2*nCompact);
    // Remote regions we want the compact number for
    List<labelHashSet> nonLocal(Pstream::nProcs());
    forAll(nonLocal, procI)
    {
        nonLocal[procI].resize((nLocalRegions-nCompact)/Pstream::nProcs());
    }

    forAll(faceRegion, faceI)
    {
        label region = faceRegion[faceI];
        if (globalRegions.isLocal(region))
        {
            Map<label>::const_iterator iter = globalToCompact.find(region);
            if (iter == globalToCompact.end())
            {
                label compactRegion = globalCompact.toGlobal
                (
                    globalToCompact.size()
                );
                globalToCompact.insert(region, compactRegion);
            }
        }
        else
        {
            nonLocal[globalRegions.whichProcID(region)].insert(region);
        }
    }

    // Now we have all the local regions compacted. Now we need to get the
    // non-local ones from the processors to whom they are local.
    // Convert the nonLocal (labelHashSets) to labelLists.

    labelListList sendNonLocal(Pstream::nProcs());
    labelList nNonLocal(Pstream::nProcs(), 0);
    forAll(sendNonLocal, procI)
    {
        sendNonLocal[procI].setSize(nonLocal[procI].size());
        forAllConstIter(labelHashSet, nonLocal[procI], iter)
        {
            sendNonLocal[procI][nNonLocal[procI]++] = iter.key();
        }
    }

    if (debug)
    {
        forAll(nNonLocal, procI)
        {
            Pout<< "    from processor " << procI
                << " want " << nNonLocal[procI]
                << " region numbers."
                << endl;
        }
        Pout<< endl;
    }


    // Get the wanted region labels into recvNonLocal
    labelListList recvNonLocal;
    Pstream::exchange<labelList, label>
    (
        sendNonLocal,
        recvNonLocal
    );

    // Now we have the wanted compact region labels that procI wants in
    // recvNonLocal[procI]. Construct corresponding list of compact
    // region labels to send back.

    labelListList sendWantedLocal(Pstream::nProcs());
    forAll(recvNonLocal, procI)
    {
        const labelList& nonLocal = recvNonLocal[procI];
        sendWantedLocal[procI].setSize(nonLocal.size());

        forAll(nonLocal, i)
        {
            sendWantedLocal[procI][i] = globalToCompact[nonLocal[i]];
        }
    }


    // Send back (into recvNonLocal)
    recvNonLocal.clear();
    Pstream::exchange<labelList, label>
    (
        sendWantedLocal,
        recvNonLocal
    );
    sendWantedLocal.clear();

    // Now recvNonLocal contains for every element in setNonLocal the
    // corresponding compact number. Insert these into the local compaction
    // map.

    forAll(recvNonLocal, procI)
    {
        const labelList& wantedRegions = sendNonLocal[procI];
        const labelList& compactRegions = recvNonLocal[procI];

        forAll(wantedRegions, i)
        {
            globalToCompact.insert(wantedRegions[i], compactRegions[i]);
        }
    }

    // Finally renumber the regions
    forAll(faceRegion, faceI)
    {
        faceRegion[faceI] = globalToCompact[faceRegion[faceI]];
    }

    return globalCompactPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::regionSplitPatch::regionSplitPatch
(
    const polyMesh& mesh,
    const indirectPrimitivePatch& pp,
    const boolList& blockedEdge,
    const bool doGlobalRegions
)
:
    MeshObject<polyMesh, Foam::TopologicalMeshObject, regionSplitPatch>(mesh),
    labelList(pp.size(), -1),
    pp_(pp)
{
    globalNumberingPtr_ = calcRegionSplit
    (
        doGlobalRegions,
        blockedEdge,      //blockedEdge
        *this
    );
}
