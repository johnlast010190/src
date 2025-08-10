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

\*---------------------------------------------------------------------------*/

#include "orientFaceZone/orientFaceZone.H"
#include "meshes/polyMesh/syncTools/syncTools.H"
#include "meshRefinement/patchFaceOrientation.H"
#include "algorithms/PatchEdgeFaceWave/PatchEdgeFaceWave.H"
#include "triSurface/orientedSurface/orientedSurface.H"
#include "algorithms/PatchEdgeFaceWave/patchEdgeFaceRegion.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(orientFaceZone, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::orientFaceZone::orientFaceZone
(
    polyMesh& mesh,
    const word& fzName
)
:
    mesh_(mesh),
    fzName_(fzName)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::orientFaceZone::orientByPoint(const point& outsidePoint)
{
    Info<< "Orienting faceZone " << fzName_
        << " such that " << outsidePoint << " is outside"
        << nl << endl;


    const faceZone& fZone = mesh_.faceZones()[fzName_];

    if (fZone.checkParallelSync())
    {
        FatalErrorInFunction
            << "Face zone " << fZone.name()
            << " is not parallel synchronised."
            << " Any coupled face also needs its coupled version to be included"
            << " and with opposite flipMap."
            << exit(FatalError);
    }

    const labelList& faceLabels = fZone;

    const indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), faceLabels),
        mesh_.points()
    );



    const PackedBoolList isMasterFace(syncTools::getMasterFaces(mesh_));


    // Data on all edges and faces
    List<patchFaceOrientation> allEdgeInfo(patch.nEdges());
    List<patchFaceOrientation> allFaceInfo(patch.size());

    // Make sure we don't walk through
    // - slaves of coupled faces
    // - non-manifold edges
    {
        const polyBoundaryMesh& bm = mesh_.boundaryMesh();

        label nProtected = 0;

        forAll(faceLabels, facei)
        {
            const label meshFacei = faceLabels[facei];
            const label patchi = bm.whichPatch(meshFacei);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFacei]
            )
            {
                // Slave side. Mark so doesn't get visited.
                allFaceInfo[facei] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }

        Info<< "Protected from visiting "
            << returnReduce(nProtected, sumOp<label>())
            << " slaves of coupled faces" << nl << endl;
    }
    {
        // Number of (master)faces per edge
        labelList nMasterFaces(patch.nEdges(), 0);

        forAll(faceLabels, facei)
        {
            const label meshFacei = faceLabels[facei];

            if (isMasterFace[meshFacei])
            {
                const labelList& fEdges = patch.faceEdges()[facei];
                forAll(fEdges, fEdgeI)
                {
                    nMasterFaces[fEdges[fEdgeI]]++;
                }
            }
        }

        syncTools::syncEdgeList
        (
            mesh_,
            patch.meshEdges(mesh_.edges(), mesh_.pointEdges()),
            nMasterFaces,
            plusEqOp<label>(),
            label(0)
        );


        label nProtected = 0;

        forAll(nMasterFaces, edgeI)
        {
            if (nMasterFaces[edgeI] > 2)
            {
                allEdgeInfo[edgeI] = orientedSurface::NOFLIP;
                nProtected++;
            }
        }

        Info<< "Protected from visiting "
            << returnReduce(nProtected, sumOp<label>())
            << " non-manifold edges" << nl << endl;
    }



    DynamicList<label> changedEdges;
    DynamicList<patchFaceOrientation> changedInfo;

    const scalar tol = PatchEdgeFaceWave
    <
        indirectPrimitivePatch,
        patchFaceOrientation
    >::propagationTol();

    int dummyTrackData;

    globalIndex globalFaces(patch.size());

    while (true)
    {
        // Pick an unset face
        label unsetFacei = labelMax;
        forAll(allFaceInfo, facei)
        {
            if (allFaceInfo[facei] == orientedSurface::UNVISITED)
            {
                unsetFacei = globalFaces.toGlobal(facei);
                break;
            }
        }

        reduce(unsetFacei, minOp<label>());

        if (unsetFacei == labelMax)
        {
            break;
        }

        label proci = globalFaces.whichProcID(unsetFacei);
        label seedFacei = globalFaces.toLocal(proci, unsetFacei);
        Info<< "Seeding from processor " << proci << " face " << seedFacei
            << endl;

        if (proci == Pstream::myProcNo())
        {
            // Determine orientation of seedFace

            vector d = outsidePoint-patch.faceCentres()[seedFacei];
            const vector& fn = patch.faceNormals()[seedFacei];

            // Set information to correct orientation
            patchFaceOrientation& faceInfo = allFaceInfo[seedFacei];
            faceInfo = orientedSurface::NOFLIP;

            if ((fn&d) < 0)
            {
                faceInfo.flip();

                Pout<< "Face " << seedFacei << " at "
                    << patch.faceCentres()[seedFacei]
                    << " with normal " << fn
                    << " needs to be flipped." << endl;
            }
            else
            {
                Pout<< "Face " << seedFacei << " at "
                    << patch.faceCentres()[seedFacei]
                    << " with normal " << fn
                    << " points in positive direction (cos = " << (fn&d)/mag(d)
                    << ")" << endl;
            }


            const labelList& fEdges = patch.faceEdges()[seedFacei];
            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                patchFaceOrientation& edgeInfo = allEdgeInfo[edgeI];

                if
                (
                    edgeInfo.updateEdge<int>
                    (
                        mesh_,
                        patch,
                        edgeI,
                        seedFacei,
                        faceInfo,
                        tol,
                        dummyTrackData
                    )
                )
                {
                    changedEdges.append(edgeI);
                    changedInfo.append(edgeInfo);
                }
            }
        }


        if (returnReduce(changedEdges.size(), sumOp<label>()) == 0)
        {
            break;
        }



        // Walk
        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            patchFaceOrientation
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );
    }


    // Push master zone info over to slave (since slave faces never visited)
    {
        const polyBoundaryMesh& bm = mesh_.boundaryMesh();

        labelList neiStatus
        (
            mesh_.nFaces()-mesh_.nInternalFaces(),
            orientedSurface::UNVISITED
        );

        forAll(faceLabels, i)
        {
            const label meshFacei = faceLabels[i];
            if (!mesh_.isInternalFace(meshFacei))
            {
                neiStatus[meshFacei-mesh_.nInternalFaces()] =
                    allFaceInfo[i].flipStatus();
            }
        }
        syncTools::swapBoundaryFaceList(mesh_, neiStatus);

        forAll(faceLabels, i)
        {
            const label meshFacei = faceLabels[i];
            const label patchi = bm.whichPatch(meshFacei);

            if
            (
                patchi != -1
             && bm[patchi].coupled()
             && !isMasterFace[meshFacei]
            )
            {
                // Slave side. Take flipped from neighbour
                label bFacei = meshFacei-mesh_.nInternalFaces();

                if (neiStatus[bFacei] == orientedSurface::NOFLIP)
                {
                    allFaceInfo[i] = orientedSurface::FLIP;
                }
                else if (neiStatus[bFacei] == orientedSurface::FLIP)
                {
                    allFaceInfo[i] = orientedSurface::NOFLIP;
                }
                else
                {
                    FatalErrorInFunction
                        << "Incorrect status for face " << meshFacei
                        << abort(FatalError);
                }
            }
        }
    }


    // Convert to flipmap and adapt faceZones

    boolList newFlipMap(allFaceInfo.size(), false);
    label nChanged = 0;
    forAll(allFaceInfo, facei)
    {
        if (allFaceInfo[facei] == orientedSurface::NOFLIP)
        {
            newFlipMap[facei] = false;
        }
        else if (allFaceInfo[facei] == orientedSurface::FLIP)
        {
            newFlipMap[facei] = true;
        }
        else
        {
            FatalErrorInFunction
                << "Problem : unvisited face " << facei
                << " centre:" << mesh_.faceCentres()[faceLabels[facei]]
                << abort(FatalError);
        }

        if (fZone.flipMap()[facei] != newFlipMap[facei])
        {
            nChanged++;
        }
    }

    reduce(nChanged, sumOp<label>());
    if (nChanged > 0)
    {
        Info<< "Flipping " << nChanged << " out of "
            << globalFaces.size() << " faces." << nl << endl;

        mesh_.faceZones()[fzName_].resetAddressing(faceLabels, newFlipMap);
        if (!mesh_.faceZones().write())
        {
            FatalErrorInFunction
                << "Failed writing faceZones" << exit(FatalError);
        }
    }
}


// ************************************************************************* //
