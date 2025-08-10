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
    (c) 2019 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "autoBlockMesh/autoBlockMesh.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "meshRefinement/meshRefinement.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddPoint.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddFace.H"
#include "polyTopoChange/polyTopoChange/addObject/polyAddCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(autoBlockMesh, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::autoBlockMesh::autoBlockMesh
(
    const refinementSurfaces& surfaces,
    const dictionary& meshDict,
    fvMesh& mesh
)
:
    mesh_(mesh),
    surfaces_(surfaces),
    coordFramePtr_(nullptr),
    blockData_(meshDict.lookup("blockData")),
    primitiveBounding_
    (
        meshDict.lookupOrDefault<Switch>("primitiveBounding",true)
    ),
    blockIso_
    (
        meshDict.lookupOrDefault<Switch>("blockIso",false)
    ),
    blockOrigin_
    (
        meshDict.lookupOrDefault("blockOrigin",vector(GREAT,GREAT,GREAT))
    )
{
    if (meshDict.found("referenceFrame"))
    {
        coordFramePtr_ = coordinateFrame::lookupNew(mesh_, meshDict);
    }
    else
    {
        dictionary coordDict(getCoordinateDict(meshDict));
        coordSysPtr_ = coordinateSystem::New("cartesian",coordDict);
    }

    const coordinateSystem& coord =
        coordFramePtr_ ? coordFramePtr_->coorSys() : *coordSysPtr_;
    if
    (
        mag(coord.e1() & coord.e2()) > 0.01
        || mag(coord.e1() & coord.e3()) > 0.01
        || mag(coord.e2() & coord.e3()) > 0.01
    )
    {
        FatalErrorInFunction
            << "Non orthogonal local coordinate system used for auto blockMesh"
            << " e1: " << coord.e1() << " e2: " << coord.e2()
            << " e3: " << coord.e1() << endl;
    }
    else
    {
        Info<<"Using local coordinates "
            << " e1: " << coord.e1() << " e2: " << coord.e2()
            << " e3: " << coord.e1() << endl;
    }

    createBlock();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary Foam::autoBlockMesh::getCoordinateDict
(
    const dictionary& meshDict
)
{
    if (meshDict.found("coordinates"))
    {
        return meshDict.subDict("coordinates");
    }
    else
    {
        dictionary defaultDict("coordinates");

        defaultDict.add("origin",vector::zero,true);
        defaultDict.add(word("rotation"), dictionary(), false);
        dictionary& rotDict = defaultDict.subDict("rotation");
        rotDict.add("type","axesRotation",true);
        rotDict.add("e1",vector(1,0,0),true);
        rotDict.add("e2",vector(0,1,0),true);

        return defaultDict;
    }
}


void Foam::autoBlockMesh::createBlock()
{
    vector greatPoint = vector(GREAT, GREAT, GREAT);

    point minPt = greatPoint;
    point maxPt =-greatPoint;

    const coordinateSystem& coord =
        coordFramePtr_ ? coordFramePtr_->coorSys() : *coordSysPtr_;

    surfaces_.getFiniteSurfacesLocalBb
    (
        coord,
        primitiveBounding_,
        minPt,
        maxPt
    );

    vector dxyz = maxPt - minPt;
    point midPt = 0.5*(maxPt+minPt);

    if (blockOrigin_ != greatPoint)
    {
        point localOrigin = coord.invTransformPoint(blockOrigin_);

        scalar maxX1 = mag(maxPt.x()-localOrigin.x());
        scalar maxX2 = mag(minPt.x()-localOrigin.x());

        scalar maxY1 = mag(maxPt.y()-localOrigin.y());
        scalar maxY2 = mag(minPt.y()-localOrigin.y());

        scalar maxZ1 = mag(maxPt.z()-localOrigin.z());
        scalar maxZ2 = mag(minPt.z()-localOrigin.z());

        vector maxBB = vector
        (
            max(maxX1,maxX2),
            max(maxY1,maxY2),
            max(maxZ1,maxZ2)
        );

        maxPt = localOrigin + maxBB;
        minPt = localOrigin - maxBB;

        dxyz = maxPt - minPt;
        midPt = 0.5*(maxPt+minPt);
    }

    Info<< nl
        << "autoBlockMesh bounding box  : "
        << coord.transformPoint(minPt) << " "
        << " " << coord.transformPoint(maxPt) << nl
        << " autoBlockMesh bounding box centre : "
        << coord.transformPoint(midPt)  << endl;

    //Calculate length scale and cell level
    scalar length = blockData_.first();
    label level = blockData_.second();
    //Convert length to level 0
    scalar cell0Len = length * (1<<level);

    labelVector n = (dxyz/cell0Len) + vector(0.5);

    vector blockSizes;
    if (blockIso_)
    {
        for (direction i = 0; i < vector::nComponents; i++)
        {
            minPt[i] = midPt[i] - 0.5*(n[i]*cell0Len);
            blockSizes[i] = cell0Len;
        }
    }
    else
    {
        if (cmptMin(n) == 0)
        {
            n = max(n, labelVector::one);
            scalar minLength = cmptMin(dxyz);
            WarningInFunction
                << "Block mesh length scale : " << length
                << " too large for minimum bounding box size : " << minLength
                << " will add a single cell in this dimension."
                << endl;
        }
        for (direction i = 0; i < vector::nComponents; i++)
        {
            blockSizes[i] = dxyz[i] / n[i];
        }
    }

    //Add single cell external to bounding box
    n += labelVector(2);

    minPt -= blockSizes;

    //Calculate global corner point
    point cornerGlobalPt = coord.transformPoint(minPt);

    SortableList<label> sorted(n.size());
    forAll(n, i)
    {
        sorted[i] = n[i];
    }
    sorted.reverseSort();
    const labelList& indices = sorted.indices();

    //Calculate slices in processor direction
    label nMax = n[indices[0]];
    label nCellsPerProc = nMax/Pstream::nProcs();
    label remain(nMax % Pstream::nProcs());

    if (Pstream::myProcNo() < remain)
    {
        nCellsPerProc++;
    }

    label nActive = 0;
    if (nCellsPerProc > 0)
    {
        nActive++;
    }
    reduce(nActive, sumOp<label>());

    //Add  default wall patch
    {
        dictionary patchInfo;
        patchInfo.add("type", wallPolyPatch::typeName);
        patchInfo.add("nFaces", 0);
        patchInfo.add("startFace", 0);
        addEmptyPatch("blockMesh", patchInfo);
    }

    label nProcStart = 1;

    DynamicList<label> nbrProcs(2);
    if (nCellsPerProc > 0 && nActive > 1)
    {
        if (Pstream::myProcNo() == 0)
        {
            nbrProcs.append(1);
        }
        else if (Pstream::myProcNo() == Pstream::nProcs() -1)
        {
            nbrProcs.append(Pstream::myProcNo()-1);
        }
        else if
        (
            nActive != Pstream::nProcs()
            && Pstream::myProcNo() == remain-1
        )
        {
            nbrProcs.append(Pstream::myProcNo()-1);
        }
        else
        {
            nbrProcs.append(Pstream::myProcNo()-1);
            nbrProcs.append(Pstream::myProcNo()+1);
        }

        forAll(nbrProcs, nbri)
        {
            word procName
            (
                processorPolyPatch::newName(Pstream::myProcNo(), nbrProcs[nbri])
            );

            dictionary patchInfo;
            patchInfo.add("type", processorPolyPatch::typeName);
            patchInfo.add("myProcNo", Pstream::myProcNo());
            patchInfo.add("neighbProcNo", nbrProcs[nbri]);
            patchInfo.add("nFaces", 0);
            patchInfo.add("startFace", 0);

            addEmptyPatch(procName, patchInfo);
        }
    }

    //Calculate start index for each parallel block and corner point
    label startIndex = Pstream::myProcNo() * (nMax/Pstream::nProcs())
        + min(Pstream::myProcNo(),remain);

    List<vector> sortedDirVec(3);
    labelVector sortedn;
    vector sortedSpacing;

    if (indices[0] == 0)
    {
        sortedDirVec[0] =  coord.e1();
        sortedDirVec[1] =  coord.e2();
        sortedDirVec[2] =  coord.e3();
        sortedn[0] = nCellsPerProc;
        sortedn[1] = n[1];
        sortedn[2] = n[2];

        sortedSpacing[0] = blockSizes[0];
        sortedSpacing[1] = blockSizes[1];
        sortedSpacing[2] = blockSizes[2];
    }
    else if (indices[0] == 1)
    {
        sortedDirVec[0] =  coord.e2();
        sortedDirVec[1] =  coord.e3();
        sortedDirVec[2] =  coord.e1();
        sortedn[0] = nCellsPerProc;
        sortedn[1] = n[2];
        sortedn[2] = n[0];

        sortedSpacing[0] = blockSizes[1];
        sortedSpacing[1] = blockSizes[2];
        sortedSpacing[2] = blockSizes[0];
    }
    else
    {
        sortedDirVec[0] =  coord.e3();
        sortedDirVec[1] =  coord.e1();
        sortedDirVec[2] =  coord.e2();
        sortedn[0] = nCellsPerProc;
        sortedn[1] = n[0];
        sortedn[2] = n[1];

        sortedSpacing[0] = blockSizes[2];
        sortedSpacing[1] = blockSizes[0];
        sortedSpacing[2] = blockSizes[1];
    }

    cornerGlobalPt += startIndex*sortedDirVec[0]*sortedSpacing[0];

    polyTopoChange meshMod(mesh_);

    if (sortedn[0] > 0)
    {
        //Add the cells
        labelList cells(sortedn[0]*sortedn[1]*sortedn[2]);

        for (label k=0; k< sortedn[2]; k++)
        {
            for (label j=0; j< sortedn[1]; j++)
            {
                for (label i=0; i< sortedn[0]; i++)
                {
                    cells[cellLabel(i,j,k,sortedn)] = meshMod.setAction
                    (
                        polyAddCell()
                    );
                }
            }
        }

        //Add the points
        labelList pts((sortedn[0]+1)*(sortedn[1]+1)*(sortedn[2]+1));

        for (label k=0; k< sortedn[2]+1; k++)
        {
            for (label j=0; j< sortedn[1]+1; j++)
            {
                for (label i=0; i< sortedn[0]+1; i++)
                {
                    point pt =
                    (
                        cornerGlobalPt
                        + sortedDirVec[0]*i*sortedSpacing[0]
                        + sortedDirVec[1]*j*sortedSpacing[1]
                        + sortedDirVec[2]*k*sortedSpacing[2]
                    );

                    pts[pointLabel(i,j,k,sortedn)] = meshMod.setAction
                    (
                        polyAddPoint(pt,-1,-1,true)
                    );
                }
            }
        }

        //Add the faces

        //First faces in i-direction

        for (label k=0; k< sortedn[2]; k++)
        {
            for (label j=0; j< sortedn[1]; j++)
            {
                for (label i=0; i< sortedn[0]+1; i++)
                {
                    face f(4);
                    f[0] = pts[pointLabel(i,j,k,sortedn)];
                    f[1] = pts[pointLabel(i,j+1,k,sortedn)];
                    f[2] = pts[pointLabel(i,j+1,k+1,sortedn)];
                    f[3] = pts[pointLabel(i,j,k+1,sortedn)];

                    label patchi = -1;
                    label own = -1;
                    label nei = -1;
                    bool flipFace = false;

                    if (i == 0)
                    {
                        own = cells[cellLabel(i,j,k,sortedn)];
                        if (nbrProcs.size() != 2)
                        {
                            if (Pstream::myProcNo() == 0)
                            {
                                patchi = 0;
                            }
                            else
                            {
                                patchi = nProcStart;
                            }
                        }
                        else
                        {
                            patchi = nProcStart;
                        }
                        flipFace = true;
                    }
                    else if (i == sortedn[0])
                    {
                        own = cells[cellLabel(i-1,j,k,sortedn)];
                        if (nbrProcs.size() != 2)
                        {
                            if (Pstream::myProcNo() == 0 && nbrProcs.size() == 1)
                            {
                                patchi = nProcStart;
                            }
                            else
                            {
                                patchi = 0;
                            }
                        }
                        else
                        {
                            patchi = nProcStart+1;
                        }
                    }
                    else
                    {
                        own = cells[cellLabel(i-1,j,k,sortedn)];
                        nei = cells[cellLabel(i,j,k,sortedn)];
                    }

                    if (flipFace)
                    {
                        f.flip();
                    }

                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            f,   // face
                            own,       // owner
                            nei,  // neighbour
                            -1,        // master point
                            -1,        // master edge
                            -1,     // master face
                            false,     // flux flip
                            patchi,        // patch for face
                            -1,     // zone for face
                            false       // face zone flip
                         )
                    );
                }
            }
        }

        //Add faces in j-direction
        for (label k=0; k< sortedn[2]; k++)
        {
            for (label i=0; i< sortedn[0]; i++)
            {
                for (label j=0; j< sortedn[1]+1; j++)
                {
                    face f(4);
                    f[0] = pts[pointLabel(i,j,k,sortedn)];
                    f[1] = pts[pointLabel(i,j,k+1,sortedn)];
                    f[2] = pts[pointLabel(i+1,j,k+1,sortedn)];
                    f[3] = pts[pointLabel(i+1,j,k,sortedn)];

                    label patchi = -1;
                    label own = -1;
                    label nei = -1;
                    bool flipFace = false;

                    if (j == 0)
                    {
                        own = cells[cellLabel(i,j,k,sortedn)];
                        patchi = 0;
                        flipFace = true;
                    }
                    else if (j == sortedn[1])
                    {
                        own = cells[cellLabel(i,j-1,k,sortedn)];
                        patchi = 0;
                    }
                    else
                    {
                        own = cells[cellLabel(i,j-1,k,sortedn)];
                        nei = cells[cellLabel(i,j,k,sortedn)];
                    }

                    if (flipFace)
                    {
                        f.flip();
                    }

                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            f,   // face
                            own,       // owner
                            nei,  // neighbour
                            -1,        // master point
                            -1,        // master edge
                            -1,     // master face
                            false,     // flux flip
                            patchi,        // patch for face
                            -1,     // zone for face
                            false       // face zone flip
                         )
                    );
                }
            }
        }

        //Add faces in k-direction
        for (label i=0; i< sortedn[0]; i++)
        {
            for (label j=0; j< sortedn[1]; j++)
            {
                for (label k=0; k< sortedn[2]+1; k++)
                {
                    face f(4);
                    f[0] = pts[pointLabel(i,j,k,sortedn)];
                    f[1] = pts[pointLabel(i+1,j,k,sortedn)];
                    f[2] = pts[pointLabel(i+1,j+1,k,sortedn)];
                    f[3] = pts[pointLabel(i,j+1,k,sortedn)];

                    label patchi = -1;
                    label own = -1;
                    label nei = -1;
                    bool flipFace = false;

                    if (k == 0)
                    {
                        own = cells[cellLabel(i,j,k,sortedn)];
                        patchi = 0;
                        flipFace = true;
                    }
                    else if (k == sortedn[2])
                    {
                        own = cells[cellLabel(i,j,k-1,sortedn)];
                        patchi = 0;
                    }
                    else
                    {
                        own = cells[cellLabel(i,j,k-1,sortedn)];
                        nei = cells[cellLabel(i,j,k,sortedn)];
                    }

                    if (flipFace)
                    {
                        f.flip();
                    }

                    meshMod.setAction
                    (
                        polyAddFace
                        (
                            f,   // face
                            own,       // owner
                            nei,  // neighbour
                            -1,        // master point
                            -1,        // master edge
                            -1,     // master face
                            false,     // flux flip
                            patchi,        // patch for face
                            -1,     // zone for face
                            false       // face zone flip
                        )
                    );
                }
            }
        }
    }

    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false);

    // Optionally inflate mesh
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes.
        mesh_.clearOut();
    }

    return;
}


void  Foam::autoBlockMesh::addEmptyPatch
(
    const word& patchName,
    const dictionary& patchInfo
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh_.boundaryMesh());
    fvBoundaryMesh& fvPatches =
        const_cast<fvBoundaryMesh&>(mesh_.boundary());

    label patchi = polyPatches.size();

    // Add polyPatch at the end
    polyPatches.setSize(patchi+1);
    polyPatches.set
    (
        patchi,
        polyPatch::New
        (
            patchName,
            patchInfo,
            patchi,
            polyPatches
        )
    );
    fvPatches.setSize(patchi+1);
    fvPatches.set
    (
        patchi,
        fvPatch::New
        (
            polyPatches[patchi],  // point to newly added polyPatch
            mesh_.boundary()
        )
    );
}


inline Foam::label Foam::autoBlockMesh::pointLabel
(
    const label i,
    const label j,
    const label k,
    const labelVector n
) const
{
    return
    (
        i
      + j*(n.x() + 1)
      + k*(n.x() + 1)*(n.y() + 1)
    );
}


inline Foam::label Foam::autoBlockMesh::cellLabel
(
    const label i,
    const label j,
    const label k,
    const labelVector n
) const
{
    return
    (
        i
      + j*n.x()
      + k*(n.x()*n.y())
    );
}

// ************************************************************************* //
