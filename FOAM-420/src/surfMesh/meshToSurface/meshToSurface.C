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

#include "meshToSurface/meshToSurface.H"
#include "MeshedSurface/MeshedSurface.H"
#include "UnsortedMeshedSurface/UnsortedMeshedSurface.H"
#include "db/Time/Time.H"
#include "meshes/polyMesh/polyMesh.H"
#include "meshes/polyMesh/polyPatches/constraint/processor/processorPolyPatch.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "meshes/primitiveMesh/primitivePatch/uindirectPrimitivePatch.H"
#include "meshes/meshTools/simpleVTKWriter.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshToSurface, 0);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshToSurface::meshToSurface
(
    const polyMesh& mesh,
    const labelHashSet& patchHash,
    const labelHashSet& zoneHash,
    const bool& consistent
)
:
    mesh_(mesh),
    patchHash_(patchHash),
    zoneHash_(zoneHash),
    consistent_(consistent)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meshToSurface::write(const fileName& outFileName) const
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();
    const faceZoneMesh& fzm = mesh_.faceZones();

    // From (name of) patch to compact 'zone' index
    HashTable<label> compactZoneID(1024);
    // Mesh face and compact zone indx
    DynamicList<label> faceLabels;
    DynamicList<bool> faceFlip;
    DynamicList<label> compactZones;

    {
        // Collect sizes. Hash on names to handle local-only patches (e.g.
        //  processor patches)
        HashTable<label> patchSize(1024);
        label nFaces = 0;
        forAllConstIter(labelHashSet, patchHash_, iter)
        {
            const polyPatch& pp = bMesh[iter.key()];
            patchSize.insert(pp.name(), pp.size());
            nFaces += pp.size();
        }

        HashTable<label> zoneSize(1024);
        forAllConstIter(labelHashSet, zoneHash_, iter)
        {
            const faceZone& pp = fzm[iter.key()];
            zoneSize.insert(pp.name(), pp.size());
            nFaces += pp.size();
        }


        Pstream::mapCombineGather(patchSize, plusEqOp<label>());
        Pstream::mapCombineGather(zoneSize, plusEqOp<label>());


        // Allocate compact numbering for all patches/faceZones
        forAllConstIter(HashTable<label>, patchSize, iter)
        {
            label sz = compactZoneID.size();
            compactZoneID.insert(iter.key(), sz);
        }

        forAllConstIter(HashTable<label>, zoneSize, iter)
        {
            label sz = compactZoneID.size();
            //Info<< "For faceZone " << iter.key() << " allocating zoneID "
            //    << sz << endl;
            compactZoneID.insert(iter.key(), sz);
        }


        Pstream::mapCombineScatter(compactZoneID);


        // Rework HashTable into labelList just for speed of conversion
        labelList patchToCompactZone(bMesh.size(), -1);
        labelList faceZoneToCompactZone(bMesh.size(), -1);
        forAllConstIter(HashTable<label>, compactZoneID, iter)
        {
            label patchi = bMesh.findPatchID(iter.key());
            if (patchi != -1)
            {
                patchToCompactZone[patchi] = iter();
            }
            else
            {
                label zoneI = fzm.findZoneID(iter.key());
                faceZoneToCompactZone[zoneI] = iter();
            }
        }


        faceLabels.setCapacity(nFaces);
        faceFlip.setCapacity(nFaces);
        compactZones.setCapacity(nFaces);

        // Collect faces on patches
        forAllConstIter(labelHashSet, patchHash_, iter)
        {
            const polyPatch& pp = bMesh[iter.key()];
            forAll(pp, i)
            {
                faceLabels.append(pp.start()+i);
                faceFlip.append(false);
                compactZones.append(patchToCompactZone[pp.index()]);
            }
        }
        // Collect faces on faceZones
        forAllConstIter(labelHashSet, zoneHash_, iter)
        {
            const faceZone& pp = fzm[iter.key()];
            forAll(pp, i)
            {
                faceLabels.append(pp[i]);
                faceFlip.append(pp.flipMap()[i]);
                compactZones.append(faceZoneToCompactZone[pp.index()]);
            }
        }
    }

    UIndirectList<face> ppFaces (mesh_.faces(), faceLabels);
    forAll(faceFlip, fI)
    {
        if (faceFlip[fI]) ppFaces[fI].flip();
    }

    // Addressing engine for all faces
    uindirectPrimitivePatch allBoundary
    (
        ppFaces,
        mesh_.points()
    );


    {
        simpleVTKWriter writeMaster
        (
            allBoundary.localFaces(),
            allBoundary.localPoints()
        );
        writeMaster.addFaceData("normal", allBoundary.faceNormals());
        writeMaster.write("writeMaster.vtk");
    }


    // Find correspondence to master points
    labelList pointToGlobal;
    labelList uniqueMeshPoints;
    autoPtr<globalIndex> globalNumbers = mesh_.globalData().mergePoints
    (
        allBoundary.meshPoints(),
        pointToGlobal,
        uniqueMeshPoints
    );

    // Gather all unique points on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()] = pointField
    (
        mesh_.points(),
        uniqueMeshPoints
    );
    Pstream::gatherList(gatheredPoints);

    // Gather all faces
    List<faceList> gatheredFaces(Pstream::nProcs());
    gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
    forAll(gatheredFaces[Pstream::myProcNo()], i)
    {
        inplaceRenumber
       (
            pointToGlobal,
            gatheredFaces[Pstream::myProcNo()][i]
       );
    }
    Pstream::gatherList(gatheredFaces);

    // Gather all ZoneIDs
    List<labelList> gatheredZones(Pstream::nProcs());
    gatheredZones[Pstream::myProcNo()] = compactZones.xfer();
    Pstream::gatherList(gatheredZones);

    // On master combine all points, faces, zones
    if (Pstream::master())
    {
        pointField allPoints = ListListOps::combine<pointField>
        (
            gatheredPoints,
            accessOp<pointField>()
        );
        gatheredPoints.clear();

        faceList allFaces = ListListOps::combine<faceList>
        (
            gatheredFaces,
            accessOp<faceList>()
        );
        gatheredFaces.clear();

        labelList allZones = ListListOps::combine<labelList>
        (
            gatheredZones,
            accessOp<labelList>()
        );
        gatheredZones.clear();

        Map<label> reorder(compactZoneID.size());

        // Zones
        surfZoneIdentifierList surfZones(compactZoneID.size());
        forAllConstIter(HashTable<label>, compactZoneID, iter)
        {
            surfZones[iter()] = surfZoneIdentifier(iter.key(), iter());
            if (!consistent_)
            {
                Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
                    << endl;
            }
            else
            {
                label patchI = bMesh.findPatchID(iter.key());
                reorder.insert(iter(),patchI);
            }
        }

        if (consistent_)
        {
            surfZoneIdentifierList reorderedSurfZones(compactZoneID.size());
            forAllConstIter(HashTable<label>, compactZoneID, iter)
            {
                label patchI = reorder(iter());
                reorderedSurfZones[patchI] = surfZoneIdentifier(iter.key(), patchI);
            }

            labelList reorderedAllZones(allZones.size());
            forAll(allZones, i)
            {
                reorderedAllZones[i] = reorder(allZones[i]);
            }
            surfZones.clear();
            surfZones = reorderedSurfZones;

            forAllConstIter(HashTable<label>, compactZoneID, iter)
            {
                    Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
                    << endl;
            }

            allZones.clear();
            allZones = reorderedAllZones;
        }

        UnsortedMeshedSurface<face> unsortedFace
        (
            xferMove(allPoints),
            xferMove(allFaces),
            xferMove(allZones),
            xferMove(surfZones)
        );

        MeshedSurface<face> sortedFace(unsortedFace);

        const Time& time = mesh_.time();

        fileName globalCasePath
        (
            outFileName.isAbsolute()
          ? outFileName
          : (
                time.processorCase()
              ? time.rootPath()/time.globalCaseName()/outFileName
              : time.path()/outFileName
            )
        );
        globalCasePath.clean();

        Info<< "Writing merged surface to " << globalCasePath << endl;

        sortedFace.write(globalCasePath);
    }
}


// ************************************************************************* //
