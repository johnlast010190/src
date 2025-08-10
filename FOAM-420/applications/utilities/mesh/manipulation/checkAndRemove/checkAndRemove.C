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
    (c) 2010-2012 Esi Ltd.
    (c) 1991-2007 OpenCFD Ltd.

Application
    checkAndRemove

Description
    Cell removal engine.

\*---------------------------------------------------------------------------*/

#include "meshRefinement/meshRefinement.H"
#include "motionSmoother/motionSmoother.H"
#include "fvMeshSubset/fvMeshSubset.H"
#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "snappyHexMeshDriver/keepData/keepData.H"
#include "algorithms/FaceCellWave/FaceCellWave.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "db/IOobjectList/IOobjectList.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "decompositionMethod/decompositionMethod.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"

using namespace Foam;


// Main program:

int main(int argc, char *argv[])
{
#include "include/addOverwriteOption.H"
#include "include/addTimeOptions.H"
#include "include/setRootCase.H"
#include "include/createTime.H"
#include "include/addOverwriteOption.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#include "include/checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

    runTime.functionObjects().off();
#include "include/createMesh.H"

    const word oldInstance = mesh.pointsInstance();
    bool overwrite = args.options().found("overwrite");

    fileName dictFile = fileName("cellRemovalDict");

    IOobject removalHeader
    (
        dictFile,
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    bool snappyHex = false;
    if (!removalHeader.typeHeaderOk<IOdictionary>(true))
    {
        dictFile = "snappyHexMeshDict";

        IOobject snappyDictHeader
        (
            dictFile,
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (!snappyDictHeader.typeHeaderOk<IOdictionary>(true))
        {
            FatalError << "No snappyHexMeshDict or cellRemovalDict present. "
                       << exit(FatalError);
        }

        snappyHex = true;
    }

    // Read setting dictionary
    IOdictionary removalDict
    (
        IOobject
        (
            dictFile,
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    Switch meshQuality(true);
    if (removalDict.found("meshQualityRemoval"))
    {
        meshQuality = readBool(removalDict.lookup("meshQualityRemoval"));
    }

    Switch multiRegion(true);
    if (removalDict.found("multiRegionRemoval"))
    {
        multiRegion = readBool(removalDict.lookup("multiRegionRemoval"));
    }

    Switch holeDetection(true);
    if (removalDict.found("holeDetection"))
    {
        holeDetection = readBool(removalDict.lookup("holeDetection"));
    }

    Switch splitMesh(false);
    if (removalDict.found("splitMesh"))
    {
        splitMesh = readBool(removalDict.lookup("splitMesh"));
    }

    Switch extrudeRepair(false);
    if (removalDict.found("extrudeRepair"))
    {
        extrudeRepair = readBool(removalDict.lookup("extrudeRepair"));
    }

    // Switch mesh quality removal off
    fileName meshDictFile = fileName("");

    if (extrudeRepair)
    {
        meshDictFile = "snappyHexMeshDict";

        IOobject autoDictHeader
        (
            meshDictFile,
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (!autoDictHeader.typeHeaderOk<IOdictionary>(true))
        {
            FatalError << "No snappyHexMeshDict required "
                       << "for extrusion repair"
                       << exit(FatalError);
        }

        multiRegion = true;
        holeDetection = false;
        meshQuality = false;
    }

    boolList flaggedCells(mesh.nCells(), true);
    int dummyTrackData = 0;

    if (multiRegion)
    {
        //- Areas to keep
        pointField keepPoints;

        if (snappyHex)
        {
            const dictionary& castDict = removalDict.subDict("castellatedMeshControls");
            if (castDict.lookup("locationInMesh"))
            {
                keepPoints = pointField(1, castDict.lookup("locationInMesh"));
            }
            else if (castDict.lookup("locationsInMesh"))
            {
                keepPoints = pointField(castDict.lookup("locationsInMesh"));
            }
        }
        else
        {
            keepPoints = pointField(removalDict.lookup("locationsInMesh"));
        }

        if (extrudeRepair)
        {
            // Read meshing dictionary
            IOdictionary extrudeDict
            (
                IOobject
                (
                    "extrudeProperties",
                    runTime.constant(),
                    runTime,
                    IOobject::MUST_READ
                )
            );

            const dictionary& extrudeModelDict = extrudeDict.subDict("fixedNormalCoeffs");
            const point axisPt(extrudeModelDict.lookup("axisPt"));
            const vector axisNormal(extrudeModelDict.lookup("axisNormal"));
            const scalar thickness(readScalar(extrudeModelDict.lookup("thickness")));

            forAll(keepPoints, i)
            {
                vector v1 = keepPoints[i]  - axisPt;
                point projPt =  keepPoints[i] - ((v1 & axisNormal) * axisNormal);
                keepPoints[i] = projPt + 0.5*thickness*axisNormal;
            }
        }

        // Cell label per point
        labelList cellLabels(keepPoints.size());
        {
            // Global calculation engine
            globalIndex globalCells(mesh.nCells());

            forAll(keepPoints, i)
            {
                const point& keepPoint = keepPoints[i];

                label localCellI = mesh.findCell(keepPoint);

                label globalCellI = -1;

                if (localCellI != -1)
                {
                    Pout<< "Found point " << keepPoint
                        << " in cell " << localCellI
                        << " on processor " << Pstream::myProcNo() << endl;
                    globalCellI = globalCells.toGlobal(localCellI);
                }

                reduceToMaster(globalCellI, maxOp<label>());

                if (globalCellI == -1)
                {
                    FatalErrorIn
                    (
                        "checkAnRemove::main"
                    )   << "Point " << keepPoint
                        << " is not inside the mesh." << nl
                        << "Bounding box of the mesh:" << mesh.bounds()
                        << exit(FatalError);
                }

                if (globalCells.isLocal(globalCellI))
                {
                    cellLabels[i] = localCellI;
                }
                else
                {
                    cellLabels[i] = -1;
                }
            }
        }


        // Initial information about connected mesh
        List<keepData> allCellInfo(mesh.nCells());

        // Initial information about connected mesh
        List<keepData> allFaceInfo(mesh.nFaces());

        // Labels of seed faces
        DynamicList<label> seedFaces(mesh.nFaces());
        //  data on seed faces
        DynamicList<keepData> seedFacesInfo(mesh.nFaces());

        if (holeDetection)
        {
            const polyBoundaryMesh& bMesh = mesh.boundaryMesh();
            const labelListList& faceEdges = mesh.faceEdges();
            boolList boundaryEdges(mesh.nEdges(), false);
            for (label patchI = 0; patchI < bMesh.size(); patchI++)
            {
                if (!isA<processorPolyPatch>(bMesh[patchI]))
                {
                    const polyPatch& patch = bMesh[patchI];
                    forAll(patch, patchFaceI)
                    {
                        label faceI = patchFaceI + patch.start();
                        forAll(faceEdges[faceI], edgeI)
                        {
                            label cEdge = faceEdges[faceI][edgeI];
                            boundaryEdges[cEdge] = true;
                        }
                    }
                }
            }

            if (splitMesh)
            {
                const faceZoneMesh& faceZones = mesh.faceZones();

                forAll(faceZones, zoneI)
                {
                    const labelList& mf = faceZones[zoneI];
                    forAll(mf, i)
                    {
                        label faceI = mf[i];
                        forAll(faceEdges[faceI], edgeI)
                        {
                            label cEdge = faceEdges[faceI][edgeI];
                            boundaryEdges[cEdge] = true;
                        }
                    }
                }
            }

            syncTools::syncEdgeList
            (
                mesh,
                boundaryEdges,
                orEqOp<bool>(),
                false
            );

            label nFaces = 1;
            if (removalDict.found("nFaces"))
            {
                nFaces = readLabel(removalDict.lookup("nFaces"));
            }
            labelList nInternalEdges(mesh.nFaces(), 0);
            forAll(mesh.neighbour(), ifI)
            {
                forAll(faceEdges[ifI], edgeI)
                {
                    label cEdge = faceEdges[ifI][edgeI];

                    if (!boundaryEdges[cEdge])
                    {
                        nInternalEdges[ifI]++;
                    }
                }
            }

            forAll(mesh.neighbour(), ifI)
            {
                if
                (
                    nInternalEdges[ifI] == 0
                    && !allFaceInfo[ifI].valid(dummyTrackData)
                 )
                {
                    seedFaces.append(ifI);
                    seedFacesInfo.append(2);
                    allFaceInfo[ifI] = seedFacesInfo[seedFacesInfo.size()-1];
                }
            }

            if (nFaces == 2)
            {
                labelList nSingleInternalEdge(mesh.nEdges(), 0);

                forAll(mesh.edges(), edgeI)
                {
                   if (!boundaryEdges[edgeI])
                   {
                       const labelList& ceFaces
                           = mesh.edgeFaces()[edgeI];

                       forAll(ceFaces, faceI)
                       {
                           if (nInternalEdges[ceFaces[faceI]] == 1)
                           {
                               nSingleInternalEdge[edgeI]++;
                           }
                       }
                   }
                }

                syncTools::syncEdgeList
                (
                    mesh,
                    nSingleInternalEdge,
                    plusEqOp<label>(),
                    labelMin
                );


               forAll(mesh.edges(), edgeI)
               {
                   if (!boundaryEdges[edgeI])
                   {
                       if (nSingleInternalEdge[edgeI] == 2)
                       {
                           const labelList& ceFaces
                               = mesh.edgeFaces()[edgeI];
                           forAll(ceFaces, faceI)
                           {
                               if (nInternalEdges[ceFaces[faceI]] == 1)
                               {
                                   if
                                   (
                                       !allFaceInfo[ceFaces[faceI]].valid
                                       (
                                           dummyTrackData
                                        )
                                    )
                                   {
                                       seedFaces.append(ceFaces[faceI]);
                                       seedFacesInfo.append(2);
                                       allFaceInfo[ceFaces[faceI]] =
                                           seedFacesInfo[seedFacesInfo.size()-1];
                                   }
                               }
                           }
                       }
                   }
               }
            }
        }

        if (splitMesh)
        {
            const faceZoneMesh& faceZones = mesh.faceZones();

            forAll(faceZones, zoneI)
            {
                const labelList& mf = faceZones[zoneI];
                forAll(mf, i)
                {
                    label faceI = mf[i];

                    if (!allFaceInfo[faceI].valid(dummyTrackData))
                    {
                        seedFaces.append(faceI);
                        seedFacesInfo.append(2);
                        allFaceInfo[faceI] = seedFacesInfo[seedFacesInfo.size()-1];
                    }
                }
            }
        }


        forAll(cellLabels, i)
        {
            label cellI = cellLabels[i];
            if (cellI != -1)
            {
                const cell& cFaces = mesh.cells()[cellI];
                forAll(cFaces, cFaceI)
                {
                    label faceI = cFaces[cFaceI];
                    if (!allFaceInfo[faceI].valid(dummyTrackData))
                    {
                        seedFaces.append(faceI);
                        seedFacesInfo.append(1);
                        allFaceInfo[faceI] =
                            seedFacesInfo[seedFacesInfo.size()-1];
                    }
                }
            }
        }

        // face-cell-face transport engine
        FaceCellWave<keepData> propagationCalc
        (
            mesh,
            allFaceInfo,
            allCellInfo
        );

        propagationCalc.setFaceInfo(seedFaces.shrink(), seedFacesInfo.shrink());
        seedFaces.clear();
        seedFacesInfo.clear();

        propagationCalc.iterate(mesh.globalData().nTotalFaces());

        flaggedCells = false;
        forAll(mesh.cells(), cellI)
        {
            if (allCellInfo[cellI].found() == 1)
            {
                flaggedCells[cellI] = true;
            }
        }
    }

    if (meshQuality)
    {
        // Checking mesh quality
        Info<< "Checking mesh ..." << endl;
        const dictionary& motionDict = removalDict.subDict("meshQualityControls");
        faceSet errorFaces(mesh, "errorFaces", mesh.nFaces()/100+100);

        motionSmoother::checkMesh(false, mesh, motionDict, errorFaces);
        Info<< "Writing set errorFaces "<<endl;
        errorFaces.write();
        forAllConstIter(labelHashSet, errorFaces, iter)
        {
            label own = mesh.faceOwner()[iter.key()];
            flaggedCells[own] = false;
            if (iter.key() < mesh.nInternalFaces())
            {
                label nei = mesh.faceNeighbour()[iter.key()];
                flaggedCells[nei] = false;
            }
        }
    }

    cellSet keptCells(mesh, "keptCells", mesh.nCells()/100);
    forAll(mesh.cells(), cellI)
    {
        if (flaggedCells[cellI] == true)
        {
            keptCells.insert(cellI);
        }
    }

    label nRemoved = mesh.nCells() - keptCells.size();

    if (returnReduce(nRemoved, sumOp<label>()) > 0)
    {
        Info<< "Writing kept cells ..." << endl;
        keptCells.write();

        // Add a new patch to the mesh to store exposed faces
        word patchName;
        if (removalDict.found("newPatch"))
        {
            patchName = word(removalDict.lookup("newPatch"));
        }
        else
        {
            patchName = "exposedFaces";
        }

        dictionary patchInfo;
        patchInfo.set("type", wallPolyPatch::typeName);
        label patchI = meshRefinement::addPatch
        (
            mesh,
            patchName,
            patchInfo
        );

        if (!overwrite)
        {
            runTime++;
        }

        // Create mesh subsetting engine
        Info<< "Subsetting mesh ..." << endl;
        fvMeshSubset subsetter(mesh);

        subsetter.setLargeCellSubset(keptCells, patchI, true);

        IOobjectList objects(mesh, runTime.timeName());

        // Perform a final balance
        // Read decomposePar dictionary
        autoPtr<mapDistributePolyMesh> map;
        if (Pstream::parRun())
        {
            Info<< "Re-balancing mesh ..." << endl;
            const scalar mergeTol =
                readScalar(removalDict.lookup("mergeTolerance"));
            const boundBox& meshBb = subsetter.subMesh().bounds();
            scalar mergeDist = mergeTol*mag(meshBb.max() - meshBb.min());

            IOdictionary decomposeDict
            (
                IOobject
                (
                    "decomposeParDict",
                    runTime.system(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            autoPtr<decompositionMethod> decomposerPtr =
                decompositionMethod::New
                (
                    decomposeDict
                );
            autoPtr<fvMeshDistribute> distributorPtr
            (
                new fvMeshDistribute
                    (
                        subsetter.subMesh(),
                        mergeDist
                    )
            );

            // Wanted distribution
            labelList distribution;
            distribution =
                decomposerPtr().decompose
                (
                    subsetter.subMesh().cellCentres()
                );
            // Do actual sending/receiving of mesh

            map = distributorPtr().distribute(distribution);
        }

        Info<< "Writing subsetted mesh to time " << runTime.timeName()
            << endl;

        if (overwrite)
        {
            subsetter.subMesh().setInstance(oldInstance);
        }

        subsetter.subMesh().write();
    }
    else
    {
        Info<< "No cells selected for deletion. No update to mesh. "<< endl;
    }
    Info<< "End\n" << endl;

    return 0;
}
