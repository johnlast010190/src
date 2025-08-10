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
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "foamMeshGenerator.H"
#include "searchableSurfaces/searchableSurfaces/searchableSurfaces.H"
#include "refinementSurfaces/refinementSurfaces.H"
#include "refinementFeatures/refinementFeatures.H"
#include "shellSurfaces/shellSurfaces.H"
#include "snappyHexMeshDriver/snapParameters/snapParameters.H"

#include "fvMeshDistribute/fvMeshDistribute.H"
#include "meshes/polyMesh/polyPatches/derived/wall/wallPolyPatch.H"
#include "sampledSetWriters/vtk/vtkSetWriter.H"
#include "motionSmoother/motionSmoother.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshes/primitiveMesh/primitivePatch/uindirectPrimitivePatch.H"
#include "meshes/polyMesh/polyPatches/constraint/symmetry/symmetryPolyPatch.H"
#include "MeshedSurface/MeshedSurface.H"
#include "db/IOstreams/IOstreams/IOmanip.H"
#include "decompositionModel.H"
#include "addHexMeshLayers/addHexMeshLayer.H"
#include "addHexInternalLayers/addHexInternalLayers.H"
#include "autoHexDualiser/autoHexDualiser.H"
#include "addHexMeshLayers/addMultiLayers.H"
#include "autoExtrude/autoExtrude.H"
#include "layerManipulate/layerManipulate.H"
#include "autoSplitCells/autoSplitCells.H"
#include "autoOptimize/autoOptimize.H"
#include "errorCellRemoval/errorCellRemoval.H"
#include "meshControl/meshControl.H"
#include "autoBlockMesh/autoBlockMesh.H"
#include "autoLayerCellsMerge/autoLayerCellsMerge.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"
#include "blockMesh/blockMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

dictionary Foam::foamMeshGenerator::finalOptimizeSetup
(
   const dictionary& meshDict
)
{
    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");
    dictionary meshOptimDict
    (
        meshDict.found("meshOptimization") ?
        meshDict.subDict("meshOptimization") :
        dictionary()
    );

    word optimType =
        meshOptimDict.lookupOrDefault<word>
        (
            "type",
            "cfMeshOptimize"
        );
    if (optimType == "foamOptimize")
    {
        if (meshOptimDict.found("foamOptimizeCoeffs"))
        {
            dictionary& coeffsDict =
                meshOptimDict.subDict("foamOptimizeCoeffs");
            if (!coeffsDict.found("meshQualityControls"))
            {
                coeffsDict.add
                (
                    "meshQualityControls",
                    motionDict,
                    true
                );
            }
        }
    }
    else if (optimType == "cfMeshOptimize")
    {
        if (!meshOptimDict.found("cfMeshOptimizeCoeffs"))
        {
            meshOptimDict.add
            (
                "cfMeshOptimizeCoeffs",
                dictionary(),
                true
             );
        }
        dictionary& coeffsDict =
            meshOptimDict.subDict("cfMeshOptimizeCoeffs");

        const Switch finalRelaxCheck
        (
            meshDict.lookupOrDefault("finalRelaxedCheck", true)
        );

        coeffsDict.add
        (
            "relaxedCheck",
            finalRelaxCheck,
            true
        );
    }

    return meshOptimDict;
}

//Perform optimization if errors still exist in the mesh
void Foam::foamMeshGenerator::cfOptimize
(
    fvMesh& mesh,
    const dictionary& meshOptimDict
)
{
    word optimType =
        meshOptimDict.lookupOrDefault<word>("type", "cfMeshOptimize");

    if (optimType != "none" && optimType == "foamOptimize")
    {
        label noFailedChecks = 0;
        if (mesh.checkCellVolumes(true)) noFailedChecks++;
        if (mesh.checkFaceOrthogonality(true)) noFailedChecks++;
        if (mesh.checkFacePyramids(true)) noFailedChecks++;

        if (noFailedChecks > 0)
        {
            //Add cfMesh optimization with default settings
            dictionary dummyOptimDict;
            dummyOptimDict.add("type","cfMeshOptimize",true);
            dummyOptimDict.add
                ("cfMeshOptimizeCoeffs",dictionary(), true);

            //optimize mesh
            {
                autoPtr<autoOptimize> optimMeshPtr
                    = autoOptimize::New(mesh, dummyOptimDict);
                optimMeshPtr->optimize();
            }
        }
    }
}


void Foam::foamMeshGenerator::extractTriSurface
(
    const polyMesh& mesh,
    const labelHashSet& includePatches,
    const fileName& outFileName
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // Collect sizes. Hash on names to handle local-only patches (e.g.
    //  processor patches)
    HashTable<label> patchSize(1024);
    label nFaces = 0;
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        patchSize.insert(pp.name(), pp.size());
        nFaces += pp.size();
    }
    Pstream::mapCombineGather(patchSize, plusEqOp<label>());


    // Allocate zone/patch for all patches
    HashTable<label> compactZoneID(1024);
    forAllConstIter(HashTable<label>, patchSize, iter)
    {
        label sz = compactZoneID.size();
        compactZoneID.insert(iter.key(), sz);
    }
    Pstream::mapCombineScatter(compactZoneID);


    // Rework HashTable into labelList just for speed of conversion
    labelList patchToCompactZone(bMesh.size(), -1);
    forAllConstIter(HashTable<label>, compactZoneID, iter)
    {
        label patchi = bMesh.findPatchID(iter.key());
        if (patchi != -1)
        {
            patchToCompactZone[patchi] = iter();
        }
    }

    // Collect faces on zones
    DynamicList<label> faceLabels(nFaces);
    DynamicList<label> compactZones(nFaces);
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        const polyPatch& pp = bMesh[iter.key()];
        forAll(pp, i)
        {
            faceLabels.append(pp.start()+i);
            compactZones.append(patchToCompactZone[pp.index()]);
        }
    }

    // Addressing engine for all faces
    uindirectPrimitivePatch allBoundary
    (
        UIndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );


    // Find correspondence to master points
    labelList pointToGlobal;
    labelList uniqueMeshPoints;
    autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
    (
        allBoundary.meshPoints(),
        pointToGlobal,
        uniqueMeshPoints
    );

    // Gather all unique points on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()] = pointField
    (
        mesh.points(),
        uniqueMeshPoints
    );
    Pstream::gatherList(gatheredPoints);

    // Gather all faces
    List<faceList> gatheredFaces(Pstream::nProcs());
    gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
    forAll(gatheredFaces[Pstream::myProcNo()], i)
    {
        inplaceRenumber(pointToGlobal, gatheredFaces[Pstream::myProcNo()][i]);
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


        // Zones
        surfZoneIdentifierList surfZones(compactZoneID.size());
        forAllConstIter(HashTable<label>, compactZoneID, iter)
        {
            surfZones[iter()] = surfZoneIdentifier(iter.key(), iter());
            Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
                << endl;
        }

        UnsortedMeshedSurface<face> unsortedFace
        (
            xferMove(allPoints),
            xferMove(allFaces),
            xferMove(allZones),
            xferMove(surfZones)
        );


        MeshedSurface<face> sortedFace(unsortedFace);

        fileName globalCasePath
        (
            runTime_.processorCase()
          ? runTime_.path()/".."/outFileName
          : runTime_.path()/outFileName
        );
        globalCasePath.clean();

        Info<< "Writing merged surface to " << globalCasePath << endl;

        sortedFace.write(globalCasePath);
    }
}


// Check writing tolerance before doing any serious work
Foam::scalar Foam::foamMeshGenerator::getMergeDistance
(
    const polyMesh& mesh,
    const boundBox& meshBb,
    const scalar mergeTol
)
{
    scalar mergeDist = mergeTol * meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    // check writing tolerance
    if (mesh.time().writeFormat() == IOstream::ASCII)
    {
#if defined(WIN64) || defined(WIN32)
        scalar writeTol = std::pow
#else
        const scalar writeTol = std::pow
#endif
        (
            scalar(10.0),
            -scalar(IOstream::defaultPrecision())
        );

#if defined(WIN32) || defined(WIN64)
        writeTol -= SMALL;
#endif

        if (mergeTol < writeTol)
        {
            FatalErrorInFunction
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << nl
                << "Your merging tolerance (" << mergeTol
                << ") is finer than this." << nl
                << "Change to binary writeFormat, "
                << "or increase the writePrecision" << endl
                << "or adjust the merge tolerance (mergeTol)."
                << exit(FatalError);
        }
    }

    return mergeDist;
}


// Write mesh and additional information
void Foam::foamMeshGenerator::writeMesh
(
    const string& msg,
    meshRefinement& meshRefiner,

    const meshRefinement::debugType debugLevel,
    const meshRefinement::writeType writeLevel,

    bool writeSurf,
    bool removeZeroSizedPatches
)
{
    fvMesh& mesh = meshRefiner.mesh();

    meshRefiner.setInstance(mesh.facesInstance());

    if (removeZeroSizedPatches)
    {
        meshRefiner.removeZeroSizedPatches(mesh);
    }

    scalar startTime = mesh.time().elapsedCpuTime();

    meshRefiner.printMeshInfo(debugLevel, msg);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    meshRefiner.write
    (
        debugLevel,
        meshRefinement::writeType(writeLevel | meshRefinement::WRITEMESH),
        "",
        writeSurf
    );
    Info<< "Wrote mesh in = "
        << mesh.time().elapsedCpuTime() - startTime << " s." << endl;
}


void Foam::foamMeshGenerator::blockMeshTriangulate
(
    const polyMesh& mesh,
    dictionary& geometryDict,
    dictionary& refineDict,
    const bool baseCheck
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    if (baseCheck)
    {
        bool validMesh = true;
        forAll(bMesh, patchi)
        {
            const polyPatch& pp = bMesh[patchi];
            if (!isA<processorPolyPatch>(pp))
            {
                label startFaceI = pp.start();
                forAll(pp, i)
                {
                    face f = mesh.faces()[startFaceI + i];
                    if (f.size() != 4)
                    {
                        validMesh = false;
                        break;
                    }
                }
            }
        }
        reduceToMaster(validMesh, andOp<bool>());

        if (!validMesh)
        {
            FatalErrorInFunction
                << "Not a valid starting blockMesh for extrude or dual method"
                << exit(FatalError);
        }
    }

    label nTris = 0;
    forAll(bMesh, patchI)
    {
        const polyPatch& pp = bMesh[patchI];

        if (!isA<processorPolyPatch>(pp))
        {
            word patchName = pp.name();
            word stlName = patchName + ".stl";

            nTris += 2*pp.size();
            List<label> faceLabels(pp.size());
            forAll(pp, i)
            {
                faceLabels[i] = pp.start() + i;
            }

            uindirectPrimitivePatch allBoundary
            (
                UIndirectList<face>(mesh.faces(), faceLabels),
                mesh.points()
            );

            // Find correspondence to master points
            labelList pointToGlobal;
            labelList uniqueMeshPoints;
            autoPtr<globalIndex> globalNumbers =
                mesh.globalData().mergePoints
                (
                    allBoundary.meshPoints(),
                    pointToGlobal,
                    uniqueMeshPoints
                );

            // Gather all unique points on master
            List<pointField> gatheredPoints(Pstream::nProcs());
            gatheredPoints[Pstream::myProcNo()] = pointField
            (
                mesh.points(),
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

                UnsortedMeshedSurface<face> unsortedFace
                (
                    xferMove(allPoints),
                    xferMove(allFaces)
                );
                MeshedSurface<face> sortedFace(unsortedFace);

                fileName globalSTLPath =
                    runTime_.rootPath()/runTime_.globalCaseName()
                    /"constant"/"triSurface"/stlName;
                sortedFace.write(globalSTLPath);
            }

            geometryDict.add(stlName, dictionary(), false);

            char defaultType[] = "triSurfaceMesh";
            geometryDict.subDict(stlName).add
            (
                word("type"),
                defaultType,
                true
            );

            geometryDict.subDict(stlName).add
            (
                word("name"),
                patchName,
                true
            );

            geometryDict.subDict(stlName).add
            (
                word("appendRegionName"),
                false,
                true
            );


            dictionary& refSurf =
                refineDict.subDict("refinementSurfaces");

            refSurf.add(patchName, dictionary(), false);
            char defaultLevel[] = "(0 0)";

            refSurf.subDict(patchName).add
            (
                word("level"),
                defaultLevel,
                true
             );
        }
    }

    Info<<"Triangulated blockMesh patches. Number of facets : "
        << returnReduce(nTris, sumOp<label>()) <<endl;
}


void Foam::foamMeshGenerator::addVDBpatches
(
    const dictionary dict,
    fvMesh& mesh
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    fvBoundaryMesh& fvPatches =
        const_cast<fvBoundaryMesh&>(mesh.boundary());

    polyPatches.setSize(6);
    fvPatches.setSize(6);

    const wordList vdbPatches({"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"});

    for (label patchi = 0; patchi < 6; patchi++)
    {
        dictionary patchInfo;

        word patchName = vdbPatches[patchi];
        word patchType = wallPolyPatch::typeName;

        if (dict.found("patches"))
        {
            if (dict.subDict("patches").found(patchName))
            {
                const dictionary& patchDict =
                    dict.subDict("patches").subDict(patchName);

                patchName = patchDict.lookupOrDefault<word>("name", patchName);
                patchType = patchDict.lookupOrDefault<word>("type", patchType);
            }
        }

        patchInfo.add("type", patchType);
        patchInfo.add("nFaces", 0);
        patchInfo.add("startFace", 0);

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

        fvPatches.set
        (
            patchi,
            fvPatch::New
            (
               polyPatches[patchi],  // point to newly added polyPatch
               mesh.boundary()
            )
        );
    }
} //addVDBpatches


Foam::dictionary Foam::foamMeshGenerator::calcReSnapNonBoundaryZoneDict
(
    const layerParameters& layerParams,
    const refinementSurfaces& surfaces,
    const dictionary& geometryDict
)
{
    dictionary reSnapZoneGeometryDict;
    if (layerParams.dualReSnapZones())
    {
        const labelList nonBoundaryNamedSurfaces =
            surfaceZonesInfo::getNonBoundaryNamedSurfaces(surfaces.surfZones());

        const labelList& surfaceGeometry = surfaces.surfaces();
        const searchableSurfaces& geometry = surfaces.geometry();

        forAll(nonBoundaryNamedSurfaces, surfi)
        {
            label geomi = surfaceGeometry[nonBoundaryNamedSurfaces[surfi]];
            word fName = geometry.fileNames()[geomi];
            reSnapZoneGeometryDict.add(fName,geometryDict.subDict(fName));
        }
    }

    return reSnapZoneGeometryDict;
}


Foam::dictionary Foam::foamMeshGenerator::calcReSnapGeometryDict
(
    const layerParameters& layerParams,
    const refinementSurfaces& surfaces,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const dictionary& geometryDict
)
{
    const labelList& patchToNLayers = layerParams.numLayers();
    const List<Switch>& reSnap = layerParams.reSnap();

    DynamicList<label> grownUpSnapPatches(patchToNLayers.size());
    forAll(patchToNLayers, patchI)
    {
        if
        (
            reSnap[patchI] && patchToNLayers[patchI] == -1
        )
        {
            grownUpSnapPatches.append(patchI);
        }
    }

    labelHashSet grownUpSnapSet(grownUpSnapPatches);
    boolList markGlobal(globalToMasterPatch.size(), false);
    forAll(globalToMasterPatch, globalI)
    {
        label masterPatchI = globalToMasterPatch[globalI];
        label slavePatchI = globalToSlavePatch[globalI];

        if (grownUpSnapSet.found(masterPatchI))
        {
            markGlobal[globalI] = true;
        }
        if (grownUpSnapSet.found(slavePatchI))
        {
            markGlobal[globalI] = true;
        }
    }

    const labelList& surfaceGeometry = surfaces.surfaces();
    const searchableSurfaces& geometry = surfaces.geometry();

    dictionary reSnapGeometryDict;
    forAll(surfaceGeometry, surfi)
    {
        label geomi = surfaceGeometry[surfi];
        const wordList& regNames = geometry.regionNames()[geomi];

        // 'Normal' surface
        forAll(regNames, i)
        {
            label globalI = surfaces.globalRegion(surfi, i);

            if (markGlobal[globalI])
            {
                word fName = geometry.fileNames()[geomi];
                reSnapGeometryDict.add(fName,geometryDict.subDict(fName));
                break;
            }
        }
    }

    return reSnapGeometryDict;
}


void Foam::foamMeshGenerator::autoCreateDecomposition
(
    const label nProcs,
    dictionary& decomposeDict
)
{
    Info<<"Calculating hierarchical coefficients automatically"<<endl;

    IOdictionary* oldDict =
        runTime_.lookupObjectRefPtr<IOdictionary>("decomposeParDict");
    if (oldDict)
    {
        // Deletes if owned
        oldDict->checkOut();
        delete oldDict;
    }
    IOdictionary* dictPtr = new IOdictionary
    (
        IOobject
        (
            "decomposeParDict",
            runTime_.system(),
            runTime_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
         )
    );

    runTime_.store(dictPtr);
    decomposeDict = *dictPtr;

    decomposeDict.add("numberOfSubdomains", Pstream::nProcs());
    decomposeDict.add("method", "hierarchical");
    decomposeDict.add("distributed", "false");
    decomposeDict.add(word("hierarchicalCoeffs"), dictionary());
    decomposeDict.subDict("hierarchicalCoeffs").add
    (
        word("order"),
        word("yxz")
    );

    decomposeDict.subDict("hierarchicalCoeffs").add
    (
        "delta",
        0.001
    );

    //Calculate hierarchical factors
    label np_3 = label(pow(nProcs, 1.0/3.0));

    label firstFactor = 0;
    for (label i = np_3; i <= nProcs; i++)
    {
        if (nProcs % i == 0)
        {
            firstFactor = i;
            break;
        }
    }
    label remainder = nProcs/firstFactor;
    label np_2 = label(pow(remainder, 0.5));

    label secondFactor = 0;
    for (label i = np_2; i <= nProcs; i++)
    {
        if (remainder % i == 0)
        {
            secondFactor = i;
            break;
        }
    }

    label thirdFactor = nProcs / firstFactor / secondFactor;

    labelVector n(firstFactor, secondFactor, thirdFactor);

    Info<<"Using hierarchical coefficients "<<n<<endl;

    decomposeDict.subDict("hierarchicalCoeffs").add("n",n);
}


void Foam::foamMeshGenerator::addProtectedCellsExtrude
(
    const fvMesh& mesh,
    volScalarField& protectedCells
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const labelList& owners = mesh.faceOwner();
    const labelList& neighbours = mesh.faceNeighbour();

    const volScalarField& layerCells =
        mesh.lookupObject<volScalarField>("layerStacks");

    forAll(layerCells, celli)
    {
        if (layerCells[celli] != -1)
        {
            protectedCells[celli] = scalar(1);
        }
    }

    labelList neiLayerCells(mesh.nFaces()-mesh.nInternalFaces());
    for
    (
        label facei = mesh.nInternalFaces();
        facei < mesh.nFaces();
        facei++
    )
    {
        neiLayerCells[facei-mesh.nInternalFaces()] =
            layerCells[owners[facei]];
    }
    syncTools::swapBoundaryFaceList(mesh, neiLayerCells);
    faceSet interFaces(mesh, "interFaces", mesh.nFaces()/100+100);
    forAll(mesh.faces(), facei)
    {
        label patchi = patches.whichPatch(facei);

        if (patchi == -1 || patches[patchi].coupled())
        {
            label own = owners[facei];
            label oLayerCell = layerCells[own];

            label nLayerCell = -1;

            if (patchi == -1)
            {
                label nei = neighbours[facei];
                nLayerCell = layerCells[nei];
            }
            else
            {
                label bfacei = facei-mesh.nInternalFaces();
                nLayerCell = neiLayerCells[bfacei];
            }
            if
            (
                (nLayerCell > -1 && oLayerCell < 0)
                || (oLayerCell > -1 && nLayerCell < 0)
            )
            {
                interFaces.insert(facei);
            }
        }
    }
    pointSet interPoints(mesh, "interPts",  mesh.nPoints()/100+100);
    forAllConstIter(pointSet, interFaces, iter)
    {
        label facei = iter.key();
        interPoints.insert(mesh.faces()[facei]);
    }

    forAllConstIter(pointSet, interPoints, iter)
    {
        label pointi = iter.key();
        const labelList& pCells = mesh.pointCells()[pointi];
        forAll(pCells, pCI)
        {
            protectedCells[pCells[pCI]] = scalar(1);
        }
    }
}


bool Foam::foamMeshGenerator::addLayers
(
    const dictionary& layerDict,
    const dictionary& motionDict,
    const meshControl& controller,
    const labelList& globalToMasterPatch,
    const labelList& globalToSlavePatch,
    const bool& overwrite,
    meshRefinement& meshRefiner,
    autoPtr<searchableSurfaces>& allGeometryPtr,
    hexReport& stats,
    decompositionMethod& decomposer,
    fvMeshDistribute& distributor,
    const label multiLayerIter
)
{
    const fvMesh& mesh = meshRefiner.mesh();
    snappyLayerDriver layerDriver
    (
        meshRefiner,
        allGeometryPtr,
        stats,
        controller,
        globalToMasterPatch,
        globalToSlavePatch
    );

    if (!overwrite)
    {
        const_cast<Time&>(mesh.time())++;
    }

    dictionary localMotionDict
    (
        layerDict.found("meshQualityControls") ?
        layerDict.subDict("meshQualityControls") :
        motionDict
    );

    return layerDriver.doLayers
    (
        layerDict,
        localMotionDict,
        decomposer,
        distributor,
        multiLayerIter
    );
}


void Foam::foamMeshGenerator::runBlockMesh
(
    const dictionary& meshDict,
    const word& regionName,
    const word& regionDir,
    autoPtr<fvMesh>& meshPtr,
    decompositionMethod& decomposer
)
{
   word defaultFacesName = "defaultFaces";
   word defaultFacesType = emptyPolyPatch::typeName;

   wordList masterPatchNames;
   PtrList<dictionary> masterBoundaryDicts;
   bool oldParRun = Pstream::parRun();

   word blockDictName("blockMeshDict");
   fileName dictPath;
   if
   (
      exists
      (
         Pstream::parRun() ?
         runTime_.path()/".."/runTime_.constant()
         /regionDir/polyMesh::meshSubDir/blockDictName
         :  runTime_.path()/runTime_.constant()
         /regionDir/polyMesh::meshSubDir/blockDictName
      )
   )
   {

      dictPath =
      (
          Pstream::parRun() ?
          runTime_.path()/".."/runTime_.constant()
          /regionDir/polyMesh::meshSubDir/blockDictName
          : runTime_.path()/runTime_.constant()
          /regionDir/polyMesh::meshSubDir/blockDictName
       );
   }
   else
   {
      dictPath =
      (
         Pstream::parRun() ?
         runTime_.path()/".."/runTime_.system()
         /regionDir/blockDictName
         : runTime_.path()/runTime_.system()
         /regionDir/blockDictName
      );
   }

   Pstream::parRun() = false;
   if (Pstream::master())
   {
       IOobject blockDictIO
       (
           dictPath,
           runTime_,
           IOobject::MUST_READ,
           IOobject::NO_WRITE,
           false,
           false
       );
       IOdictionary blockDict(blockDictIO);
       blockMesh blocks
       (
           blockDict,
           regionName,
           false //don't read pre-existing boundary files
       );

       masterPatchNames = blocks.patchNames();
       masterBoundaryDicts = blocks.patchDicts();

       meshPtr.reset
       (
           new fvMesh
           (
               Foam::IOobject
               (
                   regionName,
                   runTime_.timeName(),
                   runTime_,
                   Foam::IOobject::MUST_READ
               ),
               xferCopy(blocks.points()),
               blocks.cells(),
               blocks.patches(),
               blocks.patchNames(),
               blocks.patchDicts(),
               defaultFacesName,
               defaultFacesType,
               false
           )
       );
   }
   Pstream::parRun() = oldParRun;
   Pstream::scatter(masterPatchNames);
   Pstream::scatter(masterBoundaryDicts);
   Pstream::parRun() = false;
   if (!Pstream::master())
   {
      forAll(masterBoundaryDicts, patchi)
      {
         dictionary& pDict = masterBoundaryDicts[patchi];
         pDict.add("nFaces", label(0),true);
         pDict.add("startFace", label(0),true);
      }
      meshPtr.reset
      (
         new fvMesh
         (
             Foam::IOobject
             (
                 regionName,
                 runTime_.timeName(),
                 runTime_,
                 Foam::IOobject::MUST_READ
             ),
             xferCopy(pointField(0)),
             xferCopy(faceList(0)),
             xferCopy(cellList(0)),
             false,
             false
         )
      );
      polyBoundaryMesh& polyPatches =
         const_cast<polyBoundaryMesh&>(meshPtr().boundaryMesh());
      fvBoundaryMesh& fvPatches =
         const_cast<fvBoundaryMesh&>(meshPtr().boundary());

      forAll(masterPatchNames, i)
      {
          label patchi = polyPatches.size();
          // Add polyPatch at the end
          polyPatches.setSize(patchi+1);
          polyPatches.set
          (
              patchi,
              polyPatch::New
              (
                  masterPatchNames[i],
                  masterBoundaryDicts[i],
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
                  meshPtr().boundary()
              )
          );
      }
   }
   Pstream::parRun() = oldParRun;

   if (Pstream::parRun())
   {
       const scalar mergeDist = getMergeDistance
       (
           meshPtr(),
           meshPtr().bounds(),
           readScalar(meshDict.lookup("mergeTolerance"))
       );

       // Mesh distribution engine (uses tolerance to reconstruct meshes)
       fvMeshDistribute distributor(meshPtr(), mergeDist);

       labelList distribution =
          decomposer.decompose(meshPtr(), meshPtr().cellCentres());

       // Do actual sending/receiving of mesh
       distributor.distribute(distribution);
   }

   return;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamMeshGenerator::foamMeshGenerator(Time& time)
:
    runTime_(time)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::fvMesh> Foam::foamMeshGenerator::generateMesh
(
   const word& regionName,
   const word& regionDir,
   const fileName& dictPath,
   const fileName& decompDictFile,
   const bool& baseCheck,
   const bool& overwrite,
   const bool& writeIntermediate,
   const bool& writeConstant,
   const scalar& writeTime,
   const bool& checkGeometry,
   bool writeDict,
   bool meshWrite
)
{
    word dictName("foamHexMeshDict");

    IOobject meshHeader
    (
        dictName,
        runTime_.system(),
        regionDir,
        runTime_,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    if (!meshHeader.typeHeaderOk<IOdictionary>(true))
    {
        dictName = "snappyHexMeshDict";
    }

    IOobject dictIO
    (
        dictName,
        runTime_.system(),
        regionDir,
        runTime_,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        true,
        true
    );

    if (dictPath.size())
    {
        dictIO = IOobject
        (
            dictPath,
            runTime_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            true,
            true
        );
    }

    if (!dictIO.typeHeaderOk<IOdictionary>(true))
    {
        if (Pstream::master())
        {
            fileName failedPath = dictIO.rootPath()
                + "/"+dictIO.time().globalCaseName()
                + "/"+dictIO.instance()
                + "/"+dictIO.local();

            FatalErrorInFunction << "Cannot find meshing dictionary "
                << "foamHexMeshDict in folder : " << failedPath
                << exit(FatalError);
        }
    }

    IOdictionary meshDict(dictIO);

    //Setup overall control of mesh generator
    const meshControl controller(meshDict);

    const bool useVDB = controller.vdb() && controller.refine();
    bool autoDecompose = false;

    // Read decomposePar dictionary
    dictionary decomposeDict;
    label numSubDomains = Pstream::nProcs();
    if (Pstream::parRun() || useVDB)
    {
        // A demand-driven decompositionMethod can have issues finding
        // an alternative decomposeParDict location.

       bool checkHeaderOK = false;
       if (decompDictFile.size())
       {
           IOobject decompHeader
           (
               decompDictFile,
               runTime_,
               IOobject::MUST_READ_IF_MODIFIED,
               IOobject::NO_WRITE
           );

           if (decompHeader.typeHeaderOk<IOdictionary>(true))
           {
              checkHeaderOK = true;
           }
       }
       else
       {
           IOobject decompHeader
           (
               "decomposeParDict",
               runTime_.system(),
               regionDir,
               runTime_,
               IOobject::MUST_READ_IF_MODIFIED,
               IOobject::NO_WRITE
           );

           if (decompHeader.typeHeaderOk<IOdictionary>(true))
           {
               checkHeaderOK = true;
           }
       }

       if (checkHeaderOK)
       {
           IOdictionary* dictPtr = new IOdictionary
           (
               decompositionModel::selectIO
               (
                   IOobject
                   (
                       "decomposeParDict",
                       runTime_.system(),
                       regionDir,
                       runTime_,
                       IOobject::MUST_READ,
                       IOobject::NO_WRITE,
                       true,
                       true
                   ),
                   decompDictFile
               )
           );

           numSubDomains =
              readLabel(dictPtr->lookup("numberOfSubdomains"));

           if (numSubDomains != Pstream::nProcs())
           {
               if (useVDB)
               {
                   if (controller.snap() || controller.layers())
                   {
                      FatalErrorInFunction
                          << "The number of sub-domains "
                          << numSubDomains
                          << " differs from run settings "
                          << Pstream::nProcs()
                          << ".\nPlease switch off snap and layer in "
                          << "foamHexMeshDict."
                          << "\nAfter castellation the mesh will be "
                          << "redistributed to " << numSubDomains
                          << " domains and you can resume snap and layer"
                          << " addition using " << numSubDomains
                          << " processors"
                          << exit(FatalError);
                   }
               }
               else
               {
                   autoDecompose = true;
                   if (controller.block() == meshControl::IMPORT)
                   {
                      FatalErrorInFunction
                          << "The number of sub-domains " << numSubDomains
                          << " differs from run settings "
                          << Pstream::nProcs()
                          << " reset decomposeParDict accordingly"
                          << exit(FatalError);
                   }
               }
           }

           const word initialDecomp(dictPtr->lookup("method"));
           if
           (
               initialDecomp != "hierarchical"
               && initialDecomp != "simple" && !useVDB
           )
           {
               autoDecompose = true;
               WarningInFunction
                   << "Using decomposition method " << initialDecomp
                   << " which might result in less stable meshing. "
                   << " Automatically resetting to hierarchical " << endl;
           }

           if (!autoDecompose)
           {
               // Store it on the object registry, but to be found
               // it must also have the expected "decomposeParDict" name.
               dictPtr->rename("decomposeParDict");
               runTime_.store(dictPtr);
               decomposeDict = *dictPtr;
           }
       }
       else
       {
          autoDecompose = true;
       }

        if (autoDecompose)
        {
            autoCreateDecomposition
            (
                Pstream::nProcs(),
                decomposeDict
            );
        }
    }
    else
    {
        decomposeDict.add("method", "none");
        decomposeDict.add("numberOfSubdomains", 1);
    }

    if (controller.mode() == meshControl::DRYRUN)
    {
        Info<< "Operating in dry-run mode to detect set-up errors"
            << nl << endl;
    }

    if (!writeDict && meshDict.found("writeDict"))
    {
        writeDict = readBool(meshDict.lookup("writeDict"));
    }

    const Switch keepPatches(meshDict.lookupOrDefault("keepPatches", false));

    // format to be used for writing lines
    const word setFormat
    (
        meshDict.lookupOrDefault
        (
            "setFormat",
            vtkSetWriter<scalar>::typeName
        )
    );
    const autoPtr<writer<scalar>> setFormatter
    (
        writer<scalar>::New(setFormat)
    );

    if (writeDict)
    {
        Info<< "Writing out mesh dictionary " << endl;
        Info<< meshDict << endl;
    }

    // all surface geometry
    dictionary& geometryDict = meshDict.subDict("geometry");

    // refinement parameters
    dictionary& refineDict = meshDict.subDict("castellatedMeshControls");

    // mesh motion and mesh quality parameters
    const dictionary& motionDict = meshDict.subDict("meshQualityControls");

    // snap-to-surface parameters
    const dictionary& snapDict = meshDict.subDict("snapControls");

    // Debug
    // ~~~~~
    // Set debug level
    meshRefinement::debugType debugLevel = meshRefinement::debugType
    (
        meshDict.lookupOrDefault<label>
        (
            "debug",
            0
        )
    );
    {
        wordList flags;
        if (meshDict.readIfPresent("debugFlags", flags))
        {
            debugLevel = meshRefinement::debugType
            (
                meshRefinement::readFlags
                (
                    meshRefinement::IOdebugTypeNames,
                    flags
                )
            );
        }
    }
    if (debugLevel > 0)
    {
        meshRefinement::debug   = debugLevel;
        snappyRefineDriver::debug = debugLevel;
        snappySnapDriver::debug   = debugLevel;
        snappyLayerDriver::debug  = debugLevel;
    }

    // Set file writing level
    {
        wordList flags;
        if (meshDict.readIfPresent("writeFlags", flags))
        {
            meshRefinement::writeLevel
            (
                meshRefinement::writeType
                (
                    meshRefinement::readFlags
                    (
                        meshRefinement::IOwriteTypeNames,
                        flags
                    )
                )
            );
        }
    }

#if !defined( WIN32 ) && !defined( WIN64 )
    // for the impatient who want to see some output files:
    profiling::writeNow();
#endif

    // Parallel
    // ~~~~~~~~

    // Decomposition
    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New
        (
            decomposeDict
        )
    );
    decompositionMethod& decomposer = decomposerPtr();

    word finalDecomp =
        meshDict.lookupOrDefault<word>("finalDecomposition", "ptscotch");

    if (Pstream::parRun())
    {
        if (!decomposer.parallelAware())
        {
            FatalErrorInFunction
                << "You have selected decomposition method "
                << decomposer.typeName
                << " which is not parallel aware." << endl
                << "Please select one that is (hierarchical, ptscotch)"
                << exit(FatalError);
        }
        //Check final decomposition method valid
        decompositionMethod::dictionaryConstructorTable::iterator cstrIter =
            decompositionMethod::dictionaryConstructorTable_().
            find(finalDecomp);
        if
        (
           cstrIter == decompositionMethod::dictionaryConstructorTable_().end()
        )
        {
            cstrIter =decompositionMethod::dictionaryConstructorTable_().
                find("ptscotch");
            if
            (
               cstrIter
               != decompositionMethod::dictionaryConstructorTable_().end()
            )
            {
                WarningInFunction
                    << "Unknown finalDecomposition keyword method "
                    << finalDecomp
                    << "Reseting to ptscotch method " << nl
                    << endl;
                finalDecomp = word("ptscotch");
            }
            else
            {
               FatalErrorInFunction
                    << "Unknown finalDecomposition keyword "
                    << finalDecomp << nl << nl
                    << exit(FatalError);
            }
        }
    }

    //Read in mesh or construct from blockMeshDict
    autoPtr<fvMesh> meshPtr;
    if (controller.block() != meshControl::AUTO && !useVDB)
    {
        if (controller.block() != meshControl::IMPORT)
        {
            runBlockMesh
            (
                meshDict,
                regionName,
                regionDir,
                meshPtr,
                decomposer
            );
        }
        else
        {
            meshPtr.reset
            (
                new fvMesh
                (
                    Foam::IOobject
                    (
                        regionName,
                        runTime_.timeName(),
                        runTime_,
                        Foam::IOobject::MUST_READ
                    )
                )
            );
        }

        if
        (
            controller.algorithm() == meshControl::DUAL
         || controller.algorithm() == meshControl::EXTRUDE
        )
        {
            // Reset any symmetry patches to walls as these might cause mesh
            // motion issues for the extrude method
            polyBoundaryMesh& polyPatches =
                const_cast<polyBoundaryMesh&>(meshPtr().boundaryMesh());
            fvBoundaryMesh& fvPatches =
                const_cast<fvBoundaryMesh&>(meshPtr().boundary());
            forAll(polyPatches, patchi)
            {
                const polyPatch& pp = polyPatches[patchi];
                if (isA<symmetryPolyPatch>(pp))
                {
                    dictionary patchInfo;
                    word patchType = wallPolyPatch::typeName;
                    patchInfo.add("type", patchType);
                    patchInfo.add("nFaces", pp.size());
                    patchInfo.add("startFace", pp.start());

                    polyPatches.set
                    (
                        patchi,
                        polyPatch::New
                        (
                            pp.name(),
                            patchInfo,
                            patchi,
                            polyPatches
                        )
                    );

                    fvPatches.set
                    (
                        patchi,
                        fvPatch::New
                        (
                            polyPatches[patchi],  // point to newly added polyPatch
                            meshPtr().boundary()
                        )
                    );
                }
            }

            //create triSurface folder if not present
            if (Pstream::parRun())
            {
                fileName gPath = runTime_.rootPath()/runTime_.globalCaseName();
                if
                (
                    Pstream::master()
                    && !isDir(gPath/"constant"/"triSurface")
                )
                {
                    mkDir(gPath/"constant"/"triSurface");
                }
            }
            else
            {
                if (!isDir(runTime_.path()/"constant"/"triSurface"))
                {
                    mkDir(runTime_.path()/"constant"/"triSurface");
                }
            }
            blockMeshTriangulate
            (
                meshPtr(),
                geometryDict,
                refineDict,
                baseCheck
            );
        }
    }

    if (controller.block() != meshControl::IMPORT)
    {
        //create dummy constant folder if not present
        if
        (
            Foam::fileHandler().type() == "uncollated"
         && !isDir(runTime_.path()/"constant"/regionDir)
        )
        {
            mkDir(runTime_.path()/"constant"/regionDir);
        }
    }

    // Read geometry
    // ~~~~~~~~~~~~~
    bool singleRegionName(true);

    if (meshDict.found("appendRegionName"))
    {
        singleRegionName = !meshDict.lookup("appendRegionName");
    }
    else if (meshDict.found("singleRegionName"))
    {
        singleRegionName = meshDict.lookup("singleRegionName");
    }

    label maxRegionSize =
        meshDict.lookupOrDefault<label>("maxRegionSize", labelMax);

    List<fileName> alternativeSurfacePath =
        meshDict.lookupOrDefault<List<fileName>>
        ("alternativeSurfacePath", List<fileName>(0));

    autoPtr<searchableSurfaces> allGeometryPtr
    (
        new searchableSurfaces
        (
            IOobject
            (
                "abc",                      // dummy name
                runTime_.time().constant(),  // directory
                "triSurface",               // instance
                runTime_.time(),             // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            geometryDict,
            singleRegionName,
            alternativeSurfacePath,
            maxRegionSize
        )
    );

    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<refinementSurfaces> surfacesPtr
    (
        new refinementSurfaces
        (
            allGeometryPtr(),
            refineDict.subDict("refinementSurfaces"),
            refineDict,
            refineDict.lookupOrDefault("gapLevelIncrement", 0)
         )
    );
    refinementSurfaces& surfaces = surfacesPtr();

    Info<< "Read refinement surfaces in = "
        << runTime_.time().cpuTimeIncrement() << " s" << nl << endl;

    // Refinement parameters
    refinementParameters refineParams(refineDict);

    //Check for global geometry dict
    IOobject globalGeomHeader
    (
        "geometryDict",
        runTime_.time().system(),
        runTime_.time(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    );

    dictionary globalGeomDict;
    if (globalGeomHeader.typeHeaderOk<IOdictionary>(true))
    {
        globalGeomDict = IOdictionary(globalGeomHeader);
        if (globalGeomDict.found("surfaces"))
        {
            refineParams.addToLocationsInMesh
            (
                allGeometryPtr(),
                globalGeomDict.subDict("surfaces")
            );
        }
    }

    refineParams.addToLocationsInMesh
    (
        allGeometryPtr(),
        geometryDict
    );

    //If curvature active set curvature fields (with optional smoothing)
    surfaces.setCurvatureFields(refineParams.curvatureSmooth());

    //If distributed tri surfaces then distribute
    forAll(allGeometryPtr(), i)
    {
        List<treeBoundBox> bbs;
        autoPtr<mapDistribute> faceMap;
        autoPtr<mapDistribute> pointMap;
        allGeometryPtr()[i].distribute
        (
            bbs,
            false,          // do not keep outside triangles
            faceMap,
            pointMap
        );
    }

    // Checking only?

    if (checkGeometry)
    {
        // Extract patchInfo
        List<wordList> patchTypes(allGeometryPtr().size());

        scalar intersectionTol =
            meshDict.lookupOrDefault<scalar>("intersectionTol", -1.0);

        const PtrList<dictionary>& patchInfo = surfaces.patchInfo();
        const labelList& surfaceGeometry = surfaces.surfaces();
        forAll(surfaceGeometry, surfi)
        {
            label geomi = surfaceGeometry[surfi];
            const wordList& regNames = allGeometryPtr().regionNames()[geomi];

            patchTypes[geomi].setSize(regNames.size());
            forAll(regNames, regioni)
            {
                label globalRegioni = surfaces.globalRegion(surfi, regioni);

                if (patchInfo.set(globalRegioni))
                {
                    patchTypes[geomi][regioni] =
                        word(patchInfo[globalRegioni].lookup("type"));
                }
                else
                {
                    patchTypes[geomi][regioni] = wallPolyPatch::typeName;
                }
            }
        }

        // Write some stats
        allGeometryPtr().writeStats(patchTypes, Info);

        // Check topology
        allGeometryPtr().checkTopology(true);
        // Check geometry
        allGeometryPtr().checkGeometry
        (
            100.0,      // max size ratio
            intersectionTol,    // intersection tolerance
            autoPtr<writer<scalar>>(new vtkSetWriter<scalar>()),
            0.01,       // min triangle quality
            true
        );

        return meshPtr;
    }

    // Create mesh (automatically)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        if (controller.block() == meshControl::AUTO || useVDB)
        {
            meshPtr.reset
            (
                new fvMesh
                (
                    Foam::IOobject
                    (
                        regionName,
                        runTime_.timeName(),
                        runTime_,
                        Foam::IOobject::MUST_READ
                    ),
                    xferCopy(pointField()),
                    xferCopy(faceList()),
                    xferCopy(cellList())
                )
            );

            const Tuple2<scalar, label> blockData =
                meshDict.lookup("blockData");

            if (blockData.first() < SMALL)
            {
                FatalErrorInFunction
                    << "Length scale " << blockData.first()
                    << " is too small for generating an autoBlockMesh "
                    << exit(FatalError);
            }

            if (controller.block() == meshControl::AUTO)
            {
                if (!surfaces.checkForFiniteSurfaces())
                {
                    FatalErrorInFunction
                        << " The surface description contains no triSurfaces "
                        << "to define the automatic block mesh dimensions"
                        << exit(FatalError);
                }

                autoBlockMesh baseMesh
                (
                    surfaces,
                    meshDict,
                    meshPtr()
                );
            }
            else if (useVDB)
            {
                //add empty patches of VDBdomain
                addVDBpatches
                (
                    meshDict.subDict("VDBdomain"),
                    meshPtr()
                );
            }
        }

        // Check patches and faceZones are synchronised
        meshPtr().boundaryMesh().checkParallelSync(true);
    }

    fvMesh& mesh = meshPtr();

    const word oldInstance = mesh.pointsInstance();

    const scalar mergeDist = getMergeDistance
    (
        mesh,
        (useVDB ? allGeometryPtr().bounds() : mesh.bounds()),
        readScalar(meshDict.lookup("mergeTolerance"))
    );

    //Stop unnecessary calculations during mesh movement, etc.
    mesh.preventSolverQtyCalc();

    // Optionally read limit shells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const dictionary limitDict(refineDict.subOrEmptyDict("limitRegions"));

    if (!limitDict.empty())
    {
        Info<< "Reading limit shells." << endl;
    }

    shellSurfaces limitShells(allGeometryPtr(), limitDict, controller);

    if (!limitDict.empty())
    {
        Info<< "Read refinement shells in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }

    // Redistribute block mesh if automatically generated
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Mesh distribution engine (uses tolerance to reconstruct meshes)
    fvMeshDistribute distributor(mesh, mergeDist);

    if (controller.block() == meshControl::AUTO && Pstream::parRun())
    {
        labelList distribution =
           decomposer.decompose(mesh, mesh.cellCentres());

        // Do actual sending/receiving of mesh
        distributor.distribute(distribution);

        Random rndGen(653213);

        // Get local mesh bounding box. Single box for now.
        List<treeBoundBox> meshBb(1);
        treeBoundBox& bb = meshBb[0];
        bb = treeBoundBox(mesh.points());
        bb = bb.extend(rndGen, 1E-4);

        forAll(allGeometryPtr(), i)
        {
            autoPtr<mapDistribute> faceMap;
            autoPtr<mapDistribute> pointMap;
            allGeometryPtr()[i].distribute
            (
                meshBb,
                false,          // do not keep outside triangles
                faceMap,
                pointMap
            );

            if (faceMap.valid())
            {
                // (ab)use the instance() to signal current modification time
                allGeometryPtr()[i].instance() =
                    allGeometryPtr()[i].time().timeName();
            }

            faceMap.clear();
            pointMap.clear();
        }
    }

    // Read refinement features
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Creating refinement features." << endl;
    autoPtr<refinementFeatures> featuresPtr
    (
        new refinementFeatures
        (
            mesh,
            refineDict.lookup("features"),
            surfaces,
            refineParams,
            mergeDist // distributedTriSurfaceMesh reconstruct edge tol
         )
    );
    Info<< "Calculated refinement features in = "
        << runTime_.time().cpuTimeIncrement() << " s" << endl;

    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement shells." << endl;
    autoPtr<shellSurfaces> shellsPtr
    (
        new shellSurfaces
        (
            allGeometryPtr(),
            refineDict.subDict("refinementRegions"),
            controller,
            refineDict.lookupOrDefault<bool>("additionalInsideCheck",false)
        )
    );
    Info<< "Read refinement shells in = "
        << runTime_.time().cpuTimeIncrement() << " s" << nl << endl;

    Tuple2<bool, scalar> crackDetection;
    crackDetection.first() = false;
    crackDetection.second() = 0.05;
    Switch addedRays(false);
    scalar manualLevel0 = -GREAT;
    bool readLevel = true;

    if (controller.vdb())
    {
        if (controller.refine())
        {
            readLevel = false;

            // do not trigger hexRef8::getLevel0EdgeLength()
            manualLevel0 = GREAT;
        }
    }
    else
    {
        const Switch cacheShellLevel
        (
            meshDict.lookupOrDefault("cacheShellLevel", true)
        );


        Info<< "Setting refinement level of surface to be consistent"
            << " with shells." << endl;

        if (cacheShellLevel)
        {
            surfaces.setMinLevelFields(shellsPtr);
        }
        else
        {
            surfaces.setMinLevelFields(autoPtr<shellSurfaces>(nullptr));
        }

        Info<< "Checked shell refinement in = "
            << runTime_.time().cpuTimeIncrement() << " s" << nl << endl;


        // Refinement engine
        // ~~~~~~~~~~~~~~~~~

        Info<< nl
            << "Determining initial surface intersections" << nl
            << "-----------------------------------------" << nl
            << endl;

        if (meshDict.found("crackDetection"))
        {
           crackDetection.first() = readBool(meshDict.lookup("crackDetection"));
        }
        if (meshDict.found("crackTol"))
        {
           crackDetection.second() = readScalar(meshDict.lookup("crackTol"));
        }
        if (meshDict.found("extraRays"))
        {
           addedRays = readBool(meshDict.lookup("extraRays"));
        }

        manualLevel0 =
            meshDict.lookupOrDefault<scalar>("manualLevel0EdgeLength", -GREAT);

        if (controller.block() != meshControl::IMPORT)
        {
            readLevel = false;
        }
    }

    // Main refinement engine
    meshRefinement meshRefiner
    (
        mesh,
        mergeDist,        // tolerance used in sorting coordinates
        readLevel,        // whether to try to read level information
        overwrite,        // overwrite mesh files?
        controller,       // control of mesh generation process
        crackDetection,   // crack detection during refinement?
        addedRays,        // additional rays for crack detection?
        surfaces,         // for surface intersection refinement
        shellsPtr(),      // for volume (inside/outside) refinement
        limitShells,      // limit of volume refinement
        featuresPtr,      // for feature refinement
        manualLevel0,     // optional overwrite of level 0 edge length
        useVDB
    );

    if (!useVDB)
    {
        Info<< "Calculated surface intersections in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;

        if (controller.block() == meshControl::AUTO)
        {
            dictionary patchInfo;
            patchInfo.set("type", wallPolyPatch::typeName);

            meshRefiner.addMeshedPatch("blockMesh",patchInfo);
        }

        // Some stats
        meshRefiner.printMeshInfo(debugLevel, "Initial mesh");

        meshRefiner.write
        (
            meshRefinement::debugType
            (
                debugLevel&meshRefinement::OBJINTERSECTIONS
            ),
            meshRefinement::writeType(0),
            mesh.time().path()/meshRefiner.timeName()
        );
    }

    // Add all the cellZones and faceZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. cellZones relating to surface (faceZones added later)

    const labelList namedSurfaces
    (
        surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones())
    );

    labelList surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
    (
        surfaces.surfZones(),
        namedSurfaces,
        mesh
    );


    // 2. cellZones relating to locations

    refineParams.addCellZonesToMesh(mesh);



    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Global surface region to patch (non faceZone surface) or patches
    //  (faceZone surfaces)
    labelList globalToMasterPatch;
    labelList globalToSlavePatch;


    {
        Info<< nl
            << "Adding patches for surface regions" << nl
            << "----------------------------------" << nl
            << endl;

        // From global region number to mesh patch.
        globalToMasterPatch.setSize(surfaces.nRegions(), -1);
        globalToSlavePatch.setSize(surfaces.nRegions(), -1);

        Info<< setf(ios_base::left)
            << setw(6) << "Patch"
            << setw(20) << "Type"
            << setw(30) << "Region" << nl
            << setw(6) << "-----"
            << setw(20) << "----"
            << setw(30) << "------" << endl;

        const labelList& surfaceGeometry = surfaces.surfaces();
        const PtrList<dictionary>& surfacePatchInfo = surfaces.patchInfo();
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        forAll(surfaceGeometry, surfi)
        {
            label geomi = surfaceGeometry[surfi];

            const wordList& regNames = allGeometryPtr().regionNames()[geomi];

            Info<< surfaces.names()[surfi] << ':' << nl << nl;

            const word& fzName = surfaces.surfZones()[surfi].faceZoneName();

            if (fzName.empty())
            {
                // 'Normal' surface
                forAll(regNames, i)
                {
                    label globalRegioni = surfaces.globalRegion(surfi, i);

                    label patchi;

                    if (surfacePatchInfo.set(globalRegioni))
                    {
                        patchi = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            surfacePatchInfo[globalRegioni]
                        );
                    }
                    else
                    {
                        dictionary patchInfo;
                        patchInfo.set("type", wallPolyPatch::typeName);

                        patchi = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            patchInfo
                        );
                    }

                    Info<< setf(ios_base::left)
                        << setw(6) << patchi
                        << setw(20) << pbm[patchi].type()
                        << setw(30) << regNames[i] << nl;

                    globalToMasterPatch[globalRegioni] = patchi;
                    globalToSlavePatch[globalRegioni] = patchi;
                }
            }
            else
            {
                // Zoned surface
                forAll(regNames, i)
                {
                    label globalRegioni = surfaces.globalRegion(surfi, i);

                    // Add master side patch
                    {
                        label patchi;

                        if (surfacePatchInfo.set(globalRegioni))
                        {
                            patchi = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                surfacePatchInfo[globalRegioni]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchi = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                patchInfo
                            );
                        }

                        Info<< setf(ios_base::left)
                            << setw(6) << patchi
                            << setw(20) << pbm[patchi].type()
                            << setw(30) << regNames[i] << nl;

                        globalToMasterPatch[globalRegioni] = patchi;
                    }
                    // Add slave side patch
                    {
                        const word slaveName = regNames[i] + "_slave";
                        label patchi;

                        if (surfacePatchInfo.set(globalRegioni))
                        {
                            patchi = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                surfacePatchInfo[globalRegioni]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchi = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                patchInfo
                            );
                        }

                        Info<< setf(ios_base::left)
                            << setw(6) << patchi
                            << setw(20) << pbm[patchi].type()
                            << setw(30) << slaveName << nl;

                        globalToSlavePatch[globalRegioni] = patchi;
                    }
                }

                // For now: have single faceZone per surface. Use first
                // region in surface for patch for zoneing
                if (regNames.size())
                {
                    label globalRegioni = surfaces.globalRegion(surfi, 0);

                    meshRefiner.addFaceZone
                    (
                        fzName,
                        pbm[globalToMasterPatch[globalRegioni]].name(),
                        pbm[globalToSlavePatch[globalRegioni]].name(),
                        surfaces.surfZones()[surfi].faceType()
                    );
                }
            }

            Info<< nl;
        }
        Info<< "Added patches in = "
            << runTime_.time().cpuTimeIncrement() << " s" << nl << endl;
    }



    // Add all information for all the remaining faceZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HashTable<Pair<word>> faceZoneToPatches;
    forAll(mesh.faceZones(), zonei)
    {
        const word& fzName = mesh.faceZones()[zonei].name();

        label mpI, spI;
        surfaceZonesInfo::faceZoneType fzType;
        bool hasInfo = meshRefiner.getFaceZoneInfo(fzName, mpI, spI, fzType);

        if (!hasInfo)
        {
            // faceZone does not originate from a surface but presumably
            // from a cellZone pair instead
            string::size_type i = fzName.find("_to_");
            if (i != string::npos)
            {
                word cz0 = fzName.substr(0, i);
                word cz1 = fzName.substr(i+4, fzName.size()-i+4);
                word slaveName(cz1 + "_to_" + cz0);
                faceZoneToPatches.insert(fzName, Pair<word>(fzName, slaveName));
            }
            else
            {
                // Add as fzName + fzName_slave
                const word slaveName = fzName + "_slave";
                faceZoneToPatches.insert(fzName, Pair<word>(fzName, slaveName));
            }
        }
    }

    if (faceZoneToPatches.size())
    {
        snappyRefineDriver::addFaceZones
        (
            meshRefiner,
            refineParams,
            faceZoneToPatches
        );
    }

    // Re-do intersections on meshed boundaries since they use an extrapolated
    // other side
    if (!useVDB)
    {
        const labelList adaptPatchIDs(meshRefiner.meshedPatches());

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        label nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            nFaces += pbm[adaptPatchIDs[i]].size();
        }

        labelList faceLabels(nFaces);
        nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            const polyPatch& pp = pbm[adaptPatchIDs[i]];
            forAll(pp, i)
            {
                faceLabels[nFaces++] = pp.start()+i;
            }
        }
        meshRefiner.updateIntersections(faceLabels);
    }

    // Now do the real work -refinement -snapping -layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<snappyRefineDriver> refineDriver;
    autoPtr<snappySnapDriver> snapDriver;

    snapParameters snapParams
    (
        snapDict,
        mesh.boundaryMesh()
     );

    // create reporting object
    hexReport stats(meshRefiner);
    stats.setInitialisationTime(mesh.time().elapsedCpuTime());

    if (controller.refine())
    {
        scalar startTime = mesh.time().elapsedClockTime();

        refineDriver.reset
        (
            new snappyRefineDriver
            (
                meshRefiner,
                decomposer,
                distributor,
                globalToMasterPatch,
                globalToSlavePatch,
                controller,
                setFormatter
             )
        );

        if (controller.algorithm() == meshControl::SHELL)
        {
            refineDriver->doThinShellRefine
            (
                meshDict,
                refineParams
            );
        }
        else if (useVDB)
        {
#ifdef FOAM_USE_OPENVDB
            refineDriver->doVDBRefine
            (
                meshDict,
                refineParams,
                snapParams,
                refineParams.handleSnapProblems()
            );

            // Redistribute to wanted number of procs
            if (numSubDomains != Pstream::nProcs())
            {
                if (numSubDomains > Pstream::nProcs())
                {
                    meshRefiner.redistributeToMany
                    (
                        decomposer,
                        distributor
                    );
                }
                else
                {
                    NotImplemented;
                }

                Info<< "\nPlease restart foamHexMesh for snap and "
                    << "layer addition using " << numSubDomains
                    << " processors.\n\nEnd."
                    << endl;

                return meshPtr;
            }
#else
            FatalErrorInFunction
                << "OpenVDB has not been compiled!"
                << exit(FatalError);
#endif
        }
        else
        {
            refineDriver->doRefine
            (
                meshDict,
                refineParams,
                snapParams,
                refineParams.handleSnapProblems()
            );
        }

        if ((!overwrite && !debugLevel) || writeIntermediate)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if (meshWrite && (debugLevel || writeIntermediate))
        {
            writeMesh
            (
                "Refined mesh",
                meshRefiner,
                debugLevel,
                meshRefinement::writeLevel(),
                false,
                false
            );
        }

        Info<< "Mesh refined in = "
            << mesh.time().elapsedClockTime() - startTime << " s." << endl;

        if (useVDB)
        {
            Info<< "Finished VDB refinement in = "
                << mesh.time().elapsedClockTime() << " s." << endl;

            Info<< "End\n" << endl;
            return meshPtr;
        }
        else
        {
            Info<< "Total meshing time to this stage = "
                << mesh.time().elapsedClockTime() << " s." << endl;
        }
#if !defined( WIN32 ) && !defined( WIN64 )
        profiling::writeNow();
#endif
        stats.setRefineTime(mesh.time().elapsedCpuTime()-startTime);
    }

    //Create any fields needed during the meshing process
    //If dual mesher maintain field of layer cells (>-1)
    autoPtr<volScalarField> layerCells;

    const Switch outputLayerCells
    (
        meshDict.lookupOrDefault("outputLayerCells", false)
    );

    if (controller.algorithm() != meshControl::STANDARD)
    {
        layerCells.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "layerStacks",
                    runTime_.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    outputLayerCells ?
                    IOobject::AUTO_WRITE : IOobject::NO_WRITE,
                    true
                ),
                mesh,
                dimensionedScalar("minusone", dimless, -1),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    //If wanting to perform Adaptive Mesh Refinement (AMR)
    // add optional output field of protected cells (>-1) which cannot
    // be refined
    const Switch outputProtected
    (
        meshDict.lookupOrDefault("outputProtected", false)
    );

    const Switch smoothProtectedCells
    (
        meshDict.lookupOrDefault("smoothLayerCells", true)
    );

    autoPtr<volScalarField> protectedCells;
    if (outputProtected || !smoothProtectedCells)
    {
        protectedCells.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "protectedCells",
                    runTime_.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    outputProtected ? IOobject::AUTO_WRITE : IOobject::NO_WRITE,
                    true
                ),
                mesh,
                dimensionedScalar("minusone", dimless, -1),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        if (outputProtected && controller.algorithm() == meshControl::STANDARD)
        {
            const labelList meshedPatches = meshRefiner.meshedPatches();
            autoPtr<indirectPrimitivePatch> ppPtr
            (
                meshRefinement::makePatch
                (
                    mesh,
                    meshedPatches
                )
            );

            indirectPrimitivePatch& pp = ppPtr();
            const labelList& meshPoints = pp.meshPoints();
            forAll(meshPoints, ptI)
            {
                label meshPointI = meshPoints[ptI];
                const labelList& pCells = mesh.pointCells()[meshPointI];
                forAll(pCells, pCellI)
                {
                    label cellI = pCells[pCellI];
                    if (protectedCells()[cellI] == -1)
                    {
                        protectedCells()[cellI] = scalar(1);
                    }
                }
            }
        }
    }

    if (controller.algorithm() == meshControl::DUAL)
    {
        scalar startTime = mesh.time().elapsedCpuTime();
        const labelList meshedPatches = meshRefiner.meshedPatches();

        //Add internal layers, boundary layers and then dualise
        {
            addHexInternalLayers addInternalLayers
            (
                meshRefiner,
                decomposer,
                distributor
            );
            addInternalLayers.setRefinement();

            if ((!overwrite && !debugLevel) || writeIntermediate)
            {
                const_cast<Time&>(mesh.time())++;
            }

            if (meshWrite && (debugLevel || writeIntermediate))
            {
                writeMesh
                (
                    "Refined Internal Layers",
                    meshRefiner,
                    debugLevel,
                    meshRefinement::writeLevel(),
                    false,
                    false
                );
            }

            Info<<"Number of cells after interface refine : "
                <<returnReduceToMaster(mesh.nCells(), sumOp<label>())<<endl;

            //meshRefiner.removeOtherSide(refineParams,globalToMasterPatch);

            //remove cells that might cause problems for dualisation
            //meshRefiner.removeProblemDualisationCells(true);

            autoPtr<indirectPrimitivePatch> ppPtr
            (
                meshRefinement::makePatch
                (
                    mesh,
                    meshedPatches
                )
            );

            indirectPrimitivePatch& pp = ppPtr();

            addHexMeshLayer addRefLayers
            (
                globalToMasterPatch,
                refineParams,
                decomposer,
                distributor,
                meshRefiner,
                pp,
                0.66,//0.5
                true
            );
            addRefLayers.setRefinement();

            if ((!overwrite && !debugLevel) || writeIntermediate)
            {
                const_cast<Time&>(mesh.time())++;
            }

            if (meshWrite && (debugLevel || writeIntermediate))
            {
                writeMesh
                (
                    "Refined Boundary Layers",
                    meshRefiner,
                    debugLevel,
                    meshRefinement::writeLevel(),
                    false,
                    false
                );
            }

            Info<<"Number of cells after adding boundary cells : "
                <<returnReduceToMaster(mesh.nCells(), sumOp<label>())<<endl;

            autoHexDualiser dualiser(meshRefiner);

            //Need to duplicate non-manifold points created during dualisation
            meshRefiner.dupNonManifoldPoints();

            meshRefiner.removeDisconnectedRegions();

            //Move cells on wrong side of geometry
            meshRefiner.moveWrongSidedCells(refineParams,false);

            dictionary meshOptimDict
            (
                meshDict.found("meshOptimization") ?
                meshDict.subDict("meshOptimization") :
                dictionary()
            );

            if (meshOptimDict.found("foamOptimizeCoeffs"))
            {
                dictionary& coeffsDict =
                    meshOptimDict.subDict("foamOptimizeCoeffs");
                if (!coeffsDict.found("meshQualityControls"))
                {
                    coeffsDict.add("meshQualityControls",motionDict,true);
                }
            }

            autoPtr<autoOptimize> optimMeshPtr
                = autoOptimize::New(mesh, meshOptimDict);
            optimMeshPtr->optimize();

            //Perform optional zoning on non boundary zones
            refineDriver-> zonifyDual(refineParams,true);

            meshRefiner.balance
            (
                false,
                true, //keep zone face together
                false,
                scalarField(mesh.nCells(), 1.0),
                decomposer,
                distributor,
                false
            );
        }

        stats.setDualSetupTime(mesh.time().elapsedCpuTime()-startTime);
    }
    else if (controller.algorithm() == meshControl::EXTRUDE)
    {
        if (surfaces.cornerCellSurfaces().size() > 0)
        {
            //Move cells on wrong side of geometry if corner
            // cell selection method is used
            meshRefiner.moveWrongSidedCells(refineParams,true);
        }

        scalar startTime = mesh.time().elapsedCpuTime();

        if (snapParams.preMergeExtrude())
        {
            snapDriver.reset
            (
                new snappySnapDriver
                (
                    meshRefiner,
                    decomposer,
                    distributor,
                    meshDict,
                    controller,
                    globalToMasterPatch,
                    globalToSlavePatch
                 )
            );

            snapParams.calculateRegionAndFeatureSnapPatches
            (
                snapDict,
                meshRefiner
            );

            snapDriver->snapAndMerge(snapParams);
        }

        //Extrude zero sized layer cells on all boundary patches and expand
        autoExtrude extrude(meshRefiner,decomposer,distributor);
        extrude.extrudeZeroSizedLayer(refineParams);

        dictionary meshOptimDict
        (
            meshDict.found("meshOptimization") ?
            meshDict.subDict("meshOptimization") :
            dictionary()
        );

        if (meshOptimDict.found("foamOptimizeCoeffs"))
        {
            dictionary& coeffsDict =
                meshOptimDict.subDict("foamOptimizeCoeffs");
            if (!coeffsDict.found("meshQualityControls"))
            {
                coeffsDict.add("meshQualityControls",motionDict,true);
            }
        }

        if (meshOptimDict.found("cfMeshOptimizeCoeffs"))
        {
            meshOptimDict.subDict("cfMeshOptimizeCoeffs").add
            (
                "relaxedCheck",
                "false",
                true
             );
        }

        autoPtr<autoOptimize> optimMeshPtr
            = autoOptimize::New(mesh, meshOptimDict);
        optimMeshPtr->optimize();
        stats.setDualSetupTime(mesh.time().elapsedCpuTime()-startTime);
    }

    // initialise layer cell field
    if
    (
        controller.algorithm() == meshControl::EXTRUDE
     || controller.algorithm() == meshControl::DUAL
    )
    {
        const labelList meshedPatches = meshRefiner.meshedPatches();

        autoPtr<indirectPrimitivePatch> ppPtr
        (
            meshRefinement::makePatch
            (
                mesh,
                meshedPatches
            )
        );
        indirectPrimitivePatch& pp = ppPtr();
        scalar maxLayerIndex = gMax(layerCells());

        if (maxLayerIndex == -1)
        {
            const labelList& meshPoints = pp.meshPoints();
            label nLayerCells = 0;

            forAll(meshPoints, ptI)
            {
                label meshPointI = meshPoints[ptI];
                const labelList& pCells = mesh.pointCells()[meshPointI];
                forAll(pCells, pCellI)
                {
                    label cellI = pCells[pCellI];
                    if (layerCells()[cellI] == -1)
                    {
                        layerCells()[cellI] = nLayerCells++;
                    }
                }
            }

            globalIndex globalLayerCells(nLayerCells);

            forAll(mesh.cells(), cellI)
            {
                if (layerCells()[cellI] != -1)
                {
                    layerCells()[cellI] +=
                        globalLayerCells.offset(Pstream::myProcNo());
                }
            }
        }

        if ((!overwrite && !debugLevel) || writeIntermediate)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if (meshWrite && (debugLevel || writeIntermediate))
        {
            writeMesh
            (
                "Pre snapped mesh",
                meshRefiner,
                debugLevel,
                meshRefinement::writeLevel(),
                false,
                false
            );
        }
    }

    if (controller.snap())
    {
        scalar startTime = mesh.time().elapsedCpuTime();

#ifdef FOAM_USE_OPENVDB
        if (controller.vdb())
        {
            refineDriver.reset
            (
                new snappyRefineDriver
                (
                    meshRefiner,
                    decomposer,
                    distributor,
                    globalToMasterPatch,
                    globalToSlavePatch,
                    controller,
                    setFormatter
                 )
            );

            refineDriver->doVDBzones
            (
                meshDict,
                refineParams,
                snapParams,
                refineParams.handleSnapProblems()
            );

            Info<< "Created baffles and zones from VDB mesh in = "
                << mesh.time().elapsedCpuTime() - startTime << " s." << endl;
            Info<< "Total meshing time to this stage = "
                << mesh.time().elapsedCpuTime() << " s." << endl;
        }
#endif


        startTime = mesh.time().elapsedCpuTime();

        snapDriver.reset
        (
            new snappySnapDriver
            (
                meshRefiner,
                decomposer,
                distributor,
                meshDict,
                controller,
                globalToMasterPatch,
                globalToSlavePatch
            )
        );

        snapParams.calculateRegionAndFeatureSnapPatches
        (
            snapDict,
            meshRefiner
        );

        snapDriver->doSnap
        (
            snapParams,
            refineParams
        );

        Info<< "Finished snapping stage in = "
            << mesh.time().elapsedCpuTime() - startTime << " s." << endl;
        Info<< "Total meshing time to this stage = "
            << mesh.time().elapsedCpuTime() << " s." << endl;

#if !defined( WIN32 ) && !defined( WIN64 )
        profiling::writeNow();
#endif

        stats.setSnapTime(mesh.time().elapsedCpuTime()-startTime);


        if
        (
            controller.algorithm() == meshControl::EXTRUDE
            || controller.algorithm() == meshControl::DUAL
        )
        {
            //check for errors and perform different optimization
            dictionary meshOptimDict
            (
                meshDict.found("meshOptimization") ?
                meshDict.subDict("meshOptimization") :
                dictionary()
            );

            //If mesh still has errors try cfMesh optimization
            cfOptimize(mesh, meshOptimDict);

            //Split boundary feature cells
            if (controller.topoChanges())
            {
                // Layer addition parameters
                const polyBoundaryMesh& patches = mesh.boundaryMesh();
                const dictionary& layerDict =
                    meshDict.subDict("addLayersControls");
                layerParameters layerParams(layerDict,mesh.boundaryMesh(),true);
                const labelList& numLayers = layerParams.numLayers();
                DynamicList<label>  excludePatches(patches.size());
                forAll(patches, patchI)
                {
                    if (!patches[patchI].coupled() && numLayers[patchI] == -1)
                    {
                        excludePatches.append(patchI);
                    }
                }
                excludePatches.shrink();

                autoSplitCells splitCells
                (
                    meshRefiner,
                    excludePatches,
                    (controller.algorithm() == meshControl::EXTRUDE),
                    snapParams.preMergeExtrude()
                );
                splitCells.splitFeatureCells
                (
                    scalar(0.9659),
                    true,
                    true
                );
            }
        }

        if ((!overwrite && !debugLevel) || writeIntermediate)
        {
            const_cast<Time&>(mesh.time())++;
        }

        if (meshWrite && (debugLevel || writeIntermediate))
        {
            writeMesh
            (
                "Snapped mesh",
                meshRefiner,
                debugLevel,
                meshRefinement::writeLevel(),
                false,
                false
            );
        }
    }

    bool addedLayers(false);

    //Clear up memory before adding layers
    shellsPtr.clear();
    featuresPtr.clear();
    mesh.clearOut();

    if (controller.layers())
    {
        scalar startTime = mesh.time().elapsedCpuTime();

        if
        (
            controller.algorithm() == meshControl::DUAL
            || controller.algorithm() == meshControl::EXTRUDE
        )
        {
            // layer addition parameters
            const polyBoundaryMesh& patches = mesh.boundaryMesh();
            const dictionary& layerDict =
                meshDict.subDict("addLayersControls");
            layerParameters layerParams(layerDict,mesh.boundaryMesh(),true);

            if (controller.topoChanges())
            {
                const labelList& numLayers = layerParams.numLayers();
                DynamicList<label>  excludePatches(patches.size());
                forAll(patches, patchI)
                {
                    if (!patches[patchI].coupled() && numLayers[patchI] == -1)
                    {
                        excludePatches.append(patchI);
                    }
                }
                excludePatches.shrink();

                scalar warpage = layerParams.globalWarpedSplit();
                scalar fchwarpage = layerParams.fchWarpedSplit();

                //Split warped faces
                if (warpage > 0 || fchwarpage > 0)
                {
                    bool curvatureSplit = layerParams.curvatureSplit();
                    const List<wordList>& layerSpec = layerParams.layerSpec();
                    const labelList& patchNumLayers = layerParams.numLayers();
                    const scalarField& patchFCH = layerParams.fch();
                    const scalarField& patchExpansionRatio =
                        layerParams.expansionRatio();

                    autoSplitCells splitCells
                    (
                        meshRefiner,
                        excludePatches,
                        (controller.algorithm() == meshControl::EXTRUDE),
                        snapParams.preMergeExtrude(),
                        false
                     );

                    splitCells.splitWarpedCells
                    (
                        warpage,
                        fchwarpage,
                        layerSpec,
                        patchNumLayers,
                        patchFCH,
                        patchExpansionRatio,
                        curvatureSplit,
                        true, //perform additional convex split
                        snapParams.preMergeExtrude()
                     );

                    //optimize final mesh disabling Laplacian smoothing
                    if (controller.mode() != meshControl::DRYRUN)
                    {
                        dictionary meshOptimDict = finalOptimizeSetup(meshDict);
                        autoPtr<autoOptimize> optimMeshPtr
                            = autoOptimize::New(mesh, meshOptimDict);
                        optimMeshPtr->optimize();
                    }
                }
            }

            //Return dictionary of non-boundary zoned geometry entries
            //to re-snap to
            dictionary reSnapZoneGeometryDict = calcReSnapNonBoundaryZoneDict
            (
                layerParams,
                surfaces,
                geometryDict
            );

            //Return dictionary of geometry entries to re-snap to
            dictionary reSnapGeometryDict = calcReSnapGeometryDict
            (
                layerParams,
                surfaces,
                globalToMasterPatch,
                globalToSlavePatch,
                geometryDict
            );

            //Free up triangulated surface storage
            allGeometryPtr.clear();

            addMultiLayers addRefLayers
            (
                controller,
                globalToMasterPatch,
                refineParams,
                layerParams,
                reSnapGeometryDict,
                reSnapZoneGeometryDict,
                decomposer,
                distributor,
                meshRefiner
            );

            if (controller.algorithm() == meshControl::EXTRUDE)
            {
                addRefLayers.setRefinementExtrude();
            }
            else
            {
                addRefLayers.setRefinementDual();
            }

            addedLayers = true;

            word optimType("cfMeshOptimize");
            if (controller.mode() != meshControl::DRYRUN)
            {
                dictionary meshOptimDict = finalOptimizeSetup(meshDict);
                optimType = meshOptimDict.lookupOrDefault<word>
                (
                    "type",
                    "cfMeshOptimize"
                );
                autoPtr<autoOptimize> optimMeshPtr
                    = autoOptimize::New(mesh, meshOptimDict);
                optimMeshPtr->optimize();

                label maxMergeIter = layerParams.maxMergePostIter();
                if (maxMergeIter > 0)
                {
                    autoLayerCellsMerge autoMerge
                    (
                        meshRefiner,
                        decomposer,
                        distributor,
                        layerParams
                     );
                    autoMerge.merge(maxMergeIter);
                }

                //If mesh still has errors try alternative optimization
                cfOptimize(mesh, meshOptimDict);
            }

            layerManipulate layerManip
            (
                mesh,
                layerParams,
                meshRefiner.meshCutter().cellLevel(),
                meshRefiner.meshCutter().pointLevel(),
                meshRefiner.meshCutter().level0EdgeLength()
            );

            if (optimType == "foamOptimize")
            {
                layerManip.fitLayerPointStack();
            }

            if (outputProtected)
            {
                addProtectedCellsExtrude
                (
                    mesh,
                    protectedCells()
                );
            }

            //Write out layer info
            layerManip.writeLayerInfo();
            stats.setLayerCoverage(100.0);

            label numLayerCells = 0;
            forAll(layerCells(), cellI)
            {
                if (layerCells()[cellI] != -1)
                {
                    numLayerCells++;
                }
            }

            reduce(numLayerCells, sumOp<label>());
            stats.setNumLayerCells(numLayerCells);

            if (layerDict.found("addFaceZoneLayersControls"))
            {
                // face zone layer addition parameters
                const dictionary& zoneLayerDict =
                    layerDict.subDict("addFaceZoneLayersControls");
                addedLayers = addLayers
                (
                    zoneLayerDict,
                    motionDict,
                    controller,
                    globalToMasterPatch,
                    globalToSlavePatch,
                    overwrite,
                    meshRefiner,
                    allGeometryPtr,
                    stats,
                    decomposer,
                    distributor,
                    label(1) //Treat as multi-layer and increment count
                );
            }
        }
        else
        {
            if (meshDict.found("addMultipleLayersControls"))
            {
                const List<dictionary> multiLayerDict
                (
                    meshDict.lookup("addMultipleLayersControls")
                );

                forAll(multiLayerDict, i)
                {
                    addedLayers = addLayers
                    (
                        multiLayerDict[i],
                        motionDict,
                        controller,
                        globalToMasterPatch,
                        globalToSlavePatch,
                        overwrite,
                        meshRefiner,
                        allGeometryPtr,
                        stats,
                        decomposer,
                        distributor,
                        i //Multi layer count
                    );
                }
            }
            else
            {
                // layer addition parameters
                const dictionary& layerDict =
                    meshDict.subDict("addLayersControls");
                addedLayers = addLayers
                (
                    layerDict,
                    motionDict,
                    controller,
                    globalToMasterPatch,
                    globalToSlavePatch,
                    overwrite,
                    meshRefiner,
                    allGeometryPtr,
                    stats,
                    decomposer,
                    distributor,
                    label(0)
                );
            }

            if (meshWrite && debugLevel)
            {
                writeMesh
                (
                    "Layer mesh",
                    meshRefiner,
                    debugLevel,
                    meshRefinement::writeLevel(),
                    false,
                    false
                 );
            }

            Info<< "Finished layer stage in = "
                << mesh.time().elapsedCpuTime() - startTime << " s." << endl;
            Info<< "Total meshing time to this stage = "
                << mesh.time().elapsedCpuTime() << " s." << endl;
        }
        stats.setLayerTime(mesh.time().elapsedCpuTime()-startTime);

#if !defined( WIN32 ) && !defined( WIN64 )
        profiling::writeNow();
#endif
    }
    else if
    (
        controller.algorithm() == meshControl::DUAL
        || controller.algorithm() == meshControl::EXTRUDE
    )
    {
        //If not adding layers perform final optimisation
        dictionary meshOptimDict = finalOptimizeSetup(meshDict);
        word optimType = meshOptimDict.lookupOrDefault<word>
        (
            "type",
            "cfMeshOptimize"
        );
        autoPtr<autoOptimize> optimMeshPtr
            = autoOptimize::New(mesh, meshOptimDict);
        optimMeshPtr->optimize();
    }

    if (meshDict.found("cellRemoval"))
    {
        if
        (
            controller.algorithm() == meshControl::STANDARD
            || controller.algorithm() == meshControl::SHELL
        )
        {
            meshRefiner.updateMasterRegions();
        }

        dictionary cellRemovalDict = meshDict.subDict("cellRemoval");
        errorCellRemoval errorRemover
        (
            meshRefiner,
            cellRemovalDict
        );
        //remove cells if needed but don't update surface intersections
        errorRemover.remove(false);
    }

    //optional patch extrusion methods to final mesh
    if (meshDict.found("extrudeDict"))
    {
        const dictionary& extrudeDict = meshDict.subDict("extrudeDict");
        autoExtrude patchExtrude(meshRefiner,decomposer,distributor);
        patchExtrude.extrudeSelected(extrudeDict);
    }

    if (Pstream::parRun())
    {
        const word initialDecomp(decomposeDict.lookup("method"));
        bool constraints = false;
        if (meshDict.found("constraints"))
        {
            constraints = true;
        }

        if
        (
            (addedLayers || finalDecomp != initialDecomp || constraints)
            && (controller.mode() != meshControl::DRYRUN)
        )
        {
            decomposeDict.set("method",finalDecomp);

            if (constraints)
            {
                const dictionary& constraintDict =
                    meshDict.subDict("constraints");
                decomposeDict.set("constraints",constraintDict);
                Info<<"Adding decomposition constraints: "<<nl
                    << constraintDict
                    <<endl;
            }

            Info<< nl
                << "Doing final balancing" << nl
                << "---------------------" << nl
                << endl;

             // Decomposition
            autoPtr<decompositionMethod> decomposerFinalPtr
            (
                decompositionMethod::New
                (
                    decomposeDict
                 )
            );
            decompositionMethod& decomposerFinal = decomposerFinalPtr();

            if (!decomposerFinal.parallelAware())
            {
                FatalErrorInFunction
                    << "You have selected decomposition method "
                    << decomposerFinal.typeName
                    << " which is not parallel aware." << endl
                    << "Please select one that is (hierarchical, ptscotch)"
                    << exit(FatalError);
            }

            // Balance. No restriction on face zones and baffles.
            meshRefiner.balance
            (
                constraints,
                false,
                false,
                scalarField(mesh.nCells(), 1.0),
                decomposerFinal,
                distributor,
                false
            );
        }
    }

    if (meshDict.found("repatchRegions"))
    {
        const dictionary& regionDict = meshDict.subDict("repatchRegions");
        meshRefiner.repatchRegion(regionDict);
    }

    if (meshDict.found("surfaceToPatch"))
    {
        const dictionary& surfaceToPatchDict =
            meshDict.subDict("surfaceToPatch");
        meshRefiner.surfaceToPatch(surfaceToPatchDict);
    }

    if (meshDict.found("triangulateDict"))
    {
        const dictionary& triangulateDict = meshDict.subDict("triangulateDict");
        meshRefiner.triangulateProblemBoundaryFaces(triangulateDict);
    }

    if (meshDict.found("sourceTargetCheckDict"))
    {
        const dictionary& patchCheckDict =
           meshDict.subDict("sourceTargetCheckDict");
        meshRefiner.sourceTargetChecking(mesh,patchCheckDict);
    }

    meshRefiner.removeUnusedPoints(mesh);

    //Use built in re-numbering for final write
    polyTopoChange meshMod(mesh);
    // Change the mesh.
    autoPtr<mapPolyMesh> map = meshMod.changeMesh
    (
        mesh,
        false,      // inflate
        true,       // parallel sync
        true,       // cell ordering
        true        // point ordering
    );

    // Update fields
    mesh.updateMesh(map);

    // Move mesh (since morphing might not do this)
    if (map().hasMotionPoints())
    {
        mesh.movePoints(map().preMotionPoints());
    }
    else
    {
        // Delete mesh volumes. No other way to do this?
        mesh.clearOut();
    }
    meshRefiner.updateMesh(map, labelList(0));

    if
    (
       meshWrite &&
       ( writeIntermediate || controller.mode() != meshControl::DRYRUN )
    )
    {
        if (debugLevel || writeIntermediate)
        {
            const_cast<Time&>(mesh.time())++;
        }
        else if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }
        else if (writeConstant)
        {
            Time& runTime = const_cast<Time&>(mesh.time());
            runTime.setTime(instant(runTime.constant()), 0);
            mesh.setInstance("constant");
        }
        else
        {
            Time& runTime = const_cast<Time&>(mesh.time());
            runTime.setTime(writeTime, 0);
        }

        writeMesh
        (
            "Final mesh",
            meshRefiner,
            debugLevel,
            meshRefinement::writeLevel(),
            false
        );
    }

    bool writeTriSurf = meshDict.lookupOrDefault<Switch>
    (
       "writeTriSurf",
       false
    );

    if (writeTriSurf)
    {
#if !defined( WIN32 ) && !defined( WIN64 )
        addProfiling(surfaceSimplify, "snappyHexMesh::writeTriSurf");
#endif
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        labelHashSet includePatches(pbm.size());
        if (meshDict.found("writeTriPatches"))
        {
            includePatches = pbm.patchSet
            (
                wordReList(meshDict.lookup("writeTriPatches"))
            );
        }
        else
        {
            forAll(pbm, patchi)
            {
                const polyPatch& patch = pbm[patchi];

                if (!isA<processorPolyPatch>(patch))
                {
                    includePatches.insert(patchi);
                }
            }
        }

        fileName triOutFileName
        (
            meshDict.lookupOrDefault<fileName>
            (
               "writeTriFile",
               "constant/triSurface/simplifiedSurface.stl"
            )
        );

        extractTriSurface
        (
            mesh,
            includePatches,
            triOutFileName
        );
    }

#if !defined( WIN32 ) && !defined( WIN64 )
    profiling::writeNow();
#endif

    Info<< "Finished meshing in = "
        << runTime_.elapsedCpuTime() << " s." << endl;

    if (controller.mode() != meshControl::DRYRUN)
    {
        //Write out mesh summary report
        stats.write();
    }

    return meshPtr;
}

// ************************************************************************* //
