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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2014 blueCAPE Lda
    (c) 2015-2016 OpenCFD Ltd.
    (c) 2017-2023 Esi Ltd.

Modifications
    This file has been modified by blueCAPE's unofficial mingw patches for
    OpenFOAM.



    Modifications made:
      - Derived from the patches for blueCFD 2.1 and 2.2.
      - Changed how the tolerance is handled in 32bit Windows.

Application
    redistributePar

Group
    grpParallelUtilities

Description
    Redistributes existing decomposed mesh and fields according to the current
    settings in the decomposeParDict file.

    Must be run on maximum number of source and destination processors.
    Balances mesh and writes new mesh to new time directory.

Usage

    - redistributePar [OPTION]

    \param -region regionName \n
    Distribute named region.

    \param -allRegions \n
    Decompose all regions in regionProperties. Does not check for
    existence of processor*. Determines decomposition on merged regions
    if specified by consolidatedGroups in regionProperties.

    \param -decompose \n
    Remove any existing \a processor subdirectories and decomposes the
    geometry. Equivalent to running without processor subdirectories.

    \param -reconstruct \n
    Reconstruct geometry (like reconstructParMesh). Equivalent to setting
    numberOfSubdomains 1 in the decomposeParDict. Deprecated in FOAM v4.2.

    \param -newTimes \n
    (in combination with -reconstruct) reconstruct only new times.

    \param -cellDist \n
    Write the cell distribution as a labelList, for use with 'manual'
    decomposition method or as a volScalarField for post-processing.

    \param -cellSets \n
    Update any cellSets if present in the mesh.

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/Time/Time.H"
#include "fvMesh/fvMesh.H"
#include "fvMeshTools/fvMeshTools.H"
#include "fvMeshDistribute/fvMeshDistribute.H"
#include "decompositionMethod/decompositionMethod.H"
#include "decompositionModel.H"
#include "db/Time/timeSelector.H"
#include "db/IOstreams/Pstreams/PstreamReduceOps.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMeshDistribute/IOmapDistributePolyMesh.H"
#include "db/IOobjectList/IOobjectList.H"
#include "meshes/polyMesh/globalMeshData/globalIndex.H"
#include "loadOrCreateMesh.H"
#include "fields/fvPatchFields/constraint/processor/processorFvPatchField.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "decompositionModel.H"
#include "sets/topoSets/topoSet.H"

#include "parFvFieldReconstructor.H"
#include "parLagrangianRedistributor.H"
#include "unmappedPassiveParticleCloud.H"
#include "polyTopoChange/hexRef8/hexRef8Data.H"
#include "polyTopoChange/hexRef8/hexRef8DataList.H"
#include "meshRefinement/meshRefinement.H"
#include "regionProperties/regionProperties.H"
#include "fields/GeometricFields/pointFields/pointFields.H"
#include "sets/topoSets/cellSet.H"
#include "polyMeshAdder/polyMeshAdder.H"
#include "polyMeshAdder/faceCoupleInfo.H"
#include "primitives/ops/flipOp.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Tolerance (as fraction of the bounding box). Needs to be fairly lax since
// usually meshes get written with limited precision (6 digits)
static const scalar defaultMergeTol = 1e-6;

typedef IOList<bool> boolIOList;
defineTemplateTypeName(boolIOList);

// Get merging distance when matching face centres
scalar getMergeDistance
(
    const argList& args,
    const Time& runTime,
    const boundBox& bb
)
{
    scalar mergeTol = defaultMergeTol;
    args.optionReadIfPresent("mergeTol", mergeTol);

    scalar writeTol =
        Foam::pow(scalar(10.0), -scalar(IOstream::defaultPrecision()));

    Info<< "Merge tolerance : " << mergeTol << nl
        << "Write tolerance : " << writeTol << endl;

    #ifdef WIN32
    //this was necessary due to some crazy bug...
    writeTol -= SMALL;
    Info<< "Write tol. adj. : " << writeTol << endl;
    #endif

    if (runTime.writeFormat() == IOstream::ASCII && mergeTol < writeTol)
    {
        FatalErrorInFunction
            << "Your current settings specify ASCII writing with "
            << IOstream::defaultPrecision() << " digits precision." << endl
            << "Your merging tolerance (" << mergeTol << ") is finer than this."
            << endl
            << "Please change your writeFormat to binary"
            << " or increase the writePrecision" << endl
            << "or adjust the merge tolerance (-mergeTol)."
            << exit(FatalError);
    }

    scalar mergeDist = mergeTol * bb.mag();

    Info<< "Overall meshes bounding box : " << bb << nl
        << "Relative tolerance          : " << mergeTol << nl
        << "Absolute matching distance  : " << mergeDist << nl
        << endl;

    return mergeDist;
}


void printMeshData(PtrList<fvMesh>& meshes)
{
    labelListList localSizes(Pstream::nProcs(), labelList(meshes.size()));
    labelListList localBoundFaces(Pstream::nProcs(), labelList(meshes.size()));

    labelListList myProcPatchNeiProcNo(meshes.size());
    labelListList myProcPatchSize(meshes.size());

    forAll(meshes, i)
    {
        fvMesh& mesh = meshes[i];
        localSizes[Pstream::myProcNo()][i] = mesh.nCells();

        localBoundFaces[Pstream::myProcNo()][i] =
            mesh.nFaces() - mesh.nInternalFaces();

        labelList processorNeighbours(mesh.boundaryMesh().size());
        labelList processorSize(mesh.boundaryMesh().size());

        label nNeighbours = 0;

        forAll(mesh.boundaryMesh(), patchi)
        {
            if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchi]))
            {
                const processorPolyPatch& ppp =
                    refCast<const processorPolyPatch>
                    (
                        mesh.boundaryMesh()[patchi]
                    );
                processorNeighbours[nNeighbours] = ppp.neighbProcNo();
                processorSize[nNeighbours++] = ppp.size();
            }
        }
        processorNeighbours.setSize(nNeighbours);
        processorSize.setSize(nNeighbours);

        myProcPatchNeiProcNo[i] = processorNeighbours;
        myProcPatchSize[i] = processorSize;
    }

    Pstream::gatherList(localSizes);
    Pstream::gatherList(localBoundFaces);

    labelListListList patchNeiProcNos(Pstream::nProcs());
    labelListListList patchSizes     (Pstream::nProcs());
    patchNeiProcNos[Pstream::myProcNo()] = myProcPatchNeiProcNo;
    patchSizes     [Pstream::myProcNo()] = myProcPatchSize;
    Pstream::gatherList(patchNeiProcNos);
    Pstream::gatherList(patchSizes);

    label nFinalProcs = meshes.size();
    reduce(nFinalProcs, sumOp<label>());

    if (!Pstream::master()) return;

    labelList nMeshesInProc(Pstream::nProcs());
    labelList meshToProc(identity(nFinalProcs));
    labelListList meshNoInProc =
        fvMeshDistribute::calcMeshToProcMap
        (
            nFinalProcs,
            nMeshesInProc,
            meshToProc
        );

    label         globalCells = 0;
    labelList     localSize          (nFinalProcs);
    labelList     globalBoundaryFaces(nFinalProcs);
    labelListList patchNeiProcNo     (nFinalProcs);
    labelListList patchSize          (nFinalProcs);
    forAll(localSizes, proci)
    {
        forAll(localSizes[proci], j)
        {
            label meshi = meshNoInProc[proci][j];

            globalCells += localSizes[proci][j];

            localSize[meshi] = localSizes[proci][j];

            globalBoundaryFaces[meshi] = localBoundFaces[proci][j];

            patchNeiProcNo[meshi] = patchNeiProcNos[proci][j];

            patchSize[meshi] = patchSizes[proci][j];
        }
    }

    FOAM_ASSERT (globalCells >= 0) {
        FatalErrorInFunction
            << "Overflow : sum of sizes " << globalCells
            << " exceeds capability of label (" << labelMax
            << "). Please recompile with larger datatype for label."
            << exit(FatalError);
    }

    // Print stats

    label maxProcCells = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;
    label maxProcFaces = 0;

    for (label procI = 0; procI < nFinalProcs; procI++)
    {
        Info<< endl
            << "Processor " << procI << nl
            << "    Number of cells = " << localSize[procI]
            << endl;

        label nProcFaces = 0;

        const labelList& nei = patchNeiProcNo[procI];

        forAll(patchNeiProcNo[procI], i)
        {
            Info<< "    Number of faces shared with processor "
                << patchNeiProcNo[procI][i] << " = " << patchSize[procI][i]
                << endl;

            nProcFaces += patchSize[procI][i];
        }

        Info<< "    Number of processor patches = " << nei.size() << nl
            << "    Number of processor faces = " << nProcFaces << nl
            << "    Number of boundary faces = "
            << globalBoundaryFaces[procI]-nProcFaces
            << endl;

        maxProcCells = max(maxProcCells, localSize[procI]);
        totProcFaces += nProcFaces;
        totProcPatches += nei.size();
        maxProcPatches = max(maxProcPatches, nei.size());
        maxProcFaces = max(maxProcFaces, nProcFaces);
    }

    // Stats

    // In case of all faces on one processor, use max to avoid division by 0.
    scalar avgProcCells   = max(1, scalar(globalCells)/nFinalProcs);
    scalar avgProcPatches = max(1, scalar(totProcPatches)/nFinalProcs);
    scalar avgProcFaces   = max(1, scalar(totProcFaces)/nFinalProcs);

    if (avgProcCells > 1)
    {
        Info<< nl
            << "Number of processor faces = " << totProcFaces/2 << nl
            << "Max number of cells = " << maxProcCells
            << " (" << 100.0*(maxProcCells-avgProcCells)/avgProcCells
            << "% above average " << avgProcCells << ")" << nl
            << "Max number of processor patches = " << maxProcPatches
            << " (" << 100.0*(maxProcPatches-avgProcPatches)/avgProcPatches
            << "% above average " << avgProcPatches << ")" << nl
            << "Max number of faces between processors = " << maxProcFaces
            << " (" << 100.0*(maxProcFaces-avgProcFaces)/avgProcFaces
            << "% above average " << avgProcFaces << ")" << nl
            << endl;
    }
}

void printMeshData(const polyMesh& mesh)
{
    // Collect all data on master

    globalIndex globalCells(mesh.nCells());
    labelListList patchNeiProcNo(Pstream::nProcs());
    labelListList patchSize(Pstream::nProcs());
    const labelList& pPatches = mesh.globalData().processorPatches();
    patchNeiProcNo[Pstream::myProcNo()].setSize(pPatches.size());
    patchSize[Pstream::myProcNo()].setSize(pPatches.size());
    forAll(pPatches, i)
    {
        const processorPolyPatch& ppp = refCast<const processorPolyPatch>
        (
            mesh.boundaryMesh()[pPatches[i]]
        );
        patchNeiProcNo[Pstream::myProcNo()][i] = ppp.neighbProcNo();
        patchSize[Pstream::myProcNo()][i] = ppp.size();
    }
    Pstream::gatherList(patchNeiProcNo);
    Pstream::gatherList(patchSize);


    // Print stats

    globalIndex globalBoundaryFaces(mesh.nFaces()-mesh.nInternalFaces());

    label maxProcCells = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;
    label maxProcFaces = 0;

    for (label procI = 0; procI < Pstream::nProcs(); procI++)
    {
        Info<< endl
            << "Processor " << procI << nl
            << "    Number of cells = " << globalCells.localSize(procI)
            << endl;

        label nProcFaces = 0;

        const labelList& nei = patchNeiProcNo[procI];

        forAll(patchNeiProcNo[procI], i)
        {
            Info<< "    Number of faces shared with processor "
                << patchNeiProcNo[procI][i] << " = " << patchSize[procI][i]
                << endl;

            nProcFaces += patchSize[procI][i];
        }

        Info<< "    Number of processor patches = " << nei.size() << nl
            << "    Number of processor faces = " << nProcFaces << nl
            << "    Number of boundary faces = "
            << globalBoundaryFaces.localSize(procI)-nProcFaces << endl;

        maxProcCells = max(maxProcCells, globalCells.localSize(procI));
        totProcFaces += nProcFaces;
        totProcPatches += nei.size();
        maxProcPatches = max(maxProcPatches, nei.size());
        maxProcFaces = max(maxProcFaces, nProcFaces);
    }

    // Stats

    scalar avgProcCells = scalar(globalCells.size())/Pstream::nProcs();
    scalar avgProcPatches = scalar(totProcPatches)/Pstream::nProcs();
    scalar avgProcFaces = scalar(totProcFaces)/Pstream::nProcs();

    // In case of all faces on one processor. Just to avoid division by 0.
    if (totProcPatches == 0)
    {
        avgProcPatches = 1;
    }
    if (totProcFaces == 0)
    {
        avgProcFaces = 1;
    }

    Info<< nl
        << "Number of processor faces = " << totProcFaces/2 << nl
        << "Max number of cells = " << maxProcCells
        << " (" << 100.0*(maxProcCells-avgProcCells)/avgProcCells
        << "% above average " << avgProcCells << ")" << nl
        << "Max number of processor patches = " << maxProcPatches
        << " (" << 100.0*(maxProcPatches-avgProcPatches)/avgProcPatches
        << "% above average " << avgProcPatches << ")" << nl
        << "Max number of faces between processors = " << maxProcFaces
        << " (" << 100.0*(maxProcFaces-avgProcFaces)/avgProcFaces
        << "% above average " << avgProcFaces << ")" << nl
        << endl;
}


// Debugging: write volScalarField with decomposition for post processing.
void writeDecomposition
(
    const word& name,
    const fvMesh& mesh,
    const labelList& decomp
)
{
    // Write the decomposition as labelList for use with 'manual'
    // decomposition method.
    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            mesh.facesInstance(),  // mesh read from facesInstance
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        decomp
    );
    cellDecomposition.write();

    Info<< "Writing wanted cell distribution to volScalarField " << name
        << " for postprocessing purposes." << nl << endl;

    volScalarField procCells
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false                   // do not register
        ),
        mesh,
        dimensionedScalar(name, dimless, -1),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(procCells, cI)
    {
        procCells[cI] = decomp[cI];
    }

    procCells.correctBoundaryConditions();
    procCells.write();
}


void determineDecomposition
(
    const Time& baseRunTime,
    const fileName& decompDictFile, // location for decomposeParDict to use
    const bool decompose,       // decompose, i.e. read from undecomposed case
    const fileName& proc0CaseName,
    const UPtrList<fvMesh>& meshes,
    const bool writeCellDist,
    const scalar& dimTol,

    label& nDestProcs,
    UPtrList<labelList>& decomp
)
{

    // Read decomposeParDict (on all processors)
    const decompositionModel& method = decompositionModel::New
    (
        meshes[0],
        decompDictFile
    );

    decompositionMethod& decomposer = method.decomposer();

    if (!decomposer.parallelAware())
    {
        WarningInFunction
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which does" << endl
            << "not synchronise the decomposition across"
            << " processor patches." << endl
            << "    You might want to select a decomposition method"
            << " which is aware of this. Continuing."
            << endl;
    }

    if (Pstream::master() && decompose)
    {
        Info<< "Setting caseName to " << baseRunTime.caseName()
            << " to read decomposeParDict" << endl;
        const_cast<Time&>(meshes[0].time()).TimePaths::caseName() =
            baseRunTime.caseName();
    }

    // Read weights from all meshes
    PtrList<scalarField> cellWeights;
    if (method.found("weightField"))
    {
        cellWeights.resize(meshes.size());

        word weightName = method.lookup("weightField");

        forAll(meshes, i)
        {
            volScalarField weights
            (
                IOobject
                (
                    weightName,
                    meshes[i].time().timeName(),
                    meshes[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                meshes[i]
            );
            cellWeights.set(i, scalarField(weights.internalField()));
        }
    }

    nDestProcs = decomposer.nDomains();

    // Merge meshes for decomposer
    decomp.resize(meshes.size());
    if (meshes.size() == 1)
    {
        decomp[0] =
            decomposer.decompose
            (
                meshes[0],
                cellWeights.size() ? cellWeights[0] : scalarField()
            );
    }
    else
    {
        PtrList<labelList> separateToMergedCellMap(meshes.size());

        // Using polyMeshAdder - we could alternatively use mergePolyMesh

        // Start by creating blank mesh
        autoPtr<polyMesh> mergedMeshPtr
        (
            new polyMesh
            (
                IOobject
                (
                    "mergedMesh",
                    baseRunTime.constant(),
                    baseRunTime
                ),
                xferCopy(pointField(0)),
                xferCopy(faceList(0)),
                xferCopy(cellList(0))
            )
        );

        // Add meshes to it
        for (label i = 0; i < meshes.size(); i++)
        {
            autoPtr<mapAddedPolyMesh> mapPtr =
                polyMeshAdder::add
                (
                    mergedMeshPtr(),
                    meshes[i],
                    // Don't try to match and merge points
                    faceCoupleInfo
                    (
                        mergedMeshPtr(),
                        labelList(0),
                        meshes[i],
                        labelList(0),
                        dimTol,
                        true,
                        true,
                        false
                    )
                );

            // Modify the existing maps
            for (label j = 0; j < i; j++)
            {
                const labelUList& map(separateToMergedCellMap[j]);
                separateToMergedCellMap[j] =
                    labelList(mapPtr->oldCellMap(), map);
            }

            // Add new map
            separateToMergedCellMap.set(i, new labelList(mapPtr->addedCellMap()));
        }

        // Create merged cell weights
        scalarField mergedCellWeights;
        if (cellWeights.size())
        {
            scalarField mergedCellWeights(mergedMeshPtr->nCells());
            forAll(cellWeights, i)
            {
                forAll(cellWeights[i], cellj)
                {
                    mergedCellWeights[separateToMergedCellMap[i][cellj]] =
                        cellWeights[i][cellj];
                }
            }
        }

        // Call decompose on merged mesh
        labelList mergedDecomp =
            decomposer.decompose(mergedMeshPtr(), mergedCellWeights);

        // Redistribute the decomposition
        forAll(decomp, i)
        {
            const labelUList& map(separateToMergedCellMap[i]);
            decomp[i] = labelList(mergedDecomp, map);
        }
    }

    if (Pstream::master() && decompose)
    {
        Info<< "Restoring caseName to " << proc0CaseName << endl;
        const_cast<Time&>(meshes[0].time()).TimePaths::caseName() =
            proc0CaseName;
    }

    // Dump decomposition to volScalarField
    if (writeCellDist)
    {
        // Note: on master make sure to write to processor0
        if (decompose)
        {
            if (Pstream::master())
            {
                Info<< "Setting caseName to " << baseRunTime.caseName()
                    << " to write undecomposed cellDist" << endl;

                Time& tm = const_cast<Time&>(meshes[0].time());
                tm.TimePaths::caseName() = baseRunTime.caseName();

                forAll(meshes, i)
                {
                    writeDecomposition("cellDist", meshes[i], decomp[i]);
                }

                Info<< "Restoring caseName to " << proc0CaseName << endl;
                tm.TimePaths::caseName() = proc0CaseName;
            }
        }
        else
        {
            forAll(meshes, i)
            {
                writeDecomposition("cellDist", meshes[i], decomp[i]);
            }
        }
    }
}


void determineDecomposition
(
    const Time& baseRunTime,
    const fileName& decompDictFile, // optional location for decomposeParDict
    const bool decompose,       // decompose, i.e. read from undecomposed case
    const fileName& proc0CaseName,
    fvMesh& mesh,
    const bool writeCellDist,
    const scalar& dimTol,

    label& nDestProcs,
    labelList& decomp
)
{
    UPtrList<fvMesh> meshes(1);
    meshes.set(0, &mesh);
    UPtrList<labelList> decomps(1);
    decomps.set(0, &decomp);

    determineDecomposition
    (
        baseRunTime,
        decompDictFile,
        decompose,
        proc0CaseName,
        meshes,
        writeCellDist,
        dimTol,
        nDestProcs,
        decomps
    );
}


// Write addressing if decomposing (1 to many) or reconstructing (many to 1)
void writeProcAddressing
(
    const bool decompose,
    const fileName& meshSubDir,
    const fvMesh& mesh,
    const mapDistributePolyMesh& map
)
{
    Info<< "Writing procAddressing files to " << mesh.facesInstance()
        << endl;

    labelIOList cellMap
    (
        IOobject
        (
            "cellProcAddressing",
            mesh.facesInstance(),
            meshSubDir,
            mesh,
            IOobject::NO_READ
        ),
        0
    );

    labelIOList faceMap
    (
        IOobject
        (
            "faceProcAddressing",
            mesh.facesInstance(),
            meshSubDir,
            mesh,
            IOobject::NO_READ
        ),
        0
    );

    labelIOList pointMap
    (
        IOobject
        (
            "pointProcAddressing",
            mesh.facesInstance(),
            meshSubDir,
            mesh,
            IOobject::NO_READ
        ),
        0
    );

    labelIOList patchMap
    (
        IOobject
        (
            "boundaryProcAddressing",
            mesh.facesInstance(),
            meshSubDir,
            mesh,
            IOobject::NO_READ
        ),
        0
    );

    // Decomposing: see how cells moved from undecomposed case
    if (decompose)
    {
        cellMap = identity(map.nOldCells());
        map.distributeCellData(cellMap);

        faceMap = identity(map.nOldFaces());
        {
            const mapDistribute& faceDistMap = map.faceMap();

            if (faceDistMap.subHasFlip() || faceDistMap.constructHasFlip())
            {
                // Offset by 1
                faceMap = faceMap + 1;
            }
            // Apply face flips
            mapDistributeBase::distribute
            (
                Pstream::commsTypes::nonBlocking,
                List<labelPair>(),
                faceDistMap.constructSize(),
                faceDistMap.subMap(),
                faceDistMap.subHasFlip(),
                faceDistMap.constructMap(),
                faceDistMap.constructHasFlip(),
                faceMap,
                flipLabelOp()
            );
        }

        pointMap = identity(map.nOldPoints());
        map.distributePointData(pointMap);

        patchMap = identity(map.oldPatchSizes().size());
        const mapDistribute& patchDistMap = map.patchMap();
        // Use explicit distribute since we need to provide a null value
        // (for new patches) and this is the only call that allow us to
        // provide one ...
        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            patchDistMap.constructSize(),
            patchDistMap.subMap(),
            patchDistMap.subHasFlip(),
            patchDistMap.constructMap(),
            patchDistMap.constructHasFlip(),
            patchMap,
            eqOp<label>(),
            flipOp(),
            label(-1),
            UPstream::msgType()
        );
    }
    else    // if (nDestProcs == 1)
    {
        cellMap = identity(mesh.nCells());
        map.cellMap().reverseDistribute(map.nOldCells(), cellMap);

        faceMap = identity(mesh.nFaces());
        {
            const mapDistribute& faceDistMap = map.faceMap();

            if (faceDistMap.subHasFlip() || faceDistMap.constructHasFlip())
            {
                // Offset by 1
                faceMap = faceMap + 1;
            }

            mapDistributeBase::distribute
            (
                Pstream::commsTypes::nonBlocking,
                List<labelPair>(),
                map.nOldFaces(),
                faceDistMap.constructMap(),
                faceDistMap.constructHasFlip(),
                faceDistMap.subMap(),
                faceDistMap.subHasFlip(),
                faceMap,
                flipLabelOp()
            );
        }

        pointMap = identity(mesh.nPoints());
        map.pointMap().reverseDistribute(map.nOldPoints(), pointMap);

        const mapDistribute& patchDistMap = map.patchMap();
        patchMap = identity(mesh.boundaryMesh().size());
        patchDistMap.reverseDistribute
        (
            map.oldPatchSizes().size(),
            label(-1),
            patchMap
        );
    }

    bool cellOk = cellMap.write();
    bool faceOk = faceMap.write();
    bool pointOk = pointMap.write();
    bool patchOk = patchMap.write();

    if (!cellOk || !faceOk || !pointOk || !patchOk)
    {
        WarningInFunction
            << "Failed to write " << cellMap.objectPath()
            << ", " << faceMap.objectPath()
            << ", " << pointMap.objectPath()
            << ", " << patchMap.objectPath()
            << endl;
    }
}

template<class GeoField>
void correctCoupledBoundaryConditions(PtrList<fvMesh>& myProcMeshes)
{
    forAll(myProcMeshes, i)
    {
        fvMesh& mesh = myProcMeshes[i];

        HashTable<GeoField*> flds
        (
            mesh.objectRegistry::lookupClass<GeoField>()
        );

        forAllIter(typename HashTable<GeoField*>, flds, iter)
        {
            GeoField& fld = *iter();

            typename GeoField::Boundary& bfld = fld.boundaryFieldRef();

            forAll(bfld, patchi)
            {
                typename GeoField::Patch& pfld = bfld[patchi];

                if (pfld.patch().coupled())
                {
                    if (isA<cyclicAMIFvPatch>(pfld.patch()))
                    {
                        // Explicitly assign cell value to patch value
                        // (i.e. zeroGradient)
                        auto& f = pfld.internalField();

                        const labelUList& faceCells = pfld.patch().faceCells();

                        forAll(pfld, facei)
                        {
                            pfld[facei] = f[faceCells[facei]];
                        }
                    }
                    else
                    {
                        pfld = pfld.patchInternalField();
                    }
                }
            }
        }
    } // forAll myProcMeshes
}


// Variant of GeometricField::correctBoundaryConditions that only
// evaluates selected patch fields
template<class GeoField, class CoupledPatchType>
void correctCoupledBoundaryConditions(fvMesh& mesh)
{
    HashTable<GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllIter(typename HashTable<GeoField*>, flds, iter)
    {
        GeoField& fld = *iter();

        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();
        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            label nReq = Pstream::nRequests();

            forAll(bfld, patchi)
            {
                typename GeoField::Patch& pfld = bfld[patchi];

                if (isA<nonConformalFvPatch>(pfld.patch())) continue;

                if (pfld.patch().coupled())
                {
                    pfld.initEvaluate(Pstream::defaultCommsType);
                }
            }

            // Block for any outstanding requests
            if
            (
                Pstream::parRun()
             && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
            )
            {
                Pstream::waitRequests(nReq);
            }

            forAll(bfld, patchi)
            {
                typename GeoField::Patch& pfld = bfld[patchi];

                if (isA<nonConformalFvPatch>(pfld.patch())) continue;

                if (pfld.patch().coupled())
                {
                    pfld.evaluate(Pstream::defaultCommsType);
                }
            }
        }
        else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
        {
            const lduSchedule& patchSchedule =
                fld.mesh().globalData().patchSchedule();

            forAll(patchSchedule, patchEvali)
            {
                label patchi = patchSchedule[patchEvali].patch;
                typename GeoField::Patch& pfld = bfld[patchi];

                if (isA<nonConformalFvPatch>(pfld.patch())) continue;

                if (pfld.patch().coupled())
                {
                    if (patchSchedule[patchEvali].init)
                    {
                        pfld.initEvaluate(Pstream::commsTypes::scheduled);
                    }
                    else
                    {
                        pfld.evaluate(Pstream::commsTypes::scheduled);
                    }
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unsuported communications type "
                << Pstream::commsTypeNames[Pstream::defaultCommsType]
                << exit(FatalError);
        }
    }
}


// Inplace redistribute mesh and any fields
autoPtr<mapDistributePolyMesh> redistributeAndWrite
(
    const Time& baseRunTime,
    const scalar tolDim,
    const boolList& haveMesh,
    const fileName& meshSubDir,
    const bool doReadFields,
    const bool decompose,       // decompose, i.e. read from undecomposed case
    const bool overwrite,
    const bool updateCellSets,
    const fileName& proc0CaseName,
    const label nDestProcs,
    const labelList& decomp,
    const fileName& masterInstDir,
    fvMesh& mesh
)
{
    Time& runTime = const_cast<Time&>(mesh.time());

    //// Print some statistics
    //Info<< "Before distribution:" << endl;
    //printMeshData(mesh);


    PtrList<volScalarField> volScalarFields;
    PtrList<volVectorField> volVectorFields;
    PtrList<volSphericalTensorField> volSphereTensorFields;
    PtrList<volSymmTensorField> volSymmTensorFields;
    PtrList<volTensorField> volTensorFields;

    PtrList<surfaceScalarField> surfScalarFields;
    PtrList<surfaceVectorField> surfVectorFields;
    PtrList<surfaceSphericalTensorField> surfSphereTensorFields;
    PtrList<surfaceSymmTensorField> surfSymmTensorFields;
    PtrList<surfaceTensorField> surfTensorFields;

    PtrList<DimensionedField<scalar, volMesh>> dimScalarFields;
    PtrList<DimensionedField<vector, volMesh>> dimVectorFields;
    PtrList<DimensionedField<sphericalTensor, volMesh>> dimSphereTensorFields;
    PtrList<DimensionedField<symmTensor, volMesh>> dimSymmTensorFields;
    PtrList<DimensionedField<tensor, volMesh>> dimTensorFields;

    DynamicList<word> pointFieldNames;


    if (doReadFields)
    {
        // Create 0 sized mesh to do all the generation of zero sized
        // fields on processors that have zero sized meshes. Note that this is
        // only necessary on master but since polyMesh construction with
        // Pstream::parRun does parallel comms we have to do it on all
        // processors
        autoPtr<fvMeshSubset> subsetterPtr;

        const bool allHaveMesh = (findIndex(haveMesh, false) == -1);
        if (!allHaveMesh)
        {
            // Find last non-processor patch.
            const polyBoundaryMesh& patches = mesh.boundaryMesh();

            label nonProcI = -1;

            forAll(patches, patchI)
            {
                if (isA<processorPolyPatch>(patches[patchI]))
                {
                    break;
                }
                nonProcI++;
            }

            if (nonProcI == -1)
            {
                FatalErrorInFunction
                    << "Cannot find non-processor patch on processor "
                    << Pstream::myProcNo() << endl
                    << " Current patches:" << patches.names()
                    << abort(FatalError);
            }

            // Subset 0 cells, no parallel comms. This is used to create
            // zero-sized fields.
            subsetterPtr.reset(new fvMeshSubset(mesh));
            subsetterPtr().setLargeCellSubset(labelHashSet(0), nonProcI, false);
        }


        // Get original objects (before incrementing time!)
        if (Pstream::master() && decompose)
        {
            runTime.TimePaths::caseName() = baseRunTime.caseName();
        }
        IOobjectList objects(mesh, runTime.timeName());
        if (Pstream::master() && decompose)
        {
            runTime.TimePaths::caseName() = proc0CaseName;
        }

        Info<< "From time " << runTime.timeName()
            << " have objects:" << objects.names() << endl;

        // We don't want to map the decomposition (mapping already tested when
        // mapping the cell centre field)
        IOobjectList::iterator iter = objects.find("cellDist");
        if (iter != objects.end())
        {
            objects.erase(iter);
        }


        // volFields

        if (Pstream::master() && decompose)
        {
            runTime.TimePaths::caseName() = baseRunTime.caseName();
        }
        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volScalarFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volVectorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volSphereTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volSymmTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volTensorFields
        );


        // surfaceFields

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfScalarFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfVectorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfSphereTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfSymmTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfTensorFields
        );


        // Dimensioned internal fields
        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimScalarFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimVectorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimSphereTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimSymmTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimTensorFields
        );

        // pointFields currently not supported. Read their names so we
        // can delete them.
        {
            // Get my objects of type
            pointFieldNames.append
            (
                objects.lookupClass(pointScalarField::typeName).sortedNames()
            );
            pointFieldNames.append
            (
                objects.lookupClass(pointVectorField::typeName).sortedNames()
            );
            pointFieldNames.append
            (
                objects.lookupClass
                (
                    pointSphericalTensorField::typeName
                ).sortedNames()
            );
            pointFieldNames.append
            (
                objects.lookupClass
                (
                    pointSymmTensorField::typeName
                ).sortedNames()
            );
            pointFieldNames.append
            (
                objects.lookupClass(pointTensorField::typeName).sortedNames()
            );

            // Make sure all processors have the same set
            Pstream::scatter(pointFieldNames);
        }

        if (Pstream::master() && decompose)
        {
            runTime.TimePaths::caseName() = proc0CaseName;
        }
    }

    PtrList<cellSet> cellSets(0);
    if (updateCellSets)
    {
        // Read sets
        if (Pstream::master() && decompose)
        {
            runTime.TimePaths::caseName() = baseRunTime.caseName();
        }

        word cSetInstance =
            mesh.time().findInstance
            (
                "polyMesh/sets",
                word::null,
                IOobject::READ_IF_PRESENT
            );
        Pstream::scatter(cSetInstance);

        IOobjectList objects(mesh, cSetInstance, "polyMesh/sets");
        IOobjectList cSets(objects.lookupClass(cellSet::typeName));

        // Check that we all have all objects
        wordList objectNames = cSets.sortedNames();

        Info<< "Found cellSets in folder "
            << cSetInstance
            << ": "
            << objectNames
            << endl;

        // Get master names
        wordList masterNames(objectNames);
        Pstream::scatter(masterNames);
        cellSets.setSize(masterNames.size());

        forAll(masterNames, i)
        {
            word name = masterNames[i];
            if (cSets.found(name))
            {
                cellSets.set
                (
                    i,
                    new cellSet
                    (
                        mesh,
                        name,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    )
                );
            }
            else
            {
                //Create zero size set if not found
                cellSets.set
                (
                    i,
                    new cellSet
                    (
                        mesh,
                        name,
                        0,
                        IOobject::NO_WRITE
                    )
                );
            }
        }

        if (Pstream::master() && decompose)
        {
            runTime.TimePaths::caseName() = proc0CaseName;
        }
    }

    PtrList<Time>                  myProcTimes;
    PtrList<fvMesh>                myProcMeshes;
    PtrList<mapDistributePolyMesh> myProcDist;

    List<PtrList<volScalarField>>              myProcVolScalarFields;
    List<PtrList<volVectorField>>              myProcVolVectorFields;
    List<PtrList<volSphericalTensorField>>     myProcVolSphericalTensorFields;
    List<PtrList<volSymmTensorField>>          myProcVolSymmTensorFields;
    List<PtrList<volTensorField>>              myProcVolTensorFields;
    List<PtrList<volVector1Field>>             myProcVolVector1Fields;
    List<PtrList<volVector4Field>>             myProcVolVector4Fields;
    List<PtrList<volTensor4Field>>             myProcVolTensor4Fields;
    List<PtrList<surfaceScalarField>>          myProcSurfaceScalarFields;
    List<PtrList<surfaceVectorField>>          myProcSurfaceVectorFields;
    List<PtrList<surfaceSphericalTensorField>> myProcSurfaceSphericalTensorFields;
    List<PtrList<surfaceSymmTensorField>>      myProcSurfaceSymmTensorFields;
    List<PtrList<surfaceTensorField>>          myProcSurfaceTensorFields;

    List<PtrList<cellSet>> myProcCellSets;

    // Mesh distribution engine
    fvMeshDistribute distributor(mesh, tolDim);

    if (nDestProcs > Pstream::nProcs())
    {
        distributor.setNFinalProcs(nDestProcs);

        // Initialize empty list of meshes and fields
        // for each processor
        distributor.initMeshList
        (
            meshSubDir,
            nDestProcs,
            mesh,

            myProcTimes,
            myProcMeshes,

            myProcVolScalarFields,
            myProcVolVectorFields
            //myProcVolSphericalTensorFields,
            //myProcVolSymmTensorFields,
            //myProcVolTensorFields,
            //myProcVolVector1Fields,
            //myProcVolVector4Fields,
            //myProcVolTensor4Fields,
            //myProcSurfaceScalarFields,
            //myProcSurfaceVectorFields,
            //myProcSurfaceSphericalTensorFields,
            //myProcSurfaceSymmTensorFields,
            //myProcSurfaceTensorFields
        );

        if (updateCellSets)
        {
            List<PtrList<cellSet>> myProcCellSets0(myProcMeshes.size());

            wordList names(cellSets.size());

            forAll(myProcMeshes, i)
            {
                myProcCellSets0[i].setSize(cellSets.size());

                const fvMesh& meshi = myProcMeshes[i];

                forAll(cellSets, j)
                {
                    if (i == 0)
                    {
                        names[j] = cellSets[j].name(),
                        myProcCellSets0[i].set
                        (
                            j,
                            cellSets.set(j, nullptr)
                        );
                    }
                    else
                    {
                        myProcCellSets0[i].set
                        (
                            j,
                            new cellSet
                            (
                                meshi,
                                names[j],
                                0,
                                IOobject::NO_WRITE
                            )
                        );
                    }
                }
            } // forAll myProcMeshes

            myProcCellSets.transfer(myProcCellSets0);
        } //if updateCellSets
    } // if nDestProc > nProcs

    // Do all the distribution of mesh and fields
    autoPtr<mapDistributePolyMesh> rawMap =
        distributor.distribute
        (
            decomp,
            myProcMeshes,
            myProcDist
        );

    // Print some statistics
    Info<< "\nAfter distribution:" << endl;
    if (nDestProcs > Pstream::nProcs())
    {
        printMeshData(myProcMeshes);

        correctCoupledBoundaryConditions<volScalarField         >(myProcMeshes);
        correctCoupledBoundaryConditions<volVectorField         >(myProcMeshes);
        correctCoupledBoundaryConditions<volSphericalTensorField>(myProcMeshes);
        correctCoupledBoundaryConditions<volSymmTensorField     >(myProcMeshes);
        correctCoupledBoundaryConditions<volTensorField         >(myProcMeshes);
    }
    else
    {
        printMeshData(mesh);

        // Get other side of processor boundaries
        correctCoupledBoundaryConditions
        <
            volScalarField,
            processorFvPatchField<scalar>
        >(mesh);
        correctCoupledBoundaryConditions
        <
            volVectorField,
            processorFvPatchField<vector>
        >(mesh);
        correctCoupledBoundaryConditions
        <
            volSphericalTensorField,
            processorFvPatchField<sphericalTensor>
        >(mesh);
        correctCoupledBoundaryConditions
        <
            volSymmTensorField,
            processorFvPatchField<symmTensor>
        >(mesh);
        correctCoupledBoundaryConditions
        <
            volTensorField,
            processorFvPatchField<tensor>
        >(mesh);
        // No update surface fields
    }


    // Set the minimum write precision
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));


    if (!overwrite)
    {
        runTime++;
        mesh.setInstance(runTime.timeName());
    }
    else
    {
        mesh.setInstance(masterInstDir);
    }

    IOmapDistributePolyMesh map
    (
        IOobject
        (
            "procAddressing",
            mesh.facesInstance(),
            meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );
    map.transfer(rawMap());


    if (nDestProcs == 1)
    {
        if (Pstream::master())
        {
            Info<< "Setting caseName to " << baseRunTime.caseName()
                << " to write reconstructed mesh and fields." << endl;
            runTime.TimePaths::caseName() = baseRunTime.caseName();

            mesh.write();
            forAll(pointFieldNames, i)
            {
                IOobject io
                (
                    pointFieldNames[i],
                    runTime.timeName(),
                    mesh
                );

                fileName fieldFile(io.objectPath());
                if (topoSet::debug) DebugVar(fieldFile);
                rm(fieldFile);
            }

            // Now we've written all. Reset caseName on master
            Info<< "Restoring caseName to " << proc0CaseName << endl;
            runTime.TimePaths::caseName() = proc0CaseName;
        }
    }
    else if (nDestProcs > Pstream::nProcs())
    {
        labelListList meshNoInProc
        (
            Pstream::nProcs(),
            labelList(myProcMeshes.size(), 0)
        );
        //if (fvMeshDistribute::debug)
        {
            labelList nMeshesInProc(Pstream::nProcs());
            labelList meshToProc(identity(nDestProcs));
            meshNoInProc =
                fvMeshDistribute::calcMeshToProcMap
                (
                    nDestProcs,
                    nMeshesInProc,
                    meshToProc
                );
        }

        forAll(myProcMeshes, i)
        {
            fvMesh& meshi = myProcMeshes[i];

            Time& runTimei = const_cast<Time&>(meshi.time());

            runTimei.setTime(runTime);

            if (overwrite)
            {
                meshi.setInstance(masterInstDir);
            }
            else
            {
                meshi.setInstance(runTimei.timeName());
            }

            if (fvMeshDistribute::debug)
            {
                label meshNo = meshNoInProc[Pstream::myProcNo()][i];
                Pout<< "Writing domain "
                    << meshNo
                    << " to " << meshi.time().path()/runTimei.timeName()
                    << endl;
            }

            meshi.write();

            forAll(pointFieldNames, i)
            {
                IOobject io
                (
                    pointFieldNames[i],
                    runTimei.timeName(),
                    meshi
                );

                fileName fieldFile(io.objectPath());
                if (topoSet::debug) DebugVar(fieldFile);
                rm(fieldFile);
            }
        }

        // wait for all the processors to finish writing
        Pstream::barrier();
    }
    else
    {
        mesh.write();
        forAll(pointFieldNames, i)
        {
            IOobject io
            (
                pointFieldNames[i],
                runTime.timeName(),
                mesh
            );

            fileName fieldFile(io.objectPath());
            if (topoSet::debug) DebugVar(fieldFile);
            rm(fieldFile);
        }
    }
    Info<< "Written redistributed mesh to " << mesh.facesInstance() << nl
        << endl;


    if (decompose || nDestProcs == 1)
    {
        writeProcAddressing(decompose, meshSubDir, mesh, map);
    }


    // Refinement data
    if (nDestProcs > Pstream::nProcs())
    {
        labelList nMeshesInProc(Pstream::nProcs());
        labelList meshToProc(identity(nDestProcs));
        labelListList meshNoInProc =
            fvMeshDistribute::calcMeshToProcMap
            (
                nDestProcs,
                nMeshesInProc,
                meshToProc
            );

        bool oldParRun = UPstream::parRun();
        UPstream::parRun() = false;

        hexRef8DataList refData
        (
            myProcMeshes,
            meshToProc,
            meshNoInProc
        );

        UPstream::parRun() = oldParRun;

        //TODO refinementHistory not currently distributed
        refData.distribute(myProcDist);

        refData.write();
    }
    else
    {
        if (Pstream::master() && decompose)
        {
            runTime.TimePaths::caseName() = baseRunTime.caseName();
        }

        word hexRef8Instance =
            mesh.time().findInstance
            (
                "polyMesh",
                "cellLevel",
                IOobject::READ_IF_PRESENT
            );
        Pstream::scatter(hexRef8Instance);

        IOobject io
        (
            "dummy",
            hexRef8Instance,
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        hexRef8Data refData(io);
        if (Pstream::master() && decompose)
        {
            runTime.TimePaths::caseName() = proc0CaseName;
        }

        // Make sure all processors have valid data (since only some will
        // read)
        refData.sync(io);

        // Distribute
        refData.distribute(map);

        // Now we've read refinement data we can remove it
        word writeInstance = mesh.facesInstance();

        mesh.setInstance(hexRef8Instance);

        meshRefinement::removeFiles(mesh);

        mesh.setInstance(writeInstance);

        if (nDestProcs == 1)
        {
            if (Pstream::master())
            {
                Info<< "Setting caseName to " << baseRunTime.caseName()
                    << " to write reconstructed refinement data." << endl;
                runTime.TimePaths::caseName() = baseRunTime.caseName();

                refData.write();

                // Now we've written all. Reset caseName on master
                Info<< "Restoring caseName to " << proc0CaseName << endl;
                runTime.TimePaths::caseName() = proc0CaseName;
            }
        }
        else
        {
            refData.write();
        }
    } //hexRef8Data

    if (updateCellSets)
    {
        Info<< nl << "Distributing cellSets:";

        if (nDestProcs == 1)
        {
            topoSet::setInstance(mesh.facesInstance(), cellSets);

            if (Pstream::master())
            {
                Info<< "Setting caseName to " << baseRunTime.caseName()
                    << " to write reconstructed refinement data." << endl;
                runTime.TimePaths::caseName() = baseRunTime.caseName();

                forAll(cellSets, i)
                {
                    Info<< " " << cellSets[i].name();
                    cellSets[i].distribute(map);
                    cellSets[i].write();
                }
                Info<< endl;

                // Now we've written all. Reset caseName on master
                Info<< "Restoring caseName to " << proc0CaseName << endl;
                runTime.TimePaths::caseName() = proc0CaseName;
            }
        }
        else if (nDestProcs > Pstream::nProcs())
        {
            labelList nMeshesInProc(Pstream::nProcs());
            labelList meshToProc(identity(nDestProcs));
            labelListList meshNoInProc =
                fvMeshDistribute::calcMeshToProcMap
                (
                    nDestProcs,
                    nMeshesInProc,
                    meshToProc
                );

            const label nMyMeshes = meshNoInProc[Pstream::myProcNo()].size();

            forAll(myProcMeshes, i)
            {
                fvMesh& meshi = myProcMeshes[i];

                topoSet::setInstance
                (
                    meshi.facesInstance(),
                    myProcCellSets[i]
                );
            }

            forAll(myProcCellSets[0], csI)
            {
                Info<< " " << myProcCellSets[0][csI].name();

                PtrList<boolIOList> inSetPtrs(nMyMeshes);
                labelList           constructSizes(nMyMeshes);
                List<labelListList> subMaps(nMyMeshes);
                List<labelListList> constructMaps(nMyMeshes);
                boolList            subHasFlips(nMyMeshes);
                boolList            constructHasFlips(nMyMeshes);

                forAll(myProcDist, i)
                {
                    fvMesh& meshi = myProcMeshes[i];

                    IOobject io
                    (
                        "dummy",
                        meshi.facesInstance(),
                        polyMesh::meshSubDir,
                        meshi,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    );

                    if (i == 0)
                    {
                        inSetPtrs.set
                        (
                            i,
                            new boolIOList
                            (
                                io,
                                boolList(myProcDist[i].nOldCells(), false)
                            )
                        );

                        forAllConstIter(cellSet, myProcCellSets[i][csI], iter)
                        {
                            inSetPtrs[i][iter.key()] = true;
                        }
                    }
                    else
                    {
                        inSetPtrs.set
                        (
                            i,
                            new boolIOList
                            (
                                io,
                                boolList(myProcMeshes[i].nCells(), false)
                            )
                        );

                    }

                    constructSizes[i]    = myProcDist[i].cellMap().constructSize();
                    subMaps[i]           = myProcDist[i].cellMap().subMap();
                    constructMaps[i]     = myProcDist[i].cellMap().constructMap();
                    subHasFlips[i]       = myProcDist[i].cellMap().subHasFlip();
                    constructHasFlips[i] = myProcDist[i].cellMap().constructHasFlip();
                }

                mapDistributeBase::distribute
                (
                    constructSizes,
                    subMaps,
                    subHasFlips,
                    constructMaps,
                    constructHasFlips,
                    inSetPtrs,
                    flipOp(),
                    meshToProc,
                    meshNoInProc
                );

                forAll(myProcDist, i)
                {
                    // Count
                    label n = 0;
                    forAll(inSetPtrs[i], celli)
                    {
                        if (inSetPtrs[i][celli])
                        {
                            n++;
                        }
                    }

                    myProcCellSets[i][csI].clear();
                    myProcCellSets[i][csI].resize(n);
                    forAll(inSetPtrs[i], celli)
                    {
                        if (inSetPtrs[i][celli])
                        {
                            myProcCellSets[i][csI].insert(celli);
                        }
                    }

                    myProcCellSets[i][csI].write();
                }
            } // forAll cellSets
        }
        else
        {
            topoSet::setInstance(mesh.facesInstance(), cellSets);

            forAll(cellSets, i)
            {
                Info<< " " << cellSets[i].name();
                cellSets[i].distribute(map);
                cellSets[i].write();
            }
            Info<< endl;
        }
    }


    return autoPtr<mapDistributePolyMesh>
    (
        new mapDistributePolyMesh(map.xfer())
    );
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Field Mapping
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

autoPtr<mapDistributePolyMesh> createReconstructMap
(
    const autoPtr<fvMesh>& baseMeshPtr,
    const fvMesh& mesh,
    const labelList& cellProcAddressing,
    const labelList& faceProcAddressing,
    const labelList& pointProcAddressing,
    const labelList& boundaryProcAddressing
)
{
    // Send addressing to master
    labelListList cellAddressing(Pstream::nProcs());
    cellAddressing[Pstream::myProcNo()] = cellProcAddressing;
    Pstream::gatherList(cellAddressing);

    labelListList faceAddressing(Pstream::nProcs());
    faceAddressing[Pstream::myProcNo()] = faceProcAddressing;
    Pstream::gatherList(faceAddressing);

    labelListList pointAddressing(Pstream::nProcs());
    pointAddressing[Pstream::myProcNo()] = pointProcAddressing;
    Pstream::gatherList(pointAddressing);

    labelListList boundaryAddressing(Pstream::nProcs());
    {
        // Remove -1 entries
        DynamicList<label> patchProcAddressing(boundaryProcAddressing.size());
        forAll(boundaryProcAddressing, i)
        {
            if (boundaryProcAddressing[i] != -1)
            {
                patchProcAddressing.append(boundaryProcAddressing[i]);
            }
        }
        boundaryAddressing[Pstream::myProcNo()] = patchProcAddressing;
        Pstream::gatherList(boundaryAddressing);
    }


    autoPtr<mapDistributePolyMesh> mapPtr;

    if (baseMeshPtr.valid() && baseMeshPtr().nCells())    //baseMeshPtr.valid())
    {
        const fvMesh& baseMesh = baseMeshPtr();

        labelListList cellSubMap(Pstream::nProcs());
        cellSubMap[Pstream::masterNo()] = identity(mesh.nCells());

        mapDistribute cellMap
        (
            baseMesh.nCells(),
            cellSubMap.xfer(),
            cellAddressing.xfer()
        );

        labelListList faceSubMap(Pstream::nProcs());
        faceSubMap[Pstream::masterNo()] = identity(mesh.nFaces());

        mapDistribute faceMap
        (
            baseMesh.nFaces(),
            faceSubMap.xfer(),
            faceAddressing.xfer(),
            false,          //subHasFlip
            true            //constructHasFlip
        );

        labelListList pointSubMap(Pstream::nProcs());
        pointSubMap[Pstream::masterNo()] = identity(mesh.nPoints());

        mapDistribute pointMap
        (
            baseMesh.nPoints(),
            pointSubMap.xfer(),
            pointAddressing.xfer()
        );

        labelListList patchSubMap(Pstream::nProcs());
        // Send (filtered) patches to master
        patchSubMap[Pstream::masterNo()] =
            boundaryAddressing[Pstream::myProcNo()];

        mapDistribute patchMap
        (
            baseMesh.boundaryMesh().size(),
            patchSubMap.xfer(),
            boundaryAddressing.xfer()
        );

        const label nOldPoints = mesh.nPoints();
        const label nOldFaces = mesh.nFaces();
        const label nOldCells = mesh.nCells();

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        labelList oldPatchStarts(pbm.size());
        labelList oldPatchNMeshPoints(pbm.size());
        forAll(pbm, patchI)
        {
            oldPatchStarts[patchI] = pbm[patchI].start();
            oldPatchNMeshPoints[patchI] = pbm[patchI].nPoints();
        }

        mapPtr.reset
        (
            new mapDistributePolyMesh
            (
                nOldPoints,
                nOldFaces,
                nOldCells,
                oldPatchStarts.xfer(),
                oldPatchNMeshPoints.xfer(),
                pointMap.xfer(),
                faceMap.xfer(),
                cellMap.xfer(),
                patchMap.xfer()
            )
        );
    }
    else
    {
        labelListList cellSubMap(Pstream::nProcs());
        cellSubMap[Pstream::masterNo()] = identity(mesh.nCells());
        labelListList cellConstructMap(Pstream::nProcs());

        mapDistribute cellMap
        (
            0,
            cellSubMap.xfer(),
            cellConstructMap.xfer()
        );

        labelListList faceSubMap(Pstream::nProcs());
        faceSubMap[Pstream::masterNo()] = identity(mesh.nFaces());
        labelListList faceConstructMap(Pstream::nProcs());

        mapDistribute faceMap
        (
            0,
            faceSubMap.xfer(),
            faceConstructMap.xfer(),
            false,          //subHasFlip
            true            //constructHasFlip
        );

        labelListList pointSubMap(Pstream::nProcs());
        pointSubMap[Pstream::masterNo()] = identity(mesh.nPoints());
        labelListList pointConstructMap(Pstream::nProcs());

        mapDistribute pointMap
        (
            0,
            pointSubMap.xfer(),
            pointConstructMap.xfer()
        );

        labelListList patchSubMap(Pstream::nProcs());
        // Send (filtered) patches to master
        patchSubMap[Pstream::masterNo()] =
            boundaryAddressing[Pstream::myProcNo()];
        labelListList patchConstructMap(Pstream::nProcs());

        mapDistribute patchMap
        (
            0,
            patchSubMap.xfer(),
            patchConstructMap.xfer()
        );

        const label nOldPoints = mesh.nPoints();
        const label nOldFaces = mesh.nFaces();
        const label nOldCells = mesh.nCells();

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        labelList oldPatchStarts(pbm.size());
        labelList oldPatchNMeshPoints(pbm.size());
        forAll(pbm, patchI)
        {
            oldPatchStarts[patchI] = pbm[patchI].start();
            oldPatchNMeshPoints[patchI] = pbm[patchI].nPoints();
        }

        mapPtr.reset
        (
            new mapDistributePolyMesh
            (
                nOldPoints,
                nOldFaces,
                nOldCells,
                oldPatchStarts.xfer(),
                oldPatchNMeshPoints.xfer(),
                pointMap.xfer(),
                faceMap.xfer(),
                cellMap.xfer(),
                patchMap.xfer()
            )
        );
    }

    return mapPtr;
}


void readProcAddressing
(
    const fvMesh& mesh,
    const autoPtr<fvMesh>& baseMeshPtr,
    autoPtr<mapDistributePolyMesh>& distMap
)
{
    //IOobject io
    //(
    //    "procAddressing",
    //    mesh.facesInstance(),
    //    polyMesh::meshSubDir,
    //    mesh,
    //    IOobject::MUST_READ
    //);
    //if (io.typeHeaderOk<labelIOList>(true))
    //{
    //    Pout<< "Reading addressing from " << io.name() << " at "
    //        << mesh.facesInstance() << nl << endl;
    //    distMap.clear();
    //    distMap.reset(new IOmapDistributePolyMesh(io));
    //}
    //else
    {
        Info<< "Reading addressing from procXXXAddressing at "
            << mesh.facesInstance() << nl << endl;
        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            labelList(0)
        );
        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            labelList(0)
        );
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            labelList(0)
        );
        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            labelList(0)
        );


        if
        (
            mesh.nCells() != cellProcAddressing.size()
         || mesh.nPoints() != pointProcAddressing.size()
         || mesh.nFaces() != faceProcAddressing.size()
        )
        {
            FatalErrorInFunction
                << "Read addressing inconsistent with mesh sizes" << nl
                << "cells:" << mesh.nCells()
                << " addressing:" << cellProcAddressing.objectPath()
                << " size:" << cellProcAddressing.size() << nl
                << "faces:" << mesh.nFaces()
                << " addressing:" << faceProcAddressing.objectPath()
                << " size:" << faceProcAddressing.size() << nl
                << "points:" << mesh.nPoints()
                << " addressing:" << pointProcAddressing.objectPath()
                << " size:" << pointProcAddressing.size()
                << exit(FatalError);
        }

        distMap.clear();
        distMap = createReconstructMap
        (
            baseMeshPtr,
            mesh,
            cellProcAddressing,
            faceProcAddressing,
            pointProcAddressing,
            boundaryProcAddressing
        );
    }
}


void reconstructMeshFields
(
    const parFvFieldReconstructor& fvReconstructor,
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
)
{
    // Dimensioned fields

    fvReconstructor.reconstructFvVolumeInternalFields<scalar>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeInternalFields<vector>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeInternalFields<sphericalTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeInternalFields<symmTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeInternalFields<tensor>
    (
        objects,
        selectedFields
    );


    // volFields

    fvReconstructor.reconstructFvVolumeFields<scalar>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeFields<vector>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeFields<sphericalTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeFields<symmTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeFields<tensor>
    (
        objects,
        selectedFields
    );


    // surfaceFields

    fvReconstructor.reconstructFvSurfaceFields<scalar>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvSurfaceFields<vector>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvSurfaceFields<sphericalTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvSurfaceFields<symmTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvSurfaceFields<tensor>
    (
        objects,
        selectedFields
    );
}


void reconstructLagrangian
(
    autoPtr<parLagrangianRedistributor>& lagrangianReconstructorPtr,
    const fvMesh& baseMesh,
    const fvMesh& mesh,
    const mapDistributePolyMesh& distMap,
    const HashSet<word>& selectedLagrangianFields
)
{
    // Clouds (note: might not be present on all processors)

    wordList cloudNames;
    List<wordList> fieldNames;
    // Find all cloudNames on all processors
    parLagrangianRedistributor::findClouds(mesh, cloudNames, fieldNames);

    if (cloudNames.size())
    {
        if (!lagrangianReconstructorPtr.valid())
        {
            lagrangianReconstructorPtr.reset
            (
                new parLagrangianRedistributor
                (
                    mesh,
                    baseMesh,
                    mesh.nCells(),      // range of cell indices in clouds
                    distMap
                )
            );
        }
        const parLagrangianRedistributor& lagrangianReconstructor =
            lagrangianReconstructorPtr();

        forAll(cloudNames, i)
        {
            Info<< "Reconstructing lagrangian fields for cloud "
                << cloudNames[i] << nl << endl;

            autoPtr<mapDistributeBase> lagrangianMap =
            lagrangianReconstructor.redistributeLagrangianPositions
            (
                cloudNames[i]
            );
            IOobjectList sprayObjs
            (
                mesh,
                mesh.time().timeName(),
                cloud::prefix/cloudNames[i]
            );

            lagrangianReconstructor.redistributeLagrangianFields<label>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFieldFields<label>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFields<scalar>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFieldFields<scalar>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFields<vector>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFieldFields<vector>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFields
            <sphericalTensor>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFieldFields
            <sphericalTensor>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFields<symmTensor>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFieldFields
            <symmTensor>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFields<tensor>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeLagrangianFieldFields<tensor>
            (
                lagrangianMap,
                cloudNames[i],
                sprayObjs,
                selectedLagrangianFields
            );
        }
    }
}

// Read clouds (note: might not be present on all processors)
void readLagrangian
(
    const fvMesh& mesh,
    const wordList& cloudNames,
    const HashSet<word>& selectedLagrangianFields,
    PtrList<unmappedPassiveParticleCloud>& clouds
)
{
    (void)mesh.tetBasePtIs();

    forAll(cloudNames, i)
    {
        //Pout<< "Loading cloud " << cloudNames[i] << endl;
        clouds.set
        (
            i,
            new unmappedPassiveParticleCloud(mesh, cloudNames[i], false)
        );


        //forAllConstIter
        //(
        //    unmappedPassiveParticleCloud,
        //    clouds[i],
        //    iter
        //)
        //{
        //    Pout<< "Particle position:" << iter().position()
        //        << " cell:" << iter().cell()
        //        << " with cc:" << mesh.cellCentres()[iter().cell()]
        //        << endl;
        //}


        IOobjectList sprayObjs(clouds[i], clouds[i].time().timeName());

        //Pout<< "Found clould objects:" << sprayObjs.names() << endl;

        parLagrangianRedistributor::readLagrangianFields
        <IOField<label>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <IOField<Field<label>>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <CompactIOField<Field<label>, label>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readLagrangianFields
        <IOField<scalar>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <IOField<Field<scalar>>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <CompactIOField<Field<scalar>, scalar>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readLagrangianFields
        <IOField<vector>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <IOField<Field<vector>>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <CompactIOField<Field<vector>, vector>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readLagrangianFields
        <IOField<sphericalTensor>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <IOField<Field<sphericalTensor>>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <CompactIOField<Field<sphericalTensor>, sphericalTensor>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readLagrangianFields
        <IOField<symmTensor>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <IOField<Field<symmTensor>>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <CompactIOField<Field<symmTensor>, symmTensor>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readLagrangianFields
        <IOField<tensor>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <IOField<Field<tensor>>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readLagrangianFields
        <CompactIOField<Field<tensor>, tensor>>
        (
            clouds[i],
            sprayObjs,
            selectedLagrangianFields
        );
    }
}


void redistributeLagrangian
(
    autoPtr<parLagrangianRedistributor>& lagrangianReconstructorPtr,
    const fvMesh& mesh,
    const label nOldCells,
    const mapDistributePolyMesh& distMap,
    PtrList<unmappedPassiveParticleCloud>& clouds
)
{
    if (clouds.size())
    {
        if (!lagrangianReconstructorPtr.valid())
        {
            lagrangianReconstructorPtr.reset
            (
                new parLagrangianRedistributor
                (
                    mesh,
                    mesh,
                    nOldCells,  // range of cell indices in clouds
                    distMap
                )
            );
        }
        const parLagrangianRedistributor& distributor =
            lagrangianReconstructorPtr();

        forAll(clouds, i)
        {
            autoPtr<mapDistributeBase> lagrangianMap =
            distributor.redistributeLagrangianPositions(clouds[i]);

            distributor.redistributeStoredLagrangianFields
            <IOField<label>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <IOField<Field<label>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <CompactIOField<Field<label>, label>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredLagrangianFields
            <IOField<scalar>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <IOField<Field<scalar>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <CompactIOField<Field<scalar>, scalar>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredLagrangianFields
            <IOField<vector>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <IOField<Field<vector>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <CompactIOField<Field<vector>, vector>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredLagrangianFields
            <IOField<sphericalTensor>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <IOField<Field<sphericalTensor>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <CompactIOField<Field<sphericalTensor>, sphericalTensor>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredLagrangianFields
            <IOField<symmTensor>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <IOField<Field<symmTensor>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <CompactIOField<Field<symmTensor>, symmTensor>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredLagrangianFields
            <IOField<tensor>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <IOField<Field<tensor>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredLagrangianFields
            <CompactIOField<Field<tensor>, tensor>>
            (
                lagrangianMap,
                clouds[i]
            );
        }
    }
}


autoPtr<fvMesh> detectAndGetMesh
(
    const word& regionName,
    const fileName& meshSubDir,
    Time& runTime,
    const Time& baseRunTime,
    const scalar& masterTime,
    const word proc0CaseName,
    const bool decompose,
    const bool removeDummyMesh,

    boolList& haveMesh,
    fileName& masterInstDir
)
{
    Info<< "Using mesh subdirectory " << meshSubDir << nl << endl;

    Info<< "Setting time to that of master or undecomposed case : "
        << masterTime << endl;
    runTime.setTime(masterTime, 0);


    // Get time instance directory
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // At this point we should be able to read at least a mesh on
    // processor0. Note the changing of the processor0 casename to
    // enforce it to read/write from the undecomposed case

    if (Pstream::master())
    {
        if (decompose)
        {
            Info<< "Setting caseName to " << baseRunTime.caseName()
                << " to find undecomposed mesh" << endl;
            runTime.TimePaths::caseName() = baseRunTime.caseName();
        }

        masterInstDir = runTime.findInstance
        (
            meshSubDir,
            "faces",
            IOobject::READ_IF_PRESENT
        );

        if (decompose)
        {
            Info<< "Restoring caseName to " << proc0CaseName << endl;
            runTime.TimePaths::caseName() = proc0CaseName;
        }
    }
    Pstream::scatter(masterInstDir);

    // Check who has a mesh
    const fileName meshPath =
        runTime.path()/masterInstDir/meshSubDir/"faces";

    Info<< "Checking for mesh in " << meshPath << nl << endl;


    haveMesh = boolList(Pstream::nProcs(), false);
    haveMesh[Pstream::myProcNo()] = isFile(meshPath);
    Pstream::allGatherList(haveMesh);
    Info<< "Per processor mesh availability : " << haveMesh << endl;

    // Load mesh (or create dummy one)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (Pstream::master() && decompose)
    {
        Info<< "Setting caseName to " << baseRunTime.caseName()
            << " to read undecomposed mesh" << endl;
        runTime.TimePaths::caseName() = baseRunTime.caseName();
    }

    autoPtr<fvMesh> meshPtr = loadOrCreateMesh
    (
        IOobject
        (
            regionName,
            masterInstDir,
            runTime,
            Foam::IOobject::MUST_READ
        ),
        removeDummyMesh
    );

    if (Pstream::master() && decompose)
    {
        Info<< "Restoring caseName to " << proc0CaseName << endl;
        runTime.TimePaths::caseName() = proc0CaseName;
    }

    return meshPtr;
}


int main(int argc, char *argv[])
{
    // enable -constant ... if someone really wants it
    // enable -zeroTime to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);
    #include "include/addRegionOption.H"
    argList::addBoolOption
    (
        "allRegions",
        "operate on all regions in regionProperties"
    );
    #include "include/addOverwriteOption.H"
    argList::addBoolOption("decompose", "decompose case");
    argList::addBoolOption("reconstruct", "reconstruct case");
    argList::addOption
    (
        "mergeTol",
        "scalar",
        "specify the merge distance relative to the bounding box size "
        "(default 1e-6)"
    );
    argList::addBoolOption
    (
        "cellDist",
        "write cell distribution as a labelList - for use with 'manual' "
        "decomposition method or as a volScalarField for post-processing."
    );
    argList::addBoolOption
    (
        "newTimes",
        "only reconstruct new times (i.e. that do not exist already)"
    );

#if defined( WIN32 ) || defined( WIN64 )
//- needed to force load ptscotch and othe libs at runtime -DC
#include "include/forceLoadLibraries.H"
#endif
    argList::addBoolOption
    (
        "cellSets",
        "Maintain cellSets if present"
    );

    argList::noCheckProcessorDirectories();

    // Handle arguments
    // ~~~~~~~~~~~~~~~~
    // (replacement for setRootCase that does not abort)

    Foam::argList args(argc, argv);
    bool decompose = args.optionFound("decompose");
    const bool reconstruct = args.optionFound("reconstruct");
    bool overwrite = args.optionFound("overwrite");
    bool writeCellDist = args.optionFound("cellDist");
    bool newTimes = args.optionFound("newTimes");
    bool updateCellSets = args.optionFound("cellSets");
    bool allRegions = args.optionFound("allRegions");

    if (env("FOAM_SIGFPE") || env("FOAM_SIGFPE"))
    {
        WarningInFunction
            << "Detected floating point exception trapping (FOAM_SIGFPE)."
            << " This might give" << nl
            << "    problems when mapping fields. Switch it off in case"
            << " of problems." << endl;
    }



    const HashSet<word> selectedFields(0);
    const HashSet<word> selectedLagrangianFields(0);


    if (decompose)
    {
        Info<< "Decomposing case (like decomposePar)" << nl << endl;
        if (reconstruct)
        {
            FatalErrorInFunction
                << "Cannot specify both -decompose and -reconstruct"
                << exit(FatalError);
        }
    }
    else if (reconstruct)
    {
        Info<< "\nReconstructing case (like reconstructParMesh)" << nl << endl;

        DeprecationWarningInFunction("reconstruct", "functionality", 40200)
            << "Please use the 'reconstructParMesh' utility (for the geometry) "
            << "and 'reconstructPar' (for the fields) instead," << nl << "    "
            << "or redistribute the mesh to a small number of processors for "
            << "archiving purposes." << nl << endl;
    }


    if (decompose || reconstruct)
    {
        if (!overwrite)
        {
            WarningInFunction
                << "Working in decompose or reconstruction mode automatically"
                << " implies -overwrite" << nl << endl;
            overwrite = true;
        }
    }


    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << ": This utility can only be run parallel"
            << exit(FatalError);
    }


    if (!isDir(args.rootPath()))
    {
        FatalErrorInFunction
            << ": cannot open root directory " << args.rootPath()
            << exit(FatalError);
    }

    // Detect if running data-distributed (multiple roots)
    bool nfs = true;
    {
        List<fileName> roots(1, args.rootPath());
        combineReduce(roots, ListUniqueEqOp<fileName>());
        nfs = (roots.size() == 1);
    }

    if (!nfs)
    {
        Info<< "Detected multiple roots i.e. non-nfs running"
            << nl << endl;
    }

    if (isDir(args.path()))
    {
        if (decompose)
        {
            Info<< "Removing existing processor directories" << endl;
            rmDir(args.path());
        }
    }
    else
    {
        // Directory does not exist. If this happens on master -> decompose mode
        decompose = true;
        Info<< "No processor directories; switching on decompose mode"
            << nl << endl;
    }
    // If master changed to decompose mode make sure all nodes know about
    // it
    Pstream::scatter(decompose);


    // If running distributed we have problem of new processors not finding
    // a system/controlDict. However if we switch on the master-only reading
    // the problem becomes that the time directories are differing sizes and
    // e.g. latestTime will pick up a different time (which causes createTime.H
    // to abort). So for now make sure to have master times on all
    // processors
    {
        Info<< "Creating time directories on all processors" << nl << endl;
        instantList timeDirs;
        if (Pstream::master())
        {
            const bool oldParRun = Pstream::parRun();
            Pstream::parRun() = false;
            timeDirs = Time::findTimes(args.path(), "constant");
            Pstream::parRun() = oldParRun;
        }
        Pstream::scatter(timeDirs);
        forAll(timeDirs, i)
        {
            mkDir(args.path()/timeDirs[i].name());
        }
    }


    // Construct time
    // ~~~~~~~~~~~~~~

    #include "include/createTime.H"
    runTime.functionObjects().off();


    // Save local processor0 casename
    const fileName proc0CaseName = runTime.caseName();


    // Construct undecomposed Time
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This will read the same controlDict but might have a different
    // set of times so enforce same times

    if (!nfs)
    {
        Info<< "Creating time directories for undecomposed Time"
            << " on all processors" << nl << endl;
        instantList timeDirs;

        const fileName basePath(args.rootPath()/args.globalCaseName());

        if (Pstream::master())
        {
            const bool oldParRun = Pstream::parRun();
            Pstream::parRun() = false;
            timeDirs = Time::findTimes(basePath, "constant");
            Pstream::parRun() = oldParRun;
        }
        Pstream::scatter(timeDirs);
        forAll(timeDirs, i)
        {
            mkDir(basePath/timeDirs[i].name());
        }
    }


    Info<< "Create undecomposed database"<< nl << endl;
    Time baseRunTime
    (
        runTime.controlDict(),
        runTime.rootPath(),
        runTime.globalCaseName(),
        runTime.system(),
        runTime.constant(),
        false                   // enableFunctionObjects
    );


    HashSet<word> masterTimeDirSet;
    if (newTimes)
    {
        instantList baseTimeDirs(baseRunTime.times());
        forAll(baseTimeDirs, i)
        {
            masterTimeDirSet.insert(baseTimeDirs[i].name());
        }
    }


    // Determine any region
    word regionName = polyMesh::defaultRegion;
    fileName meshSubDir;
    hashedWordList regionNames;
    List<fileName> meshSubDirs;

    if (allRegions)
    {
        Info<< "Reconstructing for all regions in regionProperties" << nl
            << endl;
        regionProperties rp(runTime);
        forAllConstIter(HashTable<wordList>, rp, iter)
        {
            const wordList& regions = iter();
            forAll(regions, i)
            {
                if (findIndex(regionNames, regions[i]) == -1)
                {
                    regionNames.append(regions[i]);
                    meshSubDirs.append(regions[i]/polyMesh::meshSubDir);
                }
            }
        }
        if (!rp.found())
        {
            FatalErrorInFunction
                << "Region names must be specified in the regionProperties "
                << "dictionary when the -allRegions switch is used"
                << nl << exit(FatalError);
        }
    }
    else
    {
        word regionName;
        if (args.optionReadIfPresent("region", regionName))
        {
            regionNames = wordList(1, regionName);
            meshSubDirs = List<fileName>(1,regionName/polyMesh::meshSubDir);
        }
        else
        {
            regionNames = wordList(1, polyMesh::defaultRegion);
            meshSubDirs = List<fileName>(1,polyMesh::meshSubDir);
        }
    }

    dimensionedScalar origTime = runTime.startTime();
    label origTimeIndex = runTime.startTimeIndex();

    scalar masterTime;
    if (reconstruct)
    {
        forAll(regionNames, regionI)
        {
            const word& regionName = regionNames[regionI];
            const fileName meshSubDir = meshSubDirs[regionI];

            Info<< "Using mesh subdirectory " << meshSubDir << nl << endl;

            // Demand driven lagrangian mapper
            autoPtr<parLagrangianRedistributor> lagrangianReconstructorPtr;

            // use the times list from the master processor
            // and select a subset based on the command-line options
            instantList timeDirs = timeSelector::select(runTime.times(), args);
            Pstream::scatter(timeDirs);

            if (timeDirs.empty())
            {
                FatalErrorInFunction
                    << "No times selected"
                    << exit(FatalError);
            }


            // Pass1 : reconstruct mesh and addressing
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


            Info<< nl
                << "Pass1 : reconstructing mesh and addressing" << nl << endl;

            // Loop over all times
            forAll(timeDirs, timeI)
            {
                // Set time for global database
                runTime.setTime(timeDirs[timeI], timeI);
                baseRunTime.setTime(timeDirs[timeI], timeI);

                Info<< "Time = " << runTime.timeName() << endl << endl;


                // See where the mesh is
                fileName facesInstance = runTime.findInstance
                (
                    meshSubDir,
                    "faces",
                    IOobject::READ_IF_PRESENT
                );
                //Pout<< "facesInstance:" << facesInstance << endl;

                Pstream::scatter(facesInstance);

                // Check who has a mesh
                const fileName meshPath =
                    runTime.path()/facesInstance/meshSubDir/"faces";

                Info<< "Checking for mesh in " << meshPath << nl << endl;

                boolList haveMesh(Pstream::nProcs(), false);
                haveMesh[Pstream::myProcNo()] = isFile(meshPath);
                Pstream::allGatherList(haveMesh);
                Info<< "Per processor mesh availability : " << haveMesh << endl;


                // Addressing back to reconstructed mesh as xxxProcAddressing.
                // - all processors have consistent faceProcAddressing
                // - processors with no mesh don't need faceProcAddressing


                // Note: filePath searches up on processors that don't have
                //       processor if instance = constant so explicitly check
                //       found filename.
                bool haveAddressing = false;
                if (haveMesh[Pstream::myProcNo()])
                {
                    // Read faces (just to know their size)
                    faceCompactIOList faces
                    (
                        IOobject
                        (
                            "faces",
                            facesInstance,
                            meshSubDir,
                            runTime,
                            IOobject::MUST_READ
                        )
                    );

                    // Check faceProcAddressing
                    labelIOList faceProcAddressing
                    (
                        IOobject
                        (
                            "faceProcAddressing",
                            facesInstance,
                            meshSubDir,
                            runTime,
                            IOobject::READ_IF_PRESENT
                        ),
                        labelList(0)
                    );
                    if
                    (
                        faceProcAddressing.headerOk()
                     && faceProcAddressing.size() == faces.size()
                    )
                    {
                        haveAddressing = true;
                    }
                }
                else
                {
                    // Have no mesh. Don't need addressing
                    haveAddressing = true;
                }

                if (!returnReduce(haveAddressing, andOp<bool>()))
                {
                    Info<< "loading mesh from " << facesInstance << endl;
                    autoPtr<fvMesh> meshPtr = loadOrCreateMesh
                    (
                        IOobject
                        (
                            regionName,
                            facesInstance,
                            runTime,
                            Foam::IOobject::MUST_READ
                        )
                    );
                    fvMesh& mesh = meshPtr();

                    // Global matching tolerance
                    const scalar tolDim = getMergeDistance
                    (
                        args,
                        runTime,
                        mesh.bounds()
                    );


                    // Determine decomposition
                    // ~~~~~~~~~~~~~~~~~~~~~~~

                    Info<< "Reconstructing mesh for time "
                        << facesInstance << endl;

                    label nDestProcs = 1;
                    labelList finalDecomp = labelList(mesh.nCells(), 0);

                    redistributeAndWrite
                    (
                        baseRunTime,
                        tolDim,
                        haveMesh,
                        meshSubDir,
                        false,   // do not read fields
                        false,   // do not read undecomposed case on processor0
                        overwrite,
                        updateCellSets,
                        proc0CaseName,
                        nDestProcs,
                        finalDecomp,
                        facesInstance,
                        mesh
                    );
                }
            }


            // Pass2 : read mesh and addressing and reconstruct fields
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Info<< nl
                << "Pass2 : reconstructing fields" << nl << endl;

            runTime.setTime(timeDirs[0], 0);
            baseRunTime.setTime(timeDirs[0], 0);
            Info<< "Time = " << runTime.timeName() << endl << endl;


            Info<< "Reading undecomposed mesh (on master)" << endl;
            autoPtr<fvMesh> baseMeshPtr = fvMeshTools::newMesh
            (
                IOobject
                (
                    regionName,
                    baseRunTime.timeName(),
                    baseRunTime,
                    IOobject::MUST_READ
                ),
                true            // read on master only
            );

            Info<< "Reading local, decomposed mesh" << endl;
            autoPtr<fvMesh> meshPtr = loadOrCreateMesh
            (
                IOobject
                (
                    regionName,
                    baseMeshPtr().facesInstance(),
                    runTime,
                    Foam::IOobject::MUST_READ
                 )
            );
            fvMesh& mesh = meshPtr();


            // Read addressing back to base mesh
            autoPtr<mapDistributePolyMesh> distMap;
            readProcAddressing(mesh, baseMeshPtr, distMap);

            // Construct field mapper
            autoPtr<parFvFieldReconstructor> fvReconstructorPtr
            (
                new parFvFieldReconstructor
                (
                    baseMeshPtr(),
                    mesh,
                    distMap,
                    Pstream::master()       // do I need to write?
                )
            );



            // Since we start from Times[0] and not runTime.timeName() we
            // might overlook point motion in the first timestep
            // (since mesh.readUpdate() below will not be triggered). Instead
            // detect points by hand
            if (mesh.pointsInstance() != mesh.facesInstance())
            {
                Info<< "    Dected initial mesh motion; reconstructing points"
                    << nl << endl;
                fvReconstructorPtr().reconstructPoints();
            }


            // Loop over all times
            forAll(timeDirs, timeI)
            {
                if (newTimes && masterTimeDirSet.found(timeDirs[timeI].name()))
                {
                    Info<< "Skipping time " << timeDirs[timeI].name()
                        << endl << endl;
                    continue;
                }

                // Set time for global database
                runTime.setTime(timeDirs[timeI], timeI);
                baseRunTime.setTime(timeDirs[timeI], timeI);

                Info<< "Time = " << runTime.timeName() << endl << endl;


                // Check if any new meshes need to be read.
                fvMesh::readUpdateState procStat = mesh.readUpdate();

                if (procStat == fvMesh::POINTS_MOVED)
                {
                    Info<< "    Detected mesh motion; reconstructing points" << nl
                        << endl;
                    fvReconstructorPtr().reconstructPoints();
                }
                else if
                (
                    procStat == fvMesh::TOPO_CHANGE
                    || procStat == fvMesh::TOPO_PATCH_CHANGE
                )
                {
                    Info<< "    Detected topology change; reconstructing"
                        << " addressing"
                        << nl << endl;

                    if (baseMeshPtr.valid())
                    {
                        // Cannot do a baseMesh::readUpdate() since not all
                        // processors will have mesh files. So instead just
                        // recreate baseMesh
                        baseMeshPtr.clear();
                        baseMeshPtr = fvMeshTools::newMesh
                        (
                            IOobject
                            (
                                regionName,
                                baseRunTime.timeName(),
                                baseRunTime,
                                IOobject::MUST_READ
                            ),
                            true            // read on master only
                        );
                    }

                    // Re-read procXXXaddressing
                    readProcAddressing(mesh, baseMeshPtr, distMap);

                    // Reset field mapper
                    fvReconstructorPtr.reset
                    (
                        new parFvFieldReconstructor
                        (
                            baseMeshPtr(),
                            mesh,
                            distMap,
                            Pstream::master()
                        )
                    );
                    lagrangianReconstructorPtr.clear();
                }


                // Get list of objects
                IOobjectList objects(mesh, runTime.timeName());


                // Mesh fields (vol, surface, volInternal)
                reconstructMeshFields
                (
                    fvReconstructorPtr(),
                    objects,
                    selectedFields
                );

                // Clouds (note: might not be present on all processors)
                reconstructLagrangian
                (
                    lagrangianReconstructorPtr,
                    baseMeshPtr(),
                    mesh,
                    distMap,
                    selectedLagrangianFields
                );

                // If there are any "uniform" directories copy them from
                // the master processor
                if (Pstream::master())
                {
                    fileName uniformDir0 = runTime.timePath()/"uniform";
                    if (isDir(uniformDir0))
                    {
                        Info<< "Detected additional non-decomposed files in "
                            << uniformDir0 << endl;
                        cp(uniformDir0, baseRunTime.timePath());
                    }
                }
            }
            runTime.setTime(origTime,origTimeIndex);
        }
    }
    else
    {
        // Time coming from processor0 (or undecomposed if no processor0)

        if (decompose)
        {
            // Use base time. This is to handle e.g. startTime = latestTime
            // which will not do anything if there are no processor
            // directories
            masterTime = timeSelector::selectIfPresent
            (
                baseRunTime,
                args
            )[0].value();
        }
        else
        {
            masterTime = timeSelector::selectIfPresent
            (
                runTime,
                args
             )[0].value();
        }
        Pstream::scatter(masterTime);

        // Allow override of decomposeParDict location
        fileName decompDictFile;
        args.optionReadIfPresent("decomposeParDict", decompDictFile);

        HashPtrTable<labelList> predeterminedDecomps;
        HashTable<label> predeterminedNDestProcs;

        label nTargetProcs = -1;

        if (allRegions)
        {
            // Pull out the regions that are in a consolidated group
            // and calculate their decompositions first
            regionProperties rp(runTime);
            forAll(rp.consolidatedGroupNames(), groupi)
            {
                const word& groupName = rp.consolidatedGroupNames()[groupi];
                const wordList& groupRegionNames = rp[groupName];

                Info<< "Loading all meshes for consolidated group "
                    << groupName << nl << endl;

                PtrList<fvMesh> meshes(groupRegionNames.size());
                boundBox globalBb;
                forAll(groupRegionNames, groupRegioni)
                {
                    const word& regionName = groupRegionNames[groupRegioni];

                    boolList haveMesh;
                    fileName masterInstDir;

                    autoPtr<fvMesh> meshPtr =
                        detectAndGetMesh
                        (
                            regionName,
                            meshSubDirs[regionNames[regionName]],
                            runTime,
                            baseRunTime,
                            masterTime,
                            proc0CaseName,
                            decompose,
                            true,
                            haveMesh,
                            masterInstDir
                        );

                    meshes.set(groupRegioni, meshPtr);
                    globalBb += meshes[groupRegioni].bounds();
                }

                // Global matching tolerance
                const scalar tolDim = getMergeDistance
                (
                    args,
                    runTime,
                    globalBb
                );

                // Determine decomposition
                // ~~~~~~~~~~~~~~~~~~~~~~~

                Info<< "Calculating decomposition for consolidated group "
                    << groupName << nl << endl;

                label nDestProcs;
                PtrList<labelList> finalDecomps(meshes.size());
                forAll(finalDecomps, i)
                {
                    finalDecomps.set(i, new labelList);
                }
                determineDecomposition
                (
                    baseRunTime,
                    decompDictFile,
                    decompose,
                    proc0CaseName,
                    meshes,
                    writeCellDist,
                    tolDim,

                    nDestProcs,
                    finalDecomps
                );

                // Transfer pointers from finalDecomps into HashPtrTable
                forAll(finalDecomps, groupRegioni)
                {
                    predeterminedDecomps.set
                    (
                        groupRegionNames[groupRegioni],
                        finalDecomps.set(groupRegioni, nullptr).ptr()
                    );
                    predeterminedNDestProcs.set
                    (
                        groupRegionNames[groupRegioni],
                        nDestProcs
                    );
                }
            }
        }

        // Now do all regions, using the already calculated decomposition for
        // the consolidated groups
        forAll(regionNames, regionI)
        {
            const word& regionName = regionNames[regionI];
            const fileName meshSubDir = meshSubDirs[regionI];

            boolList haveMesh;
            fileName masterInstDir;

            autoPtr<fvMesh> meshPtr =
                detectAndGetMesh
                (
                    regionName,
                    meshSubDir,
                    runTime,
                    baseRunTime,
                    masterTime,
                    proc0CaseName,
                    decompose,
                    false,
                    haveMesh,
                    masterInstDir
                );

            fvMesh& mesh = meshPtr();

            label nOldCells = mesh.nCells();
            //Pout<< "Loaded mesh : nCells:" << nOldCells
            //    << " nPatches:" << mesh.boundaryMesh().size() << endl;


            // Global matching tolerance
            const scalar tolDim = getMergeDistance
            (
                args,
                runTime,
                mesh.bounds()
            );

            labelList finalDecomp;
            label nDestProcs;

            #if defined(WIN64) || defined(WIN32)
            Pstream::scatter(decompDictFile);
            #endif

            if (predeterminedDecomps.found(regionName))
            {
                Info<< "Using decomposition already determined for "
                    << "consolidated regions" << nl << endl;

                finalDecomp =
                    labelList(*predeterminedDecomps[regionName], true);
                nDestProcs = predeterminedNDestProcs[regionName];
            }
            else
            {
                // Determine decomposition
                // ~~~~~~~~~~~~~~~~~~~~~~~

                determineDecomposition
                (
                    baseRunTime,
                    decompDictFile,
                    decompose,
                    proc0CaseName,
                    mesh,
                    writeCellDist,
                    tolDim,

                    nDestProcs,
                    finalDecomp
                );
            }

            // Demand driven lagrangian mapper
            autoPtr<parLagrangianRedistributor> lagrangianReconstructorPtr;

            wordList cloudNames;
            List<wordList> fieldNames;

            // Detect lagrangian fields
            if (Pstream::master() && decompose)
            {
                runTime.TimePaths::caseName() = baseRunTime.caseName();
            }
            parLagrangianRedistributor::findClouds
            (
                mesh,
                cloudNames,
                fieldNames
            );

            // Read lagrangian fields and store on cloud (objectRegistry)
            PtrList<unmappedPassiveParticleCloud> clouds(cloudNames.size());
            readLagrangian
            (
                mesh,
                cloudNames,
                selectedLagrangianFields,
                clouds
            );
            if (Pstream::master() && decompose)
            {
                runTime.TimePaths::caseName() = proc0CaseName;
            }


            // Load fields, do all distribution (mesh and fields - but not
            // lagrangian fields; these are done later)
            autoPtr<mapDistributePolyMesh> distMap =
                redistributeAndWrite
                (
                    baseRunTime,
                    tolDim,
                    haveMesh,
                    meshSubDir,
                    true,           // read fields
                    decompose,      // decompose, i.e. read from undecomposed case
                    overwrite,
                    updateCellSets,
                    proc0CaseName,
                    nDestProcs,
                    finalDecomp,
                    masterInstDir,
                    mesh
                );


            // Redistribute any clouds
            redistributeLagrangian
            (
                lagrangianReconstructorPtr,
                mesh,
                nOldCells,
                distMap,
                clouds
            );


            // Copy any uniform data
            const fileName uniformDir("uniform");
            if (isDir(baseRunTime.timePath()/uniformDir))
            {
                Info<< "Detected additional non-decomposed files in "
                    << baseRunTime.timePath()/uniformDir << endl;
                cp
                (
                    baseRunTime.timePath()/uniformDir,
                    runTime.timePath()/uniformDir
                );
            }

            nTargetProcs = max(nTargetProcs, nDestProcs);
        }

        // Wait for all threads to get here.
        Pstream::barrier();

        if
        (
            Pstream::master()
            && overwrite && !decompose
            && nTargetProcs != -1 && nTargetProcs < Pstream::nProcs()
        )
        {
            Info<<"Deleting processor folder "<<nTargetProcs
                <<" to "<<Pstream::nProcs()-1<<endl;
            for (label proci=nTargetProcs; proci < Pstream::nProcs(); proci++)
            {
                fileName procPath
                (
                    args.path()/".."/(word("processor") + name(proci))
                );
                rmDir(procPath);
            }
        }

        runTime.setTime(origTime,origTimeIndex);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
