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
    (c) 2016 OpenCFD Ltd.
    (c) 2017-2023 Esi Ltd.

Application
    renumberMesh

Group
    grpMeshManipulationUtilities

Description
    Renumbers the cell list in order to reduce the bandwidth, reading and
    renumbering all fields from all the time directories.

    By default uses bandCompression (CuthillMcKee) but will
    read system/renumberMeshDict if -dict option is present

\*---------------------------------------------------------------------------*/

#include "global/argList/argList.H"
#include "db/IOobjectList/IOobjectList.H"
#include "fvMesh/fvMesh.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "fields/ReadFields/ReadFields.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "containers/Lists/SortableList/SortableList.H"
#include "decompositionMethod/decompositionMethod.H"
#include "renumberMethod/renumberMethod.H"
#include "fields/fvPatchFields/basic/zeroGradient/zeroGradientFvPatchFields.H"
#include "CuthillMcKeeRenumber/CuthillMcKeeRenumber.H"
#include "fvMeshSubset/fvMeshSubset.H"
#include "sets/topoSets/cellSet.H"
#include "sets/topoSets/faceSet.H"
#include "sets/topoSets/pointSet.H"
#include "processorMeshes.H"
#include "regionProperties/regionProperties.H"
#include "polyTopoChange/hexRef8/hexRef8.H"
#include "renumberUtils.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Create named field from labelList for postprocessing
tmp<volScalarField> createScalarField
(
    const fvMesh& mesh,
    const word& name,
    const labelList& elems
)
{
    tmp<volScalarField> tfld
    (
        new volScalarField
        (
            IOobject
            (
                name,
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar("zero", dimless, 0),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& fld = tfld.ref();

    forAll(fld, celli)
    {
       fld[celli] = elems[celli];
    }

    return tfld;
}


// Calculate band of matrix
label getBand(const labelList& owner, const labelList& neighbour)
{
    label band = 0;

    forAll(neighbour, facei)
    {
        label diff = neighbour[facei] - owner[facei];

        if (diff > band)
        {
            band = diff;
        }
    }
    return band;
}


// Calculate band of matrix
void getBand
(
    const bool calculateIntersect,
    const label nCells,
    const labelList& owner,
    const labelList& neighbour,
    label& bandwidth,
    scalar& profile,            // scalar to avoid overflow
    scalar& sumSqrIntersect     // scalar to avoid overflow
)
{
    labelList cellBandwidth(nCells, 0);
    scalarField nIntersect(nCells, 0.0);

    forAll(neighbour, facei)
    {
        label own = owner[facei];
        label nei = neighbour[facei];

        // Note: mag not necessary for correct (upper-triangular) ordering.
        label diff = nei-own;
        cellBandwidth[nei] = max(cellBandwidth[nei], diff);
    }

    bandwidth = max(cellBandwidth);

    // Do not use field algebra because of conversion label to scalar
    profile = 0.0;
    forAll(cellBandwidth, celli)
    {
        profile += 1.0*cellBandwidth[celli];
    }

    sumSqrIntersect = 0.0;
    if (calculateIntersect)
    {
        forAll(nIntersect, celli)
        {
            for (label colI = celli-cellBandwidth[celli]; colI <= celli; colI++)
            {
                nIntersect[colI] += 1.0;
            }
        }

        sumSqrIntersect = sum(Foam::sqr(nIntersect));
    }
}


// Determine upper-triangular face order
labelList getFaceOrder
(
    const primitiveMesh& mesh,
    const labelList& cellOrder      // New to old cell
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));

    labelList oldToNewFace(mesh.nFaces(), -1);

    label newFacei = 0;

    labelList nbr;
    labelList order;

    forAll(cellOrder, newCelli)
    {
        label oldCelli = cellOrder[newCelli];

        const cell& cFaces = mesh.cells()[oldCelli];

        // Neighbouring cells
        nbr.setSize(cFaces.size());

        forAll(cFaces, i)
        {
            label facei = cFaces[i];

            if (mesh.isInternalFace(facei))
            {
                // Internal face. Get cell on other side.
                label nbrCelli = reverseCellOrder[mesh.faceNeighbour()[facei]];
                if (nbrCelli == newCelli)
                {
                    nbrCelli = reverseCellOrder[mesh.faceOwner()[facei]];
                }

                if (newCelli < nbrCelli)
                {
                    // Celli is master
                    nbr[i] = nbrCelli;
                }
                else
                {
                    // nbrCell is master. Let it handle this face.
                    nbr[i] = -1;
                }
            }
            else
            {
                // External face. Do later.
                nbr[i] = -1;
            }
        }

        order.setSize(nbr.size());
        sortedOrder(nbr, order);

        forAll(order, i)
        {
            label index = order[i];
            if (nbr[index] != -1)
            {
                oldToNewFace[cFaces[index]] = newFacei++;
            }
        }
    }

    // Leave patch faces intact.
    for (label facei = newFacei; facei < mesh.nFaces(); facei++)
    {
        oldToNewFace[facei] = facei;
    }


    // Check done all faces.
    forAll(oldToNewFace, facei)
    {
        if (oldToNewFace[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position" << " for face " << facei
                << abort(FatalError);
        }
    }

    return invert(mesh.nFaces(), oldToNewFace);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Renumber mesh to minimise bandwidth"
    );

    #include "include/addRegionOption.H"
    argList::addBoolOption
    (
        "allRegions",
        "operate on all regions in regionProperties"
    );
    #include "include/addOverwriteOption.H"
    #include "include/addTimeOptions.H"
    #include "include/addDictOption.H"
    argList::addBoolOption
    (
        "frontWidth",
        "calculate the rms of the frontwidth"
    );
    argList::addBoolOption
    (
        "noFields",
        "do not update fields"
    );

    #include "include/setRootCase.H"
    #include "include/createTime.H"
    runTime.functionObjects().off();

    // Get times list
    instantList Times = runTime.times();

    // Set startTime and endTime depending on -time and -latestTime options
    #include "include/checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);
    const bool allRegions = args.optionFound("allRegions");

    wordList regionNames;
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
                }
            }
        }
    }
    else
    {
        word regionName;
        if (args.optionReadIfPresent("region", regionName))
        {
            regionNames = wordList(1, regionName);
        }
        else
        {
            regionNames = wordList(1, fvMesh::defaultRegion);
        }
    }

    forAll(regionNames, regionI)
    {
        runTime.setTime(Times[startTime], startTime);
        const word& regionName = regionNames[regionI];

        Info<< "\n\n Renumbering for mesh " << regionName << nl
            << endl;

        Foam::fvMesh mesh
        (
            Foam::IOobject
            (
                regionName,
                runTime.timeName(),
                runTime,
                Foam::IOobject::MUST_READ
            )
        );

        const word oldInstance = mesh.pointsInstance();

        const bool readDict = args.optionFound("dict");
        const bool doFrontWidth = args.optionFound("frontWidth");
        const bool overwrite = args.optionFound("overwrite");
        const bool fields = !args.optionFound("noFields");

        label band;
        scalar profile;
        scalar sumSqrIntersect;
        getBand
        (
            doFrontWidth,
            mesh.nCells(),
            mesh.faceOwner(),
            mesh.faceNeighbour(),
            band,
            profile,
            sumSqrIntersect
        );

        scalar totalSqrIntersect = sumSqrIntersect;
        reduce(
            std::tie(band, profile, totalSqrIntersect),
            ParallelOp<maxOp<label>, sumOp<scalar>, sumOp<scalar>>{}
        );

        scalar rmsFrontwidth = Foam::sqrt
        (
            totalSqrIntersect / mesh.globalData().nTotalCells()
        );

        Info<< "Mesh size: " << mesh.globalData().nTotalCells() << nl
            << "Before renumbering :" << nl
            << "    band           : " << band << nl
            << "    profile        : " << profile << nl;

        if (doFrontWidth)
        {
            Info<< "    rms frontwidth : " << rmsFrontwidth << nl;
        }

        Info<< endl;

        bool sortCoupledFaceCells = false;
        bool writeMaps = false;
        bool orderPoints = false;
        label blockSize = 0;

        // Construct renumberMethod
        autoPtr<IOdictionary> renumberDictPtr;
        autoPtr<renumberMethod> renumberPtr;

        if (readDict)
        {
            const word dictName("renumberMeshDict");
            #include "include/setSystemMeshDictionaryIO.H"

            Info<< "Renumber according to " << dictName << nl << endl;

            renumberDictPtr.reset(new IOdictionary(dictIO));
            const IOdictionary& renumberDict = renumberDictPtr();

            renumberPtr = renumberMethod::New(renumberDict);

            sortCoupledFaceCells = renumberDict.lookupOrDefault
            (
                "sortCoupledFaceCells",
                false
            );
            if (sortCoupledFaceCells)
            {
                Info<< "Sorting cells on coupled boundaries to be last." << nl
                    << endl;
            }

            blockSize = renumberDict.lookupOrDefault("blockSize", 0);
            if (blockSize > 0)
            {
                Info<< "Ordering cells into regions of size " << blockSize
                    << " (using decomposition);"
                    << " ordering faces into region-internal and"
                    << " region-external."
                    << nl << endl;

                if (blockSize < 0 || blockSize >= mesh.nCells())
                {
                    FatalErrorInFunction
                        << "Block size " << blockSize
                        << " should be positive integer"
                        << " and less than the number of cells in the mesh."
                        << exit(FatalError);
                }
            }

            orderPoints = renumberDict.lookupOrDefault("orderPoints", false);
            if (orderPoints)
            {
                Info<< "Ordering points into internal and boundary points."
                    << nl
                    << endl;
            }

            renumberDict.lookup("writeMaps") >> writeMaps;
            if (writeMaps)
            {
                Info<< "Writing renumber maps (new to old) to polyMesh." << nl
                    << endl;
            }
        }
        else
        {
            Info<< "Using default renumberMethod." << nl << endl;
            dictionary renumberDict;
            renumberPtr.reset(new CuthillMcKeeRenumber(renumberDict));
        }

        Info<< "Selecting renumberMethod "
            << renumberPtr().type() << nl << endl;

        // Read parallel reconstruct maps
        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
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
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            labelList(0)
        );
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                mesh.pointsInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            labelList(0)
        );

        // Read objects in time directory
        IOobjectList objects(mesh, runTime.timeName());

        if (fields) Info<< "Reading geometric fields" << nl << endl;

        PtrList<volScalarField> vsFlds;
        if (fields) ReadFields(mesh, objects, vsFlds);

        PtrList<volVectorField> vvFlds;
        if (fields) ReadFields(mesh, objects, vvFlds);

        PtrList<volSphericalTensorField> vstFlds;
        if (fields) ReadFields(mesh, objects, vstFlds);

        PtrList<volSymmTensorField> vsymtFlds;
        if (fields) ReadFields(mesh, objects, vsymtFlds);

        PtrList<volTensorField> vtFlds;
        if (fields) ReadFields(mesh, objects, vtFlds);


        PtrList<surfaceScalarField> ssFlds;
        if (fields) ReadFields(mesh, objects, ssFlds);

        PtrList<surfaceVectorField> svFlds;
        if (fields) ReadFields(mesh, objects, svFlds);

        PtrList<surfaceSphericalTensorField> sstFlds;
        if (fields) ReadFields(mesh, objects, sstFlds);

        PtrList<surfaceSymmTensorField> ssymtFlds;
        if (fields) ReadFields(mesh, objects, ssymtFlds);

        PtrList<surfaceTensorField> stFlds;
        if (fields) ReadFields(mesh, objects, stFlds);


        PtrList<pointScalarField> psFlds;
        if (fields) ReadFields(pointMesh::New(mesh), objects, psFlds);

        PtrList<pointVectorField> pvFlds;
        if (fields) ReadFields(pointMesh::New(mesh), objects, pvFlds);

        PtrList<pointSphericalTensorField> pstFlds;
        if (fields) ReadFields(pointMesh::New(mesh), objects, pstFlds);

        PtrList<pointSymmTensorField> psymtFlds;
        if (fields) ReadFields(pointMesh::New(mesh), objects, psymtFlds);

        PtrList<pointTensorField> ptFlds;
        if (fields) ReadFields(pointMesh::New(mesh), objects, ptFlds);


        // Read sets
        PtrList<cellSet> cellSets;
        PtrList<faceSet> faceSets;
        PtrList<pointSet> pointSets;
        {
            // Read sets
            IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
                ReadFields(objects, cellSets);
            ReadFields(objects, faceSets);
            ReadFields(objects, pointSets);
        }


        Info<< endl;

        // From renumbering:
        // - from new cell/face back to original cell/face
        labelList cellOrder;
        labelList faceOrder;

        if (blockSize > 0)
        {
            // Renumbering in two phases. Should be done in one so mapping of
            // fields is done correctly!

            label nBlocks = mesh.nCells()/blockSize;
            Info<< "nBlocks   = " << nBlocks << endl;

            // Read decompositionMethod dictionary
            dictionary decomposeDict(renumberDictPtr().subDict("blockCoeffs"));
            decomposeDict.set("numberOfSubdomains", nBlocks);

            bool oldParRun = UPstream::parRun();
            UPstream::parRun() = false;

            autoPtr<decompositionMethod> decomposePtr = decompositionMethod::New
            (
                decomposeDict
            );

            labelList cellToRegion
            (
                decomposePtr().decompose
                (
                    mesh,
                    mesh.cellCentres()
                )
            );

            // Restore state
            UPstream::parRun() = oldParRun;

            // For debugging: write out region
            createScalarField
            (
                mesh,
                "cellDist",
                cellToRegion
            )().write();

            Info<< nl << "Written decomposition as volScalarField to "
                << "cellDist for use in postprocessing."
                << nl << endl;


            cellOrder = regionRenumber(renumberPtr(), mesh, cellToRegion);

            // Determine new to old face order with new cell numbering
            faceOrder = getRegionFaceOrder
            (
                mesh,
                cellOrder,
                cellToRegion
            );
        }
        else
        {
            // Detemines sorted back to original cell ordering
            cellOrder = renumberPtr().renumber
            (
                mesh,
                mesh.cellCentres()
            );

            if (sortCoupledFaceCells)
            {
                // Change order so all coupled patch faceCells are at the end.
                const polyBoundaryMesh& pbm = mesh.boundaryMesh();

                // Collect all boundary cells on coupled patches
                label nBndCells = 0;
                forAll(pbm, patchi)
                {
                    if (pbm[patchi].coupled())
                    {
                        nBndCells += pbm[patchi].size();
                    }
                }

                labelList reverseCellOrder = invert(mesh.nCells(), cellOrder);

                labelList bndCells(nBndCells);
                labelList bndCellMap(nBndCells);
                nBndCells = 0;
                forAll(pbm, patchi)
                {
                    if (pbm[patchi].coupled())
                    {
                        const labelUList& faceCells = pbm[patchi].faceCells();
                        forAll(faceCells, i)
                        {
                            label celli = faceCells[i];

                            if (reverseCellOrder[celli] != -1)
                            {
                                bndCells[nBndCells] = celli;
                                bndCellMap[nBndCells++] = reverseCellOrder[celli];
                                reverseCellOrder[celli] = -1;
                            }
                        }
                    }
                }
                bndCells.setSize(nBndCells);
                bndCellMap.setSize(nBndCells);

                // Sort
                labelList order;
                sortedOrder(bndCellMap, order);

                // Redo newReverseCellOrder
                labelList newReverseCellOrder(mesh.nCells(), -1);

                label sortedI = mesh.nCells();
                forAllReverse(order, i)
                {
                    label origCelli = bndCells[order[i]];
                    newReverseCellOrder[origCelli] = --sortedI;
                }

                Info<< "Ordered all " << nBndCells
                    << " cells with a coupled face"
                    << " to the end of the cell list, starting at " << sortedI
                    << endl;

                // Compact
                sortedI = 0;
                forAll(cellOrder, newCelli)
                {
                    label origCelli = cellOrder[newCelli];
                    if (newReverseCellOrder[origCelli] == -1)
                    {
                        newReverseCellOrder[origCelli] = sortedI++;
                    }
                }

                // Update sorted back to original (unsorted) map
                cellOrder = invert(mesh.nCells(), newReverseCellOrder);
            }


            // Determine new to old face order with new cell numbering
            faceOrder = getFaceOrder
            (
                mesh,
                cellOrder      // New to old cell
            );
        }


        if (!overwrite)
        {
            runTime++;
        }


        // Change the mesh.
        autoPtr<mapPolyMesh> map = reorderMesh(mesh, cellOrder, faceOrder);


        if (orderPoints)
        {
            polyTopoChange meshMod(mesh);
            autoPtr<mapPolyMesh> pointOrderMap = meshMod.changeMesh
            (
                mesh,
                false,      // inflate
                true,       // syncParallel
                false,      // orderCells
                orderPoints // orderPoints
            );

            // Combine point reordering into map.
            const_cast<labelList&>(map().pointMap()) = UIndirectList<label>
            (
                map().pointMap(),
                pointOrderMap().pointMap()
            )();

            inplaceRenumber
            (
                pointOrderMap().reversePointMap(),
                const_cast<labelList&>(map().reversePointMap())
             );
        }


        // Update fields
        mesh.updateMesh(map);

        // Update proc maps
        if (cellProcAddressing.headerOk())
        {
            if (cellProcAddressing.size() == mesh.nCells())
            {
                Info<< "Renumbering processor cell decomposition map "
                    << cellProcAddressing.name() << endl;

                cellProcAddressing = labelList
                (
                    UIndirectList<label>(cellProcAddressing, map().cellMap())
                    );
            }
            else
            {
                Info<< "Not writing inconsistent processor cell decomposition"
                        << " map " << cellProcAddressing.filePath() << endl;
                cellProcAddressing.writeOpt() = IOobject::NO_WRITE;
            }
        }
        else
        {
            cellProcAddressing.writeOpt() = IOobject::NO_WRITE;
        }

        if (faceProcAddressing.headerOk())
        {
            if (faceProcAddressing.size() == mesh.nFaces())
            {
                Info<< "Renumbering processor face decomposition map "
                    << faceProcAddressing.name() << endl;

                faceProcAddressing = labelList
                (
                        UIndirectList<label>(faceProcAddressing, map().faceMap())
                );

                // Detect any flips.
                const labelHashSet& fff = map().flipFaceFlux();
                forAllConstIter(labelHashSet, fff, iter)
                {
                    label facei = iter.key();
                    label masterFacei = faceProcAddressing[facei];

                    faceProcAddressing[facei] = -masterFacei;

                    if (masterFacei == 0)
                    {
                        FatalErrorInFunction
                            << " masterFacei:" << masterFacei
                            << exit(FatalError);
                    }
                }
            }
            else
            {
                Info<< "Not writing inconsistent processor face decomposition"
                    << " map " << faceProcAddressing.filePath() << endl;
                faceProcAddressing.writeOpt() = IOobject::NO_WRITE;
            }
        }
        else
        {
            faceProcAddressing.writeOpt() = IOobject::NO_WRITE;
        }

        if (pointProcAddressing.headerOk())
        {
            if (pointProcAddressing.size() == mesh.nPoints())
            {
                Info<< "Renumbering processor point decomposition map "
                    << pointProcAddressing.name() << endl;

                pointProcAddressing = labelList
                (
                    UIndirectList<label>(pointProcAddressing, map().pointMap())
                );
            }
            else
            {
                Info<< "Not writing inconsistent processor point decomposition"
                    << " map " << pointProcAddressing.filePath() << endl;
                pointProcAddressing.writeOpt() = IOobject::NO_WRITE;
            }
        }
        else
        {
            pointProcAddressing.writeOpt() = IOobject::NO_WRITE;
        }

        // Move mesh (since morphing might not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }


        {
            label band;
            scalar profile;
            scalar sumSqrIntersect;
            getBand
            (
                doFrontWidth,
                mesh.nCells(),
                mesh.faceOwner(),
                mesh.faceNeighbour(),
                band,
                profile,
                sumSqrIntersect
            );

            scalar totalSqrIntersect = sumSqrIntersect;
            reduce(
                std::tie(band, profile, totalSqrIntersect),
                ParallelOp<maxOp<label>, sumOp<scalar>, sumOp<scalar>>{}
            );

            scalar rmsFrontwidth = Foam::sqrt
            (
                totalSqrIntersect / mesh.globalData().nTotalCells()
            );

            Info<< "After renumbering :" << nl
                << "    band           : " << band << nl
                << "    profile        : " << profile << nl;

            if (doFrontWidth)
            {

                Info<< "    rms frontwidth : " << rmsFrontwidth << nl;
            }

            Info<< endl;
        }

        if (orderPoints)
        {
            // Force edge calculation (since only reason that points
            // would need to be sorted)
            (void)mesh.edges();

            label nTotPoints = returnReduce
            (
                mesh.nPoints(),
                sumOp<label>()
            );
            label nTotIntPoints = returnReduce
            (
                mesh.nInternalPoints(),
                sumOp<label>()
            );

            label nTotEdges = returnReduce
            (
                mesh.nEdges(),
                sumOp<label>()
            );
            label nTotIntEdges = returnReduce
            (
                mesh.nInternalEdges(),
                sumOp<label>()
            );
            label nTotInt0Edges = returnReduce
            (
                mesh.nInternal0Edges(),
                sumOp<label>()
            );
            label nTotInt1Edges = returnReduce
            (
                mesh.nInternal1Edges(),
                sumOp<label>()
             );

            Info<< "Points:" << nl
                << "    total   : " << nTotPoints << nl
                << "    internal: " << nTotIntPoints << nl
                << "    boundary: " << nTotPoints-nTotIntPoints << nl
                << "Edges:" << nl
                << "    total   : " << nTotEdges << nl
                << "    internal: " << nTotIntEdges << nl
                << "        internal using 0 boundary points: "
                << nTotInt0Edges << nl
                << "        internal using 1 boundary points: "
                << nTotInt1Edges-nTotInt0Edges << nl
                << "        internal using 2 boundary points: "
                << nTotIntEdges-nTotInt1Edges << nl
                << "    boundary: " << nTotEdges-nTotIntEdges << nl
                << endl;
        }

        if (overwrite)
        {
            mesh.setInstance(oldInstance);
        }
        else
        {
            mesh.setInstance(runTime.timeName());
        }


        Info<< "Writing mesh to " << mesh.facesInstance() << endl;

        // Remove old procAddressing files
        processorMeshes::removeFiles(mesh);

        // Update sets
        topoSet::updateMesh(mesh.facesInstance(), map(), cellSets);
        topoSet::updateMesh(mesh.facesInstance(), map(), faceSets);
        topoSet::updateMesh(mesh.facesInstance(), map(), pointSets);

        mesh.write();

        if (writeMaps)
        {
            // For debugging: write out region
            createScalarField
            (
                mesh,
                "origCellID",
                map().cellMap()
             )().write();

            createScalarField
            (
                mesh,
                "cellID",
                identity(mesh.nCells())
             )().write();

            Info<< nl << "Written current cellID and origCellID as"
                << " volScalarField for use in postprocessing."
                << nl << endl;

            labelIOList
            (
                IOobject
                (
                    "cellMap",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                map().cellMap()
            ).write();

            labelIOList
            (
                IOobject
                (
                    "faceMap",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                map().faceMap()
            ).write();

            labelIOList
            (
                IOobject
                (
                    "pointMap",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                 ),
                map().pointMap()
            ).write();
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
