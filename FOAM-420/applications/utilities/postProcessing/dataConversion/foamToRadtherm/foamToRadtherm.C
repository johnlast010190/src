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
    foamToRadtherm

Description
    Generates Radtherm Patran neutral file

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurfaceTools/triSurfaceTools.H"

#include "db/IOstreams/IOstreams/IOmanip.H"
#include "indexedOctree/treeDataTriSurface.H"
#include "containers/Lists/ListListOps/ListListOps.H"
#include "primitives/Tuple2/Tuple2.H"

#include "cfdTools/general/include/fvCFD.H"

#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "psiThermo/psiThermo.H"

#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "polyTopoChange/polyTopoChange/modifyObject/polyModifyFace.H"

#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataPoint.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void splitBaffles
(
    const word& patchName,
    const labelList& baffleSplitFaces,
    polyMesh& mesh
)
{
    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    // Old and new patch.
    DynamicList<polyPatch*> allPatches(patches.size()+1);

    label startFaceI = mesh.nInternalFaces();

    // Copy old patches.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<directPolyPatch>(pp))
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(pp);
            if (!isA<processorPolyPatch>(dpp))
            {
                allPatches.append
                (
                    dpp.clone
                    (
                        patches,
                        patchI,
                        dpp.size(),
                        startFaceI
                     ).ptr()
                 );
                startFaceI += dpp.size();
            }
        }
    }

    word patchType("wall");
    label destPatchI = patches.findPatchID(patchName);

    if (destPatchI == -1)
    {
        destPatchI = allPatches.size();

        // Add an empty patch.
        allPatches.append
        (
            polyPatch::New
            (
                patchType,
                patchName,
                0,
                startFaceI,
                destPatchI,
                patches
             ).ptr()
         );
    }
    // Copy old processor patches.
    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];

        if (isA<directPolyPatch>(pp))
        {
            const directPolyPatch& dpp =
                refCast<const directPolyPatch>(pp);

            if (isA<processorPolyPatch>(dpp))
            {
                allPatches.append
                (
                    dpp.clone
                    (
                        patches,
                        patchI,
                        dpp.size(),
                        startFaceI
                     ).ptr()
                 );
                startFaceI += dpp.size();
            }
        }
    }

    allPatches.shrink();
    mesh.removeBoundary();
    mesh.addPatches(allPatches);

    polyTopoChange meshMod(mesh);
    label patchID = mesh.boundaryMesh().findPatchID(patchName);
    bool zoneFlip = false;

    forAll(baffleSplitFaces, i)
    {
        label meshFaceI = baffleSplitFaces[i];
        label zoneID = mesh.faceZones().whichZone(meshFaceI);

        meshMod.setAction
        (
            polyModifyFace
            (
                mesh.faces()[meshFaceI],            // face
                meshFaceI,                          // face ID
                mesh.faceOwner()[meshFaceI],        // owner
                -1,                                 // neighbour
                false,                              // flip flux
                patchID,                            // patch ID
                false,                              // remove from zone
                zoneID,                             // zone ID
                zoneFlip                            // zone flip
            )
        );
    }
    meshMod.changeMesh(mesh, false);
}


// Main program:

int main(int argc, char *argv[])
{
    argList::validOptions.insert("patches", "(patch0 .. patchN)");
    argList::validOptions.insert("file", "Radtherm file name");
    argList::validOptions.insert("compressible", "");

#include "include/addTimeOptions.H"
#include "include/setRootCase.H"
#include "include/createTime.H"
#include "include/createMesh.H"

    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#include "include/checkTimeOptions.H"
    runTime.setTime(Times[startTime], startTime);

    Info<<"Exporting time: "<<runTime.timeName()<<" to radTherm"<< nl << endl;

    //Read in Fields
#include "createFields.H"

    fileName radThermFileName("defaultRadThermFile.neu");

    IOobject radThermHeader
    (
        "radThermDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    DynamicList<Tuple2<label, vector>> baffleInfo;
    DynamicList<word > baffleNames;
    wordList patchNames;

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    if (radThermHeader.typeHeaderOk<IOdictionary>(true))
    {
        Info<<"Reading radThermDict"<<endl;
        // Read Radtherm dictionary
        IOdictionary radThermDict
        (
            IOobject
            (
                "radThermDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
         );

        if (radThermDict.found("patches"))
        {
            patchNames = wordList(radThermDict.lookup("patches"));
        }
        else if (args.optionFound("patches"))
        {
            patchNames = args.optionRead<wordList>("patches");
        }

        if (radThermDict.found("file"))
        {
            radThermFileName = fileName(radThermDict.lookup("file"));
        }
        else if (args.optionFound("file"))
        {
            radThermFileName = args.option("file");
        }

        const dictionary* bafflesDictPtr = radThermDict.subDictPtr
        (
            "baffles"
         );

        if (bafflesDictPtr)
        {
            const dictionary& bafflesDict = *bafflesDictPtr;

            forAll(bMesh, patchI)
            {
                const word& patchName = bMesh[patchI].name();
                if (bafflesDict.found(patchName))
                {
                    Info<<"Found Baffle patch: "<<patchName<<endl;
                    const dictionary& baffleDict =
                        bafflesDict.subDict(patchName);
                    pointField altDirection =
                        pointField(1, baffleDict.lookup("altDirection"));
                    baffleInfo.append
                    (
                        Tuple2<label, vector>
                        (
                            patchI,
                            altDirection[0]
                        )
                    );
                    baffleNames.append(patchName);
                }
            }
        }
    }
    else
    {
        if (args.optionFound("patches"))
        {
            patchNames = args.optionRead<wordList>("patches");
        }
        if (args.optionFound("file"))
        {
            radThermFileName = args.option("file");
        }
    }

    radThermFileName = radThermFileName.lessExt()
        + "_time"+runTime.timeName() +"."+radThermFileName.ext();

    baffleInfo.shrink();
    baffleNames.shrink();


    // Construct table of patches to include.
    labelHashSet includePatches(bMesh.size());
    if (patchNames.size())
    {
        forAll(patchNames, patchNameI)
        {
            const word& patchName = patchNames[patchNameI];

            label patchI = bMesh.findPatchID(patchName);

            if (patchI == -1)
            {
                FatalErrorIn(args.executable()) << "No such patch "
                    << patchName << endl << "Patches are " << bMesh.names()
                    << exit(FatalError);
            }
            includePatches.insert(patchI);
        }
    }
    else
    {
        forAll(bMesh, patchI)
        {
            const polyPatch& patch = bMesh[patchI];

            if (!isA<processorPolyPatch>(patch))
            {
                includePatches.insert(patchI);
            }
        }
    }
    HashSet<word> patchNamesToBaffle(baffleNames);

    Map<label> bafflePatchID(bMesh.size());
    if (baffleInfo.size())
    {
        forAll(baffleInfo, i)
        {
            label regionI = baffleInfo[i].first();
            bafflePatchID.insert(regionI, i);
        }
    }

    label nBoundaryFaces = mesh.nFaces() - mesh.nInternalFaces();
    DynamicList<scalar> patchHTC(nBoundaryFaces);
    DynamicList<scalar> patchOtherHTC(nBoundaryFaces);

    DynamicList<scalar> patchTemp(nBoundaryFaces);
    DynamicList<scalar> patchOtherTemp(nBoundaryFaces);

    DynamicList<point> patchFaceCentres(nBoundaryFaces);
    DynamicList<point> patchOtherFaceCentres(nBoundaryFaces);

    DynamicList<label>  baffleSplitFaces(nBoundaryFaces);

    forAllConstIter(labelHashSet, includePatches, iter)
    {
        label patchI = iter.key();
        const word patchName = bMesh[patchI].name();
        const bool isBaffle = patchNamesToBaffle.found(patchName);

        forAll(bMesh[patchI], j)
        {
            const label meshFaceI = bMesh[patchI].start() + j;
            const point fc = mesh.faceCentres()[meshFaceI];
            const vector fn = mesh.faceAreas()[meshFaceI];
            const scalar nwt =
                T.boundaryField()[patchI].patchInternalField()()[j];
            const scalar htc = alphaConv.boundaryField()[patchI][j];

            if (isBaffle)
            {
                const vector dir = baffleInfo[bafflePatchID[patchI]].second();

                if ((dir & fn) < 0.0)
                {
                    baffleSplitFaces.append(meshFaceI);
                    patchOtherFaceCentres.append(fc);
                    patchOtherTemp.append(nwt);
                    patchOtherHTC.append(htc);
                }
                else
                {
                    patchFaceCentres.append(fc);
                    patchTemp.append(nwt);
                    patchHTC.append(htc);
                }
            }
            else
            {
                   patchFaceCentres.append(fc);
                   patchTemp.append(nwt);
                   patchHTC.append(htc);
            }
        }
    }
    patchFaceCentres.shrink();
    patchTemp.shrink();
    patchHTC.shrink();
    patchOtherTemp.shrink();
    patchOtherHTC.shrink();
    patchOtherFaceCentres.shrink();

    baffleSplitFaces.shrink();

    if (Pstream::parRun())
    {
        // Gather all faceCentres on master
        List<pointField> gatheredPoints(Pstream::nProcs());
        // Gather all other side faceCentres on master
        List<pointField> gatheredOtherPoints(Pstream::nProcs());

        //Gather all htc on master
        List<scalarField> gatheredHTC(Pstream::nProcs());
        //Gather all temp on master
        List<scalarField> gatheredTemp(Pstream::nProcs());

        //Gather all other htc on master
        List<scalarField> gatheredOtherHTC(Pstream::nProcs());
        //Gather all other temp on master
        List<scalarField> gatheredOtherTemp(Pstream::nProcs());

        gatheredPoints[Pstream::myProcNo()].transfer(patchFaceCentres);
        Pstream::gatherList(gatheredPoints);

        gatheredOtherPoints[Pstream::myProcNo()].transfer(patchOtherFaceCentres);
        Pstream::gatherList(gatheredOtherPoints);

        gatheredHTC[Pstream::myProcNo()].transfer(patchHTC);
        Pstream::gatherList(gatheredHTC);

        gatheredTemp[Pstream::myProcNo()].transfer(patchTemp);
        Pstream::gatherList(gatheredTemp);

        gatheredOtherHTC[Pstream::myProcNo()].transfer(patchOtherHTC);
        Pstream::gatherList(gatheredOtherHTC);

        gatheredOtherTemp[Pstream::myProcNo()].transfer(patchOtherTemp);
        Pstream::gatherList(gatheredOtherTemp);

        if (Pstream::master())
        {
            /*
            // Count number of faces.
            label nFaces = 0;

            forAll(gatheredPoints, procI)
            {
                nFaces += gatheredPoints[procI].size();
            }
            */

            patchFaceCentres.clearStorage();
            patchTemp.clearStorage();
            patchHTC.clearStorage();
            patchOtherTemp.clearStorage();
            patchOtherHTC.clearStorage();
            patchOtherFaceCentres.clearStorage();

            for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {
                forAll(gatheredPoints[procI], j)
                {
                    patchFaceCentres.append(gatheredPoints[procI][j]);
                    patchTemp.append(gatheredTemp[procI][j]);
                    patchHTC.append(gatheredHTC[procI][j]);
                }
                forAll(gatheredOtherPoints[procI], j)
                {
                    patchOtherFaceCentres.append(gatheredOtherPoints[procI][j]);
                    patchOtherTemp.append(gatheredOtherTemp[procI][j]);
                    patchOtherHTC.append(gatheredOtherHTC[procI][j]);
                }

            }
            patchFaceCentres.shrink();
            patchOtherFaceCentres.shrink();

            patchTemp.shrink();
            patchOtherTemp.shrink();
            patchHTC.shrink();
            patchOtherHTC.shrink();
        }
    }

    if (returnReduce(baffleSplitFaces.size(), sumOp<label>()) != 0)
    {
        splitBaffles
        (
            "otherBaffleSide",
            baffleSplitFaces,
            mesh
         );
    }

    autoPtr<triSurface> allSurfPtr;
    if (!Pstream::parRun())
    {
        allSurfPtr.reset
        (
            new triSurface
            (
                triSurfaceTools::triangulate
                (
                    bMesh,
                    includePatches
                 )
            )
        );
    }
    else
    {
        triSurface localSurface
        (
            triSurfaceTools::triangulate
            (
                bMesh,
                includePatches
             )
         );

        // Gather all points on master
        List<pointField> gatheredPoints(Pstream::nProcs());

        gatheredPoints[Pstream::myProcNo()] = localSurface.points();

        Pstream::gatherList(gatheredPoints);

        // Gather all localSurface patches
        List<geometricSurfacePatchList> gatheredPatches(Pstream::nProcs());

        gatheredPatches[Pstream::myProcNo()] = localSurface.patches();

        Pstream::gatherList(gatheredPatches);

        // Gather all faces
        List<List<labelledTri>> gatheredFaces(Pstream::nProcs());

        gatheredFaces[Pstream::myProcNo()] = localSurface;

        Pstream::gatherList(gatheredFaces);

        if (Pstream::master())
        {
            // On master combine all points
            pointField allPoints =
                ListListOps::combine<pointField>
                (
                    gatheredPoints,
                    accessOp<pointField>()
                );

            // Count number of patches.
            label nPatches = 0;

            forAll(gatheredPatches, procI)
            {
                nPatches += gatheredPatches[procI].size();
            }

            // Count number of faces.
            label nFaces = 0;

            forAll(gatheredFaces, procI)
            {
                nFaces += gatheredFaces[procI].size();
            }

            // Loop over all processors and
            // - construct mapping from local to global patches
            // - relabel faces (both points and regions)

            label newPatchI = 0;

            // Name to new patchI
            HashTable<label> nameToIndex(2*nPatches);

            // Storage (oversized) for all patches
            geometricSurfacePatchList allPatches(nPatches);

            label newFaceI = 0;

            // Storage for all faces
            List<labelledTri> allFaces(nFaces);

            // Offset into allPoints for current processor
            label pointOffset = 0;

            for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {
                const geometricSurfacePatchList& patches =
                    gatheredPatches[procI];

                // From local patch numbering to global
                labelList localToGlobal(patches.size());

                forAll(patches, patchI)
                {
                    const geometricSurfacePatch& sp = patches[patchI];

                    if (!nameToIndex.found(sp.name()))
                    {
                        nameToIndex.insert(sp.name(), newPatchI);

                        localToGlobal[patchI] = newPatchI;

                        allPatches[newPatchI] = sp;

                        newPatchI++;
                    }
                    else
                    {
                        localToGlobal[patchI] = nameToIndex[sp.name()];
                    }
                }

                // Collect and relabel faces
                const List<labelledTri>& localFaces = gatheredFaces[procI];

                forAll(localFaces, faceI)
                {
                    const labelledTri& f = localFaces[faceI];

                    allFaces[newFaceI++] =
                        labelledTri
                        (
                            f[0] + pointOffset,
                            f[1] + pointOffset,
                            f[2] + pointOffset,
                            localToGlobal[f.region()]
                        );
                }

                pointOffset += gatheredPoints[procI].size();
            }
            allPatches.setSize(newPatchI);

            // We now have allPoints, allFaces and allPatches.
            // Construct overall (yet unmerged) surface from these.

            allSurfPtr.reset
            (
                new triSurface
                (
                    allFaces,
                    allPatches,
                    allPoints
                )
             );

            // Cleanup (which does point merge as well
            allSurfPtr().cleanup(false);
        }
    }
    Info<<"Created triangulated surface"<<endl;

    if (Pstream::master())
    {
        triSurface& allSurf = allSurfPtr();

        labelList patchFaceMap;
        const surfacePatchList surfacePatches =
            allSurf.calcPatches(patchFaceMap);

        //find surface baffle patches
        DynamicList<label> bRegions(allSurf.patches().size());
        forAll(allSurf.patches(), patchI)
        {
            word surfacePatchName = allSurf.patches()[patchI].name();
            if (patchNamesToBaffle.found(surfacePatchName))
            {
                bRegions.append(patchI);
            }
        }
        bRegions.shrink();
        labelHashSet baffleRegions(bRegions);

        if (Pstream::parRun())
        {
            radThermFileName
                =  runTime.rootPath()/runTime.caseName()/".."
                /"radTherm"/radThermFileName;
            if
            (
                !Foam::isDir
                (
                    runTime.rootPath()/runTime.caseName()/".."/"radTherm"
                )
            )
            {
                Foam::mkDir(runTime.rootPath()/runTime.caseName()
                            /".."/"radTherm");
            }
        }
        else
        {
            radThermFileName
                = runTime.rootPath()/runTime.caseName()
                /"radTherm"/radThermFileName;
            if
            (
                !Foam::isDir
                (
                    runTime.rootPath()/runTime.caseName()/"radTherm"
                 )
            )
            {
                Foam::mkDir(runTime.rootPath()/runTime.caseName()
                            /"radTherm");
            }
        }
        autoPtr<OFstream> radThermFilePtr(new OFstream(radThermFileName));

        string dateString = clock::date();
        string timeString = clock::clockTime();

        //write packet 25 (File Title)
        radThermFilePtr()<<"25"<<"       0"<<"       0"<<"       1"<<"       0"
                         <<"       0"<<"       0"<<"       0"<<"       0"<<endl;
        radThermFilePtr()<<setw(25)<<"OpenFOAM RadTherm Export "<<endl;

        //write packet 26 (Summary data)
        radThermFilePtr()<<"26"<<"       0"<<"       0"<<"       1"
                         << setw(8) << allSurf.nPoints()
                         << setw(8) << allSurf.size()
                         <<"       0"<<"       0"<<"       0"<<endl;
        radThermFilePtr() << setw(12) << dateString.c_str()
                          << setw(12) << timeString.c_str()
                          <<"   2.1"<< endl;

        //write packet 1 (Node data)
        forAll(allSurf.localPoints(), i)
        {
            const point pt =  allSurf.localPoints()[i];
            radThermFilePtr()<<" 1"<< setw(8) << i + 1 <<"       0"<<"       2"
                             <<"       0"<<"       0"<<"       0"<<"       0"
                             <<"       0"<<endl;
            radThermFilePtr()<< scientific << setw(16) << setprecision(9)
                             << pt.x() << setw(16) << pt.y()
                             << setw(16) << pt.z() <<endl;
            radThermFilePtr()<<"1G"<<"       3"<<"       0"
                             <<"       0"<<"  000000"<<endl;
        }
        //write packet 2 (Element data)
        forAll(allSurf.localFaces(), i)
        {
            const labelledTri f = allSurf.localFaces()[patchFaceMap[i]];
            label region = f.region();

            const point orient = point::zero;

            radThermFilePtr()<<" 2"<< setw(8) << i + 1 <<"       3"<<"       2"
                             <<"       0"<<"       0"<<"       0"<<"       0"
                             <<"       0"<<endl;
            radThermFilePtr()<<"       3"<<"       0"<< setw(8)
                             << region + 1 <<"       0"
                             << scientific << setw(16) << setprecision(9)
                             << orient.x()
                             << setw(16) << orient.y()
                             << setw(16) << orient.z() <<endl;
            //Radtherm requires reverse face orientation
            triFace rFace = f.reverseFace();
            forAll(rFace, fI)
            {
                radThermFilePtr()<< setw(8) << rFace[fI] + 1;
            }
            radThermFilePtr()<<endl;
        }

        //write packet 21 (Component data)
        forAll(surfacePatches, spI)
        {
            const surfacePatch& p = surfacePatches[spI];
            label kc = 1 + (2*p.size()+9)/10;

            radThermFilePtr()<<"21"<< setw(8) << spI + 1 << setw(8)
                             << 2*p.size() <<setw(8) << kc <<"       0"
                             <<"       0"<<"       0"<<"       0"
                             <<"       0"<<endl;
            radThermFilePtr()<< p.name() <<endl;

            label current = p.start();
            forAll(p, i)
            {
                radThermFilePtr()<<"       7"<< setw(8) << current + 1;
                if (((i+1) % 5) == 0 || i == p.size() -1)
                {
                    radThermFilePtr()<<endl;
                }
                current++;
            }
        }

        //construct octree
        treeBoundBox overallBb(allSurf.faceCentres());
        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        pointField pfc(patchFaceCentres.xfer());
        pointField pofc(patchOtherFaceCentres.xfer());

        indexedOctree<treeDataPoint> pointTree
        (
            treeDataPoint(pfc),
            overallBb,  // overall search domain
            8,                              // maxLevel
            10,                             // leafsize
            3.0                             // duplicity
        );

        //loop over points, find nearest and assign mapped data
        scalar spanSqr = Foam::sqr(pointTree.bb().mag());

        autoPtr<indexedOctree<treeDataPoint>> pointOtherSideTree;
        if (pofc.size())
        {
            pointOtherSideTree.reset
            (
                new indexedOctree<treeDataPoint>
                (
                    treeDataPoint(pofc),
                    overallBb,  // overall search domain
                    8,                              // maxLevel
                    10,                             // leafsize
                    3.0                             // duplicity
                 )
             );
        }

        const List<labelledTri>& localFaces = allSurf.localFaces();

        scalarField wallHTC(localFaces.size());
        DynamicList<scalar> wallHTCOther(localFaces.size());
        scalarField wallTemp(localFaces.size());
        DynamicList<scalar> wallTempOther(localFaces.size());

        forAll(localFaces, i)
        {
            label pfI = patchFaceMap[i];
            const labelledTri& tri = localFaces[pfI];
            const label region = tri.region();

            const point fc =  allSurf.faceCentres()[pfI];
            pointIndexHit info = pointTree.findNearest
            (
                fc,
                spanSqr
            );

            if (info.hit())
            {
                wallHTC[pfI] = patchHTC[info.index()] ;
                wallTemp[pfI] = patchTemp[info.index()];
            }
            else
            {
                wallHTC[pfI] = 0.0;
                wallTemp[pfI] = 273.15;
            }

            if (baffleRegions.found(region))
            {
                pointIndexHit infoOtherSide = pointOtherSideTree().findNearest
                (
                    fc,
                    spanSqr
                );

                if (infoOtherSide.hit())
                {
                    wallHTCOther.append(patchOtherHTC[infoOtherSide.index()]);
                    wallTempOther.append(patchOtherTemp[infoOtherSide.index()]);
                }
                else
                {
                    wallHTCOther.append(0.0);
                    wallTempOther.append(273.15);
                }

            }
        }
        wallTempOther.shrink();
        wallHTCOther.shrink();

        scalar defaultVal = -999.;
        //write packet 17 (HTC data)
        {
            label nBaffleFaces = 0;
            forAll(localFaces, i)
            {
                label pfI = patchFaceMap[i];
                const labelledTri& tri = localFaces[pfI];
                const label region = tri.region();

                radThermFilePtr()<<"17"<< setw(8) << i + 1
                                 << setw(8)  << region + 1
                                 << "       2" <<"       1"<<"       1"
                                 <<"       0"<<"       0"<<"       0"<<endl;

                if (baffleRegions.found(region))
                {
                    radThermFilePtr()<<"1 "<< "11000000" <<endl;
                    radThermFilePtr()<< setw(16) << setprecision(9)
                                     << wallHTC[pfI]
                                     << setw(16)
                                     << wallHTCOther[nBaffleFaces] <<endl;
                    nBaffleFaces++;
                }
                else
                {
                    radThermFilePtr()<<"1 "<< "10000000" <<endl;
                    radThermFilePtr()<< setw(16) << setprecision(9)
                                     << wallHTC[pfI]
                                     << setw(16) << defaultVal <<endl;
                }
            }
        }

        //write packet 18 (Near Wall Temperature data)
       {
            label nBaffleFaces = 0;
            forAll(localFaces, i)
            {
                label pfI = patchFaceMap[i];
                const labelledTri& tri = localFaces[pfI];
                const label region = tri.region();

                radThermFilePtr()<<"18"<< setw(8) << i + 1
                                 << setw(8)  << region + 1
                                 << "       2" <<"       1"<<"       1"
                                 <<"       0"<<"       0"<<"       0"<<endl;
                radThermFilePtr()<<"1 "<< "00000000" <<endl;

                if (baffleRegions.found(region))
                {
                    radThermFilePtr()<< setw(16) << setprecision(9)
                                     << wallTemp[pfI]
                                     << setw(16) << wallTempOther[nBaffleFaces]
                                     << setw(16) << defaultVal
                                     << setw(16) << defaultVal <<endl;
                    nBaffleFaces++;
                }
                else
                {
                    radThermFilePtr()<< setw(16) << setprecision(9)
                                     << wallTemp[pfI]
                                     << setw(16) << defaultVal
                                     << setw(16) << defaultVal
                                     << setw(16) << defaultVal <<endl;
                }
            }
        }

        //write packet 99 (End of File)
        radThermFilePtr()<<"99"<<"       0"<<"       0"<<"       1"<<"       0"
                         <<"       0"<<"       0"<<"       0"<<"       0"<<endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
