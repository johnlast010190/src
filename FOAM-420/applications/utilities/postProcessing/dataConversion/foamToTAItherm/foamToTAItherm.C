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
    foamToTAItherm

Description
    Generates TAITherm Patran neutral file with HTC and Tfluid

    foamToTAItherm -compressible -file taiResults.ntl
    -backPatches (heatshield heatshield_slave tank tank_slave)

    -file taiResults.ntl contains the target mesh
    -backPatches is used to specify patches with front & back, i.e. baffles

\*---------------------------------------------------------------------------*/

#include "triSurface/triSurfaceTools/triSurfaceTools.H"

#include "db/IOstreams/IOstreams/IOmanip.H"
#include "indexedOctree/treeDataTriSurface.H"
#include "containers/Lists/ListListOps/ListListOps.H"

#include "cfdTools/general/include/fvCFD.H"

#include "singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "turbulentTransportModels/turbulentTransportModel.H"
#include "turbulentFluidThermoModels/turbulentFluidThermoModel.H"
#include "rhoThermo/rhoThermo.H"

#include "algorithms/indexedOctree/indexedOctree.H"
#include "indexedOctree/treeDataPoint.H"

#include "fields/fvPatchFields/derived/mappedFromFile/mappedFromFileFvPatchField.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Read a column of a given width from a fixed-format NTL file
std::string readNTLtoken
(
    const string& line,
    const size_t& width,
    size_t& index
)
{
    size_t indexStart, indexEnd;

    indexStart = index;

    indexEnd = line.find(',', indexStart);
    index = indexEnd + 1;

    if (indexEnd == std::string::npos)
    {
        indexEnd = indexStart + width;
        index = indexEnd;
    }

    return line.substr(indexStart, indexEnd - indexStart);
}

int main(int argc, char *argv[])
{
    argList::validOptions.insert("patches", "list of patches (patch0 .. patchN)");
    argList::validOptions.insert("backPatches", "list of backPatches (patch0 .. patchN)");
    argList::validOptions.insert("file", "TAItherm file name");
    argList::validOptions.insert("compressible", "");

    argList::addOption
    (
        "TRef",
        "scalar",
        "Use reference temperature to calculate HTC"
    );

#include "include/addTimeOptions.H"
#include "include/setRootCase.H"
#include "include/createTime.H"
#include "include/createMesh.H"

    scalar totTimeClo = runTime.elapsedClockTime();
    scalar totTimeCPU = runTime.elapsedCpuTime();

    //Read in Fields
#include "createFieldsBuoyant.H"

    fileName TAIthermFileName("taiResults.ntl");

    IOobject TAIthermHeader
    (
        "TAIthermDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );

    wordList patchNames;
    wordList backPatchNames;

    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    if (TAIthermHeader.typeHeaderOk<IOdictionary>(true))
    {
        Info<<"Reading TAIthermDict"<<endl;
        // Read Radtherm dictionary
        IOdictionary TAIthermDict
        (
            IOobject
            (
                "TAIthermDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
         );

        if (TAIthermDict.found("patches"))
        {
            patchNames = wordList(TAIthermDict.lookup("patches"));
        }
        else if (args.optionFound("patches"))
        {
            patchNames = args.optionRead<wordList>("patches");
        }

        if (TAIthermDict.found("backPatches"))
        {
            backPatchNames = wordList(TAIthermDict.lookup("backPatches"));
        }
        else if (args.optionFound("backPatches"))
        {
            backPatchNames = args.optionRead<wordList>("backPatches");
        }


        if (TAIthermDict.found("file"))
        {
            TAIthermFileName = fileName(TAIthermDict.lookup("file"));
        }
        else if (args.optionFound("file"))
        {
            TAIthermFileName = args.option("file");
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
            TAIthermFileName = args.option("file");
        }
        if (args.optionFound("backPatches"))
        {
            backPatchNames = args.optionRead<wordList>("backPatches");
        }
    }

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

    // Construct table of patches to include as back.
    labelHashSet includeBackPatches(bMesh.size());
    if (backPatchNames.size())
    {
        forAll(backPatchNames, patchNameI)
        {
            const word& patchName = backPatchNames[patchNameI];

            label patchI = bMesh.findPatchID(patchName);

            if (patchI == -1)
            {
                FatalErrorIn(args.executable()) << "No such patch "
                    << patchName << endl << "Patches are " << bMesh.names()
                    << exit(FatalError);
            }
            includeBackPatches.insert(patchI);
        }
    }

    // ------------------------------
    // THERMAL
    // ------------------------------

    // Reading thermal mesh and creating triSurface
    Info<< "Reading thermal mesh" << endl;
    triSurface surface(TAIthermFileName);

    const Field<point>& pointsThermal = surface.points();
    //const vectorField&  normalsThermal = surface.pointNormals();
    const vectorField&  normalsThermal = surface.Sf();
    const vectorField&  centresThermal = surface.Cf();

    // Creating thermal octree
    Info<< "Constructing Octree for thermal" << endl;
    //construct octree
    //treeBoundBox overallBb(pointsThermal);
    treeBoundBox overallBb(centresThermal);
    Random rndGen(123456);
    overallBb = overallBb.extend(rndGen, 1E-4);
    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    autoPtr<indexedOctree<treeDataPoint>> thermalPointTreePtr
    (
            new indexedOctree<treeDataPoint>
            (
                treeDataPoint(centresThermal),
                //treeDataPoint(pointsThermal),
                overallBb,  // overall search domain
                8,          // maxLevel
                10,         // leafsize
                3.0         // duplicity
            )
    );

    //loop over points, find nearest and assign mapped data
    scalar thermalSpanSqr = Foam::sqr(thermalPointTreePtr->bb().mag());

    // ------------------------------
    // CFD
    // ------------------------------

    // ------------------------------
    // Creating CFD lists
    label nBoundaryFaces = mesh.nFaces() - mesh.nInternalFaces();
    DynamicList<scalar> patchHTC(nBoundaryFaces);
    DynamicList<scalar> patchTemp(nBoundaryFaces);
    DynamicList<point>  patchFaceCentres(nBoundaryFaces);

    DynamicList<scalar> patchBackHTC(nBoundaryFaces);
    DynamicList<scalar> patchBackTemp(nBoundaryFaces);
    DynamicList<point>  patchBackFaceCentres(nBoundaryFaces);

    // Create a list for front only patches
    forAllConstIter(labelHashSet, includePatches, iter)
    {
        label patchI = iter.key();
        const word patchName = bMesh[patchI].name();

        // THis should be done at a later stage to be more automatic
        ////if(isA<mappedFromFileFvPatchField<scalar>>(Tb))
        //if(isA<inletOutletFvPatchField<scalar>>(Tb))
        //{

    forAll(bMesh[patchI], j)
    {

      const label meshFaceI = bMesh[patchI].start() + j;
      const point fc = mesh.faceCentres()[meshFaceI];
      // const vector fn = mesh.faceAreas()[meshFaceI];
      const scalar nwt =
        T.boundaryField()[patchI].patchInternalField()()[j];
      const scalar htc = alphaConv.boundaryField()[patchI][j];

      pointIndexHit info = thermalPointTreePtr->findNearest
      (
       fc,
       thermalSpanSqr
      );

      if (!info.hit())
      {
        info = thermalPointTreePtr->findNearest
        (
         fc,
         GREAT*thermalSpanSqr
         );
      }

      patchFaceCentres.append(fc);
      patchTemp.append(nwt);
      patchHTC.append(htc);
    }
        //}
    }

    // Split front-and-back CFD patches into front and back and create a list
    forAllConstIter(labelHashSet, includeBackPatches, iter)
    {
        label patchI = iter.key();
        const word patchName = bMesh[patchI].name();

        forAll(bMesh[patchI], j)
        {

            const label meshFaceI = bMesh[patchI].start() + j;
            const point fc = mesh.faceCentres()[meshFaceI];
            const vector fn = mesh.faceAreas()[meshFaceI];
            const scalar nwt =
                T.boundaryField()[patchI].patchInternalField()()[j];
            const scalar htc = alphaConv.boundaryField()[patchI][j];

            pointIndexHit info = thermalPointTreePtr->findNearest
            (
                fc,
                thermalSpanSqr
            );

            if (!info.hit())
            {
                info = thermalPointTreePtr->findNearest
                (
                    fc,
                    GREAT*thermalSpanSqr
                );
            }

            if ((fn & normalsThermal[info.index()]) > 0)
            {
                patchBackFaceCentres.append(fc);
                patchBackTemp.append(nwt);
                patchBackHTC.append(htc);
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

    patchBackFaceCentres.shrink();
    patchBackTemp.shrink();
    patchBackHTC.shrink();

    thermalPointTreePtr.clear();

    //Gather everything on master

    // Gather all faceCentres on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()].transfer(patchFaceCentres);
    Pstream::gatherList(gatheredPoints);

    //Gather all htc on master
    List<scalarField> gatheredHTC(Pstream::nProcs());
    gatheredHTC[Pstream::myProcNo()].transfer(patchHTC);
    Pstream::gatherList(gatheredHTC);

    //Gather all temp on master
    List<scalarField> gatheredTemp(Pstream::nProcs());
    gatheredTemp[Pstream::myProcNo()].transfer(patchTemp);
    Pstream::gatherList(gatheredTemp);

    patchFaceCentres.clearStorage();
    patchTemp.clearStorage();
    patchHTC.clearStorage();

    // Back: Gather all faceCentres on master
    List<pointField> gatheredBackPoints(Pstream::nProcs());
    gatheredBackPoints[Pstream::myProcNo()].transfer(patchBackFaceCentres);
    Pstream::gatherList(gatheredBackPoints);

    // Back: Gather all htc on master
    List<scalarField> gatheredBackHTC(Pstream::nProcs());
    gatheredBackHTC[Pstream::myProcNo()].transfer(patchBackHTC);
    Pstream::gatherList(gatheredBackHTC);

    // Back: Gather all temp on master
    List<scalarField> gatheredBackTemp(Pstream::nProcs());
    gatheredBackTemp[Pstream::myProcNo()].transfer(patchBackTemp);
    Pstream::gatherList(gatheredBackTemp);

    patchBackFaceCentres.clearStorage();
    patchBackTemp.clearStorage();
    patchBackHTC.clearStorage();

    if (Pstream::master())
    {
        // On master combine all lists
        pointField pointsCFD =
            ListListOps::combine<pointField>
            (
                gatheredPoints,
                accessOp<pointField>()
            );

        scalarField HTC =
            ListListOps::combine<scalarField>
            (
                gatheredHTC,
                accessOp<scalarField>()
            );

        scalarField Temp =
            ListListOps::combine<scalarField>
            (
                gatheredTemp,
                accessOp<scalarField>()
            );

        pointField pointsCFDback =
            ListListOps::combine<pointField>
            (
                gatheredBackPoints,
                accessOp<pointField>()
            );

        scalarField HTCback =
            ListListOps::combine<scalarField>
            (
                gatheredBackHTC,
                accessOp<scalarField>()
            );

        scalarField TempBack =
            ListListOps::combine<scalarField>
            (
                gatheredBackTemp,
                accessOp<scalarField>()
            );

        // --- Front CFD Octree
        Info<< "Constructing Octree for CFD" << endl;

        //construct octree for Front
        treeBoundBox overallBbFront(pointsCFD);
        Random rndGenF(123456);
        overallBbFront = overallBbFront.extend(rndGenF, 1E-4);
        overallBbFront.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBbFront.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        indexedOctree<treeDataPoint> pointTreeFront
        (
            treeDataPoint(pointsCFD),
            overallBbFront,  // overall search domain
            8,          // maxLevel
            10,         // leafsize
            3.0         // duplicity
        );

        //loop over points, find nearest and assign mapped data
        scalar spanSqrFront = Foam::sqr(pointTreeFront.bb().mag());


        autoPtr<indexedOctree<treeDataPoint>> pointTreeBackPtr(nullptr);

        //loop over points, find nearest and assign mapped data
        scalar spanSqrBack = 0.;

        if (pointsCFDback.size() == 0)
        {
            Info<< "No front-and-back patches found" << endl;
        }
        else
        {
            // --- Back CFD Octree
            Info<< "Constructing Octree for CFD - Back" << endl;

            //construct octree for Front
            treeBoundBox overallBbBack(pointsCFDback);
            Random rndGenB(123456);
            overallBbBack = overallBbBack.extend(rndGenB, 1E-4);
            overallBbBack.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
            overallBbBack.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

            indexedOctree<treeDataPoint> pointTreeBackTmp
            (
                treeDataPoint(pointsCFDback),
                overallBbBack,  // overall search domain
                8,          // maxLevel
                10,         // leafsize
                3.0         // duplicity
            );

            pointTreeBackPtr = pointTreeBackTmp.clone();
            spanSqrBack = Foam::sqr(pointTreeBackPtr->bb().mag());
        }

        // --- Writing Neutral File

        Info<< "Writing Neutral File" << endl;

        IFstream is(TAIthermFileName);
        fileName outname (TAIthermFileName.substr(0,TAIthermFileName.size()-4)+"_foam.ntl");
        autoPtr<OFstream> TAIthermFilePtr(new OFstream(outname));

        string dateString = clock::date();
        string timeString = clock::clockTime();

        labelList pointIDs(pointsThermal.size());
        Map<label> indexToPoint(0);//(2*pointsThermal.size());
        labelList elementsIDs(0);
        List<vector> elementsCtrs(0);

        label pntI(0);
        label elemI(0);

        if (!is.good())
        {
            FatalErrorInFunction
                << "Cannot read file " << TAIthermFileName
                << exit(FatalError);
        }

        Info<< "Writing thermal mesh" << endl;

        while (is.good())
        {
            size_t linei = 2;
            string line;
            is.getLine(line);

            // Build inverse mapping (index to point)
            if (pntI == pointsThermal.size()
            && indexToPoint.size() == 0
            )
            {
                Info<< "Creating inverse indexing" << endl;
                indexToPoint.resize(2*pointsThermal.size());

                forAll(pointIDs, i)
                {
                    indexToPoint.insert(pointIDs[i], i);
                }
            }

            // Reading/Writing Title Card
            if (line.substr(0, 2) == "25")
            {
                TAIthermFilePtr() << line.c_str() << endl;
                TAIthermFilePtr() << "PATRAN Neutral File created by FOAM core" << endl;
                is.getLine(line);
            }
            // Reading/Writing Summary Data
            else if (line.substr(0, 2) == "26")
            {
                // Get elements count
                label index = readLabel(IStringStream(readNTLtoken(line,8,linei))());
                index = readLabel(IStringStream(readNTLtoken(line,8,linei))());
                index = readLabel(IStringStream(readNTLtoken(line,8,linei))());
                index = readLabel(IStringStream(readNTLtoken(line,8,linei))());
                index = readLabel(IStringStream(readNTLtoken(line,8,linei))());

                elementsIDs.resize(index,0);
                elementsCtrs.resize(index,vector::zero);

                TAIthermFilePtr() << line.c_str() << endl;
                is.getLine(line);
                TAIthermFilePtr() << dateString.c_str()
                                  << "    "
                                  << timeString.c_str() << endl;
            }
            // Reading/Writing Node Data
            else if (line.substr(0, 2) == " 1")
            {
                // Get ID
                label index =
                    readLabel(IStringStream(readNTLtoken(line,8,linei))());
                pointIDs[pntI]=index;
                pntI++;

                TAIthermFilePtr() << line.c_str() << endl;
                is.getLine(line);
                TAIthermFilePtr() << line.c_str() << endl;
                is.getLine(line);
                TAIthermFilePtr() << line.c_str() << endl;
            }
            // Reading/Writing Element Data
            else if (line.substr(0, 2) == " 2")
            {
                // Get ID
                label index =
                    readLabel(IStringStream(readNTLtoken(line,8,linei))());
                // Get shape
                label shape =
                    readLabel(IStringStream(readNTLtoken(line,8,linei))());

                vector elemCtr(vector::zero);

                TAIthermFilePtr() << line.c_str() << endl;
                is.getLine(line);
                TAIthermFilePtr() << line.c_str() << endl;
                is.getLine(line);
                TAIthermFilePtr() << line.c_str() << endl;

                if (shape==3)
                {
                    linei = 0;

                    label a = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                    label b = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                    label c = readLabel(IStringStream(readNTLtoken(line, 8, linei))());

                    elemCtr = ( pointsThermal[indexToPoint[a]]
                              + pointsThermal[indexToPoint[b]]
                              + pointsThermal[indexToPoint[c]]
                              )/3.;
                }
                else if (shape==4)
                {
                    linei = 0;

                    label a = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                    label b = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                    label c = readLabel(IStringStream(readNTLtoken(line, 8, linei))());
                    label d = readLabel(IStringStream(readNTLtoken(line, 8, linei))());

                    elemCtr = ( pointsThermal[indexToPoint[a]]
                              + pointsThermal[indexToPoint[b]]
                              + pointsThermal[indexToPoint[c]]
                              + pointsThermal[indexToPoint[d]]
                              )/4.;
                }
                else
                {
                    Info<< "Bar shape not implemented" << endl;
                }

                elementsIDs[elemI] = index;
                elementsCtrs[elemI] = elemCtr;
                elemI++;
            }
            // Reading/Writing Element Properties
            else if (line.substr(0, 2) == " 4")
            {
                TAIthermFilePtr() << line.c_str() << endl;
                is.getLine(line);
                TAIthermFilePtr() << line.c_str() << endl;
            }
            // Reading/Writing Unknown
            else if (line.substr(0, 2) == "98")
            {
                TAIthermFilePtr() << line.c_str() << endl;
                is.getLine(line);
                TAIthermFilePtr() << line.c_str() << endl;
            }
        }

        Info<< "Writing thermal values" << endl;

        // Write HTC and Temps
        scalar defaultValue=-999;
        forAll(elementsIDs, eI)
        {
            scalar frontHTC(0);
            scalar backHTC(0);
            scalar frontT(273.15);
            scalar backT(273.15);

            bool hasBack(false);

            pointIndexHit info = pointTreeFront.findNearest
            (
                elementsCtrs[eI],
                spanSqrFront
            );

            if (info.hit())
            {
                frontHTC = HTC[info.index()];
                frontT = Temp[info.index()];
            }


            if (pointsCFDback.size() > 0)
            {
                pointIndexHit infoB = pointTreeBackPtr->findNearest
                (
                    elementsCtrs[eI],
                    spanSqrBack
                );

                if (infoB.hit())
                {
                    backHTC = HTCback[infoB.index()];
                    backT = TempBack[infoB.index()];
                    hasBack = true;
                }
            }


            //Write HTC
            TAIthermFilePtr() << "17"
                              << setw(8) << elementsIDs[eI]
                              << "       0"
                              << "       2"
                              << "       1"
                              << "       1"
                              << "       0"
                              << "       0"
                              << "       0"
                              << endl;

            if (hasBack)
            {
                TAIthermFilePtr() << "1 "
                                  << "11000000"
                                  << endl;
            }
            else
            {
                TAIthermFilePtr() << "1 "
                                  << "10000000"
                                  << endl;
            }

            Foam::scientific(TAIthermFilePtr());

            TAIthermFilePtr() << setw(16) << setprecision(9)
                                  <<  frontHTC
                      << setw(16) << backHTC
                      << endl;

            //Write near wall temperature
            TAIthermFilePtr() << "18"
                              << setw(8) << elementsIDs[eI]
                              << "       0"
                              << "       2"
                              << "       1"
                              << "       1"
                              << "       0"
                              << "       0"
                              << "       0"
                              << endl;

             TAIthermFilePtr() << "1 "
                               << "00000000"
                               << endl;

            if (hasBack)
            {
                TAIthermFilePtr() << setw(16) << setprecision(9)
                                      << frontT
                          << setw(16) << backT
                          << setw(16) << defaultValue
                          << setw(16) << defaultValue
                          << endl;
            }
            else
            {
                TAIthermFilePtr() << setw(16) << setprecision(9)
                                      << frontT
                          << setw(16) << defaultValue
                          << setw(16) << defaultValue
                          << setw(16) << defaultValue
                          << endl;
            }
        }

        //write packet 99 (End of File)
        TAIthermFilePtr()<<"99"<<"       0"<<"       0"<<"       1"<<"       0"
                         <<"       0"<<"       0"<<"       0"<<"       0"<<endl;
    }

    totTimeClo -= runTime.elapsedClockTime();
    totTimeCPU -= runTime.elapsedCpuTime();

    Info<< "Mapped and wrote file in " << endl;
    Info<< -totTimeClo << " s Clock Time " << endl;
    Info<< -totTimeCPU << " s CPU Time " << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
