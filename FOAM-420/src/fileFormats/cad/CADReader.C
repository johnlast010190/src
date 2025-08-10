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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "CADReader.H"
#include "containers/HashTables/HashTable/HashTable.H"
#include "containers/HashTables/HashSet/HashSet.H"
#include "primitives/strings/stringOps/stringOps.H"

#include "include/TBBTimer.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileFormats
{
    defineTypeNameAndDebug(CADReader, 0);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::CADReader::CADReader
(
    const fileName& CADfileName,
    const word& nameRegionsBy,
    const bool distributed, /*= false*/
    const wordHashSet& ignoreRegions /*= wordHashSet()*/
)
:
    inFile_(CADfileName),
    nameRegionsBy_(nameRegionsBy),
    step_(false),
    iges_(false),
    regionNames_(),
    regionAreas_()
{
    //Timer t("Load CAD file");

    Info<< "Loading CAD file " << inFile_ << nl << endl;

    meshedShapes_.Clear();
    oneShape_.Nullify();

    Handle(TColStd_HSequenceOfTransient) entityList;
    Handle(Transfer_TransientProcess) tp;

    STEPControl_Reader stepReader;
    IGESControl_Reader igesReader;

    IFSelect_ReturnStatus readStat = IFSelect_RetDone;

    word ext = CADfileName.ext();

    if (ext == "step" || ext == "stp")
    {
        step_ = true;

        readStat = stepReader.ReadFile(inFile_.c_str());

        if (readStat != IFSelect_RetDone)
        {
            FatalErrorInFunction
                << "Error while opening "
                << inFile_
                << abort(FatalError);
        }

        stepReader.TransferRoots();

        if (Pstream::master())
        {
            stepReader.PrintCheckLoad(/*failsonly*/!debug, IFSelect_ItemsByEntity);
        }

        oneShape_ = stepReader.OneShape();

        const Handle(XSControl_TransferReader)& TR = stepReader.WS()->TransferReader();

        entityList = stepReader.GiveList("xst-model-all");

        tp = TR->TransientProcess();
    }
    else if (ext == "iges" || ext == "igs")
    {
        iges_ = true;

        readStat = igesReader.ReadFile(inFile_.c_str());

        if (readStat != IFSelect_RetDone)
        {
            FatalErrorInFunction
                << "Error while opening "
                << inFile_
                << abort(FatalError);
        }

        igesReader.TransferRoots();

        if (Pstream::master())
        {
            igesReader.PrintCheckLoad(/*failsonly*/!debug, IFSelect_ItemsByEntity);
        }

        oneShape_ = igesReader.OneShape();

        const Handle(XSControl_TransferReader)& TR = igesReader.WS()->TransferReader();

        entityList = igesReader.GiveList("xst-model-all");

        tp = TR->TransientProcess();

    }
    else
    {
        FatalErrorInFunction
            << "Error while loading " << inFile_
            << "\nPlease use files of types \"STEP\" or \"IGES\"."
            << abort(FatalError);
    }

    // Apply shape healing.
    // TODO propagate region names to healed shapes
    //if ( doHealing )
    //{
    //    Handle(ShapeFix_Shape) shapeHealer = new ShapeFix_Shape(oneShape_);
    //    //shapeHealer->SetPrecision ( 1e-5 );
    //    //shapeHealer->SetMaxTolerance ( maxTol );
    //    //shapeHealer->SetMinTolerance ( mintol );
    //    bool fixRes = false;
    //    try
    //    {
    //        fixRes = shapeHealer->Perform();
    //    }
    //    catch ( Standard_Failure& )
    //    {
    //        fixRes = false;
    //    }
    //    if (fixRes)
    //    {
    //        oneShape_ = shapeHealer->Shape();
    //    }
    //}

    // find name of shapes
    TopAbs_ShapeEnum sType = TopAbs_SHAPE;

    if (nameRegionsBy_ == "solid")
    {
        sType = TopAbs_SOLID;
    }
    else if (nameRegionsBy_ == "shell")
    {
        sType = TopAbs_SHELL;
    }

    TopTools_IndexedMapOfShape shapeMap;

    TopExp::MapShapes(oneShape_, sType, shapeMap);

    Info<< "Found " << shapeMap.Extent() << " "
        << nameRegionsBy_ << " shapes" << endl;

    Standard_Integer nb =  entityList->Length();

    for (Standard_Integer i = 1; i <= shapeMap.Extent(); i++)
    {
        const TopoDS_Shape& aShape = shapeMap(i);

        TCollection_AsciiString name = "";

        for (Standard_Integer i = 1; i <= nb; i++)
        {
            if (step_)
            {
                Handle(StepRepr_RepresentationItem) ent =
                    Handle(StepRepr_RepresentationItem)::DownCast(entityList->Value(i));

                if (ent.IsNull()) continue;

                TopoDS_Shape sh = TransferBRep::ShapeResult(tp, ent);

                if (sh.IsNull()) continue;

                if (!sh.IsPartner(aShape)) continue;

                if (ent->Name().IsNull()) continue;
                if (ent->Name()->IsEmpty()) continue;
                if (ent->Name()->String().IsEqual("NONE")) continue;

                name = ent->Name()->String();
            }
            else if (iges_)
            {
                Handle(IGESData_IGESEntity) ent =
                    Handle(IGESData_IGESEntity)::DownCast(entityList->Value(i));

                if (ent.IsNull()) continue;

                TopoDS_Shape sh = TransferBRep::ShapeResult(tp, ent);

                if (sh.IsNull()) continue;

                if (!sh.IsPartner(aShape)) continue;

                if (ent->HasShortLabel())
                {
                    // A short label is a non-blank 8-character string.
                    name = ent->ShortLabel()->String();
                }
                else
                {
                    name = ent->NameValue()->String();
                }
            }
        }


        if (ignoreRegions.found(name.ToCString()))
        {
            Info<< "Ignoring " << name.ToCString() <<endl;
            continue;
        }

        regionNames_.append(name.ToCString());

        if (debug)
            Info<< "name: " << name.ToCString() << endl;

        if (distributed)
        {
            GProp_GProps shapeProperties;
            BRepGProp::SurfaceProperties
            (
                aShape,
                shapeProperties,
                /*SkipShared*/true,
                /*UseTriangulation*/true
            );
            scalar area = shapeProperties.Mass();

            if (debug)
                Info<< "     area: " << area <<endl;

            regionAreas_.append(area);
        }

        meshedShapes_.Append(aShape);
    }

    meshedEntitiesNTri_.setSize(shapeMap.Extent());

    // Release memory after translation.
    stepReader.WS()->NewModel();
    igesReader.WS()->NewModel();
    Standard::Purge();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fileFormats::CADReader::distributeShapes
(
    const bool balanceArea /*= false*/
)
{
    const label nProcs = Pstream::nProcs();

    List<labelHashSet> shapesPerProc(nProcs);

    labelList cumTriPerProc(nProcs);

    scalar triPerProc = nTotalTriangles_ / nProcs;

    scalarField cumAreaPerProc(nProcs);

    scalar areaPerProc = sum(regionAreas_) / nProcs;

    TopTools_ListIteratorOfListOfShape it(meshedShapes_);

    label i = 0;

    for (; it.More(); it.Next())
    {
        label procI = i % nProcs;

        if (i < nProcs)
        {
            // Initialize with one shape per processor
            shapesPerProc[procI].insert(i);
            cumTriPerProc[procI] = meshedEntitiesNTri_[i];
            cumAreaPerProc[procI] = regionAreas_[i];
        }
        else
        {
            scalar deltaTri = abs(cumTriPerProc[procI] - triPerProc) / triPerProc;
            scalar deltaArea = abs(cumAreaPerProc[procI] - areaPerProc) / areaPerProc;

            if (balanceArea && deltaArea > 0.05)
            {
                shapesPerProc[procI].insert(i);
                cumAreaPerProc[procI] += regionAreas_[i];
            }
            else if (deltaTri > 0.05)
            {
                shapesPerProc[procI].insert(i);
                cumTriPerProc[procI] += meshedEntitiesNTri_[i];
            }
        }
        i++;
    }

    if (debug)
    {
        Pstream::gatherList(cumTriPerProc);
        Info<<"\ncoarse triangles per processor:" << nl
            <<"wanted: " << triPerProc << nl
            <<"current: " << cumTriPerProc
            <<endl;

        Pstream::gatherList(cumAreaPerProc);
        Info<<"\nshapes area per processor:" << nl
            <<"wanted: " << areaPerProc << nl
            <<"current: " << cumAreaPerProc
            <<endl;
    }

    oneShape_.Nullify();

    TopTools_ListOfShape allShapes(meshedShapes_);
    meshedShapes_.Clear();

    List<word> allRegionNames(regionNames_);
    regionNames_.clear();

    label nTransferEntity = 0;

    TopoDS_Compound aCompound;

    BRep_Builder aBuilder;

    aBuilder.MakeCompound(aCompound);

    it = TopTools_ListIteratorOfListOfShape(allShapes);

    i = 0;

    for (; it.More(); it.Next())
    {
        TopoDS_Shape aShape = it.Value();

        if (shapesPerProc[Pstream::myProcNo()][i])
        {
            aBuilder.Add(aCompound, aShape);

            meshedShapes_.Append(aShape);

            regionNames_.append(allRegionNames[i]);

            nTransferEntity++;
        }
        i++;
    }

    labelList nShapesPerProc(Pstream::nProcs(), nTransferEntity);
    Pstream::gatherList(nShapesPerProc);
    Info<<"shapesPerProc " << nShapesPerProc << endl;

    oneShape_ = aCompound;
}

void Foam::fileFormats::CADReader::distributeShapes
(
    const List<List<treeBoundBox>>& procBb
)
{
    oneShape_.Nullify();

    const treeBoundBox& bb = procBb[Pstream::myProcNo()][0];

    label nTransferEntity = 0;

    TopoDS_Compound aCompound;

    BRep_Builder aBuilder;

    aBuilder.MakeCompound(aCompound);

    TopTools_ListOfShape allShapes(meshedShapes_);
    meshedShapes_.Clear();

    List<word> allRegionNames(regionNames_);
    regionNames_.clear();

    TopTools_ListIteratorOfListOfShape it(allShapes);

    label i = 0;

    for (; it.More(); it.Next())
    {
        TopoDS_Shape aShape = it.Value();

        Bnd_Box aBndBox;

        BRepBndLib::Add(aShape, aBndBox);

        if (debug)
        {
            aBndBox.DumpJson(std::cout);
            Info<< nl;
        }

        gp_Pnt cornerMin = aBndBox.CornerMin();
        gp_Pnt cornerMax = aBndBox.CornerMax();

        boundBox bbShape
        (
            point
            (
                cornerMin.X(),
                cornerMin.Y(),
                cornerMin.Z()
            ),
            point
            (
                cornerMax.X(),
                cornerMax.Y(),
                cornerMax.Z()
            )
        );

        if (bb.contains(bbShape.midpoint()))
        {
            aBuilder.Add(aCompound, aShape);

            meshedShapes_.Append(aShape);

            regionNames_.append(allRegionNames[i]);

            nTransferEntity++;
        }
        i++;
    }

    labelList nShapesPerProc(Pstream::nProcs(), nTransferEntity);
    Pstream::gatherList(nShapesPerProc);
    Info<<"shapesPerProc " << nShapesPerProc << endl;

    oneShape_ = aCompound;
}

void Foam::fileFormats::CADReader::triangulate
(
    scalar linDeflection,
    scalar angDeflection,
    pointField& surfPoints,
    List<labelledTri>& surfTriangles,
    HashTable<label>& patchIDs
)
{
    //Timer t("CADReader::triangulate");

    Info<< "\nTriangulating... (linDeflection "
        << linDeflection
        << "; angDeflection "
        << angDeflection
        << ";)"
        << nl << endl;

    //test on F1_car.stp (Lenovo laptop AMD Ryzen 5 PRO 3500U w/ Radeon Vega Mobile Gfx 2.10 GHz 8 cores)
    ////BRepMesh_IncrementalMesh Mesh
    //    === Timer:                 CADReader::triangulate, elapsed      47.3 s (cpu   297 s)
    //    Triangles    : 5714615 in 40 region(s)
    //    Vertices     : 2991379
    ////IMeshTools_Parameters default;
    //    === Timer:                 CADReader::triangulate, elapsed      76.9 s (cpu   358 s)
    //    Triangles    : 5792321 in 40 region(s)
    //    Vertices     : 3030238
    ////IMeshTools_Parameters Delabella
    //    === Timer:                 CADReader::triangulate, elapsed      95.5 s (cpu   300 s)
    //    Triangles    : 5426105 in 40 region(s)
    //    Vertices     : 2847130

    //IMeshTools_Parameters aMeshParams;
    //aMeshParams.Deflection               = linDeflection;
    //aMeshParams.Angle                    = angDeflection;
    //aMeshParams.Relative                 = Standard_False;
    //aMeshParams.InParallel               = Standard_True;
    //aMeshParams.MinSize                  = Precision::Confusion();
    //aMeshParams.InternalVerticesMode     = Standard_True;
    //aMeshParams.ControlSurfaceDeflection = Standard_True;

    //Handle(IMeshTools_Context) aContext = new BRepMesh_Context();
    ////aContext->SetFaceDiscret (new BRepMesh_FaceDiscret (new BRepMesh_DelabellaMeshAlgoFactory()));

    //BRepMesh_IncrementalMesh aMesher;
    //aMesher.SetShape (oneShape_);
    //aMesher.ChangeParameters() = aMeshParams;
    //aMesher.Perform (aContext);

    BRepMesh_IncrementalMesh Mesh
    (
        oneShape_,
        linDeflection,
        Standard_False, //isRelative
        angDeflection,
        Standard_True   //isInParallel
    );

    if (debug)
    {
        fileName outFile =
            inFile_.lessExt()
          + "_OCCT_" + name(linDeflection)
          + "_" + name(angDeflection) + ".stl";

        Info<< "Writing OCCT triangulation to "
            << outFile
            << endl;

        StlAPI_Writer stlWriter = StlAPI_Writer();
        stlWriter.ASCIIMode() = Standard_True;
        stlWriter.Write(oneShape_, outFile.c_str());
    }

    label aNbNodes = 0;
    nTotalTriangles_ = 0;

    TopTools_ListIteratorOfListOfShape it(meshedShapes_);

    label i = 1;

    //count nodes and triangles
    for (; it.More(); it.Next())
    {
        TopoDS_Shape aShape = it.Value();

        TopTools_IndexedMapOfShape shapeMap;

        TopExp::MapShapes(aShape, TopAbs_FACE, shapeMap);

        if (debug)
        {
            Info<<"Shape " << regionNames_[i-1]
                <<" has " << shapeMap.Extent() << " FACEs" << endl;
        }

        for (Standard_Integer j = 1; j <= shapeMap.Extent(); j++)
        {
            const TopoDS_Shape& aFace = shapeMap(j);

            TopLoc_Location aLoc;

            Handle(Poly_Triangulation) aTriangulation;
            try
            {
                OCC_CATCH_SIGNALS
                aTriangulation = BRep_Tool::Triangulation(TopoDS::Face(aFace), aLoc);
            }
            catch (Standard_TypeMismatch const&)
            {
                Info<< "Can't convert to face..." <<endl;
                continue;
            }

            if (aTriangulation.IsNull())
            {
                Bnd_Box boundingBox;
                boundingBox.SetGap(0.0);
                BRepBndLib::Add(aFace, boundingBox);
                boundingBox.SetGap(0.0);

                boundBox bb;
                const point& min = bb.min();
                const point& max = bb.max();

                Standard_Real minx = min.x();
                Standard_Real miny = min.y();
                Standard_Real minz = min.z();
                Standard_Real maxx = max.x();
                Standard_Real maxy = max.y();
                Standard_Real maxz = max.z();

                boundingBox.Get(minx, miny, minz, maxx, maxy, maxz);

                //TODO write boundBox as obj for visualization!
                Info<< "WARNING: Null triangulation for face at " << bb.midpoint()
                    << " (boundBox span "<< bb.span() << ") "
                    << "Please check your CAD!"
                    << endl;
            }
            else
            {
                aNbNodes += aTriangulation->NbNodes();
                nTotalTriangles_ += aTriangulation->NbTriangles();

                meshedEntitiesNTri_[i-1] += aTriangulation->NbTriangles();
            }
        }
        i++;
    } //count nodes and triangles

    DebugInformation << "Total points : " << aNbNodes << endl;
    DebugInformation << "Total triangles : " << nTotalTriangles_ << endl;

    surfPoints.setSize(aNbNodes);
    surfTriangles.setSize(nTotalTriangles_);

    label nPatches = 0;

    // count faces missing triangulation
    Standard_Integer aNbFacesNoTri = 0;
    // fill temporary triangulation
    Standard_Integer aNodeOffset = 0;
    Standard_Integer aTriangleOffet = 0;

    TopTools_ListIteratorOfListOfShape it0(meshedShapes_);

    i = 0;

    for (; it0.More(); it0.Next())
    {
        TopoDS_Shape aShape = it0.Value();

        word patchName = regionNames_[i++];

        if (!patchIDs.found(patchName))
        {
            Info<< "Found region: " << patchName << endl;
            patchIDs.insert(patchName, nPatches);

            nPatches++;
        }

        TopTools_IndexedMapOfShape shapeMap;
        TopExp::MapShapes(aShape, TopAbs_FACE, shapeMap);

        for (Standard_Integer j = 1; j <= shapeMap.Extent(); j++)
        {
            const TopoDS_Shape& aFace = shapeMap(j);

            TopLoc_Location aLoc;

            Handle(Poly_Triangulation) aTriangulation;

            try
            {
                OCC_CATCH_SIGNALS
                aTriangulation = BRep_Tool::Triangulation(TopoDS::Face(aFace), aLoc);
            }
            catch (Standard_TypeMismatch const&)
            {
                Info<< "Can't convert to face...\n";
                ++aNbFacesNoTri;
                continue;
            }

            if (aTriangulation.IsNull())
            {
                ++aNbFacesNoTri;
                continue;
            }

            // copy nodes //TODO multithreaded
            gp_Trsf aTrsf = aLoc.Transformation();
            for (Standard_Integer aNodeIter = 1; aNodeIter <= aTriangulation->NbNodes(); aNodeIter++)
            {
                gp_Pnt aPnt = aTriangulation->Node(aNodeIter);
                aPnt.Transform(aTrsf);

                point xyz
                (
                    aPnt.X(),
                    aPnt.Y(),
                    aPnt.Z()
                );

                label pointIndex = aNodeIter + aNodeOffset - 1;
                surfPoints[pointIndex] = xyz;
            }

            // copy triangles //TODO multithreaded
            const TopAbs_Orientation anOrientation = aFace.Orientation();
            for (Standard_Integer aTriIter = 1; aTriIter <= aTriangulation->NbTriangles(); aTriIter++)
            {
                const Poly_Triangle aTri = aTriangulation->Triangle(aTriIter);

                triFace anId;

                aTri.Get(anId[0], anId[1], anId[2]);

                if (anOrientation == TopAbs_REVERSED)
                {
                    // Swap 1, 2.
                    Standard_Integer aTmpIdx = anId[1];
                    anId[1] = anId[2];
                    anId[2] = aTmpIdx;
                }

                // Update nodes according to the offset.
                anId[0] += aNodeOffset - 1;
                anId[1] += aNodeOffset - 1;
                anId[2] += aNodeOffset - 1;

                label triIndex = aTriIter + aTriangleOffet - 1;

                labelledTri aLabelledTri(anId, patchIDs[patchName]);

                surfTriangles[triIndex] = aLabelledTri;
            }

            aNodeOffset += aTriangulation->NbNodes();
            aTriangleOffet += aTriangulation->NbTriangles();
        }
    }
}
