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

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileFormats
{
    defineTypeNameAndDebug(CADReader, 0);
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::CADReader::TransferShapes()
{
    Handle(TColStd_HSequenceOfTransient) entityList =
        reader_.GiveList("xst-model-all");

    Standard_Integer nb =  entityList->Length();

    const Handle(XSControl_TransferReader)& TR = reader_.WS()->TransferReader();

    /////////////////////////
    reader_.TransferRoots();
    oneShape_ = reader_.OneShape();
    //// divide into sub shapes
    TopTools_IndexedMapOfShape shapeMap;
    TopExp::MapShapes(oneShape_, TopAbs_SOLID, shapeMap);
    Info<<"GGG Found " << shapeMap.Extent() << " SOLIDs" << endl;

    for (Standard_Integer i = 1; i <= shapeMap.Extent(); i++)
    {
        const TopoDS_Shape& aShape = shapeMap(i);

        TCollection_AsciiString name;

        //if(step)      name = this->GetNameOfShape(aShape, stepReader, shapeList);
        //else if (iges) name = this->GetNameOfShape(aShape, igesReader, shapeList);

        //if(name.IsEmpty())
        //    name += i;

        //shapesToExport.Append(aShape);
        //this->AddToNameList(name,0);
        TCollection_AsciiString result = "";
        const Handle(Transfer_TransientProcess)& tp = TR->TransientProcess();

        for (Standard_Integer i = 1; i <= nb; i++)
        {
            Handle(StepRepr_RepresentationItem) ent =
                Handle(StepRepr_RepresentationItem)::DownCast(entityList->Value(i));

            if (ent.IsNull()) continue;

            TopoDS_Shape sh = TransferBRep::ShapeResult(tp, ent);

            if (sh.IsNull()) continue;

            if (!sh.IsPartner(aShape)) continue;

            if (ent->Name().IsNull()) continue;
            result = ent->Name()->String();

            //////////
            reader_.TransferEntity(ent);

            reader_.PrintCheckTransfer
            (
                !debug, //failsonly,
                IFSelect_ItemsByEntity
            );
            meshedEntities_->Append(ent);
            //////////
        }
        Info<<"GGG result: " << result.ToCString()<<endl;
    }
    ///////////////////////

    label nTransferEntity = 0;

    //TR->BeginTransfer();


    //for (Standard_Integer i = 1; i <= nb; i++)
    //{
    //    Handle(Standard_Transient) entity = entityList->Value(i);

    //    if (entity.IsNull()) continue;

    //    if
    //    (
    //        entity->DynamicType()->SubType("IGESGeom_TrimmedSurface")
    //     || entity->DynamicType()->SubType("IGESSolid_Face")
    //     || entity->DynamicType()->SubType("StepShape_Face")
    //    )
    //    {
    //        //Info<<"GGG 0"<<endl;
    //        if (TR->TransferOne(entity) == 0) continue;
    //        //try
    //        //{
    //        //    OCC_CATCH_SIGNALS
    //        //    if (TR->TransferOne(entity) == 0) continue;
    //        //}
    //        ////catch (std::exception& ex)
    //        ////catch (Standard_TypeMismatch const&)
    //        //catch (Standard_Failure const& theFailure)
    //        //{
    //        //    Info<< "GGG catch " << theFailure.GetMessageString() <<endl;
    //        //    //throw error(ex.what());
    //        //    continue;
    //        //}
    //        //Info<<"GGG 1"<<endl;

    //        reader_.TransferEntity(entity);

    //        reader_.PrintCheckTransfer
    //        (
    //            !debug, //failsonly,
    //            IFSelect_ItemsByEntity
    //        );

    //        Handle(StepRepr_RepresentationItem) stepEntity =
    //            Handle(StepRepr_RepresentationItem)::DownCast(entity);
    //        Handle(TCollection_HAsciiString) pippoName =
    //            new TCollection_HAsciiString("pippo");
    //        stepEntity->SetName(pippoName);

    //        //meshedEntities_->Append(entity);
    //        meshedEntities_->Append(stepEntity);

    //        nTransferEntity++;
    //    } //if Face
    //}

    //if (nTransferEntity == 0)
    //{
    //    FatalErrorInFunction
    //        << "No shapes have been transferred from CAD! "
    //        << abort(FatalError);
    //}

    Info<<"Found " << nTransferEntity << " faces" << endl;

    meshedEntitiesNTri_.setSize(nTransferEntity);

    //fix the shape
    //Handle(ShapeFix_Shape) sfs = new ShapeFix_Shape;
    //sfs->Init(reader_.OneShape());
    ////sfs->SetPrecision ( 1e-5 );
    ////sfs->SetMaxTolerance ( maxTol );
    ////sfs->SetMinTolerance ( mintol );
    //sfs->Perform();

    //Handle(ShapeFix_Wireframe) SFWF = new ShapeFix_Wireframe(reader_.OneShape());
    ////sfs->SetPrecision ( Prec );
    ////sfs->SetMaxTolerance ( maxTol );
    //SFWF->FixSmallEdges();
    //SFWF->FixWireGaps();

    //oneShape_ = sfs->Shape();
    //oneShape_ = SFWF->Shape();
    oneShape_ = reader_.OneShape();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::CADReader::CADReader
(
    const fileName& CADfileName
)
:
    inFile_(CADfileName),
    step_(false),
    iges_(false)
{
    Info<< "Loading CAD file " << inFile_ << nl << endl;

    meshedEntities_ = new TColStd_HSequenceOfTransient();
    oneShape_.Nullify();
    reader_.ClearShapes();

    word ext = CADfileName.ext();

    IFSelect_ReturnStatus readStat = IFSelect_RetDone;

    if (ext == "step" || ext == "stp")
    {
        step_ = true;

        STEPControl_Reader stepReader;

        readStat = stepReader.ReadFile(inFile_.c_str());

        if (readStat != IFSelect_RetDone)
        {
            FatalErrorInFunction
                << "Error while opening "
                << inFile_
                << abort(FatalError);
        }

        if (Pstream::master())
        {
            stepReader.PrintCheckLoad(/*failsonly*/!debug, IFSelect_ItemsByEntity);
        }

        reader_ = XSControl_Reader(stepReader.WS(), /*scratch*/false);
    }
    else if (ext == "iges" || ext == "igs")
    {
        iges_ = true;

        IGESControl_Reader igesReader;

        readStat = igesReader.ReadFile(inFile_.c_str());

        if (readStat != IFSelect_RetDone)
        {
            FatalErrorInFunction
                << "Error while opening "
                << inFile_
                << abort(FatalError);
        }

        if (Pstream::master())
        {
            igesReader.PrintCheckLoad(/*failsonly*/!debug, IFSelect_ItemsByEntity);
        }

        reader_ = XSControl_Reader(igesReader.WS(), /*scratch*/false);
    }
    else
    {
        FatalErrorInFunction
            << "Error while loading " << inFile_
            << "\nPlease use files of types \"STEP\" or \"IGES\"."
            << abort(FatalError);
    }

    this->TransferShapes();
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::fileFormats::CADReader::distributeShapes()
{
    const label nProcs = Pstream::nProcs();

    List<labelHashSet> shapesPerProc(nProcs);

    labelList cumTriPerProc(nProcs);

    scalar triPerProc = nTotalTriangles_ / nProcs;

    forAll(meshedEntitiesNTri_, i)
    {
        label procI = i % nProcs;

        if (i <= nProcs)
        {
            shapesPerProc[procI].insert(i);
            cumTriPerProc[procI] = meshedEntitiesNTri_[i];
        }
        else
        {
            scalar delta = abs(cumTriPerProc[procI] - triPerProc) / triPerProc;

            if (delta > 0.05)
            {
                shapesPerProc[procI].insert(i);
                cumTriPerProc[procI] += meshedEntitiesNTri_[i];
            }
        }
    }

    if (debug)
    {
        Info<<"\ncoarse triangles per processor:" << nl
            <<"wanted: " << triPerProc << nl
            <<"current: " << cumTriPerProc
            <<endl;
    }

    oneShape_.Nullify();

    label nTransferEntity = 0;

    TopoDS_Compound aCompound;

    BRep_Builder aBuilder;

    aBuilder.MakeCompound(aCompound);

    const Handle(XSControl_TransferReader)& TR = reader_.WS()->TransferReader();

    Standard_Integer nb = meshedEntities_->Length();

    TColStd_HSequenceOfTransient allMeshedEntities = *meshedEntities_;

    meshedEntities_.reset
    (
        new TColStd_HSequenceOfTransient()
    );

    for (Standard_Integer i = 1; i <= nb; i++)
    {
        Handle(Standard_Transient) entity = allMeshedEntities.Value(i);

        TopoDS_Shape aShape = TR->ShapeResult(entity);

        if (shapesPerProc[Pstream::myProcNo()][i])
        {
            aBuilder.Add(aCompound, aShape);

            meshedEntities_->Append(entity);

            nTransferEntity++;
        }
    }

    Pout<<"nTransferEntity " << nTransferEntity << endl;

    oneShape_ = aCompound;
}

void Foam::fileFormats::CADReader::distributeShapes
(
    const List<List<treeBoundBox>>& procBb
)
{
    oneShape_.Nullify();

    const treeBoundBox& bb = procBb[Pstream::myProcNo()][0];

    gp_Pnt myProcBbMin(bb.min().x(), bb.min().y(), bb.min().z());
    gp_Pnt myProcBbMax(bb.max().x(), bb.max().y(), bb.max().z());

    Bnd_Box myProcBb(myProcBbMin, myProcBbMax);

    label nTransferEntity = 0;

    TopoDS_Compound aCompound;

    BRep_Builder aBuilder;

    aBuilder.MakeCompound(aCompound);

    const Handle(XSControl_TransferReader)& TR = reader_.WS()->TransferReader();

    Standard_Integer nb = meshedEntities_->Length();

    TColStd_HSequenceOfTransient allMeshedEntities = *meshedEntities_;

    meshedEntities_.reset
    (
        new TColStd_HSequenceOfTransient()
    );

    for (Standard_Integer i = 1; i <= nb; i++)
    {
        Handle(Standard_Transient) entity = allMeshedEntities.Value(i);

        TopoDS_Shape aShape = TR->ShapeResult(entity);

        Bnd_Box aBndBox;

        BRepBndLib::Add(aShape, aBndBox);

        if (debug)
        {
            myProcBb.DumpJson(std::cout);
            aBndBox.DumpJson(std::cout);
            Info<< nl;
        }

        if (!aBndBox.IsOut(myProcBb))
        {
            aBuilder.Add(aCompound, aShape);

            meshedEntities_->Append(entity);

            nTransferEntity++;
        }
    }

    Pout<<"nTransferEntity " << nTransferEntity << endl;

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
    Info<< "\nTriangulating... (linDeflection "
        << linDeflection
        << "; angDeflection "
        << angDeflection
        << ";)"
        << nl << endl;

    //TODO create with MeshParameters
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

    const Handle(XSControl_TransferReader)& TR =
        reader_.WS()->TransferReader();

    Handle(Transfer_TransientProcess) transProc =
        TR->TransientProcess();

    Standard_Integer nb = meshedEntities_->Length();

    //count nodes and triangles
    for (Standard_Integer i = 1; i <= nb; i++)
    {
        Handle(Standard_Transient) entity = meshedEntities_->Value(i);

        if (!transProc->IsBound(entity))
        {
            continue;
        }

        TopoDS_Shape aFace = TR->ShapeResult(entity);

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
            point& min = bb.min();
            point& max = bb.max();

            boundingBox.Get(min.x(), min.y(), min.z(), max.x(), max.y(), max.z());

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

            meshedEntitiesNTri_[i-1] = aTriangulation->NbTriangles();
        }
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

    for (Standard_Integer i = 1; i <= nb; i++)
    {
        Handle(Standard_Transient) entity = meshedEntities_->Value(i);

        Handle(TCollection_HAsciiString) hName =
            new TCollection_HAsciiString("defaultFaces");

        if (step_)
        {
            Handle(StepRepr_RepresentationItem) stepEntity =
                Handle(StepRepr_RepresentationItem)::DownCast(entity);

            if (!stepEntity->Name().IsNull())
            {
                hName = stepEntity->Name();
            }

            DebugInformation<< "STEP " << i << " " << entity->DynamicType()->Name()
                            << " " << hName->ToCString() << endl;

            if (!transProc->IsBound(stepEntity))
            {
                Info<< "WARNING: found unbound entity "
                    << hName->ToCString() << endl;
                continue;
            }
        }
        else if (iges_)
        {
            Handle(IGESData_IGESEntity) igesEntity =
                Handle(IGESData_IGESEntity)::DownCast(entity);

            if (igesEntity->HasShortLabel())
            {
                // A short label is a non-blank 8-character string.
                hName = igesEntity->ShortLabel();
            }
            else
            {
                hName = igesEntity->NameValue();
            }

            DebugInformation<< "IGES " << i << " " << entity->DynamicType()->Name()
                            << " " << hName->ToCString() << endl;

            if (!transProc->IsBound(igesEntity))
            {
                Info<< "WARNING: found unbound entity "
                    << hName->ToCString() << endl;
                continue;
            }
        }

        TopoDS_Shape aFace = TR->ShapeResult(entity);

        TopLoc_Location aLoc;

        Handle(Poly_Triangulation) aTriangulation;

        try
        {
            OCC_CATCH_SIGNALS
            aTriangulation = BRep_Tool::Triangulation(TopoDS::Face(aFace), aLoc);
        }
        catch (Standard_TypeMismatch const&)
        {
            Info<< (step_ ? "STEP " : "IGES ") << i
                << " " << entity->DynamicType()->Name()
                << " " << hName->ToCString() << endl;
            Info<< "Can't convert to face...\n";
            ++aNbFacesNoTri;
            continue;
        }

        if (aTriangulation.IsNull())
        {
            ++aNbFacesNoTri;
            continue;
        }

        word patchName(hName->ToCString());
        patchName = stringOps::trim(patchName);

        if (!patchIDs.found(patchName))
        {
            Info<< "Found region: " << patchName << endl;
            patchIDs.insert(patchName, nPatches);

            nPatches++;
        }

        // copy nodes
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

        // copy triangles
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
