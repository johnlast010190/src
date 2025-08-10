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

Class
    CADReader

Description
    Reads a CAD file in STEP or IGES format using OpenCASCADE libraries

SourceFiles
    CADReader.C

\*---------------------------------------------------------------------------*/
#ifndef CADREADER_H
#define CADREADER_H

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wnon-virtual-dtor"

#include <STEPControl_Reader.hxx>
#include <IGESControl_Reader.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
//#include <ShapeFix_Shape.hxx>
#include <XSControl_WorkSession.hxx>
#include <XSControl_TransferReader.hxx>
#include <Transfer_TransientProcess.hxx>
#include <IGESData_IGESEntity.hxx>
#include <TransferBRep.hxx>
#include <StepRepr_RepresentationItem.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>
#include <StlAPI_Writer.hxx>
#include <TopExp.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
//#include <BRepMesh_Context.hxx>
//#include <BRepMesh_FaceDiscret.hxx>
//#include <BRepMesh_DelabellaMeshAlgoFactory.hxx>

#include <BRepBndLib.hxx>
#include <Bnd_Box.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS_Compound.hxx>

#pragma GCC diagnostic pop

#include "primitives/Scalar/scalar/scalar.H"
#include "primitives/strings/fileName/fileName.H"
#include "meshes/primitiveShapes/point/pointField.H"
#include "meshes/treeBoundBox/treeBoundBox.H"
#include "meshes/meshShapes/labelledTri/labelledTri.H"
#include "containers/HashTables/HashSet/HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

/*---------------------------------------------------------------------------*\
                   Class fileFormats::CADReader Declaration
\*---------------------------------------------------------------------------*/

class CADReader
{

private:

    // Private data

        fileName inFile_;

        word nameRegionsBy_;

        bool step_, iges_;

        TopoDS_Shape oneShape_;

        TopTools_ListOfShape meshedShapes_;

        List<word> regionNames_;

        scalarField regionAreas_;

        //- store number of triangles for each meshedShapes
        labelList meshedEntitiesNTri_;

        label nTotalTriangles_;

    // Private Member Functions

public:

    //- Runtime type information
    ClassName("CADReader");

    // Constructors

        CADReader
        (
            const fileName& CADfileName,
            const word& nameRegionsBy,
            const bool distributed = false,
            const wordHashSet& ignoreRegions = wordHashSet()
        );

    // Member Functions

        void distributeShapes
        (
            const bool balanceArea = false
        );

        //- Transfers to oneShape_ only shapes inside procBb[Pstream::myProcNo()]
        void distributeShapes
        (
            const List<List<treeBoundBox>>& procBb
        );

        void triangulate
        (
            scalar linDeflection,
            scalar angDeflection,
            pointField& surfPoints,
            List<labelledTri>& surfTriangles,
            HashTable<label>& patchIDs
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fileFormats
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif // CADREADER_H
