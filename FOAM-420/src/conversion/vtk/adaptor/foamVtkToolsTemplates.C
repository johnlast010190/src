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
    (c) 2017-2019 OpenCFD Ltd.

InNamespace
    Foam::vtk::Tools

\*---------------------------------------------------------------------------*/

// OpenFOAM includes
#include "db/error/error.H"

// VTK includes
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Face>
vtkSmartPointer<vtkCellArray>
Foam::vtk::Tools::Faces(const UList<Face>& faces)
{
    label numCells = faces.size();
    label maxCellSize = 0;
    for (const auto& f : faces)
    {
        maxCellSize = std::max(maxCellSize, f.size());
    }

    auto vtkCells = vtkSmartPointer<vtkCellArray>::New();
    vtkCells->AllocateEstimate(numCells, maxCellSize);

    for (const auto& f : faces)
    {
        vtkCells->InsertNextCell(f.size());
        for (vtkIdType pId : f)
        {
            vtkCells->InsertCellPoint(pId);
        }
    }
    vtkCells->Squeeze();

    return vtkCells;
}


template<class PatchType>
vtkSmartPointer<vtkPoints>
Foam::vtk::Tools::Patch::points(const PatchType& p)
{
    // Local patch points to vtkPoints
    return Tools::Points(p.localPoints());
}


template<class PatchType>
vtkSmartPointer<vtkCellArray>
Foam::vtk::Tools::Patch::faces(const PatchType& p)
{
    return Tools::Faces(p.localFaces());
}


template<class PatchType>
vtkSmartPointer<vtkPolyData>
Foam::vtk::Tools::Patch::mesh(const PatchType& p)
{
    auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

    vtkmesh->SetPoints(points(p));
    vtkmesh->SetPolys(faces(p));

    return vtkmesh;
}


template<class Face>
vtkSmartPointer<vtkPolyData>
Foam::vtk::Tools::Patch::mesh
(
    const UList<point>& pts,
    const UList<Face>& fcs
)
{
    auto vtkmesh = vtkSmartPointer<vtkPolyData>::New();

    vtkmesh->SetPoints(Tools::Points(pts));
    vtkmesh->SetPolys(Tools::Faces(fcs));

    return vtkmesh;
}


template<class PatchType>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::Tools::Patch::faceNormals(const PatchType& p)
{
    auto array = vtkSmartPointer<vtkFloatArray>::New();

    array->SetNumberOfComponents(3);
    array->SetNumberOfTuples(p.size());

    // Unit normals for patch faces.
    // Cached values if available or loop over faces (avoid triggering cache)

    vtkIdType faceId = 0;

    if (p.hasFaceNormals())
    {
        for (const vector& n : p.faceNormals())
        {
            array->SetTuple(faceId++, n.v_);
        }
    }
    else
    {
        for (const auto& f : p)
        {
            vector normal = f.normal(p.points());
            scalar magnitude = f.mag(p.points());
            vector unitNormal = normal/magnitude;
            const vector n(unitNormal);
            array->SetTuple(faceId++, n.v_);
        }
    }

    return array;
}


template<class PatchType>
vtkSmartPointer<vtkPoints>
Foam::vtk::Tools::Patch::faceCentres(const PatchType& p)
{
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    vtkpoints->SetNumberOfPoints(p.size());

    // Use cached values if available or loop over faces
    // (avoid triggering cache)

    vtkIdType pointId = 0;

    if (p.hasFaceCentres())
    {
        for (const point& pt : p.faceCentres())
        {
            vtkpoints->SetPoint(pointId++, pt.v_);
        }
    }
    else
    {
        for (const auto& f : p)
        {
            const point pt(f.centre(p.points()));
            vtkpoints->SetPoint(pointId++, pt.v_);
        }
    }

    return vtkpoints;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//
// Low-Level conversions
//

template<class Type>
Foam::label Foam::vtk::Tools::transcribeFloatData
(
    vtkFloatArray* array,
    const UList<Type>& input,
    vtkIdType start
)
{
    const int nComp(pTraits<Type>::nComponents);

    if (array->GetNumberOfComponents() != nComp)
    {
        FatalErrorInFunction
            << "vtk array '" << array->GetName()
            << "' has mismatch in number of components for type '"
            << pTraits<Type>::typeName
            << "' : target array has " << array->GetNumberOfComponents()
            << " components instead of " << nComp
            << nl;
    }

    const vtkIdType maxSize = array->GetNumberOfTuples();
    const vtkIdType endPos = start + vtkIdType(input.size());

    if (!maxSize)
    {
        // no-op
        return 0;
    }
    else if (start < 0 || vtkIdType(start) >= maxSize)
    {
        WarningInFunction
            << "vtk array '" << array->GetName()
            << "' copy with out-of-range [0," << label(maxSize) << ")"
            << " starting at " << label(start)
            << nl;

        return 0;
    }
    else if (endPos > maxSize)
    {
        WarningInFunction
            << "vtk array '" << string(array->GetName())
            << "' copy ends out-of-range (" << label(maxSize) << ")"
            << " using sizing (start,size) = ("
            << label(start) << "," << input.size() << ")"
            << nl;

        return 0;
    }

    float scratch[pTraits<Type>::nComponents];

    for (const Type& val : input)
    {
        foamToVtkTuple(scratch, val);
        array->SetTuple(start++, scratch);
    }

    return input.size();
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::Tools::zeroField
(
    const word& name,
    const label size
)
{
    auto data = vtkSmartPointer<vtkFloatArray>::New();

    data->SetName(name.c_str());
    data->SetNumberOfComponents(static_cast<int>(pTraits<Type>::nComponents));
    data->SetNumberOfTuples(size);

    // Fill() was not available before VTK-8
    #if (VTK_MAJOR_VERSION < 8)
    for (int i = 0; i < data->GetNumberOfComponents(); ++i)
    {
        data->FillComponent(i, 0);
    }
    #else
    data->Fill(0);
    #endif

    return data;
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::Tools::convertFieldToVTK
(
    const word& name,
    const UList<Type>& fld
)
{
    auto data = vtkSmartPointer<vtkFloatArray>::New();

    data->SetName(name.c_str());
    data->SetNumberOfComponents(static_cast<int>(pTraits<Type>::nComponents));
    data->SetNumberOfTuples(fld.size());

    transcribeFloatData(data, fld);

    return data;
}


// ************************************************************************* //
