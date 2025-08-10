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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline vtkSmartPointer<vtkPoints>
Foam::vtk::vtuAdaptor::points
(
    const fvMesh& mesh
) const
{
    // Convert OpenFOAM mesh vertices to VTK
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    // Normal points
    const pointField& pts = mesh.points();

    // Additional cell centres
    const labelUList& addPoints = this->additionalIds();

    vtkpoints->SetNumberOfPoints(pts.size() + addPoints.size());

    // Normal points
    vtkIdType pointId = 0;
    for (const point& p : pts)
    {
        vtkpoints->SetPoint(pointId++, p.v_);
    }

    // Cell centres
    for (const label meshCelli : addPoints)
    {
        vtkpoints->SetPoint(pointId++, mesh.cellCentres()[meshCelli].v_);
    }

    return vtkpoints;
}


inline vtkSmartPointer<vtkPoints>
Foam::vtk::vtuAdaptor::points
(
    const fvMesh& mesh,
    const labelUList& pointMap
) const
{
    // Convert OpenFOAM mesh vertices to VTK
    auto vtkpoints = vtkSmartPointer<vtkPoints>::New();

    // Normal points
    const pointField& pts = mesh.points();

    // Additional cell centres
    const labelUList& addPoints = this->additionalIds();

    vtkpoints->SetNumberOfPoints(pointMap.size() + addPoints.size());

    // Normal points
    vtkIdType pointId = 0;
    for (const label meshPointi : pointMap)
    {
        vtkpoints->SetPoint(pointId++, pts[meshPointi].v_);
    }

    // Cell centres
    for (const label meshCelli : addPoints)
    {
        vtkpoints->SetPoint(pointId++, mesh.cellCentres()[meshCelli].v_);
    }

    return vtkpoints;
}


inline vtkSmartPointer<vtkUnstructuredGrid>
Foam::vtk::vtuAdaptor::internal
(
    const fvMesh& mesh,
    const bool decompPoly
)
{
    vtk::vtuSizing sizing(mesh, decompPoly);

    auto cellTypes = vtkSmartPointer<vtkUnsignedCharArray>::New();

    auto cells = vtkSmartPointer<vtkCellArray>::New();
    auto faces = vtkSmartPointer<vtkIdTypeArray>::New();

    auto cellLocations = vtkSmartPointer<vtkIdTypeArray>::New();
    auto faceLocations = vtkSmartPointer<vtkIdTypeArray>::New();

    UList<uint8_t> cellTypesUL =
        vtk::Tools::asUList(cellTypes, sizing.nFieldCells());

    List<vtkIdType> cellsL(sizing.sizeInternal(vtk::vtuSizing::slotType::CELLS));
    cells->AllocateEstimate(sizing.nFieldCells(), 10);

    UList<vtkIdType> cellLocationsUL =
        vtk::Tools::asUList
        (
            cellLocations,
            sizing.sizeInternal(vtk::vtuSizing::slotType::CELLS_OFFSETS)
        );

    UList<vtkIdType> facesUL =
        vtk::Tools::asUList
        (
            faces,
            sizing.sizeInternal(vtk::vtuSizing::slotType::FACES)
        );

    UList<vtkIdType> faceLocationsUL =
        vtk::Tools::asUList
        (
            faceLocations,
            sizing.sizeInternal(vtk::vtuSizing::slotType::FACES_OFFSETS)
        );


    sizing.populateInternal
    (
        mesh,
        cellTypesUL,
        cellsL,
        cellLocationsUL,
        facesUL,
        faceLocationsUL,
        static_cast<foamVtkMeshMaps&>(*this)
    );

    auto iter = cellsL.begin();
    const auto& end = cellsL.end();
    while (iter != end)
    {
        int cellSize = int(*iter);
        iter++;
        cells->InsertNextCell(cellSize);
        for (vtkIdType i = 0; i < cellSize; i++)
        {
            cells->InsertCellPoint(*iter);
            iter++;
        }
    }
    cells->Squeeze();


    auto vtkmesh = vtkSmartPointer<vtkUnstructuredGrid>::New();

    // Convert OpenFOAM mesh vertices to VTK
    // - can only do this *after* populating the decompInfo with cell-ids
    //   for any additional points (ie, mesh cell-centres)
    vtkmesh->SetPoints(this->points(mesh));

    if (!facesUL.empty())
    {
        vtkmesh->SetCells
        (
            cellTypes,
            cells,
            faceLocations,
            faces
        );
    }
    else
    {
        vtkmesh->SetCells
        (
            cellTypes,
            cells,
            nullptr,
            nullptr
        );
    }

    return vtkmesh;
}


// ************************************************************************* //
