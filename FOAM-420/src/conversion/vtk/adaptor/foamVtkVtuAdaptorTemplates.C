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

// VTK includes
#include "vtkFloatArray.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::vtuAdaptor::convertField
(
    const DimensionedField<Type, volMesh>& fld,
    const vtuAdaptor& vtuData
)
{
    const int nComp(pTraits<Type>::nComponents);
    const labelUList& cellMap = vtuData.cellMap();

    auto data = vtkSmartPointer<vtkFloatArray>::New();
    data->SetName(fld.name().c_str());
    data->SetNumberOfComponents(nComp);
    data->SetNumberOfTuples(cellMap.size());

    // DebugInformation    //     << "Convert field: " << fld.name()
    //     << " size=" << cellMap.size()
    //     << " (" << fld.size() << " + "
    //     << (cellMap.size() - fld.size())
    //     << ") nComp=" << nComp << endl;


    float scratch[pTraits<Type>::nComponents];

    vtkIdType celli = 0;
    for (const label meshCelli : cellMap)
    {
        vtk::Tools::foamToVtkTuple(scratch, fld[meshCelli]);
        data->SetTuple(celli++, scratch);
    }

    return data;
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::vtuAdaptor::convertField
(
    const GeometricField<Type, fvPatchField, volMesh>& fld,
    const vtuAdaptor& vtuData
)
{
    return convertField<Type>(fld.internalField(), vtuData);
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::vtuAdaptor::convertField
(
    const DimensionedField<Type, volMesh>& fld
) const
{
    return convertField<Type>(fld, *this);
}


template<class Type>
vtkSmartPointer<vtkFloatArray>
Foam::vtk::vtuAdaptor::convertField
(
    const GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    return convertField<Type>(fld, *this);
}


// ************************************************************************* //
