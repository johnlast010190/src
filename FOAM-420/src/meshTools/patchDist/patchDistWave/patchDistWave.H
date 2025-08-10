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
    (c) 2011-2022 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

Class
    Foam::patchDistWave

Description
    Takes a set of patches to start FaceCellWave from and computes the distance
    at patches and possibly additional transported data.

SourceFiles
    patchDistWaveTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef patchDistWave_H
#define patchDistWave_H

#include "algorithms/FaceCellWave/FaceCellWave.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class wallPoint;

namespace patchDistWave
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Set initial set of changed faces
template<class PatchPointType, class ... InitialPatchData>
void setChangedFaces
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    labelList& changedFaces,
    List<PatchPointType>& changedFacesInfo,
    const InitialPatchData& ... initialPatchData
);

//- Copy FaceCellWave values into the cell fields
template<class PatchPointType, class DataType, class DataMethod>
label getCellValues
(
    FaceCellWave<PatchPointType>& waveInfo,
    Field<DataType>& cellValues,
    DataMethod method,
    const DataType& stabiliseValue = pTraits<DataType>::zero
);

//- Wave distance data from the patches to the cells
template<class PatchPointType>
label wave
(
    const polyMesh& mesh,
    const labelHashSet& patchIDs,
    scalarField& cellDistance,
    bool correct = true
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace patchDistWave
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "patchDist/patchDistWave/patchDistWaveTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
