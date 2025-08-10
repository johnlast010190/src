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
    (c) Creative Fields, Ltd.

Authors
    Franjo Juretic (franjo.juretic@c-fields.com)

Description

\*---------------------------------------------------------------------------*/

#include "utilities/surfaceTools/meshSurfaceEdgeExtractorFUN/meshSurfaceEdgeExtractorFUN.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "include/demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// #define DEBUGSearch

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshSurfaceEngine& meshSurfaceEdgeExtractorFUN::surfaceEngine()
{
    # ifdef USE_OMP
    if (omp_in_parallel())
        FatalErrorInFunction
            << "Cannot create surface engine with a parallel region"
            << exit(FatalError);
    # endif

    if (!surfaceEnginePtr_)
        surfaceEnginePtr_ = new meshSurfaceEngine(mesh_);

    return *surfaceEnginePtr_;
}

void meshSurfaceEdgeExtractorFUN::clearOut()
{
    deleteDemandDrivenData(surfaceEnginePtr_);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh and octree
meshSurfaceEdgeExtractorFUN::meshSurfaceEdgeExtractorFUN
(
    polyMeshGen& mesh,
    const meshOctree& octree,
    const bool createWrapperSheet
)
:
    mesh_(mesh),
    meshOctree_(octree),
    surfaceEnginePtr_(nullptr),
    createWrapperSheet_(createWrapperSheet)
{
    if (Pstream::parRun())
        FatalErrorInFunction
            << "Cannot run in parallel!" << exit(FatalError);

    createBasicFundamentalSheets();

    smoothMeshSurface();

    remapBoundaryPoints();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

meshSurfaceEdgeExtractorFUN::~meshSurfaceEdgeExtractorFUN()
{
    deleteDemandDrivenData(surfaceEnginePtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
