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
    (c) 2014 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "hessianMeshOptimization/meshOptimizationMethod/fullMeshOptimization/fullMeshOptimization.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

namespace Foam
{
    defineTypeNameAndDebug(fullMeshOptimization, 0);
    addToRunTimeSelectionTable
    (
        meshOptimizationMethod,
        fullMeshOptimization,
        dictionary
    );

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //
void Foam::fullMeshOptimization::calculateActiveSet()
{
    const labelList& pS = identity(mesh_.nPoints());
    const labelList& fS = identity(mesh_.nFaces());
    const labelList& cS = identity(mesh_.nCells());
    const labelList& eS = identity(mesh_.nEdges());
    activeSet_.update(pS, fS, cS, fS, eS);
}
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fullMeshOptimization::fullMeshOptimization
(
    const dictionary& dict,
    const fvMesh& mesh,
    const meshMetric& metric
)
:
    meshOptimizationMethod(dict, mesh, metric),
    mesh_(mesh)
{
    calculateActiveSet();
}

fullMeshOptimization::~fullMeshOptimization()
{
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //
bool Foam::fullMeshOptimization::updateActiveSet()
{
    //Do nothing for now
    return false;
}

bool Foam::fullMeshOptimization::converge()
{
    return false;
}

void fullMeshOptimization::merge(const labelList& pointSet)
{
    activeSet_.updatePointSet(pointSet);
}
} /* namespace Foam */
