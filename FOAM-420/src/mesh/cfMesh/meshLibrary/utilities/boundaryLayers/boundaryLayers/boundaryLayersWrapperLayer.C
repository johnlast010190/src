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

#include "utilities/boundaryLayers/boundaryLayers/boundaryLayers.H"
#include "utilities/surfaceTools/meshSurfaceEngine/meshSurfaceEngine.H"
#include "include/demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void boundaryLayers::addWrapperLayer()
{
    createOTopologyLayers();

    if (treatedPatch_[0]) return;

    const meshSurfaceEngine& mse = surfaceEngine();

    const labelList& bPoints = mse.boundaryPoints();

    boolList treatPatches(mesh_.boundaries().size(), true);

    labelLongList newLabelForVertex(nPoints_, -1);

    pointFieldPMG& points = mesh_.points();
    points.setSize(points.size() + bPoints.size());
    forAll(bPoints, bpI)
    {
        points[nPoints_] = points[bPoints[bpI]];
        newLabelForVertex[bPoints[bpI]] = nPoints_++;
    }

    createNewFacesAndCells(treatPatches);

    forAll(treatPatches, patchI)
        if (treatPatches[patchI])
            treatedPatch_[patchI] = true;

    //- delete surface engine
    clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
