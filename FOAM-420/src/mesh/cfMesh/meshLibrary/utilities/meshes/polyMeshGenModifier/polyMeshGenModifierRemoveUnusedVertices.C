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

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "include/demandDrivenData.H"
#include "primitives/ints/lists/labelList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::removeUnusedVertices()
{
    faceListPMG& faces = mesh_.faces_;
    pointFieldPMG& points = mesh_.points_;

    boolList usePoint(points.size(), false);
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        forAll(f, pI)
            usePoint[f[pI]] = true;
    }

    labelLongList newLabel(points.size(), -1);
    label nPoints(0);
    forAll(points, pI)
        if (usePoint[pI])
            newLabel[pI] = nPoints++;

    //- remove unused points from the list
    forAll(newLabel, pI)
        if ((newLabel[pI] != -1) && (newLabel[pI] < pI))
        {
            points[newLabel[pI]] = points[pI];
        }

    points.setSize(nPoints);

    forAll(faces, faceI)
    {
        face& f = faces[faceI];

        forAll(f, pI)
            f[pI] = newLabel[f[pI]];
    }

    mesh_.updatePointSubsets(newLabel);

    mesh_.clearOut();
    this->clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
