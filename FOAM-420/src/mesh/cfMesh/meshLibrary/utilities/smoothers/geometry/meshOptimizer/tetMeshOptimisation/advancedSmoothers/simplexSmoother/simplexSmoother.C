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

#include "utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/simplexSmoother/simplexSmoother.H"
#include "meshes/primitiveShapes/tetrahedron/tetrahedron.H"
#include "utilities/meshes/partTetMesh/partTetMeshSimplex.H"

//#define DEBUGSmooth

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simplexSmoother::simplexSmoother(partTetMeshSimplex& simplex)
:
    points_(simplex.pts()),
    tets_(simplex.tets()),
    pointI_(tets_[0][3]),
    bb_()
{
    point min(VGREAT, VGREAT, VGREAT), max(-VGREAT, -VGREAT, -VGREAT);
    forAll(tets_, tetI)
    {
        const partTet& pt = tets_[tetI];
        const tetrahedron<point, point> tet
        (
            points_[pt.a()],
            points_[pt.b()],
            points_[pt.c()],
            points_[pt.d()]
        );

        min = Foam::min(min, tet.a());
        max = Foam::max(max, tet.a());

        min = Foam::min(min, tet.b());
        max = Foam::max(max, tet.b());

        min = Foam::min(min, tet.c());
        max = Foam::max(max, tet.c());
    }

    bb_.max() = max;
    bb_.min() = min;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simplexSmoother::~simplexSmoother()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Member functions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
