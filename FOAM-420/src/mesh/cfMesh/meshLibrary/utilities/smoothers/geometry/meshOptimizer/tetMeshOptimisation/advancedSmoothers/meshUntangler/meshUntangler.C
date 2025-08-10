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

#include "include/demandDrivenData.H"
#include "utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/meshUntangler/meshUntangler.H"
#include "meshes/primitiveShapes/plane/plane.H"

//#define DEBUGSmooth

#ifdef DEBUGSmooth
#include "db/Time/Time.H"
#include "db/objectRegistry/objectRegistry.H"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

meshUntangler::meshUntangler
(
    partTetMeshSimplex& simplex
)
:
    simplexSmoother(simplex)
{
}

meshUntangler::~meshUntangler()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void meshUntangler::optimizeNodePosition(const scalar /*tol*/)
{
    # ifdef DEBUGSmooth
    Info<< "Untangling point " << pointI_ << endl;
    # endif

    cutRegion cr(bb_);

    forAll(tets_, tetI)
    {
        const partTet& tet = tets_[tetI];
        vector n
        (
            (points_[tet.b()] - points_[tet.a()]) ^
            (points_[tet.c()] - points_[tet.a()])
        );

        if (mag(n) < VSMALL) continue;

        plane pl(points_[tet.a()], n);

        # ifdef DEBUGSmooth
        Info<< "tet.a() " << tet.a() << endl;
        Info<< "Cutting plane ref point " << pl.refPoint() << endl;
        Info<< "Cutting plane normal " << pl.normal() << endl;
        # endif

        cr.planeCut(pl);
    }

    if (cr.points().size())
    {
        point p(vector::zero);

        const DynList<point, 64>& pts = cr.points();
        forAll(pts, pI)
            p += pts[pI];

        p /= pts.size();

        # ifdef DEBUGSmooth
        Info<< "Corners of the feasible region " << pts << endl;
        # endif

        for (direction i=0;i<vector::nComponents;++i)
        {
            const scalar& val = p[i];
            if ((val != val) || ((val - val) != (val - val)))
                return;
        }

        points_[pointI_] = p;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
