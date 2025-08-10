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

#include "utilities/triSurfaceTools/triSurfaceDetectFeatureEdges/triSurfaceDetectFeatureEdges.H"
#include "utilities/helperFunctions/helperFunctions.H"
#include "include/demandDrivenData.H"
#include "primitives/Pair/labelPair.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void triSurfaceDetectFeatureEdges::detectFeatureEdgesAngleCriterion()
{
    const scalar tol = Foam::cos(angleTolerance_*M_PI/180.0);

    const vectorField& normals = surf_.facetNormals();

    const VRWGraph& edgeFaces = surf_.edgeFacets();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 40)
    # endif
    forAll(edgeFaces, edgeI)
    {
        const constRow eFaces = edgeFaces[edgeI];

        if (edgeFaces.sizeOfRow(edgeI) != 2)
        {
            featureEdges_[edgeI] |= 8;
            continue;
        }

        scalar cosAngle =
            (normals[eFaces[0]] & normals[eFaces[1]]) /
            (mag(normals[eFaces[0]]) * mag(normals[eFaces[1]]) + VSMALL);

        //- check the orientation of triangles at this edge
        //- check the sign of the angle if the orientation  is not consistent
        const labelledTri& tri0 = surf_[edgeFaces(edgeI, 0)];
        const labelledTri& tri1 = surf_[edgeFaces(edgeI, 1)];
        DynList<labelPair> sharedIndices;
        forAll(tri0, i)
        {
            forAll(tri1, j)
            {
                if (tri0[i] == tri1[j])
                    sharedIndices.append(labelPair(i, j));
            }
        }

        if (sharedIndices.size() == 2)
        {
            const labelPair& pair0 = sharedIndices[0];
            const labelPair& pair1 = sharedIndices[1];
            if (((pair0.first() + 1) % 3) == pair1.first())
            {
                if ((pair0.second() + 1) % 3 == pair1.second())
                    cosAngle *= -1.0;
            }
            else
            {
                if ((pair1.second() + 1) % 3 == pair0.second())
                    cosAngle *= -1.0;
            }
        }

        if (cosAngle < tol)
            featureEdges_[edgeI] |= 1;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
