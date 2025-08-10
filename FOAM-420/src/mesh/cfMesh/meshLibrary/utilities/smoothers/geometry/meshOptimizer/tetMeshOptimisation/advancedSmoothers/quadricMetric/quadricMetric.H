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

Class
    quadricMetric

Description
    Mesh smoothing using quadric metric. The metric is strictly convex

SourceFiles
    quadricMetric.C

\*---------------------------------------------------------------------------*/

#ifndef quadricMetric_H
#define quadricMetric_H

#include "utilities/smoothers/geometry/meshOptimizer/tetMeshOptimisation/advancedSmoothers/simplexSmoother/simplexSmoother.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTetMeshSimplex;

/*---------------------------------------------------------------------------*\
                           Class quadricMetric Declaration
\*---------------------------------------------------------------------------*/

class quadricMetric
:
    public simplexSmoother
{
    // Private data
        //- free vertex
        point& p_;

        //- normals of triangles forming the outer hull
        DynList<vector, 64> normals_;

        //- centres of triangles forming the outer hull
        DynList<point, 64> centres_;

    // Private member functions
        //- evaluate the value of the metric
        scalar evaluateMetric() const;

        //- evaluate metric gradients
        void evaluateGradients(vector& grad, tensor& gradGrad) const;

    public:

    // Constructor
        quadricMetric(partTetMeshSimplex& simplex);

    // Destructor
        ~quadricMetric();

    // Member functions
        //- find the best position such that all tets making a simplex
        //- have a positive volume
        void optimizeNodePosition(const scalar tol = 0.001);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
