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
    partTriMeshSimplex

Description
    A simplex which is used for smoothing purposes

SourceFiles


\*---------------------------------------------------------------------------*/

#ifndef partTriMeshSimplex_H
#define partTriMeshSimplex_H

#include "utilities/meshes/partTriMesh/partTriMesh.H"
#include "meshes/meshShapes/triFace/triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class partTriMesh;

/*---------------------------------------------------------------------------*\
                           Class partTriMeshSimplex Declaration
\*---------------------------------------------------------------------------*/

class partTriMeshSimplex
{
    // Private data
        //- points making the simplex
        DynList<point, 32> pts_;

        //- triangles making the simplex
        DynList<triFace, 32> trias_;

    public:

    // Constructors
        //- Construct from partTriMesh and point label
        partTriMeshSimplex(const partTriMesh& tm, const label pI);

    // Destructor
        ~partTriMeshSimplex();

    // Member functions
        //- return points
        inline DynList<point, 32>& pts()
        {
            return pts_;
        }

        //- return points
        inline const DynList<point, 32>& pts() const
        {
            return pts_;
        }

        //- return triangles
        inline const DynList<triFace, 32>& triangles() const
        {
            return trias_;
        }

        //- return centre point coordinates
        inline const point& centrePoint() const
        {
            return pts_[trias_[0][0]];
        }

        //- return the normal of the simplex
        vector normal() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
