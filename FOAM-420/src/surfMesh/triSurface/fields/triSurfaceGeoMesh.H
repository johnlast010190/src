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
    (c) 2011 OpenFOAM Foundation

Class
    Foam::triSurfaceGeoMesh

Description
    The triSurface GeoMesh (for holding fields).

    Similar to the volGeoMesh used for the Finite Volume discretization.

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceGeoMesh_H
#define triSurfaceGeoMesh_H

#include "meshes/GeoMesh/GeoMesh.H"
#include "triSurface/triSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class triSurfaceGeoMesh
:
    public GeoMesh<triSurface>
{

public:

    // Constructors

        //- Construct from triSurface reference
        explicit triSurfaceGeoMesh(const triSurface& mesh)
        :
            GeoMesh<triSurface>(mesh)
        {}


    // Member Functions

        //- Return size
        static label size(const triSurface& mesh)
        {
            return mesh.size();
        }

        //- Return size
        label size() const
        {
            return size(mesh_);
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
