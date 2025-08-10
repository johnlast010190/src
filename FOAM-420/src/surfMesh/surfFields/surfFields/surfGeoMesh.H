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
    Foam::surfGeoMesh

Description
    The surfMesh GeoMesh (for holding fields).

    Similar to the volGeoMesh used for the Finite Volume discretization.

\*---------------------------------------------------------------------------*/

#ifndef surfGeoMesh_H
#define surfGeoMesh_H

#include "meshes/GeoMesh/GeoMesh.H"
#include "surfMesh/surfMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class surfGeoMesh Declaration
\*---------------------------------------------------------------------------*/

class surfGeoMesh
:
    public GeoMesh<surfMesh>
{

public:

    // Constructors

        //- Construct from surfMesh reference
        explicit surfGeoMesh(const surfMesh& mesh)
        :
            GeoMesh<surfMesh>(mesh)
        {}


    // Member Functions

        //- Return size
        static label size(const surfMesh& mesh)
        {
            return mesh.nFaces();
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
