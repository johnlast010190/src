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
    (C) 2019 OpenCFD Ltd.

Class
    Foam::polySurfaceGeoMesh

Description
    The polySurface GeoMesh (for holding fields).

\*---------------------------------------------------------------------------*/

#ifndef polySurfaceGeoMesh_H
#define polySurfaceGeoMesh_H

#include "polySurface/polySurface.H"
#include "meshes/GeoMesh/GeoMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class polySurfaceGeoMesh Declaration
\*---------------------------------------------------------------------------*/

class polySurfaceGeoMesh
:
    public GeoMesh<polySurface>
{
public:

    // Constructors

        //- Construct from polySurface reference
        explicit polySurfaceGeoMesh(const polySurface& mesh)
        :
            GeoMesh<polySurface>(mesh)
        {}


    // Member Functions

        //- Return size
        static label size(const polySurface& mesh)
        {
            return mesh.faces().size();
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
