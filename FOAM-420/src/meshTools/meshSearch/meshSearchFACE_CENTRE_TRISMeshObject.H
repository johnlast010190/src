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
    (c) 2011-2015 OpenFOAM Foundation

Class
    Foam::meshSearchFACE_CENTRE_TRISMeshObject

Description
    MeshObject wrapper around meshSearch(mesh,  polyMesh::FACE_CENTRE_TRIS).

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef meshSearchFACE_CENTRE_TRISMeshObject_H
#define meshSearchFACE_CENTRE_TRISMeshObject_H

#include "meshes/MeshObject/MeshObject.H"
#include "meshSearch/meshSearch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class meshSearchFACE_CENTRE_TRISMeshObject Declaration
\*---------------------------------------------------------------------------*/

class meshSearchFACE_CENTRE_TRISMeshObject
:
    public MeshObject
    <
        polyMesh,
        GeometricMeshObject,
        meshSearchFACE_CENTRE_TRISMeshObject
    >,
    public meshSearch
{

public:

    // Declare name of the class and its debug switch
    TypeName("meshSearchFACE_CENTRE_TRISMeshObject");


    // Constructors

        //- Constructor given polyMesh
        explicit meshSearchFACE_CENTRE_TRISMeshObject(const polyMesh& mesh);

    //- Destructor
    virtual ~meshSearchFACE_CENTRE_TRISMeshObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
