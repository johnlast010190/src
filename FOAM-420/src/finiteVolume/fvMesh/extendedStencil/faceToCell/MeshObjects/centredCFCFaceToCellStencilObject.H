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
    (c) 2011-2013 OpenFOAM Foundation

Class
    Foam::centredCFCFaceToCellStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCFCFaceToCellStencilObject_H
#define centredCFCFaceToCellStencilObject_H

#include "fvMesh/extendedStencil/faceToCell/extendedCentredFaceToCellStencil.H"
#include "fvMesh/extendedStencil/faceToCell/globalIndexStencils/CFCFaceToCellStencil.H"
#include "meshes/MeshObject/MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class centredCFCFaceToCellStencilObject Declaration
\*---------------------------------------------------------------------------*/

class centredCFCFaceToCellStencilObject
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        centredCFCFaceToCellStencilObject
    >,
    public extendedCentredFaceToCellStencil
{

public:

    TypeName("centredCFCFaceToCellStencil");

    // Constructors

        //- Construct from uncompacted face stencil
        explicit centredCFCFaceToCellStencilObject
        (
            const fvMesh& mesh
        )
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                centredCFCFaceToCellStencilObject
            >(mesh),
            extendedCentredFaceToCellStencil(CFCFaceToCellStencil(mesh))
        {}


    //- Destructor
    virtual ~centredCFCFaceToCellStencilObject()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
