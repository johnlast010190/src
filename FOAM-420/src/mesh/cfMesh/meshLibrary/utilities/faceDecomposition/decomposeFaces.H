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
    decomposeFaces

Description
    Decompose selected mesh faces into triangles using midnode subdivision

SourceFiles
    decomposeFaces.C

\*---------------------------------------------------------------------------*/

#ifndef decomposeFaces_H
#define decomposeFaces_H

#include "utilities/meshes/polyMeshGenModifier/polyMeshGenModifier.H"
#include "utilities/containers/DynList/DynList.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class decomposeFaces Declaration
\*---------------------------------------------------------------------------*/

class decomposeFaces
{
    // private data
        //- reference to the mesh
        polyMeshGen& mesh_;

        //- number of points
        VRWGraph newFacesForFace_;

        //- is decomposition performed
        bool done_;

        //- disallows bitwise construct
        void operator=(const decomposeFaces&);

        //- copy constructor
        decomposeFaces(const decomposeFaces&);

public:

    // Constructors

        //- construct components
        decomposeFaces(polyMeshGen& mesh);

        //- Destructor
        ~decomposeFaces();

    // Member functions
        //- decompose selected faces into triangles using midnode subdivision
        void decomposeMeshFaces(const boolList& decomposeFace);

        //- decompose internal faces containing concave nodes
        void decomposeConcaveInternalFaces(const boolList& concaveVertex);

        //- reference to a graph containg label of new faces generated
        //- from each face
        const VRWGraph& newFacesForFace() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
