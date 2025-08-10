/*---------------------------------------------------------------------------*\
|       o        |
|    o     o     |  FOAM (R) : Open-source CFD for Enterprise
|   o   O   o    |  Version : 4.2.0
|    o     o     |  ESI Ltd. <http://esi.com/>
|       o        |
\*---------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Modified code Copyright (C) 2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
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
    (c) 2013-2016 OpenFOAM Foundation

Class
    Foam::meshStructure

Description
    Detect extruded mesh structure given a set of faces
    (uindirectPrimitivePatch).

    All indices to originating set are in terms of global faces,
    global edges, global points.

SourceFiles
    meshStructure.C

\*---------------------------------------------------------------------------*/

#ifndef meshStructure_H
#define meshStructure_H

#include "primitives/ints/lists/labelList.H"
#include "meshes/primitiveMesh/primitivePatch/uindirectPrimitivePatch.H"
#include "db/typeInfo/className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class globalIndex;

/*---------------------------------------------------------------------------*\
                         Class meshStructure Declaration
\*---------------------------------------------------------------------------*/

class meshStructure
{
    // Private data

        //- Cell to patch face
        labelList cellToPatchFaceAddressing_;

        //- Cell to layer
        labelList cellLayer_;

        //- Face to patch face
        labelList faceToPatchFaceAddressing_;

        //- Face to patch edge
        labelList faceToPatchEdgeAddressing_;

        //- Face to layer
        labelList faceLayer_;

        //- Point to patch point
        labelList pointToPatchPointAddressing_;

        //- Point to layer
        labelList pointLayer_;

        //- Is mesh structured?
        bool structured_;


   // Private Member Functions

        //- Is cell structured
        bool isStructuredCell
        (
            const polyMesh& mesh,
            const label layerI,
            const label celli
        ) const;

        //- Calculate all maps.
        void correct
        (
            const polyMesh& mesh,
            const uindirectPrimitivePatch& pp,
            const globalIndex& globalFaces,
            const globalIndex& globalEdges,
            const globalIndex& globalPoints
        );


public:

    // Declare name of the class and its debug switch
    ClassName("meshStructure");

    // Constructors

        //- Construct from mesh and faces in mesh. Any addressing to
        //  faces/edges/points on patch are global indices
        meshStructure(const polyMesh& mesh, const uindirectPrimitivePatch&);

        //- Construct from mesh and faces in mesh and global indexing
        meshStructure
        (
            const polyMesh& mesh,
            const uindirectPrimitivePatch&,
            const globalIndex& globalFaces,
            const globalIndex& globalEdges,
            const globalIndex& globalPoints
        );


    // Member Functions

        //- Is mesh structured?
        inline bool structured() const;

        //- Cell to (global) patch face
        inline const labelList& cellToPatchFaceAddressing() const;

        //- Cell to (global) patch face
        inline labelList& cellToPatchFaceAddressing();

        //- Cell to layer
        inline const labelList& cellLayer() const;

        //- Cell to layer
        inline labelList& cellLayer();

        //- Face to (global) patch face
        inline const labelList& faceToPatchFaceAddressing() const;

        //- Face to (global) patch face
        inline labelList& faceToPatchFaceAddressing();

        //- Face to (global) patch edge
        inline const labelList& faceToPatchEdgeAddressing() const;

        //- Face to (global) patch edge
        inline labelList& faceToPatchEdgeAddressing();

        //- Face to layer
        inline const labelList& faceLayer() const;

        //- Face to layer
        inline labelList& faceLayer();

        //- Point to (global) patch point
        inline const labelList& pointToPatchPointAddressing() const;

        //- Point to (global) patch point
        inline labelList& pointToPatchPointAddressing();

        //- Point to layer
        inline const labelList& pointLayer() const;

        //- Point to layer
        inline labelList& pointLayer();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "meshStructure/meshStructureI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
