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
    Foam::refinementIterator

Description
    Utility class to do iterating meshCutter until all requests satisfied.

    Needed since cell cutting can only cut cell once in one go so if
    refinement pattern is not compatible on a cell by cell basis it will
    refuse to cut.

    Parallel: communicates. All decisions done on 'reduce'd variable.

SourceFiles
    refinementIterator.C

\*---------------------------------------------------------------------------*/

#ifndef refinementIterator_H
#define refinementIterator_H

#include "meshCut/edgeVertex/edgeVertex.H"
#include "primitives/ints/lists/labelList.H"
#include "containers/HashTables/Map/Map.H"
#include "db/typeInfo/typeInfo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class refineCell;
class undoableMeshCutter;
class cellLooper;

/*---------------------------------------------------------------------------*\
                           Class refinementIterator Declaration
\*---------------------------------------------------------------------------*/

class refinementIterator
:
    public edgeVertex
{
    // Private data

        //- Reference to mesh
        polyMesh& mesh_;

        //- Reference to refinementEngine
        undoableMeshCutter& meshRefiner_;

        //- Reference to object to walk individual cells
        const cellLooper& cellWalker_;

        //- Whether to write intermediate meshes
        bool writeMesh_;


public:

    //- Runtime type information
    ClassName("refinementIterator");


    // Constructors

        //- Construct from mesh, refinementEngine and cell walking routine.
        //  If writeMesh = true increments runTime and writes intermediate
        //  meshes.
        refinementIterator
        (
            polyMesh& mesh,
            undoableMeshCutter& meshRefiner,
            const cellLooper& cellWalker,
            const bool writeMesh = false
        );


    //- Destructor
    ~refinementIterator();


    // Member Functions

        //- Try to refine cells in given direction. Constructs intermediate
        //  meshes. Returns map from old to added cells.
        Map<label> setRefinement(const List<refineCell>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
