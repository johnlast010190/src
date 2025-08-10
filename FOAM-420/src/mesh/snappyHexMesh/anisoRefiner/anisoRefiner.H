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
    (c) 2014-1014 Esi Ltd.

Class
    Foam::anisoRefiner

Description
    Anisotropic refinement of hex and split hex cells

SourceFiles
    anisoRefiner.C

\*---------------------------------------------------------------------------*/

#ifndef anisoRefiner_H
#define anisoRefiner_H

#include "primitives/ints/lists/labelList.H"
#include "primitives/Vector/labelVector/labelVector.H"
#include "meshRefinement/meshRefinement.H"
#include "polyTopoChange/polyTopoChange/polyTopoChange.H"
#include "meshCut/edgeVertex/edgeVertex.H"
#include "meshes/polyMesh/mapPolyMesh/mapDistribute/mapDistributePolyMesh.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class anisoRefiner Declaration
\*---------------------------------------------------------------------------*/

class anisoRefiner
{

      // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Per cell the anisotropic cell refinement level
        List<labelVector> cellLevel_;

        //- Per cell the anisotropic cell refinement level
        List<labelVector> pointLevel_;

      // Private Member Functions

        //- Updates refineCell so consistent 2:1 refinement. Returns local
        //  number of cells changed.
        label faceConsistentRefinement
        (
            const direction cmpt,
            PackedBoolList& isProtectedCell,
            boolList& protectedFaces,
            PackedBoolList& refineCell
        ) const;

        //- Calculate edge tolerance for splitting
        scalar edgeTol(direction cmpt) const;

        void modifyProcessorBoundaries
        (
            const PackedBoolList& refineCell,
            const direction cmpt
        );

        //- Perform pre-split of processor faces
        void splitFace
        (
            const label faceI,
            const label own,
            const label patchI,
            const direction cmpt,
            const labelList& cutPoints,
            const labelList& newEdgePts,
            const labelList& newEdgeLevels,
            const edgeVertex& ev,
            const label cellLevel,
            polyTopoChange& meshMod
        ) const;

        //- Perform refinement of cells
        void refine
        (
            const labelList& cellsToRefine,
            const direction cmpt
        );

public:

    //- Runtime type information
    ClassName("anisoRefiner");

    //- Update local numbering for mesh redistribution
    void distribute(const mapDistributePolyMesh&);

    //Anistropic refinement of cells in cmmpt direction
    label setRefinement
    (
        const direction cmpt
    );

    // Constructors
    anisoRefiner
    (
        meshRefinement& meshRefiner,
        const labelList cellLevel,
        const labelList pointLevel
    );

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
