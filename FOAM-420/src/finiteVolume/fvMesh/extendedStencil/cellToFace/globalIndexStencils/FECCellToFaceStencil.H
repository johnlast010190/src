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
    (c) 2011-2016 OpenFOAM Foundation

Class
    Foam::FECCellToFaceStencil

Description
    All cells connected via edge to face.

SourceFiles
    FECCellToFaceStencil.C

\*---------------------------------------------------------------------------*/

#ifndef FECCellToFaceStencil_H
#define FECCellToFaceStencil_H

#include "fvMesh/extendedStencil/cellToFace/globalIndexStencils/cellToFaceStencil.H"
#include "meshes/meshShapes/edge/EdgeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class FECCellToFaceStencil Declaration
\*---------------------------------------------------------------------------*/

class FECCellToFaceStencil
:
    public cellToFaceStencil
{
    // Private Member Functions

        //- Calculates per edge the neighbour data (= edgeCells)
        void calcEdgeBoundaryData
        (
            const boolList& isValidBFace,
            const labelList& boundaryEdges,
            EdgeMap<labelList>& neiGlobal
        ) const;

        //- Calculates per face the edge connected data
        //  (= cell or boundary in global numbering).
        void calcFaceStencil(labelListList& faceStencil) const;


public:

    // Constructors

        //- Construct from all cells and boundary faces
        explicit FECCellToFaceStencil(const polyMesh&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
