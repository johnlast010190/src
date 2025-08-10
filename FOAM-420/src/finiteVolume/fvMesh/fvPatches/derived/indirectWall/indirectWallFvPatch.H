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
    Foam::indirectWallFvPatch

Description
    Foam::indirectWallFvPatch

SourceFiles
    indirectWallFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef indirectWallFvPatch_H
#define indirectWallFvPatch_H

#include "fvMesh/fvPatches/derived/wall/wallFvPatch.H"
#include "meshes/polyMesh/polyPatches/derived/indirectWall/indirectWallPolyPatch.H"
#include "matrices/lduMatrix/lduAddressing/lduInterface/lduInterface.H"
#include "fvMesh/fvBoundaryMesh/fvBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class indirectWallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class indirectWallFvPatch
:
    public wallFvPatch
{

        const indirectWallPolyPatch& indirectPolyPatch_;

public:

    //- Runtime type information
    TypeName(indirectWallPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        indirectWallFvPatch(const polyPatch& patch, const fvBoundaryMesh& bm)
        :
            fvPatch(patch, bm),
            wallFvPatch(patch, bm),
            indirectPolyPatch_(refCast<const indirectWallPolyPatch>(patch))
        {}

    // Member Functions

        // Access

            label nbrPatchID() const
            {
                return indirectPolyPatch_.nbrPatchID();
            }

            const indirectWallFvPatch& neighbFvPatch() const
            {
                return refCast<const indirectWallFvPatch>
                (
                    this->boundaryMesh()[indirectPolyPatch_.nbrPatchID()]
                );
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
