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
    (c) 2023 Esi Ltd.

Class
    Foam::nonConformalOrigFvPatch

Description
    Original uncoupled non-conformal finite-volume patch, from which
    non-conformal cyclic patches will be constructed.

SourceFiles
    nonConformalOrigFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalOrigFvPatch_H
#define nonConformalOrigFvPatch_H

#include "fvMesh/fvPatches/fvPatch/fvPatch.H"
#include "meshes/polyMesh/polyPatches/derived/nonConformalOrig/nonConformalOrigPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class nonConformalOrigFvPatch Declaration
\*---------------------------------------------------------------------------*/

class nonConformalOrigFvPatch
:
    public fvPatch
{
public:

    //- Runtime type information
    TypeName(nonConformalOrigPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        nonConformalOrigFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            fvPatch(patch, bm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
