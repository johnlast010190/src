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
    Foam::nonConformalFuncs

Description
    Collection of functions used in the creation of non-conformal couples
    between non-coupled patches.

SourceFiles
    nonConformalFuncs.C

\*---------------------------------------------------------------------------*/

#ifndef nonConformalFuncs_H
#define nonConformalFuncs_H

#include "meshes/polyMesh/polyPatches/constraint/cyclic/cyclicTransform.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyPatch;
class fvMesh;


namespace nonConformalFuncs
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Helper function
void generateCouples
(
    fvMesh& mesh,
    const List<Pair<word>>& patchNames,
    const wordList& cyclicNames,
    const List<cyclicTransform>& transforms
);

//- Main driver function, to create the non-conformal couples and insert them
//  into the mesh
void createNonConformalCouples(fvMesh& mesh);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nonConformalFuncs
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
