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
    (c) 1991-2005 OpenCFD Ltd.

Description

\*---------------------------------------------------------------------------*/

#include "faMesh/faPatches/constraint/empty/emptyFaPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Patch name
defineTypeNameAndDebug(emptyFaPatch, 0);

// Add the patch constructor functions to the hash tables
addToRunTimeSelectionTable(faPatch, emptyFaPatch, dictionary);

// Over-riding the face normals return from the underlying patch
// This is the only piece of info used out of the underlying primitivePatch
// I choose to store it there because it is used in primitive patch operations
// and it should not be duplicated as before.  However, to ensure everything
// in the empty patch is sized to zero, we shall here return a regerence to
// a zero-sized field (it does not matter what the field is
//
// const vectorField& emptyFaPatch::edgeNormals() const
// {
//     return faceAreas();
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
