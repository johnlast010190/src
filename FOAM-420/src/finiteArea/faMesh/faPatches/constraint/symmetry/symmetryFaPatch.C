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
    (c) 2016 Esi Ltd.


Description

\*---------------------------------------------------------------------------*/

#include "faMesh/faPatches/constraint/symmetry/symmetryFaPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(symmetryFaPatch, 0);
addToRunTimeSelectionTable(faPatch, symmetryFaPatch, dictionary);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void symmetryFaPatch::makeCorrVecs(vectorField& cv) const
{
    // Non-orthogonal correction not allowed.  HJ, 16/Apr/2009
    cv = vector::zero;
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

// Construct from components
symmetryFaPatch::symmetryFaPatch
(
    const word& name,
    const labelList& edgeLabels,
    const label index,
    const faBoundaryMesh& bm,
    const label ngbPolyPatchIndex
)
:
    faPatch(name, edgeLabels, index, bm, ngbPolyPatchIndex)
{}

//- Construct from dictionary
symmetryFaPatch::symmetryFaPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faBoundaryMesh& bm
)
:
    faPatch(name, dict, index, bm)
{
    if (ngbPolyPatchIndex() == -1)
    {
        FatalErrorInFunction
            << "Neighbour polyPatch index is not specified for faPatch "
            << this->name()
            << exit(FatalError);
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
