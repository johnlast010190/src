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
    (c) 2021-2022 OpenFOAM Foundation
    (c) 2022 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMesh/fvPatches/constraint/nonConformalCyclic/nonConformalCyclicFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCyclicFvPatch, 0);
    addToRunTimeSelectionTable(fvPatch, nonConformalCyclicFvPatch, polyPatch);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::nonConformalCyclicFvPatch::makeWeights(scalarField& w) const
{
    nonConformalCoupledFvPatch::makeWeights
    (
        w,
        nbrPatch().Sf(),
        nbrPatch().coupledFvPatch::delta()
    );
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicFvPatch::nonConformalCyclicFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    cyclicFvPatch(patch, bm),
    nonConformalCoupledFvPatch(static_cast<const fvPatch&>(*this)),
    nonConformalCyclicPolyPatch_
    (
        refCast<const nonConformalCyclicPolyPatch>(patch)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCyclicFvPatch::~nonConformalCyclicFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nonConformalCyclicFvPatch&
Foam::nonConformalCyclicFvPatch::nbrPatch() const
{
    return refCast<const nonConformalCyclicFvPatch>(cyclicFvPatch::nbrPatch());
}


// ************************************************************************* //
