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

#include "fvMesh/fvPatches/constraint/nonConformalProcessorCyclic/nonConformalProcessorCyclicFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalProcessorCyclicFvPatch, 0);
    addToRunTimeSelectionTable
    (
        fvPatch,
        nonConformalProcessorCyclicFvPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::nonConformalProcessorCyclicFvPatch::makeWeights
(
    scalarField& w
) const
{
    if (Pstream::parRun())
    {
        nonConformalCoupledFvPatch::makeWeights
        (
            w,
          - boundaryMesh().mesh().Sf().boundaryField()[index()],
            boundaryMesh().mesh().Cf().boundaryField()[index()]
          - boundaryMesh().mesh().C().boundaryField()[index()]
        );
    }
    else
    {
        w = 1;
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalProcessorCyclicFvPatch::nonConformalProcessorCyclicFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    processorCyclicFvPatch(patch, bm),
    nonConformalCoupledFvPatch(static_cast<const fvPatch&>(*this)),
    nonConformalProcessorCyclicPolyPatch_
    (
        refCast<const nonConformalProcessorCyclicPolyPatch>(patch)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalProcessorCyclicFvPatch::~nonConformalProcessorCyclicFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField>
Foam::nonConformalProcessorCyclicFvPatch::delta() const
{
    if (Pstream::parRun())
    {
        return
            coupledFvPatch::delta
            (
                boundaryMesh().mesh().Cf().boundaryField()[index()]
              - boundaryMesh().mesh().C().boundaryField()[index()]
            );
    }
    else
    {
        return coupledFvPatch::delta();
    }
}


// ************************************************************************* //
