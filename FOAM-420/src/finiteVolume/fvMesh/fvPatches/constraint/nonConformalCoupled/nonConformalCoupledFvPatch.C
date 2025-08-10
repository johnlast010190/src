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

#include "fvMesh/fvPatches/constraint/nonConformalCoupled/nonConformalCoupledFvPatch.H"
#include "fvMesh/fvPatches/constraint/nonConformalError/nonConformalErrorFvPatch.H"
#include "db/runTimeSelection/construction/addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCoupledFvPatch, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::nonConformalCoupledFvPatch::makeWeights
(
    scalarField& w,
    const vectorField& nbrSf,
    const vectorField& nbrDelta
) const
{
    const vectorField delta(patch_.coupledFvPatch::delta());

    const scalarField nfDelta(patch_.nf() & delta);

    const scalarField nbrNfDelta((nbrSf/(mag(nbrSf) + VSMALL)) & nbrDelta);

    forAll(delta, facei)
    {
        const scalar ndoi = nfDelta[facei];
        const scalar ndni = nbrNfDelta[facei];
        const scalar ndi = ndoi + ndni;

        // !!! Note this is a different form of stabilisation compared to
        // coupledFvPatch. This is necessary to prevent negative weights on
        // very small faces. Should this also be the standard form of
        // stabilisation in coupledFvPatch and in surfaceInterpolation?

        if (ndoi > VSMALL && ndni > VSMALL)
        {
            w[facei] = ndni/ndi;
        }
        else
        {
            const scalar doi = mag(delta[facei]);
            const scalar dni = mag(nbrDelta[facei]);
            const scalar di = doi + dni;

            w[facei] = dni/di;
        }
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalCoupledFvPatch::nonConformalCoupledFvPatch
(
    const fvPatch& patch
)
:
    nonConformalFvPatch(patch),
    patch_(refCast<const coupledFvPatch>(patch)),
    nonConformalCoupledPolyPatch_
    (
        refCast<const nonConformalCoupledPolyPatch>(patch.patch())
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCoupledFvPatch::~nonConformalCoupledFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nonConformalErrorFvPatch&
Foam::nonConformalCoupledFvPatch::errorPatch() const
{
    return
        refCast<const nonConformalErrorFvPatch>
        (
            patch_.boundaryMesh()[errorPatchID()]
        );
}


// ************************************************************************* //
