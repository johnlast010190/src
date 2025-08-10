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

#include "nonConformal/polyPatches/nonConformalCoupled/nonConformalCoupledPolyPatch.H"
#include "nonConformal/polyPatches/nonConformalError/nonConformalErrorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalCoupledPolyPatch, 0);
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalCoupledPolyPatch::rename(const wordList& newNames)
{
    nonConformalPolyPatch::rename(newNames);

    errorPatchName_ = word::null;
    errorPatchID_ = -1;
}


void Foam::nonConformalCoupledPolyPatch::reorder
(
    const labelUList& oldToNewIndex
)
{
    nonConformalPolyPatch::reorder(oldToNewIndex);

    errorPatchName_ = word::null;
    errorPatchID_ = -1;
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalCoupledPolyPatch::nonConformalCoupledPolyPatch
(
    const polyPatch& patch
)
:
    nonConformalPolyPatch(patch),
    patch_(refCast<const coupledPolyPatch>(patch)),
    errorPatchName_(word::null),
    errorPatchID_(-1)
{}


Foam::nonConformalCoupledPolyPatch::nonConformalCoupledPolyPatch
(
    const polyPatch& patch,
    const word& origPatchName
)
:
    nonConformalPolyPatch(patch, origPatchName),
    patch_(refCast<const coupledPolyPatch>(patch)),
    errorPatchName_(word::null),
    errorPatchID_(-1)
{}


Foam::nonConformalCoupledPolyPatch::nonConformalCoupledPolyPatch
(
    const polyPatch& patch,
    const dictionary& dict
)
:
    nonConformalPolyPatch(patch, dict),
    patch_(refCast<const coupledPolyPatch>(patch)),
    errorPatchName_(word::null),
    errorPatchID_(-1)
{}


Foam::nonConformalCoupledPolyPatch::nonConformalCoupledPolyPatch
(
    const nonConformalCoupledPolyPatch& patch
)
:
    nonConformalPolyPatch(patch),
    patch_(refCast<const coupledPolyPatch>(patch)),
    errorPatchName_(word::null),
    errorPatchID_(-1)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalCoupledPolyPatch::~nonConformalCoupledPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::nonConformalCoupledPolyPatch::errorPatchID() const
{
    if (errorPatchID_ == -1)
    {
        forAll(patch_.boundaryMesh(), patchi)
        {
            const polyPatch& p = patch_.boundaryMesh()[patchi];

            if
            (
                isA<nonConformalErrorPolyPatch>(p)
             && refCast<const nonConformalErrorPolyPatch>(p).origPatchID()
             == origPatchID()
            )
            {
                errorPatchID_ = patchi;
                break;
            }
        }
    }

    if (errorPatchID_ == -1)
    {
        FatalErrorInFunction
            << "No error patch of type "
            << nonConformalErrorPolyPatch::typeName
            << " defined for patch " << origPatchName()
            << exit(FatalError);
    }

    return errorPatchID_;
}


const Foam::nonConformalErrorPolyPatch&
Foam::nonConformalCoupledPolyPatch::errorPatch() const
{
    return
        refCast<const nonConformalErrorPolyPatch>
        (
            patch_.boundaryMesh()[errorPatchID()]
        );
}


void Foam::nonConformalCoupledPolyPatch::write(Ostream& os) const
{
    nonConformalPolyPatch::write(os);
}


// ************************************************************************* //
