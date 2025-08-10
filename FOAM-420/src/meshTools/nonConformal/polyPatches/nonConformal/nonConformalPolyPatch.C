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
    (c) 2022-2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "nonConformal/polyPatches/nonConformal/nonConformalPolyPatch.H"
#include "AMIInterpolation/patches/cyclicAMI/cyclicAMIPolyPatch/cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalPolyPatch, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void Foam::nonConformalPolyPatch::validateSize() const
{
    if (patch_.size() != 0)
    {
        FatalErrorInFunction
            << "Patch " << patch_.name() << " has " << patch_.size()
            << " faces. Patches of type " << patch_.type()
            << " must have zero faces." << exit(FatalError);
    }
}


// * * * * * * * * * * * *  Protected Member Functions * * * * * * * * * * * //

void Foam::nonConformalPolyPatch::rename(const wordList& newNames)
{
    if (origPatchID_ != -1)
    {
        origPatchName_ = newNames[origPatchID_];
    }
    else
    {
        FatalErrorInFunction
            << "Cannot rename " << nonConformalPolyPatch::typeName
            << " without the original patch index"
            << exit(FatalError);
    }
}


void Foam::nonConformalPolyPatch::reorder(const labelUList& oldToNewIndex)
{
    if (origPatchID_ != -1)
    {
        origPatchID_ = oldToNewIndex[origPatchID_];
    }
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalPolyPatch::nonConformalPolyPatch(const polyPatch& patch)
:
    patch_(patch),
    origPatchName_(word::null),
    origPatchID_(-1)
{
    validateSize();
}


Foam::nonConformalPolyPatch::nonConformalPolyPatch
(
    const polyPatch& patch,
    const word& origPatchName
)
:
    patch_(patch),
    origPatchName_(origPatchName),
    origPatchID_(-1)
{
    validateSize();
}


Foam::nonConformalPolyPatch::nonConformalPolyPatch
(
    const polyPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    origPatchName_(dict.lookup("originalPatch")),
    origPatchID_(-1)
{
    validateSize();
}


Foam::nonConformalPolyPatch::nonConformalPolyPatch
(
    const nonConformalPolyPatch& patch
)
:
    patch_(patch.patch_),
    origPatchName_(patch.origPatchName_),
    origPatchID_(-1)
{
    validateSize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalPolyPatch::~nonConformalPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::nonConformalPolyPatch::origPatchID() const
{
    if (origPatchID_ == -1)
    {
        origPatchID_ = patch_.boundaryMesh().findPatchID(origPatchName());

        if (origPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Invalid origPatch name " << origPatchName()
                << endl << "Valid patch names are "
                << patch_.boundaryMesh().names()
                << exit(FatalError);
        }

        const polyPatch& p = patch_.boundaryMesh()[origPatchID_];

        if (p.coupled() || isA<cyclicAMIPolyPatch>(p))
        {
            FatalErrorInFunction
                << "The original patch for the " << patch_.type()
                << " patch " << patch_.name() << " is " << p.name()
                << ", which is of " << p.type() << " type." << nl
                << "This is not allowed. The original patch cannot be a coupled"
                << " or AMI patch." << nl << "Please, change the patch type to"
                << "'nonConformal' to create a non-conformal coupling,"
                << " or remove the 'nonConformalCouplesDict' to proceed with"
                << " cyclicAMI patches."
                << exit(FatalError);
        }
    }

    return origPatchID_;
}


void Foam::nonConformalPolyPatch::write(Ostream& os) const
{
    os.writeEntry("originalPatch", origPatchName_);
}


// ************************************************************************* //
