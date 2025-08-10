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

#include "fvMesh/fvPatches/constraint/nonConformal/nonConformalFvPatch.H"
#include "fields/surfaceFields/surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nonConformalFvPatch, 0);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::nonConformalFvPatch::nonConformalFvPatch
(
    const fvPatch& patch
)
:
    patch_(patch),
    nonConformalPolyPatch_
    (
        refCast<const nonConformalPolyPatch>(patch.patch())
    ),
    faceCells_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nonConformalFvPatch::~nonConformalFvPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::nonConformalFvPatch::polyFaces() const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    return
        mesh.conformal()
      ? labelList::null()
      : mesh.polyFacesBf()[patch_.index()];
}


Foam::label Foam::nonConformalFvPatch::start() const
{
    if (size())
    {
        FatalErrorInFunction
            << "The start face is not defined for a " << typeName
            << " patch with a non-zero number of faces"
            << exit(FatalError);
    }

    return patch_.patch().start();
}


const Foam::labelUList& Foam::nonConformalFvPatch::faceCells() const
{
    // !!! This needs an update mechanism, rather than re-calculating the
    // face-cells every time

    const fvMesh& mesh = patch_.boundaryMesh().mesh();

    faceCells_ = UIndirectList<label>(mesh.faceOwner(), polyFaces());

    return faceCells_;
}


// ************************************************************************* //
