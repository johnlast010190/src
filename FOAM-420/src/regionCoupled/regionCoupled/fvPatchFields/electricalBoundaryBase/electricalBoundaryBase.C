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
    (c) 2019 Esi Ltd.
    (c) 2011-2016 OpenFOAM Foundation

\*---------------------------------------------------------------------------*/

#include "electricalBoundaryBase.H"
#include "fields/volFields/volFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::electricalBoundaryBase::electricalBoundaryBase
(
    const fvPatch& patch
)
:
    patch_(patch)
{}


Foam::electricalBoundaryBase::electricalBoundaryBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    electricalBoundaryBase(patch)
{}


Foam::electricalBoundaryBase::electricalBoundaryBase
(
    const fvPatch& patch,
    const electricalBoundaryBase& base
)
:
    electricalBoundaryBase(patch)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::electricalBoundaryBase::anisotropic() const
{
    return
        patch_.boundaryMesh().mesh().foundObject<volSymmTensorField>
        (
            "electrical_sigma"
        );
}

Foam::tmp<Foam::vectorField> Foam::electricalBoundaryBase::nfSigma() const
{
    if (anisotropic())
    {
        const symmTensorField& sigmaWall =
            patch_.lookupPatchField<volSymmTensorField, scalar>
            (
                "electrical_sigma"
            );

        return patch_.nf() & sigmaWall;
    }
    else
    {
        return
            patch_.nf()*patch_.lookupPatchField<volScalarField, scalar>
            (
                "electrical_sigma"
            );
    }
}


void Foam::electricalBoundaryBase::write(Ostream& os) const
{}


// ************************************************************************* //
