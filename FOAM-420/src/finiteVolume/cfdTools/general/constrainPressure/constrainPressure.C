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
    (c) 2016 OpenFOAM Foundation
    (c) 2016 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "cfdTools/general/constrainPressure/constrainPressure.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fields/GeometricFields/geometricOneField/geometricOneField.H"
#include "fields/fvPatchFields/derived/fixedFluxPressure/fixedFluxPressureFvPatchScalarField.H"
#include "cfdTools/general/fvOptions/fvOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RhoType, class RAUType, class FvOptionType>
void Foam::constrainPressure
(
    volScalarField& p,
    const RhoType& rho,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rhorAU,
    const FvOptionType& fvOptions,
    const bool includeRhoUsf
)
{
    const fvMesh& mesh = p.mesh();

    volScalarField::Boundary& pBf = p.boundaryFieldRef();

    const volVectorField::Boundary& UBf = U.boundaryField();
    const surfaceScalarField::Boundary& phiHbyABf =
        phiHbyA.boundaryField();
    const typename RAUType::Boundary& rhorAUBf =
        rhorAU.boundaryField();
    const surfaceVectorField::Boundary& SfBf =
        mesh.Sf().boundaryField();
    const surfaceScalarField::Boundary& magSfBf =
        mesh.magSf().boundaryField();

    forAll(pBf, patchi)
    {
        fixedFluxPressureFvPatchScalarField* ffp = nullptr;
        if (isA<fixedFluxPressureFvPatchScalarField>(pBf[patchi]))
        {
            ffp =
                &refCast<fixedFluxPressureFvPatchScalarField>
                (
                    pBf[patchi]
                );
        }
        else if (isA<blendedFvPatchField<scalar>>(pBf[patchi]))
        {
            blendedFvPatchField<scalar>& bpf =
                refCast<blendedFvPatchField<scalar>>(pBf[patchi]);
            if (isA<fixedFluxPressureFvPatchScalarField>(bpf.boundaryOne()))
            {
                ffp =
                    &refCast<fixedFluxPressureFvPatchScalarField>
                    (
                        bpf.boundaryOne()
                    );
            }
            else if
            (
                isA<fixedFluxPressureFvPatchScalarField>(bpf.boundaryTwo())
            )
            {
                ffp =
                    &refCast<fixedFluxPressureFvPatchScalarField>
                    (
                        bpf.boundaryTwo()
                    );
            }
        }
        if (ffp)
        {
            scalarField sphi(SfBf[patchi] & UBf[patchi]);

            ffp->updateSnGrad
            (
                (
                    includeRhoUsf
                    ?
                    (
                        phiHbyABf[patchi]
                      - rho.boundaryField()[patchi]*
                       (
                           fvOptions.relative(sphi, patchi)
                       )
                    )
                    :
                    (
                        phiHbyABf[patchi]
                    )
                )
               /(magSfBf[patchi]*rhorAUBf[patchi])
            );
        }
    }
}


template<class RAUType>
void Foam::constrainPressure
(
    volScalarField& p,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rAU
)
{
    constrainPressure(p, rho, U, phiHbyA, rAU, NullMRF(), true);
}


template<class RAUType, class FvOptionType>
void Foam::constrainPressure
(
    volScalarField& p,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rAU,
    const FvOptionType& fvOptions
)
{
    constrainPressure(p, geometricOneField(), U, phiHbyA, rAU, fvOptions);
}


template<class RAUType>
void Foam::constrainPressure
(
    volScalarField& p,
    const volVectorField& U,
    const surfaceScalarField& phiHbyA,
    const RAUType& rAU
)
{
    constrainPressure(p, U, phiHbyA, rAU, NullMRF());
}


// ************************************************************************* //
