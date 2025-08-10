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
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "MRFSource.H"
#include "fvMesh/fvMesh.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "fvMatrices/fvMatrices.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::fv::MRFSource::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }

    const labelList& internalFaces = frameSourceFaces_.internalFaces();
    forAll(internalFaces, i)
    {
        const label facei = internalFaces[i];
        phi[facei] -= rho[facei]*getRotationalFlux(facei);
    }
    makeRelativeRhoFlux(rho.boundaryField(), phi.boundaryFieldRef());
}


template<class RhoFieldType>
void Foam::fv::MRFSource::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    FieldField<fvsPatchField, scalar>& phi
) const
{
    if (!active_)
    {
        return;
    }

    const volVectorField::Boundary& UBf =
        obr_.lookupObject<volVectorField>(UName_).boundaryField();
    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& includedFaces =
            frameSourceFaces_.includedFaces()[patchi];
        if
        (
            isA<refFvPatch>(UBf[patchi])
         && dynamic_cast<const refFvPatch&>(UBf[patchi]).inletFlux()
        )
        {
            forAll(includedFaces, i)
            {
                const label patchFacei = includedFaces[i];
                phi[patchi][patchFacei] -=
                    rho[patchi][patchFacei]
                   *getBoundaryRotationalFlux(patchi, patchFacei);
            }
        }
        else
        {
            UIndirectList<scalar>(phi[patchi], includedFaces) = 0.0;
        }
    }

    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        forAll(excludedFaces, i)
        {
            label patchFacei = excludedFaces[i];
            phi[patchi][patchFacei] -=
                rho[patchi][patchFacei]
               *getBoundaryRotationalFlux(patchi, patchFacei);
        }
    }
}


template<class RhoFieldType>
void Foam::fv::MRFSource::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    Field<scalar>& phi,
    const label patchi
) const
{
    if (!active_)
    {
        return;
    }

    const labelList& includedFaces = frameSourceFaces_.includedFaces()[patchi];
    UIndirectList<scalar>(phi, includedFaces) = 0.0;

    const labelList& excludedFaces = frameSourceFaces_.excludedFaces()[patchi];
    forAll(excludedFaces, i)
    {
        const label patchFacei = excludedFaces[i];
        phi[patchFacei] -=
            rho[patchFacei]*getBoundaryRotationalFlux(patchi, patchFacei);
    }
}


template<class RhoFieldType>
void Foam::fv::MRFSource::makeAbsoluteRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }

    forAll(frameSourceFaces_.internalFaces(), i)
    {
        const label facei = frameSourceFaces_.internalFaces()[i];
        phi[facei] += rho[facei]*getRotationalFlux(facei);
    }

    forAll(frameSourceFaces_.includedFaces(), patchi)
    {
        const labelList& patchFaces =
            frameSourceFaces_.includedFaces()[patchi];
        forAll(patchFaces, i)
        {
            const label patchFacei = patchFaces[i];
            scalarField& pphi = phi.boundaryFieldRef()[patchi];
            pphi[patchFacei] +=
                rho.boundaryField()[patchi][patchFacei]
               *getBoundaryRotationalFlux(patchi, patchFacei);
        }
    }

    forAll(frameSourceFaces_.excludedFaces(), patchi)
    {
        const labelList& excludedFaces =
            frameSourceFaces_.excludedFaces()[patchi];
        forAll(excludedFaces, i)
        {
            const label patchFacei = excludedFaces[i];
            scalarField& pphi = phi.boundaryFieldRef()[patchi];
            pphi[patchFacei] +=
                rho.boundaryField()[patchi][patchFacei]
               *getBoundaryRotationalFlux(patchi, patchFacei);
        }
    }
}


// ************************************************************************* //
