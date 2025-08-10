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
    (c) 2011-2016 OpenFOAM Foundation
    (c) 2020 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "interpolation/surfaceInterpolation/limitedSchemes/BRICS/BRICS.H"
#include "finiteVolume/fvc/fvc.H"
#include "fields/volFields/volFields.H"
#include "fields/surfaceFields/surfaceFields.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "meshes/primitiveShapes/point/point.H"
#include <set>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::BRICS::limiter
(
    const scalar cdWeight,
    const scalar faceFlux,
    const scalar& phiP,
    const scalar& phiN,
    const vector& gradcP,
    const vector& gradcN,
    const scalar Cof,
    const vector d,
    const scalar blendingFactor
) const
{
    // Calculate BRICS weight
    scalar w =
        weight
        (
            faceFlux,
            phiP,
            phiN,
            gradcP,
            gradcN,
            Cof,
            d,
            blendingFactor
        );

    // Calculate BRICS face value
    scalar phif = w*phiP + (1 - w)*phiN;

    // Calculate UD and CD face value
    scalar phiU = faceFlux >= 0 ? phiP : phiN;
    scalar phiCD = cdWeight*phiP + (1 - cdWeight)*phiN;

    // Calculate the effective limiter for the QUICK interpolation
    scalar CLimiter = (phif - phiU)/stabilise(phiCD - phiU, SMALL);

    Info<< "Foam::scalar Foam::BRICS::limiter is called ...\n" << endl;
    // Limit the limiter between upwind and downwind
    return max(min(CLimiter, 2), 0);
}


Foam::tmp<Foam::surfaceScalarField> Foam::BRICS::limiter
(
    const volScalarField& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tLimiter
    (
        new surfaceScalarField
        (
            IOobject
            (
                type() + "Limiter(" + phi.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimless
        )
    );
    surfaceScalarField& lim = tLimiter.ref();

    volVectorField gradc( fvc::grad(phi) );

    const surfaceScalarField faceFlux = this->faceFlux_;

    volScalarField outFlux( fvc::surfaceSum(mag(faceFlux)) );
    outFlux.primitiveFieldRef() /= mesh.V();

    // this calculation is only valid for a divergence-free flux field!!!
    // otherwise: need to sum up all fluxes leaving the centre cell
    surfaceScalarField Cof
    (
        0.5 * mesh.time().deltaT()
       *upwind<scalar>(mesh, faceFlux).interpolate
        (
            outFlux
        )
    );

    calculateTheta();
    calculateBlendingFactor();

    const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const vectorField& C = mesh.C();

    scalarField& pLim = lim.primitiveFieldRef();

    forAll(pLim, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector d = C[nei] - C[own];

        pLim[faceI] = limiter
        (
            CDweights[faceI],
            this->faceFlux_[faceI],
            phi[own],
            phi[nei],
            gradc[own],
            gradc[nei],
            Cof[faceI],
            d,
            blendingFactor_[faceI]
        );
    }

    surfaceScalarField::Boundary& bLim = lim.boundaryFieldRef();
    forAll(bLim, patchi)
    {
        scalarField& pLim = bLim[patchi];

        if (bLim[patchi].coupled())
        {
            const scalarField& pCDweights = CDweights.boundaryField()[patchi];

            const scalarField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchi];

            scalarField pphiP
            (
                phi.boundaryField()[patchi].patchInternalField()
            );

            scalarField pphiN
            (
                phi.boundaryField()[patchi].patchNeighbourField()
            );

            vectorField pGradcP
            (
                gradc.boundaryField()[patchi].patchInternalField()
            );

            vectorField pGradcN
            (
                gradc.boundaryField()[patchi].patchNeighbourField()
            );

            const scalarField& pCof = Cof.boundaryField()[patchi];
            const scalarField& pblendingFactor = blendingFactor_.boundaryField()[patchi];

            // Build the d-vectors
            // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
            vectorField pd( bLim[patchi].patch().delta() );

            forAll(pLim, faceI)
            {
                pLim[faceI] = limiter
                (
                    pCDweights[faceI],
                    pFaceFlux[faceI],
                    pphiP[faceI],
                    pphiN[faceI],
                    pGradcP[faceI],
                    pGradcN[faceI],
                    pCof[faceI],
                    pd[faceI],
                    pblendingFactor[faceI]
                );
            }
        }
        else
        {
            pLim = 1.0;
        }
    }

    return tLimiter;
}


Foam::scalar Foam::BRICS::weight
(
    const scalar faceFlux,
    const scalar& phiP,
    const scalar& phiN,
    const vector& gradcP,
    const vector& gradcN,
    const scalar Cof,
    const vector d,
    const scalar blendingFactor
) const
{
    // Calculate upwind value, faceFlux C tilde and do a stabilisation
    scalar phict = 0;
    scalar phiupw = 0;

    if (faceFlux > 0)
    {
        phiupw = phiN - 2*(gradcP & d);
        phiupw = max(min(phiupw, 1.0), 0.0);

        if ((phiN - phiupw) > 0)
        {
            phict = (phiP - phiupw)/(phiN - phiupw + SMALL);
        }
        else
        {
            phict = (phiP - phiupw)/(phiN - phiupw - SMALL);
        }
    }
    else
    {
        phiupw = phiP + 2*(gradcN & d);
        phiupw = max(min(phiupw, 1.0), 0.0);

        if ((phiP - phiupw) > 0)
        {
            phict = (phiN - phiupw)/(phiP - phiupw + SMALL);
        }
        else
        {
            phict = (phiN - phiupw)/(phiP - phiupw - SMALL);
        }
    }

    // BRICS weights
    scalar pCo = 0;
    scalar alphapCo = 1;
    scalar betaBICS = 0;
    scalar phifupwt = 0;
    scalar phifGDS = 0;
    scalar phifBICSstar = 0;
    scalar phifBRICS = 0;
    scalar weight;

    if (Cof > 0.3)
    {
        alphapCo = (Cof-0.3)/(Foam::exp(Cof-0.3) - 1.);
    }
    pCo = alphapCo * p_IGDS_ + (1.-alphapCo) * p_GDS_;
    betaBICS = beta_IGDS_ - (beta_GDS_ - beta_IGDS_)/(p_GDS_ - p_IGDS_)*p_IGDS_
        + (beta_GDS_ - beta_IGDS_)/(p_GDS_ - p_IGDS_)*pCo;

    // Calculate UDS-Scheme
    phifupwt = phict;

    // Calculate GAMMA-Scheme
    if (phict >= 0 && phict < beta_GDS_)
    {
        phifGDS = - Foam::pow(phict, 2.)/(2.*beta_GDS_) + (1.+1./(2.*beta_GDS_))*phict;
    }
    else if (phict >= beta_GDS_ && phict <= 1)
    {
        phifGDS = 0.5*phict + 0.5;
    }
    else
    {
        phifGDS = phict;
    }

    // Calculate BICS-Star-Scheme
    if (phict >= 0 && phict < betaBICS)
    {
        phifBICSstar = - Foam::pow(phict, 2.)*(1.-pCo)/Foam::pow(betaBICS, 2.) + (pCo+2.*(1.-pCo)/betaBICS)*phict;
    }
    else if (phict >= betaBICS && phict <= 1)
    {
        phifBICSstar = pCo*phict + (1.-pCo);
    }
    else
    {
        phifBICSstar = phict;
    }

    // Calculate the weighting factors for BRICS
    if (phict >= 1 || phict <= 0)  // use upwind
    {
        phifBRICS = phifupwt;
        weight = 0;
    }
    else
    {
        phifBRICS = phifBICSstar * blendingFactor + phifGDS * (1-blendingFactor);
        weight = (phifBRICS - phict)/(1 - phict);
    }

    if (faceFlux > 0)
    {
        return 1 - weight;
    }
    else
    {
        return weight;
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::BRICS::calcWeights
(
    const volScalarField& phi
) const
{
    const fvMesh& mesh = this->mesh();

    tmp<surfaceScalarField> tWeightingFactors
    (
        new surfaceScalarField(mesh.surfaceInterpolation::weights())
    );
    surfaceScalarField& weightingFactors = tWeightingFactors.ref();

    volVectorField gradc( fvc::grad(phi) );

    const surfaceScalarField faceFlux = this->faceFlux_;

    volScalarField outFlux( fvc::surfaceSum(mag(faceFlux)) );
    outFlux.primitiveFieldRef() /= mesh.V();

    //- fix mass conservation problems with coupled boundaries (parallel runs)
    //  i.e. correct calculation of local Courant number at coupled boundaries
    if (Pstream::parRun())
    {
        outFlux.correctBoundaryConditions();
    }

    surfaceScalarField Cof
    (
        0.5 * mesh.time().deltaT()
       *upwind<scalar>(mesh, faceFlux).interpolate
        (
            outFlux
        )
    );

    calculateTheta();
    calculateBlendingFactor();

    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();

    const vectorField& C = mesh.C();

    scalarField& w = weightingFactors.primitiveFieldRef();

    forAll(w, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector d = C[nei] - C[own];

        w[faceI] = weight
        (
            this->faceFlux_[faceI],
            phi[own],
            phi[nei],
            gradc[own],
            gradc[nei],
            Cof[faceI],
            d,
            blendingFactor_[faceI]
        );
    }

    surfaceScalarField::Boundary& bWeights =
        weightingFactors.boundaryFieldRef();

    forAll(bWeights, patchi)
    {
        scalarField& pWeights = bWeights[patchi];

        if (bWeights[patchi].coupled())
        {
            const scalarField& pFaceFlux =
                this->faceFlux_.boundaryField()[patchi];

            scalarField pphiP
            (
                phi.boundaryField()[patchi].patchInternalField()
            );

            scalarField pphiN
            (
                phi.boundaryField()[patchi].patchNeighbourField()
            );

            vectorField pGradcP
            (
                gradc.boundaryField()[patchi].patchInternalField()
            );

            vectorField pGradcN
            (
                gradc.boundaryField()[patchi].patchNeighbourField()
            );

            const scalarField& pCof = Cof.boundaryField()[patchi];
            const scalarField& pblendingFactor = blendingFactor_.boundaryField()[patchi];

            // Build the d-vectors
            // Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
            vectorField pd( bWeights[patchi].patch().delta() );

            forAll(pWeights, faceI)
            {
                pWeights[faceI] = weight
                (
                    pFaceFlux[faceI],
                    pphiP[faceI],
                    pphiN[faceI],
                    pGradcP[faceI],
                    pGradcN[faceI],
                    pCof[faceI],
                    pd[faceI],
                    pblendingFactor[faceI]
                );
            }
        }
        else
        {
            pWeights = 1.0;
        }
    }

    return tWeightingFactors;
}


Foam::tmp<Foam::surfaceScalarField> Foam::BRICS::correction
(
    const volScalarField& vf
) const
{
    const fvMesh& mesh = this->mesh();
    const surfaceScalarField faceFlux = this->faceFlux_;

    // Correction = full interpolation minus upwinded part
    tmp<surfaceScalarField> tcorr
    (
        new surfaceScalarField
        (
            IOobject
            (
                type() + "correction(" + vf.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            vf.dimensions()
        )
    );

    surfaceScalarField& corr = tcorr.ref();

    upwind<scalar> upwindScheme(mesh, faceFlux);

    corr = surfaceInterpolationScheme<scalar>::interpolate(vf, calcWeights(vf)) - upwindScheme.interpolate(vf);

    return tcorr;
}


void Foam::BRICS::calculateTheta() const
{
    // Blending only in narrow band around interface
    const fvMesh& mesh = this->mesh();

    const volScalarField& vsf = mesh.objectRegistry::lookupObject<const volScalarField> (blendingField_);
    const volVectorField vsfGrad(fvc::grad(vsf,"alphaS"));
    blendingFactor_ =
    (
        linearInterpolate(
            vsfGrad/(mag(vsfGrad)+deltaN_)
        )
        &
        mesh.Sf()/mesh.magSf()
    );

    //surfaceScalarField ssf = linearInterpolate(vsf);
    //blendingFactor_ *= pos(ssf-0.001)*pos(0.999-ssf);
}


void Foam::BRICS::calculateBlendingFactor() const
{
    // limit blending factor
    blendingFactor_ = Foam::mag(blendingFactor_);
    blendingFactor_ = max(min(blendingFactor_, scalar(1.)), scalar(0.));

    // Compute Courant number
    const fvMesh& mesh = this->mesh();
    const surfaceScalarField faceFlux = this->faceFlux_;

    volScalarField outFlux( fvc::surfaceSum(mag(faceFlux)) );
    outFlux.primitiveFieldRef() /= mesh.V();
    if (Pstream::parRun())
    {
        outFlux.correctBoundaryConditions();
    }
    surfaceScalarField Cof
    (
        0.5 * mesh.time().deltaT()
       *upwind<scalar>(mesh, faceFlux).interpolate
        (
            outFlux
        )*dimensionedScalar("fac", dimless/dimVolume, 1.)
    );

    // Blending factor for BRICS:
    blendingFactor_ = Foam::pow(blendingFactor_, Cof);
}


namespace Foam
{
defineTypeNameAndDebug(BRICS, 0);

surfaceInterpolationScheme<scalar>::addMeshConstructorToTable<BRICS>
    addBRICSshConstructorToTable_;

surfaceInterpolationScheme<scalar>::addMeshFluxConstructorToTable<BRICS>
    addBRICSMeshFluxConstructorToTable_;

limitedSurfaceInterpolationScheme<scalar>::addMeshConstructorToTable<BRICS>
    addBRICSMeshConstructorToLimitedTable_;

limitedSurfaceInterpolationScheme<scalar>::
addMeshFluxConstructorToTable<BRICS>
    addBRICSMeshFluxConstructorToLimitedTable_;
}

// ************************************************************************* //
