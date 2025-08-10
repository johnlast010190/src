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
    (c) 2011-2017 OpenFOAM Foundation
    (c) 2023 Esi Ltd.

\*---------------------------------------------------------------------------*/

#include "fvMatrices/solvers/MULES/MULES.H"
#include "interpolation/surfaceInterpolation/limitedSchemes/upwind/upwind.H"
#include "finiteVolume/fvc/fvcSurfaceIntegrate.H"
#include "finiteVolume/ddtSchemes/localEulerDdtScheme/localEulerDdtScheme.H"
#include "fields/surfaceFields/slicedSurfaceFields.H"
#include "fvMesh/fvPatches/constraint/wedge/wedgeFvPatch.H"
#include "finiteVolume/fvc/fvcCellFaceFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su
)
{
    Info<< "MULES: Solving for " << psi.name() << endl;

    const fvMesh& mesh = psi.mesh();

    scalarField& psiIf = psi;
    const scalarField& psi0 = psi.oldTime();

    psiIf = 0.0;
    psiIf = fvc::surfaceIntegrate(phiPsi).ref();

    if (mesh.moving())
    {
        psiIf =
        (
            mesh.Vsc0()().field()*rho.oldTime().field()
           *psi0*rDeltaT/mesh.Vsc()().field()
          + Su.field()
          - psiIf
        )/(rho.field()*rDeltaT - Sp.field());
    }
    else
    {
        psiIf =
        (
            rho.oldTime().field()*psi0*rDeltaT
          + Su.field()
          - psiIf
        )/(rho.field()*rDeltaT - Sp.field());
    }

    psi.correctBoundaryConditions();
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su
)
{
    const fvMesh& mesh = psi.mesh();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);
        explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();
        explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
    }
}


template<class RhoType, class SpType, class SuType>
void Foam::MULES::explicitSolve
(
    const RhoType& rho,
    volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin
)
{
    const fvMesh& mesh = psi.mesh();

    psi.correctBoundaryConditions();

    if (fv::localEulerDdt::enabled(mesh))
    {
        const volScalarField& rDeltaT = fv::localEulerDdt::localRDeltaT(mesh);

        limit
        (
            rDeltaT,
            rho,
            psi,
            phi,
            phiPsi,
            Sp,
            Su,
            psiMax,
            psiMin,
            false
        );

        explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
    }
    else
    {
        const scalar rDeltaT = 1.0/mesh.time().deltaTValue();

        limit
        (
            rDeltaT,
            rho,
            psi,
            phi,
            phiPsi,
            Sp,
            Su,
            psiMax,
            psiMin,
            false
        );

        explicitSolve(rDeltaT, rho, psi, phiPsi, Sp, Su);
    }
}


template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::limiter
(
    surfaceScalarField& lambda,
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phiBD,
    const surfaceScalarField& phiCorr,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin
)
{
    const scalarField& psiIf = psi;
    const volScalarField::Boundary& psiBf = psi.boundaryField();

    const fvMesh& mesh = psi.mesh();

    const dictionary& MULEScontrols = mesh.solution().solverDict(psi.name());

    const label nLimiterIter
    (
        MULEScontrols.lookupOrDefault<label>("nLimiterIter", 3)
    );

    const scalar smoothLimiter
    (
        MULEScontrols.lookupOrDefault<scalar>("smoothLimiter", 0)
    );

    const scalar extremaCoeff
    (
        MULEScontrols.lookupOrDefault<scalar>("extremaCoeff", 0)
    );

    const scalarField& psi0 = psi.oldTime();

    const labelUList& owner = mesh.owner();
    const labelUList& neighb = mesh.neighbour();
    tmp<volScalarField::Internal> tVsc = mesh.Vsc();
    const scalarField& V = tVsc();

    const scalarField& phiCorrIf = phiCorr;
    const surfaceScalarField::Boundary& phiCorrBf =
        phiCorr.boundaryField();

    scalarField& lambdaIf = lambda;
    surfaceScalarField::Boundary& lambdaBf = lambda.boundaryFieldRef();

    scalarField psiMaxn(psiIf.size(), psiMin);
    scalarField psiMinn(psiIf.size(), psiMax);

    scalarField sumPhiBD(psiIf.size(), 0.0);

    scalarField sumPhip(psiIf.size(), 0.0);
    scalarField mSumPhim(psiIf.size(), 0.0);

    surfaceScalarField mpsiOwn(fvc::ownerField(psi));
    surfaceScalarField mpsiNei(fvc::neighbourField(psi));
    fvc::applyFaceMaskTo(mpsiOwn);
    fvc::applyFaceMaskTo(mpsiNei);
    // Add and sync internal and boundary values for min/max
    fvc::sumIndirectAndInternalFaces(mpsiOwn, mpsiNei);
    tmp<surfaceScalarField> tmphiBD = fvc::applyFaceMask(phiBD);
    const surfaceScalarField& mphiBD = tmphiBD();
    const surfaceScalarField::Boundary& mphiBDBf = mphiBD.boundaryField();
    tmp<surfaceScalarField> tmphiCorr = fvc::applyFaceMask(phiCorr);
    const surfaceScalarField& mphiCorr = tmphiCorr();
    const surfaceScalarField::Boundary& mphiCorrBf = mphiCorr.boundaryField();

    forAll(phiCorrIf, facei)
    {
        const label own = owner[facei];
        const label nei = neighb[facei];

        psiMaxn[own] = max(psiMaxn[own], mpsiNei[facei]);
        psiMinn[own] = min(psiMinn[own], mpsiNei[facei]);

        psiMaxn[nei] = max(psiMaxn[nei], mpsiOwn[facei]);
        psiMinn[nei] = min(psiMinn[nei], mpsiOwn[facei]);

        sumPhiBD[own] += mphiBD[facei];
        sumPhiBD[nei] -= mphiBD[facei];

        const scalar phiCorrf = mphiCorr[facei];

        if (phiCorrf > 0)
        {
            sumPhip[own] += phiCorrf;
            mSumPhim[nei] += phiCorrf;
        }
        else
        {
            mSumPhim[own] -= phiCorrf;
            sumPhip[nei] -= phiCorrf;
        }
    }

        forAll(phiCorrBf, patchi)
    {
        const fvPatchScalarField& psiPf = psiBf[patchi];
        const fvsPatchScalarField& mpsiNeiPf = mpsiNei.boundaryField()[patchi];
        const scalarField& mphiBDPf = mphiBDBf[patchi];
        const scalarField& mphiCorrPf = mphiCorrBf[patchi];

        const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

        if (psiPf.coupled())
        {
            forAll(mphiCorrPf, pFacei)
            {
                const label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], mpsiNeiPf[pFacei]);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], mpsiNeiPf[pFacei]);
            }
        }
        else if (psiPf.fixesValue())
        {
            forAll(mphiCorrPf, pFacei)
            {
                const label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], mpsiNeiPf[pFacei]);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], mpsiNeiPf[pFacei]);
            }
        }
        else
        {
            forAll(mphiCorrPf, pFacei)
            {
                const label pfCelli = pFaceCells[pFacei];

                psiMaxn[pfCelli] = max(psiMaxn[pfCelli], psiMax);
                psiMinn[pfCelli] = min(psiMinn[pfCelli], psiMin);
            }
        }

        forAll(mphiCorrPf, pFacei)
        {
            const label pfCelli = pFaceCells[pFacei];

            sumPhiBD[pfCelli] += mphiBDPf[pFacei];

            const scalar phiCorrf = mphiCorrPf[pFacei];

            if (phiCorrf > 0)
            {
                sumPhip[pfCelli] += phiCorrf;
            }
            else
            {
                mSumPhim[pfCelli] -= phiCorrf;
            }
        }
    }

    psiMaxn = min(psiMaxn + extremaCoeff*(psiMax - psiMin), psiMax);
    psiMinn = max(psiMinn - extremaCoeff*(psiMax - psiMin), psiMin);

    if (smoothLimiter > SMALL)
    {
        psiMaxn =
            min(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMaxn, psiMax);
        psiMinn =
            max(smoothLimiter*psiIf + (1.0 - smoothLimiter)*psiMinn, psiMin);
    }

    if (mesh.moving())
    {
        tmp<volScalarField::Internal> V0 = mesh.Vsc0();

        psiMaxn =
            V
           *(
               (rho.field()*rDeltaT - Sp.field())*psiMaxn
             - Su.field()
            )
          - (V0().field()*rDeltaT)*rho.oldTime().field()*psi0
          + sumPhiBD;

        psiMinn =
            V
           *(
               Su.field()
             - (rho.field()*rDeltaT - Sp.field())*psiMinn
            )
          + (V0().field()*rDeltaT)*rho.oldTime().field()*psi0
          - sumPhiBD;
    }
    else
    {
        psiMaxn =
            V
           *(
               (rho.field()*rDeltaT - Sp.field())*psiMaxn
             - Su.field()
             - (rho.oldTime().field()*rDeltaT)*psi0
            )
          + sumPhiBD;

        psiMinn =
            V
           *(
               Su.field()
             - (rho.field()*rDeltaT - Sp.field())*psiMinn
             + (rho.oldTime().field()*rDeltaT)*psi0
            )
          - sumPhiBD;
    }

    volScalarField sumlPhip("sumlPhip", psi.db(), mesh, phiCorr.dimensions());
    volScalarField mSumlPhim("sumlPhim", psi.db(), mesh, phiCorr.dimensions());

    for (int j=0; j<nLimiterIter; j++)
    {
        sumlPhip = dimensionedScalar("sumlPhip", phiCorr.dimensions(), 0.0);
        mSumlPhim = dimensionedScalar("sumlPhim", phiCorr.dimensions(), 0.0);

        forAll(lambdaIf, facei)
        {
            const label own = owner[facei];
            const label nei = neighb[facei];

            scalar lambdaPhiCorrf = lambdaIf[facei]*mphiCorr[facei];

            if (lambdaPhiCorrf > 0)
            {
                sumlPhip[own] += lambdaPhiCorrf;
                mSumlPhim[nei] += lambdaPhiCorrf;
            }
            else
            {
                mSumlPhim[own] -= lambdaPhiCorrf;
                sumlPhip[nei] -= lambdaPhiCorrf;
            }
        }

        forAll(lambdaBf, patchi)
        {
            scalarField& lambdaPf = lambdaBf[patchi];
            const scalarField& mphiCorrfPf = mphiCorrBf[patchi];

            const labelList& pFaceCells = mesh.boundary()[patchi].faceCells();

            forAll(lambdaPf, pFacei)
            {
                const label pfCelli = pFaceCells[pFacei];
                const scalar lambdaPhiCorrf =
                    lambdaPf[pFacei]*mphiCorrfPf[pFacei];

                if (lambdaPhiCorrf > 0)
                {
                    sumlPhip[pfCelli] += lambdaPhiCorrf;
                }
                else
                {
                    mSumlPhim[pfCelli] -= lambdaPhiCorrf;
                }
            }
        }

        forAll(sumlPhip, celli)
        {
            sumlPhip[celli] =
                max(min
                (
                    (sumlPhip[celli] + psiMaxn[celli])
                   /(mSumPhim[celli] - ROOTVSMALL),
                    1.0), 0.0
                );

            mSumlPhim[celli] =
                max(min
                (
                    (mSumlPhim[celli] + psiMinn[celli])
                   /(sumPhip[celli] + ROOTVSMALL),
                    1.0), 0.0
                );
        }

        const volScalarField& lambdam = sumlPhip;
        const volScalarField& lambdap = mSumlPhim;

        surfaceScalarField mlambdamOwn(fvc::ownerField(lambdam));
        surfaceScalarField mlambdamNei(fvc::neighbourField(lambdam));
        surfaceScalarField mlambdapOwn(fvc::ownerField(lambdap));
        surfaceScalarField mlambdapNei(fvc::neighbourField(lambdap));

        fvc::applyFaceMaskTo(mlambdamOwn);
        fvc::applyFaceMaskTo(mlambdamNei);
        fvc::applyFaceMaskTo(mlambdapOwn);
        fvc::applyFaceMaskTo(mlambdapNei);

        // Add and sync internal and boundary values for min/max
        fvc::sumIndirectAndInternalFaces(mlambdamOwn, mlambdamNei);
        fvc::sumIndirectAndInternalFaces(mlambdapOwn, mlambdapNei);

        forAll(lambdaIf, facei)
        {
            if (phiCorrIf[facei] > 0)
            {
                lambdaIf[facei] = min
                (
                    lambdaIf[facei],
                    min(mlambdapOwn[facei], mlambdamNei[facei])
                );
            }
            else
            {
                lambdaIf[facei] = min
                (
                    lambdaIf[facei],
                    min(mlambdamOwn[facei], mlambdapNei[facei])
                );
            }
        }

        forAll(lambdaBf, patchi)
        {
            fvsPatchScalarField& lambdaPf = lambdaBf[patchi];
            const scalarField& phiCorrfPf = phiCorrBf[patchi];
            const fvPatchScalarField& psiPf = psiBf[patchi];
            const scalarField& pmlambdapOwn =
                mlambdapOwn.boundaryField()[patchi];
            const scalarField& pmlambdamOwn =
                mlambdamOwn.boundaryField()[patchi];

            if (isA<wedgeFvPatch>(mesh.boundary()[patchi]))
            {
                lambdaPf = 0;
            }
            else if (psiPf.coupled())
            {
                forAll(lambdaPf, pFacei)
                {
                    if (phiCorrfPf[pFacei] > 0)
                    {
                        lambdaPf[pFacei] =
                            min(lambdaPf[pFacei], pmlambdapOwn[pFacei]);
                    }
                    else
                    {
                        lambdaPf[pFacei] =
                            min(lambdaPf[pFacei], pmlambdamOwn[pFacei]);
                    }
                }
            }
        }

        // Take minimum value of limiter across coupled patches
        surfaceScalarField::Boundary lambdaNbrBf
        (
            surfaceScalarField::Internal::null(),
            lambdaBf.boundaryNeighbourField()
        );
        forAll(lambdaBf, patchi)
        {
            fvsPatchScalarField& lambdaPf = lambdaBf[patchi];
            const fvsPatchScalarField& lambdaNbrPf = lambdaNbrBf[patchi];
            if (lambdaPf.coupled())
            {
                lambdaPf = min(lambdaPf, lambdaNbrPf);
            }
        }
    }
}


template<class RdeltaTType, class RhoType, class SpType, class SuType>
void Foam::MULES::limit
(
    const RdeltaTType& rDeltaT,
    const RhoType& rho,
    const volScalarField& psi,
    const surfaceScalarField& phi,
    surfaceScalarField& phiPsi,
    const SpType& Sp,
    const SuType& Su,
    const scalar psiMax,
    const scalar psiMin,
    const bool returnCorr
)
{
    const fvMesh& mesh = psi.mesh();

    surfaceScalarField phiBD(upwind<scalar>(psi.mesh(), phi).flux(psi));

    surfaceScalarField::Boundary& phiBDBf = phiBD.boundaryFieldRef();
    const surfaceScalarField::Boundary& phiPsiBf = phiPsi.boundaryField();

    forAll(phiBDBf, patchi)
    {
        fvsPatchScalarField& phiBDPf = phiBDBf[patchi];

        if (!phiBDPf.coupled())
        {
            phiBDPf = phiPsiBf[patchi];
        }
    }

    surfaceScalarField& phiCorr = phiPsi;
    phiCorr -= phiBD;

    surfaceScalarField lambda
    (
        IOobject
        (
            "lambda",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        mesh,
        dimensionedScalar("lambda", dimless, 1)
    );

    limiter
    (
        lambda,
        rDeltaT,
        rho,
        psi,
        phiBD,
        phiCorr,
        Sp,
        Su,
        psiMax,
        psiMin
    );

    if (returnCorr)
    {
        phiCorr *= lambda;
    }
    else
    {
        phiPsi = phiBD + lambda*phiCorr;
    }
}


template<class SurfaceScalarFieldList>
void Foam::MULES::limitSum(SurfaceScalarFieldList& phiPsiCorrs)
{
    {
        UPtrList<scalarField> phiPsiCorrsInternal(phiPsiCorrs.size());
        forAll(phiPsiCorrs, phasei)
        {
            phiPsiCorrsInternal.set(phasei, &phiPsiCorrs[phasei]);
        }

        limitSum(phiPsiCorrsInternal);
    }

    const surfaceScalarField::Boundary& bfld =
        phiPsiCorrs[0].boundaryField();

    forAll(bfld, patchi)
    {
        if (bfld[patchi].coupled())
        {
            UPtrList<scalarField> phiPsiCorrsPatch(phiPsiCorrs.size());
            forAll(phiPsiCorrs, phasei)
            {
                phiPsiCorrsPatch.set
                (
                    phasei,
                    &phiPsiCorrs[phasei].boundaryFieldRef()[patchi]
                );
            }

            limitSum(phiPsiCorrsPatch);
        }
    }
}


// ************************************************************************* //
