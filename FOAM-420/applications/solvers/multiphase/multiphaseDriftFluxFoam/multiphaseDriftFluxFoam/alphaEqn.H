{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    if (MULESCorr)
    {
        fvScalarMatrix alphaEqn
        (
            fv::EulerDdtScheme<scalar>(mesh).fvmDdt(alpha)
          + fv::gaussConvectionScheme<scalar>
            (
                mesh,
                phi,
                upwind<scalar>(mesh, phi)
            ).fvmDiv(phi, alpha)
        );

        solve(alphaEqn);

        Info<< alpha.name() << " volume fraction = "
            << alpha.weightedAverage(mesh.Vsc()).value()
            << "  Min(" << alpha.name() << ") = " << min(alpha).value()
            << "  Max(" << alpha.name() << ") = " << max(alpha).value()
            << endl;

        tmp<surfaceScalarField> talphaPhiUD(alphaEqn.flux());
        alphaPhi = talphaPhiUD();

        if (alphaApplyPrevCorr && talphaPhiCorr0[phasei].valid())
        {
            Info<< "Applying the previous iteration correction flux" << endl;

            MULES::correct
            (
                alpha,
                alphaPhi,
                talphaPhiCorr0[phasei].ref(),
                mixture.alphaMax(),
                0
            );

            alphaPhi += talphaPhiCorr0[phasei];
        }

        // Cache the upwind-flux
        talphaPhiCorr0[phasei] = talphaPhiUD;
    }

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phir(fvc::flux(alpha.Udm()));
        #include "computeDriftFlux.H"

        tmp<surfaceScalarField> talphaPhiUn
        (
            fvc::flux
            (
                phi,
                alpha,
                alphaScheme
            )
          + fvc::flux
            (
                phir,
                alpha,
                alpharScheme
            )
        );

        if (MULESCorr)
        {
            tmp<surfaceScalarField> talphaPhiCorr(talphaPhiUn() - alphaPhi);
            volScalarField alpha10("alpha10", alpha);

            MULES::correct
            (
                alpha,
                talphaPhiUn(),
                talphaPhiCorr.ref(),
                mixture.alphaMax(),
                0
            );

            // Under-relax the correction for all but the 1st corrector
            if (aCorr == 0)
            {
                alphaPhi += talphaPhiCorr();
            }
            else
            {
                volScalarField& alphai = alpha;
                alphai = 0.5*alphai + 0.5*alpha10;
                alphaPhi += 0.5*talphaPhiCorr();
            }
        }
        else
        {
            alphaPhi = talphaPhiUn;

            MULES::explicitSolve
            (
                alpha,
                phi,
                alphaPhi,
                mixture.alphaMax(),
                0
            );
        }
    }

    if (alphaApplyPrevCorr && MULESCorr)
    {
        talphaPhiCorr0[phasei] = alphaPhi - talphaPhiCorr0[phasei];
    }

    Info<< alpha.name() << " volume fraction = "
        << alpha.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha.name() << ") = " << min(alpha).value()
        << "  Max(" << alpha.name() << ") = " << max(alpha).value()
        << endl;
}
