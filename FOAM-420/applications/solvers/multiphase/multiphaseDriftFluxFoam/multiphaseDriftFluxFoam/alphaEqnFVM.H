{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // include turbulent diffusion
    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phir(fvc::flux(alpha.Udm()));
        #include "computeDriftFlux.H"

        fvScalarMatrix alphaEqn
        (
            fvm::ddt(alpha)
          + fvm::div(phi, alpha, alphaScheme)
          + fvm::div(phir, alpha, alpharScheme)
        );

        if (alphaDiffusion && alpha.name() != mixture.gasPhase())
        {
            // set diffusion at boundaries to zero
            volScalarField diffCoeff("nut",turbulence->nut());
            forAll(diffCoeff.boundaryField(), pI)
            {
                if (!diffCoeff.boundaryField()[pI].coupled())
                {
                    forAll(diffCoeff.boundaryField()[pI], fI)
                    {
                        diffCoeff.boundaryFieldRef()[pI][fI] = 0.0;
                    }
                }
            }
            if (mixture.gasPhase() != "none")
            {
                const volScalarField& alphaGas =
                    mesh.lookupObject<volScalarField>("alpha." + mixture.gasPhase());

                alphaEqn -=
                    fvm::laplacian(diffCoeff*rSct*pos0(0.3 - alphaGas), alpha);
            }
            else
            {
                alphaEqn -= fvm::laplacian(diffCoeff*rSct, alpha);
            }
        }

        alphaEqn.relax();
        fvOptions.constrain(alphaEqn);
        solve(alphaEqn);
        fvOptions.correct(alpha);

        alpha.max(0); // prevent phase going negative

        tmp<surfaceScalarField> talphaPhiUn(alphaEqn.flux());
        alphaPhi = talphaPhiUn;
    }

    Info<< alpha.name() << " volume fraction = "
        << alpha.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha.name() << ") = " << min(alpha).value()
        << "  Max(" << alpha.name() << ") = " << max(alpha).value()
        << endl;
}
