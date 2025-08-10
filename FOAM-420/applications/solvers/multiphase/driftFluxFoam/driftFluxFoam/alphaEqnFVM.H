{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    // include turbulent diffusion
    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        fvScalarMatrix alpha1Eqn
        (
            fvm::ddt(alpha1)
          + fvm::div(phi, alpha1, alphaScheme)
          + fvm::div(phir, alpha1, alpharScheme)
          - fvm::laplacian(turbulence->nut(), alpha1)
        );

        alpha1Eqn.relax();

        solve(alpha1Eqn);

        tmp<surfaceScalarField> talphaPhiUn(alpha1Eqn.flux());
        alphaPhi = talphaPhiUn();
    }

    alpha2 = 1.0 - alpha1;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}
