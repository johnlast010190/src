{
    surfaceScalarField phir(fvc::flux(UrelModel.Urel()/rho)*rho2);

    surfaceScalarField phirvol(fvc::flux(UrelModel.Urel()));


    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        volumetricRelativeFlux = fvc::flux
        (
            fvc::flux(phirvol, scalar(1) - alpha1, alpharScheme),
            alpha1,
            alpharScheme
        );

        alphaPhi =
        (
            // divergence-free volumetric flux
            fvc::flux
            (
                phiAlphaAdvection,
                alpha1,
                alphaScheme
            )
            // relative flux
          + volumetricRelativeFlux
        );

        MULES::explicitSolve(geometricOneField(), alpha1, phi, alphaPhi, zeroField(), zeroField(), 1, 0);
    }

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(alpha1) = " << min(alpha1).value()
        << "  Max(alpha1) = " << max(alpha1).value()
        << endl;
}
