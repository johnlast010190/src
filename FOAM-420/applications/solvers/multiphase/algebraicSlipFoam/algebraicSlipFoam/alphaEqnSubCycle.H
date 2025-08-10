
volScalarField divU(fvc::div(phi));
word alphaScheme("div(phi,alpha)");
word alpharScheme("div(phirb,alpha)");

surfaceScalarField sumVolumetricRelativeFlux(phi*0);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();
    surfaceScalarField alphaPhiSum(0.0*alphaPhi);

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        alphaPhiSum += (runTime.deltaT()/totalDeltaT)*alphaPhi;
        sumVolumetricRelativeFlux += (runTime.deltaT()/totalDeltaT)*volumetricRelativeFlux;
    }

    alphaPhi = alphaPhiSum;
    volumetricRelativeFlux = sumVolumetricRelativeFlux;
}
else
{
    #include "alphaEqn.H"
}

// Apply the diffusion term separately to allow implicit solution
// and boundedness of the explicit advection
/*
{
    fvScalarMatrix alpha1Eqn
    (
        fvm::ddt(alpha1) - fvc::ddt(alpha1)
      - fvm::laplacian(turbulence->nut(), alpha1)
    );

    alpha1Eqn.solve(mesh.solution().solver("alpha1Diffusion"));

    alphaPhi += alpha1Eqn.flux();
    alpha2 = 1.0 - alpha1;

    Info<< "Phase-1 volume fraction = "
        << alpha1.weightedAverage(mesh.Vsc()).value()
        << "  Min(" << alpha1.name() << ") = " << min(alpha1).value()
        << "  Max(" << alpha1.name() << ") = " << max(alpha1).value()
        << endl;
}

*/
mixture.correct();

rho.forceAssign(alpha1*rho1 + (scalar(1) - alpha1)*rho2);

// compute mass-based relative flux
limitedRelativeFlux = volumetricRelativeFlux/linearInterpolate(rho)*rho2;
rhoPhi = (alphaPhi - limitedRelativeFlux)*(rho1 - rho2) + phi*rho2;
