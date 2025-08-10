{
    // store additional drift flux component from mixture
    tmp<volVectorField> Hvol = mixture.Hvol();

    // iterator over all phases
    PtrDictionary<phase>::iterator iter = phases.begin();

    // Reset rhoPhi to zero
    rhoPhi = dimensionedScalar("0", dimensionSet(1, 0, -1, 0, 0), 0);

    // phase counter
    label phasei=0;

    // phase loop (disperse phases only)
    for (iter = phases.begin(); iter != phases.end(); ++iter)
    {
        if (iter().name() != mixture.continuousPhase())
        {
            phase& alpha = iter();

            surfaceScalarField alphaPhi
            (
                IOobject
                (
                    "alphaPhi",
                    runTime.timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("0", phi.dimensions(), 0)
            );

            if (nAlphaSubCycles > 1)
            {
                dimensionedScalar totalDeltaT = runTime.deltaT();
                surfaceScalarField alphaPhiSum
                (
                    IOobject
                    (
                        "alphaPhiSum",
                        runTime.timeName(),
                        mesh
                    ),
                    mesh,
                    dimensionedScalar("0", phi.dimensions(), 0)
                );

                for
                (
                    subCycle<volScalarField> alphaSubCycle(alpha, nAlphaSubCycles);
                    !(++alphaSubCycle).end();
                )
                {
                    if (MULESadvection)
                    {
                        #include "alphaEqn.H"
                    }
                    else
                    {
                        #include "alphaEqnFVM.H"
                    }
                    alphaPhiSum += (runTime.deltaT()/totalDeltaT)*alphaPhi;
                }

                alphaPhi = alphaPhiSum;
            }
            else
            {
                if (MULESadvection)
                {
                    #include "alphaEqn.H"
                }
                else
                {
                    #include "alphaEqnFVM.H"
                }
            }

            // Apply the diffusion term separately to allow implicit solution
            // and boundedness of the explicit advection
            if (MULESadvection && alphaDiffusion)
            {
                fvScalarMatrix alphaEqn
                (
                    fvm::ddt(alpha) - fvc::ddt(alpha)
                  - fvm::laplacian(turbulence->nut()*rSct, alpha)
                );

                alphaEqn.solve(mesh.solution().solver("alphaDiffusion"));

                alphaPhi += alphaEqn.flux();
                Info<< alpha.name() << " volume fraction = "
                    << alpha.weightedAverage(mesh.Vsc()).value()
                    << "  Min(" << alpha.name() << ") = " << min(alpha).value()
                    << "  Max(" << alpha.name() << ") = " << max(alpha).value()
                    << endl;
            }

            rhoPhi += alphaPhi*(alpha.rho() - mixture.phasec().rho());
        }

        phasei++;
    }

    // compute water volume fraction
    volScalarField& alphac = mixture.alphac();
    alphac = scalar(1);
    for (iter = phases.begin(); iter != phases.end(); ++iter)
    {
        if (iter().name() != mixture.continuousPhase())
        {
            alphac -= iter();
        }
    }
    alphac.max(0); // prevent phase going negative

    // add water contribution
    rhoPhi += phi*mixture.phasec().rho();
    rho = mixture.rho();
}
