/*---------------------------------------------------------------------------*\
| Modified 2010-2012 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));
    tUEqn.clear();

    surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(rho*HbyA)
    );

    tmp<surfaceScalarField> rhof(fvc::interpolate(rho));
    fvOptions.makeRelative(rhof, phiHbyA);

    if (simple.modifiedMomentumInterp())
    {
        fvOptions.makeAbsolute(rhof, phi);

        phiHbyA
            += (1 - mesh.solution().equationRelaxationFactor(U.name()))
            *(phi - (fvc::interpolate(rho.prevIter()*U.prevIter())&mesh.Sf()));
    }

    bool closedVolume = adjustPhi(phiHbyA, U, p_rgh);

    phiHbyA += phig;

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, fvOptions);

    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rhorAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference
        (
            pressureControl.refCell(),
            pressureControl.refValue()
        );

        p_rghEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            U = HbyA + rAU*fvc::reconstruct((phig - p_rghEqn.flux())/rhorAUf);
            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    p = p_rgh + rho*gh;

    if (pressureControl.limit(p))
    {
        p.correctBoundaryConditions();
    }

    #include "cfdTools/incompressible/continuityErrs.H"

    // For closed-volume cases adjust the pressure level
    // to obey overall mass continuity
    if (closedVolume)
    {
        if (thermo.incompressible())
        {
            p += dimensionedScalar
            (
                "p",
                p.dimensions(),
                pressureControl.refValue()
              - getRefCellValue(p, pressureControl.refCell())
            );
        }
        else
        {
            p += (initialMass - fvc::domainIntegrate(psi*p))
                /fvc::domainIntegrate(psi);

        }
        p_rgh = p - rho*gh;
    }

    rho = thermo.rho();
    rho.relax();
    Info<< "rho max/min : " << max(rho).value() << " " << min(rho).value()
        << endl;
}
