/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
{
    volScalarField rAU("rAU", 1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
    volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

    surfaceScalarField phiHbyA
    (
        "phiHbyA",
        fvc::flux(HbyA)
      + fvOptions.zeroFilter(rAUf*fvc::ddtCorr(U, phi))
    );

    fvOptions.makeRelative(phiHbyA);

    autoPtr<volVectorField> Fb;
    autoPtr<surfaceScalarField> phig;

    if (cellBuoyancy)
    {
        dimensionedScalar dimlessOne("1", dimless, 1);
        Fb.reset
        (
            new volVectorField(g*(rhok - dimlessOne))

        );

        Fb->correctBoundaryConditions();

        phiHbyA += (rAUf*(fvc::interpolate(Fb())) & mesh.Sf());
    }
    else
    {
        phig.reset
        (
            new surfaceScalarField
            (-rAUf*ghf*fvc::snGrad(rhok)*mesh.magSf())
        );

        phiHbyA += phig();
    }

    // Update the pressure BCs to ensure flux consistency
    constrainPressure(p_rgh, U, phiHbyA, rAUf, fvOptions);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix p_rghEqn
        (
            fvm::laplacian(rAUf, p_rgh) == fvc::div(phiHbyA)
        );

        p_rghEqn.setReference(pRefCell, getRefCellValue(p_rgh, pRefCell));

        p_rghEqn.solve(mesh.solution().solver(p_rgh.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            // Calculate the conservative fluxes
            phi = phiHbyA - p_rghEqn.flux();

            // Explicitly relax pressure for momentum corrector
            p_rgh.relax();

            // Correct the momentum source with the pressure gradient flux
            // calculated from the relaxed pressure
            if (cellBuoyancy)
            {
                U = HbyA - rAU*fvc::grad(p_rgh)
                    + rAU*fvc::reconstructSmooth(Fb());
            }
            else
            {
                U = HbyA + rAU*fvc::reconstruct
                (
                    (phig() - p_rghEqn.flux())/rAUf
                );
            }

            U.correctBoundaryConditions();
            fvOptions.correct(U);
        }
    }

    #include "cfdTools/incompressible/continuityErrs.H"

    p = p_rgh + rhok*gh;

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rhok*gh;
    }
}
