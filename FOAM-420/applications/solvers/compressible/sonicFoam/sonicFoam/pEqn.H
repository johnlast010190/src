rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phid
(
    "phid",
    fvc::interpolate(psi)
   *(
        fvc::flux(HbyA)
      + rhorAUf*fvc::ddtCorr(rho, U, phi)/fvc::interpolate(rho)
    )
);

{
    tmp<surfaceScalarField> psif(fvc::interpolate(psi));
    fvOptions.makeRelative(psif, phid);
}

// Non-orthogonal pressure corrector loop
while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi, p)
      + fvm::div(phid, p) + fvc::div(phid*thermo.pRef())
      - fvm::laplacian(rhorAUf, p)
     ==
        fvOptions(psi, p, rho.name())
    );

    pEqn.solve();

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phid*thermo.pRef() + pEqn.flux();
    }
}

#include "cfdTools/compressible/rhoEqn.H"
#include "cfdTools/compressible/compressibleContinuityErrs.H"

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);
