/*---------------------------------------------------------------------------*\
| Modified 2010-2017 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
if (!pimple.SIMPLErho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::flux(rho*HbyA)
);
if (ddtPhiCorr)
{
    phiHbyA += rhorAUf*fvc::ddtCorr(rho, U, phi);
}


tmp<surfaceScalarField> rhof(fvc::interpolate(rho));
fvOptions.makeRelative(rhof, phiHbyA);


// Update the pressure BCs to ensure flux consistency
constrainPressure(p, rho, U, phiHbyA, rhorAUf, fvOptions);

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        (fvc::interpolate(psi)/rhof())*phiHbyA
    );

    phiHbyA -= fvc::interpolate(psi*p)*phiHbyA/rhof();

    // For correct behaviour this therm has to be after constrain pressure
    // It compensates the fvc::ddt(rho) that produce additional flux (mesh)
    // Even for constant rho
    fvc::makeRelative(phiHbyA, rho, U);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
          + fvm::div(phid, p)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        // Relax the pressure equation to ensure diagonal-dominance
        pEqn.relax();

        pEqn.solve(mesh.solution().solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}
else
{
    // For correct behaviour this therm has to be after constrain pressure
    // It compensates the fvc::ddt(rho) that produce additional flux (mesh)
    // Even for constant rho
    fvc::makeRelative(phiHbyA, rho, U);

    while (pimple.correctNonOrthogonal())
    {
        // Pressure corrector
        fvScalarMatrix pEqn
        (
            fvc::ddt(rho) + psi*correction(fvm::ddt(p))
          + fvc::div(phiHbyA)
          - fvm::laplacian(rhorAUf, p)
         ==
            fvOptions(psi, p, rho.name())
        );

        pEqn.solve(mesh.solution().solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi = phiHbyA + pEqn.flux();
        }
    }
}

#include "cfdTools/compressible/rhoEqn.H"
#include "cfdTools/compressible/compressibleContinuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

U = HbyA - rAU*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);
K = 0.5*magSqr(U);

if (pressureControl.limit(p))
{
    p.correctBoundaryConditions();
}

thermo.correctRho(psi*p - psip0, rhoMin, rhoMax);
rho = thermo.rho();

{
    rhoUf = fvc::interpolate(rho*U);
    surfaceVectorField n(mesh.Sf()/mesh.magSf());
    rhoUf += n*(fvc::absolute(phi, rho, U)/mesh.magSf() - (n & rhoUf));
}

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}
