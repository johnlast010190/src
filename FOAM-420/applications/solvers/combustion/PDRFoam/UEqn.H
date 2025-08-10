/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        betav*fvm::ddt(rho, U) + fvm::div(phi, U)
      + fvOptions.MRFDDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        betav*rho*g
      + betav*fvOptions(rho, U)
    );

    fvOptions.constrain(UEqn);

    volSymmTensorField invA(inv(I*UEqn.A() + drag->Dcu()));

    if (pimple.momentumPredictor())
    {
        U = invA & (UEqn.H() - betav*fvc::grad(p));
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
