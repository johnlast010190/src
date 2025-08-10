/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
    // Solve the Momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + fvOptions.MRFDDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
