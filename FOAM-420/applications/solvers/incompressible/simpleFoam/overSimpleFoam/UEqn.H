    // Momentum predictor

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + fvOptions.MRFDDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve(UEqn == -cellMask*fvc::grad(p));

    fvOptions.correct(U);
