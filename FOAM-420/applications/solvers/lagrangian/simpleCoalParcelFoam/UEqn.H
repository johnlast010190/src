
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + fvOptions.MRFDDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho()*g
      + parcels.SU(U)
      + fvOptions(rho, U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
