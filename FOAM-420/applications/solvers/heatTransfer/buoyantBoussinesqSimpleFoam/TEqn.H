/*---------------------------------------------------------------------------*\
| Modified 2010-2015 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
{
    volScalarField kappaEff
    (
        "kappaEff",
        turbulence->alphaEff()
    );

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
     ==
        fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    radiation.correct();

    rhok = updateRho(laminarTransport, T)();
    rhok.correctBoundaryConditions();
    rhoEff = turbulence->rho()() * rhok;
}
