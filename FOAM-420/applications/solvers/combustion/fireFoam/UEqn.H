/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(phi, U)
      + fvOptions.MRFDDt(rho, U)
      + turbulence->divDevRhoReff(U)
     ==
        parcels.SU(U)
      + fvOptions(rho, U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
          ==
            fvc::reconstruct
            (
                (
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
