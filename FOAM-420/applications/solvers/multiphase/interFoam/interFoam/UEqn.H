/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + fvOptions.MRFDDt(rho, U)
      + turbulence->divDevRhoReff(rho, U)
     ==
        fvOptions(rho, U)
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
                    mixture.surfaceTensionForce()
                  - ghf*fvc::snGrad(rho)
                  - fvc::snGrad(p_rgh)
                ) * mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }
    else
    {
        UEqn.boundaryManipulate(U.boundaryFieldRef());
    }
