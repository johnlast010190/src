/*---------------------------------------------------------------------------*\
| Modified 2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
    // Solve the Momentum equation

    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
      + fvOptions.MRFDDt(rho, U)
      //+ fvc::ddt(mixture.Hdm())
      //+ fvc::div(mixture.tauDrift())
      + turbulence->divDevRhoReff(U)
     ==
        fvOptions(rho, U)
    );

    if (momentumDDtHdm)
    {
        UEqn += fvc::ddt(mixture.Hdm());
    }
    if (momentumDivTauDm)
    {
        UEqn += fvc::div(mixture.tauDrift());
    }

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
    }
    else
    {
        UEqn.boundaryManipulate(U.boundaryFieldRef());
    }
