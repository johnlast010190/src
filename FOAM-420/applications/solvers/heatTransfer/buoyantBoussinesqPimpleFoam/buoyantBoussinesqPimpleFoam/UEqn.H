/*---------------------------------------------------------------------------*\
| Modified 2010-2017 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/

    // Solve the momentum equation

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::ddt(U) + fvm::div(phi, U)
      + fvOptions.MRFDDt(U)
      + turbulence->divDevReff(U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        if (cellBuoyancy)
        {
            dimensionedScalar dimlessOne("1", dimless, 1);
            solve
            (
                UEqn == -fvc::grad(p_rgh)
                        + fvc::reconstructSmooth(g*(rhok - dimlessOne))
                        //smoothed cell average body force
            );
        }
        else
        {
            solve
            (
                UEqn
                ==
                fvc::reconstruct
                (
                    (
                      - ghf*fvc::snGrad(rhok)
                      - fvc::snGrad(p_rgh)
                    )*mesh.magSf()
                )
            );
        }
    }
    fvOptions.correct(U);
