{
    // Correct p using Newton's method to solve p-rho(p)g.h = p_rgh
    Info<< "Correcting and writing field p from p_rgh\n" << endl;
    label iter = 0;
    while (1)
    {
        tmp<volScalarField> deltap =
            (p-thermo.rho()*gh-p_rgh) /
            (1.0-thermo.psi()*gh);
        p -= deltap();
        thermo.correct();
        dimensionedScalar smallp("smallp", p.dimensions(), SMALL);
        if (max(mag(deltap/stabilise(p, smallp))).value() < 1e-5)
        {
            break;
        }
        if (iter++ == 100)
        {
            FatalErrorInFunction
                << "Unable to determine p from p_rgh; maximum iterations reached"
                << endl
                << exit(FatalError);
        }
    }
    rho = thermo.rho();

    // Rewrite correct p to replace temporary one written earlier
    p.write();
}

    pressureControl pressureControl
    (
        p,
        p_rgh,
        rho,
        simple.dict(),
        thermo.incompressible()
    );
{
    if (thermo.incompressible() && p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pressureControl.refValue()
          - getRefCellValue(p, pressureControl.refCell())
        );
    }

    // Force p_rgh to be consistent with p
    p_rgh = p - rho*gh;

    mesh.schemes().setFluxRequired(p_rgh.name());
}
