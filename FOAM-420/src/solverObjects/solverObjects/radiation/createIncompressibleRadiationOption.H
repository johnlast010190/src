// Create radiation solver object for use in legacy solvers, in place of
// radiationModel
    autoPtr<fv::option> radiationOption
    (
        new fv::SolverOption<fv::radiationSolver>
        (
            "radiationSolver",
            fv::radiationSolver::typeName,
            dictionary(),
            T.db()
        )
    );

    fv::SolverOption<fv::radiationSolver>& radiation =
        refCast<fv::SolverOption<fv::radiationSolver>>(radiationOption());

    // The fvOption is deselected if there is no radiation, but the model
    // remains stored so the reference above is valid either way
    fvOptions.addOption(radiationOption);

    Info<< "Added radiation solver " << radiation.name() << endl;
