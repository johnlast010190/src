Info<< "Constructing momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(phases.size());

{
    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransfer());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        const volScalarField& alpha = phase;
        tmp<volScalarField> trho(phase.rho());
        const volScalarField& rho = trho();
        volVectorField& U = phase.U();

        UEqns.set
        (
            phasei,
            new fvVectorMatrix
            (
                phase.UEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvOptions(alpha, rho, U)
            )
        );

        UEqns[phasei].relax();
        fvOptions.constrain(UEqns[phasei]);
        fvOptions.correct(U);
    }
}
