{
    autoPtr<phaseSystem::massTransferTable>
        massTransferPtr(fluid.massTransfer());

    phaseSystem::massTransferTable&
        massTransfer(massTransferPtr());

    forAll(phases, phasei)
    {
        phaseModel& phase = phases[phasei];

        PtrList<volScalarField>& Y = phase.Y();
        const volScalarField& alpha = phase;

        tmp<volScalarField> trho(phase.rho());
        const volScalarField& rho = trho();

        forAll(Y, i)
        {
            tmp<fvScalarMatrix> YiEqn(phase.YiEqn(Y[i]));

            if (YiEqn.valid())
            {
                YiEqn =
                (
                    YiEqn
                  ==
                    *massTransfer[Y[i].name()]
                   + fvOptions(alpha, rho, Y[i])
                );

                YiEqn->relax();
                YiEqn->solve(mesh.solution().solver("Yi"));
            }
        }
    }

    fluid.massTransfer(); // updates interfacial mass flow rates
}
