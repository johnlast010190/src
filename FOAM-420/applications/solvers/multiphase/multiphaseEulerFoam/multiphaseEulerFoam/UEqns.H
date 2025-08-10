
PtrList<fvVectorMatrix> UEqns(fluid.phases().size());
autoPtr<multiphaseSystem::dragCoeffFields> dragCoeffs(fluid.dragCoeffs());

int phasei = 0;
forAllIter(PtrDictionary<phaseModel>, fluid.phases(), iter)
{
    phaseModel& phase = iter();
    const volScalarField& alpha = phase;
    volVectorField& U = phase.U();

    volScalarField nuEff(turbulence->nut() + iter().nu());
    // Store nu for porosity source
    volScalarField nu
    (
        IOobject
        (
            "nu."+phase.name(),
            runTime.timeName(),
            mesh
        ),
        mesh,
        iter().nu(),
        calculatedFvPatchScalarField::typeName
    );

    UEqns.set
    (
        phasei,
        new fvVectorMatrix
        (
            fvm::ddt(alpha, U)
          + fvOptions.MRFDDt
            (
                alpha*(1 + (1/phase.rho())*fluid.Cvm(phase)),
                U
            )
          + fvm::div(phase.alphaPhi(), U)

          + (alpha/phase.rho())*fluid.Cvm(phase)*
            (
                fvm::ddt(U)
              + fvm::div(phase.phi(), U)
              - fvm::Sp(fvc::div(phase.phi()), U)
            )

          - fvm::laplacian(alpha*nuEff, U)
          - fvc::div
            (
                alpha*(nuEff*dev(T(fvc::grad(U))) /*- ((2.0/3.0)*I)*k*/),
                "div(Rc)"
            )
         ==
        //- fvm::Sp(fluid.dragCoeff(phase, dragCoeffs())/phase.rho(), U)
        //- (alpha*phase.rho())*fluid.lift(phase)
          //+
            (alpha/phase.rho())*fluid.Svm(phase)
          - fvm::Sp
            (
                slamDampCoeff
               *max
                (
                    mag(U()) - maxSlamVelocity,
                    dimensionedScalar("U0", dimVelocity, 0)
                )
               /pow(mesh.V(), 1.0/3.0),
                U
            )
          + alpha*fvOptions(U)
        )
    );
    //UEqns[phasei].relax();

    fvOptions.constrain(UEqns[phasei]);

    phasei++;
}
