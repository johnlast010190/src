Info<< "Constructing face momentum equations" << endl;

fvVectorMatrix U1Eqn(U1, rho1.dimensions()*U1.dimensions()*dimVol/dimTime);
fvVectorMatrix U2Eqn(U2, rho2.dimensions()*U2.dimensions()*dimVol/dimTime);

{
    volScalarField Vm(fluid.Vm());

    fvVectorMatrix UgradU1
    (
        fvm::div(phi1, U1) - fvm::Sp(fvc::div(phi1), U1)
      + fvOptions.MRFDDt(U1)
    );

    fvVectorMatrix UgradU2
    (
        fvm::div(phi2, U2) - fvm::Sp(fvc::div(phi2), U2)
      + fvOptions.MRFDDt(U2)
    );

    {
        U1Eqn =
        (
            fvm::div(alphaRhoPhi1, U1) - fvm::Sp(fvc::div(alphaRhoPhi1), U1)
          + fvOptions.MRFDDt(alpha1*rho1, U1)
          + phase1.turbulence().divDevRhoReff(U1)
          + Vm*(UgradU1 - (UgradU2 & U2))
          - fvOptions(alpha1, rho1, U1)
        );
        U1Eqn.relax();
        fvOptions.constrain(U1Eqn);
        U1.correctBoundaryConditions();
        fvOptions.correct(U1);
    }

    {
        U2Eqn =
        (
            fvm::div(alphaRhoPhi2, U2) - fvm::Sp(fvc::div(alphaRhoPhi2), U2)
          + fvOptions.MRFDDt(alpha2*rho2, U2)
          + phase2.turbulence().divDevRhoReff(U2)
          + Vm*(UgradU2 - (UgradU1 & U1))
          - fvOptions(alpha2, rho2, U2)
        );
        U2Eqn.relax();
        fvOptions.constrain(U2Eqn);
        U2.correctBoundaryConditions();
        fvOptions.correct(U2);
    }
}
