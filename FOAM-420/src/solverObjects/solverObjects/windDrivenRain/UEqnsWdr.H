volTensorField Rca(-nutrain()*(T(fvc::grad(Urain(phase_no)))));

if (solveTD_)
{
    Rca = Rca + (2.0/3.0)*sqr(Ctrain(phase_no))*I*k() - (2.0/3.0)*I*tr(Rca);
}

Rca.correctBoundaryConditions();

surfaceScalarField phiRa
(
    - fvc::interpolate(nutrain())
    *mesh_.magSf()*fvc::snGrad(alpharain(phase_no))/fvc::interpolate(alpharain(phase_no) + stabDivAlpha_ + VSMALL)
);

fvVectorMatrix UrainEqn
(
    fvm::ddt(Urain(phase_no))
    + fvm::div(phirain(phase_no), Urain(phase_no),"div(phirain,Urain)")
    - fvm::Sp(fvc::div(phirain(phase_no)),Urain(phase_no))

    - fvm::laplacian(nutrain(), Urain(phase_no))
    + fvc::div(Rca)

    + fvm::div(phiRa, Urain(phase_no),"div(phirain,Urain)")
    - fvm::Sp(fvc::div(phiRa), Urain(phase_no))
    + (fvc::grad(alpharain(phase_no))/(fvc::average(alpharain(phase_no)) + stabDivAlpha_ + VSMALL) & Rca)
    ==
    g
    + ((3* mua_ * CdRe())/(4 * rhop_ * dp * dp))*U()
    - fvm::Sp(((3* mua_ * CdRe())/(4 * rhop_ * dp * dp)),Urain(phase_no))
);
UrainEqn.relax();

UrainEqn.solve();

Urain(phase_no).correctBoundaryConditions();

phirain(phase_no) = fvc::interpolate(Urain(phase_no)) & mesh_.Sf();
