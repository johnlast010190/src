dimensionedScalar rhod = phases_[phase_no].rhod();
volScalarField rhoc = phases_[phase_no].rhoc();

fvVectorMatrix UdEqn
(
    fvm::ddt(Ud(phase_no))
  + fvm::div(phid(phase_no), Ud(phase_no),"div(phid,Ud)")
  + fvm::SuSp(-fvc::div(phid(phase_no)),Ud(phase_no))
  ==
    fvOptions(Ud(phase_no))
  + g*(1.-rhoc/rhod)
  + Kd* U()
  - fvm::Sp(Kd, Ud(phase_no))
  + forcesL_[phase_no].F()
  - forcesTD_[phase_no].F()
);

if (debug)
{
    volVectorField Fd("Fd" + Foam::word(phases_[phase_no].name()), Kd * (U() - Ud(phase_no)));
    volVectorField Fl("Fl" + Foam::word(phases_[phase_no].name()), forcesL_[phase_no].F());
    volVectorField Ftd("Ftd" + Foam::word(phases_[phase_no].name()), forcesTD_[phase_no].F());

    if (runTime.outputTime())
    {
        Fd.write();
        Fl.write();
        Ftd.write();
    }
}

if (solveTD_)
{
    volTensorField Rca("Rca",-*nutd*(T(fvc::grad(Ud(phase_no)))));
    Rca = Rca + (2.0/3.0)*sqr(*Ctd)*I*k() - (2.0/3.0)*I*tr(Rca);
    Rca.correctBoundaryConditions();

    const scalar stabAlpha = stabDivAlpha_*gMax(alphad(phase_no)) + VSMALL;

    surfaceScalarField phiRa
    (
        - fvc::interpolate(*nutd)
        *mesh_.magSf()*fvc::snGrad(alphad(phase_no))/fvc::interpolate(alphad(phase_no) + stabAlpha)
    );

    UdEqn += (
        - fvm::laplacian(*nutd, Ud(phase_no),"laplacian(nutd,Ud)")
        + fvc::div(Rca)

        + fvm::div(phiRa, Ud(phase_no),"div(phiRa,Ud)")
        - fvm::Sp(fvc::div(phiRa), Ud(phase_no))
        + (fvc::grad(alphad(phase_no))/(fvc::average(alphad(phase_no)) + stabAlpha) & Rca)
    );

    if (debug)
    {
        volVectorField turbDisp
        (
            "turbDisp",
            - fvc::laplacian(*nutd, Ud(phase_no))
            + fvc::div(Rca)

            + fvc::div(phiRa, Ud(phase_no),"div(phiRa,Ud)")
            - fvc::div(phiRa) * Ud(phase_no)
            + (fvc::grad(alphad(phase_no))/(fvc::average(alphad(phase_no)) + stabAlpha) & Rca)
        );

        if (runTime.outputTime())
        {
            turbDisp.write();
        }
    }
}

UdEqn.relax();
fvOptions.constrain(UdEqn);
UdEqn.solve();
fvOptions.correct(Ud(phase_no));

Ud(phase_no).correctBoundaryConditions();
