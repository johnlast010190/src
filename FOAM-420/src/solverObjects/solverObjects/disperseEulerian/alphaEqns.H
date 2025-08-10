fv::options& fvOptions = this->fvOptions();

fvScalarMatrix alphaEqn
(
    fvm::ddt(alphad(phase_no))
  + fvm::div(phid(phase_no),alphad(phase_no),"div(phid,alphad)")
  ==
    fvOptions(alphad(phase_no))
);

if (solveTD_ && alphaTD_)
{
    alphaEqn += -fvm::laplacian((*nutd),alphad(phase_no),"laplacian(nutd,alphad)");
}
alphaEqn.relax();

fvOptions.constrain(alphaEqn);
alphaEqn.solve();
fvOptions.correct(alphad(phase_no));

alphad(phase_no).correctBoundaryConditions();


//patch reports
if (reportOn_)
{
    Info<<"Time integration for reports is: "<<timeIntOn_<<endl;

    const scalar currentTime = runTime.time().value();
    const surfaceScalarField alpha_flux = alphaEqn.flux()();
    const scalar dT = runTime.deltaTValue();
    const surfaceScalarField& magSf = mesh_.magSf();

    forAll(sMass_[phase_no].boundaryField(), patchI)
    {
    if (!sMass_[phase_no].boundaryField()[patchI].coupled())
    {
        if (timeIntOn_)
        {
        if (currentTime>startTime_)
        {
            // this is [Kg/m2]
            sMass_[phase_no].boundaryFieldRef()[patchI] +=
            (alpha_flux.boundaryField()[patchI]/
             magSf.boundaryField()[patchI])*dT*rhop.value();
        }
        }
        else
        {
        // this is [Kg/m2s]
        sMass_[phase_no].boundaryFieldRef()[patchI] =
            (alpha_flux.boundaryField()[patchI]/
             magSf.boundaryField()[patchI])*rhop.value();
        }
    }
    }

    forAll(gMass_[phase_no].boundaryField(), patchI)
    {
        if (!gMass_[phase_no].boundaryField()[patchI].coupled())
        {
            if (timeIntOn_)
            {
            if (currentTime>startTime_)
            {
                // this is [Kg]
                gMass_[phase_no].boundaryFieldRef()[patchI] +=
                alpha_flux.boundaryField()[patchI]*
                dT*rhop.value();

            if (mesh_.time().outputTime()) // to limit the output in the log file
            {
            scalar massOnPatch =gSum(gMass_[phase_no].boundaryField()[patchI]);
            Info<<"Phase"<<phase_no<<" on patch "<<mesh_.boundaryMesh()[patchI].name()
                <<" is "<<massOnPatch<<" [Kg]."<<endl;
            }
            }
            }
            else
            {
            // this is [Kg/s]
            gMass_[phase_no].boundaryFieldRef()[patchI] =
                alpha_flux.boundaryField()[patchI]*rhop.value();

        if (mesh_.time().outputTime()) // to limit the output in the log file
        {
            scalar massOnPatch =gSum(gMass_[phase_no].boundaryField()[patchI]);
            Info<<"Phase"<<phase_no<<" on patch "<<mesh_.boundaryMesh()[patchI].name()
            <<" is "<<massOnPatch<<" [Kg/s]."<<endl;
        }
            }
        }
    }
}//end if report
