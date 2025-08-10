
dimensionedScalar dp ("dp", dimensionSet(0,1,0,0,0,0,0), phases_[phase_no][0]);

volScalarField magUr( mag(U() - Urain(phase_no)) );

Re() = (magUr*dp*rhoa_)/mua_;
computeCdRe();
CdRe().correctBoundaryConditions();

if (solveTD_)
{
    // stabilize for division by zero
    dimensionedScalar stab("stab", epsilon().dimensions(), VSMALL);

    Ctrain(phase_no) = sqrt( 0.2*(k()/(epsilon() + stab)) / (0.2*(k()/(epsilon() + stab)) + ((4*rhop_*dp*dp)/(3*mua_*(CdRe() + VSMALL)))) );
}

nutrain() = nut()*sqr(Ctrain(phase_no));
