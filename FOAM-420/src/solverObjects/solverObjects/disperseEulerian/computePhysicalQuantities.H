
// get total drag force coefficient
volScalarField Kd(forcesD_[phase_no].K());

if (solveTD_)
{
    // stabilize for division by zero
    dimensionedScalar stabE("stabE", epsilon()().dimensions(), VSMALL);
    dimensionedScalar stabU("stabU", dimVelocity, VSMALL);
    dimensionedScalar stabK("stabK", Kd.dimensions(), VSMALL);

    // relate disperse phase turbulence to continuous phase
    // turbulence based on time scale ratio
    volScalarField tTurb( 0.2*(k()/(epsilon() + stabE)) );   // standard in windDrivenRain
    volScalarField tFluid( 1./(Kd + stabK) );                // standard in windDrivenRain

    // model taken from Lagrangian library StochasticDispersionRAS.H
    if (stochasticDispersionModel_)
    {
        const scalar cps = 0.16432;
        tTurb = min(k()/(epsilon() + stabE), cps*pow(k(), 1.5)/(epsilon() + stabE)/(phases_[phase_no].magUr() + stabU));
    }

    Ctd = new volScalarField
    (
        sqrt( tTurb/(tTurb + tFluid) )
    );

    nutd = new volScalarField
    (
        "nutd",
        nut() * sqr(*Ctd)
    );
}
