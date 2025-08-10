// Solving the energy equation

{
    fvScalarMatrix TEqn
    (
        fvm::ddt(cp, T)
      + fvm::div(phi*fvc::interpolate(cp), T)
      + hs*4.0*exp(-pow(4.0*(T-Tmelt)/(Tl-Ts),2))/Foam::sqrt(PI)/(Tl-Ts)*fvm::ddt(T)
      + hs*4.0*exp(-pow(4.0*(T-Tmelt)/(Tl-Ts),2))/Foam::sqrt(PI)/(Tl-Ts)*(U & fvc::grad(T))
      - fvm::laplacian(lambda/rho, T)
    );

    TEqn.relax();
    TEqn.solve();

    alpha = 0.5*Foam::erf(4.0*(T-Tmelt)/(Tl-Ts))+scalar(0.5);

    cp = alpha*cpL+(1.0-alpha)*cpS;
    lambda = alpha*lambdaL+(1.0-alpha)*lambdaS;
    nu = alpha*nuL+(1.0-alpha)*nuS;
    rhok = 1.0 - Foam::max(beta*(T - Tl),scalar(0.0));
    DC = DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs);
}
