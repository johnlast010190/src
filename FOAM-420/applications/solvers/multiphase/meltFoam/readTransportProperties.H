   IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    dimensionedScalar PI = constant::mathematical::pi;

    // solid -> phase 1
    // liquid -> phase 2

    // Reading density rho
    dimensionedScalar rho("rho", transportProperties.lookup("rho"));

    // Reading thermal conductivity lambda
    dimensionedScalar lambdaS
    (
        "lambdaS", transportProperties.lookup("lambdaS")
    );
    dimensionedScalar lambdaL
    (
        "lambdaL", transportProperties.lookup("lambdaL")
    );

    // Reading heat capacity cp
    dimensionedScalar cpS("cpS", transportProperties.lookup("cpS"));
    dimensionedScalar cpL("cpL", transportProperties.lookup("cpL"));

    // Reading kinematic viscosity
    dimensionedScalar nuS("nuS", transportProperties.lookup("nuS"));
    dimensionedScalar nuL("nuL", transportProperties.lookup("nuL"));

    // Reading latent heat of fusion hs
    dimensionedScalar hs("hs", transportProperties.lookup("hs"));

    // Reading solid bound of melting temperature Ts
    dimensionedScalar Ts("Ts", transportProperties.lookup("Ts"));

    // Reading liquid bound of melting temperature Tl
    dimensionedScalar Tl("Tl", transportProperties.lookup("Tl"));

    // Reading volume expansion factor beta
    dimensionedScalar beta("beta", transportProperties.lookup("beta"));

    // Reading large Darcy-type source term constant DCl
    dimensionedScalar DCl("DCl", transportProperties.lookup("DCl"));

    // Reading small Darcy-type source term constant DCs
    dimensionedScalar DCs("DCs", transportProperties.lookup("DCs"));

