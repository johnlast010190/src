Info<< "\nReading vibroacousticProperties" << endl;
IOdictionary vibroacousticProperties
(
    IOobject
    (
        "vibroacousticProperties",
        runTime.constant(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

dimensionedScalar c0 = dimensionedScalar(
                                              "c0",
                                              dimensionSet(0,1,-1,0,0,0,0),
                                              vibroacousticProperties.lookupOrDefault<scalar>("c0",343)
                                        );

dimensionedScalar EModul = dimensionedScalar(
                                                  "EModul",
                                                  dimensionSet(1,-1,-2,0,0,0,0),
                                                  vibroacousticProperties.lookupOrDefault<scalar>("EModul",7.e10)
                                              );

dimensionedScalar PoissonRatio = dimensionedScalar(
                                                    "PoissonRatio",
                                                    dimless,
                                                    vibroacousticProperties.lookupOrDefault<scalar>("PoissonRatio",0.3298)
                                                    );

dimensionedScalar damp0 = dimensionedScalar(
                                                "damp0",
                                                dimless,
                                                vibroacousticProperties.lookupOrDefault<scalar>("damp0",.0)
                                            );

dimensionedScalar damp1 = dimensionedScalar(
                                                "damp1",
                                                dimensionSet(0,0,-1,0,0,0,0),
                                                vibroacousticProperties.lookupOrDefault<scalar>("damp1",.0)
                                            );

dimensionedScalar damp2 = dimensionedScalar(
                                                "damp2",
                                                dimensionSet(0,0,1,0,0,0,0),
                                                vibroacousticProperties.lookupOrDefault<scalar>("damp2",.0)
                                            );

dimensionedScalar plateThickness = dimensionedScalar(
                                                        "plateThickness",
                                                        dimensionSet(0,1,0,0,0,0,0),
                                                        vibroacousticProperties.lookupOrDefault<scalar>("plateThickness",0.00315)
                                                        );

dimensionedScalar plateRho = dimensionedScalar(
                                                "plateRho",
                                                dimensionSet(1,-3,0,0,0,0,0),
                                                vibroacousticProperties.lookupOrDefault<scalar>("plateRho",2700.)
                                                );
dimensionedScalar c02 = c0*c0;
dimensionedScalar plateB         = ( EModul*pow(plateThickness,3) ) / (12.*(1.-sqr(PoissonRatio)) );
dimensionedScalar plateCrossMass = plateRho*plateThickness;
dimensionedScalar plateD         = plateB/plateCrossMass;


scalar freq = 1000;
scalar propSpeed = Foam::sqrt(2*constant::mathematical::pi*freq)
                 * Foam::sqrt(Foam::sqrt(plateD.value()));


Info<< "Propagating speed for " << freq << " Hz :" <<  propSpeed << endl;


