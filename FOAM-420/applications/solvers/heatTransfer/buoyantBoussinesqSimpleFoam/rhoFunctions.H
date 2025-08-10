tmp<volScalarField> updateRho
(
    const singlePhaseTransportModel& laminarTransport,
    const volScalarField& T
)
{

    // Reference temperature [K]
    dimensionedScalar TRef
    (
        "TRef",
        dimTemperature,
        gSum(T.primitiveField()*T.mesh().V())/gSum(T.mesh().V())
    );

    if (laminarTransport.found("TRef"))
    {
        TRef = dimensionedScalar
        (
            "TRef",
            dimTemperature,
            laminarTransport.lookup("TRef")
        );
    }

    // Thermal expansion coefficient [1/K]
    dimensionedScalar beta
    (
        laminarTransport.found("beta")
      ? dimensionedScalar
        (
            "beta",
            dimless/dimTemperature,
            laminarTransport.lookup("beta")
        )
      : 1/TRef
    );


    word rhoType =
        laminarTransport.lookupOrDefault<word>("rhoType", "Beta");

    tmp<volScalarField> trhok;

    if (rhoType == "Beta")
    {
        trhok = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "trhok",
                    T.mesh().time().timeName(),
                    T.mesh()
                ),
                1.0/(1 + beta*(T - TRef))
            )
        );
    }
    else if (rhoType == "BetaLinear")
    {
        trhok = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "trhok",
                    T.mesh().time().timeName(),
                    T.mesh()
                ),
                1.0 - beta*(T - TRef)
            )
        );
    }
    else if (rhoType == "Thiesen")
    {
        dimensionedScalar a1
        (
            dimensionedScalar::lookupOrDefault
            (
                "a1",
                laminarTransport,
                dimTemperature,
                -3.983035
            )
        );
        dimensionedScalar a2
        (
            dimensionedScalar::lookupOrDefault
            (
                "a2",
                laminarTransport,
                dimTemperature,
                301.797
            )
        );
        dimensionedScalar a3
        (
            dimensionedScalar::lookupOrDefault
            (
                "a3",
                laminarTransport,
                sqr(dimTemperature),
                522528.9
            )
        );
        dimensionedScalar a4
        (
            dimensionedScalar::lookupOrDefault
            (
                "a4",
                laminarTransport,
                dimTemperature,
                69.34881
            )
        );
        dimensionedScalar a5
        (
            dimensionedScalar::lookupOrDefault
            (
                "a5",
                laminarTransport,
                dimDensity,
                999.74950
            )
        );

        dimensionedScalar T0("T0", dimTemperature, 273.15);
        dimensionedScalar rho("rho", dimDensity, laminarTransport.lookup("rho"));

        trhok = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "trhok",
                    T.mesh().time().timeName(),
                    T.mesh()
                ),
                a5*(1.0 - (sqr(T - T0 + a1) * (T - T0 + a2))
                /(a3*(T - T0 + a4)))/rho
            )
        );
    }
    else if (rhoType == "BokilBewtra")
    {
        const dictionary& dict = laminarTransport.subDict("BokilBewtraCoeffs");
        word XName = word(dict.lookup("XName"));
        dimensionedScalar rho("rho", dimDensity, laminarTransport.lookup("rho"));

        if (T.mesh().foundObject<volScalarField>(XName))
        {
            dimensionedScalar Xc = dimensionedScalar::lookupOrDefault
            (
                "Xc",
                dict,
                dimDensity,
                0.7
            );
            const volScalarField& X = T.mesh().lookupObject<volScalarField>(XName);
            dimensionedScalar fac = dimensionedScalar::lookupOrDefault
            (
                "fac",
                dict,
                dimDensity/X.dimensions(),
                1.0
            );
            dimensionedScalar rhoSolid
            (
                "rhoSolid",
                dimDensity,
                laminarTransport.lookup("rhoSolid")
            );

            volScalarField volumeFraction
            (
                min
                (
                    max(X*fac/rhoSolid, scalar(0.0)),
                    scalar(1.0)
                )
            );

            trhok = tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "trhok",
                        T.mesh().time().timeName(),
                        T.mesh()
                    ),
                    (1.-volumeFraction) + volumeFraction * rhoSolid/rho
                )
            );
        }
        else
        {
            trhok = tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "trhok",
                        T.mesh().time().timeName(),
                        T.mesh()
                    ),
                    T.mesh(),
                    dimensionedScalar("trhok", dimless, 1.0)
                )
            );
        }
    }
    else
    {
        FatalErrorInFunction
            << "Invalid density calculation method: " << rhoType << "." << nl
            << "Valid methods are: Beta, BetaLinear, Thiesen (water)." << endl;
    }

    return trhok;
}
