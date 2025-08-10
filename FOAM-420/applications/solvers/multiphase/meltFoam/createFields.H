    // Reading fields
    Info<< "Reading field T\n" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field alpha\n" << endl;
    volScalarField alpha
    (
        IOobject
        (
            "alpha",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    #include "cfdTools/incompressible/createPhi.H"

    // Reading transport properties
    Info<< "Reading thermophysical properties\n" << endl;
    #include "readTransportProperties.H"

    // Calculating fit-parameters for phase change function
    Info<< "Calculating phase change properties\n" << endl;

    dimensionedScalar Tmelt
    (
        "Tmelt",
        (Tl+Ts)/2.0
    );

    // Kinematic density for buoyancy force
    volScalarField rhok
    (
        IOobject
        (
            "rhok",
            runTime.timeName(),
            mesh
        ),
        1.0 - Foam::max(beta*(T - Tl), scalar(0.0))
    );

    // D'arcy-type source term field
    volScalarField DC
    (
        IOobject
        (
            "DC",
            runTime.timeName(),
            mesh
        ),
        DCl*Foam::pow(1.0-alpha,2)/(Foam::pow(alpha,3)+DCs)
    );

    // Thermal conductivity field
    volScalarField lambda
    (
        IOobject
        (
            "lambda",
            runTime.timeName(),
            mesh
        ),
        alpha*lambdaL+(1.0-alpha)*lambdaS
    );

    // Heat capacity field
    volScalarField cp
    (
        IOobject
        (
            "cp",
            runTime.timeName(),
            mesh
        ),
        alpha*cpL+(1.0-alpha)*cpS
    );

    // Kinematic viscosity field
    volScalarField nu
    (
        IOobject
        (
            "nu",
            runTime.timeName(),
            mesh
        ),
        alpha*nuL+(1.0-alpha)*nuS
    );

    Info<< "Calculating field g.h\n" << endl;
    volScalarField gh("gh", g & mesh.C());
    surfaceScalarField ghf("ghf", g & mesh.Cf());

    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh + rhok*gh
    );

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        mesh.solution().dict().subDict("PIMPLE"),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
