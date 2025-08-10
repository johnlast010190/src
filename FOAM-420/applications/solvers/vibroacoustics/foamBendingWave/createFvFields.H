Info<< "\nCreating Fv Fields" << endl;

volScalarField wV
(
    IOobject
    (
        "wV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("w",dimensionSet(0,1,0,0,0,0,0),0)
);

volScalarField loadV
(
    IOobject
    (
        "loadV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("load",dimensionSet(1,-1,-2,0,0,0,0),0)
);

volScalarField magSf
(
    IOobject
    (
        "magSf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("load",dimless,0)
);

volScalarField crit
(
    IOobject
    (
        "crit",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("load",dimless,0)
);

/*
volScalarField wLapV
(
    IOobject
    (
        "wLapV",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wLap",dimensionSet(0,-3,0,0,0,0,0),0)
);
*/
