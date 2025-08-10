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

volScalarField pInt
(
    IOobject
    (
        "pInt",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
    //dimensionedScalar("pInt",dimensionSet(1,-1,-2,0,0,0,0),0)
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
        //IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("wLap",dimensionSet(0,-3,0,0,0,0,0),0)
);
*/
