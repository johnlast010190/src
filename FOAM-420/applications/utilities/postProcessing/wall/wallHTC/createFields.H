/*---------------------------------------------------------------------------*\
| Modified 2010-2016 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/

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



// Read radiative heat-flux if available
volScalarField Qr
(
    IOobject
    (
        "Qr",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Qr", dimMass/pow3(dimTime), 0.0)
);

// Read gravity field if present
autoPtr<uniformDimensionedVectorField> g;
IOobject gIOObj
(
    "g",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE
);
if (gIOObj.headerOk())
{
    Info<< "\nReading g" << endl;
    g.set(new uniformDimensionedVectorField(gIOObj));
}
