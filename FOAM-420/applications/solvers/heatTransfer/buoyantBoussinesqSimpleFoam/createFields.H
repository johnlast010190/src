/*---------------------------------------------------------------------------*\
| Modified 2010-2017 Copyright (C) Esi Ltd                                  |
\*---------------------------------------------------------------------------*/
Switch cellBuoyancy = true;
if (mesh.solution().dict().subDict("SIMPLE").found("buoyancyFormulation"))
{
    word buoyancyFormulation
    (
        mesh.solution().dict().subDict("SIMPLE").lookup("buoyancyFormulation")
    );

    if (buoyancyFormulation == "face")
    {
        cellBuoyancy = false;
    }
    else if (buoyancyFormulation == "cell")
    {
        cellBuoyancy = true;
    }
    else
    {
        FatalError << "Unknown buoyancy formulation " << buoyancyFormulation
            << " specified. Valid types are 'cell' and 'face'" << nl
            << exit(FatalError);
    }
}

Info<< "Reading thermophysical properties\n" << endl;

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

#include "cfdTools/general/include/readGravitationalAcceleration.H"
#include "cfdTools/general/include/readhRef.H"
#include "cfdTools/general/include/gh.H"

#include "cfdTools/incompressible/createPhi.H"

#include "readTransportProperties.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
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
    updateRho(laminarTransport, T)()
);

// effective density for post-processing
volScalarField rhoEff
(
    IOobject
    (
        "rhoEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    turbulence->rho()() * rhok
);
rhoEff.write();


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
    simple.dict(),
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

mesh.schemes().setFluxRequired(p_rgh.name());

p.write();

#include "cfdTools/general/include/createFvOptions.H"
#include "radiation/createIncompressibleRadiationOption.H"