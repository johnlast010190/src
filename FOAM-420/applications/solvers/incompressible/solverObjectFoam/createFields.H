
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

singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

tmp<volScalarField> rho(turbulence->rho());

#include "cfdTools/general/include/createFvOptions.H"

IOdictionary fvOptionsDict
(
    IOobject
    (
        "fvOptions",
        mesh.time().system(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const word solutionMode = fvOptionsDict.lookupOrDefault<word>("solutionMode", Foam::word("pimple"));
