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

// Creating e based thermo
autoPtr<twoPhaseMixtureEThermo> thermo;
thermo.set(new twoPhaseMixtureEThermo(U, phi));

// Create mixture and
Info<< "Creating temperaturePhaseChangeTwoPhaseMixture\n" << endl;
autoPtr<temperaturePhaseChangeTwoPhaseMixture> mixture =
    temperaturePhaseChangeTwoPhaseMixture::New(thermo(), mesh);


// Correct e from T and alpha
//thermo->correct();

volScalarField& alpha1(thermo->alpha1());
volScalarField& alpha2(thermo->alpha2());

const dimensionedScalar& rho1 = thermo->rho1();
const dimensionedScalar& rho2 = thermo->rho2();

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2,
    alpha1.boundaryField().types()
);
rho.oldTime();

// Construct interface from alpha1 distribution
interfaceProperties interface
(
    alpha1,
    U,
    thermo->transportPropertiesDict()
);

#include "cfdTools/general/include/readGravitationalAcceleration.H"
#include "cfdTools/general/include/readhRef.H"
#include "cfdTools/general/include/gh.H"

// Construct incompressible turbulence model
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, thermo())
);


volScalarField& p = thermo->p();

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
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
    p_rgh = p - rho*gh;
}

mesh.schemes().setFluxRequired(p_rgh.name());
mesh.schemes().setFluxRequired(alpha1.name());

// Turbulent Prandtl number
dimensionedScalar Prt("Prt", dimless, thermo->transportPropertiesDict());

volScalarField kappaEff
(
    IOobject
    (
        "kappaEff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo->kappa()
);

Info<< "Creating field pDivU\n" << endl;
volScalarField pDivU
(
    IOobject
    (
        "pDivU",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pDivU", p.dimensions()/dimTime, 0)
);

#include "cfdTools/general/include/createFvOptions.H"

bool reportMdot(mesh.solution().solverDict(alpha1.name()).lookupOrDefault<Switch>("reportMdot", false));

tmp<volScalarField> tmDotAlphaApprox;

if (reportMdot)
{
    Info<< "Creating field mDotAlphaApprox\n" << endl;
    tmDotAlphaApprox = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "mDotAlphaApprox",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar("mdot", dimDensity/dimTime, 0.0)
        )
    );
}
